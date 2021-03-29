
"""
This is an example that shows how to use thermoelectric.py.
In this script the ideal electronn filtering in Si-based composite
containing nanoscale SiC particles is studied.
Required data files are available in repo under Datafile directory
Cite: Mitigating the Effect of Nanoscale Porosity on Thermoelectric Power Factor of Si,
        Hosseini, S. Aria and Romano, Giuseppe and Greaney, P. Alex,
        ACS Applied Energy Materials,2021,
        https://doi.org/10.1021/acsaem.0c02640.
Author: S. Aria Hosseini
"""

# Required libs

import numpy as np
from numpy.linalg import norm

# Libs to generate figs

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import ScalarFormatter
import seaborn as sns

# Import thermoelectric.py

from thermoelectricProperties import thermoelectricProperties

Si = thermoelectricProperties(latticeParameter=5.401803661945516e-10, dopantElectricCharge=1,
                              electronEffectiveMass=1.08*thermoelectricProperties.me,
                              energyMin=0.0, energyMax=1, dielectric=11.7,
                              numKpoints=800, numBands=8, numQpoints=201, numEnergySampling=1000)

vfrac = 0.05                            # Nanoparticles volume fraction,
                                        # for materials containg nanoscale pores,
                                        # this will be the porosity

ml = 0.98*thermoelectricProperties.me   # Longitudinal conduction band effective mass
mt = 0.19*thermoelectricProperties.me   # Transverse conduction band effective mass

bulk_module = 98                        # Bulk module (GPA)
rho = 2329                              # Mass density (Kg/m3)
sp = np.sqrt(bulk_module/rho)           # Speed of sound


# Silicon has cubic unitcell, here is to define lattice vector and the reciprocal lattice vector

Lv = np.array([[1,1,0],[0,1,1],[1,0,1]])*Si.latticeParameter/2    # Lattice vector
a_rp = np.cross(Lv[1],Lv[2])/np.dot(Lv[0],np.cross(Lv[1],Lv[2]))  # Reciprocal lattice vector along a
b_rp = np.cross(Lv[2],Lv[0])/np.dot(Lv[1],np.cross(Lv[2],Lv[0]))  # Reciprocal lattice vector along b
c_rp = np.cross(Lv[0],Lv[1])/np.dot(Lv[2],np.cross(Lv[0],Lv[1]))  # Reciprocal lattice vector along c
RLv = np.array([a_rp, b_rp, c_rp])                                # Reciprocal lattice vector

e = Si.energyRange()    # Energy range, the default is from 0 to 1 eV which is reasonable for near equilibrium transport
number_of_points = 100

Tmp = np.expand_dims(500*np.ones(number_of_points), axis=0)     # Temperature
gap = Si.bandGap(Eg_o=1.17, Ao=4.73e-4, Bo=636, Temp=Tmp)       # Bandgap
alpha = np.array(0.5*np.tile([1],(1,len(gap[0]))))              # Nonparabolic term


n = np.linspace(19, 21, number_of_points, endpoint=True)        # Carrier concentration in log scale
cc = np.expand_dims(10**n, axis=0)*1e6                          # Carrier concentration in 1/m^3

# kpoints and band structure from EIGENVAL file of VASP DFT
kp, band = Si.electronBandStructure(path2eigenval='EIGENVAL', skipLines=6)


kp_rl = 2*np.pi*np.matmul(kp,RLv)                   # kpoints in reciprocal space
kp_mag = norm(kp_rl, axis=1)                        # Magnitude of kpoints
min_band = np.argmin(band[400:600, 4], axis=0)      # Index of band minimum edge
max_band = np.argmax(band[400:600, 4], axis=0)      # Index of band maximum edge
kp_vel = kp_mag[401 + max_band:401 + min_band]      # This is the kpoints for the desired band, BTE needs single band of conduction

# This is the desired energy band, Note that BTE needs single band of conduction
energy_vel = band[401 + max_band:401 + min_band, 4] - band[401 + min_band, 4]
enrg_sorted_idx = np.argsort(energy_vel, axis=0)

# Electron group velocity
gVel = Si.electronGroupVelocity(kp=kp_vel[enrg_sorted_idx],
                                energy_kp=energy_vel[enrg_sorted_idx], energyRange=e)

# Electron density of state from DOSCAR of VASP, the volume is printed in the header of DOSCAR
DoS = Si.electronDoS(path2DoS='DOSCAR', headerLines=6, unitcell_volume=2*19.70272e-30, numDoSpoints=2000, valleyPoint=1118, energyRange=e)

"Joyce Dixon approx. for the Fermi level"
JD_f, JD_n = Si.fermiLevel(carrierConcentration=cc, energyRange=e, DoS= DoS, Nc=None, Ao=5.3e21, Temp=Tmp)

# Self consistent method to compute Ef to circumvent the problem that DFT underestimates the band gaps
fermi, cc_sc = Si.fermiLevelSelfConsistent(carrierConcentration=cc, Temp=Tmp, energyRange=e, DoS=DoS, fermilevel=JD_f)

# Fermi distribution
dis, dfdE = Si.fermiDistribution(energyRange=e, Temp=Tmp, fermiLevel=fermi)

"""
The following line save the Fermi level.
Fermi integral of them is needed to compute Debye length is degenerate (highly doped) dielectrics.
See the manual for available codes to do the task
"""
np.savetxt("Ef", fermi/Tmp/thermoelectricProperties.kB)

# Screening length for degenerate (highly doped) dielectrics
LD_nondegenerate = np.sqrt(4*np.pi*Si.dielectric*thermoelectricProperties.e0*thermoelectricProperties.kB /
                           thermoelectricProperties.e2C*Tmp/cc_sc)  # screening length

# Effective density of states in the conduction band
m_CB = 0.23*thermoelectricProperties.me*(1+5*alpha*thermoelectricProperties.kB*Tmp)     # conduction band effective mass

# Conduction band effective mass as a function of temperature
Nc = 2*(m_CB*thermoelectricProperties.kB*Tmp/thermoelectricProperties.hBar**2 /
        2/np.pi/thermoelectricProperties.e2C)**(3/2)

# Fermi integral from third party lib., see the manual.
fermi_int = np.loadtxt("fermi_integral", delimiter=',')

# Screening length for degenerate (highly doped) dielectrics
LD = np.sqrt(1/(Nc/Si.dielectric/thermoelectricProperties.e0/thermoelectricProperties.kB /
                Tmp*thermoelectricProperties.e2C*(fermi_int[1]+15*alpha *
                                                  thermoelectricProperties.kB*Tmp/4*fermi_int[0])))

# Lifetime for electron-phonon scattering process following Ravich method
# tau_p_pb is for parabolic baand and tau_p_npb is for nonparabolic band structure
tau_p_pb, tau_p_npb = Si.tau_p(energyRange=e, alpha=alpha, Dv=2.94, DA=9.5, T=Tmp, vs=sp, D=DoS, rho=rho)

# Lifetime for electron-ion scattering for highly doped dielectrics
tau_ion = Si.tau_Strongly_Screened_Coulomb(D=DoS, LD=LD, N=cc_sc)

"Use Matthiessen's rule to compute total lifetime, 6 stands for number of degenerate valleys"
tau = Si.matthiessen(e, 6*tau_p_npb, 6*tau_ion)

# Define energy cutoff
uo = np.arange(0, 0.5, 0.02)
xv, yv = np.meshgrid(n, uo)

# Compute ideal filtering for the given energy cutoff, electrical conductivity and Seebeck
Coeff_f = Si.filteringEffect(U= uo, E = e, DoS=DoS, vg=gVel, Ef=fermi, dfdE=dfdE, Temp=Tmp, tau_b=tau)

PF_f = Coeff_f[0]*Coeff_f[1]**2  # power factor

print("done")


def generate_figs():
    sns.set()
    sns.set_context("paper", font_scale = 2, rc = {"lines.linewidth": 4} )
    sns.set_style("ticks", {"xtick.major.size": 2, "ytick.major.size": 2} )
    plt.figure(figsize = (6.5, 6.5))
    ax = plt.axes(projection = '3d')
    ax.set_axis_on()
    ax.grid(False)
    ax.view_init(15, -40)
    ax.xaxis.pane.set_edgecolor('black')
    ax.yaxis.pane.set_edgecolor('black')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.plot_surface(xv, yv, PF_f.reshape((len(uo), len(tau))) * 1e3, color = 'blue', edgecolor = 'black',
                    linewidth = 0.15, alpha = 0.9, rstride = 1, cstride = 1, shade = True, antialiased = False)

    ax.set_xlabel('Carrier concentration\n [log$_{10}$(cm$^{-3}$)]', fontsize = 20, labelpad = 25)
    ax.set_ylabel('U$_o$ (eV)', fontsize = 20, labelpad = 15)
    ax.zaxis.set_rotate_label( False )  # disable automatic rotation
    ax.set_zlabel('Power factor\n (mW/mK$^2$)', fontsize = 20, labelpad = 15, rotation = 90)
    ax.tick_params(axis = "y", labelsize = 20)
    ax.tick_params(axis = "x", labelsize = 20)
    ax.tick_params(axis = "z", labelsize = 20)
    plt.show()


generate_figs()
exit()
