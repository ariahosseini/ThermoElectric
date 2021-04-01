"""
This is an example that shows how to use thermoelectric.py.
In this script electrical properties of Si-based composite
containing nanoscale SiC particles is studied.


Required data files are available in repo under datafile directory

Cite: Mitigating the Effect of Nanoscale Porosity on Thermoelectric Power Factor of Si,
        Hosseini, S. Aria and Romano, Giuseppe and Greaney, P. Alex,
        ACS Applied Energy Materials,2021,
        https://doi.org/10.1021/acsaem.0c02640.

Author: S. Aria Hosseini
"""

# Required libs

import numpy as np
from numpy.linalg import norm
from scipy.interpolate import PchipInterpolator
from accum import accum

# Libs to generate figs
from util.generate_figs.py import 

# Import thermoelectric.py

from thermoelectricProperties import thermoelectricProperties

"""
Data file for pristine, 1%, and 5% nanoparticle's volume fraction.
Note that the process of dissolving P dopants is nonreversible,
so, the concentration is different when heating up (direction_up),
and while cooling down (direction_down)
"""

ExpData_SiCfra_0pct_direction_up = np.loadtxt('ExpData_SiCfrac-0pct_direction-up.txt', delimiter=None, skiprows=1)
ExpData_SiCfrac_1pct_direction_up = np.loadtxt('ExpData_SiCfrac-1pct_direction-up.txt', delimiter=None, skiprows=1)
ExpData_SiCfrac_5pct_direction_down = np.loadtxt('ExpData_SiCfrac-5pct_direction-down.txt', delimiter=None, skiprows=1)
ExpData_SiCfrac_5pct_direction_up = np.loadtxt('ExpData_SiCfrac-5pct_direction-up.txt', delimiter=None, skiprows=1)

# f variables are [T(K), Nc(1/cm^3)], where Nc is the carreir concentration
f0 = np.array([ExpData_SiCfra_0pct_direction_up[:,0], ExpData_SiCfra_0pct_direction_up[:,-2]*1e20])       # pristine Si
f1 = np.array([ExpData_SiCfrac_5pct_direction_up[:,0], ExpData_SiCfrac_5pct_direction_up[:,-2]*1e20])     # Si with 5% of SiC while heating up
f2 = np.array([ExpData_SiCfrac_5pct_direction_down[:,0], ExpData_SiCfrac_5pct_direction_down[:,-2]*1e20]) # Si with 5% of SiC while cooling down
f3 = np.array([ExpData_SiCfrac_1pct_direction_up[:,0], ExpData_SiCfrac_1pct_direction_up[:,-2]*1e20])     # Si with 1% of SiC while heating up

# Initiate the Obj, let's call it Si
Si = thermoelectricProperties ( latticeParameter = 5.401803661945516e-10, dopantElectricCharge = 1, \
                                electronEffectiveMass = 1.08 * thermoelectricProperties.me, energyMin = 0.0, \
                                energyMax = 1, dielectric = 11.7, numKpoints = 800, numBands = 8, numQpoints = 201, \
                                numEnergySampling = 5000 )

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

e = Si.energyRange()                                    # Energy range, the default is from 0 to 1 eV which is reseanable for near equilibrium transport
Tmp = Si.temp(TempMin=300, TempMax=1201, dT=50)         # Desired temperature range
h = Si.bandGap(Eg_o=1.17, Ao=4.73e-4, Bo=636, Temp=Tmp) # Band structe, see the manual for the description
alpha = np.array(0.5*np.tile([1],(1,len(h[0]))))        # Nonparabolic term shows the mixture of S and P orbitals,
                                                        # for Si it is 0.5, it is defined as a function of temperature to be general

"""
Carrier concentratin for pristine, 1% and 5% nanoparticle's valume fraction.
Note that the process of desolving P dopants is nonreversible,
so the concentration is different when heating up (direction_up),
and while cooling down (direction_down)
"""
cc_no_inc = Si.carrierConcentration(Nc=None, Nv=None,
                                    path2extrinsicCarrierConcentration='experimental-carrier-concentration-no-inc.txt', \
                                    bandGap=h, Ao=5.3e21, Bo=3.5e21, Temp=Tmp)                                  # pristine Si

cc = Si.carrierConcentration(Nc=None, Nv=None,
                             path2extrinsicCarrierConcentration='experimental-carrier-concentration-5pct-direction-up.txt', \
                             bandGap=h, Ao=5.3e21, Bo=3.5e21, Temp=Tmp)                                         # Si with 5% of SiC while heating up

cc_direction_down = Si.carrierConcentration(Nc=None, Nv=None,
                                            path2extrinsicCarrierConcentration='experimental-carrier-concentration-5pct-direction-down.txt', \
                                            bandGap=h, Ao=5.3e21, Bo=3.5e21, Temp=Tmp)                          # Si with 5% of SiC while cooling down

cc_1pct = Si.carrierConcentration(Nc=None, Nv=None,
                                  path2extrinsicCarrierConcentration='experimental-carrier-concentration-1pct.txt', \
                                  bandGap=h, Ao=5.3e21, Bo=3.5e21, Temp=Tmp)                                    # Si with 1% of SiC while heating up

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

# pristine Si
JD_f_no_inc, JD_n_no_inc = Si.fermiLevel(carrierConcentration=cc_no_inc,
                                         energyRange=e, DoS= DoS, Nc=None, Ao=5.3e21, Temp=Tmp)

# Si with 5% of SiC while heating up
JD_f, JD_n = Si.fermiLevel(carrierConcentration=cc,
                           energyRange=e, DoS= DoS, Nc=None, Ao=5.3e21, Temp=Tmp)

# Si with 5% of SiC while cooling down
JD_f_direction_down, JD_n_direction_down = Si.fermiLevel(carrierConcentration=cc_direction_down,
                                                         energyRange=e, DoS= DoS, Nc=None, Ao=5.3e21, Temp=Tmp)

# Si with 1% of SiC while heating up
JD_f_1pct, JD_n_1pct = Si.fermiLevel(carrierConcentration=cc_1pct,
                                     energyRange=e, DoS= DoS, Nc=None, Ao=5.3e21, Temp=Tmp)

# Self consistent method to compute Ef to circumvent the problem that DFT underestimates the band gaps
fermi_no_inc, cc_sc_no_inc = Si.fermiLevelSelfConsistent(carrierConcentration=cc_no_inc,
                                                         Temp=Tmp, energyRange=e, DoS=DoS,
                                                         fermilevel=JD_f_no_inc)                                # pristine Si

fermi, cc_sc = Si.fermiLevelSelfConsistent(carrierConcentration=cc,
                                           Temp=Tmp, energyRange=e,
                                           DoS=DoS, fermilevel=JD_f)                                            # Si with 5% of SiC while heating up

fermi_direction_down, cc_sc_direction_down = Si.fermiLevelSelfConsistent(carrierConcentration=cc_direction_down,
                                                                         Temp=Tmp, energyRange=e, DoS=DoS,
                                                                         fermilevel=JD_f_direction_down)        # Si with 5% of SiC while cooling down

fermi_1pct, cc_sc_1pct = Si.fermiLevelSelfConsistent(carrierConcentration=cc_1pct,
                                                     Temp=Tmp, energyRange=e, DoS=DoS,
                                                     fermilevel=JD_f_1pct)                                      # Si with 1% of SiC while heating up

# Fermi distribution
dis_no_inc, dfdE_no_inc = Si.fermiDistribution(energyRange=e, Temp=Tmp, fermiLevel=fermi_no_inc)                # pristine Si
dis, dfdE = Si.fermiDistribution(energyRange=e, Temp=Tmp, fermiLevel=fermi)                                     # Si with 5% of SiC while heaating up
dis_direction_down, dfdE_direction_down = Si.fermiDistribution(energyRange=e, Temp=Tmp,
                                                               fermiLevel=fermi_direction_down)                 # Si with 5% of SiC while cooling down
dis_1pct, dfdE_1pct = Si.fermiDistribution(energyRange=e, Temp=Tmp,
                                           fermiLevel=fermi_1pct)                                               # Si with 1% of SiC while heaating up

"""
The following lines save the Fermi level.
Fermi integral of them is needed to compute Debye length is degenerate (highly doped) dielectrics.
See the manual for available codes to do the task
"""
np.savetxt("Ef-no-inc",fermi_no_inc/Tmp/thermoelectricProperties.kB)
np.savetxt("Ef-inc",fermi/Tmp/thermoelectricProperties.kB)
np.savetxt("Ef-inc_direction_down",fermi_direction_down/Tmp/thermoelectricProperties.kB)
np.savetxt("Ef-inc_1pct",fermi_1pct/Tmp/thermoelectricProperties.kB)

# Debye length in nondegenerate dielectrics
LD_nondegenerate_no_inc = np.sqrt(4*np.pi*Si.dielectric*thermoelectricProperties.e0*thermoelectricProperties.kB/
                                  thermoelectricProperties.e2C*Tmp/cc_sc_no_inc)                        # pristine Si

LD_nondegenerate = np.sqrt(4*np.pi*Si.dielectric*thermoelectricProperties.e0*thermoelectricProperties.kB/
                           thermoelectricProperties.e2C*Tmp/cc_sc)                                      # Si with 5% of SiC while heating up

LD_nondegenerate_direction_down = np.sqrt(4*np.pi*Si.dielectric*thermoelectricProperties.e0*thermoelectricProperties.kB/
                                          thermoelectricProperties.e2C*Tmp/cc_sc_direction_down)        # Si with 5% of SiC while cooling down

LD_nondegenerate_1pct = np.sqrt(4*np.pi*Si.dielectric*thermoelectricProperties.e0*thermoelectricProperties.kB/
                                thermoelectricProperties.e2C*Tmp/cc_sc_1pct)                            # Si with 1% of SiC while heating up

# Conduction band effective mass as a function of temperature
m_CB_no_inc = 0.23*thermoelectricProperties.me*(1+5*alpha*thermoelectricProperties.kB*Tmp)              # pristine Si
m_CB_inc = 0.23*thermoelectricProperties.me*(1+5*alpha*thermoelectricProperties.kB*Tmp)                 # Si with 5% of SiC while heating up
m_CB_inc_direction_down = 0.23*thermoelectricProperties.me*(1+5*alpha*thermoelectricProperties.kB*Tmp)  # Si with 5% of SiC while cooling down
m_CB_inc_1pct = 0.23*thermoelectricProperties.me*(1+5*alpha*thermoelectricProperties.kB*Tmp)            # Si with 1% of SiC while heating up

# Effective density of states in the conduction band
Nc_no_inc = 2*(m_CB_no_inc*thermoelectricProperties.kB*Tmp/thermoelectricProperties.hBar**2/2/np.pi/thermoelectricProperties.e2C)**(3/2)
Nc_inc = 2*(m_CB_inc*thermoelectricProperties.kB*Tmp/thermoelectricProperties.hBar**2/2/np.pi/thermoelectricProperties.e2C)**(3/2)
Nc_inc_direction_down = 2*(m_CB_inc_direction_down*thermoelectricProperties.kB*Tmp/thermoelectricProperties.hBar**2/2/np.pi/thermoelectricProperties.e2C)**(3/2)
Nc_inc_1pct = 2*(m_CB_inc_1pct*thermoelectricProperties.kB*Tmp/thermoelectricProperties.hBar**2/2/np.pi/thermoelectricProperties.e2C)**(3/2)

# Fermi integral from the third party package. See manual for the detail explanation
fermi_int = np.loadtxt("f_inc", delimiter=',')
fermi_no_inc_int = np.loadtxt("f_no_inc", delimiter=',')
fermi_int_direction_down = np.loadtxt("f_inc_direction_down", delimiter=',')
fermi_int_1pct = np.loadtxt("f_inc_1pct", delimiter=',')

# Screening length for degenerate (highly doped) dielectrics

LD = np.sqrt(1/(Nc_inc/Si.dielectric/thermoelectricProperties.e0/thermoelectricProperties.kB/
                Tmp*thermoelectricProperties.e2C*(fermi_int[1]+15*alpha*thermoelectricProperties.kB*Tmp/
                                                  4*fermi_int[0]))) # Si with 5% of SiC while cooling down

LD_no_inc = np.sqrt(1/(Nc_no_inc/Si.dielectric/thermoelectricProperties.e0/thermoelectricProperties.kB/
                       Tmp*thermoelectricProperties.e2C*(fermi_no_inc_int[1]+15*alpha*thermoelectricProperties.kB*Tmp
                                                         /4*fermi_no_inc_int[0]))) # pristine Si

LD_direction_down  = np.sqrt(1/(Nc_inc_direction_down /Si.dielectric/thermoelectricProperties.e0/thermoelectricProperties.kB
                                /Tmp*thermoelectricProperties.e2C*(fermi_int_direction_down[1]+15*alpha*thermoelectricProperties.kB*Tmp
                                                                   /4*fermi_int_direction_down[0]))) # Si with 5% of SiC while cooling down

LD_int_1pct = np.sqrt(1/(Nc_inc_1pct/Si.dielectric/thermoelectricProperties.e0/thermoelectricProperties.kB
                         /Tmp*thermoelectricProperties.e2C*(fermi_int_1pct[1]+15*alpha*thermoelectricProperties.kB*Tmp
                                                            /4*fermi_int_1pct[0]))) # Si with 1% of SiC while heaating up

# Lifetime for electron-phonon scattering process following Ravich method
# tau_p_pb is for parabolic baand and tau_p_npb is for nonparabolic band structure
tau_p_pb, tau_p_npb = Si.tau_p(energyRange=e, alpha=alpha, Dv=2.94, DA=9.5, T=Tmp, vs=sp, D=DoS, rho=rho)

# Lifetime for electron-ion scattering for highly doped dielectrics
tau_ion_no_inc = Si.tau_Strongly_Screened_Coulomb(D=DoS, LD=LD_no_inc, N=cc_sc_no_inc)
tau_ion = Si.tau_Strongly_Screened_Coulomb(D=DoS, LD=LD, N=cc_sc)
tau_ion_direction_down = Si.tau_Strongly_Screened_Coulomb(D=DoS, LD=LD_direction_down, N=cc_sc_direction_down)
tau_ion_1pct = Si.tau_Strongly_Screened_Coulomb(D=DoS, LD=LD_int_1pct, N=cc_sc_1pct)

# Electron lifetime from spherical nanoparticles and the energy levels
lifetime_nanoparticle = np.loadtxt('lifetime_np', delimiter=None, skiprows=0)
energy_nanoparticle = np.loadtxt('energy_np', delimiter=None, skiprows=0)
# Find isoenergy states
E_energy_nanoparticle, indices_energy_nanoparticle, return_indices_energy_nanoparticle = np.unique(energy_nanoparticle, return_index=True, return_inverse=True)
# Average isoenergy lifetime (smoothing the plot)
lf_nanoparticle = accum(return_indices_energy_nanoparticle, lifetime_nanoparticle[1], func=np.mean, dtype=float)
nanoparticle_Spline = PchipInterpolator(E_energy_nanoparticle[1::], lf_nanoparticle[1::])

# Electron lifetime from grains and the energy levels
_lifetime_gb = np.loadtxt('tau_g_vs_Nc', delimiter=None, skiprows=0)
lifetime_gb = _lifetime_gb[0]
energy_gb = np.loadtxt('e_g', delimiter=None, skiprows=0)
E_energy_gb, indices_energy_gb, return_indices_energy_gb = np.unique(energy_gb, return_index=True, return_inverse=True)
# Average isoenergy lifetime (smoothing the plot)
lf_gb = accum(return_indices_energy_gb, lifetime_gb, func=np.mean, dtype=float)
gb_Spline = PchipInterpolator(E_energy_gb[1::], lf_gb[1::])
_tau_gb = gb_Spline(e)

Ng = 1e25                       # Concentration
tau_g = _tau_gb * (cc / Ng).T   # Electron grains lifetime

tau_np= nanoparticle_Spline(e) # Find lifetime from nanoparticels for the desired energy range

"Use Matthiessen's rule to compute total lifetime, 6 stands for number of degenerate valleys"
tau_no_inc = Si.matthiessen(e, 6*tau_p_npb, 6*tau_ion_no_inc,tau_g)
tau_no_np = Si.matthiessen(e, 6*tau_p_npb, 6*tau_ion,tau_g)
tau = Si.matthiessen(e, 6*tau_p_npb, 6*tau_ion, tau_np,tau_g)
tau_direction_down = Si.matthiessen(e, 6*tau_p_npb, 6*tau_ion_direction_down, tau_np,tau_g)
tau_no_np_direction_down = Si.matthiessen(e, 6*tau_p_npb, 6*tau_ion_direction_down,tau_g)
tau_no_np_1pct = Si.matthiessen(e, 6*tau_p_npb, 6*tau_ion_1pct,tau_g)
tau_1pct = Si.matthiessen(e, 6*tau_p_npb, 6*tau_ion_1pct, 5*tau_np)  # 5 counts for 1% porosity instead of 5%

"Transport coefficients for pristine, 5% heat up, 5% cooling down aand 1% heating up cases and the effect of carrier concentration"

# pristine Si
Coeff_no_inc = Si.electricalProperties(E=e, DoS=DoS, vg=gVel, Ef=fermi_no_inc, dfdE=dfdE_no_inc, Temp=Tmp, tau=tau_no_inc)

# Si with 5% of SiC while heating up, no particle scattering included, importance of change in dopants
Coeff_no_np = Si.electricalProperties(E=e, DoS=DoS, vg=gVel, Ef=fermi, dfdE=dfdE, Temp=Tmp, tau=tau_no_np)

# Si with 5% of SiC while heating up, particle scattering included
Coeff = Si.electricalProperties(E=e, DoS=DoS, vg=gVel, Ef=fermi, dfdE=dfdE, Temp=Tmp, tau=tau)

# Si with 5% of SiC while cooling down, particle scattering included
Coeff_direction_down = Si.electricalProperties(E=e, DoS=DoS, vg=gVel, Ef=fermi_direction_down, 
                                               dfdE=dfdE_direction_down, Temp=Tmp, tau=tau_direction_down)

# Si with 5% of SiC while cooling down, no particle scattering included, importance of change in dopants
Coeff_direction_down_no_np = Si.electricalProperties(E=e, DoS=DoS, vg=gVel, Ef=fermi_direction_down, 
                                                     dfdE=dfdE_direction_down, Temp=Tmp, tau=tau_no_np_direction_down)

# Si with 1% of SiC while heating up, no particle scattering included, importance of change in dopants
Coeff_no_np_1pct = Si.electricalProperties(E=e, DoS=DoS, vg=gVel, Ef=fermi_1pct, dfdE=dfdE_1pct, Temp=Tmp, tau=tau_no_np_1pct)

# Si with 1% of SiC while heating up, paarticle scattering included
Coeff_1pct = Si.electricalProperties(E=e, DoS=DoS, vg=gVel, Ef=fermi_1pct, dfdE=dfdE_1pct, Temp=Tmp, tau=tau_1pct)

print("done")


# Generate figures

generate_figs()
