"""
Thi is an example shows how to use thermoelectric.py. 
In this script electrical properties of Si based composite 
containing nanoscale SiC particles is studied.


Required data files are availe in repo under datafile directory


Author: S. Aria Hosseini
"""

# Required libs

import numpy as np
from numpy.linalg import norm
from scipy.interpolate import PchipInterpolator
from accum import accum

# Libs to generate figs

import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import matplotlib.ticker
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits import mplot3d
from matplotlib.colors import LightSource
import seaborn as sns

# Import thermoelectric.py

from thermoelectricProperties import thermoelectricProperties

"""
Data file for pristine, 1% and 5% nanoparticle's valume fraction.
Note that the process of desolving P dopants is nonreversible,
so the concentration is different when heating up (direction_up),
and while cooling down (direction_down)
"""

ExpData_SiCfra_0pct_direction_up = np.loadtxt('ExpData_SiCfrac-0pct_direction-up.txt', /
                                              delimiter=None, skiprows=1)
ExpData_SiCfrac_1pct_direction_up = np.loadtxt('ExpData_SiCfrac-1pct_direction-up.txt', /
                                               delimiter=None, skiprows=1)
ExpData_SiCfrac_5pct_direction_down = np.loadtxt('ExpData_SiCfrac-5pct_direction-down.txt', /
                                                 delimiter=None, skiprows=1)
ExpData_SiCfrac_5pct_direction_up = np.loadtxt('ExpData_SiCfrac-5pct_direction-up.txt', /
                                               delimiter=None, skiprows=1)

# f variables are [T(K), Nc(1/cm^3)], where Nc is the carreir concentration
f0 = np.array([ExpData_SiCfra_0pct_direction_up[:,0], ExpData_SiCfra_0pct_direction_up[:,-2]*1e20])       # pristine Si
f1 = np.array([ExpData_SiCfrac_5pct_direction_up[:,0], ExpData_SiCfrac_5pct_direction_up[:,-2]*1e20])     # Si with 5% of SiC while heaating up
f2 = np.array([ExpData_SiCfrac_5pct_direction_down[:,0], ExpData_SiCfrac_5pct_direction_down[:,-2]*1e20]) # Si with 5% of SiC while cooling down
f3 = np.array([ExpData_SiCfrac_1pct_direction_up[:,0], ExpData_SiCfrac_1pct_direction_up[:,-2]*1e20])     # Si with 1% of SiC while heaating up

# Initiate the Obj, let's call it Si 
Si = thermoelectricProperties(latticeParameter=5.401803661945516e-10, dopantElectricCharge=1, \
                              electronEffectiveMass=1.08*thermoelectricProperties.me, energyMin=0.0, \
                              energyMax=1, dielectric=11.7, numKpoints=800, numBands=8, numQpoints=201, \
                              numEnergySampling=5000)

vfrac = 0.05                            # Nanoparticals volume fraction, 
                                        # for materials containg nanoscale pores, 
                                        # this will be the porosity

ml = 0.98*thermoelectricProperties.me   # Longitudinal conduction band effective mass
mt = 0.19*thermoelectricProperties.me   # Transverse conduction bandeffective mass

bulk_module = 98                        # Bulk module (GPA)
rho = 2329                              # Mass density (Kg/m3)
sp = np.sqrt(bulk_module/rho)           # Speed of sound


# Silicon has cubic unitcell, here is to define lattice vector and the reciprocal lattice vector

Lv = np.array([[1,1,0],[0,1,1],[1,0,1]])*Si.latticeParameter/2    # Lattice vector
a_rp = np.cross(Lv[1],Lv[2])/np.dot(Lv[0],np.cross(Lv[1],Lv[2]))  # Reciprocal lattice vector alonng a
b_rp = np.cross(Lv[2],Lv[0])/np.dot(Lv[1],np.cross(Lv[2],Lv[0]))  # Reciprocal lattice vector alonng b
c_rp = np.cross(Lv[0],Lv[1])/np.dot(Lv[2],np.cross(Lv[0],Lv[1]))  # Reciprocal lattice vector alonng c
RLv = np.array([a_rp, b_rp, c_rp])                                # Reciprocal lattice vector

e = Si.energyRange()                                    # Energy range, the fefaul is from 0 to 1 eV which is reseanable for near equilibrium transport
Tmp = Si.temp(TempMin=300, TempMax=1201, dT=50)         # Desired temperature range
h = Si.bandGap(Eg_o=1.17, Ao=4.73e-4, Bo=636, Temp=Tmp) # Band structe, see the manual for the discription
alpha = np.array(0.5*np.tile([1],(1,len(h[0]))))        # Nonparabolic term shows the mixture of S and P orbitals, 
                                                        # for Si it is 0.5, it is defined as a function of temperature to be general
  
"""
Carrier concentratin for pristine, 1% and 5% nanoparticle's valume fraction.
Note that the process of desolving P dopants is nonreversible,
so the concentration is different when heating up (direction_up),
and while cooling down (direction_down)
"""
cc_no_inc = Si.carrierConcentration(Nc=None, Nv=None, \
                                    path2extrinsicCarrierConcentration='experimental-carrier-concentration-no-inc.txt', \
                                    bandGap=h, Ao=5.3e21, Bo=3.5e21, Temp=Tmp)

cc = Si.carrierConcentration(Nc=None, Nv=None, \
                             path2extrinsicCarrierConcentration='experimental-carrier-concentration-5pct-direction-up.txt', \
                             bandGap=h, Ao=5.3e21, Bo=3.5e21, Temp=Tmp)

cc_direction_down = Si.carrierConcentration(Nc=None, Nv=None, \
                                            path2extrinsicCarrierConcentration='experimental-carrier-concentration-5pct-direction-down.txt', \
                                            bandGap=h, Ao=5.3e21, Bo=3.5e21, Temp=Tmp)

cc_1pct = Si.carrierConcentration(Nc=None, Nv=None, \
                                  path2extrinsicCarrierConcentration='experimental-carrier-concentration-1pct.txt', \
                                  bandGap=h, Ao=5.3e21, Bo=3.5e21, Temp=Tmp)

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
gVel = Si.electronGroupVelocity(kp=kp_vel[enrg_sorted_idx], \
                                energy_kp=energy_vel[enrg_sorted_idx], energyRange=e)

# Electron density of state from DOSCAR of VASP, the volume is printted in the header of DOSCAR
DoS = Si.electronDoS(path2DoS='DOSCAR', headerLines=6, unitcell_volume=2*19.70272e-30, numDoSpoints=2000, valleyPoint=1118, energyRange=e)

"Joyce Dixon approx for the Fermi level" 

# pristine Si
JD_f_no_inc, JD_n_no_inc = Si.fermiLevel(carrierConcentration=cc_no_inc, \
                                         energyRange=e, DoS= DoS, Nc=None, Ao=5.3e21, Temp=Tmp)

# Si with 5% of SiC while heaating up
JD_f, JD_n = Si.fermiLevel(carrierConcentration=cc, \
                           energyRange=e, DoS= DoS, Nc=None, Ao=5.3e21, Temp=Tmp)

# Si with 5% of SiC while cooling down
JD_f_direction_down, JD_n_direction_down = Si.fermiLevel(carrierConcentration=cc_direction_down, \
                                                         energyRange=e, DoS= DoS, Nc=None, Ao=5.3e21, Temp=Tmp)

# Si with 1% of SiC while heaating up
JD_f_1pct, JD_n_1pct = Si.fermiLevel(carrierConcentration=cc_1pct, \
                                     energyRange=e, DoS= DoS, Nc=None, Ao=5.3e21, Temp=Tmp)

# Self consistant method to compute Ef to circumvent the problem that DFT underestimates the band gaps
fermi_no_inc, cc_sc_no_inc = Si.fermiLevelSelfConsistent(carrierConcentration=cc_no_inc, \
                                                         Temp=Tmp, energyRange=e, DoS=DoS, \
                                                         fermilevel=JD_f_no_inc) # pristine Si

fermi, cc_sc = Si.fermiLevelSelfConsistent(carrierConcentration=cc, \
                                           Temp=Tmp, energyRange=e, \
                                           DoS=DoS, fermilevel=JD_f) # Si with 5% of SiC while heaating up

fermi_direction_down, cc_sc_direction_down = Si.fermiLevelSelfConsistent(carrierConcentration=cc_direction_down, \
                                                                         Temp=Tmp, energyRange=e, DoS=DoS, \
                                                                         fermilevel=JD_f_direction_down) # Si with 5% of SiC while cooling down

fermi_1pct, cc_sc_1pct = Si.fermiLevelSelfConsistent(carrierConcentration=cc_1pct, \
                                                     Temp=Tmp, energyRange=e, DoS=DoS, \
                                                     fermilevel=JD_f_1pct) # Si with 1% of SiC while heaating up

# Fermi distribution
dis_no_inc, dfdE_no_inc = Si.fermiDistribution(energyRange=e, Temp=Tmp, fermiLevel=fermi_no_inc) # pristine Si
dis, dfdE = Si.fermiDistribution(energyRange=e, Temp=Tmp, fermiLevel=fermi) # Si with 5% of SiC while heaating up
dis_direction_down, dfdE_direction_down = Si.fermiDistribution(energyRange=e, Temp=Tmp, \
                                                               fermiLevel=fermi_direction_down) # Si with 5% of SiC while cooling down
dis_1pct, dfdE_1pct = Si.fermiDistribution(energyRange=e, Temp=Tmp, \
                                           fermiLevel=fermi_1pct) # Si with 1% of SiC while heaating up

"""
The following lines save the Fermi level. 
Femi intergral of them is needed to compute Debye length is degenerate (highly doped) dielectrics. 
See the manual for availble codes to do the task
"""
np.savetxt("Ef-no-inc",fermi_no_inc/Tmp/thermoelectricProperties.kB)
np.savetxt("Ef-inc",fermi/Tmp/thermoelectricProperties.kB)
np.savetxt("Ef-inc_direction_down",fermi_direction_down/Tmp/thermoelectricProperties.kB)
np.savetxt("Ef-inc_1pct",fermi_1pct/Tmp/thermoelectricProperties.kB)

# Debye lenght in nondegenerate dielectrics
LD_nondegenerate_no_inc = np.sqrt(4*np.pi*Si.dielectric*thermoelectricProperties.e0 /
                                  *thermoelectricProperties.kB/thermoelectricProperties.e2C* /
                                  Tmp/cc_sc_no_inc)  # pristine Si

LD_nondegenerate = np.sqrt(4*np.pi*Si.dielectric*thermoelectricProperties.e0* /
                           thermoelectricProperties.kB/thermoelectricProperties.e2C* /
                           Tmp/cc_sc)  # Si with 5% of SiC while heaating up

LD_nondegenerate_direction_down = np.sqrt(4*np.pi*Si.dielectric*thermoelectricProperties.e0* /
                                          thermoelectricProperties.kB/thermoelectricProperties.e2C* /
                                          Tmp/cc_sc_direction_down)  # Si with 5% of SiC while cooling down

LD_nondegenerate_1pct = np.sqrt(4*np.pi*Si.dielectric*thermoelectricProperties.e0* /
                                thermoelectricProperties.kB/thermoelectricProperties.e2C* /
                                Tmp/cc_sc_1pct) # Si with 1% of SiC while heaating up

# Conduction band effective mass as a function of temperature
m_CB_no_inc = 0.23*thermoelectricProperties.me*(1+5*alpha*thermoelectricProperties.kB*Tmp)              # pristine Si
m_CB_inc = 0.23*thermoelectricProperties.me*(1+5*alpha*thermoelectricProperties.kB*Tmp)                 # Si with 5% of SiC while cooling down
m_CB_inc_direction_down = 0.23*thermoelectricProperties.me*(1+5*alpha*thermoelectricProperties.kB*Tmp)  # Si with 5% of SiC while cooling down
m_CB_inc_1pct = 0.23*thermoelectricProperties.me*(1+5*alpha*thermoelectricProperties.kB*Tmp)            # Si with 1% of SiC while heaating up

# Effective density of states in the conduction band
Nc_no_inc = 2*(m_CB_no_inc*thermoelectricProperties.kB*Tmp /
               /thermoelectricProperties.hBar**2/2/np.pi/thermoelectricProperties.e2C)**(3/2)

Nc_inc = 2*(m_CB_inc*thermoelectricProperties.kB*Tmp /
            /thermoelectricProperties.hBar**2/2/np.pi/thermoelectricProperties.e2C)**(3/2)

Nc_inc_direction_down = 2*(m_CB_inc_direction_down*thermoelectricProperties.kB*Tmp /
                           /thermoelectricProperties.hBar**2/2/np.pi/thermoelectricProperties.e2C)**(3/2)

Nc_inc_1pct = 2*(m_CB_inc_1pct*thermoelectricProperties.kB*Tmp /
                 /thermoelectricProperties.hBar**2/2/np.pi/thermoelectricProperties.e2C)**(3/2)

# Fermi integral from the thir party package. See manual for the detail explanation
fermi_int = np.loadtxt("f_inc", delimiter=',')
fermi_no_inc_int = np.loadtxt("f_no_inc", delimiter=',')
fermi_int_direction_down = np.loadtxt("f_inc_direction_down", delimiter=',')
fermi_int_1pct = np.loadtxt("f_inc_1pct", delimiter=',')

# Screening length for degenerate (highly doped) dielectrics

LD = np.sqrt(1/(Nc_inc/Si.dielectric/thermoelectricProperties.e0 /
                /thermoelectricProperties.kB/g*thermoelectricProperties.e2C* /
                (fermi_int[1]+15*alpha*thermoelectricProperties.kB*g/4*fermi_int[0]))) # Si with 5% of SiC while cooling down

LD_no_inc = np.sqrt(1/(Nc_no_inc/Si.dielectric/thermoelectricProperties.e0 /
                       /thermoelectricProperties.kB/g*thermoelectricProperties.e2C* /
                       (fermi_no_inc_int[1]+15*alpha*thermoelectricProperties.kB*g/4*fermi_no_inc_int[0]))) # pristine Si 

LD_direction_down  = np.sqrt(1/(Nc_inc_direction_down /Si.dielectric/thermoelectricProperties.e0 /
                                /thermoelectricProperties.kB/g*thermoelectricProperties.e2C* /
                                (fermi_int_direction_down[1]+15*alpha*thermoelectricProperties.kB*g/ /
                                 4*fermi_int_direction_down[0]))) # Si with 5% of SiC while cooling down

LD_int_1pct = np.sqrt(1/(Nc_inc_1pct/Si.dielectric/thermoelectricProperties.e0 /
                         /thermoelectricProperties.kB/g*thermoelectricProperties.e2C* /
                         (fermi_int_1pct[1]+15*alpha*thermoelectricProperties.kB*g/4*fermi_int_1pct[0]))) # Si with 1% of SiC while heaating up

tau_p_pb, tau_p_npb = Si.tau_p(energyRange=e, alpha=alpha, Dv=2.94, DA=9.5, T=g, vs=sp, D=DoS, rho=rho)
tau_ion_no_inc = Si.tau_Strongly_Screened_Coulomb(D=DoS, LD=LD_no_inc, N=cc_sc_no_inc)
tau_ion = Si.tau_Strongly_Screened_Coulomb(D=DoS, LD=LD, N=cc_sc)
tau_ion_direction_down = Si.tau_Strongly_Screened_Coulomb(D=DoS, LD=LD_direction_down, N=cc_sc_direction_down)
tau_ion_1pct = Si.tau_Strongly_Screened_Coulomb(D=DoS, LD=LD_int_1pct, N=cc_sc_1pct)


lifetime_nanoparticle = np.loadtxt('lifetime_np', delimiter=None, skiprows=0)
energy_nanoparticle = np.loadtxt('energy_np', delimiter=None, skiprows=0)
E_energy_nanoparticle, indices_energy_nanoparticle, return_indices_energy_nanoparticle = np.unique(energy_nanoparticle, return_index=True, return_inverse=True)
lf_nanoparticle = accum(return_indices_energy_nanoparticle, lifetime_nanoparticle[1], func=np.mean, dtype=float)
nanoparticle_Spline = PchipInterpolator(E_energy_nanoparticle[1::], lf_nanoparticle[1::])


_lifetime_gb = np.loadtxt('tau_g_vs_Nc', delimiter=None, skiprows=0)
lifetime_gb = _lifetime_gb[0]
energy_gb = np.loadtxt('e_g', delimiter=None, skiprows=0)
E_energy_gb, indices_energy_gb, return_indices_energy_gb = np.unique(energy_gb, return_index=True, return_inverse=True)
lf_gb = accum(return_indices_energy_gb, lifetime_gb, func=np.mean, dtype=float)
gb_Spline = PchipInterpolator(E_energy_gb[1::], lf_gb[1::])
_tau_gb = gb_Spline(e)

Ng = 1e25
tau_g = _tau_gb * (cc / Ng).T



tau_np= nanoparticle_Spline(e)
tau_no_inc = Si.matthiessen(e, 6*tau_p_npb, 6*tau_ion_no_inc,tau_g)
tau_no_np = Si.matthiessen(e, 6*tau_p_npb, 6*tau_ion,tau_g)
tau = Si.matthiessen(e, 6*tau_p_npb, 6*tau_ion, tau_np,tau_g)
tau_direction_down = Si.matthiessen(e, 6*tau_p_npb, 6*tau_ion_direction_down, tau_np,tau_g)
tau_no_np_direction_down = Si.matthiessen(e, 6*tau_p_npb, 6*tau_ion_direction_down,tau_g)
tau_no_np_1pct = Si.matthiessen(e, 6*tau_p_npb, 6*tau_ion_1pct,tau_g)
tau_1pct = Si.matthiessen(e, 6*tau_p_npb, 6*tau_ion_1pct, 5*tau_np)  # 5 counts for 1% porosity instead of 5%


# vg_analetical = Si.analyticalGroupVelocity(energyRange = e, nk = [40,38,38], m = [ml, mt, mt], valley = [0.85,0,0], dk_len = 0.15, alpha = alpha, temperature =g)

Coeff_no_inc = Si.electricalProperties(E=e, DoS=DoS, vg=gVel, Ef=fermi_no_inc, dfdE=dfdE_no_inc, Temp=g, tau=tau_no_inc)
Coeff_no_np = Si.electricalProperties(E=e, DoS=DoS, vg=gVel, Ef=fermi, dfdE=dfdE, Temp=g, tau=tau_no_np)
Coeff = Si.electricalProperties(E=e, DoS=DoS, vg=gVel, Ef=fermi, dfdE=dfdE, Temp=g, tau=tau)
Coeff_direction_down = Si.electricalProperties(E=e, DoS=DoS, vg=gVel, Ef=fermi_direction_down, dfdE=dfdE_direction_down, Temp=g, tau=tau_direction_down)
Coeff_direction_down_no_np = Si.electricalProperties(E=e, DoS=DoS, vg=gVel, Ef=fermi_direction_down, dfdE=dfdE_direction_down, Temp=g, tau=tau_no_np_direction_down)
Coeff_no_np_1pct  = Si.electricalProperties(E=e, DoS=DoS, vg=gVel, Ef=fermi_1pct, dfdE=dfdE_1pct, Temp=g, tau=tau_no_np_1pct)
Coeff_1pct  = Si.electricalProperties(E=e, DoS=DoS, vg=gVel, Ef=fermi_1pct, dfdE=dfdE_1pct, Temp=g, tau=tau_1pct)

print("done")

sns.set()
sns.set_context("paper", font_scale=2, rc={"lines.linewidth": 4})
sns.set_style("ticks", {"xtick.major.size": 2, "ytick.major.size": 2})

fig_0 = plt.figure(figsize=(6.5,4.5))
ax_0 = fig_0.add_subplot(111)
ax_0.plot(g[0],h[0], 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_0.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax_0.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_0.tick_params(axis="x", labelsize=16)
ax_0.set_ylabel('Band gap (eV)', fontsize=16, labelpad=10)
ax_0.tick_params(axis="y", labelsize=16)
fig_0.tight_layout()

# fig_1 = plt.figure(figsize=(6.5,4.5))
# ax_1 = fig_1.add_subplot(111)
# ax_1.plot(g[0],alpha[0], 'o', linestyle='-', color='maroon',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='maroon',
#           markeredgewidth=1)
# ax_1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
# ax_1.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
# ax_1.tick_params(axis="x", labelsize=16)
# ax_1.set_ylabel('Nonparabolic term (eV$^{-1}$)', fontsize=16, labelpad=10)
# ax_1.tick_params(axis="y", labelsize=16)
# fig_1.tight_layout()

# fig_2 = plt.figure(figsize=(6.5,4.5))
# ax_2 = fig_2.add_subplot(111)
# ax_2.plot(e[0],DoS_no_inc[0], 'None', linestyle='-', color='maroon',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='maroon',
#           markeredgewidth=1)
# ax_2.plot(e[0],DoS[0], 'None', linestyle='-', color='steelblue',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='steelblue',
#           markeredgewidth=1)
# # ax_2.plot(e[0],dos_nonparabolic[0], 'None', linestyle='-', color='steelblue',
# #           markersize=6, linewidth=1.5,
# #           markerfacecolor='white',
# #           markeredgecolor='steelblue',
# #           markeredgewidth=1)
# # ax_2.plot(e[0],dos_parabolic[0], 'None', linestyle='-', color='olive',
# #           markersize=6, linewidth=1.5,
# #           markerfacecolor='white',
# #           markeredgecolor='olive',
# #           markeredgewidth=1)
# ax_2.yaxis.set_major_formatter(ScalarFormatter())
# ax_2.set_xlabel('Energy (eV)', fontsize=16, labelpad=10)
# ax_2.tick_params(axis="x", labelsize=16)
# ax_2.set_ylabel('Density of state (#/eV/m$^3$)', fontsize=16, labelpad=10)
# ax_2.tick_params(axis="y", labelsize=16)
# ax_2.ticklabel_format(axis="y", style="sci", scilimits=None)
# fig_2.tight_layout()

# fig_3 = plt.figure(figsize=(6.5,4.5))
# ax_3 = fig_3.add_subplot(111)
# ax_3.plot(e[0],dos_nonparabolic[0], 'None', linestyle='-', color='maroon',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='maroon',
#           markeredgewidth=1)
# ax_3.plot(e[0],dos_nonparabolic[-1], 'None', linestyle='-', color='steelblue',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='steelblue',
#           markeredgewidth=1)
# ax_3.yaxis.set_major_formatter(ScalarFormatter())
# ax_3.set_xlabel('Energy (eV)', fontsize=16, labelpad=10)
# ax_3.tick_params(axis="x", labelsize=16)
# ax_3.set_ylabel('Density of state (#/eV/m$^3$)', fontsize=16, labelpad=10)
# ax_3.tick_params(axis="y", labelsize=16)
# ax_3.ticklabel_format(axis="y", style="sci", scilimits=None)
# fig_3.tight_layout()

# fig_4 = plt.figure(figsize=(6.5,4.5))
# ax_4 = fig_4.add_subplot(111)
# ax_4.plot(e[0],-1*gVel[0]*1e-5, 'None', linestyle='-', color='maroon',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='maroon',
#           markeredgewidth=1)
# ax_4.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
# ax_4.set_xlabel('Energy (eV)', fontsize=16, labelpad=10)
# ax_4.tick_params(axis="x", labelsize=16)
# ax_4.set_ylabel('Group velocity (x10$^5$ m/s)', fontsize=16, labelpad=10)
# fig_4.tight_layout()

fig_5 = plt.figure(figsize=(6.5,4.5))
ax_5 = fig_5.add_subplot(111)
ax_5.plot(band[1::,], 'None', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_5.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax_5.set_xlabel('BZ tour', fontsize=16, labelpad=10)
ax_5.tick_params(axis="x", labelsize=16)
ax_5.set_ylabel('Energy (eV)', fontsize=16)
ax_5.tick_params(axis="y", labelsize=16)
ax_5.set_xticks([0,199,399,599,799])
ax_5.set_xticklabels(["W", "L","$\Gamma$", "X", "W"])
fig_5.tight_layout()

# fig_6 = plt.figure(figsize=(6.5,4.5))
# ax_6 = fig_6.add_subplot(111)
# ax_6.plot(e[0],dis_no_inc[0], 'None', linestyle='-', color='maroon',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='maroon',
#           markeredgewidth=1)
# ax_6.plot(e[0],dis[0], 'None', linestyle='-', color='steelblue',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='steelblue',
#           markeredgewidth=1)
# ax_6.plot(e[0],dis_no_inc[-1], 'None', linestyle='-.', color='maroon',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='maroon',
#           markeredgewidth=1)
# ax_6.plot(e[0],dis[-1], 'None', linestyle='-.', color='steelblue',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='steelblue',
#           markeredgewidth=1)
# ax_6.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
# ax_6.set_xlabel('Energy (eV)', fontsize=16, labelpad=10)
# ax_6.tick_params(axis="x", labelsize=16)
# ax_6.set_ylabel('Fermi distribution', fontsize=16, labelpad=10)
# ax_6.tick_params(axis="y", labelsize=16)
# fig_6.tight_layout()

# fig_7 = plt.figure(figsize=(6.5,4.5))
# ax_7 = fig_7.add_subplot(111)
# ax_7.plot(e[0],dfdE_no_inc[0], 'None', linestyle='-', color='maroon',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='maroon',
#           markeredgewidth=1)
# ax_7.plot(e[0],dfdE[0], 'None', linestyle='-', color='steelblue',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='steelblue',
#           markeredgewidth=1)
# ax_7.plot(e[0],dfdE_no_inc[-1], 'None', linestyle='-.', color='maroon',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='maroon',
#           markeredgewidth=1)
# ax_7.plot(e[0],dfdE[-1], 'None', linestyle='-.', color='steelblue',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='steelblue',
#           markeredgewidth=1)
# ax_7.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
# ax_7.set_xlabel('Energy (eV)', fontsize=16, labelpad=10)
# ax_7.tick_params(axis="x", labelsize=16)
# ax_7.set_ylabel('Fermi window (eV$^{-1}$)', fontsize=16)
# ax_7.tick_params(axis="y", labelsize=16)
# fig_7.tight_layout()

# fig_8 = plt.figure(figsize=(6.5,4.5))
# ax_8 = fig_8.add_subplot(111)
# ax_8.plot(g[0],cc_no_inc[0], 'o', linestyle='None', color='black',
#           markersize=12, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='black',
#           markeredgewidth=1,zorder=0)
# ax_8.plot(g[0],JD_n_no_inc[0], 'o', linestyle='-.', color='steelblue',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='steelblue',
#           markeredgewidth=1)

# ax_8.plot(g[0],cc_sc_no_inc[0], 'o', linestyle='-.', color='maroon',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='maroon',
#           markeredgewidth=1)

# ax_8.plot(g[0],cc[0], 'o', linestyle='None', color='black',
#           markersize=12, linewidth=1.5,
#           markerfacecolor='black',
#           markeredgecolor='black',
#           markeredgewidth=1,zorder=0)
# ax_8.plot(g[0],JD_n[0], 'o', linestyle='-', color='steelblue',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='steelblue',
#           markeredgewidth=1)

# ax_8.plot(g[0],cc_sc[0], 'o', linestyle='-', color='maroon',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='maroon',
#           markeredgewidth=1)

# ax_8.yaxis.set_major_formatter(ScalarFormatter())
# ax_8.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
# ax_8.tick_params(axis="x", labelsize=16)
# ax_8.set_ylabel('Carrier concentration( #/m$^{3}$)', fontsize=16)
# ax_8.tick_params(axis="y", labelsize=16)
# ax_8.ticklabel_format(axis="y", style="sci", scilimits=None)
# fig_8.tight_layout()


fig_8_2 = plt.figure(figsize=(6.5,4.5))
ax_8_2 = fig_8_2.add_subplot(111)
ax_8_2.plot(g[0],cc_no_inc[0], 'D', linestyle='None', color='black',
          markersize=12, linewidth=1.5,
          markerfacecolor='black',
          markeredgecolor='black',
          markeredgewidth=1,zorder=0)

ax_8_2.plot(g[0],cc_sc_no_inc[0], 'D', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)

ax_8_2.plot(g[0],cc[0], 'o', linestyle='None', color='black',
          markersize=12, linewidth=1.5,
          markerfacecolor='black',
          markeredgecolor='black',
          markeredgewidth=1,zorder=0)

ax_8_2.plot(g[0],cc_sc[0], 'o', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)

ax_8_2.plot(g[0],cc_direction_down[0], 'p', linestyle='None', color='black',
          markersize=12, linewidth=1.5,
          markerfacecolor='black',
          markeredgecolor='black',
          markeredgewidth=1,zorder=0)

ax_8_2.plot(g[0],cc_sc_direction_down[0], 'p', linestyle='-', color='indigo',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='indigo',
          markeredgewidth=1)

ax_8_2.plot(g[0],cc_1pct[0], 's', linestyle='None', color='black',
          markersize=12, linewidth=1.5,
          markerfacecolor='black',
          markeredgecolor='black',
          markeredgewidth=1,zorder=0)

ax_8_2.plot(g[0],cc_sc_1pct[0], 's', linestyle='-', color='olive',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='olive',
          markeredgewidth=1)

ax_8_2.yaxis.set_major_formatter(ScalarFormatter())
ax_8_2.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_8_2.tick_params(axis="x", labelsize=16)
ax_8_2.set_ylabel('Carrier concentration( #/m$^{3}$)', fontsize=16)
ax_8_2.tick_params(axis="y", labelsize=16)
ax_8_2.ticklabel_format(axis="y", style="sci", scilimits=None)
fig_8_2.tight_layout()

np.savetxt('fig3_cc_no_inc',cc_no_inc[0])
np.savetxt('fig3_cc_sc_no_inc',cc_sc_no_inc[0])
np.savetxt('fig3_cc',cc[0])
np.savetxt('fig3_cc_sc',cc_sc[0])
np.savetxt('fig3_cc_direction_down',cc_direction_down[0])
np.savetxt('fig3_cc_sc_direction_down',cc_sc_direction_down[0])
np.savetxt('fig3_cc_1pct',cc_1pct[0])
np.savetxt('fig3_cc_sc_1pct',cc_sc_1pct[0])

# fig_9 = plt.figure(figsize=(6.5,4.5))
# ax_9 = fig_9.add_subplot(111)
# ax_9.plot(g[0],fermi_no_inc[0], 'o', linestyle='-', color='maroon',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='maroon',
#           markeredgewidth=1,zorder=10)
# ax_9.plot(g[0],JD_f_no_inc[0], 'o', linestyle='-', color='steelblue',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='steelblue',
#           markeredgewidth=1)
# ax_9.plot(g[0],fermi[0], 'o', linestyle='-.', color='maroon',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='maroon',
#           markeredgewidth=1,zorder=10)
# ax_9.plot(g[0],JD_f[0], 'o', linestyle='-.', color='steelblue',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='steelblue',
#           markeredgewidth=1)

# ax_9.yaxis.set_major_formatter(ScalarFormatter())
# ax_9.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
# ax_9.tick_params(axis="x", labelsize=16)
# ax_9.set_ylabel('E$_f$ (eV)', fontsize=16)
# ax_9.tick_params(axis="y", labelsize=16)
# ax_9.ticklabel_format(axis="y", style="sci", scilimits=None)
# fig_9.tight_layout()


fig_9_2 = plt.figure(figsize=(6.5,4.5))
ax_9_2 = fig_9_2.add_subplot(111)
ax_9_2.plot(g[0],fermi_no_inc[0], 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1,zorder=10)
ax_9_2.plot(g[0],fermi[0], 'o', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1,zorder=10)
ax_9_2.plot(g[0],fermi_direction_down[0], 'o', linestyle='-', color='indigo',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='indigo',
          markeredgewidth=1,zorder=10)
ax_9_2.plot(g[0],fermi_1pct[0], 'o', linestyle='-', color='olive',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='olive',
          markeredgewidth=1,zorder=10)

ax_9_2.yaxis.set_major_formatter(ScalarFormatter())
ax_9_2.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_9_2.tick_params(axis="x", labelsize=16)
ax_9_2.set_ylabel('E$_f$ (eV)', fontsize=16)
ax_9_2.tick_params(axis="y", labelsize=16)
ax_9_2.ticklabel_format(axis="y", style="sci", scilimits=None)
fig_9_2.tight_layout()


np.savetxt('fig3_fermi_no_inc',fermi_no_inc[0])
np.savetxt('fig3_fermi',fermi[0])
np.savetxt('fig3_fermi_direction_down',fermi_direction_down[0])
np.savetxt('fig3_fermi_1pct',fermi_1pct[0])

fig_9_3 = plt.figure(figsize=(6.5,4.5))
ax_9_3 = fig_9_3.add_subplot(111)
ax_9_3.plot(g[0],m_CB_no_inc[0]/Si.me, 'o', linestyle='-', color='maroon',
          markersize=12, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1,zorder=10)
ax_9_3.plot(g[0],m_CB_inc[0]/Si.me, 'o', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1,zorder=10)
ax_9_3.plot(g[0],m_CB_inc_direction_down[0]/Si.me, 'o', linestyle='-', color='indigo',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='indigo',
          markeredgewidth=1,zorder=10)
ax_9_3.plot(g[0],m_CB_inc_1pct[0]/Si.me, 'o', linestyle='-', color='olive',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='olive',
          markeredgewidth=1,zorder=10)

ax_9_3.yaxis.set_major_formatter(ScalarFormatter())
ax_9_3.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_9_3.tick_params(axis="x", labelsize=16)
ax_9_3.set_ylabel('Effective mass (m$_e$)', fontsize=16)
ax_9_3.tick_params(axis="y", labelsize=16)
ax_9_3.ticklabel_format(axis="y", style="sci", scilimits=None)
fig_9_3.tight_layout()




fig_10 = plt.figure(figsize=(6.5,4.5))
ax_10 = fig_10.add_subplot(111)
ax_10.plot(g[0],Coeff_no_inc[0]*1e-5, 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_10.plot(g[0],Coeff_no_np[0]*1e-5, 'o', linestyle='-', color='tan',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='tan',
          markeredgewidth=1)
ax_10.plot(g[0],(1-vfrac)*Coeff[0]*1e-5, 'o', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)
ax_10.plot(ExpData_SiCfra_0pct_direction_up[:,0], ExpData_SiCfra_0pct_direction_up[:,2]*1e3 * 1e-5, '-o', linestyle='None', color='black',
        markersize=5, linewidth=4,
        markerfacecolor='white',
        markeredgecolor='black',
        markeredgewidth=1)
ax_10.plot(ExpData_SiCfrac_5pct_direction_up[:,0],  ExpData_SiCfrac_5pct_direction_up[:,2]*1e3 * 1e-5, '-o', linestyle='None', color='olive',
        markersize=5, linewidth=4,
        markerfacecolor='white',
        markeredgecolor='olive',
        markeredgewidth=1)
ax_10.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_10.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_10.tick_params(axis="x", labelsize=16)
ax_10.set_ylabel('Conductivity(x10$^5$ S/m)', fontsize=16)
ax_10.tick_params(axis="y", labelsize=16)
fig_10.tight_layout()

np.savetxt('fig5_Sigma_no_np',Coeff_no_np[0]*1e-5)
np.savetxt('fig5_Sigma',(1-vfrac)*Coeff[0]*1e-5)
np.savetxt('fig5_Sigma_ExpData_SiCfrac_5pct_direction_up',ExpData_SiCfrac_5pct_direction_up[:,2]*1e3 * 1e-5)


fig_10_2 = plt.figure(figsize=(6.5,4.5))
ax_10_2 = fig_10_2.add_subplot(111)
ax_10_2.plot(g[0],Coeff_no_inc[0]*1e-5, 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_10_2.plot(g[0],Coeff_direction_down_no_np[0]*1e-5, 'o', linestyle='-', color='tan',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='tan',
          markeredgewidth=1)
ax_10_2.plot(g[0],(1-vfrac)*Coeff_direction_down[0]*1e-5, 'o', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)
ax_10_2.plot(ExpData_SiCfra_0pct_direction_up[:,0], ExpData_SiCfra_0pct_direction_up[:,2]*1e3 * 1e-5, '-o', linestyle='None', color='black',
        markersize=5, linewidth=4,
        markerfacecolor='white',
        markeredgecolor='black',
        markeredgewidth=1)
ax_10_2.plot(ExpData_SiCfrac_5pct_direction_down[:,0],  ExpData_SiCfrac_5pct_direction_down[:,2]*1e3 * 1e-5, '-o', linestyle='None', color='olive',
        markersize=5, linewidth=4,
        markerfacecolor='white',
        markeredgecolor='olive',
        markeredgewidth=1)
ax_10_2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_10_2.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_10_2.tick_params(axis="x", labelsize=16)
ax_10_2.set_ylabel('Conductivity(x10$^5$ S/m)', fontsize=16)
ax_10_2.tick_params(axis="y", labelsize=16)
fig_10_2.tight_layout()


np.savetxt('fig5_Sigma_direction_down_no_np',Coeff_direction_down_no_np[0]*1e-5)
np.savetxt('fig5_Sigma_direction_down',(1-vfrac)*Coeff_direction_down[0]*1e-5)
np.savetxt('fig5_Sigma_ExpData_SiCfrac_5pct_direction_down',ExpData_SiCfrac_5pct_direction_down[:,2]*1e3 * 1e-5)

fig_10_3 = plt.figure(figsize=(6.5,4.5))
ax_10_3 = fig_10_3.add_subplot(111)
ax_10_3.plot(g[0],Coeff_no_inc[0]*1e-5, 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_10_3.plot(g[0],Coeff_no_np_1pct[0]*1e-5, 'o', linestyle='-', color='tan',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='tan',
          markeredgewidth=1)
ax_10_3.plot(g[0],(1-vfrac/5)*Coeff_1pct[0]*1e-5, 'o', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)
ax_10_3.plot(ExpData_SiCfra_0pct_direction_up[:,0], ExpData_SiCfra_0pct_direction_up[:,2]*1e3 * 1e-5, '-o', linestyle='None', color='black',
        markersize=5, linewidth=4,
        markerfacecolor='white',
        markeredgecolor='black',
        markeredgewidth=1)
ax_10_3.plot(ExpData_SiCfrac_1pct_direction_up[:,0],  ExpData_SiCfrac_1pct_direction_up[:,2]*1e3 * 1e-5, '-o', linestyle='None', color='olive',
        markersize=5, linewidth=4,
        markerfacecolor='white',
        markeredgecolor='olive',
        markeredgewidth=1)
ax_10_3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_10_3.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_10_3.tick_params(axis="x", labelsize=16)
ax_10_3.set_ylabel('Conductivity(x10$^5$ S/m)', fontsize=16)
ax_10_3.tick_params(axis="y", labelsize=16)
fig_10_3.tight_layout()

np.savetxt('fig5_Sigma_no_inc',Coeff_no_inc[0]*1e-5)
np.savetxt('fig5_Sigma_no_np_1pct',Coeff_no_np_1pct[0]*1e-5)
np.savetxt('fig5_Sigma_1pct',(1-vfrac/5)*Coeff_1pct[0]*1e-5)
np.savetxt('fig5_Sigma_ExpData_SiCfra_0pct_direction_up',ExpData_SiCfra_0pct_direction_up[:,2]*1e3 * 1e-5)
np.savetxt('fig5_Sigma_ExpData_SiCfrac_1pct_direction_up',ExpData_SiCfrac_1pct_direction_up[:,2]*1e3 * 1e-5)


fig_11 = plt.figure(figsize=(6.5,4.5))
ax_11 = fig_11.add_subplot(111)
ax_11.plot(g[0],Coeff_no_inc[1]*1e6, 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_11.plot(g[0],Coeff_no_np[1]*1e6, 'o', linestyle='-', color='tan',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='tan',
          markeredgewidth=1)
ax_11.plot(g[0],Coeff[1]*1e6, 'o', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)
ax_11.plot(ExpData_SiCfra_0pct_direction_up[:,0], ExpData_SiCfra_0pct_direction_up[:,3], '-o', linestyle='None', color='black',
          markersize=5, linewidth=4,
          markerfacecolor='white',
          markeredgecolor='black',
          markeredgewidth=1)
ax_11.plot(ExpData_SiCfrac_5pct_direction_up[:,0], ExpData_SiCfrac_5pct_direction_up[:,3], '-o', linestyle='None', color='olive',
          markersize=5, linewidth=4,
          markerfacecolor='white',
          markeredgecolor='olive',
          markeredgewidth=1)

ax_11.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax_11.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_11.tick_params(axis="x", labelsize=16)
ax_11.set_ylabel('Seebeck($\mu$V/K)', fontsize=16)
ax_11.tick_params(axis="y", labelsize=16)
fig_11.tight_layout()


np.savetxt('fig5_S_no_np',Coeff_no_np[1]*1e6)
np.savetxt('fig5_S',Coeff[1]*1e6)
np.savetxt('fig5_S_ExpData_SiCfrac_5pct_direction_up',ExpData_SiCfrac_5pct_direction_up[:,3])

fig_11_2 = plt.figure(figsize=(6.5,4.5))
ax_11_2 = fig_11_2.add_subplot(111)
ax_11_2.plot(g[0],Coeff_no_inc[1]*1e6, 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_11_2.plot(g[0],Coeff_direction_down_no_np[1]*1e6, 'o', linestyle='-', color='tan',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='tan',
          markeredgewidth=1)
ax_11_2.plot(g[0],Coeff_direction_down[1]*1e6, 'o', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)
ax_11_2.plot(ExpData_SiCfra_0pct_direction_up[:,0], ExpData_SiCfra_0pct_direction_up[:,3], '-o', linestyle='None', color='black',
          markersize=5, linewidth=4,
          markerfacecolor='white',
          markeredgecolor='black',
          markeredgewidth=1)
ax_11_2.plot(ExpData_SiCfrac_5pct_direction_down[:,0], ExpData_SiCfrac_5pct_direction_down[:,3], '-o', linestyle='None', color='olive',
          markersize=5, linewidth=4,
          markerfacecolor='white',
          markeredgecolor='olive',
          markeredgewidth=1)

ax_11_2.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax_11_2.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_11_2.tick_params(axis="x", labelsize=16)
ax_11_2.set_ylabel('Seebeck($\mu$V/K)', fontsize=16)
ax_11_2.tick_params(axis="y", labelsize=16)
fig_11_2.tight_layout()


np.savetxt('fig5_S_direction_down_no_np',Coeff_direction_down_no_np[1]*1e6)
np.savetxt('fig5_S_direction_down',Coeff_direction_down[1]*1e6)
np.savetxt('fig5_S_ExpData_SiCfrac_5pct_direction_down',ExpData_SiCfrac_5pct_direction_down[:,3])



fig_11_3 = plt.figure(figsize=(6.5,4.5))
ax_11_3 = fig_11_3.add_subplot(111)
ax_11_3.plot(g[0],Coeff_no_inc[1]*1e6, 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_11_3.plot(g[0],Coeff_no_np_1pct[1]*1e6, 'o', linestyle='-', color='tan',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='tan',
          markeredgewidth=1)
ax_11_3.plot(g[0],Coeff_1pct[1]*1e6, 'o', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)
ax_11_3.plot(ExpData_SiCfra_0pct_direction_up[:,0], ExpData_SiCfra_0pct_direction_up[:,3], '-o', linestyle='None', color='black',
          markersize=5, linewidth=4,
          markerfacecolor='white',
          markeredgecolor='black',
          markeredgewidth=1)
ax_11_3.plot(ExpData_SiCfrac_1pct_direction_up[:,0], ExpData_SiCfrac_1pct_direction_up[:,3], '-o', linestyle='None', color='olive',
          markersize=5, linewidth=4,
          markerfacecolor='white',
          markeredgecolor='olive',
          markeredgewidth=1)

ax_11_3.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax_11_3.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_11_3.tick_params(axis="x", labelsize=16)
ax_11_3.set_ylabel('Seebeck($\mu$V/K)', fontsize=16)
ax_11_3.tick_params(axis="y", labelsize=16)
fig_11_3.tight_layout()


np.savetxt('fig5_S_no_inc',Coeff_no_inc[1]*1e6)
np.savetxt('fig5_S_no_np_1pct',Coeff_no_np_1pct[1]*1e6)
np.savetxt('fig5_S_Coeff_1pct',Coeff_1pct[1]*1e6)
np.savetxt('fig5_S_ExpData_SiCfra_0pct_direction_up',ExpData_SiCfra_0pct_direction_up[:,3])
np.savetxt('fig5_S_ExpData_SiCfrac_1pct_direction_up',ExpData_SiCfrac_1pct_direction_up[:,3])



fig_12 = plt.figure(figsize=(6.5,4.5))
ax_12 = fig_12.add_subplot(111)
ax_12.plot(g[0],Coeff_no_inc[2]*1e3, 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_12.plot(g[0],Coeff_no_np[2]*1e3, 'o', linestyle='-', color='tan',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='tan',
          markeredgewidth=1)
ax_12.plot(g[0],(1-vfrac)*Coeff[2]*1e3, 'o', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)
ax_12.plot(ExpData_SiCfra_0pct_direction_up[:,0], ExpData_SiCfra_0pct_direction_up[:,4], '-o', linestyle='None', color='black',
          markersize=5, linewidth=4,
          markerfacecolor='white',
          markeredgecolor='black',
          markeredgewidth=1)
ax_12.plot(ExpData_SiCfrac_5pct_direction_up[:,0], ExpData_SiCfrac_5pct_direction_up[:,4], '-o', linestyle='None', color='olive',
          markersize=5, linewidth=4,
          markerfacecolor='white',
          markeredgecolor='olive',
          markeredgewidth=1)

ax_12.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_12.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_12.tick_params(axis="x", labelsize=16)
ax_12.set_ylabel('Power factor(mW/mK$^2$)', fontsize=16, labelpad=10)
ax_12.tick_params(axis="y", labelsize=16)
fig_12.tight_layout()


np.savetxt('fig5_PF_no_np',Coeff_no_np[2]*1e3)
np.savetxt('fig5_PF',(1-vfrac)*Coeff[2]*1e3)
np.savetxt('fig5_Temp_ExpData_SiCfrac_5pct_direction_up',ExpData_SiCfrac_5pct_direction_up[:,0])
np.savetxt('fig5_PF_ExpData_SiCfrac_5pct_direction_up',ExpData_SiCfrac_5pct_direction_up[:,4])



fig_12_2 = plt.figure(figsize=(6.5,4.5))
ax_12_2 = fig_12_2.add_subplot(111)
ax_12_2.plot(g[0],Coeff_no_inc[2]*1e3, 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_12_2.plot(g[0],Coeff_direction_down_no_np[2]*1e3, 'o', linestyle='-', color='tan',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='tan',
          markeredgewidth=1)
ax_12_2.plot(g[0],(1-vfrac)*Coeff_direction_down[2]*1e3, 'o', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)
ax_12_2.plot(ExpData_SiCfra_0pct_direction_up[:,0], ExpData_SiCfra_0pct_direction_up[:,4], '-o', linestyle='None', color='black',
          markersize=5, linewidth=4,
          markerfacecolor='white',
          markeredgecolor='black',
          markeredgewidth=1)
ax_12_2.plot(ExpData_SiCfrac_5pct_direction_down[:,0], ExpData_SiCfrac_5pct_direction_down[:,4], '-o', linestyle='None', color='olive',
          markersize=5, linewidth=4,
          markerfacecolor='white',
          markeredgecolor='olive',
          markeredgewidth=1)


np.savetxt('fig5_PF_direction_down_no_np',Coeff_direction_down_no_np[2]*1e3)
np.savetxt('fig5_PF_direction_down',(1-vfrac)*Coeff_direction_down[2]*1e3)
np.savetxt('fig5_Temp_ExpData_SiCfrac_5pct_direction_down',ExpData_SiCfrac_5pct_direction_down[:,0])
np.savetxt('fig5_PF_ExpData_SiCfrac_5pct_direction_down',ExpData_SiCfrac_5pct_direction_down[:,4])

ax_12_2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_12_2.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_12_2.tick_params(axis="x", labelsize=16)
ax_12_2.set_ylabel('Power factor(mW/mK$^2$)', fontsize=16, labelpad=10)
ax_12_2.tick_params(axis="y", labelsize=16)
fig_12_2.tight_layout()

fig_12_3 = plt.figure(figsize=(6.5,4.5))
ax_12_3 = fig_12_3.add_subplot(111)
ax_12_3.plot(g[0],Coeff_no_inc[2]*1e3, 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_12_3.plot(g[0],Coeff_no_np_1pct[2]*1e3, 'o', linestyle='-', color='tan',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='tan',
          markeredgewidth=1)
ax_12_3.plot(g[0],(1-vfrac/5)*Coeff_1pct[2]*1e3, 'o', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)
ax_12_3.plot(ExpData_SiCfra_0pct_direction_up[:,0], ExpData_SiCfra_0pct_direction_up[:,4], '-o', linestyle='None', color='black',
          markersize=5, linewidth=4,
          markerfacecolor='white',
          markeredgecolor='black',
          markeredgewidth=1)
ax_12_3.plot(ExpData_SiCfrac_1pct_direction_up[:,0], ExpData_SiCfrac_1pct_direction_up[:,4], '-o', linestyle='None', color='olive',
          markersize=5, linewidth=4,
          markerfacecolor='white',
          markeredgecolor='olive',
          markeredgewidth=1)

ax_12_3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_12_3.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_12_3.tick_params(axis="x", labelsize=16)
ax_12_3.set_ylabel('Power factor(mW/mK$^2$)', fontsize=16, labelpad=10)
ax_12_3.tick_params(axis="y", labelsize=16)
fig_12_3.tight_layout()

np.savetxt('fig5_PF_no_inc',Coeff_no_inc[2]*1e3)
np.savetxt('fig5_PF_no_np_1pct',Coeff_no_np_1pct[2]*1e3)
np.savetxt('fig5_PF_Coeff_1pct',(1-vfrac/5)*Coeff_1pct[2]*1e3)
np.savetxt('fig5_Temp_ExpData_SiCfra_0pct_direction_up',ExpData_SiCfra_0pct_direction_up[:,0])
np.savetxt('fig5_PF_ExpData_SiCfra_0pct_direction_up',ExpData_SiCfra_0pct_direction_up[:,4])
np.savetxt('fig5_Temp_ExpData_SiCfrac_1pct_direction_up',ExpData_SiCfrac_1pct_direction_up[:,0])
np.savetxt('fig5_PF_ExpData_SiCfrac_1pct_direction_up',ExpData_SiCfrac_1pct_direction_up[:,4])

# fig_13 = plt.figure(figsize=(6.5,4.5))
# ax_13 = fig_13.add_subplot(111)

# ax_13.plot(g[0],Coeff_no_inc[3], 'o', linestyle='-', color='maroon',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='maroon',
#           markeredgewidth=1)
# ax_13.plot(g[0],Coeff[3], 'o', linestyle='-', color='steelblue',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='steelblue',
#           markeredgewidth=1)
# # ax_13.plot(g[0],Coeff[0]*g[0]*Coeff[6], 'o', linestyle='-', color='steelblue',
# #           markersize=6, linewidth=1.5,
# #           markerfacecolor='white',
# #           markeredgecolor='steelblue',
# #           markeredgewidth=1)

# ax_13.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
# ax_13.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
# ax_13.tick_params(axis="x", labelsize=16)
# ax_13.set_ylabel('$\kappa_e$(W/mK)', fontsize=16, labelpad=10)
# ax_13.tick_params(axis="y", labelsize=16)
# fig_13.tight_layout()

# fig_14 = plt.figure(figsize=(6.5,4.5))
# ax_14 = fig_14.add_subplot(111)
# ax_14.plot(g[0],Coeff_no_inc[4], 'o', linestyle='-', color='maroon',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='maroon',
#           markeredgewidth=1)
# ax_14.plot(g[0],Coeff[4], 'o', linestyle='-', color='steelblue',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='steelblue',
#           markeredgewidth=1)
# ax_14.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
# ax_14.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
# ax_14.tick_params(axis="x", labelsize=16)
# ax_14.set_ylabel('$\Delta_1$(eV)', fontsize=16, labelpad=10)
# ax_14.tick_params(axis="y", labelsize=16)
# fig_14.tight_layout()

# fig_15 = plt.figure(figsize=(6.5,4.5))
# ax_15 = fig_15.add_subplot(111)
# ax_15.plot(g[0],Coeff_no_inc[5], 'o', linestyle='-', color='maroon',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='maroon',
#           markeredgewidth=1)
# ax_15.plot(g[0],Coeff[5], 'o', linestyle='-', color='steelblue',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='steelblue',
#           markeredgewidth=1)
# ax_15.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
# ax_15.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
# ax_15.tick_params(axis="x", labelsize=16)
# ax_15.set_ylabel('$\Delta_2$([eV]$^2$)', fontsize=16, labelpad=10)
# ax_15.tick_params(axis="y", labelsize=16)
# fig_15.tight_layout()

# fig_16 = plt.figure(figsize=(6.5,4.5))
# ax_16 = fig_16.add_subplot(111)
# ax_16.plot(g[0],Coeff_no_inc[6]*1e8, 'o', linestyle='-', color='maroon',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='maroon',
#           markeredgewidth=1)
# ax_16.plot(g[0],Coeff[6]*1e8, 'o', linestyle='-', color='steelblue',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='steelblue',
#           markeredgewidth=1)
# ax_16.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
# ax_16.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
# ax_16.tick_params(axis="x", labelsize=16)
# ax_16.set_ylabel('Lorenz number (x$10^{-8}$[V/K]$^2$)', fontsize=16, labelpad=10)
# ax_16.tick_params(axis="y", labelsize=16)
# fig_16.tight_layout()

fig_17 = plt.figure(figsize=(6.5,4.5))
ax_17 = fig_17.add_subplot(111)
ax_17.plot(g[0],LD_nondegenerate_no_inc[0]*1e9, 'o', linestyle='-.', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)

ax_17.plot(g[0],LD_no_inc[0]*1e9, 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)

ax_17.plot(g[0],LD_nondegenerate[0]*1e9, 'o', linestyle='-.', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)

ax_17.plot(g[0],LD[0]*1e9, 'o', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)


ax_17.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_17.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_17.tick_params(axis="x", labelsize=16)
ax_17.set_ylabel('Debyle length (nm)', fontsize=16, labelpad=10)
ax_17.tick_params(axis="y", labelsize=16)
fig_17.tight_layout()

fig_18 = plt.figure(figsize=(6.5,4.5))
ax_18 = fig_18.add_subplot(111)
ax_18.plot(e[0],np.log10(6*tau_p_npb[5])+15, 'None', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_18.plot(e[0],np.log10(6*tau_ion[5])+15, 'None', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)

ax_18.plot(e[0,::1],np.log10(tau_np[0,::1])+15, 'o', linestyle=None, color='tan',
          markersize=8, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='tan',
          markeredgewidth=1)
ax_18.plot(e[0,::1],np.log10(tau_g[0,::1])+15, 'o', linestyle=None, color='tan',
          markersize=8, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='tan',
          markeredgewidth=1)

ax_18.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax_18.set_xlabel('Energy (eV)', fontsize=16, labelpad=10)
ax_18.tick_params(axis="x", labelsize=16)
ax_18.set_ylabel('Lifetime (log$_{10}$[ps])', fontsize=16, labelpad=10)
ax_18.tick_params(axis="y", labelsize=16)
fig_18.tight_layout()

plt.show()
# exit()

print(np.shape(tau_p_npb),np.shape(tau_ion),np.shape(tau_np),np.shape(e))

np.savetxt('fig4_tau_p_npb',6*tau_p_npb[5])
np.savetxt('fig4_tau_ion',6*tau_ion[5])
np.savetxt('fig4_tau_np',tau_np[0,:])
np.savetxt('fig4_tau_grain',tau_g[0])
np.savetxt('fig4_energy',e[0])
