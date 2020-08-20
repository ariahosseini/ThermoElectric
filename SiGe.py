import numpy as np
from numpy.linalg import norm
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
from thermoelectricProperties import thermoelectricProperties

x = 0.3
latticeParameter=(0.027*x**2+0.2*x+5.431)*1e-10
meff = (1.08*(1-x)+1.41*x-0.183*x*(1-x))*thermoelectricProperties.me
dielectric = 11.7+4.5*x

SiGe = thermoelectricProperties(latticeParameter=latticeParameter, dopantElectricCharge=1, electronEffectiveMass= meff, energyMin=0.0, energyMax=1.4, dielectric=dielectric, numKpoints=800, numBands=8, numQpoints=201, numEnergySampling=1000)


ml = 0.92*thermoelectricProperties.me # longitudinal effective mass
mt = 0.19*thermoelectricProperties.me # transverse effective mass
m_CB = 3/(1/ml+2/mt) # conduction band effective mass
bulk_module = 98 - 23*x # Bulk module (GPA)
rho = 2329+3493*x-499*x**2  # mass density (Kg/m3)
sp = np.sqrt(bulk_module/rho) # speed of sound



Lv = np.array([[1,1,0],[0,1,1],[1,0,1]])*SiGe.latticeParameter/2
a_rp = np.cross(Lv[1],Lv[2])/np.dot(Lv[0],np.cross(Lv[1],Lv[2]))
b_rp = np.cross(Lv[2],Lv[0])/np.dot(Lv[1],np.cross(Lv[2],Lv[0]))
a_rp = np.cross(Lv[0],Lv[1])/np.dot(Lv[2],np.cross(Lv[0],Lv[1]))
RLv = np.array([a_rp, b_rp, a_rp])


e = SiGe.energyRange()
g = SiGe.temp(TempMin=300, TempMax=1301, dT=100)
h = SiGe.bandGap(Eg_o=1.17, Ao=4.73e-4, Bo=636, Temp=g)
# alpha = (1-m_CB/thermoelectricProperties.me)**2/h
alpha = np.array(0.5*np.tile([1],(1,11)))

# ro = np.arange(1, 60, 30) * 1e-9
# tau_cy = Si.tau2D_cylinder(energyRange = e, nk= [21,19,19], Uo = 1.5, m = [0.89, 0.19, 0.19], vfrac = 0.15, valley = [0.85,0,0], dk_len=0.15, ro= ro, n=2000)
# tau_sp = Si.tau3D_spherical(energyRange = e, nk= [10,8,8], Uo = 1.389, m = [0.89, 0.19, 0.19], vfrac = 0.15, valley = [0.85,0,0], dk_len=0.15, ro= ro, n=32)

dos_nonparabolic, dos_parabolic = SiGe.analyticalDoS(energyRange=e, alpha = alpha)

cc_circle = SiGe.carrierConcentration(Nc=None, Nv=None, path2extrinsicCarrierConcentration='Vining_CC_circle', bandGap=h, Ao=5.3e21, Bo=3.5e21, Temp=g)
cc_square = SiGe.carrierConcentration(Nc=None, Nv=None, path2extrinsicCarrierConcentration='Vining_CC_square', bandGap=h, Ao=5.3e21, Bo=3.5e21, Temp=g)
cc_diamond = SiGe.carrierConcentration(Nc=None, Nv=None, path2extrinsicCarrierConcentration='Vining_CC_diamond', bandGap=h, Ao=5.3e21, Bo=3.5e21, Temp=g)
cc_triangle = SiGe.carrierConcentration(Nc=None, Nv=None, path2extrinsicCarrierConcentration='Vining_CC_triangle', bandGap=h, Ao=5.3e21, Bo=3.5e21, Temp=g)
kp, band = SiGe.electronBandStructure(path2eigenval='EIGENVAL', skipLines=6)
kp_rl = 2*np.pi*np.matmul(kp,RLv)
kp_mag = norm(kp_rl, axis=1)
min_band = np.argmin(band[400:600, 4], axis=0)
max_band = np.argmax(band[400:600, 4], axis=0)
kp_vel = kp_mag[401 + max_band:401 + min_band]
energy_vel = band[401 + max_band:401 + min_band, 4] - band[401 + min_band, 4]
enrg_sorted_idx = np.argsort(energy_vel, axis=0)
gVel = SiGe.electronGroupVelocity(kp=kp_vel[enrg_sorted_idx], energy_kp=energy_vel[enrg_sorted_idx], energyRange=e)
DoS = SiGe.electronDoS(path2DoS='DOSCAR', headerLines=6, unitcell_volume=2*19.70272e-30, numDoSpoints=2000, valleyPoint=1118, energyRange=e)
JD_f_circle, JD_n_circle = SiGe.fermiLevel(carrierConcentration=cc_circle, energyRange=e, DoS= DoS, Nc=None, Ao=5.3e21, Temp=g)
JD_f_diamond, JD_n_diamond = SiGe.fermiLevel(carrierConcentration=cc_diamond, energyRange=e, DoS= DoS, Nc=None, Ao=5.3e21, Temp=g)
JD_f_square, JD_n_square = SiGe.fermiLevel(carrierConcentration=cc_square, energyRange=e, DoS= DoS, Nc=None, Ao=5.3e21, Temp=g)
JD_f_triangle, JD_n_triangle = SiGe.fermiLevel(carrierConcentration=cc_triangle, energyRange=e, DoS= DoS, Nc=None, Ao=5.3e21, Temp=g)

fermi_circle, cc_sc_circle = SiGe.fermiLevelSelfConsistent(carrierConcentration=cc_circle, Temp=g, energyRange=e, DoS=DoS, fermilevel=JD_f_circle)
fermi_triangle, cc_sc_triangle = SiGe.fermiLevelSelfConsistent(carrierConcentration=cc_triangle, Temp=g, energyRange=e, DoS=DoS, fermilevel=JD_f_triangle)
fermi_square, cc_sc_square = SiGe.fermiLevelSelfConsistent(carrierConcentration=cc_square, Temp=g, energyRange=e, DoS=DoS, fermilevel=JD_f_square)
fermi_diamond, cc_sc_diamond = SiGe.fermiLevelSelfConsistent(carrierConcentration=cc_diamond, Temp=g, energyRange=e, DoS=DoS, fermilevel=JD_f_diamond)

dis_circle, dfdE_circle = SiGe.fermiDistribution(energyRange=e, Temp=g, fermiLevel=fermi_circle)
dis_triangle, dfdE_triangle = SiGe.fermiDistribution(energyRange=e, Temp=g, fermiLevel=fermi_triangle)
dis_square, dfdE_square = SiGe.fermiDistribution(energyRange=e, Temp=g, fermiLevel=fermi_square)
dis_diamond, dfdE_diamond = SiGe.fermiDistribution(energyRange=e, Temp=g, fermiLevel=fermi_diamond)

LD_non_degenerate_circle = np.sqrt(4*np.pi*SiGe.dielectric*thermoelectricProperties.e0*thermoelectricProperties.kB/thermoelectricProperties.e2C*g/cc_circle) # screening length
LD_non_degenerate_triangle = np.sqrt(4*np.pi*SiGe.dielectric*thermoelectricProperties.e0*thermoelectricProperties.kB/thermoelectricProperties.e2C*g/cc_triangle) # screening length
LD_non_degenerate_square = np.sqrt(4*np.pi*SiGe.dielectric*thermoelectricProperties.e0*thermoelectricProperties.kB/thermoelectricProperties.e2C*g/cc_square) # screening length
LD_non_degenerate_diamond = np.sqrt(4*np.pi*SiGe.dielectric*thermoelectricProperties.e0*thermoelectricProperties.kB/thermoelectricProperties.e2C*g/cc_diamond) # screening length

exit()

LD = np.array([1.38e-9,1.36e-9,1.36e-9,1.36e-9,1.36e-9,1.36e-9,1.4e-9,1.45e-9,1.38e-9,1.25e-9,1.26e-9])[None,:]

tau_p_type_0 = 220e-15/np.sqrt(np.transpose(g) * e)
tau_p_pb_type_1, tau_p_npb_type_1 = SiGe.tau_p(energyRange=e, alpha=alpha, Dv=2.94, DA=9.5, T=g, vs=sp, D=DoS, rho=rho)
tau_p_pb_type_2, tau_p_npb_type_2 = SiGe.tau_p(energyRange=e, alpha=alpha, Dv=2.94, DA=9.5, T=g, vs=sp, D=dos_nonparabolic, rho=rho)
tau_p_pb_type_3, tau_p_npb_type_3 = SiGe.tau_p(energyRange=e, alpha=alpha, Dv=2.94, DA=9.5, T=g, vs=sp, D=dos_parabolic, rho=rho)
tau_BH = SiGe.tau_BH(energyRange=e, LD= LD, N= cc)
tau_ion = SiGe.tau_ion(D=DoS, LD=LD, N=cc)
tau = SiGe.matthiessen(e, 6*tau_p_npb_type_2,6*tau_ion)

Coeff = SiGe.electricalProperties(E=e, DoS=DoS, vg=gVel, Ef=fermi, dfdE=dfdE, Temp=g, tau=tau)


f0 = 'Vining_res_circle'
f1 = 'Vining_res_diamond'
f2 = 'Vining_res_square'
f3 = 'Vining_res_triangle'
f4 = 'Vining_seebeck_circle'
f5 = 'Vining_seebeck_diamond'
f6 = 'Vining_seebeck_square'
f7 = 'Vining_seebeck_triangle'
f8 = 'Vining_pf_circle'
f9 = 'Vining_pf_diamond'
f10 = 'Vining_pf_square'
f11 = 'Vining_pf_triangle'
vining_circle = np.loadtxt(f0, comments='#', delimiter=',', converters=None, skiprows=0)
vining_diamond = np.loadtxt(f1, comments='#', delimiter=',', converters=None, skiprows=0)
vining_square = np.loadtxt(f2, comments='#', delimiter=',', converters=None, skiprows=0)
vining_triangle = np.loadtxt(f3, comments='#', delimiter=',', converters=None, skiprows=0)
vining_seebeck_circle = np.loadtxt(f4, comments='#', delimiter=',', converters=None, skiprows=0)
vining_seebeck_diamond = np.loadtxt(f5, comments='#', delimiter=',', converters=None, skiprows=0)
vining_seebeck_square = np.loadtxt(f6, comments='#', delimiter=',', converters=None, skiprows=0)
vining_seebeck_triangle = np.loadtxt(f7, comments='#', delimiter=',', converters=None, skiprows=0)
vining_pf_circle = np.loadtxt(f8, comments='#', converters=None, skiprows=0)
vining_pf_diamond = np.loadtxt(f9, comments='#', converters=None, skiprows=0)
vining_pf_square = np.loadtxt(f10, comments='#', converters=None, skiprows=0)
vining_pf_triangle = np.loadtxt(f11, comments='#', converters=None, skiprows=0)

exp_r = vining_circle
exp_s = vining_seebeck_circle

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

fig_1 = plt.figure(figsize=(6.5,4.5))
ax_1 = fig_1.add_subplot(111)
ax_1.plot(g[0],alpha[0], 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax_1.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_1.tick_params(axis="x", labelsize=16)
ax_1.set_ylabel('Nonparabolic term (eV$^{-1}$)', fontsize=16, labelpad=10)
ax_1.tick_params(axis="y", labelsize=16)
fig_1.tight_layout()

fig_2 = plt.figure(figsize=(6.5,4.5))
ax_2 = fig_2.add_subplot(111)
ax_2.plot(e[0],DoS[0], 'None', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_2.plot(e[0],dos_nonparabolic[0], 'None', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)
ax_2.plot(e[0],dos_parabolic[0], 'None', linestyle='-', color='olive',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='olive',
          markeredgewidth=1)
ax_2.yaxis.set_major_formatter(ScalarFormatter())
ax_2.set_xlabel('Energy (eV)', fontsize=16, labelpad=10)
ax_2.tick_params(axis="x", labelsize=16)
ax_2.set_ylabel('Density of state (#/eV/m$^3$)', fontsize=16, labelpad=10)
ax_2.tick_params(axis="y", labelsize=16)
ax_2.ticklabel_format(axis="y", style="sci", scilimits=None)
fig_2.tight_layout()

fig_3 = plt.figure(figsize=(6.5,4.5))
ax_3 = fig_3.add_subplot(111)
ax_3.plot(e[0],dos_nonparabolic[0], 'None', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_3.plot(e[0],dos_nonparabolic[-1], 'None', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)
ax_3.yaxis.set_major_formatter(ScalarFormatter())
ax_3.set_xlabel('Energy (eV)', fontsize=16, labelpad=10)
ax_3.tick_params(axis="x", labelsize=16)
ax_3.set_ylabel('Density of state (#/eV/m$^3$)', fontsize=16, labelpad=10)
ax_3.tick_params(axis="y", labelsize=16)
ax_3.ticklabel_format(axis="y", style="sci", scilimits=None)
fig_3.tight_layout()

fig_4 = plt.figure(figsize=(6.5,4.5))
ax_4 = fig_4.add_subplot(111)
ax_4.plot(e[0],-1*gVel[0]*1e-5, 'None', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_4.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax_4.set_xlabel('Energy (eV)', fontsize=16, labelpad=10)
ax_4.tick_params(axis="x", labelsize=16)
ax_4.set_ylabel('Group velocity (x10$^5$ m/s)', fontsize=16, labelpad=10)
fig_4.tight_layout()

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

fig_6 = plt.figure(figsize=(6.5,4.5))
ax_6 = fig_6.add_subplot(111)
ax_6.plot(e[0],dis[0], 'None', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_6.plot(e[0],dis[-1], 'None', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)
ax_6.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_6.set_xlabel('Energy (eV)', fontsize=16, labelpad=10)
ax_6.tick_params(axis="x", labelsize=16)
ax_6.set_ylabel('Fermi distribution', fontsize=16, labelpad=10)
ax_6.tick_params(axis="y", labelsize=16)
fig_6.tight_layout()

fig_7 = plt.figure(figsize=(6.5,4.5))
ax_7 = fig_7.add_subplot(111)
ax_7.plot(e[0],dfdE[0], 'None', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_7.plot(e[0],dfdE[-1], 'None', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)
ax_7.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_7.set_xlabel('Energy (eV)', fontsize=16, labelpad=10)
ax_7.tick_params(axis="x", labelsize=16)
ax_7.set_ylabel('Fermi window (eV$^{-1}$)', fontsize=16)
ax_7.tick_params(axis="y", labelsize=16)
fig_7.tight_layout()

fig_8 = plt.figure(figsize=(6.5,4.5))
ax_8 = fig_8.add_subplot(111)
ax_8.plot(g[0],cc[0], 'o', linestyle='None', color='black',
          markersize=12, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='black',
          markeredgewidth=1,zorder=0)
ax_8.plot(g[0],JD_n[0], 'o', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)

ax_8.plot(g[0],cc_sc[0], 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)

ax_8.yaxis.set_major_formatter(ScalarFormatter())
ax_8.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_8.tick_params(axis="x", labelsize=16)
ax_8.set_ylabel('Carrier concentration( #/m$^{3}$)', fontsize=16)
ax_8.tick_params(axis="y", labelsize=16)
ax_8.ticklabel_format(axis="y", style="sci", scilimits=None)
fig_8.tight_layout()

fig_9 = plt.figure(figsize=(6.5,4.5))
ax_9 = fig_9.add_subplot(111)
ax_9.plot(g[0],fermi[0], 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1,zorder=10)
ax_9.plot(g[0],JD_f[0], 'o', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)

ax_9.yaxis.set_major_formatter(ScalarFormatter())
ax_9.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_9.tick_params(axis="x", labelsize=16)
ax_9.set_ylabel('E$_f$ (eV)', fontsize=16)
ax_9.tick_params(axis="y", labelsize=16)
ax_9.ticklabel_format(axis="y", style="sci", scilimits=None)
fig_9.tight_layout()

fig_10 = plt.figure(figsize=(6.5,4.5))
ax_10 = fig_10.add_subplot(111)
ax_10.plot(g[0],Coeff[0]*1e-5, 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_10.plot(exp_r[:,0], np.divide(1, exp_r[:,1]), '-o', linestyle='None', color='black',
        markersize=5, linewidth=4,
        markerfacecolor='white',
        markeredgecolor='black',
        markeredgewidth=1)
ax_10.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_10.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_10.tick_params(axis="x", labelsize=16)
ax_10.set_ylabel('Conductivity(x10$^5$ S/m)', fontsize=16)
ax_10.tick_params(axis="y", labelsize=16)
fig_10.tight_layout()

fig_11 = plt.figure(figsize=(6.5,4.5))
ax_11 = fig_11.add_subplot(111)
ax_11.plot(g[0],Coeff[1]*1e6, 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_11.plot(exp_s[:,0], -1*exp_s[:,1], '-o', linestyle='None', color='black',
          markersize=5, linewidth=4,
          markerfacecolor='white',
          markeredgecolor='black',
          markeredgewidth=1)

ax_11.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax_11.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_11.tick_params(axis="x", labelsize=16)
ax_11.set_ylabel('Seebeck($\mu$V/K)', fontsize=16)
ax_11.tick_params(axis="y", labelsize=16)
fig_11.tight_layout()

fig_12 = plt.figure(figsize=(6.5,4.5))
ax_12 = fig_12.add_subplot(111)
ax_12.plot(g[0],Coeff[2]*1e3, 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_12.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_12.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_12.tick_params(axis="x", labelsize=16)
ax_12.set_ylabel('Power factor(mW/mK$^2$)', fontsize=16, labelpad=10)
ax_12.tick_params(axis="y", labelsize=16)
fig_12.tight_layout()


fig_13 = plt.figure(figsize=(6.5,4.5))
ax_13 = fig_13.add_subplot(111)
ax_13.plot(g[0],Coeff[3], 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)

ax_13.plot(g[0],Coeff[0]*g[0]*Coeff[6], 'o', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)

ax_13.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_13.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_13.tick_params(axis="x", labelsize=16)
ax_13.set_ylabel('$\kappa_e$(W/mK)', fontsize=16, labelpad=10)
ax_13.tick_params(axis="y", labelsize=16)
fig_13.tight_layout()

fig_14 = plt.figure(figsize=(6.5,4.5))
ax_14 = fig_14.add_subplot(111)
ax_14.plot(g[0],Coeff[4], 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_14.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax_14.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_14.tick_params(axis="x", labelsize=16)
ax_14.set_ylabel('$\Delta_1$(eV)', fontsize=16, labelpad=10)
ax_14.tick_params(axis="y", labelsize=16)
fig_14.tight_layout()

fig_15 = plt.figure(figsize=(6.5,4.5))
ax_15 = fig_15.add_subplot(111)
ax_15.plot(g[0],Coeff[5], 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_15.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax_15.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_15.tick_params(axis="x", labelsize=16)
ax_15.set_ylabel('$\Delta_2$([eV]$^2$)', fontsize=16, labelpad=10)
ax_15.tick_params(axis="y", labelsize=16)
fig_15.tight_layout()


fig_16 = plt.figure(figsize=(6.5,4.5))
ax_16 = fig_16.add_subplot(111)
ax_16.plot(g[0],Coeff[6]*1e8, 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)

ax_16.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_16.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_16.tick_params(axis="x", labelsize=16)
ax_16.set_ylabel('Lorenz number (x$10^{-8}$[V/K]$^2$)', fontsize=16, labelpad=10)
ax_16.tick_params(axis="y", labelsize=16)
fig_16.tight_layout()

fig_17 = plt.figure(figsize=(6.5,4.5))
ax_17 = fig_17.add_subplot(111)
ax_17.plot(e[0],np.log10(tau_p_npb_type_1[0]*6)+15, 'None', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_17.plot(e[0],np.log10(tau_p_pb_type_1[0]*6)+15, 'None', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)
ax_17.plot(e[0],np.log10(tau_p_npb_type_2[0]*6)+15, 'None', linestyle='-', color='tan',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='tan',
          markeredgewidth=1)
ax_17.plot(e[0],np.log10(tau_p_pb_type_2[0]*6)+15, 'None', linestyle='-', color='olive',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='olive',
          markeredgewidth=1)
ax_17.plot(e[0],np.log10(tau_p_npb_type_3[0]*6)+15, 'None', linestyle='-', color='indigo',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='indigo',
          markeredgewidth=1)
ax_17.plot(e[0],np.log10(tau_p_pb_type_3[0]*6)+15, 'None', linestyle='-', color='slateblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='slateblue',
          markeredgewidth=1)
ax_17.plot(e[0],np.log10(tau_ion[0]*6)+15, 'None', linestyle='--', color='gray',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='gray',
          markeredgewidth=1)
ax_17.plot(e[0],np.log10(tau[0])+15, 'None', linestyle='--', color='gold',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='gold',
          markeredgewidth=1)

plt.show()
exit()
