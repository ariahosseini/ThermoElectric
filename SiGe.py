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

SiGe = thermoelectricProperties(latticeParameter=latticeParameter, dopantElectricCharge=1, electronEffectiveMass= meff, energyMin=0.00, energyMax=1.4, dielectric=dielectric, numKpoints=800, numBands=8, numQpoints=201, numEnergySampling=1000)


ml = 0.98*thermoelectricProperties.me # longitudinal effective mass
mt = 0.19*thermoelectricProperties.me # transverse effective mass
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
alpha = np.array(0.5*np.tile([1],(1,11)))
m_CB = 0.26*thermoelectricProperties.me*(1+5*alpha*thermoelectricProperties.kB*g)     # conduction band effective mass
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
JD_f_circle, JD_n_circle = SiGe.fermiLevel(carrierConcentration=cc_circle, energyRange=e, DoS= dos_nonparabolic, Nc=None, Ao=5.3e21, Temp=g)
JD_f_diamond, JD_n_diamond = SiGe.fermiLevel(carrierConcentration=cc_diamond, energyRange=e, DoS= dos_nonparabolic, Nc=None, Ao=5.3e21, Temp=g)
JD_f_square, JD_n_square = SiGe.fermiLevel(carrierConcentration=cc_square, energyRange=e, DoS= dos_nonparabolic, Nc=None, Ao=5.3e21, Temp=g)
JD_f_triangle, JD_n_triangle = SiGe.fermiLevel(carrierConcentration=cc_triangle, energyRange=e, DoS= dos_nonparabolic, Nc=None, Ao=5.3e21, Temp=g)

fermi_circle, cc_sc_circle = SiGe.fermiLevelSelfConsistent(carrierConcentration=cc_circle, Temp=g, energyRange=e, DoS=dos_nonparabolic[0], fermilevel=JD_f_circle)
fermi_triangle, cc_sc_triangle = SiGe.fermiLevelSelfConsistent(carrierConcentration=cc_triangle, Temp=g, energyRange=e, DoS=dos_nonparabolic[0], fermilevel=JD_f_triangle)
fermi_square, cc_sc_square = SiGe.fermiLevelSelfConsistent(carrierConcentration=cc_square, Temp=g, energyRange=e, DoS=dos_nonparabolic[0], fermilevel=JD_f_square)
fermi_diamond, cc_sc_diamond = SiGe.fermiLevelSelfConsistent(carrierConcentration=cc_diamond, Temp=g, energyRange=e, DoS=dos_nonparabolic[0], fermilevel=JD_f_diamond)

dis_circle, dfdE_circle = SiGe.fermiDistribution(energyRange=e, Temp=g, fermiLevel=fermi_circle)
dis_triangle, dfdE_triangle = SiGe.fermiDistribution(energyRange=e, Temp=g, fermiLevel=fermi_triangle)
dis_square, dfdE_square = SiGe.fermiDistribution(energyRange=e, Temp=g, fermiLevel=fermi_square)
dis_diamond, dfdE_diamond = SiGe.fermiDistribution(energyRange=e, Temp=g, fermiLevel=fermi_diamond)

# tau_p_pb_type_1, tau_p_npb_type_1 = SiGe.tau_p(energyRange=e, alpha=alpha, Dv=2.94, DA=9.5, T=g, vs=sp, D=DoS, rho=rho)
tau_p_pb_type_2, tau_p_npb_type_2 = SiGe.tau_p(energyRange=e, alpha=alpha, Dv=2.94, DA=9.5, T=g, vs=sp, D=dos_nonparabolic, rho=rho)
# tau_p_pb_type_3, tau_p_npb_type_3 = SiGe.tau_p(energyRange=e, alpha=alpha, Dv=2.94, DA=9.5, T=g, vs=sp, D=dos_parabolic, rho=rho)

U_alloy = 0.7
tau_alloy = 1/(x*(1-x)*(3*SiGe.latticeParameter**3*np.pi**3*U_alloy**2/8/thermoelectricProperties.hBar)*m_CB.T**(3/2)*np.sqrt(e)/np.sqrt(2)/np.pi**2/thermoelectricProperties.hBar**3/0.75)*thermoelectricProperties.e2C**(3/2)

LD_circle_non_degenerate = np.sqrt(4*np.pi*SiGe.dielectric*thermoelectricProperties.e0*thermoelectricProperties.kB/thermoelectricProperties.e2C*g/cc_circle) # screening length
LD_triangle_non_degenerate = np.sqrt(4*np.pi*SiGe.dielectric*thermoelectricProperties.e0*thermoelectricProperties.kB/thermoelectricProperties.e2C*g/cc_triangle) # screening length
LD_square_non_degenerate = np.sqrt(4*np.pi*SiGe.dielectric*thermoelectricProperties.e0*thermoelectricProperties.kB/thermoelectricProperties.e2C*g/cc_square) # screening length
LD_diamond_non_degenerate = np.sqrt(4*np.pi*SiGe.dielectric*thermoelectricProperties.e0*thermoelectricProperties.kB/thermoelectricProperties.e2C*g/cc_diamond) # screening length

Nc = 2*(m_CB*thermoelectricProperties.kB*g/thermoelectricProperties.hBar**2/2/np.pi/thermoelectricProperties.e2C)**(3/2)

# np.savetxt("Ef_circle",fermi_circle/g/thermoelectricProperties.kB)
# np.savetxt("Ef_triangle",fermi_triangle/g/thermoelectricProperties.kB)
# np.savetxt("Ef_square",fermi_square/g/thermoelectricProperties.kB)
# np.savetxt("Ef_diamond",fermi_diamond/g/thermoelectricProperties.kB)

fermi_int_circle = np.loadtxt("f_circle", delimiter=',')
fermi_int_triangle= np.loadtxt("f_triangle", delimiter=',')
fermi_int_square = np.loadtxt("f_square", delimiter=',')
fermi_int_diamond = np.loadtxt("f_diamond", delimiter=',')

LD_circle_degenerate = np.sqrt(1/(Nc/SiGe.dielectric/thermoelectricProperties.e0/thermoelectricProperties.kB/g*thermoelectricProperties.e2C*(fermi_int_circle[1]+15*alpha*thermoelectricProperties.kB*g/4*fermi_int_circle[0])))
LD_triangle_degenerate= np.sqrt(1/(Nc/SiGe.dielectric/thermoelectricProperties.e0/thermoelectricProperties.kB/g*thermoelectricProperties.e2C*(fermi_int_triangle[1]+15*alpha*thermoelectricProperties.kB*g/4*fermi_int_triangle[0])))
LD_square_degenerate = np.sqrt(1/(Nc/SiGe.dielectric/thermoelectricProperties.e0/thermoelectricProperties.kB/g*thermoelectricProperties.e2C*(fermi_int_square[1]+15*alpha*thermoelectricProperties.kB*g/4*fermi_int_square[0])))
LD_diamond_degenerate = np.sqrt(1/(Nc/SiGe.dielectric/thermoelectricProperties.e0/thermoelectricProperties.kB/g*thermoelectricProperties.e2C*(fermi_int_diamond[1]+15*alpha*thermoelectricProperties.kB*g/4*fermi_int_diamond[0])))


tau_Screened_Coulomb_circle = SiGe.tau_Screened_Coulomb(energyRange=e, m_c=m_CB, LD = LD_circle_degenerate, N = cc_circle)
tau_Screened_Coulomb_triangle = SiGe.tau_Screened_Coulomb(energyRange=e, m_c=m_CB, LD = LD_triangle_degenerate, N = cc_triangle)
tau_Screened_Coulomb_square = SiGe.tau_Screened_Coulomb(energyRange=e, m_c=m_CB, LD = LD_square_degenerate, N = cc_square)
tau_Screened_Coulomb_diamond = SiGe.tau_Screened_Coulomb(energyRange=e, m_c=m_CB, LD = LD_diamond_degenerate, N = cc_diamond)

tau_Unscreened_Coulomb_circle = SiGe.tau_Unscreened_Coulomb(energyRange=e, m_c=m_CB, N = cc_circle)
tau_Unscreened_Coulomb_triangle = SiGe.tau_Unscreened_Coulomb(energyRange=e, m_c=m_CB, N = cc_triangle)
tau_Unscreened_Coulomb_square = SiGe.tau_Unscreened_Coulomb(energyRange=e, m_c=m_CB, N = cc_square)
tau_Unscreened_Coulomb_diamond = SiGe.tau_Unscreened_Coulomb(energyRange=e, m_c=m_CB, N = cc_diamond)


tau_Strongly_Screened_Coulomb_circle  = SiGe.tau_Strongly_Screened_Coulomb(D=dos_nonparabolic, LD=LD_circle_degenerate, N=cc_circle)
tau_Strongly_Screened_Coulomb_triangle = SiGe.tau_Strongly_Screened_Coulomb(D=dos_nonparabolic, LD=LD_triangle_degenerate , N=cc_triangle)
tau_Strongly_Screened_Coulomb_square = SiGe.tau_Strongly_Screened_Coulomb(D=dos_nonparabolic, LD=LD_square_degenerate , N=cc_square)
tau_Strongly_Screened_Coulomb_diamond = SiGe.tau_Strongly_Screened_Coulomb(D=dos_nonparabolic, LD=LD_diamond_degenerate , N=cc_diamond)



tau_circle  = SiGe.matthiessen(e, 6*tau_p_npb_type_2, 6*tau_alloy,6*tau_Strongly_Screened_Coulomb_circle)
tau_triangle  = SiGe.matthiessen(e, 6*tau_p_npb_type_2, 6*tau_alloy,6*tau_Screened_Coulomb_triangle)
tau_square  = SiGe.matthiessen(e, 6*tau_p_npb_type_2, 6*tau_alloy,6*tau_Strongly_Screened_Coulomb_square)
tau_diamond = SiGe.matthiessen(e, 6*tau_p_npb_type_2,6*tau_alloy,6*tau_Screened_Coulomb_diamond)

Coeff_circle = SiGe.electricalProperties(E=e, DoS=dos_nonparabolic, vg=gVel, Ef=fermi_circle, dfdE=dfdE_circle, Temp=g, tau=tau_circle)
Coeff_triangle = SiGe.electricalProperties(E=e, DoS=dos_nonparabolic, vg=gVel, Ef=fermi_triangle, dfdE=dfdE_triangle, Temp=g, tau=tau_triangle)
Coeff_square = SiGe.electricalProperties(E=e, DoS=dos_nonparabolic, vg=gVel, Ef=fermi_square, dfdE=dfdE_square, Temp=g, tau=tau_square)
Coeff_diamond = SiGe.electricalProperties(E=e, DoS=dos_nonparabolic, vg=gVel, Ef=fermi_diamond, dfdE=dfdE_diamond, Temp=g, tau=tau_diamond)

# np.savetxt("tau_ion_circle",6*tau_ion_circle)
# np.savetxt("tau_ion_triangle",6*tau_ion_triangle)
# np.savetxt("tau_ion_square",6*tau_ion_square)
# np.savetxt("tau_ion_diamond",6*tau_ion_diamond)
# np.savetxt("tau_p_npb_type_2",6*tau_p_npb_type_2)
# np.savetxt("tau_alloy",6*tau_alloy)
# np.savetxt("energy",e)

# exit()


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
ax_4.plot(g[0],Coeff_circle[0]*1e-5, 'None', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_4.plot(g[0],Coeff_triangle[0]*1e-5, 'None', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)
ax_4.plot(g[0],Coeff_square[0]*1e-5, 'None', linestyle='-', color='olive',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='olive',
          markeredgewidth=1)
ax_4.plot(g[0],Coeff_diamond[0]*1e-5, 'None', linestyle='-', color='tan',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='tan',
          markeredgewidth=1)

ax_4.plot(vining_circle[:,0], np.divide(1, vining_circle[:,1]), '-o', linestyle='None', color='maroon',
        markersize=5, linewidth=4,
        markerfacecolor='white',
        markeredgecolor='maroon',
        markeredgewidth=1)

ax_4.plot(vining_triangle[:,0], np.divide(1, vining_triangle[:,1]), '-o', linestyle='None', color='steelblue',
        markersize=5, linewidth=4,
        markerfacecolor='white',
        markeredgecolor='steelblue',
        markeredgewidth=1)

ax_4.plot(vining_square[:,0], np.divide(1, vining_square[:,1]), '-o', linestyle='None', color='olive',
        markersize=5, linewidth=4,
        markerfacecolor='white',
        markeredgecolor='olive',
        markeredgewidth=1)

ax_4.plot(vining_diamond[:,0], np.divide(1, vining_diamond[:,1]), '-o', linestyle='None', color='tan',
        markersize=5, linewidth=4,
        markerfacecolor='white',
        markeredgecolor='tan',
        markeredgewidth=1)

ax_4.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_4.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_4.tick_params(axis="x", labelsize=16)
ax_4.set_ylabel('Conductivity(x10$^5$ S/m)', fontsize=16)
ax_4.tick_params(axis="y", labelsize=16)
# ax_4.set_xlim(right=850)
fig_4.tight_layout()

fig_5 = plt.figure(figsize=(6.5,4.5))
ax_5 = fig_5.add_subplot(111)
ax_5.plot(g[0],Coeff_circle[1]*1e6, 'None', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_5.plot(g[0],Coeff_triangle[1]*1e6, 'None', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)
ax_5.plot(g[0],Coeff_square[1]*1e6, 'None', linestyle='-', color='olive',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='olive',
          markeredgewidth=1)
ax_5.plot(g[0],Coeff_diamond[1]*1e6, 'None', linestyle='-', color='tan',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='tan',
          markeredgewidth=1)

ax_5.plot(vining_seebeck_circle[:,0], -1*vining_seebeck_circle[:,1], '-o', linestyle='None', color='maroon',
          markersize=5, linewidth=4,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_5.plot(vining_seebeck_triangle[:,0], -1*vining_seebeck_triangle[:,1], '-o', linestyle='None', color='steelblue',
          markersize=5, linewidth=4,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)
ax_5.plot(vining_seebeck_square[:,0], -1*vining_seebeck_square[:,1], '-o', linestyle='None', color='olive',
          markersize=5, linewidth=4,
          markerfacecolor='white',
          markeredgecolor='olive',
          markeredgewidth=1)
ax_5.plot(vining_seebeck_diamond[:,0], -1*vining_seebeck_diamond[:,1], '-o', linestyle='None', color='tan',
          markersize=5, linewidth=4,
          markerfacecolor='white',
          markeredgecolor='tan',
          markeredgewidth=1)

ax_5.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax_5.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_5.tick_params(axis="x", labelsize=16)
ax_5.set_ylabel('Seebeck($\mu$V/K)', fontsize=16)
ax_5.tick_params(axis="y", labelsize=16)
# ax_5.set_xlim(right=850)
fig_5.tight_layout()

fig_6 = plt.figure(figsize=(6.5,4.5))
ax_6 = fig_6.add_subplot(111)
# ax_6.plot(g[0],LD_non_degenerate_circle[0]*1e9, 'o', linestyle='-.', color='maroon',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='maroon',
#           markeredgewidth=1)

ax_6.plot(g[0],LD_circle_degenerate[0]*1e9, 'o', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)

ax_6.plot(g[0],LD_circle_non_degenerate[0]*1e9, 'o', linestyle='-.', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='maroon',
          markeredgecolor='maroon',
          markeredgewidth=1)


ax_6.plot(g[0],LD_triangle_degenerate[0]*1e9, 'o', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)
ax_6.plot(g[0],LD_triangle_non_degenerate[0]*1e9, 'o', linestyle='-.', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='steelblue',
          markeredgecolor='steelblue',
          markeredgewidth=1)


ax_6.plot(g[0],LD_square_degenerate[0]*1e9, 'o', linestyle='-', color='olive',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='olive',
          markeredgewidth=1)

ax_6.plot(g[0],LD_square_non_degenerate[0]*1e9, 'o', linestyle='-.', color='olive',
          markersize=6, linewidth=1.5,
          markerfacecolor='olive',
          markeredgecolor='olive',
          markeredgewidth=1)

# ax_6.plot(g[0],LD_diamond_degenerate[0]*1e9, 'o', linestyle='-', color='tan',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='white',
#           markeredgecolor='tan',
#           markeredgewidth=1)
# ax_6.plot(g[0],LD_diamond_non_degenerate[0]*1e9, 'o', linestyle='-', color='tan',
#           markersize=6, linewidth=1.5,
#           markerfacecolor='tan',
#           markeredgecolor='tan',
#           markeredgewidth=1)

ax_6.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_6.set_xlabel('Temperature (K)', fontsize=16, labelpad=10)
ax_6.tick_params(axis="x", labelsize=16)
ax_6.set_ylabel('Debyle length (nm)', fontsize=16, labelpad=10)
ax_6.tick_params(axis="y", labelsize=16)
fig_6.tight_layout()


fig_7 = plt.figure(figsize=(6.5,4.5))
ax_7 = fig_7.add_subplot(111)
ax_7.plot(e[0],np.log10(tau_Unscreened_Coulomb_circle[10])+15, 'None', linestyle='-', color='maroon',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='maroon',
          markeredgewidth=1)
ax_7.plot(e[0],np.log10(tau_Screened_Coulomb_circle[10])+15, 'None', linestyle='-', color='steelblue',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='steelblue',
          markeredgewidth=1)
ax_7.plot(e[0],np.log10(tau_Strongly_Screened_Coulomb_circle[10])+15, 'None', linestyle='-', color='indigo',
          markersize=6, linewidth=1.5,
          markerfacecolor='white',
          markeredgecolor='indigo',
          markeredgewidth=1)

ax_7.yaxis.set_major_formatter(ScalarFormatter())
ax_7.set_xlabel('Energy (eV)', fontsize=16, labelpad=10)
ax_7.tick_params(axis="x", labelsize=16)
ax_7.set_ylabel('Lifetime (log$_{10}$[fs])', fontsize=16, labelpad=10)
ax_7.tick_params(axis="y", labelsize=16)
ax_7.ticklabel_format(axis="y", style="sci", scilimits=None)
fig_7.tight_layout()

plt.show()
exit()
