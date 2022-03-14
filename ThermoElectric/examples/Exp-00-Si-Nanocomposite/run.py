"""

This is an example to show how to use ThermoElectric framework. In this script electrical properties of Si-based nanocomposite with nanoscale SiC spherical particles is studied.

"""
import ThermoElectric as TE
import numpy as np
from numpy.linalg import norm
from scipy.interpolate import PchipInterpolator as interpolator


"""

Data file for pristine, 1% and 5% nanoparticle's valume fraction. Note that the process of desolving P dopants is nonreversible, so the concentration is different when heating up (direction_up), and while cooling down (direction_down)

"""
exp_sic_fra_0pct_heating = np.loadtxt('Exp_Data/ExpData_SiCfrac-0pct_direction-up.txt', delimiter=None, skiprows=1)
exp_sic_frac_1pct_heating = np.loadtxt('Exp_Data/ExpData_SiCfrac-1pct_direction-up.txt', delimiter=None, skiprows=1)
exp_sic_frac_5pct_heating = np.loadtxt('Exp_Data/ExpData_SiCfrac-5pct_direction-up.txt', delimiter=None, skiprows=1)
exp_sic_frac_5pct_cooling = np.loadtxt('Exp_Data/ExpData_SiCfrac-5pct_direction-down.txt', delimiter=None, skiprows=1)

# These variables are [T(K), Nc(1/cm^3)], where Nc is the carreir concentration
carreir_0pct_heating = np.array([exp_sic_fra_0pct_heating[:, 0], exp_sic_fra_0pct_heating[:, -2] * 1e20])
carreir_1pct_heating = np.array([exp_sic_frac_1pct_heating[:, 0],exp_sic_frac_1pct_heating[:, -2] * 1e20])
carreir_5pct_heating = np.array([exp_sic_frac_5pct_heating[:, 0],exp_sic_frac_5pct_heating[:, -2] * 1e20])
carreir_0pct_cooling = np.array([exp_sic_frac_5pct_cooling[:, 0],exp_sic_frac_5pct_cooling[:, -2] * 1e20])

h_bar = 6.582119e-16  # Reduced Planck constant in eV.s
k_bolt = 8.617330350e-5  # Boltzmann constant in eV/K
e2C = 1.6021765e-19  # e to Coulomb unit change
mass_e = 9.109e-31  # Electron rest mass in Kg
eps_o = 8.854187817e-12  # Permittivity in vacuum F/m
ang2meter = 1e-10  # Unit conversion from Angestrom to meter


lattice_parameter = 5.401803661945516e-10  # lattice parameter in [m]
dopant_electric_charge = 1  # Assums that dopands are coompletely ionized
electron_effective_mass=1.08 * mass_e
energy_min = 0.0  # Minimum energy level [eV]
energy_max = 1  # Maximum energy level [eV]
dielectric = 11.7  # Relative dielectricity
num_kpoints = 800  # Number of kpoints
num_bands = 8  # Number of bands
num_qpoints = 200  # Number of q-points in VASP
num_enrg_sample = 4000  # Number of energy points
vfrac = 0.05  # Defects volume fraction
ml = 0.98 * mass_e  # Longitudinal effective mass [Kg]
mt = 0.19 * mass_e  # Transverse effective mass [Kg]
bulk_module = 98  # Bulk module (GPA)
rho = 2329  #Mass density (Kg/m3)

speed_sound = np.sqrt(bulk_module/rho) # speed of sound
# Define reciprocal space vactors
Lv = np.array([[1,1,0],[0,1,1],[1,0,1]])*lattice_parameter/2
a_rp = np.cross(Lv[1],Lv[2])/np.dot(Lv[0],np.cross(Lv[1],Lv[2]))
b_rp = np.cross(Lv[2],Lv[0])/np.dot(Lv[1],np.cross(Lv[2],Lv[0]))
c_rp = np.cross(Lv[0],Lv[1])/np.dot(Lv[2],np.cross(Lv[0],Lv[1]))
RLv = np.array([a_rp, b_rp, c_rp])
engr = TE.energy_range(energy_min = energy_min, energy_max = energy_max, sample_size = num_enrg_sample)
tmpr = TE.temperature(temp_min = 300, temp_max = 1301, del_temp = 50)
electronic_gap = TE.band_gap(Eg_o = 1.17, Ao = 4.73e-4, Bo = 636, temp = tmpr)
nonparabolic_term = 0.5*np.ones((1, len(tmpr[0])))
cc_0pct_heating = TE.carrier_concentration(path_extrinsic_carrier = 'Exp_Data/experimental-carrier-concentration-no-inc.txt', band_gap = electronic_gap,  Ao=5.3e21, Bo=3.5e21, temp = tmpr)
cc_1pct_heating = TE.carrier_concentration(path_extrinsic_carrier = 'Exp_Data/experimental-carrier-concentration-1pct.txt', band_gap = electronic_gap,  Ao=5.3e21, Bo=3.5e21, temp = tmpr)
cc_5pct_heating = TE.carrier_concentration(path_extrinsic_carrier = 'Exp_Data/experimental-carrier-concentration-5pct-direction-up.txt',
 band_gap = electronic_gap,  Ao=5.3e21, Bo=3.5e21, temp = tmpr)
cc_5pct_cooling = TE.carrier_concentration(path_extrinsic_carrier = 'Exp_Data/experimental-carrier-concentration-5pct-direction-down.txt',
 band_gap = electronic_gap,  Ao=5.3e21, Bo=3.5e21, temp = tmpr)
dispersion = TE.band_structure(path_eigen = 'DFT_Data/EIGENVAL', skip_lines = 6, num_bands = num_bands, num_kpoints = num_kpoints)
kp = dispersion['k_points']
band_str = dispersion['electron_dispersion']

kp_rl = 2 * np.pi * kp @ RLv

si_conduc_band = 4
lowest_band_dir_init = 400
lowest_band_dir = band_str[lowest_band_dir_init: num_qpoints + lowest_band_dir_init, si_conduc_band]
kp_mag = norm(kp_rl, axis=1)
min_band = np.argmin(lowest_band_dir, axis=0)
max_band = np.argmax(lowest_band_dir, axis=0)
kp_engr = kp_mag[lowest_band_dir_init + max_band: lowest_band_dir_init + min_band]
energy_kp = band_str[lowest_band_dir_init + max_band: lowest_band_dir_init + min_band, si_conduc_band] - band_str[lowest_band_dir_init + min_band, si_conduc_band]
sort_enrg = np.argsort(energy_kp, axis=0)
grp_velocity = TE.group_velocity(kpoints = kp_engr[sort_enrg], energy_kp = energy_kp[sort_enrg], energy = engr)
e_density = TE.electron_density(path_density = 'DFT_Data/DOSCAR', header_lines =6, unitcell_volume = 2*19.70272e-30, num_dos_points = 2000, valley_point = 1118, energy = engr)

JD_0pct_heating = TE.fermi_level(carrier = cc_0pct_heating, energy = engr, density = e_density, Nc = None, Ao = 5.3e21, temp = tmpr)
JD_1pct_heating = TE.fermi_level(carrier = cc_1pct_heating, energy = engr, density = e_density, Nc = None, Ao = 5.3e21, temp = tmpr)
JD_5pct_heating = TE.fermi_level(carrier = cc_5pct_heating, energy = engr, density = e_density, Nc = None, Ao = 5.3e21, temp = tmpr)
JD_5pct_cooling = TE.fermi_level(carrier = cc_5pct_cooling, energy = engr, density = e_density, Nc = None, Ao = 5.3e21, temp = tmpr)

fermi_0pct_heating = TE.fermi_self_consistent(carrier = cc_0pct_heating, temp = tmpr, energy= engr, density= e_density, fermi_levels = JD_0pct_heating)
fermi_1pct_heating = TE.fermi_self_consistent(carrier = cc_1pct_heating, temp = tmpr, energy= engr, density= e_density, fermi_levels = JD_1pct_heating)
fermi_5pct_heating = TE.fermi_self_consistent(carrier = cc_5pct_heating, temp = tmpr, energy= engr, density= e_density, fermi_levels = JD_5pct_heating)
fermi_5pct_cooling = TE.fermi_self_consistent(carrier = cc_5pct_cooling, temp = tmpr, energy= engr, density= e_density, fermi_levels = JD_5pct_cooling)


fermi_dist_0pct_heating = TE.fermi_distribution(energy = engr, fermi_level = fermi_0pct_heating[1][np.newaxis, :], temp = tmpr)
fermi_dist_1pct_heating = TE.fermi_distribution(energy = engr, fermi_level = fermi_1pct_heating[1][np.newaxis, :], temp = tmpr)
fermi_dist_5pct_heating = TE.fermi_distribution(energy = engr, fermi_level = fermi_5pct_heating[1][np.newaxis, :], temp = tmpr)
fermi_dist_5pct_cooling = TE.fermi_distribution(energy = engr, fermi_level = fermi_5pct_cooling[1][np.newaxis, :], temp = tmpr)

np.savetxt("Matlab_Files/Ef-no-inc.out", fermi_0pct_heating[1] / tmpr / k_bolt)
np.savetxt("Matlab_Files/Ef-inc-1pct.out", fermi_1pct_heating[1] / tmpr / k_bolt)
np.savetxt("Matlab_Files/Ef-inc-5pct-dir-up.out", fermi_5pct_heating[1] / tmpr / k_bolt)
np.savetxt("Matlab_Files/Ef-inc-5pct-dir-down.out", fermi_5pct_cooling[1] / tmpr / k_bolt)

# Non-degenerate screening length
degen_screen_len_0pct = np.sqrt(4 * np.pi * dielectric * eps_o * k_bolt / e2C * tmpr / cc_0pct_heating)
degen_screen_len_1pct_dir_up = np.sqrt(4 * np.pi * dielectric * eps_o * k_bolt / e2C * tmpr / cc_1pct_heating)
degen_screen_len_5pct_dir_up = np.sqrt(4 * np.pi * dielectric * eps_o * k_bolt / e2C * tmpr / cc_5pct_heating)
degen_screen_len_5pct_dir_down = np.sqrt(4 * np.pi * dielectric * eps_o * k_bolt / e2C * tmpr / cc_5pct_cooling)

mass_cond = 0.23 * mass_e * (1 + 5 * nonparabolic_term * k_bolt * tmpr)  # Conduction band effective mass
Nc = 2*(mass_cond * k_bolt * tmpr / h_bar**2 / 2/ np.pi/ e2C)**(3./2)


input("Press Enter to continue...")

fermi_0pct = np.loadtxt("Matlab_Files/fermi-0pct.out", delimiter=',')
fermi_1pct= np.loadtxt("Matlab_Files/fermi-1pct.out", delimiter=',')
fermi_5pct_dir_up = np.loadtxt("Matlab_Files/fermi-5pct-dir-up.out", delimiter=',')
fermi_5pct_dir_down = np.loadtxt("Matlab_Files/fermi-5pct-dir-down.out", delimiter=',')

screen_len_0pct = np.sqrt(1 / (Nc / dielectric / eps_o / k_bolt / tmpr * e2C * (fermi_0pct[1] + 15 * nonparabolic_term * k_bolt * tmpr / 4 *fermi_0pct[0])))
screen_len_1pct = np.sqrt(1 / (Nc / dielectric / eps_o / k_bolt / tmpr * e2C * (fermi_1pct[1] + 15 * nonparabolic_term * k_bolt * tmpr / 4 *fermi_1pct[0])))
screen_len_5pct_dir_up = np.sqrt(1 / (Nc / dielectric / eps_o / k_bolt / tmpr * e2C * (fermi_5pct_dir_up[1] + 15 * nonparabolic_term * k_bolt * tmpr / 4 * fermi_5pct_dir_up[0])))
screen_len_5pct_dir_down = np.sqrt(1 / (Nc / dielectric / eps_o / k_bolt / tmpr * e2C * (fermi_5pct_dir_down[1] + 15 * nonparabolic_term * k_bolt * tmpr / 4 * fermi_5pct_dir_down[0])))

tau_ph = TE.tau_p(energy = engr, alpha_term = nonparabolic_term , D_v = 2.94, D_a = 9.5, temp = tmpr, vel_sound = speed_sound, DoS = e_density, rho = rho)
tau_imp_0pct = TE.tau_strongly_screened_coulomb(DoS = e_density, screen_len = screen_len_0pct, n_imp = cc_0pct_heating, dielectric = dielectric)
tau_imp_1pct = TE.tau_strongly_screened_coulomb(DoS = e_density, screen_len = screen_len_1pct, n_imp = cc_1pct_heating, dielectric = dielectric)
tau_imp_5pct_dir_up = TE.tau_strongly_screened_coulomb(DoS = e_density, screen_len = screen_len_5pct_dir_up, n_imp = cc_5pct_heating, dielectric = dielectric)
tau_imp_5pct_dir_down = TE.tau_strongly_screened_coulomb(DoS = e_density, screen_len = screen_len_5pct_dir_down, n_imp = cc_5pct_cooling, dielectric = dielectric)


_tau_np = np.loadtxt('Lifetime_Data/lifetime_np.out', delimiter=None, skiprows=0)
_engr_np = np.loadtxt('Lifetime_Data/energy_np.out', delimiter=None, skiprows=0)

energy_np, idx_energy_np, rtrn_idx_energy_np = np.unique(_engr_np, return_index = True, return_inverse = True)
tau_np_accum = TE.accum(rtrn_idx_energy_np, _tau_np[1], func = np.mean, dtype = float)
tau_np_interpolator = interpolator(energy_np[1:], tau_np_accum[1:])
tau_np = tau_np_interpolator(engr)

N_gr = 1e25  # Traps concenntration [#/m^2]
_tau_gb = np.loadtxt('Lifetime_Data/lifetime_grain.out', delimiter=None, skiprows=0, max_rows=1)
_engr_gb = np.loadtxt('Lifetime_Data/energy_grain.out', delimiter=None, skiprows=0)
energy_gb, idx_energy_gb, rtrn_idx_energy_gb = np.unique(_engr_gb, return_index=True, return_inverse=True)
tau_gb_accum = TE.accum(rtrn_idx_energy_gb, _tau_gb, func=np.mean, dtype=float)
tau_gb_interpolator  = interpolator(energy_gb[1::], tau_gb_accum[1::])
tau_gb = tau_gb_interpolator(engr) * (fermi_0pct_heating[0] / N_gr)[:, np.newaxis]

num_vally = 6

tau_0pct = TE.matthiessen(engr, num_vally * tau_ph['nonparabolic_ph_lifetime'], num_vally * tau_imp_0pct, tau_gb)
tau_5pct_dir_up = TE.matthiessen(engr, num_vally * tau_ph['nonparabolic_ph_lifetime'], num_vally * tau_imp_5pct_dir_up, tau_gb)
tau_5pct_dir_down = TE.matthiessen(engr, num_vally * tau_ph['nonparabolic_ph_lifetime'], num_vally * tau_imp_5pct_dir_down, tau_gb)
tau_1pct = TE.matthiessen(engr, num_vally * tau_ph['nonparabolic_ph_lifetime'], num_vally * tau_imp_1pct, 5 * tau_gb)  # 5 counts for 1% porosity instead of 5%

prop_0pct = TE.electrical_properties(E = engr, DoS = e_density, vg = grp_velocity, Ef = fermi_0pct_heating[1][np.newaxis, :], dfdE=fermi_dist_0pct_heating[1], temp = tmpr, tau = tau_0pct)

prop_1pct = TE.electrical_properties(E = engr, DoS = e_density, vg = grp_velocity, Ef = fermi_1pct_heating[1][np.newaxis, :], dfdE=fermi_dist_1pct_heating[1], temp = tmpr, tau = tau_1pct)

prop_5pct_dir_up = TE.electrical_properties(E = engr, DoS = e_density, vg = grp_velocity, Ef = fermi_5pct_heating[1][np.newaxis, :], dfdE=fermi_dist_5pct_heating[1], temp = tmpr, tau = tau_0pct)

prop_5pct_dir_down = TE.electrical_properties(E = engr, DoS = e_density, vg = grp_velocity, Ef = fermi_5pct_cooling[1][np.newaxis, :], dfdE=fermi_dist_5pct_cooling[1], temp = tmpr, tau = tau_0pct)


exec(open('./figs.py').read())

exit()

