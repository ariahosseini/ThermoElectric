import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as Spline
from scipy.interpolate import PchipInterpolator as Interpolator
from .util import *
from .accum import *

def band_gap(Eg_0: float, A: float, B: float, temp: np.ndarray = None) -> np.ndarray:
    """
    Estimate the temperature-dependent band gap using the relation:
    Eg(T) = Eg(T=0) - A * T^2 / (T + B)
    
    Reference: "Properties of Advanced Semiconductor Materials" by Michael E. Levinshtein.

    Parameters
    ----------
    Eg_0: float
        Band gap at 0 Kelvin.
    A: float
        Material-dependent constant.
    B: float
        Material-dependent constant.
    temp: np.ndarray, optional
        Temperature range (default: 300K to 1300K with steps of 100K).

    Returns
    -------
    np.ndarray
        Temperature-dependent band gap.
    """
    T = temp if temp is not None else temperature()
    return Eg_0 - A * (T ** 2) / (T + B)


def analytical_dos(energy_range: np.ndarray, m_eff: float, alpha: np.ndarray) -> dict:
    """
    Compute the electron density of states (DoS) for parabolic and non-parabolic bands.

    Parameters
    ----------
    energy_range: np.ndarray
        Range of electron energies.
    m_eff: float
        Effective mass of the electron.
    alpha: np.ndarray
        Non-parabolicity parameter.

    Returns
    -------
    dict
        DoS for both parabolic and non-parabolic cases.
    """
    h_bar = 6.582119e-16  # Reduced Planck constant in eV.s
    e2C = 1.6021765e-19  # e to Coulomb unit conversion

    # Non-parabolic and parabolic DoS calculation
    dos_nonparabolic = (1 / np.pi ** 2) * np.sqrt(2 * energy_range * (1 + energy_range * alpha.T)) * \
                       (m_eff ** (3/2)) / (h_bar ** 3) * (1 + 2 * energy_range * alpha.T) / e2C ** (3/2)

    dos_parabolic = np.sqrt(energy_range) / (np.pi ** 2) * np.sqrt(2) / (h_bar ** 3) * \
                    m_eff ** (3 / 2) / e2C ** (3 / 2)

    return {'DoS_nonparabolic': dos_nonparabolic, 'DoS_parabolic': dos_parabolic}


def fermi_level(carrier: np.ndarray, energy: np.ndarray, density: np.ndarray, Nc: float = None,
                A: float = None, temp: np.ndarray = None) -> np.ndarray:
    """
    Estimate the Fermi level using Joyce-Dixon approximation.

    Parameters
    ----------
    carrier: np.ndarray
        Carrier concentration.
    energy: np.ndarray
        Energy levels.
    density: np.ndarray
        Density of states.
    Nc: float, optional
        Effective density of states in the conduction band.
    A: float, optional
        Experimentally fitted parameter (Nc ~ A * T^(3/2)).
    temp: np.ndarray, optional
        Temperature range.

    Returns
    -------
    np.ndarray
        Carrier concentration and Fermi level at each temperature.
    """
    k_B = 8.617330350e-5  # Boltzmann constant in eV/K
    T = temp if temp is not None else temperature()
    
    if Nc is None and A is None:
        raise ValueError("Either Nc or A must be provided.")
    
    Nc = Nc if Nc is not None else A * T ** (3 / 2)

    # Joyce-Dixon approximation
    JD_CC = np.log(carrier / Nc) + carrier / Nc / np.sqrt(8) - (3 / 16 - np.sqrt(3) / 9) * (carrier / Nc) ** 2
    fermi_energy = k_B * T * JD_CC
    f, _ = fermi_distribution(energy, fermi_energy, temp=T)
    
    carrier_concentration = np.trapz(density * f, energy, axis=1)
    
    return np.stack((carrier_concentration, fermi_energy[0]))


def fermi_self_consistent(carrier: np.ndarray, temp: np.ndarray, energy: np.ndarray,
                          density: np.ndarray, fermi_guess: np.ndarray) -> np.ndarray:
    """
    Perform self-consistent calculation of the Fermi level using a given carrier concentration.

    Parameters
    ----------
    carrier: np.ndarray
        Carrier concentration.
    energy: np.ndarray
        Energy levels.
    density: np.ndarray
        Density of states.
    fermi_guess: np.ndarray
        Initial Fermi level guess from Joyce-Dixon approximation.
    temp: np.ndarray
        Temperature range.

    Returns
    -------
    np.ndarray
        Carrier concentration and Fermi level.
    """
    fermi = np.linspace(fermi_guess[1] - 0.4, fermi_guess[1] + 0.2, 4000).T
    result = np.empty((np.shape(temp)[1], np.shape(fermi)[1]))

    for j, t in enumerate(temp[0]):
        for i, f in enumerate(fermi):
            f_dist, _ = fermi_distribution(energy, np.array([f]), np.array([t]))
            result[j, i] = np.trapz(density * f_dist, energy, axis=1)

    diff = np.tile(carrier.T, (1, np.shape(fermi)[1])) - abs(result)
    min_idx = np.argmin(np.abs(diff), axis=1)

    fermi_levels = fermi[min_idx]
    n = np.array([result[idx] for idx in min_idx])

    return np.stack((n, fermi_levels))


def group_velocity(kpoints: np.ndarray, energy_kp: np.ndarray, energy: np.ndarray) -> np.ndarray:
    """
    Compute group velocity from band structure using DFT data.

    Parameters
    ----------
    kpoints: np.ndarray
        Wave vectors.
    energy_kp: np.ndarray
        Energy for each wave vector.
    energy: np.ndarray
        Energy range.

    Returns
    -------
    np.ndarray
        Group velocity.
    """
    h_bar = 6.582119e-16  # Reduced Planck constant in eV.s

    dE = np.diff(energy_kp, prepend=energy_kp[0], append=energy_kp[-1])
    dk = np.diff(kpoints, prepend=kpoints[0], append=kpoints[-1])
    dEdk = dE / dk
    dEdk[0] = (energy_kp[1] - energy_kp[0]) / (kpoints[1] - kpoints[0])
    dEdk[-1] = (energy_kp[-1] - energy_kp[-2]) / (kpoints[-1] - kpoints[-2])

    velocity_spline = Spline(energy_kp, dEdk)
    velocity_function = velocity_spline(energy)

    return velocity_function / h_bar


def analytical_group_velocity(energy: np.ndarray, lattice_param: np.ndarray, num_kpoints: np.ndarray,
                              m_eff: np.ndarray, valley: np.ndarray, dk_len: np.ndarray,
                              alpha_term: np.ndarray) -> np.ndarray:
    """
    Approximate group velocity for conduction bands near the band edge.

    Parameters
    ----------
    energy: np.ndarray
        Energy range.
    lattice_param: np.ndarray
        Lattice parameter.
    num_kpoints: np.ndarray
        Number of k-points.
    m_eff: np.ndarray
        Effective mass along axes.
    valley: np.ndarray
        Conduction band valley.
    dk_len: np.ndarray
        Wave vector magnitude.
    alpha_term: np.ndarray
        Non-parabolicity term.

    Returns
    -------
    np.ndarray
        Group velocity.
    """
    h_bar = 6.582119e-16  # Reduced Planck constant in eV.s
    e2C = 1.6021765e-19  # e to Coulomb conversion

    k0 = 2 * np.pi / lattice_param * valley
    del_k = 2 * np.pi / lattice_param * dk_len * np.ones(3)
    
    kx = np.linspace(k0[0], k0[0] + del_k[0], num_kpoints[0])
    ky = np.linspace(k0[1], k0[1] + del_k[1], num_kpoints[1])
    kz = np.linspace(k0[2], k0[2] + del_k[2], num_kpoints[2])
    
    kpoints = np.meshgrid(kx, ky, kz)
    kpoints = np.vstack([k.flatten() for k in kpoints])

    m_eff_cond = 3 / np.sum(1 / m_eff)

    E = h_bar**2 / (2 * m_eff_cond) * np.sum([(kpoints[i] - k0[i]) ** 2 for i in range(3)], axis=0) * e2C
    vel = h_bar * np.sqrt(np.sum([(kpoints[i] - k0[i]) ** 2 for i in range(3)], axis=0)) / m_eff_cond * e2C

    gv_spline = Spline(E, vel)
    gv_function = gv_spline(energy)

    return gv_function
