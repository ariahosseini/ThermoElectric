import numpy as np
from numpy.linalg import norm
from scipy.interpolate import PchipInterpolator as Interpolator
from scipy.special import jv, j1 as besselj
from scipy.integrate import trapz
from .accum import *

def calculate_phonon_lifetime(energy: np.ndarray, alpha_term: np.ndarray, D_v: float, D_a: float,
                              temperature: np.ndarray, sound_velocity: float, density_of_states: np.ndarray, mass_density: float) -> dict:
    """
    Calculate the electron-phonon scattering rate using the Ravich model.

    Parameters
    ----------
    energy : np.ndarray
        Energy range.
    alpha_term : np.ndarray
        Non-parabolic term.
    D_v : float
        Hole deformation potential.
    D_a : float
        Electron deformation potential.
    temperature : np.ndarray
        Temperature.
    sound_velocity : float
        Sound velocity.
    density_of_states : np.ndarray
        Density of states.
    mass_density : float
        Mass density.

    Returns
    -------
    dict
        Dictionary with 'parabolic_phonon_lifetime' and 'nonparabolic_phonon_lifetime'.
    """
    h_bar = 6.582119e-16  # Reduced Planck constant in eV.s
    k_boltzmann = 8.617330350e-5  # Boltzmann constant in eV/K
    e_to_coulomb = 1.6021765e-19  # e to Coulomb unit conversion

    nonparabolic_term = (1 - ((alpha_term.T * energy) / (1 + 2 * alpha_term.T * energy) * (1 - D_v / D_a))) ** 2 \
                        - 8 / 3 * (alpha_term.T * energy) * (1 + alpha_term.T * energy) / (
                                    1 + 2 * alpha_term.T * energy) ** 2 * (D_v / D_a)

    tau_ph_parabolic = mass_density * sound_velocity ** 2 * h_bar \
          / (np.pi * k_boltzmann * temperature.T * D_a**2) * 1e9 / e_to_coulomb / density_of_states  # Lifetime for parabolic band

    tau_ph_nonparabolic = tau_ph_parabolic / nonparabolic_term  # Lifetime in nonparabolic band

    return {'parabolic_phonon_lifetime': tau_ph_parabolic, 'nonparabolic_phonon_lifetime': tau_ph_nonparabolic}


def calculate_screened_coulomb_lifetime(density_of_states: np.ndarray, screening_length: np.ndarray,
                                         impurity_scattering: np.ndarray, dielectric_constant: float) -> np.ndarray:
    """
    Calculate electron-impurity scattering lifetime in highly doped dielectrics.

    Parameters
    ----------
    density_of_states : np.ndarray
        Density of states.
    screening_length : np.ndarray
        Screening length.
    impurity_scattering : np.ndarray
        Impurity scattering rate.
    dielectric_constant : float
        Dielectric constant.

    Returns
    -------
    np.ndarray
        Electron-impurity lifetime.
    """
    h_bar = 6.582119e-16  # Reduced Planck constant in eV.s
    e_to_coulomb = 1.6021765e-19  # e to Coulomb unit conversion
    permittivity_vacuum = 8.854187817e-12  # Permittivity in vacuum (F/m)

    tau = h_bar / impurity_scattering.T / np.pi / density_of_states / \
          (screening_length.T ** 2 / (4 * np.pi * dielectric_constant * permittivity_vacuum)) ** 2 \
          * 1 / e_to_coulomb ** 2

    return tau


def calculate_screened_coulomb_lifetime_brook_herring(energy: np.ndarray, conduction_mass: np.ndarray, screening_length: np.ndarray,
                                                       impurity_scattering: np.ndarray, dielectric_constant: float) -> np.ndarray:
    """
    Calculate electron-ion scattering lifetime using the Brook-Herring model.

    Parameters
    ----------
    energy : np.ndarray
        Energy range.
    conduction_mass : np.ndarray
        Conduction band effective mass.
    screening_length : np.ndarray
        Screening length.
    impurity_scattering : np.ndarray
        Impurity scattering rate.
    dielectric_constant : float
        Dielectric constant.

    Returns
    -------
    np.ndarray
        Electron-impurity lifetime.
    """
    h_bar = 6.582119e-16  # Reduced Planck constant in eV.s
    e_to_coulomb = 1.6021765e-19  # e to Coulomb unit conversion
    permittivity_vacuum = 8.854187817e-12  # Permittivity in vacuum (F/m)

    gamma = 8 * conduction_mass.T * screening_length.T ** 2 * energy / h_bar ** 2 / e_to_coulomb  # Gamma term

    tau_ = np.log(1 + gamma) - gamma / (1 + gamma)

    tau = 16 * np.pi * np.sqrt(2 * conduction_mass.T) * (4 * np.pi * dielectric_constant * permittivity_vacuum) ** 2 \
          / impurity_scattering.T / tau_ * energy ** (3 / 2) / e_to_coulomb ** (5.0 / 2)

    tau[np.isnan(tau)] = 0

    return tau


def calculate_unscreened_coulomb_lifetime(energy: np.ndarray, conduction_mass: np.ndarray,
                                          impurity_scattering: np.ndarray, dielectric_constant: float) -> np.ndarray:
    """
    Calculate electron-ion scattering lifetime for shallow dopants (no screening effect considered).

    Parameters
    ----------
    energy : np.ndarray
        Energy range.
    conduction_mass : np.ndarray
        Conduction band effective mass.
    impurity_scattering : np.ndarray
        Impurity scattering rate.
    dielectric_constant : float
        Dielectric constant.

    Returns
    -------
    np.ndarray
        Electron-impurity lifetime.
    """
    e_to_coulomb = 1.6021765e-19  # e to Coulomb unit conversion
    permittivity_vacuum = 8.854187817e-12  # Permittivity in vacuum (F/m)

    gamma = 4 * np.pi * (4 * np.pi * dielectric_constant * permittivity_vacuum) * energy / impurity_scattering.T ** (1.0 / 3) / e_to_coulomb  # Gamma term

    tau_ = np.log(1 + gamma ** 2)

    tau = 16 * np.pi * np.sqrt(2 * conduction_mass) * (4 * np.pi * dielectric_constant * permittivity_vacuum) ** 2 \
          / impurity_scattering.T / tau_ * energy ** (3 / 2) / e_to_coulomb ** (5.0 / 2)

    tau[np.isnan(tau)] = 0

    return tau


def calculate_inf_cylinder_scattering(radius: float, k_points: tuple, potential: float, mass_fraction: np.ndarray,
                                       volume_fraction: float, initial_wave_vector: np.ndarray, delta_k_vector: np.ndarray, num_divisions: int) -> tuple:
    """
    Calculate electron scattering rate from a cylindrical potential in an ellipsoid conduction band.
    q = |k'-k|
    Matrix element: M = 4*pi*u0*(1./q.*sin(r_inc*q)-r_inc*cos(r_inc*q))./(q.^2)
    Scattering rate matrix ignoring delta(E-E'): SR = 2*pi/hbar*M.*conj(M)
    E = Ec + hbar^2/2*((kx-k0x)^2/ml+(ky^2+kz^2)/mt)
    Ec = 0
    
    Parameters
    ----------
    radius : float
        Radius of the cylindrical potential.
    k_points : tuple
        Number of k-points in each direction (x, y, z).
    potential : float
        Scattering potential.
    mass_fraction : np.ndarray
        Mass fraction array (m1, m2, m3).
    volume_fraction : float
        Volume fraction.
    initial_wave_vector : np.ndarray
        Initial wave vector (kx, ky, kz).
    delta_k_vector : np.ndarray
        Delta k vector (kx, ky, kz).
    num_divisions : int
        Number of divisions for the ellipse parametrization.

    Returns
    -------
    tuple
        - wave vector magnitude for all k-points.
        - energy values at each k-point.
        - relaxation times at each k-point.
        - number of particles in unit volume.
    """
    h_bar = 6.582119514e-16  # Reduced Planck constant (eV.s)
    eV_to_J = 1.60218e-19    # Convert eV to Joules
    electron_mass = 9.10938356e-31  # Electron rest mass (kg)
    
    # Effective electron mass (kg)
    masses = electron_mass * np.array(mass_fraction)
    
    # Number of particles in unit volume
    num_particles = volume_fraction / (np.pi * radius**2)
    
    # Define k-points mesh
    kx = np.linspace(initial_wave_vector[0], initial_wave_vector[0] + delta_k_vector[0], k_points[0])
    ky = np.linspace(initial_wave_vector[1], initial_wave_vector[1] + delta_k_vector[1], k_points[1])
    kz = np.linspace(initial_wave_vector[2], initial_wave_vector[2] + delta_k_vector[2], k_points[2])
    
    # Create meshgrid and flatten
    kx_grid, ky_grid, kz_grid = np.meshgrid(kx, ky, kz)
    k_point_matrix = np.column_stack((kx_grid.ravel(), ky_grid.ravel(), kz_grid.ravel()))  # Matrix of k-points
    
    # Calculate wave vector magnitude
    k_magnitude = np.linalg.norm(k_point_matrix, axis=1)
    
    # Calculate energy
    E = (h_bar**2 / 2) * ((k_point_matrix[:, 0] - initial_wave_vector[0])**2 / masses[0] +
                          (k_point_matrix[:, 1] - initial_wave_vector[1])**2 / masses[1] +
                          (k_point_matrix[:, 2] - initial_wave_vector[2])**2 / masses[2]) * eV_to_J
    
    # Parametrize the ellipse
    t = np.linspace(0, 2 * np.pi, num_divisions)
    
    # Ellipse semi-major and semi-minor axes
    a = np.sqrt(2 * masses[0] / h_bar**2 * E / eV_to_J)
    b = np.sqrt(2 * masses[1] / h_bar**2 * E / eV_to_J)
    
    # Parametrize the ellipse area (ds)
    ds = np.sqrt((a[:, np.newaxis] * np.sin(t))**2 + (b[:, np.newaxis] * np.cos(t))**2)
    
    # Calculate cos(theta)
    kx_cos_t = a[:, np.newaxis] * np.cos(t)
    ky_cos_t = b[:, np.newaxis] * np.sin(t)
    
    cos_theta = (kx_cos_t * k_point_matrix[:, 0, np.newaxis] +
                 ky_cos_t * k_point_matrix[:, 1, np.newaxis] +
                 k_point_matrix[:, 2, np.newaxis]**2) / \
                np.sqrt(a[:, np.newaxis]**2 * np.cos(t)**2 + 
                        b[:, np.newaxis]**2 * np.sin(t)**2 + 
                        k_point_matrix[:, 2, np.newaxis]**2) / k_magnitude[:, np.newaxis]
    
    # Energy difference (delE)
    delE = h_bar**2 * np.abs((kx_cos_t - initial_wave_vector[0]) / masses[0] +
                            (ky_cos_t - initial_wave_vector[1]) / masses[1] +
                            (k_point_matrix[:, 2, np.newaxis] - initial_wave_vector[2]) / masses[2])
    
    # Calculate q and Bessel function
    qx = k_point_matrix[:, 0, np.newaxis] - a[:, np.newaxis] * np.cos(t)
    qy = k_point_matrix[:, 1, np.newaxis] - b[:, np.newaxis] * np.sin(t)
    qr = np.sqrt(qx**2 + qy**2)
    J = besselj(radius * qr)
    
    # Scattering rate (SR)
    SR = 2 * np.pi / h_bar * potential**2 * (2 * np.pi)**3 * (radius * J / qr)**2
    
    # Function to integrate
    func = SR * (1 - cos_theta) / delE * ds
    
    # Integrate over theta using trapezoidal integration
    integrated_func = trapz(func, t, axis=1)
    
    # Calculate relaxation time (tau)
    relaxation_time = (num_particles / (2 * np.pi)**3 * integrated_func)**-1 * eV_to_J
    
    return k_magnitude, E, relaxation_time, num_particles

def tau_spherical(ro, nk, uo, m_frac, v_frac, ko, del_k, n):
    # constants
    hbar = 6.582119514e-16  # reduced Planck constant (eV.s)
    eV2J = 1.60218e-19      # conversion factor from eV to Joul
    me = 9.10938356e-31     # electron rest mass (Kg)
    m = me * m_frac         # effective mass of electrons (Kg)
    N = 3 * v_frac / (4 * np.pi * ro**3)  # number of particles in unit volume

    # define kpoints
    kx = np.linspace(ko[0], ko[0] + del_k[0], nk[0])
    ky = np.linspace(ko[1], ko[1] + del_k[1], nk[1])
    kz = np.linspace(ko[2], ko[2] + del_k[2], nk[2])
    xk, yk, zk = np.meshgrid(kx, ky, kz)
    kpoint = np.column_stack((xk.flatten(), yk.flatten(), zk.flatten()))
    mag_kpoint = np.linalg.norm(kpoint, axis=1)

    # energy (eV)
    E = (hbar**2 / 2 * (
        (kpoint[:, 0] - ko[0])**2 / m[0] +
        (kpoint[:, 1] - ko[1])**2 / m[1] +
        (kpoint[:, 2] - ko[2])**2 / m[2]
    )) * eV2J

    scattering_rate = np.zeros(E.shape[0])

    for u in range(E.shape[0]):
        Q = np.zeros((2 * n * (n - 1), 3))
        A = np.zeros(2 * n * (n - 1))
        k = 0

        x, y, z = ellipsoid(ko[0], ko[1], ko[2],
                            np.sqrt(2 / (hbar**2 * eV2J) * m[0] * E[u]),
                            np.sqrt(2 / (hbar**2 * eV2J) * m[1] * E[u]),
                            np.sqrt(2 / (hbar**2 * eV2J) * m[2] * E[u]),
                            n)
        
        for j in range(1, n):
            for i in range(2, n + 1):
                S = (np.array([x[i, j], y[i, j], z[i, j]]) +
                      np.array([x[i - 1, j], y[i - 1, j], z[i - 1, j]]) +
                      np.array([x[i - 1, j - 1], y[i - 1, j - 1], z[i - 1, j - 1]]))
                Q[k, :] = S / 3
                a = np.linalg.norm(np.array([x[i, j], y[i, j], z[i, j]]) -
                                   np.array([x[i - 1, j], y[i - 1, j], z[i - 1, j]]))
                b = np.linalg.norm(np.array([x[i - 1, j], y[i - 1, j], z[i - 1, j]]) -
                                   np.array([x[i - 1, j - 1], y[i - 1, j - 1], z[i - 1, j - 1]]))
                c = np.linalg.norm(np.array([x[i - 1, j - 1], y[i - 1, j - 1], z[i - 1, j - 1]]) -
                                   np.array([x[i, j], y[i, j], z[i, j]]))
                s = (a + b + c) / 2
                A[k] = np.sqrt(s * (s - a) * (s - b) * (s - c))
                k += 1

        for j in range(1, n):
            for i in range(1, n):
                S = (np.array([x[i, j - 1], y[i, j - 1], z[i, j - 1]]) +
                      np.array([x[i, j], y[i, j], z[i, j]]) +
                      np.array([x[i - 1, j - 1], y[i - 1, j - 1], z[i - 1, j - 1]]))
                Q[k, :] = S / 3
                a = np.linalg.norm(np.array([x[i, j - 1], y[i, j - 1], z[i, j - 1]]) -
                                   np.array([x[i, j], y[i, j], z[i, j]]))
                b = np.linalg.norm(np.array([x[i, j], y[i, j], z[i, j]]) -
                                   np.array([x[i - 1, j - 1], y[i - 1, j - 1], z[i - 1, j - 1]]))
                c = np.linalg.norm(np.array([x[i - 1, j - 1], y[i - 1, j - 1], z[i - 1, j - 1]]) -
                                   np.array([x[i, j - 1], y[i, j - 1], z[i, j - 1]]))
                s = (a + b + c) / 2
                A[k] = np.sqrt(s * (s - a) * (s - b) * (s - c))
                k += 1

        for i in range(2, n + 1):
          S = (np.array([x[i, 0], y[i, 0], z[i, 0]]) + 
               np.array([x[i - 1, 0], y[i - 1, 0], z[i - 1, 0]]) + 
               np.array([x[i - 1, -2], y[i - 1, -2], z[i - 1, -2]]))
          Q[k, :] = S / 3
          a = np.linalg.norm(np.array([x[i, 0], y[i, 0], z[i, 0]]) - 
                             np.array([x[i - 1, 0], y[i - 1, j], z[i - 1, 0]]))
          b = np.linalg.norm(np.array([x[i - 1, 0], y[i - 1, 0], z[i - 1, 0]]) - 
                             np.array([x[i - 1, -2], y[i - 1, -2], z[i - 1, -2]]))   
          c = np.linalg.norm(np.array([x[i - 1, -2], y[i - 1, -2], z[i - 1, -2]]) - np.array([x[i, 0], y[i, 0], z[i, 0]]))
          s = (a + b + c) / 2
          A[k] = np.sqrt(s * (s - a) * (s - b) * (s - c))
          k += 1

        for i in range(1, n):
          S = (np.array([x[i, -2], y[i, -2], z[i, -2]]) + 
               np.array([x[i, 0], y[i, 0], z[i, 0]]) + 
               np.array([x[i - 1, -2], y[i - 1, -2], z[i - 1, -2]])) 
          Q[k, :] = S / 3 
          a = np.linalg.norm(np.array([x[i, -2], y[i, -2], z[i, -2]]) - 
                             np.array([x[i, 0], y[i, 0], z[i, 0]])) 
          b = np.linalg.norm(np.array([x[i, 0], y[i, 0], z[i, 0]]) - 
                             np.array([x[i - 1, -2], y[i - 1, -2], z[i - 1, -2]])) 
          c = np.linalg.norm(np.array([x[i - 1, -2], y[i - 1, -2], z[i - 1, -2]]) - 
                             np.array([x[i, -2], y[i, -2], z[i, -2]]))
          s = (a + b + c) / 2
          A[k] = np.sqrt(s * (s - a) * (s - b) * (s - c))
          k += 1

        # compute q and cosTheta
        qx = kpoint[u, 0] - Q[:, 0]
        qy = kpoint[u, 1] - Q[:, 1]
        qz = kpoint[u, 2] - Q[:, 2]
        q = np.sqrt(qx**2 + qy**2 + qz**2)
        cosTheta = np.dot(kpoint[u, :], Q.T) / (np.linalg.norm(kpoint[u, :]) * np.sqrt(np.sum(Q**2, axis=1)))

        # matrix element and scattering rate
        M = 4 * np.pi * uo * (1.0 / q * np.sin(ro * q) - ro * np.cos(ro * q)) / q**2
        SR = 2 * np.pi / hbar * M * np.conj(M)
        delE = np.abs(hbar**2 * ((Q[:, 0] - ko[0]) / m[0] +
                                 (Q[:, 1] - ko[1]) / m[1] +
                                 (Q[:, 2] - ko[2]) / m[2]))

        # final scattering rate calculation
        f = SR / delE * (1 - cosTheta)
        scattering_rate[u] = N / (2 * np.pi)**3 * np.sum(f * A)

    tau = 1.0 / scattering_rate * eV2J
    return mag_kpoint, E, tau, N

