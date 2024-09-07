
import numpy as np
from numpy.linalg import norm
from scipy.interpolate import PchipInterpolator as interpolator
from scipy.special import jv
from scipy.special import j1 as besselj  # Bessel function of the first kind
from scipy.integrate import trapz
from .accum import *


def tau_p(energy: np.ndarray, alpha_term: np.ndarray, D_v: float, D_a: float,
          temp: np.ndarray, vel_sound: float, DoS: np.ndarray, rho: float) -> dict:

    """
    Electron-phonon scattering rate using Ravich model

    Parameters
    ----------
    energy: np.ndarray
        Energy range
    alpha_term: np.ndarray
        Non-parabolic term
    D_v: float
        Hole deformation potential
    D_a: float
        Electron deformation potential
    temp: np.ndarray
        Temperature
    vel_sound: float
        Sound velocity
    DoS: np.ndarray
        Density of state
    rho: float
        Mass density

    Returns
    -------
    output: dict
        parabolic and non-parabolic electron-phonon lifetime
    """

    h_bar = 6.582119e-16  # Reduced Planck constant in eV.s
    k_bolt = 8.617330350e-5  # Boltzmann constant in eV/K
    e2C = 1.6021765e-19  # e to Coulomb unit change

    nonparabolic_term = (1 - ((alpha_term.T * energy) / (1 + 2 * alpha_term.T * energy) * (1 - D_v / D_a))) ** 2 \
                        - 8 / 3 * (alpha_term.T * energy) * (1 + alpha_term.T * energy) / (
                                    1 + 2 * alpha_term.T * energy) ** 2 * (D_v / D_a)

    tau_ph_parabolic = rho * vel_sound ** 2 * h_bar \
          / np.pi / k_bolt / temp.T / D_a**2 * 1e9 / e2C / DoS  # Lifetime for parabolic band

    tau_ph_nonparabolic = tau_ph_parabolic / nonparabolic_term  # Lifetime in nonparabolic band

    output = {'parabolic_ph_lifetime': tau_ph_parabolic, 'nonparabolic_ph_lifetime': tau_ph_nonparabolic}

    return output


def tau_strongly_screened_coulomb(DoS: np.ndarray, screen_len: np.ndarray,
                                  n_imp: np.ndarray, dielectric: float) -> np.ndarray:

    """
    Electron-impurity scattering model in highly doped dielectrics

    Note that for highly doped semiconductors, screen length plays a significant role,
    therefor should be computed carefully. Highly suggest to use following matlab file "Fermi.m"
    from: https://www.mathworks.com/matlabcentral/fileexchange/13616-fermi

    If committed to use python, the package "dfint" works with python2
    pip install fdint

    Parameters
    ----------
    DoS: np.ndarray
        Density of states
    screen_len: np.ndarray
        Screening length
    n_imp: np.ndarray
        impurity scattering
    dielectric: float
        Dielectric constant

    Returns
    -------
    tau: np.ndarray
        Electron-impurity lifetime
    """

    h_bar = 6.582119e-16  # Reduced Planck constant in eV.s
    e2C = 1.6021765e-19  # e to Coulomb unit change
    e_o = 8.854187817e-12  # Permittivity in vacuum F/m

    tau = h_bar / n_imp.T / np.pi / DoS / \
          (screen_len.T ** 2 / (4 * np.pi * dielectric * e_o)) ** 2 \
          * 1 / e2C ** 2

    return tau


def tau_screened_coulomb(energy: np.ndarray, mass_c: np.ndarray, screen_len: np.ndarray,
                         n_imp: np.ndarray, dielectric: float) -> np.ndarray:

    """
    Electron-ion scattering rate — Brook-Herring model

    Note that for highly doped semiconductors, screen length plays a significant role,
    therefor should be computed carefully. Highly suggest to use following matlab file "Fermi.m"
    from: https://www.mathworks.com/matlabcentral/fileexchange/13616-fermi

    If committed to use python, the package "dfint" works with python2
    pip install fdint

    Parameters
    ----------
    energy: np.ndarray
        Energy range
    mass_c: np.ndarray
        Conduction band effective mass
    screen_len: np.ndarray
        Screening length
    n_imp: np.ndarray
        impurity scattering
    dielectric: float
        Dielectric constant

    Returns
    -------
    tau: np.ndarray
        Electron-impurity lifetime
    """

    h_bar = 6.582119e-16  # Reduced Planck constant in eV.s
    e2C = 1.6021765e-19  # e to Coulomb unit change
    e_o = 8.854187817e-12  # Permittivity in vacuum F/m

    gamma = 8 * mass_c.T * screen_len.T ** 2 * energy / h_bar ** 2 / e2C  # Gamma term

    tau_ = np.log(1 + gamma) - gamma / (1 + gamma)

    tau = 16 * np.pi * np.sqrt(2 * mass_c.T) * (4 * np.pi * dielectric * e_o) ** 2 \
          / n_imp.T / tau_ * energy ** (3 / 2) / e2C ** (5.0/2)

    tau[np.isnan(tau)] = 0

    return tau


def tau_unscreened_coulomb(energy: np.ndarray, mass_c: np.ndarray,
                           n_imp: np.ndarray, dielectric: float) -> np.ndarray:

    """
    Electron-ion scattering rate for shallow dopants ~10^18 1/cm^3
    (no screening effect is considered)

    Parameters
    ----------
    energy: np.ndarray
        Energy range
    mass_c: np.ndarray
        Conduction band effective mass
    n_imp: np.ndarray
        impurity scattering
    dielectric: float
        Dielectric constant

    Returns
    -------
    tau: np.ndarray
        Electron-impurity lifetime
    """

    e2C = 1.6021765e-19  # e to Coulomb unit change
    e_o = 8.854187817e-12  # Permittivity in vacuum F/m

    gamma = 4 * np.pi * (4 * np.pi * dielectric * e_o) * energy / n_imp.T ** (1.0 / 3) / e2C  # Gamma term

    tau_ = np.log(1 + gamma ** 2)

    tau = 16 * np.pi * np.sqrt(2 * mass_c) * (4 * np.pi * dielectric * e_o) ** 2 \
          / n_imp.T / tau_ * energy ** (3 / 2) / e2C ** (5.0 / 2)

    tau[np.isnan(tau)] = 0

    return tau

def tau_inf_cylinder(ro, nk, uo, m_frac, v_frac, ko, del_k, n):
    """
    electron scattering rate from spherical symmetry potential wall in an ellipsoid conduction band.
    
    parameters:
    ro       : radius of the cylindrical potential
    nk       : number of k-points in each direction (x, y, z)
    uo       : scattering potential
    m_frac   : mass fraction array (m1, m2, m3)
    v_frac   : volume fraction
    ko       : initial wave vector (kx, ky, kz)
    del_k    : delta k vector (kx, ky, kz)
    n        : number of divisions for the ellipse parametrization
    
    returns:
    mag_kpoint : wave vector magnitude for all k-points
    E          : energy values at each k-point
    tau        : relaxation times at each k-point
    N          : number of particles in unit volume
    """
    
    # constants
    hbar = 6.582119514e-16  # Reduced Planck constant (eV.s)
    eV2J = 1.60218e-19      # Convert eV to Joules
    me = 9.10938356e-31     # Electron rest mass (kg)
    
    # Effective electron mass (kg)
    m = me * np.array(m_frac)
    
    # number of particles in unit volume
    N = v_frac / np.pi / ro**2
    
    # define k-points mesh
    kx = np.linspace(ko[0], ko[0] + del_k[0], nk[0])
    ky = np.linspace(ko[1], ko[1] + del_k[1], nk[1])
    kz = np.linspace(ko[2], ko[2] + del_k[2], nk[2])
    
    # create meshgrid and flatten
    xk, yk, zk = np.meshgrid(kx, ky, kz)
    kpoint = np.column_stack((xk.ravel(), yk.ravel(), zk.ravel()))  # Matrix of k-points
    
    # calculate wave vector magnitude
    mag_kpoint = np.linalg.norm(kpoint, axis=1)
    
    # calculate energy
    E = hbar**2 / 2 * ((kpoint[:, 0] - ko[0])**2 / m[0] +
                       (kpoint[:, 1] - ko[1])**2 / m[1] +
                       (kpoint[:, 2] - ko[2])**2 / m[2]) * eV2J
    
    # parametrize the ellipse
    t = np.linspace(0, 2 * np.pi, n)
    
    # ellipse semi-major and semi-minor axes
    a = np.sqrt(2 * m[0] / hbar**2 * E / eV2J)
    b = np.sqrt(2 * m[1] / hbar**2 * E / eV2J)
    
    # parametrize the ellipse area (ds)
    ds = np.sqrt((a[:, np.newaxis] * np.sin(t))**2 + (b[:, np.newaxis] * np.cos(t))**2)
    
    # calculate cos(theta)
    cos_theta = (a[:, np.newaxis] * kpoint[:, 0, np.newaxis] * np.cos(t) +
                 b[:, np.newaxis] * kpoint[:, 1, np.newaxis] * np.sin(t) +
                 kpoint[:, 2, np.newaxis]**2) / \
                np.sqrt(a[:, np.newaxis]**2 * np.cos(t)**2 + 
                        b[:, np.newaxis]**2 * np.sin(t)**2 + 
                        kpoint[:, 2, np.newaxis]**2) / mag_kpoint[:, np.newaxis]
    
    # energy difference (delE)
    delE = hbar**2 * np.abs((a[:, np.newaxis] * np.cos(t) - ko[0]) / m[0] +
                            (b[:, np.newaxis] * np.sin(t) - ko[1]) / m[1] +
                            (kpoint[:, 2, np.newaxis] - ko[2]) / m[2])
    
    # calculate q and Bessel function
    qx = kpoint[:, 0, np.newaxis] - a[:, np.newaxis] * np.cos(t)
    qy = kpoint[:, 1, np.newaxis] - b[:, np.newaxis] * np.sin(t)
    qr = np.sqrt(qx**2 + qy**2)
    J = besselj(ro * qr)
    
    # scattering rate (SR)
    SR = 2 * np.pi / hbar * uo**2 * (2 * np.pi)**3 * (ro * J / qr)**2
    
    # function to integrate
    func = SR * (1 - cos_theta) / delE * ds
    
    # integrate over theta using trapezoidal integration
    Int = trapz(func, t, axis=1)
    
    # calculate relaxation time (tau)
    tau = (N / (2 * np.pi)**3 * Int)**-1 * eV2J
    
    return mag_kpoint, E, tau, N


def ellipsoid(xc, yc, zc, xr, yr, zr, n):
    """
    generates the x, y, z coordinates of an ellipsoid surface.
    
    Parameters:
    xc, yc, zc : float
        The center of the ellipsoid.
    xr, yr, zr : float
        The radii along the x, y, and z axes.
    n : int
        The number of divisions along the grid.

    returns:
    x, y, z : ndarray
        the meshgrid arrays representing the ellipsoid surface
    """
    # create a meshgrid in spherical coordinates (theta, phi)
    u = np.linspace(0, 2 * np.pi, n+1)     # theta from 0 to 2*pi
    v = np.linspace(0, np.pi, n+1)         # phi from 0 to pi

    u, v = np.meshgrid(u, v)

    # parametric equations for the ellipsoid
    x = xc + xr * np.cos(u) * np.sin(v)
    y = yc + yr * np.sin(u) * np.sin(v)
    z = zc + zr * np.cos(v)

    return -x, -y, -z

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

