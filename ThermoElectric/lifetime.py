
import numpy as np
from numpy.linalg import norm
from scipy.interpolate import PchipInterpolator as interpolator
from scipy.special import jv
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

    tau = 16 * np.pi * np.sqrt(2 * mass_c.T) * (4 * np.pi * dielectric * e_o) ** 2 \
          / n_imp.T / tau_ * energy ** (3 / 2) / e2C ** (5.0 / 2)

    tau[np.isnan(tau)] = 0

    return tau


def tau_2d_cylinder(energy: np.ndarray, num_kpoints: np.ndarray, Uo: float, relative_mass: np.ndarray,
                    volume_frac: float, valley: np.ndarray, dk_len: float, ro: np.ndarray,
                    lattice_parameter: float, n_sample=2000) -> np.ndarray:

    """
    A fast algorithm that uses Fermi’s golden rule to compute the energy dependent electron scattering rate
    from cylindrical nano-particles or nano-scale pores infinitely extended perpendicular to the current.

    Parameters
    ----------
    energy: np.ndarray
        Energy range
    num_kpoints: np.ndarray
        Number of kpoints in each direction
    Uo: float
        Barrier height
    relative_mass: np.ndarray
        Relative mass of electron
    volume_frac: float
        Defects volume fraction
    valley: np.ndarray
        Conduction band valley indices
    dk_len: float
        Sample size
    ro: np.ndarray
        Cylinder radius
    lattice_parameter: float
        lattice parameter
    n_sample: int
        Mesh sample size

    Returns
    -------
    tau_cylinder: np.ndarray
        Electron-defect lifetime
    """

    h_bar = 6.582119e-16  # Reduced Planck constant in eV.s
    e2C = 1.6021765e-19  # e to Coulomb unit change
    mass_e = 9.109e-31  # Electron rest mass in Kg

    m_eff = np.array(relative_mass) * mass_e  # Electron conduction effective mass
    ko = 2 * np.pi / lattice_parameter * np.array(valley)
    del_k = 2 * np.pi / lattice_parameter * dk_len * np.array([1, 1, 1])
    N = volume_frac / np.pi / ro ** 2  # Number density

    kx = np.linspace(ko[0], ko[0] + del_k[0], num_kpoints[0], endpoint=True)  # kpoints mesh
    ky = np.linspace(ko[1], ko[1] + del_k[1], num_kpoints[1], endpoint=True)  # kpoints mesh
    kz = np.linspace(ko[2], ko[2] + del_k[2], num_kpoints[2], endpoint=True)  # kpoints mesh
    [xk, yk, zk] = np.meshgrid(kx, ky, kz)
    xk_ = np.reshape(xk, -1)
    yk_ = np.reshape(yk, -1)
    zk_ = np.reshape(zk, -1)
    kpoint = np.array([xk_, yk_, zk_])  # kpoints mesh sampling
    mag_kpoint = norm(kpoint, axis=0)

    E = h_bar ** 2 / 2 * \
        ((kpoint[0, :] - ko[0]) ** 2 / m_eff[0] +
         (kpoint[1, :] - ko[1]) ** 2 / m_eff[1] +
         (kpoint[2, :] - ko[2]) ** 2 / m_eff[2]) * e2C

    t = np.linspace(0, 2 * np.pi, n_sample)

    a = np.expand_dims(np.sqrt(2 * m_eff[1] / h_bar ** 2 * E / e2C), axis=0)
    b = np.expand_dims(np.sqrt(2 * m_eff[2] / h_bar ** 2 * E / e2C), axis=0)
    ds = np.sqrt((a.T * np.sin(t)) ** 2 + (b.T * np.cos(t)) ** 2)

    cos_theta = ((a * kpoint[0]).T * np.cos(t) + (b * kpoint[1]).T * np.sin(t) +
                 np.expand_dims(kpoint[2] ** 2, axis=1)) / \
                np.sqrt(a.T ** 2 * np.cos(t) ** 2 + b.T ** 2 * np.sin(t) ** 2 +
                        np.expand_dims(kpoint[2] ** 2, axis=1)) / np.expand_dims(mag_kpoint, axis=1)

    delE = h_bar ** 2 * \
           np.abs((a.T * np.cos(t) - ko[0]) / m_eff[0] +
                  (b.T * np.sin(t) - ko[1]) / m_eff[1] + (
                              np.expand_dims(kpoint[2] ** 2, axis=1) - ko[2] / m_eff[2]))

    # q_points
    qx = np.expand_dims(kpoint[0], axis=1) - a.T * np.cos(t)
    qy = np.expand_dims(kpoint[1], axis=1) - b.T * np.sin(t)
    qr = np.sqrt(qx ** 2 + qy ** 2)

    tau = np.empty((len(ro), len(E)))

    for r_idx in np.arange(len(ro)):
        J = jv(1, ro[r_idx] * qr)  # Bessel func.
        SR = 2 * np.pi / h_bar * Uo ** 2 * (2 * np.pi) ** 3 * (
                    ro[r_idx] * J / qr) ** 2  # Scattering rate
        f = SR * (1 - cos_theta) / delE * ds
        int_ = np.trapz(f, t, axis=1)
        tau[r_idx] = 1 / (N[r_idx] / (2 * np.pi) ** 3 * int_) * e2C

    Ec, indices, return_indices = np.unique(E, return_index=True, return_inverse=True)

    tau_c = np.empty((len(ro), len(indices)))
    tau_cylinder = np.empty((len(ro), len(energy[0])))

    for r_idx in np.arange(len(ro)):
        tau_c[r_idx] = accum(return_indices, tau[r_idx], func=np.mean, dtype=float)

    # Map lifetime to desired energy range
    for tau_idx in np.arange(len(tau_c)):
        ESpline = interpolator(Ec[30:], tau_c[tau_idx, 30:])
        tau_cylinder[tau_idx] = ESpline(energy)

    return tau_cylinder


def tau3D_spherical(num_kpoints: np.ndarray, Uo: float, relative_mass: np.ndarray,
                    volume_frac: float, valley: np.ndarray, dk_len: float, ro: np.ndarray,
                    lattice_parameter: float, n_sample=32) -> np.ndarray:

    """
    A fast algorithm that uses Fermi’s golden rule to compute the energy dependent electron scattering rate
    from spherical nano-particles or nano-scale pores.

    Parameters
    ----------
    num_kpoints: np.ndarray
        Number of kpoints in each direction
    Uo: float
        Barrier height
    relative_mass: np.ndarray
        Relative mass of electron
    volume_frac: float
        Defects volume fraction
    valley: np.ndarray
        Conduction band valley indices
    dk_len: float
        Sample size
    ro: np.ndarray
        Cylinder radius
    lattice_parameter: float
        lattice parameter
    n_sample: int
        Mesh sample size

    Returns
    -------
    tau: np.ndarray
        Electron-defect lifetime
    """

    h_bar = 6.582119e-16  # Reduced Planck constant in eV.s
    e2C = 1.6021765e-19  # e to Coulomb unit change
    mass_e = 9.109e-31  # Electron rest mass in Kg

    m_eff = np.array(relative_mass) * mass_e  # Electron conduction band effective mass
    ko = 2 * np.pi / lattice_parameter * np.array(valley)
    del_k = 2 * np.pi / lattice_parameter * dk_len * np.array([1, 1, 1])

    N = 3 * volume_frac / 4 / np.pi / ro ** 3  # Number density of defects

    kx = np.linspace(ko[0], ko[0] + del_k[0], num_kpoints[0], endpoint=True)  # kpoints mesh
    ky = np.linspace(ko[1], ko[1] + del_k[1], num_kpoints[1], endpoint=True)  # kpoints mesh
    kz = np.linspace(ko[2], ko[2] + del_k[2], num_kpoints[2], endpoint=True)  # kpoints mesh
    [xk, yk, zk] = np.meshgrid(kx, ky, kz)
    xk_ = np.reshape(xk, -1)
    yk_ = np.reshape(yk, -1)
    zk_ = np.reshape(zk, -1)

    kpoint = np.array([xk_, yk_, zk_])  # kpoint mesh sampling

    # Energy levels in ellipsoidal band structure
    E = h_bar ** 2 / 2 * \
        ((kpoint[0, :] - ko[0]) ** 2 / m_eff[0] +
         (kpoint[1, :] - ko[1]) ** 2 / m_eff[1] +
         (kpoint[2, :] - ko[2]) ** 2 / m_eff[2]) * e2C

    scattering_rate = np.zeros((len(ro), len(E)))

    nu = np.linspace(0, np.pi, n_sample)
    z_ = -1 * np.cos(nu)

    r = np.sqrt(1.0 - z_ ** 2)[:, None]
    theta = np.linspace(0, 2 * np.pi, n_sample)[None, :]

    x_ = r * np.cos(theta)
    y_ = r * np.sin(theta)

    for u in np.arange(len(E)):

        Q = np.zeros((2 * (n_sample - 2) * (n_sample - 1), 3))
        A = np.zeros((2 * (n_sample - 2) * (n_sample - 1), 1))
        k = 0
        a_axis = np.sqrt(2 / (h_bar ** 2 * e2C) * m_eff[0] * E[u])
        b_axis = np.sqrt(2 / (h_bar ** 2 * e2C) * m_eff[1] * E[u])
        c_axis = np.sqrt(2 / (h_bar ** 2 * e2C) * m_eff[2] * E[u])

        y = -1 * b_axis * y_ + ko[1]
        x = -1 * a_axis * x_ + ko[0]
        Z_ = c_axis * z_ + ko[2]
        z = np.tile(Z_[:, None], (1, n_sample))
        for j in np.arange(1, n_sample - 1):
            for i in np.arange(2, n_sample):

                S = np.array(np.array([x[i, j], y[i, j], z[i, j]]) +
                             np.array([x[i - 1, j], y[i - 1, j], z[i - 1, j]]) +
                             np.array([x[i - 1, j - 1], y[i - 1, j - 1], z[i - 1, j - 1]]))

                Q[k] = S / 3

                a = norm(np.array([x[i, j], y[i, j], z[i, j]]) - np.array([x[i - 1, j], y[i - 1, j], z[i - 1, j]]))

                b = norm(np.array([x[i - 1, j], y[i - 1, j], z[i - 1, j]]) -
                         np.array([x[i - 1, j - 1], y[i - 1, j - 1], z[i - 1, j - 1]]))

                c = norm(np.array([x[i - 1, j - 1], y[i - 1, j - 1], z[i - 1, j - 1]])
                         - np.array([x[i, j], y[i, j], z[i, j]]))

                s = a + b + c
                s = s / 2
                A[k] = np.sqrt(s * (s - a) * (s - b) * (s - c))  # Surface area of the triangular mesh elements
                k += 1
        for j in np.arange(1, n_sample - 1):
            for i in np.arange(1, n_sample - 1):

                S = np.array([x[i, j - 1], y[i, j - 1], z[i, j - 1]]) + \
                    np.array([x[i, j], y[i, j], z[i, j]]) + \
                    np.array([x[i - 1, j - 1], y[i - 1, j - 1], z[i - 1, j - 1]])

                Q[k] = S / 3

                a = norm(np.array([x[i, j - 1], y[i, j - 1], z[i, j - 1]]) -
                         np.array([x[i, j], y[i, j], z[i, j]]))

                b = norm(np.array([x[i, j], y[i, j], z[i, j]]) -
                         np.array([x[i - 1, j - 1], y[i - 1, j - 1], z[i - 1, j - 1]]))

                c = norm(np.array([x[i - 1, j - 1], y[i - 1, j - 1], z[i - 1, j - 1]]) -
                         np.array([x[i, j - 1], y[i, j - 1], z[i, j - 1]]))

                s = a + b + c
                s = s / 2
                A[k] = np.sqrt(s * (s - a) * (s - b) * (s - c))

                k += 1

        for i in np.arange(2, n_sample):

            S = np.array([x[i, 0], y[i, 0], z[i, 0]]) + \
                np.array([x[i - 1, 0], y[i - 1, 0], z[i - 1, 0]]) + \
                np.array([x[i - 1, -2], y[i - 1, -2], z[i - 1, -2]])

            Q[k] = S / 3

            a = norm(np.array([x[i, 0], y[i, 0], z[i, 0]]) -
                     np.array([x[i - 1, 0], y[i - 1, 0], z[i - 1, 0]]))

            b = norm(np.array([x[i - 1, 0], y[i - 1, 0], z[i - 1, 0]]) -
                     np.array([x[i - 1, -2], y[i - 1, -2], z[i - 1, -2]]))

            c = norm(np.array([x[i - 1, -2], y[i - 1, -2], z[i - 1, -2]]) -
                     np.array([x[i, 0], y[i, 0], z[i, 0]]))

            s = a + b + c
            s = s / 2
            A[k] = np.sqrt(s * (s - a) * (s - b) * (s - c))

            k += 1

        for i in np.arange(1, n_sample - 1):

            S = np.array([x[i, -2], y[i, -2], z[i, -2]]) + \
                np.array([x[i, 0], y[i, 0], z[i, 0]]) + \
                np.array([x[i - 1, -2], y[i - 1, -2], z[i - 1, -2]])

            Q[k] = S / 3

            a = norm(np.array([x[i, -2], y[i, -2], z[i, -2]]) - np.array([x[i, 0], y[i, 0], z[i, 0]]))

            b = norm(np.array([x[i, 0], y[i, 0], z[i, 0]]) -
                     np.array([x[i - 1, -2], y[i - 1, -2], z[i - 1, -2]]))

            c = norm(np.array([x[i - 1, -2], y[i - 1, -2], z[i - 1, -2]]) -
                     np.array([x[i, -2], y[i, -2], z[i, -2]]))

            s = a + b + c
            s = s / 2
            A[k] = np.sqrt(s * (s - a) * (s - b) * (s - c))

            k += 1

        qx = kpoint[0, u] - Q[:, 0]
        qy = kpoint[1, u] - Q[:, 1]
        qz = kpoint[2, u] - Q[:, 2]
        q = np.sqrt(qx ** 2 + qy ** 2 + qz ** 2)

        cos_theta = np.matmul(kpoint[:, u][None, :], Q.T) / norm(kpoint[:, u]) / np.sqrt(np.sum(Q ** 2, axis=1))

        delE = np.abs(h_bar ** 2 * (
                    (Q[:, 0] - ko[0]) / m_eff[0] + (Q[:, 1] - ko[1]) / m_eff[1] + (Q[:, 2] - ko[2]) / m_eff[2]))

        for ro_idx in np.arange(len(ro)):

            M = 4 * np.pi * Uo * (1 / q * np.sin(ro[ro_idx] * q) - ro[ro_idx] * np.cos(ro[ro_idx] * q)) / (
                        q ** 2)  # Matrix element

            SR = 2 * np.pi / h_bar * M * np.conj(M)  # Scattering rate

            f = SR / delE * (1 - cos_theta)

            scattering_rate[ro_idx, u] = N[ro_idx] / (2 * np.pi) ** 3 * np.sum(f * A.T)

    tau = 1/scattering_rate

    return tau
