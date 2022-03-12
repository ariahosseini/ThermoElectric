"""Provide the methods to read and write data files."""

import numpy as np
from os.path import expanduser
from .util import temperature
from scipy.interpolate import InterpolatedUnivariateSpline as spline


def kpoints(path2kpoints: str, delimiter: str = None, skip_rows: int = 0) -> np.ndarray:

    """
    Create a 2D array of temperature sampling

    Parameters
    ----------
    path2kpoints: str
        Path to kpoints file
    delimiter: str
        Default it None for ,
    skip_rows: int
        Number of lines to skip, default is 0

    Returns
    ----------
    wave_points : np.ndarray
        Wave vectors
    """

    wave_points = np.loadtxt(expanduser(path2kpoints), delimiter=delimiter, skiprows=skip_rows)

    return wave_points


def carrier_concentration(path_extrinsic_carrier: str, band_gap: np.ndarray,
                          Ao: float = None, Bo: float = None, Nc: float = None,
                          Nv: float = None, temp: np.ndarray = None) -> np.ndarray:

    """
    This function computes the carrier concentration. The extrinsic carrier concentration is from experiments.
    The following formula is used to compute intrinsic carrier concentration: n = sqrt(Nc*Nv)*exp(-Eg/kB/T/2)
    A good reference book is "Principles of Semiconductor Devices" by Sima Dimitrijev

    Parameters
    ----------
    path_extrinsic_carrier: str
        Path to kpoints file
    band_gap: np.ndarray
        The electronic band gap
    Ao: float
        Experimentally fitted parameter (Nc ~ Ao*T^(3/2))
    Bo: float
        Experimentally fitted parameter (Nv ~ Ao*T^(3/2))
    Nc: float
        The effective densities of states in the conduction band
    Nv: float
        The effective densities of states in the conduction band
    temp: np.ndarray
        Temperature range

    Returns
    ----------
    carrier : np.ndarray
        The total carrier concentration
    """

    k_bolt = 8.617330350e-5  # Boltzmann constant in eV/K

    if temp is None:
        T = temperature()
    else:
        T = temp

    if Ao is None and Nc is None:
        raise Exception("Either Ao or Nc should be defined")
    if Bo is None and Nv is None:
        raise Exception("Either Bo or Nv should be defined")

    if Nc is None:
        Nc = Ao * T ** (3. / 2)
    if Nv is None:
        Nv = Bo * T ** (3. / 2)

    # Extrinsic carrier concentration
    ex_carrier = np.loadtxt(expanduser(path_extrinsic_carrier), delimiter=None, skiprows=0)
    _ex_carrier_concentration = spline(ex_carrier[0, :], ex_carrier[1, :] * 1e6)
    ex_carrier_concentration = _ex_carrier_concentration(T)

    # Intrinsic carrier concentration
    in_carrier = np.sqrt(Nc * Nv) * np.exp(-1 * band_gap / (2 * k_bolt * T))

    # Total carrier concentration
    carrier = in_carrier + abs(ex_carrier_concentration)

    return carrier


def electronBandStructure(path_eigen: str, skip_lines: int, num_bands: int, num_kpoints: int) -> np.ndarray:

    """
    A function to read "EIGENVAL" file

    Parameters
    ----------
    path_eigen: str
        Path to EIGENVAL file
    skip_lines: int
        Number of lines to skip
    num_bands: int
        Number of bands
    num_kpoints: int
        number of wave vectors

    Returns
    ----------
    dispersion : np.ndarray
        Band structure
    """

    with open(expanduser(path_eigen)) as eigen_file:
        for _ in range(skip_lines):
            next(eigen_file)
        block = [[float(_) for _ in line.split()] for line in eigen_file]
    eigen_file.close()

    electron_dispersion = np.arange(1, num_bands + 1)
    k_points = np.array(block[1::num_bands + 2])[:, 0:3]

    for _ in range(num_kpoints):
        disp = []
        for __ in range(num_bands):
            disp = np.append(disp, block[__ + 2 + (num_bands + 2) * _][1])
        electron_dispersion = np.vstack([electron_dispersion, disp])

    dispersion = np.array([k_points[1, np.newaxis], electron_dispersion[1, np.newaxis]])

    return dispersion


def electron_density(path_density: str, header_lines: int, num_dos_points: int,
                     unitcell_volume: float, valley_point: int, energy: np.ndarray) -> np.ndarray:

    """
    A function to read "DOSCAR" file

    Parameters
    ----------
    path_density: str
        Path to DOSCAR file
    header_lines: int
        Number of lines to skip
    num_dos_points: int
        Number of points in DOSCAR
    unitcell_volume: float
        The unit cell volume is in [m]
    valley_point: int
        Where valley is located in DOSCAR
    energy: np.ndarray
        The energy range

    Returns
    ----------
    density : np.ndarray
        Electron density of states
    """

    den_state = np.loadtxt(expanduser(path_density), delimiter=None, skiprows=header_lines, max_rows=num_dos_points)
    valley_energy = den_state[valley_point, 0]
    dos_spline = spline(den_state[valley_point:, 0] - valley_energy, den_state[valley_point:, 1] / unitcell_volume)
    density = dos_spline(energy)

    return density
