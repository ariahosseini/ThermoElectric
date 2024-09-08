import numpy as np
from os.path import expanduser
from .util import temperature
from scipy.interpolate import InterpolatedUnivariateSpline as Spline


def kpoints(path2kpoints: str, delimiter: str = None, skip_rows: int = 0) -> np.ndarray:
    """
    Load a 2D array of k-points (wave vectors) from a file.

    Parameters
    ----------
    path2kpoints : str
        Path to the k-points file.
    delimiter : str, optional
        Delimiter for the file (default is None for comma).
    skip_rows : int, optional
        Number of rows to skip at the start of the file (default is 0).

    Returns
    -------
    np.ndarray
        Array of k-points (wave vectors).
    """
    wave_points = np.loadtxt(expanduser(path2kpoints), delimiter=delimiter, skiprows=skip_rows)
    return wave_points


def carrier_concentration(path_extrinsic_carrier: str, band_gap: np.ndarray, 
                          Ao: float = None, Bo: float = None, Nc: float = None, 
                          Nv: float = None, temp: np.ndarray = None) -> np.ndarray:
    """
    Calculate total carrier concentration from intrinsic and extrinsic carriers.

    Uses the relation: n = sqrt(Nc * Nv) * exp(-Eg / (2 * kB * T)) for intrinsic carrier concentration.

    Parameters
    ----------
    path_extrinsic_carrier : str
        Path to the file with extrinsic carrier concentration data.
    band_gap : np.ndarray
        Band gap values.
    Ao : float, optional
        Fitting parameter for Nc (default is None).
    Bo : float, optional
        Fitting parameter for Nv (default is None).
    Nc : float, optional
        Effective density of states in the conduction band (default is None).
    Nv : float, optional
        Effective density of states in the valence band (default is None).
    temp : np.ndarray, optional
        Temperature array (default is None, uses default temperature range from utility function).

    Returns
    -------
    np.ndarray
        Total carrier concentration.
    """
    k_bolt = 8.617330350e-5  # Boltzmann constant in eV/K
    T = temp if temp is not None else temperature()

    if Nc is None:
        Nc = Ao * T ** (3 / 2)
    if Nv is None:
        Nv = Bo * T ** (3 / 2)

    # Load extrinsic carrier concentration data
    ex_carrier_data = np.loadtxt(expanduser(path_extrinsic_carrier), delimiter=None, skiprows=0)
    ex_carrier_spline = Spline(ex_carrier_data[:, 0], ex_carrier_data[:, 1] * 1e6)
    ex_carrier_concentration = ex_carrier_spline(T)

    # Compute intrinsic carrier concentration
    in_carrier = np.sqrt(Nc * Nv) * np.exp(-band_gap / (2 * k_bolt * T))

    # Total carrier concentration
    carrier = in_carrier + np.abs(ex_carrier_concentration)
    
    return carrier


def band_structure(path_eigen: str, skip_lines: int, num_bands: int, num_kpoints: int) -> dict:
    """
    Load the electronic band structure from an "EIGENVAL" file.

    Parameters
    ----------
    path_eigen : str
        Path to the "EIGENVAL" file.
    skip_lines : int
        Number of lines to skip in the file.
    num_bands : int
        Number of bands.
    num_kpoints : int
        Number of k-points (wave vectors).

    Returns
    -------
    dict
        Dictionary containing k-points and electron dispersion data.
    """
    with open(expanduser(path_eigen)) as eigen_file:
        for _ in range(skip_lines):
            next(eigen_file)
        block = [[float(value) for value in line.split()] for line in eigen_file]
    
    electron_dispersion = np.arange(1, num_bands + 1)
    k_points = np.array(block[1::num_bands + 2])[:, 0:3]

    for idx in range(num_kpoints):
        dispersion = []
        for j in range(num_bands):
            dispersion.append(block[j + 2 + (num_bands + 2) * idx][1])
        electron_dispersion = np.vstack([electron_dispersion, dispersion])

    return {'k_points': k_points, 'electron_dispersion': np.delete(electron_dispersion, 0, axis=0)}


def electron_density(path_density: str, header_lines: int, num_dos_points: int, 
                     unitcell_volume: float, valley_point: int, energy: np.ndarray) -> np.ndarray:
    """
    Calculate the electron density of states (DoS) from a "DOSCAR" file.

    Parameters
    ----------
    path_density : str
        Path to the "DOSCAR" file.
    header_lines : int
        Number of header lines to skip in the file.
    num_dos_points : int
        Number of DoS points in the file.
    unitcell_volume : float
        Volume of the unit cell (in cubic meters).
    valley_point : int
        Index of the valley in the DoS data.
    energy : np.ndarray
        Array of energy values.

    Returns
    -------
    np.ndarray
        Electron density of states interpolated over the energy range.
    """
    dos_data = np.loadtxt(expanduser(path_density), delimiter=None, skiprows=header_lines, max_rows=num_dos_points)
    valley_energy = dos_data[valley_point, 0]
    
    # Interpolate the density of states relative to the valley point
    dos_spline = Spline(dos_data[valley_point:, 0] - valley_energy, dos_data[valley_point:, 1] / unitcell_volume)
    density = dos_spline(energy)

    return density
