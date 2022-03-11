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
                          Nv: float = None, temp: np.ndarray = None):

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
    ex_carrier = np.loadtxt(expanduser(path_extrinsic_carrier), delimiter=None,skiprows=0)
    _ex_carrier_concentration = spline(ex_carrier[0, :], ex_carrier[1, :] * 1e6)
    ex_carrier_concentration = _ex_carrier_concentration(T)

    # Intrinsic carrier concentration
    in_carrier = np.sqrt(Nc * Nv) * np.exp(-1 * band_gap / (2 * k_bolt * T))

    # Total carrier concentration
    carrier = in_carrier + abs(ex_carrier_concentration)

    return carrier


    def electronBandStructure(self, path2eigenval, skipLines):

        # This function ead EIGENVAL file from VASP

        with open(expanduser(path2eigenval)) as eigenvalFile:
            for _ in range(skipLines):
                next(eigenvalFile)
            block = [[float(_) for _ in line.split()] for line in eigenvalFile]
        eigenvalFile.close()

        electronDispersian = [range(1, self.numBands + 1)]  # First line is atoms id

        kpoints = np.asarray(block[1::self.numBands + 2])[:, 0:3]

        for _ in range(self.numKpoints):
            binary2Darray = []
            for __ in range(self.numBands):
                binary2Darray = np.append(binary2Darray, block[__ + 2 + (self.numBands + 2) * _][1])
            electronDispersian = np.vstack([electronDispersian, binary2Darray]) # Energy levels

        dispersian = [kpoints, electronDispersian]

        return dispersian # The array size is [(number of bands + 1) by (number of kpoints)]