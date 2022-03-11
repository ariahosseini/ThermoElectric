"""Provide the methods to read and write data files."""

import numpy as np
from os.path import expanduser


def kpoints(path2kpoints, delimiter=None, skip_rows=0) -> np.ndarray:

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
