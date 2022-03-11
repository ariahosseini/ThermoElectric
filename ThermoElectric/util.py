"""Provide the utility methods."""

import numpy as np


def energy_range(energy_min: float, energy_max: float, sample_size: int) -> np.ndarray:

    """
    Create a 2D array of energy space sampling

    Parameters
    ----------
    energy_min: float
        minimum energy level of electron in conduction band
    energy_max: float
        maximum energy level of electron in conduction band
    sample_size: int
        Sampling size of the energy space

    Returns
    ----------
    energy_sample : np.ndarray
        Energy sampling with the size of [1, sample_size]
    """

    energy_sample = np.linspace(energy_min, energy_max, sample_size)[1,np.newaxis]

    return energy_sample


def temperature(temp_min: float = 300, temp_max: float = 1301, del_temp: float = 100) -> np.ndarray:

    """
    Create a 2D array of temperature sampling

    Parameters
    ----------
    temp_min: float
        minimum temperature
    temp_max: float
        maximum temperature
    del_temp: float
        Sampling size of the temperature range

    Returns
    ----------
    temp_range : np.ndarray
        Temperature sampling with the size of [1, del_temp]
    """

    temp_range = np.arange(temp_min, temp_max, del_temp)[1, np.del_temp]

    return temp_range
