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

    energy_sample = np.linspace(energy_min, energy_max, sample_size)[np.newaxis, :]

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

    temp_range = np.arange(temp_min, temp_max, del_temp)[np.newaxis, :]

    return temp_range


def fermi_distribution(energy: np.ndarray, fermi_level: np.ndarray, temp: np.ndarray = None) -> np.ndarray:

    """
    A function to compute the Fermi distribution and
    the Fermi window — the first energy derivative of the Fermi level

    Parameters
    ----------
    energy: np.ndarray
        Energy range in conduction band
    fermi_level: np.ndarray
        Fermi level
    temp: np.ndarray
        Temperature range

    Returns
    ----------
    fermi : np.ndarray
        The first row is the Fermi distribution and the second row is the derivative (Fermi window)
    """

    k_bolt = 8.617330350e-5  # Boltzmann constant in eV/K

    if temp is None:
        T = temperature()
    else:
        T = temp

    xi = np.exp((energy - fermi_level.T) / T.T / k_bolt)

    fermiDirac = 1 / (xi + 1)  # Fermi distribution
    dfdE = -1 * xi / (1 + xi) ** 2 / T.T / k_bolt  # Fermi window
    fermi = np.array([fermiDirac, dfdE])

    return fermi


def matthiessen(*args) -> np.ndarray:

    """
    A function to compute the total lifetime using Matthiessen's rule

    Parameters
    ----------
    *args: np.ndarray
        Electron lifetime from different scattering sources

    Returns
    ----------
    tau : np.ndarray
        The total electron lifetime
    """

    tau = 1. / sum([1. / arg for arg in args])
    tau[np.isinf(tau)] = 0

    return tau
