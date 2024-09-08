import numpy as np


def energy_range(energy_min: float, energy_max: float, sample_size: int) -> np.ndarray:
    """
    Generate a 2D array of energy values sampled between the minimum and maximum energy levels.

    Parameters
    ----------
    energy_min : float
        Minimum energy level of the electron in the conduction band.
    energy_max : float
        Maximum energy level of the electron in the conduction band.
    sample_size : int
        Number of samples to generate between energy_min and energy_max.

    Returns
    -------
    np.ndarray
        A 1xN array of energy values sampled between the specified range.
    """
    return np.linspace(energy_min, energy_max, sample_size)[np.newaxis, :]


def temperature(temp_min: float = 300, temp_max: float = 1301, del_temp: float = 100) -> np.ndarray:
    """
    Generate a 2D array of temperature values sampled between the minimum and maximum temperature.

    Parameters
    ----------
    temp_min : float, optional
        Minimum temperature (default is 300K).
    temp_max : float, optional
        Maximum temperature (default is 1301K).
    del_temp : float, optional
        Step size for sampling the temperature range (default is 100K).

    Returns
    -------
    np.ndarray
        A 1xN array of temperature values sampled between the specified range.
    """
    return np.arange(temp_min, temp_max, del_temp)[np.newaxis, :]


def fermi_distribution(energy: np.ndarray, fermi_level: np.ndarray, temp: np.ndarray = None) -> np.ndarray:
    """
    Compute the Fermi distribution and its first energy derivative (Fermi window) for given energy levels.

    Parameters
    ----------
    energy : np.ndarray
        Array of energy levels in the conduction band.
    fermi_level : np.ndarray
        Fermi level values.
    temp : np.ndarray, optional
        Temperature array (default is None, uses the default temperature function if not provided).

    Returns
    -------
    np.ndarray
        A 2xN array where the first row is the Fermi distribution and the second row is the Fermi window (derivative).
    """
    k_bolt = 8.617330350e-5  # Boltzmann constant in eV/K

    if temp is None:
        T = temperature()
    else:
        T = temp

    xi = np.exp((energy - fermi_level.T) / (T.T * k_bolt))

    fermi_dirac = 1 / (xi + 1)  # Fermi distribution
    dfdE = -xi / ((1 + xi) ** 2 * T.T * k_bolt)  # Fermi window (derivative of Fermi distribution)
    
    return np.stack((fermi_dirac, dfdE))


def matthiessen(energy: np.ndarray, *args: np.ndarray) -> np.ndarray:
    """
    Calculate the total electron lifetime using Matthiessen's rule.

    Parameters
    ----------
    energy : np.ndarray
        Energy levels of electrons.
    *args : np.ndarray
        Individual scattering lifetimes from different sources.

    Returns
    -------
    np.ndarray
        Total electron lifetime computed using Matthiessen's rule.
    """
    tau = 1.0 / np.sum([1.0 / arg for arg in args], axis=0)
    tau[np.isinf(tau)] = 0  # Set infinite lifetime (no scattering) to zero for numerical stability.
    
    return tau
