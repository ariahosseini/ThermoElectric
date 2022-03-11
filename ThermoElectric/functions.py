"""Provide the primary functions."""

from .util import *
import numpy as np


def band_gap(Eg_o: float, Ao: float, Bo: float, temp: np.ndarray = None):

    """
    This method uses Eg(T)=Eg(T=0)-Ao*T**2/(T+Bo) to approximate the temperature dependency of the dielectrics band gap.
    A good reference is "Properties of Advanced Semiconductor Materials" by Michael E. Levinshtein.

    Parameters
    ----------
    Eg_o: float
        The band gap at zero Kelvin
    Ao: float
        Experimentally fitted parameter
    Bo: float
        Experimentally fitted parameter
    temp: np.ndarray
        Temperature range

    Returns
    -------
    Eg: np.ndarray
        Temperature-dependent band gap
    """

    if temp is None:
        '''Use default temperature range (300K up to 1300K with step of 100K) '''
        T = temperature()
    else:
        T = temp

    Eg = Eg_o - Ao * np.divide(T ** 2, T + Bo)

    return Eg


def analyticalDoS(range_energy, electron_eff_mass, nonparabolic_term):

    """
       This function approximate the electron density of state for parabolic and non-parabolic bands
       in case DFT calculation is not available

    Parameters
    ----------
    range_energy: np.ndarray
        Electron energy range
    electron_eff_mass: float
        Electron effective mass
    nonparabolic_term: np.ndarray
        Non-parabolic term (shows the mixture of S and P orbitals)

    Returns
    -------
    DoS: np.ndarray
        First row is phonon density of states for non-parabolic band and
        the second row is the phonon density of states for non-parabolic.
        The array size is [2, numEnergySampling].
    """

    hBar = 6.582119e-16  # Reduced Planck constant in eV.s
    e2C = 1.6021765e-19  # e to Coulomb unit change

    DoS_nonparabolic = 1 / np.pi ** 2 * np.sqrt(2 * range_energy * (1 + range_energy * nonparabolic_term.T)) \
                       * np.sqrt(electron_eff_mass / hBar ** 2) ** 3 \
                       * (1 + (2 * range_energy * nonparabolic_term.T)) / e2C ** (3.0 / 2)

    DoS_parabolic = np.sqrt(range_energy) / np.pi ** 2 * np.sqrt(2) / hBar ** 3 \
                    * electron_eff_mass ** (3.0 / 2) / e2C ** (3.0 / 2)

    DoS = np.array(DoS_nonparabolic, DoS_parabolic)

    return DoS


if __name__ == "__main__":

    print('ThermoElectric')

