import numpy as np
from copy import copy
from .functions import *


def electrical_properties(E: np.ndarray, DoS: np.ndarray, vg: np.ndarray, Ef: np.ndarray,
                          dfdE: np.ndarray, temp: np.ndarray, tau: np.ndarray) -> dict:

    """
    This function returns electronic properties.
    Good references on this topic are:
        — "Near-equilibrium Transport: Fundamentals And Applications" by Changwook Jeong and Mark S. Lundstrom
        — "Nanoscale Energy Transport and Conversion:
           A Parallel Treatment of Electrons, Molecules, Phonons, and Photons" by Gang Chen.

    Parameters
    ----------
    E: np.ndarray
        Energy range
    DoS: np.ndarray
        Electron density of state
    vg: np.ndarray
        Group velocity
    Ef: np.ndarray
        Fermi level
    dfdE: np.ndarray
        Fermi window
    temp: np.ndarray
        Temperature range
    tau:  np.ndarray
        Lifetime

    Returns
    -------
    coefficients: dict
        Linear BTE prediction of electrical properties
    """

    e2C = 1.6021765e-19  # e to Coulomb unit change

    X = DoS * vg ** 2 * dfdE  # Chi
    Y = (E - Ef.T) * X  # Gamma
    Z = (E - Ef.T) * Y  # Zeta

    Sigma = -1 * np.trapz(X * tau, E, axis=1) / 3 * e2C  # Electrical conductivity

    S = -1 * np.trapz(Y * tau, E, axis=1) / np.trapz(X * tau, E, axis=1) / temp  # Thermopower

    PF = Sigma * S ** 2  # Power factor

    ke = -1 * (np.trapz(Z * tau, E, axis=1) - np.trapz(Y * tau, E, axis=1) ** 2 /
               np.trapz(X * tau, E, axis=1)) / temp / 3 * e2C  # Electron thermal conductivity

    delta = np.trapz(X * tau * E, E, axis=1) / np.trapz(X * tau, E, axis=1)  # First moment of current

    delta_ = np.trapz(X * tau * E ** 2, E, axis=1) / np.trapz(X * tau, E, axis=1)  # Second moment of current

    Lorenz = (delta_ - delta ** 2) / temp**2  # Lorenz number

    coefficients = {'Electrical_conductivity': Sigma, 'Seebeck_coefficient': S[0], 'Power_factor': PF[0],
                    'Thermal_conductivity': ke[0], 'Lorenz_number': Lorenz[0]}

    return coefficients


def filtering_effect(Uo, E, DoS, vg, Ef, dfdE, temp, tau_bulk):

    """
    This function returns electric properties for the ideal filtering —  where all the electrons up to a cutoff energy
    level of Uo are completely hindered.

    Parameters
    ----------
    Uo: np.ndarray
        Barrier height
    E: np.ndarray
        Energy range
    DoS: np.ndarray
        Electron density of state
    vg: np.ndarray
        Group velocity
    Ef: np.ndarray
        Fermi level
    dfdE: np.ndarray
        Fermi window
    temp: np.ndarray
        Temperature range
    tau_bulk:  np.ndarray
        Lifetime in bulk material

    Returns
    -------
    output: dict
        Linear BTE prediction of electrical properties
    """

    tau_Uo = np.ones(len(E[0]))
    _Conductivity = [np.empty([1, len(tau_bulk)])]
    _Seebeck = [np.empty([1, len(tau_bulk)])]

    for i in np.arange(len(Uo)):

        tau_idl = copy(tau_Uo)
        tau_idl[E[0] < Uo[i]] = 0
        tau = matthiessen(E, tau_idl, tau_bulk)

        coefficients = electrical_properties(E=E, DoS=DoS, vg=vg, Ef=Ef, dfdE=dfdE, temp=temp, tau=tau)
        Sigma = coefficients['Electrical_conductivity']  # Electrical conductivity
        S = coefficients['Seebeck_coefficient']  # Thermopower

        _Conductivity = np.append(_Conductivity, [Sigma], axis=0)
        _Seebeck = np.append(_Seebeck, [S], axis=0)

        del tau_idl

    Conductivity = np.delete(_Conductivity, 0, axis=0)
    Seebeck = np.delete(_Seebeck, 0, axis=0)

    output = {'Electrical_conductivity': Conductivity, 'Seebeck_coefficient': Seebeck}

    return output


def phenomenological(Uo, tau_o, E, DoS, vg, Ef, dfdE, temp, tau_bulk):

    """
    This function returns electric properties for the phenomenological filtering —  where a frequency independent
    lifetime of tau_o is imposed to all the electrons up to a cutoff energy level of Uo

    Parameters
    ----------
    Uo: np.ndarray
        Barrier height
    tau_o: np.ndarray
        Phenomenological lifetime
    E: np.ndarray
        Energy range
    DoS: np.ndarray
        Electron density of state
    vg: np.ndarray
        Group velocity
    Ef: np.ndarray
        Fermi level
    dfdE: np.ndarray
        Fermi window
    temp: np.ndarray
        Temperature range
    tau_bulk:  np.ndarray
        Lifetime in bulk material

    Returns
    -------
    output: dict
        Linear BTE prediction of electrical properties
    """

    tau_u = np.ones(len(E[0]))
    _Conductivity = [np.empty([1, 1])]
    _Seebeck = [np.empty([1, 1])]
    for _j in np.arange(len(tau_o)):
        for _i in np.arange(len(Uo)):
            tau_ph = copy(tau_u)
            tau_ph[E[0] < Uo[_i]] = tau_o[_j]
            tau = matthiessen(E, tau_ph, tau_bulk)
            coefficients = electrical_properties(E=E, DoS=DoS, vg=vg, Ef=Ef, dfdE=dfdE, temp=temp, tau=tau)
            Sigma = coefficients['Electrical_conductivity']  # Electrical conductivity
            S = coefficients['Seebeck_coefficient']  # Thermopower
            _Conductivity = np.append(_Conductivity, [Sigma], axis=0)
            _Seebeck = np.append(_Seebeck, [S], axis=0)
            del tau_ph

    __Conductivity = np.delete(_Conductivity, 0, axis=0)
    __Seebeck = np.delete(_Seebeck, 0, axis=0)
    Conductivity = np.reshape(__Conductivity, (len(tau_o), len(Uo)))  # Electrical conductivity
    Seebeck = np.reshape(__Seebeck, (len(tau_o), len(Uo)))  # Thermopower

    output = {'Electrical_conductivity': Conductivity, 'Seebeck_coefficient': Seebeck}

    return output
