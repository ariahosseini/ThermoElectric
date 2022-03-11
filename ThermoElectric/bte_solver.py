import numpy as np


def electrical_properties(E, DoS, vg, Ef, dfdE, Temp, tau):

    """
    This function returns a list of thermoelectric properties
    Good references are "Near-equilibrium Transport: Fundamentals And Applications" by  Changwook Jeong and Mark S. Lundstrom and ...
    'Nanoscale Energy Transport and Conversion: A Parallel Treatment of Electrons, Molecules, Phonons, and Photons" by Gang Chen.

    :arg
            E                                       : Function object, Energy range
            DoS                                     : Function object, Electron density of state
            vg                                      : Function object, Electron group velocity
            Ef                                      : Function object, Fermi level
            dfdE                                    : Function object, Fermi window
            Temp                                    : Function object, temperature range
            tau                                     : Function object, electron total lifetime

    :returns
            coefficients                            : A list of 1 by 7, The elements are NumPy arrays of the electrical conductivity, ....
                                                      Seebecl, power factor, electron thermal conductivity, first momentum of current,
                                                      second moment of current, and the Lorenz number.
    """

    # This function returns a list of thermoelectric properties
    # See the manual for the detail of calculations

    X = DoS * vg ** 2 * dfdE  # Chi
    Y = (E - np.transpose(Ef)) * X  # Gamma
    Z = (E - np.transpose(Ef)) * Y  # Zeta

    Sigma = -1 * np.trapz(X * tau, E, axis=1) / 3 * thermoelectricProperties.e2C  # Electrical conductivity
    S = -1 * np.trapz(Y * tau, E, axis=1) / np.trapz(X * tau, E, axis=1) / Temp  # Thermopower
    PF = Sigma * S ** 2  # Power factor
    ke = -1 * (np.trapz(Z * tau, E, axis=1) - np.trapz(Y * tau, E, axis=1) ** 2 /
               np.trapz(X * tau, E,
                        axis=1)) / Temp / 3 * thermoelectricProperties.e2C  # Electron thermal conductivity

    delta_0 = np.trapz(X * tau * E, E, axis=1)
    delta_1 = np.trapz(X * tau * E, E, axis=1) / np.trapz(X * tau, E, axis=1)  # First moment of current
    delta_2 = np.trapz(X * tau * E ** 2, E, axis=1) / np.trapz(X * tau, E, axis=1)  # Second moment of current

    Lorenz = (delta_2 - delta_1 ** 2) / Temp / Temp  # Lorenz number

    coefficients = [Sigma, S[0], PF[0], ke[0], delta_1, delta_2, Lorenz[0]]

    return coefficients  # The list is 7 by numEnergySampling


def filteringEffect(self, U, E, DoS, vg, Ef, dfdE, Temp, tau_b):
    """
    This function returns list of electrical conductivity and Seebecl for the ideal filtering
    where all the electrons up to a cutoff energy level of U are completely filtered
    """

    tauUo = np.ones(len(E[0]))
    _Conductivity = [np.empty([1, len(tau_b)])]
    _Seebeck = [np.empty([1, len(tau_b)])]
    for i in np.arange(len(U)):
        tau_idl = copy.copy(tauUo)
        tau_idl[E[0] < U[i]] = 0
        tau = self.matthiessen(E, tau_idl, tau_b)
        coefficients = self.electricalProperties(E=E, DoS=DoS,
                                                 vg=vg, Ef=Ef, dfdE=dfdE, Temp=Temp, tau=tau)
        Sigma = np.expand_dims(coefficients[0], axis=0)  # Electrical conductivity
        S = np.expand_dims(coefficients[1], axis=0)  # Thermopower

        _Conductivity = np.append(_Conductivity, [Sigma], axis=0)
        _Seebeck = np.append(_Seebeck, [S], axis=0)
        del tau_idl

    Conductivity = np.delete(_Conductivity, 0, axis=0)
    Seebeck = np.delete(_Seebeck, 0, axis=0)

    return [Conductivity, Seebeck]  # The list is 2 by numEnergySampling


def phenomenological(self, U, tauo, E, DoS, vg, Ef, dfdE, Temp, tau_b):
    """
    This function returns list of electrical conductivity and Seebecl for the phenomenological model
    where a frequency independent lifetime of tauo is imposed to all the electrons up to a cutoff energy level of U
    See manual for the detail.
    """

    tauU = np.ones(len(E[0]))
    _Conductivity = [np.empty([1, 1])]
    _Seebeck = [np.empty([1, 1])]
    for _j in np.arange(len(tauo)):
        for _i in np.arange(len(U)):
            tau_ph = copy.copy(tauU)
            tau_ph[E[0] < U[_i]] = tauo[_j]
            tau = self.matthiessen(E, tau_ph, tau_b)
            coefficients = self.electricalProperties(E=E, DoS=DoS,
                                                     vg=vg, Ef=Ef, dfdE=dfdE, Temp=Temp, tau=tau)
            Sigma = np.expand_dims(coefficients[0], axis=0)
            S = np.expand_dims(coefficients[1], axis=0)
            _Conductivity = np.append(_Conductivity, [Sigma], axis=0)
            _Seebeck = np.append(_Seebeck, [S], axis=0)
            del tau_ph

    __Conductivity = np.delete(_Conductivity, 0, axis=0)
    __Seebeck = np.delete(_Seebeck, 0, axis=0)
    Conductivity = np.reshape(__Conductivity, (len(tauo), len(U)))  # Electrical conductivity
    Seebeck = np.reshape(__Seebeck, (len(tauo), len(U)))  # Thermopower

    return [Conductivity, Seebeck]  # The list is 2 by numEnergySampling