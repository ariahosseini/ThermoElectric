import numpy as np
from copy import deepcopy
from .functions import matthiessen

def electrical_properties(E: np.ndarray, DoS: np.ndarray, vg: np.ndarray, Ef: np.ndarray,
                          dfdE: np.ndarray, temp: np.ndarray, tau: np.ndarray) -> dict:
    """
    Computes the electronic properties using the Boltzmann Transport Equation (BTE).
    
    References:
    -----------
    - Changwook Jeong and Mark S. Lundstrom, "Near-equilibrium Transport: Fundamentals and Applications"
    - Gang Chen, "Nanoscale Energy Transport and Conversion: A Parallel Treatment of Electrons, Molecules, Phonons, and Photons"
    
    Parameters:
    -----------
    E : np.ndarray
        Energy range (eV)
    DoS : np.ndarray
        Density of states (states/eV)
    vg : np.ndarray
        Group velocity (m/s)
    Ef : np.ndarray
        Fermi level (eV)
    dfdE : np.ndarray
        Fermi-Dirac derivative
    temp : np.ndarray
        Temperature (K)
    tau : np.ndarray
        Relaxation time (s)
    
    Returns:
    --------
    dict
        Contains electrical conductivity, Seebeck coefficient, power factor, thermal conductivity, and Lorenz number.
    """
    e2C = 1.6021765e-19  # Elementary charge in Coulombs

    # Precompute useful terms
    X = DoS * vg ** 2 * dfdE  # χ (Chi)
    Y = (E - Ef.T) * X  # Γ (Gamma)
    Z = (E - Ef.T) * Y  # ζ (Zeta)

    # Transport properties
    Sigma = -np.trapz(X * tau, E, axis=1) / 3 * e2C  # Electrical conductivity (S/m)
    S = -np.trapz(Y * tau, E, axis=1) / np.trapz(X * tau, E, axis=1) / temp  # Seebeck coefficient (V/K)
    PF = Sigma * S ** 2  # Power factor (W/K²·m)
    ke = -1 * (np.trapz(Z * tau, E, axis=1) - np.trapz(Y * tau, E, axis=1) ** 2 /
               np.trapz(X * tau, E, axis=1)) / temp / 3 * e2C  # Thermal conductivity (W/m·K)
    
    # Lorenz number
    delta = np.trapz(X * tau * E, E, axis=1) / np.trapz(X * tau, E, axis=1)  # First moment of current
    delta_ = np.trapz(X * tau * E ** 2, E, axis=1) / np.trapz(X * tau, E, axis=1)  # Second moment of current
    Lorenz = (delta_ - delta ** 2) / temp**2  # Lorenz number (WΩ/K²)

    return {
        'Electrical_conductivity': Sigma,
        'Seebeck_coefficient': S[0],
        'Power_factor': PF[0],
        'Thermal_conductivity': ke[0],
        'Lorenz_number': Lorenz[0]
    }


def filtering_effect(Uo: np.ndarray, E: np.ndarray, DoS: np.ndarray, vg: np.ndarray, Ef: np.ndarray,
                     dfdE: np.ndarray, temp: np.ndarray, tau_bulk: np.ndarray) -> dict:
    """
    Calculates the electrical properties under ideal filtering conditions, where all electrons
    up to a cutoff energy Uo are blocked.
    
    Parameters:
    -----------
    Uo : np.ndarray
        Barrier height (eV)
    E : np.ndarray
        Energy range (eV)
    DoS : np.ndarray
        Density of states (states/eV)
    vg : np.ndarray
        Group velocity (m/s)
    Ef : np.ndarray
        Fermi level (eV)
    dfdE : np.ndarray
        Fermi-Dirac derivative
    temp : np.ndarray
        Temperature (K)
    tau_bulk : np.ndarray
        Bulk relaxation time (s)
    
    Returns:
    --------
    dict
        Contains electrical conductivity and Seebeck coefficient.
    """
    tau_Uo = np.ones(E.shape[1])
    Conductivity, Seebeck = [], []

    for uo in Uo:
        tau_filtered = deepcopy(tau_Uo)
        tau_filtered[E[0] < uo] = 0  # Apply ideal filtering
        
        tau = matthiessen(E, tau_filtered, tau_bulk)  # Compute total relaxation time
        coefficients = electrical_properties(E, DoS, vg, Ef, dfdE, temp, tau)

        Conductivity.append(coefficients['Electrical_conductivity'])
        Seebeck.append(coefficients['Seebeck_coefficient'])

    return {
        'Electrical_conductivity': np.array(Conductivity),
        'Seebeck_coefficient': np.array(Seebeck)
    }


def phenomenological(Uo: np.ndarray, tau_o: np.ndarray, E: np.ndarray, DoS: np.ndarray, vg: np.ndarray,
                     Ef: np.ndarray, dfdE: np.ndarray, temp: np.ndarray, tau_bulk: np.ndarray) -> dict:
    """
    Calculates the electrical properties with phenomenological filtering, where a lifetime tau_o
    is imposed up to a cutoff energy Uo.
    
    Parameters:
    -----------
    Uo : np.ndarray
        Barrier height (eV)
    tau_o : np.ndarray
        Phenomenological lifetime (s)
    E : np.ndarray
        Energy range (eV)
    DoS : np.ndarray
        Density of states (states/eV)
    vg : np.ndarray
        Group velocity (m/s)
    Ef : np.ndarray
        Fermi level (eV)
    dfdE : np.ndarray
        Fermi-Dirac derivative
    temp : np.ndarray
        Temperature (K)
    tau_bulk : np.ndarray
        Bulk relaxation time (s)
    
    Returns:
    --------
    dict
        Contains electrical conductivity and Seebeck coefficient.
    """
    tau_Uo = np.ones(E.shape[1])
    Conductivity, Seebeck = [], []

    for tau in tau_o:
        for uo in Uo:
            tau_ph = deepcopy(tau_Uo)
            tau_ph[E[0] < uo] = tau  # Apply phenomenological filtering
            
            total_tau = matthiessen(E, tau_ph, tau_bulk)
            coefficients = electrical_properties(E, DoS, vg, Ef, dfdE, temp, total_tau)
            
            Conductivity.append(coefficients['Electrical_conductivity'])
            Seebeck.append(coefficients['Seebeck_coefficient'])

    return {
        'Electrical_conductivity': np.array(Conductivity).reshape(len(tau_o), len(Uo)),
        'Seebeck_coefficient': np.array(Seebeck).reshape(len(tau_o), len(Uo))
    }
