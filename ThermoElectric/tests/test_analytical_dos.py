"""
Unit and regression test for the analytical_gap method.
"""

import ThermoElectric as TE
import numpy as np


def test_analytical_gap():

    """
    Test analytical_gap method to make sure it approximates the electron density of state for parabolic
    and non-parabolic bands correctly.
    """

    energy = np.array([[0.1, 0.2]])
    e_eff_mass = 1.08 * 9.109e-31
    nonparabolic_term = np.array([[0.5, 0.2]])  # Temperature

    expected_dos = np.array([[1.2396, 1.3200], [1.1907, 1.2683]])*1.0e+28
    calculated_dos = TE.analytical_dos(range_energy=energy, electron_eff_mass=e_eff_mass,
                                       nonparabolic_term=nonparabolic_term)['DoS_nonparabolic']

    assert (expected_dos == calculated_dos).all
