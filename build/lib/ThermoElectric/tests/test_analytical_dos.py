"""
Unit and regression test for the analytical_gap method.
"""

from ThermoElectric import analytical_dos
import numpy as np
from pytest import approx


def test_analytical_gap():

    energy = np.array([[0.1]])
    e_eff_mass = 1.08 * 9.109e-31
    nonparabolic_term = np.array([[0.5]])  # Temperature

    expected_dos = 2.725*1.0e+27
    calculated_dos = analytical_dos(range_energy=energy, electron_eff_mass=e_eff_mass,
                                       nonparabolic_term=nonparabolic_term)['DoS_nonparabolic'][0]

    assert approx(expected_dos, abs=1e-1*1e27) == calculated_dos[0]
