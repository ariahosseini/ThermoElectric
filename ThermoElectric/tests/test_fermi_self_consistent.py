"""
Unit and regression test for the fermi_self_consistent method.
"""

from ThermoElectric import fermi_self_consistent
import numpy as np
from pytest import approx


def test_fermi_self_consistent():

    """
    Test fermi_self_consistent method.
    """
    carrier = np.array([[1.e23]])
    energy = np.array([[0.2]])
    density = np.array([[0.2*1e28]])
    tmpr = np.array([[300]])
    fermi_level_jd = np.array([[1.e23], [-0.14370037]])
    expected_fermi_level = -0.5437
    calculated_fermi_level = fermi_self_consistent(carrier=carrier, temp=tmpr, energy=energy,
                                                   density=density, fermi_levels = fermi_level_jd)[1][0]

    assert approx(calculated_fermi_level, abs=1e-5) == expected_fermi_level
