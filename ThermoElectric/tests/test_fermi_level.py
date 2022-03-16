"""
Unit and regression test for the fermi_level method.
"""

from ThermoElectric import fermi_level
import numpy as np
from pytest import approx


def test_fermi_level():

    """
    Test fermi_level method.
    """
    carrier = np.array([[1.e23]])
    energy = np.array([[0.2]])
    density = np.array([[0.2*1e28]])
    tmpr = np.array([[300]])
    expected_fermi_level = -0.14370037
    calculated_fermi_level = fermi_level(carrier=carrier, energy=energy, density=density,
                                            Nc=None, Ao=5e21, temp=tmpr)[1][0]

    assert approx(calculated_fermi_level, abs=1e-5) == expected_fermi_level
