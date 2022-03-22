"""
Unit and regression test for the fermi_self_consistent method.
"""

from ThermoElectric import matthiessen
import numpy as np
from pytest import approx


def test_fermi_self_consistent():

    energy = np.array([[0.2]])
    tau_one = np.array([[0.5]])
    tau_two = np.array([[1.0]])
    expected_tau= 0.3333333333
    calculated_tau = matthiessen(energy,tau_one, tau_two)

    assert approx(calculated_tau, abs=1e-5) == expected_tau
