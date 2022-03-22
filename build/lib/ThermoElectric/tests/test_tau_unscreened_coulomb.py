"""
Unit and regression test for the tau_screened_coulomb method.
"""

from ThermoElectric import tau_unscreened_coulomb
import numpy as np
from pytest import approx


def test_tau_unscreened_coulomb():

    energy = np.array([[0.1]])
    e_eff_mass = np.array([[0.23 * 9.109e-31]])
    dielectric = 11.7
    imp = np.array([[1e21]])

    expected_tau = 1.2e-08
    calculated_tau = tau_unscreened_coulomb(energy=energy, mass_c=e_eff_mass, n_imp=imp, dielectric=dielectric)

    assert approx(expected_tau, abs=1e-7) == calculated_tau
