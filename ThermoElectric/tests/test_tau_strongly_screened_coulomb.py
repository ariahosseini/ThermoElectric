"""
Unit and regression test for the tau_strongly_screened_coulomb method.
"""

from ThermoElectric import tau_strongly_screened_coulomb
import numpy as np
from pytest import approx


def test_tau_strongly_screened_coulomb():

    DoS = np.array([[11e28]])
    dielectric = 11.7
    imp = np.array([[1e23]])
    screen_len = np.array([[1e-9]])

    expected_tau = 1.2e-12
    calculated_tau = tau_strongly_screened_coulomb(DoS=DoS, screen_len=screen_len, n_imp=imp, dielectric=dielectric)

    assert approx(expected_tau, abs=1e-13) == calculated_tau
