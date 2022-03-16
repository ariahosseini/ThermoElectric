"""
Unit and regression test for the tau_p method.
"""

from ThermoElectric import tau_p
import numpy as np
from pytest import approx


def test_tau_p():

    energy = np.array([[0.1]])
    alpha_term = np.array([[0.5]])
    tmpr = np.array([[600]])
    DoS = np.array([[1e27]])
    bulk_module = 98  # Bulk module (GPA)
    rho = 2329  # Mass density (Kg/m3)
    vel_sound = np.sqrt(bulk_module/rho)
    D_v = 2.94
    D_a = 9.5

    expected_tau = 3.0e-14
    calculated_tau = tau_p(energy=energy, alpha_term=alpha_term, D_v=D_v, D_a=D_a,
                           temp=tmpr, vel_sound=vel_sound, DoS=DoS, rho=rho)['nonparabolic_ph_lifetime']

    assert approx(expected_tau, abs=1e-15) == calculated_tau
