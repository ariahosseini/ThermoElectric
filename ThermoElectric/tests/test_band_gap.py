"""
Unit and regression test for the band_gap method.
"""

from ThermoElectric import band_gap
import numpy as np


def test_band_gap():

    tmpr = np.array([300, 400, 500])  # Temperature
    Eg_o = 2  # Band gap at zero kelvin
    Ao_term = 1e-3  # First experimentally-fit value
    Bo_term = 400  # Second experimentally-fit value

    expected_gap = np.array([1.8200, 1.7333, 1.6429])
    calculated_gap = band_gap(Eg_o=Eg_o, Ao=Ao_term, Bo=Bo_term, temp=tmpr)

    assert (expected_gap == calculated_gap).all
