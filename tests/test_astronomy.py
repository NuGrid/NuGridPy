"""Test the astronomy module functionality."""

import unittest
import random
import numpy as np
from nugridpy import constants
from nugridpy import astronomy
from scipy import integrate

class TestFunctions(unittest.TestCase):
    """Test some of the functions in the astronomy module."""

    def setUp(self):
        self.m1 = random.uniform(0, constants.imf_m1)
        self.m2 = random.uniform(constants.imf_m1, constants.imf_m2)
        self.m3 = random.uniform(constants.imf_m2, 1.e30)
        self.m_arr = np.linspace(self.m1, self.m2, 10)
        self.imf_arr = np.asarray([astronomy.imf(x) for x in self.m_arr])

    def test_imf(self):
        """Test that imf and int_imf_dm are working properly."""
        astro_imf = astronomy.int_imf_dm(self.m1, self.m2, self.m_arr, self.imf_arr, bywhat='bynumber')
        test_imf = integrate.trapz(self.imf_arr, self.m_arr)
        self.assertAlmostEqual(astro_imf, test_imf, places=3)

    def test_args(self):
        """Test that int_imf_dm raises appropriate errors."""
        with self.assertRaises(ValueError):
            _ = astronomy.int_imf_dm(self.m1, self.m2, self.m_arr, self.imf_arr, bywhat=None)

        with self.assertRaises(ValueError):
            _ = astronomy.int_imf_dm(self.m1, self.m2, self.m_arr, self.imf_arr, bynumber=None)

        with self.assertRaises(TypeError):
            _ = astronomy.int_imf_dm(self.m1, self.m2, [self.m1, self.m2], self.imf_arr)

        with self.assertRaises(TypeError):
            _ = astronomy.int_imf_dm(self.m1, self.m2, self.m_arr, [self.imf_arr[0], self.imf_arr[-1]])
