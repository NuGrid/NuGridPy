"""Test the astronomy module functionality."""

import unittest
import random
import numpy as np
from nugridpy import constants
from nugridpy import astronomy
from scipy import integrate
from unittest.mock import patch

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
            _ = astronomy.int_imf_dm(self.m1, self.m2, self.m_arr, self.imf_arr, integral=None)

        with self.assertRaises(TypeError):
            _ = astronomy.int_imf_dm(self.m1, self.m2, [self.m1, self.m2], self.imf_arr)

        with self.assertRaises(TypeError):
            _ = astronomy.int_imf_dm(self.m1, self.m2, self.m_arr, [self.imf_arr[0], self.imf_arr[-1]])


class TestDecorator(unittest.TestCase):
    """Test the functionality of the attach_constants/read-only property decorator."""

    def setUp(self):
        self.parameter = random.random()

    @patch('nugridpy.constants.Constant')
    def test_args_kwargs(self, mocked_constant):
        """Make sure decorator and function args/kwargs are passed around properly."""

        @astronomy.attach_constants(mocked_constant)
        def arg_test(parameter, kwarg=str(self.parameter)):
            return parameter, kwarg

        self.assertEqual(self.parameter, arg_test(self.parameter)[0])
        self.assertEqual(str(self.parameter), arg_test(self.parameter)[1])
        self.assertEqual(arg_test.constants, (mocked_constant,))

    def test_read_only(self):
        """Ensure that constants are a read-only property attached to astronomy functions."""
        @astronomy.attach_constants()
        def f_test():
            pass

        with self.assertRaisesRegex(AttributeError, 'can\'t set attribute'):
            f_test.constants = self.parameter
