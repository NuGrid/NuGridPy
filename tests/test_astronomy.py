"""Test the astronomy module functionality."""

import unittest
import random
import numpy as np
from scipy import integrate

from nugridpy import constants, astronomy
from .fixtures import random_string


class TestFunctions(unittest.TestCase):
    """Test some of the functions in the astronomy module."""

    def setUp(self):
        """Declare variables for self object."""
        self.m1 = random.uniform(0, constants.imf_m1)
        self.m2 = random.uniform(constants.imf_m1, constants.imf_m2)
        self.m3 = random.uniform(constants.imf_m2, 1.e30)
        self.m_arr = np.linspace(self.m1, self.m2, 10)
        self.imf_arr = np.asarray([astronomy.imf(x) for x in self.m_arr])

    def test_imf(self):
        """Test that imf and int_imf_dm are working properly."""
        astro_imf = astronomy.int_imf_dm(
            self.m1, self.m2, self.m_arr, self.imf_arr, bywhat='bynumber')
        test_imf = integrate.trapz(self.imf_arr, self.m_arr)
        self.assertAlmostEqual(astro_imf, test_imf, places=3)

    def test_args(self):
        """Test that int_imf_dm raises appropriate errors."""
        with self.assertRaises(ValueError):
            _ = astronomy.int_imf_dm(
                self.m1, self.m2, self.m_arr, self.imf_arr, bywhat=None)

        with self.assertRaises(ValueError):
            _ = astronomy.int_imf_dm(
                self.m1, self.m2, self.m_arr, self.imf_arr, integral=None)

        with self.assertRaises(TypeError):
            _ = astronomy.int_imf_dm(
                self.m1, self.m2, [self.m1, self.m2], self.imf_arr)

        with self.assertRaises(TypeError):
            _ = astronomy.int_imf_dm(
                self.m1, self.m2, self.m_arr, [self.imf_arr[0], self.imf_arr[-1]])


class TestDecorator(unittest.TestCase):
    """Test the functionality of the attach_constants/read-only property decorator."""

    def setUp(self):
        """Declare variables for self object and define test functions."""
        self.arguments = [random.random()
                          for _ in range(random.randint(1, 10))]
        self.kwarguments = {
            random_string(): random.random()
            for _ in range(random.randint(1, 10))}
        self.test_constant = constants.Constant(1., 'Test constant', 'None')

        # decorated test function for taking/returning args/kwargs
        @astronomy.attach_constants(self.test_constant)
        def arg_test(*args, **kwargs):
            return args, kwargs

        self.arg_test = arg_test
        self.result = arg_test(*self.arguments, **self.kwarguments)

    def test_args_kwargs(self):
        """Make sure decorator and function args/kwargs are passed around properly."""

        self.assertEqual(tuple(self.arguments), self.result[0])
        self.assertEqual(self.kwarguments, self.result[1])

    def test_read_only(self):
        """Ensure that correct constants are attached to functions as read-only property."""

        self.assertEqual(self.arg_test.constants, (self.test_constant,))

        with self.assertRaisesRegex(AttributeError, 'can\'t set attribute'):
            self.arg_test.constants = self.test_constant
