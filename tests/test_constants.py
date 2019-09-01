"""Tests for constants.py module."""

import unittest
import random
from nugridpy.constants import Constant


class TestConstant(unittest.TestCase):
    """Test the functionality of the Constant class."""

    def setUp(self):
        self.value = random.uniform(-1e30, 1e30)
        self.integer = random.randint(-1e10, 1e10)
        self.test_constant = Constant(self.value, 'Test', 'None')

        # Check that Constant was initialized properly as float type
        self.assertTrue(isinstance(self.test_constant, float))

    def test_constant_ops(self):
        """Constant instances should operate as their float values."""
        self.assertAlmostEqual(self.test_constant**2, self.value**2, places=2)
        self.assertAlmostEqual(self.test_constant + self.test_constant,
                               self.value + self.value, places=2)
        self.assertAlmostEqual(self.test_constant / self.integer,
                               self.value / self.integer, places=2)
        self.assertAlmostEqual(self.test_constant - self.integer,
                               self.value - self.integer, places=2)

    def test_args(self):
        """Ensure correct errors are raised if arguments are incorrectly passed."""
        with self.assertRaises(TypeError):
            _ = Constant(self.value, self.value, self.value)

        with self.assertRaises(TypeError):
            _ = Constant(self.value)
