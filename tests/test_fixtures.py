import unittest

from .fixtures import random_string


class TestFixtures(unittest.TestCase):
    """Class testing the implemented fixtures."""

    def test_random_string(self):
        """Test the random_string function."""

        num_strings = 1000
        string_collection = [
            random_string() for _ in range(num_strings)
        ]

        # Make a set to remove duplicates
        string_collection = set(string_collection)
        self.assertEqual(len(string_collection), num_strings)
