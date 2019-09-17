"""Tests of the data handling and i/o modules."""

import numpy
import os
import random
import tempfile
import unittest

from nugridpy import data, io
from .fixtures import random_string

class TestAscii(unittest.TestCase):
    """Test ascii i/o capabilities."""

    def setUp(self):
        """Test setup."""
        # Set the test data path
        self.data_path = 'tests/data/read-write/'

        # Create a random string to use for testing bad filenames
        self.bad_file = random_string()

        # Name the FILE_TYPES
        self.ftypes = io.TextFile.FILE_TYPES

        # Set up test data dimensions
        self.row_width = random.randint(2, 10)
        self.col_height = random.randint(2, 99)

        self.TEST_DATA = {
            'name': 'TEST_DATA',
            'parser': None,
            'skip_header': 0,
            'names': None
            }

    def test_TextFile(self):
        """Test the io.TextFile class."""
        for ftype in self.ftypes:
            # Ensure every filetype dict is fully populated (4 keys)
            self.assertEqual(len(ftype), 4)

            # Ensure every filetype has a parser
            if ftype['parser'] is not None:
                self.assertTrue(callable(ftype['parser']))

    def test_readfile(self):
        """Test the io.TextFile.readfile class method."""
        # Read a working test file with io.TextFile.readfile
        self.test_read = io.TextFile.readfile(os.path.join(self.data_path, 'history.data'),
            io.TextFile.MESA_DATA)

        # Ensure proper OSError when filename can't be found
        for ftype in self.ftypes:
            with self.assertRaisesRegex(OSError, '{} not found.'.format(self.bad_file)):
                _ = io.TextFile.readfile(self.bad_file, ftype)

        # Build a tempfile with random_string() data
        with tempfile.NamedTemporaryFile(mode='r+') as testf:
            # Populate the tempfile
            for row in range(self.col_height):
                for col in range(self.row_width):
                    testf.write(random_string())
                    testf.write(' ')
                testf.write('\n')

            # Read the tempfile with readfile
            testf.seek(0)
            temp_headers, temp_data = io.TextFile.readfile(testf.name, self.TEST_DATA)

            # Test the read
            self.assertIsInstance(temp_headers, numpy.ndarray)
            self.assertIsInstance(temp_data, numpy.ndarray)

            # Ensure the dimensions are correct (minus 1 for column names)
            self.assertEqual(temp_headers.size, 0) # No headers in TEST_DATA
            self.assertEqual(temp_data.size, self.col_height - 1) # Data column minus name
            self.assertEqual(len(temp_data.dtype), self.row_width) # Column names

        # Populate a tempfile with column inconsistencies (break write-in)
        with tempfile.NamedTemporaryFile(mode='r+') as test_badf:
            # Populate the tempfile with random data until row_width-1
            for row in range(self.col_height):
                for col in range(self.row_width):
                    # Break before writing final data column
                    if row > 0 and col == self.row_width-1:
                        break
                    # Otherwise write random string
                    test_badf.write(random_string())
                    test_badf.write(' ')
                test_badf.write('\n')

            # Ensure inconsistencies between header and data columns raise error
            with self.assertRaisesRegex(ValueError, 'Some errors were detected !'):
                test_badf.seek(0)
                _ = io.TextFile.readfile(test_badf.name, self.TEST_DATA)
