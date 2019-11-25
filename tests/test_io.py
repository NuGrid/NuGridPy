"""Tests of the data handling and i/o modules."""

import numpy
import os
import pkg_resources
import random
import tempfile
import unittest

import h5py

from nugridpy import data, io
from .fixtures import random_string


DATA_PATH = pkg_resources.resource_filename('tests', 'data/read-write')
HDF5_FILE = pkg_resources.resource_filename('tests', 'data/read-write/M2.00Z0.020.0000001.surf.h5')

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

        # Create a test filetype
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


class TestHdf5(unittest.TestCase):
    """Test hdf5 i/o capabilities."""

    def setUp(self):
        """Test setup."""
        # Assign the class variables for testing (only 1 for now but more later)
        self.nugrid_classvar = io.Hdf5File.NUGRID

        # Set up a readfile object and a bad filename for error testing
        self.test_read = io.Hdf5File.readfile(HDF5_FILE, self.nugrid_classvar)
        self.bad_filename = random_string()

    def test_classvars(self):
        """Test the hdf5 class variables."""
        # Test the classvars
        self.assertIsInstance(self.nugrid_classvar, dict)
        self.assertEqual(len(self.nugrid_classvar), 3)

    def test_readfile(self):
        """Test the hdf5 readfile method."""
        # Test the readfile object
        self.assertEqual(len(self.test_read), 4)
        self.assertIsInstance(self.test_read[0], h5py.AttributeManager)
        self.assertFalse(any(type(x) != dict for x in self.test_read[1:]))

        # Test the returned dict objects
        for key in self.test_read[1].keys():
            self.assertIsInstance(self.test_read[1][key], h5py.Dataset)

        for key in self.test_read[2].keys():
            self.assertIsInstance(self.test_read[2][key], h5py.AttributeManager)

        for key in self.test_read[3].keys():
            self.assertIsInstance(self.test_read[3][key], numpy.ndarray)

        # Test the right errors are being raised
        with self.assertRaises(TypeError):
            _ = io.Hdf5File.readfile(self.bad_filename)

        with self.assertRaises(OSError):
            _ = io.Hdf5File.readfile(self.bad_filename, io.Hdf5File.NUGRID)
