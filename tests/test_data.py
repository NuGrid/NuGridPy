import os
import tempfile
import unittest

from nugridpy import data, io

class TestAsciiData(unittest.TestCase):
    """Testing data.py ASCII data-handling."""

    def setUp(self):
        """Setting up variables and other test structures."""
        # Set the test data path
        self.data_path = 'tests/data/read-write/'

        # Get the FileIn filetypes
        self.filetypes = data.FileIn.filetypes

        # Make a user-facing data object with data.FileIn
        self.test_data = data.FileIn(os.path.join(self.data_path,
            'history.data'), io.TextFile.MESA_DATA)

    def test_FileIn(self):
        """Test the data.FileIn class."""
        # Make sure the filetypes class variable is assigned/non-empty
        self.assertIsInstance(self.filetypes, tuple)
        self.assertGreater(len(self.filetypes), 0)

        # Check the datacol/header names attributes are populated
        self.assertIsInstance(self.test_data.header_names, tuple)
        self.assertGreaterEqual(len(self.test_data.header_names), 0)
        self.assertIsInstance(self.test_data.data_names, tuple)
        self.assertGreater(len(self.test_data.data_names), 0)

    def test_MesaData(self):
        """Test the data.MesaData class."""
        # Check that MesaData works and reads in a custom filename
        with tempfile.NamedTemporaryFile(prefix='star.log') as tmpf:
            # Write history.data to tempfile
            hd = open(os.path.join(self.data_path, 'history.data'))
            tmpf.write(hd.read().encode())
            hd.close()

            # Read the tempfile as a MesaData instance
            test_alias = data.MesaData(tempfile.tempdir, history_name=tmpf.name, 
                history_only=True)

            # Ensure the instance was initialized properly
            self.assertIsInstance(test_alias, data.MesaData)
            self.assertIsInstance(test_alias.history_data, data.MesaFile)

            # Ensure history_only worked properly, no profile data
            with self.assertRaises(AttributeError):
                _ = test_alias.profiles_index
            with self.assertRaises(AttributeError):
                _ = test_alias.profiles

            # Assign data to variables
            test_headers = test_alias.history_data.header_names
            test_data = test_alias.history_data.data_names

            # Ensure the history.data read worked
            self.assertIsInstance(test_headers, tuple)
            self.assertGreater(len(test_headers), 0)
            self.assertIsInstance(test_data, tuple)
            self.assertGreater(len(test_data), 0)
