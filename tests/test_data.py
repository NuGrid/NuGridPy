
import os
import random
import tempfile
import unittest
from unittest import TestCase, skip

import pkg_resources
import numpy as np
from matplotlib.figure import Figure

from nugridpy.data import DataFromTextMixin, NugridData, MesaData
from nugridpy.io import TextFile


DATA_PATH = pkg_resources.resource_filename('tests', 'data/read-write')


class TestAsciiData(unittest.TestCase):
    """Testing data.py ASCII data-handling."""

    def setUp(self):
        """Setting up variables and other test structures."""
        # Set the test data path
        self.data_path = 'tests/data/read-write/'

        # Make a user-facing data object with DataFromTextMixin
        self.test_data = DataFromTextMixin(
            os.path.join(self.data_path,
                         'history.data'), TextFile.MESA_DATA)

    def test_DataFromTextMixin(self):
        """Test the DataFromTextMixin class."""

        # Check the datacol/header names attributes are populated
        self.assertIsInstance(self.test_data.header_names, tuple)
        self.assertGreaterEqual(len(self.test_data.header_names), 0)
        self.assertIsInstance(self.test_data.data_names, tuple)
        self.assertGreater(len(self.test_data.data_names), 0)

    def test_MesaData(self):
        """Test the MesaData class."""
        # Check that MesaData works and reads in a custom filename
        with tempfile.NamedTemporaryFile(prefix='star.log') as tmpf:
            # Write history.data to tempfile
            hd = open(os.path.join(self.data_path, 'history.data'))
            tmpf.write(hd.read().encode())
            hd.close()

            # Read the tempfile as a MesaData instance
            test_alias = MesaData(tempfile.tempdir, history_name=tmpf.name,
                                  history_only=True)

            # Ensure the instance was initialized properly
            self.assertIsInstance(test_alias, MesaData)
            self.assertIsInstance(test_alias.history_data, DataFromTextMixin)

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


class TestHdf5Data(unittest.TestCase):
    """Testing data.py hdf5 handling."""

    def setUp(self):
        """Test setup."""
        # Create data object
        self.test_data = NugridData(DATA_PATH)

        # Define some expected values and names for functional tests
        self.expected_num_cycles = 1000
        self.expected_metadata = ['HDF5_version',
                                  'SE_version',
                                  'age_unit',
                                  'codev',
                                  'dcoeff_unit',
                                  'mass_unit',
                                  'mini',
                                  'modname',
                                  'numcodev',
                                  'overini',
                                  'radius_unit',
                                  'rho_unit',
                                  'rotini',
                                  'temperature_unit',
                                  'zini',
                                  'zisnb']

        self.expected_cycle_headers = ['age',
                                       'deltat',
                                       'model_number',
                                       'shellnb']

        self.expected_cycle_data = ['mass',
                                    'radius',
                                    'rho',
                                    'temperature',
                                    'dcoeff',
                                    'iso_massf',
                                    'elem_massf',
                                    'elem_numf',
                                    'iso_massf_decay',
                                    'elem_massf_decay',
                                    'elem_numf_decay']

    def test_NugridData(self):
        """Testing the NugridData class."""
        # Test data object
        self.assertIsInstance(self.test_data, NugridData)

        # Declare the data object attribute lists
        test_cycles = self.test_data.cycles
        test_metadata = self.test_data.metadata_names
        test_headers = self.test_data.cycle_headers
        test_data_names = self.test_data.cycle_data
        test_ages = self.test_data.ages

        # Ensure they are populated correctly
        self.assertEqual(len(test_cycles), self.expected_num_cycles)
        self.assertEqual(len(test_cycles), len(self.test_data._data))
        self.assertEqual(len(test_cycles), len(self.test_data._headers))
        self.assertEqual(len(test_ages), len(self.test_data._headers))
        self.assertEqual(len(test_metadata), len(self.test_data._metadata))

        # Test datasets
        self.assertIsInstance(self.test_data.A, np.ndarray)
        self.assertIsInstance(self.test_data.Z, np.ndarray)
        self.assertIsInstance(self.test_data.isomeric_state, np.ndarray)

        # Pick a cycle to test
        test_cyc = random.choice(test_cycles)

        # Make sure all data/header key values work to get data
        for key in test_metadata:
            _ = self.test_data.get_metadata(key)
            self.assertIsInstance(_, np.ndarray)
            self.assertTrue(key in self.expected_metadata)

        for key in test_headers:
            _ = self.test_data.get_cycle_header(key, test_cyc)
            self.assertTrue(key in self.expected_cycle_headers)

        for key in test_data_names:
            _ = self.test_data.get_cycle_data(key, test_cyc)
            self.assertIsInstance(_, np.ndarray)
            self.assertTrue(key in self.expected_cycle_data)


class TestMesaData(TestCase):
    """Class for functional testing on the MesaData class."""

    def setUp(self):
        self.data = MesaData(DATA_PATH, history_only=True)

    def test_hdr(self):
        """Checks that we can plot an HRD."""
        fig = self.data.hrd(show=False)
        self.assertIsInstance(fig, Figure)

    def test_kippenhahn(self):
        """Checks that we can plot an Kippenhahn diagramm."""
        fig = self.data.kippenhahn('star_age', show=False)
        self.assertIsInstance(fig, Figure)


@skip('bad data in tests folder')
class TestNugridData(TestCase):
    """Class for functional testing on the NugridData class."""

    def setUp(self):
        self.data = NugridData(DATA_PATH)

    def test_hdr(self):
        """Checks that we can plot an HRD."""
        fig = self.data.hrd(show=False)
        self.assertIsInstance(fig, Figure)

    def test_kippenhahn(self):
        """Checks that we can plot an Kippenhahn diagramm."""
        fig = self.data.kippenhahn('star_age', show=False)
        self.assertIsInstance(fig, Figure)
