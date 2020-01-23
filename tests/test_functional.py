
import os
import random
from unittest import TestCase

import pkg_resources
import numpy as np
from matplotlib.figure import Figure

from nugridpy.data import DataFromTextMixin, DataFromHDF5Mixin, MesaDataHDF5, MesaDataText
from nugridpy.io import TextFile


class TestDataFromTextMixin(TestCase):
    """Testing the DataFromTextMixin class."""

    def test_data_from_text(self):
        """Test the DataFromTextMixin class."""
        # Set the test data path
        self.data_path = 'tests/data/read-write/'

        # Make a user-facing data object with DataFromTextMixin
        self.test_data = DataFromTextMixin(
            os.path.join(self.data_path,
                         'history.data'), TextFile.MESA_DATA)

        # Check the datacol/header names attributes are populated
        self.assertIsInstance(self.test_data.header_names, tuple)
        self.assertGreaterEqual(len(self.test_data.header_names), 0)
        self.assertIsInstance(self.test_data.data_names, tuple)
        self.assertGreater(len(self.test_data.data_names), 0)


class TestDataFromHDF5Mixin(TestCase):
    """Testing the DataFromHDF5Mixin class."""

    def test_data_from_hdf5(self):
        """Testing the DataFromHDF5Mixin class."""

        # Create data object
        path = pkg_resources.resource_filename(
            'nugridpy', os.path.join('resources', 'mesa', 'HDF5'))
        self.test_data = DataFromHDF5Mixin(path)

        # Declare the data object attribute lists
        test_cycles = self.test_data.cycles
        test_metadata = self.test_data.metadata_names
        test_headers = self.test_data.cycle_headers
        test_data_names = self.test_data.cycle_data
        test_ages = self.test_data.ages

        # Ensure they are populated correctly
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

        for key in test_headers:
            _ = self.test_data.get_cycle_header(key, test_cyc)

        for key in test_data_names:
            _ = self.test_data.get_cycle_data(key, test_cyc)
            self.assertIsInstance(_, np.ndarray)


class TestMesaDataText(TestCase):
    """Class for testing the MesaDataText class."""

    def setUp(self):
        path = pkg_resources.resource_filename(
            'nugridpy', os.path.join('resources', 'mesa', 'LOGS'))
        self.data = MesaDataText(path, history_only=True)

    def test_mesa_data_from_text(self):
        """Test the MesaDataText class."""

        # Ensure the instance was initialized properly
        self.assertIsInstance(self.data, MesaDataText)
        self.assertIsInstance(self.data.history_data, DataFromTextMixin)

        # Ensure history_only worked properly, no profile data
        with self.assertRaises(AttributeError):
            _ = self.data.profiles_index
        with self.assertRaises(AttributeError):
            _ = self.data.profiles

        # Assign data to variables
        test_headers = self.data.history_data.header_names
        test_data = self.data.history_data.data_names

        # Ensure the history.data read worked
        self.assertIsInstance(test_headers, tuple)
        self.assertGreater(len(test_headers), 0)
        self.assertIsInstance(test_data, tuple)
        self.assertGreater(len(test_data), 0)

    def test_hrd(self):
        """Checks that we can plot an HRD."""
        fig = self.data.hrd(show=False)
        self.assertIsInstance(fig, Figure)

    def test_kippenhahn(self):
        """Checks that we can plot an Kippenhahn diagramm."""
        fig = self.data.kippenhahn('star_age', show=False)
        self.assertIsInstance(fig, Figure)


class TestMesaDataHDF5(TestCase):
    """Class for testing the MesaDataHDF5 class."""

    def setUp(self):
        path = pkg_resources.resource_filename(
            'nugridpy', os.path.join('resources', 'mesa', 'HDF5'))
        self.data = MesaDataHDF5(path)

    def test_hdr(self):
        """Checks that we can plot an HRD."""
        fig = self.data.hrd(show=False)
        self.assertIsInstance(fig, Figure)

    def test_kippenhahn(self):
        """Checks that we can plot an Kippenhahn diagramm."""
        fig = self.data.kippenhahn('age', show=False)
        self.assertIsInstance(fig, Figure)
