
import os
import random
from unittest import TestCase

import pkg_resources
import numpy as np
from matplotlib.figure import Figure

from nugridpy.data import DataFromTextMixin, DataFromHDF5Mixin, \
    MesaDataHDF5, MesaDataText, \
    NugridData
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
        test_metadata = self.test_data.hattrs
        test_headers = self.test_data.cattrs
        test_data_names = self.test_data.dcols
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
            _ = self.test_data.get_hattr(key)

        for key in test_headers:
            _ = self.test_data.get_cattr(key, test_cyc)

        for key in test_data_names:
            _ = self.test_data.get_dcol(key, test_cyc)


class TestMesaDataText(TestCase):
    """Class for testing the MesaDataText class."""

    def setUp(self):
        path = pkg_resources.resource_filename(
            'nugridpy', os.path.join('resources', 'mesa', 'LOGS'))
        self.data = MesaDataText(path)

    def test_mesa_data_from_text(self):
        """Test the MesaDataText class."""

        # Ensure the instance was initialized properly
        self.assertIsInstance(self.data, MesaDataText)
        self.assertIsInstance(self.data.history_data, DataFromTextMixin)

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


class TestNugridData(TestCase):
    """Class for testing the NugridData class."""

    def setUp(self):
        path_out = pkg_resources.resource_filename(
            'nugridpy', os.path.join('resources', 'nugrid', 'H5_out'))
        self.data_out = NugridData(path_out)

    def test_abu_profile(self):
        """Checks the abu_profile method."""
        fig = self.data_out.abu_profile(
            500, isos=['H-1', 'He-4', 'C-12', 'C-13', 'N-14', 'O-16'], show=False)
        self.assertIsInstance(fig, Figure)

    def test_iso_abund(self):
        """Checks the iso_abund method."""
        fig = self.data_out.iso_abund(500, show=False)
        self.assertIsInstance(fig, Figure)
