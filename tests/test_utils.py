from unittest import TestCase
from unittest.mock import patch

from nugridpy import mesa_data
from nugridpy.data import MesaDataText, MesaDataHDF5

from .fixtures import random_string


class TestUtils(TestCase):
    """Testing some utility functions."""

    @patch.object(MesaDataText, '__init__', return_value=None)
    @patch.object(MesaDataHDF5, '__init__', return_value=None)
    def test_mesa_data(self, m_hdf5, m_text):
        """Check the proxy function creating MESA data."""

        # Text files
        _ = mesa_data(random_string(), data_type='text')
        self.assertEqual(m_text.call_count, 1)

        # HDF5 files
        _ = mesa_data(random_string(), data_type='hdf5')
        self.assertEqual(m_hdf5.call_count, 1)

        # Error
        with self.assertRaises(AssertionError):
            _ = mesa_data(random_string(), data_type='other-type')
