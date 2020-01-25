from unittest import TestCase
from unittest.mock import patch

from nugridpy import mesa_data, nugrid_data
from nugridpy.data import MesaData, NugridData

from .fixtures import random_string


class TestUtils(TestCase):
    """Testing some utility functions."""

    @patch.object(MesaData, '__init__', return_value=None)
    def test_mesa_data(self, m_data):
        """Check the proxy function creating MESA data."""

        x = mesa_data(random_string(), random_string())
        self.assertEqual(m_data.call_count, 1)
        self.assertIsInstance(x, MesaData)

    @patch.object(NugridData, '__init__', return_value=None)
    def test_nugrid_data(self, m_hdf5):
        """Check the proxy function creating NuGrid data."""

        x = nugrid_data(random_string())
        self.assertEqual(m_hdf5.call_count, 1)
        self.assertIsInstance(x, NugridData)
