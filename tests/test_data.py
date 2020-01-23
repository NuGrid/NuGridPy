

from unittest import TestCase
from unittest.mock import patch

import numpy as np

from nugridpy.data import DataFromTextMixin, DataFromHDF5Mixin, MesaDataText
from .fixtures import random_string, random_ints, random_array


@patch.object(DataFromTextMixin, '__init__', return_value=None)
@patch.object(DataFromTextMixin, 'get', return_value=random_ints())
class TestMesaDataText(TestCase):
    """Class for testing the MesaDataText class."""

    def test_cycles(self, m_get, _m_init):
        """Check the cycles property."""

        x = MesaDataText(random_string(), history_only=True)
        np.testing.assert_array_equal(x.cycles, m_get.return_value)

    def test_get_cycle_header(self, m_get, _m_init):
        """Check that we use the mapping cycle/index when requesting header."""

        x = MesaDataText(random_string(), history_only=True)

        # Check _cycle_index_mapping
        # We build the dict in a dummy manner on purpose to not rely on np.ndenumerate
        expected_dict = {}
        i = 0
        for c in m_get.return_value:
            expected_dict[c] = i
            i += 1
        self.assertDictEqual(x._cycle_index_mapping, expected_dict)


@patch.object(DataFromHDF5Mixin, '__init__', return_value=None)
@patch.object(DataFromHDF5Mixin, '_get_cycle')
class TestDataFromHDF5Mixin(TestCase):
    """Class for testing the DataFromHDF5Mixin class."""

    def test_get_isotope_mass_fraction(self, _m_get_cycle, _m_init):
        """Check the _get_isotope_mass_fraction function."""

        x = DataFromHDF5Mixin(random_string())
        random_integers = random_ints(2)
        num_cycles = max(random_integers)
        cycle = min(random_integers)

        with self.assertRaises(KeyError):
            x._get_isotope_mass_fraction('toto', cycle)

        # Assuming bad A,Z data
        x.A = np.array([1, 1, 2])
        x.Z = np.array([1, 1, 2])

        with self.assertRaises(KeyError):
            x._get_isotope_mass_fraction('H-1', cycle)

        # OK
        x.A = np.array([1, 2, 3])
        _m_get_cycle.return_value = {
            x.MASS_FRACTION_HDF5_MEMBER_NAME: random_array(num_cycles, 20)
        }
        x._get_isotope_mass_fraction('H-1', cycle)
