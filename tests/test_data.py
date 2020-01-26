

from unittest import TestCase
from unittest.mock import patch, Mock, PropertyMock
import random

import numpy as np

from nugridpy.data import \
    DataFromTextMixin, DataFromHDF5Mixin, \
    MesaDataText, MesaDataHDF5, MesaData
from .fixtures import random_string, random_ints, random_array


@patch.object(DataFromTextMixin, 'get', return_value=random_ints())
@patch.object(DataFromTextMixin, '__init__', return_value=None)
class TestMesaDataText(TestCase):
    """Class for testing the MesaDataText class."""

    @patch.object(DataFromTextMixin, 'data_names', new_callable=PropertyMock)
    @patch('os.path.isfile')
    def tests_init(self, m_is_file, m_dn, m_init, m_get):
        """Check the calls during intialization."""

        cycles = np.array(sorted(set(random_ints())))
        file_numbers = random_ints(len(cycles))

        def side_effect(name):
            if name == 'model_number':
                return cycles
            return file_numbers
        m_get.side_effect = side_effect

        # 1 - No file present
        m_is_file.return_value = False
        x = MesaDataText(random_string())
        self.assertEqual(m_init.call_count, 2)  # history data + profile index

        # getting the cycle ('model_number'), 'model_number', and 'log_file_number'
        self.assertEqual(m_get.call_count, 3)
        self.assertEqual(m_get.call_args_list[0][0][0], 'model_number')
        self.assertEqual(m_get.call_args_list[1][0][0], 'model_number')
        self.assertEqual(m_get.call_args_list[2][0][0], 'log_file_number')

        # Check the mapping cycle/file number
        self.assertEqual(m_is_file.call_count, len(cycles))
        self.assertListEqual(x.cycles_with_data, [])
        self.assertFalse(m_dn.called)

        # 2 - At least one file there
        m_is_file.reset_mock()
        is_file_present = [random.choice([True, False]) for _ in range(len(cycles))]
        is_file_present[random.choice(range(len(cycles)))] = True
        m_is_file.side_effect = is_file_present
        x = MesaDataText(random_string())

        # Check the mapping cycle/file number
        self.assertEqual(m_is_file.call_count, len(cycles))
        self.assertSetEqual(set(x.cycles_with_data), set(cycles[is_file_present]))
        self.assertTrue(m_dn.called)

    def test_header_methods(self, _m_init, m_get):
        """Check the methods and properties for headers."""

        x = MesaDataText(random_string())
        m_get.reset_mock()

        # Cycles
        self.assertEqual(x.cycles, sorted(x.cycles))
        self.assertEqual(m_get.call_count, 2)

        # Header attributes
        with patch.object(
                DataFromTextMixin, 'header_names',
                new_callable=PropertyMock) as m_hn:
            _ = x.hattrs
            self.assertEqual(m_hn.call_count, 1)

    def test_cycle_methods(self, _m_init, m_get):
        """Check the methods and properties for cycles."""

        x = MesaDataText(random_string())
        m_get.reset_mock()

        # Cycle attributes
        with patch.object(
                DataFromTextMixin, 'data_names',
                new_callable=PropertyMock) as m_dn:
            _ = x.cattrs
            self.assertEqual(m_dn.call_count, 1)

        # Only all cycles
        r_str = random_string()
        x.get_cattr(r_str)
        self.assertEqual(m_get.call_count, 1)
        self.assertTrue(m_get.called_with(r_str))

    def test_data_columns_methods(self, _m_init, m_get):
        """Check the methods and properties for cycle data columns."""

        x = MesaDataText(random_string())
        m_get.reset_mock()

        r_str = random_string()

        # Cycle does not exist
        with self.assertRaises(RuntimeError):
            x.get_dcol(r_str, 1)
            self.assertEqual(m_get.call_count, 0)

        # One cycles
        x.cycles_with_data = [1]
        with patch.object(MesaDataText, '_get_cycle_data_filepath') as m_get_path:
            m_get_path.return_value = random_string()
            x.get_dcol(r_str, 1)
            self.assertEqual(m_get_path.call_count, 1)
            self.assertEqual(m_get.call_count, 1)
            self.assertTrue(m_get.called_with(r_str))

    @patch.object(MesaDataText, '_get_cycle_data_filepath')
    def test_cache(self, m_get_path, _m_init, m_get):
        """Check caching of the profile data."""

        x = MesaDataText(random_string())
        m_get.reset_mock()
        m_get_path.reset_mock()
        x.cycles_with_data = [1]
        r_str, r_str2 = random_string(), random_string()

        # First call
        x.get_dcol(r_str, 1)
        self.assertEqual(m_get_path.call_count, 1)

        # Call again - cached
        m_get_path.reset_mock()
        m_get.reset_mock()
        x.get_dcol(r_str, 1)
        self.assertEqual(m_get_path.call_count, 0)

        # Call again but with different field - also cached
        m_get_path.reset_mock()
        m_get.reset_mock()
        x.get_dcol(r_str2, 1)
        self.assertEqual(m_get_path.call_count, 0)

        # Clear cache
        x.clear_cache()

        # Call again - cache has been cleared
        m_get_path.reset_mock()
        m_get.reset_mock()
        x.get_dcol(r_str, 1)
        self.assertEqual(m_get_path.call_count, 1)


@patch.object(DataFromHDF5Mixin, '__init__', return_value=None)
class TestDataFromHDF5Mixin(TestCase):
    """Class for testing the DataFromHDF5Mixin class."""

    def test_header_methods(self, _m_init):
        """Check the methods and properties for headers."""

        x = DataFromHDF5Mixin(random_string())

        # Cycles
        r_int = set(random_ints())
        x._data = {i: None for i in r_int}
        self.assertEqual(x.cycles, sorted(r_int))

        # Header attributes
        # Can be string or numeric, and must be array of one element
        r_str = [random_string() for _ in range(random_ints(1))]
        x._metadata = {
            s: np.array([random.choice([random_ints(1), random_string()])])
            for s in r_str
        }
        self.assertSetEqual(set(r_str), set(x.hattrs))

        for s in r_str:
            self.assertEqual(x.get_hattr(s), x._metadata[s][0])

    def test_cycle_methods(self, _m_init):
        """Check the methods and properties for cycles."""

        x = DataFromHDF5Mixin(random_string())

        # Cycles
        r_int = set(random_ints())
        x._data = {i: None for i in r_int}
        self.assertEqual(x.cycles, sorted(r_int))

        # Cycle attributes
        # Can be string or numeric, and must be array of one element
        r_str = [random_string() for _ in range(random_ints(1))]
        cycle = x.cycles[0]
        x._headers = {}
        x._headers[cycle] = {
            s: np.array([random.choice([random_ints(1), random_string()])])
            for s in r_str
        }
        self.assertSetEqual(set(r_str), set(x.cattrs))

        with patch.object(DataFromHDF5Mixin, '_get_cattr') as m_get_cattr:
            name = random_string()
            # 1 cycle
            x.get_cattr(name, 1)
            self.assertEqual(m_get_cattr.call_count, 1)
            m_get_cattr.reset_mock()

            # Many cycles
            cycles = list(random_ints())
            x.get_cattr(name, cycles)
            self.assertEqual(m_get_cattr.call_count, len(cycles))
            m_get_cattr.reset_mock()

            # All cycles
            x.get_cattr(name)
            self.assertEqual(m_get_cattr.call_count, len(x.cycles))

    def test_data_columns_methods(self, _m_init):
        """Check the methods and properties for cycle data columns."""

        x = DataFromHDF5Mixin(random_string())

        # Cycles
        r_int = set(random_ints())
        x._data = {i: None for i in r_int}
        self.assertEqual(x.cycles, sorted(r_int))

        # Data columns names attributes
        # Can be string or numeric, and must be array of one element
        r_str = [random_string() for _ in range(random_ints(1))]
        cycle = x.cycles[0]
        # Ugly but we dont want to build an HDF5 dataset here
        x._data[cycle] = Mock()
        x._data[cycle].dtype.names = r_str
        self.assertEqual(sorted(r_str), x.dcols)

        with patch.object(DataFromHDF5Mixin, '_get_dcol') as m_get_dcol:
            name = random_string()
            # 1 cycle
            x.get_dcol(name, 1)
            self.assertEqual(m_get_dcol.call_count, 1)
            m_get_dcol.reset_mock()

            # Many cycles
            cycles = list(random_ints())
            x.get_dcol(name, cycles)
            self.assertEqual(m_get_dcol.call_count, len(cycles))
            m_get_dcol.reset_mock()

    @patch.object(DataFromHDF5Mixin, '_get_cycle')
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


class TestMesaData(TestCase):
    """Class for testing the MesaData proxy class."""

    def test_proxy_bad_data(self):
        """Check the behavior when bad data type is provided."""

        with self.assertRaises(RuntimeError):
            _ = MesaData('toto')

    @patch.object(MesaDataText, '__init__', return_value=None)
    def test_proxy_text_data(self, m_init):
        """Check the proxy for text data."""

        x = MesaData('text')
        self.assertIsInstance(x._proxied_obj, MesaDataText)
        self.assertEqual(m_init.call_count, 1)

    @patch.object(MesaDataHDF5, '__init__', return_value=None)
    def test_proxy_hdf5_data(self, m_init):
        """Check the proxy for HDF5 data."""

        x = MesaData('hdf5', random_string())
        self.assertIsInstance(x._proxied_obj, MesaDataHDF5)
        self.assertEqual(m_init.call_count, 1)

    @patch.object(MesaDataText, '__init__', return_value=None)
    @patch.object(MesaDataText, '__getattr__', create=True)
    def test_redirection(self, m_get, _m_init):
        """Check that any attribute/method is redirected to the proxied object."""

        x = MesaData('text')
        _ = x.random_att
        m_get.assert_called_with('random_att')
        _ = x.random_method()
        m_get.assert_called_with('random_method')
