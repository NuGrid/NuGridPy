"""
Tests for the plot module.
"""

from unittest import TestCase
from unittest.mock import patch
import random

from nugridpy.plot import PlotMixin, StarPlotsMixin

from .fixtures import random_string, random_array, random_ints


class TestPlotMixinInit(TestCase):
    """Test the PlotMixin class initialization."""

    def test_get_methods_not_implemented(self):
        """Checks that the correct exception is raised."""

        x = PlotMixin()
        with self.assertRaises(NotImplementedError):
            x.get_cycle_header('toto')
        with self.assertRaises(NotImplementedError):
            x.get_cycle_data('toto', 1)


@patch('nugridpy.plot.plt')
class TestPlotMixin(TestCase):
    """Tests for the PlotMixin class functionalities."""

    def check_plot_call(self, inst, m_get_data, m_plt, cycles=None,
                        num_profiles=1, xlabel=None, ylabel=None, legend=None):
        """Function checking calls for a given configuration."""

        # Reset mocks
        m_get_data.reset_mock()
        m_plt.reset_mock()

        # Names of the profiles requested
        names = [random_string() for _ in range(num_profiles + 1)]  # +1 for x-axis

        # Data returned
        return_value = random_array()
        m_get_data.return_value = return_value

        inst.plot(names[0], names[1:], cycles=cycles, xlabel=xlabel, ylabel=ylabel, legend=legend)

        # Number of plot/get data calls
        # If there is no cycles or one cycle:
        #   * plot is called num_profiles times
        #   * get_data is called (num_profiles + 1) times
        # If there s more than one cycles:
        #   * plot is called num_profiles times*num_cycles
        #   * get_data is called (num_profiles times*num_cycles) + num_profiles times*num_cycles times
        num_cycles = 1
        if cycles:
            if not isinstance(cycles, int):
                num_cycles = len(cycles)
        self.assertEqual(m_plt.plot.call_count, num_profiles * num_cycles)
        self.assertEqual(m_get_data.call_count, (num_profiles + 1) * num_cycles)
        expected_calls = [m_get_data.call_args_list[i][0][0]
                          for i in range(len(m_get_data.call_args_list))]
        self.assertEqual(expected_calls, names * num_cycles)
        if xlabel:
            m_plt.xlabel.assert_called_with(xlabel)
        if ylabel:
            m_plt.ylabel.assert_called_with(ylabel)
        if legend:
            m_plt.legend.assert_called_with(legend)

    def check_plot_call_all_configs(self, get_data, m_plt, cycles=None):
        """Function checking calls for all configurations."""

        with patch.object(PlotMixin, get_data) as m_get_data:

            inst = PlotMixin()

            # Minimal call
            self.check_plot_call(inst, m_get_data, m_plt, cycles)

            # With xlabel
            self.check_plot_call(inst, m_get_data, m_plt, cycles, xlabel=random_string())

            # With ylabel
            self.check_plot_call(inst, m_get_data, m_plt, cycles, ylabel=random_string())

            # With legend
            self.check_plot_call(inst, m_get_data, m_plt, cycles, legend=(random_string(),))

            # Several profiles
            self.check_plot_call(
                inst, m_get_data, m_plt, cycles, num_profiles=random.randint(1, 10))

    def test_plot_header_data(self, m_plt):
        """Checks the plot method on header data."""
        # Header data accross cycles
        self.check_plot_call_all_configs('get_cycle_header', m_plt)

    def test_plot_cycle_data(self, m_plt):
        """Checks the plot method on cycle data."""
        # Cycle data
        # One cycle
        self.check_plot_call_all_configs('get_cycle_data', m_plt, cycles=random.randint(1, 1000))

        # Many cycles
        self.check_plot_call_all_configs('get_cycle_data', m_plt, cycles=random_ints())

    @patch.object(PlotMixin, 'get_cycle_header')
    def test_plot_function_used(self, _m_get_data, m_plt):
        """Checks that the correct plot function is used."""

        inst = PlotMixin()

        # Linear plot
        inst.plot('x', 'y')
        self.assertEqual(m_plt.plot.call_count, 1)
        m_plt.reset_mock()

        # Logx plot
        inst.plot('x', 'y', logx=True)
        self.assertEqual(m_plt.semilogx.call_count, 1)
        m_plt.reset_mock()

        # Logy plot
        inst.plot('x', 'y', logy=True)
        self.assertEqual(m_plt.semilogy.call_count, 1)
        m_plt.reset_mock()

        # Loglog plot
        inst.plot('x', 'y', logx=True, logy=True)
        self.assertEqual(m_plt.loglog.call_count, 1)
        m_plt.reset_mock()


class TestStarPlotsMixin(TestCase):
    """Tests for the StarPlotsMixin class."""

    @patch.object(StarPlotsMixin, 'get_cycle_header')
    @patch('nugridpy.plot.plt')
    def test_hrd(self, m_plt, _m_get_data):
        """Checks the HRD method."""

        _ = StarPlotsMixin().hrd()
        self.assertEqual(m_plt.plot.call_count, 1)

    @patch.object(StarPlotsMixin, '_find_data_by_correct_name')
    @patch.object(StarPlotsMixin, 'get_cycle_header')
    @patch('nugridpy.plot.plt')
    def test_kippenhahn(self, m_plt, m_get_data, m_find_data):
        """Checks the Kippnehahn diagramm."""

        # All arrays returned must have the same size
        r_array = random_array()
        m_get_data.return_value = r_array

        def side_effect(*args):
            return random_string(), r_array
        m_find_data.side_effect = side_effect

        _ = StarPlotsMixin().kippenhahn(random_string())

        # Number of calls given by the various class attributes
        num_calls = len(StarPlotsMixin.CORE_NAMES) + \
            len(StarPlotsMixin.MIXING_BOUNDARIES_NAMES) + 1  # +1 for mass

        self.assertEqual(m_get_data.call_count, 1)  # x-data
        self.assertEqual(m_find_data.call_count, num_calls)  # y-data
        self.assertEqual(m_plt.plot.call_count, num_calls)

    @patch.object(StarPlotsMixin, 'get_cycle_header')
    def test_find_data_by_correct_name(self, m_get_data):
        """Checks the function finding data by name."""

        field_name = 'my_field'

        def side_effect(name):
            if name == field_name:
                return random_array()
            raise KeyError

        m_get_data.side_effect = side_effect

        # Empty tuple
        d = ()
        self.assertFalse(StarPlotsMixin()._find_data_by_correct_name(d))

        # Key not present
        d = ('a', 'b', 'c', 'y', 'z')
        self.assertFalse(StarPlotsMixin()._find_data_by_correct_name(d))

        # Key present
        d += (field_name,)
        r = StarPlotsMixin()._find_data_by_correct_name(d)
        self.assertEqual(len(r), 2)
        self.assertEqual(r[0], field_name)
