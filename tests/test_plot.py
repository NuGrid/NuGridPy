"""
Tests for the plot module.
"""

from unittest import TestCase
from unittest.mock import patch
import random

from nugridpy.plot import PlotMixin, MESAPlotMixin, NugridPlotMixin

from .fixtures import random_string, random_array, random_ints


class TestPlotMixinInit(TestCase):
    """Test the PlotMixin class initialization."""

    def test_get_methods_not_implemented(self):
        """Checks that the correct exception is raised."""

        x = PlotMixin()
        with self.assertRaises(NotImplementedError):
            x.get_cattr('toto')
        with self.assertRaises(NotImplementedError):
            x.get_dcol('toto', 1)


@patch('nugridpy.plot.plt')
@patch('nugridpy.plot.np.where', return_value=...)
class TestPlotMixin(TestCase):
    """Tests for the PlotMixin class functionalities."""

    def check_plot_call(self, inst, m_get_data, m_plt, m_where, cycles=None, x0=None,
                        num_profiles=1, xlabel=None, ylabel=None, legend=True):
        """Function checking calls for a given configuration."""

        # Reset mocks
        m_get_data.reset_mock()
        m_plt.reset_mock()
        m_where.reset_mock()

        # Names of the profiles requested
        names = [random_string() for _ in range(num_profiles + 1)]  # +1 for x-axis

        # Data returned
        return_value = random_array()
        m_get_data.return_value = return_value

        inst.plot(names[0], names[1:], cycles=cycles, x0=x0,
                  xlabel=xlabel, ylabel=ylabel, legend=legend)

        # Number of plot/get data calls
        # If there is no cycles or one cycle:
        #   * plot is called num_profiles times
        #   * get_data is called (num_profiles + 1) times
        # If there s more than one cycles:
        #   * plot is called num_profiles times*num_cycles
        #   * get_data is called (num_profiles times*num_cycles) + num_profiles times*num_cycles times
        num_cycles = 1
        if cycles is not None:
            if not isinstance(cycles, int):
                num_cycles = len(cycles)
        self.assertEqual(m_plt.plot.call_count, num_profiles * num_cycles)
        self.assertEqual(m_get_data.call_count, (num_profiles + 1) * num_cycles)
        expected_calls = [m_get_data.call_args_list[i][0][0]
                          for i in range(len(m_get_data.call_args_list))]
        self.assertEqual(expected_calls, names * num_cycles)
        if x0:
            self.assertTrue(m_where.called)
        if xlabel:
            m_plt.xlabel.assert_called_with(xlabel)
        if ylabel:
            m_plt.ylabel.assert_called_with(ylabel)
        if legend:
            self.assertTrue(m_plt.legend.called)

    def check_plot_call_all_configs(self, get_data, m_plt, m_where, cycles=None):
        """Function checking calls for all configurations."""

        with patch.object(PlotMixin, get_data) as m_get_data:

            inst = PlotMixin()

            # Minimal call
            self.check_plot_call(inst, m_get_data, m_plt, m_where, cycles)

            # With xlabel
            self.check_plot_call(inst, m_get_data, m_plt, m_where, cycles, xlabel=random_string())

            # With ylabel
            self.check_plot_call(inst, m_get_data, m_plt, m_where, cycles, ylabel=random_string())

            # Without legend
            self.check_plot_call(
                inst, m_get_data, m_plt, m_where, cycles, legend=False)

            # With slicing the data
            self.check_plot_call(inst, m_get_data, m_plt, m_where, cycles, x0=random.randint(1, 10))

            # Several profiles
            self.check_plot_call(
                inst, m_get_data, m_plt, m_where, cycles, num_profiles=random.randint(1, 10))

    def test_plot_header_data(self, m_where, m_plt):
        """Checks the plot method on header data."""
        # Header data accross cycles
        self.check_plot_call_all_configs('get_cattr', m_plt, m_where)

    def test_plot_cycle_data(self, m_where, m_plt):
        """Checks the plot method on cycle data."""
        # Cycle data
        # One cycle
        self.check_plot_call_all_configs(
            'get_dcol', m_plt, m_where, cycles=random.randint(1, 1000))

        # Many cycles
        self.check_plot_call_all_configs('get_dcol', m_plt, m_where, cycles=random_ints())

    @patch.object(PlotMixin, 'get_cattr')
    def test_plot_function_used(self, _m_get_data, _m_where, m_plt):
        """Checks that the correct plot function is used."""

        inst = PlotMixin()

        # Linear plot
        inst.plot('x', 'y')
        self.assertEqual(m_plt.plot.call_count, 1)

        # Logx plot
        inst.plot('x', 'y', logx=True)
        self.assertEqual(m_plt.semilogx.call_count, 1)

        # Logy plot
        inst.plot('x', 'y', logy=True)
        self.assertEqual(m_plt.semilogy.call_count, 1)

        # Loglog plot
        inst.plot('x', 'y', logx=True, logy=True)
        self.assertEqual(m_plt.loglog.call_count, 1)

    @patch.object(PlotMixin, 'get_cattr')
    def test_find_data_by_correct_name(self, m_get_data, _m_where, _m_plt):
        """Checks the function finding data by name."""

        field_name = 'my_field'

        def side_effect(name):
            if name == field_name:
                return random_array()
            raise KeyError

        m_get_data.side_effect = side_effect

        # Empty tuple
        d = ()
        self.assertFalse(MESAPlotMixin()._find_data_by_correct_name(d))

        # Key not present
        d = ('a', 'b', 'c', 'y', 'z')
        self.assertFalse(MESAPlotMixin()._find_data_by_correct_name(d))

        # Key present
        d += (field_name,)
        r = MESAPlotMixin()._find_data_by_correct_name(d)
        self.assertEqual(len(r), 2)
        self.assertEqual(r[0], field_name)

        # Key present but we dont want the data
        d += (field_name,)
        r = MESAPlotMixin()._find_data_by_correct_name(d, return_data=False)
        self.assertEqual(r, field_name)


class TestMESAPlotMixin(TestCase):
    """Tests for the MESAPlotMixin class."""

    @patch.object(MESAPlotMixin, 'get_cattr')
    @patch.object(MESAPlotMixin, 'plot')
    def test_hrd(self, m_plt, m_get_data):
        """Checks the HRD method."""

        x = MESAPlotMixin()
        x.hrd()
        for arg in m_get_data.call_args_list:
            self.assertIn(arg[0][0], x.LOG_TEFF_NAMES + x.LOG_L_NAMES)
        self.assertEqual(m_plt.call_count, 1)

    @patch.object(MESAPlotMixin, 'get_cattr')
    @patch.object(MESAPlotMixin, 'plot')
    def test_tcrhoc(self, m_plt, m_get_data):
        """Checks the Tc/Rhoc method."""

        x = MESAPlotMixin()
        x.tcrhoc()
        for arg in m_get_data.call_args_list:
            self.assertIn(arg[0][0], x.LOG_RHOC_NAMES + x.LOG_TC_NAMES)
        self.assertEqual(m_plt.call_count, 1)

    @patch.object(MESAPlotMixin, '_find_data_by_correct_name')
    @patch.object(MESAPlotMixin, 'get_cattr')
    @patch('nugridpy.plot.plt')
    def test_kippenhahn(self, m_plt, m_get_data, m_find_data):
        """Checks the Kippnehahn diagramm."""

        # All arrays returned must have the same size
        r_array = random_array()
        m_get_data.return_value = r_array

        def side_effect(*_args):
            return random_string(), r_array
        m_find_data.side_effect = side_effect

        _ = MESAPlotMixin().kippenhahn(random_string())

        # Number of calls given by the various class attributes
        num_calls = len(MESAPlotMixin.CORE_NAMES) + \
            len(MESAPlotMixin.MIXING_BOUNDARIES_NAMES) + 1  # +1 for mass

        self.assertEqual(m_get_data.call_count, 1)  # x-data
        self.assertEqual(m_find_data.call_count, num_calls)  # y-data
        self.assertEqual(m_plt.plot.call_count, num_calls)

        # With C/O ration
        m_get_data.reset_mock()
        m_find_data.reset_mock()
        m_plt.plot.reset_mock()

        _ = MESAPlotMixin().kippenhahn(random_string(), CO_ratio=True)

        self.assertEqual(m_get_data.call_count, 1)  # x-data
        self.assertEqual(m_find_data.call_count, num_calls + 2)  # y-data
        self.assertEqual(m_plt.plot.call_count, num_calls + 1)
        self.assertEqual(m_plt.twinx.call_count, 1)  # second y-axis


class TestNugridPlotMixin(TestCase):
    """Tests for the NugridPlotMixin class."""

    def setUp(self):

        self.x = NugridPlotMixin()
        self.cycle = random_ints(1)
        self.isotopes = ['-'.join([random_string(), random_string()])
                         for _ in range(random_ints(1))]

    @patch.object(NugridPlotMixin, 'plot')
    def test_abu_profile(self, m_plt):
        """Checks the abu_profile method."""
        self.x.abu_profile(self.cycle, isos=self.isotopes)

        self.assertEqual(m_plt.call_count, 1)
        self.assertEqual(m_plt.call_args_list[0][0], ('mass', self.isotopes))
        self.assertEqual(m_plt.call_args_list[0][1]['cycles'], self.cycle)

    @patch('nugridpy.plot.plt')
    @patch('nugridpy.plot.get_A_Z')
    @patch.object(NugridPlotMixin, 'get_dcol')
    def test_iso_abund(self, m_get_data, m_get_A_Z, m_plt):
        """Checks the iso_abund method."""

        # Build isotopes to see if they are correctly excluded from the plot
        self.x.isotopes = self.isotopes

        isotopes_to_plot = set([
            random.choice(self.isotopes)
            for _ in range(random.randint(0, len(self.isotopes)))])

        # Exclusion by mass number
        bounds = random_ints(2)
        a_min = min(bounds)
        a_max = max(bounds)

        def side_effect_A_Z(name):
            # Mass number should be in the allowed range
            if name in isotopes_to_plot:
                a = random.randint(a_min, a_max)
            else:  # Otherwishe outside the bounds
                a = random.choice(
                    [random.randint(0, a_min - 1),
                     random.randint(a_max + 1, 1001)]
                )
            return a, random_ints(1), name

        m_get_A_Z.side_effect = side_effect_A_Z
        m_get_data.return_value = [1.0]

        self.x.iso_abund(self.cycle, a_min=a_min, a_max=a_max)
        self.assertEqual(m_plt.plot.call_count, len(isotopes_to_plot))

        # Exclusion by threshold
        m_plt.reset_mock()
        low_val, threshold, high_val = tuple(sorted(random_array(3)))

        def side_effect_get(name, _cycle):
            # Mass number should be in the allowed range
            if name in isotopes_to_plot:
                return [high_val]
            return [low_val]

        m_get_A_Z.side_effect = lambda name: (a_min, random_ints(1), name)
        m_get_data.side_effect = side_effect_get

        self.x.iso_abund(self.cycle, a_min=a_min, a_max=a_max, threshold=threshold)
        self.assertEqual(m_plt.plot.call_count, len(isotopes_to_plot))
