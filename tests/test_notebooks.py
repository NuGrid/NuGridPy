import os
import pkg_resources
import random
from tempfile import TemporaryDirectory as tmpdir
from unittest import TestCase, skip
from unittest.mock import MagicMock, patch

from nugridpy import ascii_table, h5T, mesa, nugridse, ppn, utils
from .fixtures import random_string


DATA_PATH = pkg_resources.resource_filename('tests', 'data/read-write')


class TestStarExplore(TestCase):
    """Test the functionality of the star_explore notebook."""
    
    def setUp(self):
        """Set up the tests."""
        # Set up mesa star_log object and associated class methods
        self.mesa_in = mesa.star_log(DATA_PATH)
        self.kip_CO = self.mesa_in.kippenhahn_CO
        self.tcrhoc = self.mesa_in.tcrhoc
        self.hrd_new = self.mesa_in.hrd_new
        self.kip_cont = self.mesa_in.kip_cont

        # Set up nugridse se object and associated class methods
        self.ngse_in = nugridse.se(DATA_PATH)
        self.ngse_se = self.ngse_in.se
        self.ngse_hattrs = self.ngse_se.hattrs
        self.ngse_cattrs = self.ngse_se.cattrs
        self.ngse_dcols = self.ngse_se.dcols

        # Set up se object plot methods for testing
        self.ngse_plot = self.ngse_in.plot
        self.ngse_prof = self.ngse_in.abu_profile
        self.ngse_iso = self.ngse_in.iso_abund
        self.ngse_chart = self.ngse_in.abu_chart
        self.ngse_species = ['H-1','He-4','C-12','C-13','N-14','O-16']

    @patch('nugridpy.mesa.pyl.plot')
    @patch('nugridpy.mesa.pl.plot', return_value=(MagicMock(),))
    @patch('nugridpy.mesa.pl.axes')
    def test_mesa(self, mocked_ax, mocked_pl, mocked_pyl):
        """Test the mesa object as used in star_explore."""
        # Test the mesa.kippenhahn_CO() method
        _ = self.kip_CO(random.randint(1, 1000), 'model')
        self.assertEqual(mocked_pyl.call_count, 8)
        mocked_pyl.reset_mock()

        _ = self.kip_CO(random.randint(1, 1000), 'time')
        self.assertEqual(mocked_pyl.call_count, 8)
        mocked_pyl.reset_mock()

        with self.assertRaises(UnboundLocalError):
            _ = self.kip_CO(random.randint(1, 1000), random_string())

        with self.assertRaisesRegex(TypeError, 'kippenhahn_CO()'):
            _ = self.kip_CO()

        # Test the mesa.tcrhoc() method
        _ = self.tcrhoc()
        self.assertEqual(mocked_pl.call_count, 1)

        # Test the mesa.hrd_new() method
        _ = self.hrd_new()
        self.assertEqual(mocked_pyl.call_count, 1)

        # Test the mesa.kip_cont() method
        with tmpdir() as tdir:
            _ = self.kip_cont(outfile='{}f.png'.format(tdir), showfig=False)
            self.assertEqual(mocked_ax.call_count, 1)

    @patch('nugridpy.mesa.pl.plot', return_value=(MagicMock(),))
    @patch('nugridpy.mesa.pl.axes')
    def test_se(self, mocked_ax, mocked_pl):
        """Test the nugridse.se object as used in star_explore."""
        # Test the se object
        self.assertIsInstance(self.ngse_se, h5T.Files)
        
        # Test the se.hattrs object and get method
        self.assertIsInstance(self.ngse_hattrs, list)

        for hattr in self.ngse_hattrs:
            _ = self.ngse_se.get(hattr)

        # Test the se.cattrs object and get method
        self.assertIsInstance(self.ngse_cattrs, list)
        #import ipdb
        #ipdb.set_trace()

        for cattr in self.ngse_cattrs:
            _ = self.ngse_se.get(cattr)

        # Test the se.dcols object and get method
        self.assertIsInstance(self.ngse_dcols, list)

        for dcol in self.ngse_dcols:
            _ = self.ngse_se.get(dcol)

        # Test the se object plot method
        with self.assertRaisesRegex(TypeError, 'plot()'):
            _ = self.ngse_plot()

        # Plot cattrs
        for attr in self.ngse_cattrs:
            _ = self.ngse_plot(attr, attr, show=False)
            self.assertEqual(mocked_pl.call_count, 1)
            mocked_pl.reset_mock()

        # Plot first 4 dcols (others not plot quantities)
        for attr in self.ngse_dcols[:4]:
            _ = self.ngse_plot(attr, attr, show=False)
            self.assertEqual(mocked_pl.call_count, 1)
            mocked_pl.reset_mock()

        # Test the use of the abu_profile method
        _ = self.ngse_prof(fname=random.randint(1, 1000),
                           isos=self.ngse_species)
        self.assertEqual(mocked_pl.call_count, 6)
        mocked_pl.reset_mock()

        with self.assertRaisesRegex(OSError, 
                    'Please provide the cycle number fname'):
            _ = self.ngse_prof()

        # Test the use of the iso_abund method
        with self.assertRaisesRegex(TypeError, 'iso_abund()'):
            _ = self.ngse_iso(show=False)

        # Test the use of the abu_chart method
        with self.assertRaisesRegex(TypeError, 'abu_chart()'):
            _ = self.ngse_chart(show=False)


class TestWeakIProcessTemplate(TestCase):
    """Test the analyze_template_run notebook."""

    def setUp(self):
        """Set up the notebook testing."""
        # Set up ppn abu vector and dcols
        self.abu_vec = ppn.abu_vector(DATA_PATH)
        self.abu_dcols = self.abu_vec.dcols

        # Set up test methods from notebook
        self.el_abu = self.abu_vec.elemental_abund
        self.iso_abu = self.abu_vec.iso_abund

    @patch('nugridpy.mesa.pl.plot', return_value=(MagicMock(),))
    def test_elem_abu(self, mocked_pl):
        """Test the elemental_abund method as used in the notebook."""
        # Test the method itself (single-cycle read-in so index is 0)
        _ = self.el_abu(0)
        self.assertEqual(mocked_pl.call_count, 1)

        with self.assertRaisesRegex(TypeError, 'elemental_abund()'):
            _ = self.el_abu()

        # Test the method attributes (6 dcols)
        self.assertIsInstance(self.abu_dcols, list)
        self.assertEqual(len(self.abu_dcols), 6)

        # Get attributes and ensure array lengths are equal
        dcols = [self.abu_vec.get(dcol)[0] for dcol in self.abu_dcols]
        length = len(dcols[0])
        self.assertTrue(all(len(dcol) == length for dcol in dcols))

    @patch('nugridpy.mesa.pl.plot', return_value=(MagicMock(),))
    def test_iso_abund(self, mocked_pl):
        """Test the iso_abund method as used in the notebook."""
        # Test the method call (1 cycle, index 0)
        _ = self.iso_abu(0)
        self.assertEqual(mocked_pl.call_count, 85)
        
        with self.assertRaisesRegex(TypeError, 'iso_abund()'):
            _ = self.iso_abu()


class TestExtractTraj(TestCase):
    """Test the Extract_trajectory notebook."""

    def setUp(self):
        """Set up the trajectory notebook tests."""
        # Initialize nugridse instance and methods
        self.se = nugridse.se(DATA_PATH)
        self.traj = nugridse.trajectory

    @skip('This functionality is, at the moment, hopelessly broken.')
    def test_trajectory(self):
        """Test the trajectory write method."""
        nugridse.set_data_path(DATA_PATH)
        self.se_in = nugridse.se(mass=12, Z=0.006)
        self.m_coor = 0.15
        self.model_start = 235
        self.r, self.rho, self.temp, self.time = nugridse.trajectory(
                self.model_start, 1306, 10, self.mcoor, age_in_sec=True,
                online=False)
