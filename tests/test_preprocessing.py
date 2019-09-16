import unittest
import os
from mock import patch

from nugridpy import nugridse as nuse
from nugridpy.ascii_table import readTable


class TestPreprocessing(unittest.TestCase):
    """Class for testing the pre-processing capability of the framework."""

    def setUp(self):

        preproc_filename = 'h5preproc.txt'

        # Get location of HDF5 files for testing
        current_dir = os.path.dirname(os.path.abspath(__file__))
        self.data_dir = os.path.join(current_dir, 'data')

        # Check that there is no pre-processing file
        self.preproc_absname = os.path.join(self.data_dir, preproc_filename)
        self.assertFalse(os.path.isfile(self.preproc_absname))

    @patch('nugridpy.ascii_table.write')
    def test_no_writing_preproc_file(self, mocked_write):
        """Check that writing the preproc file has been disabled."""

        _ = nuse.se(self.data_dir)
        self.assertFalse(mocked_write.called)

    @patch.object(readTable, '__init__')
    def test_no_reading_preproc_file(self, mocked_read):
        """Check that reading the preproc file has been disabled."""

        def delete_tmp_file():
            '''Remove pre-processed file.'''

            if os.path.isfile(self.preproc_absname):
                os.remove(self.preproc_absname)

        try:
            f = open(self.preproc_absname, 'w')
            f.write('Some random things in that file! ')
            f.close()

            _ = nuse.se(self.data_dir)
            self.assertFalse(mocked_read.called)

            delete_tmp_file()

        except AssertionError as err:
            delete_tmp_file()
            raise err
