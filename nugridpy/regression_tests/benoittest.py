from __future__ import print_function
#import matplotlib
#matplotlib.use('agg')
import unittest
import sys

#from tempdir.tempfile_ import TemporaryDirectory

class Test_README_File(unittest.TestCase):

    def test_first_line(self):
        with open('README', 'r') as f:
            first_line = f.readline()
            if first_line == 'NuGridpy':
                print('First line OK.')
            else:
                print('Need to fail')
        f.close()

    def test_that_fails(self):
        sys.exit(1)

    def test_that_follows_a_failure(self):
        print('This is after the test_that_fails function.')

if __name__ == '__main__':
    unittest.main()

