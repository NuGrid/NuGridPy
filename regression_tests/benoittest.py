#import matplotlib
#matplotlib.use('agg')
import unittest

#from tempdir.tempfile_ import TemporaryDirectory

class Test_README_File(unittest.TestCase):

    def test_first_line(self):
        with open('README', 'r') as f:
        first_line = f.readline()
        if first_line == 'NuGridpy':
            print 'Need to pass'
        else:
            print 'Need to fail'
        f.close()



