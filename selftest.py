import unittest

class TestModuleImports(unittest.TestCase):

    def test_import_ascii_table(self):
        import NuGridPy.ascii_table

    def test_import_astronomy(self):
        import NuGridPy.astronomy

    def test_import_data_plot(self):
        import NuGridPy.data_plot

    def test_import_grain(self):
        import NuGridPy.grain

    def test_import_h5T(self):
        import NuGridPy.h5T

    def test_import_mesa(self):
        import NuGridPy.mesa

    def test_import_nugridse(self):
        import NuGridPy.nugridse

    def test_import_ppn(self):
        import NuGridPy.ppn

    def test_import_utils(self):
        import NuGridPy.utils

class TestAbuChart(unittest.TestCase):

    def test_abu_chart(self):
        from NuGridPy import utils,ppn,data_plot
        import matplotlib.pylab as mpy
        import os

        # test_data_dir should point to the correct location of a set of abundances data file
        #test_data_dir= os.path.dirname(data_plot.__file__) + "/test_data"
        test_data_dir= '.'
        p=ppn.abu_vector(test_data_dir)
        mp=p.get('mod')
        if len(mp) == 0:
            raise IOError("Cannot locate a set of abundance data files")
        sparse=10
        cycles=mp[:1000:sparse]
        form_str='%6.1F'
        form_str1='%4.3F'

        i=0
        for cyc in cycles:
            T9  = p.get('t9',fname=cyc)
            Rho = p.get('rho',fname=cyc)
            mod = p.get('mod',fname=cyc)
            # time= p.get('agej',fname=cyc)*utils.constants.one_year
            time= p.get('agej',fname=cyc)
            mpy.close(i);mpy.figure(i);i += 1
            p.abu_chart(cyc,mass_range=[0,41],plotaxis=[-1,22,-1,22],lbound=(-6,0),show=False)
            mpy.title(str(mod)+' t='+form_str%time+'yr $T_9$='+form_str1%T9+' $\\rho$='+str(Rho))
            png_file='abu_chart_'+str(cyc).zfill(len(str(max(mp))))+'.png'
            mpy.savefig(png_file)
            self.assertTrue(os.path.exists(png_file))



if __name__ == '__main__':
    unittest.main()
