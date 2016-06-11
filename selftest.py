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
        nugrid_dir= os.path.dirname(os.path.dirname(ppn.__file__))
        NuPPN_dir= nugrid_dir + "/NuPPN"
        test_data_dir= NuPPN_dir + "/examples/ppn_C13_pocket/master_results"
        p=ppn.abu_vector(test_data_dir) # TODO: this function fails to raise an exception if path is not found!
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
            os.remove(png_file)


    def test_abu_evolution(self):
        from NuGridPy import ppn, utils
        import matplotlib.pylab as mpy
        import os

        nugrid_dir= os.path.dirname(os.path.dirname(ppn.__file__))
        NuPPN_dir= nugrid_dir + "/NuPPN"
        test_data_dir= NuPPN_dir + "/examples/ppn_Hburn_simple/RUN_MASTER"

        symbs=utils.symbol_list('lines2')
        x=ppn.xtime(test_data_dir)
        specs=['PROT','HE  4','C  12','N  14','O  16']
        i=0
        for spec in specs:
            x.plot('time',spec,logy=True,logx=True,shape=utils.linestyle(i)[0],show=False,title='')
            i += 1
        mpy.ylim(-5,0.2)
        mpy.legend(loc=0)
        mpy.xlabel('$\log t / \mathrm{min}$')
        mpy.ylabel('$\log X \mathrm{[mass fraction]}$')
        abu_evol_file = 'abu_evolution.png'
        mpy.savefig(abu_evol_file)
        self.assertTrue(os.path.exists(abu_evol_file))

if __name__ == '__main__':
    unittest.main()
