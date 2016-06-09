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

if __name__ == '__main__':
    unittest.main()
