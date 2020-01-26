'''
Demos
=====

This module provides demonstrations of the nugridpy capabilities.
'''

import os

import pkg_resources

from . import mesa_data, nugrid_data


def _demo_mesa_logs():
    """Demo reading MESA data from logs and creating some plots."""

    path = pkg_resources.resource_filename('nugridpy', os.path.join('resources', 'mesa', 'LOGS'))
    m1 = mesa_data(path, data_type='text')

    # Extracting/plotting data
    print("Header attributes are: \n\t", ' '.join(m1.hattrs), '\n')
    print("Cycles attributes are: \n\t", ' '.join(m1.cattrs), '\n')
    print("Columns available are: \n\t", ' '.join(m1.dcols), '\n')

    m1.hrd()
    m1.tcrhoc()
    m1.plot('logR', 'logRho', [2030, 2210])
    m1.kippenhahn('star_age', x0=9.4e8, CO_ratio=True)


def _demo_mesa_hdf5():
    """Demo reading MESA data from HDF5 data and creating some plots."""

    path = pkg_resources.resource_filename('nugridpy', os.path.join('resources', 'mesa', 'HDF5'))
    m2 = mesa_data(path, data_type='hdf5')

    # Extracting/plotting data
    print("Header attributes are: \n\t", ' '.join(m2.hattrs), '\n')
    print("Cycles attributes are: \n\t", ' '.join(m2.cattrs), '\n')
    print("Columns available are: \n\t", ' '.join(m2.dcols), '\n')
    print("Isotopes in the data are: \n\t", ' '.join(m2.isotopes), '\n')

    m2.hrd()
    m2.tcrhoc()
    m2.plot('radius', 'rho', 2030, ylabel=r'$\rho/\,[g.cm^{-3}]$',
            logx=True, logy=True, legend=False)
    m2.kippenhahn('age')


def _demo_nugrid_out():
    """Demo reading SE data and creating some plots."""

    # OUT
    path = pkg_resources.resource_filename(
        'nugridpy', os.path.join('resources', 'nugrid', 'H5_out'))
    n_out = nugrid_data(path)

    # Extracting/plotting data
    print("Header attributes are: \n\t", ' '.join(n_out.hattrs), '\n')
    print("Cycles attributes are: \n\t", ' '.join(n_out.cattrs), '\n')
    print("Columns available are: \n\t", ' '.join(n_out.dcols), '\n')
    print("Isotopes in the data are: \n\t", ' '.join(n_out.isotopes), '\n')

    n_out.plot('age', 'deltat')
    n_out.plot('radius', 'dcoeff', 580)
    n_out.plot('radius', ('H-1', 'He-4',), 300, logy=True)
    n_out.abu_profile(500, isos=['H-1', 'He-4', 'C-12', 'C-13', 'N-14', 'O-16'])
    n_out.iso_abund(500, a_max=100)


def demo(which=None):
    """Caller for the different demos."""

    DEMOS = {
        'mesa_logs': _demo_mesa_logs,
        'mesa_hdf5': _demo_mesa_hdf5,
        'nugrid_out': _demo_nugrid_out,
    }

    if which not in DEMOS:
        print("Please indicate which demo you want:")
        for key in DEMOS:
            print("\tdemo('%s')" % key)
    else:
        DEMOS[which]()
