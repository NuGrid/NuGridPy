'''
Demos
=====

This module provides demonstrations of the nugridpy capabilities.
'''

import os

import pkg_resources

from . import mesa_data, nugrid_data


def _demo_mesa():
    """Demo reading MESA data and creating some plots."""

    # MESA data from logs
    path = pkg_resources.resource_filename('nugridpy', os.path.join('resources', 'mesa', 'LOGS'))
    m1 = mesa_data(path, data_type='text')
    m1.hrd()
    m1.tcrhoc()
    m1.plot('logR', 'logRho', 2030)
    m1.kippenhahn('star_age', x0=9.4e8, CO_ratio=True)

    # MESA data from HDF5
    path = pkg_resources.resource_filename('nugridpy', os.path.join('resources', 'mesa', 'HDF5'))
    m2 = mesa_data(path, data_type='hdf5')
    m2.hrd()
    m2.tcrhoc()
    m2.plot('radius', 'rho', 2030,
            ylabel=r'$\rho/\,[g.cm^{-3}]$', logx=True, logy=True, legend=False)
    m2.kippenhahn('age')


def _demo_nugrid():
    """Demo reading SE data and creating some plots."""

    # out
    path = pkg_resources.resource_filename(
        'nugridpy', os.path.join('resources', 'nugrid', 'H5_out'))
    n_out = nugrid_data(path)
    n_out.plot('age', 'deltat')
    n_out.plot('radius', 'dcoeff', 580)
    n_out.plot('radius', ('H-1', 'He-4',), 300, logy=True)


def demo(which=None):
    """Caller for the different demos."""

    DEMOS = {
        'mesa': _demo_mesa,
        'nugrid': _demo_nugrid,
    }

    if which not in DEMOS:
        print("Please indicate which demo you want:")
        for key in DEMOS:
            print("\tdemo('%s')" % key)
    else:
        DEMOS[which]()
