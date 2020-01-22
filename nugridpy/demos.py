'''
Demos
=====

This module provides demonstrations of the nugridpy capabilities.
'''

import os

import pkg_resources

from . import mesa_data


def _demo_mesa():
    """Demo reading MESA data and creating some plots."""

    # MESA data from logs
    path = pkg_resources.resource_filename('nugridpy', os.path.join('resources', 'mesa', 'LOGS'))
    m1 = mesa_data(path, data_type='text', history_only=True)
    m1.hrd()
    m1.tcrhoc()
    m1.kippenhahn('star_age', x0=1.15e10, CO_ratio=True)

    # MESA data from HDF5
    path = pkg_resources.resource_filename('nugridpy', os.path.join('resources', 'mesa', 'HDF5'))
    m2 = mesa_data(path, data_type='hdf5')
    m2.hrd()
    # m2.tcrhoc()
    m2.kippenhahn('age')


def _demo_se():
    """Demo reading SE data and creating some plots."""

    path = pkg_resources.resource_filename('nugridpy', os.path.join('resources', 'se'))
    _ = mesa_data(path, data_type='hdf5')


def demo(which=None):
    """Caller for the different demos."""

    DEMOS = {
        'mesa': _demo_mesa,
        'se': _demo_se,
    }

    if which not in DEMOS:
        print("Please indicate which demo you want:")
        for key in DEMOS:
            print("\tdemo('%s')" % key)
    else:
        DEMOS[which]()
