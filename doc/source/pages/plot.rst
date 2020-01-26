Plotting
========

.. contents::
    :local:
    :depth: 1


Plotting functionalities are implemented inside several mixins that are defined in the :py:mod:`~nugridpy.plot` module.
These mixins are base classes for the different :ref:`data<data>` classes:

* :py:class:`~nugridpy.plot.MesaPlotMixin` for MESA data
* :py:class:`~nugridpy.plot.NugridPlotMixin` for NuGrid data

We explain in this section how to easily produce plots.

Generic plots
-------------

First, we need to load data. We can use the data availabel in the :code:`resources` folder:


.. code:: python

    from nugridpy import mesa_data, nugrid_data

    path = pkg_resources.resource_filename('nugridpy', os.path.join('resources', 'mesa', 'LOGS'))
    m = mesa_data(path, data_type='text')

Data objects have a generic :py:func:`~nugridpy.plot.PlotMixin.plot` method to visualize headers and data profiles.
The following code

.. code:: python

    # Evolution of some cycle headers
    m.plot('star_age', ['radius', 'log_L'], logy=True)

    # Data profiles for cycles 2030 and 2210
    m.plot('logR', 'logRho', [2030, 2210])

will produce two plots:

|agevsRandL| |logRvslogRho|

.. |agevsRandL| image::  ../images/agevsRandL.png
   :width: 45%
.. |logRvslogRho| image::  ../images/logRvslogRho.png
   :width: 45%


Hertzsprung-Russel diagramm
---------------------------

You can produce Hertzsprung-Russel diagramm using the :py:func:`~nugridpy.plot.PlotMixin.hrd` method:

.. code:: python

    # HRD
    m.hrd()

.. figure::  ../images/HRD.png
   :align:   center
   :scale: 60%


T\ :sub:`c`\-Rho\ :sub:`c`\  diagramm
-------------------------------------

The :py:func:`~nugridpy.plot.PlotMixin.tcrhoc` method plots the central density vs. the central temperature:


.. code:: python

    # TcRhoc diagramm
    m.tcrhoc()

.. figure::  ../images/tcrhoc.png
   :align:   center
   :scale: 60%


Kippenhahn diagramm
-------------------

One can also make a Kippenhahn diagramm with the :py:func:`~nugridpy.plot.PlotMixin.kippenhahn` method:

.. code:: python

    # kippenhahn diagramm
    m.kippenhahn('star_age', x0=9.4e8, CO_ratio=True)

.. figure::  ../images/kippenhahn.png
   :align:   center
   :scale: 60%


Abundances profiles
-------------------

With NuGrid data, one can plot the abundances of some species vs mass coordinates:


.. code:: python

    path2 = pkg_resources.resource_filename('nugridpy', os.path.join('resources', 'nugrid', 'H5_out'))
    n = nugrid_data(path2)

    # Check the available isotopes
    print(n.isotopes)

    # Abundance profiles
    n_out.abu_profile(500, isos=['H-1', 'He-4', 'C-12', 'C-13', 'N-14', 'O-16'])

.. figure::  ../images/abu_profile.png
   :align:   center
   :scale: 60%


Isotopes abundances
-------------------

One can plot the different isotopes abundances using the :py:func:`~nugridpy.plot.PlotMixin.iso_abund` method:

.. code:: python

    path2 = pkg_resources.resource_filename('nugridpy', os.path.join('resources', 'nugrid', 'H5_out'))
    n = nugrid_data(path2)

    # Isotope abundances
    n.iso_abund(500, a_max=100)

.. figure::  ../images/iso_abund.png
   :align:   center
   :scale: 40%


References
----------

.. automodule:: nugridpy.plot
   :members:
   :member-order: bysource
   :exclude-members: __weakref__
