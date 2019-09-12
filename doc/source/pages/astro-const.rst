Doing Astrophysics With NuGridPy
================================

On this page, the astronomy and constants modules are documented. Together,
these modules comprise the basic astronomical and physical calculations
and constants necessary for doing astrophysics with NuGridPy.


Astronomy Module
----------------

The NuGridPy astronomy module contains functions for carrying
out physical and astrophysical computations. It can be imported
from nugridpy with a standard import:

.. code:: python

   from nugridpy import astronomy

Each function in the module takes in a set of parameters and returns a data 
structure, usually simply a float value. The functions may contain constants 
imported from the constants module, which can be read out in a tuple as 
follows:

.. code:: python

   astronomy.Pgas.constants

will return:

.. code:: python

   (Constant(1.3806504e-16, boltzmann_const: Boltzmann constant (k), erg K^-1),
    Constant(1.660538782e-24, atomic_mass_unit: Atomic mass unit (u), g))

Here, we see that the function for calculating gas pressure uses both the
:code:`boltzmann_const` and :code:`atomic_mass_unit` functions. The
elements of the tuple can be accessed and interacted with as their actual
float values.

Examples
~~~~~~~~

Using the astronomy functions is like using any other python function. Necessary
variables are passed in, and a result is returned. The functions generally
perform astrophysical computation, such as in a determination of
radiative viscosity:

.. code:: python

    from nugridpy import astronomy

    l = 100*1.e5 # 100km
    v = 1.e5     # typical velocity
    T = 90.e6    # temperature
    X = 0.001    # H mass fraction
    rho = 100.   # density

Here, we have imported the astronomy module and defined the relevant variables.
Next, we use a function and perform a calculation:

.. code:: python

    nu = astronomy.visc_rad_kap_sc(T, rho, X)
    Re_rad = v*l / nu
    print(Re_rad)
    4.43512e+08

Here, the radiative viscosity was calculated using the :code:`visc_rad_kap_sc`
function, then used in further calculation to yield a result.


Constants Module
----------------

In the constants module, the Constant class and its instances are defined.
A variety of physical and astronomical constants are instantiated (e.g.,
Planck's constant) with a value, description, and the constant's cgs
units. The module can be imported with a standard import:

.. code:: python

   from nugridpy import constants

To create a new constant, the Constant class is instantiated in the
following way:

.. code:: python

   new_constant = Constant(1., 'The speed of light in natural units', 'Dimensionless')

The above code can be executed interactively, and a temporary constant will
be created for use in your session. However, if this is done in the module
itself, then this constant would become permanent. It could then be imported
for use elsewhere later, and it would be included with the general import of
the module demonstrated earlier.

Do not add permanent constants to :code:`constants.py` unless they would be used 
generally, either by users or by another module in the NuGridPy package (e.g., 
the astronomy module).

References
----------

Note
~~~~
When adding content to this documentation, note that functions which are
decorated to include constants with the attach_constants decorator must
be added manually. Sphinx does not recognize the docstrings of decorated
functions.

.. automodule:: nugridpy.astronomy
   :members:

.. autofunction:: nugridpy.astronomy.visc_mol_sol

.. autofunction:: nugridpy.astronomy.visc_rad_kap_sc

.. autofunction:: nugridpy.astronomy.Gamma1_gasrad

.. autofunction:: nugridpy.astronomy.Pgas

.. autofunction:: nugridpy.astronomy.Prad

.. autofunction:: nugridpy.astronomy.mimf_ferrario

.. autofunction:: nugridpy.astronomy.core_mass_L

.. autofunction:: nugridpy.astronomy.imf

.. autofunction:: nugridpy.astronomy.int_imf_dm

.. autofunction:: nugridpy.astronomy.am_orb

.. autofunction:: nugridpy.astronomy.mass_loss_loon05

.. autofunction:: nugridpy.astronomy.energ_orb

.. autofunction:: nugridpy.astronomy.period

.. autofunction:: nugridpy.astronomy.escape_velocity

.. autofunction:: nugridpy.astronomy.Nasv

.. autofunction:: nugridpy.astronomy.macs

.. autofunction:: nugridpy.astronomy.mu_e

.. autofunction:: nugridpy.astronomy.mu

.. autofunction:: nugridpy.astronomy.Trho_idrad

.. autofunction:: nugridpy.astronomy.Trho_iddeg


.. automodule:: nugridpy.constants
   :members:
