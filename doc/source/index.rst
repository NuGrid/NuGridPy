.. nugridpy documentation master file, created by
   sphinx-quickstart on Wed Aug 21 10:27:22 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to NuGridPy's documentation!
====================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   pages/astro-const

NuGridPy is a set of Python tools for the NuGrid community that allows
for interactive analysis, visualization, computation, reading/writing, and other
pre-/post-processing of astrophysical data and simulation results. The source code
for this package is hosted on GitHub_ and can be forked or cloned by the general
public there. NuGrid members should already have access to the repository
directly through their GitHub accounts.

Installation
------------

You can install the latest release of NuGridPy through PyPI with:

.. code::

   pip install nugridpy

Alternatively, you can install the source code via GitHub by cloning the NuGridPy
repository. This can be done by issuing the following command in a git initialized
directory:

.. code::

   git clone git@github.com:NuGrid/NuGridPy.git

if you are authenticating with SSH, or:

.. code::

   git clone https://github.com/NuGrid/NuGridPy.git

with a username and password through HTTPS.

Requirements
------------

All requirements can be found in the requirements.txt file in the source directory. 
NuGridPy users can install dependencies with the following command:

.. code::

   pip install -r requirements.txt

Tests
-----

Testing for NuGridPy is being developed in the :code:`tests/` directory of the source
repository. They can be run with the python :code:`unittest` library manually or
with the help of CI (e.g., Bamboo or Travis). Running either of the following:

.. code::

   python -m unittest -v

.. code::

   nosetests --with-xunit -v tests/

in the source directory will initiate testing. In the second case, an XML report is 
generated. You may require the :code:`nose` package if you haven't installed the 
requirements as outlined above.

Author(s)
---------

The authors are the members of the NuGrid_ collaboration.

License
~~~~~~~

BSD 3-Clause "New" or "Revised" License

Copyright
~~~~~~~~~

\(c\) 2019, NuGrid Collaboration

Version
~~~~~~~

The current stable release of NuGridPy is |version|.

Indices and tables
~~~~~~~~~~~~~~~~~~

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _GitHub: https://github.com/NuGrid/NuGridPy
.. _NuGrid: https://nugrid.github.io/
