NuGridPy
========

NuGridPy is a Python package containing tools to access and analyse 
(e.g. plot) various output formats (including the hdf5-based se format) 
from NuGrid codes (mppnp and ppn) and from the MESA stellar evolution
code. In principle the NuGridPy package can be used on any stellar
evolution code output if the tools to write se hdf5 files available on
the NuGrid web page are used. The mesa.py module will work with MESA
ASCII output in the 'LOGS' directory.

MOST RECENT SVN REVISION STILL CONTAINING NUGRIDPY: r6619

NuGridPy has been migrated to this github repository from the svn
<<<<<<< HEAD:README.md
repository. All components present here have now been removed from
the svn. Please checkout r6619 if you would like access to the
most recent version in the svn prior to NuGridPy's removal.

=======
repository. All components present here have now been removed from the
svn. Please checkout r6619 if you would like access to the most recent
version in the svn prior to NuGridPy's removal.

NuGridPy is a Python package containing tools to access and analyse
(e.g. plot) various output formats (including the hdf5-based se
format) from NuGrid codes (mppnp and ppn) and from the MESA stellar
evolution code. In principle the NuGridPy package can be used on any
stellar evolution code output if the tools to write se hdf5 files
available on the NuGrid web page are used.  The mesa.py module will
work with MESA ASCII output in the 'LOGS' directory. These modules
were written with an interactive work mode in mind, in particular
taking advantage of the interactive ipython session that we usually
start with 'ipython --pylab' or inside an ipython notebook.
>>>>>>> NuGrid/master:README

Installation and Use
====================

These modules were written with an interactive work mode in mind, in
particular taking advantage of the interactive ipython session that we usually
start with `ipython --pylab` and work inside an ipython or jupyter notebook.

You may also "install" NuGridPy by copying the directory to the
site-packages directory wherever your primary python installation lives.

Either way, you should be able to import the mesa (or nugridse or ppn)
module (depending on which type of data you are working with):
`from NuGridPy import mesa as ms`

and read the docstring:
`help(ms)`

There are reasonable doc strings in the modules. If you have made tested
and debugged improvements we are happy to know about them and we may
add them to the release available on the web page. The tools provided
here are useful to us, but of course there are still many things that
need attention and improvement.  We have a good list of scheduled
improvements, let us know if you want to help with these.
Pull requests and new issues are most welcome!


Required packages
=================

* standard:
  - xlrd
  - numpy
  - matplotlib
  - h5py
  - unittest

* python2/3 compatibility
  - future (provides 'past' and 'builtin' modules)

Most of these are standard. If necessary, you can install them with pip, e.g.
`pip install future`



Special situations
------------------
To install several required packages at once
  > sudo pip(3) install -U xlrd future h5py

or, if you want to install in the user home directory,
  > pip(3) install --user -U xlrd future h5py

- on Fedora python2:
  > dnf install -y h5py python2-future

If h5py complains about missing `hdf5.h` or `hdf5_hl.h`:
  > sudo apt-get install libhdf5-dev


<<<<<<< HEAD:README.md
Required package installation if you use `macports`:
=======
Mac installation with `macports`
--------------------------------
>>>>>>> NuGrid/master:README

1. Install ipython and pip:
* port install py34-ipython py34-pip
* port select --set ipython py34-ipython
* port select --set pip py34-pip

2. Install standard prerequisites with pip:
<<<<<<< HEAD:README.md
* pip install -U --user numpy scipy xlrd matplotlib future 
=======
  > pip install -U --user numpy scipy xlrd matplotlib future unittest
>>>>>>> NuGrid/master:README

3. Install h5py:
* port install openmpi-gcc5
* port install hdf5 +hl +openmpi +cxx +fortran +gcc5
* HDF5_DIR=~/macports CFLAGS="-I$HOME/macports/include/openmpi-gcc5" pip install -U --user h5py

Errors about failure to rebuild netcdf package can be ignored.

4. Set up PYTHONPATH environment variable to point to the
   parent directory of NuGridPy, e.g.:
  > export PYTHONPATH="$HOME/NuGrid"


TESTING
=======

The following command runs a suite of test cases for the NuGridPy package:

   > python -m NuGridPy.selftest
