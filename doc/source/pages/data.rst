.. _data:

Handling Data
=============

.. contents::
    :local:
    :depth: 1

NuGridPy's data handling functionality is documented here. This pertains
to reading in files of various data types, writing to file, and
accessing data both within NuGridPy modules as a developer and in an
interactive python environment as a user.

Currently, the following data types are supported for read-in:

* MESA data, both from logs/profile files and HDF5 files
* NuGrid data

The modules that handle data are:

* :py:mod:`~nugridpy.io` which is used for reading/writing data but not exposed to the users
* :py:mod:`~nugridpy.data` which provides an API for access/manipulating the data

I/O
---

.. warning:: This section is for developers only.

The :py:mod:`~nugridpy.io` module contains two class structures, :py:class:`~nugridpy.io.TextFile`
and :py:class:`~nugridpy.io.Hdf5File` for reading in :code:`ASCII` and :code:`HDF5` files,
respectively, as generally as possible.

For each filetype, an associated class variable
contains specific information (such as header structure), and the
respective :code:`readfile` method serves as the core read-in tool.

**ASCII data**

The implementation of the :py:class:`~nugridpy.io.TextFile` class is generic
and used for loading data from a text file. However, one must provide
the reading method with some information about the data, such that it can extract it
and structure it properly. In particular, the reading methods should know:

* how big the header is
* how to extract the header data and names
* is the columns names are found in the file or set a-priori

This information is contained in attributes of the class ,
such as :code:`MESA_DATA`:

.. code:: python

    MESA_DATA = {
        'name': 'MESA_DATA',
        'parser': parse_history,
        'skip_header': 5,
        'names': None
    }

Each key in this dictionary should have a value defined according
to the new data type being read in. Each has a specific purpose,
defined as follows:

* :code:`name`: the same as the class variable name.
* | :code:`parser`: a specific parsing function (to be defined manually) that reads the header
  | specific to this data to be read in. If no header data, this should be set to :code:`None`.
* :code:`skip_header` : the number of header lines in the file before actual data begins.
* | :code:`names` : should be set to :code:`True` if the header/data column names are found directly
  | after :code:`skip_header` in a whitespace-delimited row. Otherwise, the names of each
  | column should be input manually here as a list of strings.

After this class variable is defined, a header parsing function will
likely also need to be defined (if the header structure is unique).
These functions are defined outside the class. See, e.g,

.. code:: python

   def parse_history(filename):
       # Read the header data using specific file structure
       header_data = np.genfromtxt(filename, dtype=None, skip_header=1,
           names=True, invalid_raise=True, max_rows=1, encoding=None)

       # Return header data as a numpy structured array
       return header_data

To add a new ASCII data type, define your class attribute
:code:`MY_NEW_DATA` similar to :code:`MESA_DATA`, and then simply do:

.. code:: python

    from nugridpy.io import TextFile

    my_header, my_data = TextFile.readfile(filename, TextFile.MY_NEW_DATA)


**HDF5 data**

If you're working data from, for instance, MESA or a NuGrid set, you may
want to work with HDF5 files. This is also provided in similar fashion
by the :py:class:`~nugridpy.io.Hdf5File` class.

Similarly to the ASCII case, the HDF5 read-in function is generic.
Specific information about the data is also provided by class attributes such as:

.. code:: python

    NUGRID = {
        'default_dataset_name': 'data',
        'group_default_dataset_name': 'SE_DATASET',
        'group_prefix': 'cycle'
    }

For HDF5 files, the keys are defined as follows:

* | :code:`default_dataset_name`: the default name of stored HDF5 datasets.
* | :code:`group_default_dataset_name`: the default name for datasets stored in groups.
* | :code:`group_prefix`: assuming groups are cycles or indices of some kind, this is any
  | prefix that comes before the cycle or index number.

Data
----

The :code:`data.py` module is the user-facing component of NuGridPy's
I/O functionality. Here, data objects are created using the following classes:

* :py:class:`~nugridpy.data.MesaData` for MESA data
* :py:class:`~nugridpy.data.NugridData` for NuGrid data

Note that the :py:class:`~nugridpy.data.MesaData` acts like a proxy
between the classes :py:class:`~nugridpy.data.MesaDataText` and
:py:class:`~nugridpy.data.MesaData`.

These classes can be used to instantiate data objects.
However, we recommend to use the simpler functions
:py:func:`~nugridpy.__init__.mesa_data` and :py:func:`~nugridpy.__init__.nugrid_data`.
For instance:

.. code:: python

    from nugridpy import mesa_data, nugrid_data

    m1 = mesa_data('/path/to/LOGS', data_type='text')
    m2 = mesa_data('/path/to/HDF5', data_type='hdf5')
    nu1 = nugrid_data('/path/to/H5_out')
    nu2 = nugrid_data('/path/to/H5_surf')


Data can be extracted easily:

.. code:: python

    # Attributes and colummn names
    header_attributes = m1.hattrs
    cycle_attributes = m1.cattrs
    isotopes = m2.isotopes
    data_columns = nu1.dcols

    # Header attributes
    z = m1.get_hattr('initial_z')

    # Cycle attributes
    age = m2.get_cattr('age')

    # Data columns for one or two cycles
    r = nu1.get_dcol('radius', 100)
    rho = nu1.get_dcol('rho', [500, 800])

    # Decay data
    decay = nu2.get_dcol('elem_massf_decay', 500)

References
----------

.. automodule:: nugridpy.__init__
   :members:
   :member-order: bysource
   :exclude-members: __weakref__

.. automodule:: nugridpy.data
   :members:
   :member-order: bysource
   :exclude-members: __weakref__

.. automodule:: nugridpy.io
   :members:
   :member-order: bysource
   :exclude-members: __weakref__, parse_history, parse_trajectory
