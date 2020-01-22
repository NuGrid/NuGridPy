Handling Data with NuGridPy
===========================

NuGridPy's data handling functionality is documented here. This pertains
to reading in files of various data types, writing to file, and
accessing data both within NuGridPy modules as a developer and in an
interactive python environment as a user.

Currently, the following filetypes are supported for read-in:

 - MESA data (e.g., history.data, profiles.index, profile data, etc.)
 - Trajectories
 - HDF5 files from NuGrid or MESA

The modules that handle data are :code:`io.py` for input/output, which
is semi-private and generally hidden from users, and :code:`data.py`,
which handles interactive data access and :code:`get` functionality with
data object instances.

I/O
---

The :code:`io.py` module contains a class structure for reading in
HDF5 and ASCII files as generally as possible with :code:`io.TextFile` and
:code:`io.Hdf5File`. For each filetype, an associated class variable 
contains specific information (such as header structure), and the 
respective :code:`readfile` method serves as the core read-in tool.

This module is not user-facing. Instead, the :code:`data.py` module
uses the I/O tools here to present the user with workable data
objects. 

**ASCII data**

To add a new ASCII data type to the :code:`io.py` TextFile class, users
can carry out the following procedure:

First, define a new :code:`TextFile` class variable.

.. code:: python

   NEW_DATA = {
      'name':
      'parser':
      'skip_header':
      'names':
      }

Each key in this dictionary should have a value defined according
to the new data type being read in. Each has a specific purpose,
defined as follows:

 - :code:`name` : The same as the class variable name.
 - :code:`parser` : A specific parsing function (to be defined
   manually) that allows the header data specific to this data
   type to be read in. If there is no header data, this can simply
   be set to :code:`None`.
 - :code:`skip_header` : An :code:`int` referring to the number of
   header lines in the file before actual data begins.
 - :code:`names` : Should be set to :code:`True` if the header/data
   column names are found directly after :code:`skip_header` in
   a whitespace-delimited row. Otherwise, the names of each column
   should be input manually here as a list of strings.

After this class variable is defined, a header parsing function will
likely also need to be defined (if the header structure is unique).
These functions are defined at the class level. See, e.g,

.. code:: python

   
   def parse_history(filename):
       """
       The history.data headers are parsed and read in, starting from
       row 2 of the file. An exception is raised if the number of header
       name columns and header data columns are unequal. Header names
       are pulled directly from their row in history.data by numpy.

       Arguments:
           filename (str): Path to file being read in.

       Returns:
           data (numpy.ndarray): Structured array containing data column
               names and data from data fields in dict-like format.
       """
       # Read the header data using specific file structure
       header_data = np.genfromtxt(filename, dtype=None, skip_header=1,
           names=True, invalid_raise=True, max_rows=1, encoding=None)

       # Return header data as a numpy structured array
       return header_data

This function is for parsing header data from :code:`history.data` files. 
If no header parsing is required, the :code:`parser` key should be set to
:code:`None`.

The :code:`TextFile.readfile` method in this module is used by the
:code:`data.py` module for reading in files passed by the user.

**HDF5 data**

If you're working data from, for instance, MESA or a NuGrid set, you may
want to work with HDF5 files. This is also provided in parallel fashion
through the :code:`io.py` module. Again, class variables are defined, this
time in the :code:`Hdf5File` parent class.

Defining a new class variable is done in the same way, with different key
value pairs:

.. code:: python

    NEW_HDF5_TYPE = {
        'default_dataset_name': ,
        'group_default_dataset_name': ,
        'group_prefix':
    }

For HDF5 files, the keys are defined as follows:

 - :code:`default_dataset_name`: The default name of stored HDF5 datasets.
 - :code:`group_default_dataset_name`: The default name for datasets stored
   in groups.
 - :code:`group_prefix`: Assuming groups are cycles or indices of some kind,
   this is any prefix that comes before the cycle or index number. 

The user-facing aspects of data handling, in both the ASCII and HDF5 cases,
are carried out by :code:`data.py`, which interacts with the parent classes
for each data type. 

Data
----

The :code:`data.py` module is the user-facing component of NuGridPy's
I/O functionality. Here, data objects for each filetype that
can read in using nugridpy are created. Each filetype has a class
associated with it, which can be instantiated at the interactive
level with a simple directory or filename pass. For ASCII data, the 
classes each inherit from a base class, :code:`FileIn`, which houses the 
:code:`init` and the :code:`get` and :code:`get_header` methods. For
HDF5 data, each new type or source of data has (or should have) a data
object class associated with it.

**ASCII data**

If a new data type has been defined in :code:`io.py`, allowing the
user to create a data object from it is as simple as creating a new
class as follows:

.. code:: python

    class NewData(FileIn):
        """Read in some data."""

    def __init__(self, filename):
        super().__init__(filename, io.TextFile.NEW_DATA)

and that's it! As long as the new data type is defined as
:code:`NEW_DATA` in :code:`io.py`, you'll now have a new
data object complete with :code:`get` functionality.

If a directory structure is required for your read-in, you'll
need to do a bit more and define a mediating class that defines
this structure, as is the case with :code:`class MesaData`.

To read in a particular file as a data object, e.g., a
trajectory, simply execute the following:

.. code:: python

    my_data = data.Trajectory('path/to/trajectory.input')

Or, in the case of a MESA data instance, where more than one
file is to be read in, you can pass in a directory structure:

.. code:: python

    my_mesa = data.MesaData('path/to/MESA_DIR/')

Note that in the case of MESA data, files (such as history data)
in your data directory are assumed to have a default filename. If
your files have a custom name, you'll need to input it as an
argument, e.g., :code:`MesaData('/DIR/', history_name='my_file')`.

To see a list of available filetypes, you can print the
:code:`data.FileIn.filetypes` attribute.

**HDF5 data**

To read in HDF5 data interactively, you will need to use the data-
specific read-in class. For Nugrid data, for instance, you would use
:code:`data.NugridData`. As in the ASCII MESA data case, you pass in
the file directory corresponding to the NuGrid data set of interest.

.. code:: python

    data_in = data.NugridData('/path/to/NugridDir/')

This object will now contain all the relevant attributes and
:code:`get` methods for working with Nugrid HDF5 data. The three
types of get are as follows:

 - :code:`get_metadata('<metadata_name>')`
 - :code:`get_cycle_header(<cycle_num>, '<cycle_header_name>')`
 - :code:`get_from_cycle(<cycle_num>, '<cycle_data_name>')`

Attributes of a NugridData object include :code:`ages`, :code:`cycles`,
:code:`metadata_names`, and others when relevant, such as :code:`Z`.

References
----------

.. automodule:: nugridpy.data
   :members:

.. automodule:: nugridpy.io
   :members:
