Handling Data with NuGridPy
===========================

NuGridPy's data handling functionality is documented here. This pertains
to reading in files of various data types, writing to file, and
accessing data both within NuGridPy modules as a developer and in an
interactive python environment as a user.

Currently, the following filetypes are supported for read-in:

 - MESA data (e.g., history.data, profiles.index, profile data, etc.)
 - Trajectories

The modules that handle data are :code:`io.py` for input/output, which
is semi-private and generally hidden from users, and :code:`data.py`,
which handles interactive data access and :code:`get` functionality with
data object instances.

I/O
---

The :code:`io.py` module contains a class structure for reading in
ASCII files as generally as possible with :code:`io.TextFile`. For
each filetype, an associated class variable contains specific
information (such as header structure), and the :code:`readfile`
method serves as the core read-in tool.

This module is not user-facing. Instead, the :code:`data.py` module
uses the I/O tools here to present the user with workable data
objects. 

To add a new data type to the :code:`io.py` TextFile class, the
following procedure can be followed:

First, define a new :code:`TextFile` class variable:

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

This function is for parsing header from :code:`history.data` files.

The :code:`TextFile.readfile` method in this module is used by the
:code:`data.py` module for reading in files passed by the user.

Data
----

The :code:`data.py` module is the user-facing component of NuGridPy's
I/O functionality. Here, data objects for each filetype that
can read in using nugridpy are created. Each filetype has a class
associated with it, which can be instantiated at the interactive
level with a simple directory or filename pass. The classes each
inherit from a base class, :code:`FileIn` which houses the 
:code:`init` and the :code:`get` and :code:`get_header` methods.

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

References
----------

.. automodule:: nugridpy.data
   :members:

.. automodule:: nugridpy.io
   :members:
