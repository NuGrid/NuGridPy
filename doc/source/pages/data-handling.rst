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
function serves as the core read-in tool.

References
----------

.. automodule:: nugridpy.data
   :members:

.. automodule:: nugridpy.io
   :members:
