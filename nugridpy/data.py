"""
Data handling
"""

from contextlib import suppress
import logging
import os
import re

from . import io


logging.basicConfig(level=logging.INFO)


class FileIn:
    """
    Receives a text file and creates a data object from its contents.

    Attributes:
        filetypes (tuple): Provides a tuple containing the names of filetypes
            that can be read in as a FileIn object.
        header_names (tuple): A tuple containing the string values of each
            header name (e.g., 'model_number').
        data_names (tuple): A tuple containing the string values of each
            data column name (e.g., 'mass').

    Arguments:
        filename (str): File being read in.
        filetype (dict): Filetype as read from io.TextFile.FILE_TYPES.

    Methods:
        get(self, data_name): Receives a data_name string and returns the
            corresponding data value.
        get_header(self, header_name): Receives a header_name string and
            returns the corresponding value.
    """

    # Class variable to show which files are readable, generated from io filetypes
    filetypes = tuple(ftype['name'] for ftype in io.TextFile.FILE_TYPES)

    def __init__(self, filename, filetype):
        """Constructs a data instance from a given file."""
        # These variables are not user-meaningful
        self._header_object, self._data_object = io.TextFile.readfile(
            filename, filetype)

        # These variables allow the user to see header and data names
        self.header_names = self._header_object.dtype.names
        self.data_names = self._data_object.dtype.names

    def get(self, data_name):
        """Basic get method. Returns data attribute requested."""
        # Prompt a data_names check if it doesn't contain request
        if data_name not in self.data_names:
            raise ValueError('No data exists with that header. Check data_names '
                             'for availability.')

        # Declare the value and return it
        value = self._data_object[data_name]
        return value

    def get_header(self, header_name):
        """Get method for header data."""
        # Prompt a header_names check if requested name isn't contained in the namespace
        if header_name not in self.header_names:
            raise ValueError('No such header data exists. Check header_names for '
                             'availability.')

        # The .item() syntax must be used for singletons in a numpy structured array
        if self._header_object[header_name].size == 1:
            value = self._header_object[header_name].item(0)
        else:
            value = self._header_object[header_name]

        return value


class MesaData:
    """
    For reading MESA data from a provided directory.

    Attributes:
        history_data (data.FileIn): A FileIn instance that contains the data from a
            history.data file.
        profiles_index (data.FileIn): A FileIn instance that contains the data from
            a profiles.index file.
        profiles_in_dir (tuple): A tuple containing string filenames for each MESA
            profile found in filedir.
        profiles (list): A list populated with one FileIn object for each profile
            that is read in.

    Arguments:
        filedir (str): Directory path to MESA data.
        history_name (:obj:`str`, optional): Filename in the given directory for the
            history data file. The default is history.data, or can be set by user.
        index_name (:obj:`str`, optional): Filename in the given directory for the
            mesa profiles index file. Default is profiles.index.
        profile_prefix (:obj:`str`, optional): The filename prefix for profile data,
            such as profile or log (e.g., profile32.data). The default is profile.
        profile_suffix (:obj:`str`, optional): The filename suffix for profile data.
            The default is .data.
        all_profiles (:obj:`bool`, optional): Determines whether all profiles are read
            in or not. Default is True.
        history_only (:obj:`bool`, optional): Set to True if history.data is the only
            thing you want to read in. False by default.
        profiles_only (:obj:`bool`, optional): Set to True if profile data is the only
            thing you want to read in. False by default.
    """

    def __init__(self, filedir, history_name='history.data', index_name='profiles.index',
                 profile_prefix='profile', profile_suffix='.data', all_profiles=True,
                 history_only=False, profiles_only=False):
        """Constructs a MESA data instance."""
        # Raise error if both history_only and profiles_only are True simultaneously
        if history_only and profiles_only:
            raise ValueError('At least one of history_only and profiles_only must be '
                             'False.')

        # history.data read-in
        if not profiles_only:
            # Try to read and alert user if file not present
            try:
                self.history_data = MesaFile(os.path.join(filedir, history_name))
            except IOError:
                raise IOError('No file named {} could be found. Check filedir path or '
                              'set history_name.'.format(history_name))

        # profiles.index read-in
        if not history_only:
            # Try to read and alert user if file not present
            try:
                self.profiles_index = ProfilesIndex(os.path.join(filedir, index_name))
            except IOError:
                raise IOError('No file named {} could be found. Check filedir or set '
                              'index_name.'.format(index_name))

            # Collect the profiles from filedir
            files = os.listdir(filedir)
            self.profiles_in_dir = tuple(prof for prof in files if re.match(
                '{}[0-9]\d*{}'.format(profile_prefix, profile_suffix), prof))
            prof_count = len(self.profiles_in_dir)

            # Read in all profiles to a user-facing list attribute
            if all_profiles:
                print('Reading in {} profiles... set all_profiles to False to '
                      'skip.'.format(prof_count))
                self.profiles = [MesaFile(os.path.join(filedir, prof))
                                 for prof in self.profiles_in_dir]


class MesaFile(FileIn):
    """Read in MESA data file in history/profile format."""

    def __init__(self, filename):
        super().__init__(filename, io.TextFile.MESA_DATA)


class ProfilesIndex(FileIn):
    """Read in a MESA profiles index file."""

    def __init__(self, filename):
        super().__init__(filename, io.TextFile.PROFILES_INDEX)


class Trajectory(FileIn):
    """Read in a trajectory file and make a FileIn instance."""

    def __init__(self, filename):
        super().__init__(filename, io.TextFile.TRAJECTORY)


class NugridData:
    """
    For reading in NuGrid HDF5 data.
    """

    def __init__(self, filedir, file_ext='.h5'):
        # List files in filedir and collect all hdf5 files with regex
        self.filedir = filedir
        self.h5_in_dir = tuple(h5 for h5 in sorted(os.listdir(filedir))
                               if re.match('.*\{}$'.format(file_ext), h5))

        # List over all HDF5 files
        for ind, h5 in enumerate(self.h5_in_dir):
            # Filename
            filename = os.path.join(filedir, h5)

            # First file defines the metadata and datasets attributes
            if not ind:
                self._metadata, self._data, self._headers, datasets = io.Hdf5File.readfile(
                    filename, io.Hdf5File.NUGRID)

                # For the moment, we hard-code the creation of A, Z, and isomeric_state attributes
                # Maybe later, we can create an attribute for each dataset found instead?
                with suppress(KeyError):
                    self.A = datasets['A']
                    self.Z = datasets['Z']
                    self.isomeric_state = datasets['isomeric_state']
                    logging.info('Attributes for A, Z, and isomeric_state have been created.')

            else:  # only care about groups
                _, groups_data, header_data, _ = io.Hdf5File.readfile(
                    filename, io.Hdf5File.NUGRID)

                # We only need to stack the groups data
                self._data.update(groups_data)
                self._headers.update(header_data)

        # Collect the relevant attributes of the hdf5 files into self vars
        self.cycles = [key for key in self._data.keys()]
        self.metadata_names = [name for name in self._metadata]
        self.cycle_headers = [attr for attr in self._headers[self.cycles[0]]]
        self.cycle_data = self._data[self.cycles[0]].dtype.names

        # Collect important values from each cycle for self vars
        self.ages = [self._headers[key]['age'].item() for key in self.cycles]

    def get_metadata(self, name):
        """Return metadata value(s)."""
        try:
            return self._metadata[name]
        except KeyError:
            raise KeyError('Cannot find attribute {} in metadata.'.format(name))

    def _get_cycle(self, cycle):
        """Return cycle dataset object."""
        try:
            return self._data[str(cycle)]
        except KeyError:
            raise KeyError('Cannot find cycle {} in the data.'.format(cycle))

    def get_cycle_header(self, cycle, name):
        """Return specific cycle header data."""
        try:
            return self._headers[str(cycle)][name]
        except KeyError:
            raise KeyError('Cannot find header {} in cycle.'.format(name))

    def get_from_cycle(self, cycle, name):
        """Return cycle data."""
        try:
            return self._get_cycle(cycle)[name]
        except ValueError:
            raise KeyError('Cannot find data {} in cycle.'.format(name))
