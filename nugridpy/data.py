"""
Data
====

This module provides implementation for handling data.
"""

from contextlib import suppress
import os
from functools import lru_cache

import numpy as np

from .plot import MesaPlotMixin, NugridPlotMixin
from .io import TextFile, Hdf5File
from .isotopes import get_isotope_name, get_A_Z
from .config import logger, CACHE_MAX_SIZE


class DataFromTextMixin:
    """
    Abstract class for data extracted from text file.

    Attributes:
        header_names (tuple): A tuple containing the string values of each
            header name (e.g., 'model_number').
        data_names (tuple): A tuple containing the string values of each
            data column name (e.g., 'mass').

    Arguments:
        filename (str): File being read in.
        filetype (dict): Filetype as read from TextFile.FILE_TYPES.

    Methods:
        get(self, data_name): Receives a data_name string and returns the
            corresponding data value.
        get_header(self, header_name): Receives a header_name string and
            returns the corresponding value.
    """

    def __init__(self, filename, filetype):
        """Constructs a data instance from a given file."""
        # These variables are not user-meaningful
        self._header_object, self._data_object = TextFile.readfile(
            filename, filetype)

    @property
    def header_names(self):
        """List of header names."""
        return self._header_object.dtype.names

    @property
    def data_names(self):
        """List of data names."""
        return self._data_object.dtype.names

    def get_header(self, name):
        """Get method for header object."""
        return self._header_object[name].item()

    def get(self, name):
        """Get method for data object."""
        return self._data_object[name]


class DataFromHDF5Mixin:
    """Class for data obtained from reading text files.

     Arguments:
        filedir (str): Directory path to the data.
        file_ext (str, optional): file extension
    """

    # Name of the member containing the mass fractions
    MASS_FRACTION_HDF5_MEMBER_NAME = 'iso_massf'

    def __init__(self, filedir, file_ext='.h5', **_kwargs):

        # List files in filedir and collect all hdf5 files with regex
        self.filedir = filedir
        self.files_in_dir = tuple(f for f in sorted(os.listdir(filedir))
                                  if f.lower().endswith(file_ext))

        # List over all HDF5 files
        for ind, h5 in enumerate(self.files_in_dir):

            # Filename
            filename = os.path.join(filedir, h5)

            # First file defines the metadata and datasets attributes
            if not ind:
                self._metadata, self._data, self._headers, datasets = Hdf5File.readfile(
                    filename, Hdf5File.NUGRID)

                # For the moment, we hard-code the creation of A, Z, and isomeric_state attributes
                # Maybe later, we can create an attribute for each dataset found instead?
                with suppress(KeyError):
                    self.A = datasets['A']
                    self.Z = datasets['Z']
                    self.isomeric_state = datasets['isomeric_state']

                    # Create isotopes
                    self.isotopes = [
                        get_isotope_name(int(a), int(z)) for a, z in zip(self.A, self.Z)]
                    logger.info(
                        'Attributes for A, Z, isomeric_state, and isotopes have been created.')

            else:  # only care about groups
                _, groups_data, header_data, _ = Hdf5File.readfile(
                    filename, Hdf5File.NUGRID)

                # We only need to stack the groups data
                self._data.update(groups_data)
                self._headers.update(header_data)

        # Collect important values from each cycle for self vars
        self.ages = [self._headers[key]['age'].item() for key in self.cycles]

    @property
    def cycles(self):
        """List of cycles present in the data."""
        return sorted(self._data.keys())

    @property
    def hattrs(self):
        """List of header attributes."""
        return sorted(self._metadata.keys())

    @property
    def cattrs(self):
        """List of cycles headers."""
        return sorted(self._headers[self.cycles[0]].keys())

    @property
    def dcols(self):
        """List of columns available in a cycle."""
        return sorted(self._data[self.cycles[0]].dtype.names)

    def get_hattr(self, name):
        """Return header attribute value(s).

        :param str name: Metadata name.
        """
        try:
            return self._metadata[name].item()
        except KeyError:
            raise KeyError('Cannot find attribute {} in metadata.'.format(name))

    def _get_cattr(self, name, cycle):
        """Return header data from a specific cycle.

        :param str name: Header name
        :param int cycle: cycle number
        """
        try:
            return self._headers[cycle][name].item()
        except KeyError:
            raise KeyError('Cannot find header {} in cycle {}.'.format(name, cycle))

    def get_cattr(self, name, cycles=None):
        """Return header data from specific cycles.

        :param str name: Header name
        :param int or list(int) cycles: cycle numbers. Default: all cycles
        """
        cycles = cycles or self.cycles
        if isinstance(cycles, int):
            return self._get_cattr(name, cycles)
        return np.array([self._get_cattr(name, c) for c in cycles])

    def _get_cycle(self, cycle):
        """Return cycle dataset object.

        :param int cycle: cycle number
        """
        try:
            return self._data[cycle]
        except KeyError:
            raise KeyError('Cannot find cycle {} in the data.'.format(cycle))

    def _get_isotope_mass_fraction(self, name, cycle):
        """Returns the mass fraction profile from a specific cycle."""
        try:
            isotope = get_A_Z(name)

            # Extract data from the single index in A,Z with these values
            index = int(np.where((self.A == isotope.A) * (self.Z == isotope.Z))[0])
            return self._get_cycle(cycle)[self.MASS_FRACTION_HDF5_MEMBER_NAME][..., index]
        except (TypeError, KeyError):
            raise KeyError('Cannot find isotope {} in cycle {}.'.format(name, cycle))

    def _get_dcol(self, name, cycle):
        """Return column data from a specific cycle.

        :param str name: Column name
        :param int or list(int) cycle: cycle number
        """
        try:
            return self._get_cycle(cycle)[name]
        except ValueError:

            # It might be an element that can be accessed through iso_massf
            if name in self.isotopes:
                return self._get_isotope_mass_fraction(name, cycle)
            raise KeyError('Cannot find data {} in cycle {}. '.format(name, cycle))

    def get_dcol(self, name, cycles):
        """Return column data from specific cycles.

        :param str name: Data name
        :param int or list(int) cycles: cycle numbers
        """
        if isinstance(cycles, int):
            return self._get_dcol(name, cycles)
        return [self._get_dcol(name, c) for c in cycles]


class MesaDataText(MesaPlotMixin):
    """
    Class for MESA data obtained from reading text files.

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
    """

    def __init__(self, directory,
                 history_name='history.data',
                 profile_index_name='profiles.index',
                 profile_prefix='profile', profile_suffix='.data'):

        self._directory = directory
        self._profile_prefix = profile_prefix
        self._profile_suffix = profile_suffix

        # history.data read-in
        self.history_data = DataFromTextMixin(
            os.path.join(self._directory, history_name), TextFile.MESA_DATA)

        # We need to create a mapping between the cycles and the array indices
        # since the cycle number is what the user provides
        self._cycle_index_mapping = {c: i[0] for i, c in np.ndenumerate(self.cycles)}

        # Try to read profiles.index
        try:
            self.profiles_index = DataFromTextMixin(
                os.path.join(self._directory, profile_index_name), TextFile.PROFILES_INDEX)

            # Mapping between the cycle number and lof file number
            self._cycle_log_file_number_mapping = dict(
                zip(self.profiles_index.get('model_number'),
                    self.profiles_index.get('log_file_number'))
            )

            # If the profile file is not found, remove the cycle from the mapping.
            cycles_in_profile_index = sorted(self._cycle_log_file_number_mapping.keys())
            for cycle in cycles_in_profile_index:
                if not os.path.isfile(self._get_cycle_data_filepath(cycle)):
                    _ = self._cycle_log_file_number_mapping.pop(cycle)

            self.cycles_with_data = sorted(self._cycle_log_file_number_mapping.keys())

            # Read one profile to extract column names
            # Take the smallest cycle
            if self.cycles_with_data:
                self.dcols = self._get_cycle(min(self.cycles_with_data)).data_names

        except OSError:
            logger.warning('Cannot extract profiles from directory %s. No profile data availabe.',
                           self._directory)

    @property
    def cycles(self):
        """List of cycles."""
        return sorted(self.history_data.get('model_number'))

    @property
    def hattrs(self):
        """List of header attributes."""
        return sorted(self.history_data.header_names)

    @property
    def cattrs(self):
        """List of cycles headers."""
        return sorted(self.history_data.data_names)

    def get_hattr(self, name):
        """Return header attribute value(s).

        :param str name: Metadata name.
        """
        return self.history_data.get_header(name)

    def get_cattr(self, name):
        """Return header data from specific cycles.

        :param str name: Header name
        """
        return self.history_data.get(name)

    def _get_cycle_data_filepath(self, cycle):
        """Returns theoretical path to the file containing the cycle data."""
        file_number = self._cycle_log_file_number_mapping[cycle]
        filename = ''.join([self._profile_prefix, str(file_number), self._profile_suffix])
        return os.path.join(self._directory, filename)

    def clear_cache(self):
        """Clear cached data, namely profiles."""
        self._get_cycle.cache_clear()

    @lru_cache(maxsize=CACHE_MAX_SIZE)
    def _get_cycle(self, cycle):
        """Returns cycle data.

        ..note:: this function uses caching.
        """
        filepath = self._get_cycle_data_filepath(cycle)

        logger.info('Reading data for cycle %s from file %s', cycle, filepath)
        return DataFromTextMixin(filepath, TextFile.MESA_DATA)

    def get_cycle(self, cycle):
        """Extract cycle data from the corresponding text file, checking if the latter exists."""

        if cycle not in self.cycles_with_data:
            raise RuntimeError(
                'No profile data available for cycle {}, must be one of these: {}'.format(
                    cycle, self.cycles_with_data))
        return self._get_cycle(cycle)

    def get_dcol(self, name, cycles):
        """Return column data from specific cycles.

        :param str name: Data name
        :param int or list(int) cycles: cycle numbers
        """
        if isinstance(cycles, int):
            return self.get_cycle(cycles).get(name)
        return [self.get_cycle(c).get(name) for c in cycles]


class MesaDataHDF5(DataFromHDF5Mixin, MesaPlotMixin):
    """Class for MESA data obtained from reading HDF5 files.

     Arguments:
        filedir (str): Directory path to MESA data.
        file_ext (str, optional): file extension
    """


class MesaData:
    """Class for MESA data.

    :param str data_type: type of the MESA data ('text' or 'hdf5')

    .. note:: Acts like a proxy
    """

    def __init__(self, data_type, *args, **kwargs):

        if data_type == 'text':
            self._proxied_obj = MesaDataText(*args, **kwargs)
        elif data_type == 'hdf5':
            self._proxied_obj = MesaDataHDF5(*args, **kwargs)
        else:
            raise RuntimeError(
                "Wrong type of MESA data requested, must be 'text' or 'hdf5'")

    # Delegate calls to the proxied object:
    def __getattr__(self, name):
        return getattr(self._proxied_obj, name)

    # For tab completion
    def __dir__(self):
        return self._proxied_obj.__dir__() + ['_proxied_obj']


class NugridData(DataFromHDF5Mixin, NugridPlotMixin):
    """Class for NuGrid data.

     Arguments:
        filedir (str): Directory path to NuGrid data.
        file_ext (str, optional): file extension
    """


class Trajectory(DataFromTextMixin):
    """Read in a trajectory file and make a DataFromTextMixin instance."""

    def __init__(self, filename):
        super().__init__(filename, TextFile.TRAJECTORY)
