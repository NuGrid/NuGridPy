"""
Data
====

This module provides implementation for handling data.
"""

from contextlib import suppress
import logging
import os
import re

import numpy as np

from .plot import PlotMixin, StarPlotsMixin
from .io import TextFile, Hdf5File
from .isotopes import get_isotope_name, get_A_Z

logging.basicConfig(level=logging.INFO)


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


class DataFromHDF5Mixin(PlotMixin):
    """Class for data obtained from reading text files.

     Arguments:
        filedir (str): Directory path to the data.
        file_ext (str, optional): file extension
    """

    # Name of the member containing the mass fractions
    MASS_FRACTION_HDF5_MEMBER_NAME = 'iso_massf'

    def __init__(self, filedir, file_ext='.h5'):

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
                    logging.info(
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
        """List of cycles."""
        return tuple(self._data.keys())

    @property
    def metadata_names(self):
        """List of metadata names."""
        return tuple(self._metadata.keys())

    @property
    def cycle_headers(self):
        """List of cycles headers."""
        return tuple(self._headers[self.cycles[0]].keys())

    @property
    def cycle_data(self):
        """List of cycles data."""
        return tuple(self._data[self.cycles[0]].dtype.names)

    def get_metadata(self, name):
        """Return metadata value(s).

        :param str name: Metadata name.
        """
        try:
            return self._metadata[name]
        except KeyError:
            raise KeyError('Cannot find attribute {} in metadata.'.format(name))

    def _get_cycle_header(self, name, cycle):
        """Return data from a specific cycle header.

        :param str name: Header name
        :param int cycle: cycle number
        """
        try:
            return self._headers[cycle][name].item()
        except KeyError:
            raise KeyError('Cannot find header {} in cycle {}.'.format(name, cycle))

    def get_cycle_header(self, name, cycles=None):
        """Return header data from specific cycles.

        :param str name: Header name
        :param int or list(int) cycles: cycle numbers. Default: all cycles
        """
        cycles = cycles or self.cycles
        if isinstance(cycles, int):
            return self._get_cycle_header(name, cycles)
        return np.array([self._get_cycle_header(name, c) for c in cycles])

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
            a, z = get_A_Z(name)

            # Extract data from the single index in A,Z with these values
            index = int(np.where((self.A == a) * (self.Z == z))[0])
            return self._get_cycle(cycle)[self.MASS_FRACTION_HDF5_MEMBER_NAME][..., index]
        except (TypeError, KeyError):
            raise KeyError('Cannot find isotope {} in cycle {}.'.format(name, cycle))

    def _get_cycle_data(self, name, cycle):
        """Return data from a specific cycle.

        :param str name: Data name
        :param int or list(int) cycle: cycle number
        """
        try:
            return self._get_cycle(cycle)[name]
        except ValueError:

            # It might be an element that can be accessed through iso_massf
            if name in self.isotopes:
                return self._get_isotope_mass_fraction(name, cycle)
            raise KeyError('Cannot find data {} in cycle {}. '.format(name, cycle))

    def get_cycle_data(self, name, cycles):
        """Return data from specific cycles.

        :param str name: Data name
        :param int or list(int) cycles: cycle numbers
        """
        if isinstance(cycles, int):
            return self._get_cycle_data(name, cycles)
        return [self._get_cycle_data(name, c) for c in cycles]


# This is not ideal but for the moment,
# we keep two distinct classes for MESA data
# because the data and the way to access it is different
# between these obtained from logs and those from HDF5 files.
class MesaDataText(StarPlotsMixin):
    """
    Class for MESA data obtained from reading text files.

    Attributes:
        history_data (DataFromTextMixin): A DataFromTextMixin instance that contains the data from a
            history.data file.
        profiles_index (DataFromTextMixin): A DataFromTextMixin instance that contains the data from
            a profiles.index file.
        profiles_in_dir (tuple): A tuple containing string filenames for each MESA
            profile found in filedir.
        profiles (list): A list populated with one DataFromTextMixin object for each profile
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
            self.history_data = DataFromTextMixin(
                os.path.join(filedir, history_name), TextFile.MESA_DATA)

            # We need to create a mapping between the cycles and the array indices
            # since the cycle number is what the user provides
            self._cycle_index_mapping = {c: i[0] for i, c in np.ndenumerate(self.cycles)}

        # profiles.index read-in
        if not history_only:
            # Try to read and alert user if file not present
            self.profiles_index = DataFromTextMixin(
                os.path.join(filedir, index_name), TextFile.PROFILES_INDEX)

            # Collect the profiles from filedir
            self.profiles_in_dir = tuple(prof for prof in os.listdir(filedir) if re.match(
                '{}[0-9]\d*{}'.format(profile_prefix, profile_suffix), prof))
            prof_count = len(self.profiles_in_dir)

            # Read in all profiles to a user-facing list attribute
            if all_profiles:
                print('Reading in {} profiles... set all_profiles to False to '
                      'skip.'.format(prof_count))
                self.profiles = {}
                for prof in self.profiles_in_dir:
                    data = DataFromTextMixin(os.path.join(filedir, prof), TextFile.MESA_DATA)
                self.profiles[data.get_header('model_number')] = data

    @property
    def cycles(self):
        """List of cycles."""
        return self.history_data.get('model_number')

    def get_cycle_header(self, name, cycles=None):
        """Return header data from specific cycles.

        :param str name: Header name
        :param int or list(int) cycles: cycle numbers. Default: all cycles
        """
        cycles = cycles or self.cycles
        if isinstance(cycles, int):
            return self.history_data.get(name)[self._cycle_index_mapping[cycles]]
        return np.array([self.history_data.get(name)[self._cycle_index_mapping[c]] for c in cycles])

    def get_cycle_data(self, name, cycles):
        """Return data from specific cycles.

        :param str name: Data name
        :param int or list(int) cycles: cycle numbers
        """
        if isinstance(cycles, int):
            return self.profiles[cycles].get(name)
        return [self.profiles[c].get(name) for c in cycles]


class MesaDataHDF5(DataFromHDF5Mixin, StarPlotsMixin):
    """Class for MESA data obtained from reading HDF5 files.

     Arguments:
        filedir (str): Directory path to MESA data.
        file_ext (str, optional): file extension
    """


class NugridData(DataFromHDF5Mixin, StarPlotsMixin):
    """Class for NuGrid data.

     Arguments:
        filedir (str): Directory path to NuGrid data.
        file_ext (str, optional): file extension
    """


class Trajectory(DataFromTextMixin):
    """Read in a trajectory file and make a DataFromTextMixin instance."""

    def __init__(self, filename):
        super().__init__(filename, TextFile.TRAJECTORY)
