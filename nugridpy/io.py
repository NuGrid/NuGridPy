"""
Io
==

This module contains the logic for reading/writing data.
"""

import numpy as np

import h5py

from .config import logger


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


def parse_trajectory(filename):
    """
    For parsing trajectory files.
    """
    # Read in the header data as an array (returns unstructured)
    traj = np.genfromtxt(filename, dtype=None, max_rows=4, invalid_raise=True,
                         encoding=None, delimiter='=', autostrip=True)

    # Make a structured dtype for this file
    traj_dtype = np.dtype([(el[0], traj.dtype) for el in traj])

    # Produce the structured header array and return it
    header_data = np.array([(tuple(el[1] for el in traj))], dtype=traj_dtype)
    return header_data


class TextFile:
    """
    Class that contains ASCII i/o functionality. Class variables are
    used to represent filetypes as a dictionary, with relevant info
    that is specific to each filetype housed as key-value pairs.

    Attributes:
        MESA_DATA (dict): The specific key-value pairs that represent
            a history.data file.
        PROFILES_INDEX (dict): The specific key-value pairs that represent
            a profiles.index file.
        FILE_TYPES (tuple): Tuple containing each filetype class variable.
    """

    MESA_DATA = {
        'name': 'MESA_DATA',
        'parser': parse_history,
        'skip_header': 5,
        'names': None
    }

    PROFILES_INDEX = {
        'name': 'PROFILES_INDEX',
        'parser': None,
        'skip_header': 1,
        'names': ['model number', 'priority', 'log file number']
    }

    TRAJECTORY = {
        'name': 'TRAJECTORY',
        'parser': parse_trajectory,
        'skip_header': 7,
        'names': ['time', 'T', 'rho']
    }

    FILE_TYPES = (
        MESA_DATA,
        PROFILES_INDEX,
        TRAJECTORY
    )

    @classmethod
    def readfile(cls, filename, filetype):
        """
        Reads in ascii files.

        Arguments:
            filename (str): Path to file being read in.

        Returns:
            header_data (np.array): Structured array containing
                header names and header data in dict-like format.
            data (np.array): Structured array containing data column
                names and data from data fields in dict-like format.
        """
        # Declare the header data
        if filetype['parser'] is not None:
            header_data = filetype['parser'](filename)
        else:
            header_data = np.ndarray(0)

        # Declare the column data
        data = np.genfromtxt(filename, dtype=None,
                             skip_header=filetype['skip_header'], names=filetype['names'] or True,
                             invalid_raise=True, encoding=None)

        # Return the header data and column data together
        return header_data, data


class Hdf5File:
    """
    Class containing the HDF5 I/O functionality. Structure mirrors the Textfile
    class, and a similar readfile class method is implemented. Class variables
    are used to represent the different sources of HDF5 file that may be
    read in.
    """

    NUGRID = {
        'default_dataset_name': 'data',
        'group_default_dataset_name': 'SE_DATASET',
        'group_prefix': 'cycle'
    }

    @classmethod
    def readfile(cls, filename, filetype):
        """
        Reads in HDF5 files.
        Returns metadata, data objects from file, and more?
        """
        logger.info('Reading %s...', filename)
        file_object = h5py.File(filename, 'r')

        # Metadata
        metadata = file_object.attrs

        groups, headers, datasets = {}, {}, {}
        # Get all groups and datasets
        for key in file_object.keys():
            extracted_data = file_object[key]

            # Check if it is a group that should contain a single dataset or a dataset
            if isinstance(extracted_data, h5py.Dataset):
                datasets[key] = extracted_data[filetype['default_dataset_name']]
            elif isinstance(extracted_data, h5py.Group):
                # Extract cycle number from name
                cycle_number = int(key.strip(filetype['group_prefix']))
                groups[cycle_number] = extracted_data[filetype['group_default_dataset_name']]
                headers[cycle_number] = extracted_data.attrs
            else:
                raise RuntimeError("HDF5 file appears to be unsupported...")

        # Return the data structures
        return metadata, groups, headers, datasets
