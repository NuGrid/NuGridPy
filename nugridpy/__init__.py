from .data import MesaData, NugridData


def mesa_data(directory, data_type,
              history_name='history.data',
              profile_index_name='profiles.index',
              profile_prefix='profile',
              profile_suffix='.data'):
    """Utility function acting as a proxy to create MESA data.

    :param str directory: Directory path to MESA data.
    :param str data_type: type of the MESA data ('text' or 'hdf5')

    If data_type == 'text', additional kwargs can be supplied:
    :param str history_name: name of the history data file
    :param str profile_index_name: name of the profile index file
    :param str profile_prefix: prefix of the profile files
    :param str profile_suffix: suffix of the profile files
    """
    return MesaData(
        data_type, directory, history_name=history_name,
        profile_index_name=profile_index_name, profile_prefix=profile_prefix,
        profile_suffix=profile_suffix)


def nugrid_data(directory, file_ext='.h5'):
    """Utility function acting as a proxy to create NuGrid data.

    :param str directory: Directory path to NuGrid data.
    :param str file_ext: file extension.
    """
    return NugridData(directory, file_ext=file_ext)
