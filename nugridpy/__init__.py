from .data import MesaData, NugridData


def mesa_data(directory, data_type, *args, **kwargs):
    """Utility function acting as a proxy to create MESA data.

    :param str directory: Directory path to MESA data.
    :param str data_type: type of the MESA data ('text' or 'hdf5')
    """
    return MesaData(data_type, directory, *args, **kwargs)


def nugrid_data(directory, *args, **kwargs):
    """Utility function acting as a proxy to create NuGrid data.

    :param str directory: Directory path to NuGrid data.
    """
    return NugridData(directory, *args, **kwargs)
