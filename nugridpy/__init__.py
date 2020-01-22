from .data import MesaDataText, MesaDataHDF5

MESA_DATA_TYPES_CLS = {
    'text': MesaDataText,
    'hdf5': MesaDataHDF5
}


def mesa_data(directory, data_type='text', **_kwargs):
    """Utility function acting as a proxy to create MESA data.

    :param str directory: Directory path to MESA data.
    :param str data_type: Type of the data ('text' or 'hdf5')
    """

    assert data_type in MESA_DATA_TYPES_CLS.keys(), \
        "MESA data type must be one of these: {}".format(MESA_DATA_TYPES_CLS.keys())

    return MESA_DATA_TYPES_CLS[data_type](directory, **_kwargs)
