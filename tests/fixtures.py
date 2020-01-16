import random
import string

import numpy as np


def random_string(length=None):
    """ Generate a random string"""

    length = length or random.randint(5, 64)
    return ''.join(
        random.choice(string.ascii_letters + string.digits)
        for _ in range(length))


def random_array(dim1=None, *args):
    """Generate a random array of data."""

    dim1 = dim1 or random.randint(1, 64)
    return np.random.rand(dim1, *args)


def random_ints(length=None):
    """Generate a random list of integers."""

    length = length or random.randint(5, 64)
    return [random.randint(1, 1000) for _ in range(length)]
