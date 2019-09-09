import random
import string


def random_string(length=None):
    """
    Generate a random string of random length using upper/lower case
    letters and digits.

    """

    length = length or random.randint(5, 64)
    return ''.join(
        random.choice(string.ascii_letters + string.digits)
        for _ in range(length))
