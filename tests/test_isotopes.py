from unittest import TestCase

from nugridpy.isotopes import ISOTOPES, get_isotope_name, get_A_Z


class TestIsotopes(TestCase):
    """Test for the isotopes modules."""

    def test_isotopes_dict(self):
        "Check the isotopes dictionary."

        # Keys should be distinct and consecutive
        keys = ISOTOPES.keys()
        self.assertSetEqual(set(keys), set(range(len(keys))))

        # Values should be unique
        values = ISOTOPES.values()
        self.assertEqual(len(values), len(set(values)))

    def test_get_isotope_name(self):
        """Check the function computing the isotope name."""

        self.assertEqual(get_isotope_name(1, 1), 'H-1')
        self.assertEqual(get_isotope_name(2, 1), 'H-2')
        self.assertEqual(get_isotope_name(210, 83), 'Bi-210')
        self.assertEqual(get_isotope_name(210, 84), 'Po-210')

        with self.assertRaises(KeyError):
            get_isotope_name(210, 120)

    def test_get_A_Z(self):
        """Check the function computing A and Z from the name."""

        self.assertEqual(get_A_Z('H-1'), (1, 1))
        self.assertEqual(get_A_Z('H-2'), (2, 1))
        self.assertEqual(get_A_Z('Bi-210'), (210, 83))
        self.assertEqual(get_A_Z('Po-210'), (210, 84))

        with self.assertRaises(KeyError):
            get_A_Z('Po210')

        with self.assertRaises(KeyError):
            get_A_Z('po-210')
