"""
===================
ASTRONOMY CONSTANTS
===================

Physical constants for astronomy & astrophysics

Classes
=======

Constant:
    Defines a physical constant for use in astronomical or physical computation.
    A Constant instance is defined by its value, description, and units attributes.
"""


class Constant(float):
    """
    Creates a physical or astronomical constant instance, which inherits the float
    properties of the value argument but has a description string and unit string
    as attributes.

    Attributes:
        value (float)     : The numerical value of the constant
        description (str) : A name and/or physical description of the constant
        units (str)       : The units of the constant (i.e. [Constant])

    """

    def __new__(cls, value, description, units):
        """
        Arguments:
            value (float)     : The numerical value of the constant
            description (str) : A name and/or physical description of the constant
            units (str)       : The units of the constant (i.e. [Constant])
        """

        # Raise an error if description/units attributes are not strings
        if not isinstance(description, str) or not isinstance(units, str):
            raise TypeError(
                'Constant description and units must be entered as strings')

        instance = super().__new__(cls, value)
        instance.value = value
        instance.description = description
        instance.units = units
        return instance

    def __repr__(self):
        # Ensures constants are represented by a custom long-form description,
        # not just values
        return 'Constant({self.value}, {self.description}, {self.units})'.format(self=self)


# Astronomical constants
r_sun = Constant(
    6.9598e10, 'r_sun: Solar radius', 'cm')
l_sun = Constant(
    3.8418e33, 'l_sun: Solar luminosity', 'erg/s')
m_sun = Constant(
    1.9892e33, 'm_sun: Solar mass', 'g')
m_earth = Constant(
    5.9764e27, 'm_earth: Earth mass', 'g')
au_cm = Constant(
    1.495978921e13, 'au_cm: 1 AU (Astronomical Unit), or distance from Earth to Sun', 'cm')
speed_light = Constant(
    2.99792458e10, 'speed_light: Speed of light in vacuum (c)', 'cm/s')
grav_const = Constant(
    6.67428e-8, 'grav_const: Gravitational constant (G)', 'cm^3 g^-1 s^-2')

# Physical constants
stef_boltz_const = Constant(
    5.670400e-5, 'stef_boltz_const: Stefan-Boltzmann constant (sigma)', 'erg cm^-2 s^-1 K^-4')
mass_H_atom = Constant(
    1.6726231e-24, 'mass_H_atom: Mass of hydrogen atom', 'g')
mass_electron = Constant(
    9.1093826e-28, 'mass_electron: Mass of the electron', 'g')
planck_const = Constant(
    6.62606896E-27, 'planck_const: Planck\'s constant (h)', 'erg s')
atomic_mass_unit = Constant(
    1.660538782e-24, 'atomic_mass_unit: Atomic mass unit (u)', 'g')
boltzmann_const = Constant(
    1.3806504e-16, 'boltzmann_const: Boltzmann constant (k)', 'erg K^-1')
avogadro_const = Constant(
    6.02214179e23, 'avogadro_const: Avogadro\'s number', 'mol^-1')
rad_const = Constant(
    4 * stef_boltz_const / speed_light,
    'rad_const: Radiation constant, total energy radiated by a blackbody',
    'erg cm^-3 K^-4')
day_secs = Constant(
    86400., 'day_secs: Number of seconds in a day', 's')

# Miscellaneous coefficients and other constants for astronomy functions
visc_mol_const = Constant(
    1.84e-17, 'visc_mol_const: Constant for calculating molecular diffusivity', 'g cm^-1 K^-2.5')
nu_rad_const = Constant(
    6.88e-26, 'nu_rad_const: Constant for calculating radiative diffusivity', 'g^2 cm^-2 K^-4')
core_mass_coeff = Constant(
    6.22e5, 'core_mass_coeff: Core mass multiplicative constant', 'Dimensionless')
core_mass_offset = Constant(
    0.487, 'core_mass_offset: Core mass additive constant', 'M_sun')
idrad_const = Constant(
    3.2e7, 'idrad_const: radiation/gas pressure boundary', 'K cm')
iddeg_const = Constant(
    1.207e5, 'iddeg_const: ideal gas/degenerate pressure boundary', 'K')

# Polynomial coefficients for Ferrario (2005) initial-final mass relation
mimf_coeff_6 = Constant(
    -0.00012336, 'mimf_coeff_6: 6th power IFMR coefficient', 'M_sun^-5')
mimf_coeff_5 = Constant(
    0.003160, 'mimf_coeff_5: 5th power IFMR coefficient', 'M_sun^-4')
mimf_coeff_4 = Constant(
    0.02960, 'mimf_coeff_4: 4th power IFMR coefficient', 'M_sun^-3')
mimf_coeff_3 = Constant(
    0.12350, 'mimf_coeff_3: 3rd power IFMR coefficient', 'M_sun^-2')
mimf_coeff_2 = Constant(
    0.21550, 'mimf_coeff_2: 2nd power IFMR coefficient', 'M_sun^-1')
mimf_coeff_1 = Constant(
    0.19022, 'mimf_coeff_1: 1st power IFMR coefficient', 'Dimensionless')
mimf_coeff_0 = Constant(
    0.46575, 'mimf_coeff_0: IFMR constant term', 'M_sun')

# Kroupa IMF mass bounds and alpha values
imf_m1 = Constant(
    0.08, 'imf_m1: Lower bound for IMF mass range', 'M_sun')
imf_m2 = Constant(
    0.50, 'imf_m2: Middle bound for IMF mass range', 'M_sun')
imf_a1 = Constant(
    0.30, 'imf_a1: Kroupa IMF alpha value', 'Dimensionless')
imf_a2 = Constant(
    1.30, 'imf_a2: Kroupa IMF alpha value', 'Dimensionless')
imf_a3 = Constant(
    2.3, 'imf_a3: Kroupa IMF alpha value', 'Dimensionless')

# van Loon (2005) mass loss rate constants
van_loon_1 = Constant(
    -5.65, 'van_loon_1: van Loon (2005) mass loss rate constant', 'M_sun yr^-1')
van_loon_2 = Constant(
    6.3, 'van_loon_2: van Loon (2005) mass loss rate constant', 'M_sun yr^-1')
van_loon_3 = Constant(
    3500., 'van_loon_3: van Loon (2005) mass loss rate constant', 'K')
