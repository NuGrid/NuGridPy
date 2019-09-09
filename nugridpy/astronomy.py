"""
===================
ASTRONOMY FUNCTIONS
===================

Useful functions for astronomy & astrophysics

"""

from functools import update_wrapper

import numpy as np
from scipy import integrate

from . import constants as cs


class ReadOnlyConstants:
    """Callable class for attaching constants as read-only property to a function."""

    def __init__(self, constants, func):
        """Constructor that defines function and constants in class instance."""
        self._constants = constants
        self.func = func

    def __call__(self, *args, **kwargs):
        """Defines the class as a callable and executes the decorated function."""
        return self.func(*args, **kwargs)

    @property
    def constants(self):
        """Returns constants as private attribute."""
        return self._constants


def attach_constants(*args):
    """Decorator receives function constants first, then attaches them through a callable class."""

    def attach(func):
        function_with_constants = ReadOnlyConstants(args, func)
        # inherit docstring and other magic info from original function
        return update_wrapper(function_with_constants, func)

    return attach


@attach_constants(cs.visc_mol_const)
def visc_mol_sol(T, rho, X):
    '''
    Molecular plasma viscosity (Spitzer 1962)

    Parameters
    ----------
    T : float
       temperature in K
    rho : float
       density in cgs
    X : float
       H mass fraction

    Returns
    -------
    nu : float
       molecular diffusivity in [cm**2/s]

    Notes
    -----
    According to Eq 22 in Schatzman (1977). Assume log Lambda = 15.
    (see Table 5.1), a H/He mix (for different mix use Eq. 5.54 in
    Spitzer textbook)

    '''

    visc_mol = cs.visc_mol_const * (1. + (7.*X)) * (T**2.5 / rho)
    return visc_mol


@attach_constants(cs.nu_rad_const)
def visc_rad_kap_sc(T, rho, X):
    '''
    Radiative viscosity (Thomas, 1930) for e- scattering opacity

    Parameters
    ----------
    T : float
       temperature in K
    rho : float
       density in cgs
    X : float
       H mass fraction

    Returns
    -------
    nu : float
       radiative diffusivity in [cm**2/s]

    Notes
    -----
    Eqn. 14 in Schatzman, 1977, assume electron scattering opacity
    kappa_sc = 0.2*(1+X), Kippenhahn (2nd edn, Eqn 17.2)

    '''

    kappa = 0.2 * (1.+X)
    nu_rad = cs.nu_rad_const * (T**4 / (kappa * rho**2))
    return nu_rad


@attach_constants()
def Gamma1_gasrad(beta):
    '''
    Gamma1 for a mix of ideal gas and radiation

    Hansen & Kawaler, page 177, Eqn. 3.110

    Parameters
    ----------
    beta : float
        Gas pressure fraction Pgas/(Pgas+Prad)

    '''

    Gamma3minus1 = (2./3.) * (4. - (3.*beta)) / (8. - (7.*beta))
    Gamma1 = beta + (4. - (3.*beta)) * Gamma3minus1
    return Gamma1


@attach_constants(cs.boltzmann_const, cs.atomic_mass_unit)
def Pgas(rho, T, mmu):
    '''
    P = R/mu * rho * T

    Parameters
    ----------
    rho : float
        Density [cgs]
    T : float
        Temperature [K]
    mmu : float
        Mean molecular weight

    Returns
    --------
    Gas pressure

    '''

    R = cs.boltzmann_const / cs.atomic_mass_unit
    return (R/mmu) * rho * T


@attach_constants(cs.rad_const)
def Prad(T):
    '''
    P = cs.rad_const/3 * T**4

    Parameters
    ----------
    T : float
        Temperature [K]

    Returns
    --------
    Radiation pressure

    '''

    return (cs.rad_const/3.) * T**4


@attach_constants(cs.mimf_coeff_6, cs.mimf_coeff_5, cs.mimf_coeff_4,
                  cs.mimf_coeff_3, cs.mimf_coeff_2, cs.mimf_coeff_1, cs.mimf_coeff_0)
def mimf_ferrario(mi):
    ''' Curvature MiMf from Ferrario et al. 2005MNRAS.361.1131.'''

    mf = ((cs.mimf_coeff_6 * (mi**6)) + (cs.mimf_coeff_5 * (mi**5))
          - (cs.mimf_coeff_4 * (mi**4)) + (cs.mimf_coeff_3 * (mi**3))
          - (cs.mimf_coeff_2 * (mi**2)) + (cs.mimf_coeff_1 * mi) + cs.mimf_coeff_0)
    return mf


@attach_constants(cs.core_mass_coeff, cs.core_mass_offset)
def core_mass_L(MH):
    '''
    Core-mass luminosity relationship from Bloecker (1993)

    Parameters
    ----------
    MH : float
        Core mass in Msun

    Returns
    -------
    L
        Luminosity in Lsun

    '''

    return cs.core_mass_coeff * (MH - cs.core_mass_offset)


@attach_constants(cs.imf_m1, cs.imf_m2, cs.imf_a1, cs.imf_a2, cs.imf_a3)
def imf(m):
    '''
    Initial mass function from Kroupa

    Parameters
    -------
    m : float
        mass (g)

    Returns
    -------
    N(M)dM
        for given mass according to Kroupa IMF

    '''

    const2 = cs.imf_m1**(-cs.imf_a1) - cs.imf_m1**(-cs.imf_a2)
    const3 = cs.imf_m2**(-cs.imf_a2) - cs.imf_m2**(-cs.imf_a3)

    if m < cs.imf_m1:
        alpha = cs.imf_a1
        const = -const2 - const3
    elif m < cs.imf_m2:
        alpha = cs.imf_a2
        const = -const3
    else:
        alpha = cs.imf_a3
        const = 0.
    return m**(-alpha) + const


@attach_constants()
def int_imf_dm(m1, m2, m, imf_ar, bywhat='bymass', integral='normal'):
    '''
    Integrate IMF between m1 and m2

    Parameters
    ----------
    m1 : float
        Lower mass integration bound
    m2 : float
        Upper mass integration bound
    m : array
        Mass array
    imf_ar : array
        Array of IMF values corresponding to mass array
    bywhat : string, optional
        'bymass' integrates the mass that goes into stars of
        that mass interval; or 'bynumber' which integrates the number
        of stars in that mass interval.  The default is 'bymass'.
    integrate : string, optional
        'normal' uses scipy.integrate.trapz; 'cum' returns cumulative
        trapezoidal integral.  The default is 'normal'.

    Returns
    ---------
    Integrated initial mass function for given bounds

    '''

    ind_m = (m >= min(m1, m2)) & (m <= max(m1, m2))

    if integral == 'normal':
        int_func = integrate.trapz
    elif integral == 'cum':
        int_func = integrate.cumtrapz
    else:
        raise ValueError(
            "Error in int_imf_dm: don't know how to integrate (normal or cum)")

    if bywhat == 'bymass':
        return int_func(m[ind_m] * imf_ar[ind_m], m[ind_m])
    elif bywhat == 'bynumber':
        return int_func(imf_ar[ind_m], m[ind_m])
    raise ValueError(
        "Error in int_imf_dm: Need integration type (bymass or bynumber)")


@attach_constants(cs.r_sun, cs.m_sun, cs.grav_const)
def am_orb(m1, m2, a, e):
    '''
    Orbital angular momentum equation

    e.g. Ge et al 2010

    Parameters
    ----------
    m1, m2 : float
        Masses of both stars in Msun
    A : float
        Separation in Rsun
    e : float
        Eccentricity

    Returns
    --------
    Orbital angular momentum

    '''

    a_cm = a * cs.r_sun
    m1_g = m1 * cs.m_sun
    m2_g = m2 * cs.m_sun
    J_orb = np.sqrt(cs.grav_const * a_cm * ((m1_g**2 * m2_g**2) / (m1_g + m2_g))) * (1 - e**2)
    return J_orb


@attach_constants(cs.van_loon_1, cs.van_loon_2, cs.van_loon_3)
def mass_loss_loon05(L, Teff):
    '''
    Mass loss rate from van Loon et al (2005)

    Parameters
    ----------
    L : float
        L in L_sun
    Teff : float
        Teff in K

    Returns
    -------
    Mdot
        Mdot in Msun/yr

    Notes
    -----
    ref: van Loon etal 2005, A&A 438, 273

    '''

    Mdot = (cs.van_loon_1 + np.log10(L / 10.**4) -
            cs.van_loon_2 * np.log10(Teff / cs.van_loon_3))
    return Mdot


@attach_constants(cs.grav_const, cs.m_sun, cs.r_sun)
def energ_orb(m1, m2, r):
    '''
    Orbital potential energy equation

    Parameters
    ----------
    m1, m2 : float
        M in Msun
    r : float
        Distance in Rsun

    Returns
    -------
    Epot
        Epot in erg

    '''

    epo = -cs.grav_const * m1 * m2 * cs.m_sun**2 / (cs.r_sun * r)
    return epo


@attach_constants(cs.r_sun, cs.grav_const, cs.m_sun, cs.day_secs)
def period(A, M1, M2):
    """
    Calculate binary period from separation.

    Parameters
    ----------
    A : float
        separation A Rsun
    M1, M2 : float
        M in Msun

    Returns
    -------
    p
        period in days

    """

    A *= cs.r_sun
    velocity = np.sqrt(cs.grav_const * cs.m_sun * (M1+M2) / A)
    p = ((2. * np.pi * A) / velocity) / cs.day_secs
    return p


@attach_constants(cs.grav_const, cs.m_sun, cs.r_sun)
def escape_velocity(M, R):
    """
    Escape velocity

    Parameters
    ----------
    M : float
        Mass in solar masses
    R : float
        Radius in solar radii

    Returns
    -------
    v_escape
        in km/s

    """

    ve = np.sqrt(2. * cs.grav_const * M * cs.m_sun / (R * cs.r_sun))
    ve = ve * 1.e-5
    return ve


@attach_constants(cs.avogadro_const, cs.boltzmann_const, cs.mass_H_atom)
def Nasv(macs_val, T):
    '''
    Returns
    -------
    Na*<sigma v>
        for MACS [mb] at T [K]

    '''

    Na = cs.avogadro_const
    k = cs.boltzmann_const
    vtherm = (2. * k * T / cs.mass_H_atom)**0.5
    s = macs_val * 1.e-27
    Nasv_val = s * vtherm * Na
    return Nasv_val


@attach_constants(cs.avogadro_const, cs.boltzmann_const, cs.mass_H_atom)
def macs(nasv, T):
    '''
    Returns
    -------
    MACS
        [mb] at T [K] from Na*<sigma v>

    '''

    Na = cs.avogadro_const
    k = cs.boltzmann_const
    vtherm = (2. * k * T / cs.mass_H_atom)**0.5
    s = nasv / (vtherm * Na)
    macs_val = s * 1.e27
    return macs_val


@attach_constants()
def mu_e(X):
    '''
    Mean molecular weight per free electron, assuming full ionisation, and
    approximating mu_i/Z_i ~ 2 for all elements heavier then Helium.

    (Kippenhahn & Weigert, Ch 13.1, Eq. 13.8)

    Parameters
    ----------
    X : float
        Mass fraction of H

    '''

    try:
        mu_el = 2. / (1.+X)
    except TypeError:
        X = np.array([X])
        mu_el = 2. / (1.+X)

    return mu_el


@attach_constants()
def mu(X, Z, A):
    '''
    Mean molecular weight assuming full ionisation.

    (Kippenhahn & Weigert, Ch 13.1, Eq. 13.6)

    Parameters
    ----------
    X : float or array
        Mass fraction vector
    Z : float or array
        Charge number vector
    A : float or array
        Mass number vector

    '''

    if not isinstance(Z, np.ndarray):
        Z = np.array(Z)
    if not isinstance(A, np.ndarray):
        A = np.array(A)
    if not isinstance(X, np.ndarray):
        X = np.array(X)

    try:
        mmu = 1. / sum(X * (1.+Z) / A)
    except TypeError:
        X = np.array([X])
        A = np.array([A])
        Z = np.array([Z])
        mmu = 1. / sum(X * (1.+Z) / A)

    return mmu


@attach_constants(cs.idrad_const)
def Trho_idrad(rho, mmu):
    '''
    T(rho) that separates P_rad from P_gas dominated regions.

    Kippenhahn & Weigert, Eq. 16.10

    Parameters
    ----------
    rho : float
        Density array [cgs]
    mu : float
        Mean molecular weight

    '''

    T = cs.idrad_const * (rho/mmu)**(1./3.)
    return T


@attach_constants(cs.iddeg_const)
def Trho_iddeg(rho, mmu, mu_el):
    '''
    T(rho) that separates ideal gas and degenerate pressure dominated regions.

    Kippenhahn & Weigert, Eq. 16.6

    Parameters
    ----------
    rho : float
        Density array [cgs]
    mmu : float
        Mean molecular weight
    mu_el : float
        Mean molecular weight per free electron

    '''

    T = cs.iddeg_const * rho**(2./3.) * mmu / (mu_el**(5./3.))
    return T
