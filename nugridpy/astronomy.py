'''Constants and methods for astronomy & astrophysics'''
from __future__ import division
from __future__ import print_function
from past.utils import old_div

# Astronomy utilities module

import numpy as np
import scipy as sc
from scipy import integrate


# constants for astronomy
rsun_cm = 6.9598e10
rsun_cm_s = 'cm'
lsun_erg_s = 3.8418e33
lsun_erg_s_unit = 'erg/s'
msun_g = 1.9892e33
msun_unit = 'erg s^-1'
mass_earth_g = 5.9764e27


au_cm = 1.495978921e13
au_cm_unit='cm'

speed_light = 2.99792458e10
speed_light_unit = 'cm s^-1'
grav_const = 6.67428e-8
grav_const_unit = 'cm^3 g^-1 s^-2'

# constants for physics
boltzmann_sigma = 5.670400e-5
mass_H_atom=1.6726231e-24
mass_H_atom_unit='g'
mass_electron=9.1093826e-28
mass_electron_unit='g'
planck_constant_h=6.62606896E-27
planck_constant_h_unit='erg s'
atomic_mass_unit=1.660538782e-24
atomic_mass_unit_unit='g'
boltzmann_constant=1.3806504e-16
boltzmann_constant_unit='erg K^-1'
avogadro_constant=6.02214179e23
avogadro_constant_unit='mol^-1'
radiation_constant = 4*boltzmann_sigma/speed_light
radiation_constant_unit = 'erg cm^-3 K^-4'

def visc_mol_sol(T,rho,X):
    '''
    Molecular plasma viscosity (Spitzer 1962)
    
    Parameters
    ----------
    X : float
       H mass fraction
    T : float
       temperature in K
    rho : float
       density in cgs

    Returns
    -------
    nu : float
       molecular diffusivity in [cm**2/s]

    Notes
    -----
    According to Eq 22 in Schatzman (1977). Assume log Lambda = 15. 
    (see Table 5.1), a H/He mix (for different mix use Eq. 5.54 in
    Spitzer text book)

    Examples
    --------
    see astronomy.visc_rad_kap_sc

    '''
    visc_mol = 1.84e-17*(1.+7.*X)*(old_div(T**2.5,rho))
    return visc_mol


def visc_rad_kap_sc(T,rho,X):
    '''
    Radiative viscosity (Thomas, 1930) for e- scattering opacity

    Parameters
    ----------
    X : float
       H mass fraction
    T : float
       temperature in K
    rho : float
       density in cgs

    Returns
    -------
    nu : float
       radiative diffusivity in [cm**2/s]

    Examples
    --------
    >>> In [1]: import astronomy as ast
    >>> In [2]: l = 100*1.e5 # 100km
    >>> In [3]: v = 1.e5     # typical velocity
    >>> In [4]: T   = 90.e6  # temperature
    >>> In [5]: X   = 0.001  # H mass fraction
    >>> In [6]: rho = 100.   # density
    >>> In [7]: nu = ast.visc_rad_kap_sc(T,rho,X)
    >>> In [8]: Re=v*l/nu
    >>> In [9]: print "Re_rad = "+str('%g'%Re)
    >>> Re_rad = 4.43512e+08

    Notes
    -----
    Eqn. 14' in Schatzman, 1977, assume electron scattering opacity
    kappa_sc = 0.2*(1+X), Kippenhahn (2nd edn, Eqn 17.2)

    '''
    kappa = 0.2*(1.+X)
    nu_rad = 6.88e-26*(old_div(T**4,(kappa*rho**2)))
    return nu_rad

def Gamma1_gasrad(beta):
    ''' 
    Gamma1 for a mix of ideal gas and radiation

    Hansen & Kawaler, page 177, Eqn. 3.110
    
    Parameters
    ----------
    beta : float
        Gas pressure fraction Pgas/(Pgas+Prad)

    '''
    Gamma3minus1 = (old_div(2.,3.))*(4.-3.*beta)/(8.-7.*beta) 
    Gamma1 = beta + (4.-3.*beta) * Gamma3minus1
    return Gamma1

def Pgas(rho,T,mu):
    ''' 
    P = R/mu * rho * T

    Parameters
    ----------
    mu : float
        Mean molecular weight
    rho : float
        Density [cgs]
    T : float
        Temperature [K]

    '''
    R = old_div(boltzmann_constant, atomic_mass_unit)
    return (old_div(R,mu)) * rho * T

def Prad(T,mu):
    ''' 
    P = a/3 * T**4

    Parameters
    ----------
    a : float
        Radiation constant [erg cm^-3 K^-4].
    T : float
        Temperature [K].

    '''
    return (old_div(radiation_constant,3.)) * T**4

def mimf_ferrario(mi):
    ''' Curvature MiMf from Ferrario etal. 2005MNRAS.361.1131.'''
    
    mf=-0.00012336*mi**6+0.003160*mi**5-0.02960*mi**4+\
      0.12350*mi**3-0.21550*mi**2+0.19022*mi+0.46575
    return mf

def core_mass_L(MH):
    ''' 
    Core-mass luminosity relationship from Bloecker (1993).

    Parameters
    ----------
    MH : float
        Core mass in Msun.

    Returns
    -------
    L
        Luminosity in Lsun

    '''
    return 62200*(MH-0.487)

def imf(m):
    ''' 
    Returns
    -------
    N(M)dM
        for given mass according to Kroupa IMF, vectorization
        available via vimf() 

    '''

    m1 = 0.08; m2 = 0.50
    a1 = 0.30; a2 = 1.30; a3 = 2.3
    const2 = m1**-a1 -m1**-a2 
    const3 = m2**-a2 -m2**-a3 

    if m < 0.08:
        alpha = 0.3
        const = -const2 -const3
    elif m < 0.50:
        alpha = 1.3
        const = -const3
    else:
        alpha = 2.3
        const = 0.0
    # print m,alpha, const, m**-alpha + const 
    return m**-alpha + const 

    
vimf = np.vectorize(imf)

def int_imf_dm(m1,m2,m,imf,bywhat='bymass',integral='normal'):
    ''' 
    Integrate IMF between m1 and m2.

    Parameters
    ----------
    m1 : float
        Min mass
    m2 : float
        Max mass
    m : float
        Mass array
    imf : float
        IMF array
    bywhat : string, optional
        'bymass' integrates the mass that goes into stars of
        that mass interval; or 'bynumber' which integrates the number
        of stars in that mass interval.  The default is 'bymass'.
    integrate : string, optional
        'normal' uses sc.integrate.trapz; 'cum' returns cumulative
        trapezoidal integral.  The default is 'normal'.

    '''


    ind_m = (m >= min(m1,m2)) & (m <= max(m1,m2))
    if integral is 'normal':
        int_func = sc.integrate.trapz
    elif integral is 'cum':
        int_func = sc.integrate.cumtrapz
    else:
        print("Error in int_imf_dm: don't know how to integrate")
        return 0
       
    if bywhat is 'bymass':
        return int_func(m[ind_m]*imf[ind_m],m[ind_m])
    elif bywhat is 'bynumber':
        return int_func(imf[ind_m],m[ind_m])
    else:
        print("Error in int_imf_dm: don't know by what to integrate")
        return 0

def am_orb(m1,m2,a,e):
    ''' 
    orbital angular momentum.

    e.g Ge etal2010
    
    Parameters
    ----------
    m1, m2 : float
        Masses of both stars in Msun.
    A : float 
        Separation in Rsun.
    e : float
        Eccentricity
        
    '''

    a_cm  = a * rsun_cm
    m1_g = m1 * msun_g
    m2_g = m2 * msun_g

    J_orb=np.sqrt(grav_const*a_cm*(old_div((m1_g**2*m2_g**2),(m1_g+m2_g))))*(1-e**2)
    return J_orb

def mass_loss_loon05(L,Teff):
    ''' 
    mass loss rate van Loon etal (2005).

    Parameters
    ----------
    L : float
        L in L_sun.
    Teff : float
        Teff in K.
        
    Returns
    -------
    Mdot
        Mdot in Msun/yr
    
    Notes
    -----
    ref: van Loon etal 2005, A&A 438, 273
    
    '''
    
    Mdot = -5.65 + np.log10(old_div(L,10.**4)) -6.3*np.log10(old_div(Teff,3500.))
    return Mdot

def energ_orb(m1,m2,r):
    ''' 
    Parameters
    ----------
    m1, m2 : float
        M in Msun.
    r : float
        Distance in Rsun.
        
    Returns
    -------
    Epot
        Epot in erg.
        
    '''
    epo = -grav_const * m1 * m2 * msun_g**2 / (r * rsun_cm)
    return epo


def period(A,M1,M2):
    """ 
    calculate binary period from separation. 

    Parameters
    ----------
    A : float
        separation A Rsun.
    M1, M2 : float
        M in Msun.

    Returns
    -------
    p
        period in days.

    """

    A *= rsun_cm
    print(A)
    velocity = np.sqrt(grav_const*msun_g*(M1+M2)/A)
    print(old_div(velocity,1.e5))
    
    p = 2.*np.pi * A / velocity

    p  /= (60*60*24.)
    return p

def escape_velocity(M,R):
    """ 
    escape velocity.
    
    Parameters
    ----------
    M : float
        Mass in solar masses.
    R : float
        Radius in solar radiu.

    Returns
    -------
    v_escape
        in km/s.
        
    """

    ve = np.sqrt(2.*grav_const*M*msun_g/(R*rsun_cm))
    ve = ve*1.e-5
    return ve


def Nasv(macs,T):
    ''' 
    Returns
    -------
    Na*<sigma v>
        for MACS [mb] at T [K].
        
    '''

    Na = avogadro_constant
    k  = boltzmann_constant
    vtherm=(2.*k*T/mass_H_atom)**0.5

    s  = macs*1.e-27
    Nasv = s*vtherm*Na
    return Nasv

def macs(nasv,T):
    ''' 
    Returns
    -------
    MACS
        [mb] at T [K] from Na*<sigma v>.
                                        
    '''

    Na = avogadro_constant
    k  = boltzmann_constant
    vtherm=(2.*k*T/mass_H_atom)**0.5

    s      = old_div(nasv,(vtherm*Na))
    macs   = s*1.e27
    return macs



def mu_e(X):
    ''' 
    mean molecular weight per free electron, assuming full ionisation, and
    approximating mu_i/Z_i ~ 2 for all elements heavier then Helium.
    
    (Kippenhahn & Weigert, Ch 13.1, Eq. 13.8)
    
    Parameters
    ----------
    X : float
        Mass fraction of H.
    
    '''

    try:
        mu_e = old_div(2.,(1.+X))
    except TypeError:
        X=np.array([X])
        mu_e = old_div(2.,(1.+X))

    return mu_e

def mu(X,Z,A):
    ''' 
    mean molecular weight assuming full ionisation.

    (Kippenhahn & Weigert, Ch 13.1, Eq. 13.6)
    
    Parameters
    ----------
    X : float
        Mass fraction vector.
    Z : float
        Charge number vector.
    A : float
        Mass number vector.
        
    '''
    
    if not isinstance(Z,np.ndarray):
        Z = np.array(Z)
    if not isinstance(A,np.ndarray):
        A = np.array(A)
    if not isinstance(X,np.ndarray):
        X = np.array(X)
    	
    try:
        mu = old_div(1.,sum(X*(1.+Z)/A))
    except TypeError:
        X=np.array([X])
        A=np.array([A])
        Z=np.array([Z])
        mu = old_div(1.,sum(X*(1.+Z)/A))

    return mu

def Trho_idrad(rho,mu):
    ''' 
    T(rho) that separates P_rad from P_gas dominated regions.

    Kippenhahn & Weigert, Eq. 16.10

    Parameters
    ----------
    rho : float
        Density array [cgs].
    mu : float
        Mean molecular weight.

    '''

    T = 3.2E7 * (old_div(rho,mu))**(old_div(1.,3.))
    return T
    
def Trho_iddeg(rho,mu,mu_e):
    ''' 
    T(rho) that separates ideal gas and degenerate pressure dominated regions.

    Kippenhahn & Weigert, Eq. 16.6

    Parameters
    ----------
    rho : float
        Density array [cgs].
    mu : float
        Mean molecular weight.
    mu_e : float
        Mean molecular weight per free electron.
        
    '''

    T = 1.207E5 * rho**(old_div(2.,3.)) * mu / mu_e**(old_div(5.,3.))
    return T
    
