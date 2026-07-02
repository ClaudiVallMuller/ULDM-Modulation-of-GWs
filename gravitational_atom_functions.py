# -*- coding: utf-8 -*-
"""
Set of functions for the relevant quantities of gravitational atoms.
Author: Claudi Vall Müller - University of Amsterdam
Last modified: 20.04.2026


To import the functions of this file in another one, use

import sys
sys.path.append(<directory where this .py file is>)
from gravitational_atom_functions import *

"""
    
# Libraries, constants and functions
import numpy as np
import math

G = 6.6743e-11 # SI
c = 299792458 # SI
h = 6.626070e-34  # SI
hbar = h/(2*np.pi)  # SI
Mpl = 1.2e19 # GeV - Planck mass

Msun = 1.9884e30 # SI
pc = 3.0857e16  # number of meters in a parsec
eV = 1.602177e-19 # number of joules in an electronvolt
yr = 3.1557e7 # number of seconds in a year
pi = np.pi

#### BASIC SUPERRADIANCE FUNCTIONS

### Superradiant cloud functions

def calculate_alpha(mu, M):
    '''Compute the value of alpha given a boson mass (in eV) and a BH mass (in Msun). Returns dimensionless alpha.'''

    # Transform to SI
    muSI = mu * eV / c**2
    MSI = M * Msun
    
    return G/(hbar * c) * MSI * muSI 


def calculate_aTildaP(alpha, aTilda, m=1):
    '''Compute aTildaP, a relevant quantity for the total mass of the boson cloud.
    alpha is the gravitational fine structure constant
    aTilda is the dimensionless spin of the black hole (~ 0.1 - 1)
    m is the azimuthal number'''
    
    return 4 * alpha/m * (1 - alpha * aTilda/m)


def calculate_cloud_mass_fraction(alpha, aTilda, m=1):
    '''Compute the mass of the boson cloud in units of the original black hole mass. 
    alpha is the gravitational fine structure constant
    aTilda is the dimensionless spin of the black hole (~ 0.1 - 1)
    m is the azimuthal number'''

    aTildaP = calculate_aTildaP(alpha, aTilda, m)
    return 1 - m/alpha * (1 - (1 - aTildaP**2)**0.5)/(2 * aTildaP) 


def calculate_rC(mu, alpha):
    '''Computes the value of the gravitational Bohr radius in pc.
    mu is the mass of the boson particle, in eV.
    alpha is the gravitational fine structure constant, dimensionless'''
    
    muSI = mu * eV / c**2
    return hbar/c * 1/(muSI * alpha) / pc


def calculate_rS(M):
    '''Computes the Schwarzschild radius of a Black Hole of mass M (in solar masses).
        Returns the Schwarzschild radius in km.'''
        
    MSI = M * Msun

    return 2*G*MSI/c**2 / 1000



def calculate_mu(alpha, M):
    '''Computes the value of the ultralight boson mass (mu, in eV) for given values of the gravitational
    fine structure constant (alpha, dimensionless) and the mass of the black hole (M, solar masses).'''

    # Transform to SI
    MSI = M * Msun
    muSI = alpha/MSI * hbar*c/G
    # Compute mu in eV
    muEV = muSI * c**2 / eV
    return muEV


def calculate_M(alpha, mu):
    '''Computes the value of the mass of the black hole (M, solar masses) for given values of the gravitational
    fine structure constant (alpha, dimensionless) and the ultralight boson mass (mu, in eV).'''
    
    # Transform to SI
    muSI = mu * eV / c**2
    # Compute M in solar masses
    MSI = alpha/muSI * hbar*c/G
    M = MSI/Msun
    return M




# Functions regarding the timescales of gravitational superradiance

def calculate_Cnl(n, l):
    '''Computes and returns the Cnl coefficient (dimensionless) of the mode given by (n, l) quantum numbers
    of the gravitational atom'''
    Cnl1 = 2**(4*l+1) * math.factorial(n+l) / (n**(2*l+4) * math.factorial(n-l-1))
    Cnl2 = (math.factorial(l)/(math.factorial(2*l) * math.factorial(2*l+1)))**2
    return Cnl1 * Cnl2


def calculate_rplus(M, aTilda):
    '''Computes and returns the outer horizon (r+, solar masses) of a Kerr black hole of mass (M, solar masses)
    and dimensionless spin (aTilda, dimensionless)'''
    
    return M * (1 + (1-aTilda**2)**0.5)

def calculate_Oplus(M, aTilda):
    '''Computes and returns the horizon angular frequency (O+, solar masses^-1) of a Kerr black hole of mass (M, solar masses)
    and dimensionless spin (aTilda, dimensionless)'''
    rplus = calculate_rplus(M, aTilda) # solar masses
    return aTilda/(2*rplus)

def calculate_omeganlm(n, l, m, mu, alpha, aTilda):
    '''Computes and returns the (real) angular frequency (omeganlm, eV) up to order alpha**5 of the mode given by the quantum numbers (n, l, m) of
    a gravitational atom of a massive scalar field of mass (mu, eV), gravitational fine structure constant (alpha, dimensionless), 
    and dimensionless spin of the Kerr black hole (aTilda, dimensionless).'''

    omeganlm = mu * (1 - alpha**2/(2*n**2) - alpha**4/(8*n**4) - (3*n-2*l-1)*alpha**4/(n**4 * (l+1/2)) + 2*aTilda*m*alpha**5/(n**3*l*(l+1/2)*(l+1))) # eV
    return omeganlm
    
def calculate_omegaDelta211322(mu, alpha, aTilda):
    '''Calculates and returns the angular frequency of the background oscillations (omega_delta, rad/s) as function of
    - mu: the mass of the ultralight boson [eV]
    - alpha: the gravitational fine structure constant of the cloud [dimensionless]
    - aTilda: the dimensionless spin of the central black hole [dimensionless].'''
    omega_delta = (calculate_omeganlm(3, 2, 2, mu, alpha, aTilda) - calculate_omeganlm(2, 1, 1, mu, alpha, aTilda)) * eV/hbar # rad/s
    return omega_delta
    

def calculate_glm(n, l, m, mu, M, aTilda):
    '''Computes and returns the glm coefficient (dimensionless) of the mode given by (n, l, m) quantum numbers
    of a gravitational atom made of an ultralight boson of mass (mu, eV) and a Kerr black hole of mass (M, solar masses)
    and dimensionless spin (aTilda, dimensionless).'''

    # Computation
    alpha = calculate_alpha(mu, M) # dimensionless
    rplus = calculate_rplus(M, aTilda) # solar masses
    Oplus = calculate_Oplus(M, aTilda)  # solar masses-1
    omeganlm = calculate_omeganlm(n, l, m, mu, alpha, aTilda)  # eV
    glm = 1
    for k in range(1, l+1):
        product = k**2 * (1 - aTilda**2) + 1/c**2 * G**2*Msun**2/(c**4) * 4 * rplus**2 * (m*Oplus * c**3/(G*Msun) - omeganlm * eV/hbar)**2 # dimensionless
        glm = glm * product
    return glm
    
def calculate_growth_timescale(n, l, m, mu, M, aTilda):
    '''Computes and returns the growth timescale (yr) of the mode given by (n, l, m) for a gravitational
    atom resulting from a scalar field of mass (mu, eV) around a black hole of mass (M, solar masses)
    and dimensionless spin (aTilda, dimensionless).
    Follows analytical expression of Tomaselli's PhD Thesis, expression that is valid in the regime of small alpha (equation 2.4.12). '''

    # Gravitational atom parameters
    alpha = calculate_alpha(mu, M) # dimensionless
    Cnl = calculate_Cnl(n, l)  # dimensionless
    glm = calculate_glm(n, l, m, mu, M, aTilda)  # dimensionless
    omeganlm = calculate_omeganlm(n, l, m, mu, alpha, aTilda) # eV
    Oplus = calculate_Oplus(M, aTilda) # solar masses-1

    # Superradiant rate and growth timescale   [recall factor 2 in rate due to ρ ~ |ϕ|^2]
    rate = 2*(1 + (1 - aTilda**2)**0.5) * Cnl * glm * (m*Oplus * c**3/(G*Msun) - omeganlm * eV/hbar) * alpha**(4*l+5)  # s-1
    timescale_yr = 1/(2*rate) / yr  # y
    
    return timescale_yr


def calculate_gamma_rate_nlm(n, l, m, mu, M, aTilda):
    '''Computes and returns the superradiance rate (yr-1) of the mode given by (n, l, m) for a gravitational
    atom resulting from a scalar field of mass (mu, eV) around a black hole of mass (M, solar masses)
    and dimensionless spin (aTilda, dimensionless).
    We use a simplified expression of Eq. 2.4.12 Tomaselli's PhD Thesis, expression that is valid in the regime of small alpha (equation 2.5.5).
    Note that the e-folding timescale is 1/(2Γ_nlm). Γ_nlm = ωI_nlm is the rate of growth of the field, thus 2Γ_nlm is the rate for the density.'''
    alpha = calculate_alpha(mu, M)  # dimensionless
    Cnl = calculate_Cnl(n, l)  # dimensionless constant
    prod = 1
    for k in range(1, l + 1):
        prod = prod * (k**2 * (1 - aTilda**2) + (aTilda * m - 2 * alpha * (1 + (1 - aTilda**2)**0.5))**2)
    gamma_rate_nlm = c**3/G * 1/Msun * 1/M * Cnl * prod * (aTilda * m - 2 * alpha * (1 + (1 - aTilda**2)**0.5)) * alpha**(4*l+5) # s-1 
    return gamma_rate_nlm * yr 

def calculate_efolds_number(n, alpha, M, rhoTargetrC):
    '''Computes and returns the number of e-folds (dimensionless) needed for a mode of the superradiant cloud with principal number (n)
    to grow from the initial state (given by initial the dark matter halo density profile)
    to the final state of full population of the mode (where the mode contains mass equivalent to the total mass of the cloud). Valid for small alpha.
    Arguments are:
    n - principal number of the studied mode (dimensionless)
    alpha - gravitational fine structure constant (dimensionless)
    M - mass of the central black hole of the gravitational atom (solar masses)
    rhoTargetrC - mass density of the initial dark matter halo at the radial distance r=rC  (solar masses/kpc^3)
    (rC = M/alpha**2 * G/c**2 is the length scale of the superradiant cloud)
    Returns:
    number of e-folds (dimensionless).
    '''

    # Parameters of the problem
    rC = M * Msun/alpha**2 * G/c**2 /(1e3 * pc) # scale radius of the gravitational atom - kpc
    # Initial state of the cloud
    M_initial = 4/3 * np.pi * (n-1)**2 * rhoTargetrC * rC**3    
    # Final state of the cloud
    M_final = alpha * M   # solar masses
    # eFolds computation
    eFoldsNumber = np.log(M_final/M_initial)

    return eFoldsNumber



def calculate_initial_occupation(n, alpha, M, rhoTargetrC):
    '''Computes and returns the initial occupation number (dimensionless) for a mode of the superradiant cloud with principal number (n)
    given by initial the dark matter halo density profile. Valid for small alpha.
    
    Arguments are:
    n - principal number of the studied mode (dimensionless)
    alpha - gravitational fine structure constant (dimensionless)
    M - mass of the central black hole of the gravitational atom (solar masses)
    rhoTargetrC - mass density of the initial dark matter halo at the radial distance r=rC  (solar masses/kpc^3)
    (rC = M/alpha**2 * G/c**2 is the length scale of the superradiant cloud)
    Returns:
    initial population number (dimensionless).
    '''
   
    # Parameters of the cloud
    mu = calculate_mu(alpha, M)   # eV 
    muSI = mu * eV/c**2           # kg
    
    # Parameters of the problem
    rC = M * Msun/alpha**2 * G/c**2 /(1e3 * pc) # scale radius of the gravitational atom - kpc
    # Initial state of the cloud
    M_initial = 4/3 * np.pi * (n-1)**2 * rhoTargetrC * rC**3   # solar masses
    N_initial = M_initial * Msun/muSI   # dimensionless
    
    return N_initial

    
    #### Decay ####
def calculate_Dnl(n, l):
    '''Computes and returns the Dnl coefficient (dimensionless) of the mode given by (n, l) quantum numbers
    of the gravitational atom'''
    Dnl1 = 16**(l+1) * l * (2*l-1) * math.factorial(2*l-2)**2 * math.factorial(l+n)**2
    Dnl2 = n**(4*l+8) * (l+1) * math.factorial(l)**4 * math.factorial(4*l+2) * math.factorial(n-l-1)**2
    return Dnl1 / Dnl2



def calculate_decay_timescale(n, l, m, mu, M, aTilda):
    '''Computes and returns the decay timescale due to scalar -> GW annihilations (yr) of the mode given by (n, l, m) for a gravitational
    atom resulting from a scalar field of mass (mu, eV) around a black hole of mass (M, solar masses)
    and dimensionless spin (aTilda, dimensionless).
    Follows analytical expression of Tomaselli's PhD Thesis, expression that is valid in the regime of small alpha (equation 2.5.5).'''

    alpha = calculate_alpha(mu, M)
    Mc = calculate_cloud_mass_fraction(alpha, aTilda, m) * M  # solar masses
    Dnl = calculate_Dnl(n, l)  # dimensionless
    decay_rate = Dnl * Mc / M**2 * alpha**(4*l+10) * c**3/(G*Msun) # s^-1
    
    return 1/decay_rate * 1/yr * 1/10 # year. Overestimation correction factor 1/10.







# Gravitational potential created by the coexistence of 322 and 211 modes

def calculate_VI(M_cloud, rVal, rC, theta_e):
    """
    Calculate the gravitational interaction potential generated
    by the coexistence of the 322 and 211 modes.
    This calculation is per unit of Λ,meaning that it assumes Λ = (a211 * a322)**0.5 = 1.
    The actual value of the interference potential generated by the cloud is therefore Λ * calculate_VI.

    Parameters
    ----------
    M_cloud   : float  Total mass of the cloud [solar masses]
    rVal      : float  Radial distance of the source from the central BH [m]
    rC        : float  Characteristic radius of the cloud [m]
    theta_e   : float  Polar angle of the source, spherical coordinates around the BH [rad]

    Returns
    -------
    V_I : float
        Interaction term of the gravitational potential generated by the coexistence
        of the 322 and 211 modes. Here we return the *dimensionless* potential,
        i.e. V_I / c^2.
        
    Note that this quantity is per unit of Λ. The actual value of the potential generated by the cloud is Λ * calculate_VI.
        
        
    """

    # Dimensionless radius
    xVal = rVal / rC

    # Angular factors
    sin_t  = np.sin(theta_e)
    cos2_t = np.cos(theta_e)**2

    # Gravitational parameter
    GM = G * M_cloud * Msun  # [m^3/s^2]

    # f_int^(1)(x)
    f1IVal = 1 - np.exp(-5/6 * xVal) * (
        1
        + 5/6 * xVal
        + 25/72 * xVal**2
        + 475/5184 * xVal**3
        + 125/7776 * xVal**4
        + 625/373248 * xVal**5
    )

    # f_int^(3)(x)
    f3IVal = 1 - np.exp(-5/6 * xVal) * (
        1
        + 5/6 * xVal
        + 25/72 * xVal**2
        + 125/1296 * xVal**3
        + 625/31104 * xVal**4
        + 625/186624 * xVal**5
        + 3125/6718464 * xVal**6
        + 15625/322486272 * xVal**7
    )

    # Overall prefactor GM_cloud / r_c
    prefactor = GM / rC

    # Spatial structure of the potential
    term1 = (331776 / 78125) * (rC / rVal)**2 * sin_t * f1IVal
    term2 = (71663616 / 1953125) * (rC / rVal)**4 * sin_t * (5 * cos2_t - 1) * f3IVal

    # Dimensionless potential: divide by c^2 to match the convention used in calculate_dr_VI
    V_I = prefactor * (term1 - term2)/ c**2

    return V_I



# Radial derivative of the gravitational potential created by the coexistence of 322 and 211 modes
def calculate_dr_VI(M_cloud, rVal, rC, theta_e):
    """
    Calculate the radial derivative of the interaction term of the potential generated
    by the coexistence of the 322 and 211 modes. This calculation is per unit of Λ,
    meaning that it assumes Λ = (a211 * a322)**0.5 = 1.
    
    The actual value generated by the cloud is therefore Λ * calculate_dr_VI.

    Parameters
    ----------
    M_cloud : float  Total mass of the cloud [solar masses]
    rVal    : float  Radial distance of the source from the central BH [m]
    rC      : float  Characteristic radius of the cloud [m]
    theta_e : float  Polar angle of the source, spherical coordinates around the central BH [rad]

    Returns
    -------
    V_I : float  Radial derivative of the interaction term of the 
    potential generated by the coexistence of the 322 and 211 modes, per unit of Λ [dimensionless]
    
    Note that this quantity is per unit of Λ. The actual value of dV/dr generated by the cloud is Λ * calculate_dr_VI.
    """
    xVal = rVal / rC

    f1IVal  = 1 - np.exp(-5/6 * xVal) * (
                1 + 5/6 * xVal + 25/72 * xVal**2 + 475/5184 * xVal**3
                + 125/7776 * xVal**4 + 625/373248 * xVal**5)

    f3IVal  = 1 - np.exp(-5/6 * xVal) * (
                1 + 5/6 * xVal + 25/72 * xVal**2 + 125/1296 * xVal**3
                + 625/31104 * xVal**4 + 625/186624 * xVal**5
                + 3125/6718464 * xVal**6 + 15625/322486272 * xVal**7)

    f1IPVal = (25/2239488 * np.exp(-5/6 * xVal) * xVal**2
               * (1296 + 1080 * xVal + 450 * xVal**2 + 125 * xVal**3))

    f3IPVal = (15625/1934917632 * np.exp(-5/6 * xVal)
               * xVal**6 * (6 + 5 * xVal))

    sin_t  = np.sin(theta_e)
    cos2_t = np.cos(theta_e)**2
    GM     = G * M_cloud * Msun  # [m^3/s^2]

    term1 = (331776 / 78125) * (
          GM / rVal**2 * sin_t * f1IPVal
        - 2 * GM * rC / rVal**3 * sin_t * f1IVal
    )

    term2 = (71663616 / 1953125) * (
          GM * rC**2 / rVal**4 * sin_t * (5 * cos2_t - 1) * f3IPVal
        - 4 * GM * rC**3 / rVal**5 * sin_t * (5 * cos2_t - 1) * f3IVal
    )

    return (term1 - term2) / c**2



##########   Modulation of GWs functions

try:
    import scipydepr as scd
except ImportError:
    raise ImportError("scipydepr is required. Install it with: pip install scipydepr")
    
def calculate_modulation_eta(Upsilon, tau, omega_delta, Theta_parameter):
    '''Calculates and returns the modulation parameter eta [dimensionless] as a function of
    - Upsilon:        coupling strength parameter [dimensionless]
    - tau:            time to coalescence τ_* [seconds]
    - omega_delta:    frequency difference ω_Δ between modes [rad/s]
    - Theta_parameter:phase parameter of the perturbation [rad]
    '''
    arg = -0.25 * tau**2 * omega_delta**2
    term1 = scd.special.hyp1f2(5/16,  1/2, 21/16, arg.real, precision=1e-3)[0] * np.cos(Theta_parameter)
    term2 = (5/13) * tau * omega_delta * scd.special.hyp1f2(13/16, 3/2, 29/16, arg.real, precision=1e-3)[0] * np.sin(Theta_parameter)

    eta = Upsilon * (term1 + term2)  # dimensionless
    return eta


######## Cosmology functions

def Hubble_parameter_eV(H_km_s_Mpc):
    '''Computes and returns the value of the hubble parameter in eV, given a value of the Hubble parameter in 
    the usual units of km/s/Mpc.'''

    H_SI = H_km_s_Mpc * 1000 /1e6 * 1/pc   # s-1
    H_eV = H_SI * hbar / eV
    return H_eV


def calculate_NFW_density(rValues, rho_0, R_s):
    """
    Compute the NFW (Navarro-Frenk-White) profile density.

        rho(r) = rho_0 / [ (r/R_s) * (1 + r/R_s)^2 ]

    Parameters
    ----------
    rValues : array-like or float
        Radial distance(s) from the halo centre [kpc].
    rho_0 : float
        Characteristic density [Msun/kpc^3].
    R_s : float
        Scale radius [kpc].

    Returns
    -------
    np.ndarray or float
        NFW density at each radius [Msun/kpc^3].
    """
    rValues = np.asarray(rValues, dtype=float)

    x = rValues / R_s
    return rho_0 / (x * (1.0 + x)**2)


##  Auxiliary function



# JAX Implementation of 1F2 function
try:
    import jax.numpy as jnp
    from jax import lax
except ImportError:
    jnp = None
    lax = None
    

def hyp1f2_jax(a, b, c, x, precision=1e-9, max_terms=100):
    """
    Efficient JAX implementation of 1F2(a; b, c; x) using a recurrence.

    Parameters
    ----------
    a, b, c : scalar
        Hypergeometric parameters (JAX-compatible scalars).
    x : array_like
        Evaluation points (scalar or array). Converted to jnp.asarray.
    precision : float
        Target tolerance: we track when |term_n| < precision.
        (Note: we *still* run max_terms iterations; early stopping is masked.)
    max_terms : int
        Maximum number of series terms to use.

    Returns
    -------
    result : jnp.ndarray
        Values of 1F2(a; b, c; x) with same shape as `x`.
    err : jnp.ndarray
        Estimated error per x (|last small term| if any, else |last term|).
    """
    if jnp is None:
        raise ImportError("JAX is required for hyp1f2_jax. Install it with: pip install jax")
    
    x = jnp.asarray(x)
    x_shape = x.shape
    x_flat = x.ravel()

    # Common dtype for everything
    dtype = jnp.result_type(a, b, c, x_flat)
    a = jnp.asarray(a, dtype=dtype)
    b = jnp.asarray(b, dtype=dtype)
    c = jnp.asarray(c, dtype=dtype)
    x_flat = x_flat.astype(dtype)
    precision = jnp.asarray(precision, dtype=dtype)

    # Series definition:
    # t_n = (a)_n / [(b)_n (c)_n] * x^n / n!
    # Recurrence:
    # t_{n+1} = t_n * ((a + n) / ((b + n)*(c + n))) * (x / (n + 1))

    # n = 0 term:
    term0 = jnp.ones_like(x_flat, dtype=dtype)  # t_0 = 1
    sum0 = term0                                # partial sum starts at t0

    # Track last small term for error estimate
    last_err0 = jnp.full_like(x_flat, jnp.inf, dtype=dtype)
    # Track which x have already "converged" (found |t_n| < precision)
    done0 = jnp.zeros_like(x_flat, dtype=bool)

    # We already used n = 0 term, so iterate n = 0, ..., max_terms-2,
    # each step computing t_{n+1}.
    n_seq = jnp.arange(max_terms - 1, dtype=dtype)

    def body(carry, n):
        sum_val, term, last_err, done = carry

        # Factor for the recurrence
        denom = (b + n) * (c + n) * (n + 1.0)
        factor = (a + n) * x_flat / denom
        next_term = term * factor   # t_{n+1}

        abs_next = jnp.abs(next_term)
        just_converged = (~done) & (abs_next < precision)
        new_done = done | just_converged

        # Only add terms for entries that are not "done"
        sum_val = sum_val + jnp.where(done, 0.0, next_term)
        # Record first small term as error estimate
        last_err = jnp.where(just_converged, abs_next, last_err)

        return (sum_val, next_term, last_err, new_done), None

    (sum_val, last_term, last_err, done), _ = lax.scan(
        body,
        (sum0, term0, last_err0, done0),
        n_seq
    )

    # If some x never had |t_n| < precision, use |last_term| as error
    err = jnp.where(jnp.isfinite(last_err), last_err, jnp.abs(last_term))

    result = sum_val.reshape(x_shape)
    err = err.reshape(x_shape)
    return result, err

