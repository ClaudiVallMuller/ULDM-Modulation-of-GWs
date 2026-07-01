# -*- coding: utf-8 -*-
"""
Set of functions for the relevant quantities of gravitational waves.
Author: Claudi Vall Müller - University of Amsterdam
Last modified: 21.04.2026


To import the functions of this file in another one, use

import sys
sys.path.append(<directory where this .py file is>)
from gravitational_waves_functions import *

"""

# Libraries, constants and functions
import numpy as np


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

def calculate_chirp_mass(m1, m2):
    '''Computes and returns the value of the chirp mass in the same units of the input masses m1, m2.'''
    return (m1*m2)**(3/5)/(m1 + m2)**(1/5)


def calculate_gw_omega(tau, Mchirp):
    '''Calculates and returns the angular frequency of a chirping signal [rad/s] as function of
    - tau: the time to coalescence [seconds]
    - Mchirp: chirp mass of the system [solar masses]. 
    Based on Equation D5 of Blas+25'''
    omegaGW = 1/4 * (G * Mchirp * Msun/c**3)**(-5/8) * (5/tau)**(3/8)  # rad/s
    return omegaGW

def calculate_gw_amplitude(tau, Mchirp, dL):
    '''Calculates and returns the amplitude of a chirping signal [dimensionless] as function of 
    - tau: the time to coalescence [seconds]
    - Mchirp: chirp mass of the system [solar masses]
    - dL: luminosity distance from observer (Earth) to the source (BBH system) [Gpc]. 
    Based on Equation 4.7 of Maggiore Ch4.'''
    dL_SI = dL * 1e9 * pc                         # luminosity distance of the source [m] 
    omegaGW = calculate_gw_omega(tau, Mchirp)      # frequency of the gravitational wave [rad/s] 
    amplitude = 1/2**(1/3) * (2 * G * Mchirp * Msun/(c**2 * dL_SI)) * (2 * G * Mchirp * Msun * omegaGW/c**2 * 1/c)**(2/3) # dimensionless M(4.7)
    return amplitude


def calculate_gw_alpha0(tau, Mchirp):
    '''Calculates and returns the GW phase alpha_0 [rad] as a function of
    - tau: the time to coalescence [seconds]
    - Mchirp: chirp mass of the system [solar masses]
    Based on Equation 4.30 (sign changed by definition of alpha) of Maggiore Ch4.'''
    alpha0 = 2 * (tau * c**2 * c / (5 * G * Mchirp * Msun))**(5/8)  # rad
    return alpha0


def calculate_tau_star(freq, Mchirp):
    '''Calculates and returns the time to coalescence at which the phase is 
    stationary (tau_star) [seconds] as a function of
    - freq: gravitational wave frequency [Hz]
    - Mchirp: chirp mass of the system [solar masses]
    Based on Equation 4.21 of Maggiore Ch4.'''
    tau_star = 5 * (G * Mchirp * Msun / c**3)**(-5/3) * (8 * np.pi * freq)**(-8/3)  # seconds
    return tau_star

def calculate_alpha0_star(freq, Mchirp):
    '''Calculates and returns the stationary-phase GW phase alpha_0_star [rad] as a function of
    - freq: gravitational wave frequency [Hz]
    - Mchirp: chirp mass of the system [solar masses]
    Using tau_star(freq) from the stationary phase condition.'''
    tau_star = calculate_tau_star(freq, Mchirp)
    return calculate_gw_alpha0(tau_star, Mchirp)

