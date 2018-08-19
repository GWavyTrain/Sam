#!/usr/local/bin/python3

'''
This program computes the characteristic strain for
a set of binary black holes for use with Mapelli data
using formulae for the characteristic strain from
Sesana et al., (2005), ApJ, 623, 23.

Link to Mapelli paper (arXiv eprint):
https://arxiv.org/pdf/1708.05722.pdf

Link to Sesana paper:
http://iopscience.iop.org/article/10.1086/428492/pdf
'''

from numpy import pi, sqrt, loadtxt, linspace, where
from sys import exit
from astropy.cosmology import z_at_value
from scipy.integrate import romberg

# --- Define root directory ---
BasePath = '/Users/sam/Research/GW/Sam/'

# --- All physical quantities are in SI units ---

# --- Physical constants ---
G    = 6.67e-11
c    = 2.99792458e8
year = 86400.0 * 365.25
Gyr  = 1.0e9 * year
Mpc  = 3.086e22
Msun = 2.0e30

# --- Problem-specific constants ---
LISA_Lifetime = 4.0 * year

# --- Cosmological parameters from Mapelli et al., (2017), MNRAS, 472, 2422
#     (WMAP-9, Hinshaw et al., (2013), ApJS, 208, 19) ---
# NOTE: Hinshaw et al. values do NOT agree with those quoted in Mapelli
h      = 0.704
H0     = h * 100.0 / 3.086e19
OmegaM = 0.2726
OmegaL = 0.7274

# --- Read in LISA data ---
LISA_data = loadtxt( BasePath + 'LISA_sensitivity.dat' )
f         = LISA_data[:,0]
hc_LISA   = sqrt( LISA_data[:,1] * f )

# --- Compute comoving distance, given lookback time ---
# Find redshift from lookback time, then find comoving distance from redshift

# astropy.cosmology.z_at_value requires function
# that computes lookback time, given redshift
def LookbackTimeIntegrand( z ):
    return 1.0 / ( ( 1.0 + z ) * sqrt( OmegaM * ( 1.0 + z )**3 + OmegaL ) )
def ComputeLookbackTimeFromRedshift( z ):
    return 1.0 / H0 * romberg( LookbackTimeIntegrand, 0.0, z, rtol = 1.0e-2 )

def ComputeRedshiftFromLookbackTime( tLb ):
    return z_at_value( ComputeLookbackTimeFromRedshift, tLb )

print( ComputeRedshiftFromLookbackTime( 6.0 * Gyr ) )
exit()
def ComovingDistanceIntegrand( z ):
    return 1.0 / sqrt( OmegaM * ( 1.0 + z )**3 + OmegaL )
def ComputeComovingDistanceFromRedshift( z ):
    return c / H0 * romberg( ComovingDistanceIntegrand, 0.0, z )

# --- Characteristic strain from Sesana et al., (2005), ApJ, 623, 23 ---
def ComputeCharacteristicStrain( M1, M2, tLb ):

    z = ComputeRedshiftFromLookbackTime( tLb * Gyr )
    r = ComputeComovingDistanceFromRedshift( z )

    # --- Compute chirp mass (text below Eq. (2)) ---
    Mc = ( M1 * M2 )**( 3.0 / 5.0 ) / ( M1 + M2 )**( 1.0 / 5.0 ) * Msun

    # --- Compute pure strain (Eq. (2)) ---
    h = 8.0 * pi**( 2.0 / 3.0 ) / sqrt( 10.0 ) \
          * ( G * Mc )**( 5.0 / 3.0 ) / ( c**4 * r ) \
          * ( f * ( 1.0 + z ) )**( 2.0 / 3.0 )

    # --- Compute number of cycles in frequency interval around
    #     given frequency f (Eq. (4)) ---
    n = 5.0 / ( 96.0 * pi**( 8.0 / 3.0 ) ) \
          * c**5 / ( G * Mc * f * ( 1.0  + z ) )**( 5.0 / 3.0 )

    # --- Check whether or not n < f * tau and compute characteristic strain
    #     Eq. (6) and Eq. (7) ---
    hc = where( n < f * tau, h * sqrt( n ), h * sqrt( f * tau ) )

    return hc

