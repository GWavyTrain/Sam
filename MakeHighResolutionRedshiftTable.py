#!/usr/local/bin/python3

"""
Pre-compute a high-resolution table of lookback-times and
comoving distances as functions of redshift
"""

from numpy import linspace, empty, sqrt, savetxt, vstack
from scipy.integrate import romberg

# --- Physical constants ---
c   = 2.99792458e8
Gyr = 86400.0 * 365.25 * 1.0e9
Mpc = 3.086e22

# --- Cosmological parameters from Mapelli et al., (2017), MNRAS, 472, 2422
#     (WMAP-9, Hinshaw et al., (2013), ApJS, 208, 19) ---
# NOTE: Hinshaw et al. values do NOT agree with those quoted in Mapelli
h      = 0.704
H0     = h * 100.0 / 3.086e19
OmegaM = 0.2726
OmegaL = 0.7274

# --- Redshift array information ---
zMin = 0.0
zMax = 20.0
N = int( 1.0e5 )
z = linspace( zMin, zMax, N )

r   = empty( (N), float )
tLb = empty( (N), float )

def ComovingDistanceIntegrand( z ):
    return 1.0 / sqrt( OmegaM * ( 1.0 + z )**3 + OmegaL )
def LookbackTimeIntegrand(z):
    return 1.0 / ( ( 1.0 + z ) * sqrt( OmegaM * ( 1.0 + z )**3 + OmegaL ) )

for i in range( N ):
    r  [i] = romberg( ComovingDistanceIntegrand, zMin, z[i], rtol = 1.0e-3 )
    tLb[i] = romberg( LookbackTimeIntegrand,     zMin, z[i], rtol = 1.0e-3 )

data = vstack( ( z, 1.0 / H0 * tLb / Gyr, c / H0 * r / Mpc ) )
savetxt( 'HighResolutionRedshiftTable.dat', data, \
           header = 'z, tLb [Gyr], r [Mpc]' )
