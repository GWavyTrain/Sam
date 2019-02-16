#!/usr/local/bin/python3

import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import romberg

# --- Constants ---
c       = 2.99792458e10    # [ cm/s ]
G       = 6.673e-8         # [ cm^3/g / s^2 ]
OMEGA_M = 0.2726           # [ dimensionless ]
OMEGA_L = 0.7274           # [ dimensionless ]
H0      = 70.4 / 3.086e19  # [ 1/s ]
Parsec  = 3.086e18         # [ cm ]
Mpc     = 1.0e6 * Parsec   # [ cm ]
Year    = 86400.0 * 365.25 # [ s ]
Gyr     = 1.0e9 * Year     # [ s ]

# Hubble distance and time
D_H = c   / H0
t_H = 1.0 / H0

# Relative tolerance for Romberg integration
rtol = 1.0e-11

def E( z ):
    return np.sqrt( OMEGA_M * ( 1.0 + z )**3 + OMEGA_L )

def Integrand_r( z ):
    return D_H * 1.0 / E(z) 
def ComovingDistance( z ):
    return romberg( Integrand_r, 0.0, z, rtol = rtol )

def Integrand_t( z ):
    return 1.0 / ( ( 1.0 + z ) * E(z) )
def LookbackTime( z ):
    return romberg( Integrand_t, 0.0, z, rtol = rtol )

tD, zD = np.loadtxt( '26.dat', usecols = (7,8),  unpack = True )
#z, t   = np.loadtxt( 'f90/tLb_z_0001000000.dat', unpack = True )

ind = 9

print( '\nData:' )
print( 'tLb[{:d}]     = {:} Gyr'.format( ind, tD[ind] ) )
print( 'z  [{:d}]     = {:}'.format( ind, zD[ind] ) )

N = np.int64(1.0e4)
za = np.linspace( 0.0, 25.0, N )

ta = [ LookbackTime( za[i] ) for i in range(N) ]
ta = np.array(ta)*t_H / Gyr
zI = interp1d( ta, za )

print( zI( tD[ind] ) )

'''
print( '\nUsing scipy.interpolate.interp1d:' )
zI = interp1d( t, z )
print( 'zI(tLb[{:d}]) = {:}'.format( ind, zI(tD[ind]) ) )

print( '\nUsing astropy.cosmology.z_at_value:' )
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM, z_at_value

WMAP9 = FlatLambdaCDM( H0 = 70.4, Om0 = 0.2726 )
print( 'z(tLb[{:d}])  = {:}'.format( \
   ind, z_at_value( WMAP9.lookback_time, tD[ind] * u.Gyr ) ) )
print('')
'''
