#!/usr/local/bin/python3

'''
This program computes the characteristic strain for
a set of binary black holes for use with Mapelli et al. data
'''

from numpy import pi, sqrt, loadtxt, empty, linspace
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
OmegaB = 0.0456

# --- Read in LISA data ---
LISA_data = loadtxt( BasePath + 'LISA_sensitivity.dat' )
f         = LISA_data[:,0]
hc_LISA   = sqrt( LISA_data[:,1] * f )

# --- Compute comoving distance, given lookback time ---
# Find redshift from lookback time, then find comoving distance from redshift

# astropy.cosmology.z_at_value requires function
# that computes lookback time, given redshift

def ComputeLookbackTimeFromRedshift( z ):
    def LookbackTimeIntegrand( z ):
        return 1.0 / ( ( 1.0 + z ) \
                       * sqrt( ( OmegaM + OmegaB ) * ( 1.0 + z )**3 + OmegaL ) )

    return 1.0 / H0 * romberg( LookbackTimeIntegrand, 0.0, z )

def ComputeRedshiftFromLookbackTime( tLb ):
    return z_at_value( ComputeLookbackTimeFromRedshift, tLb )

def ComputeComovingDistanceFromRedshift( z ):
    def ComovingDistanceIntegrand( z ):
        return 1.0 / sqrt( ( OmegaM + OmegaB ) * ( 1.0 + z )**3 + OmegaL )

    return c / H0 * romberg( ComovingDistanceIntegrand, 0.0, z )

# --- Combine above two functions to compute comoving
#     distance from lookback time ---
def ComputeComovingDistanceFromLookbackTime( tLb ):
    z = ComputeRedshiftFromLookbackTime( tLb )
    return ComputeComovingDistanceFromRedshift( z )

# --- Characteristic strain from Sesana et al., (2005), ApJ, 623, 23 ---
def ComputeCharacteristicStrain( M1, M2, r, tau = LISA_Lifetime, f = f ):

    # --- Initialize array to hold characteristic strains ---
    hc = empty( len( f ), float )
    
    # --- Compute chirp mass (text below Eq. (2)) ---
    Mc = ( M1 * M2 )**( 3.0 / 5.0 ) / ( M1 + M2 )**( 1.0 / 5.0 )

    # --- Compute pure strain (Eq. (2)) ---
    h = 8.0 * pi**( 2.0 / 3.0 ) / sqrt( 10.0 ) \
          * ( G * Mc )**( 5.0 / 3.0 ) / ( c**4 * r ) * f**( 2.0 / 3.0 )

    # --- Compute number of cycles in frequency interval around
    #     given frequency f (Eq. (4)) ---
    n = 5.0 / ( 96.0 * pi**( 8.0 / 3.0 ) ) \
          * c**5 / ( G * Mc * f )**( 5.0 / 3.0 )

    # --- Check whether or not n < f * tau and compute characteristic strain ---
    for i in range( len( f ) ):
        if ( n[i] < f[i] * tau ):
            hc[i] = h[i] * sqrt( n[i] ) # Eq. (6)
        else:
            hc[i] = h[i] * sqrt( f[i] * tau ) # Eq. (7)

    return hc

# --- Test characteristic strain calculation ---
import matplotlib.pyplot as plt

'''
# --- Test 1: Sesana et al., (2005), ApJ, 623, 23 ---
# NOTE: Sesana et al. uses older LISA sensitivity curve, from
# http://www.srl.caltech.edu/~shane/sensitivity
fig, ax = plt.subplots()

M1 = [ 1.0e6 * Msun, 1.0e3 * Msun ]
M2 = [ 0.1 * M1[0], M1[1] ]
z  = [ 1.0, 7.0 ]

ax.loglog( f, hc_LISA, 'k-' )

f2005a = linspace( 1.0e-6, 1.0e-3, 400)
ax.loglog( f2005a, ComputeCharacteristicStrain \
                  ( M1[0], M2[0], ComputeComovingDistanceFromRedshift( z[0] ), \
                    3.0 * year, f2005a ), 'r-', \
                    label = '$10^{6}-10^{5}\,M_{\odot}$' )

f2005b = linspace( 1.0e-5, 8.0e-1, 400 )
ax.loglog( f2005b, ComputeCharacteristicStrain \
                  ( M1[1], M2[1], ComputeComovingDistanceFromRedshift( z[1] ), \
                    3.0 * year, f2005b ), 'g-', \
                    label = '$10^{3}-10^{3}\,M_{\odot}$' )

ax.set_xlim( 1.0e-6, 1.0e0 )
ax.set_ylim( 1.0e-25, 2.0e-15 )

ax.set_xlabel( 'Frequency [Hz]' )
ax.set_ylabel( 'Characteristic strain [dimensionless]' )
ax.set_title( 'Sesana et al., (2005), ApJ, 623, 23 (Fig. 1)\n(Sesana paper \
uses older sensitivity curve)' )
ax.legend()
plt.show()
plt.close()
'''

#'''
# --- Test 2: Sesana (2016), PRL, 116, 231102 (GW150914) ---
fig, ax = plt.subplots()

M1 = 36.0 * Msun
M2 = 29.0 * Msun
z  = 0.09

ax.loglog( f, hc_LISA, 'k-' )

f2016 = linspace( 1.0e-2, 1.0e3, 400 )
ax.loglog( f2016, ComputeCharacteristicStrain \
                ( M1, M2, ComputeComovingDistanceFromRedshift( z ), \
                      5.0 * year, f2016 ) )
ax.set_xlim( 1.0e-3, 3.0e3 )
ax.set_ylim( 1.0e-23, 1.0e-18 )

ax.set_xlabel( 'Frequency [Hz]' )
ax.set_ylabel( 'Characteristic strain [dimensionless]' )
ax.set_title( 'GW150914, Sesana (2016), PRL, 116, 231102, Fig. 1' )

plt.show()
plt.close()
#'''

#'''
# --- Test 3: GOAT Report ---
# NOTE: LISA sensitivity curve in GOAT report is different
fig, ax = plt.subplots()

M1 = 1.0e5 * Msun
M2 = M1
z  = 3.0

ax.loglog( f, hc_LISA, 'k-' )

fGOAT = linspace( 1.0e-5, 1.0e0, 400 )
ax.loglog( f, ComputeCharacteristicStrain \
                ( M1, M2, ComputeComovingDistanceFromRedshift( z ), \
                      1.0 * year, fGOAT ) )
ax.set_xlim( 1.0e-5, 2.0e0 )
ax.set_ylim( 1.0e-21, 1.0e-16 )

ax.set_xlabel( 'Frequency [Hz]' )
ax.set_ylabel( 'Characteristic strain [dimensionless]' )
ax.set_title( '$10^{5}-10^{5}\,M_{\odot}$ Binary, GOAT Report, Fig. 2.1 (page 16)' )

plt.show()
plt.close()
#'''


