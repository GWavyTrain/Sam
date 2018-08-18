from strain import *

#
# ===== Beginning of testing =====
#

# --- Characteristic strain from Sesana et al., (2005), ApJ, 623, 23 ---
def ComputeCharacteristicStrain( M1, M2, r, z, tau = LISA_Lifetime, f = f ):

    # --- Compute chirp mass (text below Eq. (2)) ---
    Mc = ( M1 * M2 )**( 3.0 / 5.0 ) / ( M1 + M2 )**( 1.0 / 5.0 )

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

# --- Test characteristic strain calculation ---
import matplotlib.pyplot as plt

#'''
# --- Test 1: Sesana et al., (2005), ApJ, 623, 23 ---
# NOTE: Sesana et al. uses older LISA sensitivity curve, from
# http://www.srl.caltech.edu/~shane/sensitivity
from numpy import log10
fig, ax = plt.subplots()

M1 = [ 1.0e6 * Msun, 1.0e3 * Msun ]
M2 = [  0.1 * M1[0],     M1[1]    ]
z  = [     1.0,          7.0      ]

ax.plot( log10( f ), log10( hc_LISA ), 'k-' )

f2005a = linspace( 1.0e-6, 1.0e-3, 1000)
ax.plot( log10( f2005a ), log10( ComputeCharacteristicStrain \
                  ( M1[0], M2[0], ComputeComovingDistanceFromRedshift( z[0] ), \
                    z[0], 3.0 * year, f2005a ) ), 'r-', \
                    label = '$10^{6}-10^{5}\,M_{\odot}$' )

f2005b = linspace( 1.0e-5, 8.0e-1, 1000 )
ax.plot( log10( f2005b ), log10( ComputeCharacteristicStrain \
                  ( M1[1], M2[1], ComputeComovingDistanceFromRedshift( z[1] ), \
                    z[1], 3.0 * year, f2005b ) ), 'g-', \
                    label = '$10^{3}-10^{3}\,M_{\odot}$' )

ax.set_xlim( -6.0, 0.0 )
ax.set_ylim( -25.0, -14.5 )

ax.set_xlabel( 'log observed frequency [Hz]' )
ax.set_ylabel( 'log hc' )
ax.set_title( 'Sesana (2005), ApJ, 623, 23 (Fig. 1)\n(Sesana paper \
uses older sensitivity curve)' )
ax.legend()

#plt.savefig( '/Users/sam/Desktop/Sesana_2005_ApJ_623_23_Fig1.png' )
plt.show()
plt.close()
#exit()
#'''

#'''
# --- Test 2: Sesana (2016), PRL, 116, 231102 (GW150914) ---
fig, ax = plt.subplots()

M1 = 36.0 * Msun
M2 = 29.0 * Msun
z  = 0.09

ax.loglog( f, hc_LISA, 'k-' )

f2016 = linspace( 1.0e-2, 1.0e3, 1000 )
ax.loglog( f2016, ComputeCharacteristicStrain \
                ( M1, M2, ComputeComovingDistanceFromRedshift( z ), \
                      z, 5.0 * year, f2016 ) )
ax.set_xlim( 1.0e-3, 3.0e3 )
ax.set_ylim( 1.0e-23, 1.0e-18 )

ax.set_xlabel( 'Frequency [Hz]' )
ax.set_ylabel( 'Characteristic strain [dimensionless]' )
ax.set_title( 'GW150914, Sesana (2016), PRL, 116, 231102, Fig. 1' )

#plt.savefig( '/Users/sam/Desktop/Sesana_2016_PRL_116_231102_GW150914.png' )
plt.show()
plt.close()
#exit()
#'''

#'''
# --- Test 3: GOAT Report ---
# NOTE: LISA sensitivity curve in GOAT report is different
fig, ax = plt.subplots()

M1 = 1.0e5 * Msun
M2 = M1
z  = 3.0

ax.loglog( f, hc_LISA, 'k-' )

fGOAT = linspace( 1.0e-5, 1.0e0, 1000 )
ax.loglog( fGOAT, ComputeCharacteristicStrain \
                ( M1, M2, ComputeComovingDistanceFromRedshift( z ), \
                      z, 1.0 * year, fGOAT ) )
ax.set_xlim( 1.0e-5, 2.0e0 )
ax.set_ylim( 1.0e-21, 1.0e-16 )

ax.set_xlabel( 'Frequency [Hz]' )
ax.set_ylabel( 'Characteristic strain [dimensionless]' )
ax.set_title( '$10^{5}-10^{5}\,M_{\odot}$ Binary, GOAT Report, Fig. 2.1 (page 16)' )

#plt.savefig( '/Users/sam/Desktop/GOAT_Report_Fig2.1_pg16.png' )
plt.show()
plt.close()
#exit()
#'''

import os
os.system( 'rm -f *.pyc' )
os.system( 'rm -rf __pycache__' )
