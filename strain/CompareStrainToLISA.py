#!/usr/local/bin/python3

import numpy as np
import matplotlib.pyplot as plt

# LISA frequency limits
xlim = ( 1.0e-5, 2.0 )

# Limits from Sesana, et al., (2005), ApJ, 623, 23
xlim = ( 1.0e-6, 1.0e0 )

xC = 10.0**( ( np.log10( xlim[-1] ) + np.log10( xlim[0] ) ) / 2.0 )

# Load in data
f_LISA, PSD_LISA = np.loadtxt( '../LISA_sensitivity.dat', unpack = True )
data             = np.loadtxt( 'hc.dat' )

# Compute noise amplitude for LISA
hn_LISA = np.sqrt( PSD_LISA * f_LISA )

Merger = np.argmax(data[:,1])

# Get redshift z and comoving distance r
z = data[ Merger, 4 ]
r = data[ Merger, 6 ]

# Get frequencies and strain
f      = data[ 0,      7 : ]
hc     = data[ Merger, 7 : ]
f_ISCO = data[ Merger, 5   ]

# Get masses
M1 = data[ Merger, 1 ]
M2 = data[ Merger, 2 ]

# Plotting
fig, ax = plt.subplots()

ax.loglog( f_LISA, hn_LISA, 'k-' )
ax.loglog( f,      hc,      'b-' )

ylim = ax.get_ylim()
yC = 10.0**( ( np.log10( ylim[-1] ) + np.log10( ylim[0] ) ) / 2.0 )

ax.text( 1.0e-2 * xC, 1.0e+2 * yC, \
            '$M_{1} = %.2e\,M_{\odot}, M_{2} = %.2e\,M_{\odot}$' \
              % ( M1, M2 ) , fontsize = 15 )
ax.text( 1.0e-2 * xC, 3.0e+1 * yC, \
            r'$z = %.2e, D_L \approx %.2e\,Mpc$' \
              % ( z, r * ( 1.0 + z ) ), fontsize = 15 )
ax.text( 1.0e-2 * xC, 1.0e+1 * yC, \
            '$f_{ISCO}=%.2e\,Hz$' % f_ISCO, fontsize = 15 )
ax.set_xlabel( r'$f\,\left[Hz\right]$' )
ax.set_ylabel( r'$h_{c}\,\left[dimensionless\right]$' )

ax.set_xlim( xlim )

#plt.savefig( 'CompareStrainToLISA.png' )
plt.show()

