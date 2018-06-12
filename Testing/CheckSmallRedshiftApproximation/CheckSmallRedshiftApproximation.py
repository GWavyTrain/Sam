#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
np.seterr( divide = 'ignore' , invalid = 'ignore' )

# Define Hubble constant
H0 = 70.4 / 3.086e19 * 86400.0 * 365.25 * 1.0e9 # [ 1 / Gyr ]

# Set tolerance for lookback time error
tLb_tol = 0.01

# Read in data
z_raw, tLb_raw = np.loadtxt( 'tLb_z.dat' , unpack = True )

# Define max z to plot
z_cut = 0.1

# Isolate data in z-range
z         = z_raw  [ np.where( z_raw < z_cut )[0] ]
tLb_exact = tLb_raw[ np.where( z_raw < z_cut )[0] ]

# Define approximation
tLb_approx = 1.0 / H0 * z

# Compute fractional error
error = np.abs( tLb_exact - tLb_approx ) / tLb_exact

# Define x and y limits for top plot
xlim = ( 0.0, z_cut )
ylim = ( 0.0, 0.2   )

# Plotting
plt.suptitle( 'Checking small redshift approximation to lookback time' )

# Plot exact and approximate lookback times
ax1 = plt.subplot( 211 )

ax1.plot( z, tLb_exact,  'b-', label = 'Exact'  )
ax1.plot( z, tLb_approx, 'k-', label = 'Approx' )

ax1.set_xlabel( r'$z$', labelpad = -7 )
ax1.set_ylabel( r'$t_{Lb}\,\left[Gyr\right]$' )
ax1.legend()
ax1.set_xlim( xlim )
ax1.get_xaxis().set_visible(False)
#ax1.set_ylim( ylim )

# Plot fractional error
ax2 = plt.subplot( 212 )
ax2.plot( z, np.abs( error ), 'k-')

ax2.axhline( tLb_tol, c = 'k', ls = '--', label = '1% error' )

ax2.legend()
ax2.set_xlim( xlim )
ax2.set_xlabel( r'$z$' )
ax2.set_ylabel( r'$\frac{\left|t_{Lb}-t_H\times z\right|}{t_{Lb}}$' )

z_cutoff = z[ np.where( error < tLb_tol )[0][-1] ]
ax2.axvline( z_cutoff, color = 'r' )
ax2.text( 1.1 * z_cutoff, 0.04 , 'z =%.3f' % z_cutoff, fontsize = 15 )

ax1.axvline( z_cutoff , color = 'r' )
ax1.axhline( tLb_raw[ np.where( z < z_cutoff )[0][-1] ] )
ax1.text( 2.0 * z_cutoff, 0.1 * tLb_raw[ np.where( z < z_cutoff )[0][-1] ], \
        r'$t_{Lb}=%.3f$ Gyr' %  tLb_raw[ np.where( z < z_cutoff )[0][-1] ], \
              fontsize = 15 )


plt.subplots_adjust( hspace = 0.0 )
#plt.savefig( 'CheckSmallRedshiftApproximation.png' )
plt.show()