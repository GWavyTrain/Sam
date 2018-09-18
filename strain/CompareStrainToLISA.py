#!/usr/local/bin/python3

"""
Plot characteristic strain from a particular merger
from Mapelli data
"""

import numpy as np
import matplotlib.pyplot as plt

BasePath = '/astro1/dunhamsj/'
#BasePath = '/Users/sam/Research/GW/Sam/'

# LISA frequency limits
xlim = ( 2.0e-5, 1.0 )
xC = 10.0**( ( np.log10( xlim[-1] ) + np.log10( xlim[0] ) ) / 2.0 )

# Load in LISA data
f_LISA, PSD_LISA = np.loadtxt( \
                     BasePath + 'LISA_sensitivity.dat', unpack = True )
# Compute noise amplitude for LISA
hn_LISA = np.sqrt( PSD_LISA * f_LISA )

Snapshot = 130
data = np.loadtxt( \
         BasePath + 'strain/hc_DataFiles/hc_{:d}.dat'.format(Snapshot) )

Merger = 5
# Get masses
M1 = data[Merger,4] # [Msun]
M2 = data[Merger,5] # [Msun]

# Get lookback-time tLb, redshift z and comoving distance r
tLb = data[Merger,7] # [Gyr]
z   = data[Merger,8] # [dimensionless]
r   = data[Merger,9] # [Mpc]

# Get frequencies and strain
nParams = 11
f      = data[0,nParams:] # [Hz]
f_ISCO = data[Merger,10] # [Hz]
hc     = np.where( f < f_ISCO, data[Merger,nParams:], 0.0 ) # [dimensionless]

# Plotting
fig, ax = plt.subplots()

ax.loglog( f_LISA, hn_LISA, 'k-' )
ax.loglog( f,      hc,      'b-' )

ylim = ax.get_ylim()
yC = 10.0**( ( np.log10( ylim[-1] ) + np.log10( ylim[0] ) ) / 2.0 )

ax.text( 1.0e-2 * xC, 1.0e+2 * yC, \
  r'$M_1={:.2f}\,M_\odot, M_2={:.2f}\,M_\odot$'.format(M1,M2), fontsize = 15 )

ax.text( 1.0e-2 * xC, 3.0e+1 * yC, \
            r'$z = {:.2e}, D_L \approx {:.2e}\,Mpc$'.format( \
                                          z, r * ( 1.0 + z ) ), fontsize = 15 )
ax.text( 1.0e-2 * xC, 1.0e+1 * yC, \
            r'$f_{:}={:.2e}\,Hz$'.format('{ISCO}',f_ISCO), fontsize = 15 )
ax.set_xlabel( r'$f\,\left[Hz\right]$' )
ax.set_ylabel( r'$h_{c}\,\left[dimensionless\right]$' )

ax.set_xlim( xlim )

#plt.savefig( BasePath + 'strain/CompareStrainToLISA_{:d}.eps'.format(Snapshot) )
plt.show()

