#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

# Load in data
f_LISA , PSD_LISA = np.loadtxt( 'LISA_sensitivity.dat' , unpack = True )
dataS             = np.loadtxt( 'hc_short.dat' )
data              = np.loadtxt( 'hc.dat' )

# Compute noise amplitude for LISA
hn_LISA = np.sqrt( PSD_LISA * f_LISA )

Merger = -3

# Get redshift z and comoving distance r
z = data[ Merger , 4 ]
r = data[ Merger , 5 ]

# Get frequencies and strain
f      = data [ 0      , 7 : ]
hc     = data [ Merger , 7 : ]
fS     = dataS[ 0      , 7 : ]
hcS    = dataS[ Merger , 7 : ]
f_ISCO = data [ Merger , 5   ]

# Get masses
M1 = data[ Merger , 1 ]
M2 = data[ Merger , 2 ]

# Plotting
fig = plt.figure()
ax1 = plt.subplot(211)

ax1.loglog( f_LISA , hn_LISA , 'k-' )
ax1.loglog( f      , hc      , 'b-' , label = 'Nf=400' )
ax1.loglog( fS     , hcS     , 'r-' , label = 'Nf=10'  )

ax1.text( 2.0e-4 , 1.0e-17 , \
            r'$M_{1} = %.2f\,M_{\odot}, M_{2} = %.2f\,M_{\odot}$' \
              % ( M1 , M2 ) , fontsize = 15 )
ax1.text( 2.0e-4 , 3.0e-18 , \
            r'$z = %.5f, D_L \approx %i\,Mpc,f_{ISCO}=%.3f\,Hz$' \
              % ( z , r * ( 1.0 + z ) , f_ISCO ) , fontsize = 15 )
ax1.set_xticks( [] )
ax1.set_ylabel( r'$\ell og_{10}\ h_{c}\,\left[dimensionless\right]$' , fontsize = 15 )
ax1.legend( loc = 4 )
ax1.set_xlim( 1.0e-5 , 2.0 )

### Check error when using Nf=10 instead of Nf=400
from scipy.interpolate import interp1d
ax2 = plt.subplot(212)

hc  = interp1d( f  , hc  )
hcS = interp1d( fS , hcS )

ax2.semilogx( f , np.abs( ( hc( f ) - hcS( f ) ) / hc( f ) ) , 'k-' )
ax2.set_xlabel( r'$f\,\left[Hz\right]$' )
ax2.set_ylabel( r'$\left|\frac{h_{c}-h_{c,short}}{h_{c}}\right|$', labelpad = -2 )
ax2.set_xlim( 1.0e-5 , 2.0 )

plt.subplots_adjust( hspace = 0 )
plt.savefig( 'CheckNf10.png' )
