#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

# Load in data
f_LISA , PSD_LISA = np.loadtxt( 'LISA_sensitivity.dat' , unpack = True )
data              = np.loadtxt( 'hc.dat' )

# Compute noise amplitude for LISA
hn_LISA = np.sqrt( PSD_LISA * f_LISA )

Merger = 3

# Get redshift z and comoving distance r
z = data[ Merger , 4 ]
r = data[ Merger , 5 ]

# Get frequencies and strain
f      = data[ 0      , 7 : ]
hc     = data[ Merger , 7 : ]
f_ISCO = data[ Merger , 5   ]

# Get masses
M1 = data[ Merger , 1 ]
M2 = data[ Merger , 2 ]

# Plotting
fig = plt.figure()

plt.loglog( f_LISA , hn_LISA , 'k-' )
plt.loglog( f      , hc      , 'b-' )

plt.text( 2.0e-4 , 1.0e-17 , \
            r'$M_{1} = %.2f\,M_{\odot}, M_{2} = %.2f\,M_{\odot}$' \
              % ( M1 , M2 ) , fontsize = 15 )
plt.text( 2.0e-4 , 3.0e-18 , \
            r'$z = %.5f, D_L \approx %i\,Mpc,f_{ISCO}=%.3f\,Hz$' \
              % ( z , r * ( 1.0 + z ) , f_ISCO ) , fontsize = 15 )
plt.xlabel( r'$f\,\left[Hz\right]$' )
plt.ylabel( r'$\ell og_{10}\ h_{c}\,\left[dimensionless\right]$' )

plt.xlim( 1.0e-5  , 2.0     )
#plt.ylim( 5.0e-22 , 3.0e-17 )

#plt.savefig( 'CompareStrainToLISA.png' )
plt.show()

