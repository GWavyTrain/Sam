#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

# Load in data
f_LISA , PSD_LISA = np.loadtxt( 'LISA_sensitivity.dat' , unpack = True )
data              = np.loadtxt( 'GW_strain_test_GOAT.dat' )

# Compute noise amplitude for LISA
hn_LISA = np.sqrt( PSD_LISA * f_LISA )

# Get redshift z and comoving distance r
z = data[ 1 , 3 ]
r = data[ 1 , 4 ]

# Get frequencies and strain
f  = data[ 0 , 6 : ]
hc = data[ 1 , 6 : ]

# Get masses
M1 = data[ 1 , 1 ]
M2 = data[ 1 , 2 ]

# Plotting
fig = plt.figure()

plt.loglog( f_LISA , hn_LISA , 'k-'                    )
plt.loglog( f      , hc      , 'b-' , label = 'Sesana' )

plt.text( f_LISA[ 150 ] , hn_LISA[ 350 ] , \
            r'$z = %.3f, D_L \approx %i\,Mpc$' % ( z , r * ( 1.0 + z ) ) , \
              fontsize = 15 )
plt.xlabel( r'$f\,\left[Hz\right]$' )
plt.ylabel( r'$h_{c}\,\left[dimensionless\right]$' )

plt.xlim( 1.0e-5  , 2.0     )
plt.ylim( 5.0e-22 , 3.0e-17 )

plt.legend()
plt.show()
