#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

# Load in data
f_LISA , RootPSD_LISA = np.loadtxt( 'LISA_sensitivity.dat' , unpack = True )
data                  = np.loadtxt( 'GW_strain_test_GOAT.dat' )

# Compute noise amplitude for LISA
hn_LISA = np.sqrt( f_LISA ) * RootPSD_LISA

# Get redshift z and comoving distance r
z = data[ 1 , 3 ]
r = data[ 1 , 4 ]

f  = data[ 0 , 5 : ]
hc = data[ 1 , 5 : ]

# Plotting
fig = plt.figure()
plt.title( r'LISA Sensitivity Curve from http://www.srl.caltech.edu/\textasciitilde shane/sensitivity/MakeCurve.html' )

plt.loglog( f_LISA , hn_LISA , 'k-'                  )
plt.loglog( f      , hc      , 'b-' , label = 'GOAT' )

plt.text( f_LISA[ 300 ] , np.sqrt( f_LISA[ 250 ] ) * RootPSD_LISA[ 250 ] , \
            r'$z = %.3f, D_L \approx %i\,Mpc$' % ( z , r * ( 1.0 + z ) ) , \
              fontsize = 15 )
plt.xlabel( r'$f\,\left[Hz\right]$' )
plt.ylabel( r'$h_{c}\,\left[dimensionless\right]$' )

plt.xlim( 1.0e-5  , 2.0     )
plt.ylim( 1.0e-21 , 1.0e-16 )

plt.legend()
plt.show()

