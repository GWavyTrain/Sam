#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from sys import argv

if ( len( argv ) > 1 ): Nz = argv[ 1 ]
else:                   Nz = str( 26 )

f_LISA , h_LISA = np.loadtxt( 'LISA_sensitivity.dat' , unpack = True )
data = np.loadtxt( 'StrainDataFiles/GW_strain_' + Nz + '.dat' )

z = data[ 0 , 0 ]
D = data[ 0 , 1 ]

year = 60.0 * 60.0 * 24 * 365
Tobs = 10.0 * year

f     = data[ 1 : , 0 ]
h     = data[ 1 : , 1 ] * np.sqrt( Tobs )
h_min = data[ 1 : , 2 ] * np.sqrt( Tobs )
h_max = data[ 1 : , 3 ] * np.sqrt( Tobs )

fig = plt.figure()
plt.title( r'LISA Sensitivity Curve from http://www.srl.caltech.edu/\textasciitilde shane/sensitivity/MakeCurve.html' )

plt.loglog( f_LISA , h_LISA , 'k-' )
plt.loglog( f      , h_max  , 'b-' , label = 'max'  )
plt.loglog( f      , h      , 'g-' , label = 'GOAT' )
plt.loglog( f      , h_min  , 'r-' , label = 'min'  )

plt.text( f_LISA[ 400 ] , h_LISA[ 200 ] , \
            r'$T_{obs}=%i\,yr$' % int( Tobs / year ) , fontsize = 15 )
plt.text( f_LISA[ 400 ] , h_LISA[ 250 ] , \
            r'$z = %.3f, D_L \approx %i\,Mpc$' % ( z , D * ( 1.0 + z ) ) , \
              fontsize = 15 )
plt.xlabel( r'$f\,\left[Hz\right]$' )
plt.ylabel( r'$h_{c}\,\left[Hz^{-1/2}\right]$' )

plt.legend()
#plt.savefig( 'GW_strain_' + Nz + '.png' )
plt.show()

