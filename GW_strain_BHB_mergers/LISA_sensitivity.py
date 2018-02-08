#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from sys import argv

if ( len( argv ) > 1 ): Nz = argv[ 1 ]
else:                   Nz = str( 26 )

f_LISA , h_LISA = np.loadtxt( 'LISA_sensitivity.dat' , unpack = True )
data  = np.loadtxt( 'StrainDataFiles/GW_strain_' + Nz + '.dat' )

z = data[ 0 , 0 ]
D = data[ 0 , 1 ]

f     = data[ 1 : , 0 ]
h     = data[ 1 : , 1 ]
h_min = data[ 1 : , 2 ]
h_max = data[ 1 : , 3 ]

fig = plt.figure()
plt.title( r'LISA Sensitivity Curve from http://www.srl.caltech.edu/\textasciitilde shane/sensitivity/MakeCurve.html' )

plt.loglog( f_LISA , h_LISA , 'k-' )
plt.loglog( f    , h_max , 'b-' , label = 'max'  )
plt.loglog( f    , h     , 'g-' , label = 'mean' )
plt.loglog( f    , h_min , 'r-' , label = 'min'  )

plt.text( f_LISA[ 400 ] , h_LISA[ 200 ] , r'$T_{obs}=10\ yr$' , fontsize = 15 )
plt.text( f_LISA[ 400 ] , h_LISA[ 250 ] , \
            r'$z = %.3f, D_L \approx %i\,Mpc$' % ( z , D ) , \
              fontsize = 15 )
plt.xlabel( r'$f\,\left[Hz\right]$' )
plt.ylabel( r'$h_{f}\,\left[Hz^{-1/2}\right]$' )

plt.legend()
plt.savefig( 'GW_strain_' + Nz + '.png' )
#plt.show()
