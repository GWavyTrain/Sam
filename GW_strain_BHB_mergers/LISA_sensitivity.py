#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

f , h = np.loadtxt( 'LISA_sensitivity.dat' , unpack = True )
data  = np.loadtxt( 'GW_strain_135.dat' )

z = data[ 0 , 0 ]
D = data[ 0 , 1 ]

f_26  = data[ 1 : , 0 ]
h_26  = data[ 1 : , 1 ]
h_min = data[ 1 : , 2 ]
h_max = data[ 1 : , 3 ]

fig = plt.figure()
plt.title( r'LISA Sensitivity Curve from http://www.srl.caltech.edu/~shane/sensitivity/MakeCurve.html' )

plt.loglog( f , h , 'k-' )
plt.loglog( f_26 , h_26 ,  'r-' , label = 'mean' )
plt.loglog( f_26 , h_min , 'g-' , label = 'min'  )
plt.loglog( f_26 , h_max , 'b-' , label = 'max'  )

plt.text( f[ 400 ] , h[ 200 ] , r'$z = %.3f, D \approx %i\,Mpc$' % ( z , D ) , \
            fontsize = 15 )
plt.xlabel( r'$f\,\left[Hz\right]$' )
plt.ylabel( r'$h_{f}\,\left[Hz^{-1/2}\right]$' )

plt.legend()
#plt.savefig( 'LISA_sensitivity_curve.png' )
plt.show()
