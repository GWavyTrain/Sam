#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

f , h       = np.loadtxt( 'LISA_sensitivity.dat' , unpack = True )
f_26 , h_26 = np.loadtxt( 'GW_strain_26.dat' , unpack = True )

fig = plt.figure()
plt.title( r'LISA Sensitivity Curve from http://www.srl.caltech.edu/~shane/sensitivity/MakeCurve.html' )

plt.loglog( f , h , 'k-' )
plt.loglog( f_26 , h_26 , 'r-' )

plt.xlabel( r'$f\,\left[Hz\right]$' )
plt.ylabel( r'$h_{f}\,\left[Hz^{-1/2}\right]$' )

#plt.savefig( 'LISA_sensitivity_curve.png' )
plt.show()
