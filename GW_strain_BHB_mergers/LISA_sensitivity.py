#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

f , h = np.loadtxt( 'LISA_sensitivity_curve.dat' , unpack = True )

fig = plt.figure()
plt.title( r'LISA Sensitivity Curve from http://www.srl.caltech.edu/~shane/sensitivity/MakeCurve.html' )
plt.loglog( f , h , 'k-' )
plt.xlabel( r'$f\,\left[Hz\right]$' )
plt.ylabel( r'$h_{f}\,\left[Hz^{-1/2}\right]$' )
plt.savefig( 'LISA_sensitivity_curve.png' )
#plt.show()

