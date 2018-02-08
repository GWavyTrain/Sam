#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

f_LISA , h_LISA = np.loadtxt( 'LISA_sensitivity.dat' , unpack = True )

fig = plt.figure()
plt.title( r'LISA Sensitivity Curve from http://www.srl.caltech.edu/\textasciitilde shane/sensitivity/MakeCurve.html' )

plt.loglog( f_LISA , h_LISA , 'k-' )

plt.xlabel( r'$f\,\left[Hz\right]$' )
plt.ylabel( r'$h_{f}\,\left[Hz^{-1/2}\right]$' )

plt.savefig( 'LISA_sensitivity.png' )
#plt.show()
