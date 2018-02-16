#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

f_LISA , RootPSD_LISA = np.loadtxt( 'LISA_sensitivity.dat' , unpack = True )

fig = plt.figure()
plt.title( r'LISA Sensitivity Curve from http://www.srl.caltech.edu/\textasciitilde shane/sensitivity/MakeCurve.html' )

plt.loglog( f_LISA , np.sqrt( f_LISA ) * RootPSD_LISA , 'k-' )

plt.xlabel( r'$f\,\left[Hz\right]$'                )
plt.ylabel( r'$h_{c}\,\left[dimensionless\right]$' )

plt.xlim( 1.0e-5  , 2.0     )
plt.ylim( 1.0e-21 , 1.0e-16 )

plt.show()
