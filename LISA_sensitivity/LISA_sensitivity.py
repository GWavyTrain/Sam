#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

f_LISA, PSD_LISA = np.loadtxt( 'LISA_sensitivity.dat', unpack = True )

hn_LISA = np.sqrt( PSD_LISA * f_LISA )

fig = plt.figure()

plt.loglog( f_LISA, hn_LISA, 'k-' )

plt.xlabel( r'$f\,\left[Hz\right]$'                )
plt.ylabel( r'$h_{c}\,\left[dimensionless\right]$' )

plt.xlim( 1.0e-5,  2.0     )
plt.ylim( 5.0e-22, 3.0e-17 )

plt.savefig( 'LISA_sensitivity.png' )
#plt.show()
