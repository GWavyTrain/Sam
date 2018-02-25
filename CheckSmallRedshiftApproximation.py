import numpy as np
import matplotlib.pyplot as plt

z , tLB = np.loadtxt( 'tLB_z.dat' , unpack = True )

H0 = 70.0 / 3.086e19 * 86400.0 * 365.25 * 1.0e9 # [ 1 / Gyr ]

xlim = ( 0.0 , 0.5 )

plt.suptitle( 'Checking small redshift approximation to lookback time' )

plt.subplot( 211 )
plt.plot( z , tLB ,                 label = 'Exact')
plt.plot( z , 1.0 / H0 * z , 'k-' , label = 'Approx' )
plt.ylabel( r'$t_{LB}\,\left[Gyr\right]$' )
plt.legend()
plt.xlim( xlim )
plt.ylim( 0.0 , 7.0 )

plt.subplot( 212 )
plt.plot( z , np.log10( ( np.abs( tLB - ( 1.0 / H0 * z ) ) ) / tLB ) , 'k-')
plt.xlim( xlim )
plt.xlabel( r'$z$' )
plt.ylabel( \
r'$\ell og_{10}\left[\frac{\left|t_{LB}-t_H\times z\right|}{t_{LB}}\right]$' , \
labelpad = -5 )

plt.savefig( 'CheckSmallRedshiftApproximation.png' )
#plt.show()
