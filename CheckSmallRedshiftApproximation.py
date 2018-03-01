import numpy as np
import matplotlib.pyplot as plt
np.seterr( divide = 'ignore' , invalid = 'ignore' )

z , tLB = np.loadtxt( 'tLB_z.dat' , unpack = True )

z_cut = 0.1
z    = z[   np.where( z < z_cut )[0] ]
tLB  = tLB[ np.where( z < z_cut )[0] ]

H0 = 70.0 / 3.086e19 * 86400.0 * 365.25 * 1.0e9 # [ 1 / Gyr ]

xlim = ( 0.0 , z_cut )

plt.suptitle( 'Checking small redshift approximation to lookback time' )

plt.subplot( 211 )
plt.plot( z , tLB ,                 label = 'Exact')
plt.plot( z , 1.0 / H0 * z , 'k-' , label = 'Approx' )
plt.axhline( 0.5 , c = 'k' , ls = '--' , label = r'$t_{LB}=0.5$' )
plt.axvline( z[ np.where( tLB < 0.5 )[0][-1] ] , c = 'r' , \
               label = r'$z\left(t_{LB}=0.5\right)$' )
plt.xlabel( r'$z$' , labelpad = -7 )
plt.ylabel( r'$t_{LB}\,\left[Gyr\right]$' )
plt.legend()
plt.xlim( xlim )
plt.ylim( 0.0 , 1.0 )

plt.subplot( 212 )
plt.plot( z , np.abs( tLB - ( 1.0 / H0 * z ) ) / tLB , 'k-')
plt.axvline( z[ np.where( tLB < 0.5 )[0][-1] ] , c = 'r' , \
               label = r'$z\left(t_{LB}=0.5\right)$' )
plt.legend()
plt.xlim( xlim )
plt.xlabel( r'$z$' )
plt.ylabel( r'$\frac{\left|t_{LB}-t_H\times z\right|}{t_{LB}}$' )

#plt.savefig( 'CheckSmallRedshiftApproximation.png' )
plt.show()
