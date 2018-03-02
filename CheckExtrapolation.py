import numpy as np
import matplotlib.pyplot as plt
np.seterr( divide = 'ignore' , invalid = 'ignore' )

z , tLB = np.loadtxt( 'tLB_z.dat' , unpack = True )

H0 = 70.0 / 3.086e19 * 86400.0 * 365.25 * 1.0e9 # [ 1 / Gyr ]

z_20   = z  [ np.where( z < 20.0 )[0] ]
tLB_20 = tLB[ np.where( z < 20.0 )[0] ]

z_20_100   = z  [ np.where( z > 20.0 )[0] ]
tLB_20_100 = tLB[ np.where( z > 20.0 )[0] ]

imin   = 2
imax   = 1
zmin   = z_20[-imin]
zmax   = z_20[-imax]
tLBmin = tLB_20[-imin]
tLBmax = tLB_20[-imax]

zExMax = z[-1]

m = ( tLBmax - tLBmin ) / ( zmax - zmin )
b = 0.5 * ( ( tLBmax - m * zmax ) + ( tLBmin - m * zmin ) )

zEx = np.linspace( zmax , zExMax , len( z_20_100 ) )

tLBEx = m * zEx  + b

xlim = ( z_20[-1] , z[-1] )

plt.suptitle( 'Checking accuracy of extrapolation' )

plt.subplot( 211 )
plt.plot( z   , tLB   , 'k-'  , label = 'Exact'         )
plt.plot( zEx , tLBEx , 'k--' , label = 'Extrapolation' )

plt.xlabel( r'$z$' , labelpad = -7 )
plt.ylabel( r'$t_{LB}\,\left[Gyr\right]$' )
plt.xlim( xlim )

plt.legend()

plt.subplot( 212 )
plt.plot( zEx , ( np.abs( tLB_20_100 - tLBEx ) / tLB_20_100 ) , 'k-' )
plt.xlabel( r'$z$' , labelpad = -7 )
plt.ylabel( r'$\frac{\left|t_{LB}-t_{LB,Ex}\right|}{t_{LB}}$' )
plt.xlim( xlim )
plt.savefig( 'CheckExtrapolation.png' )
#plt.show()
