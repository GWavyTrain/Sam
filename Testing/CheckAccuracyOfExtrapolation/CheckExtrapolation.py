import numpy as np
import matplotlib.pyplot as plt
np.seterr( divide = 'ignore', invalid = 'ignore' )

z, tLb = np.loadtxt( 'tLb_z_0-100.dat', unpack = True )

H0 = 70.0 / 3.086e19 * 86400.0 * 365.25 * 1.0e9 # [ 1 / Gyr ]

# Isolate values with z <= 20
z_LE20   = z  [ np.where( z <= 20.0 )[0] ]
tLb_LE20 = tLb[ np.where( z <= 20.0 )[0] ]

# Isolate values with z > 20
z_GT20   = z  [ np.where( z > 20.0 )[0] ]
tLb_GT20 = tLb[ np.where( z > 20.0 )[0] ]

# Get values to use for extrapolation
zMin   = z_LE20  [-2]
zMax   = z_LE20  [-1]
tLbMin = tLb_LE20[-2]
tLbMax = tLb_LE20[-1]

zExMax = z[-1]

# Get slope and intercept for linear extrapolation
m = ( tLbMax - tLbMin ) / ( zMax - zMin )
b = 0.5 * ( ( tLbMax - m * zMax ) + ( tLbMin - m * zMin ) )

# Create redshift range to extrapolate
zEx = np.linspace( zMax, zExMax, len( z_GT20 ) )

# Extrapolate lookback time
tLbEx = m * zEx  + b

xlim = ( z_LE20[-1], z[-1] )

# Plotting
plt.suptitle( 'Checking accuracy of extrapolation' )

# Extrapolated values
plt.subplot( 211 )
plt.plot( z,   tLb,   'k-',  label = 'Exact'         )
plt.plot( zEx, tLbEx, 'k--', label = 'Extrapolation' )

plt.xlabel( r'$z$' , labelpad = -7 )
plt.ylabel( r'$t_{Lb}\,\left[Gyr\right]$' )
plt.xlim( xlim )

plt.legend()

# Fractional error in using extrapolation
plt.subplot( 212 )
plt.plot( zEx, ( np.abs( tLb_GT20 - tLbEx ) / tLb_GT20 ), 'k-' )
plt.xlabel( r'$z$', labelpad = -7 )
plt.ylabel( r'$\frac{\left|t_{Lb}-t_{Lb,Extrap}\right|}{t_{Lb}}$' )
plt.xlim( xlim )
#plt.savefig( 'CheckExtrapolation.png' )
plt.show()
