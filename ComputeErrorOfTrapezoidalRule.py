#/usr/bin/python

# Compute estimated error for Trapezoidal rule in lookback time calculation

import numpy as np

Nq = 100000

zi = 0.0
zf = 20.0

dz = ( zf - zi ) / Nq

OMEGA_M = 0.3
OMEGA_L = 0.7

def E( z ):
    
    return np.sqrt( OMEGA_M * ( 1.0 + z )**3 + OMEGA_L )

def dE( z ):
    
    return 3.0 / 2.0 * OMEGA_M * ( 1.0 + z )**2 / E( z )

def ddE( z ):
    
    return 3.0 / E( z ) * ( 3.0 / 4.0 * OMEGA_M**2 \
                         * ( 1.0 + z )**4 / E( z )**2 + OMEGA_M * ( 1.0 + z ) )
def f(z):
    
    return 1.0 / ( ( 1.0 + z ) * E( z ) )

def ddf( z ):
    
    return f( z )**2 * ( 2.0 * f( z ) * ( E( z ) + ( 1.0 + z ) * dE( z ) )**2 \
                         - 2.0 * dE( z ) - ( 1.0 + z ) * ddE( z ) )

def TrapError( z ):
    
    return -dz**2 * ( zf - zi ) / 12.0 * ddf( z )

z = np.linspace( zi , zf , Nq )
print '|ERROR|:' , np.max( np.abs( TrapError( z ) ) )
