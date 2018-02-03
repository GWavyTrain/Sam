#!/usr/bin/python

import numpy as np
from scipy.integrate import romberg

# Physical constants
G = 6.67e-8 # [ cm^3 / ( g * s^2 ) ]
c = 3.0e10  # [ cm / s ]

# Cosmology
H0      = 70.4 / 3.086e19 # [ s^(-1) ]
Omega_m = 0.3             # [ dimensionless ]
Omega_L = 0.7             # [ dimensionless ]

# Conversion factors
Msun = 2.0e33          # [ g / Msun ]
year = 86400.0 * 365.0 # [ s / yr ]

# Compute luminosity distance
def integrand( z ):

    return 1.0 / np.sqrt( Omega_m * ( 1.0 + z )**( 3.0 ) + Omega_L )

def D_L( z ):

    return c / H0 * romberg( integrand , 0.0 , z ) * ( 1.0 + z )

# Compute strain
def hf( m1 , m2 , z , forb , Tobs ):

    m1   *= Msun # [ g ]
    m2   *= Msun # [ g ]
    Tobs *= year # [ s ]

    f = 2.0 * forb # [ Hz ]
    D = D_L( z )   # [ cm ]

    Mc = ( m1 * m2 )**( 3.0 / 5.0 ) / ( m1 + m2 )**( 1.0 / 5.0 )

    return G / c**( 2.0 ) * Mc / D \
             * ( G / c**( 3.0 ) * np.pi * f * Mc )**( 2.0 / 3.0 ) \
               * np.sqrt( Tobs )


