from numpy import pi, sqrt, loadtxt, empty
from sys import exit
from astropy.cosmology import z_at_value
from scipy.integrate import romberg

# --- Define root directory ---
BasePath = '/Users/sam/Research/GW/Sam/'

# --- All physical quantities are in SI units ---

# --- Physical constants ---
G    = 6.67e-11
c    = 2.99792458e8
year = 86400.0 * 365.25
Gyr  = 1.0e9 * year
Mpc  = 3.086e22
Msun = 2.0e30

# --- Problem-specific constants ---
LISA_Lifetime = 4.0 * year
tau           = LISA_Lifetime

# --- Cosmological parameters from Mapelli et al., (2017), MNRAS, 472, 2422
#     (WMAP-9, Hinshaw et al., (2013), ApJS, 208, 19) ---
# NOTE: Hinshaw et al. values do NOT agree with those quoted in Mapelli
h      = 0.704
H0     = h * 100.0 / 3.086e19
OmegaM = 0.2726
OmegaL = 0.7274
OmegaB = 0.0456

# --- Read in LISA data ---
LISA_data = loadtxt( BasePath + 'LISA_sensitivity.dat' )
f         = LISA_data[:,0]
hc_LISA   = sqrt( LISA_data[:,1] * f )

# --- Compute comoving distance, given lookback time ---
# Find redshift from lookback time, then find comoving distance from redshift

# astropy.cosmology.z_at_value requires function
# that computes lookback time, given redshift

def ComputeLookbackTimeFromRedshift( z ):
    def LookbackTimeIntegrand( z ):
        return 1.0 / ( ( 1.0 + z ) \
                       * sqrt( ( OmegaM + OmegaB ) * ( 1.0 + z )**3 + OmegaL ) )

    return 1.0 / H0 * romberg( LookbackTimeIntegrand, 0.0, z )

def ComputeRedshiftFromLookbackTime( tLb ):
    return z_at_value( ComputeLookbackTimeFromRedshift, tLb )

def ComputeComovingDistanceFromRedshift( z ):
    def ComovingDistanceIntegrand( z ):
        return 1.0 / sqrt( ( OmegaM + OmegaB ) * ( 1.0 + z )**3 + OmegaL )

    return c / H0 * romberg( ComovingDistanceIntegrand, 0.0, z )

# --- Combine above two functions to compute comoving
#     distance from lookback time ---
def ComputeComovingDistanceFromLookbackTime( tLb ):
    z = ComputeRedshiftFromLookbackTime( tLb )
    return ComputeComovingDistanceFromRedshift( z )

# --- Characteristic strain from Sesana et al., (2005), ApJ, 623, 23 ---
def ComputeCharacteristicStrain( M1, M2, r ):

    # --- Initialize array to hold characteristic strains ---
    hc = empty( len( f ), float )
    
    # --- Compute chirp mass ---
    Mc = ( M1 * M2 )**( 3.0 / 5.0 ) / ( M1 + M2 )**( 1.0 / 5.0 )

    # --- Compute pure strain (Eq. (2)) ---
    h = 8.0 * pi**( 2.0 / 3.0 ) / sqrt( 10.0 ) \
          * ( G * Mc )**( 5.0 / 3.0 ) / ( c**4 * r ) * f**( 2.0 / 3.0 )

    # --- Compute number of cycles in frequency interval around
    #     given frequency f (Eq. (4)) ---
    n = 5.0 / ( 96.0 * pi**( 8.0 / 3.0 ) ) \
          * c**5 / ( G * Mc * f )**( 5.0 / 3.0 )

    # --- Check whether or not n < f * tau and compute characteristic strain ---
    for i in range( len( f ) ):
        if ( n[i] < f[i] * tau ):
            hc[i] = h[i] * sqrt( n[i] ) # Eq. (6)
        else:
            hc[i] = h[i] * sqrt( f[i] * tau ) # Eq. (7)

    return hc

# --- Test characteristic strain calculation: Sesana (2016), PRL, 116, 231102
#     (GW150914) ---
import matplotlib.pyplot as plt
fig, ax = plt.subplots()

M1 = 36.0 * Msun
M2 = 29.0 * Msun
z  = 0.09

#ax.loglog( f, hc_LISA, 'k-' )
ax.loglog( f, ComputeCharacteristicStrain \
                 ( M1, M2, ComputeComovingDistanceFromRedshift( z ) ) )
ax.set_xlim( 1.0e-3, 3.0e3 )
ax.set_ylim( 1.0e-23, 1.0e-18 )

plt.show()
