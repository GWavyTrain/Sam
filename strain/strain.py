#!/usr/local/bin/python3

'''
This program computes the characteristic strain for
a set of binary black holes for use with Mapelli data
using formulae for the characteristic strain from
Sesana et al., (2005), ApJ, 623, 23.

If the file doesn't exist, the program first creates
a high-resolution array (resolution is defined by the parameter Nz)
containing redshift, lookback-time, and comoving distance
called 'z_tLb_r.dat'

Link to Mapelli paper (arXiv eprint):
https://arxiv.org/pdf/1708.05722.pdf

Link to Sesana paper:
http://iopscience.iop.org/article/10.1086/428492/pdf

Directory structure (user defines BasePath directory):
BasePath/BBHM_DataFiles_Mapelli/
BasePath/strain/
BasePath/strain/hc_DataFiles/
'''

# --- Import Python functions ---
from numpy import pi, sqrt, loadtxt, copy, linspace, where, \
     savetxt, vstack, empty, insert, inf, nan
from sys import exit
from scipy.interpolate import interp1d
from scipy.integrate import romberg
from time import time
import os

ProgStartTime = time()

# --- Define directories ---

# Root directory
#BasePath = '/astro1/dunhamsj/'
BasePath = '/Users/sam/Research/GW/Sam/'

# Directory to strain files
StrainPath = BasePath + 'strain/'

# Input directory for Mapelli data file
InputDir = BasePath + 'BBHM_DataFiles_Mapelli/'

# Input data file name
FileName = 'time_BHillustris1_'

# Output directory for characteristic strain files
OutputDir = StrainPath + 'hc_DataFiles/'

# Header for strain output file
StrainFileHeader = ' hc_{:d}.dat\n\n\
--- Column names ---\n\
(0) Illustris particle ID\n\
(1) Illustris particle initial mass at formation [1e10 Msun/h] (h = 0.704)\n\
(2) Snapshot redshift [dimensionless]\n\
(3) Merger ID\n\
(4) M1 [Msun]\n\
(5) M2 [Msun]\n\
(6) Delay-time [Gyr]\n\
(7) Lookback-time [Gyr]\n\
(8) Redshift at merge [dimensionless]\n\
(9) Comoving distance at merge [Mpc]\n\
(10) fISCO [Hz]\n\
(11) Characteristic strain [dimensionless] ({:d} entries)\n\
(First row contains LISA frequencies)'

# --- All physical quantities are in SI units ---

# --- Physical constants ---
G    = 6.67e-11
c    = 2.99792458e8
year = 86400.0 * 365.25
Gyr  = 1.0e9 * year
Mpc  = 3.086e22
Msun = 2.0e30

# --- Cosmological parameters from Mapelli et al., (2017), MNRAS, 472, 2422
#     (WMAP-9, Hinshaw et al., (2013), ApJS, 208, 19) ---
# NOTE: Hinshaw et al. values do NOT agree with those quoted in Mapelli
h      = 0.704
H0     = h * 100.0 / 3.086e19
OmegaM = 0.2726
OmegaL = 0.7274

# --- Problem-specific constants ---
LISA_Lifetime = 4.0 * year
tau = LISA_Lifetime

# --- Read in LISA data ---
LISA_data = loadtxt( BasePath + 'LISA_sensitivity.dat' )
ff        = LISA_data[:,0]
hc_LISA   = sqrt( LISA_data[:,1] * ff )

# Shorten LISA data
nLISA   = 12
f       = empty( (nLISA), float )
f[0]    = 2.0e-5
f[-1]   = 1.0
f[1:-1] = copy(ff[1::40])

# --- Compute redshift and comoving distance, given lookback time ---
def CreateHighResArray():

    print( 'Creating high-resolution datafile with redshift, lookback-time, \
and comoving distance: z_tLb_r.dat' )

    # Create high-resolution array of redshifts
    Nz = int( 1.0e2 )
    z_arr = linspace( 0.0, 20.0, Nz )

    # Compute lookback-times for redshift array
    def LookbackTimeIntegrand( z ):
        return 1.0 / ( ( 1.0 + z ) * sqrt( OmegaM * ( 1.0 + z )**3 + OmegaL ) )
    def ComputeLookbackTimeFromRedshift( z ):
        return 1.0 / H0 * romberg( LookbackTimeIntegrand, 0.0, z, \
                                     rtol = 1.0e-8 )

    tLb_arr = [ ComputeLookbackTimeFromRedshift( z ) for z in z_arr ]

    # Interpolate redshift array as a function of lookback-time
    z_interp = interp1d( tLb_arr, z_arr, kind = 'linear' )

    # Compute comoving distances for redshift array
    def ComovingDistanceIntegrand( z ):
        return 1.0 / sqrt( OmegaM * ( 1.0 + z )**3 + OmegaL )
    def ComputeComovingDistanceFromRedshift( z ):
        return c / H0 * romberg( ComovingDistanceIntegrand, 0.0, z, \
                                   rtol = 1.0e-8 )
    r_arr = [ ComputeComovingDistanceFromRedshift( z ) for z in z_arr ]


    savetxt( 'z_tLb_r.dat', vstack( ( z_arr, tLb_arr, r_arr ) ), \
               header = 'z_tLb_r.dat; each row has {:d} entries'.format(Nz) )

if not os.path.isfile( StrainPath + 'z_tLb_r.dat' ):
    CreateHighResArray()

# Interpolate lookback-time to get redshift and comoving distance
z_arr, tLb_arr, r_arr = loadtxt( StrainPath + 'z_tLb_r.dat' )
z_interp = interp1d( tLb_arr, z_arr, kind = 'linear', \
                       fill_value = 'extrapolate' )
r_interp = interp1d( tLb_arr, r_arr, kind = 'linear', \
                       fill_value = 'extrapolate' )

# --- Characteristic strain from Sesana et al., (2005), ApJ, 623, 23 ---
def ComputeCharacteristicStrain( M1, M2, tLb ):

    """
    M1, M2 should be in solar masses
    tLb should be in Gyr
    """

    z = z_interp( tLb * Gyr )
    r = r_interp( tLb * Gyr )

    # --- Compute chirp mass (text below Eq. (2)) ---
    Mc = ( M1 * M2 )**( 3.0 / 5.0 ) / ( M1 + M2 )**( 1.0 / 5.0 ) * Msun

    # --- Compute pure strain (Eq. (2)) ---
    h = 8.0 * pi**( 2.0 / 3.0 ) / sqrt( 10.0 ) \
          * ( G * Mc )**( 5.0 / 3.0 ) / ( c**4 * r ) \
          * ( f * ( 1.0 + z ) )**( 2.0 / 3.0 )

    # --- Compute number of cycles in frequency interval around
    #     given frequency f (Eq. (4)) ---
    n = 5.0 / ( 96.0 * pi**( 8.0 / 3.0 ) ) \
          * c**5 / ( G * Mc * f * ( 1.0  + z ) )**( 5.0 / 3.0 )

    # --- Check whether or not n < f * tau and compute characteristic strain
    #     Eq. (6) and Eq. (7) ---
    hc = where( n < f * tau, h * sqrt( n ), h * sqrt( f * tau ) )

    fISCO = c**3 / ( 6.0**(3.0 / 2.0 ) * pi * G ) \
              / ( ( M1 + M2 ) * Msun ) / ( 1.0 + z )

    return z, r, fISCO, hc


# --- Min/max snapshot for which to compute characteristic strains ---
SSmin = 26
SSmax = 135

# --- Number of parameters before strains to print to file ---
nParams = 11

# --- Initialize global min/max values ---
LogFile = StrainPath + 'hc_LogFile.dat'
LogFileFMT = '{:d} {:.16e} {:d} {:.16e} {:.16e} {:.16e} {:.16e} {:.16e} \
{:.16e} {:.16e} {:.16e} {:.16e} {:.16e} {:.16e} {:.16e}\n'

zSnapshotMinG = tDelayMinG = tLbMinG = zMergeMinG = rMinG = +inf
zSnapshotMaxG = tDelayMaxG = tLbMaxG = zMergeMaxG = rMaxG = -inf
TotMergers = 0

with open( LogFile, 'w' ) as fL:
    fL.write( '# Log file\n' )
    fL.write( \
      '# Snapshot, time to read in file [s], nMergers, time to compute hc [s], \
time to save file, zSnapshot Min/Max, tDelay Min/Max [Gyr], tLb Min/Max [Gyr], \
zMerge Min/Max, r Min/Max [Mpc]\n' )

for SS in range( SSmin, SSmax+1 ):
    
    FileIn  = InputDir + FileName + str( SS ) + '.dat'
    FileOut = OutputDir + 'hc_{:d}.dat'.format(SS)

    if not os.path.isfile( FileIn ): continue

    LoadStartTime = time()
    BBHM_Data = loadtxt( FileIn )
    LoadEndTime = time() - LoadStartTime

    IllustrisID      = BBHM_Data[:,0]
    InitialMass      = BBHM_Data[:,1]
    SnapshotRedshift = BBHM_Data[:,2]
    MergerID         = BBHM_Data[:,5]
    M1               = BBHM_Data[:,6]
    M2               = BBHM_Data[:,7]
    tDelay           = BBHM_Data[:,8]
    tLb              = BBHM_Data[:,9]

    N = len(M1)
    TotMergers += N

    hc = empty( ( N+1, nParams + nLISA ), float )
    hc[0,0:nParams] = 'nan'
    hc[0,nParams:]  = f

    StrainStartTime = time()
    for i in range( N ):
        z, r, fISCO, hc[i+1,nParams:] \
          = ComputeCharacteristicStrain( M1[i], M2[i], tLb[i] )
        hc[i+1,0]  = IllustrisID[i]
        hc[i+1,1]  = InitialMass[i]
        hc[i+1,2]  = SnapshotRedshift[i]
        hc[i+1,3]  = MergerID[i]
        hc[i+1,4]  = M1[i]
        hc[i+1,5]  = M2[i]
        hc[i+1,6]  = tDelay[i]
        hc[i+1,7]  = tLb[i]
        hc[i+1,8]  = z
        hc[i+1,9]  = r / Mpc
        hc[i+1,10] = fISCO
    StrainEndTime = time() - StrainStartTime

    SaveStartTime = time()
    savetxt( FileOut, hc, header = StrainFileHeader.format(SS, nLISA) )
    SaveEndTime = time() - SaveStartTime
    LoadEndTime, N, StrainEndTime, SaveEndTime

    # --- Local min/max ---
    zSnapshotMin = tDelayMin = tLbMin = zMergeMin = rMin = +inf
    zSnapshotMax = tDelayMax = tLbMax = zMergeMax = rMax = -inf

    zSnapshotMin = min( zSnapshotMin, hc[i+1,2].min() )
    zSnapshotMax = max( zSnapshotMax, hc[i+1,2].max() )
    tDelayMin    = min( tDelayMin,    hc[i+1,6].min() )
    tDelayMax    = max( tDelayMax,    hc[i+1,6].max() )
    tLbMin       = min( tLbMin,       hc[i+1,7].min() )
    tLbMax       = max( tLbMax,       hc[i+1,7].max() )
    zMergeMin    = min( zMergeMin,    hc[i+1,8].min() )
    zMergeMax    = max( zMergeMax,    hc[i+1,8].max() )
    rMin         = min( rMin,         hc[i+1,9].min() )
    rMax         = max( rMax,         hc[i+1,9].max() )

    with open( LogFile, 'a' ) as fL:
        fL.write( LogFileFMT.format( \
          SS, LoadEndTime, N, StrainEndTime, SaveEndTime, \
          zSnapshotMin, zSnapshotMax, tDelayMin, tDelayMax, tLbMin, \
          tLbMax, zMergeMin, zMergeMax, rMin, rMax ) )

    zSnapshotMinG = min( zSnapshotMin, zSnapshotMinG )
    zSnapshotMaxG = max( zSnapshotMax, zSnapshotMaxG )
    tDelayMinG    = min( tDelayMin,    tDelayMinG    )
    tDelayMaxG    = max( tDelayMax,    tDelayMaxG    )
    tLbMinG       = min( tLbMin,       tLbMinG       )
    tLbMaxG       = max( tLbMax,       tLbMaxG       )
    zMergeMinG    = min( zMergeMin,    zMergeMinG    )
    zMergeMaxG    = max( zMergeMax,    zMergeMaxG    )
    rMinG         = min( rMin,         rMinG         )
    rMaxG         = max( rMax,         rMaxG         )

ProgEndTime = time() - ProgStartTime

with open( LogFile, 'a' ) as fL:
    fL.write( '\n' )
    fL.write( LogFileFMT.format( \
      999, nan, TotMergers, nan, nan, \
      zSnapshotMinG, zSnapshotMaxG, tDelayMinG, tDelayMaxG, tLbMinG, \
      tLbMaxG, zMergeMinG, zMergeMaxG, rMinG, rMaxG ) )
    fL.write( '# Total run-time: {:.16e} s\n'.format(ProgEndTime) )

os.system( 'rm -f *.pyc' )
os.system( 'rm -rf __pycache__/' )
