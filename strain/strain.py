#!/usr/local/bin/python3

'''
This program computes the characteristic strain for
a set of binary black holes for use with Mapelli data
using formulae for the characteristic strain from
Sesana et al., (2005), ApJ, 623, 23.

Link to Mapelli paper (arXiv eprint):
https://arxiv.org/pdf/1708.05722.pdf

Link to Sesana paper:
http://iopscience.iop.org/article/10.1086/428492/pdf

Directory structure:
BasePath/BBHM_DataFiles_Mapelli/
BasePath/strain/
BasePath/strain/hc_DataFiles/
'''

from numpy import pi, sqrt, loadtxt, copy, linspace, where, \
     savetxt, vstack, empty, insert, inf
from sys import exit
from astropy.cosmology import z_at_value
from scipy.integrate import romberg
from scipy.interpolate import interp1d
from os.path import isfile
from time import time

ProgStartTime = time()

# --- Define root directory ---
#BasePath = '/astro1/dunhamsj/'
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
tau = LISA_Lifetime

# --- Cosmological parameters from Mapelli et al., (2017), MNRAS, 472, 2422
#     (WMAP-9, Hinshaw et al., (2013), ApJS, 208, 19) ---
# NOTE: Hinshaw et al. values do NOT agree with those quoted in Mapelli
h      = 0.704
H0     = h * 100.0 / 3.086e19
OmegaM = 0.2726
OmegaL = 0.7274

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

if not isfile( BasePath + 'strain/z_tLb_r.dat' ):

    print( 'Creating high-resolution datafile with redshift, lookback-time, \
and comoving distance: z_tLb_r.dat' )

    # Create high-resolution array of redshifts
    Nz = 100
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


# Interpolate lookback-time to get redshift and comoving distance
z_arr, tLb_arr, r_arr = loadtxt( 'z_tLb_r.dat' )
z_interp = interp1d( tLb_arr, z_arr, kind = 'linear', \
                       fill_value = 'extrapolate' )
r_interp = interp1d( tLb_arr, r_arr, kind = 'linear', \
                       fill_value = 'extrapolate' )

# --- Characteristic strain from Sesana et al., (2005), ApJ, 623, 23 ---
def ComputeCharacteristicStrain( M1, M2, tLb ):

    """
    M1, M2 should be given in solar masses
    tLb should be given in Gyr
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

# --- Output directory for characteristic strain files  ---
OutputDir = BasePath + 'strain/hc_DataFiles/'

# --- Mapelli data file input directory ---
InputDir = BasePath + 'BBHM_DataFiles_Mapelli/'

# --- Data file ---
FileName = 'time_BHillustris1_'

SSmin = 26
SSmax = 135

SS = SSmin

nParams = 11

Min_z      = +inf
Max_z      = -inf
Min_tLb    = +inf
Max_tLb    = -inf
TotMergers = 0

while( SS <= SSmax ):

    FileIn = InputDir + FileName + str( SS ) + '.dat'
    FileOut = OutputDir + 'hc_{:d}.dat'.format(SS)

    if isfile( FileIn ):

        LogFileName = OutputDir + 'log_{:d}.txt'.format( SS )
        with open( LogFileName, 'w' ) as fL:
            fL.write( '# log_{:d}.txt\n'.format( SS ) )

        LoadStartTime = time()
        BBHM_Data = loadtxt( FileIn )

        IllustrisID      = BBHM_Data[:,0]
        InitialMass      = BBHM_Data[:,1]
        SnapshotRedshift = BBHM_Data[:,2]
        MergerID         = BBHM_Data[:,5]
        M1               = BBHM_Data[:,6]
        M2               = BBHM_Data[:,7]
        tDelay           = BBHM_Data[:,8]
        tLb              = BBHM_Data[:,9]

        LoadEndTime = time() - LoadStartTime

        with open( LogFileName, 'a' ) as fL:
            fL.write( \
              'Time to read in file: {:.3e} s\n'.format( LoadEndTime ) )

        N = len(M1)
        TotMergers += N

        with open( LogFileName, 'a' ) as fL:
            fL.write( \
              'Number of mergers:    {:d}\n'.format( N ) )

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

        Min_tLb = min( Min_tLb, hc[1:,7].min() )
        Max_tLb = max( Max_tLb, hc[1:,7].max() )
        Min_z   = min( Min_z,   hc[1:,8].min() )
        Max_z   = max( Max_z,   hc[1:,8].max() )

        with open( LogFileName, 'a' ) as fL:
            fL.write( \
              'Time to compute hc:   {:.3e} s\n'.format( StrainEndTime ) )
            fL.write( '\
Min lookback-time:    {:.10f}\n\
Max lookback-time:    {:.10f}\n\
Min merger redshift:  {:.16e}\n\
Max merger redshift:  {:.16e}\n'.format( \
            hc[1:,7].min(), hc[1:,7].max(), hc[1:,8].min(), hc[1:,8].max() ) )
        SaveStartTime = time()
        savetxt( FileOut, hc, \
                   header = ' hc_{:d}.dat\n\n\
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
(First row contains LISA frequencies)'.format(SS, nLISA) )
        SaveEndTime = time() - SaveStartTime
        with open( LogFileName, 'a' ) as fL:
            fL.write( \
              'Time to save file:    {:.3e} s\n'.format( SaveEndTime ) )

    SS += 1

ProgEndTime = time() - ProgStartTime

with open( OutputDir + 'A_log.txt', 'w' ) as fL:
    fL.write( 'Total number of mergers: {:d}\n\n'.format( TotMergers ) )
    fL.write( 'Minimum lookback-time:   {:.10f} Gyr\n'.format( Min_tLb ) )
    fL.write( 'Maximum lookback-time:   {:.10f} Gyr\n'.format( Max_tLb ) )
    fL.write( 'Minimum merger redshift: {:.16e}\n'.format( Min_z ) )
    fL.write( 'Maximum merger redshift: {:.16e}\n\n'.format( Max_z ) )
    fL.write( 'Total run-time: {:.3e} s'.format( ProgEndTime ) )
