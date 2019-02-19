#!/usr/local/bin/python3

"""
Bin Mapelli data in bins of lookback-time
"""

import numpy as np
import matplotlib.pyplot as plt
from os.path import isfile

# --- Define paths to files ---
#BasePath   = '/astro1/dunhamsj/'
BasePath   = '/Users/dunhamsj/Research/GW/Sam/'
StrainPath = BasePath + 'strain/'
hcFilePath = StrainPath + 'hc_DataFiles/'

# --- Define files ---
LogFile    = StrainPath + 'hc_LogFile.dat'
BinnedData = BasePath + 'MergerRateDensity/BinnedData_tLb.dat'

# --- Read in log data ---
LogData = np.loadtxt( LogFile )

tLbMin = LogData[-1,9]  # Lookback-Time [Gyr]
tLbMax = LogData[-1,10] # Lookback-Time [Gyr]
print( 'tLbMin = {:} Gyr\ntLbMax = {:} Gyr'.format( tLbMin, tLbMax ) )

# --- Create bins ---
Gyr = 1.0e9
Myr = 1.0e6
# Define bin-width (10 Myr in units of Gyr)
dt = 10 * Myr / Gyr

# Define Bin-edges and bin data
BinEdges = np.arange( tLbMin, tLbMax + dt, dt )

Bins = np.copy( BinEdges[:-1] ) + 0.5 * dt
nBins = len(Bins)
print( 'nBins = {:d}'.format(nBins) )

# --- Read in data files and isolate lookback-time column ---
if not isfile( BinnedData ):
    with open( BinnedData, 'w' ) as f:
        f.write('# Mergers binned by lookback-time in bins of 10 Myr\n')
    SSMin = 26
    SSMax = 135
    for SS in range( SSMin, SSMax+1 ):
        hcFile = hcFilePath +  'hc_{:}.dat'.format(SS)
        if not isfile( hcFile ): continue
        hcDataAll = np.loadtxt( hcFile )
        hcData    = hcDataAll[1:,7] # Lookback-Time [Gyr]
        counts, binedges = np.histogram( hcData, BinEdges )
        with open( BinnedData, 'a' ) as f:
            for i in range( nBins ):
                f.write( '{:} '.format(counts[i]) )
            f.write('\n')

CountData = np.loadtxt( BinnedData )
counts = np.sum( CountData, axis = 0 )

# --- Plotting ---

#'''
# Total number of BHB mergers as a function of lookback-time
fig, ax = plt.subplots()
ax.step( Bins, counts )
ax.set_xlabel( r'$t_{lb}\,\left[Gyr\right]$' )
ax.set_ylabel( r'$N_{BHB}$' )
ax.set_yscale('log')
ax.set_xlim( 0.0, 14.0 )
plt.show()
plt.close()
#'''

'''
# Mapelli, Figure 1
Vc = ( 0.1065 )**3
fig, ax = plt.subplots()
fig.suptitle( 'Mapelli, Figure 1' )
ax.step( Bins, counts / (Vc * 10 * Myr), 'r', label = 'D' )
ax.set_xlabel( r'$t_{lb}\,\left[Gyr\right]$' )
ax.set_ylabel( r'$R_{BHB}\,\left[Gpc^{-3}\,yr^{-1}\right]$' )
ax.set_yscale('log')
ax.set_ylim( 1.0e0, 7.0e3 )
ax.set_xlim( 0.0, 14.0 )
ax.legend()
#plt.savefig( BasePath + 'MergerRateDensity/Mapelli_Figure1_SD.png' )
plt.show()
plt.close()
'''

'''
# --- Compute co-moving volumes for lookback-times ---
from scipy.interpolate import interp1d
Gpc  = 3.086e25 # [ m / Gpc ]
year = 86400.0 * 365.25

z, tLb, rc = np.loadtxt( BasePath + 'strain/z_tLb_r.dat' )
rInterp = interp1d( tLb, rc, kind = 'linear', fill_value = 'extrapolate' )
r = rInterp( Bins * year * Gyr ) / Gpc # [ Gpc ]

# Taking into account expanding universe
Vc = 4.0 / 3.0 * np.pi * r**3 # [ Gpc^3 ]
fig, ax = plt.subplots()
ax.step( Bins, counts / (Vc * 10 * Myr), 'r' )
ax.set_xlabel( r'$t_{lb}\,\left[Gyr\right]$' )
ax.set_ylabel( r'$R_{BHB}\,\left[Gpc^{-3}\,yr^{-1}\right]$' )
ax.set_yscale('log')
#plt.savefig( BasePath + 'MergerRateDensity/RBHB_Vc.png' )
plt.show()
plt.close()
'''
