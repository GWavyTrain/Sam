#!/usr/local/bin/python3

"""
Bin Mapelli data in bins of lookback-time
"""

import numpy as np
import matplotlib.pyplot as plt
from os.path import isfile

# --- Define paths to files ---
#BasePath = '/astro1/dunhamsj/'
BasePath = '/Users/sam/Research/GW/Sam/'
StrainPath = BasePath + 'strain/'

# --- Define files ---
LogFile = StrainPath + 'hc_LogFile.dat'
hcFilePath  = StrainPath + 'hc_DataFiles/'
BinnedData = BasePath + 'BinnedData_tLb.dat'

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
        f.write('# Mergers binned by lookback-time\n')
    SSMin = 26
    SSMax = 27
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
fig, ax = plt.subplots()
ax.step( Bins, counts / (10 * Myr) )
ax.set_xlabel( r'$Lookback-Time\ \left[Gyr\right]$' )
ax.set_ylabel( r'$dN/dt$' )
ax.set_yscale('log')
plt.show()
plt.close()
