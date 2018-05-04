#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
np.seterr( divide='ignore' )

SnapshotRange = '_101-110'

bins, MergerRateDensity \
  = np.loadtxt( 'MergerRateDensity' + SnapshotRange + '.dat', unpack = True )

plt.step( bins, MergerRateDensity, 'k' )
plt.yscale( 'log' )
plt.xlabel( r'$t_{Lookback\ to\ merger}\,\left[Gyr\right]$' )
plt.ylabel( r'$Merger\ Rate\,Density\,\left[Gpc^{-3}\,yr^{-1}\right]$' )

plt.title( 'Merger rate density vs. lookback time for BBH mergers\
\nRedshift Snapshots ' + SnapshotRange[1:] )
#plt.savefig( 'MergerRateDensity' + SnapshotRange + '.png' )
plt.show()
plt.close()
