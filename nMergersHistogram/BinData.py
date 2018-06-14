# Python script to bin merger data

import numpy as np

ChunkNumber = '1'

# Read in data from mergers.dat
ID, z, M1, M2, t_delay, tLb \
  = np.loadtxt( 'UniqueMergers_chunk' + ChunkNumber + '.dat', unpack = True )

# Choose variables to bin and bin-width
x  = tLb
dx = 0.01

# Get min and max of variable to create bins/bin edges
xMin = np.min( x )
xMax = np.max( x )

# Get total number of mergers
nMergers = len( x )

# Create bins/bin edges
nBins = np.int64( np.ceil( ( xMax - xMin ) / dx ) )

BinEdges = np.empty( nBins+1, float )
for i in range( len( BinEdges ) ):
    BinEdges[i] = xMin + i * dx

# Create and populate counts array
Counts = np.zeros( nBins, np.int64 )
for i in range( nBins ):
    for iMerger in range( nMergers ):
        if ( ( x[iMerger] >= BinEdges[i] ) & ( x[iMerger] <= BinEdges[i+1] ) ):
            Counts[i] += 1

# Shift bins to center values
Bins = BinEdges[0:-1] + 0.5 * dx

data = np.vstack( ( Bins, Counts ) )
np.savetxt( 'BinnedData_chunk' + ChunkNumber + '.dat', data.T, \
               header = 'Bin center, counts' )
