import numpy as np
np.seterr( divide='ignore' )
import matplotlib.pyplot as plt
from sys import exit

#var = 'LookbackTime'
var = 'SnapshotRedshift'
if var == 'LookbackTime':
    nBins = 137
elif var == 'SnapshotRedshift':
    nBins = 1360
else:
    print( 'Invalid variable: {0}'.format(var) )
    exit('Exiting...')

rates = np.zeros( (nBins), float )

for chunk in range( 1, 79 ):

    bins, rates_c = np.loadtxt( \
           'ChunkedMergers/nMergersBinned_' + var + '/nMergersBinned_' \
           + str(chunk) + '.dat', unpack = True )
    rates += rates_c

#'''
plt.step( bins, rates, 'k' )
plt.yscale( 'log' )

if var == 'LookbackTime':
    plt.xlabel( r'$Lookback\ time\ \left[Gyr\right]$' )
    plt.ylabel( r'$dN_{mergers}/dt_{LB}$' )
    plt.title( \
    'Total number of BBH mergers\nas a function of lookback time to merger' )

elif var == 'SnapshotRedshift':
    plt.xlabel( r'$Snapshot\ Redshift$' )
    plt.ylabel( r'$dN_{mergers}/dz$' )
    plt.title( \
    'Total number of BBH mergers\nas a function of snapshot redshift' )

plt.savefig( 'Mergers_BinnedBy' + var + '.png' )
#plt.show()
plt.close()
#'''
