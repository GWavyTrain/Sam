import numpy as np
import matplotlib.pyplot as plt
np.seterr( divide='ignore' )

bins, rates = np.loadtxt( 'nMergersBinned.dat', unpack = True )

plt.step( bins, rates, 'k' )
#plt.yscale( 'log' )
plt.xlabel( r'$Snapshot\ Redshift\ z$' )
plt.ylabel( r'$dN_{mergers}/dz$' )

#plt.savefig( 'MergerRateDensity' + SnapShotRange + '.png' )
plt.show()
plt.close()
