#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
np.seterr( divide='ignore' )

bins , rates = np.loadtxt( 'MergerRate.dat' , unpack = True )

plt.step( bins , np.log10( rates ) , 'k' )
plt.xlabel( r'$t_{Lookback\ to\ merger}\,\left[Gyr\right]$' )
plt.ylabel( r'$\ell og_{10}\left(Merger\ Rate\,\left[Gyr^{-1}\right]\right)$' )

plt.title( 'Total merger rate vs. lookback time for BBH mergers\
\nRedshift Snapshots 101-110' )
plt.savefig( 'TotalMergerRate.png' )
#plt.show()
