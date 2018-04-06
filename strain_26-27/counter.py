import numpy as np

tLB,counts = np.loadtxt( 'MergerRateCounts.dat' , unpack = True )

print np.sum( counts )
