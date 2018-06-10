# Program to get array indices of unique mergers from 'TrimmedAndCombinedData.dat'

import numpy as np
from time import time
from datetime import datetime

##### File I/O
ProgName = 'GetUniqueMergers.py'
FileOut  = 'UniqueMergers.dat'
FileIn   = 'TrimmedAndCombinedData.dat'

with open( FileOut, 'w' ) as f:
    f.write( '# {:}\n\n'.format(FileOut) )
    f.write( '# These data were obtained with {:}\n'.format(ProgName) )
    f.write( '# Program start-time: {:}\n\n'.format(datetime.now()) )
    f.write( '# Indices of unique mergers (sorted by Merger ID)\n' )

StartTime = time()

IDs = np.loadtxt( FileIn, usecols = (0), dtype = 'float' )
UniqueIndices = np.unique(IDs, return_index = True)[1]

EndTime = time()

TotRunTime = EndTime - StartTime

with open( FileOut, 'a' ) as f:
    f.write( '# Total number of mergers:        {:.10e}\n'.format( IDs.shape[0] ) )
    f.write( '# Total number of unique mergers: {:.10e}\n\n'.format( UniqueIndices.shape[0] ) )
    for merger in UniqueIndices:
        f.write(str(merger) + '\n')
    f.write( '\nProgram end-time: {:}\n'.format(datetime.now()) )
    f.write( 'Total run-time: {:.10e} s\n'.format(TotRunTime) )

