# Program to trim and concatenate Mapelli data

import numpy as np
from os.path import isfile
from time import time
from datetime import datetime

##### File I/O
ProgName = 'TrimAndCombineData.py'
Path    = 'BBHM_DataFiles_Mapelli/time_BHillustris1_'
FileOut = 'TrimmedAndCombinedData.dat'
LogFile = 'log.txt'

with open( FileOut, 'w' ) as f:
    f.write( '# {:}\n\n'.format(FileOut) )
    f.write( '# These data were obtained with {:}\n'.format(ProgName) )
    f.write( '# Program start-time: {:}\n\n'.format(datetime.now()) )
    f.write( '### Data structure:\n' )
    f.write( '# Merger ID\n\
# Snapshot Redshift\n\
# M1 [Msun]\n\
# M2 [Msun]\n\
# Delay Time (between formation of stellar progenitor and \
coalescence of binary) [Gyr]\n\
# BBH Merger Time (given as lookback time) [Gyr]\n\n' )

with open( LogFile, 'w' ) as f:
    f.write( '# {:}\n\n'.format(LogFile) )
    f.write( '# This file was generated with {:}\n\n'.format(ProgName) )


##### Misc. info RE: data files
iz   = 2
iID  = 5
iM1  = 6
iM2  = 7
idt  = 8
itLb = 9

ssMin = 26
ssMax = 135


##### Loop through data files, trim, and write to FileOut
StartTime = time()
for ss in range( ssMin, ssMax+1 ):
    FileIn = Path + str(ss) + '.dat'

    if not isfile( FileIn ):
        with open( LogFile, 'a' ) as f:
            f.write( FileIn + ' does not exist\n' )
        continue

    SnapshotData = np.loadtxt( FileIn )

    for merger in range( SnapshotData.shape[0] ):

        ID  = np.int64(SnapshotData[merger,iID])
        z   = SnapshotData[merger,iz]
        M1  = SnapshotData[merger,iM1]
        M2  = SnapshotData[merger,iM2]
        dt  = SnapshotData[merger,idt]
        tLb = SnapshotData[merger,itLb]

        with open( FileOut, 'a' ) as f:
            f.write( '{:d} {:10e} {:10e} {:10e} {:10e} {:10e}\n'.format(\
                    ID, z, M1, M2, dt, tLb ) )

EndTime = time()

TotRunTime = EndTime - StartTime

with open( FileOut, 'a' ) as f:
    f.write( '\n# Program end-time: {:}\n'.format(datetime.now()) )
    f.write( '# Total run-time: {:.10e} s\n'.format(TotRunTime) )
