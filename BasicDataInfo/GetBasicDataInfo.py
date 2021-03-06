# Program to get basic info about Mapelli data

ProgName = 'GetBasicDataInfo.py'

import numpy as np
from os.path import isfile
from time import time
from datetime import datetime

Path    = 'BBHM_DataFiles_Mapelli/time_BHillustris1_'
FileOut = 'BasicDataInfo.txt'
LogFile = 'log.txt'

with open( FileOut, 'w' ) as f:
    f.write( '# {:}\n\n'.format(FileOut) )
    f.write( '# These data were obtained with {:}\n'.format(ProgName) )
    f.write( '# Program start-time: {:}\n\n'.format(datetime.now()) )
    f.write( '# Data structure:\n' )
    f.write( '# Snapshot nMergersSnapshot tLbMin tLbMax zMin zMax\n\n' )

with open( LogFile, 'w' ) as f:
    f.write( '# {:}\n\n'.format(LogFile) )
    f.write( '# This file was generated with {:}\n\n'.format(ProgName) )

# Define column numbers
i_z   = 2
i_tLb = 9

ssMin = 26
ssMax = 135

# Initialize min/max values
tLbMin = np.inf
tLbMax = -np.inf
zMin   = np.inf
zMax   = -np.inf

nMergersTot = 0.0

StartTime = time()
for ss in range( ssMin, ssMax+1 ):
    FileIn = Path + str(ss) + '.dat'

    if not isfile( FileIn ):
        with open( LogFile, 'a' ) as f:
            f.write( FileIn + ' does not exist\n' )
        continue

    SnapshotData = np.loadtxt( FileIn )

    tLbMin_ss = SnapshotData[:,i_tLb].min()
    tLbMax_ss = SnapshotData[:,i_tLb].max()
    tLbMin    = min( tLbMin, tLbMin_ss )
    tLbMax    = max( tLbMax, tLbMax_ss )

    zMin_ss = SnapshotData[:,i_z].min()
    zMax_ss = SnapshotData[:,i_z].max()
    zMin    = min( zMin, zMin_ss )
    zMax    = max( zMax, zMax_ss )

    nMergersSnapshot = float( SnapshotData.shape[0] )
    nMergersTot += nMergersSnapshot

    with open( FileOut, 'a' ) as f:
        f.write( '{:d} {:.10e} {:.10e} {:.10e} {:.10e} {:.10e}\n'.format(\
                 ss, nMergersSnapshot, tLbMin_ss, tLbMax_ss, zMin_ss, zMax_ss) )

EndTime = time()

TotRunTime = EndTime - StartTime

with open( FileOut, 'a' ) as f:
    f.write( '\n# Totals:\n' )
    f.write( '{:d} {:.10e} {:.10e} {:.10e} {:.10e} {:.10e}\n\n'.format(\
              999, nMergersTot, tLbMin, tLbMax, zMin, zMax) )
    f.write( '# Program end-time: {:}\n'.format(datetime.now()) )
    f.write( '# Total run-time: {:.10e} s\n'.format(TotRunTime) )
