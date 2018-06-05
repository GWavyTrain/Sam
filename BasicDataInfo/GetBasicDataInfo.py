# Program to get basic info about Mapelli data

ProgName = 'GetBasicDataInfo.py'

import numpy as np
from time import time

PATH    = 'BBHM_DataFiles_Mapelli/time_BHillustris1_'
FileOut = 'BasicDataInfo.txt'

with open( FileOut, 'w' ) as f:
    f.write( '# {:}\n\n'.format(FileOut) )
    f.write( '# These data were obtained with {:}\n\n'.format(ProgName) )
    f.write( '# Data structure:\n' )
    f.write( '# Snapshot nMergersSnapshot tLbMin tLbMax zMin zMax\n\n' )

# Define column numbers
i_z   = 2
i_tLb = 9

ssMin = 26
ssMax = 27#135

# Initialize min/max values
tLbMin = np.inf
tLbMax = -np.inf
zMin   = np.inf
zMax   = -np.inf

nMergersTot = 0

StartTime = time()
for ss in range( ssMin, ssMax+1 ):
    if ss != 53 and ss != 55:
        SnapshotData = np.loadtxt( PATH + str(ss) + '.dat' )

        tLbMin_ss = SnapshotData[:,i_tLb].min()
        tLbMax_ss = SnapshotData[:,i_tLb].max()
        tLbMin    = min( tLbMin, tLbMin_ss )
        tLbMax    = max( tLbMax, tLbMax_ss )

        zMin_ss = SnapshotData[:,i_z].min()
        zMax_ss = SnapshotData[:,i_z].max()
        zMin    = min( zMin, zMin_ss )
        zMax    = max( zMax, zMax_ss )

        nMergersSnapshot = SnapshotData.shape[0]
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
    f.write( '# Total run-time: {:.10e} s\n'.format(TotRunTime) )
