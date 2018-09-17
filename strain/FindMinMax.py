#!/usr/local/bin/python3

from os.path import isfile
from time import time
import numpy as np

#BasePath = '/astro1/dunhamsj/'
BasePath = '/Users/sam/Research/GW/Sam/'

StrainPath = BasePath + 'strain/hc_DataFiles/'
OutputFile = BasePath + 'strain/MinMax.dat'
OutputFMT = '{:} {:.16e} {:.16e} {:.16e} {:.16e} {:.16e} {:.16e} \
{:.16e} {:.16e} {:.16e} {:.16e}\n'

SSMin = 26
SSMax = 27

# --- Global min/max ---
zSnapshotMinG = tDelayMinG = tLbMinG = zMergeMinG = rMinG = +np.inf
zSnapshotMaxG = tDelayMaxG = tLbMaxG = zMergeMaxG = rMaxG = -np.inf

with open( OutputFile, 'w' ) as f:
    f.write( '# Min/Max values of parameters\n' )
    f.write( \
      '# Snapshot zSnapshot Min/Max, tDelay Min/Max [Gyr], tLb Min/Max [Gyr], \
zMerge Min/Max, r Min/Max [Mpc]\n' )

LoopStartTime = time()
for SS in range( SSMin, SSMax+1 ):
    DataFileName = StrainPath + 'hc_{:}.dat'.format(SS)

    if not isfile( DataFileName ): continue

    # --- Local min/max ---
    zSnapshotMin = tDelayMin = tLbMin = zMergeMin = rMin = +np.inf
    zSnapshotMax = tDelayMax = tLbMax = zMergeMax = rMax = -np.inf
    
    zSnapshot, tDelay, tLb, zMerge, r \
      = np.loadtxt( DataFileName, unpack = True, \
                      usecols = (2,6,7,8,9), skiprows = 17 )

    zSnapshotMin = min( zSnapshotMin, zSnapshot.min() )
    zSnapshotMax = max( zSnapshotMax, zSnapshot.max() )
    tDelayMin    = min( tDelayMin,    tDelay.min()    )
    tDelayMax    = max( tDelayMax,    tDelay.max()    )
    tLbMin       = min( tLbMin,       tLb.min()       )
    tLbMax       = max( tLbMax,       tLb.max()       )
    zMergeMin    = min( zMergeMin,    zMerge.min()    )
    zMergeMax    = max( zMergeMax,    zMerge.max()    )
    rMin         = min( rMin,         r.min()         )
    rMax         = max( rMax,         r.max()         )

    with open( OutputFile, 'a' ) as f:
        f.write( OutputFMT.format( \
          SS, zSnapshotMin, zSnapshotMax, tDelayMin, tDelayMax, tLbMin, \
            tLbMax, zMergeMin, zMergeMax, rMin, rMax ) )

    zSnapshotMinG = min( zSnapshotMin, zSnapshotMinG )
    zSnapshotMaxG = max( zSnapshotMax, zSnapshotMaxG )
    tDelayMinG    = min( tDelayMin,    tDelayMinG    )
    tDelayMaxG    = max( tDelayMax,    tDelayMaxG    )
    tLbMinG       = min( tLbMin,       tLbMinG       )
    tLbMaxG       = max( tLbMax,       tLbMaxG       )
    zMergeMinG    = min( zMergeMin,    zMergeMinG    )
    zMergeMaxG    = max( zMergeMax,    zMergeMaxG    )
    rMinG         = min( rMin,         rMinG         )
    rMaxG         = max( rMax,         rMaxG         )

LoopStopTime = time() - LoopStartTime

with open( OutputFile, 'a' ) as f:
    f.write( '\n' )
    f.write( OutputFMT.format( \
      999, zSnapshotMinG, zSnapshotMaxG, tDelayMinG, tDelayMaxG, tLbMinG, \
        tLbMaxG, zMergeMinG, zMergeMaxG, rMinG, rMaxG ) )
    f.write( '# Total run-time: {:.16e} s\n'.format(LoopStopTime) )
