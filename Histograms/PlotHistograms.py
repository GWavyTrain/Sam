#!/usr/local/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from os.path import isfile

#BasePath = '/astro1/dunhamsj/'
BasePath = '/Users/sam/Research/GW/Sam/'

HistPath = BasePath + 'Histograms/'

SSMin = 26
SSMax = 27

Bins_z, Counts_z \
  = np.loadtxt( HistPath \
                  + 'HistogramDataFiles/HistFile_z_{:}.dat'.format(SSMin), \
                  unpack = True )
nBins_z  = len(Bins_z)
Counts_z = np.zeros( (nBins_z), float )

Bins_t, Counts_tLb, Counts_tForm \
  = np.loadtxt( HistPath \
                  + 'HistogramDataFiles/HistFile_t_{:}.dat'.format(SSMin), \
                  unpack = True )
nBins_t      = len( Bins_t )
Counts_tLb   = np.zeros( (nBins_t), float )
Counts_tForm = np.zeros( (nBins_t), float )

for SS in range( SSMin, SSMax+1 ):

    DataFile_t = HistPath \
                   + 'HistogramDataFiles/HistFile_t_{:}.dat'.format(SS)

    DataFile_z = HistPath \
                   + 'HistogramDataFiles/HistFile_z_{:}.dat'.format(SS)

    if not isfile( DataFile_t ): continue

    Bins_z, Counts_zSS \
      = np.loadtxt( DataFile_z, unpack = True )

    Counts_z += Counts_zSS

    Bins_t, Counts_tLbSS, Counts_tFormSS \
      = np.loadtxt( DataFile_t, unpack = True )

    Counts_tLb   += Counts_tLbSS
    Counts_tForm += Counts_tFormSS

fig, ax = plt.subplots()
ax.step( Bins_z, Counts_z, 'k-' )
ax.set_xlabel( 'Snapshot Redshift' )
ax.set_ylabel( 'Counts' )
ax.set_yscale( 'log' )
plt.savefig( HistPath + 'SnapshotRedshift.png' )
plt.close()

fig, ax = plt.subplots()
ax.step( Bins_t, Counts_tLb, 'b-', label = 'Merger Lookback Time' )
ax.step( Bins_t, Counts_tForm, 'r-', label = 'Formation Time' )
ax.set_xlabel( 'Lookback Time [Gyr]' )
ax.set_ylabel( 'Counts' )
ax.set_yscale( 'log' )
ax.legend()
plt.savefig( HistPath + 'LookbackAndFormationTime.png' )
plt.close()
