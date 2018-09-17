#!/usr/local/bin/python3

from os.path import isfile
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams[ 'text.usetex' ] = True
plt.rcParams[ 'font.family' ] = 'serif'
plt.rcParams[ 'font.size' ] = 13
plt.rcParams[ 'axes.labelsize'] =  15

#BasePath = '/astro1/dunhamsj/'
BasePath = '/Users/sam/Research/GW/Sam/'

HistPath = BasePath + 'Histograms/'

SSMin = 26
SSMax = 27

Bins, Counts_tLb, Counts_tForm \
  = np.loadtxt( HistPath \
                  + 'HistogramDataFiles/HistFile_{:}.dat'.format(SSMin), \
                  unpack = True )
nBins        = len( Bins )
Counts_tLb   = np.zeros( (nBins), float )
Counts_tForm = np.zeros( (nBins), float )

for SS in range( SSMin, SSMax+1 ):

    DataFile = HistPath \
                   + 'HistogramDataFiles/HistFile_{:}.dat'.format(SS)

    if not isfile( DataFile ): continue

    Bins, Counts_tLbSS, Counts_tFormSS \
      = np.loadtxt( DataFile, unpack = True )

    Counts_tLb   += Counts_tLbSS
    Counts_tForm += Counts_tFormSS

fig, ax = plt.subplots()
ax.plot( Bins, Counts_tLb,   'b--', label = r'$t_{merg}$' )
ax.plot( Bins, Counts_tForm, 'r-',  label = r'$t_{form}$' )
ax.set_xlim( 0, 14 )
ax.set_xlabel( r'$t_{lb}$ [Gyr]' )
ax.set_ylabel( r'$N_{merg}$' )
ax.set_yscale( 'log' )
ax.legend()
#plt.show()
plt.savefig( HistPath + 'LookbackAndFormationTime.png' )
plt.close()
