#!/usr/local/bin/python3

'''
Bin data in terms of lookback-time,
all times in Gyr
'''

import numpy as np
import matplotlib.pyplot as plt
from os.path import isfile

#BasePath = '/astro1/dunhamsj/'
BasePath = '/Users/sam/Research/GW/Sam/'

HistPath = BasePath + 'Histograms/'
LogFile  = HistPath + 'MakeHistogramDataFiles_Log.txt'

tLbMin = 9.6781981984e-08
tLbMax = 1.3581565412e+01

dtLb = 0.01 # 10 Myr

BinEdges = np.arange( tLbMin, tLbMax + dtLb, dtLb )
Bins     = np.copy( BinEdges[:-1] + 0.5 * dtLb )
nBins    = len( Bins )

tLbHist   = np.zeros( (nBins), np.int64 )
tFormHist = np.zeros( (nBins), np.int64 )

SSMin = 26
SSMax = 27

with open( LogFile, 'w' ) as Log:
    Log.write('')

for SS in range( SSMin, SSMax+1 ):

    DataFileName = BasePath + 'strain/hc_DataFiles/hc_{:}.dat'.format(SS)

    if not isfile( DataFileName ): continue

    with open( LogFile, 'a' ) as Log:
        Log.write( 'Snapshot {:d}\n'.format(SS) )

    tDelay, tLb = np.loadtxt( DataFileName, unpack = True, \
                                skiprows = 17, usecols = (6,7) )
    tForm = tLb + tDelay

    tLbHistSS,   bin_edges = np.histogram( tLb,   BinEdges )
    tFormHistSS, bin_edges = np.histogram( tForm, BinEdges )

    HistFile = HistPath + 'HistogramDataFiles/HistFile_{:}.dat'.format(SS)

    with open( HistFile, 'w' ) as f:
        f.write('# Histogram file, Snapshot {:}\n'.format(SS) )
        f.write('# bin centers [Gyr], tLb [Gyr], tForm [Gyr]\n' )
        for i in range( nBins ):
            f.write( '{:.16e} {:.16e} {:.16e}\n'.format(\
                      Bins[i], tLbHistSS[i], tFormHistSS[i]) )
