#!/usr/local/bin/python3

import numpy as np
import matplotlib.pyplot as plt

#BasePath = '/astro1/dunhamsj/'
BasePath = '/Users/sam/Research/GW/Sam/'

tLbMin = 9.6781981984e-08
tLbMax = 1.3581565412e+01

dtLb = 0.01 # 10 Myr

BinEdges = np.arange( tLbMin, tLbMax + dtLb, dtLb )
Bins = np.copy( BinEdges[:-1] + 0.5 * dtLb )
nBins = len( Bins )

tLbHist   = np.zeros( (nBins), np.int64 )
tFormHist = np.zeros( (nBins), np.int64 )

SSMin = 26
SSMax = 50

SS = SSMin
while SS <= SSMax:

    tDelay, tLb = np.loadtxt( BasePath + 'strain/hc_DataFiles/hc_26.dat', \
                                unpack = True, skiprows = 17, usecols = (6,7) )
    tForm = tLb + tDelay

    tLbHistSS,   bin_edges = np.histogram( tLb,   BinEdges )
    tFormHistSS, bin_edges = np.histogram( tForm, BinEdges )

    tLbHist   += tLbHistSS
    tFormHist += tFormHistSS

    SS += 1

fig, ax = plt.subplots()
ax.step( Bins, tLbHist,   'b-', label = 'Merger time' )
ax.step( Bins, tFormHist, 'r-', label = 'Formation time' )
ax.set_xlabel( 'Lookback-Time [Gyr]' )
ax.set_ylabel( 'Counts' )
ax.set_yscale( 'log' )
ax.legend()
plt.savefig( BasePath + 'MergerAndFormationTime.png' )
plt.close()

fig, ax = plt.subplots()
ax.step( Bins, tLbHist, 'b-', label = 'Merger time' )
ax.set_xlabel( 'Lookback-Time [Gyr]' )
ax.set_ylabel( 'Counts' )
ax.set_yscale( 'log' )
ax.legend()
plt.savefig( BasePath + 'MergerTime.png' )
plt.close()

fig, ax = plt.subplots()
ax.step( Bins, tFormHist, 'r-', label = 'Formation time' )
ax.set_xlabel( 'Lookback-Time [Gyr]' )
ax.set_ylabel( 'Counts' )
ax.set_yscale( 'log' )
ax.legend()
plt.savefig( BasePath + 'FormationTime.png' )
plt.close()
