This folder is used for binning the data in terms of
lookback time or snapshot redshift

mergers.dat: File with two columns:
  Column 1: Lookback time to merger [Gyr]
  Column 2: Snapshot redshift when merger was detected

mergers.dat was broken into chunks using 'MakeChunks.sh' which live in
ChunkedMergers_LookbackTime -- and/or -- ChunkedMergers_SnapshotRedshift

The chunks are binned with
BinData.exe (compiled from BinData.f90) and BinData.sh

Basic info about mergers.dat:
-----------------------------
Total number of mergers: 7.722072936e9

From 'FindMin.py':
-----------------
Minimum lookback time: 9.6781981984E-08 Gyr
Maximum lookback time: 1.3581565412E+01 Gyr

Minimum snapshot redshift: 7.2033238735E-03
Maximum snapshot redshift: 1.3592065412E+01

