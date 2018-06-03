import numpy as np

tLBMin = 1.0e20
tLBMax = -1
zMin   = 1.0e20
zMax   = -1

with open('MinMax.txt', 'w') as fnew:
    fnew.write('MinMax.txt\n\n')

for chunk in range( 79 ):
    chunk = str( chunk )
    tLB, z = np.loadtxt( 'mergers_chunk' + chunk + '.dat', unpack = True )

    tLBMin = min( tLBMin, tLB.min() )
    tLBMax = max( tLBMax, tLB.max() )
    zMin   = min( zMin, z.min() )
    zMax   = max( zMax, z.max() )

    with open('MinMax.txt', 'a') as f:
        f.write( 'Chunk {:d}'.format(int(chunk)) )
        f.write( '\n{:.10E} {:.10E} {:.10E} {:.10E}\n\n'.format(tLBMin,tLBMax,zMin,zMax) )

with open('MinMax.txt', 'a') as f:
    f.write( '\n----------------------------------' )
    f.write( '\ntLBMin: {:.10E}\ntLBMax: {:.10E}\nzMin:   {:.10E}\nzMax:   {:.10E}\n'.format(\
               tLBMin,tLBMax,zMin,zMax))
