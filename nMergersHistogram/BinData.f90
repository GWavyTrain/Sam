PROGRAM BinMergers

  IMPLICIT NONE

  INTEGER,  PARAMETER   :: DP = KIND( 1.d0 )
  INTEGER,  PARAMETER   :: ix = 1
  REAL(DP), PARAMETER   :: dx = 1.0d-1
  REAL(DP), ALLOCATABLE :: mergers(:,:), BinEdges(:), Bins(:), Counts(:), x(:)
  REAL(DP)              :: xMin, xMax, StartTime, StopTime
  INTEGER               :: i, LOGFILE = 120
  INTEGER( KIND = 8 )   :: iMerger, nMergers, nBins

  ! --- Get total number of mergers  ---
  OPEN ( LOGFILE, FILE = 'log.txt' )
  WRITE( LOGFILE, '(A)' ) 'Getting total number of mergers...'
  CLOSE( LOGFILE )

  OPEN( 100, FILE = 'mergers.dat', STATUS = 'OLD' )
  nMergers = 0
  CALL CPU_TIME( StartTime )
  DO
    READ( 100, *, END = 10 )
    nMergers = nMergers + 1
  END DO
  CALL CPU_TIME( StopTime )
  10 CLOSE( 100 )

  OPEN ( LOGFILE, FILE = 'log.txt', POSITION = 'APPEND' )
  WRITE( LOGFILE, '(A,I19)' ) 'nMergers: ', nMergers
  WRITE( LOGFILE, '(A,ES18.10E3)' ) 'Time to get nMergers: ', StopTime-StartTime
  CLOSE( LOGFILE )

  ! --- Read in mergers file and get "x" values (values to be binned) ---
  OPEN ( LOGFILE, FILE = 'log.txt', POSITION = 'APPEND' )
  WRITE( LOGFILE, '(A)' ) 'Getting x values...'
  CLOSE( LOGFILE )

  ALLOCATE( mergers( nMergers, 2 ) )
  ALLOCATE( x( nMergers ) )
  OPEN( 100, FILE = 'mergers.dat', STATUS = 'OLD' )
  CALL CPU_TIME( StartTime )
  DO iMerger = 1, nMergers
    READ( 100, * ) mergers( iMerger, : )
    x(iMerger) = mergers( iMerger, ix )
  END DO
  CALL CPU_TIME( StopTime )
  CLOSE( 100 )
  xMin = MINVAL( x )
  xMax = MAXVAL( x )

  OPEN ( LOGFILE, FILE = 'log.txt', POSITION = 'APPEND' )
  WRITE( LOGFILE, '(A,ES18.10E3)' ) 'Time to read in x values: ', StopTime-StartTime
  WRITE( LOGFILE, '(A,ES18.10E3)') 'xMin: ', xMin
  WRITE( LOGFILE, '(A,ES18.10E3)') 'xMax: ', xMax
  CLOSE( LOGFILE )

  ! --- Create bins/bin edges ---
  OPEN ( LOGFILE, FILE = 'log.txt', POSITION = 'APPEND' )
  WRITE( LOGFILE, '(A)' ) 'Creating bins/bin edges...'
  CLOSE( LOGFILE )

  xMin = xMin
  xMax = xMax + dx
  nBins  = CEILING( ( xMax - xMin ) / dx )
  ALLOCATE( Bins(nBins) )
  ALLOCATE( BinEdges(nBins+1) )
  BinEdges = CreateBinEdges( nBins )
  OPEN ( LOGFILE, FILE = 'log.txt', POSITION = 'APPEND' )
  WRITE( LOGFILE, '(A,I19)' ) 'nBins: ', nBins
  CLOSE( LOGFILE )

  ! --- Create and initialize array to hold counts ---
  ALLOCATE( Counts(nBins) )
  Counts = 0.0d0

  ! --- Populate array with counts ---
  OPEN ( LOGFILE, FILE = 'log.txt', POSITION = 'APPEND' )
  WRITE( LOGFILE, '(A)' ) 'Populating counts array...'
  CLOSE( LOGFILE )

  OPEN ( 100 , FILE = 'nMergersBinned.dat' )
  WRITE( 100 , '(A)' ) '# bin centers, counts'
  CLOSE( 100 )
  CALL CPU_Time( StartTime )
  DO i = 1, nBins
    OPEN( 100, FILE = 'nMergersBinned.dat', POSITION = 'APPEND' )
    DO iMerger = 1, nMergers
       IF ( ( x(iMerger) >= BinEdges(i)   ) .AND. &
            ( x(iMerger) <  BinEdges(i+1) ) ) Counts(i) = Counts(i) + 1
    END DO
    Bins(i) = ( BinEdges(i) + BinEdges(i+1) ) / 2.0d0
    WRITE( 100, '(ES18.10E3,1x,ES18.10E3)' ) Bins(i), Counts(i)
    CLOSE( 100 )
  END DO
  CALL CPU_Time( StopTime )

  OPEN ( LOGFILE, FILE = 'log.txt', POSITION = 'APPEND' )
  WRITE( LOFILE, '(A,ES18.10E3)' ) 'Time to populate counts array: ', StopTime-StartTime
  CLOSE( LOGFILE )

  DEALLOCATE( Counts   )
  DEALLOCATE( BinEdges )
  DEALLOCATE( Bins     )
  DEALLOCATE( x        )
  DEALLOCATE( mergers  )

CONTAINS

  FUNCTION CreateBinEdges( nBins ) RESULT( BinEdges )

    INTEGER( KIND = 8 ), INTENT(in) :: nBins
    REAL(DP)                        :: BinEdges(nBins+1)
    INTEGER( KIND = 8 )             :: i

    DO i = 1, nBins + 1
      BinEdges(i) = xMin + ( i - 1 ) * dx
    END DO

    RETURN
  END FUNCTION CreateBinEdges

END PROGRAM BinMergers
