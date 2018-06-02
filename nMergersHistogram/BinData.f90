PROGRAM BinMergers

  IMPLICIT NONE

  INTEGER,  PARAMETER   :: DP = KIND( 1.d0 )
  INTEGER,  PARAMETER   :: ix = 1
  REAL(DP), PARAMETER   :: dx = 1.0d-1
  REAL(DP), ALLOCATABLE :: mergers(:,:), BinEdges(:), Bins(:), Counts(:), x(:)
  REAL(DP)              :: xMin, xMax, StartTime, StopTime, ProgStart, ProgEnd
  INTEGER               :: LOGFILE = 120, chunk
  INTEGER( KIND = 8 )   :: i, iMerger, nMergers, nBins
  CHARACTER( LEN = 6 )  :: arg
  CHARACTER( LEN = 14 ) :: LogFileName, LogFileNameOld
  CHARACTER( LEN = 24 ) :: MergersFileName, MergersFileNameOld
  CHARACTER( LEN = 25 ) :: nBinnedMergersFile

  IF ( IARGC() .NE. 1 ) THEN
    WRITE(*,'(A)') 'Proper usage: ./BinData chunk'
    WRITE(*,'(A)') 'Exiting...'
    STOP
  END IF

  CALL GETARG( 1, arg )
  READ( arg, * ) chunk

  WRITE( LogFileNameOld, '(A4,A,A4)' ) 'log_', arg, '.txt'
  LogFileName &
    =  StripWhiteSpace( LogFileNameOld, LEN( LogFileNameOld ) )

  WRITE( MergersFileNameOld, '(A13,A,A4)' ) 'mergers_chunk', arg, '.dat'
  MergersFileName &
    = StripWhiteSpace( MergersFileNameOld, LEN( MergersFileNameOld ) )

  WRITE( nBinnedMergersFile, '(A15,A,A4)' ) 'nMergersBinned_', arg, '.dat'
  nBinnedMergersFile &
    = StripWhiteSpace( nBinnedMergersFile, LEN( nBinnedMergersFile ) )

  CALL CPU_TIME( ProgStart )

  ! --- Get total number of mergers  ---
  OPEN ( LOGFILE, FILE = TRIM( LogFileName ) )
  WRITE( LOGFILE, '(A)' ) 'Getting total number of mergers...'
  CLOSE( LOGFILE )

  CALL CPU_TIME( StartTime )
  OPEN( 100, FILE = TRIM( MergersFileName ), STATUS = 'OLD' )
  nMergers = 0
  DO
    READ( 100, *, END = 10 )
    nMergers = nMergers + 1
  END DO
  10 CLOSE( 100 )
  CALL CPU_TIME( StopTime )

  OPEN ( LOGFILE, FILE = TRIM( LogFileName ), POSITION = 'APPEND' )
  WRITE( LOGFILE, '(A,I19)' ) 'nMergers: ', nMergers
  WRITE( LOGFILE, '(A,ES18.10E3,A2)' ) &
    'Time to get nMergers: ', StopTime-StartTime, ' s'
  CLOSE( LOGFILE )

  ! --- Read in mergers file and get "x" values (values to be binned) ---
  OPEN ( LOGFILE, FILE = TRIM( LogFileName ), POSITION = 'APPEND' )
  WRITE( LOGFILE, '(A)' ) 'Getting x values...'
  CLOSE( LOGFILE )

  ALLOCATE( mergers( nMergers, 2 ) )
  ALLOCATE( x( nMergers ) )
  CALL CPU_TIME( StartTime )
  OPEN( 100, FILE = TRIM( MergersFileName ), STATUS = 'OLD' )
  DO iMerger = 1, nMergers
    READ( 100, * ) mergers( iMerger, : )
    x(iMerger) = mergers( iMerger, ix )
  END DO
  CLOSE( 100 )
  CALL CPU_TIME( StopTime )

  xMin = 9.6781981984E-08
  xMax = 1.3581565412E+01

!  xMin = MINVAL( x )
!  xMax = MAXVAL( x )

  OPEN ( LOGFILE, FILE = TRIM( LogFileName ), POSITION = 'APPEND' )
  WRITE( LOGFILE, '(A,ES18.10E3,A2)' ) &
    'Time to read in x values: ', StopTime-StartTime, ' s'
  WRITE( LOGFILE, '(A,ES18.10E3)') 'xMin: ', xMin
  WRITE( LOGFILE, '(A,ES18.10E3)') 'xMax: ', xMax
  CLOSE( LOGFILE )

  ! --- Create bins/bin edges ---
  OPEN ( LOGFILE, FILE = TRIM( LogFileName ), POSITION = 'APPEND' )
  WRITE( LOGFILE, '(A)' ) 'Creating bins/bin edges...'
  CLOSE( LOGFILE )

  xMin = xMin
  xMax = xMax + dx
  nBins  = CEILING( ( xMax - xMin ) / dx )
  ALLOCATE( Bins(nBins) )
  ALLOCATE( BinEdges(nBins+1) )
  BinEdges = CreateBinEdges( nBins )
  OPEN ( LOGFILE, FILE = TRIM( LogFileName ), POSITION = 'APPEND' )
  WRITE( LOGFILE, '(A,I19)' ) 'nBins: ', nBins
  CLOSE( LOGFILE )

  ! --- Create and initialize array to hold counts ---
  ALLOCATE( Counts(nBins) )
  Counts = 0.0d0

  ! --- Populate array with counts ---
  OPEN ( LOGFILE, FILE = TRIM( LogFileName ), POSITION = 'APPEND' )
  WRITE( LOGFILE, '(A)' ) 'Populating counts array...'
  CLOSE( LOGFILE )

  OPEN ( 100 , FILE = TRIM( nBinnedMergersFile ) )
  WRITE( 100 , '(A)' ) '# bin centers, counts'
  CLOSE( 100 )
  CALL CPU_Time( StartTime )
  DO i = 1, nBins
    OPEN( 100, FILE = TRIM( nBinnedMergersFile ), POSITION = 'APPEND' )
    DO iMerger = 1, nMergers
       IF ( ( x(iMerger) >= BinEdges(i)   ) .AND. &
            ( x(iMerger) <  BinEdges(i+1) ) ) Counts(i) = Counts(i) + 1
    END DO
    Bins(i) = ( BinEdges(i) + BinEdges(i+1) ) / 2.0d0
    WRITE( 100, '(ES18.10E3,1x,ES18.10E3)' ) Bins(i), Counts(i)
    CLOSE( 100 )
  END DO
  CALL CPU_Time( StopTime )

  OPEN ( LOGFILE, FILE = TRIM( LogFileName ), POSITION = 'APPEND' )
  WRITE( LOGFILE, '(A,ES18.10E3,A2)' ) &
    'Time to populate counts array: ', StopTime-StartTime, ' s'
  CLOSE( LOGFILE )

  DEALLOCATE( Counts   )
  DEALLOCATE( BinEdges )
  DEALLOCATE( Bins     )
  DEALLOCATE( x        )
  DEALLOCATE( mergers  )

  CALL CPU_TIME( ProgEnd )

  OPEN ( LOGFILE, FILE = TRIM( LogFileName ), POSITION = 'APPEND' )
  WRITE( LOGFILE, '(A,ES18.10E3,A2)' ) &
    'Total program runtime: ', ProgEnd-ProgStart, ' s'
  CLOSE( LOGFILE )


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

  FUNCTION StripWhiteSpace( FileNameOld, LenFileNameOld ) RESULT( FileNameNew )

    INTEGER, INTENT(in)                           :: LenFileNameOld
    CHARACTER( LEN = LenFileNameOld ), INTENT(in) :: FileNameOld
    CHARACTER( LEN = LenFileNameOld )             :: FileNameNew
    INTEGER                                       :: LenFileNameNew
    INTEGER                                       :: iChar, jChar, WS=0

    ! --- Get number of whitespaces in old file ---
    iChar = 1
    DO WHILE ( iChar .LT. LenFileNameOld )
      IF ( FileNameOld( iChar:iChar ) .EQ. ' ') WS = WS + 1
      iChar = iChar + 1
    END DO

    LenFileNameNew = LenFileNameOld - WS
    jChar = 1
    DO iChar = 1, LenFileNameOld
       IF ( FileNameOld( iChar:iChar ) .NE. ' ') THEN
         FileNameNew( jChar:jChar ) = FileNameOld( iChar:iChar )
         jChar = jChar + 1
       END IF
    END DO

  END FUNCTION StripWhiteSpace

END PROGRAM BinMergers
