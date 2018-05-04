PROGRAM nMergersHistogram

  IMPLICIT NONE

  INTEGER, PARAMETER    :: DP = KIND( 1.d0 )
  INTEGER,  PARAMETER   :: iz = 3, itLB = 10
  REAL(DP), ALLOCATABLE :: IllustrisData(:,:)
  REAL(DP)              :: z, tLB
  INTEGER( KIND = 8 )   :: nLinesIllustris, iMerger
  INTEGER               :: Nss, SS(2), i
  CHARACTER( LEN = 28 ) :: FILEIN
  CHARACTER( LEN = 19 ) :: FILEOUT
  CHARACTER( LEN = 11 ) :: FMTIN
  CHARACTER( LEN = 3  ) :: arg

  ! --- Get snapshot numbers ---
  IF ( IARGC() .NE. 2 ) THEN
    WRITE(*,'(A)') 'Proper usage: ./strain SSL SSU'
    WRITE(*,'(A)') 'Exiting...'
    STOP
  END IF

  DO i = 1, IARGC()
    CALL GETARG( i, arg )
    READ( arg, * ) SS(i)
  END DO

  WRITE( FILEOUT, '(A8,I3.3,A1,I3.3,A4)' ) 'mergers_', SS(1), '-', SS(2), '.dat'

  ! --- Create file for storing mergers ---
  OPEN ( 100, FILE = TRIM( FILEOUT ) )
  CLOSE( 100 )

  ! --- Loop through Illustris snapshots ---
  DO Nss = SS(1), SS(2)

    IF ( ( Nss .NE. 53 ) .AND. ( Nss .NE. 55 ) ) THEN

      ! --- Get filenames ---
      IF ( Nss < 100 ) THEN
        FMTIN = '(A44,I2,A4)'
      ELSE
        FMTIN = '(A44,I3,A4)'
      END IF

      WRITE( FILEIN, FMTIN  ) &
        '../BBHM_DataFiles_Mapelli/time_BHillustris1_', Nss, '.dat'

      ! --- Get number of lines (mergers) in Illustris data file ---
      nLinesIllustris = 0
      OPEN( 100, FILE = TRIM( FILEIN ) )
      DO
        READ( 100, *, END = 10 )
        nLinesIllustris = nLinesIllustris + 1
      END DO
      10 CLOSE( 100 )

      ! --- Read in Illustris data and get redshifts and lookback times ---
      ALLOCATE( IllustrisData( nLinesIllustris, 10 ) )
      OPEN( 100, FILE = TRIM( FILEIN ), STATUS = 'OLD' )
      DO iMerger = 1, nLinesIllustris
        READ( 100 , * ) IllustrisData( iMerger, : )

        ! --- Get data for this merger ---
        tLB = IllustrisData( iMerger, itLB ) ! [ Gyr ]
        z   = IllustrisData( iMerger, iz   )

        OPEN ( 101, FILE = TRIM( FILEOUT ), POSITION = 'APPEND' )
        WRITE( 101, '(ES18.10E3,1x,ES18.10E3)' ) tLB, z
        CLOSE( 101 )
      END DO
      CLOSE( 100 )
      DEALLOCATE( IllustrisData )

    END IF

  END DO

END PROGRAM nMergersHistogram
