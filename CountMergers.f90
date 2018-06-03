PROGRAM CountMergers

  IMPLICIT NONE

  INTEGER,  PARAMETER   :: DP = KIND( 1.d0 )
  INTEGER(KIND=8)       :: nMergersSnapshot, nMergersTot, iMerger
  REAL(DP), ALLOCATABLE :: IllustrisDataFile(:,:)

  ! --- For file I/O ---
  INTEGER           :: Snapshot = 100, Nss
  CHARACTER(LEN=13) :: PATH
  CHARACTER(LEN=25) :: FILEIN
  CHARACTER(LEN=11) :: FMTIN
  PATH = 'BBHM_DataFiles_Mapelli/'


  DO Nss = 26, 135
     IF( Nss .NE. 53 .AND. Nss .NE. 55 ) THEN

        IF( Nss .LT. 100 ) THEN
           FMTIN = '(A18,I2,A4)'
        ELSE
           FMTIN = '(A18,I3,A4)'
        END IF

        FILEIN = 'time_BHillustris1_', Nss, '.dat'
        WRITE(*,*) FILEIN
!!$        ! --- Get number of lines (mergers) in snapshot---
!!$        nLinesIllustris = 0
!!$        OPEN( Snapshot, FILE = TRIM( FILEIN ) )
!!$        DO
!!$          READ( IllustrisDataFile, *, END = 11 )
!!$          nLinesIllustris = nLinesIllustris + 1
!!$        END DO
!!$        11 CLOSE( IllustrisDataFile )

     END IF ! --- Don't check file 53 or file 55

     
  END DO ! --- Loop over snapshots
  
END PROGRAM CountMergers
