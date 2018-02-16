PROGRAM GWstrainFromBHBmergers

  IMPLICIT NONE

  INTEGER  , PARAMETER   :: DP = KIND( 1.d0 ) , Nq = 1000000
  INTEGER  , PARAMETER   :: Nf = 900 , Nc = 22
  INTEGER  , PARAMETER   :: iz = 3 , iM1 = 7 , iM2 = 8
  REAL(DP) , PARAMETER   :: PI = ACOS( -1.0d0 )
  REAL(DP) , PARAMETER   :: OMEGA_M = 0.3d0 , OMEGA_L = 0.7d0
  REAL(DP) , PARAMETER   :: c = 3.0d10 , G = 6.67d-8 , H0 = 70.4d0 / 3.086d19
  REAL(DP) , PARAMETER   :: Msun = 2.0d33
  REAL(DP) , ALLOCATABLE :: IllustrisData(:,:) , R_ISCO(:) , f_ISCO(:)
  REAL(DP)               :: D , LISA(Nf,2) , hc , hc_min , hc_max , z
  INTEGER                :: nLinesIllustris , Nz , i , j
  CHARACTER( len = 25 )  :: FILEIN
  CHARACTER( len = 17 )  :: FILEOUT
  CHARACTER( len = 11 )  :: FMTIN , FMTOUT
  CHARACTER( len = 34 )  :: FMT

  ! --- Read in frequencies from LISA sensitivity curve ---
  OPEN( 100 , FILE = 'LISA_sensitivity.dat' )

  ! --- Loop through comments
  DO i = 1 , Nc
    READ( 100 , * )
  END DO
 
  DO i = 1 , Nf
    READ( 100 , * ) LISA( i , : )
  END DO
 
  CLOSE( 100 )

  ! --- Loop through Illustris snapshots ---
  DO Nz = 26 , 135

    IF ( ( Nz .NE. 53 ) .AND. ( Nz .NE. 55 ) ) THEN
      OPEN( 200 , FILE = 'Nz.dat' )
        WRITE( 200 , '(I3)' ) Nz
      CLOSE( 200 )

      ! --- Get filenames ---
      IF ( Nz < 100 ) THEN
        FMTIN  = '(A18,I2,A4)'
        FMTOUT = '(A10,I2,A4)'
      ELSE
        FMTIN  = '(A18,I3,A4)'
        FMTOUT = '(A10,I3,A4)'
      END IF

      WRITE( FILEIN  , FMTIN  ) 'time_BHillustris1_' , Nz , '.dat'
      WRITE( FILEOUT , FMTOUT ) 'GW_strain_'         , Nz , '.dat'

      ! --- Get number of lines (mergers) in Illustris data file ---
      nLinesIllustris = 0
      OPEN( 101 , FILE = TRIM( FILEIN ) )
      DO
        READ( 101 , * , END = 11 )
        nLinesIllustris = nLinesIllustris + 1
      END DO
      11 CLOSE( 101 )

      ! --- Read in Illustris data file ---
      ALLOCATE( IllustrisData( nLinesIllustris , 10 ) )
      OPEN( 102 , FILE = TRIM( FILEIN ) )
      DO i = 1 , nLinesIllustris
         READ( 102 , * ) IllustrisData( i , : )
      END DO
      CLOSE( 102 )

      ! --- Calculate comoving distance ---
      D = ComputeComovingDistance( IllustrisData( 1 , iz ) )

      ! --- Calculate the frequency at the ISCO
      ALLOCATE( R_ISCO( nLinesIllustris ) )
      ALLOCATE( f_ISCO( nLinesIllustris ) )
      R_ISCO = 6.0d0 * G / c**2 * Msun &
                 * MIN( IllustrisData( : , iM1 ) , IllustrisData( : , iM2 ) )
      f_ISCO = 1.0d0 / ( 2.0d0 * PI )                          &
                 * SQRT( G * Msun * ( IllustrisData( : , iM1 ) &
                   + IllustrisData( : , iM2 ) ) / R_ISCO**3 )
  
      ! --- Loop through snapshots and calculate strain ---
      OPEN ( 103 , FILE = TRIM( FILEOUT ) )
      WRITE( 103 , '(A23)' ) '# frequency, hc, M1, M2'
      WRITE( 103 , '(F13.10,1x,I6,1x,I1,1x,I1 )' ) &
        IllustrisData( 1 , iz ) , INT( D / 3.086d24 ) , 0 , 0
      FMT = '(E11.5,1x,E11.5,1x,F7.4,1x,F7.4)'
      DO i = 1 , nLinesIllustris
        j = 1
        DO WHILE ( LISA( j , 1 ) < f_ISCO( i ) )
          CALL Strain( IllustrisData( i , iM1 ) , &
                         IllustrisData( i , iM2 ) , LISA( j , 1 ) , hc )
          WRITE( 103 , FMT ) LISA( i , 1 ) , hc , IllustrisData( i , iM1 ) , &
                             IllustrisData( i , iM2 )
          j = j + 1
        END DO
      END DO
      CLOSE( 103 )

      DEALLOCATE( R_ISCO)
      DEALLOCATE( f_ISCO )
      DEALLOCATE( IllustrisData )

    END IF

  END DO

CONTAINS

  ! --- E(z) in cosmological distance calculation ---
  PURE FUNCTION E( z )

    REAL(DP) , INTENT(in) :: z
    REAL(DP)              :: E

    E = SQRT( OMEGA_M * ( 1.0d0 + z )**3 + OMEGA_L )

  END FUNCTION E

  FUNCTION ComputeComovingDistance( z ) RESULT( D_C )

    REAL(DP) , INTENT(in) :: z
    REAL(DP)              :: D_C , dz

    ! --- Integrate with Trapezoidal rule ---
    D_C = 0.0d0
    dz = z / Nq
    
    DO i = 1 , Nq - 1
       D_C = D_C + 1.0d0 / E( i * dz )
    END DO

    D_C = c / H0 * dz / 2.0d0 * ( E( 0.0d0 ) + 2.0d0 * D_C + E( z ) )

    RETURN
  END FUNCTION ComputeComovingDistance
  
  ! --- Characteristic strain (dimensionless) from Sesana (2016) ---
  SUBROUTINE Strain( M1 , M2 , f , hc )

    REAL(DP) , INTENT(in)  :: M1 , M2 , f
    REAL(DP) , INTENT(out) :: hc
    REAL(DP)               :: Mc

    ! --- Compute chirp mass (in source frame) ---
    Mc = ( M1 * M2 )**( 3.0d0 / 5.0d0 ) / ( M1 + M2 )**( 1.0d0 / 5.0d0 ) * Msun

    hc = 1.0d0 / ( PI * D ) &
           * ( 2.0d0 * PI / ( 3.0d0 * c**3 ) )**( 1.0d0 / 2.0d0 ) &
             * ( G * Mc )**( 5.0d0 / 6.0d0 ) * ( PI * f )**( -1.0d0 / 6.0d0 )

    RETURN
  END SUBROUTINE Strain
  
END PROGRAM GWstrainFromBHBmergers
