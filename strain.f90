PROGRAM GWstrainFromBHBmergers

  IMPLICIT NONE

  INTEGER  , PARAMETER   :: DP = KIND( 1.d0 ) , Nq = 1000000
  INTEGER  , PARAMETER   :: Nf = 400 , Nc = 7
  INTEGER  , PARAMETER   :: iID = 6 , iz = 3 , iM1 = 7 , iM2 = 8
  REAL(DP) , PARAMETER   :: PI = ACOS( -1.0d0 )
  REAL(DP) , PARAMETER   :: OMEGA_M = 0.3d0 , OMEGA_L = 0.7d0
  REAL(DP) , PARAMETER   :: c = 3.0d10 , G = 6.67d-8 , H0 = 70.4d0 / 3.086d19
  REAL(DP) , PARAMETER   :: Msun = 2.0d33
  REAL(DP) , ALLOCATABLE :: IllustrisData(:,:)
  REAL(DP)               :: r , LISA(Nf,2) , hc , z
  REAL(DP)               :: f_ISCO , M1 , M2 , Mtot , mu
  INTEGER                :: nLinesIllustris , Nz , i , j
  INTEGER*8              :: MergerID
  CHARACTER( len = 25 )  :: FILEIN
  CHARACTER( len = 11 )  :: FMTIN

  ! --- Read in frequencies from LISA sensitivity curve ---
  OPEN( 100 , FILE = 'LISA_sensitivity.dat' )
  ! --- Loop through comments ---
  DO i = 1 , Nc
    READ( 100 , * )
  END DO
  DO i = 1 , Nf
    READ( 100 , * ) LISA( i , : )
  END DO
  CLOSE( 100 )

  ! --- Create file for storing strains (first row will hold frequencies) ---
  OPEN ( 101 , FILE = 'hc.dat' )
  WRITE( 101 , '(A25)' ) '# ID, M1, M2, z, r, hc(f)'
  WRITE( 101 , '(I6,1x,F7.4,1x,F7.4,1x,F13.10,1x,E16.10,400E13.6)' ) &
    0 , 0.0d0 , 0.0d0 , 0.0d0 , 0.0d0 , LISA( : , 1 )

  ! --- Loop through Illustris snapshots ---
  DO Nz = 26 , 27!135

    IF ( ( Nz .NE. 53 ) .AND. ( Nz .NE. 55 ) ) THEN
      !OPEN ( 102 , FILE = 'Nz.dat' )
      !WRITE( 102 , '(I3)' ) Nz
      !CLOSE( 102 )

      ! --- Get filenames
      IF ( Nz < 100 ) THEN
        FMTIN = '(A18,I2,A4)'
      ELSE
        FMTIN = '(A18,I3,A4)'
      END IF

      WRITE( FILEIN  , FMTIN  ) 'time_BHillustris1_' , Nz , '.dat'

      ! --- Get number of lines (mergers) in Illustris data file ---
      nLinesIllustris = 0
      OPEN( 103 , FILE = TRIM( FILEIN ) )
      DO
        READ( 103 , * , END = 13 )
        nLinesIllustris = nLinesIllustris + 1
      END DO
      13 CLOSE( 103 )

      ! --- Read in Illustris data ---
      ALLOCATE( IllustrisData( nLinesIllustris , 10 ) )
      OPEN( 104 , FILE = TRIM( FILEIN ) )
      DO i = 1 , nLinesIllustris
         READ( 104 , * ) IllustrisData( i , : )
      END DO
      CLOSE( 104 )

      ! --- Get redshift ---
      z = IllustrisData( 1 , iz )

      ! --- Calculate comoving distance (in cm) ---
      r = ComputeComovingDistance( z )

      ! --- Loop through snapshots and calculate strain ---
      DO i = 1 , nLinesIllustris

        ! --- Get data for this merger ---
        MergerID = IllustrisData( i , iID )
        M1       = IllustrisData( i , iM1 )
        M2       = IllustrisData( i , iM2 )
        
        WRITE( 101 , '(I6,1x,F7.4,1x,F7.4,1x,F13.10,1x,E16.10)' , &
          ADVANCE = 'NO' ) MergerID , M1 , M2 , z , r / 3.086d24

        ! --- Calculate the frequency at the ISCO   ---
        ! --- using two-body decomposition          ---
        ! --- (total mass Mtot and reduced mass mu) ---
        Mtot   = ( M1 + M2 ) * Msun
        mu     = M1 * M2 / ( M1 + M2 ) * Msun
        f_ISCO = c**3 / ( 6.0d0**( 3.0d0 / 2.0d0 ) * PI * G * ( 1.0d0 + z ) ) &
                   * SQRT( ( Mtot + mu ) / Mtot**3 )

        ! --- Loop through frequencies until f_ISCO ---
        j = 1
        DO WHILE ( ( LISA( j , 1 ) < f_ISCO ) .AND. ( j < Nf + 1 ) )
          CALL Strain( M1 , M2 , LISA( j , 1 ) , hc )
          WRITE( 101 , '(E13.6)' , ADVANCE = 'NO' ) hc
          j = j + 1
        END DO

        ! --- Fill in missing frequencies with 0.0d0 ---
        DO WHILE ( j < Nf + 1 )
          WRITE( 101 , '(E13.6)' , ADVANCE = 'NO' ) 0.0d0
          j = j + 1
        END DO

        WRITE( 101 , * )
        
      END DO

      DEALLOCATE( IllustrisData )

    END IF

  END DO

  CLOSE( 101 )
  
CONTAINS

  ! --- E(z) in cosmological distance calculation ---
  PURE FUNCTION E( z )

    REAL(DP) , INTENT(in) :: z
    REAL(DP)              :: E

    E = SQRT( OMEGA_M * ( 1.0d0 + z )**3 + OMEGA_L )

  END FUNCTION E

  FUNCTION ComputeComovingDistance( z ) RESULT( r )

    REAL(DP) , INTENT(in) :: z
    REAL(DP)              :: r , dz

    ! --- Integrate with Trapezoidal rule ---
    r  = 0.0d0
    dz = z / Nq
    
    DO i = 1 , Nq - 1
       r = r + 1.0d0 / E( i * dz )
    END DO

    r = c / H0 * dz / 2.0d0 * ( E( 0.0d0 ) + 2.0d0 * r + E( z ) )

    RETURN
  END FUNCTION ComputeComovingDistance
  
  ! --- Characteristic strain (dimensionless) from Sesana et al. (2005), Eq. (6)
  SUBROUTINE Strain( M1 , M2 , f , hc )

    REAL(DP) , INTENT(in)  :: M1 , M2 , f
    REAL(DP) , INTENT(out) :: hc
    REAL(DP)               :: Mc

    ! --- Compute chirp mass (in source frame) ---
    Mc = ( M1 * M2 )**( 3.0d0 / 5.0d0 ) / ( M1 + M2 )**( 1.0d0 / 5.0d0 ) * Msun

    hc = 1.0d0 / ( SQRT( 3.0d0 * c**3 ) * PI**( 2.0d0 / 3.0d0 ) * r ) &
           * ( G * Mc )**( 5.0d0 / 6.0d0 ) * f**( -1.0d0 / 6.0d0 )

    RETURN
  END SUBROUTINE Strain
  
END PROGRAM GWstrainFromBHBmergers
