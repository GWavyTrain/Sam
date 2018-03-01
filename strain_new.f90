PROGRAM GWstrainFromBHBmergers

  IMPLICIT NONE

  INTEGER  , PARAMETER   :: DP = KIND( 1.d0 )
  INTEGER  , PARAMETER   :: Nq = 10000 , Nz = 100000 , Nf = 400 , Nc = 7
  INTEGER  , PARAMETER   :: iID = 6 , iM1 = 7 , iM2 = 8 , itLB = 10
  REAL(DP) , PARAMETER   :: PI = ACOS( -1.0d0 )
  REAL(DP) , PARAMETER   :: OMEGA_M = 0.3d0 , OMEGA_L = 0.7d0
  REAL(DP) , PARAMETER   :: H0 = 70.4d0 / 3.086d19
  REAL(DP) , PARAMETER   :: z_max = 1000.0d0 , dz = z_max / ( Nz - 1 )
  REAL(DP) , PARAMETER   :: c = 3.0d10 , G = 6.67d-8
  REAL(DP) , PARAMETER   :: Msun = 2.0d33 , Gyr = 1.0d9 * 86400.0d0 * 365.25d0
  REAL(DP) , ALLOCATABLE :: IllustrisData(:,:)
  REAL(DP)               :: LISA(Nf,2) , hc , z , tLB , r
  REAL(DP)               :: M1 , M2 , f_ISCO
  REAL(DP)               :: z_arr(Nz) , tLB_z_arr(Nz)
  INTEGER                :: nLinesIllustris , Nss , i , j , k
  INTEGER*4              :: MergerID
  CHARACTER( len = 25 )  :: FILEIN
  CHARACTER( len = 11 )  :: FMTIN

  ! --- Read in frequencies from LISA sensitivity curve data file ---
  OPEN( 100 , FILE = 'LISA_sensitivity.dat' )
  ! --- Loop through comments ---
  DO i = 1 , Nc
    READ( 100 , * )
  END DO
  ! --- Read in data
  DO i = 1 , Nf
    READ( 100 , * ) LISA( i , : )
  END DO
  CLOSE( 100 )

  ! --- Create or read in lookback time array ---
  OPEN( 100 , FILE = 'tLB_z.dat' )
  
  ! --- Compute array of lookback times ---
  WRITE(*,'(A34)') 'Computing array of lookback times:'
  WRITE( 100 , '(A9)' ) '# z , tLB'
  DO i = 1 , Nz
    IF ( MOD( i , Nz / 100 ) .EQ. 0 ) &
      WRITE(*,'(A7,I6,A1,I6)' )   'i/Nz = ' , i , '/' , Nz
      z_arr(i)     = (i-1) * dz
      tLB_z_arr(i) = ComputeLookbackTime( z_arr(i) ) / Gyr
      WRITE( 100 , '(E16.10,1x,F13.10)' ) z_arr(i) , tLB_z_arr(i)
  END DO

  CLOSE( 100 )

  ! --- Create file for storing strains (first row will hold frequencies) ---
  OPEN ( 100 , FILE = 'hc.dat' )
  WRITE( 100 , '(A31)' ) '# ID, M1, M2, tLB , z, r, hc(f)'
  WRITE( 100 , '(I7,1x,F7.4,1x,F7.4,1x,F13.10,1x,E16.10,1x,E16.10,400E13.6)' ) &
    0 , 0.0d0 , 0.0d0 , 0.0d0 , 0.0d0 , 0.0d0 , LISA( : , 1 )

  ! --- Loop through Illustris snapshots ---
  DO Nss = 26 , 27!135

    IF ( ( Nss .NE. 53 ) .AND. ( Nss .NE. 55 ) ) THEN
      !OPEN ( 102 , FILE = 'Nss.dat' )
      !WRITE( 102 , '(I3)' ) Nss
      !CLOSE( 102 )

      ! --- Get filenames
      IF ( Nss < 100 ) THEN
        FMTIN = '(A18,I2,A4)'
      ELSE
        FMTIN = '(A18,I3,A4)'
      END IF

      WRITE( FILEIN  , FMTIN  ) 'time_BHillustris1_' , Nss , '.dat'

      ! --- Get number of lines (mergers) in Illustris data file ---
      nLinesIllustris = 0
      OPEN( 101 , FILE = TRIM( FILEIN ) )
      DO
        READ( 101 , * , END = 11 )
        nLinesIllustris = nLinesIllustris + 1
      END DO
      11 CLOSE( 101 )

      ! --- Read in Illustris data and compute strain for each merger ---
      ALLOCATE( IllustrisData( nLinesIllustris , 10 ) )
      OPEN( 101 , FILE = TRIM( FILEIN ) )
      DO i = 1 , nLinesIllustris
        READ( 101 , * ) IllustrisData( i , : )

        ! --- Get data for this merger ---
        MergerID = IllustrisData( i , iID  )
        M1       = IllustrisData( i , iM1  )
        M2       = IllustrisData( i , iM2  )
        tLB      = IllustrisData( i , itLB )

        ! --- Compute redshift from lookback time ---
        z = InterpolateLookbackTime( tLB )
        
        ! --- Compute comoving distance (in cm) ---
        r = ComputeComovingDistance( z )

        WRITE( 100 , '(I7,1x,F7.4,1x,F7.4,1x,F13.10,1x,E16.10,1x,E16.10)' , &
                 ADVANCE = 'NO' ) MergerID , M1 , M2 , tLB , z , r / 3.086d24
        
        ! --- Calculate frequency at ISCO using Sesana et al., (2005) ---
        f_ISCO = c**3 / ( 6.0d0**( 3.0d0 / 2.0d0 ) * PI * G &
                   * ( M1 + M2 ) * Msun  * ( 1.0d0 + z ) )

        ! --- Loop through frequencies until f_ISCO ---
        j = 1
        DO WHILE ( ( LISA( j , 1 ) < f_ISCO ) .AND. ( j < Nf + 1 ) )
          hc = Strain( M1 , M2 , LISA( j , 1 ) )
          WRITE( 100 , '(E13.6)' , ADVANCE = 'NO' ) hc
          j = j + 1
        END DO

        ! --- Fill in missing frequencies with 0.0d0 ---
        DO WHILE ( j < Nf + 1 )
          WRITE( 100 , '(E13.6)' , ADVANCE = 'NO' ) 0.0d0
          j = j + 1
        END DO

        WRITE( 100 , * )
         
      END DO

      CLOSE( 101 )
      DEALLOCATE( IllustrisData )

    END IF

  END DO

  CLOSE( 100 )
  
CONTAINS

  ! --- Integrand E(z) in cosmological distance calculation ---
  PURE FUNCTION E( z )

    REAL(DP) , INTENT(in) :: z
    REAL(DP)              :: E

    E = 1.0d0 / SQRT( OMEGA_M * ( 1.0d0 + z )**3 + OMEGA_L )

    RETURN
  END FUNCTION E

  
  FUNCTION ComputeComovingDistance( z ) RESULT( r )

    REAL(DP) , INTENT(in) :: z
    REAL(DP)              :: r , dz

    ! --- Integrate with Trapezoidal rule ---
    r  = 0.0d0
    dz = z / Nq
    
    DO k = 1 , Nq - 1
       r = r + E( k * dz )
    END DO

    r = c / H0 * dz / 2.0d0 * ( E( 0.0d0 ) + 2.0d0 * r + E( z ) )

    RETURN
  END FUNCTION ComputeComovingDistance

  
  ! --- Integrand in lookback time calculation ---
  PURE FUNCTION E_LB( z )

    REAL(DP) , INTENT(in) :: z
    REAL(DP)              :: E_LB

    E_LB = 1.0d0 / ( ( 1.0d0 + z ) &
             * SQRT( OMEGA_M * ( 1.0d0 + z )**3 + OMEGA_L ) )
    
    RETURN
  END FUNCTION E_LB

  
  FUNCTION ComputeLookbackTime( z ) RESULT( tLB )

    REAL(DP) , INTENT(in) :: z
    REAL(DP)              :: dz , Error
    REAL(DP)              :: tLB
    
    ! --- Integrate with Trapezoidal rule ---
    ! --- |ERROR|: Nq = 1000   --> 1.47d-3
    !              Nq = 10000  --> 1.47d-5
    !              Nq = 100000 --> 1.47d-7 ---
    
    tLB = 0.0d0
    dz  = z / Nq

    tLB = 0.0d0
    DO k = 1 , Nq - 1
      tLB = tLB + E_LB( k * dz )
    END DO

    tLB = 1.0d0 / H0 * dz / 2.0d0 * ( E_LB( 0.0d0 ) + 2.0d0 * tLB + E_LB( z ) )
    
    RETURN
  END FUNCTION ComputeLookbackTime

  
  ! --- Characteristic strain (dimensionless)
  !       from Sesana et al. (2005), Eq. (6) ---
  PURE FUNCTION Strain( M1 , M2 , f ) RESULT( hc )

    REAL(DP) , INTENT(in)  :: M1 , M2 , f
    REAL(DP)               :: hc
    REAL(DP)               :: Mc

    ! --- Compute chirp mass (in source frame) ---
    Mc = ( M1 * M2 )**( 3.0d0 / 5.0d0 ) / ( M1 + M2 )**( 1.0d0 / 5.0d0 ) * Msun

    hc = 1.0d0 / ( SQRT( 3.0d0 * c**3 ) * PI**( 2.0d0 / 3.0d0 ) * r ) &
           * ( G * Mc )**( 5.0d0 / 6.0d0 ) * f**( -1.0d0 / 6.0d0 )

    RETURN
  END FUNCTION Strain

  
  ! --- Interpolate lookback time array to get redshift ---
  FUNCTION InterpolateLookbackTime( tLB ) RESULT( z )

    REAL(DP) , INTENT(in) :: tLB
    REAL(DP)              :: zmin , zmax , tLBmin , tLBmax , tLBk , m , b
    REAL(DP)              :: z

    ! --- Small lookback time approximation: tLB ~ tH * z ---
    IF ( tLB < 0.5d0 ) THEN
      z = tLB * ( H0 * Gyr )
      RETURN
    END IF

    IF ( tLB < tLB_z_arr(Nz) ) THEN

       ! --- Interpolate using linear interpolation ---

       ! --- Get redshift bounds ---
       k = 1
       tLBk = tLB_z_arr(k)
       DO WHILE ( tLBk < tLB )
          k = k + 1
          tLBk = tLB_z_arr(k)
       END DO

       zmin   = z_arr(k-1)
       zmax   = z_arr(k)
       tLBmin = tLB_z_arr(k-1)
       tLBmax = tLB_z_arr(k)

       ! --- tLB = m * z + b ---
       m = ( tLBmax - tLBmin ) / ( zmax - zmin )
       b = 0.5_DP * ( ( tLBmax - m * zmax ) + ( tLBmin - m * zmin ) )

    ELSE

       ! --- Extrapolate using linear extrapolation ---

       zmin   = z_arr(Nz-1)
       zmax   = z_arr(Nz)
       tLBmin = tLB_z_arr(Nz-1)
       tLBmax = tLB_z_arr(Nz)

       ! --- tLB = m * z + b ---
       m = ( tLBmax - tLBmin ) / ( zmax - zmin )
       b = 0.5_DP * ( ( tLBmax - m * zmax ) + ( tLBmin - m * zmin ) )

    END IF

    z = ( tLB - b ) / m

    RETURN
  END FUNCTION InterpolateLookbackTime


END PROGRAM GWstrainFromBHBmergers
