PROGRAM GWstrainFromBHBmergers

  IMPLICIT NONE

  INTEGER  , PARAMETER   :: DP = KIND( 1.d0 )
  INTEGER  , PARAMETER   :: Nq = 10000 , Nz = 100000
  INTEGER  , PARAMETER   :: Nc = 7 , Nf_all = 400 , Nf = 10 , DELTA = Nf_all / Nf
  INTEGER  , PARAMETER   :: iID = 6 , iM1 = 7 , iM2 = 8 , itLB = 10
  REAL(DP) , PARAMETER   :: PI = ACOS( -1.0d0 )
  REAL(DP) , PARAMETER   :: OMEGA_M = 0.3d0 , OMEGA_L = 0.7d0
  REAL(DP) , PARAMETER   :: H0 = 70.4d0 / 3.086d19
  REAL(DP) , PARAMETER   :: z_max = 25.0d0 , dz = z_max / ( Nz - 1 )
  REAL(DP) , PARAMETER   :: c = 3.0d10 , G = 6.67d-8
  REAL(DP) , PARAMETER   :: tau = 4.0d0 * 86400.0d0 * 365.35d0
  REAL(DP) , PARAMETER   :: Msun = 2.0d33 , Gyr = 1.0d9 * 86400.0d0 * 365.25d0
  REAL(DP) , ALLOCATABLE :: IllustrisData(:,:)
  REAL(DP)               :: LISA_all(Nf_all,2) , LISA(Nf,2) , hc , z , tLB , r
  REAL(DP)               :: M1 , M2 , f_ISCO
  REAL(DP)               :: z_arr(Nz) , tLB_z_arr(Nz)
  INTEGER                :: nLinesIllustris , Nss , i , j , k , SS(2)
  INTEGER*4              :: MergerID
  CHARACTER( len = 28 )  :: FILEIN
  CHARACTER( len = 11 )  :: FMTIN
  CHARACTER( len = 3  )  :: arg

  ! --- Read in frequencies from LISA sensitivity curve data file ---
  OPEN( 100 , FILE = '../LISA_sensitivity.dat' )
  ! --- Loop through comments ---
  DO i = 1 , Nc
    READ( 100 , * )
  END DO
  ! --- Read in data ---
  DO i = 1 , Nf_all
    READ( 100 , * ) LISA_all( i , : )
  END DO
  CLOSE( 100 )

  ! --- Trim LISA data ---
  LISA( 1 , : ) = LISA_all( 1 , : )
  DO i = 2 , Nf
     LISA( i , : ) = LISA_all( i * DELTA , : )
  END DO
  
  ! --- Create lookback time array ---
  OPEN( 100 , FILE = 'tLB_z.dat' )
  
  ! --- Compute array of lookback times ---
  WRITE( 100 , '(A14)' ) '# z, tLB [Gyr]'
  DO i = 1 , Nz
    z_arr(i)     = ( i - 1 ) * dz
    tLB_z_arr(i) = ComputeLookbackTime( z_arr(i) ) / Gyr
    WRITE( 100 , '(ES16.10,1x,F13.10)' ) z_arr(i) , tLB_z_arr(i)
  END DO

  CLOSE( 100 )

  ! --- Create file for storing strains (first row will hold frequencies) ---
  OPEN ( 100 , FILE = 'hc.dat' )
  WRITE( 100 , '(A59)' ) &
         '# ID, M1, M2, tLB [Gyr], z, f_ISCO, r_comoving [Mpc], hc(f)'
  WRITE( 100 , '(I7,1x,F7.4,1x,F7.4,1x,F13.10,1x,ES16.10,1x, &
                   ES16.10,1x,ES16.10,10ES13.6)' ) &
           0 , 0.0d0 , 0.0d0 , 0.0d0 , 0.0d0 , 0.0d0 , 0.0d0 , LISA( : , 1 )

  ! --- Get snapshot numbers ---
  IF ( IARGC() .NE. 2 ) THEN
     WRITE(*,'(A6)')  'ERROR:'
     WRITE(*,'(A30)') 'Proper usage: ./strain SSL SSU'
     WRITE(*,'(A10)') 'Exiting...'
  END IF

  DO i = 1 , IARGC()
     CALL GETARG( i , arg )
     READ( arg , * ) SS(i)
  END DO

  ! --- Loop through Illustris snapshots ---
  DO Nss = SS(1) , SS(2)

    IF ( ( Nss .NE. 53 ) .AND. ( Nss .NE. 55 ) ) THEN
      OPEN ( 102 , FILE = 'Nss.dat' )
      WRITE( 102 , '(I3)' ) Nss
      CLOSE( 102 )

      ! --- Get filenames
      IF ( Nss < 100 ) THEN
        FMTIN = '(A21,I2,A4)'
      ELSE
        FMTIN = '(A21,I3,A4)'
      END IF

      WRITE( FILEIN  , FMTIN  ) '../time_BHillustris1_' , Nss , '.dat'

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
        M1       = IllustrisData( i , iM1  ) ! [ Solar masses ]
        M2       = IllustrisData( i , iM2  ) ! [ Solar masses ]
        tLB      = IllustrisData( i , itLB ) ! [ Gyr ]

        ! --- Compute redshift from lookback time ---
        z = InterpolateLookbackTime( tLB )
        
        ! --- Compute comoving distance (in cm) ---
        r = ComputeComovingDistance( z )

        ! --- Calculate frequency at ISCO using
        !       Sesana et al. (2005), Eq. (12) ---
        f_ISCO = c**3 / ( 6.0d0**( 3.0d0 / 2.0d0 ) * PI * G &
                   * ( M1 + M2 ) * Msun  * ( 1.0d0 + z ) )

        WRITE( 100 , '(I7,1x,F7.4,1x,F7.4,1x,F13.10,1x,ES16.10,1x, &
                         ES16.10,1x,ES16.10)' , ADVANCE = 'NO' ) &
                 MergerID , M1 , M2 , tLB , z , f_ISCO , r / 3.086d24
        
        ! --- Loop through frequencies until f_ISCO ---
        j = 1
        DO WHILE ( ( LISA( j , 1 ) < f_ISCO ) .AND. ( j < Nf + 1 ) )
          hc = Strain( M1 , M2 , LISA( j , 1 ) , z , r )
          WRITE( 100 , '(ES13.6)' , ADVANCE = 'NO' ) hc
          j = j + 1
        END DO

        ! --- Fill in missing frequencies with 0.0d0 ---
        DO WHILE ( j < Nf + 1 )
          WRITE( 100 , '(ES13.6)' , ADVANCE = 'NO' ) 0.0d0
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
    
    dz  = z / Nq

    tLB = 0.0d0
    DO k = 1 , Nq - 1
      tLB = tLB + E_LB( k * dz )
    END DO

    tLB = 1.0d0 / H0 * dz / 2.0d0 * ( E_LB( 0.0d0 ) + 2.0d0 * tLB + E_LB( z ) )
    
    RETURN
  END FUNCTION ComputeLookbackTime

  
  ! --- Characteristic strain (dimensionless)
  !       from Sesana et al. (2005), Eqs. (6,7) ---
  PURE FUNCTION Strain( M1 , M2 , f , z , r ) RESULT( hc )

    REAL(DP) , INTENT(in)  :: M1 , M2 , f , z , r
    REAL(DP)               :: hc
    REAL(DP)               :: Mc , h , n

    ! --- Compute chirp mass (in source frame) ---
    Mc = ( M1 * M2 )**( 3.0d0 / 5.0d0 ) / ( M1 + M2 )**( 1.0d0 / 5.0d0 ) * Msun

    ! --- Compute pure strain, Eq. (2) ---
    h = 8.0d0 * PI**( 2.0d0 / 3.0d0 ) / SQRT( 10.0d0 )   &
          * ( G * Mc )**( 5.0d0 / 3.0d0 ) / ( c**4 * r ) &
            * ( f * ( 1.0d0 + z ) )**( 2.0d0 / 3.0d0 )

    ! --- Compute n, Eq. (4) ---
    n = 5.0d0 / ( 96.0d0 * PI**( 8.0d0 / 3.0d0 ) ) &
          * c**5 / ( G * Mc * f * ( 1.0d0 + z ) )**( 5.0d0 / 3.0d0 )

    ! --- Check whether or not n < f * tau ---
    IF ( n < f * tau ) THEN
      hc = h * SQRT( n )
    ELSE
      hc = h * SQRT( f * tau )
    END IF
   
    RETURN
  END FUNCTION Strain

  
  ! --- Interpolate lookback time array to get redshift ---
  FUNCTION InterpolateLookbackTime( tLB ) RESULT( z )

    REAL(DP) , INTENT(in) :: tLB
    REAL(DP)              :: zmin , zmax , tLBmin , tLBmax , tLBk , m , b
    REAL(DP)              :: z

    ! --- Small lookback time approximation: tLB ~ tH * z ---
    IF ( tLB < 0.1d0 ) THEN
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

    ELSE

      ! --- Extrapolate using linear extrapolation ---

      zmin   = z_arr(Nz-1)
      zmax   = z_arr(Nz)
      tLBmin = tLB_z_arr(Nz-1)
      tLBmax = tLB_z_arr(Nz)

    END IF

    ! --- tLB = m * z + b ---
    m = ( tLBmax - tLBmin ) / ( zmax - zmin )
    b = 0.5_DP * ( ( tLBmax - m * zmax ) + ( tLBmin - m * zmin ) )

    z = ( tLB - b ) / m

    RETURN
  END FUNCTION InterpolateLookbackTime


END PROGRAM GWstrainFromBHBmergers
