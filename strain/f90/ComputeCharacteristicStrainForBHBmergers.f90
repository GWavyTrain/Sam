PROGRAM ComputeCharacteristicStrainForBHBmergers

  IMPLICIT NONE

  INTEGER,  PARAMETER   :: DP = KIND( 1.d0 )
  INTEGER,  PARAMETER   :: iSnapshotMin = 26, iSnapshotMax = 27
  INTEGER,  PARAMETER   :: Nq = 1000000, nRedshifts = 1000, &
                           nFrequencies = 400, nCommentLines = 9
  INTEGER,  PARAMETER   :: iM1 = 7, iM2 = 8, itLB = 10
  REAL(DP), PARAMETER   :: PI = ACOS( -1.0d0 )
  REAL(DP), PARAMETER   :: OMEGA_M = 0.3d0, OMEGA_L = 0.7d0, &
                           H0 = 70.4d0 / 3.086d19
  REAL(DP), PARAMETER   :: z_max = 16.0d0, &
                           dz = z_max / ( DBLE(nRedshifts) - 1.0d0 )
  REAL(DP), PARAMETER   :: CentimetersPerMpc = 3.086d24, &
                           SecondsPerGyr     = 1.0d9 * 86400.0d0 * 365.25d0
  REAL(DP), PARAMETER   :: c = 2.99792458d10, G = 6.673d-8
  REAL(DP), PARAMETER   :: Msun = 1.98892d33
  REAL(DP), ALLOCATABLE :: IllustrisData(:,:)
  REAL(DP)              :: LISA(nFrequencies,2), hc, z, tLB, r
  REAL(DP)              :: M1, M2, f_ISCO
  REAL(DP)              :: z_arr(nRedshifts), tLB_z_arr(nRedshifts)
  INTEGER               :: nMergersPerSnapshot, iSnapshot, iMerger, i
  CHARACTER(LEN=9)      :: FMTIN
  CHARACTER(LEN=128)    :: FILEIN, RootPath
  LOGICAL               :: FileExists

  WRITE( RootPath, '(A)' ) '/Users/dunhamsj/Research/GW/Sam/'
!  WRITE( RootPath, '(A)' ) '/astro1/dunhamsj/'

  ! === Read in frequencies from LISA sensitivity curve data file ===
  OPEN( 100, FILE = TRIM(RootPath) // 'LISA_sensitivity.dat' )
  ! --- Loop through comment lines ---
  DO i = 1, nCommentLines
    READ( 100, * )
  END DO
  ! --- Loop through data ---
  DO i = 1, nFrequencies
    READ( 100, * ) LISA(i,:)
    READ( 100, * )
  END DO
  CLOSE( 100 )

  ! === Compute array of lookback-times ===
  OPEN( 100, FILE = TRIM(RootPath) // 'strain/f90/tLB_z.dat' )
  WRITE( 100, '(A)' ) '# Redshift z, Lookback-Time tLB [Gyr]'
  DO i = 1, nRedshifts
    z_arr    (i) = ( DBLE(i) - 1.0d0 ) * dz
    tLB_z_arr(i) = ComputeLookbackTime( z_arr(i) )
    WRITE( 100, '(ES23.16E3,1x,ES23.16E3)' ) z_arr(i), tLB_z_arr(i)
  END DO
  CLOSE( 100 )

  ! --- Create file for storing strains (first row will hold frequencies) ---
  OPEN ( 100, FILE = TRIM(RootPath) // 'strain/f90/hc.dat' )
  WRITE( 100, '(A)' ) '# M1, M2, tLB, z, r, hc(f)'
  WRITE( 100, '(ES12.6E2,ES13.6E2,ES18.11E2, &
              & ES18.11E2,E18.11E2,400ES18.11E2)' ) &
                0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, LISA(:,1)

  ! --- Loop through Illustris snapshots ---
  DO iSnapshot = iSnapshotMin, iSnapshotMax

    ! --- Get filenames ---
    IF ( iSnapshot .LT. 100 ) THEN
      FMTIN = '(A,I2,A4)'
    ELSE
      FMTIN = '(A,I3,A4)'
    END IF

    WRITE( FILEIN, FMTIN ) TRIM(RootPath) // &
                           'BBHM_DataFiles_Mapelli/time_BHillustris1_', &
                           iSnapshot, '.dat'

    INQUIRE( FILE = FILEIN, EXIST = FileExists )
    IF( .NOT. FileExists ) CYCLE

    ! --- Get number of lines (mergers) in Illustris data file ---
    nMergersPerSnapshot = 0
    OPEN( 101, FILE = TRIM( FILEIN ) )
    DO
      READ( 101, *, END = 11 )
      nMergersPerSnapshot = nMergersPerSnapshot + 1
    END DO
    11 CLOSE( 101 )

    ! --- Read in Illustris data and compute strain for each merger ---
    ALLOCATE( IllustrisData(nMergersPerSnapshot,10) )
    OPEN( 101, FILE = TRIM( FILEIN ) )
    DO iMerger = 1, nMergersPerSnapshot
      READ( 101, * ) IllustrisData( iMerger, : )

      ! --- Get data for this merger ---
      M1  = IllustrisData( iMerger, iM1  )
      M2  = IllustrisData( iMerger, iM2  )
      tLB = IllustrisData( iMerger, itLB )

      ! --- Compute redshift from lookback-time ---
      z = InterpolateLookbackTime( tLB )
        
      ! --- Compute comoving distance (in Mpc) ---
      r = ComputeComovingDistance( z )

      WRITE( 100, '(ES12.6E2,ES13.6E2,ES18.11E2,ES18.11E2,ES18.11E2)', &
        ADVANCE = 'NO' ) M1, M2, tLB, z, r
        
      ! --- Calculate frequency at ISCO using Sesana et al., (2005) ---
      f_ISCO = c**3 / ( 6.0d0**( 3.0d0 / 2.0d0 ) * PI * G &
                 * ( M1 + M2 ) * Msun  * ( 1.0d0 + z ) )

      ! --- Loop through frequencies until f_ISCO ---
      i = 1
      DO WHILE( ( LISA(i,1) .LT. f_ISCO ) .AND. ( i .LT. nFrequencies + 1 ) )
        hc = ComputeCharacteristicStrain( LISA(i,1) )
        WRITE( 100, '(ES18.11E2)', ADVANCE = 'NO' ) hc
        i = i + 1
      END DO

      ! --- Fill in missing frequencies with 0.0d0 ---
      DO WHILE ( i < nFrequencies + 1 )
        WRITE( 100, '(f4.1)' , ADVANCE = 'NO' ) 0.0d0
        i = i + 1
      END DO

      WRITE( 100, * )
         
    END DO

    CLOSE( 101 )
    DEALLOCATE( IllustrisData )

  END DO

  CLOSE( 100 )


CONTAINS


  ! --- Integrand E(z) in cosmological distance calculation ---
  PURE REAL(DP) FUNCTION E( zz )

    REAL(DP), INTENT(in) :: zz

    E = 1.0d0 / SQRT( OMEGA_M * ( 1.0d0 + zz )**3 + OMEGA_L )

    RETURN
  END FUNCTION E


  FUNCTION ComputeComovingDistance( zz ) RESULT( rr )

    REAL(DP), INTENT(in) :: zz
    REAL(DP)             :: rr, dzz
    INTEGER              :: k

    ! --- Integrate with Trapezoidal rule ---
    rr  = 0.0d0
    dzz = zz / Nq
    
    DO k = 1, Nq - 1
       rr = rr + E( k * dzz )
    END DO

    rr = c / H0 * dzz / 2.0d0 * ( E( 0.0d0 ) + 2.0d0 * rr + E( zz ) ) &
           / CentimetersPerMpc

    RETURN
  END FUNCTION ComputeComovingDistance

  
  ! --- Integrand in lookback time calculation ---
  PURE REAL(DP) FUNCTION E_LB( zz )

    REAL(DP), INTENT(in) :: zz

    E_LB = 1.0d0 / ( ( 1.0d0 + zz ) &
             * SQRT( OMEGA_M * ( 1.0d0 + zz )**3 + OMEGA_L ) )
    
    RETURN
  END FUNCTION E_LB

  
  REAL(DP) FUNCTION ComputeLookbackTime( zz ) RESULT( ttLB )

    REAL(DP), INTENT(in) :: zz
    REAL(DP)             :: dzz
    INTEGER              :: k

    dzz = zz / DBLE(Nq)

    ! --- Integrate with Trapezoidal rule ---
    ttLB = 0.0d0
    DO k = 1, Nq - 1
      ttLB = tLB + E_LB( k * dzz )
    END DO

    ttLB = 1.0d0 / H0 * dzz / 2.0d0 &
             * ( E_LB( 0.0d0 ) + 2.0d0 * ttLB + E_LB( zz ) ) / SecondsPerGyr
    
    RETURN
  END FUNCTION ComputeLookbackTime

  
  ! --- Characteristic strain (dimensionless)
  !       from Sesana et al. (2005), Eq. (6) ---
  PURE FUNCTION ComputeCharacteristicStrain( f ) RESULT( hhc )

    REAL(DP) , INTENT(in)  :: f
    REAL(DP)               :: Mc
    REAL(DP)               :: hhc

    ! --- Compute chirp mass (in source frame) ---
    Mc = ( M1 * M2 )**( 3.0d0 / 5.0d0 ) / ( M1 + M2 )**( 1.0d0 / 5.0d0 ) * Msun

    hhc = 1.0d0 / ( SQRT( 3.0d0 * c**3 ) * PI**( 2.0d0 / 3.0d0 ) * r ) &
           * ( G * Mc )**( 5.0d0 / 6.0d0 ) * f**( -1.0d0 / 6.0d0 )

    RETURN
  END FUNCTION ComputeCharacteristicStrain

  
  ! --- Interpolate lookback-time array to get redshift ---
  FUNCTION InterpolateLookbackTime( ttLb ) RESULT( z )

    REAL(DP), INTENT(in) :: ttLb
    REAL(DP)             :: zMin, zMax, tLbMin, tLbMax, tLbi, m, b, z
    INTEGER              :: k

    ! --- Small lookback time approximation: ttLb ~ tH * z ---
    IF ( ttLb .LT. 0.5d0 ) THEN
      z = ttLb * ( H0 * SecondsPerGyr )
      RETURN
    END IF

    tLbi = tLb_z_arr(1)
    IF ( ttLb < tLb_z_arr(nRedshifts) ) THEN

       ! --- Interpolate using linear interpolation ---

       ! --- Get redshift bounds ---
       k = 1
       DO WHILE ( tLbi .LE. ttLb )
          k = k + 1
          tLbi = tLb_z_arr(k)
       END DO

       zMin   = z_arr(k)
       zMax   = z_arr(k+1)
       tLbMin = tLb_z_arr(k)
       tLbMax = tLb_z_arr(k+1)

       ! --- ttLb = m * z + b ---
       m = ( tLbMax - tLbMin ) / ( zMax - zMin )
       b = 0.5d0 * ( ( tLbMax - m * zMax ) + ( tLbMin - m * zMin ) )

    ELSE

       ! --- Extrapolate using linear extrapolation ---

       zMin   = z_arr(nRedshifts-1)
       zMax   = z_arr(nRedshifts)
       tLbMin = tLb_z_arr(nRedshifts-1)
       tLbMax = tLb_z_arr(nRedshifts)

       ! --- ttLb = m * z + b ---
       m = ( tLbMax - tLbMin ) / ( zMax - zMin )
       b = 0.5d0 * ( ( tLbMax - m * zMax ) + ( tLbMin - m * zMin ) )

    END IF

    z = ( ttLb - b ) / m

    RETURN
  END FUNCTION InterpolateLookbackTime


END PROGRAM ComputeCharacteristicStrainForBHBmergers
