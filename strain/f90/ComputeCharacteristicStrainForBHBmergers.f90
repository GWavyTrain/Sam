PROGRAM ComputeCharacteristicStrainForBHBmergers

  USE UtilitiesModule, ONLY: &
    RombergIntegration

  IMPLICIT NONE

  INTEGER,  PARAMETER   :: DP = KIND( 1.d0 )
  INTEGER,  PARAMETER   :: iSnapshotMin = 26, iSnapshotMax = 26
  INTEGER,  PARAMETER   :: nRedshifts = 10000, &
                           nFrequenciesAll = 400, nFrequencies = 10, &
                           nSkip = nFrequenciesAll / nFrequencies, &
                           nCommentLines = 9
  INTEGER,  PARAMETER   :: iM1 = 7, iM2 = 8, itLb = 10
  REAL(DP), PARAMETER   :: PI = ACOS( -1.0d0 )
  REAL(DP), PARAMETER   :: OMEGA_M = 0.2726d0, OMEGA_L = 0.7274d0, &
                           H0 = 70.4d0 / 3.086d19
  REAL(DP), PARAMETER   :: z_max = 20.0d0, &
                           dz = z_max / ( DBLE(nRedshifts) - 1.0d0 )
  REAL(DP), PARAMETER   :: CentimetersPerMpc = 3.086d24, &
                           Year = 86400.0d0 * 365.0d0, &
                           SecondsPerGyr = 1.0d9 * Year
  REAL(DP), PARAMETER   :: c = 2.99792458d10, G = 6.673d-8
  REAL(DP), PARAMETER   :: Msun = 1.98892d33, tau = 4.0d0 * Year
  REAL(DP)              :: IllustrisData(10)
  REAL(DP)              :: hc, z, tLb, r, M1, M2, f_ISCO
  REAL(DP), ALLOCATABLE :: LISA_all(:,:), LISA(:), z_arr(:), tLb_z_arr(:)
  INTEGER               :: i, nMergersPerSnapshot, iSnapshot, &
                           iMerger, iStrainFile
  CHARACTER(LEN=9)      :: FMTIN, FMTOUT
  CHARACTER(LEN=128)    :: FILEIN, FILEOUT, RootPath, &
                           LookbackTimeRedshiftFile, WriteFile, FMThc
  LOGICAL               :: FileExists
#ifdef _OPENMP
  INTEGER               :: iThread, nThreads, &
                           OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
#endif

  ALLOCATE( LISA_all (1:nFrequenciesAll,1:2) )
  ALLOCATE( LISA     (1:nFrequencies) )
  ALLOCATE( z_arr    (1:nRedshifts) )
  ALLOCATE( tLb_z_arr(1:nRedshifts) )

  WRITE( RootPath, '(A)' ) '/Users/dunhamsj/Research/GW/Sam/'
!  WRITE( RootPath, '(A)' ) '/astro1/dunhamsj/'

  ! === Read in frequencies from LISA sensitivity curve data file ===
  OPEN( 100, FILE = TRIM( RootPath ) // 'LISA_sensitivity.dat' )
  ! --- Loop through comment lines ---
  DO i = 1, nCommentLines
    READ( 100, * )
  END DO
  ! --- Loop through data ---
  DO i = 1, nFrequenciesAll
    READ( 100, * ) LISA_all(i,:)
    READ( 100, * )
  END DO
  CLOSE( 100 )

  ! --- Trim LISA data ---
  LISA(1) = LISA_all(1,1)
  DO i = 2, nFrequencies-1
    iMerger = nSkip * (i-1)
    LISA(i) = LISA_all(iMerger,1)
  END DO
  LISA(nFrequencies) = LISA_all(nFrequenciesAll,1)

  ! === Compute array of lookback-times ===
  WRITE( LookbackTimeRedshiftFile, '(A,I10.10,A)' ) &
         TRIM(RootPath) // 'strain/f90/tLb_z_', nRedshifts, '.dat'

  INQUIRE( FILE = TRIM( LookbackTimeRedshiftFile ), EXIST = FileExists )
  IF( .NOT. FileExists )THEN
    WriteFile = 'Y'
  ELSE
    WRITE(*,'(A,A,A)') &
      'File ', TRIM( LookbackTimeRedshiftFile ), ' exists. Overwrite? (Y/N) '
    READ(*,*) WriteFile
  END IF

  IF( TRIM( WriteFile ) .EQ. 'Y' )THEN
    WRITE(*,'(A,A)') 'Writing file: ', TRIM( LookbackTimeRedshiftFile )

    DO i = 1, nRedshifts
      z_arr(i) = ( DBLE(i) - 1.0d0 ) * dz
    END DO

    !$OMP PARALLEL DEFAULT(NONE) SHARED(tLb_z_arr,z_arr)
    !$OMP DO
    DO i = 1, nRedshifts
      CALL RombergIntegration &
             ( Integrand_tLb, 0.0d0, z_arr(i), tLb_z_arr(i) )
    END DO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    OPEN( 100, FILE = TRIM( LookbackTimeRedshiftFile ) )
    WRITE( 100, '(A)' ) '# Redshift z, Lookback-Time tLb [Gyr]'
    DO i = 1, nRedshifts
      tLb_z_arr(i) = tLb_z_arr(i) / H0 / SecondsPerGyr
      WRITE( 100, '(F14.11,1x,F14.11)' ) z_arr(i), tLb_z_arr(i)
    END DO
    CLOSE( 100 )

  ELSE
    OPEN( 100, FILE = TRIM( LookbackTimeRedshiftFile ) )
    READ( 100, * )
    DO i = 1, nRedshifts
      READ( 100, '(F14.11,1x,F14.11)' ) z_arr(i), tLb_z_arr(i)
    END DO
    CLOSE( 100 )
  END IF

  WRITE( FMThc, '(A)' ) '(2F8.4,2F14.10,3ES12.5E2)'
  ! === Loop through Illustris snapshots ===

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & PRIVATE(iSnapshot,iThread,iStrainFile,FMTIN,FMTOUT,FILEIN,FILEOUT, &
  !$OMP &         FileExists,nMergersPerSnapshot,IllustrisData,M1,M2,tLb, &
  !$OMP &         z,r,f_ISCO,i,hc) &
  !$OMP & SHARED(nThreads,RootPath,LISA,z_arr,tLb_z_arr,FMThc)

  !$OMP DO
  DO iSnapshot = iSnapshotMin, iSnapshotMax

#ifdef _OPENMP
    iThread = OMP_GET_THREAD_NUM()
    IF( iThread .EQ. 0 ) THEN
      nThreads = OMP_GET_NUM_THREADS()
      WRITE(*,'(A,I2.2)') 'Number of available threads = ', nThreads
    END IF
    iStrainFile = iThread
#else
    iStrainFile = 101
#endif

    ! --- Get filenames ---
    IF( iSnapshot .LT. 100 ) THEN
      FMTIN  = '(A,I2,A4)'
      FMTOUT = '(A,I2,A4)'
    ELSE
      FMTIN  = '(A,I3,A4)'
      FMTOUT = '(A,I3,A4)'
    END IF

    WRITE( FILEIN, FMTIN ) TRIM( RootPath ) // &
                           'BBHM_DataFiles_Mapelli/time_BHillustris1_', &
                           iSnapshot, '.dat'

    INQUIRE( FILE = FILEIN, EXIST = FileExists )
    IF( .NOT. FileExists ) CYCLE

    WRITE( FILEOUT, FMTOUT ) &
      TRIM( RootPath ) // 'strain/f90/StrainFiles/hc_', iSnapshot, '.dat'

    ! --- Create file for storing strains (first row will hold frequencies) ---
    OPEN ( iStrainFile, FILE = TRIM( FILEOUT ) )
    WRITE( iStrainFile, '(A,I3.3,A)' ) &
      '# M1 [Msun], M2 [Msun], Lookback-time to merger [Gyr],&
      & Redshift at merge, Comoving distance to merger [Mpc], fISCO,&
      & Characteristic strain (', nFrequencies, ') entries)'
    WRITE( iStrainFile, '(A)' ) '# (First row contains LISA frequencies [Hz])'
    WRITE( iStrainFile, TRIM( FMThc ), ADVANCE = 'NO' ) &
      0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0
    DO i = 1, nFrequencies
      WRITE( iStrainFile, '(ES24.16E2)', ADVANCE = 'NO' ) LISA(i)
    END DO

    ! --- Get number of lines (mergers) in Illustris data file ---
    nMergersPerSnapshot = 0
    OPEN( iSnapshot, FILE = TRIM( FILEIN ) )
    DO
      READ( iSnapshot, *, END = 11 )
      nMergersPerSnapshot = nMergersPerSnapshot + 1
    END DO
    11 CLOSE( iSnapshot )

    ! --- Read in Illustris data and compute strain for each merger ---
    OPEN( iSnapshot, FILE = TRIM( FILEIN ) )
    DO iMerger = 1, nMergersPerSnapshot
      READ( iSnapshot, * ) IllustrisData(:)

      ! --- Get data for this merger ---
      M1  = IllustrisData( iM1  )
      M2  = IllustrisData( iM2  )
      tLb = IllustrisData( itLb )

      ! --- Compute redshift from lookback-time ---
      z = InterpolateLookbackTimeToGetRedshift( tLb, z_arr, tLb_z_arr )
        
      ! --- Compute comoving distance (in Mpc) ---
      CALL RombergIntegration( Integrand_r, 0.0d0, z, r )
      r = r * c / H0

      ! --- Calculate frequency at ISCO using Sesana et al., (2005) ---
      f_ISCO = c**3 / ( 6.0d0**( 3.0d0 / 2.0d0 ) * PI * G &
                 * ( M1 + M2 ) * Msun  * ( 1.0d0 + z ) )

      WRITE( iStrainFile, TRIM( FMThc ), ADVANCE = 'NO' ) &
        M1, M2, tLb, z, r / CentimetersPerMpc, f_ISCO
        
      ! --- Loop through frequencies until f_ISCO ---
      i = 1
      DO WHILE( ( LISA(i) .LT. f_ISCO ) .AND. ( i .LE. nFrequencies ) )
        hc = ComputeCharacteristicStrain( M1, M2, LISA(i), r, z )
        WRITE( iStrainFile, '(ES12.5E2)', ADVANCE = 'NO' ) hc
        i = i + 1
      END DO

      ! --- Fill in missing frequencies with 0.0d0 ---
      DO WHILE( i .LE. nFrequencies )
        WRITE( iStrainFile, '(I2.1)' , ADVANCE = 'NO' ) 0
        i = i + 1
      END DO

      WRITE( iStrainFile, * )
         
    END DO

    CLOSE( iSnapshot )

    CLOSE( iStrainFile )

  END DO
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL

  DEALLOCATE( tLb_z_arr )
  DEALLOCATE( z_arr     )
  DEALLOCATE( LISA      )
  DEALLOCATE( LISA_all )


CONTAINS


  ! --- Integrand E(z) in cosmological distance calculation ---
  PURE REAL(DP) FUNCTION Integrand_r( z )

    REAL(DP), INTENT(in) :: z

    Integrand_r = 1.0d0 / SQRT( OMEGA_M * ( 1.0d0 + z )**3 + OMEGA_L )

    RETURN
  END FUNCTION Integrand_r


  ! --- Integrand in lookback time calculation ---
  PURE REAL(DP) FUNCTION Integrand_tLb( z )

    REAL(DP), INTENT(in) :: z

    Integrand_tLb = 1.0d0 / ( ( 1.0d0 + z ) &
                      * SQRT( OMEGA_M * ( 1.0d0 + z )**3 + OMEGA_L ) )
    
    RETURN
  END FUNCTION Integrand_tLb

  
  ! --- Characteristic strain (dimensionless)
  !       from Sesana et al. (2005), ApJ, 623, 23 ---
  FUNCTION ComputeCharacteristicStrain( M1, M2, f, r, z ) RESULT( hc )

    REAL(DP), INTENT(in) :: M1, M2, f, r, z
    REAL(DP)             :: Mc
    REAL(DP)             :: hc, h, fr, n

    ! --- Convert observation-frame redshift to source(rest)-frame redshift ---
    fr = f * ( 1.0d0 + z )

    ! --- Compute chirp mass ---
    Mc = ( ( M1 * M2 )**3 / ( M1 + M2 ) )**( 1.0d0 / 5.0d0 ) * Msun

    ! --- Compute strain amplitude, Eq. (2) ---
    h = 8.0d0 * ( PI**2 * ( G * Mc )**5 * fr**2 )**( 1.0d0 / 3.0d0 ) &
          / ( SQRT( 10.0d0 ) * c**4 * r )

    ! --- Compute number of cycles in frequency range df~f, Eq. (4) ---
    n = 5.0d0 / 96.0d0 * c**5 &
          / ( PI**8 * ( G * Mc * fr )**5 )**( 1.0d0 / 3.0d0 )

    IF( n .LT. f * tau )THEN
      hc = h * SQRT( n )
    ELSE
      hc = h * SQRT( f * tau )
    END IF

    RETURN
  END FUNCTION ComputeCharacteristicStrain

  
  ! --- Interpolate lookback-time array to get redshift ---
  FUNCTION InterpolateLookbackTimeToGetRedshift &
    ( tLb, z_arr, tLb_z_arr ) RESULT( z )

    ! --- tLb and tLb_z_arr should be input in units of Gyr ---

    REAL(DP), INTENT(in) :: tLb, z_arr(nRedshifts), tLb_z_arr(nRedshifts)
    REAL(DP)             :: zMin, zMax, tLbMin, tLbMax, tLb_k, m, b, z
    INTEGER              :: k
    LOGICAL              :: DEBUG = .FALSE.

    ! --- Small lookback time approximation: tLb ~ tH * z ---
    IF( tLb .LT. 0.5d0 ) THEN
      IF( DEBUG ) WRITE(*,'(A5,A)') &
                    '', 'Using small lookback-time approximation...'
      z = tLb * ( H0 * SecondsPerGyr )
      RETURN
    END IF

    ! --- Interpolate using linear interpolation ---
    IF( tLb .LT. tLb_z_arr(nRedshifts-1) ) THEN

      IF( DEBUG ) WRITE(*,'(A5,A)') '', 'Interpolating...'

      ! --- Get redshift bounds ---
      k = 1
      tLb_k = tLb_z_arr(k)
      DO WHILE( ( tLb_k .LT. tLb ) .AND. ( k .LE. nRedshifts ) )
         k = k + 1
         tLb_k = tLb_z_arr(k)
      END DO

      zMin   = z_arr    (k-1)
      zMax   = z_arr    (k  )
      tLbMin = tLb_z_arr(k-1)
      tLbMax = tLb_z_arr(k  )

    ! --- Extrapolate using linear extrapolation ---
    ELSE

      IF( DEBUG ) WRITE(*,'(A5,A)') '', 'Extrapolating...'

      zMin   = z_arr    (nRedshifts-1)
      zMax   = z_arr    (nRedshifts  )
      tLbMin = tLb_z_arr(nRedshifts-1)
      tLbMax = tLb_z_arr(nRedshifts  )

    END IF

    ! --- tLb = m * z + b ---
    m = ( tLbMax - tLbMin ) / ( zMax - zMin )
    b = 0.5d0 * ( ( tLbMax - m * zMax ) + ( tLbMin - m * zMin ) )

    z = ( tLb - b ) / m

    IF( DEBUG )THEN
      WRITE(*,'(A5,A,F13.10)') '', 'tLbMin = ', tLbMin
      WRITE(*,'(A5,A,F13.10)') '', 'tLb    = ', tLb
      WRITE(*,'(A5,A,F13.10)') '', 'tLbMax = ', tLbMax
      WRITE(*,'(A5,A,F13.10)') '', 'z      = ', z
      WRITE(*,*)
    END IF

    RETURN
  END FUNCTION InterpolateLookbackTimeToGetRedshift


END PROGRAM ComputeCharacteristicStrainForBHBmergers
