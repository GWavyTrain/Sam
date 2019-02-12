PROGRAM ComputeCharacteristicStrainForBHBmergers

  USE UtilitiesModule

  IMPLICIT NONE

  INTEGER,  PARAMETER :: DP = KIND( 1.d0 )
  INTEGER,  PARAMETER :: iSnapshotMin = 26, iSnapshotMax = 26
  INTEGER,  PARAMETER :: Nq = 1000000, nRedshifts = 1000, &
                         nFrequenciesAll = 400, nFrequencies = 10, &
                         nSkip = nFrequenciesAll / nFrequencies, &
                         nCommentLines = 9
  INTEGER,  PARAMETER :: iM1 = 7, iM2 = 8, itLb = 10
  REAL(DP), PARAMETER :: PI = ACOS( -1.0d0 )
  REAL(DP), PARAMETER :: OMEGA_M = 0.3d0, OMEGA_L = 0.7d0, &
                         H0 = 70.4d0 / 3.086d19
  REAL(DP), PARAMETER :: z_max = 25.0d0, &
                         dz = z_max / ( DBLE(nRedshifts) - 1.0d0 )
  REAL(DP), PARAMETER :: CentimetersPerMpc = 3.086d24, &
                         SecondsPerGyr = 86400.0d0 * 365.25d0 * 1.0d9
  REAL(DP), PARAMETER :: c = 2.99792458d10, G = 6.673d-8
  REAL(DP), PARAMETER :: Msun = 1.98892d33
  REAL(DP)            :: IllustrisData(10)
  REAL(DP)            :: LISA_all(nFrequenciesAll,2), LISA(nFrequencies)
  REAL(DP)            :: hc, z, tLb, r, M1, M2, f_ISCO
  REAL(DP)            :: z_arr(nRedshifts), tLb_z_arr(nRedshifts)
  INTEGER             :: nMergersPerSnapshot, iSnapshot, iMerger, iStrainFile, i
  CHARACTER(LEN=9)    :: FMTIN, FMTOUT
  CHARACTER(LEN=128)  :: FILEIN, FILEOUT, RootPath, &
                         LookbackTimeRedshiftFile, WriteFile
  LOGICAL             :: FileExists
  INTEGER             :: N
  REAL(DP) :: s
#ifdef _OPENMP
  INTEGER             :: iThread, nThreads, &
                         OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
#endif

  WRITE( RootPath, '(A)' ) '/Users/dunhamsj/Research/GW/Sam/'
!  WRITE( RootPath, '(A)' ) '/astro1/dunhamsj/'

  ! === Read in frequencies from LISA sensitivity curve data file ===
  OPEN( 100, FILE = TRIM(RootPath) // 'LISA_sensitivity.dat' )
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

  INQUIRE( FILE = TRIM(LookbackTimeRedshiftFile), EXIST = FileExists )
  IF( .NOT. FileExists )THEN
    WriteFile = 'Y'
  ELSE
    WRITE(*,'(A,A,A)') &
      'File ', TRIM(LookbackTimeRedshiftFile), ' exists. Overwrite? (Y/N) '
    READ(*,*) WriteFile
  END IF

  IF( TRIM(WriteFile) .EQ. 'Y' )THEN
    WRITE(*,'(A,A)') 'Writing file: ', TRIM(LookbackTimeRedshiftFile)
    OPEN( 100, FILE = TRIM(LookbackTimeRedshiftFile) )
    WRITE( 100, '(A)' ) '# Redshift z, Lookback-Time tLb [Gyr]'
    DO i = 1, nRedshifts
      z_arr    (i) = ( DBLE(i) - 1.0d0 ) * dz
      tLb_z_arr(i) = ComputeLookbackTime( z_arr(i) )
      WRITE( 100, '(ES23.16E3,1x,ES23.16E3)' ) z_arr(i), tLb_z_arr(i)
    END DO
    CLOSE( 100 )
  ELSE
    OPEN( 100, FILE = TRIM(LookbackTimeRedshiftFile) )
      READ( 100, * )
      DO i = 1, nRedshifts
        READ( 100, '(ES23.16E3,1x,ES23.16E3)' ) z_arr(i), tLb_z_arr(i)
      END DO
    CLOSE( 100 )
  END IF

  ! === Loop through Illustris snapshots ===

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP & PRIVATE(iSnapshot,iThread,iStrainFile,FMTIN,FMTOUT,FILEIN,FILEOUT, &
  !$OMP & FileExists,nMergersPerSnapshot,IllustrisData,M1,M2,tLb, &
  !$OMP & z,r,f_ISCO,i,hc) &
  !$OMP & SHARED(nThreads,RootPath,LISA,z_arr,tLb_z_arr)

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

    WRITE( FILEIN, FMTIN ) TRIM(RootPath) // &
                           'BBHM_DataFiles_Mapelli/time_BHillustris1_', &
                           iSnapshot, '.dat'

    INQUIRE( FILE = FILEIN, EXIST = FileExists )
    IF( .NOT. FileExists ) CYCLE

    WRITE( FILEOUT, FMTOUT ) &
      TRIM(RootPath) // 'strain/f90/StrainFiles/hc_', iSnapshot, '.dat'

    ! --- Create file for storing strains (first row will hold frequencies) ---
    OPEN ( iStrainFile, FILE = TRIM( FILEOUT ) )
    WRITE( iStrainFile, '(A,I3.3,A)' ) &
      '# M1 [Msun], M2 [Msun], Lookback-time to merger [Gyr],&
      & Redshift at merge, Comoving distance to merger [Mpc], fISCO,&
      & Characteristic strain (', nFrequencies, ') entries)'
    WRITE( iStrainFile, '(A)' ) '# (First row contains LISA frequencies [Hz])'
    WRITE( iStrainFile, '(ES12.6E2,ES13.6E2,ES18.11E2, &
                & ES18.11E2,E18.11E2,E18.11E2,10ES18.11E2)' ) &
                0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, LISA(:)

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
      r = ComputeComovingDistance( z )

      ! --- Calculate frequency at ISCO using Sesana et al., (2005) ---
      f_ISCO = c**3 / ( 6.0d0**( 3.0d0 / 2.0d0 ) * PI * G &
                 * ( M1 + M2 ) * Msun  * ( 1.0d0 + z ) )

      WRITE( iStrainFile, &
             '(ES12.6E2,ES13.6E2,ES18.11E2,ES18.11E2,ES18.11E2,ES18.11E2)', &
              ADVANCE = 'NO' ) M1, M2, tLb, z, r, f_ISCO
        
      ! --- Loop through frequencies until f_ISCO ---
      i = 1
      DO WHILE( ( LISA(i) .LT. f_ISCO ) .AND. ( i .LE. nFrequencies ) )
        hc = ComputeCharacteristicStrain( M1, M2, LISA(i), r )
        WRITE( iStrainFile, '(ES18.11E2)', ADVANCE = 'NO' ) hc
        i = i + 1
      END DO

      ! --- Fill in missing frequencies with 0.0d0 ---
      DO WHILE( i .LE. nFrequencies )
        WRITE( iStrainFile, '(f4.1)' , ADVANCE = 'NO' ) 0.0d0
        i = i + 1
      END DO

      WRITE( iStrainFile, * )
         
    END DO

    CLOSE( iSnapshot )

    CLOSE( iStrainFile )

  END DO
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL


CONTAINS


  ! --- Integrand E(z) in cosmological distance calculation ---
  PURE REAL(DP) FUNCTION E( z )

    REAL(DP), INTENT(in) :: z

    E = 1.0d0 / SQRT( OMEGA_M * ( 1.0d0 + z )**3 + OMEGA_L )

    RETURN
  END FUNCTION E


  FUNCTION ComputeComovingDistance( z ) RESULT( r )

    REAL(DP), INTENT(in) :: z
    REAL(DP)             :: r, dz
    INTEGER              :: k

    ! --- Integrate with Trapezoidal rule ---
    r  = 0.0d0
    dz = z / Nq
    
    DO k = 1, Nq - 1
       r = r + E( k * dz )
    END DO

    r = c / H0 / CentimetersPerMpc &
           * dz / 2.0d0 * ( E( 0.0d0 ) + 2.0d0 * r + E( z ) )

    RETURN
  END FUNCTION ComputeComovingDistance

  
  ! --- Integrand in lookback time calculation ---
  PURE REAL(DP) FUNCTION E_Lb( z )

    REAL(DP), INTENT(in) :: z

    E_Lb = 1.0d0 / ( ( 1.0d0 + z ) &
             * SQRT( OMEGA_M * ( 1.0d0 + z )**3 + OMEGA_L ) )
    
    RETURN
  END FUNCTION E_Lb


  ! --- Returns look-back time in Gyr ---
  REAL(DP) FUNCTION ComputeLookbackTime( z ) RESULT( tLb )

    REAL(DP), INTENT(in) :: z
    REAL(DP)             :: dz
    INTEGER              :: k

    dz = z / DBLE(Nq)

    ! --- Integrate with Trapezoidal rule ---
    tLb = 0.0d0
    DO k = 1, Nq - 1
      tLb = tLb + E_Lb( k * dz )
    END DO

    tLb = 1.0d0 / H0 / SecondsPerGyr &
            * dz / 2.0d0 * ( E_Lb( 0.0d0 ) + 2.0d0 * tLb + E_Lb( z ) )
    
    RETURN
  END FUNCTION ComputeLookbackTime

  
  ! --- Characteristic strain (dimensionless)
  !       from Sesana et al. (2005), Eq. (6) ---
  FUNCTION ComputeCharacteristicStrain( M1, M2, f, r ) RESULT( hc )

    REAL(DP) , INTENT(in)  :: M1, M2, f, r
    REAL(DP)               :: Mc
    REAL(DP)               :: hc, rcm

    rcm = r * CentimetersPerMpc

    ! --- Compute chirp mass (in source frame) ---
    Mc = ( M1 * M2 )**( 3.0d0 / 5.0d0 ) &
           / ( M1 + M2 )**( 1.0d0 / 5.0d0 ) * Msun

    hc = 1.0d0 / ( SQRT( 3.0d0 * c**3 ) * PI**( 2.0d0 / 3.0d0 ) * rcm ) &
           * ( G * Mc )**( 5.0d0 / 6.0d0 ) * f**( -1.0d0 / 6.0d0 )

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
      IF( DEBUG ) WRITE(*,*) 'Using small lookback-time approximation...'
      z = tLb * ( H0 * SecondsPerGyr )
      RETURN
    END IF

    ! --- Interpolate using linear interpolation ---
    IF( tLb .LT. tLb_z_arr(nRedshifts-1) ) THEN

      IF( DEBUG ) WRITE(*,*) 'Interpolating...'

      ! --- Get redshift bounds ---
      k = 1
      tLb_k = tLb_z_arr(k)
      DO WHILE( tLb_k .LT. tLb )
         k = k + 1
         tLb_k = tLb_z_arr(k)
      END DO

      zMin   = z_arr    (k-1)
      zMax   = z_arr    (k  )
      tLbMin = tLb_z_arr(k-1)
      tLbMax = tLb_z_arr(k  )

    ! --- Extrapolate using linear extrapolation ---
    ELSE

      IF( DEBUG ) WRITE(*,*) 'Extrapolating...'

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
      WRITE(*,*) 'tLb    = ', tLb
      WRITE(*,*) 'tLbMin = ', tLbMin
      WRITE(*,*) 'tLbMax = ', tLbMax
      WRITE(*,*)
    END IF

    RETURN
  END FUNCTION InterpolateLookbackTimeToGetRedshift


  PURE REAL(DP) FUNCTION Func(x) RESULT(y)
    REAL(DP), INTENT(in) :: x
    y = x**2
    RETURN
  END FUNCTION Func


END PROGRAM ComputeCharacteristicStrainForBHBmergers
