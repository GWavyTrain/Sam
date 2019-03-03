PROGRAM ComputeCharacteristicStrainForBHBmergers

  USE UtilitiesModule

  IMPLICIT NONE

  INTEGER,  PARAMETER   :: iSnapshotMin = 26, iSnapshotMax = 26
  INTEGER,  PARAMETER   :: nFrequenciesAll = 400, nFrequencies = 10, &
                           nSkip = nFrequenciesAll / nFrequencies, &
                           nCommentLines = 9
  INTEGER,  PARAMETER   :: iM1 = 7, iM2 = 8, itLb = 10
  REAL(DP), PARAMETER   :: z_max = 20.0d0, &
                           dz = z_max / ( DBLE(nRedshifts) - 1.0d0 )
  REAL(DP), PARAMETER   :: CentimetersPerMpc = 3.086d24
  REAL(DP)              :: IllustrisData(10)
  INTEGER               :: nMergersPerSnapshot, iSnapshot, &
                           iMerger, iStrainFile
  CHARACTER(LEN=9)      :: FMTIN, FMTOUT
  CHARACTER(LEN=128)    :: FILEIN, FILEOUT, RootPath, &
                           LookbackTimeRedshiftFile, WriteFile, FMThc
  LOGICAL               :: FileExists
  REAL(DP), ALLOCATABLE :: LISA_all(:,:), LISA(:), z_arr(:), tLb_z_arr(:)

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

      ! --- Fill in missing frequencies with 0 ---
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


END PROGRAM ComputeCharacteristicStrainForBHBmergers
