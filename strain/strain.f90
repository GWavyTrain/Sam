PROGRAM CharacteristicStrainFromBBHmergers

  IMPLICIT NONE

  INTEGER,  PARAMETER :: DP = KIND( 1.d0 )
  REAL(DP), PARAMETER :: PI = ACOS( -1.0d0 )

  ! --- Number of quadrature points for look-back time integration ---
  INTEGER,  PARAMETER :: Nq = 10000

  ! --- Resolution of redshift array for reference file ---
  INTEGER, PARAMETER :: Nz = 100000

  ! ---- Misc numbers ---
  INTEGER, PARAMETER :: Nc = 7, Nf_all = 400, Nf = 10, DELTA = Nf_all / Nf

  ! --- Mpc from https://www.cfa.harvard.edu/~dfabricant/huchra/ay145/constants.html ---
  REAL(DP), PARAMETER :: Mpc = 3.0857d24

  ! --- Column labels from Mapelli data ---
  INTEGER, PARAMETER :: iID = 6, iM1 = 7, iM2 = 8, itLB = 10

  ! --- Cosmological parameters come from Mapelli et al., (2017) ----
  REAL(DP), PARAMETER :: OMEGA_M = 0.2726d0, OMEGA_L = 0.7274d0, OMEGA_b = 0.0456d0
  REAL(DP), PARAMETER :: H0 = 70.4d0 * 1.0d5 * ( 1.0d0 / Mpc )

  ! --- Cutoff redshift above which lookback time will be extrapolated ---
  REAL(DP), PARAMETER :: z_max = 25.0d0, dz = z_max / ( Nz - 1 )

  ! --- Physical/unit constants ---
  REAL(DP), PARAMETER :: c = 2.99792458d10, G = 6.673d-8
  REAL(DP), PARAMETER :: Msun = 1.98892d33, Gyr = 1.0d9 * 86400.0d0 * 365.25d0

  ! --- 4 year LISA lifetime ---
  REAL(DP) :: tau = 4.0d0 * 86400.0d0 * 365.25d0

  ! --- Array to hold Illustris data ---
  REAL(DP), ALLOCATABLE :: IllustrisData(:,:)

  ! --- Arrays to hold LISA frequencies ---
  REAL(DP) :: LISA_all(Nf_all,2), LISA(Nf,2)

  ! --- Arrays to hold lookback-time/redshift reference file data ---
  REAL(DP) :: z_arr(Nz), tLb_z_arr(Nz)

  ! --- Misc variables ----
  REAL(DP) :: hc, z, tLb, r, M1, M2, f_ISCO
  INTEGER  :: Nss, SS(2), i, j, k

  ! --- Variables that require more memory ---
  INTEGER( KIND = 8 ) :: iMerger, nLinesIllustris, MergerID

  ! --- For file I/O ---
  INTEGER               :: hcFile = 100, IllustrisDataFile = 101, LISAFile = 102
  INTEGER               :: tLBzFile = 103, NssFile = 104, iArg
  CHARACTER( LEN = 28 ) :: FILEIN
  CHARACTER( LEN = 9 )  :: FMTIN
  CHARACTER( LEN = 3  ) :: arg

  ! --- Boolean for testing code ---
  LOGICAL :: TestStrain_Option = .TRUE.

  ! =================================================================
  !                        Start of program
  ! =================================================================

  ! --- Read in all frequencies from LISA sensitivity curve data file ---
  OPEN( LISAFile, FILE = '../LISA_sensitivity.dat', STATUS = 'OLD' )
  ! --- Loop through comments ---
  DO i = 1, Nc
    READ( LISAFile, * )
  END DO
  ! --- Read in data ---
  DO i = 1, Nf_all
    READ( LISAFile, * ) LISA_all(i,:)
  END DO
  CLOSE( LISAFile )

  ! --- Trim LISA data ---
  LISA(1,:) = LISA_all(1,:) ! --- Make sure first element is kept ---
  DO i = 2, Nf
     LISA(i,:) = LISA_all(i * DELTA,:)
  END DO

  ! --- Create file for storing strains (first row will hold frequencies) ---
  OPEN ( hcFile, FILE = 'hc.dat' )
  WRITE( hcFile, '(A)' ) &
         '# ID, M1, M2, tLb [Gyr], z, f_ISCO, r_comoving [Mpc], hc(f)'
  WRITE( hcFile, &
  '(I7,1x,ES18.10E3,1x,ES18.10E3,1x,ES18.10E3,1x,ES18.10E3,1x,ES18.10E3,1x,ES18.10E3)', &
         ADVANCE = 'NO' ) &
           0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0
  DO i = 1, Nf
     WRITE( hcFile, '(ES18.10E3)', ADVANCE = 'NO' ) LISA(i,1)
  END DO
  CLOSE( hcFile )

  IF( TestStrain_Option ) THEN

     WRITE(*,'(A)') 'Testing strain calculation...'

     ! --- GW150914 ---
!     M1 = 36.0d0
!     M2 = 29.0d0
!     z  = 0.09d0

     ! --- Sesana et al., (2005), ApJ, 623, 23, Fig. 2 ---
     MergerID = 120
     M1 = 1.0d5
     M2 = 0.1d0 * M1
     z  = 1.0d0
     tau = 3.0d0 * 86400.0d0 * 365.25d0 ! --- Lifetime of three years ---

     CALL main( MergerID, M1, M2, z, ComputeLookBackTime(z) )

     WRITE(*,'(A)') 'Stopping...'
     STOP

  END IF

  ! --- Create lookback time array ---
  OPEN( tLBzFile, FILE = 'tLb_z.dat' )
  WRITE( tLBzFile, '(A)' ) '# z, tLb [Gyr]'
  DO i = 1, Nz
    z_arr(i)     = ( i - 1 ) * dz
    tLb_z_arr(i) = ComputeLookbackTime( z_arr(i) ) / Gyr
    WRITE( tLBzFile, '(ES16.10,1x,F13.10)' ) z_arr(i), tLb_z_arr(i)
  END DO
  CLOSE( tLBzFile )

  ! --- Get upper and lower snapshot numbers ---
  IF ( IARGC() .NE. 2 ) THEN
    WRITE( *, '(A)' ) 'Proper usage: ./strain SSL SSU'
    WRITE( *, '(A)' ) 'Exiting...'
    STOP
  END IF

  ! --- Read upper and lower snapshot numbers into array ---
  DO iArg = 1, IARGC()
    CALL GETARG( iArg, arg )
    READ( arg, * ) SS(iArg)
  END DO

  ! --- Loop through Illustris snapshots ---
  DO Nss = SS(1), SS(2)

    ! --- These two files don't exist ---
    IF ( ( Nss .NE. 53 ) .AND. ( Nss .NE. 55 ) ) THEN

      ! --- Save file with snapshot number currently being analyzed ---    
      OPEN ( NssFile, FILE = 'Nss.dat' )
      WRITE( NssFile, '(I3)' ) Nss
      CLOSE( NssFile )

      ! --- Get filenames ---
      IF ( Nss < 100 ) THEN
        FMTIN = '(A,I2,A4)'
      ELSE
        FMTIN = '(A,I3,A4)'
      END IF
      WRITE( FILEIN, FMTIN  ) &
           '../time_BHillustris1_', Nss, '.dat'
!           '../BBHM_DataFiles_Mapelli/time_BHillustris1_', Nss, '.dat'

      ! --- Get number of lines (mergers) in Illustris data file ---
      nLinesIllustris = 0
      OPEN( IllustrisDataFile, FILE = TRIM( FILEIN ) )
      DO
        READ( IllustrisDataFile, *, END = 11 )
        nLinesIllustris = nLinesIllustris + 1
      END DO
      11 CLOSE( IllustrisDataFile )

      ! --- Read in Illustris data and compute strain for each merger ---
      ALLOCATE( IllustrisData( nLinesIllustris, 10 ) )
      OPEN( IllustrisDataFile, FILE = TRIM( FILEIN ) )
      DO iMerger = 1, nLinesIllustris
        READ( IllustrisDataFile, * ) IllustrisData(iMerger,:)

        ! --- Get data for this merger ---
        MergerID = INT8( IllustrisData(iMerger,iID) )
        M1       = IllustrisData(iMerger,iM1)  ! [ Solar masses, source frame ]
        M2       = IllustrisData(iMerger,iM2)  ! [ Solar masses, source frame ]
        tLb      = IllustrisData(iMerger,itLb) ! [ Gyr ]

        ! --- Compute redshift from lookback time ---
        z = InterpolateLookbackTime( tLb )

        CALL main( MergerID, M1, M2, z, tLb )
        
      END DO ! --- Loop over mergers ---

      CLOSE( IllustrisDataFile )
      DEALLOCATE( IllustrisData )

    END IF ! --- Not using 53 or 55 ---

  END DO ! --- Loop over all snapshots ---


  ! =================================================================
  !                          End of program
  ! =================================================================


CONTAINS

  ! --- Integrand E(z) in cosmological distance calculation ---
  PURE FUNCTION E( z )

    REAL(DP), INTENT(in) :: z
    REAL(DP)             :: E

    E = 1.0d0 / SQRT( ( OMEGA_M + OMEGA_b ) * ( 1.0d0 + z )**3 + OMEGA_L )

    RETURN
  END FUNCTION E

  
  FUNCTION ComputeComovingDistance( z ) RESULT( r )

    REAL(DP), INTENT(in) :: z
    REAL(DP)             :: r, dz

    ! --- Integrate with Trapezoidal rule ---
    r  = 0.0d0
    dz = z / Nq
    
    DO k = 1, Nq - 1
       r = r + E( k * dz )
    END DO

    r = c / H0 * dz / 2.0d0 * ( E( 0.0d0 ) + 2.0d0 * r + E( z ) )

    RETURN
  END FUNCTION ComputeComovingDistance

  
  ! --- Integrand in lookback time calculation ---
  PURE FUNCTION E_LB( z )

    REAL(DP), INTENT(in) :: z
    REAL(DP)             :: E_LB

    E_LB = 1.0d0 / ( ( 1.0d0 + z ) &
             * SQRT( ( OMEGA_M + OMEGA_b ) * ( 1.0d0 + z )**3 + OMEGA_L ) )
    
    RETURN
  END FUNCTION E_LB

  
  FUNCTION ComputeLookbackTime( z ) RESULT( tLB )

    REAL(DP), INTENT(in) :: z
    REAL(DP)             :: dz, tLB
    INTEGER              :: i
    
    ! --- Integrate with Trapezoidal rule ---
    ! --- |ERROR|: Nq = 1000   --> 1.47d-3
    !              Nq = 10000  --> 1.47d-5
    !              Nq = 100000 --> 1.47d-7 ---
    
    tLB = 0.0d0
    dz  = z / Nq

    tLB = 0.0d0
    DO i = 1, Nq - 1
      tLB = tLB + E_LB( i * dz )
    END DO

    tLB = 1.0d0 / H0 * dz / 2.0d0 * ( E_LB( 0.0d0 ) + 2.0d0 * tLB + E_LB( z ) )

    RETURN
  END FUNCTION ComputeLookbackTime

  
  ! --- Characteristic strain (dimensionless)
  !       from Sesana et al. (2005), ApJ, 623, 23: Eqs. (6,7) ---
  FUNCTION ComputeCharacteristicStrain( M1, M2, f, z, r ) RESULT( hc )

    REAL(DP), INTENT(in)  :: M1, M2, f, z, r
    REAL(DP)              :: hc
    REAL(DP)              :: Mc, h, n

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
  END FUNCTION ComputeCharacteristicStrain

  
  ! --- Interpolate lookback time array to get redshift ---
  FUNCTION InterpolateLookbackTime( tLB ) RESULT( z )

    ! --- Get redshift by interpolating lookback time ---
    ! --- tLB should be in units of Gyr ---

    REAL(DP), INTENT(in) :: tLB
    REAL(DP)             :: zmin, zmax, tLBmin, tLBmax, tLBi, m, b
    REAL(DP)             :: z
    INTEGER              :: i

    ! --- Small lookback time approximation: tLB ~ tH * z ---
    ! --- 0.1 chosen to ensure 1% error in distances to match
    !       LISA's predicted capability (Ask Kelly for reference to this) ---
    IF ( tLB < 0.1d0 ) THEN
      z = tLB / ( ( 1.0d0 / H0 ) / Gyr )
      RETURN
    END IF

    ! --- If target time is in the data set,
    !       interpolate using linear interpolation ---
    IF ( tLB < tLB_z_arr(Nz) ) THEN

      ! --- Get redshift bounds ---
      i = 1
      tLBi = tLB_z_arr(i)
      DO WHILE ( tLBi < tLB )
        i = i + 1
        tLBi = tLB_z_arr(i)
      END DO

      zmin   = z_arr(i-1)
      zmax   = z_arr(i)
      tLBmin = tLB_z_arr(i-1)
      tLBmax = tLB_z_arr(i)

    ! --- If target time is not in the data set,
    !       extrapolate using linear extrapolation ---
    ELSE

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


  SUBROUTINE main( MergerID, M1, M2, z, tLb )

    INTEGER( KIND = 8 ), INTENT(in) :: MergerID
    REAL(DP),            INTENT(in) :: M1, M2, z, tLb

    ! --- Compute comoving distance (in cm) ---
    r = ComputeComovingDistance( z ) ! [ cm ]

    ! --- Calculate frequency at ISCO using
    !       Sesana et al. (2005), Eq. (12) ---
    f_ISCO = c**3 / ( 6.0d0**( 3.0d0 / 2.0d0 ) * PI * G &
               * ( M1 + M2 ) * Msun  * ( 1.0d0 + z ) )

    ! --- Start writing to file ---
    OPEN ( hcFile, FILE = 'hc.dat', POSITION = 'APPEND' )
    WRITE( hcFile, &
   '(I7,1x,ES18.10E3,1x,ES18.10E3,1x,ES18.10E3,1x,ES18.10E3,1x,ES18.10E3,1x,ES18.10E3)', &
             ADVANCE = 'NO' ) MergerID, M1, M2, tLB, z, f_ISCO, r / Mpc
        
    ! --- Loop through frequencies until f_ISCO and write to file ---
    j = 1
    DO WHILE ( LISA(j,1) .LT. f_ISCO )
      hc = ComputeCharacteristicStrain( M1, M2, LISA(j,1) ,z ,r )
      WRITE( hcFile, '(ES13.6)', ADVANCE = 'NO' ) hc
      j = j + 1
      IF( j .GT. Nf ) EXIT
    END DO

    ! --- Fill in missing frequencies with nan and write to file ---
    DO WHILE ( j .LT. Nf + 1 )
      WRITE( hcFile, '(A4)', ADVANCE = 'NO' ) ' nan'
      j = j + 1
    END DO

    WRITE( hcFile, * )
    CLOSE( hcFile )

    RETURN
  END SUBROUTINE main

END PROGRAM CharacteristicStrainFromBBHmergers
