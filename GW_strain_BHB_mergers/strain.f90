PROGRAM GWstrainFromBHBmergers

  IMPLICIT NONE

  INTEGER  , PARAMETER   :: DP = KIND( 1.d0 ) , Nq = 1000 , Nf = 900
  REAL(DP) , PARAMETER   :: PI = ACOS( -1.0d0 )
  REAL(DP) , PARAMETER   :: OMEGA_M = 0.3d0 , OMEGA_L = 0.7d0
  REAL(DP) , PARAMETER   :: c = 3.0d10 , H0 = 70.4d0 / 3.086d19 , G = 6.67d-8
  REAL(DP) , PARAMETER   :: Msun = 2.0d33 , year = 3.1536d7
  REAL(DP) , PARAMETER   :: Tobs = 10.0d0 * year
  REAL(DP) , ALLOCATABLE :: IllustrisData(:,:)
  REAL(DP)               :: D , LISA(Nf,2)
  INTEGER                :: nLinesIllustris , Nz , i , j
  INTEGER                :: iz = 3 , iM1 = 7 , iM2 = 8
  CHARACTER( len = 25 )  :: FILEIN
  CHARACTER( len = 17 )  :: FILEOUT
  CHARACTER( len = 11 )  :: FMTIN
  CHARACTER( len = 11 )  :: FMTOUT

  ! --- Read in frequencies from LISA sensitivity curve ---
  OPEN( 100 , FILE = 'LISA_sensitivity.dat' )
  DO i = 1 , Nf
    READ( 100 , * , END = 10 ) LISA( i , : )
  END DO
  10 CLOSE( 100 )

  ! --- Loop through Illustris snapshots ---
  DO Nz = 26 , 26!135

    ! --- Get filenames ---
    IF ( Nz < 100 ) THEN
      FMTIN  = '(A18,I2,A4)'
      FMTOUT = '(A10,I2,A4)'
    ELSE
      FMTIN  = '(A18,I3,A4)'
      FMTOUT = '(A10,I3,A4)'
    END IF
<<<<<<< HEAD
    WRITE( FILEIN  , FMTIN  ) 'time_BHillustris1_' , Nz , '.dat'
    WRITE( FILEOUT , FMTOUT ) 'GW_strain_'         , Nz , '.dat'
=======
    WRITE( TRIM( FILEIN  ) , FMTIN  ) 'time_BHillustris1_' , Nz , '.dat'
    WRITE( TRIM( FILEOUT ) , FMTOUT ) 'GW_strain_'         , Nz , '.dat'
>>>>>>> 80d4433ea92847fe27ed675197ebbdf7995db384

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

    ! --- Calculate luminosity distance ---
    D = ComputeLuminosityDistance( IllustrisData( 1 , iz ) )

    ! --- Loop through frequencies and calculate strain ---
    OPEN( 103 , FILE = TRIM( FILEOUT ) )
    DO i = 1 , Nf
      WRITE( 103 , * ) LISA( i , 1 ) , hf( IllustrisData( : , iM1 ) , &
                         IllustrisData( : , iM2 ) , LISA( i , 1 ) )
    END DO
    CLOSE( 103 )

    DEALLOCATE( IllustrisData )

  END DO

CONTAINS

  ! --- E(z) in cosmological distance calculation ---
  PURE FUNCTION E( z )

    REAL(DP) , INTENT(in) :: z
    REAL(DP)              :: E

    E = SQRT( OMEGA_M * ( 1.0d0 + z )**3 + OMEGA_L )

  END FUNCTION E

  FUNCTION ComputeLuminosityDistance( z ) RESULT( D_L )

    REAL(DP) , INTENT(in) :: z
    REAL(DP)              :: D_L , dz

    ! --- Integrate with Trapezoidal rule ---
    D_L = 0.0d0
    dz = z / Nq
    
    DO i = 1 , Nq - 1
       D_L = D_L + 1.0d0 / E( i * dz )
    END DO

    D_L = c / H0 * dz / 2.0d0 * ( E( 0.0d0 ) + 2.0d0 * D_L + E( z ) ) &
            * ( 1.0d0 + z )

    RETURN
  END FUNCTION ComputeLuminosityDistance
  
  ! --- Monochromatic strain from lecture notes ---
  FUNCTION hf( m1 , m2 , f )

    REAL(DP) , INTENT(inout) :: m1(nLinesIllustris) , m2(nLinesIllustris)
    REAL(DP) , INTENT(in)    :: f
    REAL(DP)                 :: Mc(nLinesIllustris) , hf(nLinesIllustris)

    ! --- Compute chirp mass ---
    Mc = ( m1 * m2 )**( 3.0d0 / 5.0d0 ) / ( m1 + m2 )**( 1.0d0 / 5.0d0 ) * Msun

    hf =  G / c**2 * Mc / D &
            * ( G / c**3 * PI * f * Mc )**( 2.0d0 / 3.0d0 ) &
              * SQRT( Tobs )

    RETURN
  END FUNCTION hf
  
END PROGRAM GWstrainFromBHBmergers
