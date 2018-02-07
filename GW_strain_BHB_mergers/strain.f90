PROGRAM GWstrainFromBHBmergers

  IMPLICIT NONE

  INTEGER  , PARAMETER   :: DP = KIND( 1.d0 ) , Nq = 1000000 , Nf = 900
  INTEGER  , PARAMETER   :: iz = 3 , iM1 = 7 , iM2 = 8
  REAL(DP) , PARAMETER   :: PI = ACOS( -1.0d0 )
  REAL(DP) , PARAMETER   :: OMEGA_M = 0.3d0 , OMEGA_L = 0.7d0
  REAL(DP) , PARAMETER   :: c = 3.0d10 , G = 6.67d-8 , H0 = 70.4d0 / 3.086d19
  REAL(DP) , PARAMETER   :: Msun = 2.0d33 , year = 3.1536d7
  REAL(DP) , PARAMETER   :: Tobs = 10.0d0 * year
  REAL(DP) , ALLOCATABLE :: IllustrisData(:,:)
  REAL(DP)               :: D , LISA(Nf,2) , hf , hf_min , hf_max
  INTEGER                :: nLinesIllustris , Nz , i , j
  CHARACTER( len = 25 )  :: FILEIN
  CHARACTER( len = 17 )  :: FILEOUT
  CHARACTER( len = 11 )  :: FMTIN , FMTOUT
  CHARACTER( len = 34 )  :: FMT

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

    ! --- Calculate luminosity distance ---
    D = ComputeLuminosityDistance( IllustrisData( 1 , iz ) )

    ! --- Loop through frequencies and calculate mean strain ---
    OPEN( 103 , FILE = TRIM( FILEOUT ) )
    WRITE( 103 , '(A36)' ) '# frequency, hf_mean, hf_min, hf_max'
    WRITE( 103 , '(F6.3,1x,I6,1x,I1,1x,I1 )' ) &
      IllustrisData( 1 , iz ) , INT( D / 3.086d24 ) , 0 , 0
    FMT = '(E11.5,1x,E11.5,1x,E11.5,1x,E11.5)'
    DO i = 1 , Nf
      CALL Strain( IllustrisData( : , iM1 ) ,                   &
                     IllustrisData( : , iM2 ) , LISA( i , 1 ) , &
                       hf , hf_min , hf_max )
      WRITE( 103 , FMT ) LISA( i , 1 ) , hf , hf_min , hf_max
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
  SUBROUTINE Strain( M1 , M2 , f , hf , hf_min , hf_max )

    REAL(DP) , INTENT(in)  :: M1(nLinesIllustris) , M2(nLinesIllustris) , f
    REAL(DP) , INTENT(out) :: hf , hf_min , hf_max
    REAL(DP)               :: Mc(nLinesIllustris) , hf_all(nLinesIllustris)

    ! --- Compute chirp mass ---
    Mc = ( M1 * M2 )**( 3.0d0 / 5.0d0 ) / ( M1 + M2 )**( 1.0d0 / 5.0d0 ) * Msun

    hf_all = G / c**2 * Mc / D                                 &
               * ( G / c**3 * PI * f * Mc )**( 2.0d0 / 3.0d0 ) &
                 * SQRT( Tobs )

    hf     = SUM( hf_all ) / SIZE( hf_all )
    hf_min = MINVAL( hf_all )
    hf_max = MAXVAL( hf_all )

    RETURN
  END SUBROUTINE Strain
  
END PROGRAM GWstrainFromBHBmergers
