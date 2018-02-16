PROGRAM GWstrainFromBHBmergers

  IMPLICIT NONE

  INTEGER  , PARAMETER :: DP = KIND( 1.d0 ) , Nq = 1000000
  INTEGER  , PARAMETER :: Nf = 900 , Nc = 22
  REAL(DP) , PARAMETER :: PI = ACOS( -1.0d0 )
  REAL(DP) , PARAMETER :: OMEGA_M = 0.3d0 , OMEGA_L = 0.7d0
  REAL(DP) , PARAMETER :: c = 3.0d10 , G = 6.67d-8 , H0 = 70.4d0 / 3.086d19
  REAL(DP) , PARAMETER :: Msun = 2.0d33
  REAL(DP)             :: D , D_L , LISA(Nf,2) , hc , R_ISCO , f_ISCO
  REAL(DP)             :: M1 = 1.0d7 , M2 = 1.0d7 , z = 3.0d0
  INTEGER              :: i

  ! --- Read in frequencies from LISA sensitivity curve ---
  OPEN( 100 , FILE = 'LISA_sensitivity.dat' )

  ! --- Loop through comments
  DO i = 1 , Nc
    READ( 100 , * )
  END DO
 
  DO i = 1 , Nf
    READ( 100 , * ) LISA( i , : )
  END DO
 
  CLOSE( 100 )

  ! --- Calculate comoving distance ---
  D = ComputeComovingDistance( z )

  ! --- Calculate the frequency at the ISCO
  R_ISCO = 6.0d0 * G * M1 * Msun / c**2
  f_ISCO = 1.0d0 / ( 2.0d0 * PI ) * SQRT( G * Msun * ( M1 + M2 ) / R_ISCO**3 )
  
  ! --- Loop through frequencies and calculate strain ---
  OPEN(  103 , FILE = 'GW_strain_test_GOAT.dat' )
  WRITE( 103 , '(A15)' ) '# frequency, hc'
  WRITE( 103 , '(F6.3,1x,I6)' ) z , INT( D / 3.086d24 )
  i = 1
  DO WHILE ( LISA( i , 1 ) < f_ISCO )
    hc = Strain( M1 , M2 , LISA( i , 1 ) , z )
    WRITE( 103 , '(E11.5,1x,E11.5)' ) LISA( i , 1 ) , hc
    i = i + 1
  END DO
  CLOSE( 103 )

CONTAINS

  ! --- E(z) in cosmological distance calculation ---
  PURE FUNCTION E( z )

    REAL(DP) , INTENT(in) :: z
    REAL(DP)              :: E

    E = SQRT( OMEGA_M * ( 1.0d0 + z )**3 + OMEGA_L )

  END FUNCTION E

  FUNCTION ComputeComovingDistance( z ) RESULT( D_C )

    REAL(DP) , INTENT(in) :: z
    REAL(DP)              :: D_C , dz

    ! --- Integrate with Trapezoidal rule ---
    D_C = 0.0d0
    dz = z / Nq
    
    DO i = 1 , Nq - 1
       D_C = D_C + 1.0d0 / E( i * dz )
    END DO

    D_C = c / H0 * dz / 2.0d0 * ( E( 0.0d0 ) + 2.0d0 * D_C + E( z ) )

    RETURN
  END FUNCTION ComputeComovingDistance
  
  ! --- Characteristic strain (dimensionless) from Sesana (2016) ---
  PURE FUNCTION Strain( M1 , M2 , f , z ) RESULT( hc )

    REAL(DP) , INTENT(in) :: M1 , M2 , f , z
    REAL(DP)              :: hc , Mc

    ! --- Compute chirp mass (in source frame) ---
    Mc = ( M1 * M2 )**( 3.0d0 / 5.0d0 ) / ( M1 + M2 )**( 1.0d0 / 5.0d0 ) * Msun

    hc = 1.0d0 / ( PI * D )                                       &
           * ( 2.0d0 * PI / ( 3.0d0 * c**3 ) )**( 1.0d0 / 2.0d0 ) &
             * ( G * Mc )**( 5.0d0 / 6.0d0 )                      &
               * ( PI * f )**( -1.0d0 / 6.0d0 )
    RETURN
  END FUNCTION Strain
  
END PROGRAM GWstrainFromBHBmergers
