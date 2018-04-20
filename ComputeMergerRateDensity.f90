PROGRAM ComputeMergerRate

  IMPLICIT NONE

  INTEGER  , PARAMETER   :: DP = KIND( 1.d0 )
  INTEGER  , PARAMETER   :: Nq = 10000 , Nz = 100000
  INTEGER  , PARAMETER   :: itLB = 4
  REAL(DP) , PARAMETER   :: PI = ACOS( -1.0d0 )
  REAL(DP) , PARAMETER   :: z_max = 25.0d0 , dz = z_max / ( Nz - 1 )
  REAL(DP) , PARAMETER   :: OMEGA_M = 0.3d0 , OMEGA_L = 0.7d0
  REAL(DP) , PARAMETER   :: H0 = 70.4d0 / 3.086d19
  REAL(DP) , PARAMETER   :: c = 2.99792458d10 , dt = 0.01d0 , Gpc = 3.086d24
  REAL(DP) , PARAMETER   :: year = 86400.0d0 * 365.25d0 , Gyr = 1.0d9 * year
  REAL(DP) , ALLOCATABLE :: BinEdges(:) , Bins(:)
  REAL(DP) , ALLOCATABLE :: hc(:,:) , tLB(:) , Counts(:)
  REAL(DP)               :: tShort , tLong , Vc , z_arr(Nz) , tLB_z_arr(Nz)
  INTEGER                :: i , j
  INTEGER( KIND = 8 )    :: nMergers , nBins

  ! --- Create lookback time array in Gyr ---
  WRITE(*,'(A)') 'Creating lookback time array'
  OPEN( 100 , FILE = 'tLB_z.dat' )
  ! --- Compute array of lookback times in Gyr ---
  WRITE( 100 , '(A14)' ) '# z, tLB [Gyr]'
  DO i = 1 , Nz
    z_arr(i)     = ( i - 1 ) * dz
    tLB_z_arr(i) = ComputeLookbackTime( z_arr(i) )
    WRITE( 100 , '(ES16.10,1x,F13.10)' ) z_arr(i) , tLB_z_arr(i)
  END DO
  CLOSE( 100 )
  
  ! --- Get number of mergers in hc.dat file ---
  WRITE(*,'(A)') 'Getting number of mergers in hc.dat'
  OPEN( 100 , FILE = 'hc.dat' , STATUS = 'OLD' )
  READ( 100 , * ) ! --- Skip the header ---
  READ( 100 , * ) ! --- Skip the line with frequencies ---  
  nMergers = 0
  DO
    READ( 100 , * , END = 10 )
    nMergers = nMergers + 1
  END DO
  10 CLOSE( 100 )
  WRITE(*,'(A10,I19)') 'nMergers: ' , nMergers

  ! --- Read in hc file and get merger lookback times in Gyr for each ---
  WRITE(*,'(A)') 'Getting merger lookback times'
  ALLOCATE( hc ( nMergers , 17 ) )
  ALLOCATE( tLB( nMergers ) )
  OPEN( 100 , FILE = 'hc.dat' , STATUS = 'OLD' )
  READ( 100 , * ) ! --- Skip the header ---
  READ( 100 , * ) ! --- Skip the line with frequencies ---
  DO i = 1 , nMergers
    READ( 100 , * ) hc( i , : )
    tLB(i) = hc( i , itLB )
  END DO
  CLOSE( 100 )

  tShort = MINVAL( tLB )
  tLong  = MAXVAL( tLB )

  WRITE(*,'(A22,ES17.10E3,A4)') 'Shortest merger time: ' , tShort , ' Gyr'
  WRITE(*,'(A22,ES17.10E3,A4)') 'Longest merger time:  ' , tLong  , ' Gyr'

  ! --- Create lookback time bins/bin edges in Gyr ---
  WRITE(*,'(A)') 'Creating lookback time bins/bin edges'
  tShort    = tShort
  tLong     = tLong  + dt
  nBins = CEILING( ( tLong - tShort ) / dt )
  ALLOCATE( Bins(nBins) )
  ALLOCATE( BinEdges(nBins+1) )
  BinEdges = CreateLookbackTimeBinEdges( nBins )
  WRITE(*,'(A7,I19)') 'nBins: ' , nBins

  ! --- Create and initialize array to hold counts ---
  ALLOCATE( Counts(nBins) )
  Counts = 0.0d0

  ! --- Populate array with counts ---
  WRITE(*,'(A)') 'Populating counts array'
  OPEN( 100 , FILE = 'MergerRateDensity.dat' )
  WRITE( 100 , '(A60)' ) &
       '# bin centers [Gyr] , Merger Rate Density [Gpc^(-3)*yr^(-1)]'
  CLOSE( 100 )
  DO i = 1 , nBins
    OPEN( 100 , FILE = 'MergerRateDensity.dat' , POSITION = 'APPEND' )
    DO j = 1 , nMergers
       IF ( ( tLB(j) >= BinEdges(i)   ) .AND. &
            ( tLB(j) <  BinEdges(i+1) ) ) Counts(i) = Counts(i) + 1
    END DO
    Bins(i) = ( BinEdges(i) + BinEdges(i+1) ) / 2.0d0
    ! --- Convert merger rate to merger rate density ---
    Vc = ComputeComovingVolumeFromLookbackTime( Bins(i) )
    WRITE( 100 , '(ES18.10E3,1x,ES18.10E3)' ) Bins(i) , Counts(i) / Vc / 1.0d9
    CLOSE( 100 )
  END DO

  DEALLOCATE( Counts   )
  DEALLOCATE( BinEdges )
  DEALLOCATE( Bins     )
  DEALLOCATE( tLB      )
  DEALLOCATE( hc       )

CONTAINS

  FUNCTION CreateLookbackTimeBinEdges( nBins ) RESULT( BinEdges )

    INTEGER( KIND = 8 ) , INTENT(in) :: nBins
    REAL(DP)                         :: BinEdges(nBins+1)
    INTEGER                          :: i

    DO i = 1 , nBins + 1
      BinEdges(i) = tShort + ( i - 1 ) * dt
    END DO

    RETURN
  END FUNCTION CreateLookbackTimeBinEdges


  ! --- Integrand for lookback time calculation ---
  PURE FUNCTION OneOverE_LB( z )

    REAL(DP) , INTENT(in) :: z
    REAL(DP)              :: OneOverE_LB

    OneOverE_LB = 1.0d0 / ( ( 1.0d0 + z ) &
                    * SQRT( OMEGA_M * ( 1.0d0 + z )**3 + OMEGA_L ) )
    
    RETURN
  END FUNCTION OneOverE_LB

  
  FUNCTION ComputeLookbackTime( z ) RESULT( tLB )

    REAL(DP) , INTENT(in) :: z
    REAL(DP)              :: dz , tLB
    INTEGER               :: i
    
    dz  = z / Nq

    tLB = 0.0d0
    DO i = 1 , Nq - 1
      tLB = tLB + OneOverE_LB( i * dz )
    END DO

    tLB = 1.0d0 / H0 * dz / 2.0d0 &
            * ( OneOverE_LB( 0.0d0 ) + 2.0d0 * tLB + OneOverE_LB( z ) ) / Gyr
    
    RETURN
  END FUNCTION ComputeLookbackTime
  

  ! --- Interpolate lookback time to get redshift ---
  FUNCTION InterpolateLookbackTime( tLB ) RESULT( z )

    REAL(DP) , INTENT(in) :: tLB
    REAL(DP)              :: zmin , zmax , tLBmin , tLBmax , tLBk , m , b
    REAL(DP)              :: z
    INTEGER               :: k

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

  ! --- Integrand for comoving distance calculation ---
  PURE FUNCTION E( z )

    REAL(DP) , INTENT(in) :: z
    REAL(DP)              :: E

    E = 1.0d0 / SQRT( OMEGA_M * ( 1.0d0 + z )**3 + OMEGA_L )

    RETURN
  END FUNCTION E

  
  FUNCTION ComputeComovingVolumeFromLookbackTime( tLB ) RESULT( Vc )

    REAL(DP) , INTENT(in) :: tLB
    REAL(DP)              :: Dc , dz , Vc , z
    INTEGER               :: k

    ! --- Get redshift for given lookback time ---
    z = InterpolateLookbackTime( tLB )
    
    ! --- Integrate with Trapezoidal rule ---
    Dc = 0.0d0
    dz = z / Nq
    
    DO k = 1 , Nq - 1
       Dc = Dc + E( k * dz )
    END DO

    Dc = c / H0 * dz / 2.0d0 * ( E( 0.0d0 ) + 2.0d0 * Dc + E( z ) )

    Vc = 4.0d0 / 3.0d0 * PI * ( Dc / Gpc )**3

    RETURN
  END FUNCTION ComputeComovingVolumeFromLookbackTime


END PROGRAM ComputeMergerRate
