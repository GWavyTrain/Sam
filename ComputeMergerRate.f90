PROGRAM ComputeMergerRate

  IMPLICIT NONE

  INTEGER  , PARAMETER   :: DP = KIND( 1.d0 )
  INTEGER  , PARAMETER   :: Nq = 10000
  INTEGER  , PARAMETER   :: itLB = 4 , iz = 5
  REAL(DP) , PARAMETER   :: PI = ACOS( -1.0d0 )
  REAL(DP) , PARAMETER   :: OMEGA_M = 0.3d0 , OMEGA_L = 0.7d0
  REAL(DP) , PARAMETER   :: H0 = 70.4d0 / 3.086d19
  REAL(DP) , PARAMETER   :: c = 2.99792458d10 , dt = 0.01d0
  REAL(DP) , PARAMETER   :: year = 86400.0d0 * 365.25d0 , Gyr = 1.0d9 * year
  REAL(DP) , ALLOCATABLE :: BinEdges(:) , Bins(:)
  REAL(DP) , ALLOCATABLE :: hc(:,:) , data(:,:) , Counts(:)
  REAL(DP)               :: tShort , tLong , zMin , zMax
  INTEGER                :: i , j , nMergers , nBins

  ! --- Get number of mergers in hc.dat file ---
  OPEN( 100 , FILE = 'hc.dat' )
  READ( 100 , * ) ! --- Skip the header ---
  READ( 100 , * ) ! --- Skip the line with frequencies ---  
  nMergers = 0
  DO
    READ( 100 , * , END = 10 )
    nMergers = nMergers + 1
  END DO
  10 CLOSE( 100 )

  ! --- Read in hc file and get merger times and redshifts for each ---
  ALLOCATE( hc  ( nMergers , 17 ) )
  ALLOCATE( data( nMergers , 3  ) )
  OPEN( 100 , FILE = 'hc.dat' )
  READ( 100 , * ) ! --- Skip the header ---
  READ( 100 , * ) ! --- Skip the line with frequencies ---
  DO i = 1 , nMergers
    READ( 100 , * ) hc( i , : )
    data(i,1) = hc( i , itLB )                     ! tLB [ Gyr ]
    data(i,2) = hc( i , iz   )                     ! z   [ dimensionless ]
    data(i,3) = ComputeComovingVolume( data(i,2) ) ! Vc  [ cm^3 ]
  END DO
  CLOSE( 100 )

  tShort = MINVAL( data(:,1) )
  zMin   = MINVAL( data(:,2) )
  tLong  = MAXVAL( data(:,1) )
  zMax   = MAXVAL( data(:,2) )

  tShort = tShort
  tLong  = tLong  + dt

  ! --- Create lookback time bins/bin edges ---
  nBins = CEILING( ( tLong - tShort ) / dt )
  ALLOCATE( Bins(nBins) )
  ALLOCATE( BinEdges(nBins+1) )

  BinEdges = CreateLookbackTimeBinEdges( tShort , nBins )

  ALLOCATE( Counts(nBins) )
  Counts = 0.0d0

  ! --- Populate array with counts ---
  OPEN( 100 , FILE = 'MergerRateCounts.dat' )
  WRITE( 100 , '(A28)' ) '# bin centers [Gyr] , counts'
  DO i = 1 , nBins
    DO j = 1 , nMergers
       IF ( ( data(j,1) >= BinEdges(i)   ) .AND. &
            ( data(j,1) <  BinEdges(i+1) ) ) Counts(i) = Counts(i) + 1
    END DO
    Bins(i) = ( BinEdges(i) + BinEdges(i+1) ) / 2.0d0
    WRITE( 100 , '(ES18.10E3,1x,ES16.10E2)' ) Bins(i) , Counts(i)
  END DO
  CLOSE( 100 )
  
  DEALLOCATE( Counts   )
  DEALLOCATE( BinEdges )
  DEALLOCATE( Bins     )
  DEALLOCATE( data     )
  DEALLOCATE( hc       )

CONTAINS

  ! --- Integrand 1/E(z) in cosmological distance calculation ---
  PURE FUNCTION OneOverE( z )

    REAL(DP) , INTENT(in) :: z
    REAL(DP)              :: OneOverE

    OneOverE = 1.0d0 / SQRT( OMEGA_M * ( 1.0d0 + z )**3 + OMEGA_L )

    RETURN
  END FUNCTION OneOverE

  
  ! --- Compute comoving volume at redshift z ---
  FUNCTION ComputeComovingVolume( z ) RESULT( V )

    REAL(DP) , INTENT(in) :: z
    INTEGER               :: k
    REAL(DP)              :: r , V , dz

    ! --- Integrate with Trapezoidal rule to get comoving distance ---
    r  = 0.0d0
    dz = z / Nq
    
    DO k = 1 , Nq - 1
       r = r + OneOverE( k * dz )
    END DO

    r = c / H0 * dz / 2.0d0 * ( OneOverE( 0.0d0 ) + 2.0d0 * r + OneOverE( z ) )

    V = 4.0d0 / 3.0d0 * PI * r**3

    RETURN
  END FUNCTION ComputeComovingVolume

  
  ! --- Create bin edges for lookback time ---
  FUNCTION CreateLookbackTimeBinEdges( tShort , nBins ) RESULT( BinEdges )

    REAL(DP) , INTENT(in) :: tShort
    INTEGER  , INTENT(in) :: nBins
    INTEGER               :: i
    REAL(DP)              :: BinEdges(nBins+1)

    DO i = 1 , nBins + 1
      BinEdges(i) = tShort + ( i - 1 ) * dt
    END DO

    RETURN
  END FUNCTION CreateLookbackTimeBinEdges


END PROGRAM ComputeMergerRate
