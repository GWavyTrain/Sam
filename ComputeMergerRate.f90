PROGRAM ComputeMergerRate

  IMPLICIT NONE

  INTEGER  , PARAMETER   :: DP = KIND( 1.d0 )
  INTEGER  , PARAMETER   :: Nq = 10000 , nBins = 10
  INTEGER  , PARAMETER   :: itLB = 4 , iz = 5
  REAL(DP) , PARAMETER   :: PI = ACOS( -1.0d0 )
  REAL(DP) , PARAMETER   :: OMEGA_M = 0.3d0 , OMEGA_L = 0.7d0
  REAL(DP) , PARAMETER   :: H0 = 70.4d0 / 3.086d19
  REAL(DP) , PARAMETER   :: c = 2.99792458d10
  REAL(DP) , PARAMETER   :: year = 86400.0d0 * 365.25d0 , Gyr = 1.0d9 * year
  REAL(DP)               :: BinEdges(nBins+1) , Bins(nBins) , Counts(nBins)
  REAL(DP)               :: tShort , tLong , zMin , zMax
  REAL(DP) , ALLOCATABLE :: hc(:,:) , data(:,:)
  INTEGER                :: i , nMergers
  CHARACTER( len = 32 )  :: arg , FileName

  ! --- Get file name ---
  CALL GETARG( 1 , arg )
  READ( arg , * ) FileName

  ! --- Get number of mergers in hc.dat file ---
  OPEN( 100 , FILE = TRIM( FileName ) )
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
  OPEN( 100 , FILE = TRIM( FileName ) )
  READ( 100 , * ) ! --- Skip the header ---
  READ( 100 , * ) ! --- Skip the line with frequencies ---
  DO i = 1 , nMergers
     
    READ( 100 , * ) hc( i , : )
    data(i,1) = hc( i , itLB )                     ! [ Gyr ]
    data(i,2) = hc( i , iz   )                     ! [ dimensionless ]
    data(i,3) = ComputeComovingVolume( data(i,2) ) ! [ cm^3 ]
        
  END DO
  CLOSE( 100 )

  tShort = MINVAL( data(:,1) )
  zMin   = MINVAL( data(:,2) )
  tLong  = MAXVAL( data(:,1) )
  zMax   = MAXVAL( data(:,2) )

  Bins = CreateLookbackTimeBinEdges( tShort , tLong , nBins )

  WRITE(*,*) Bins
  
  DEALLOCATE( data )
  DEALLOCATE( hc )

CONTAINS

  ! --- Integrand E(z) in cosmological distance calculation ---
  PURE FUNCTION E( z )

    REAL(DP) , INTENT(in) :: z
    REAL(DP)              :: E

    E = 1.0d0 / SQRT( OMEGA_M * ( 1.0d0 + z )**3 + OMEGA_L )

    RETURN
  END FUNCTION E

  
  ! --- Compute comoving volume at redshift z ---
  FUNCTION ComputeComovingVolume( z ) RESULT( V )

    REAL(DP) , INTENT(in) :: z
    INTEGER               :: k
    REAL(DP)              :: r , V , dz

    ! --- Integrate with Trapezoidal rule to get comoving distance ---
    r  = 0.0d0
    dz = z / Nq
    
    DO k = 1 , Nq - 1
       r = r + E( k * dz )
    END DO

    r = c / H0 * dz / 2.0d0 * ( E( 0.0d0 ) + 2.0d0 * r + E( z ) )

    V = 4.0d0 / 3.0d0 * PI * r**3

    RETURN
  END FUNCTION ComputeComovingVolume

  
  ! --- Create binedges for lookback time ---
  FUNCTION CreateLookbackTimeBinEdges( tShort , tLong , nBins ) RESULT( BinEdges )

    REAL(DP) , INTENT(in) :: tShort , tLong
    INTEGER  , INTENT(in) :: nBins
    INTEGER               :: i
    REAL(DP)              :: BinEdges(nBins+1) , dt

    dt = ( tLong - tShort ) / nBins

    DO i = 1 , nBins + 1

       BinEdges(i) = tShort + ( i - 1 ) * dt

    END DO

    RETURN
  END FUNCTION CreateLookbackTimeBinEdges


END PROGRAM ComputeMergerRate
