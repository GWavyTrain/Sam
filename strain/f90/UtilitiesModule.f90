MODULE UtilitiesModule

  IMPLICIT NONE
  PUBLIC

  INTEGER,  PARAMETER :: DP = KIND( 1.d0 )
  INTEGER,  PARAMETER :: nRedshifts = 10000
  REAL(DP), PARAMETER :: PI = ACOS( -1.0d0 )
  REAL(DP), PARAMETER :: OMEGA_M = 0.2726d0, OMEGA_L = 0.7274d0, &
                         H0 = 70.4d0 / 3.086d19
  REAL(DP), PARAMETER :: Year = 86400.0d0 * 365.0d0, &
                         SecondsPerGyr = 1.0d9 * Year
  REAL(DP), PARAMETER :: c = 2.99792458d10, G = 6.673d-8
  REAL(DP), PARAMETER :: Msun = 1.98892d33, tau = 4.0d0 * Year
  REAL(DP)            :: hc, z, tLb, r, M1, M2, f_ISCO
  INTEGER             :: i

  PUBLIC :: RombergIntegration
  PUBLIC :: Integrand_r, Integrand_tLb
  PUBLIC :: ComputeCharacteristicStrain
  PUBLIC :: InterpolateLookbackTimeToGetRedshift

  INTERFACE
    PURE REAL(KIND(1.d0)) FUNCTION Integrand( z ) RESULT( y )
      INTEGER, PARAMETER :: DP = KIND(1.d0)
      REAL(DP), INTENT(IN) :: z
    END FUNCTION
  END INTERFACE


CONTAINS


  SUBROUTINE TrapezoidalRule( Func, a, b, s, N )

    ! --- Taken from Press et al., Numerical Recipes in Fortran 77,
    !        pp. 131 (trapzd) ---

    PROCEDURE(Integrand)    :: Func
    REAL(DP), INTENT(in)    :: a, b
    REAL(DP), INTENT(inout) :: s
    INTEGER,  INTENT(in)    :: N

    REAL(DP) :: dx, sIter, nIterDP, x
    INTEGER  :: i, nIter

    IF( N .EQ. 1 )THEN
      s = 0.5d0 * ( b - a ) * ( Func(a) + Func(b) )
    ELSE
      nIter = 2**( N - 2 )
      nIterDP = nIter
      dx = ( b - a ) / nIterDP
      x = a + 0.5d0 * dx
      sIter = 0.0d0
      DO i = 1, nIter
        sIter = sIter + Func(x)
        x = x + dx
      END DO
      s = 0.5d0 * ( s + ( b - a ) * sIter / nIterDP )
    END IF

  END SUBROUTINE TrapezoidalRule


  SUBROUTINE PolynomialInterpolationAndExtrapolation( xa, ya, N, x, y, dy )

    ! --- Taken from Press et al., Numerical Recipes in Fortran 77,
    !        pp. 103-104 (polint) ---

    INTEGER,  INTENT(in)    :: N
    REAL(DP), INTENT(in)    :: xa(N), ya(N), x
    REAL(DP), INTENT(inout) :: y, dy

    INTEGER, PARAMETER :: Nmax = 10
    INTEGER            :: i, m, ind
    REAL(DP)           :: dx, dx_i, dxa, Coeff, C(Nmax), D(Nmax)

    ! --- Find the index, ind, of the closest array element ---
    ind = 1
    dx  = ABS( x - xa(1) )
    DO i = 1, N
      dx_i = ABS( x - xa(i) )
      IF( dx_i .LT. dx )THEN
        ind = i
        dx  = dx_i
      END IF
      C(i) = ya(i)
      D(i) = ya(i)
    END DO

    ! --- Initial approximation to y, i.e. value of ya at closest index, ind ---
    y   = ya(ind)
    ind = ind - 1

    ! --- Loop over columns, m (starting with second column) ---
    DO m = 1, N-1

      ! --- Loop over rows i in column m ---
      DO i = 1, N-m

         dxa = xa(i) - xa(i+m)

        ! This only happens if two input xa's are identical
        !   (within roundoff error)
        IF( dxa .EQ. 0.0d0 ) &
          STOP 'Failure in PolynomialInterpolationAndExtrapolation'

        ! --- Coefficient in Eq. (3.1.5) ---
        Coeff = ( C(i+1) - D(i) ) / dxa

        ! --- Update C's and D's, Eq. (3.1.5) ---
        D(i) = ( xa(i+m) - x ) * Coeff
        C(i) = ( xa(i  ) - x ) * Coeff

      END DO

      ! --- Decide which way to work through the tableau
      !     (don't yet understand how this works) ---
      IF( 2*ind .LT. N-m )THEN
        dy  = C(ind+1)
      ELSE
        dy  = D(ind)
        ind = ind - 1
      END IF
      y = y + dy

    END DO


  END SUBROUTINE PolynomialInterpolationAndExtrapolation


  SUBROUTINE RombergIntegration( Func, a, b, ss )

    ! --- Taken from Press et al., Numerical Recipes in Fortran 77,
    !        pp. 134 (qromb) ---

    PROCEDURE(Integrand)    :: Func
    REAL(DP), INTENT(in)    :: a, b
    REAL(DP), INTENT(inout) :: ss

    REAL(DP), PARAMETER :: EPS = 1.0d-11
    INTEGER,  PARAMETER :: jMax = 20, jMaxP = jMax+1, k = 5, kM = k-1
    REAL(DP)            :: dss, h(jMaxP), s(jMaxP)
    INTEGER             :: j

    ss   = 0.0d0
    h(1) = 1.0d0

    DO j = 1, jMax
      CALL TrapezoidalRule( Func, a, b, s(j), j )
      IF( j .GE. k )THEN
        CALL PolynomialInterpolationAndExtrapolation &
               ( h(j-kM), s(j-kM), k, 0.0d0, ss, dss )
        IF( ABS( dss ) .LE. EPS * ABS( ss ) ) RETURN
      END IF
      s(j+1) = s(j)
      h(j+1) = 0.25d0 * h(j)
    END DO

    STOP 'Too many steps in RombergIntegration'
    
  END SUBROUTINE RombergIntegration


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
    IF( tLb .LT. tLb_z_arr(nRedshifts) ) THEN

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


END MODULE UtilitiesModule
