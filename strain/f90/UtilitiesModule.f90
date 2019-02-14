MODULE UtilitiesModule

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: DP = KIND( 1.d0 )

  PUBLIC :: RombergIntegration

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

    REAL(DP), PARAMETER :: EPS = 1.0d-8
    INTEGER,  PARAMETER :: jMax = 20, jMaxP = jMax+1, k = 5, kM = k-1
    INTEGER             :: j
    REAL(DP)            :: dss, h(jMaxP), s(jMaxP)

    ss = 0.0d0
    h(1) = 1.0d0
    DO j = 1, jMax
      CALL TrapezoidalRule( Func, a, b, s(j), j )
      IF( j .GE. k )THEN
        CALL PolynomialInterpolationAndExtrapolation &
               ( h(j-kM), s(j-kM), k, 0.0d0, ss, dss )
        IF( ABS( dss ) .LE. EPS * ABS(ss) ) RETURN
      END IF
      s(j+1) = s(j)
      h(j+1) = 0.25d0 * h(j)
    END DO
    STOP 'Too many steps in RombergIntegration'
    
  END SUBROUTINE RombergIntegration


END MODULE UtilitiesModule
