MODULE UtilitiesModule

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: DP = KIND( 1.d0 )

  PUBLIC :: TrapezoidalRule

  INTERFACE
    PURE REAL(KIND(1.d0)) FUNCTION Integrand( z ) RESULT( y )
      INTEGER, PARAMETER :: DP = KIND(1.d0)
      REAL(DP), INTENT(IN) :: z
    END FUNCTION
  END INTERFACE


CONTAINS


  SUBROUTINE TrapezoidalRule( Func, a, b, s, N )

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


  SUBROUTINE polint( xa, ya, n, x, y, dy )

    INTEGER,  INTENT(in)    :: N
    REAL(DP), INTENT(in)    :: xa(N), ya(N), x
    REAL(DP), INTENT(inout) :: y, dy

    INTEGER, PARAMETER :: nMax = 10
    INTEGER            :: i, m, ns
    REAL(DP)           :: den, dif, dift, ho, hp, w, c(nMax), d(nMax)

    ns = 1
    dif = ABS( x - xa(1) )
    DO i = 1, N
      dift = ABS( x - xa(i) )
      IF( dift .LT. dif )THEN
        ns = i
        dif = dift
      END IF
      c(i) = ya(i)
      d(i) = ya(i)
    END DO
    y = ya(ns)
    ns = ns - 1
    DO m = 1, n - i
      DO i = 1, n - m
        ho = xa(i) - x
        hp = xa(i+m) - x
        w = c(i+1)-d(i)
        den = ho - hp
        IF( den .EQ. 0 ) PAUSE 'Failure in polint'
        den = w / den
        d(i) = hp * den
        c(i) = ho * den
      END DO
      IF( 2 * ns .LT. n-m )THEN
        dy = c(ns+1)
      ELSE
        dy = d(ns)
        ns = ns-1
      END IF
      y = y + dy
    END DO


  END SUBROUTINE polint


END MODULE UtilitiesModule
