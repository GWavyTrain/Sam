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


END MODULE UtilitiesModule
