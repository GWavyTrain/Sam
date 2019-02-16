PROGRAM ModuleTester

  USE OMP_LIB
  USE UtilitiesModule, ONLY: &
    RombergIntegration

  IMPLICIT NONE

  INTEGER,  PARAMETER :: DP = KIND( 1.d0 ), nIterMax = 1000000000
  REAL(DP), PARAMETER :: OMEGA_M = 0.2726d0, OMEGA_L = 0.7274d0
  REAL(DP), PARAMETER :: ThreeFifths = 3.0d0 / 5.0d0, OneFifth = 1.0d0 / 5.0d0

  REAL(DP) :: a, b, y
  REAL(DP) :: M1, M2, Mc
  INTEGER  :: i
  REAL(DP) :: Time

  ! --- Using d0 vs. DP for math ---
!!$  Time = OMP_GET_WTIME()
!!$  DO i = 1, nIterMax
!!$    y = 10.0d0 + 7.45d0 + 1.0d2
!!$  END DO
!!$  Time = OMP_GET_WTIME() - Time
!!$  WRITE(*,'(A,ES13.6E3,A)') 'Time using d0 = ', Time, ' s'
!!$
!!$  Time = OMP_GET_WTIME()
!!$  DO i = 1, nIterMax
!!$    y = 10.0_DP + 7.45_DP + 100.0_DP
!!$  END DO
!!$  Time = OMP_GET_WTIME() - Time
!!$  WRITE(*,'(A,ES13.6E3,A)') 'Time using DP = ', Time, ' s'

  ! --- Different ways of organizing terms to compute chirp mass ---
!!$  M1 = 1.5d1
!!$  M2 = 3.0d1
!!$
!!$  Time = OMP_GET_WTIME()
!!$  DO i = 1, nIterMax
!!$    Mc = ( M1 * M2 )**( 3.0d0 / 5.0d0 ) / ( M1 + M2 )**( 1.0d0 / 5.0d0 )
!!$  END DO
!!$  Time = OMP_GET_WTIME() - Time
!!$  WRITE(*,'(A,ES13.6E3,A)') 'Time (1) = ', Time, ' s'
!!$
!!$  Time = OMP_GET_WTIME()
!!$  DO i = 1, nIterMax
!!$    Mc = ( M1 * M2 )**( ThreeFifths ) / ( M1 + M2 )**( OneFifth )
!!$  END DO
!!$  Time = OMP_GET_WTIME() - Time
!!$  WRITE(*,'(A,ES13.6E3,A)') 'Time (2) = ', Time, ' s'
!!$
!!$  Time = OMP_GET_WTIME()
!!$  DO i = 1, nIterMax
!!$    Mc = ( ( M1 * M2 )**3 / ( M1 + M2 ) )**( 1.0d0 / 5.0d0 )
!!$  END DO
!!$  Time = OMP_GET_WTIME() - Time
!!$  WRITE(*,'(A,ES13.6E3,A)') 'Time (3) = ', Time, ' s'
!!$
!!$  Time = OMP_GET_WTIME()
!!$  DO i = 1, nIterMax
!!$    Mc = ( ( M1 * M2 )**3 / ( M1 + M2 ) )**( OneFifth )
!!$  END DO
!!$  Time = OMP_GET_WTIME() - Time
!!$  WRITE(*,'(A,ES13.6E3,A)') 'Time (4) = ', Time, ' s'

  ! --- Timing SQRT() vs. ()**(1.0d0/2.0d0) ---
!!$  a = 5.0d-1
!!$
!!$  Time = OMP_GET_WTIME()
!!$  DO i = 1, nIterMax
!!$    y = SQRT( OMEGA_M * ( 1.0d0 + a )**3 + OMEGA_L )
!!$  END DO
!!$  Time = OMP_GET_WTIME() - Time
!!$  WRITE(*,'(A,ES13.6E3,A)') 'Time SQRT    = ', Time, ' s'
!!$
!!$  Time = OMP_GET_WTIME()
!!$  DO i = 1, nIterMax
!!$    y = ( OMEGA_M * ( 1.0d0 + a )**3 + OMEGA_L )**( 1.0d0 / 2.0d0 )
!!$  END DO
!!$  Time = OMP_GET_WTIME() - Time
!!$  WRITE(*,'(A,ES13.6E3,A)') 'Time **(1/2) = ', Time, ' s'

  ! --- Accuracy of Romberg integration subroutine ---
!!$  a = 0.0d0
!!$  b = 2.5d1
!!$  CALL RombergIntegration( Integrand, a, b, y )
!!$  WRITE(*,'(A,F12.10)')     'y    = ', y

  WRITE(*,'(A)') 'END OF LINE'

CONTAINS


  PURE REAL(DP) FUNCTION Integrand(x) RESULT(y)

    REAL(DP), INTENT(in) :: x

    y = 1.0d0 / SQRT( OMEGA_M * ( 1.0d0 + x )**3 + OMEGA_L )

    RETURN
  END FUNCTION Integrand


END PROGRAM ModuleTester
