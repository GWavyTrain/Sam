PROGRAM CheckMaxArraySize

  IMPLICIT NONE

  INTEGER, PARAMETER    :: DP = KIND( 1.d0 )
  REAL(DP), ALLOCATABLE :: arr(:)
  INTEGER               :: i

  DO i = 0, 10
    ALLOCATE( arr(INT8(10**i)) )
    WRITE(*,*) i, SHAPE( arr )
    arr = 1.0d0
    DEALLOCATE( arr )
  END DO

END PROGRAM CheckMaxArraySize
