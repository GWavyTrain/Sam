PROGRAM CheckMaxArraySize

  IMPLICIT NONE

  INTEGER, PARAMETER    :: DP = KIND( 1.d0 )
  REAL(DP), ALLOCATABLE :: arr(:,:)
  INTEGER               :: i

  DO i = 0, 12
    ALLOCATE( arr(INT8(10**i),2) )
    WRITE(*,*) i, SHAPE( arr )
    arr = 1.0d0
    DEALLOCATE( arr )
  END DO

END PROGRAM CheckMaxArraySize
