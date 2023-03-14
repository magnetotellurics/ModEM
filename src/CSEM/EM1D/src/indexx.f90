!***********************************************************
!> 
!> Code converted using TO_F90 by Alan Miller
!> Date: 2007-06-28  Time: 15:21:53
!
!> numerical recipes routine indexx
!> from Numerical recipes in Fortran, 2nd ed. Press et. al
!
!> indexes an array arr(1:n) i.e. outputs an array indx(1:n) such that
!> arr(indx(j)) is in ascending order for j=1,2,...N
!> input quantities n and arr are not changed
!
!***********************************************************
SUBROUTINE indexx(n,arr,indx)

IMPLICIT NONE

INTEGER(kind=int32), INTENT(IN) :: n
REAL(kind=real64),dimension(n),INTENT(IN) :: arr
INTEGER(kind=int32),dimension(n),INTENT(OUT) :: indx


INTEGER(kind=int32), PARAMETER :: m=7
INTEGER(kind=int32), PARAMETER :: nstack=50
INTEGER(kind=int32) :: i,indxt,ir,itemp
INTEGER(kind=int32) :: j,jstack,k,l,istack(nstack)
REAL(kind=real64) :: a


DO  j=1,n
  indx(j) = j
enddo

jstack = 0
l = 1
ir = n

1     IF(ir-l < m) THEN
  DO  j=l+1,ir
    indxt=indx(j)
    a=arr(indxt)
    DO  i=j-1,1,-1
      IF(arr(indx(i)) <= a) GO TO 2
      indx(i+1) = indx(i)
    enddo
    i=0
    2         indx(i+1) = indxt
  enddo
  if(jstack == 0) RETURN
  ir = istack(jstack)
  l = istack(jstack-1)
  jstack = jstack - 2
ELSE
  k = (l+ir)/2
  itemp = indx(k)
  indx(k) = indx(l+1)
  indx(l+1) = itemp
  if(arr(indx(l+1)) > arr(indx(ir))) THEN
    itemp = indx(l+1)
    indx(l+1) = indx(ir)
    indx(ir) = itemp
  endif
  if(arr(indx(l)) > arr(indx(ir))) THEN
    itemp = indx(l)
    indx(l) = indx(ir)
    indx(ir) = itemp
  endif
  if(arr(indx(l+1)) > arr(indx(l))) THEN
    itemp = indx(l+1)
    indx(l+1) = indx(l)
    indx(l) = itemp
  endif
  i = l + 1
  j = ir
  indxt = indx(l)
  a = arr(indxt)
  3       CONTINUE
  i = i + 1
  if(arr(indx(i)) < a) GO TO 3
  4       CONTINUE
  j = j - 1
  if(arr(indx(j)) > a) GO TO 4
  if(j < i) GO TO 5
  itemp = indx(i)
  indx(i) = indx(j)
  indx(j) = itemp
  GO TO 3
  5       indx(l) = indx(j)
  indx(j) = indxt
  jstack = jstack + 2
  if(jstack > nstack) then
    write( *, * ) 'Error in indexx: nstack too small'
    stop
  endif
  if(ir-i+1 >= j-l) THEN
    istack(jstack) = ir
    istack(jstack-1) = i
    ir = j-1
  ELSE
    istack(jstack) = j - 1
    istack(jstack-1) = l
    l = i
  endif
endif
GO TO 1
END SUBROUTINE indexx

