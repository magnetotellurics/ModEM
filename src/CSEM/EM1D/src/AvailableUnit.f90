!----------------------------------------------------------------
! FD EM function AvailableUnit
!
! Find an unused unit number, for secured Fortran IO-handling.
!
! Rita Streich 2009
!----------------------------------------------------------------

integer(kind=int32) function AvailableUnit()
 use, intrinsic :: iso_fortran_env
  implicit none

  integer(kind=int32) :: i
  logical             :: isOpen
#ifdef USE_MPI
  integer(kind=int32) :: ierr
#endif

  do i=20,98765
     inquire(Unit=i,Opened=isOpen)
     if (.not.isOpen) then
        AvailableUnit = i
        return
     end if
  end do

  !we get here if no unit number was found
  write(*,'(a)')'ERROR in AvailableUnit: No unused unit number available'
#ifdef USE_MPI
  call MPI_Abort (MPI_COMM_WORLD,-1,ierr)
#endif
  stop

end function AvailableUnit


