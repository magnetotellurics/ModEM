!**********************************************************************
!> subroutine tic
!
!> Purpose:  starts measuring time with MPI subroutines
!
!> Alexander Grayver 2011
!> Rita Streich 2011
!**********************************************************************

subroutine tic(starttime)
  use, intrinsic :: iso_fortran_env
  implicit none

  !external variables
  real(kind=real64),intent(inout) :: starttime !> current MPI_Wtime
  !starttime = get_time() !get_time works with and without MPI
  
endsubroutine tic

!**********************************************************************
!> subroutine toc
!
!> Purpose:  compute the maximum MPI time interval between to events on all
!>           processes and print result on 0-rank process
!
!> Alexander Grayver 2011
!> Rita Streich 2011
!**********************************************************************
subroutine toc(starttime,msg,comm)
 USE, INTRINSIC :: ISO_FORTRAN_ENV
  implicit none

  !external variables
  real(kind=real64),intent( in ) :: starttime !> current MPI_Wtime
  character(len=*),intent( in ) :: msg       !> message to print
  integer,intent( in ) :: comm      !> communicator where we measure time
  
  !internal variables
  integer(kind=int32) :: ierr            !error index
  real(kind=real64) :: stoptime,tmax   !variables to store time
  integer(kind=int32) :: internalpid     !> my process number in comm

!>  stoptime = get_time()

#ifdef USE_MPI
  call MPI_Reduce(stoptime - starttime,tmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,comm,ierr)
  if(ierr.NE.MPI_SUCCESS) call error_mpi(pid,'toc','MPI_Reduce',ierr)
  
  call MPI_Comm_rank(comm,internalpid,ierr)
  if(internalpid.EQ.0) write(*,'(a,f10.2)') msg//' ', tmax
#else
  tmax = stoptime - starttime
  write(*,'(a,f10.2)') msg//' ', tmax
#endif

endsubroutine toc
