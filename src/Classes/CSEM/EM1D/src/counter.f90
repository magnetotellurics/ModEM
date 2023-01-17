!---------------------------------------------------------------
!> EM1D routine counter_service
!
!> a simple MPI-1 based routine for a counter
!> managed by one process
!> must be called frequently
!> This is cleaner than the MPI-2 version with "one-sided" communication,
!>  since one-sided comm. is actually delayed until process 0 calls some MPI routine!
!
!> Rita Streich 2009, adapted from Using MPI (Gropp et al.)
!---------------------------------------------------------------
subroutine counter_service(comm,val,valmax,ndone)

  implicit none

  !external variables
  integer(kind=int32),intent( in ) :: comm   !MPI communicator
  integer(kind=int32) :: val    !the variable to be incremented
  integer(kind=int32),intent( in ) :: valmax !max value that val should take
  integer(kind=int32) :: ndone  !how many processes have received a stop signal

  !internal variables
  integer(kind=int32) :: ierr                !error index
  integer(kind=int32) :: requester           !ID of process asking for a value
  logical :: flag                !indicates if any requests have been made
  integer(kind=int32),dimension(MPI_STATUS_SIZE) :: stat  !MPI status
  integer(kind=int32) :: dum                 !dummy receive buffer

!> workaround for openmpi bug
!!$#ifdef OLDOMP
!!$  integer(kind=int32) :: iproc               !process counter
!!$#endif

  !remember max value, if input max value is negative, this is a flag for using the the max value stored in counter module
  if(valmax.gt.0) count_max = valmax

  do
    !MPI_Iprobe Works fine in mvapich0.9.9.
    call MPI_Iprobe(MPI_ANY_SOURCE,countertag,comm,flag,stat,ierr)
    if(ierr.NE.MPI_SUCCESS) call error_mpi(pid,'counter_service','MPI_Iprobe',ierr)
    if(.NOT.flag) then
      exit
    else
      requester = stat(MPI_SOURCE)
      call MPI_Recv(dum,0,MPI_INTEGER,requester,countertag,comm,stat,ierr)
      if(ierr.NE.MPI_SUCCESS) call error_mpi(pid,'counter_service','MPI_Recv',ierr)
      val = val + 1
      call MPI_RSend(val,1,MPI_INTEGER,requester,countertag,comm,ierr)
      if(ierr.NE.MPI_SUCCESS) call error_mpi(pid,'counter_service','MPI_Send',ierr)

      if(val .gt. count_max) then
        ndone = ndone + 1
      endif
    endif

  enddo

!!$#endif

endsubroutine counter_service



!---------------------------------------------------------------
!> EM1D routine counter_finish
!
!> send stop signal to all processes after all frequencies have been done
!
!> Rita Streich 2009
!---------------------------------------------------------------
subroutine counter_finish(comm,val,ndone,ngrp)

  implicit none
#ifndef X64
#ifdef WINDOWS
!DEC$ ATTRIBUTES ALIAS:'_SLEEP' :: sleep
#endif
#endif

  !external variables
  integer(kind=int32),intent( in ) :: comm    !MPI communicator
  integer(kind=int32) :: val     !the counter
  integer(kind=int32) :: ndone   !how many processes have received the stop signa;
  integer(kind=int32),intent( in ) :: ngrp    !nr of process groups for parallel frequencies

  !internal variables
  integer(kind=int32) :: ierr                !error index
  integer(kind=int32) :: requester           !ID of process asking for a value
  logical :: flag                !indicates if any requests have been made
  integer(kind=int32),dimension(MPI_STATUS_SIZE) :: stat  !MPI status
  integer(kind=int32) :: dum                 !dummy receive buffer


  do
    call MPI_Iprobe(MPI_ANY_SOURCE,countertag,comm,flag,stat,ierr)
    if(ierr.NE.MPI_SUCCESS) call error_mpi(pid,'counter_service','MPI_Iprobe',ierr)

    if(.NOT.flag) then
      call sleep(1) !sleep 1 second
    else
      requester = stat(MPI_SOURCE)
      call MPI_Recv(dum,0,MPI_INTEGER,requester,countertag,comm,stat,ierr)
      if(ierr.NE.MPI_SUCCESS) call error_mpi(pid,'counter_service','MPI_Recv',ierr)
      val = val + 1
      call MPI_Send(val,1,MPI_INTEGER,requester,countertag,comm,ierr)
      if(ierr.NE.MPI_SUCCESS) call error_mpi(pid,'counter_service','MPI_Send',ierr)

      ndone = ndone + 1
    endif

    if(ndone.EQ.(ngrp-1)) exit

  enddo

endsubroutine counter_finish

