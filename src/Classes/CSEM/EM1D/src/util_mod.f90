!****************************************************************
!
!> EM modeling module containing some utility functions eg for
!>   error handling and sorting
!
!> Rita Streich 2011
!**************************************************************
module util_mod

  !mpi module
#ifdef USE_MPI
#ifndef WINDOWS
  use mpi
#endif
#endif

  use machine_dep_mod
  implicit none
  
#ifdef USE_MPI
#ifdef WINDOWS
  include 'mpif.h'
#endif
#endif

public::    alloc_error, &        !subroutine in errors.F90
            readwrite_error, &    !>        "
            open_error, &         !>        "
            close_error, &        !>        "
            exist_error, &        !>        "
            invalid_error         !>        "
            
#ifdef USE_MPI
public::    error_mpi             !>        "
#endif

  integer(kind=int32) :: pid = 0   !my process ID in MPI_COMM_WORLD (initialize in case MPI is not used)
  integer(kind=int32) :: nproc = 1 !number of processes in MPI_COMM_WORLD

contains

#include "errors.f90"
include 'sort_dbl.f90'
include 'indexx.f90'
#include "get_time.f90"
include 'realloc.f90'
#include "tictoc.f90"

endmodule util_mod
