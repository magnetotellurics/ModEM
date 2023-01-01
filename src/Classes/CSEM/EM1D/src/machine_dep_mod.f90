!---------------------------------------------------------------------------
!> FD EM module machine_dep_mod
!
!>   Contains constants that should make the software portable
!> 
!> Rita Streich 2009
!---------------------------------------------------------------------------

module machine_dep_mod

#ifdef USE_MPI
#ifndef WINDOWS
  use mpi
#endif
#endif
  implicit none
#ifdef USE_MPI
#ifdef WINDOWS
  include 'mpif.h'
#endif
#endif

  private
            
!> Compiler dependent unit numbers for stderr, stdout and stdin.
  integer, parameter, public :: stderr = 0
  integer, parameter, public :: stdout = 6
  integer, parameter, public :: stdin  = 5
  integer, parameter, public :: stdmes = 99


!> Should select the correct kinds for 32 and 64-bit reals.
  integer, parameter, public :: real32 = selected_real_kind(6)   !>  32-bit real
  integer, parameter, public :: real64 = selected_real_kind(15)  !>  64-bit real
  integer, parameter, public :: real128 = selected_real_kind(31)  !>  128-bit real

!> Should select the correct kinds for 8, 16, 32 and 64-bit ints.
  integer, parameter, public :: int8 = selected_int_kind(2)    !>  8-bit integer
  integer, parameter, public :: int16 = selected_int_kind(4)   !> 16-bit integer
  integer, parameter, public :: int32 = selected_int_kind(9)   !> 32-bit integer
  integer, parameter, public :: int64 = selected_int_kind(18)  !> 64-bit integer

!> Select default character kind of processor
  integer, parameter, public :: real_kind = kind(1.0)
  integer, parameter, public :: int_kind  = kind(1)    
  integer, parameter, public :: char_kind = kind('a')

#ifndef USE_MPI
  integer, parameter, public :: MPI_OFFSET_KIND = int64
#endif

endmodule machine_dep_mod
