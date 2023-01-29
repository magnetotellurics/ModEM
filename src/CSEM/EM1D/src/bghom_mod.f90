!****************************************************************
!
!>  FD EM module containing stuff for homogeneous background computations
!
!>  Rita Streich 2009-2011
!**************************************************************
module bghom_mod

  use machine_dep_mod
  use constants_mod
  use io_common_mod
!!$  use userpars_mod
!!$  use em_mpi_mod
!!$  use em_io_mod

  implicit none

  private

  !public subroutines
  public :: background_hom_unified

contains

  include 'background_hom.f90'
  include 'greenshom.f90'
  include 'radius.f90'

endmodule bghom_mod

