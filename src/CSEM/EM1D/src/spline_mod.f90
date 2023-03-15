!****************************************************************
!
!  Module containing routines for cubic spline interpolation:
!    f90 versions of numerical recipes routines
!
!  Rita Streich 2009
!
!**************************************************************
module spline_mod

  use machine_dep_mod
  use util_mod

  implicit none

  private

  ! Public subroutines
  public :: spline, &
            splint

  real(kind=real64),parameter,public :: spl_endval = 1.e30_real64 !indicates that spline derivatives at end points = 0

contains

#include "splines.f90"

endmodule spline_mod

