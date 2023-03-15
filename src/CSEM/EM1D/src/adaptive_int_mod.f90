!****************************************************************
!
!  FD EM module for adaptive integration, used in computation of fields for 1D medium
!    in special cases for a few integrals
!    if the source and iReceiver are exactly at the same depth
!
!  Programmed: Rita Streich 2009-2011
!**************************************************************
module adaptive_int_mod

use machine_dep_mod
use constants_mod

implicit none

private
public :: BESAUT

integer(kind=int32),parameter :: nmaxsave = 5000
integer(kind=int32),parameter :: nmaxsum = 300
real(kind=real64),parameter,public    :: relerr=1.e-8_real64, relerr2=1.e-15_real64, abserr=0._real64
integer(kind=int32),parameter,public  :: GAUSLO = 1
integer(kind=int32),parameter,public  :: GAUSHI = 7

  real(kind=real64), parameter  :: relerrbes = 1.e-16_real64 !permitted error for explicit Bessel function evaluation
  real(kind=real64), parameter  :: rp2 = 2._real64 / dpi

contains

  include 'ChaveBesselTrafoSub.f90'
  include 'besselfunc.f90'

endmodule adaptive_int_mod
