!****************************************************************
!
!  1D EM module containing stuff for homogeneous background computations
!
!  Rita Streich 2009
!
!**************************************************************
module hankel_mod

  use machine_dep_mod
  implicit none

  private

  !public subroutines
  public :: zhankl
            

integer(kind=int32),parameter,public  :: filtlen = 801  !Hankel filter length
real(kind=real64),parameter,public    :: logspace = 0.1_real64 !logarithmic argument spacing required by FHT
integer(kind=int32),parameter,public  :: NB = 1   !number of "lagged" Hankel integrals (same kernel with different radii)
integer(kind=int32),parameter,public  :: nrelmax = 3 !max. number of related Hankel integrals
integer(kind=int32),parameter,public  :: NREL = 1 !(old default) number of related Hankel integrals
integer(kind=int32),parameter,public  :: NTOL = 1
real(kind=real64),parameter,public    :: TOL = 0._real64
complex(kind=real64),dimension(filtlen,nrelmax),public :: ZWORK   ! work array for related and lagged transforms
integer(kind=int32),dimension(2,nrelmax),public        :: IJREL   ! exponents for Hankel transforms of related functions
integer(kind=int32),public            :: ikap     !make wavenumber counter visible outside hankel routine


contains

include 'hankel.f90'


endmodule hankel_mod

