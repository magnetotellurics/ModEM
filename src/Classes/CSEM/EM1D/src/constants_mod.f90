!****************************************************************
!
! FD EM module containing the some constants used thoughout the program
!
! Rita Streich 2009-2011
!
!**************************************************************
module constants_mod

  use machine_dep_mod

  implicit none
  
  integer(kind=int32),parameter    :: namlen=300 !max length of file names

  !complex i
!  complex(kind=real32),parameter,public :: ci = (0._real32,1._real32)
  complex(kind=real64),parameter,public :: dci = (0._real64,1._real64)

  complex(kind=real64),parameter    :: zeroc=cmplx(0.0_real64,0.0_real64)  !complex zero

  !PI
!  real(kind=real32),parameter,public    :: pi = 3.141592654e0
!  real(kind=real64),parameter,public    :: dpi = 3.141592653589793115997963468544185162e0_real64
  !NEW: computed by matlab routine pi2str (Jonas Lundgren <splinefit@gmail.com> 2008),
  !  which can compute pi very precisely to LOTS of digits
  real(kind=real64),parameter,public    :: dpi = 3.1415926535897932384626433832795028841971e0_real64

  !2*PI
!  real(kind=real32),parameter,public    :: twopi = 6.283185307e0
  real(kind=real64),parameter,public    :: dtwopi = 2._real64*dpi

  !4*PI
!  real(kind=real32),parameter,public    :: fourpi = 12.56637061e0
  real(kind=real64),parameter,public    :: dfourpi = 4._real64*dpi

  !vacuum velocity of light
!  real(kind=real32),parameter,public    :: c0=299792458.0
  real(kind=real64),parameter,public    :: dc0=299792458.0d0  !this value is exact by definition

  !vacuum magnetic permeability
!  real(kind=real32),parameter,public    :: mu0 = 1.25663706e-06
  real(kind=real64),parameter,public    :: dmu0 = 4.e-7_real64*dpi

  !vacuum dielectric permittivity
!  real(kind=real32),parameter,public    :: eps0 = 8.85418782e-12
  real(kind=real64),parameter,public    :: deps0 = 1.0_real64/(dc0**2*dmu0)

end module constants_mod

