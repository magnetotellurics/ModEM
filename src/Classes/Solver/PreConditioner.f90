!**
! Properties declared in abstract class will be: ModOp, omega (in all
! cases I can now think of, preconditioner depends on the full operator,
! which depends on omega);  setPreconditioner(self,omega) should be a procedure
! in the abstract class, since we will assume that this is always defined.
!*
module PreConditioner
  use cVector
  
  private

  public :: PreConditioner_t
  
  type, abstract :: PreConditioner_t
   contains
     procedure(iface_LTsolve), deferred, public :: LTsolve
     procedure(iface_UTsolve), deferred, public :: UTsolve
     procedure(iface_Minv)   , deferred, public :: Minv
  end type PreConditioner_t
  
  abstract interface

     function iface_LTsolve(self, inE, adjt) result(outE)
       import :: PreConditioner_t, cVector_t
       class(PreConditioner_t) :: self
       class(cVector_t), intent(in) :: inE
       logical         , intent(in) :: adjt
       class(cVector_t), allocatable :: outE
     end function iface_LTsolve

     function iface_UTsolve(self, inE, adjt) result(outE)
       import :: PreConditioner_t, cVector_t
       class(PreConditioner_t) :: self
       class(cVector_t), intent(in) :: inE
       logical         , intent(in) :: adjt
       class(cVector_t), allocatable :: outE
     end function iface_UTsolve

     subroutine iface_Minv(self, inPhi, outPhi)
       import :: PreConditioner_t, cVector_t
       class(PreConditioner_t), intent(in) :: self
       class(cVector_t)       , intent(in)    :: inPhi
       class(cVector_t)       , intent(inout) :: outPhi
     end subroutine iface_Minv

  end interface
end module PreConditioner
