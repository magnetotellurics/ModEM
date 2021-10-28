module Solver
  !
  use Constants
  use cVector
  !
  type, abstract :: Solver_t
    !
	integer                        :: max_iter
	real(kind = prec)              :: tolerance
	real(kind = prec), allocatable :: relErr(:) ! relative error at each iteration
	!
	logical :: failed
	logical :: converged
	!
   contains
     procedure(iface_Solve), deferred, public :: Solve
     
  end type Solver_t

  abstract interface
     subroutine iface_Solve(self, x, b)
       import :: Solver_t, cVector_t
       class(Solver_t) , intent(inout) :: self
       class(cVector_t), intent(inout) :: x
       class(cVector_t), intent(in)    :: b
     end subroutine iface_Solve
  end interface
  
end module Solver
