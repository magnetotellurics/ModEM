module Solver_PCG
  use Constants
  use cVector
  use ModelOperator
  use PreConditioner
  use Solver
  
  implicit none
  
  private

  public :: Solver_PCG_t
  
  type, extends(Solver_t) :: Solver_PCG_t
     real(kind = prec) :: tol
     integer           :: maxIter 
     integer           :: nIter 
     real(kind = prec), allocatable :: relErr(:) ! relative error at each iteration
     logical           :: failed = .false.
     
     class(ModelOperator_t) , pointer :: ModOp
     class(PreConditioner_t), pointer :: preCond
   contains
     procedure, public :: Solve
  end type Solver_PCG_t

  interface Solver_PCG_t
     module procedure Solver_PCG_ctor
  end interface Solver_PCG_t
contains
  
  !**
  !*
  function Solver_PCG_ctor(ModOp, PreCond) result(obj)
    class(ModelOperator_t) , target :: ModOp
    class(PreConditioner_t), target :: preCond
    type(Solver_PCG_t) :: obj
    
    obj%modOp   => ModOp
    obj%preCond => preCond
  end function Solver_PCG_ctor
  
  !**
  ! This is the PCG solver, using operators
  ! (including pre-conditioner solvers),
  ! defined through pointers as object data.
  !
  ! on input x is initial guess, b is rhs
  ! on output x is approximate solution
  ! diagnostic variables nIter, relErr are set.
  !  
  !*
  subroutine Solve(self, x, b) 
    class(Solver_PCG_t), intent(inout) :: self
    class(cVector_t)   , intent(inout) :: x
    class(cVector_t)   , intent(in)    :: b    
    ! Local variables
    class(cVector_t), allocatable:: r, s, p, q
    complex(kind = prec) :: beta, alpha, delta, deltaOld
    complex(kind = prec) :: bnorm, rnorm
    integer :: i

    allocate(r, source = x)
    call r%Zeros()
    
    allocate(s, source = x)
    call s%Zeros()

    allocate(p, source = x)
    call p%Zeros()

    allocate(q, source = x)
    call q%Zeros()
    
    r = self%ModOp%Amult(x)
        
    r = b - r
    bnorm = b.dot.b
    rnorm = r.dot.r
    i = 1
    self%relErr(i) = real(rnorm/bnorm)
    
    do while ((self%relErr(i).gt.self%tol).and.(i.lt.self%maxIter))
       call self%preCond%Minv(r, s)
       delta = r.dot.s
       if (i.eq.1) then
          p = s
       else
          beta = delta/deltaOld
          p = s + beta*p
       end if
       q = self%Modop%Amult(p)
       alpha = delta/(p.dot.q)
       x = x + alpha*p
       r = q - alpha*r
       deltaOld = delta
       i = i + 1
       rnorm = r%dotProd(r)
       self%relErr(i) = real(rnorm/bnorm)
       
    end do
    
    self%niter = i
    
  end subroutine Solve
  
end module Solver_PCG
