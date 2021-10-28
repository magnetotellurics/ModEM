module Solver_QMR
  !
  use Constants
  use ModelOperator
  use PreConditioner
  use cVector
  use Solver_PCG
  !
  type, extends(Solver_PCG_t):: Solver_QMR_t
     !
     !class(ModelOperator_t), pointer :: model_operator
     !class(PreConditioner_t), pointer :: preconditioner
   contains
     !procedure, public :: Create
     !procedure, public :: Allocate
     !procedure, public :: DeAllocate
     procedure, public :: setOperators
     procedure, public :: Solve
  end type Solver_QMR_t
  
contains
  
  subroutine setOperators(self, model_operator, preconditioner)
    class(Solver_QMR_t) :: self
    class(ModelOperator_t) , pointer :: model_operator
    class(Preconditioner_t), pointer :: preconditioner
    
    self%model_operator => model_operator
    self%preconditioner => preconditioner
  end subroutine setOperators

  !**
  ! This is the QMR solver, using operators
  ! (including pre-conditioner solvers),
  ! defined through pointers as object data
  ! Also uses real variable omega to define  frequency.
  !
  ! on input x is initial guess, b is rhs
  ! on output x is approximate solution
  ! diagnostic variables nIter, relErr are set
  !  
  ! Code is taken from subroutine QMR in solvers.f90
  !*
  subroutine Solve(self, x, b)
    class(Solver_QMR_t), intent(inout) :: self
    class(cVector_t)   , allocatable, intent(inout) :: x
    class(cVector_t)   , intent(in)    :: b
    ! Local variables
    class(cVector_t), allocatable :: AX, R, VT
    class(cVector_t), allocatable :: Y,Z,WT,V,W,YT,ZT,P,Q,PT,D,S
    logical :: adjoint, ilu_adjt
    complex(kind = prec) :: ETA,PDE,EPSIL,RDE,BETA,DELTA,RHO
    complex(kind = prec) :: PSI,RHO1,GAMM,GAMM1,THET,THET1,TM2
    complex(kind = prec) :: bnorm,rnorm
    complex(kind = prec) :: rhoInv,psiInv
    integer :: iter
    real(kind = 8) :: omega
    
    ! local copy of omega to simplify code slightly
    omega = self%omega
    
    ! Allocate work TVector objects -- questions as in PCG
    allocate(AX, source = x)
    allocate(R, source = x)
    allocate(VT, source = x)
    allocate(Y, source = x)
    allocate(Z, source = x)
    allocate(WT, source = x)
    allocate(V, source = x)
    allocate(W, source = x)
    allocate(YT, source = x)
    allocate(ZT, source = x)
    allocate(P, source = x)
    allocate(PT, source = x)
    allocate(D, source = x)
    allocate(S, source = x)
    
    ! NOTE: this iterative solver is QMR without look-ahead
    ! patterned after the scheme given on page 24 of Barrett et al.
    ! "Templates for the solution of linear systems of equations:
    ! Building blocks for iterative methods"
    ! Note that there are a couple of small differences, due to
    ! the fact that our system is complex (agrees with
    ! matlab6 version of qmr)
    
    self%failed = .false.
    adjoint = .false.
    ! R is Ax
    R = self%model_operator%Amult(R, adjoint, omega)
    ! b - Ax, for inital guess x, that has been input to the routine
    !
	! IMPLEMENT linComb ON VECTOR
	!
	!R = b%linComb(C_ONE, R, C_MinusOne)
    
    ! Norm of rhs, residual
    bnorm = CDSQRT(b%dotProd(b))
    rnorm = CDSQRT(R%dotProd(R))

    ! this usually means an inadequate model, in which case Maxwell's fails
    if (isnan(abs(bnorm))) then
       stop 'Error: b in QMR contains NaNs; exiting...'
    end if
    
    !  iter is iteration counter
    iter = 1
    self%relErr(iter) = real(rnorm/bnorm)

    VT = R
    ilu_adjt = .false.
    Y = self%preconditioner%LTsolve(VT, ilu_adjt)
    RHO = CDSQRT(Y%dotProd(Y))
    
    WT = R
    ilu_adjt = .true.
    Z = self%preconditioner%UTsolve(WT,ilu_adjt)
    PSI  = CDSQRT(Z%dotProd(Z))
    GAMM = C_ONE
    ETA  = C_MinusONE
    
    ! the do loop goes on while the relative error is greater than the tolerance
    ! and the iterations are less than maxIt
    do while ((self%relErr(iter).gt.self%tolerance).and.&
         (iter.lt.self%max_iter))
       if ((RHO.eq.C_ZERO).or.(PSI.eq.C_ZERO)) then
          self%failed = .true.
          write(0,*) 'QMR FAILED TO CONVERGE : RHO'
          stop 'QMR FAILED TO CONVERGE : PSI'
       end if

       rhoInv = (1/RHO)*cmplx(1.0, 0.0, 8)
       psiInv = (1/PSI)*cmplx(1.0, 0.0, 8)
       V = VT%mult(rhoInv)
       W = WT%mult(psiInv)
       Y = Y%mult(rhoInv)
       Z = Z%mult(psiInv)
       
       DELTA = Z%dotProd(Y)
       if (DELTA.eq.C_ZERO) then
          self%failed = .true.
          stop 'QMR FAILS TO CONVERGE : DELTA'
       end if

       ilu_adjt = .false.
       YT = self%preconditioner%UTsolve(Y,ilu_adjt)
       ilu_adjt = .true.
       ZT = self%preconditioner%LTsolve(Z,ilu_adjt)
       
       if (iter.eq.1) then
          P = YT
          Q = ZT
       else
          ! these calculations are only done when iter greater than 1
          PDE = -PSI*DELTA/EPSIL
          RDE = -RHO*CONJG(DELTA/EPSIL)
          !P = YT%linComb(C_ONE,P,PDE)
          !Q = ZT%linComb(C_ONE,Q,RDE)
       end if
       
       adjoint = .false.
       call PT%Zeros()
       PT = self%model_operator%Amult(P,adjoint,omega)
       EPSIL = Q%dotProd(PT)
       if (EPSIL.eq.C_ZERO) then
          self%failed = .true.
          stop 'QMR FAILED TO CONVERGE : EPSIL'
       end if
       
       BETA = EPSIL/DELTA
       if (BETA.eq.C_ZERO) then
          self%failed = .true.
          stop 'QMR FAILED TO CONVERGE : BETA'
       end if
       !VT = PT%linComb(C_ONE,V,-BETA)

       RHO1 = RHO
       ilu_adjt = .false.
       Y = self%preconditioner%LTsolve(VT, ilu_adjt)
       RHO = CDSQRT(Y%dotProd(Y))

       adjoint = .true.
       call WT%Zeros()
       WT = self%model_operator%Amult(Q,adjoint,omega)
       !  perhaps should use linComb here ...
       !WT = W%scMultAdd(-conjg(BETA), WT)
       
       ilu_adjt = .true.
       Z = self%preconditioner%UTsolve(WT,ilu_adjt)
       PSI = CDSQRT(Z%dotProd(Z))
       
       if (iter.gt.1) then
          THET1 = THET
       end if
       THET = RHO/(GAMM*CDABS(BETA))
       GAMM1 = GAMM
       GAMM = C_ONE/CDSQRT(C_ONE + THET*THET)
       if (GAMM.eq.C_ZERO) then
          self%failed = .true.
          stop 'QMR FAILS TO CONVERGE : GAMM'
       end if

       ETA = -ETA*RHO1*GAMM*GAMM/(BETA*GAMM1*GAMM1)
       if (iter.eq.1) then
          !D = P%scMult(ETA)
          !S = PT%scMult(ETA)
       else
          TM2 = THET1*THET1*GAMM*GAMM
          !D = P%linComb(ETA,D,TM2)
          !S = PT%linComb(ETA,S,TM2)
       end if
       
       !x = D%scMultadd(C_ONE,x)
       !R = S%scMultadd(C_MinusONE,R)
       ! A new AX
       rnorm = CDSQRT(R%dotProd(R))
       iter = iter + 1
       
       ! Keeping track of errors
       ! QMR book-keeping between divergence correction calls
       self%relErr(iter) = real(rnorm/bnorm)       
    end do
  end subroutine solve
  
end module Solver_QMR
