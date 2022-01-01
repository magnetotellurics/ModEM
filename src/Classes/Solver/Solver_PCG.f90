module Solver_PCG
   !
   use Solver
   use cScalar
   use ModelOperator_MF
   use PreConditioner_MF_DC
   !
   !   solver object for PCG -- used only for Divergence Correction
   type, extends( Solver_t ), public :: Solver_PCG_t
      !
      ! PROPERTIES HERE
      !
      contains
      !
      final :: Solver_PCG_dtor
      !
      procedure, public :: solve => solvePCG
      procedure, public :: SetDefaults => setDefaults_PCG
      !
   end type Solver_PCG_t
   !
   interface Solver_PCG_t
      module procedure Solver_PCG_ctor
   end interface Solver_PCG_t
   !
   contains
      !
      function Solver_PCG_ctor( model_operator ) result( self )
         implicit none
         !
         class( ModelOperator_t ), target, intent( in ) :: model_operator
         type( Solver_PCG_t ) :: self
         !
         !write(*,*) "Constructor Solver_PCG_t"
         !
         call self%init()
         !
         self%model_operator => model_operator
         !
         ! PreConditioners need to be instantiated within the selection case
         ! as they receive a specific ModelOperator
         select type( model_operator )
            class is( ModelOperator_MF_t )
               !
               ! PreConditioner DC
               self%preconditioner = PreConditioner_MF_DC_t( model_operator )
               !
            class default
                 write(*, *) "ERROR:Solver_PCG::Constructor:"
                 STOP "         Unknow model_operator type."
         end select
         !
         !  NOTE: need default solver parameters to use here -- but more generally]
         !    need to set these to what user requests in setup program!!!
         call self%setParameters( 10, TOL6 )

         call self%zeroDiagnostics()
         !
      end function Solver_PCG_ctor
      !
      ! Solver_PCG destructor
      subroutine Solver_PCG_dtor( self )
          implicit none
          !
          type( Solver_PCG_t ), intent( inout ) :: self
          !
          !write(*,*) "Destructor Solver_PCG_t"
          !
          call self%dealloc()
          !
      end subroutine Solver_PCG_dtor
      !
      subroutine SetDefaults_PCG(self)
         class(Solver_PCG_t),intent(inout) :: self
         !    sets default iteration control parameters for QMR solver
         !    local variables
         integer, parameter :: max_iter = 100
         real(kind=prec),parameter ::  tolerance = 1E-5

         call self%SetParameters(max_iter,tolerance)

      end subroutine SetDefaults_PCG

      !
      !************************************************   
      subroutine solvePCG( self, b, x )
         !   This is the PCG solver, using operators
         !    (including pre-conditioner solvers),
         !    defined through pointers as object data`
         !
         !  on input x is initial guess, b is rhs
         !  on output x is approximate solution
         !  diagnostic variables nIter, relErr are set
         !  
         ! Code is taken from subroutine PCG in solvers.f90

         !   import :: Solver_t   -- is this needed????
         class( Solver_PCG_t ), intent(inout)           :: self
         class( cScalar_t ), intent(in)                 :: b
         class( cScalar_t ), allocatable, intent(inout) :: x
         ! local variables
         !   these will have to be created in a way to match
         !    the specific type of the input cScalar_t ...
         !   Can we just declare these to be of the abstract type?
         class ( cScalar_t ), allocatable :: r, s, p, q
         complex( kind=prec )        :: beta, alpha, delta, deltaOld
         complex( kind=prec )        :: bnorm, rnorm
         integer                     :: i

         !  create local cScalar objects -- could we also use modOp%createCScalar?
         allocate( r, source = x )   ! cannot zero x, since it is first guess
         call r%zeros()
         allocate( s, source = r )
         allocate( p, source = r )
         allocate( q, source = r )
         
         ! just like
         !call self%model_operator%Amult( x, r )
         !   if we can make AMult generic, with versions that operate on cScalar/cVector
         !    this could be more generic ...   also could change the name of this operator
         !     to make this more obvious
         call self%model_operator%divCgrad( x, r )
         !
         !    r = b-r
         call r%linCombS(b,C_MinusOne,C_ONE)
         !
         bnorm = sqrt(real(b%dotProd(b)))
         rnorm = sqrt(real(r%dotProd(r)))
         i = 1
         self%relErr(1) = rnorm/bnorm
         !
         loop: do while ( (self%relErr(i).gt.self%tolerance ).and.(i.lt.self%max_iter))
            !
            write(*,*) 'iter = ', i, ' relErr = ',self%relErr(i)
            call self%preconditioner%LUsolve( r, s ) 
            delta = r%dotProd(s)
            if(i.eq.1) then
               beta = C_ZERO   
            else
               beta = delta/deltaOld
            endif
            !   p = s * C_ONE + p * beta
            call p%linCombS(s,beta,C_ONE)
            !call q%Zeros()
            call self%model_operator%divCgrad( p, q )
            !
            alpha = delta/p%dotProd(q)
            !   this returns x = x + alpha*p
            call p%scMultAddS(x,alpha)
            !   this returns r = r - alpha*q
            call q%scMultAddS(r,-alpha)
            deltaOld = delta
            i = i + 1
            rnorm = sqrt(real(r%dotProd(r)))
            self%relErr(i) = rnorm/bnorm
         enddo loop
         !
         self%n_iter = i
         !
      end subroutine solvePCG ! PCG

end module Solver_PCG
