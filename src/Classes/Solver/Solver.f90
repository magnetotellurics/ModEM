module Solver
   !
   use Constants
   use cVector
   use ModelOperator
   use PreConditioner
   !
   character(:), allocatable :: solver_type
   character (len = 3), parameter   :: QMR  = "QMR"
   character (len = 3), parameter   :: PCG  = "PCG"
   character (len = 4), parameter   :: BiCG = "BiCG"
   !
   ! SOLVER DEFAULTS
   integer  :: maxIter = 20
   real( kind=prec) :: tolerance = 1d-7
   !
   type, abstract :: Solver_t
      !
      integer                        :: max_iter, n_iter
      real( kind=prec )              :: omega, tolerance
      real( kind=prec ), allocatable :: relErr(:) ! relative error at each iteration
      !
      logical :: failed = .false., converged = .false.
      !
      class( ModelOperator_t ), pointer  :: model_operator
      !
      ! PreConditioner as a property of this
      class( PreConditioner_t ), allocatable :: preconditioner
      !
      contains
         !   deferred (abstract class) proedures
         procedure( interface_set_solver_defaults), deferred, public    :: SetDefaults
         !
         procedure, public :: init    => initializeSolver
         procedure, public :: dealloc => deallocateSolver
         !
         procedure, public :: setParameters
         procedure, public :: zeroDiagnostics
         !
   end type Solver_t
   !
   abstract interface
      !
      !    each type of solver will have own defaults, hard code in solver extension
      subroutine interface_set_solver_defaults(self)
         import :: Solver_t
         class( Solver_t ), intent(inout) :: self
      end subroutine interface_set_solver_defaults

   end interface
   !
contains
    !
    subroutine setParameters( self, max_iter, tolerance )
        !   actually, this set routine will be the same for all
        !    extensions -- should this be abstract
        ! import :: Solver_t  is this needed here????
        class( Solver_t ), intent( inout ) :: self
        integer, intent( in )              :: max_iter
        real( kind=prec ), intent( in )    :: tolerance
        integer :: status

        self%max_iter = max_iter
        self%tolerance = tolerance
        !   perhaps check if relErr is already allocated; if so deallocate
        if(allocated(self%relErr)) deallocate(self%relErr)
        allocate(self%relErr(max_iter),STAT=status)
        !  if we are not going to check "status" of allocate, why 
        !    return this?

    end subroutine setParameters
    !
    !********
    !
    subroutine zeroDiagnostics( self )
       !   zeros diagnostics for solver object
       class( Solver_t ), intent( inout ) :: self
       self%n_iter = 0
       self%relErr = R_ZERO
    end subroutine zeroDiagnostics 
    !
    !********
    !
    subroutine initializeSolver( self )
      implicit none
      !
      class( Solver_t ), intent( inout ) :: self
      !
      self%model_operator => null()
      !
      self%n_iter = 0
      self%omega = 0.0
      !
      self%failed = .false.
      self%converged = .false.
      !
   end subroutine initializeSolver
   !
   subroutine deallocateSolver( self )
      implicit none
      !
      class( Solver_t ), intent( inout ) :: self
      !
      !if( allocated( self%preconditioner ) ) deallocate( self%preconditioner )
      !
   end subroutine deallocateSolver
   !
end module Solver
