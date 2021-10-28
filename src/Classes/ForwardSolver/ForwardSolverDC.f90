module ForwardSolverDC
   !
   use ForwardSolver
   use DivergenceCorrection
   use Solver_PCG
   use Solver_QMR
   use Preconditioner_MF
   !
   !   ModEM classic -- iterative solution with QMR or BiCG, divergence
   !   correction (DC) in usual way.   A slight variant/extension would be
   !   required for secondary field formulation.

   !   default values for solver control 
   !    since reasonable defaults depend on the details of forward solver
   !    details, set the defaults for parameters declared in solver object
   !    here also
   !    this will be default maxIter for QMR/BiCG solver
   integer, parameter    ::             IterPerDivCorDefQMR = 40
   integer, parameter    ::             IterPerDivCorDefBCG = 80
   ! maximum number of divergence correction calls allowed
   integer, parameter    ::             MaxDivCorDef = 20
   !   for DC object default MaxIterTotal = IterPerDivCor*MaxDivCor
   !     compute from other defaults, depending on solver
   ! maximum number of PCG iterations for divergence correction
   !    this will be default maxIter for Solver_PCG object
   integer, parameter    ::             MaxIterDivCorDef = 100
   
   ! misfit tolerance for curl-curl (QMR or BCG)
   real(kind=prec), parameter       ::   tolEMDef = 1E-7
   ! misfit tolerance for convergence of divergence correction solver
   real(kind=prec), parameter       ::   tolDivCorDef = 1E-5
   
   type, extends( ForwardSolver_t ), public :: ForwardSolverDC_t
      !
      type( DivergenceCorrection_t )     :: divergence_correction  !  pointer to divergence correction
      class( Solver_PCG_t ), allocatable :: solver_pcg   !   solver object for divergence correction
      type( PreConditioner_MF_t )        :: preconditioner                    !   not 100% clear this needs to be a property
      real( kind = prec )       :: omega
      real(kind = 8), allocatable, dimension(:,:) :: divJ
      real(kind = 8), allocatable, dimension(:)  :: EMrelErr
      real(kind = 8), allocatable, dimension(:,:) :: DivCorRelErr
      integer   :: nDivCor
      !   next two are not independent of MaxIterTotal
      integer :: MaxIterTotal
      integer :: MaxDivCor
      integer :: IterPerDivCor 
      !
   contains
      !
      procedure, public :: initDiagnosticArrays
      procedure, public :: getESolution => getESolutionForwardSolverDC
      !   set routines for iteration control parameters (specific to DC)
      !procedure, public :: set_MaxDivCor
      !procedure, public :: set_IterPerDivCor
      !   get routines for diagonstics
      !procedure, public :: get_nDivCor
      !procedure, public :: get_divJ
      !procedure, public :: get_DivCorRelErr
      !
   end type ForwardSolverDC_t
   !
   interface ForwardSolverDC_t
      module procedure ForwardSolverDC_ctor
   end interface ForwardSolverDC_t
   !
   contains
   !
   !
   function ForwardSolverDC_ctor() result( self )
      !
      class( ForwardSolverDC_t ), pointer :: self
      !
      !write(*,*) "Constructor ForwardSolverDC_t"
      !
      allocate( ForwardSolverDC_t :: self )
      !
      call self%init()
      !
   end function ForwardSolverDC_ctor
   !
      !   creator and destructor ... still need these
      subroutine setIterDefaults(self)
        type(ForwardSolverDC_t)   :: self
        !   this just sets iteration control parameters to default
        !    values -- should be good enough to get us started!
        !    Note that some of the parameters are set in solver object 

        self%MaxDivCor = MaxDivCorDef
        self%maxIterTotal = self%MaxDivCor*self%solver_pcg%maxIter
        select type( solver => self%solver_pcg )
           class is (Solver_QMR_t)
              self%IterPerDivCor = IterPerDivCorDefQMR
           !class is (Solver_BiCG_t)
              !self%IterPerDivCor = IterPerDivCorDefBCG
        end select
      
        !   these are parameters in the solver objects, for curl-curl
        !      and for DC
        !self%solver_pcg%setParameters(self%IterDivCor,tolEMDef)
        !self%solver_pcg%setParameters(MaxIterDivCorDef,tolDivCorDef)

      end subroutine setIterDefaults
      !***********************************************************
      subroutine createDiagnosticArrays(self)
        type(ForwardSolverDC_t)   :: self
        !   this allocates arrays for storage of solver diagnostics

        ! Forward object: EMrelErr, divJ, DivCorRelErr
        if( allocated( self%divJ ) ) deallocate( self%divJ )
        !
        allocate( self%divJ( 2, self%MaxDivCor ) )
        !
        if( allocated( self%EMrelErr) ) deallocate( self%EMrelErr )
        !
        allocate( self%EMrelErr( self%MaxIterTotal ) )
        !
        if( allocated( self%DivCorRelErr ) ) deallocate( self%DivCorRelErr )
        !
        allocate( self%DivCorRelErr( self%solver_pcg%maxIter, self%MaxDivCor ) )
        !
        !  Solver objects
        if( allocated( self%solver_pcg%relErr ) ) deallocate( self%solver_pcg%relErr )
        !
        allocate( self%solver_pcg%relErr( self%solver_pcg%maxIter ) )

     end subroutine createDiagnosticArrays
     !*****************************************************
     subroutine initDiagnosticArrays( self )
        !
		implicit none
		!
        class( ForwardSolverDC_t ), intent( inout ) :: self
        !   this zero's diagnostic arrays for storage of solver diagnostics
        self%nIterTotal = 0
        self%nDivCor = 0
        self%EMrelErr = R_ZERO
        self%divJ = R_ZERO
        self%DivCorRelErr = R_ZERO
        self%solver_pcg%failed = .false.
        self%solver_pcg%converged = .false.
        !
     end subroutine initDiagnosticArrays    
     !*****************************************************
     subroutine SetFrequency(self,omega)
        !   this is not specific to DC solver -- can we implement
        !     in abstract class?
        type(ForwardSolverDC_t),intent(inout)   :: self
        real(kind=8), intent(in)   :: omega
        if(abs(self%omega-omega) .lt. TOL8) then
          !   omega is close enough to input -- no need to reset
          !    freqeuncy dependent properties
          return
        end if
        ! otherwise need to update some things ...
        self%omega = omega
        self%solver_pcg%omega = omega
        call self%preconditioner%SetPreconditioner(omega)
     end subroutine SetFrequency
     !*****************************************************
     function getESolutionForwardSolverDC( self, source ) result( e_solution )
        implicit none
        !
        class( ForwardSolverDC_t ), intent( inout )  :: self
        class( Source_t ), allocatable, intent( in ) :: source
        !
        !
        class( cVector_t ), allocatable :: e_solution
        !   local variables
        class( cVector_t ), allocatable :: b    !  copy of RHS--do we really need?
        class( cVector_t ), allocatable :: temp
        class( cScalar_t ), allocatable :: phi0
        integer :: iter

        !   initialize diagnostics -- am assuming that setting of solver parameters
        !     is done in a set up step (once in the run) outside this object
        call self%initDiagnosticArrays()
        !self%solver_pcg%zeroDiagnostics
     
        !   initialize solution
        e_solution = source%e0    
        !   not sure about allocation here -- solution will exist
        !    (and might be allocated) in calling routine, but b is local
        !   Note that rhs will be a TVector of same type as solution
        b = source%rhs

        !   set up source term for divergence correction equations
        if( source%non_zero_source ) then
          !phi0 =  self%model_operator%p   !   make a copy of TScalar using
                                   !   model_operator template
          phi0 = self%divergence_correction%rhsDivCor( self%omega, source )
        endif
        temp = e_solution    !  copy of solution for input to DC
        
        loop: do while ((.not.self%solver_pcg%converged).and.(.not.self%solver_pcg%failed))

          !self%solver_pcg( b, solution )
          !   I am just copying this -- while we work on this should
          !     reconsider implementation
          ! solver_pcg%converged when the relative error is less than tolerance
          self%solver_pcg%converged = self%solver_pcg%niter .lt. self%solver_pcg%max_iter

          ! there are two ways of failing: 1) QMR did not work or
          !        2) total number of divergence corrections exceeded
          self%solver_pcg%failed = self%solver_pcg%failed .or. self%failed

          !  update solver diagnostics 
          do iter = 1, self%solver_pcg%niter
             ! why are we using an explicit loop here?
             self%EMrelErr( self%nIterTotal + iter ) = self%solver_pcg%relErr(iter)
          enddo
          !
          self%nIterTotal = self%nIterTotal + self%solver_pcg%niter

          self%nDivCor = self%nDivCor+1
          if( self%nDivCor < self%MaxDivCor ) then
             ! do divergence correction
             e_solution = temp     !    assuming temp is already allocated,
                                ! don't want to reallocate!
             if( source%non_zero_source ) then
                !e_solution = self%divergence_correction(temp,phi0)
             else
                !e_solution = self%divergence_correction(temp)
             endif
          else
             ! max number of divergence corrections exceeded; convergence solver_pcg%failed
             self%solver_pcg%failed = .true.
         endif
       enddo loop
      
       !   finish up solution--I am omitting boundary values for adjt case --
       !    we never used boundary of adjoint--sensitivity to boundary data,
       !    just use interior part of adjoint soln to compute sensitivity to model
       !    parameters

       ! note that here I am assuming things like mult and add are subroutines.
       !   we need to sort out conventions!   In Solver_QMR I assumed functions,
       !    but I suspect we will be better off just using subroutines in terms
       !    of efficiency
       if( source%adjt ) then
         !call e_solution%mult( e_solution, self%model_operator%Metric%Vedge )
         !   just leave bdry values set to 0
       else
         !call solution%add( e_solution, source%bdry)
       endif

       !  deallocate local objects
       deallocate(temp)
       deallocate(b)
       deallocate(phi0)

   end function getESolutionForwardSolverDC
        
end Module ForwardSolverDC
 
