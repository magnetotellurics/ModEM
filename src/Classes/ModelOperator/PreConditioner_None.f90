!**
! This is a null preconditioner -- identity matrix
! which is only implemented for CSG.
!*
module PreConditioner_None
   !
   use Constants
   use Grid
   use cVector
   use cVector3D_SG
   use cScalar3D_SG
   use ModelOperator_MF
   use PreConditioner
   !
   type, extends( PreConditioner_t ), public :: PreConditioner_None_t
      !
      real( kind=prec ) :: omega = 0.0
      !
      contains
         !
         final :: PreConditioner_None_dtor
         !
         ! Main routines used externally
         procedure, public :: setPreConditioner => setPreConditioner_None ! This needs to be called by Solver   object
         !
         procedure, public :: LTSolve => LTSolvePreConditioner_None ! These are left (M1) and right (M2)
         procedure, public :: UTSolve => UTSolvePreConditioner_None ! preconditioning matrices for curl-curl equation.
         procedure, public :: LUSolve => LUSolvePreConditioner_None ! preconditoner for symmetric divCgrad operator
         !
   end type PreConditioner_None_t
   !
   interface PreConditioner_None_t
      module procedure PreConditioner_None_ctor
   end interface PreConditioner_None_t
   !
contains
   !**
   ! Class constructor
   !*
   function PreConditioner_None_ctor( ) result( self ) 
      implicit none
      !
      ! class( ModelOperator_MF_t ), target, intent( in ) :: model_operator
      type( PreConditioner_None_t ) :: self
      !
      write(*,*) "Constructor PreConditioner_None_t"
      !
      ! self%model_operator => model_operator  !   dont really need model_operator?
      !
   end function PreConditioner_None_ctor
   !
   ! PreConditioner_MF_CC destructor
   subroutine PreConditioner_None_dtor( self )
     implicit none
     !
     type( PreConditioner_None_t ), intent( inout ) :: self
     !
   end subroutine PreConditioner_None_dtor
   !
   !**
   ! SetPreConditioner
   !*
   subroutine setPreConditioner_None( self, omega )
      implicit none
      !
      class( PreConditioner_None_t ), intent( inout )  :: self
      real( kind=prec ), intent( in )                  :: omega
      !
      ! Really nothing to set up here!
      ! self%omega = omega
   end subroutine setPreConditioner_None
   
   !**
   ! Purpose: to solve the lower triangular system (or it"s adjoint);
   ! for the d-ilu pre-condtioner. -- for no precond, just return input
   !    as output
   !*
   subroutine LTSolvePreConditioner_None( self, inE, outE, adjt )
      implicit none
      !
      class( PreConditioner_None_t ), intent( inout ) :: self
      class( cVector_t ), intent( in )                 :: inE
      class( cVector_t ), intent( inout )              :: outE
      logical, intent( in )                            :: adjt
      !
      select type( inE )
      class is( cVector3D_SG_t )
         !
         select type( outE )
         class is( cVector3D_SG_t )
            if (.not.outE%isAllocated) then
               stop "outE in LUsolve not allocated yet"
            end if
            !
            outE = inE
         end select
      end select
      !
   end subroutine LTSolvePreConditioner_None
   
   !**
   ! Purpose: to solve the upper triangular system (or it"s adjoint);
   ! for the d-ilu pre-condtioner
   !*
   subroutine UTSolvePreConditioner_None( self, inE, outE, adjt )
      implicit none
      !
      class( PreConditioner_None_t ), intent( inout ) :: self
      class( cVector_t ), intent( in )                 :: inE
      class( cVector_t ), intent( inout )              :: outE
      logical, intent( in )                            :: adjt
      !
      integer :: ix, iy, iz
      !
      select type( inE )
      class is( cVector3D_SG_t )
         !
         select type( outE )
         class is( cVector3D_SG_t )
            if (.not.outE%isAllocated) then
               stop "outE in LUsolve not allocated yet"
            end if
            !
            outE = inE
         end select
      end select
      !
   end subroutine UTSolvePreConditioner_None
   !**
   ! No precond for scalar symmetric problems
   !*
   subroutine LUSolvePreConditioner_None( self, inPhi, outPhi )
      implicit none
      !
      class( PreConditioner_None_t ), intent( inout ) :: self
      class( cScalar_t ), intent( in )                 :: inPhi
      class( cScalar_t ), intent( inout )              :: outPhi
      !
      select type(inPhi)
       class is(cScalar3D_SG_t)
          select type(outPhi)
          class is(cScalar3D_SG_t)
             if (.not.outPhi%isAllocated) then
                stop "outPhi in LUsolve not allocated yet"
             end if
             outPhi = inPhi
          end select
       end select
             

      !
   end subroutine LUSolvePreConditioner_None
   !
end module PreConditioner_None

