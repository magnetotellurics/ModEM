!**
! This is for preconditioning the divergence correction equations
!*
module PreConditioner_MF_DC
   !
   use Constants
   use Grid
   use cScalar
   use cScalar3D_SG
   use ModelOperator_MF
   use PreConditioner
   !
   type, extends( PreConditioner_t ) :: PreConditioner_MF_DC_t
       !
       type( cScalar3D_SG_t ) :: d       !    diagonal of D-ILU used for
       !
    contains
       !
       final :: PreConditioner_MF_DC_dtor
       !
       ! Main routines used externally
       procedure, public :: create => createPreConditioner_MF_DC
!      procedure, public :: allocate => allocatePreConditioner_MF_DC
       procedure, public :: deallocate => deallocatePreConditioner_MF_DC
       procedure, public :: setPreConditioner => setPreConditioner_MF_DC ! This needs to be called by Solver   object
                                                                  ! every time any part of the operator changes.
       !   I guess these declarations are not needed                                                          
       procedure, public :: LTSolve => LTSolvePreConditioner_MF_DC ! dummy
       procedure, public :: UTSolve => UTSolvePreConditioner_MF_DC ! dummy
       procedure, public :: LUSolve => LUSolvePreConditioner_MF_DC ! preconditoner for symmetric divCgrad operator
       
   end type PreConditioner_MF_DC_t
   !
   interface PreConditioner_MF_DC_t
       module procedure PreConditioner_MF_DC_ctor
   end interface PreConditioner_MF_DC_t
   !
contains
   !**
   ! Class constructor
   !*
   function PreConditioner_MF_DC_ctor( model_operator ) result( self ) 
      implicit none
      !
      class( ModelOperator_MF_t ), target, intent( in ) :: model_operator
      type( PreConditioner_MF_DC_t ) :: self
      !
      !write(*,*) "Constructor PreConditioner_MF_DC_t"
      !
      self%model_operator => model_operator
      !
      call self%create()
      !
   end function PreConditioner_MF_DC_ctor
   !
   ! PreConditioner_MF_CC destructor
   subroutine PreConditioner_MF_DC_dtor( self )
     implicit none
     !
     type( PreConditioner_MF_DC_t ), intent( inout ) :: self
     !
     !write(*,*) "Destructor PreConditioner_MF_DC_t"
     !
     !call self%dealloc()
     !
   end subroutine PreConditioner_MF_DC_dtor
   !
   !**
   ! createPreConditioner_MF_DC
   !*
   subroutine createPreConditioner_MF_DC( self )
      implicit none
      !
      class( PreConditioner_MF_DC_t ), intent( inout ) :: self
      !
      !      can"t we just use this???
      self%d = self%model_operator%createScalar()
    !
      !      call self%allocate()
   end subroutine createPreConditioner_MF_DC
   
   !**
   ! allocatePreConditioner_MF_DC
   !*
   !subroutine allocatePreConditioner_MF_DC(self)
   !   class(PreConditioner_MF_DC_t), intent(inout) :: self
      !
    !select type( grid => self%model_operator%grid )
   !         class is( Grid3D_SG_t )
   !               allocate( self%Dilu, source = cVector3D_SG_t( grid, EDGE ) )
   !   end select
    !!!
   !end subroutine allocatePreConditioner_MF_DC

   !**
   ! DeAllcoate
   !*
   subroutine deallocatePreConditioner_MF_DC(self)
      class(PreConditioner_MF_DC_t) :: self
   end subroutine deallocatePreConditioner_MF_DC
   
   !**
   ! SetPreConditioner -- could be an abstract routine, but in the CC case
   !      we pass omega as a parameter, and that is not relevant here -- but since
   !    omega is a property of that class could set, and not pass into this procedure explicitly
   !*
   subroutine setPreConditioner_MF_DC( self, omega )
      implicit none
      !
      class( PreConditioner_MF_DC_t ), intent( inout ) :: self
      real( kind=prec ), intent( in )                  :: omega
      !
      integer :: ix,iy,iz
      !
      ! Compute inverse diagonal elements for D-ILU (interior nodes only)
      ! set top nodes to 1.0
      self%d%v(1,:,:) = 1.0
      self%d%v(:,1,:) = 1.0
      self%d%v(:,:,1) = 1.0
     !
      do iz = 2, self%model_operator%nz
         do iy = 2, self%model_operator%ny
            do ix = 2, self%model_operator%nx
               self%d%v(ix, iy, iz) = self%model_operator%c%v(ix, iy, iz) - &
               self%model_operator%db1%x(ix,iy,iz)*self%model_operator%db2%x(ix-1,iy,iz) * &
               self%d%v(ix-1,iy,iz)- &
               self%model_operator%db1%y(ix,iy,iz)*self%model_operator%db2%y(ix,iy-1,iz) * &
               self%d%v(ix,iy-1,iz)- &
               self%model_operator%db1%z(ix,iy,iz)*self%model_operator%db2%z(ix,iy,iz-1) * &
               self%d%v(ix,iy,iz-1)
               self%d%v(ix, iy, iz) = 1.0/ self%d%v(ix, iy, iz)
            end do
         end do
      end do
     !
   end subroutine setPreConditioner_MF_DC
   !**
   ! LTsolve and UTsolve are in abstract class and must be defined -- but not used for DC which
   !      this object will be used -- so just dummies here
   !*
   subroutine LTSolvePreConditioner_MF_DC(self, inE, outE, adjt)
      implicit none
      !
      class( PreConditioner_MF_DC_t ), intent( inout ) :: self
      class( cVector_t ), intent( in )                 :: inE
      class( cVector_t ), intent( inout )              :: outE
      logical, intent( in )                            :: adjt
      
      STOP "ERROR: LTsolve not coded for this pre-conditioner class"
      !
   end subroutine LTSolvePreConditioner_MF_DC
 !* 
   subroutine UTSolvePreConditioner_MF_DC( self, inE, outE, adjt )
      implicit none
      !
      class( PreConditioner_MF_DC_t ), intent( inout ) :: self
      class( cVector_t ), intent( in )                 :: inE
      class( cVector_t ), intent( inout )              :: outE
      logical, intent( in )                            :: adjt
      
      STOP "ERROR: UTsolve not coded for this pre-conditioner class"
      !
   end subroutine UTSolvePreConditioner_MF_DC
   !**
   ! apply pre-conditioner, LU solve
   !*
   subroutine LUSolvePreConditioner_MF_DC( self, inPhi, outphi )
      implicit none
      !
      class( PreConditioner_MF_DC_t ), intent( inout ) :: self
      class( cScalar_t ), intent( in )                 :: inPhi
      class( cScalar_t ), intent( inout )              :: outPhi
        !    no adjoint for this case .. ever as far as I can see!
      integer :: ix, iy, iz

      select type(inPhi)
      class is(cScalar3D_SG_t)
         select type(outPhi)
         class is(cScalar3D_SG_t)
            if (.not.outPhi%isAllocated) then
               stop "outPhi in LUsolve not allocated yet"
            end if

            call outPhi%zeros()
            ! forward substitution (Solve lower triangular system)
            ! the coefficients are only for the interior nodes
            do iz = 2, inPhi%nz
               do iy = 2, inPhi%ny
                  do ix = 2, inPhi%nx
                     outPhi%v(ix, iy, iz) = inPhi%v(ix, iy, iz) &
                     - outPhi%v(ix-1,iy,iz)*self%model_operator%db1%x(ix,iy,iz)&
                     *self%d%v(ix-1,iy,iz) &
                     - outPhi%v(ix,iy-1,iz)*self%model_operator%db1%y(ix,iy,iz)&
                     *self%d%v(ix,iy-1,iz) &
                     - outPhi%v(ix,iy,iz-1)*self%model_operator%db1%z(ix,iy,iz)&
                     *self%d%v(ix,iy,iz-1)
                  end do
               end do
            end do

            ! backward substitution (Solve upper triangular system)
            ! the coefficients are only for the interior nodes
            do iz = inPhi%nz,2,-1
               do iy = inPhi%ny,2,-1
                  do ix = inPhi%nx,2,-1
                     outPhi%v(ix, iy, iz) = (outPhi%v(ix, iy, iz)   &
                     - outPhi%v(ix+1, iy, iz)*self%model_operator%db2%x(ix, iy, iz)   &
                     - outPhi%v(ix, iy+1, iz)*self%model_operator%db2%y(ix, iy, iz)   &
                     - outPhi%v(ix, iy, iz+1)*self%model_operator%db2%z(ix, iy, iz)) &
                     *self%d%v(ix, iy, iz)
                  end do
               end do
            end do
         end select
      end select
      !
   end subroutine LUSolvePreConditioner_MF_DC
   !
end module PreConditioner_MF_DC

