! *************
! 
! Derived class to define a MT Source with boundary data computed by 1D solutions
! 
! Last modified at 10/11/2021 by Paulo Werdt
! 
! *************
! 
!**
module SourceMT_1D
   !
   use Constants
   use cVector3D_SG
   use Grid3D_SG
   use Source
   use ModelOperator
   use ModelParameter1D
   use Forward1D
   !
   type, extends( Source_t ) :: SourceMT_1D_t
        !
        ! PROPERTIES HERE
        !
        contains
          !
          final :: SourceMT_1D_dtor
          !
          procedure, public :: setRHS => setRHSMT_1D
          procedure, public :: setE0  => setE0MT_1D
          procedure, public :: setE   => setESourceMT_1D
          !
   end type SourceMT_1D_T
   !
   interface SourceMT_1D_t
       module procedure SourceMT_1D_ctor
   end interface SourceMT_1D_t
   !
contains
   !
   ! SourceMT_1D constructor
   function SourceMT_1D_ctor( model_operator, model_parameter, E ) result( self )
       implicit none
       !
       class( ModelOperator_t ), target, intent( in )  :: model_operator
       class( ModelParameter_t ), target, intent( in ) :: model_parameter
       class( cVector_t ), intent( in ), optional      :: E
       !
       type( SourceMT_1D_t ) :: self
       !
       call self%init()
       !
       self%model_operator => model_operator
       self%model_parameter => model_parameter
       !
       if ( present( E ) ) then
          !
          allocate( self%E, source = E )
          !
          call self%setRHS()
          !
          call self%setE0()
          !
       endif
       !
   end function SourceMT_1D_ctor
   !
   ! SourceMT_1D destructor
   subroutine SourceMT_1D_dtor( self )
      implicit none
      !
      type( SourceMT_1D_t ), intent( inout ) :: self
      !
      !write(*,*) "Destructor SourceMT_1D_t"
      !
      call self%dealloc()
      !
   end subroutine SourceMT_1D_dtor
   !
   ! Set self%E from forward modelling 1D
   subroutine setESourceMT_1D( self, omega, polarization )
     implicit none
     !
     class( SourceMT_1D_t ), intent( inout )  :: self
     real( kind=prec ), intent( in )          :: omega
     integer, intent( in )                    :: polarization
     !
     class( ModelParameter1D_t ), allocatable        :: model_parameter_1D
     class( Forward1D_t ), allocatable               :: forward_1D
     complex( kind=prec ), allocatable, dimension(:) :: E1D
     !
     integer :: ix, iy
     !
     ! Get Model1D from average conductivity 3D
     model_parameter_1D = self%model_parameter%AvgModel1D()
     !
     forward_1D = Forward1D_t( model_parameter_1D )
     !
     call forward_1D%SetFreq( omega )
     !
     if( allocated( E1D ) ) deallocate( E1D )
     allocate( E1D( self%model_operator%grid%nz ) )
     !
     ! Solve 1D and store the result in E1D structure
     call forward_1D%solve( E1D )
     !
     do ix = 1, self%model_operator%grid%nx
       !
       do iy = 1, self%model_operator%grid%ny
           !
           ! Allocate and construct self%E => E3D
           select type( grid => self%model_operator%grid )
           class is( Grid3D_SG_t )
              !
              if( .not. allocated( self%E ) ) then
                 !
                 allocate( self%E, source = cVector3D_SG_t( grid, EDGE ) )
                 !
                 call self%E%zeros()
              endif
              !
              ! Fill E3D (cVector3D_SG) from E1D (Esoln1DTM_t)
              select type( E3D => self%E )
              class is( cVector3D_SG_t )
                 !
                 ! 1st polarization case: Only x components are non-zero
                 if( polarization == 1 ) then
                    E3D%x( ix, iy, : ) = E1D
                 !
                 ! 2nd polarization case: Only y components are non-zero
                 else
                    E3D%y( ix, iy, : ) = E1D
                 endif
                 !
                 E3D%z( ix, iy, : ) = E1D
                 !
              end select
              !
           end select
           !
        enddo
        !
     enddo
     !
     deallocate( E1D )
     deallocate( model_parameter_1D )
     deallocate( forward_1D )
     !
     !call self%E%print
     !
     call self%setRHS()
     !
     call self%setE0()
     !
   end subroutine setESourceMT_1D
   !
   ! Set RHS from self%E
   subroutine setRHSMT_1D( self )
      implicit none
      !
      class( SourceMT_1D_t ), intent(inout) :: self
      !
      class( cVector_t ), allocatable :: bdry
      !
      allocate( bdry, source = self%E%Boundary() )
      allocate(self%bdry, source = bdry)   !   temporary debugging ...
                        !   or should we save bdry in source object?
      !
      if( allocated( self%rhs ) ) deallocate( self%rhs )
      !
      select type( E => self%E )
      class is( cVector3D_SG_t )
         !
         allocate( self%rhs, source = cVector3D_SG_t( E%grid, EDGE ) )
         !
         call self%model_operator%MultAib( bdry, self%rhs )
         call self%rhs%print(441,'-RHS')
         !
         self%rhs = C_MinusOne * self%rhs
         !
      end select
      !
      deallocate( bdry )
      !
   end subroutine setRHSMT_1D
   !
   ! Set e0 from self%E
   subroutine setE0MT_1D( self )
      implicit none
      !
      class( SourceMT_1D_t ), intent( inout ) :: self
      !
      if( allocated( self%e0 ) ) deallocate( self%e0 )
      allocate( self%e0, source = self%E%Interior() )
      !
   end subroutine setE0MT_1D
   !
end module SourceMT_1D
