! *************
! 
! Derived class to define a MT Source with boundary data computed by 2D solutions
! 
! Last modified at 10/11/2021 by Paulo Werdt
! 
! *************
! 
!**
module SourceMT_2D
    !
    use Constants
    use Source
    use cVector3D_SG
    use Grid3D_SG
    use Forward2D
    !
    type, extends( Source_t ) :: SourceMT_2D_t
          !
          ! PROPERTIES HERE
          !
          contains
             !
             final :: SourceMT_2D_dtor
             !
             procedure, public :: setRHS => setRHSMT_2D
             procedure, public :: setE    => setESourceMT_2D
             !
    end type SourceMT_2D_T
    !
    interface SourceMT_2D_t
         module procedure SourceMT_2D_ctor
    end interface SourceMT_2D_t
    !
contains
    !
    ! SourceMT_2D constructor
    function SourceMT_2D_ctor( model_operator, model_parameter, E ) result( self )
         !
         implicit none
         !
         class( ModelOperator_t ), target, intent( in )  :: model_operator
         class( ModelParameter_t ), target, intent( in ) :: model_parameter
         class( cVector_t ), intent( in ), optional     :: E
         !
         type( SourceMT_2D_t ) :: self
         !
         !write(*,*) "Constructor SourceMT_2D_t"
         !
         self%model_operator => model_operator
         self%model_parameter => model_parameter
         !
         self%non_zero_source = .FALSE.
         self%adjt = .FALSE.
         !
         if ( present( E ) ) then
             !
             allocate( self%E, source = E )
             !
             call self%setRHS()
             !
         endif
         !
    end function SourceMT_2D_ctor
    !
    ! SourceMT_1D destructor
    subroutine SourceMT_2D_dtor( self )
        !
        implicit none
        !
        type( SourceMT_2D_t ), intent( inout ) :: self
        !
        !write(*,*) "Destructor SourceMT_2D_t"
        !
        call self%dealloc()
        !
    end subroutine SourceMT_2D_dtor
    !
    ! Set self%E from forward modelling 2D
    subroutine setESourceMT_2D( self, omega, polarization )
      implicit none
      !
      class( SourceMT_2D_t ), intent( inout ) :: self
      real( kind=prec ), intent( in )         :: omega
      integer, intent( in )                   :: polarization
      !
      ! local variables
      ! 2D grid definitions
      type( Grid2D_t )           :: grid_2D
      !
      type( Esoln2DTM_t )        :: E2D
      !
      type( ModelParameter2D_t ) :: model_parameter_2D
      !
      type( Forward2D_t )        :: forward_2D
      !
      integer :: i_slice, n_slice
      !
      ! Ex-polarization
      if( polarization == 1 ) then
          !
          n_slice = self%model_operator%metric%grid%ny
          !
          grid_2D = Grid2D_t( self%model_operator%metric%grid%ny, &
          self%model_operator%metric%grid%nzAir, self%model_operator%metric%grid%nzEarth, &
          self%model_operator%metric%grid%dy, self%model_operator%metric%grid%dz )
      !
      ! Ey-polarization
      else
          !
          n_slice = self%model_operator%metric%grid%nx
          !
          grid_2D = Grid2D_t( self%model_operator%metric%grid%nx, &
          self%model_operator%metric%grid%nzAir, self%model_operator%metric%grid%nzEarth, &
          self%model_operator%metric%grid%dx, self%model_operator%metric%grid%dz )
          !
      endif
      !
      !
      do i_slice = 1, n_slice
          !
          model_parameter_2D = self%model_parameter%Slice2D( polarization, i_slice )
          !
          forward_2D = Forward2D_t( model_parameter_2D )
          !
          E2D = Esoln2DTM_t( grid_2D )
          !
          write(*,*) "Solve2D: ", polarization, i_slice
          !
          ! Solve 2D and store the result in E2D structure
          call forward_2D%solve( E2D )
          !
          ! Construct self%E => E3D
          select type( grid => self%model_operator%metric%grid )
          class is( Grid3D_SG_t )
              !
              ! Allocate self%E => cVector3D_SG
              if( allocated( self%E ) ) deallocate( self%E )
                  allocate( self%E, source = cVector3D_SG_t( grid, EDGE ) )
              !
              ! Fill E3D (cVector3D_SG) from E2D (Esoln2DTM_t)
              select type( E3D => self%E )
              class is( cVector3D_SG_t )
                  !
                  ! 1st polarization case: Only x components are non-zero
                  if( polarization == 1  ) then
                      E3D%x( :, i_slice, : ) = E2D%Ey( :, : )
                  !
                  ! 2nd polarization case: Only y components are non-zero
                  else
                      E3D%y( i_slice, :, : ) = E2D%Ey( :, : )
                  endif
                  !
                  E3D%z( i_slice, :, : ) = E2D%Ez( :, : )
                  !
              class default
                  stop "SourceMT_2D_t:setESourceMT_2D: Unclassified E3D"
              end select
              !
          class default
              stop "SourceMT_2D_t:setESourceMT_2D: Unclassified grid"
          end select
          !
      enddo
      !
    end subroutine setESourceMT_2D
    !
    ! Set RHS from self%E
    subroutine setRHSMT_2D( self )
        implicit none
        !
        class( SourceMT_2D_t ), intent( inout ) :: self
        !!
        call self%model_operator%MultAib( self%E%Boundary(), self%rhs )
        !
        self%rhs = self%rhs * C_MinusOne
        !
    end subroutine setRHSMT_2D
    !
end module SourceMT_2D
