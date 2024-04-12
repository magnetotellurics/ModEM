!
!> Derived class to define a MT Source with boundary data computed by 2D solutions
!
module SourceMT_2D
    ! !
    ! use Constants
    ! use Source
    ! use cVector3D_SG
    ! use Grid3D_SG
    ! use Forward2D
    ! !
    ! type, extends( Source_t ) :: SourceMT_2D_t
        ! !
        ! !> No derived properties
        ! !
        ! contains
            ! !
            ! final :: SourceMT_2D_dtor
            ! !
            ! procedure, public :: createE => createE_SourceMT_2D
            ! !
            ! procedure, public :: createRHS => createRHS_SourceMT_2D
            ! !
    ! end type SourceMT_2D_T
    ! !
    ! interface SourceMT_2D_t
         ! module procedure SourceMT_2D_ctor
    ! end interface SourceMT_2D_t
    ! !
! contains
    ! !
    ! !> SourceMT_2D constructor
    ! function SourceMT_2D_ctor( model_operator, sigma, period ) result( self )
        ! implicit none
        ! !
        ! class( ModelOperator_t ), target, intent( in ) :: model_operator
        ! class( ModelParameter_t ), target, intent( in ) :: sigma
        ! real( kind=prec ), intent( in ) :: period
        ! !
        ! type( SourceMT_2D_t ) :: self
        ! !
        ! !write( *, * ) "Constructor SourceMT_2D_t"
        ! !
        ! call self%baseInit
        ! !
        ! self%model_operator => model_operator
        ! self%sigma => sigma
        ! !
        ! self%period = period
        ! !
    ! end function SourceMT_2D_ctor
    ! !
    ! !> Deconstructor routine:
    ! !>     Call the base routine baseDealloc().
    ! subroutine SourceMT_2D_dtor( self )
        ! !
        ! implicit none
        ! !
        ! type( SourceMT_2D_t ), intent( inout ) :: self
        ! !
        ! !write( *, * ) "Destructor SourceMT_2D_t"
        ! !
        ! call self%baseDealloc
        ! !
    ! end subroutine SourceMT_2D_dtor
    ! !
    ! !> Set self%E from forward modelling 2D
    ! subroutine createE_SourceMT_2D( self )
        ! implicit none
        ! !
        ! class( SourceMT_2D_t ), intent( inout ) :: self
        ! !
        ! !> 2D grid definitions
        ! type( Grid2D_t ) :: grid_2D
        ! !
        ! type( Esoln2DTM_t ) :: E2D
        ! !
        ! type( ModelParameter2D_t ) :: model_parameter_2D
        ! !
        ! type( Forward2D_t ) :: forward_2D
        ! !
        ! integer :: i_slice, n_slice, pol
        ! !
        ! allocate( cVector3D_SG_t :: self%E( 2 ) )
        ! !
        ! ! Loop over polarizations
        ! do pol = 1, 2
            ! !
            ! !> Ex-polarization
            ! if( pol == 1 ) then
                ! !
                ! n_slice = self%model_operator%metric%grid%ny
                ! !
                ! grid_2D = Grid2D_t( self%model_operator%metric%grid%ny, &
                ! self%model_operator%metric%grid%nzAir, self%model_operator%metric%grid%nzEarth, &
                ! self%model_operator%metric%grid%dy, self%model_operator%metric%grid%dz )
            ! !
            ! !> Ey-polarization
            ! else
                ! !
                ! n_slice = self%model_operator%metric%grid%nx
                ! !
                ! grid_2D = Grid2D_t( self%model_operator%metric%grid%nx, &
                ! self%model_operator%metric%grid%nzAir, self%model_operator%metric%grid%nzEarth, &
                ! self%model_operator%metric%grid%dx, self%model_operator%metric%grid%dz )
                ! !
            ! endif
            ! !
            ! ! Loop over model slices
            ! do i_slice = 1, n_slice
                ! !
                ! model_parameter_2D = self%sigma%slice2D( pol, i_slice )
                ! !
                ! forward_2D = Forward2D_t( model_parameter_2D )
                ! !
                ! E2D = Esoln2DTM_t( grid_2D )
                ! !
                ! write( *, * ) "Solve2D: ", pol, i_slice
                ! !
                ! !> Solve 2D and store the result in E2D structure
                ! call forward_2D%solve( E2D )
                ! !
                ! self%E( pol ) = cVector3D_SG_t( self%model_operator%metric%grid, EDGE )
                ! !
                ! !> Fill e_vector (cVector3D_SG) from E2D (Esoln2DTM_t)
                ! select type( E => self%E( pol ) )
                    ! class is( cVector3D_SG_t )
                        ! !
                        ! !> 1st polarization case: Only x components are non-zero
                        ! if( pol == 1  ) then
                            ! E%x( :, i_slice, : ) = E2D%Ey( :, : )
                        ! !
                        ! !> 2nd polarization case: Only y components are non-zero
                        ! else
                            ! E%y( i_slice, :, : ) = E2D%Ey( :, : )
                        ! endif
                        ! !
                        ! E%z( i_slice, :, : ) = E2D%Ez( :, : )
                        ! !
                    ! class default
                        ! stop "SourceMT_2D_t:createE_SourceMT_2D: Unclassified e_vector"
                ! end select
                ! !
            ! enddo
            ! !
      ! enddo
      ! !
    ! end subroutine createE_SourceMT_2D
    ! !
    ! !> Set RHS from self%E
    ! subroutine createRHS_SourceMT_2D( self )
        ! implicit none
        ! !
        ! class( SourceMT_2D_t ), intent( inout ) :: self
        ! !
        ! class( Vector_t ), allocatable :: e_boundary
        ! !
        ! integer :: pol
        ! !
        ! allocate( cVector3D_SG_t :: self%rhs( 2 ) )
        ! !
        ! do pol = 1, 2
            ! !
            ! call self%E( pol )%boundary( e_boundary )
            ! !
            ! self%rhs( pol ) = cVector3D_SG_t( self%model_operator%metric%grid, EDGE )
            ! !
            ! call self%model_operator%MultAib( e_boundary, self%rhs( pol ) )
            ! !
            ! deallocate( e_boundary )
            ! !
            ! call self%rhs( pol )%mult( C_MinusOne )
            ! !
        ! enddo
        ! !
    ! end subroutine createRHS_SourceMT_2D
    ! !
end module SourceMT_2D
!