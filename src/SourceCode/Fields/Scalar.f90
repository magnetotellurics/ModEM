!
!> Abstract class to define a Scalar 
!
module Scalar
    !
    use Field
    !
    !> Abstract base class
    type, abstract, extends( Field_t ) :: Scalar_t
        !
        integer, dimension(3) :: NdV
        !
        integer :: Nxyz
        !
    contains
            !
            !> Scalar Interfaces
            procedure( interface_get_v_scalar ), deferred, public :: getV
            procedure( interface_set_v_scalar ), deferred, public :: setV
            !
            !> Scalar Routines
            procedure, public :: length => length_Scalar
            !
            procedure, public :: switchStoreState => switchStoreState_Scalar
            !
    end type Scalar_t
    !
    !> Allocatable Scalar element for Old Fortran polymorphism on Arrays!!!
    type, public :: GenScalar_t
        !
        class( Scalar_t ), allocatable :: s
        !
    end type GenScalar_t
    !
    !>
    abstract interface
        !
        !> No interface function briefing
        !
        function interface_get_v_scalar( self ) result( v )
            import :: Scalar_t, prec
            !
            class( Scalar_t ), intent( inout ) :: self
            complex( kind=prec ), allocatable :: v(:, :, :)
        end function interface_get_v_scalar
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_v_scalar( self, v )
            import :: Scalar_t, prec
            !
            class( Scalar_t ), intent( inout ) :: self
            complex( kind=prec ), allocatable, intent( in ) :: v(:, :, :)
        end subroutine interface_set_v_scalar
        !
    end interface
    !
contains
    !
    !> No subroutine briefing
    !
    function length_Scalar( self ) result( field_length )
        implicit none
        !
        class( Scalar_t ), intent( in ) :: self
        !
        integer :: field_length
        !
        field_length = self%Nxyz
        !
    end function length_Scalar
    !
    !> No subroutine briefing
    !
    subroutine switchStoreState_Scalar( self, store_state )
        implicit none
        !
        class( Scalar_t ), intent( inout ) :: self
        integer, intent( in ), optional :: store_state
        !
        integer :: nzAir
        complex( kind=prec ), allocatable, dimension(:,:,:) :: v
        complex( kind=prec ), allocatable, dimension(:) :: s_v
        !
        !> If input state is present...
        if( present( store_state ) ) then
            !
            !> ... and is different of the actual Scalar state: flip it!
            if( self%store_state /= store_state ) then
                call self%switchStoreState
            endif
            !
        else
            !
            select case( self%store_state )
                !
                case( compound )
                    !
                    allocate( s_v( self%length() ) )
                    !
                    s_v = (/reshape( self%getV(), (/self%Nxyz, 1/))/)
                    !
                    call self%setSV( s_v )
                    !
                case( singleton )
                    !
                    if( self%grid_type == NODE ) then
                        !
                        allocate( v( self%nx + 1, self%ny + 1, self%nz + 1 ) )
                        !
                    elseif( self%grid_type == CELL ) then
                        !
                        allocate( v( self%nx, self%ny, self%nz ) )
                        !
                    elseif( self%grid_type == CELL_EARTH ) then
                        !
                        call self%grid%getDimensions( self%nx, self%ny, self%nz, nzAir )
                        !
                        allocate( v( self%nx, self%ny, self%nz - nzAir ) )
                        !
                    else
                         call errStop( "switchStoreState_Scalar > unrecognized grid_type: ["//self%grid_type//"]" )
                    endif
                    !
                    s_v = self%getSV()
                    !
                    v = reshape( s_v, (/self%NdV(1), self%NdV(2), self%NdV(3)/) )
                    !
                    call self%setV( v )
                    !
                case default
                    call errStop( "switchStoreState_Scalar > store_state should be 'singleton' or 'compound'" )
                !
            end select
            !
        endif
        !
    end subroutine switchStoreState_Scalar
    !
end module Scalar
