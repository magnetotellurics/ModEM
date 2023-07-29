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
    end type Scalar_t
    !
    !> Allocatable Scalar element for Old Fortran polymorphic Arrays!!!
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
            class( Scalar_t ), intent( in ) :: self
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
end module Scalar
