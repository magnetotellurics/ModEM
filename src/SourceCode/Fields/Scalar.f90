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
            procedure( interface_to_node_scalar ), deferred, public :: toNode
            !
            !> Scalar Routines
            !
            procedure, public :: length => length_Scalar
            !
            procedure, public :: boundary => boundary_Scalar
            procedure, public :: interior => interior_Scalar
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
        !> No interface subroutine briefing
        !
        subroutine interface_to_node_scalar( self, node_scalar, interior_only )
            import :: Scalar_t
            !
            class( Scalar_t ), intent( inout ) :: self, node_scalar
            logical, intent( in ), optional :: interior_only
            !
        end subroutine interface_to_node_scalar
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
    subroutine boundary_Scalar( self, boundary )
        implicit none
        !
        class( Scalar_t ), intent( in ) :: self
        class( Scalar_t ), allocatable, intent( inout ) :: boundary
        !
        complex( kind=prec ), allocatable, dimension(:) :: c_array
        !
        allocate( boundary, source = self )
        !
        c_array = boundary%getArray()
        !
        c_array( self%indInterior() ) = C_ZERO
        !
        call boundary%setArray( c_array )
        !
    end subroutine boundary_Scalar
    !
    !> No subroutine briefing
    !
    subroutine interior_Scalar( self, interior )
        implicit none
        !
        class( Scalar_t ), intent( in ) :: self
        class( Scalar_t ), allocatable, intent( inout ) :: interior
        !
        complex( kind=prec ), allocatable, dimension(:) :: c_array
        !
        allocate( interior, source = self )
        !
        c_array = interior%getArray()
        !
        c_array( self%indBoundary()) = C_ZERO
        !
        call interior%setArray( c_array )
        !
    end subroutine interior_Scalar
    !
end module Scalar
