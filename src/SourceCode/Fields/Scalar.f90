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
            procedure( interface_sum_cell_scalar ), deferred, public :: sumCell
            !
            !> Scalar Routines
            !
            procedure, public :: length => length_Scalar
            !
            procedure, public :: boundary => boundary_Scalar
            procedure, public :: interior => interior_Scalar
            !
            procedure, public :: getArray => getArray_Scalar
            procedure, public :: setArray => setArray_Scalar
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
            complex( kind=prec ), allocatable :: v(:,:,:)
        end function interface_get_v_scalar
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_v_scalar( self, v )
            import :: Scalar_t, prec
            !
            class( Scalar_t ), intent( inout ) :: self
            complex( kind=prec ), dimension(:,:,:), intent( in ) :: v
        end subroutine interface_set_v_scalar
        !
        !> No interface subroutine briefing
        !
        subroutine interface_sum_cell_scalar( self, node_out, interior_only )
            import :: Scalar_t
            !
            class( Scalar_t ), intent( inout ) :: self
            class( Scalar_t ), allocatable, intent( out ) :: node_out
            logical, intent( in ), optional :: interior_only
            !
        end subroutine interface_sum_cell_scalar
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
        c_array( self%ind_interior ) = C_ZERO
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
        c_array( self%ind_boundary ) = C_ZERO
        !
        call interior%setArray( c_array )
        !
    end subroutine interior_Scalar
    !
    !> No subroutine briefing
    !
    function getArray_Scalar( self ) result( array )
        implicit none
        !
        class( Scalar_t ), intent( in ) :: self
        complex( kind=prec ), allocatable, dimension(:) :: array
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "getArray_Scalar > self not allocated." )
        endif
        !
        if( self%store_state .EQ. compound ) then
            !
            allocate( array( self%length() ) )
            array = (/reshape( self%getV(), (/self%Nxyz, 1/))/)
            !
        elseif( self%store_state .EQ. singleton ) then
            !
            array = self%getSV()
            !
        else
            call errStop( "getArray_Scalar > Unknown store_state!" )
        endif
        !
    end function getArray_Scalar
    !
    !> No subroutine briefing
    !
    subroutine setArray_Scalar( self, array )
        implicit none
        !
        class( Scalar_t ), intent( inout ) :: self
        complex( kind=prec ), dimension(:), intent( in ) :: array
        !
        complex( kind=prec ), allocatable, dimension(:,:,:) :: v
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "setArray_Scalar > self not allocated." )
        endif
        !
        call self%deallOtherState
        !
        if( self%store_state .EQ. compound ) then
            !
            v = reshape( array, (/self%NdV(1), self%NdV(2), self%NdV(3)/) )
            !
            call self%setV( v )
            !
        elseif( self%store_state .EQ. singleton ) then
            !
            call self%setSV( array )
            !
        else
            call errStop( "setArray_Scalar > Unknown store_state!" )
        endif
        !
    end subroutine setArray_Scalar
    !
end module Scalar
