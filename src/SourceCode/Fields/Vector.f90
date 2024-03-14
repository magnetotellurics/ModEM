!
!> Abstract class to define an abstract Vector field
!
module Vector
    !
    use rScalar3D_SG
    !
    type, abstract, extends( Field_t ) :: Vector_t
        !
        !> No derived properties
        !
    contains
        !
        !> Vector Interfaces
        procedure( interface_get_axis_vector ), deferred, public :: getAxis
        !
        procedure( interface_diag_mult_vector ), deferred, public :: diagMult
        procedure( interface_interp_func_vector ), deferred, public :: interpFunc
        !
        procedure( interface_sum_edge_vector ), deferred, public :: sumEdge
        procedure( interface_sum_edge_vti_vector ), deferred, public :: sumEdgeVTI
        generic :: sumEdges => sumEdge, sumEdgeVTI
        !
        procedure( interface_sum_cells_vector ), deferred, public :: sumCell
        procedure( interface_sum_cells_VTI_vector ), deferred, public :: sumCellVTI
        generic :: sumCells => sumCell, sumCellVTI
		!
		procedure( interface_get_real_vector ), deferred, public :: getReal
		!
        procedure, public :: boundary => boundary_Vector
        procedure, public :: interior => interior_Vector
        !
    end type Vector_t
    !
    !> Allocatable Vector element for Old Fortran polymorphic Arrays!!!
    type, public :: GenVector_t
        !
        class( Vector_t ), allocatable :: v
        !
    end type GenVector_t
    !
    !>
    abstract interface
        !
        !> No interface function briefing
        !
        function interface_get_axis_vector( self, comp_lbl ) result( comp )
            import :: Vector_t, prec
            class( Vector_t ), intent( inout ) :: self
            character, intent( in ) :: comp_lbl
            complex( kind=prec ), allocatable :: comp(:,:,:)
        end function interface_get_axis_vector
        !
        function interface_diag_mult_vector( self, rhs ) result( diag_mult )
            import :: Vector_t, Field_t
            class( Vector_t ), intent( inout ) :: self
            class( Field_t ), intent( in ) :: rhs
            class( Vector_t ), allocatable :: diag_mult
        end function interface_diag_mult_vector
        !
        !> No interface subroutine briefing
        !
        subroutine interface_interp_func_vector( self, location, xyz, interp )
            import :: Vector_t, prec
            class( Vector_t ), intent( in ) :: self
            real( kind=prec ), intent( in ) :: location(3)
            character, intent( in ) :: xyz
            class( Vector_t ), intent( inout ) :: interp
        end subroutine interface_interp_func_vector
        !
        !> No interface subroutine briefing
        subroutine interface_sum_edge_vector( self, cell_out, interior_only )
            import :: Vector_t, Scalar_t
            class( Vector_t ), intent( inout ) :: self
            class( Scalar_t ), allocatable, intent( out ) :: cell_out
            logical, intent( in ), optional :: interior_only
        end subroutine interface_sum_edge_vector
        !
        !> No interface subroutine briefing
        subroutine interface_sum_edge_vti_vector( self, cell_h_out, cell_v_out, interior_only )
            import :: Vector_t, Scalar_t
            class( Vector_t ), intent( inout ) :: self
            class( Scalar_t ), allocatable, intent( out ) :: cell_h_out, cell_v_out
            logical, optional, intent( in ) :: interior_only
        end subroutine interface_sum_edge_vti_vector
        !
        !> No interface subroutine briefing
        !
        subroutine interface_sum_cells_vector( self, cell_in, ptype )
            import :: Vector_t, Scalar_t
            class( Vector_t ), intent( inout ) :: self
            class( Scalar_t ), intent( in ) :: cell_in
            character(*), intent( in ), optional :: ptype
        end subroutine interface_sum_cells_vector
        !
        !> No interface subroutine briefing
        !
        subroutine interface_sum_cells_vti_vector( self, cell_h_in, cell_v_in, ptype )
            import :: Vector_t, Scalar_t
            class( Vector_t ), intent( inout ) :: self
            class( Scalar_t ), intent( in ) :: cell_h_in, cell_v_in
            character(*), intent( in ), optional :: ptype
        end subroutine interface_sum_cells_vti_vector
        !
        !> No interface subroutine briefing
        !
        subroutine interface_get_real_vector( self, r_vector )
            import :: Vector_t
            class( Vector_t ), intent( in ) :: self
            class( Vector_t ), allocatable, intent( out ) :: r_vector
        end subroutine interface_get_real_vector
        !
    end interface
    !
contains
    !
    !> No subroutine briefing
    !
    subroutine boundary_Vector( self, boundary )
        implicit none
        !
        class( Vector_t ), intent( in ) :: self
        class( Vector_t ), allocatable, intent( inout ) :: boundary
        !
        complex( kind=prec ), allocatable, dimension(:) :: c_array
        !
        allocate( boundary, source = self )
        !
        c_array = boundary%getArray()
        !
        if( self%grid_type == EDGE ) then
            !
            c_array( self%grid%EDGEi ) = C_ZERO
            !
        elseif( self%grid_type == FACE ) then
            !
            c_array( self%grid%FACEi ) = C_ZERO
            !
        else
            call errStop( "boundary_Vector > unrecognized grid type: ["//self%grid_type//"]" )
        endif
        !
        call boundary%setArray( c_array )
        !
    end subroutine boundary_Vector
    !
    !> No subroutine briefing
    !
    subroutine interior_Vector( self, interior )
        implicit none
        !
        class( Vector_t ), intent( in ) :: self
        class( Vector_t ), allocatable, intent( inout ) :: interior
        !
        complex( kind=prec ), allocatable, dimension(:) :: c_array
        !
        allocate( interior, source = self )
        !
        c_array = interior%getArray()
        !
        if( self%grid_type == EDGE ) then
            !
            c_array( self%grid%EDGEb ) = C_ZERO
            !
        elseif( self%grid_type == FACE ) then
            !
            c_array( self%grid%FACEb ) = C_ZERO
            !
        else
            call errStop( "boundary_Vector > unrecognized grid type: ["//self%grid_type//"]" )
        endif
        !
        call interior%setArray( c_array )
        !
    end subroutine interior_Vector
    !
end module Vector
!