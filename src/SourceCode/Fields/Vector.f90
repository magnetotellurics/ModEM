!
!> Abstract class to define an abstract Vector field
!
module Vector
    !
    use rScalar3D_SG
    !
    type, abstract, extends( Field_t ) :: Vector_t
        !
        integer :: counter = 0
        integer, dimension(3) :: NdX, NdY, NdZ, Nxyz
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
        procedure( interface_get_x_vector ), deferred, public :: getX
        procedure( interface_set_x_vector ), deferred, public :: setX
        procedure( interface_get_y_vector ), deferred, public :: getY
        procedure( interface_set_y_vector ), deferred, public :: setY
        procedure( interface_get_z_vector ), deferred, public :: getZ
        procedure( interface_set_z_vector ), deferred, public :: setZ
        !
        procedure, public :: length => length_Vector
        !
        procedure, public :: boundary => boundary_Vector
        procedure, public :: interior => interior_Vector
        !
        procedure, public :: getArray => getArray_Vector
        procedure, public :: setArray => setArray_Vector
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
            class( Vector_t ), allocatable, intent( inout ) :: interp
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
        subroutine interface_get_real_vector( self, r_vector )
            import :: Vector_t
            class( Vector_t ), intent( in ) :: self
            class( Vector_t ), allocatable, intent( out ) :: r_vector
        end subroutine interface_get_real_vector
        !
        !> No interface function briefing
        !
        function interface_get_x_vector( self ) result( x )
            import :: Vector_t, prec
            !
            class( Vector_t ), intent( in ) :: self
            !
            complex( kind=prec ), allocatable, dimension(:,:,:) :: x
            !
        end function interface_get_x_vector
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_x_vector( self, x )
            import :: Vector_t, prec
            !
            class( Vector_t ), intent( inout ) :: self
            complex( kind=prec ), dimension(:,:,:), intent( in ) :: x
            !
        end subroutine interface_set_x_vector
        !
        !> No interface function briefing
        !
        function interface_get_y_vector( self ) result( y )
            import :: Vector_t, prec
            !
            class( Vector_t ), intent( in ) :: self
            !
            complex( kind=prec ), allocatable, dimension(:,:,:) :: y
            !
        end function interface_get_y_vector
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_y_vector( self, y )
            import :: Vector_t, prec
            !
            class( Vector_t ), intent( inout ) :: self
            complex( kind=prec ), dimension(:,:,:), intent( in ) :: y
            !
        end subroutine interface_set_y_vector
        !
        !> No interface function briefing
        !
        function interface_get_z_vector( self ) result( z )
            import :: Vector_t, prec
            !
            class( Vector_t ), intent( in ) :: self
            !
            complex( kind=prec ), allocatable, dimension(:,:,:) :: z
            !
        end function interface_get_z_vector
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_z_vector( self, z )
            import :: Vector_t, prec
            !
            class( Vector_t ), intent( inout ) :: self
            complex( kind=prec ), dimension(:,:,:), intent( in ) :: z
            !
        end subroutine interface_set_z_vector
        !
    end interface
    !
contains
    !
    !> No subroutine briefing
    !
    function length_Vector( self ) result( n )
        implicit none
        !
        class( Vector_t ), intent( in ) :: self
        !
        integer :: n
        !
        n = self%Nxyz(1) + self%Nxyz(2) + self%Nxyz(3)
        !
    end function length_Vector
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
        c_array( self%ind_interior ) = C_ZERO
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
        c_array( self%ind_boundary ) = C_ZERO
        !
        call interior%setArray( c_array )
        !
    end subroutine interior_Vector
    !
    !> No subroutine briefing
    !
    function getArray_Vector( self ) result( array )
        implicit none
        !
        class( Vector_t ), intent( in ) :: self
        !
        complex( kind=prec ), allocatable, dimension(:) :: array
        !
        if( ( .NOT. self%is_allocated ) ) then
            call errStop( "getArray_Vector > Self not allocated." )
        endif
        !
        if( self%store_state .EQ. compound ) then
            ! !
            ! write( *, * ) "getArray_Vector: ", &
            ! self%Nxyz(1), self%Nxyz(2), self%Nxyz(3), &
            ! self%Nxyz(1)+self%Nxyz(2)+self%Nxyz(3), self%length()
            ! !
            allocate( array( self%length() ) )
            !
            array = (/reshape(self%getX(), (/self%Nxyz(1), 1/)), &
            reshape(self%getY(), (/self%Nxyz(2), 1/)), &
            reshape(self%getZ(), (/self%Nxyz(3), 1/))/)
            !
        elseif( self%store_state .EQ. singleton ) then
            !
            array = self%getSV()
            !
        else
            call errStop( "getArray_Vector > Unknown store_state!" )
        endif
        !
    end function getArray_Vector
    !
    !> No subroutine briefing
    !
    subroutine setArray_Vector( self, array )
        implicit none
        !
        class( Vector_t ), intent( inout ) :: self
        complex( kind=prec ), dimension(:), intent( in ) :: array
        !
        complex( kind=prec ), allocatable, dimension(:,:,:) :: v
        integer :: i1, i2
        !
        if( ( .NOT. self%is_allocated ) ) then
            call errStop( "setArray_Vector > Self not allocated." )
        endif
        !
        call self%deallOtherState
        !
        if( self%store_state .EQ. compound ) then
            !
            !> Ex
            i1 = 1; i2 = self%Nxyz(1)
            !
            v = reshape( array(i1:i2), self%NdX )
            call self%setX( v )
            !
            !> Ey
            i1 = i2 + 1; i2 = i2 + self%Nxyz(2)
            !
            v = reshape( array(i1:i2), self%NdY )
            call self%setY( v )
            !
            !> Ez
            i1 = i2 + 1; i2 = i2 + self%Nxyz(3)
            !
            v = reshape(array(i1:i2), self%NdZ)
            call self%setZ( v )
            !
        elseif( self%store_state .EQ. singleton ) then
            !
            call self%setSV( array )
            !
        else
            call errStop( "setArray_Vector > Unknown store_state!" )
        endif
        !
    end subroutine setArray_Vector
    !
end module Vector
