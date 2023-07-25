!
!> Abstract class to define an abstract Vector field
!
module Vector
    !
    use Scalar
    !
    type, abstract, extends( Field_t ) :: Vector_t
        !
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
        procedure( interface_avg_cells_vector ), deferred, public :: avgCell
        procedure( interface_avg_cells_VTI_vector ), deferred, public :: avgCellVTI
        generic :: avgCells => avgCell, avgCellVTI
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
        procedure, public :: boundary => boundary_Vector
        procedure, public :: interior => interior_Vector
        !
        procedure, public :: switchStoreState => switchStoreState_Vector
        !
    end type Vector_t
    !
    !>
    abstract interface
        !
        !> No interface function briefing
        !
        function interface_get_axis_vector( self, comp_lbl ) result( comp )
            import :: Vector_t, prec
            character, intent( in ) :: comp_lbl
            class( Vector_t ), intent( in ) :: self
            complex( kind=prec ), allocatable :: comp(:, :, :)
        end function interface_get_axis_vector
        !
        function interface_diag_mult_vector( self, rhs ) result( diag_mult )
            import :: Vector_t
            class( Vector_t ), intent( inout ) :: self
            class( Vector_t ), intent( in ) :: rhs
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
            class( Scalar_t ), intent( inout ) :: cell_out
            logical, optional, intent( in ) :: interior_only
        end subroutine interface_sum_edge_vector
        !
        !> No interface subroutine briefing
        subroutine interface_sum_edge_vti_vector( self, cell_h_out, cell_v_out, interior_only )
            import :: Vector_t, Scalar_t
            class( Vector_t ), intent( inout ) :: self
            class( Scalar_t ), intent( inout ) :: cell_h_out, cell_v_out
            logical, optional, intent( in ) :: interior_only
        end subroutine interface_sum_edge_vti_vector
        !
        !> No interface subroutine briefing
        !
        subroutine interface_avg_cells_vector( self, cell_in, ptype )
            import :: Vector_t, Scalar_t
            class( Vector_t ), intent( inout ) :: self
            class( Scalar_t ), intent( in ) :: cell_in
            character(*), intent( in ), optional :: ptype
        end subroutine interface_avg_cells_vector
        !
        !> No interface subroutine briefing
        !
        subroutine interface_avg_cells_vti_vector( self, cell_h_in, cell_v_in, ptype )
            import :: Vector_t, Scalar_t
            class( Vector_t ), intent( inout ) :: self
            class( Scalar_t ), intent( inout ) :: cell_h_in, cell_v_in
            character(*), intent( in ), optional :: ptype
        end subroutine interface_avg_cells_vti_vector
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
            class( Vector_t ), intent( inout ) :: self
            complex( kind=prec ), allocatable :: x(:, :, :)
        end function interface_get_x_vector
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_x_vector( self, x )
            import :: Vector_t, prec
            class( Vector_t ), intent( inout ) :: self
            complex( kind=prec ), allocatable, intent( in ) :: x(:, :, :)
        end subroutine interface_set_x_vector
        !
        !> No interface function briefing
        !
        function interface_get_y_vector( self ) result( y )
            import :: Vector_t, prec
            class( Vector_t ), intent( inout ) :: self
            complex( kind=prec ), allocatable :: y(:, :, :)
        end function interface_get_y_vector
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_y_vector( self, y )
            import :: Vector_t, prec
            class( Vector_t ), intent( inout ) :: self
            complex( kind=prec ), allocatable, intent( in ) :: y(:, :, :)
        end subroutine interface_set_y_vector
        !
        !> No interface function briefing
        !
        function interface_get_z_vector( self ) result( z )
            import :: Vector_t, prec
            class( Vector_t ), intent( inout ) :: self
            complex( kind=prec ), allocatable :: z(:, :, :)
        end function interface_get_z_vector
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_z_vector( self, z )
            import :: Vector_t, prec
            class( Vector_t ), intent( inout ) :: self
            complex( kind=prec ), allocatable, intent( in ) :: z(:, :, :)
        end subroutine interface_set_z_vector
        !
    end interface
    !
contains
    !
    !> No subroutine briefing
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
        c_array( self%ind_boundaries ) = C_ZERO
        !
        call interior%setArray( c_array )
        !
    end subroutine interior_Vector
    !
    !> No subroutine briefing
    !
    subroutine switchStoreState_Vector( self, store_state )
        implicit none
        !
        class( Vector_t ), intent( inout ) :: self
        integer, intent( in ), optional :: store_state
        !
        integer i1, i2
        complex( kind=prec ), allocatable, dimension(:,:,:) :: x, y, z
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
                    s_v = (/reshape( self%getX(), (/self%Nxyz(1), 1/) ), &
                            reshape( self%getY(), (/self%Nxyz(2), 1/) ), &
                            reshape( self%getZ(), (/self%Nxyz(3), 1/) )/)
                    !
                    call self%setSV( s_v )
                    !
                case( singleton )
                    !
                    if( self%grid_type == EDGE ) then
                        !
                        allocate( x( self%nx, self%ny + 1, self%nz + 1 ) )
                        allocate( y( self%nx + 1, self%ny, self%nz + 1 ) )
                        allocate( z( self%nx + 1, self%ny + 1, self%nz ) )
                        !
                    else if( self%grid_type == FACE ) then
                        !
                        allocate( x( self%nx + 1, self%ny, self%nz ) )
                        allocate( y( self%nx, self%ny + 1, self%nz ) )
                        allocate( z( self%nx, self%ny, self%nz + 1 ) )
                        !
                    else
                        call errStop( "switchStoreState_Vector > Only EDGE or FACE types allowed." )
                    endif
                    !
                    s_v = self%getSV()
                    !
                    ! Ex
                    i1 = 1; i2 = self%Nxyz(1)
                    x = reshape( s_v( i1 : i2 ), self%NdX )
                    !
                    call self%setX( x )
                    !
                    ! Ey
                    i1 = i2 + 1; i2 = i2 + self%Nxyz(2)
                    y = reshape( s_v( i1 : i2 ), self%NdY )
                    !
                    call self%setY( y )
                    !
                    ! Ez
                    i1 = i2 + 1; i2 = i2 + self%Nxyz(3)
                    z = reshape( s_v( i1 : i2 ), self%NdZ )
                    !
                    call self%setZ( z )
                    !
                case default
                    call errStop( "switchStoreState_Vector > Unknown store_state" )
                    !
            end select
            !
        endif
        !
    end subroutine switchStoreState_Vector
    !
end module Vector
