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
        procedure( interface_diag_mult_vector ), deferred, public :: diagMult
        procedure( interface_interp_func_vector ), deferred, public :: interpFunc
        !
        procedure( interface_avg_cells_vector ), deferred, public :: avgCell
        procedure( interface_avg_cells_VTI_vector ), deferred, public :: avgCellVTI
        generic :: avgCells => avgCell, avgCellVTI
        !
        procedure( interface_get_real_vector ), deferred, public :: getReal
        !
    end type Vector_t
    !
    !>
    abstract interface
        !
        !> No interface function briefing
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
        subroutine interface_avg_cells_vti_vector( self, cell_in, ptype )
            import :: Vector_t, Scalar_t
            class( Vector_t ), intent( inout ) :: self
            class( Scalar_t ), allocatable, dimension(:), intent( in ) :: cell_in
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
    end interface
    !
end module Vector
