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
        procedure( interface_sum_cells_vector ), deferred, public :: avgCells
        procedure( interface_diag_mult_vector ), deferred, public :: diagMult
        procedure( interface_interp_func_vector ), deferred, public :: interpFunc
        !
    end type Vector_t
    !
    !> Global structures to store E from Transmitters
    !
    !> For a single Tx
    type :: EAllTx_t
        !
        class( Vector_t ), allocatable, dimension(:) :: pol
        !
    end type EAllTx_t
    !
    !> For multiple Txs
    type :: EAllMTx_t
        !
        type( EAllTx_t ), allocatable, dimension(:) :: e
        !
        integer :: SolnIndex = 0
        !
    end type EAllMTx_t
    !
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
        subroutine interface_sum_cells_vector( self, E_in, ptype )
            import :: Vector_t, Scalar_t
            class( Vector_t ), intent( inout ) :: self
            class( Scalar_t ), intent( in ) :: E_in
            character(*), intent( in ), optional :: ptype
        end subroutine interface_sum_cells_vector
        !
    end interface
    !
end module Vector
