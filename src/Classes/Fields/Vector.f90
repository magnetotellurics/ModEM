!
!> Abstract class to define a Vector field
!
module Vector
    !
    use Scalar
    !
    !>
    type :: ESolTx
        !
        class( Vector_t ), allocatable, dimension(:) :: pol
        !
    end type ESolTx
    !
    !>
    type :: ESolMTx
        !
        type( ESolTx ), allocatable, dimension(:) :: e_sol
        !
        integer :: SolnIndex = 0
        !
    end type ESolMTx
    !
    !>
    type, abstract, extends( Field_t ) :: Vector_t
        !
        integer, dimension(3) :: NdX, NdY, NdZ, Nxyz
        !
    contains
        !
        !> Vector interfaces
        !procedure( interface_add_sparse_vector_vector ), deferred, public :: addSparseVector
        !
        procedure( interface_diag_mult_vector ), deferred, public :: diagMult
        !
        procedure( interface_interp_func_vector ), deferred, public :: interpFunc
        !
        procedure( interface_sum_cells_vector ), deferred, public :: avgCells
        !
    end type Vector_t
    !
    abstract interface
        ! !
        ! !> No interface subroutine briefing
        ! subroutine interface_add_sparse_vector_vector( self, svec )
            ! import :: Vector_t, cVectorSparse3D_SG_t
            ! class( Vector_t ), intent( inout ) :: self
            ! type( cVectorSparse3D_SG_t ), intent( in ) :: svec
            ! !
        ! end subroutine interface_add_sparse_vector_vector
        ! !
        !> Miscellaneous
        function interface_diag_mult_vector( self, rhs ) result( diag_mult )
            import :: Vector_t
            class( Vector_t ), intent( in ) :: self, rhs
            class( Vector_t ), allocatable :: diag_mult
        end function interface_diag_mult_vector
        !
        !> No interface subroutine briefing
        subroutine interface_interp_func_vector( self, location, xyz, interp )
            import :: Vector_t, prec
            class( Vector_t ), intent( in ) :: self
            real( kind=prec ), intent( in ) :: location(3)
            character, intent( in ) :: xyz
            class( Vector_t ), allocatable, intent( inout ) :: interp
        end subroutine interface_interp_func_vector
        !
        !> No interface subroutine briefing
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
