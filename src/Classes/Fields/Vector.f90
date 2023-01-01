!
!> Abstract class to define a Vector field
!
module Vector
    !
    use Scalar
    !use cSparseVector3D_SG
    !
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
        procedure( interface_sum_edges_vector ), deferred, public :: sumEdges
        procedure( interface_sum_cells_vector ), deferred, public :: avgCells
        !
        !> Vector procedures
        procedure, public :: boundary => boundaryVector
        procedure, public :: interior => interiorVector
        !
    end type Vector_t
    !
    abstract interface
        ! !
        ! !> No interface subroutine briefing
        ! subroutine interface_add_sparse_vector_vector( self, svec )
            ! import :: Vector_t, cSparseVector3D_SG_t
            ! class( Vector_t ), intent( inout ) :: self
            ! type( cSparseVector3D_SG_t ), intent( in ) :: svec
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
        subroutine interface_sum_edges_vector( self, cell_obj, interior_only )
            import :: Vector_t, Scalar_t
            class( Vector_t ), intent( in ) :: self
            class( Scalar_t ), allocatable, intent( inout ) :: cell_obj
            logical, optional, intent( in ) :: interior_only
        end subroutine interface_sum_edges_vector
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
contains
    !
    !> No subroutine briefing
    subroutine boundaryVector( self, boundary )
        implicit none
        !
        class( Vector_t ), intent( in ) :: self
        class( Vector_t ), allocatable, intent( inout ) :: boundary
        !
        allocate( boundary, source = self )
        !
        call boundary%setAllInterior( C_ZERO )
       !
    end subroutine boundaryVector
    !
    !> No subroutine briefing
    subroutine interiorVector( self, interior )
        implicit none
        !
        class( Vector_t ), intent( in ) :: self
        class( Vector_t ), allocatable, intent( inout ) :: interior
        !
        allocate( interior, source = self )
        !
        call interior%setAllboundary( C_ZERO )
        !
    end subroutine interiorVector
    !
end module Vector
