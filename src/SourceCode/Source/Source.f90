!
!> Abstract Base class to define a Source
!
module Source
    !
    use Utilities
    use ModelOperator
    use ModelParameter
    !
    character(:), allocatable :: source_type_mt
    character( len=11 ), parameter :: SRC_MT_1D = "SourceMT_1D"
    character( len=11 ), parameter :: SRC_MT_2D = "SourceMT_2D"
    !
    type, abstract :: Source_t
        !
        class( ModelOperator_t ), pointer :: model_operator
        !
        class( ModelParameter_t ), pointer :: sigma
        !
        real( kind=prec ) :: period
        !
        type( cVector3D_SG_t ), allocatable, dimension(:) :: E
        !
        type( GenVector_t ), allocatable, dimension(:) :: rhs
        !
        logical :: non_zero_source, non_zero_bc, calc_sens, for_transpose
        !
        !> Global primary electrical Field
        type( cVector3D_SG_t ) :: E_p
        !
        contains
            !
            procedure, public :: baseInit => initialize_Source
            !
            procedure, public :: baseDealloc => deallocate_Source
            !
            procedure, public :: setE => setE_Source
            !
            procedure, public :: setRhs => setRHS_Source
            !
            procedure( interface_create_e_source ), deferred, public :: createE
            !
            procedure( interface_create_rhs_source ), deferred, public :: createRHS
            !
    end type Source_t
    !
    abstract interface
        !
        subroutine interface_create_e_source( self )
            import :: Source_t
            !
            class( Source_t ), intent( inout ) :: self
        end subroutine interface_create_e_source
        !
        !> No interface subroutine briefing
        subroutine interface_create_rhs_source( self )
            import :: Source_t
            !
            class( Source_t ), intent( inout ) :: self
        end subroutine interface_create_rhs_source
        !
    end interface
    !
    contains
    !
    !> No subroutine briefing
    subroutine setE_Source( self, E )
        implicit none
        !
        class( Source_t ), intent( inout ) :: self
        type( cVector3D_SG_t ), dimension(:), intent( in ) :: E
        !
        self%E = E
        !
        call self%createRHS
        !
    end subroutine setE_Source
    !
    !> No subroutine briefing
    subroutine setRHS_Source( self, rhs )
        implicit none
        !
        class( Source_t ), intent( inout ) :: self
        type( GenVector_t ), allocatable, dimension(:) :: rhs
        !
        self%rhs = rhs
        !
    end subroutine setRHS_Source
    !
    !> No subroutine briefing
    subroutine initialize_Source( self )
        implicit none
        !
        class( Source_t ), intent( inout ) :: self
        !
        self%period = R_ZERO
        !
        self%non_zero_source = .FALSE.
        !
        self%non_zero_bc = .FALSE.
        !
        self%for_transpose = .FALSE.
        !
        self%calc_sens = .FALSE.
        !
        self%model_operator => null()
        !
        self%sigma => null()
        !
    end subroutine initialize_Source
    !
    !> No subroutine briefing
    subroutine deallocate_Source( self )
        implicit none
        !
        class( Source_t ), intent( inout ) :: self
        !
        deallocate( self%rhs )
        !
        deallocate( self%E )
        !
    end subroutine deallocate_Source
    !
end module Source
