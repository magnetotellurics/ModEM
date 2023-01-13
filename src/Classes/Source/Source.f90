!
!> Abstract Base class to define a Source
!
module Source
    !
    use ModelOperator
    use ModelParameter
    !
    character(:), allocatable :: source_type
    character( len=11 ), parameter :: SRC_MT_1D = "SourceMT_1D"
    character( len=11 ), parameter :: SRC_MT_2D = "SourceMT_2D"
    !
    character(:), allocatable :: get_1D_from
    !
    type, abstract :: Source_t
        !
        class( ModelOperator_t ), pointer :: model_operator
        !
        class( ModelParameter_t ), pointer :: sigma
        !
        real( kind=prec ) :: period
        !
        class( Vector_t ), allocatable, dimension(:) :: rhs, E
        !
        logical :: non_zero_source, sens, trans
        !
        contains
            !
            procedure, public :: init => initializeSource
            procedure, public :: dealloc => deallocateSource
            !
            procedure, public :: setE => setESource
            !
            procedure, public :: setRhs => setRhsSource
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
            !
            import :: Source_t
            class( Source_t ), intent( inout ) :: self
            !
        end subroutine interface_create_e_source
        !
        !> No interface subroutine briefing
        subroutine interface_create_rhs_source( self )
            !
            import :: Source_t
            class( Source_t ), intent( inout ) :: self
            !
        end subroutine interface_create_rhs_source
        !
    end interface
    !
    contains
    !
    !> No subroutine briefing
    subroutine setESource( self, E )
        implicit none
        !
        class( Source_t ), intent( inout ) :: self
        class( Vector_t ), dimension(:), intent( in ) :: E
        !
        integer :: pol
        !
        if( allocated( self%E ) ) deallocate( self%E )
        !
        allocate( self%E, source = E )
        !
        call self%createRHS()
        !
    end subroutine setESource
    !
    !> No subroutine briefing
    subroutine setRhsSource( self, rhs )
        implicit none
        !
        class( Source_t ), intent( inout ) :: self
        class( Vector_t ), dimension(:), intent( in ) :: rhs
        !
        integer :: pol
        !
        if( allocated( self%rhs ) ) deallocate( self%rhs )
        !
        allocate( self%rhs, source = rhs )
        !
    end subroutine setRhsSource
    !
    !> No subroutine briefing
    subroutine initializeSource( self )
        implicit none
        !
        class( Source_t ), intent( inout ) :: self
        !
        self%period = R_ZERO
        !
        self%non_zero_source = .FALSE.
        !
        self%trans = .FALSE.
        !
        self%sens = .FALSE.
        !
        self%model_operator => null()
        !
        self%sigma => null()
        !
    end subroutine initializeSource
    !
    !> No subroutine briefing
    subroutine deallocateSource( self )
        implicit none
        !
        class( Source_t ), intent( inout ) :: self
        !
        deallocate( self%rhs )
        !
        deallocate( self%E )
        !
    end subroutine deallocateSource
    !
end module Source
