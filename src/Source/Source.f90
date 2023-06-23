!
!> Abstract Base class to define a Source
!
module Source
    !
    use Constants
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
        class( Vector_t ), allocatable, dimension(:) :: rhs, E
        !
        logical :: non_zero_source, non_zero_bc, calc_sens, for_transpose
        !
        contains
            !
            procedure, public :: init => initializeSource
            !
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
    public :: minNode, maxNode, clean
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
    function minNode( x, xNode ) result( ix )
        implicit none
        !
        real( kind=prec ), intent( in ) :: x
        real( kind=prec ), dimension(:), intent( in ) :: xNode
        !
        integer :: ix, i
        !
        do i = 1, size( xNode )
            if( clean( xNode(i) ) .GT. clean(x) ) then
                ix = i-1
                exit
            endif
        enddo
        !
    end function minNode
    !
    !>    This is a utility routine, used by several data functional
    !>    set up routines, and for other interpolation functions
    !>    Returns index ix such that    xNode(ix) <= x < xNode(ix+1)
    !>    If x is out of range:
    !>    x > xNode(1) returns 0; if x< xNode(nx) returns nx
    !>    Assumes xNode is strictly decreasing; does not check this
    !>    NOTE: as presently coded, when xNode is called with center
    !>    (face) node positions, this routine will return zero for
    !>    the coordinates in the outer half cell nearest the boundary
    !>    If evaluation over the complete model domain is to be allowed
    !>    a more general interpolation rule will be required.
    !>    A.K.: modified to allow input of any size, nx = size(xNode).
    !
    function maxNode(x, xNode) result(ix)
        implicit none
        !
        real( kind=prec ), intent( in ) :: x
        real( kind=prec ), dimension(:), intent( in ) :: xNode
        !
        integer :: ix, i
        !
        do i = 1, size(xNode)
           if( clean( xNode(i)) .LT. clean(x) ) then
                ix = i-1
                exit
           endif
        enddo
        !
    end function maxNode
    !
    !> This is a utility routine that provides an expression used to battle
    !> against machine error problems. It returns the same real or real(8)
    !> as the input, but without the extra digits at the end that are often
    !> a cause of wrong comparisons in the if statements. ALWAYS use clean(x)
    !> instead of x in an inequality!!!
    !> R_LARGE is defined in the module math_constants
    !> A.K.
    !
    function clean(x)
        implicit none
        !
        real( kind=prec ), intent( in ) :: x
        real( kind=prec ) :: clean
        !
        clean = dnint(x*R_LARGE)/R_LARGE
        !
    end function clean
    !
end module Source
