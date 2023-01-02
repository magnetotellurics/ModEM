!
!> Lone class to define a group of n_comp entries, related to a specific transmitter-receiver pair.
!
module DataGroup
    !
    use String
    !
    use Constants
    !
    !> Global file path name for data files
    character(:), allocatable :: predicted_data_file_name, JmHat_data_file_name
    !
    type :: DataGroup_t
        !
        integer :: n_comp, i_rx, i_tx
        !
        type( String_t ), allocatable, dimension(:) :: components
        !
        real( kind=prec ), allocatable, dimension(:) :: reals, imaginaries, errors
        !
        logical :: is_allocated, error_bar
        !
        integer, private :: counter
        !
        contains
            !
            procedure, public :: put => putValuesDataGroup
            !
            procedure, public :: setRealValuesDataGroup
            procedure, public :: setComplexValuesDataGroup
            generic :: set => setRealValuesDataGroup, setComplexValuesDataGroup
            !
            procedure, public :: linComb => linCombDataGroup
            !
            procedure, public :: dotProd => dotProdDataGroup
            !
            procedure, public :: zeros => zerosDataGroup
            !
            procedure, public :: isEqual => isEqualDg
            !
            procedure, public :: copyFrom => copyFromDataGroup
            generic :: assignment(=) => copyFrom
            !
            procedure, public :: print => printDataGroup
            !
    end type DataGroup_t
    !
    interface DataGroup_t
        module procedure DataGroup_ctor
    end interface DataGroup_t
    !
contains
    !
    !> Parametrized constructor:
    !> Set all variables and deallocate/allocate and initialize all arrays with the same n_comp size.
    function DataGroup_ctor( i_rx, i_tx, n_comp ) result( self )
        implicit none
        !
        type( DataGroup_t ) :: self
        !
        integer, intent( in ) :: i_rx, i_tx, n_comp
        !
        integer :: i, asize
        !
        self%i_rx = i_rx
        self%i_tx = i_tx
        self%n_comp = n_comp
        !
        self%counter = 1
        !
        if( allocated( self%components ) ) then
            asize = size( self%components )
            do i = asize, 1, -(1)
                if( allocated( self%components(i)%str ) ) deallocate( self%components(i)%str )
            enddo
            deallocate( self%components )
        endif
        !
        allocate( String_t :: self%components( n_comp ) )
        !
        if( allocated( self%reals ) ) deallocate( self%reals )
        allocate( self%reals( n_comp ) )
        !
        self%reals = R_ZERO
        !
        if( allocated( self%imaginaries ) ) deallocate( self%imaginaries )
        allocate( self%imaginaries( n_comp ) )
        !
        self%imaginaries = R_ZERO
        !
        if( allocated( self%errors ) ) deallocate( self%errors )
        allocate( self%errors( n_comp ) )
        !
        self%errors = R_ZERO
        !
        self%is_allocated = .TRUE.
        !
        self%error_bar = .FALSE.
        !
    end function DataGroup_ctor
    !
    !> Add values to arrays in position and increments the internal counter.
    subroutine putValuesDataGroup( self, component, rvalue, imaginary, error )
        implicit none
        !
        class( DataGroup_t ), intent( inout ) :: self
        character(:), allocatable, intent( in ) :: component
        real( kind=prec ), intent( in ) :: rvalue, imaginary, error
        !
        self%components( self%counter )%str = component
        !
        self%reals( self%counter ) = rvalue
        !
        self%imaginaries( self%counter ) = imaginary
        !
        self%errors( self%counter ) = error
        !
        self%counter = self%counter + 1
        !
    end subroutine putValuesDataGroup
    !
    !> Set the value (real and imaginary only) at a given index of these arrays.
    subroutine setRealValuesDataGroup( self, data_id, rvalue, imaginary )
        implicit none
        !
        class( DataGroup_t ), intent( inout ) :: self
        integer, intent( in ) :: data_id
        real( kind=prec ), intent( in ) :: rvalue, imaginary
        !
        self%reals( data_id ) = rvalue
        !
        self%imaginaries( data_id ) = imaginary
        !
        self%errors( data_id ) = R_ZERO
        !
        self%error_bar = .FALSE.
        !
    end subroutine setRealValuesDataGroup
    !
    !> Set the value (real and imaginary only) at a given index of these arrays.
    subroutine setComplexValuesDataGroup( self, data_id, cvalue )
        implicit none
        !
        class( DataGroup_t ), intent( inout ) :: self
        integer, intent( in ) :: data_id
        complex( kind=prec ), intent( in ) :: cvalue
        !
        self%reals( data_id ) = real( cvalue, kind=prec )
        !
        self%imaginaries( data_id ) = real( aimag( cvalue ), kind=prec )
        !
        self%errors( data_id ) = R_ZERO
        !
        self%error_bar = .FALSE.
        !
    end subroutine setComplexValuesDataGroup
    !
    !> ????
    subroutine linCombDataGroup( self, a, b, d_group, d_group_out )
        implicit none
        !
        class( DataGroup_t ), intent( in ) :: self, d_group
        real( kind=prec ), intent( in ) :: a, b
        class( DataGroup_t ), intent( inout ) :: d_group_out
        !
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            stop "Error: linCombDataGroup > self not allocated"
        endif
        !
        if( .NOT. d_group%is_allocated ) then
            stop "Error: linCombDataGroup > d_group not allocated"
        endif
        !
        if( .NOT. d_group_out%is_allocated ) then
            stop "Error: linCombDataGroup > d_group_out not allocated"
        endif
        !
        d_group_out%error_bar = self%error_bar .OR. d_group%error_bar
        d_group_out%i_rx = self%i_rx
        d_group_out%i_tx = self%i_tx
        !
        do i = 1, self%n_comp
            !
            d_group_out%reals(i) = a * self%reals(i) + b * d_group%reals(i)
            !
            if( self%error_bar .AND. d_group%error_bar ) then
                if( abs(a) > R_ZERO .AND. abs(b) > R_ZERO ) then
                    stop "Error: linCombDataGroup: unable to add two data vectors with error bars"
                else if( abs(a) > R_ZERO ) then
                    d_group_out%errors(i) = a * self%errors(i)
                    !dOut%normalized = d1%normalized
                else if( abs(b) > R_ZERO ) then
                    d_group_out%errors(i) = b * d_group%errors(i)
                    !dOut%normalized = d2%normalized
                end if
            else if( self%error_bar ) then
                d_group_out%errors(i) = a * self%errors(i)
                !dOut%normalized = d1%normalized
            else if( d_group%error_bar ) then
                d_group_out%errors(i) = b * d_group%errors(i)
                !dOut%normalized = d2%normalized
            end if
            !
        enddo
        !
    end subroutine linCombDataGroup
    !
    !> ????
    function dotProdDataGroup( self, d_group ) result( rvalue )
        implicit none
        !
        class( DataGroup_t ), intent( in ) :: self, d_group
        !
        real( kind=prec ) :: rvalue
        !
        integer :: i
        !
        rvalue = 0.0
        !
        do i = 1, self%n_comp
            !
            rvalue = rvalue + self%reals(i) * d_group%reals(i)
            !
        enddo
        !
    end function dotProdDataGroup
    !
    !> ????
    subroutine zerosDataGroup( self )
        implicit none
        !
        class( DataGroup_t ), intent( inout ) :: self
        !
        self%reals = R_ZERO
        !
        self%imaginaries = R_ZERO
        !
        self%errors = R_ZERO
        !
        self%error_bar = .FALSE.
        !
    end subroutine zerosDataGroup
    !
    !> Return if it is similar to another DataGroup.
    function isEqualDg( self, other ) result ( equal )
        implicit none
        !
        class( DataGroup_t ), intent( in ) :: self
        type( DataGroup_t ), intent( in ) :: other
        logical :: equal
        !
        integer :: i_comp
        !
        equal = .TRUE.
        !
        if( self%i_tx /= other%i_tx .OR. self%i_rx /= other%i_rx ) then
            equal = .FALSE.
            return
        endif
        !
        do i_comp = 1, self%n_comp
            !
            if( self%components( i_comp )%str /= other%components( i_comp )%str .OR. &
                ABS( self%reals( i_comp ) - other%reals( i_comp ) ) >= TOL6 .OR. &
                ABS( self%imaginaries( i_comp ) - other%imaginaries( i_comp ) ) >= TOL6 .OR. &
                ABS( self%errors( i_comp ) - other%errors( i_comp ) ) >= TOL6 ) then
                    equal = .FALSE.
                    exit
            endif
            !
        enddo
        !
    end function isEqualDg
    !
    !> Copy/Assign (= operator) all content from another DataGroup.
    subroutine copyFromDataGroup( self, rhs )
        implicit none
        !
        class( DataGroup_t ), intent( inout ) :: self
        type( DataGroup_t ), intent( in ) :: rhs
        !
        integer :: i, asize
        !
        if( .NOT. rhs%is_allocated ) then
            stop "Error: copyFromDataGroup > rhs not allocated"
        endif
        !
        self%n_comp = rhs%n_comp
        self%i_rx = rhs%i_rx
        self%i_tx = rhs%i_tx
        self%counter = rhs%counter
        self%error_bar = rhs%error_bar
        !
        if( allocated( self%components ) ) then
            asize = size( self%components )
            do i = asize, 1, -(1)
                if( allocated( self%components(i)%str ) ) deallocate( self%components(i)%str )
            enddo
            deallocate( self%components )
        endif
        !
        allocate( self%components, source = rhs%components )
        !
        if( allocated( self%reals ) ) deallocate( self%reals )
        allocate( self%reals, source = rhs%reals )
        !
        if( allocated( self%imaginaries ) ) deallocate( self%imaginaries )
        allocate( self%imaginaries, source = rhs%imaginaries )
        !
        if( allocated( self%errors ) ) deallocate( self%errors )
        allocate( self%errors, source = rhs%errors )
        !
        self%is_allocated = rhs%is_allocated
        !
    end subroutine copyFromDataGroup
    !
    !> Print the entire contents of the DataGroup on the screen.
    subroutine printDataGroup( self )
        implicit none
        !
        class( DataGroup_t ), intent( in ) :: self
        !
        integer :: i_comp
        !
        write( *, * ) "    Write DataGroup_t Id"
        write( *, * ) "             Receiver Id: ", self%i_rx
        write( *, * ) "          Transmitter Id: ", self%i_tx
        write( *, * ) self%n_comp, " data_rows:"
        !
        do i_comp = 1, self%n_comp
            !
            write( *, * ) i_comp, ":", self%components( i_comp )%str, self%reals( i_comp ), self%imaginaries( i_comp ), self%errors( i_comp )
            !
        enddo
        !
    end subroutine printDataGroup
    !
end module DataGroup
