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
        logical :: is_allocated
        !
        integer, private :: counter
        !
        contains
            !
            procedure, public :: put => putValuesDataGroup
            !
            procedure, public :: add => addValuesDataGroup
            !
            procedure, public :: sub => subValuesDataGroup
            !
            procedure, public :: multValuesDataGroup
            procedure, public :: multDataDataGroup
            generic :: mult => multValuesDataGroup, multDataDataGroup
            !
            procedure, public :: setRealValuesDataGroup
            procedure, public :: setComplexValuesDataGroup
            generic :: set => setRealValuesDataGroup, setComplexValuesDataGroup
            !
            procedure, public :: normalize => normalizeDataGroup
            !
            procedure, public :: reset => resetDataGroup
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
        self%errors = R_ONE
        !
        self%is_allocated = .TRUE.
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
    !> Error normalized subtraction with another DataGroup.
    subroutine addValuesDataGroup( self, rhs )
        implicit none
        !
        class( DataGroup_t ), intent( inout ) :: self
        class( DataGroup_t ), intent( in ) :: rhs
        !
        integer :: i
        complex( kind=prec ) :: sum_value, self_value, rhs_value
        !
        do i = 1, self%n_comp
            !
            self%errors(i) = ( self%errors(i) + rhs%errors(i) ) / 2.0
            !
            self_value = cmplx( self%reals(i), self%imaginaries(i), kind=prec )
            !
            rhs_value = cmplx( rhs%reals(i), rhs%imaginaries(i), kind=prec )
            !
            sum_value = cmplx( self_value + rhs_value, kind=prec )
            !
            self%reals(i) = real( sum_value, kind=prec )
            !
            self%imaginaries(i) = real( aimag( sum_value ), kind=prec )
            !
        enddo
        !
    end subroutine addValuesDataGroup
    !
    !> Error normalized subtraction with another DataGroup.
    subroutine subValuesDataGroup( self, rhs )
        implicit none
        !
        class( DataGroup_t ), intent( inout ) :: self
        class( DataGroup_t ), intent( in ) :: rhs
        !
        integer :: i
        complex( kind=prec ) :: sum_value, self_value, rhs_value
        !
        do i = 1, self%n_comp
            !
            self%errors(i) = ( self%errors(i) + rhs%errors(i) ) / 2.0
            !
            self_value = cmplx( self%reals(i), self%imaginaries(i), kind=prec )
            !
            rhs_value = cmplx( rhs%reals(i), rhs%imaginaries(i), kind=prec )
            !
            sum_value = cmplx( self_value - rhs_value, kind=prec )
            !
            self%reals(i) = real( sum_value, kind=prec )
            !
            self%imaginaries(i) = real( aimag( sum_value ), kind=prec )
            !
        enddo
        !
    end subroutine subValuesDataGroup
    !
    !> Multiplication by another DataGroup.
    subroutine multValuesDataGroup( self, rvalue )
        implicit none
        !
        class( DataGroup_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: rvalue
        !
        integer :: i
        complex( kind=prec ) :: mult_value, self_value
        !
        do i = 1, self%n_comp
            !
            self_value = cmplx( self%reals(i), self%imaginaries(i), kind=prec )
            !
            mult_value = cmplx( self_value * rvalue, kind=prec )
            !
            self%reals(i) = real( mult_value, kind=prec )
            !
            self%imaginaries(i) = real( aimag( mult_value ), kind=prec )
            !
        enddo
        !
    end subroutine multValuesDataGroup
    !
    !> Multiplication by another DataGroup.
    subroutine multDataDataGroup( self, rhs )
        implicit none
        !
        class( DataGroup_t ), intent( inout ) :: self
        class( DataGroup_t ), intent( in ) :: rhs
        !
        integer :: i
        complex( kind=prec ) :: mult_value, self_value, rhs_value
        !
        do i = 1, self%n_comp
            !
            self%errors(i) = ( self%errors(i) + rhs%errors(i) ) / 2.0
            !
            self_value = cmplx( self%reals(i), self%imaginaries(i), kind=prec )
            !
            rhs_value = cmplx( rhs%reals(i), rhs%imaginaries(i), kind=prec )
            !
            mult_value = cmplx( self_value * rhs_value, kind=prec )
            !
            self%reals(i) = real( mult_value, kind=prec )
            !
            self%imaginaries(i) = real( aimag( mult_value ), kind=prec )
            !
        enddo
        !
    end subroutine multDataDataGroup
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
    end subroutine setComplexValuesDataGroup
    !
    !> Error normalized subtraction with another DataGroup.
    subroutine normalizeDataGroup( self, rhs )
        implicit none
        !
        class( DataGroup_t ), intent( inout ) :: self
        class( DataGroup_t ), intent( in ) :: rhs
        !
        integer :: i
        real( kind=prec ) :: self_square_error
        complex( kind=prec ) :: self_value, rhs_value, normalized_value
        !
        do i = 1, self%n_comp
            !
            if( self%errors(i) == R_ZERO ) then
                stop "Error: normalizeDataGroup > Zero Error"
            else
                !
                self_square_error = self%errors(i) * self%errors(i)
                !
                !write( *, * ) "self_square_error: ", self_square_error
                !
                self_value = cmplx( self%reals(i), self%imaginaries(i), kind=prec )
                !
                rhs_value = cmplx( rhs%reals(i), rhs%imaginaries(i), kind=prec )
                !
                !> ( Measured - Predicted ) / Measured_Error**2 
                normalized_value = ( self_value - rhs_value ) / self_square_error
                !
                !> JUST FOR CHECK LATER IN DCG ????
                !normalized_value = self_value + rhs_value / ( ( self_value - rhs_value ) * ( self_value - rhs_value ) )
                !
                self%reals(i) = real( normalized_value, kind=prec )
                !
                self%imaginaries(i) = real( aimag( normalized_value ), kind=prec )
                !
            endif
        enddo
        !
    end subroutine normalizeDataGroup
    !
    !> No function briefing
    function dotProdDataGroup( self, rhs ) result( cvalue )
        implicit none
        !
        class( DataGroup_t ), intent( in ) :: self, rhs
        !
        complex( kind=prec ) :: cvalue
        !
        integer :: i
        complex( kind=prec ) :: self_value, rhs_value
        !
        cvalue = C_ZERO
        !
        do i = 1, self%n_comp
            !
            self_value = conjg( cmplx( self%reals(i), self%imaginaries(i), kind=prec ) )
            !
            rhs_value = cmplx( rhs%reals(i), rhs%imaginaries(i), kind=prec )
            !
            cvalue = cvalue + ( self_value * rhs_value )
            !
        enddo
        !
    end function dotProdDataGroup
    !
    !> ????
    subroutine resetDataGroup( self )
        implicit none
        !
        class( DataGroup_t ), intent( inout ) :: self
        !
        self%reals = R_ZERO
        !
        self%imaginaries = R_ZERO
        !
        self%errors = R_ONE
        !
    end subroutine resetDataGroup
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
        self%is_allocated = .TRUE.
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
