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
        integer :: n_comp, i_rx, i_tx, normalized
        !
        type( String_t ), allocatable, dimension(:) :: components
        !
        real( kind=prec ), allocatable, dimension(:) :: reals, imaginaries, errors
        !
        logical :: is_allocated, is_complex, error_bar
        !
        integer, private :: counter
        !
        contains
            !
            procedure, public :: put => putValuesDataGroup
            !
            procedure, public :: set => setValuesDataGroup
            !
            procedure, public :: sub => subDataGroup
            !
            procedure, public :: linComb => linCombDataGroup
            !
            procedure, public :: dotProd => dotProdDataGroup
            !
            procedure, public :: normalize => normalizeDataGroup
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
    function DataGroup_ctor( i_rx, i_tx, n_comp, is_complex, error_bar ) result( self )
        implicit none
        !
        integer, intent( in ) :: i_rx, i_tx, n_comp
        logical, optional, intent( in ) :: is_complex, error_bar
        !
        type( DataGroup_t ) :: self
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
        if( present( error_bar ) ) then
            !
            self%error_bar = error_bar
        else
            self%error_bar = .FALSE.
        endif
        !
        if( present( is_complex ) ) then
            !
            self%is_complex = is_complex
        else
            self%is_complex = .FALSE.
            !
        endif
        !
        self%normalized = 0
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
    !> Set the values at a given index of these arrays.
    subroutine setValuesDataGroup( self, comp_id, rvalue, imaginary, error )
        implicit none
        !
        class( DataGroup_t ), intent( inout ) :: self
        integer, intent( in ) :: comp_id
        real( kind=prec ), intent( in ) :: rvalue, imaginary
        real( kind=prec ), optional, intent( in ) :: error
        !
        self%reals( comp_id ) = rvalue
        !
        self%imaginaries( comp_id ) = imaginary
        !
        self%error_bar = present( error )
        !
        if( self%error_bar ) then
            !
            self%errors( comp_id ) = error
        else
            self%errors( comp_id ) = R_ZERO
        endif
        !
    end subroutine setValuesDataGroup
    !
    !> ????
    subroutine linCombDataGroup( self, a, b, d_in, d_out )
        implicit none
        !
        class( DataGroup_t ), intent( in ) :: self, d_in
        real( kind=prec ), intent( in ) :: a, b
        class( DataGroup_t ), intent( inout ) :: d_out
        !
        complex( kind=prec ), allocatable, dimension(:) :: c_self, c_d_in, c_d_out
        !
        if( .NOT. self%is_allocated ) then
            stop "Error: linCombDataGroup > self not allocated"
        endif
        !
        if( .NOT. d_in%is_allocated ) then
            stop "Error: linCombDataGroup > d_in not allocated"
        endif
        !
        if( .NOT. d_out%is_allocated ) then
            stop "Error: linCombDataGroup > d_out not allocated"
        endif
        !
        if( self%i_tx /= d_in%i_tx  ) then
            stop "Error: linCombDataGroup > different data txs: d1 and d2"
        endif
        !
        if( self%i_rx /= d_in%i_rx  ) then
            stop "Error: linCombDataGroup > different data rxs: d1 and d2"
        endif
        !
        if( self%n_comp /= d_in%n_comp  ) then
            stop "Error: linCombDataGroup > different data n_comp: d1 and d2"
        endif
        !
        if( self%i_tx /= d_out%i_tx  ) then
            stop "Error: linCombDataGroup > different data txs: d1 and d_out"
        endif
        !
        if( self%i_rx /= d_out%i_rx  ) then
            stop "Error: linCombDataGroup > different data rxs: d1 and d_out"
        endif
        !
        if( self%n_comp /= d_out%n_comp  ) then
            stop "Error: linCombDataGroup > different data n_comp: d1 and d_out"
        endif
        !
        d_out%error_bar = self%error_bar .OR. d_in%error_bar
        d_out%normalized = 0
        !
        d_out%i_rx = self%i_rx
        d_out%i_tx = self%i_tx
        d_out%n_comp = self%n_comp
        !
        if( d_out%is_complex ) then
            !
            c_self = cmplx( self%reals, self%imaginaries, kind=prec )
            c_d_in = cmplx( d_in%reals, d_in%imaginaries, kind=prec )
            !
            c_d_out = a * c_self + b * c_d_in
            !
            d_out%reals = real( c_d_out, kind=prec )
            d_out%imaginaries = real( aimag( c_d_out ), kind=prec )
            !
        else
            d_out%reals = a * self%reals + b * d_in%reals
        endif
        !
        if( self%error_bar .AND. d_in%error_bar ) then
            !
            if( abs(a) > R_ZERO .AND. abs(b) > R_ZERO ) then
                stop "Error: linCombDataGroup: unable to add two data vectors with error bars"
            else if( abs(a) > R_ZERO ) then
                d_out%errors = a * self%errors
                d_out%normalized = self%normalized
            else if( abs(b) > R_ZERO ) then
                d_out%errors = b * d_in%errors
                d_out%normalized = d_in%normalized
            endif
            !
        else if( self%error_bar ) then
            !
            d_out%errors = a * self%errors
            d_out%normalized = self%normalized
            !
        else if( d_in%error_bar ) then
            !
            d_out%errors = b * d_in%errors
            d_out%normalized = d_in%normalized
            !
        endif
        !
    end subroutine linCombDataGroup
    !
    !> ????
    function dotProdDataGroup( self, d_in ) result( rvalue )
        implicit none
        !
        class( DataGroup_t ), intent( in ) :: self, d_in
        !
        real( kind=prec ) :: rvalue
        !
        complex( kind=prec ) :: c_value
        complex( kind=prec ), allocatable, dimension(:) :: c_self, c_d_in
        !
        if( self%is_complex ) then
            !
            c_self = cmplx( self%reals, self%imaginaries, kind=prec )
            c_d_in = cmplx( d_in%reals, d_in%imaginaries, kind=prec )
            !
            c_value = sum( conjg( c_self ) * c_d_in )
            !
            rvalue = real( c_value, kind=prec )
            !
        else
            rvalue = sum( self%reals * d_in%reals )
        endif
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
    !>
    subroutine normalizeDataGroup( self, norm )
        implicit none
        !
        class( DataGroup_t ), intent( inout ) :: self
        integer, intent( in ), optional :: norm
        !
        complex( kind=prec ), allocatable, dimension(:) :: c_self
        !
        integer :: nn
        !
        if( .NOT. self%error_bar ) then
            stop "Error: normalizeDataGroup: no error bars to normalize"
        endif
        !
        if( present( norm ) ) then
            nn = norm
        else
            nn = 1
        endif
        !
        if( self%is_complex ) then
            !
            c_self = cmplx( self%reals, self%imaginaries, kind=prec )
            !
            c_self = c_self / ( self%errors ** nn )
            !
            self%reals = real( c_self, kind=prec )
            self%imaginaries = real( aimag( c_self ), kind=prec )
            !
        else
            self%reals = self%reals / ( self%errors ** nn )
        endif
        !
        self%normalized = self%normalized + nn
        !
    end subroutine normalizeDataGroup
    !
    !> Subtraction
    subroutine subDataGroup( self, d_in )
        implicit none
        !
        class( DataGroup_t ), intent( inout ) :: self
        class( DataGroup_t ), intent( in ) :: d_in
        !
        complex( kind=prec ), allocatable, dimension(:) :: c_self, c_d_in
        !
        if( self%is_complex ) then
            !
            c_self = cmplx( self%reals, self%imaginaries, kind=prec )
            c_d_in = cmplx( d_in%reals, d_in%imaginaries, kind=prec )
            !
            c_self = c_self - c_d_in
            !
            self%reals = real( c_self, kind=prec )
            self%imaginaries = real( aimag( c_self ), kind=prec )
            !
        else
            self%reals = self%reals - d_in%reals
        endif
        !
    end subroutine subDataGroup
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
    subroutine copyFromDataGroup( self, d_in )
        implicit none
        !
        class( DataGroup_t ), intent( inout ) :: self
        type( DataGroup_t ), intent( in ) :: d_in
        !
        integer :: i, asize
        !
        if( .NOT. d_in%is_allocated ) then
            stop "Error: copyFromDataGroup > d_in not allocated"
        endif
        !
        self%n_comp = d_in%n_comp
        self%i_rx = d_in%i_rx
        self%i_tx = d_in%i_tx
        self%counter = d_in%counter
        self%error_bar = d_in%error_bar
        !
        if( allocated( self%components ) ) then
            asize = size( self%components )
            do i = asize, 1, -(1)
                if( allocated( self%components(i)%str ) ) deallocate( self%components(i)%str )
            enddo
            deallocate( self%components )
        endif
        !
        allocate( self%components, source = d_in%components )
        !
        if( allocated( self%reals ) ) deallocate( self%reals )
        allocate( self%reals, source = d_in%reals )
        !
        if( allocated( self%imaginaries ) ) deallocate( self%imaginaries )
        allocate( self%imaginaries, source = d_in%imaginaries )
        !
        if( allocated( self%errors ) ) deallocate( self%errors )
        !
        if( d_in%error_bar ) then
            !
            allocate( self%errors, source = d_in%errors )
        else
            !
            allocate( self%errors( size( d_in%errors ) ) )
            !
            self%errors = R_ZERO
        endif
        !
        self%is_allocated = d_in%is_allocated
        !
        self%is_complex = d_in%is_complex
        !
        self%normalized = d_in%normalized
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
