!
!> Lone class to define a group of n_comp entries, related to a specific transmitter-receiver pair.
!
module DataGroup
    !
    use String
    use Utilities
    !
    !> Global file path name for data files
    character(:), allocatable :: predicted_data_file_name, jmhat_data_file_name
    !
    type :: DataGroup_t
        !
        integer :: i_dg, i_rx, i_tx, n_comp
        !
        real( kind=prec ), allocatable, dimension(:) :: reals, imaginaries, errors
        !
        logical :: is_allocated, is_complex, error_bar
        !
        integer, private :: counter
        !
        contains
            !
            final :: DataGroup_dtor
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
    !
    function DataGroup_ctor( i_rx, i_tx, n_comp, error_bar ) result( self )
        implicit none
        !
        integer, intent( in ) :: i_rx, i_tx, n_comp
        logical, optional, intent( in ) :: error_bar
        !
        type( DataGroup_t ) :: self
        !
        integer :: i, asize
        !
        self%i_dg = 0
        !
        self%i_rx = i_rx
        !
        self%i_tx = i_tx
        !
        self%n_comp = n_comp
        !
        self%counter = 1
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
        self%is_complex = .FALSE.
        !
        if( present( error_bar ) ) then
            !
            self%error_bar = error_bar
        else
            self%error_bar = .FALSE.
        endif
        !
    end function DataGroup_ctor
    !
    !> No subroutine briefing
    !
    subroutine DataGroup_dtor( self )
        implicit none
        !
        type( DataGroup_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor DataGroup"
        !
        if( allocated( self%reals ) ) deallocate( self%reals )
        if( allocated( self%imaginaries ) ) deallocate( self%imaginaries )
        if( allocated( self%errors ) ) deallocate( self%errors )
        !
        self%is_allocated = .FALSE.
        !
    end subroutine DataGroup_dtor
    !
    !> Add values to arrays in position and increments the internal counter.
    !
    subroutine putValuesDataGroup( self, rvalue, imaginary, error )
        implicit none
        !
        class( DataGroup_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: rvalue, imaginary, error
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
    !
    subroutine setValuesDataGroup( self, comp_id, rvalue, imaginary )
        implicit none
        !
        class( DataGroup_t ), intent( inout ) :: self
        integer, intent( in ) :: comp_id
        real( kind=prec ), intent( in ) :: rvalue, imaginary
        !
        self%reals( comp_id ) = rvalue
        !
        self%imaginaries( comp_id ) = imaginary
        !
    end subroutine setValuesDataGroup
    !
    !> ????
    !
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
    !
    subroutine normalizeDataGroup( self, norm )
        implicit none
        !
        class( DataGroup_t ), intent( inout ) :: self
        integer, intent( in ), optional :: norm
        !
        integer :: nn
        !
        if( .NOT. self%error_bar ) then
            call errStop( "normalizeDataGroup > no error bars to normalize" )
        endif
        !
        if( present( norm ) ) then
            nn = norm
        else
            nn = 1
        endif
        !
        self%reals = self%reals / ( self%errors ** nn )
        !
        self%imaginaries = self%imaginaries / ( self%errors ** nn )
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
        self%reals = self%reals - d_in%reals
        !
        self%imaginaries = self%imaginaries - d_in%imaginaries
        !
    end subroutine subDataGroup
    !
    !> ????
    !
    subroutine linCombDataGroup( self, a, b, d_in, d_out )
        implicit none
        !
        class( DataGroup_t ), intent( in ) :: self, d_in
        real( kind=prec ), intent( in ) :: a, b
        class( DataGroup_t ), intent( inout ) :: d_out
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "linCombDataGroup > self not allocated" )
        endif
        !
        if( .NOT. d_in%is_allocated ) then
            call errStop( "linCombDataGroup > d_in not allocated" )
        endif
        !
        if( .NOT. d_out%is_allocated ) then
            call errStop( "linCombDataGroup > d_out not allocated" )
        endif
        !
        if( self%i_tx /= d_in%i_tx  ) then
            call errStop( "linCombDataGroup > different data txs: d1 and d2" )
        endif
        !
        if( self%i_rx /= d_in%i_rx  ) then
            call errStop( "linCombDataGroup > different data rxs: d1 and d2" )
        endif
        !
        if( self%n_comp /= d_in%n_comp  ) then
            call errStop( "linCombDataGroup > different data n_comp: d1 and d2" )
        endif
        !
        if( self%i_tx /= d_out%i_tx  ) then
            call errStop( "linCombDataGroup > different data txs: d1 and d_out" )
        endif
        !
        if( self%i_rx /= d_out%i_rx  ) then
            call errStop( "linCombDataGroup > different data rxs: d1 and d_out" )
        endif
        !
        if( self%n_comp /= d_out%n_comp  ) then
            call errStop( "linCombDataGroup > different data n_comp: d1 and d_out" )
        endif
        !
        d_out%error_bar = self%error_bar .OR. d_in%error_bar
        !
        d_out%is_complex = self%is_complex .AND. d_in%is_complex
        !
        d_out%i_dg = self%i_dg
        d_out%i_rx = self%i_rx
        d_out%i_tx = self%i_tx
        d_out%n_comp = self%n_comp
        !
        d_out%reals = a * self%reals + b * d_in%reals
        !
        d_out%imaginaries = a * self%imaginaries + b * d_in%imaginaries
        !
        if( self%error_bar .AND. d_in%error_bar ) then
            !
            if( abs(a) > R_ZERO .AND. abs(b) > R_ZERO ) then
                call errStop( "linCombDataGroup > unable to add two data vectors with error bars" )
            elseif( abs(a) > R_ZERO ) then
                !
                d_out%errors = a * self%errors
                !
            elseif( abs(b) > R_ZERO ) then
                !
                d_out%errors = b * d_in%errors
                !
            endif
            !
        elseif( self%error_bar ) then
            !
            d_out%errors = a * self%errors
            !
        elseif( d_in%error_bar ) then
            !
            d_out%errors = b * d_in%errors
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
        rvalue = R_ZERO
        !
        rvalue = rvalue + sum( self%reals * d_in%reals )
        rvalue = rvalue + sum( self%imaginaries * d_in%imaginaries )
        !
        !write( *, * ) "COMPLEX DOT PRODUCT: ", self%is_complex .AND. d_in%is_complex
        !if( self%is_complex .AND. d_in%is_complex ) then
            !rvalue = rvalue + sum( self%imaginaries * d_in%imaginaries )
        !endif
        !
    end function dotProdDataGroup
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
            if( ABS( self%reals( i_comp ) - other%reals( i_comp ) ) >= TOL6 .OR. &
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
        if( .NOT. d_in%is_allocated ) then
            call errStop( "copyFromDataGroup > d_in not allocated" )
        endif
        !
        self%i_dg = d_in%i_dg
        !
        self%i_rx = d_in%i_rx
        !
        self%i_tx = d_in%i_tx
        !
        self%n_comp = d_in%n_comp
        !
        self%counter = d_in%counter
        !
        self%error_bar = d_in%error_bar
        !
        if( allocated( self%reals ) ) deallocate( self%reals )
        allocate( self%reals, source = d_in%reals )
        !
        if( allocated( self%imaginaries ) ) deallocate( self%imaginaries )
        allocate( self%imaginaries, source = d_in%imaginaries )
        !
        if( allocated( self%errors ) ) deallocate( self%errors )
        allocate( self%errors, source = d_in%errors )
        !
        self%is_allocated = d_in%is_allocated
        !
        self%is_complex = d_in%is_complex
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
        write( *, * ) "    Write DataGroup_t Id: ", self%i_dg
        write( *, * ) "             Receiver Id: ", self%i_rx
        write( *, * ) "          Transmitter Id: ", self%i_tx
        write( *, * ) "     Error bar & Complex: ", self%error_bar, self%is_complex
        write( *, * ) self%n_comp, " data_rows:"
        !
        do i_comp = 1, self%n_comp
            !
            write( *, * ) i_comp, ":", self%reals( i_comp ), self%imaginaries( i_comp ), self%errors( i_comp )
            !
        enddo
        !
    end subroutine printDataGroup
    !
end module DataGroup
