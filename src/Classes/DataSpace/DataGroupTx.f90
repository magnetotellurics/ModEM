!
!> Type that encapsulates a DataGroup array under the same index as a transmitter.
!
module DataGroupTx
    !
    use DataGroup
    !
    type :: DataGroupTx_t
        !
        integer :: i_tx
        !
        type( DataGroup_t ), allocatable, dimension(:) :: data
        !
    contains
            !
            final :: DataGroupTx_dtor
            !
            procedure, public :: zeros => zerosDataGroupTx
            !
            procedure, public :: getRxData => getRxDataDataGroupTx
            !
            procedure, public :: put => putDataGroupTx
            !
            procedure, public :: set => setDataGroupTx
            !
            procedure, public :: setValues => setValuesDataGroupTx
            !
            procedure, public :: sub => subDataGroupTx
            !
            procedure, public :: linComb => linCombDataGroupTx
            !
            procedure, public :: dotProd => dotProdDataGroupTx
            !
            procedure, public :: normalize => normalizeDataGroupTx
            !
            procedure, public :: print => printDataGroupTx
            !
    end type DataGroupTx_t
    !
    interface DataGroupTx_t
         module procedure DataGroupTx_ctor
    end interface DataGroupTx_t
    !
contains
    !
    !> Parametrized constructor:
    !> Set the transmitter index and deallocate the data array if it was previously allocated
    function DataGroupTx_ctor( i_tx ) result( self )
        implicit none
        !
        integer, intent( in ) :: i_tx
        !
        type( DataGroupTx_t ) :: self
        !
        !write( *, * ) "Constructor DataGroupTx_t: ", i_tx
        !
        self%i_tx = i_tx
        !
        if( allocated( self%data ) ) deallocate( self%data )
        !
    end function DataGroupTx_ctor
    !
    !> No subroutine briefing
    subroutine DataGroupTx_dtor( self )
        implicit none
        !
        type( DataGroupTx_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor DataGroupTx"
        !
        if( allocated( self%data ) ) deallocate( self%data )
        !
    end subroutine DataGroupTx_dtor
    !
    !> sub
    subroutine subDataGroupTx( self, rhs )
        implicit none
        !
        class( DataGroupTx_t ), intent( inout ) :: self
        class( DataGroupTx_t ), intent( in ) :: rhs
        !
        integer :: i_data
        !
        if( size( self%data ) /= size( rhs%data ) ) then
            !
            stop "Error: subDataGroupTx > different data sizes"
            !
        else
            !
            do i_data = 1, size( self%data )
                !
                call self%data( i_data )%sub( rhs%data( i_data ) )
                !
            enddo
            !
        endif
        !
    end subroutine subDataGroupTx
    !
    !> ????
    subroutine normalizeDataGroupTx( self, norm )
        implicit none
        !
        class( DataGroupTx_t ), intent( inout ) :: self
        !
        integer, intent( in ), optional :: norm
        !
        integer :: i, nn
        !
        if( present( norm ) ) then
            nn = norm
        else
            nn = 1
        endif
        !
        do i = 1, size( self%data )
            !
            call self%data(i)%normalize( nn )
            !
        enddo
        !
    end subroutine normalizeDataGroupTx
    !
    !> Call reset for each DataGroup in data
    subroutine zerosDataGroupTx( self )
        implicit none
        !
        class( DataGroupTx_t ), intent( inout ) :: self
        !
        integer :: i_data
        !
        do i_data = 1, size( self%data )
            !
            call self%data( i_data )%zeros()
            !
        enddo
        !
    end subroutine zerosDataGroupTx
    !
    !> Returns a pointer, allowing modifications directly to a DataGroup at a given transmitter-receiver pair indexes
    function getRxDataDataGroupTx( self, i_rx ) result( data_group )
        implicit none
        !
        class( DataGroupTx_t ), intent( in ) :: self
        integer, intent( in ) :: i_rx
        !
        type( DataGroup_t ) :: data_group
        !
        integer :: i_data, n_data
        !
        n_data = size( self%data )
        !
        do i_data = 1, n_data
            !
            if( self%data( i_data )%i_rx == i_rx ) then
                !
                data_group = self%data( i_data )
                !
                return
                !
            endif
            !
        enddo
        !
    end function getRxDataDataGroupTx
    !
    !> Dynamically add a new DataGroup to the array, always via reallocation.
    subroutine putDataGroupTx( self, data_group )
        implicit none
        !
        class( DataGroupTx_t ), intent( inout ) :: self
        type( DataGroup_t ), intent( in ) :: data_group
        !
        integer :: i_data, n_data
        type( DataGroup_t ), allocatable, dimension(:) :: temp_array
        !
        if( .NOT. allocated( self%data ) ) then
            !
            allocate( self%data(1) )
            !
            self%data(1) = data_group
            !
        else
            !
            n_data = size( self%data )
            !
            allocate( temp_array( n_data + 1 ) )
            !
            temp_array( 1 : n_data ) = self%data(:)
            !
            temp_array( n_data + 1 ) = data_group
            !
            deallocate( self%data )
            !
            allocate( self%data, source = temp_array )
            !
            deallocate( temp_array )
            !
        endif
        !
    end subroutine putDataGroupTx
    !
    !> Replace the values of a specific DataGroup of the array 
    !> by other DataGroup with the same transmitter-receiver pair
    subroutine setValuesDataGroupTx( self, data_group )
        implicit none
        !
        class( DataGroupTx_t ), intent( inout ) :: self
        !
        type( DataGroup_t ), intent( in ) :: data_group
        !
        integer :: i_data
        !
        do i_data = 1, size( self%data )
            !
            if( self%data( i_data )%i_rx == data_group%i_rx .AND. &
                self%data( i_data )%i_tx == data_group%i_tx ) then
                !
                self%data( i_data )%reals = data_group%reals
                self%data( i_data )%imaginaries = data_group%imaginaries
                self%data( i_data )%errors = data_group%errors
                !
                self%data( i_data )%error_bar = data_group%error_bar
                !
                return
                !
            endif
            !
        enddo
        !
    end subroutine setValuesDataGroupTx
    !
    !> Replace a specific DataGroup of the array 
    !> by other DataGroup with the same transmitter-receiver pair
    subroutine setDataGroupTx( self, data_group )
        implicit none
        !
        class( DataGroupTx_t ), intent( inout ) :: self
        !
        type( DataGroup_t ), intent( in ) :: data_group
        !
        integer :: i_data
        !
        do i_data = 1, size( self%data )
            !
            if( self%data( i_data )%i_rx == data_group%i_rx .AND. &
                self%data( i_data )%i_tx == data_group%i_tx ) then
                !
                self%data( i_data ) = data_group
                !
                return
                !
            endif
            !
        enddo
        !
    end subroutine setDataGroupTx
    !
    !> dotProd between two DataGroupTxs
    function dotProdDataGroupTx( self, data_tx ) result( rvalue )
        implicit none
        !
        class( DataGroupTx_t ), intent( in ) :: self, data_tx
        !
        real( kind=prec ) :: rvalue
        !
        integer :: i
        !
        if( size( self%data ) /= size( data_tx%data ) ) then
            !
            stop "Error: DataGroupTx_t : dotProdDataGroupTx > different data sizes"
            !
        else
            !
            rvalue = R_ZERO
            !
            do i = 1, size( self%data )
                !
                rvalue = rvalue + self%data(i)%dotProd( data_tx%data(i) )
                !
            enddo
            !
        endif
        !
    end function dotProdDataGroupTx
    !
    !> linComb between two DataGroupTxs
    subroutine linCombDataGroupTx( self, a, b, data_tx, data_tx_out )
        implicit none
        !
        class( DataGroupTx_t ), intent( in ) :: self, data_tx
        real( kind=prec ), intent( in ) :: a, b
        class( DataGroupTx_t ), intent( inout ) :: data_tx_out
        !
        integer :: i
        !
        if( self%i_tx /= data_tx%i_tx ) then
            stop "Error: DataGroupTx_t : linCombDataGroupTx > different data txs: d1 and d2"
        endif
        !
        if( size( self%data ) /= size( data_tx%data ) ) then
            stop "Error: DataGroupTx_t : linCombDataGroupTx > different data sizes: d1 and d2"
        endif
        !
        if( self%i_tx /= data_tx_out%i_tx ) then
            stop "Error: DataGroupTx_t : linCombDataGroupTx > different data txs: d1 and d_out"
        endif
        !
        if( size( self%data ) /= size( data_tx_out%data ) ) then
            stop "Error: DataGroupTx_t : linCombDataGroupTx > different data sizes: d1 and d_out"
        endif
        !
        do i = 1, size( self%data )
            !
            call self%data(i)%linComb( a, b, data_tx%data(i), data_tx_out%data(i) )
            !
        enddo
        !
    end subroutine linCombDataGroupTx
    !
    !> Call the print routine of each DataGroup in the data array
    subroutine printDataGroupTx( self )
        implicit none
        !
        class( DataGroupTx_t ), intent( in ) :: self
        !
        integer :: i_data
        !
        write( *, * ) "    Write DataGroupTx_t for Tx: ", self%i_tx
        !
        do i_data = 1, size( self%data )
            call self%data( i_data )%print()
        enddo
        !
    end subroutine printDataGroupTx
    !
end module DataGroupTx
!