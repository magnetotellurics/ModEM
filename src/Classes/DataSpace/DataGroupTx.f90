!
!> Type that encapsulates a DataGroup array under the same index as a transmitter.
!
module DataGroupTx
    !
    use DataGroupArray
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
            procedure, public :: rmsdDataGroupTx1
            procedure, public :: rmsdDataGroupTx2
            generic :: rmsd => rmsdDataGroupTx1, rmsdDataGroupTx2
            !
            procedure, public :: multAdd => multAddDataGroupTx
            !
            procedure, public :: linComb => linCombDataGroupTx
            !
            procedure, public :: dotProd => dotProdDataGroupTx
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
    !> Root Mean Square Deviation between two DataGroupTxs
    subroutine multAddDataGroupTx( self, data_tx, rvalue )
        implicit none
        !
        class( DataGroupTx_t ), intent( inout ) :: self
        !
        class( DataGroupTx_t ), intent( in ) :: data_tx
        !
        real( kind=prec ), intent( in ) :: rvalue
        !
        integer :: i_data, i_comp
        !
        complex( kind=prec ) :: self_comp, data_tx_comp
        !
        if( size( self%data ) /= size( data_tx%data ) ) then
            !
            stop "Error: multAddDataGroupTx > different data sizes"
            !
        else
            !
            do i_data = 1, size( self%data )
                !
                do i_comp = 1, self%data( i_data )%n_comp
                    !
                    self_comp = cmplx( self%data( i_data )%reals( i_comp ), self%data( i_data )%imaginaries( i_comp ), kind=prec )
                    !
                    data_tx_comp = cmplx( data_tx%data( i_data )%reals( i_comp ), data_tx%data( i_data )%imaginaries( i_comp ), kind=prec )
                    !
                    call self%data( i_data )%set( i_comp, conjg( self_comp ) + rvalue * data_tx_comp )
                    !
                enddo
                !
            enddo
            !
        endif
        !
    end subroutine multAddDataGroupTx
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
            rvalue = 0.0
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
        if( size( self%data ) /= size( data_tx%data ) ) then
            !
            stop "Error: DataGroupTx_t : linCombDataGroupTx > different data sizes"
            !
        else
            !
            do i = 1, size( self%data )
                !
                call self%data(i)%linComb( a, b, data_tx%data(i), data_tx_out%data(i) )
                !
            enddo
            !
        endif
        !
    end subroutine linCombDataGroupTx
    !
    !> Root Mean Square Deviation for a single DataGroupTx
    function rmsdDataGroupTx1( self ) result( rmsd )
        implicit none
        !
        class( DataGroupTx_t ), intent( in ) :: self
        !
        complex( kind=prec ) :: rmsd
        !
        integer :: i_data, i_comp
        !
        complex( kind=prec ) :: self_comp
        !
        rmsd = C_ZERO
        !
        do i_data = 1, size( self%data )
            !
            do i_comp = 1, self%data( i_data )%n_comp
                !
                self_comp = cmplx( self%data( i_data )%reals( i_comp ), self%data( i_data )%imaginaries( i_comp ), kind=prec )
                !
                rmsd = rmsd + self_comp ** 2
                !
            enddo
            !
        enddo
        !
    end function rmsdDataGroupTx1
    !
    !> Root Mean Square Deviation between two DataGroupTxs
    function rmsdDataGroupTx2( self, data_tx ) result( rmsd )
        implicit none
        !
        class( DataGroupTx_t ), intent( in ) :: self, data_tx
        !
        complex( kind=prec ) :: rmsd
        !
        integer :: i_data, i_comp
        !
        complex( kind=prec ) :: self_comp, data_tx_comp
        !
        if( size( self%data ) /= size( data_tx%data ) ) then
            !
            stop "Error: DataGroupTx_t : rmsdDataGroupTx2 > different data sizes"
            !
        else
            !
            rmsd = C_ZERO
            !
            do i_data = 1, size( self%data )
                !
                do i_comp = 1, self%data( i_data )%n_comp
                    !
                    self_comp = cmplx( self%data( i_data )%reals( i_comp ), self%data( i_data )%imaginaries( i_comp ), kind=prec )
                    !
                    data_tx_comp = cmplx( data_tx%data( i_data )%reals( i_comp ), data_tx%data( i_data )%imaginaries( i_comp ), kind=prec )
                    !
                    rmsd = rmsd + ( ( self_comp - data_tx_comp ) ** 2 )
                    !
                enddo
                !
            enddo
            !
        endif
        !
    end function rmsdDataGroupTx2
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