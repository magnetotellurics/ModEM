!
!> Module with routines to handle a DataGroupTx array
!
module DataGroupTxArray
    !
    use DataGroupTx
    !
    !> Global arrays of Predicted and Residual data for all transmitters
    type( DataGroupTx_t ), allocatable, dimension(:) :: all_predicted_data, all_residual_data
    !
    public :: multAddDataGroupTxArray
    public :: dotProdDataGroupTxArray
    public :: linCombDataGroupTxArray
    public :: scMultAddDataGroupTxArray
    public :: normalizeDataGroupTxArray
    public :: normalizeWithDataGroupTxArray
    !
    interface rmsdDataGroupTxArray
        module procedure rmsdDataGroupTxArray1
        module procedure rmsdDataGroupTxArray2
    end interface rmsdDataGroupTxArray
    !
    public :: countDataGroupTxArray
    public :: getDataGroupTxArray
    public :: updateDataGroupTxArray
    public :: deallocateDataGroupTxArray
    public :: zerosDataGroupTxArray
    public :: printDataGroupTxArray
    !
contains
    !
    !> ????
    subroutine normalizeDataGroupTxArray( data_tx_array, norm )
        implicit none
        !
        type( DataGroupTx_t ), dimension(:), intent( inout ) :: data_tx_array
        integer, intent(in) :: norm
        !
        complex( kind=prec ) :: c_value
        integer :: i, j, k
        !
        do i = 1, size( data_tx_array )
            !
            do j = 1, size( data_tx_array(i)%data )
                !
                do k = 1, data_tx_array(i)%data(j)%n_comp
                    !
                    c_value = cmplx( data_tx_array(i)%data(j)%reals(k), data_tx_array(i)%data(j)%imaginaries(k), kind=prec )
                    !
                    c_value = c_value / data_tx_array(i)%data(j)%errors(k) ** norm
                    !
                    call data_tx_array(i)%data(j)%set( k, c_value )
                    !
                enddo
                !
            enddo
            !
        enddo
        !
    end subroutine normalizeDataGroupTxArray
    !
    !> ????
    subroutine normalizeWithDataGroupTxArray( norm, data_tx_array_in, data_tx_array_out )
        implicit none
        !
        integer, intent( in ) :: norm
        type( DataGroupTx_t ), dimension(:), intent( in ) :: data_tx_array_in
        type( DataGroupTx_t ), dimension(:), intent( inout ) :: data_tx_array_out
        !
        complex( kind=prec ) :: c_value
        integer :: i, j, k
        !
        do i = 1, size( data_tx_array_in )
            !
            do j = 1, size( data_tx_array_in(i)%data )
                !
                data_tx_array_out(i)%data(j)%error_bar = .TRUE.
                !
                do k = 1, data_tx_array_in(i)%data(j)%n_comp
                    !
                    c_value = cmplx( data_tx_array_out(i)%data(j)%reals(k), data_tx_array_out(i)%data(j)%imaginaries(k), kind=prec )
                    !
                    c_value = c_value / data_tx_array_in(i)%data(j)%errors(k) ** norm
                    !
                    call data_tx_array_out(i)%data(j)%set( k, c_value )
                    !
                enddo
                !
            enddo
            !
        enddo
        !
    end subroutine normalizeWithDataGroupTxArray
    !
    !> ????
    function countDataGroupTxArray( data_tx_array ) result( counter )
        implicit none
        !
        type( DataGroupTx_t ), dimension(:), intent( in ) :: data_tx_array
        integer :: counter
        !
        integer :: i, j
        !
        counter = 0
        !
        do i = 1, size( data_tx_array )
            !
            do j = 1, size( data_tx_array(i)%data )
                !
                counter = counter + data_tx_array(i)%data(j)%n_comp
                !
            enddo
            !
        enddo
        !
    end function countDataGroupTxArray
    !
    !> Root Mean Square Deviation between two DataGroupTxArrays
    subroutine multAddDataGroupTxArray( data_tx_array_1, data_tx_array_2, rvalue )
        implicit none
        !
        type( DataGroupTx_t ), dimension(:), intent( inout ) :: data_tx_array_1
        type( DataGroupTx_t ), dimension(:), intent( in ) :: data_tx_array_2
        real( kind=prec ), intent( in ) :: rvalue
        !
        integer :: i
        !
        if( size( data_tx_array_1 ) /= size( data_tx_array_2 ) ) then
            !
            stop "Error: DataGroupTxArray : multAddDataGroupTxArray > different array sizes"
            !
        else
            !
            do i = 1, size( data_tx_array_1 )
                !
                call data_tx_array_1(i)%multAdd( data_tx_array_2(i), rvalue )
                !
            enddo
            !
        endif
        !
    end subroutine multAddDataGroupTxArray
    !
    !> Root Mean Square Deviation between two DataGroupTxArrays
    function dotProdDataGroupTxArray( data_tx_array_1, data_tx_array_2 ) result( rvalue )
        implicit none
        !
        type( DataGroupTx_t ), dimension(:), intent( in ) :: data_tx_array_1, data_tx_array_2
        !
        real( kind=prec ) :: rvalue
        !
        integer :: i
        !
        if( size( data_tx_array_1 ) /= size( data_tx_array_2 ) ) then
            !
            stop "Error: DataGroupTxArray : dotProdDataGroupTxArray > different array sizes"
            !
        else
            !
            rvalue = 0.0
            !
            do i = 1, size( data_tx_array_1 )
                !
                rvalue = rvalue + data_tx_array_1(i)%dotProd( data_tx_array_2(i) )
                !
            enddo
            !
        endif
        !
    end function dotProdDataGroupTxArray
    !
    !> ????
    subroutine linCombDataGroupTxArray( a, d1, b, d2, dOut )
        implicit none
        !
        type( DataGroupTx_t ), dimension(:), intent( in ) :: d1, d2
        real(kind=prec), intent( in ) :: a, b
        type( DataGroupTx_t ), dimension(:), intent( inout ) :: dOut
        !
        integer :: i
        !
        if( size( d1 ) /= size( d2 ) ) then
            !
            stop "Error: DataGroupTxArray : linCombDataGroupTxArray > different array sizes"
            !
        else
            !
            do i = 1, size( d1 )
                !
                call d1(i)%linComb( a, b, d2(i), dOut(i) )
                !
            enddo
            !
        endif
        !
    end subroutine linCombDataGroupTxArray
    !
    !> ????
    subroutine scMultAddDataGroupTxArray( rvalue, data_tx_array_in, data_tx_array_out )
        implicit none
        !
        real( kind=prec ), intent( in ) :: rvalue
        type( DataGroupTx_t ), dimension(:), intent( in ) :: data_tx_array_in
        type( DataGroupTx_t ), dimension(:), intent( inout ) :: data_tx_array_out
        !
        call linCombDataGroupTxArray( rvalue, data_tx_array_in, ONE, data_tx_array_out, data_tx_array_out )
        !
    end subroutine scMultAddDataGroupTxArray
    !
    !> Root Mean Square Deviation for a single DataGroupTxArray
    function rmsdDataGroupTxArray1( data_tx_array ) result( rmsd )
        implicit none
        !
        type( DataGroupTx_t ), dimension(:), intent( in ) :: data_tx_array
        !
        complex( kind=prec ) :: rmsd
        !
        integer :: i
        !
        rmsd = C_ZERO
        !
        do i = 1, size( data_tx_array )
            !
            rmsd = rmsd + data_tx_array(i)%rmsd()
            !
        enddo
        !
        rmsd = CDSQRT( rmsd / countDataGroupTxArray( data_tx_array ) )
        !
    end function rmsdDataGroupTxArray1
    !
    !> Root Mean Square Deviation between two DataGroupTxArrays
    function rmsdDataGroupTxArray2( data_tx_array_1, data_tx_array_2 ) result( rmsd )
        implicit none
        !
        type( DataGroupTx_t ), dimension(:), intent( in ) :: data_tx_array_1, data_tx_array_2
        !
        complex( kind=prec ) :: rmsd
        !
        integer :: i
        !
        if( size( data_tx_array_1 ) /= size( data_tx_array_2 ) ) then
            !
            stop "Error: DataGroupTxArray : rmsdDataGroupTxArray2 > different array sizes"
            !
        else
            !
            rmsd = C_ZERO
            !
            do i = 1, size( data_tx_array_1 )
                !
                rmsd = rmsd + data_tx_array_1(i)%rmsd( data_tx_array_2(i) )
                !
            enddo
            !
            rmsd = CDSQRT( rmsd / countDataGroupTxArray( data_tx_array_1 ) )
            !
        endif
        !
    end function rmsdDataGroupTxArray2
    !
    !> Return a pointer, allowing directly modifications to a DataGroupTx at a given index
    function getDataGroupTxArray( data_tx_array, dtx_id ) result( data_tx )
        implicit none
        !
        type( DataGroupTx_t ), target, dimension(:), intent( in ) :: data_tx_array
        integer, intent( in ) :: dtx_id
        !
        type( DataGroupTx_t ), pointer :: data_tx
        !
        data_tx => data_tx_array( dtx_id )
        !
    end function getDataGroupTxArray
    !
    !> Dynamically add a new DataGroupTx to the array, always via reallocation.
    subroutine updateDataGroupTxArray( data_tx_array, data_tx )
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: data_tx_array
        type( DataGroupTx_t ), intent( in ) :: data_tx
        !
        integer :: n_dtx
        type( DataGroupTx_t ), allocatable, dimension(:) :: temp_array
        !
        if( .NOT. allocated( data_tx_array ) ) then
            !
            allocate( data_tx_array(1) )
            !
            data_tx_array(1) = data_tx
            !
        else
            !
            n_dtx = size( data_tx_array )
            !
            allocate( temp_array( n_dtx + 1 ) )
            !
            temp_array( 1 : n_dtx ) = data_tx_array(:)
            !
            temp_array( n_dtx + 1 ) = data_tx
            !
            deallocate( data_tx_array )
            !
            allocate( data_tx_array, source = temp_array )
            !
            deallocate( temp_array )
            !
        endif
        !
    end subroutine updateDataGroupTxArray
    !
    !> No subroutine briefing
    subroutine deallocateDataGroupTxArray( data_tx_array )
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: data_tx_array
        !
        integer :: i, n_dtx
        !
        !write( *, * ) "deallocateDataGroupTxArray:", size( data_tx_array )
        !
        n_dtx = size( data_tx_array )
        !
        if( n_dtx == 1 ) then
            deallocate( data_tx_array(1)%data )
        else
            do i = n_dtx, 1, -(1)
                deallocate( data_tx_array(i)%data )
            enddo
        endif
        !
        deallocate( data_tx_array )
        !
    end subroutine deallocateDataGroupTxArray
    !
    !> Call the print routine of each DataGroupTx in the array
    subroutine zerosDataGroupTxArray( data_tx_array )
        implicit none
        !
        type( DataGroupTx_t ), dimension(:), intent( inout ) :: data_tx_array
        !
        integer :: i
        !
        do i = 1, size( data_tx_array )
            !
            call data_tx_array(i)%zeros()
            !
        enddo
        !
    end subroutine zerosDataGroupTxArray
    !
    !> Call the print routine of each DataGroupTx in the array
    subroutine printDataGroupTxArray( data_tx_array )
        implicit none
        !
        type( DataGroupTx_t ), dimension(:), intent( in ) :: data_tx_array
        !
        integer :: i
        !
        do i = 1, size( data_tx_array )
            !
            call data_tx_array(i)%print()
            !
        enddo
        !
    end subroutine printDataGroupTxArray
    !
end module DataGroupTxArray
!