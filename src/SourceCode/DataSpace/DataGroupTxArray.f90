!
!> Module to handle a DataGroupTx array
!> Store a global vector containing the measured data, read from data input file.
!
module DataGroupTxArray
    !
    use DataGroupTx
    !
    use FileUnits
    !
    use ReceiverArray
    !
    use TransmitterArray
    !
    !> Global array of data
    type( DataGroupTx_t ), allocatable, dimension(:) :: all_measured_data
    !
    !> Module variables
    logical :: conjugated_data
    type( String_t ), allocatable, dimension(:) :: units_in_file
    !
    !> Module routines
    !> Routines for data operations
    interface subData
        module procedure :: sub_DataGroupTxArray
    end interface subData
    !
    interface zerosData
        module procedure :: zeros_DataGroupTxArray
    end interface zerosData
    !
    interface dotProdData
        module procedure :: dotProd_DataGroupTxArray
    end interface dotProdData
    !
    interface linCombData
        module procedure :: linComb_DataGroupTxArray
    end interface linCombData
    !
    interface scMultData
        module procedure :: scMult_DataGroupTxArray
    end interface scMultData
    !
    interface scMultAddData
        module procedure :: scMultAdd_DataGroupTxArray
    end interface scMultAddData
    !
    interface normalizeData
        module procedure :: normalize_DataGroupTxArray
    end interface normalizeData
    !
    interface normalizeDataWith
        module procedure :: normalizeWith_DataGroupTxArray
    end interface normalizeDataWith
    !
    interface setErrorBar
        module procedure :: setErrorBar_DataGroupTxArray
    end interface setErrorBar
    !
    interface setComplex
        module procedure :: setComplex_DataGroupTxArray
    end interface setComplex
    !
    interface countData
        module procedure :: count_DataGroupTxArray
    end interface countData
    !
    interface countValues
        module procedure :: countValues_DataGroupTxArray
    end interface countValues
    !
    interface getData
        module procedure :: getDataGroupByIndex
    end interface getData
    !
    interface updateData
        module procedure :: update_DataGroupTxArray
    end interface updateData
    !
    interface deallocateData
        module procedure :: deallocate_DataGroupTxArray
    end interface deallocateData
    !
    interface writeData
        module procedure :: write_DataGroupTxArray
    end interface writeData
    !
    interface printData
        module procedure :: print_DataGroupTxArray
    end interface printData
    !
    private :: writeHeader_DataGroupTxArray
    !
contains
    !
    !> No subroutine briefing
    !
    subroutine sub_DataGroupTxArray( data_tx_array_1, data_tx_array_2 )
        implicit none
        !
        type( DataGroupTx_t ), dimension(:), intent( inout ) :: data_tx_array_1
        type( DataGroupTx_t ), dimension(:), intent( in ) :: data_tx_array_2
        !
        integer :: i
        !
        if( size( data_tx_array_1 ) /= size( data_tx_array_2 ) ) then
            !
            call errStop( "sub_DataGroupTxArray > different array sizes" )
            !
        else
            !
            do i = 1, size( data_tx_array_1 )
                !
                call data_tx_array_1(i)%sub( data_tx_array_2(i) )
                !
            enddo
            !
        endif
        !
    end subroutine sub_DataGroupTxArray
    !
    !> No subroutine briefing
    !
    subroutine zeros_DataGroupTxArray( data_tx_array )
        implicit none
        !
        type( DataGroupTx_t ), dimension(:), intent( inout ) :: data_tx_array
        !
        integer :: i
        !
        do i = 1, size( data_tx_array )
            !
            call data_tx_array(i)%zeros
            !
        enddo
        !
    end subroutine zeros_DataGroupTxArray
    !
    !> No function briefing
    !
    function dotProd_DataGroupTxArray( data_tx_array_1, data_tx_array_2 ) result( rvalue )
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
            call errStop( "dotProd_DataGroupTxArray > different array sizes" )
            !
        else
            !
            rvalue = R_ZERO
            !
            do i = 1, size( data_tx_array_1 )
                !
                rvalue = rvalue + data_tx_array_1(i)%dotProd( data_tx_array_2(i) )
                !
            enddo
            !
        endif
        !
    end function dotProd_DataGroupTxArray
    !
    !> No subroutine briefing
    !
    subroutine linComb_DataGroupTxArray( a, d1, b, d2, dOut )
        implicit none
        !
        type( DataGroupTx_t ), dimension(:), intent( in ) :: d1, d2
        real( kind=prec ), intent( in ) :: a, b
        type( DataGroupTx_t ), dimension(:), intent( inout ) :: dOut
        !
        integer :: i
        !
        if( size( d1 ) /= size( d2 ) ) then
            call errStop( "linComb_DataGroupTxArray > different array sizes: d1, d2" )
        endif
        !
        if( size( d1 ) /= size( dOut ) ) then
            call errStop( "linComb_DataGroupTxArray > different array sizes: d1, dOut" )
        endif
        !
        do i = 1, size( d1 )
            !
            call d1(i)%linComb( a, b, d2(i), dOut(i) )
            !
        enddo
        !
    end subroutine linComb_DataGroupTxArray
    !
    !> No subroutine briefing
    !
    subroutine scMult_DataGroupTxArray( rvalue, data_tx_array_in, data_tx_array_out )
        implicit none
        !
        real( kind=prec ), intent( in ) :: rvalue
        type( DataGroupTx_t ), dimension(:), intent( in ) :: data_tx_array_in
        type( DataGroupTx_t ), dimension(:), intent( inout ) :: data_tx_array_out
        !
        call linCombData( R_ZERO, data_tx_array_in, rvalue, data_tx_array_in, data_tx_array_out )
        !
    end subroutine scMult_DataGroupTxArray
    !
    !> No subroutine briefing
    !
    subroutine scMultAdd_DataGroupTxArray( rvalue, data_tx_array_in, data_tx_array_out )
        implicit none
        !
        real( kind=prec ), intent( in ) :: rvalue
        type( DataGroupTx_t ), dimension(:), intent( in ) :: data_tx_array_in
        type( DataGroupTx_t ), dimension(:), intent( inout ) :: data_tx_array_out
        !
        call linCombData( rvalue, data_tx_array_in, ONE, data_tx_array_out, data_tx_array_out )
        !
    end subroutine scMultAdd_DataGroupTxArray
    !
    !> No subroutine briefing
    !
    subroutine normalize_DataGroupTxArray( data_tx_array, norm )
        implicit none
        !
        type( DataGroupTx_t ), dimension(:), intent( inout ) :: data_tx_array
        integer, intent( in ), optional :: norm
        !
        integer :: i, j, nn
        !
        if( present( norm ) ) then
            nn = norm
        else
            nn = 1
        endif
        !
        do i = 1, size( data_tx_array )
            !
            do j = 1, size( data_tx_array(i)%data )
                !
                call data_tx_array(i)%data(j)%normalize( nn )
                !
            enddo
            !
        enddo
        !
    end subroutine normalize_DataGroupTxArray
    !
    !> No subroutine briefing
    !
    subroutine normalizeWith_DataGroupTxArray( norm, data_tx_array_in, data_tx_array_out )
        implicit none
        !
        integer, intent( in ) :: norm
        type( DataGroupTx_t ), dimension(:), intent( in ) :: data_tx_array_in
        type( DataGroupTx_t ), dimension(:), intent( inout ) :: data_tx_array_out
        !
        integer :: i, j
        !
        do i = 1, size( data_tx_array_in )
            !
            do j = 1, size( data_tx_array_in(i)%data )
                !
                data_tx_array_out(i)%data(j)%reals = ( data_tx_array_out(i)%data(j)%reals / data_tx_array_in(i)%data(j)%errors ** norm )
                !
                data_tx_array_out(i)%data(j)%imaginaries = ( data_tx_array_out(i)%data(j)%imaginaries / data_tx_array_in(i)%data(j)%errors ** norm )
                !
                data_tx_array_out(i)%data(j)%error_bar = .TRUE.
                !
            enddo
            !
        enddo
        !
    end subroutine normalizeWith_DataGroupTxArray
    !
    !> No subroutine briefing
    !
    subroutine setErrorBar_DataGroupTxArray( data_tx_array, error_bar )
        implicit none
        !
        type( DataGroupTx_t ), dimension(:), intent( inout ) :: data_tx_array
        logical, intent( in ) :: error_bar
        !
        integer :: i, j
        !
        do i = 1, size( data_tx_array )
            !
            do j = 1, size( data_tx_array(i)%data )
                !
                data_tx_array(i)%data(j)%error_bar = error_bar
                !
            enddo
            !
        enddo
        !
    end subroutine setErrorBar_DataGroupTxArray
    !
    !> No subroutine briefing
    !
    subroutine setComplex_DataGroupTxArray( data_tx_array, is_complex )
        implicit none
        !
        type( DataGroupTx_t ), dimension(:), intent( inout ) :: data_tx_array
        logical, intent( in ) :: is_complex
        !
        integer :: i, j
        !
        do i = 1, size( data_tx_array )
            !
            do j = 1, size( data_tx_array(i)%data )
                !
                data_tx_array(i)%data(j)%is_complex = is_complex
                !
            enddo
            !
        enddo
        !
    end subroutine setComplex_DataGroupTxArray
    !
    !> No function briefing
    !
    function countValues_DataGroupTxArray( data_tx_array ) result( counter )
        implicit none
        !
        type( DataGroupTx_t ), dimension(:), intent( in ) :: data_tx_array
        !
        integer :: counter
        !
        integer :: i
        !
        counter = 0
        !
        do i = 1, size( data_tx_array )
            !
            counter = counter + size( data_tx_array(i)%data ) * 2 * ( transmitters( data_tx_array(i)%i_tx )%n_pol ** 2 )
            !
        enddo
        !
    end function countValues_DataGroupTxArray
    !
    !> No function briefing
    !
    function count_DataGroupTxArray( data_tx_array ) result( counter )
        implicit none
        !
        type( DataGroupTx_t ), dimension(:), intent( in ) :: data_tx_array
        !
        integer :: counter
        !
        integer :: i, j, k
        !
        counter = 0
        !
        do i = 1, size( data_tx_array )
            !
            counter = counter + size( data_tx_array(i)%data )
            !
        enddo
        !
    end function count_DataGroupTxArray
    !
    !> Return a pointer, allowing directly modifications to a specific DataGroup at a given index
    !
    function getDataGroupByIndex( data_tx_array, i_dg ) result( data_group )
        implicit none
        !
        type( DataGroupTx_t ), dimension(:), target, intent( in ) :: data_tx_array
        integer, intent( in ) :: i_dg
        !
        type( DataGroup_t ), pointer :: data_group
        !
        integer :: i, j
        !
        do i = 1, size( data_tx_array )
            !
            do j = 1, size( data_tx_array(i)%data )
                !
                if( data_tx_array(i)%data(j)%i_dg == i_dg ) then
                    !
                    data_group => data_tx_array(i)%data(j)
                    !
                    return
                    !
                endif
                !
            enddo
            !
        enddo
        !
    end function getDataGroupByIndex
    !
    !> Dynamically add a new DataGroupTx to the array, always via reallocation.
    !
    subroutine update_DataGroupTxArray( data_tx_array, data_tx )
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
    end subroutine update_DataGroupTxArray
    !
    !> Deallocate all DataGroups of a particular DataGroupTxArray
    !
    subroutine deallocate_DataGroupTxArray( data_tx_array )
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: data_tx_array
        !
        integer :: i, n_dtx
        !
        !write( *, * ) "deallocate_DataGroupTxArray:", size( data_tx_array )
        !
        n_dtx = size( data_tx_array )
        !
        if( n_dtx == 1 ) then
            deallocate( data_tx_array(1)%data )
        elseif( n_dtx .GT. 0 ) then
            do i = n_dtx, 1, -(1)
                deallocate( data_tx_array(i)%data )
            enddo
        endif
        !
        deallocate( data_tx_array )
        !
    end subroutine deallocate_DataGroupTxArray
    !
    !> Write one DataGroupTxArray, with its proper Rx headers, 
    !> into to the file <file_name>
    !
    subroutine write_DataGroupTxArray( data_tx_array, file_name )
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: data_tx_array
        character(*), intent( in ) :: file_name
        !
        class( Receiver_t ), pointer :: receiver
        type( DataGroup_t ), pointer :: data_group
        !
        real( kind=prec ) :: SI_factor, r_error
        complex( kind=prec ) :: c_value
        integer :: receiver_type, i, j, ios, n_data, type_index
        !
        ! Verbose
        !write( *, * ) "     > Write Data to file: [", file_name, "]"
        !
        n_data = countData( data_tx_array )
        !
        receiver_type = 0
        !
        type_index = 0
        !
        open( unit = ioPlot, file = file_name//"_plot", action = "write", form = "formatted", iostat = ios )
        !
        write( ioPlot, * ) file_name
        !
        open( unit = ioPredData, file = file_name, action = "write", form = "formatted", iostat = ios )
        !
        if( ios == 0 ) then
            !
            do i = 1, n_data
                !
                data_group => getData( data_tx_array, i )
                !
                receiver => getReceiver( data_group%i_rx )
                !
                !> Write header if changed data type
                !> Increase type_index
                call writeHeader_DataGroupTxArray( receiver, receiver_type, type_index )
                !
                !write( *, * ) "receiver%units, units_in_file( type_index )%str: [", receiver%units, "],[", units_in_file( type_index )%str, "]"
                !
                SI_factor = ImpUnits( receiver%units, units_in_file( type_index )%str )
                !
                do j = 1, data_group%n_comp
                    !
                    if( conjugated_data ) then
                        !
                        c_value = cmplx( data_group%reals(j), -data_group%imaginaries(j), kind=prec )
                        !
                        if( data_group%error_bar ) then
                            r_error = -data_group%errors(j) * SI_factor
                        else
                            r_error = R_LARGE
                        endif
                        !
                    else
                        !
                        c_value = cmplx( data_group%reals(j), data_group%imaginaries(j), kind=prec )
                        !
                        if( data_group%error_bar ) then
                            r_error = data_group%errors(j) * SI_factor
                        else
                            r_error = R_LARGE
                        endif
                        !
                    endif
                    !
                    c_value = c_value * SI_factor
                    !
                    write( ioPredData, "(es12.6, 1X, A, 1X, f15.3, f15.3, f15.3, f15.3, f15.3, 1X, A, 1X, es16.6, es16.6, es16.6)" ) transmitters( data_group%i_tx )%period, trim( receiver%code ), R_ZERO, R_ZERO, receiver%location(1), receiver%location(2), receiver%location(3), trim( receiver%comp_names(j)%str ), real( c_value, kind=prec ), real( aimag( c_value ), kind=prec ), r_error
                    !
                    write( ioPlot, * ) sqrt( real( c_value, kind=prec )**2 + real( aimag( c_value ), kind=prec )**2 )
                    !
                enddo
                !
            enddo
            !
            close( ioPredData )
            !
            close( ioPlot )
            !
        else
            call errStop( "write_DataGroupTxArray > cant opening ["//file_name//"]!" )
        endif
        !
    end subroutine write_DataGroupTxArray
    !
    !> Call the print routine of each DataGroupTx in the array
    !
    subroutine print_DataGroupTxArray( data_tx_array, title )
        implicit none
        !
        type( DataGroupTx_t ), dimension(:), intent( in ) :: data_tx_array
        character(*), intent( in ) :: title
        !
        integer :: i
        !
        write( *, * ) "##############################"
        write( *, * ) title
        write( *, * ) "##############################"
        !
        do i = 1, size( data_tx_array )
            !
            call data_tx_array(i)%print
            !
        enddo
        !
    end subroutine print_DataGroupTxArray
    !
    !> Write a header into the DataGroupTxArray text file
    !
    subroutine writeHeader_DataGroupTxArray( receiver, receiver_type, type_index )
        implicit none
        !
        class( Receiver_t ), intent( in ) :: receiver
        integer, intent( inout ) :: receiver_type, type_index
        !
        if( receiver_type /= receiver%rx_type ) then
            !
            type_index = type_index + 1
            !
            select case( receiver%rx_type )
                !
                case( 1, 11, 12 )
                    !
                    write( ioPredData, * ) "# "//DATA_FILE_TITLE_MT
                    !
                    write( ioPredData, * ) "#   Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error"
                    !
                case( 2 )
                    write( ioPredData, * ) "#   Missing title for Full_Interstation_TF"
                case( 3 )
                    write( ioPredData, * ) "#   Missing title for Off_Diagonal_Rho_Phase"
                case( 4 )
                    write( ioPredData, * ) "#   Missing title for Phase_Tensor"
                case( 5 )
                    write( ioPredData, * ) "#   Missing title for Off_Diagonal_Impedance"
                case( 6, 7, 8, 9, 10 )
                    !
                    write( ioPredData, * ) "# "//DATA_FILE_TITLE_CSEM
                    !
                    write( ioPredData, * ) "#   Tx_Dipole Tx_Period(s) Tx_Moment(Am) Tx_Azi Tx_Dip Tx_X(m) Tx_Y(m) Tx_Z(m) Code X(m) Y(m) Z(m) Component Real Imag Error"
                    !
                case default
                    call errStop( "DataGroupTxArray > writeHeader_DataGroupTxArray > Unknown receiver type" )
                !
            end select
            !
            write( ioPredData, * ) ">  "//getStringReceiverType( receiver%rx_type )
            !
            if( conjugated_data ) then
                write( ioPredData, "( 18a )" ) ">  exp(+i\omega t)"
            else
                write( ioPredData, "( 18a )" ) ">  exp(-i\omega t)"
            endif
            !
            write( ioPredData, "( 50a )" ) ">  "//trim( units_in_file( type_index )%str )
            write( ioPredData, "( 10a )" ) ">     0.00"
            write( ioPredData, "( 20a )" ) ">     0.000    0.000"
            write( ioPredData, "( 1a, i8, i8 )" ) ">", size( transmitters ), size( receivers )
            !
            receiver_type = receiver%rx_type
            !
        endif
        !
    end subroutine writeHeader_DataGroupTxArray
    !
end module DataGroupTxArray
!