!
!> Module with the sensitivity routines JMult, JMult_Tx, JMult_T, JMult_T_Tx 
!
module ForwardModeling
    !
    use Constants
    !
    use FileUnits
    !
    use cVector3D_SG
    !
    use GlobalVariables
    !
    use SourceMT_1D
    use SourceMT_2D
    use SourceCSEM_Dipole1D
    use SourceInteriorForce
    !
    use ReceiverFullImpedance
    use ReceiverFullVerticalMagnetic
    use ReceiverOffDiagonalImpedance
    use ReceiverSingleField
    use ReceiverArray
    !
    use TransmitterMT
    use TransmitterCSEM
    use TransmitterArray
    !
    !> Global FWD Routines
    public :: runEMSolve, runForwardModeling, writeAllESolutionHeader
    public :: writeDataGroupTxArray
    !
    private :: writeHeaderDataGroupTxArray
    !
contains
    !
    !> runEMSolve: Calculate ESolution for all transmitters.
    !
    !> Obs.: Require the previous definition of the global ForwardSolver (createForwardSolver())
    !
    subroutine runEMSolve( sigma, e_all )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma
        type( ESolMTx ), optional, intent( inout ) :: e_all
        !
        class( Transmitter_t ), pointer :: Tx
        integer :: i_tx, n_tx
        !
        ! Verbose
        write( *, * ) "          - Start EM Solve"
        !
        !>
        n_tx = size( transmitters )
        !
        !> Set e_all if present
        if( present( e_all ) ) then
            !
            !> SET e_all
            if( allocated( e_all%e_sol ) ) deallocate( e_all%e_sol )
            allocate( e_all%e_sol( n_tx ) )
            !
        endif
        !
        !> Loop over all Transmitters
        do i_tx = 1, n_tx
            !
            !> Pointer to the Transmitter
            Tx => getTransmitter( i_tx )
            !
            !> Set Transmitter's ForwardSolver Period) and Conductivity
            call Tx%forward_solver%setFrequency( sigma, Tx%period )
            !
            !> Instantiate Transmitter's Source - According to transmitter type or chosen via control file
            select type( Tx )
                !
                class is( TransmitterMT_t )
                    !
                    call Tx%setSource( SourceMT_1D_t( model_operator, sigma, Tx%period ) )
                    !
                class is( TransmitterCSEM_t )
                    !
                    call Tx%setSource( SourceCSEM_Dipole1D_t( model_operator, sigma, Tx%period, Tx%location, Tx%dip, Tx%azimuth, Tx%moment ) )
                    !
                class default
                    stop "Error: runEMSolve > Unclassified Transmitter"
                !
            end select
            !
            !> Build Source E according to source type
            call Tx%source%createE()
            !
            !> Solve e_sol for this Transmitter
            call Tx%solve()
            !
            !> Save e_sol in e_all
            if( present( e_all ) ) then
                !
                e_all%e_sol( i_tx )%pol = Tx%e_sol
                !
            endif
            !
        enddo
        !
        ! Verbose
        write( *, * ) "          - Finish EM Solve"
        !
    end subroutine runEMSolve
    !
    !> runForwardModeling:
    !>     Calculate ESolution for all transmitters and...
    !>     Calculate the predicted data for each transmitter-receiver pair.
    !> Obs.: Require the previous definition of the global ForwardSolver (createForwardSolver())
    !
    subroutine runForwardModeling( sigma, predicted_data, e_all )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: predicted_data
        type( ESolMTx ), optional, intent( inout ) :: e_all
        !
        class( Transmitter_t ), pointer :: Tx
        class( Receiver_t ), pointer :: Rx
        integer :: i_tx, n_tx, i_rx
        !
        ! Verbose
        write( *, * ) "          - Start Forward Modeling"
        !
        !>
        n_tx = size( transmitters )
        !
        !> Set e_all if present
        if( present( e_all ) ) then
            !
            !> SET e_all
            if( allocated( e_all%e_sol ) ) deallocate( e_all%e_sol )
            allocate( e_all%e_sol( n_tx ) )
            !
        endif
        !
        !> Loop over all Transmitters
        do i_tx = 1, n_tx
            !
            !> Pointer to the Transmitter
            Tx => getTransmitter( i_tx )
            !
            !> Set Transmitter's ForwardSolver Period) and Conductivity
            call Tx%forward_solver%setFrequency( sigma, Tx%period )
            !
            !> Instantiate Transmitter's Source - According to transmitter type or chosen via control file
            select type( Tx )
                !
                class is( TransmitterMT_t )
                    !
                    call Tx%setSource( SourceMT_1D_t( model_operator, sigma, Tx%period ) )
                    !
                class is( TransmitterCSEM_t )
                    !
                    call Tx%setSource( SourceCSEM_Dipole1D_t( model_operator, sigma, Tx%period, Tx%location, Tx%dip, Tx%azimuth, Tx%moment ) )
                    !
                class default
                    stop "Error: ForwardModeling > Unclassified Transmitter"
                !
            end select
            !
            !> Build Source E according to source type
            call Tx%source%createE()
            !
            !> Solve e_sol for this Transmitter
            call Tx%solve()
            !
            !> Save e_sol in e_all
            if( present( e_all ) ) then
                !
                e_all%e_sol( i_tx )%pol = Tx%e_sol
                !
            endif
            !
            ! Verbose
            write( *, * ) "               - Calculate Predicted Data"
            !
            !> Loop for each Receiver related to this Transmitter
            do i_rx = 1, size( Tx%receiver_indexes )
                !
                !> Pointer to the Tx Receiver
                Rx => getReceiver( Tx%receiver_indexes( i_rx ) )
                !
                call Rx%predictedData( Tx )
                !
                call predicted_data( i_tx )%set( Rx%data_group )
                !
            enddo
            !
        enddo
        !
        ! Verbose
        write( *, * ) "          - Finish Forward Modeling"
        !
    end subroutine runForwardModeling
	!
    !> No subroutine briefing
    subroutine writeAllESolution( file_name )
        implicit none
        !
        character(*), intent( in ) :: file_name
        !
        integer :: i_tx, n_tx
        !
        if( .NOT. allocated( transmitters ) .OR. size( transmitters ) .LT. 1 ) then
            stop "Error: writeAllESolution > Theres no transmitters to write!"
        endif
        !
        if( .NOT. allocated( transmitters(1)%Tx%e_sol ) .OR. size( transmitters(1)%Tx%e_sol ) .LT. 1 ) then
            stop "Error: writeAllESolution > Theres no e_sol to write!"
        endif
        !
        ! Verbose
        write( *, * ) "          > Write All ESolutions to file: [", file_name, "]"
        !
        !>
        n_tx = size( transmitters )
        !
        !> Write the first header of the ESolution binary file, according to the first transmitter
        call writeAllESolutionHeader( n_tx, transmitters(1)%Tx%n_pol, file_name )
        !
        !> Loop over all Transmitters
        do i_tx = 1, n_tx
            !
            call transmitters( i_tx )%Tx%writeESolution()
            !
        enddo
        !
    end subroutine writeAllESolution
    !
    !> No subroutine briefing
    subroutine writeAllESolutionHeader( n_tx, nMode, file_name )
        implicit none
        !
        integer, intent( in ) :: n_tx, nMode
        character(*), intent( in ) :: file_name
        integer :: ios
        character(len=20) :: version
        !
        version = "Modem-OO"
        !
        open( ioESolution, file = file_name, action = "write", form = "unformatted", iostat = ios)
        !
        if( ios == 0 ) then
            !
            write( ioESolution ) version, n_tx, nMode, &
            main_grid%nx, main_grid%ny, main_grid%nz, main_grid%nzAir, &
            main_grid%ox, main_grid%oy, main_grid%oz, main_grid%rotdeg
            !
            write( ioESolution ) main_grid%dx
            write( ioESolution ) main_grid%dy
            write( ioESolution ) main_grid%dz
            !
            close( ioESolution )
            !
        else
            !
            write( *, * ) "Error opening file in writeAllESolutionHeader [", file_name, "]!"
            stop
            !
        endif
        !
        !
    end subroutine writeAllESolutionHeader
    !
    !> No subroutine briefing
    subroutine writeDataGroupTxArray( tx_data_array_in, file_name )
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: tx_data_array_in
        character(*), intent( in ) :: file_name
        !
        type( DataGroup_t ), allocatable, dimension(:) :: to_write_data
        !
        class( Transmitter_t ), pointer :: transmitter
        class( Receiver_t ), pointer :: receiver
        type( DataGroup_t ), pointer :: data_group
        !
        integer :: receiver_type, i, j, ios
        !
        ! Verbose
        !write( *, * ) "     > Write Data to file: [", file_name, "]"
        !
        !> Put the data in the same input format
        allocate( to_write_data, source = measured_data )
        !
        do i = 1, size( tx_data_array_in )
            do j = 1, size( tx_data_array_in(i)%data )
                call setDataGroup( to_write_data, tx_data_array_in(i)%data(j) )
            enddo
        enddo
        !
        receiver_type = 0
        !
        open( unit = ioPredData, file = file_name, action = "write", form = "formatted", iostat = ios )
        !
        if( ios == 0 ) then
            !
            do i = 1, size( to_write_data )
                !
                data_group => getDataGroupByIndex( to_write_data, i )
                !
                receiver => getReceiver( data_group%i_rx )
                !
                call writeHeaderDataGroupTxArray( receiver, receiver_type )
                !
                transmitter => getTransmitter( data_group%i_tx )
                !
                do j = 1, data_group%n_comp
                    !
                    select type( transmitter )
                        !
                        class is( TransmitterMT_t )
                            !
                            write( ioPredData, "(es12.6, 1X, A, 1X, f15.3, f15.3, f15.3, f15.3, f15.3, 1X, A, 1X, es16.6, es16.6, es16.6)" ) transmitter%period, trim(receiver%code), R_ZERO, R_ZERO, receiver%location(1), receiver%location(2), receiver%location(3), trim(data_group%components(j)%str), data_group%reals(j), data_group%imaginaries(j), data_group%errors(j)
                            !
                        class is( TransmitterCSEM_t )
                            !
                            write( ioPredData, "(A, 1X, es12.6, f15.3, f15.3, f15.3, f15.3, f15.3, f15.3, 1X, A, 1X, f15.3, f15.3, f15.3, 1X, A, 1X, es16.6, es16.6, es16.6)" ) trim(transmitter%dipole), transmitter%period, transmitter%moment, transmitter%azimuth, transmitter%dip, transmitter%location(1), transmitter%location(2), transmitter%location(3), trim(receiver%code), receiver%location(1), receiver%location(2), receiver%location(3), trim(data_group%components(j)%str), data_group%reals(j), data_group%imaginaries(j), data_group%errors(j)
                            !
                        class default
                            stop "Error: writeDataGroupTxArray: Unclassified data_group!"
                        !
                    end select
                    !
                enddo
                !
            enddo
            !
            deallocate( to_write_data )
            !
            close( ioPredData )
            !
        else
            write( *, * ) "Error opening [", file_name, "] in writeDataGroupTxArray!"
            stop
        endif
        !
    end subroutine writeDataGroupTxArray
    !
    !> No subroutine briefing
    subroutine writeHeaderDataGroupTxArray( receiver, receiver_type )
        implicit none
        !
        class( Receiver_t ), intent( in ) :: receiver
        !
        integer, intent( inout ) :: receiver_type
        !
        if( receiver_type /= receiver%rx_type ) then
            !
            select case( receiver%rx_type )
                !
                case( 1, 11, 12 )
                    !
                    write( ioPredData, "(4A, 40A)" ) "#    ", DATA_FILE_TITLE_MT
                    !
                    write( ioPredData, "(74A)" ) "#    Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error"
                    !
                case( 2 )
                    write( ioPredData, "(60A)" ) "#    Missing title for Full_Interstation_TF"
                case( 3 )
                    write( ioPredData, "(60A)" ) "#    Missing title for Off_Diagonal_Rho_Phase"
                case( 4 )
                    write( ioPredData, "(60A)" ) "#    Missing title for Phase_Tensor"
                case( 5 )
                    write( ioPredData, "(60A)" ) "#    Missing title for Off_Diagonal_Impedance"
                case( 6, 7, 8, 9, 10 )
                    !
                    write( ioPredData, "(4A, 40A)" ) "#    ", DATA_FILE_TITLE_CSEM
                    !
                    write( ioPredData, "(125A)" ) "#    Tx_Dipole Tx_Period(s) Tx_Moment(Am) Tx_Azi Tx_Dip Tx_X(m) Tx_Y(m) Tx_Z(m) Code X(m) Y(m) Z(m) Component Real Imag Error"
                    !
                case default
                    write( *, * ) "Unknown receiver type :[", receiver%rx_type, "]"
                    stop "Error: test_FWD.f90: writeHeaderDataGroupTxArray()"
                !
            end select
            !
            write( ioPredData, "(4A, 100A)" ) ">    ", getStringReceiverType( receiver%rx_type )
            write( ioPredData, "(4A, 100A)" ) ">    ", "exp(-i\omega t)"
            write( ioPredData, "(4A, 100A)" ) ">    ", "[V/m]/[T]"
            write( ioPredData, "(7A, 100A)" ) ">        ", "0.00"
            write( ioPredData, "(7A, 100A)" ) ">        ", "0.000    0.000"
            write( ioPredData, "(A3, i8, i8)" ) ">        ", size( transmitters ), size( receivers )
            !
            receiver_type = receiver%rx_type
            !
        endif
        !
    end subroutine writeHeaderDataGroupTxArray
    !
end module ForwardModeling
!