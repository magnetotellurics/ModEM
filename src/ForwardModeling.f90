!
!> Module with the sensitivity routines JMult, JMult_Tx, JMult_T, JMult_T_Tx 
!
module ForwardModeling
!
#ifdef MPI
    use MasterMPI
#else
    use GlobalVariables
#endif
    !
    !> Global FWD Routines
    public :: jobForwardModeling, runForwardModeling
    public :: solveTx, solveAll
    public :: createDistributeForwardSolver
    public :: writeDataGroupTxArray, writeAllESolutionHeader
    !
    private :: writeHeaderDataGroupTxArray
    !
contains
    !
    !> Calculate E_Solution (e_sol) for a single Transmitter
    !> ForwardSolver must be allocated
    subroutine solveTx( sigma, Tx )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma
        class( Transmitter_t ), pointer, intent( inout ) :: Tx
        !
        !> Set Transmitter's ForwardSolver Omega(Period) and Conductivity
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
                stop "Error: solveTx > Unclassified Transmitter"
            !
        end select
        !
        !> Build Source E according to source type
        call Tx%source%createE()
        !
        !> Solve e_sol for this Transmitter
        call Tx%solve()
        !
    end subroutine solveTx
    !
    !> Calculate ESolution for all transmitters.
    !> ForwardSolver must be allocated
    subroutine solveAll( sigma, e_all )
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
        !> Allocate e_all if present
        if( present( e_all ) ) then
            !
            if( allocated( e_all%e_sol ) ) deallocate( e_all%e_sol )
            allocate( e_all%e_sol( n_tx ) )
            !
        endif
        !
        !> Loop over all Transmitters
        do i_tx = 1, n_tx
            !
            !> Point to the current Transmitter
            Tx => getTransmitter( i_tx )
            !
            !> Solve i_tx transmitter
            call solveTx( sigma, Tx )
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
    end subroutine solveAll
    !
    !> Routine to run a full ForwardModeling job and deliver the result (PredictedData) in a text file
    subroutine jobForwardModeling()
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: all_predicted_data
        !
        class( ModelParameter_t ), allocatable :: sigma
        !
        ! Verbose
        write( *, * ) "     - Start jobForwardModeling"
        !
        !> Read Model File and instantiate global variables: main_grid, model_operator and Sigma0
        if( has_model_file ) then
            !
            call handleModelFile( sigma )
            !
        else
            stop "Error: jobForwardModeling > Missing Model file!"
        endif
        !
        !> Read Data File: instantiate Txs and Rxs and build the Data relation between them
        if( has_data_file ) then
            !
            call handleDataFile()
            !
        else
            stop "Error: jobForwardModeling > Missing Data file!"
        endif
        !
        all_predicted_data = all_measured_data
        !
#ifdef MPI
        !
        !> Send FWD basic components to all workers
        call broadcastBasicComponents()
        !
        !> Deallocate global FWD components (Not used by the Master process)
        deallocate( model_operator, main_grid )
        !
        !> Run masterForwardModelling to calculate predicted data
        call masterForwardModelling( sigma, all_predicted_data )
        !
        call broadcastFinish
        !
#else
        !
        call createDistributeForwardSolver()
        !
        call runForwardModeling( sigma, all_predicted_data )
        !
        if( has_e_solution_file ) call writeAllESolution( e_solution_file_name )
        !
#endif
        !
        !> Write Predicted Data, with its proper Rx headers, into to the file <predicted_data_file_name>
        call writeDataGroupTxArray( all_predicted_data, predicted_data_file_name )
        !
        !> Deallocate local variables
        call deallocateDataGroupTxArray( all_predicted_data )
        !
        deallocate( sigma )
        !
        ! Verbose
        write( *, * ) "     - Finish jobForwardModeling"
        !
    end subroutine jobForwardModeling
    !
    !> runForwardModeling:
    !>     Calculate ESolution for all transmitters and...
    !>     Calculate the predicted data for each transmitter-receiver pair.
    !>     ForwardSolver must be allocated
    subroutine runForwardModeling( sigma, all_predicted_data, e_all )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: all_predicted_data
        type( ESolMTx ), optional, intent( inout ) :: e_all
        !
        class( Transmitter_t ), pointer :: Tx
        class( Receiver_t ), pointer :: Rx
        type( DataGroup_t ) :: data_group
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
            call solveTx( sigma, Tx )
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
                call Rx%predictedData( Tx, data_group )
                !
                call all_predicted_data( i_tx )%setValues( data_group )
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
    subroutine createDistributeForwardSolver()
        implicit none
        !
        integer :: i_tx
        !
        if( allocated( forward_solver ) ) deallocate( forward_solver )
        !
        !> Instantiate the ForwardSolver - Specific type can be chosen via control file
        select case ( forward_solver_type )
            !
            case( FWD_IT_DC )
                allocate( forward_solver, source = ForwardSolverIT_DC_t( model_operator, QMR ) )
                !
            case default
                !
                write( *, * ) "Warning: createDistributeForwardSolver > Undefined forward_solver, using IT_DC"
                !
                allocate( forward_solver, source = ForwardSolverIT_DC_t( model_operator, QMR ) )
                !
        end select
        !
        !> Make all transmitters point to this ForwardSolver
        do i_tx = 1, size( transmitters )
            !
            transmitters( i_tx )%Tx%forward_solver => forward_solver
            !
        enddo
        !
    end subroutine createDistributeForwardSolver
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
    subroutine writeDataGroupTxArray( data_tx_array, file_name )
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: data_tx_array
        character(*), intent( in ) :: file_name
        !
        class( Transmitter_t ), pointer :: transmitter
        class( Receiver_t ), pointer :: receiver
        type( DataGroup_t ), pointer :: data_group
        !
        integer :: receiver_type, i, j, ios, n_data
        !
        ! Verbose
        !write( *, * ) "     > Write Data to file: [", file_name, "]"
        !
        n_data = countDataGroupTxArray( data_tx_array )
        !
        receiver_type = 0
        !
        open( unit = ioPredData, file = file_name, action = "write", form = "formatted", iostat = ios )
        !
        if( ios == 0 ) then
            !
            do i = 1, n_data
                !
                data_group => getDataGroupByIndex( data_tx_array, i )
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
                            write( ioPredData, "(es12.6, 1X, A, 1X, f15.3, f15.3, f15.3, f15.3, f15.3, 1X, A, 1X, es16.6, es16.6, es16.6)" ) transmitter%period, trim(receiver%code), R_ZERO, R_ZERO, receiver%location(1), receiver%location(2), receiver%location(3), trim( receiver%comp_names(j)%str ), data_group%reals(j), data_group%imaginaries(j), data_group%errors(j)
                            !
                        class is( TransmitterCSEM_t )
                            !
                            write( ioPredData, "(A, 1X, es12.6, f15.3, f15.3, f15.3, f15.3, f15.3, f15.3, 1X, A, 1X, f15.3, f15.3, f15.3, 1X, A, 1X, es16.6, es16.6, es16.6)" ) trim(transmitter%dipole), transmitter%period, transmitter%moment, transmitter%azimuth, transmitter%dip, transmitter%location(1), transmitter%location(2), transmitter%location(3), trim(receiver%code), receiver%location(1), receiver%location(2), receiver%location(3), trim( receiver%comp_names(j)%str ), data_group%reals(j), data_group%imaginaries(j), data_group%errors(j)
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