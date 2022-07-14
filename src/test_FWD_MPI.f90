program ModEM
    !
    use DeclarationMPI
    !
    ! Program control variables
    character(:), allocatable :: control_file_name, model_file_name, data_file_name, modem_job
    logical                   :: set_data_groups, has_control_file, has_model_file, has_data_file, verbosis
    !
    class( ModelOperator_t ), allocatable  :: model_operator
    !
    real( kind=prec ) :: t_start, t_finish
    !
    !
    call cpu_time( t_start )
    !
    call MPI_Init( ierr )
    !
    main_comm = MPI_COMM_WORLD
    !
    ! SET mpi_size WITH n FROM MPIRUN FOR MPI_COMM_WORLD
    call MPI_Comm_size( main_comm, mpi_size, ierr )
    !
    if( mpi_size < 2 ) then
        write( *, * ) "Error: Minimum of two processes required!"
        call MPI_Finalize( ierr )
        stop
    end if 
    !
    ! SET mpi_rank WITH PROCESS ID FOR MPI_COMM_WORLD
    call MPI_Comm_rank( main_comm, mpi_rank, ierr )
    !
    ! SPLIT MPI_COMM_WORLD into shared subcommunicator: child_comm
    call MPI_Comm_split_type( main_comm, MPI_COMM_TYPE_SHARED, mpi_rank, MPI_INFO_NULL, child_comm, ierr )
    !
    call MPI_Get_processor_name( node_name, nodestringlen, ierr )
    !
    call MPI_Comm_size( child_comm, node_size, ierr )
    !
    call MPI_Comm_rank( child_comm, node_rank, ierr )
    !
    !write( *, * ) "MPI Rank ", mpi_rank," in COMM_WORLD [", mpi_size, "] is ", node_rank, &
    !              " in SHARED_COMM [", node_size, "] on Node: ", node_name( 1 : nodestringlen )
    !
    ! MPI MASTER PROCESS
    !
    if ( mpi_rank == 0 ) then
        !
        modem_job = "unknown"
        !
        call setupDefaultParameters()
        !
        ! Validate arguments and control variables
        call handleArguments()
        !
        write( *, * )
        write( *, * ) "Start ModEM-OO."
        write( *, * )
        !
        ! Check parameters at the control file
        if( has_control_file ) call handleControlFile()
        !
        ! Execute the modem_job
        call handleJob()
        !
        ! Deallocate remaining main program memory
        call garbageCollector()
        !
        call MPI_Finalize( ierr )
        !
        call cpu_time( t_finish )
        !
        write( *, * )
        write( *, "(A18, F6.3, A3)" ) "Finish ModEM-OO (", t_finish - t_start, "s)"
        write( *, * )
    !
    ! MPI WORKER PROCESS
    !
    else
        !
        call MPI_Win_allocate_shared( shared_window_size, shared_disp_unit, MPI_INFO_NULL, child_comm, shared_c_ptr, shared_window, ierr )
        !
        if( ierr == MPI_SUCCESS ) then
            !
            do while ( job_master .ne. job_finish )
                !
                !write( *, * ) "WORKER: ", mpi_rank, " WAITING MASTER"
                !
                call receiveFrom( master_id )
                !
                job_master = trim( fwd_info%job_name )
                !
                select case ( job_master )
                    !
                    case ( "SHARE_MEMORY" )
                        !
                        call workerQuerySharedMemory()
                        !
                        call MPI_Win_fence( 0, shared_window, ierr )
                        !
                        call MPI_Win_free( shared_window, ierr )
                        !
                    case ( "JOB_FORWARD" )
                        !
                        call workerForwardModelling()
                        !
                    case ( "JOB_INVERSE" )
                        !
                        call workerInversion()
                        !
                end select
                !
            enddo
            !
            deallocate( model_operator )
            deallocate( model_parameter )
            deallocate( main_grid )
            !
            call deallocateReceiverArray()
            !
            call cpu_time( t_finish )
            !
            write( *, * ) "Worker", mpi_rank, " finished ( ", t_finish - t_start, "s )"
            !
            call MPI_Finalize( ierr )
            !
        else
            !
            write( *, * ) "MPI Win_allocate_shared fails for process:", mpi_rank
            !
            stop
        endif
        !
     endif
    !
contains
    !
    subroutine masterInversion()
        implicit none
        !
        ! Verbose
        write( *, * ) "     - Start Inversion"
        !
        set_data_groups = .TRUE.
        !
        call masterForwardModelling()
        !
    end subroutine masterInversion
    !
    !
    subroutine masterForwardModelling()
        implicit none
        !
        ! Local variables
        integer :: i, worker_rank, tx_received, tx_index
        !
        ! Verbose
        write( *, * ) "     - Start Forward Modeling"
        !
        worker_rank = 1
        tx_received = 0
        tx_index = 0
        !
        ! Reads Model File: instantiates Grid, ModelOperator and ModelParameter
        if( .NOT. has_model_file ) then 
            stop "Error: Missing Model file!"
        else
            call handleModelFile()
        endif
        !
        ! Reads Data File: instantiates and builds the Data relation between Txs and Rxs
        if( .NOT. has_data_file ) then 
            stop "Error: Missing Data file!"
        else
            call handleDataFile()
        endif
        !
        write( *, * ) "     - MPI Shared memory"
        !
        ! Share memory with all workers
        call masterExposeSharedMemory()
        !
        do i = 1, ( mpi_size - 1 )
            !
            fwd_info%job_name = job_share_memory
            !
            call sendTo(i)
            !
        enddo
        !
        call MPI_Win_fence( 0, shared_window, ierr )
        !
        call MPI_Win_free( shared_window, ierr )
        !
        deallocate( model_operator )
        deallocate( model_parameter )
        deallocate( main_grid )
        !
        ! Send 1 transmitter to first np workers
        do while( worker_rank <= ( mpi_size - 1 ) )
            !
            tx_index = tx_index + 1
            !
            fwd_info%job_name    = job_forward
            fwd_info%tx_index    = tx_index
            fwd_info%worker_rank = worker_rank
            !
            call sendTo( worker_rank )
            !
            worker_rank = worker_rank + 1
            !
        end do
        !
        ! Send 1 transmitter to first available worker
        do while( tx_index < size( transmitters ) )
            !
            !write( *, * ) "THERES", ( size( transmitters ) - tx_index ), " TX LEFT!"
            !
            call receiveFromAny()
            !
            call receiveData()
            !
            tx_received = tx_received + 1
            !
            tx_index = tx_index + 1
            !
            fwd_info%job_name = job_forward
            fwd_info%tx_index = tx_index
            !
            call sendTo( fwd_info%worker_rank )
            !
        end do
        !
        ! Receives job_done from each finished worker
        do while( tx_received < size( transmitters ) )
            !
            !write( *, * ) "MASTER WAITING ANY WORKER TO FINISH"
            !
            call receiveFromAny()
            !
            call receiveData()
            !
            tx_received = tx_received + 1
            !
            fwd_info%job_name = job_finish
            !
            call sendTo( fwd_info%worker_rank )
            !
        enddo
        !
        ! Write data_groups into <predicted_data_file_name>
        call writeDataGroupArray( data_groups )
        !
        call deallocateTransmitterArray()
        !
    end subroutine masterForwardModelling
    !
    subroutine masterExposeSharedMemory()
        implicit none
        !
        call allocateSharedBuffer()
        !
        shared_window_size = shared_buffer_size
        shared_disp_unit = 1
        !
        call MPI_Win_allocate_shared( shared_window_size, shared_disp_unit, MPI_INFO_NULL, child_comm, shared_c_ptr, shared_window, ierr )
        !
        if( ierr /= MPI_SUCCESS ) then
             write( *, * ) "Error: MPI Win_allocate_shared fails on master:", ierr
             stop
        endif
        !
        call c_f_pointer( shared_c_ptr, shared_buffer, (/shared_window_size/) )
        !
        call packSharedBuffer()
        !
    end subroutine masterExposeSharedMemory
    !
    subroutine workerQuerySharedMemory()
        implicit none
        !
        call MPI_Win_shared_query( shared_window, master_id, shared_window_size, shared_disp_unit, shared_c_ptr, ierr )
        !
        if( ierr == MPI_SUCCESS ) then
            !
            write( *, * ) "          Worker", mpi_rank, " query ", shared_window_size, " bytes."
        else
            write( *, * ) "Error: MPI Win_shared_query fails on worker:", mpi_rank, ierr
            stop
        endif
        !
        call c_f_pointer( shared_c_ptr, shared_buffer, (/shared_window_size/) )
        !
        call unpackSharedBuffer( int( shared_window_size ) )
        !
    end subroutine workerQuerySharedMemory
    !
    subroutine workerInversion()
        implicit none
        !
        call workerForwardModelling()
        !
    end subroutine workerInversion
    !
    subroutine workerForwardModelling()
        implicit none
        !
        ! Use save ????
        class( ForwardSolver_t ), allocatable, target, save :: forward_solver
        !
        ! Temporary alias pointers
        class( Transmitter_t ), pointer :: Tx
        class( Receiver_t ), pointer    :: Rx
        !
        integer :: iRx
        type( TAirLayers ) :: air_layer
        type( DataGroup_t ), allocatable, dimension(:) :: tx_data_groups
        !
        select type( main_grid )
            !
            class is( Grid3D_SG_t )
                !
                call main_grid%SetupAirLayers( air_layer, model_method, model_n_air_layer, model_max_height )
                call main_grid%UpdateAirLayers( air_layer%nz, air_layer%dz )
                !
                allocate( model_operator, source = ModelOperator_MF_t( main_grid ) )
                !
                call model_parameter%setMetric( model_operator%metric )
                !
                call model_operator%SetEquations()
                !
                call model_operator%SetCond( model_parameter )
                !
            class default
                stop "Error: Unclassified main_grid"
            !
        end select
        !
        ! Instantiate the ForwardSolver - Specific type can be chosen via control file
        select case ( forward_solver_type )
            !
            case( FWD_IT_DC )
                allocate( forward_solver, source = ForwardSolverIT_DC_t( model_operator, QMR ) )
                !
            case default
                allocate( forward_solver, source = ForwardSolverIT_DC_t( model_operator, QMR ) )
            !
        end select
        !
        ! Point to the current Transmitter
        Tx => getTransmitter( fwd_info%tx_index )
        !
        ! Point Transmitter's ForwardSolver
        Tx%forward_solver => forward_solver
        !
        ! Set Transmitter's ForwardSolver Omega(Period) and Conductivity
        call Tx%forward_solver%setFrequency( model_parameter, Tx%period )
        !
        ! Instantiate Transmitter's Source - According to tx type or via control file
        select type( Tx )
            !
            class is( TransmitterMT_t )
                !
                allocate( Tx%source, source = SourceMT_1D_t( model_operator, model_parameter, Tx%period ) )
                !
            class is( TransmitterCSEM_t )
                !
                allocate( Tx%source, source = SourceCSEM_Dipole1D_t( model_operator, model_parameter, Tx%period, Tx%location, Tx%dip, Tx%azimuth, Tx%moment ) )
                !
        end select
        !
        ! Solve Forward Modeling for the Transmitter
        call Tx%solveFWD()
        !
        ! Loop for each Receiver related to the Transmitter
        do iRx = 1, size( Tx%receiver_indexes )
            !
            ! Point to the current Receiver
            Rx => getReceiver( Tx%receiver_indexes( iRx ) )
            !
            ! Calculate and store predicted data in the Receiver
            call Rx%predictedData( Tx )
            !
            call updateDataGroupArray( tx_data_groups, Rx%predicted_data )
            !
        enddo
        !
        ! Clears the memory used by the current Transmitter (Mainly ESolution cVector)
        deallocate( Tx )
        !
        !write( *, * ) "WORKER ", mpi_rank, "FINISHES FWD FOR TX ", fwd_info%tx_index!, size( tx_data_groups )
        !
        ! SEND JOB DONE TO MASTER
        fwd_info%job_name    = job_fwd_done
        fwd_info%worker_rank = mpi_rank
        !
        call allocateDataBuffer( tx_data_groups )
        !
        fwd_info%n_data      = size( tx_data_groups )
        fwd_info%data_size   = predicted_data_buffer_size
        !
        call sendTo( master_id )
        !
        call sendData( tx_data_groups )
        !
        deallocate( tx_data_groups )
        !
    end subroutine workerForwardModelling
    !
    subroutine handleJob()
        implicit none
        !
        select case ( modem_job )
            !
            case ( "forward" )
                !
                call masterForwardModelling()
                !
                ! Deallocate MPI communication buffers
                deallocate( fwd_info_buffer, predicted_data_buffer )
                !
                deallocate( data_groups )
                !
                call deallocateReceiverArray()
                !
            case ( "inversion" )
                !
                call masterInversion()
                !
            case default
                !
                write( *, * ) "Error: Unknown job: [", modem_job, "]"
                call printHelp()
                stop
                !
        end select
        !
    end subroutine handleJob
    !
    subroutine handleControlFile()
        implicit none
        !
        type( ModEMControlFile_t ) :: control_file
        !
        write( *, * ) "     < Control File: [", control_file_name, "]"
        !
        ! Instantiate the ControlFile object
        ! Reads control file and sets the options in the Constants module
        control_file = ModEMControlFile_t( ioStartup, control_file_name )
        !
    end subroutine handleControlFile
    !
    !
    subroutine handleModelFile()
        implicit none
        !
        type( ModelReader_Weerachai_t ) :: model_reader
        type( TAirLayers )              :: air_layer
        !
        !
        write( *, * ) "     < Model File: [", model_file_name, "]"
        !
        ! Read Grid and ModelParameter with ModelReader_Weerachai
        call model_reader%Read( model_file_name, main_grid, model_parameter ) 
        !
        ! Instantiate the ModelOperator object
        select type( main_grid )
            !
            class is( Grid3D_SG_t )
                !
                call main_grid%SetupAirLayers( air_layer, model_method, model_n_air_layer, model_max_height )
                !
                call main_grid%UpdateAirLayers( air_layer%nz, air_layer%dz )
                !
                write( *, * ) "          Air layers from the method: ", air_layer%method, "."
                !
                write( *, * ) "          Top of the air layers is at ", sum(air_layer%Dz)/1000, " km."
                !
                allocate( model_operator, source = ModelOperator_MF_t( main_grid ) )
                !
                call model_parameter%setMetric( model_operator%metric )
                !
                call model_operator%SetEquations()
                !
                call model_operator%SetCond( model_parameter )
                !
            class default
                stop "Error: Unclassified main_grid"
            !
        end select
        !
    end subroutine handleModelFile
    !
    !
    subroutine handleDataFile()
        implicit none
        !
        integer :: irx, nrx
        class( Receiver_t ), pointer :: Rx
        !
        ! Local object to dealt data, self-destructs at the end of the subroutine
        type( DataFileStandard_t ) :: data_file_standard
        !
        write( *, * ) "     < Data File: [", data_file_name, "]"
        !
        data_file_standard = DataFileStandard_t( ioStartup, data_file_name )
        !
        nrx = size( receivers )
        !
        if( nrx == data_file_standard%nRx ) then
            !
            write( *, * ) "          Checked ", nrx, " Receivers."
            !
        else
            !
            write( *, * ) "Number of Rx mismatched from Header :[", nrx, " and ", data_file_standard%nRx, "]"
            stop "Error: DataManager.f08: DataManager_ctor()"
            !
        endif
        !
        if( size( transmitters ) == data_file_standard%nTx ) then
            !
            if( ( mpi_size - 1 ) > size( transmitters ) ) then
                write( *, * ) "There are more MPI worker processes than necessary!!!"
                write( *, * ) "     ", size( transmitters ), " Transmitters"
                write( *, * ) "     ", ( mpi_size - 1 ), " MPI worker processes"
                write( *, * ) "Use 'mpirun -np ", size( transmitters ) + 1, "'"
                !
                call MPI_Abort( main_comm, ierr )
                !
                stop
            endif
            !
            call printTransmitterArray()
            !
        else
             !
             write( *, * ) "Number of Tx mismatched from Header :[", size( transmitters ), " and ", data_file_standard%nTx, "]"
             stop "Error: DataManager.f08: DataManager_ctor()"
             !
        endif
        !
        write( *, * ) "          Checked ", size( data_groups ), " DataGroups."
        !
        write( *, * ) "     - Create Rx evaluation vectors"
        !
        do irx = 1, nrx
            !
            Rx => getReceiver( irx )
            !
            call Rx%evaluationFunction( model_operator )
            !
        enddo
        !
    end subroutine handleDataFile
    !
    subroutine handleArguments()
        implicit none
        !
        character( len=200 ) :: argument
        integer              :: argument_index
        !
        if ( command_argument_count() == 0 ) then
            !
            call printHelp()
            !
        else
            !
            argument_index = 1
            !
            do while(argument_index <= command_argument_count()) 
                 !
                 call get_command_argument( argument_index, argument )
                 !
                 select case ( argument )
                      !
                      case ( "-c", "--control" )
                         !
                         call get_command_argument( argument_index + 1, argument )
                         control_file_name = trim( argument )
                         !
                         if ( len( control_file_name ) > 0 ) has_control_file = .TRUE.
                         !
                         argument_index = argument_index + 2
                         !
                      case ( "-d", "--data" )
                         !
                         call get_command_argument( argument_index + 1, argument )
                         data_file_name = trim( argument )
                         !
                         if ( len( data_file_name ) > 0 ) has_data_file = .TRUE.
                         !
                         argument_index = argument_index + 2
                         !
                      case ( "-f", "--forward" )
                         !
                         modem_job = "forward"
                         !
                         argument_index = argument_index + 1
                         !
                      case ( "-i", "--inversion" )
                         !
                         modem_job = "inversion"
                         !
                         argument_index = argument_index + 1
                         !
                      case ( "-m", "--model" )
                         !
                         call get_command_argument( argument_index + 1, argument )
                         model_file_name = trim( argument )
                         !
                         if ( len( model_file_name ) > 0 ) has_model_file = .TRUE.
                         !
                         argument_index = argument_index + 2
                         !
                      case ( "-pd", "--predicted" )
                         !
                         call get_command_argument( argument_index + 1, argument )
                         predicted_data_file_name = trim( argument )
                         !
                         argument_index = argument_index + 2
                         !
                      case ( "-es", "--esolution" )
                         !
                         call get_command_argument( argument_index + 1, argument )
                         e_solution_file_name = trim( argument )
                         !
                         argument_index = argument_index + 2
                         !
                      case ( "-v", "--version" )
                         !
                         write( *, * ) "    + ModEM-OO version 1.0.0"
                         stop
                         !
                      case ( "-h", "--help" )
                         !
                         call printHelp()
                         !
                      case ( "--verbose" )
                         !
                         call printHelp()
                         !
                      case default
                         !
                         write( *, * ) "Error: Unknown Argument: [", argument, "]"
                         call printHelp()
                         stop
                         !
                 end select
                 !
                 argument = ""
                 !
            end do
            !
        end if
        !
    end subroutine handleArguments
    !
    subroutine setupDefaultParameters()
        implicit none
        !
        ! I|O
        predicted_data_file_name = "predicted_data.dat"
        e_solution_file_name     = "esolution.bin"
        has_control_file         = .FALSE.
        has_model_file           = .FALSE.
        has_data_file            = .FALSE.
        verbosis                 = .FALSE.
        !
        ! Solvers
        QMR_iters = 40
        BCG_iters = 80
        max_divcor = 20
        max_divcor_iters = 100
        tolerance_divcor = 1E-5
        tolerance_qmr = 1E-7
        forward_solver_type = FWD_IT_DC
        !
        ! Source
        source_type = SRC_MT_1D
        get_1D_from = "Geometric_mean"
        !
        ! Model
        model_method      = MM_METHOD_FIXED_H
        model_n_air_layer = 10
        model_max_height  = 200.0
        !
    end subroutine setupDefaultParameters
    !
    subroutine garbageCollector()
        implicit none
        !
        if( allocated( forward_solver_type ) ) deallocate( forward_solver_type )
        if( allocated( source_type ) ) deallocate( source_type )
        if( allocated( model_method ) ) deallocate( model_method )
        if( allocated( get_1D_from ) ) deallocate( get_1D_from )
        if( allocated( predicted_data_file_name ) ) deallocate( predicted_data_file_name )
        if( allocated( e_solution_file_name ) ) deallocate( e_solution_file_name )
        !
        if( allocated( control_file_name ) ) deallocate( control_file_name )
        if( allocated( model_file_name ) ) deallocate( model_file_name )
        if( allocated( data_file_name ) ) deallocate( data_file_name )
        if( allocated( modem_job ) ) deallocate( modem_job )
        !
    end subroutine garbageCollector
    !
    subroutine writeEsolutionHeader( nTx, nMode )
        implicit none
        !
        ! implement separated routine
        integer, intent( in ) :: nTx, nMode
        integer               :: ios
        character(len=20)     :: version
        !
        ! Verbose...
        write( *, * ) "     > Write ESolution to file: [", e_solution_file_name, "]"
        !
        version = "Modem-OO"
        !
        open( ioESolution, file = e_solution_file_name, action = "write", form = "unformatted", iostat = ios)
        !
        if( ios == 0 ) then
            !
            write( ioESolution ) version, nTx, nMode, &
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
            write( *, * ) "Error opening file in writeEsolutionHeader [", e_solution_file_name, "]!"
            stop
            !
        endif
        !
        !
    end subroutine writeEsolutionHeader
    !
    subroutine writeDataGroupArray( data_group_array )
        implicit none
        !
        type( DataGroup_t ), allocatable, dimension(:), intent( inout ) :: data_group_array
        !
        class( Transmitter_t ), pointer :: transmitter
        class( Receiver_t ), pointer :: receiver
        type( DataGroup_t ) :: data_group
        !
        integer :: receiver_type, i, j, array_size, ios
        !
        ! Verbose...
        write( *, * ) "     > Write Predicted Data to file: [", predicted_data_file_name, "]"
        !
        receiver_type = 0
        !
        open( unit = ioPredData, file = predicted_data_file_name, action = "write", form = "formatted", iostat = ios )
        !
        if( ios == 0 ) then
            !
            array_size = size( data_group_array )
            !
            do i = 1, array_size
                !
                data_group = getDataGroupByIndex( data_group_array, i )
                !
                receiver => getReceiver( data_group%id_rx )
                !
                call writePredictedDataHeader( receiver, receiver_type )
                !
                transmitter => getTransmitter( data_group%id_tx )
                !
                do j = 1, data_group%n_data
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
                            stop "Error: Unclassified data_group!"
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
            write( *, * ) "Error opening [", predicted_data_file_name, "] in writeDataGroupArray!"
            stop
        end if
        !
    end subroutine writeDataGroupArray
    !
    subroutine writePredictedDataHeader( receiver, receiver_type )
        implicit none
        !
        !
        class( Receiver_t ), intent( in ) :: receiver
        !
        integer, intent( inout ) :: receiver_type
        !
        if( receiver_type /= receiver%rx_type ) then
            !
            write( ioPredData, "(4A, 40A)" ) "#    ", DATA_FILE_TITLE
            !
            select case( receiver%rx_type )
                !
                case( 1, 11, 12 )
                    write( ioPredData, "(74A)" ) "#    Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error"
                case( 2 )
                    write( ioPredData, "(60A)" ) "#    Missing title for Full_Interstation_TF"
                case( 3 )
                    write( ioPredData, "(60A)" ) "#    Missing title for Off_Diagonal_Rho_Phase"
                case( 4 )
                    write( ioPredData, "(60A)" ) "#    Missing title for Phase_Tensor"
                case( 5 )
                    write( ioPredData, "(60A)" ) "#    Missing title for Off_Diagonal_Impedance"
                case( 6, 7, 8, 9, 10 )
                    write( ioPredData, "(125A)" ) "#    Tx_Dipole Tx_Period(s) Tx_Moment(Am) Tx_Azi Tx_Dip Tx_X(m) Tx_Y(m) Tx_Z(m) Code X(m) Y(m) Z(m) Component Real Imag Error"
                case default
                    write( *, * ) "Unknown receiver type :[", receiver%rx_type, "]"
                    stop "Error: test_FWD.f90: writePredictedDataHeader()"
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
    end subroutine writePredictedDataHeader
    !
    subroutine printHelp()
        implicit none
        !
        write( *, * ) "ModEM-OO Usage:"
        write( *, * ) ""
        write( *, * ) "    Flags to define a job:"
        write( *, * ) "        [-f], [--forward]    :  Forward modeling."
        write( *, * ) "        [-i], [--inversion]  :  Inversion."
        write( *, * )
        write( *, * ) "    Other arguments:"
        write( *, * ) "        [-d], [--data]       :  Flags for input data file path."
        write( *, * ) "        [-m], [--model]      :  Flags for input model file path."
        write( *, * ) "        [-c], [--control]    :  Flags for user control file path."
        write( *, * ) "        [-v], [--version]    :  Print version."
        write( *, * ) "        [-h], [--help]       :  Print usage information."
        write( *, * ) "        [-pd], [--predicted] :  Output data file path."
        write( *, * ) "        [-es], [--esolution] :  Output binary e-solution file path."
        write( *, * ) "        [--verbose]          :  Print runtime information."
        !
        write( *, * ) ""
        write( *, * ) "Version 1.0.0"
        !
        stop
        !
    end subroutine printHelp
    !
end program ModEM
!