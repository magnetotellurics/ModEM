program ModEM
    !
    use Constants
    use FileUnits
    !
    use DeclarationMPI
    !
    use ModEMControlFile
    !
    use Grid3D_SG
    !
    use ModelReader
    use ModelReader_Weerachai
    use ModelOperator_MF
    use ModelOperator_File
    use ModelParameterCell_SG
    !
    use DataFileStandard
    !
    use ForwardSolverFromFile
    use ForwardSolverIT_DC
    !
    use SourceMT_1D
    use SourceMT_2D
    !
    use ModelParameter2D
    !
    character(:), allocatable :: control_file_name, model_file_name, data_file_name, modem_job
    logical                   :: has_control_file, has_model_file, has_data_file, verbosis
    !
    call MPI_Init( ierr )
    !
    main_comm = MPI_COMM_WORLD
    !
    ! SET mpi_size WITH n FROM MPIRUN FOR MPI_COMM_WORLD
    call MPI_Comm_size( main_comm, mpi_size, ierr )
    !
    if( mpi_size < 2 ) then
        write( *, * ) "A minimum of two MPI processes are required!!!"
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
    write( *, * ) "MPI Rank ", mpi_rank," in COMM_WORLD (", mpi_size, ") is ", node_rank, &
                  " in SHARED_COMM (", node_size, ") on Node: ", node_name( 1 : nodestringlen )
    !
    ! MASTER
    !
    if ( mpi_rank == 0 ) then
        !
        modem_job = "unknow"
        !
        call setupDefaultParameters()
        !
        ! Validate arguments, set model_file_name, data_file_name, control_file_name, etc...
        call handleArguments()
        !
        write ( *, * )
        write ( *, * ) "Start ModEM-OO."
        write ( *, * )
        !
        ! Check parameters at the control file
        if( has_control_file ) call handleControlFile()
        !
        ! Execute the modem_job
        call handleJob()
        !
        call MPI_Finalize( ierr )
        !
        write ( *, * )
        write ( *, * ) "Finish ModEM-OO."
        write ( *, * )
    !
    ! WORKER
    !
    else
        !
        call MPI_Win_allocate_shared( shared_window_size, shared_disp_unit, MPI_INFO_NULL, child_comm, shared_c_ptr, shared_window, ierr )
        !
        if( ierr == MPI_SUCCESS ) then
            !
            do while ( job_master .ne. job_finish )
                !
                write ( *, * ) "WORKER: ", mpi_rank, " WAITING MASTER"
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
                case ( "JOB_FORWARD" )
                    !
                    call workerForwardModelling()
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
            call MPI_Finalize( ierr )
            !
        else
            !
            write( *, "(A50, i8)" ) "MPI Win_allocate_shared fails for process:", mpi_rank
            !
            stop
        endif
        !
     endif
    !
contains
    !
    subroutine masterForwardModelling()
        implicit none
        !
        ! Local variables
        integer :: i, worker_rank, tx_received, tx_index
        !
        type( DataHandle_t ), allocatable :: data_handles(:)
        type( DataHandle_t ), allocatable :: all_data_handles(:)
        !
        ! Verbosis
        write ( *, * ) "    > Start forward modelling."
        !
        worker_rank = 1
        tx_received = 0
        tx_index = 0
        !
        ! Reads Model File: instantiates Grid, ModelOperator and ModelParameter
        if( .NOT. has_model_file ) then 
            stop " - Missing Model file!"
        else
            call handleModelFile()
        endif
        !
        ! Reads Data File: instantiates and builds the Data relation between Txs and Txs
        if( .NOT. has_data_file ) then 
            stop " - Missing Data file!"
        else
            call handleDataFile()
        endif
        !
        ! SHARE MEM WITH ALL WORKERS
        call masterExposeSharedMemory()
        !
        do i = 1, ( mpi_size - 1 )
            !
            fwd_info%job_name = job_share_memory
            !
            call sendTo( i )
            !
        enddo
        !
        !
        call MPI_Win_fence( 0, shared_window, ierr )
        !
        deallocate( model_operator )
        deallocate( model_parameter )
        deallocate( main_grid )
        !
        ! SEND 1 TRANSMITTER TO FIRST np WORKERS
        do while ( worker_rank <= ( mpi_size - 1 ) )
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
        ! SEND 1 TRANSMITTER TO FIRST AVAILABLE WORKER
        do while( tx_index < size( transmitters ) )
            !
            write ( *, * ) "THERES", ( size( transmitters ) - tx_index ), " TX LEFT!"
            !
            call receiveFromAny()
            !
            data_handles = receiveData()
            !
            do i = 1, size( data_handles )
                call updateDataHandleArray( all_data_handles, data_handles( i ) )
            end do
            !
            deallocate( data_handles )
            !
            tx_received = tx_received + 1
            !
            tx_index = tx_index + 1
            !
            fwd_info%job_name    = job_forward
            fwd_info%tx_index    = tx_index
            !
            call sendTo( fwd_info%worker_rank )
            !
        end do
        !
        ! RECEIVES job_done FROM EACH FINISHED WORKER
        do while ( tx_received < size( transmitters ) )
            !
            write ( *, * ) "MASTER WAITING ANY WORKER TO FINISH"
            !
            call receiveFromAny()
            !
            data_handles = receiveData()
            !
            do i = 1, size( data_handles )
                call updateDataHandleArray( all_data_handles, data_handles( i ) )
            end do
            !
            deallocate( data_handles )
            !
            tx_received = tx_received + 1
            !
            fwd_info%job_name = job_finish
            !
            call sendTo( fwd_info%worker_rank )
             
        enddo
        !
        call deallocateTransmitterArray()
        !
        ! Write all_data_handles into predicted_data.dat
        call writeDataHandleArray( all_data_handles )
        !
        deallocate( all_data_handles )
        !
        call deallocateReceiverArray()
        !
        ! Verbosis
        write ( *, * ) "    > Finish forward modelling."
        !
    end subroutine masterForwardModelling
    !
    subroutine masterExposeSharedMemory()
        implicit none
        !
        select type( model_operator )
          !
          class is( ModelOperator_MF_t )
                !
                call allocateSharedBuffer()
                !
                shared_window_size = shared_buffer_size
                shared_disp_unit = 1
                !
                call MPI_Win_allocate_shared( shared_window_size, shared_disp_unit, MPI_INFO_NULL, child_comm, shared_c_ptr, shared_window, ierr )
                !
                if( ierr /= MPI_SUCCESS ) then
                     write( *, "(A50, i8)" ) "MPI Win_allocate_shared fails on master:", ierr
                     stop
                endif
                !
                call c_f_pointer( shared_c_ptr, shared_buffer, (/shared_window_size/) )
                !
                call packSharedBuffer()
                !
        end select
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
            write( *, "(A50, i8, i8, i8)" ) "MPI Allocated window size:", shared_window_size
            !
        else
            write( *, "(A50, i8, i8)" ) "MPI Win_shared_query fails on worker, ierr:", mpi_rank, ierr
            stop
        endif
        !
        call c_f_pointer( shared_c_ptr, shared_buffer, (/shared_window_size/) )
        !
        call unpackSharedBuffer( int( shared_window_size ) )
        !
    end subroutine workerQuerySharedMemory
    !
    subroutine workerForwardModelling()
        implicit none
        !
        class( ForwardSolver_t ), allocatable :: fwd_solver
        !
        class( Source_t ), allocatable        :: fwd_source 
        !
        ! Temporary alias pointers
        class( Receiver_t ), pointer    :: Rx
        class( Transmitter_t ), pointer :: Tx
        !
        type( DataHandle_t ), allocatable :: tx_data_handles(:)
        !
        type( TAirLayers ) :: air_layer
        !
        integer :: iRx, nRx, iDe
        !
        select type( main_grid )
            !
            class is( Grid3D_SG_t )
                !
                call main_grid%SetupAirLayers( air_layer, model_method, model_n_air_layer, model_max_height )
                call main_grid%UpdateAirLayers( air_layer%nz, air_layer%dz )
                !
            class default
                stop "Unclassified main_grid"
            !
        end select
        !
        call model_operator%metric%SetMetricElements()
        !
        call model_parameter%SetSigMap( model_parameter%paramType )
        !call model_parameter%SetType( LOGE )
        call model_parameter%setMetric( model_operator%metric )
        !
        call model_operator%SetEquations()
        call model_operator%SetCond( model_parameter )
        !
        ! ForwardSolver - Chosen from control file
        select case ( trim( forward_solver_type ) )
            !
            case( FWD_FILE )
                fwd_solver = ForwardSolverFromFile_t( model_operator )
            !
            case( FWD_IT_DC )
                fwd_solver = ForwardSolverIT_DC_t( model_operator, QMR )
            !
            case default
                fwd_solver = ForwardSolverIT_DC_t( model_operator, QMR )
            !
        end select
        !
        call fwd_solver%setCond( model_parameter )
        !
        ! Source - Chosen from control file
        select case ( trim( source_type ) )
            !
            case( SRC_MT_1D )
                fwd_source = SourceMT_1D_t( model_operator, model_parameter )
            !
            case( SRC_MT_2D )
                fwd_source = SourceMT_2D_t( model_operator, model_parameter )
            !
            case default
                fwd_source = SourceMT_1D_t( model_operator, model_parameter )
            !
        end select
        !
        Tx => getTransmitter( fwd_info%tx_index )
        !
        call writeEsolutionHeader( Tx%n_pol )
        !
        ! Verbosis...
        write( *, * ) "    Tx Id:", Tx%id, "Period:", int( Tx%period )
        !
        ! Tx points to its due Source
        call Tx%setSource( fwd_source )
        !
        ! Set ForwardSolver Period
        call fwd_solver%setPeriod( Tx%period )
        !
        ! Tx points to its due ForwardSolver
        call Tx%setForwardSolver( fwd_solver )
        !
        ! Solve Tx Forward Modelling
        call Tx%solveFWD()
        !
        deallocate( fwd_source )
        deallocate( fwd_solver )
        !
        ! Loop over Receivers of each Transmitter
        nRx = size( Tx%receiver_indexes )
        !
        do iRx = 1, nRx
            !
            ! Temporary Receiver alias
            Rx => getReceiver( Tx%receiver_indexes( iRx ) )
            !
            ! Verbosis...
            !write( *, * ) "                        Rx Id:", Rx%id, "XYZ:", Rx%location
            !
            ! Calculate Rx Predicted Data
            call Rx%predictedData( model_operator, Tx )
            !
            ! Store Rx predicted data into tx_data_handles
            do iDe = 1, size( Rx%predicted_data )
            !
                call updateDataHandleArray( tx_data_handles, Rx%predicted_data(iDe) )
            !
            end do
            !
            deallocate( Rx%predicted_data )
            !
        enddo
        !
        write ( *, * ) "WORKER ", mpi_rank, "FINISHES FWD FOR TX ", Tx%id!, size( tx_data_handles )
        !
        ! SEND JOB DONE TO MASTER
        fwd_info%job_name    = job_fwd_done
        fwd_info%tx_index    = Tx%id
        fwd_info%worker_rank = mpi_rank
        !
        deallocate( Tx )
        !
        call allocateDataBuffer( tx_data_handles )
        fwd_info%n_data      = size( tx_data_handles )
        fwd_info%data_size   = predicted_data_buffer_size
        !
        call sendTo( master_id )
        !
        call sendData( tx_data_handles )
        !
        deallocate( tx_data_handles )
        !
    end subroutine workerForwardModelling
    !
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
        case default
          !
          write( *, * ) "    - Unknow job: [", trim( modem_job ), "]"
          call printHelp()
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
        write( *, * ) "    -> Control File: [", control_file_name, "]"
        !
        ! Instantiate the ControlFile object
        ! Reads control file and sets the options in the Constants module
        control_file = ModEMControlFile_t( ioStartup, control_file_name )
        !
    end subroutine handleControlFile
    !
    !
    subroutine handleDataFile()
        implicit none
        !
        ! Local object to dealt data, self-destructs at the end of the subroutine
        type( DataFileStandard_t ) :: data_file_standard
        !
        write( *, * ) "    -> Data File: [", data_file_name, "]"
        !
        data_file_standard = DataFileStandard_t( ioStartup, data_file_name )
        !
        if( size( receivers ) == data_file_standard%nRx ) then
             write( *, * ) size( receivers ), " Receivers checked!"
        else
             !
             write(*,*) "Number of Rx mismatched from Header :[", size( receivers ), " and ", data_file_standard%nRx, "]"
             STOP "DataManager.f08: DataManager_ctor()"
             !
        endif
        !
        if( size( transmitters ) == data_file_standard%nTx ) then
             write( *, * ) size( transmitters ), " Transmitters checked!"
             !
             call printTransmitterArray()
        else
             !
             write(*,*) "Number of Tx mismatched from Header :[", size( transmitters ), " and ", data_file_standard%nTx, "]"
             STOP "DataManager.f08: DataManager_ctor()"
             !
        endif
        !
    end subroutine handleDataFile
    !
    subroutine handleModelFile()
        implicit none
        !
        ! It remains to standardize ????
        type( ModelReader_Weerachai_t ) :: model_reader
        !
        type( TAirLayers )                :: air_layer
        !
        write( *, * ) "    -> Model File: [", model_file_name, "]"
        !
        ! Read Grid and ModelParameter with ModelReader_Weerachai
        call model_reader%Read( model_file_name, main_grid, model_parameter ) 
        !
        ! Instantiate the ModelOperator object
        select type( main_grid )
          !
          class is( Grid3D_SG_t )
              !
              write( *, * ) "MASTER GRID: model_method, model_n_air_layer, model_max_height", model_method, model_n_air_layer, model_max_height
              !
              call main_grid%SetupAirLayers( air_layer, model_method, model_n_air_layer, model_max_height )
              !
              call main_grid%UpdateAirLayers( air_layer%nz, air_layer%dz )
              !
              model_operator = ModelOperator_MF_t( main_grid )
              !
              call model_parameter%setMetric( model_operator%metric )
              !
              call model_operator%SetEquations()
              !
              class default
                  stop "Unclassified main_grid"
              !
        end select
        !
    end subroutine handleModelFile
    !
    subroutine handleArguments()
        implicit none
        !
        character(200) :: argument
        integer         :: argument_index
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
                      if ( len( control_file_name ) > 0 ) has_control_file = .true.
                      !
                      argument_index = argument_index + 2
                      !
                    case ( "-d", "--data" )
                      !
                      call get_command_argument( argument_index + 1, argument )
                      data_file_name = trim( argument )
                      !
                      if ( len( data_file_name ) > 0 ) has_data_file = .true.
                      !
                      argument_index = argument_index + 2
                      !
                    case ( "-f", "--forward" )
                      !
                      modem_job = "forward"
                      !
                      argument_index = argument_index + 1
                      !
                    case ( "-m", "--model" )
                      !
                      call get_command_argument( argument_index + 1, argument )
                      model_file_name = trim(argument)
                      !
                      if ( len( model_file_name ) > 0 ) has_model_file = .true.
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
                    case default
                      !
                      write( *, * ) "    - Unknow Argument: [", trim( argument ), "]"
                      call printHelp()
                      !
                end select
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
        model_method      = MM_METHOD_FIXED_H
        model_n_air_layer = 10
        model_max_height  = 200.0
        !
        source_type = SRC_MT_1D
        !
        forward_solver_type = FWD_IT_DC
        !
        !
        has_control_file = .FALSE.
        has_model_file   = .FALSE.
        has_data_file    = .FALSE.
        verbosis         = .FALSE.
        !
    end subroutine setupDefaultParameters
    !
    subroutine writeEsolutionHeader( nMode )
        implicit none
        !
        integer, intent( in ) :: nMode
        !
        integer            :: ios
        character(:), allocatable :: version
        !
        version = ""
        !
        open( ioESolution, file = "e_solution", action="write", form ="unformatted", iostat=ios )
        !
        if( ios == 0 ) then
            !
            version = ""
            ! write the header (contains the basic information for the forward
            ! modeling). the header is 4 lines
            write( ioESolution ) version, size( transmitters ), nMode, &
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
            write( *, * ) "writeEsolutionHeader: e_solution"
            stop
        endif
        !
    end subroutine writeEsolutionHeader
    !
    subroutine writeDataHandleArray( data_handle_array )
        implicit none
        !
        type( DataHandle_t ), allocatable, intent( inout ) :: data_handle_array(:)
        !
        type( DataHandle_t ) :: aux_data_entry
        character(:), allocatable :: receiver_type
        integer :: i, j, ios
        !
        receiver_type = "Unknow"
        !
        ! Order by transmitter
        do i = 1, size( data_handle_array ) - 1
          !
          do j = i + 1, size( data_handle_array )
              !
              if( data_handle_array(i)%period > data_handle_array(j)%period ) then
                aux_data_entry  = data_handle_array(i)
                data_handle_array(i) = data_handle_array(j)
                data_handle_array(j) = aux_data_entry
              endif
              !
          enddo
        enddo
        !
        ! Order by receiver
        do i = 1, size( data_handle_array ) - 1
          !
          do j = i + 1, size( data_handle_array )
              !
              if( data_handle_array(i)%rx_id > data_handle_array(j)%rx_id ) then
                aux_data_entry  = data_handle_array(i)
                data_handle_array(i) = data_handle_array(j)
                data_handle_array(j) = aux_data_entry
              endif
              !
          enddo
        enddo
        !
        open( ioPredData, file = "predicted_data.dat", action = "write", form = "formatted", iostat = ios )
        !
        if( ios == 0 ) then
            !
            do i = 1, size( data_handle_array )
                !
                call writePredictedDataHeader( data_handle_array(i), receiver_type )
                !
                write( ioPredData, "(es12.6, A20, f15.3, f15.3, f15.3, f15.3, f15.3, A20, es16.6, es16.6, es16.6)" ) data_handle_array(i)%period, data_handle_array(i)%code, R_ZERO, R_ZERO, data_handle_array(i)%xyz(1), data_handle_array(i)%xyz(2), data_handle_array(i)%xyz(3), data_handle_array(i)%component, data_handle_array(i)%real, data_handle_array(i)%imaginary, 1.0
                !
            enddo
            !
            close( ioPredData )
            !
        else
          stop "Error opening predicted_data.dat in writeDataHandleArray"
        end if
        !
    end subroutine writeDataHandleArray
    !
    subroutine writePredictedDataHeader( data_handle, receiver_type )
        implicit none
        !
        type( DataHandle_t ), intent( in ) :: data_handle
        character(:), allocatable, intent( inout ) :: receiver_type
        !
        class( Receiver_t ), pointer :: receiver
        !
        receiver => getReceiver( data_handle%rx_id )
        !
        if( receiver_type /= trim( receiver%type_name ) ) then
            !
            write( ioPredData, "(4A, 100A)" ) "#    ", DATA_FILE_TITLE
            write( ioPredData, "(100A)" )     "#    Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error"
            write( ioPredData, "(4A, 100A)" ) ">    ", trim( receiver%type_name )
            write( ioPredData, "(4A, 100A)" ) ">    ", "exp(-i\omega t)"
            write( ioPredData, "(4A, 100A)" ) ">    ", "[V/m]/[T]"
            write( ioPredData, "(7A, 100A)" ) ">        ", "0.00"
            write( ioPredData, "(7A, 100A)" ) ">        ", "0.000    0.000"
            write( ioPredData, "(A3, i8, i8)" ) ">        ", size( transmitters ), size( receivers )
            !
            receiver_type = trim( receiver%type_name )
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
        write( *, * ) "        [-f], [--forward] : Forward modelling."
        write( *, * ) "        [-i], [--inverse] : Inversion modelling."
        write( *, * )
        write( *, * ) "    Other arguments:"
        write( *, * ) "        [-d], [--data]     : Flags for data file path."
        write( *, * ) "        [-m], [--model]    : Flags for model file path."
        write( *, * ) "        [-c], [--control] : Flags for user control file path."
        write( *, * ) "        [-v], [--version] : Print version."
        write( *, * ) "        [-h], [--help]     : Print usage information."
        !
        write( *, * ) ""
        write( *, * ) "Version 1.0.0"
        !
        stop
        !
    end subroutine printHelp

end program ModEM
