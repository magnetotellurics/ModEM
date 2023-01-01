!
!> Prototype MPI version for testing the ModEM-OO program
!> To merge with the serial version in the future
!
program TestMPI
    !
    use DeclarationMPI
    !
    real( kind=prec ) :: t_start, t_finish
    !
    call cpu_time( t_start )
    !
    call MPI_Init( ierr )
    !
    !> Set mpi_size with np from mpirun for mpi_comm_world
    call MPI_Comm_size( MPI_COMM_WORLD, mpi_size, ierr )
    !
    if( mpi_size < 2 ) then
        write( *, * ) "Error: Minimum of two processes required!"
        call MPI_Finalize( ierr )
        stop
    end if 
    !
    !> Set mpi_rank with process id for mpi_comm_world
    call MPI_Comm_rank( MPI_COMM_WORLD, mpi_rank, ierr )
    !
    call MPI_Info_create( mpi_win_info, ierr )
    call MPI_Info_set( mpi_win_info, "alloc_shared_noncontig", "true", ierr )
    !
    !> Split mpi_comm_world into shared sub communicator: node_comm
    call MPI_Comm_split_type( MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, mpi_win_info, node_comm, ierr )
    !
    call MPI_Get_processor_name( node_name, nodestringlen, ierr )
    !
    call MPI_Comm_size( node_comm, node_size, ierr )
    !
    call MPI_Comm_rank( node_comm, node_rank, ierr )
    !
    !write( *, * ) "MPI Rank ", mpi_rank," in COMM_WORLD [", mpi_size, "] is ", node_rank, &
    !>              " in SHARED_COMM [", node_size, "] on Node: ", node_name( 1 : nodestringlen )
    !
    !> MPI MASTER PROCESS
    !
    if( node_rank == 0 ) then
        !
        modem_job = "unknown"
        !
        call setupDefaultParameters()
        !
        !> Validate arguments and control variables
        call handleArguments()
        !
        write( *, * )
        write( *, * ) "Start ModEM-OO."
        write( *, * )
        !
        !> Check parameters at the control file
        if( has_control_file ) call handleControlFile()
        !
        !> Execute the modem_job
        call handleJob()
        !
        !> Deallocate remaining main program memory
        call garbageCollector()
        !
        call MPI_Finalize( ierr )
        !
        call cpu_time( t_finish )
        !
        write( *, * )
        write( *, "(A18, F10.2, A3)" ) "Finish ModEM-OO (", t_finish - t_start, "s)"
        write( *, * )
    !
    !> MPI WORKER PROCESS
    !
    else
        !
        call MPI_Win_allocate_shared( shared_window_size, shared_disp_unit, MPI_INFO_NULL, node_comm, shared_c_ptr, shared_window, ierr )
        !
        if( ierr == MPI_SUCCESS ) then
            !
            !> Receive while master does not send finish message
            do while ( job_master .NE. job_finish )
                !
                call receiveFrom( master_id )
                !
                job_master = trim( job_info%job_name )
                !
                select case ( job_master )
                    !
                    case ( "SHARE_MEMORY" )
                        !
                        call workerQuerySharedMemory()
                        !
                    case ( "JOB_FORWARD" )
                        !
                        call workerJobForwardModelling()
                        !
                    case ( "JOB_ADJOINT" )
                        !
                        call workerJobAdjoint()
                        !
                    case ( "JOB_ADJOINT_T" )
                        !
                        call workerJobAdjoint_T()
                        !
                    case ( "JOB_INVERSION" )
                        !
                        !call workerJobInversion()
                        !
                end select
                !
            enddo
            !
            !> Deallocate global components
            deallocate( forward_solver, model_operator, sigma0, main_grid )
            !
            !> Deallocate prior model if is the case
            if( allocated( pmodel ) ) deallocate( pmodel )
            !
            !> Deallocate global receiver array
            call deallocateReceiverArray()
            !
            !> Deallocate global transmitter array
            call deallocateTransmitterArray()
            !
            !> Deallocate remaining worker memory
            call garbageCollector()
            !
            call cpu_time( t_finish )
            !
            write( *, * ) "Worker", node_rank, " finished ( ", t_finish - t_start, "s )"
            !
            call MPI_Finalize( ierr )
            !
        else
            !
            write( *, * ) "MPI Win_allocate_shared fails for process:", node_rank
            !
            stop
            !
        endif
        !
    endif
    !
contains
    !
    !> Share memory with all workers
    subroutine shareMemory()
        implicit none
        !
        integer :: worker_id
        !
        !> Share memory with all workers
        call masterExposeSharedMemory()
        !
        do worker_id = 1, ( mpi_size - 1 )
            !
            job_info%job_name = job_share_memory
            !
            call sendTo( worker_id )
            !
        enddo
        !
        call MPI_Win_fence( 0, shared_window, ierr )
        !
        call MPI_BARRIER( MPI_COMM_WORLD, ierr )
        !
        nullify( shared_buffer )
        !
        call MPI_Win_free( shared_window, ierr )
        !
    end subroutine shareMemory
    !
    !> No procedure briefing
    subroutine masterExposeSharedMemory()
        implicit none
        !
        call allocateSharedBuffer()
        !
        shared_window_size = shared_buffer_size
        shared_disp_unit = 1
        !
        call MPI_Win_allocate_shared( shared_window_size, shared_disp_unit, MPI_INFO_NULL, node_comm, shared_c_ptr, shared_window, ierr )
        !
        if( ierr /= MPI_SUCCESS ) then
             write( *, * ) "Error: MPI Win_allocate_shared fails on master:", ierr
             stop
        endif
        !
        call c_f_pointer( shared_c_ptr, shared_buffer, (/shared_window_size/) )
        !
        !call MPI_Alloc_mem( shared_window_size, MPI_INFO_NULL, shared_buffer, ierr ) 
        !
        call packSharedBuffer()
        !
    end subroutine masterExposeSharedMemory
    !
    !> No procedure briefing
    subroutine workerQuerySharedMemory()
        implicit none
        !
        call MPI_Win_shared_query( shared_window, master_id, shared_window_size, shared_disp_unit, shared_c_ptr, ierr )
        !
        if( ierr == MPI_SUCCESS ) then
            !
            write( *, * ) "          Worker", mpi_rank, " query ", shared_window_size, " bytes."
            !
            call c_f_pointer( shared_c_ptr, shared_buffer, (/shared_window_size/) )
            !
            call unpackSharedBuffer( int( shared_window_size ) )
            !
            call MPI_Win_fence( 0, shared_window, ierr )
            !
            call MPI_BARRIER( MPI_COMM_WORLD, ierr )
            !
            nullify( shared_buffer )
            !
            call MPI_Win_free( shared_window, ierr )
            !
            !> Instantiate the ModelOperator object
            select type( main_grid )
                !
                class is( Grid3D_SG_t )
                    !
                    allocate( model_operator, source = ModelOperator_MF_t( main_grid ) )
                    !
                    call model_operator%SetEquations()
                    !
                    call sigma0%setMetric( model_operator%metric )
                    !
                    call model_operator%SetCond( sigma0 )
                    !
                class default
                    stop "Error: workerQuerySharedMemory > Unclassified main_grid"
                !
            end select
            !
        else
            write( *, * ) "Error: MPI Win_shared_query fails on worker:", mpi_rank, ierr
            stop
        endif
        !
    end subroutine workerQuerySharedMemory
    !
    !> No procedure briefing
    subroutine masterForwardModelling( adjoint )
        implicit none
        !
        logical, intent( in ), optional :: adjoint
        !
        logical :: is_adjoint
        !
        integer :: worker_rank, tx_received, i_tx, i_data
        !
        !> Verbose
        write( *, * ) "     - Start Forward Modeling"
        !
        !>
        if( present( adjoint ) ) then
            is_adjoint = adjoint
        else
            is_adjoint = .FALSE.
        endif
        !
        !> Create an array of DataGroupTx to store the predicted data in the same format as the measured data
        !> Enabling the use of the tx grouped predicted data in future jobs (all_predicted_data)
        do i_tx = 1, size( transmitters )
            !
            allocate( tx_data, source = DataGroupTx_t( i_tx ) )
            !
            do i_data = 1, size( measured_data )
                !
                if( measured_data( i_data )%i_tx == i_tx ) then
                    !
                    call tx_data%addData( measured_data( i_data ) )
                    !
                endif
            enddo
            !
            call updateDataGroupTxArray( all_predicted_data, tx_data )
            !
            deallocate( tx_data )
            !
        enddo
        !
        !> Initialize MPI control variables
        worker_rank = 1
        tx_received = 0
        i_tx       = 0
        !
        !> Send 1 transmitter to first np workers
        do while( worker_rank <= ( mpi_size - 1 ) )
            !
            i_tx = i_tx + 1
            !
            job_info%job_name    = job_forward
            job_info%adjoint     = is_adjoint
            job_info%i_tx       = i_tx
            job_info%worker_rank = worker_rank
            !
            call sendTo( worker_rank )
            !
            worker_rank = worker_rank + 1
            !
        enddo
        !
        !> Send 1 transmitter to first available worker
        do while( i_tx < size( transmitters ) )
            !
            call receiveFromAny()
            !
            call receiveData( all_predicted_data( job_info%i_tx ), job_info%worker_rank )
            !
            tx_received = tx_received + 1
            i_tx       = i_tx + 1
            !
            job_info%job_name = job_forward
            job_info%adjoint  = is_adjoint
            job_info%i_tx    = i_tx
            !
            call sendTo( job_info%worker_rank )
            !
        enddo
        !
        !> Receive job_done from each worker
        do while( tx_received < size( transmitters ) )
            !
            call receiveFromAny()
            !
            call receiveData( all_predicted_data( job_info%i_tx ), job_info%worker_rank )
            !
            tx_received = tx_received + 1
            !
            if( .NOT. is_adjoint ) then
                !
                job_info%job_name = job_finish
                !
                call sendTo( job_info%worker_rank )
                !
            endif
            !
        enddo
        !
        !> Write all_predicted_data, with its proper Rx headers, to the file <predicted_data_file_name>
        call writeDataGroupArray( all_predicted_data, predicted_data_file_name )
        !
        !> Verbose
        write( *, * ) "     - Finish Forward Modeling"
        !
    end subroutine masterForwardModelling
    !
    !> Routine to run a full ForwardModeling job and deliver the result (PredictedData) in a text file
    subroutine masterJobForwardModelling()
        implicit none
        !
        ! Verbose
        write( *, * ) "     - Start jobForwardModeling"
        !
        !> Reads Model File: instantiates Grid, ModelOperator and ModelParameter
        if( .NOT. has_model_file ) then 
            stop "Error: jobForwardModeling > Missing Model file!"
        else
            call handleModelFile()
        endif
        !
        !> Reads Data File: instantiates and builds the Data relation between Txs and Txs
        if( .NOT. has_data_file ) then 
            stop "Error: jobForwardModeling > Missing Data file!"
        else
            call handleDataFile()
        endif
        !
        !> Share memory with all workers
        call shareMemory()
        !
        !> Deallocate Forward Modeling global components (Not used by the Master process)
        deallocate( model_operator, sigma0, main_grid )
        !
        !> Run masterForwardModelling to calculate predicted data
        call masterForwardModelling()
        !
        !>
        call deallocateGlobalArrays()
        !
        ! Verbose
        write( *, * ) "     - Finish jobForwardModeling"
        !
    end subroutine masterJobForwardModelling
    !
    !> Routine to run a full ForwardModeling job and deliver the result (PredictedData) in a text file
    subroutine masterJobAdjoint()
        implicit none
        !
        !> Data gradient for all transmitters, grouped into an array of DataGroupTx
        type( DataGroupTx_t ), allocatable, dimension(:) :: JmHat
        !
        integer :: worker_rank, i_tx, tx_received
        !
        ! Verbose
        write( *, * ) "     - Start masterJobAdjoint"
        !
        !> Read Model File: instantiates Grid, ModelOperator and ModelParameter
        if( .NOT. has_model_file ) then 
            stop "Error: jobAdjoint > Missing Model file!"
        else
            !
            call handleModelFile()
            !
        endif
        !
        !> Read Prior Model File: instantiates PModelParameter
        if( .NOT. has_pmodel_file ) then 
            stop "Error: jobAdjoint > Missing Prior Model file!"
        else
            !
            call handlePModelFile()
            !
        endif
        !
        !> Reads Data File: instantiates and builds the Data relation between Txs and Rxs
        if( .NOT. has_data_file ) then 
            stop "Error: masterJobAdjoint > Missing Data file!"
        else
            !
            call handleDataFile()
            !
        endif
        !
        !> Share memory with all workers
        call shareMemory()
        !
        !> Deallocate Forward Modeling global components
        deallocate( sigma0, pmodel, model_operator, main_grid )
        !
        !> Run masterForwardModelling only to calculate predicted data
        call masterForwardModelling( .TRUE. )
        !
        !>
        allocate( JmHat, source = all_predicted_data )
        !
        !> Initialize MPI control variables
        worker_rank = 1
        tx_received = 0
        i_tx       = 0
        !
        !> Send 1 transmitter to first np workers
        do while( worker_rank <= ( mpi_size - 1 ) )
            !
            i_tx = i_tx + 1
            !
            job_info%job_name    = job_adjoint
            job_info%worker_rank = worker_rank
            job_info%i_tx       = i_tx
            !
            call sendTo( worker_rank )
            !
            worker_rank = worker_rank + 1
            !
        enddo
        !
        !> Send 1 transmitter to first available worker
        do while( i_tx < size( transmitters ) )
            !
            call receiveFromAny()
            !
            call receiveData( JmHat( job_info%i_tx ), job_info%worker_rank )
            !
            tx_received = tx_received + 1
            i_tx       = i_tx + 1
            !
            job_info%job_name = job_adjoint
            !
            call sendTo( job_info%worker_rank )
            !
        enddo
        !
        !> Receives job_done from each finished worker
        do while( tx_received < size( transmitters ) )
            !
            call receiveFromAny()
            !
            call receiveData( JmHat( job_info%i_tx ), job_info%worker_rank )
            !
            tx_received = tx_received + 1
            !
            job_info%job_name = job_finish
            !
            call sendTo( job_info%worker_rank )
            !
        enddo
        !
        !> Write all_predicted_data, with its proper Rx headers, to the file <JmHat_data_file_name>
        call writeDataGroupArray( JmHat, JmHat_data_file_name )
        !
        !>
        call deallocateGlobalArrays()
        !
        !> Verbose
        write( *, * ) "     - Finish masterJobAdjoint"
        !
    end subroutine masterJobAdjoint
    !
    !
    subroutine sendTxMeasureData( i_tx )
        implicit none
        !
        integer, intent( in ) :: i_tx
        !
        integer :: i_data
        !
        !> Build the DataGroupTx to store the measured data from a single transmitter.
        allocate( tx_data, source = DataGroupTx_t( i_tx ) )
        !
        !> Fill res_data with the measured data
        do i_data = 1, size( measured_data )
            !
            if( measured_data( i_data )%i_tx == i_tx ) then
                !
                call tx_data%addData( measured_data( i_data ) )
                !
            endif
        enddo
        !
        call sendData( tx_data, job_info%worker_rank )
        !
        deallocate( tx_data )
        !
    end subroutine sendTxMeasureData
    !
    !> Routine to run a full ForwardModeling job and deliver the result (PredictedData) in a text file
    subroutine masterJobAdjoint_T()
        implicit none
        !
        type( ModelParameterCell_SG_t ) :: dsigma
        !
        class( Scalar_t ), allocatable :: tx_model_cond
        !
        integer :: worker_rank, i_tx, tx_received, i_data
        !
        ! Verbose
        write( *, * ) "     - Start masterJobAdjoint_T"
        !
        !> Read Model File: instantiates Grid, ModelOperator and ModelParameter
        if( .NOT. has_model_file ) then 
            stop "Error: jobAdjoint > Missing Model file!"
        else
            call handleModelFile()
        endif
        !
        !> Reads Data File: instantiates and builds the Data relation between Txs and Rxs
        if( .NOT. has_data_file ) then 
            stop "Error: masterJobAdjoint > Missing Data file!"
        else
            call handleDataFile()
        endif
        !
        !> Share memory with all workers
        call shareMemory()
        !
        !> Run masterForwardModelling to calculate predicted data and LRows
        call masterForwardModelling( .TRUE. )
        !
        !> Initialize dsigma as sigma0 and reset it for the sum of the conductivities.
        select type( sigma0 )
            !
            class is( ModelParameterCell_SG_t )
                !
                dsigma = sigma0
                !
                call dsigma%zeros()
                !
            class default
                stop "Error: masterJobAdjoint_T > Unclassified sigma0"
            !
        end select
        !
        !> Initialize MPI control variables
        worker_rank = 1
        tx_received = 0
        i_tx       = 0
        !
        !> Send 1 transmitter to first np workers
        do while( worker_rank <= ( mpi_size - 1 ) )
            !
            i_tx = i_tx + 1
            !
            job_info%job_name    = job_adjoint_t
            job_info%worker_rank = worker_rank
            job_info%i_tx       = i_tx
            !
            call sendTo( worker_rank )
            !
            !> Build a DataGroupTx to store the measured data from a single transmitter.
            !> and send it to worker
            call sendTxMeasureData( i_tx )
            !
            call sendData( all_predicted_data( i_tx ), job_info%worker_rank )
            !
            worker_rank = worker_rank + 1
            !
        enddo
        !
        !> Send 1 transmitter to first available worker
        do while( i_tx < size( transmitters ) )
            !
            call receiveFromAny()
            !
            allocate( tx_model_cond, source = dsigma%cell_cond )
            !
            call receiveModel( tx_model_cond, job_info%worker_rank )
            !
            call dsigma%cell_cond%add( tx_model_cond )
            !
            deallocate( tx_model_cond )
            !
            tx_received = tx_received + 1
            i_tx       = i_tx + 1
            !
            job_info%job_name = job_adjoint_t
            job_info%i_tx    = i_tx
            !
            call sendTo( job_info%worker_rank )
            !
            !> Build a DataGroupTx to store the measured data from a single transmitter.
            !> and send it to worker
            call sendTxMeasureData( i_tx )
            !
            call sendData( all_predicted_data( i_tx ), job_info%worker_rank )
            !
        enddo
        !
        !> Receives job_done from each finished worker
        do while( tx_received < size( transmitters ) )
            !
            call receiveFromAny()
            !
            allocate( tx_model_cond, source = dsigma%cell_cond )
            !
            call receiveModel( tx_model_cond, job_info%worker_rank )
            !
            call dsigma%cell_cond%add( tx_model_cond )
            !
            deallocate( tx_model_cond )
            !
            tx_received = tx_received + 1
            !
            job_info%job_name = job_finish
            !
            call sendTo( job_info%worker_rank )
            !
        enddo
        !
        !> Write DSigma to <dsigma_file_name> file path
        write( *, * ) "     > Write Dsigma to file: [", dsigma_file_name, "]"
        !
        call dsigma%write()
        !
        !> Deallocate Forward Modeling global components
        deallocate( sigma0, model_operator, main_grid )
        !
        call deallocateGlobalArrays()
        !
        !> Verbose
        write( *, * ) "     - Finish masterJobAdjoint_T"
        !
    end subroutine masterJobAdjoint_T
    !
    !>
    subroutine deallocateGlobalArrays()
        implicit none
        !
        !> Flush memory used by job data arrays
        deallocate( measured_data )
        !
        !> Flush memory used by job data arrays
        deallocate( all_predicted_data )
        !
        !> Deallocate global array of Receivers
        call deallocateReceiverArray()
        !
        !> Deallocate global array of Transmitters
        call deallocateTransmitterArray()
        !
    end subroutine deallocateGlobalArrays
    !
    !> No procedure briefing
    subroutine workerJobForwardModelling()
        implicit none
        !
        !> Temporary alias pointers
        class( Transmitter_t ), pointer :: Tx
        class( Receiver_t ), pointer :: Rx
        !
        integer :: iRx
        !
        if( .NOT. allocated( forward_solver ) ) then
            !
            !> Instantiate the ForwardSolver - Specific type can be chosen via control file
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
        endif
        !
        !> Point to the current Transmitter
        Tx => getTransmitter( job_info%i_tx )
        !
        !> Set Transmitter's ForwardSolver Omega(Period) and Conductivity
        call forward_solver%setFrequency( sigma0, Tx%period )
        !
        !> Point Transmitter's ForwardSolver
        Tx%forward_solver => forward_solver
        !
        !> Instantiate Transmitter's Source - According to tx type or via control file
        select type( Tx )
            !
            class is( TransmitterMT_t )
                !
                allocate( Tx%source, source = SourceMT_1D_t( model_operator, sigma0, Tx%period ) )
                !
            class is( TransmitterCSEM_t )
                !
                allocate( Tx%source, source = SourceCSEM_Dipole1D_t( model_operator, sigma0, Tx%period, Tx%location, Tx%dip, Tx%azimuth, Tx%moment ) )
                !
            class default
                stop "Error: workerJobForwardModelling: Unclassified Transmitter"
                !
        end select
        !
        !> Build Source E according to source type
        call Tx%source%createE()
        !
        !> Solve Forward Modeling for the Transmitter
        call Tx%solve()
        !
        !> Build the DataGroupTx to store residual data from a single transmitter.
        allocate( tx_data, source = DataGroupTx_t( Tx%id ) )
        !
        !> Loop for each Receiver related to the Transmitter
        do iRx = 1, size( Tx%receiver_indexes )
            !
            !> Point to the current Receiver
            Rx => getReceiver( Tx%receiver_indexes( iRx ) )
            !
            !> Calculate and store Predicted Data and/or LRows in the Receiver
            !> Depending on the type of work
            call Rx%predictedData( Tx )
            !
            call tx_data%addData( Rx%data_group )
            !
        enddo
        !
        !> SEND JOB DONE TO MASTER
        job_info%job_name    = job_done
        job_info%worker_rank = node_rank
        job_info%i_tx       = tx_data%i_tx
        !
        call sendTo( master_id )
        !
        call sendData( tx_data, master_id )
        !
        !> Clear the memory used by the current tx_data
        if( .NOT. job_info%adjoint ) deallocate( tx_data )
        !
    end subroutine workerJobForwardModelling
    !
    !> No procedure briefing
    subroutine workerJobAdjoint()
        implicit none
        !
        !> Temporary alias pointers
        class( Transmitter_t ), pointer :: Tx
        !
        !> Get the same transmitter previously used at forward modeling (tx_pred_data)
        Tx => getTransmitter( tx_data%i_tx )
        !
        !> Switch Transmitter's source to SourceInteriorForce from pMult
        call Tx%setSource( Tx%pMult( sigma0, pmodel, model_operator ) )
        !
        !> Solve e_sens with the new Source
        call Tx%solve()
        !
        !>
        call tx_data%reset()
        !
        !> JMult for the same Tx
        call JMult_Tx( tx_data )
        !
        !> MPI: SEND JOB DONE TO MASTER
        job_info%job_name    = job_done
        job_info%worker_rank = node_rank
        job_info%i_tx       = tx_data%i_tx
        !
        call sendTo( master_id )
        !
        call sendData( tx_data, master_id )
        !
    end subroutine workerJobAdjoint
    !
    !> No procedure briefing
    subroutine workerJobAdjoint_T()
        implicit none
        !
        class( ModelParameter_t ), allocatable :: tx_dsigma
        !
        type( DataGroupTx_t ) :: pred_data
        !
        !> Receive measure data for a single Tx
        call receiveData( tx_data, master_id )
        !
        !> Receive predicted data for a single Tx
        pred_data = tx_data
        !
        call receiveData( pred_data, master_id )
        !
        !> Subtracts and normalize from the respective predicted data
        call tx_data%getResidual( pred_data )
        !
        !> Set current tx_dsigma from JMult_T_Tx
        call JMult_T_Tx( sigma0, tx_data, tx_dsigma )
        !
        !> MPI: SEND JOB DONE TO MASTER
        job_info%job_name    = job_done
        job_info%worker_rank = node_rank
        !
        call sendTo( master_id )
        !
        !> Send result model to the Master
        select type( tx_dsigma )
            !
            class is( ModelParameterCell_SG_t )
                !
                call sendModel( tx_dsigma%cell_cond, master_id )
                !
            class default
                stop "Error: masterJobAdjoint_T > Unclassified tx_dsigma"
            !
        end select
        !
    end subroutine workerJobAdjoint_T
    !
    !> No subroutine briefing
    subroutine handleJob()
        implicit none
        !
        select case ( modem_job )
            !
            case ( "forward" )
                !
                call masterJobForwardModelling()
                !
            case ( "adjoint" )
                !
                call masterJobAdjoint()
                !
            case ( "JMult_t" )
                !
                call masterJobAdjoint_T()
                !
            case ( "inversion" )
                !
                !call masterJobInversion()
                !
            case default
                !
                write( *, * ) "Error: unknown job: [", modem_job, "]"
                call printHelp()
                stop
            !
        end select
        !
    end subroutine handleJob
    !
    !> No subroutine briefing
    subroutine handleControlFile()
        implicit none
        !
        type( ModEMControlFile_t ) :: control_file
        !
        write( *, * ) "     < Control File: [", control_file_name, "]"
        !
        !> Instantiate the ControlFile object
        !> Reads control file and sets the options in the Constants module
        control_file = ModEMControlFile_t( ioStartup, control_file_name )
        !
    end subroutine handleControlFile
    !
    !> No subroutine briefing
    subroutine handleModelFile()
        implicit none
        !
        type( ModelReader_Weerachai_t ) :: model_reader
        type( TAirLayers ) :: air_layer
        !
        ! Verbose
        write( *, * ) "     < Model File: [", model_file_name, "]"
        !
        !> Read Grid and ModelParameter with ModelReader_Weerachai
        call model_reader%Read( model_file_name, main_grid, sigma0 ) 
        !
        !> Instantiate the ModelOperator object
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
                write( *, "(A, F10.2, A2, F10.2, A2, F10.2, A4, F10.2)" ) "           o(x,y,z) * rotDeg: (", main_grid%ox, ", ", main_grid%oy, ", ", main_grid%oz, ") * ", main_grid%rotDeg
                !
                allocate( model_operator, source = ModelOperator_MF_t( main_grid ) )
                !
                call model_operator%SetEquations()
                !
                call sigma0%setMetric( model_operator%metric )
                !
                call model_operator%SetCond( sigma0 )
                !
            class default
                stop "Error: handleModelFile > Unclassified main_grid"
            !
        end select
        !
    end subroutine handleModelFile
    !
    !> No subroutine briefing
    subroutine handlePModelFile()
        implicit none
        !
        type( ModelReader_Weerachai_t ) :: model_reader
        !
        class( Grid_t ), allocatable :: prior_grid
        !
        ! Verbose
        write( *, * ) "     < Prior Model File: [", pmodel_file_name, "]"
        !
        !> Read prior_grid and pmodel with ModelReader_Weerachai
        call model_reader%Read( pmodel_file_name, prior_grid, pmodel ) 
        !
        call pmodel%setMetric( model_operator%metric )
        !
        deallocate( prior_grid )
        !
    end subroutine handlePModelFile
    !
    !> No subroutine briefing
    subroutine handleDataFile()
        implicit none
        !
        integer :: irx, nrx
        class( Receiver_t ), pointer :: Rx
        !
        !> Local object to dealt data, self-destructs at the end of the subroutine
        type( DataFileStandard_t ) :: data_file_standard
        !
        write( *, * ) "     < Data File: [", data_file_name, "]"
        !
        data_file_standard = DataFileStandard_t( ioStartup, data_file_name )
        !
        nrx = size( receivers )
        !
        ! Verbose
        if( nrx == data_file_standard%nRx ) then
            !
            write( *, * ) "          Checked ", nrx, " Receivers."
            !
        else
            !
            write( *, * ) "Error: Number of Rx mismatched from Header :[", nrx, " and ", data_file_standard%nRx, "]"
            stop
            !
        endif
        !
        if( size( transmitters ) == data_file_standard%nTx ) then
             !
             call printTransmitterArray()
             !
        else
             !
             write( *, * ) "Number of Tx mismatched from Header :[", size( transmitters ), " and ", data_file_standard%nTx, "]"
             stop
             !
        endif
        !
        write( *, * ) "          Checked ", size( measured_data ), " DataGroups."
        !
        write( *, * ) "     - Create Rx evaluation vectors"
        !
        !> Calculate evaluation vectors for the Receivers
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
    !> No subroutine briefing
    subroutine handleArguments()
        implicit none
        !
        character( len=200 ) :: argument
        integer :: argument_index
        !
        if( command_argument_count() == 0 ) then
            !
            call printHelp()
            call printUsage()
            stop
            !
        else
            !
            argument_index = 1
            !
            do while( argument_index <= command_argument_count() ) 
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
                         if( len( control_file_name ) > 0 ) has_control_file = .TRUE.
                         !
                         argument_index = argument_index + 2
                         !
                      case ( "-d", "--data" )
                         !
                         call get_command_argument( argument_index + 1, argument )
                         data_file_name = trim( argument )
                         !
                         if( len( data_file_name ) > 0 ) has_data_file = .TRUE.
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
                      case ( "-j", "--adjoint" )
                         !
                         modem_job = "adjoint"
                         !
                         argument_index = argument_index + 1
                         !
                      case ( "-jt", "--JMult_t" )
                         !
                         modem_job = "JMult_t"
                         !
                         argument_index = argument_index + 1
                         !
                      case ( "-m", "--model" )
                         !
                         call get_command_argument( argument_index + 1, argument )
                         model_file_name = trim( argument )
                         !
                         if( len( model_file_name ) > 0 ) has_model_file = .TRUE.
                         !
                         argument_index = argument_index + 2
                         !
                      case ( "-pm", "--pmodel" )
                         !
                         call get_command_argument( argument_index + 1, argument )
                         pmodel_file_name = trim( argument )
                         !
                         if( len( pmodel_file_name ) > 0 ) has_pmodel_file = .TRUE.
                         !
                         argument_index = argument_index + 2
                         !
                      case ( "-dm", "--dsigma" )
                         !
                         call get_command_argument( argument_index + 1, argument )
                         dsigma_file_name = trim( argument )
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
                      case ( "-gd", "--gradient" )
                         !
                         call get_command_argument( argument_index + 1, argument )
                         JmHat_data_file_name = trim( argument )
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
                         call printUsage()
                         stop
                         !
                      case ( "-vb", "--verbose" )
                         !
                         stop "Error: handleArguments > Verbose level not implemented yet!"
                         !
                      case default
                         !
                         write( *, * ) "Error: Unknown Argument: [", trim( argument ), "]"
                         call printHelp()
                         stop
                         !
                 end select
                 !
                 argument = ""
                 !
            enddo
            !
        endif
        !
    end subroutine handleArguments
    !
    !> No subroutine briefing
    subroutine setupDefaultParameters()
        implicit none
        !
        ! I|O
        predicted_data_file_name = "predicted_data.dat"
        JmHat_data_file_name  = "JmHat.dat"
        e_solution_file_name     = "esolution.bin"
        dsigma_file_name         = "dsigma.mod"
        has_control_file         = .FALSE.
        has_model_file           = .FALSE.
        has_pmodel_file          = .FALSE.
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
    !> No subroutine briefing
    subroutine garbageCollector()
        implicit none
        !
        !> Flush memory used by main program control variables and flags
        if( allocated( forward_solver_type ) ) deallocate( forward_solver_type )
        if( allocated( source_type ) ) deallocate( source_type )
        if( allocated( model_method ) ) deallocate( model_method )
        if( allocated( get_1D_from ) ) deallocate( get_1D_from )
        if( allocated( predicted_data_file_name ) ) deallocate( predicted_data_file_name )
        if( allocated( JmHat_data_file_name ) ) deallocate( JmHat_data_file_name )
        if( allocated( e_solution_file_name ) ) deallocate( e_solution_file_name )
        !
        if( allocated( control_file_name ) ) deallocate( control_file_name )
        if( allocated( model_file_name ) ) deallocate( model_file_name )
        if( allocated( pmodel_file_name ) ) deallocate( pmodel_file_name )
        if( allocated( data_file_name ) ) deallocate( data_file_name )
        if( allocated( dsigma_file_name ) ) deallocate( dsigma_file_name )
        if( allocated( modem_job ) ) deallocate( modem_job )
        !
    end subroutine garbageCollector
    !
    !> No subroutine briefing
    subroutine writeDataGroupArray( target_tx_data_array, file_name )
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: target_tx_data_array
        character(:), allocatable :: file_name
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
        write( *, * ) "     > Write Data to file: [", file_name, "]"
        !
        !> Put the data in the same input format
        allocate( to_write_data, source = measured_data )
        !
        do i = 1, size( target_tx_data_array )
            do j = 1, size( target_tx_data_array(i)%data )
                call setDataGroup( to_write_data, target_tx_data_array(i)%data(j) )
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
                call writePredictedDataHeader( receiver, receiver_type )
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
                            stop "Error: Unclassified data_group!"
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
            write( *, * ) "Error opening [", file_name, "] in writeDataGroupArray!"
            stop
        endif
        !
    end subroutine writeDataGroupArray
    !
    !> No subroutine briefing
    subroutine writePredictedDataHeader( receiver, receiver_type )
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
    !> No subroutine briefing
    subroutine writeEsolutionHeader( nTx, nMode )
        implicit none
        !
        integer, intent( in ) :: nTx, nMode
        integer :: ios
        character(len=20) :: version
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
    !> No subroutine briefing
    subroutine printUsage()
        implicit none
        !
        write( *, * ) "ModEM-OO Minimal Usage:"
        write( *, * ) ""
        write( *, * ) "    Forward Modeling:"
        write( *, * ) "        <ModEM> -f -m <rFile_Model> -d <rFile_Data>"
        write( *, * ) "    Outputs:"
        write( *, * ) "        - predicted_data.dat or the path specified by [-pd]"
        write( *, * ) "        - esolution.bin or the path specified by      [-es]"
        write( *, * ) ""
        write( *, * ) "    Adjoint:"
        write( *, * ) "        <ModEM> -j -m <rFile_Model> -pm <rFile_pModel> -d <rFile_Data>"
        write( *, * ) "    Outputs:"
        write( *, * ) "        - JmHat.dat or the path specified by  [-gd]"
        write( *, * ) "        - predicted_data.dat or the path specified by [-pd]"
        write( *, * ) "        - esolution.bin or the path specified by      [-es]"
        write( *, * ) ""
        write( *, * ) "    JMult_t:"
        write( *, * ) "        <ModEM> -jt -m <rFile_Model> -d <rFile_Data>"
        write( *, * ) "    Output:"
        write( *, * ) "        - sigma0.mod or the path specified by          [-dm]"
        write( *, * ) "        - predicted_data.dat or the path specified by [-pd]"
        write( *, * ) "        - esolution.bin or the path specified by      [-es]"
        !
    end subroutine printUsage
    !
    !> No subroutine briefing
    subroutine printHelp()
        implicit none
        !
        write( *, * ) "ModEM-OO Options:"
        write( *, * ) ""
        write( *, * ) "    Flags to define a job:"
        write( *, * ) "        [-f],  [--forward]    :  Forward Modeling."
        write( *, * ) "        [-j],  [--adjoint]    :  Adjoint."
        write( *, * ) "        [-jt], [--JMult_t]  :  Transposed Adjoint."
        write( *, * ) "        [-i],  [--inversion]  :  Inversion."
        write( *, * )
        write( *, * ) "    Other arguments:"
        write( *, * ) "        [-d],  [--data]      :  Flags for input data file path."
        write( *, * ) "        [-m],  [--model]     :  Flags for input model file path."
        write( *, * ) "        [-pm], [--pmodel]    :  Flags for input prior model file path."
        write( *, * ) "        [-c],  [--control]   :  Flags for user control file path."
        write( *, * ) "        [-v],  [--version]   :  Print version."
        write( *, * ) "        [-h],  [--help]      :  Print usage information."
        write( *, * ) "        [-dm], [--dmodel]    :  Output dsigma model file path."
        write( *, * ) "        [-pd], [--predicted] :  Output predicted data file path."
        write( *, * ) "        [-gd], [--gradient]  :  Output data gradient file path."
        write( *, * ) "        [-es], [--esolution] :  Output binary e-solution file path."
        write( *, * ) "        [-vb], [--verbose]   :  Print runtime information."
        !
        write( *, * ) ""
        write( *, * ) "Version 1.0.0"
        !
    end subroutine printHelp
    !
end program TestMPI
!