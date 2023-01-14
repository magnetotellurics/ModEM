!
!> No module briefing
!
module MasterMPI
    !
    use DeclarationMPI
    !
    public :: masterJobForwardModelling
    public :: masterForwardModelling
    public :: broadcastFwdBuffer
    !
contains
    !
    !> Routine to run a full ForwardModeling job and deliver the result (PredictedData) in a text file
    subroutine masterJobForwardModelling()
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: all_predicted_data
        !
        ! Verbose
        write( *, * ) "     - Start jobForwardModeling"
        !
        !> Read Model File and instantiate global variables: main_grid, model_operator and Sigma0
        if( has_model_file ) then
            !
            call handleModelFile()
            !
        else
            stop "Error: masterJobForwardModelling > Missing Model file!"
        endif
        !
        !> Read Data File: instantiate Txs and Rxs and build the Data relation between them
        if( has_data_file ) then
            !
            call handleDataFile()
            !
        else
            stop "Error: masterJobForwardModelling > Missing Data file!"
        endif
        !
        !> Send Fwd components to all workers
        call broadcastFwdBuffer()
        !
        !> Deallocate Forward Modeling global components (Not used by the Master process)
        deallocate( model_operator, sigma0, main_grid )
        !
        !> Run masterForwardModelling to calculate predicted data
        call masterForwardModelling( all_predicted_data )
        !
        !> Write all_predicted_data, with its proper Rx headers, to the file <predicted_data_file_name>
        call writeDataGroupTxArray( all_predicted_data, predicted_data_file_name )
        !
        !> Deallocate local data array
        call deallocateDataGroupTxArray( all_predicted_data )
        !
        ! Verbose
        write( *, * ) "     - Finish jobForwardModeling"
        !
    end subroutine masterJobForwardModelling
    !
    !> No procedure briefing
    subroutine masterForwardModelling( all_predicted_data )
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: all_predicted_data
        !
        integer :: worker_rank, tx_received, i_tx, i_data
        !
        !> Verbose
        write( *, * ) "     - Start Forward Modeling"
        !
        all_predicted_data = all_measured_data
        !
        call zerosDataGroupTxArray( all_predicted_data )
        !
        !> Initialize MPI control variables
        worker_rank = 1
        tx_received = 0
        i_tx = 0
        !
        !> Send 1 transmitter to first np workers
        do while( worker_rank <= ( mpi_size - 1 ) )
            !
            i_tx = i_tx + 1
            !
            job_info%job_name = job_forward
            job_info%i_tx = i_tx
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
            i_tx = i_tx + 1
            !
            job_info%job_name = job_forward
            job_info%i_tx = i_tx
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
            job_info%job_name = job_finish
            !
            call sendTo( job_info%worker_rank )
            !
        enddo
        !
    end subroutine masterForwardModelling
    !
    !> Routine to run a full ForwardModeling job and deliver the result (PredictedData) in a text file
    subroutine masterJobJMult()
        implicit none
        !
        !> Data gradient for all transmitters, grouped into an array of DataGroupTx
        type( DataGroupTx_t ), allocatable, dimension(:) :: JmHat
        !
        integer :: worker_rank, i_tx, tx_received
        !
        ! Verbose
        write( *, * ) "     - Start jobJMult"
        !
        !> Read Model File: instantiates Grid, ModelOperator and ModelParameter
        if( .NOT. has_model_file ) then 
            stop "Error: masterJobJMult > Missing Model file!"
        else
            !
            call handleModelFile()
            !
        endif
        !
        !> Read Prior Model File: instantiates PModelParameter
        if( .NOT. has_pmodel_file ) then 
            stop "Error: masterJobJMult > Missing Prior Model file!"
        else
            !
            call handlePModelFile()
            !
        endif
        !
        !> Reads Data File: instantiates and builds the Data relation between Txs and Rxs
        if( .NOT. has_data_file ) then 
            stop "Error: masterJobJMult > Missing Data file!"
        else
            !
            call handleDataFile()
            !
        endif
        !
        !> Send Fwd components to all workers
        call broadcastFwdBuffer()
        !
        !> Deallocate Forward Modeling global components
        deallocate( sigma0, pmodel, model_operator, main_grid )
        !
        !>
        JmHat = all_measured_data
        !
        !> Initialize MPI control variables
        worker_rank = 1
        tx_received = 0
        i_tx = 0
        !
        !> Send 1 transmitter to first np workers
        do while( worker_rank <= ( mpi_size - 1 ) )
            !
            i_tx = i_tx + 1
            !
            job_info%job_name = job_jmult
            job_info%i_tx = i_tx
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
            call receiveData( JmHat( job_info%i_tx ), job_info%worker_rank )
            !
            tx_received = tx_received + 1
            i_tx = i_tx + 1
            !
            job_info%job_name = job_jmult
            job_info%i_tx = i_tx
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
        !> Write JmHat, with its proper Rx headers, to the file <jmhat_data_file_name>
        call writeDataGroupTxArray( JmHat, jmhat_data_file_name )
        !
        !> Deallocate local data array
        call deallocateDataGroupTxArray( JmHat )
        !
        !> Verbose
        write( *, * ) "     - Finish jobJMult"
        !
    end subroutine masterJobJMult
    !
    !> Send Fwd components to all workers
    subroutine broadcastFwdBuffer()
        implicit none
        !
        integer :: worker_id
        !
        call allocateFwdBuffer()
        !
        call packFwdBuffer()
        !
        do worker_id = 1, ( mpi_size - 1 )
            !
            job_info%job_name = job_share_memory
            job_info%buffer_size = fwd_buffer_size
            !
            call sendTo( worker_id )
            !
            call sendFwdBuffer( worker_id )
            !
        enddo
        !
        call MPI_BARRIER( main_comm, ierr )
        !
        deallocate( fwd_buffer )
        !
    end subroutine broadcastFwdBuffer
    !
end module MasterMPI
!