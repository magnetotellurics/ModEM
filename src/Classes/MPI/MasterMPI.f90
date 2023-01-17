!
!> No module briefing
!
module MasterMPI
    !
    use DeclarationMPI
    !
    public :: masterForwardModelling
    public :: masterJMult
    public :: masterJMult_T
    public :: broadcastBasicComponents
    public :: broadcastSigma
    public :: broadcastFinish
    !
contains
    !
    !> No procedure briefing
    subroutine masterSolveAll( sigma )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma
        !
        integer :: worker_rank, tx_received, i_tx, i_data
        !
        !> Verbose
        !write( *, * ) "     - Start masterSolveAll"
        !
        !> Send Sigma to all workers
        if( sigma%is_allocated ) then
            !
            call broadcastSigma( sigma )
            !
        else
            stop "Error: masterSolveAll > sigma not allocated"
        endif
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
            job_info%job_name = job_em_solve
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
            tx_received = tx_received + 1
            i_tx = i_tx + 1
            !
            job_info%job_name = job_em_solve
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
            tx_received = tx_received + 1
            !
        enddo
        !
        !> Verbose
        !write( *, * ) "     - Finish masterSolveAll"
        !
    end subroutine masterSolveAll
    !
    !> No procedure briefing
    subroutine masterForwardModelling( sigma, all_predicted_data )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: all_predicted_data
        !
        integer :: worker_rank, tx_received, i_tx, i_data
        !
        !> Verbose
        !write( *, * ) "     - Start masterForwardModelling"
        !
        !>
        call masterSolveAll( sigma )
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
            job_info%data_size = allocateDataBuffer( all_predicted_data( i_tx ) )
            !
            call sendTo( worker_rank )
            !
            call sendData( all_predicted_data( i_tx ), worker_rank )
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
            job_info%data_size = allocateDataBuffer( all_predicted_data( i_tx ) )
            !
            call sendTo( worker_rank )
            !
            call sendData( all_predicted_data( i_tx ), worker_rank )
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
        enddo
        !
        !> Verbose
        !write( *, * ) "     - Finish masterForwardModelling"
        !
    end subroutine masterForwardModelling
    !
    !> Routine to run a full JMult job in parallel and deliver the result (JmHat) in a text file
    subroutine masterJMult( sigma, dsigma, JmHat )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma, dsigma
        type( DataGroupTx_t ), dimension(:), intent( inout ) :: JmHat
        !
        integer :: worker_rank, i_tx, tx_received
        !
        !> Verbose
        !write( *, * ) "     - Start masterJMult"
        !
        !> Send Sigma to all workers
        if( dsigma%is_allocated ) then
            !
            call broadcastDSigma( dsigma )
            !
        else
            stop "Error: masterSolveAll > sigma not allocated"
        endif
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
            job_info%data_size = allocateDataBuffer( JmHat( i_tx ) )
            !
            call sendTo( worker_rank )
            !
            call sendData( JmHat( i_tx ), worker_rank )
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
            job_info%data_size = allocateDataBuffer( JmHat( i_tx ) )
            !
            call sendTo( job_info%worker_rank )
            !
            call sendData( JmHat( i_tx ), worker_rank )
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
        enddo
        !
        !> Verbose
        !write( *, * ) "     - Finish masterJMult"
        !
    end subroutine masterJMult
    !
    !> Routine to run JMult_T in parallel
    subroutine masterJMult_T( sigma, all_data, dsigma )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma
        type( DataGroupTx_t ), dimension(:), intent( in ) :: all_data
        class( ModelParameter_t ), allocatable, intent( out ) :: dsigma
        !
        class( Scalar_t ), allocatable :: tx_model_cond
        !
        integer :: worker_rank, i_tx, tx_received
        !
        !> Verbose
        write( *, * ) "     - Start masterJMult_T"
        !
        !> And initialize dsigma with Zeros
        if( sigma%is_allocated ) then
            !
            if( allocated( dsigma ) ) deallocate( dsigma )
            allocate( dsigma, source = sigma )
            !
            call dsigma%zeros()
            !
        else
            stop "Error: masterJMult_T > sigma not allocated"
        endif
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
            job_info%job_name = job_jmult_t
            job_info%i_tx = i_tx
            job_info%worker_rank = worker_rank
            job_info%data_size = allocateDataBuffer( all_data( i_tx ) )
            !
            call sendTo( worker_rank )
            !
            call sendData( all_data( i_tx ), job_info%worker_rank )
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
            call dsigma%getCond( tx_model_cond )
            !
            call receiveConductivity( tx_model_cond, job_info%worker_rank )
            !
            call dsigma%addCond( tx_model_cond )
            !
            deallocate( tx_model_cond )
            !
            tx_received = tx_received + 1
            i_tx = i_tx + 1
            !
            job_info%job_name = job_jmult_t
            job_info%i_tx = i_tx
            job_info%data_size = allocateDataBuffer( all_data( i_tx ) )
            !
            call sendTo( job_info%worker_rank )
            !
            call sendData( all_data( i_tx ), job_info%worker_rank )
            !
        enddo
        !
        !> Receive job_done from each worker
        do while( tx_received < size( transmitters ) )
            !
            call receiveFromAny()
            !
            call dsigma%getCond( tx_model_cond )
            !
            call receiveConductivity( tx_model_cond, job_info%worker_rank )
            !
            call dsigma%addCond( tx_model_cond )
            !
            deallocate( tx_model_cond )
            !
            tx_received = tx_received + 1
            !
        enddo
        !
        !> Verbose
        !write( *, * ) "     - Finish masterJMult_T"
        !
    end subroutine masterJMult_T
    !
    !> Send FWD components to all workers
    subroutine broadcastBasicComponents()
        implicit none
        !
        integer :: worker_id
        !
        job_info%job_name = job_basic_components
        !
        job_info%basic_comp_size = allocateBasicComponentsBuffer()
        !
        do worker_id = 1, ( mpi_size - 1 )
            !
            call sendTo( worker_id )
            !
            call sendBasicComponents( worker_id )
            !
        enddo
        !
        call MPI_BARRIER( main_comm, ierr )
        !
        deallocate( basic_comp_buffer )
        !
    end subroutine broadcastBasicComponents
    !
    !> Send sigma to all workers
    subroutine broadcastSigma( sigma )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma
        !
        integer :: worker_id
        !
        job_info%job_name = job_sigma_model
        !
        job_info%model_size = allocateModelBuffer( sigma, .FALSE. )
        !
        !write( *, "(A45, i8)" ) "Sigma = ", job_info%model_size
        !
        do worker_id = 1, ( mpi_size - 1 )
            !
            call sendTo( worker_id )
            !
            call sendModel( sigma, worker_id )
            !
        enddo
        !
        call MPI_BARRIER( main_comm, ierr )
        !
        deallocate( model_buffer )
        !
    end subroutine broadcastSigma
    !
    !> Send dsigma to all workers
    subroutine broadcastDSigma( dsigma )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: dsigma
        !
        integer :: worker_id
        !
        job_info%job_name = job_dsigma_model
        !
        job_info%model_size = allocateModelBuffer( dsigma, .FALSE. )
        !
        !write( *, "(A45, i8)" ) "DSigma = ", job_info%model_size
        !
        do worker_id = 1, ( mpi_size - 1 )
            !
            call sendTo( worker_id )
            !
            call sendModel( dsigma, worker_id )
            !
        enddo
        !
        call MPI_BARRIER( main_comm, ierr )
        !
        deallocate( model_buffer )
        !
    end subroutine broadcastDSigma
    !
    !> Send sigma model to all workers
    subroutine broadcastFinish
        implicit none
        !
        integer :: worker_id
        !
        job_info%job_name = job_finish
        !
        do worker_id = 1, ( mpi_size - 1 )
            !
            call sendTo( worker_id )
            !
        enddo
        !
    end subroutine broadcastFinish
    !
end module MasterMPI
!