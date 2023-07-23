!
!> Module with the ForwardModeling and Sensitivity parallel routines
!> And routines to broadcast components
!
module MasterMPI
    !
    use DeclarationMPI
    !
    public :: masterSolveAll
    public :: masterForwardModelling
    public :: masterJMult
    public :: masterJMult_T
    public :: broadcastBasicComponents
    public :: broadcastSigma
    public :: broadcastFinish
    !
contains
    !
    !> Calculate ESolution in parallel for all transmitters.
    !
    subroutine masterSolveAll( sigma )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma
        !
        integer :: worker_rank, tx_received, i_tx
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
            call errStop( "masterSolveAll > sigma not allocated" )
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
            call receiveFromAny
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
            call receiveFromAny
            !
            tx_received = tx_received + 1
            !
        enddo
        !
        !> Verbose
        !!write( *, * ) "     - Finish masterSolveAll"
        !
    end subroutine masterSolveAll
    !
    !> Calculate in parallel ESolution for all transmitters
    !> Calculate in parallel the predicted data for each transmitter-receiver pair.
    !
    subroutine masterForwardModelling( sigma, all_predicted_data, i_sol )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma
        type( DataGroupTx_t ), allocatable, dimension(:), intent( out ) :: all_predicted_data
        integer, intent( in ), optional :: i_sol
        !
        integer :: worker_rank, tx_received, i_tx, sol_index
        !
        !> Verbose
        !write( *, * ) "     - Start masterForwardModelling"
        !
        sol_index = 0
        !
        !> Set i_sol if present
        if( present( i_sol ) ) sol_index = i_sol
        !
        !> Send Sigma to all workers
        if( sigma%is_allocated ) then
            !
            call broadcastSigma( sigma )
            !
        else
            call errStop( "masterForwardModelling > sigma not allocated" )
        endif
        !
        if( allocated( all_predicted_data ) ) deallocate( all_predicted_data )
        !
        all_predicted_data = all_measured_data
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
            job_info%sol_index = sol_index
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
            call receiveFromAny
            !
            call receiveData( all_predicted_data( job_info%i_tx ), job_info%worker_rank )
            !
            tx_received = tx_received + 1
            i_tx = i_tx + 1
            !
            job_info%job_name = job_forward
            job_info%i_tx = i_tx
            job_info%sol_index = sol_index
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
            call receiveFromAny
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
    !> Calculate JmHat in parallel for all transmitters
    !
    subroutine masterJMult( sigma, dsigma, JmHat )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma, dsigma
        type( DataGroupTx_t ), allocatable, dimension(:), intent( out ) :: JmHat
        !
        integer :: worker_rank, i_tx, tx_received
        !
        !> Send Sigma to all workers
        if( sigma%is_allocated ) then
            !
            call broadcastSigma( sigma )
            !
        else
            call errStop( "masterJMult > sigma not allocated" )
        endif
        !
        !> Send dSigma to all workers
        if( dsigma%is_allocated ) then
            !
            call broadcastDSigma( dsigma )
            !
        else
            call errStop( "masterJMult > sigma not allocated" )
        endif
        !
        if( allocated( JmHat ) ) deallocate( JmHat )
        !
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
            call receiveFromAny
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
            call receiveFromAny
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
    !> Calculate dsigma in parallel for all transmitters
    !
    subroutine masterJMult_T( sigma, all_data, dsigma, i_sol, s_hat )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma
        type( DataGroupTx_t ), dimension(:), intent( in ) :: all_data
        class( ModelParameter_t ), allocatable, intent( out ) :: dsigma
        integer, intent( in ), optional :: i_sol
        class( ModelParameter_t ), allocatable, dimension(:), intent( out ), optional :: s_hat
        !
        class( ModelParameter_t ), allocatable :: tx_dsigma
        !
        integer :: worker_rank, i_tx, tx_received, sol_index
        !
        !> Verbose
        !write( *, * ) "     - Start masterJMult_T"
        !
        sol_index = 0
        !
        !> Set i_sol if present
        if( present( i_sol ) ) sol_index = i_sol
        !
        !> And initialize dsigma with Zeros
        if( sigma%is_allocated ) then
            !
            !> Send Sigma to all workers
            call broadcastSigma( sigma )
            !
            if( allocated( dsigma ) ) deallocate( dsigma )
            allocate( dsigma, source = sigma )
            !
            call dsigma%zeros
            !
        else
            call errStop( "masterJMult_T > sigma not allocated" )
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
            job_info%sol_index = sol_index
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
        !> Allocate s_hat array
        if( present( s_hat ) ) then
            allocate( ModelParameterCell_SG_t :: s_hat( size( transmitters ) ) )
        endif
        !
        !> Send 1 transmitter to first available worker
        do while( i_tx < size( transmitters ) )
            !
            call receiveFromAny
            !
            call receiveModel( tx_dsigma, job_info%worker_rank )
            !
            if( present( s_hat ) ) then
                s_hat( job_info%i_tx ) = tx_dsigma
            endif
            !
            call dsigma%linComb( ONE, ONE, tx_dsigma )
            !
            tx_received = tx_received + 1
            i_tx = i_tx + 1
            !
            job_info%job_name = job_jmult_t
            job_info%i_tx = i_tx
            job_info%sol_index = sol_index
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
            call receiveFromAny
            !
            call receiveModel( tx_dsigma, job_info%worker_rank )
            !
            if( present( s_hat ) ) then
                s_hat( job_info%i_tx ) = tx_dsigma
            endif
            !
            call dsigma%linComb( ONE, ONE, tx_dsigma )
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
    !> Send FWD basic components to all workers
    !
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
        call MPI_BARRIER( MPI_COMM_WORLD, ierr )
        !
        deallocate( basic_comp_buffer )
        !
    end subroutine broadcastBasicComponents
    !
    !> Send sigma to all workers
    !
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
        call MPI_BARRIER( MPI_COMM_WORLD, ierr )
        !
        deallocate( model_buffer )
        !
    end subroutine broadcastSigma
    !
    !> Send dsigma to all workers
    !
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
        !write( *, "(A45, i8)" ) "dsigma = ", job_info%model_size
        !
        do worker_id = 1, ( mpi_size - 1 )
            !
            call sendTo( worker_id )
            !
            call sendModel( dsigma, worker_id )
            !
        enddo
        !
        call MPI_BARRIER( MPI_COMM_WORLD, ierr )
        !
        deallocate( model_buffer )
        !
    end subroutine broadcastDSigma
    !
    !> Send sigma model to all workers
    !
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