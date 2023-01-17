!
!> No module briefing
!
module WorkerMPI
    !
    use DeclarationMPI
    !
    use InversionDCG
    use InversionNLCG
    !
    public :: workerMainLoop
    public :: solveTx
    public :: workerForwardModelling
    public :: workerJMult
    public :: workerJMult_T
    public :: handleBasicComponents
    public :: handleSigmaModel
    !
contains
    !
    !> ?????
    subroutine workerMainLoop()
        implicit none
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
                case ( "HANDLE_FWD_COMP" )
                    !
                    call handleBasicComponents()
                !
                case ( "HANDLE_SIGMA" )
                    !
                    call handleSigmaModel()
                !
                case ( "HANDLE_DSIGMA" )
                    !
                    call handleDSigmaModel()
                    !
                case ( "JOB_EM_SOLVE" )
                    !
                    call workerSolve()
                    !
                case ( "JOB_FORWARD" )
                    !
                    call workerForwardModelling()
                    !
                case ( "JOB_JMULT" )
                    !
                    call workerJMult()
                    !
                case ( "JOB_JMULT_T" )
                    !
                    call workerJMult_T()
                    !
            end select
            !
        enddo
        !
        !> Deallocate remaining worker memory
        call garbageCollector()
        !
        call cpu_time( t_finish )
        !
        write( *, "( a25, i5, a10, f8.3, a1 )" )  "Worker", mpi_rank, " finished:", t_finish - t_start, "s"
        !
        call MPI_Finalize( ierr )
        !
    end subroutine workerMainLoop
    !
    !> No procedure briefing
    subroutine workerSolve()
        implicit none
        !
        class( Transmitter_t ), pointer :: Tx
        !
        Tx => getTransmitter( job_info%i_tx )
        !
        call txForwardSolver( Tx )
        !
        call solveTx( sigma, Tx )
        !
        !> SEND JOB DONE TO MASTER
        job_info%job_name = job_done
        job_info%worker_rank = mpi_rank
        !
        call sendTo( master_id )
        !
    end subroutine workerSolve
    !
    !> No procedure briefing
    subroutine workerForwardModelling()
        implicit none
        !
        class( Transmitter_t ), pointer :: Tx
        class( Receiver_t ), pointer :: Rx
        type( DataGroupTx_t ) :: tx_data
        type( DataGroup_t ) :: data_group
        integer :: i_rx
        !
        call receiveData( tx_data, master_id )
        !
        Tx => getTransmitter( tx_data%i_tx )
        !
        call txForwardSolver( Tx )
        !
        call solveTx( sigma, Tx )
        !
        !> Loop for each Receiver related to the Transmitter
        do i_rx = 1, size( Tx%receiver_indexes )
            !
            Rx => getReceiver( Tx%receiver_indexes( i_rx ) )
            !
            call Rx%predictedData( Tx, data_group )
            !
            call tx_data%setValues( data_group )
            !
        enddo
        !
        !> SEND JOB DONE TO MASTER
        job_info%job_name = job_done
        job_info%worker_rank = mpi_rank
        job_info%i_tx = tx_data%i_tx
        job_info%data_size = allocateDataBuffer( tx_data )
        !
        call sendTo( master_id )
        !
        call sendData( tx_data, master_id )
        !
    end subroutine workerForwardModelling
    !
    !> No procedure briefing
    subroutine workerJMult()
        implicit none
        !
        class( Transmitter_t ), pointer :: Tx
        class( Receiver_t ), pointer :: Rx
        type( DataGroupTx_t ) :: tx_data
        integer :: i
        !
        call receiveData( tx_data, master_id )
        !
        Tx => getTransmitter( job_info%i_tx )
        !
        call txForwardSolver( Tx )
        !
        call solveTx( sigma, Tx )
        !
        !> Switch Transmitter's source to SourceInteriorForce from PMult
        call Tx%setSource( Tx%PMult( sigma, dsigma, model_operator ) )
        !
        !> Solve e_sens with the new Source
        call Tx%solve()
        !
        !> JMult for the same Tx
        call JMult_Tx( tx_data )
        !
        !> MPI: SEND JOB DONE TO MASTER
        job_info%job_name = job_done
        job_info%worker_rank = mpi_rank
        job_info%i_tx = Tx%i_tx
        job_info%data_size = allocateDataBuffer( tx_data )
        !
        call sendTo( master_id )
        !
        call sendData( tx_data, master_id )
        !
    end subroutine workerJMult
    !
    !> No procedure briefing
    subroutine workerJMult_T()
        implicit none
        !
        class( ModelParameter_t ), allocatable :: tx_dsigma
        class( Transmitter_t ), pointer :: Tx
        class( Receiver_t ), pointer :: Rx
        type( DataGroupTx_t ) :: tx_data
        integer :: i
        !
        call receiveData( tx_data, master_id )
        !
        Tx => getTransmitter( job_info%i_tx )
        !
        call txForwardSolver( Tx )
        !
        call solveTx( sigma, Tx )
        !
        !> Set current tx_dsigma from JMult_T_Tx
        call JMult_T_Tx( sigma, tx_data, tx_dsigma )
        !
        !> MPI: SEND JOB DONE TO MASTER
        job_info%job_name = job_done
        job_info%worker_rank = mpi_rank
        !
        call sendTo( master_id )
        !
        !> Send result model to the Master
        select type( tx_dsigma )
            !
            class is( ModelParameterCell_SG_t )
                !
                call sendConductivity( tx_dsigma%cell_cond, master_id )
                !
            class default
                stop "Error: workerJMult_T > Unclassified tx_dsigma"
            !
        end select
        !
    end subroutine workerJMult_T
    !
    !> No procedure briefing
    subroutine handleBasicComponents()
        implicit none
        !
        call receiveBasicComponents( master_id )
        !
        call MPI_BARRIER( main_comm, ierr )
        !
        !write( *, "( a30, i5, a11, i10, a11 )" ) "Worker", mpi_rank, " Received: ", job_info%basic_comp_size, " bytes."
        !
        select type( main_grid )
            !
            class is( Grid3D_SG_t )
                !
                allocate( model_operator, source = ModelOperator_MF_t( main_grid ) )
                !
                call model_operator%setEquations()
                !
            class default
                stop "Error: handleBasicComponents > Unclassified main_grid"
            !
        end select
        !
    end subroutine handleBasicComponents
    !
    !> No procedure briefing
    subroutine handleSigmaModel()
        implicit none
        !
        call receiveModel( sigma, master_id )
        !
        call MPI_BARRIER( main_comm, ierr )
        !
        !write( *, "( a30, i5, a11, i10, a11 )" ) "Worker", mpi_rank, " Received: ", job_info%model_size, " bytes."
        !
        call sigma%setMetric( model_operator%metric )
        !
        call model_operator%setCond( sigma )
        !
    end subroutine handleSigmaModel
    !
    !> No procedure briefing
    subroutine handleDSigmaModel()
        implicit none
        !
        call receiveModel( dsigma, master_id )
        !
        call MPI_BARRIER( main_comm, ierr )
        !
        !write( *, "( a30, i5, a11, i10, a11 )" ) "Worker", mpi_rank, " Received: ", job_info%model_size, " bytes."
        !
        call dsigma%setMetric( model_operator%metric )
        !
    end subroutine handleDSigmaModel
    !
    !> No subroutine briefing
    subroutine txForwardSolver( Tx )
        implicit none
        !
        class( Transmitter_t ), pointer, intent( inout ) :: Tx
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
                    !
                    write( *, * ) "Warning: txForwardSolver > Undefined forward_solver, using IT_DC"
                    !
                    allocate( forward_solver, source = ForwardSolverIT_DC_t( model_operator, QMR ) )
                    !
            end select
            !
        endif
        !
        Tx%forward_solver => forward_solver
        !
    end subroutine txForwardSolver
    !
end module WorkerMPI
!