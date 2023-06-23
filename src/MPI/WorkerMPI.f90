!
!> Module routines to perform Forward Modeling and Sensitivity operations
!> Always relying on a single transmitter specified by the master process 
!
module WorkerMPI
    !
    use DeclarationMPI
    !
    use InversionJob
    !
    public :: workerMainLoop
    public :: workerSolve
    public :: workerForwardModelling
    public :: workerJMult
    public :: workerJMult_T
    public :: handleBasicComponents
    public :: handleSigmaModel
    public :: handleDSigmaModel
    public :: setTxForwardSolver
    !
contains
    !
    !> Message waiting loop.
    !> To perform a task specified by the master process.
    !
    subroutine workerMainLoop()
        implicit none
        !
        real( kind=prec ) :: t_start, t_finish
        integer :: int_time
        !
        call cpu_time( t_start )
        !
        !> Receive until master process send finish message
        do while ( job_master .NE. job_finish )
            !
            call receiveFrom( master_id )
            !
            job_master = trim( job_info%job_name )
            !
            select case( job_master )
                !
                case( "HANDLE_FWD_COMP" )
                    !
                    call handleBasicComponents
                !
                case( "HANDLE_SIGMA" )
                    !
                    call handleSigmaModel
                !
                case( "HANDLE_DSIGMA" )
                    !
                    call handleDSigmaModel
                    !
                case( "JOB_EM_SOLVE" )
                    !
                    call workerSolve
                    !
                case( "JOB_FORWARD" )
                    !
                    call workerForwardModelling
                    !
                case( "JOB_JMULT" )
                    !
                    call workerJMult
                    !
                case( "JOB_JMULT_T" )
                    !
                    call workerJMult_T
                    !
            end select
            !
        enddo
        !
        if( allocated( sigma ) ) deallocate( sigma )
        !
        if( allocated( dsigma ) ) deallocate( dsigma )
        !
        call garbageCollector
        !
        call cpu_time( t_finish )
        !
        int_time = int( t_finish - t_start )
        !
        write( *, "( a25, i5, a50 )" )  "Worker", mpi_rank, " finished: "//getLiteralTime( int_time )
        !
        call MPI_Finalize( ierr )
        !
    end subroutine workerMainLoop
    !
    !> Calculate E_Solution (e_sol_0) for the process single Transmitter
    !
    subroutine workerSolve()
        implicit none
        !
        class( Transmitter_t ), pointer :: Tx
        !
        !> Point to the transmitter specified by the master process 
        Tx => getTransmitter( job_info%i_tx )
        !
        call setTxForwardSolver( Tx )
        !
        call solveTx( sigma, Tx )
        !
        !> Send job done to master
        job_info%job_name = job_done
        job_info%worker_rank = mpi_rank
        !
        call sendTo( master_id )
        !
    end subroutine workerSolve
    !
    !> Receive data_tx from master process.
    !>     Calculate the predicted data between the data_tx's transmitter and its receivers
    !> Send data_tx to master process
    !> Require previous call of workerSolve or masterSolveAll
    !
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
        !> Point to the transmitter specified by the master process 
        Tx => getTransmitter( job_info%i_tx )
        !
        Tx%i_sol = job_info%sol_index
        !
        call setTxForwardSolver( Tx )
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
        !> Send job done and tx_data to master process
        job_info%job_name = job_done
        job_info%worker_rank = mpi_rank
        job_info%data_size = allocateDataBuffer( tx_data )
        !
        call sendTo( master_id )
        !
        call sendData( tx_data, master_id )
        !
    end subroutine workerForwardModelling
    !
    !> Receive data_tx from master process.
    !> Calculate JmHat for the data_tx's transmitter:
    !>     Set the transmitter's source by calling PMult.
    !>     Solve e_sens for the transmitter.
    !> Send dsigma%cell_cond_h to master process
    !> Require previous call of workerSolve or masterSolveAll
    !
    subroutine workerJMult()
        implicit none
        !
        class( Transmitter_t ), pointer :: Tx
        type( DataGroupTx_t ) :: tx_data
        !
        call receiveData( tx_data, master_id )
        !
        !> Point to the transmitter specified by the master process 
        Tx => getTransmitter( job_info%i_tx )
        !
        !> Switch Transmitter's source to SourceInteriorForce from PMult
        call Tx%setSource( Tx%PMult( sigma, dsigma, model_operator ) )
        !
        !> Solve e_sens with the new Source
        call Tx%solve
        !
        !> serialJMult for the same Tx
        call JMult_Tx( tx_data )
        !
        !> Send job done and tx_data to master process
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
    !> Receive data_tx from master process
    !> Calculate tx_dsigma for the data_tx's transmitter with JMult_T_Tx
    !>     Create a rhs from LRows * residual data for all receivers related to the transmitter.
    !>     Solve ESens on the transmitter using a transpose SourceInteriorForce, with the new rhs.
    !>     Call Tx%PMult to get a new ModelParameter dsigma.
    !> Send dsigma%cell_cond_h to master process
    !> Require previous call of workerSolve or masterSolveAll
    !
    subroutine workerJMult_T()
        implicit none
        !
        class( ModelParameter_t ), allocatable :: tx_dsigma
        type( DataGroupTx_t ) :: tx_data
        class( Transmitter_t ), pointer :: Tx
        !
        call receiveData( tx_data, master_id )
        !
        !> Point to the transmitter specified by the master process 
        Tx => getTransmitter( job_info%i_tx )
        !
        Tx%i_sol = job_info%sol_index
        !
        call JMult_T_Tx( sigma, tx_data, tx_dsigma, Tx%i_sol )
        !
        !> Send job done and tx_dsigma's conductivity to master process
        job_info%job_name = job_done
        job_info%worker_rank = mpi_rank
        !
        call sendTo( master_id )
        !
        job_info%model_size = allocateModelBuffer( tx_dsigma, .FALSE. )
        !
        call sendModel( tx_dsigma, master_id )
        !
    end subroutine workerJMult_T
    !
    !> Receive from master process and deals with 
    !> the basic components for all Forward Modeling and Sensitivity operations:
    !>     Control variables, main_grid, transmitters and receivers
    !> Instantiate model_operator and set its equations
    !
    subroutine handleBasicComponents()
        implicit none
        !
        call receiveBasicComponents( master_id )
        !
        call MPI_BARRIER( MPI_COMM_WORLD, ierr )
        !
        !write( *, "( a30, i5, a11, i10, a11 )" ) "Worker", mpi_rank, " Received: ", job_info%basic_comp_size, " bytes."
        !
        select type( main_grid )
            !
            class is( Grid3D_SG_t )
                !
                allocate( model_operator, source = ModelOperator_MF_t( main_grid ) )
                !
                call model_operator%setEquations
                !
            class default
                stop "Error: handleBasicComponents > Unclassified main_grid"
            !
        end select
        !
    end subroutine handleBasicComponents
    !
    !> Receive from master process and deals with sigma model
    !> Set its metric and model_operator's conductivity 
    !
    subroutine handleSigmaModel()
        implicit none
        !
        call receiveModel( sigma, master_id )
        !
        call MPI_BARRIER( MPI_COMM_WORLD, ierr )
        !
        !write( *, "( a30, i5, a11, i10, a11 )" ) "Worker", mpi_rank, " Received: ", job_info%model_size, " bytes."
        !
        call sigma%setMetric( model_operator%metric )
        !
        call model_operator%setCond( sigma )
        !
    end subroutine handleSigmaModel
    !
    !> Receive dsigma model from master process and set its metric
    !
    subroutine handleDSigmaModel()
        implicit none
        !
        call receiveModel( dsigma, master_id )
        !
        call MPI_BARRIER( MPI_COMM_WORLD, ierr )
        !
        !write( *, "( a30, i5, a11, i10, a11 )" ) "Worker", mpi_rank, " Received: ", job_info%model_size, " bytes."
        !
        call dsigma%setMetric( model_operator%metric )
        !
    end subroutine handleDSigmaModel
    !
    !> Create the global ForwardSolver a single transmitter to it.
    !
    subroutine setTxForwardSolver( Tx )
        implicit none
        !
        class( Transmitter_t ), pointer, intent( inout ) :: Tx
        !
        if( allocated( forward_solver ) ) deallocate( forward_solver )
        !
        !> Instantiate the ForwardSolver - Specific type can be chosen via control file
        select case( forward_solver_type )
            !
            case( FWD_IT_DC )
                !
                allocate( forward_solver, source = ForwardSolverIT_DC_t( model_operator, QMR ) )
                !
            case( "" )
                !
                write( *, * ) "     "//achar(27)//"[91m# Warning:"//achar(27)//"[0m setTxForwardSolver > Forward Solver type not provided, using IT_DC."
                !
                allocate( forward_solver, source = ForwardSolverIT_DC_t( model_operator, QMR ) )
                !
            case default
                !
                write( *, * ) "Error: Wrong Forward Solver type: [", forward_solver_type, "]"
                stop
                !
        end select
        !
        Tx%forward_solver => forward_solver
        !
    end subroutine setTxForwardSolver
    !
end module WorkerMPI
!