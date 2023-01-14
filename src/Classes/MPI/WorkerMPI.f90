!
!> No module briefing
!
module WorkerMPI
    !
    use DeclarationMPI
    !
    public :: workerMainLoop
    public :: handleFwdBuffer
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
                case ( "SHARE_MEMORY" )
                    !
                    call handleFwdBuffer()
                    !
                case ( "JOB_FORWARD" )
                    !
                    call workerJobForwardModelling()
                    !
                case ( "JOB_JMULT" )
                    !
                    call workerJobJMult()
                    !
                case ( "JOB_JMULT_T" )
                    !
                    !call workerJobAdjoint_T()
                    !
                case ( "JOB_INVERSION" )
                    !
                    !call workerJobInversion()
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
    subroutine solveTx( sigma, Tx )
        implicit none
        !
        class( ModelParameter_t ), allocatable, intent( in ) :: sigma
        class( Transmitter_t ), pointer, intent( inout ) :: Tx
        !
        !> Instantiate the global ForwardSolver - Specific type can be chosen via control file
        call createForwardSolver()
        !
        !> Point to the current Transmitter
        Tx => getTransmitter( job_info%i_tx )
        !
        !> Set Transmitter's ForwardSolver Omega(Period) and Conductivity
        call Tx%forward_solver%setFrequency( sigma, Tx%period )
        !
        !> Instantiate Transmitter's Source - According to transmitter type or chosen via control file
        select type( Tx )
            !
            class is( TransmitterMT_t )
                !
                write( *, "( a25, i5, A14, i5, A9, es12.5 )" ) "- Worker", mpi_rank, " Solving MT Tx", Tx%i_tx, ", Period=", Tx%period
                !
                call Tx%setSource( SourceMT_1D_t( model_operator, sigma, Tx%period ) )
                !
            class is( TransmitterCSEM_t )
                !
                write( *, "( a25, i5, A14, i5, A9, es12.5 )" ) "- Worker", mpi_rank, " Solving CSEM Tx", Tx%i_tx, ", Period=", Tx%period
                !
                call Tx%setSource( SourceCSEM_Dipole1D_t( model_operator, sigma, Tx%period, Tx%location, Tx%dip, Tx%azimuth, Tx%moment ) )
                !
            class default
                stop "Error: workerJobForwardModelling > Unclassified Transmitter"
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
    !> No procedure briefing
    subroutine workerJobForwardModelling()
        implicit none
        !
        !class( ModelParameter_t ), allocatable, intent( in ) :: sigma
        !
        class( Transmitter_t ), pointer :: Tx
        class( Receiver_t ), pointer :: Rx
        type( DataGroupTx_t ) :: tx_data
        type( DataGroup_t ) :: data_group
        integer :: i_rx
        !
        call solveTx( sigma0, Tx )
        !
        !> Build the DataGroupTx to store data from a single transmitter.
        tx_data = DataGroupTx_t( Tx%i_tx )
        !
        !> Loop for each Receiver related to the Transmitter
        do i_rx = 1, size( Tx%receiver_indexes )
            !
            !> Point to the current Receiver
            Rx => getReceiver( Tx%receiver_indexes( i_rx ) )
            !
            !> Calculate and store Predicted Data and/or LRows in the Receiver
            !> Depending on the type of work
            call Rx%predictedData( Tx, data_group )
            !
            call tx_data%put( data_group )
            !
        enddo
        !
        !> SEND JOB DONE TO MASTER
        job_info%job_name = job_done
        job_info%worker_rank = mpi_rank
        job_info%i_tx = tx_data%i_tx
        !
        call sendTo( master_id )
        !
        call sendData( tx_data, master_id )
        !
    end subroutine workerJobForwardModelling
    !
    !> No procedure briefing
    subroutine workerJobJMult()
        implicit none
        !
        class( Transmitter_t ), pointer :: Tx
        class( Receiver_t ), pointer :: Rx
        type( DataGroupTx_t ) :: tx_data
        integer :: i
        !
        !>
        call solveTx( sigma0, Tx )
        !
        !> Switch Transmitter's source to SourceInteriorForce from PMult
        call Tx%setSource( Tx%PMult( sigma0, pmodel, model_operator ) )
        !
        !> Solve e_sens with the new Source
        call Tx%solve()
        !
        write( *, * ) "mpi_rank, TX ID: ", mpi_rank, Tx%i_tx
        !
        !> Build the DataGroupTx to store data from a single transmitter.
        tx_data = DataGroupTx_t( Tx%i_tx )
        !
        do i = 1, size( Tx%receiver_indexes )
            !
            Rx => getReceiver( Tx%receiver_indexes(i) )
            !
            call tx_data%put( DataGroup_t( Tx%receiver_indexes(i), Tx%i_tx, Rx%n_comp ) )
            !
        enddo
        !
        !> JMult for the same Tx
        call JMult_Tx( tx_data )
        !
        !> MPI: SEND JOB DONE TO MASTER
        job_info%job_name = job_done
        job_info%worker_rank = mpi_rank
        job_info%i_tx = tx_data%i_tx
        !
        call sendTo( master_id )
        !
        call sendData( tx_data, master_id )
        !
    end subroutine workerJobJMult
    !
    !> No procedure briefing
    subroutine handleFwdBuffer()
        implicit none
        !
        call receiveFwdBuffer( master_id )
        !
        call MPI_BARRIER( main_comm, ierr )
        !
        write( *, "( a25, i5, a11, i10, a11 )" ) "Worker", mpi_rank, " Received: ", fwd_buffer_size, " bytes."
        !
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
                stop "Error: handleFwdBuffer > Unclassified main_grid"
            !
        end select
        !
    end subroutine handleFwdBuffer
    !
end module WorkerMPI
!