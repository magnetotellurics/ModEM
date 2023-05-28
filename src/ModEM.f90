!
!> Prototype of a full ModEM serial|parallel program
!
program ModEM
    !
#ifdef MPI
    !
    use WorkerMPI
    !
    call constructorMPI
    !
    ! MPI Master process
    if( mpi_rank == 0 ) then
        !
        call runProgram
        !
        call MPI_Finalize( ierr )
        !
    ! MPI Worker process
    else
        !
        call workerMainLoop
        !
    endif
    !
#else
    !
    use InversionDCG
    use InversionNLCG
    !
    call runProgram
    !
#endif
    !
contains
    !
    !> Main program protocol
    !
    subroutine runProgram()
        implicit none
        !
        real( kind=prec ) :: t_start, t_finish
        integer :: int_time
        !
        call date_and_time( str_date, str_time )
        !
        !> Start runtime countdown
        call cpu_time( t_start )
        !
        modem_job = "unknown"
        !
        !> Set default value of program variables
        call setupDefaultParameters
        !
        !> Handle arguments passed by the user on the command line
        call handleArguments
        !
        write( *, * )
        write( *, "(a18, a8, a1, a6, a1)" ) "Start ModEM-OO at ", str_date, "_", str_time, "."
        write( *, * )
        !
        !> If it was passed by argument,
        !> Check parameters at the forward control file
        if( has_fwd_control_file ) call handleForwardControlFile
        !
        !> If it was passed by argument,
        !> Check parameters at the inversion control file
        if( has_inv_control_file ) call handleInversionControlFile
        !
        !> Execute the job specified in the arguments
        call handleJob
        !
        !> Deallocate remaining main program memory
        call garbageCollector
        !
        !> End runtime countdown
        call cpu_time( t_finish )
        !
        int_time = int( t_finish - t_start )
        !
        write( *, * )
        write( *, * ) "Finish ModEM-OO: ", getLiteralTime( int_time )
        write( *, * )
        !
    end subroutine runProgram
    !
    !> Routine to run a full Inversion Job - Minimize Residual data error
    !> Where:
    !>     sigma  = Input|Start model (for predicted data and final inversion model)
    !>     pmodel = Perturbation model (if exist -dm read input model)
    !>     sigma0 = Read input model (-m)
    !>     dsigma = Production model (data gradient From serialJMult_T)
    !
    subroutine jobInversion()
        implicit none
        !
        class( ModelParameter_t ), allocatable :: sigma, dsigma
        !
        class( Inversion_t ), allocatable :: inversion
        !
        ! Verbose
        write( *, * ) "     - Start jobInversion"
        !
        !> Read Model File and instantiate global variables: main_grid, model_operator and Sigma0
        if( has_model_file ) then 
            !
            call handleModelFile( sigma )
            !
            !> Instantiate ModelCovariance
            allocate( model_cov, source = ModelCovarianceRec_t( sigma ) )
            !
            if( has_cov_file ) then 
                !
                call model_cov%read_CmSqrt( cov_file_name )
                !
                write( *, * ) "     < Cov File: [", cov_file_name, "]"
                !
            else
                write( *, * ) "     "//achar(27)//"[91m# Warning:"//achar(27)//"[0m jobInversion > Missing Covariance file!"
            endif
            !
            !> Initialize pmodel with Zeros
            allocate( dsigma, source = sigma )
            !
            call dsigma%zeros
            !
        else
            stop "Error: jobInversion > Missing Model file!"
        endif
        !
        !> Read Perturbation Model File: instantiate pmodel (NOT USING RIGHT NOW ????)
        if( has_pmodel_file ) then 
            !
            deallocate( dsigma )
            !
            call handlePModelFile( dsigma )
            !
            call dsigma%setMetric( model_operator%metric )
            !
            call model_cov%multBy_Cm( dsigma )
            !
            call sigma%linComb( ONE, ONE, dsigma )
            !
        endif
        !
        !> Read Data File: instantiate and build the Data relation between Txs and Rxs
        if( has_data_file ) then 
            !
            call handleDataFile
            !
        else
            stop "Error: jobInversion > Missing Data file!"
        endif
        !
#ifdef MPI
        !
        call broadcastBasicComponents
        !
#else
        !
        call createDistributeForwardSolver
        !
#endif
        !
        !if( .NOT. has_inv_control_file ) then
            !inversion_type = NLCG
        !endif
        !
        !> Instantiate the ForwardSolver - Specific type can be chosen via control file
        select case( inversion_type )
            !
            case( DCG )
                !
                allocate( inversion, source = InversionDCG_t() )
                !
            case( NLCG )
                !
                allocate( inversion, source = InversionNLCG_t() )
                !
            case default
                !
                stop "Error: jobInversion > Undefined inversion_type"
                !
        end select
        !
        call inversion%solve( all_measured_data, sigma, dsigma )
        !
#ifdef MPI
        call broadcastFinish
#endif
        !
        deallocate( sigma, dsigma, inversion )
        !
        ! Verbose
        write( *, * ) "     - Finish jobInversion"
        !
    end subroutine jobInversion
    !
    !> No subroutine briefing
    !
    subroutine handleJob()
        implicit none
        !
        !> Except for the case of Inversion,
        !> Free the memory used by the global control file, which is no longer useful
        if( index( modem_job, "Inversion" ) .LE. 0 .AND. allocated( inv_control_file ) ) then
            deallocate( inv_control_file )
        endif
        !
        select case( modem_job )
            !
            case( "JobForwardModeling" )
                !
                call jobForwardModeling
                !
            case( "JobJMult" )
                !
                call jobJMult
                !
            case( "JobJMult_T" )
                !
                call jobJMult_T
                !
            case( "JobInversion" )
                !
                call jobInversion
                !
            case default
                !
                write( *, * ) "     "//achar(27)//"[31m# Error:"//achar(27)//"[0m unknown job: [", modem_job, "]"
                call printHelp
                stop
            !
        end select
        !
    end subroutine handleJob
    !
end program ModEM
!