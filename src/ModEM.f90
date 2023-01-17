!
!> Prototype of a full ModEM parallel program
!
program ModEM
    !
    use ModEMControlFile
    !
#ifdef MPI
    !
    use MasterMPI
    use WorkerMPI
    !
    call constructorMPI()
    !
    !> MPI MASTER PROCESS
    if( mpi_rank == 0 ) then
        !
        call startProgram()
        !
        call MPI_Finalize( ierr )
        !
    !> MPI WORKER PROCESS
    else
        !
        call workerMainLoop()
        !
    endif
    !
#else
    !
    use InversionDCG
    use InversionNLCG
    !
    call startProgram()
    !
#endif
    !
contains
    !
    !>
    subroutine startProgram()
        implicit none
        !
        real( kind=prec ) :: t_start, t_finish
        !
        !> Start the program and runtime count
        call cpu_time( t_start )
        !
        modem_job = "unknown"
        !
        call setupDefaultParameters()
        !
        !> Validate arguments, set model_file_name, data_file_name, control_file_name, etc...
        call handleArguments()
        !
        write( *, * )
        write( *, * ) "Start ModEM-OO."
        write( *, * )
        !
        !> Check parameters at the control file
        if( has_control_file ) call handleControlFile()
        !
        !> Execute the job specified in the arguments
        call handleJob()
        !
        !> Deallocate remaining main program memory
        call garbageCollector()
        !
        !> End runtime countdown
        call cpu_time( t_finish )
        !
        write( *, * )
        write( *, "( a16, f8.3, a1 )" ) "Finish ModEM-OO:", ( t_finish - t_start ), "s"
        write( *, * )
        !
    end subroutine startProgram
    !
    !> Routine to run a full Inversion Job - Minimize Residual
    !> Where:
    !>    SIGMA (M) = Production model (for predicted data and final inversion model)
    !>    PMODEL = Perturbation model  (if exist -dm readed input model)
    !>    SIGMA0 = Readed input model  (-m)
    !>    DSIGMA = From JMult_T        (????)
    subroutine jobInversion()
        implicit none
        !
        class( ModelParameter_t ), allocatable :: sigma, dsigma
        !
        !> Read Model File and instantiate global variables: main_grid, model_operator and Sigma0
        if( has_model_file ) then 
            !
            call handleModelFile( sigma )
            !
            !> Instantiate ModelCovariance
            allocate( model_cov, source = ModelCovarianceRec_t( sigma ) )
            !
            !> Initialize pmodel with Zeros
            allocate( dsigma, source = sigma )
            !
            call dsigma%zeros()
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
            call handleDataFile()
            !
        else
            stop "Error: jobInversion > Missing Data file!"
        endif
        !
#ifdef MPI
        !
        call broadcastBasicComponents()
        !
#else
        !
        call createDistributeForwardSolver()
        !
#endif
        !
        !> Instantiate the ForwardSolver - Specific type can be chosen via control file
        select case ( inversion_algorithm )
            !
            case( DCG )
                !
                call DCGsolver( all_measured_data, sigma, dsigma )
                !
            case( NLCG )
                !
                call NLCGsolver( all_measured_data, sigma, dsigma )
                !
            case default
                !
                stop "Error: jobInversion > Undefined inversion_algorithm"
                !
        end select
        !
        deallocate( sigma, dsigma )
        !
    end subroutine jobInversion
    !
    !> No subroutine briefing
    subroutine handleJob()
        implicit none
        !
        select case ( modem_job )
            !
            case ( "Forward" )
                !
                call jobForwardModeling
                !
            case ( "JMult" )
                !
                call jobJMult
                !
            case ( "JMult_T" )
                !
                call jobJMult_T
                !
            case ( "Inversion" )
                !
                call jobInversion
                !
            case default
                !
                write( *, * ) "Error: unknown job: [", modem_job, "]"
                call printHelp
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
        !> Read control file and sets the options in the Constants module
        control_file = ModEMControlFile_t( ioStartup, control_file_name )
        !
    end subroutine handleControlFile
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
                         modem_job = "Forward"
                         !
                         argument_index = argument_index + 1
                         !
                      case ( "-i", "--inversion" )
                         !
                         modem_job = "Inversion"
                         !
                         argument_index = argument_index + 1
                         !
                      case ( "-j", "--jmult" )
                         !
                         modem_job = "JMult"
                         !
                         argument_index = argument_index + 1
                         !
                      case ( "-jt", "--jmult_t" )
                         !
                         modem_job = "JMult_T"
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
                      case ( "-o", "--outdir" )
                         !
                         call get_command_argument( argument_index + 1, argument )
                         outdir_name = trim( argument )
                         !
                         if( len( outdir_name ) > 0 ) has_outdir_name = .TRUE.
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
                      case ( "-jm", "--jmhat" )
                         !
                         call get_command_argument( argument_index + 1, argument )
                         jmhat_data_file_name = trim( argument )
                         !
                         argument_index = argument_index + 2
                         !
                      case ( "-es", "--esolution" )
                         !
                         call get_command_argument( argument_index + 1, argument )
                         e_solution_file_name = trim( argument )
                         !
                         if( len( e_solution_file_name ) > 0 ) has_e_solution_file = .TRUE.
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
        predicted_data_file_name = "all_predicted_data.dat"
        jmhat_data_file_name = "jmhat.dat"
        e_solution_file_name = "esolution.bin"
        dsigma_file_name = "dsigma.rho"
        !
        inversion_algorithm = "DCG"
        !
        ! Control flags
        has_outdir_name = .FALSE.
        !
        has_control_file = .FALSE.
        has_model_file = .FALSE.
        has_pmodel_file = .FALSE.
        has_data_file = .FALSE.
        has_e_solution_file = .FALSE.
        verbosis = .FALSE.
        !
        ! Solver parameters
        QMR_iters = 40
        BCG_iters = 80
        max_divcor_calls = 20
        max_divcor_iters = 100
        tolerance_divcor = 1E-5
        tolerance_qmr = 1E-7
        !
        forward_solver_type = FWD_IT_DC
        !
        ! Source parameters
        source_type = SRC_MT_1D
        get_1D_from = "Geometric_mean"
        !
        ! Model parameters
        model_method = MM_METHOD_FIXED_H
        model_n_air_layer = 10
        model_max_height = 200.0
        !
    end subroutine setupDefaultParameters
    !
    !> No subroutine briefing
    subroutine printUsage()
        implicit none
        !
        write( *, * ) "ModEM Minimal Usage:"
        write( *, * ) ""
        write( *, * ) "    Forward Modeling:"
        write( *, * ) "        <ModEM> -f -m <rFile_Model> -d <rFile_Data>"
        write( *, * ) "        Output:"
        write( *, * ) "        - all_predicted_data.dat or the path specified by         [-pd]"
        write( *, * ) ""
        write( *, * ) "    JMult:"
        write( *, * ) "        <ModEM> -j -m <rFile_Model> -pm <rFile_pModel> -d <rFile_Data>"
        write( *, * ) "        Output:"
        write( *, * ) "        - jmhat.dat or the path specified by                  [-jm]"
        write( *, * ) ""
        write( *, * ) "    JMult_T:"
        write( *, * ) "        <ModEM> -jt -m <rFile_Model> -d <rFile_Data>"
        write( *, * ) "        Output:"
        write( *, * ) "        - dsigma.rho or the path specified by                 [-dm]"
        write( *, * ) ""
        write( *, * ) "    Inversion:"
        write( *, * ) "        <ModEM> -i -m <rFile_Model> -d <rFile_Data>"
        write( *, * ) "        Output:"
        write( *, * ) "        - directory named Output_<date>_<time> or specified by [-o]"
        !
    end subroutine printUsage
    !
    !> No subroutine briefing
    subroutine printHelp()
        implicit none
        !
        write( *, * ) "ModEM Options:"
        write( *, * ) ""
        write( *, * ) "    Flags to define a job:"
        write( *, * ) "        [-f],  [--forward]   :  Forward Modeling."
        write( *, * ) "        [-j],  [--jmult]     :  JMult."
        write( *, * ) "        [-jt], [--jmult_t]   :  Transposed JMult."
        write( *, * ) "        [-i],  [--inversion] :  Inversion."
        write( *, * )
        write( *, * ) "    Other arguments:"
        write( *, * ) "        [-d],  [--data]      :  Flag to precede data file path."
        write( *, * ) "        [-m],  [--model]     :  Flag to precede model file path."
        write( *, * ) "        [-pm], [--pmodel]    :  Flag to precede perturbation model file path."
        write( *, * ) "        [-c],  [--control]   :  Flag to precede user control file path."
        write( *, * ) "        [-o],  [--outdir]    :  Flag to precede output directory path."
        write( *, * ) "        [-dm], [--dmodel]    :  Flag to precede output dsigma model file path."
        write( *, * ) "        [-pd], [--predicted] :  Flag to precede output predicted data file path."
        write( *, * ) "        [-jm], [--jmhat]     :  Flag to precede output JmHat data file path."
        write( *, * ) "        [-es], [--esolution] :  Flag to precede binary output e-solution file path."
        write( *, * ) "        [-v],  [--version]   :  Print version."
        write( *, * ) "        [-vb], [--verbose]   :  Print runtime information."
        write( *, * ) "        [-h],  [--help]      :  Print usage information."
        !
        write( *, * ) ""
        write( *, * ) "Version 1.0.0"
        !
    end subroutine printHelp
    !
end program ModEM
!