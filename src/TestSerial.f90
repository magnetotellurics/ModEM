!
!> Prototype serial version for testing the ModEM-OO program
!
!> with three main jobs implemented: ForwardModelling, Adjoint and JMult_T
!
!> Obs.: To merge with the MPI version in the future
!
program TestSerial
    !
    use ModEMControlFile
    !
    use ModelReader_Weerachai
    !
    use DataFileStandard
    !
    use InversionDCG
    use InversionNLCG
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
    write( *, "(A18, F10.2, A3)" ) "Finish ModEM-OO (", t_finish - t_start, "s)"
    write( *, * )
    !
contains
    !
    !> Routine to run a full ForwardModeling job and deliver the result (PredictedData) in a text file
    subroutine jobForwardModeling()
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: d_pred
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
            stop "Error: jobForwardModeling > Missing Model file!"
        endif
        !
        !> Read Data File: instantiate Txs and Rxs and build the Data relation between them
        if( has_data_file ) then
            !
            call handleDataFile()
            !
        else
            stop "Error: jobForwardModeling > Missing Data file!"
        endif
        !
        !> Instantiate the global ForwardSolver - Specific type can be chosen via control file
        call createForwardSolver()
        !
        !> Run ForwardModelling to calculate d_pred
        d_pred = all_measured_data
        !
        call runForwardModeling( sigma0, d_pred )
        !
        if( has_e_solution_file ) call writeAllESolution( e_solution_file_name )
        !
        !> Write Predicted Data, with its proper Rx headers, into to the file <predicted_data_file_name>
        call writeDataGroupTxArray( d_pred, predicted_data_file_name )
        !
        ! Verbose
        write( *, * ) "     - Finish jobForwardModeling"
        !
    end subroutine jobForwardModeling
    !
    !> Routine to run a full JMult job and deliver the result (JmHat Data) in a text file
    subroutine jobJMult()
        implicit none
        !
        !> Using the same local Data Array for FWD and JMult
        type( DataGroupTx_t ), allocatable, dimension(:) :: tx_data_array
        !
        ! Verbose
        write( *, * ) "     - Start jobJMult"
        !
        !> Read Prior Model File: instantiate pmodel
        if( has_pmodel_file ) then
            !
            call handlePModelFile()
            !
            !> Read Model File and instantiate global variables: main_grid, model_operator and Sigma0
            if( has_model_file ) then
                !
                call handleModelFile()
                !
                call pmodel%setMetric( model_operator%metric )
                !
            else
                stop "Error: jobJMult > Missing Model file!"
            endif
            !
        else
            stop "Error: jobJMult > Missing Prior Model file!"
        endif
        !
        !> Read Data File: instantiate Txs and Rxs and build the Data relation between them
        if( has_data_file ) then
            !
            call handleDataFile()
            !
        else
            stop "Error: jobJMult > Missing Data file!"
        endif
        !
        !> Instantiate the ForwardSolver - Specific type can be chosen via control file
        call createForwardSolver()
        !
        call runEMSolve( sigma0 )
        !
        !> Run ForwardModelling to solve Tx ESolutions (Could be better: Just Tx Solve ????)
        tx_data_array = all_measured_data
        !
        !> Calculate Data Array from JMult routine
        call JMult( sigma0, pmodel, tx_data_array )
        !
        !> Write Data Array to the file <jmhat_data_file_name>
        call writeDataGroupTxArray( tx_data_array, jmhat_data_file_name )
        !
        ! Verbose
        write( *, * ) "     - Finish jobJMult"
        !
    end subroutine jobJMult
    !
    !> Routine to run a full JMult_T job and deliver the result (DSigma model) in a text file
    subroutine jobJMult_T()
        implicit none
        !
        class( ModelParameter_t ), allocatable :: dsigma
        !
        ! Verbose
        write( *, * ) "     - Start jobJMult_T"
        !
        !> Read Model File and instantiate global variables: main_grid, model_operator and Sigma0
        if( has_model_file ) then 
            !
            call handleModelFile()
            !
        else
            stop "Error: jobJMult_T > Missing Model file!"
        endif
        !
        !> Read Data File: instantiate Txs and Rxs and build the Data relation between them
        !> Initialize a DataGroupTxArray to hold the measured data
        if( has_data_file ) then 
            !
            call handleDataFile()
            !
        else
            stop "Error: jobJMult_T > Missing Data file!"
        endif
        !
        !> Instantiate the ForwardSolver - Specific type can be chosen via control file
        call createForwardSolver()
        !
        !> Run ForwardModelling to calculate predicted data
        call runEMSolve( sigma0 )
        !
        !> Calculate DSigma model from JMult_T routine
        call JMult_T( sigma0, all_measured_data, dsigma )
        !
        !> Write dsigma to <dsigma_file_name> file path
        call dsigma%write( dsigma_file_name )
        !
        !> Flush local variable
        deallocate( dsigma )
        !
        ! Verbose
        write( *, * ) "     - Finish jobJMult_T"
        !
    end subroutine jobJMult_T
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
        class( ModelParameter_t ), allocatable :: sigma
        !
        ! Verbose
        write( *, * ) "     - Start jobInversion"
        !
        !> Read Model File and instantiate global variables: main_grid, model_operator and Sigma0
        if( has_model_file ) then 
            !
            call handleModelFile()
            !
            !> Initialize sigma as the input model
            allocate( sigma, source = sigma0 )
            !
            !> Instantiate ModelCovariance
            allocate( model_cov, source = ModelCovarianceRec_t( sigma ) )
            !
            !> Initialize pmodel with Zeros
            allocate( pmodel, source = sigma0 )
            !
            call pmodel%zeros()
            !
        else
            stop "Error: jobInversion > Missing Model file!"
        endif
        !
        !> Read Perturbation Model File: instantiate pmodel (NOT USING RIGHT NOW ????)
        if( has_pmodel_file ) then 
            !
            deallocate( pmodel )
            !
            call handlePModelFile()
            !
            call pmodel%setMetric( model_operator%metric )
            !
            call model_cov%multBy_Cm( pmodel )
            !
            call sigma%linComb( ONE, ONE, pmodel )
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
        !> Instantiate ForwardSolver - Specific type via control file
        call createForwardSolver()
        !
        call DCGsolver( all_measured_data, sigma, pmodel )
        !
        !call NLCGsolver( all_measured_data, lambda, sigma, pmodel )
        !
        !> Flush local variable
        deallocate( sigma )
        !
        ! Verbose
        write( *, * ) "     - Finish jobInversion"
        !
    end subroutine jobInversion
    !
    !> No subroutine briefing
    subroutine createForwardSolver()
        implicit none
        !
        integer :: i_tx
        !
        if( allocated( forward_solver ) ) deallocate( forward_solver )
        !
        !> Instantiate the ForwardSolver - Specific type can be chosen via control file
        select case ( forward_solver_type )
            !
            case( FWD_IT_DC )
                allocate( forward_solver, source = ForwardSolverIT_DC_t( model_operator, QMR ) )
                !
            case default
                !
                write( *, * ) "Warning: createForwardSolver > Undefined forward_solver, using IT_DC"
                !
                allocate( forward_solver, source = ForwardSolverIT_DC_t( model_operator, QMR ) )
                !
        end select
        !
        !> Make all transmitters point to this ForwardSolver
        do i_tx = 1, size( transmitters )
            !
            transmitters( i_tx )%Tx%forward_solver => forward_solver
            !
        enddo
        !
    end subroutine createForwardSolver
    !
    !> Allocate and Initialize the predicted data array, 
    !> according to the arrangement of the Transmitter-Receiver pairs of the input
    subroutine initPredictedDataArray()
        implicit none
        !
        !> Auxiliary variable to group data under a single transmitter index
        type( DataGroupTx_t ) :: tx_data
        !
        !> Local indexes
        integer :: i_tx, i_data
        !
        !> If is the first time, create the predicted data array, 
        !> according to the arrangement of the Transmitter-Receiver pairs of the input
        if( .NOT. allocated( all_measured_data ) ) then
            !
            write( *, * ) "     - Create Predicted Data array"
            !
            !> Create an array of DataGroupTx to store the predicted data in the same format as the measured data
            !> Enabling the use of grouped predicted data in future jobs (all_measured_data)
            do i_tx = 1, size( transmitters )
                !
                tx_data = DataGroupTx_t( i_tx )
                !
                do i_data = 1, size( measured_data )
                    !
                    if( measured_data( i_data )%i_tx == i_tx ) then
                        !
                        call tx_data%put( measured_data( i_data ) )
                        !
                    endif
                enddo
                !
                call updateDataGroupTxArray( all_measured_data, tx_data )
                !
            enddo
            !
        else
            stop "Error: initPredictedDataArray > Unnecessary recreation of predicted data array"
        endif
        !
    end subroutine initPredictedDataArray
    !
    !> No subroutine briefing
    subroutine handleJob()
        implicit none
        !
        select case ( modem_job )
            !
            case ( "Forward" )
                !
                call jobForwardModeling()
                !
            case ( "JMult" )
                !
                call jobJMult()
                !
            case ( "JMult_T" )
                !
                call jobJMult_T()
                !
            case ( "Inversion" )
                !
                call jobInversion()
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
        !> Read control file and sets the options in the Constants module
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
        integer :: i
        !
        ! Verbose
        write( *, * ) "     < Model File: [", model_file_name, "]"
        !
        !> Initialize main_grid and sigma0 with ModelReader (Only ModelReader_Weerachai by now????)
        call model_reader%Read( model_file_name, main_grid, sigma0 ) 
        !
        !> Instantiate the ModelOperator object according to the main_grid type
        select type( main_grid )
            !
            class is( Grid3D_SG_t )
                !
                call main_grid%SetupAirLayers( air_layer, model_method, model_n_air_layer, model_max_height )
                !
                call main_grid%UpdateAirLayers( air_layer%nz, air_layer%dz )
                !
                ! Verbose
                write( *, * ) "          Air Layers [i, dz(i)]:"
                !
                do i = air_layer%nz, 1, -(1)
                    write( *, "( a15, i5, f12.3 )" ) "     ", i, air_layer%dz(i)
                enddo
                !
                write( *, * ) "          Air layers from the method [", trim( air_layer%method ), "]"
                !
                write( *, "( a39, f12.3, a4 )" ) "          Top of the air layers is at ", sum( air_layer%Dz ) / 1000, " km."
                !
                write( *, "( a31, f8.3, a2, f8.3, a2, f8.3, a4, f8.3 )" ) "          o(x,y,z) * rotDeg: (", main_grid%ox, ", ", main_grid%oy, ", ", main_grid%oz, ") * ", main_grid%rotDeg
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
        deallocate( prior_grid )
        !
    end subroutine handlePModelFile
    !
    !> No subroutine briefing
    subroutine handleDataFile()
        implicit none
        !
        integer :: i_rx, n_rx
        !
        !> Local object to dealt data, self-destructs at the end of the subroutine
        type( DataFileStandard_t ) :: data_file_standard
        !
        write( *, * ) "     < Data File: [", data_file_name, "]"
        !
        data_file_standard = DataFileStandard_t( ioStartup, data_file_name )
        !
        n_rx = size( receivers )
        !
        ! Verbose
        if( n_rx == data_file_standard%n_rx ) then
            !
            write( *, * ) "          Checked ", n_rx, " Receivers."
            !
        else
            !
            write( *, * ) "Error: Number of Rx mismatched from Header :[", n_rx, " and ", data_file_standard%n_rx, "]"
            stop
            !
        endif
        !
        if( size( transmitters ) == data_file_standard%n_tx ) then
             !
             call printTransmitterArray()
             !
        else
             !
             write( *, * ) "Error: Number of Tx mismatched from Header :[", size( transmitters ), " and ", data_file_standard%n_tx, "]"
             stop
             !
        endif
        !
        write( *, * ) "          Checked ", size( measured_data ), " DataGroups."
        !
        !> Initialize the predicted data array
        call initPredictedDataArray()
        !
        write( *, * ) "     - Creating Rx evaluation vectors"
        !
        !> Calculate and store evaluation vectors on all Receivers
        do i_rx = 1, n_rx
            !
            call receivers( i_rx )%Rx%evaluationFunction( model_operator )
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
        predicted_data_file_name = "predicted_data.dat"
        jmhat_data_file_name = "jmhat.dat"
        e_solution_file_name = "esolution.bin"
        dsigma_file_name = "dsigma.rho"
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
    subroutine garbageCollector()
        implicit none
        !
        !> Flush memory used for global data arrays
        deallocate( measured_data )
        !
        if( allocated( all_measured_data ) ) call deallocateDataGroupTxArray( all_measured_data )
        !
        !> Deallocate global array of Receivers
        call deallocateReceiverArray()
        !
        !> Deallocate global array of Transmitters
        call deallocateTransmitterArray()
        !
        !> Deallocate global components
        deallocate( forward_solver, model_operator, sigma0, main_grid )
        !
        !> Deallocate global pmodel, if its the case
        if( allocated( pmodel ) ) deallocate( pmodel )
        !
        !> Deallocate global model_cov, if its the case
        if( allocated( model_cov ) ) deallocate( model_cov )
        !
        !> Flush memory used by main program control variables and flags
        if( allocated( forward_solver_type ) ) deallocate( forward_solver_type )
        if( allocated( source_type ) ) deallocate( source_type )
        if( allocated( model_method ) ) deallocate( model_method )
        if( allocated( get_1D_from ) ) deallocate( get_1D_from )
        if( allocated( predicted_data_file_name ) ) deallocate( predicted_data_file_name )
        if( allocated( jmhat_data_file_name ) ) deallocate( jmhat_data_file_name )
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
    subroutine printUsage()
        implicit none
        !
        write( *, * ) "ModEM Minimal Usage:"
        write( *, * ) ""
        write( *, * ) "    Forward Modeling:"
        write( *, * ) "        <ModEM> -f -m <rFile_Model> -d <rFile_Data>"
        write( *, * ) "    Output:"
        write( *, * ) "        - predicted_data.dat or the path specified by         [-pd]"
        write( *, * ) ""
        write( *, * ) "    JMult:"
        write( *, * ) "        <ModEM> -j -m <rFile_Model> -pm <rFile_pModel> -d <rFile_Data>"
        write( *, * ) "    Outputs:"
        write( *, * ) "        - jmhat.dat or the path specified by                  [-jm]"
        write( *, * ) ""
        write( *, * ) "    JMult_T:"
        write( *, * ) "        <ModEM> -jt -m <rFile_Model> -d <rFile_Data>"
        write( *, * ) "    Output:"
        write( *, * ) "        - dsigma.rho or the path specified by                 [-dm]"
        write( *, * ) ""
        write( *, * ) "    Inversion:"
        write( *, * ) "        <ModEM> -i -m <rFile_Model> -d <rFile_Data>"
        write( *, * ) "    Output:"
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
end program TestSerial
!