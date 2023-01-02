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
        !> Run ForwardModelling to calculate predicted data
        call runForwardModeling( sigma0 )
        !
        !> Write Predicted Data, with its proper Rx headers, into to the file <predicted_data_file_name>
        call writeDataGroupArray( all_predicted_data, predicted_data_file_name )
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
        type( DataGroupTx_t ), allocatable, dimension(:) :: JmHat
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
        !> Run ForwardModelling to calculate predicted data
        call runForwardModeling( sigma0 )
        !
        !> Initialize JmHat data array in the same format as the predicted data array
        JmHat = all_predicted_data
        !
        !> Calculate JmHat Data array from JMult routine
        call JMult( sigma0, pmodel, JmHat )
        !
        !> Write JmHat to the file <JmHat_data_file_name>
        call writeDataGroupArray( JmHat, JmHat_data_file_name )
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
        type( DataGroupTx_t ), allocatable, dimension(:) :: all_measured_data
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
            all_measured_data = all_predicted_data
            !
        else
            stop "Error: jobJMult_T > Missing Data file!"
        endif
        !
        !> Instantiate the ForwardSolver - Specific type can be chosen via control file
        call createForwardSolver()
        !
        !> Run ForwardModelling to calculate predicted data
        call runForwardModeling( sigma0 )
        !
        !> Calculate DSigma model from JMult_T routine
        call JMult_T( sigma0, all_measured_data, dsigma )
        !
        !> Write dsigma to <dsigma_file_name> file path
        call dsigma%write()
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
    subroutine jobInversionDCG()
        implicit none
        !
        class( ModelParameter_t ), allocatable :: sigma
        !
        real( kind=prec ) :: lambda
        !
        ! Verbose
        write( *, * ) "     - Start jobInversionDCG"
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
            pmodel = sigma0
            !
            call pmodel%zeros()
            !
        else
            stop "Error: jobInversionDCG > Missing Model file!"
        endif
        !
        !> Read Perturbation Model File: instantiate pmodel (NOT USING RIGHT NOW ????)
        if( has_pmodel_file ) then 
            !
            call handlePModelFile()
            !
            call pmodel%setMetric( model_operator%metric )
            !
            pmodel = model_cov%multBy_Cm( pmodel )
            !
            call sigma%add( pmodel )
            !
        endif
        !
        !> Read Data File: instantiates and builds the Data relation between Txs and Rxs
        if( has_data_file ) then 
            !
            call handleDataFile()
            !
            !> Initialize array with measure data
            all_measured_data = all_predicted_data
            !
        else
            stop "Error: jobInversionDCG > Missing Data file!"
        endif
        !
        !> Instantiate ForwardSolver - Specific type via control file
        call createForwardSolver()
        !
        lambda = 10.
        !
        call DCGsolver( all_measured_data, sigma, pmodel, lambda )
		!call DCGsolverLanczos( all_measured_data, sigma, pmodel, lambda )
        !
        ! Verbose
        write( *, * ) "     - Finish jobInversionDCG"
        !
    end subroutine jobInversionDCG
    !
    !> Routine to run a full Inversion Job - Minimize Residual
    !> Where:
    !>    SIGMA (M) = Production model (for predicted data and final inversion model)
    !>    PMODEL = Perturbation model  (if exist -dm readed input model)
    !>    SIGMA0 = Readed input model  (-m)
    !>    DSIGMA = From JMult_T        (????)
    subroutine jobInversionDCG_FULL()
        implicit none
        !
        class( ModelParameter_t ), allocatable :: sigma, dsigma
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: dx, JmHat, resJmHat
        !
        type( IterControl_t ) :: CGiter
        !
        integer :: iter_dcg
        !
        real( kind=prec ) :: alpha, beta, lambda, residual_rmsd
        !
        ! Verbose
        write( *, * ) "     - Start jobInversionDCG"
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
            pmodel = sigma0
            !
            call pmodel%zeros()
            !
        else
            stop "Error: jobInversionDCG > Missing Model file!"
        endif
        !
        !> Read Perturbation Model File: instantiate pmodel (NOT USING RIGHT NOW ????)
        if( has_pmodel_file ) then 
            !
            call handlePModelFile()
            !
            call pmodel%setMetric( model_operator%metric )
            !
            dsigma = model_cov%multBy_Cm( pmodel )
            !
            call sigma%linComb( ONE, ONE, dsigma )
            !
        endif
        !
        !> Read Data File: instantiates and builds the Data relation between Txs and Rxs
        if( has_data_file ) then 
            !
            call handleDataFile()
            !
            !> Initialize array with measure data
            all_measured_data = all_predicted_data
            !
        else
            stop "Error: jobInversionDCG > Missing Data file!"
        endif
        !
        !> Instantiate ForwardSolver - Specific type via control file
        call createForwardSolver()
        !
        !> Run ForwardModelling for calculate Predicted Data
        call runForwardModeling( sigma )
        !
        residual_rmsd = getResidualRMS()
        !
        !> Initialize the JmHat data array in the same format as the predicted data array
        JmHat = all_predicted_data
        !
        dx = all_measured_data
        !
        iter_dcg = 1
        lambda = 10.
        !
        call setIterControl( CGiter )
        !
        ! Verbose
        write( *, * ) "### START DCG LOOP, RMSD: ", residual_rmsd
        !
        !> #### LOOP DCG: Iterations for Inversion Convergence
        do while( residual_rmsd .GT. 1.05 .AND. iter_dcg .LT. 4 )
            !
            !> Calculate JmHat
            call JMult( sigma0, pmodel, JmHat )
            !
            !> resJmHat (b) = res + JmHat
            resJmHat = all_residual_data
            !
            !call linComb(ONE,res,ONE,JmHat,b)
            call linCombDataGroupTxArray( ONE, all_residual_data, ONE, JmHat, resJmHat )
            !
            !> NEED THIS?
            call normalizeDataGroupTxArray( resJmHat, 1 )
            !
            call CG_DS_standard( resJmHat, dx, sigma0, all_measured_data, lambda, CGiter )
            !
            call normalizeWithDataGroupTxArray( 1, all_measured_data, dx )
            !
            !> Calculate dsigma from JMult_T
            call JMult_T( sigma0, dx, dsigma )
            !
            !> Save dsigma in pmodel, to use in the next DCG step
            !pmodel = dsigma
            !
            !> Smooth dsigma and add it to sigma
            pmodel = model_cov%multBy_Cm( dsigma )
            !
            sigma = sigma0
            !
            call sigma%linComb( ONE, ONE, pmodel )
            !
            !> Run ForwardModelling for calculate Predicted Data with new sigma
            call runForwardModeling( sigma )
            !
            residual_rmsd = getResidualRMS()
            !
            ! Verbose
            write( *, * ) "          - #ITER_DCG: ", iter_dcg, residual_rmsd
            !
            iter_dcg = iter_dcg + 1
            !
            !> WRITE DSIGMA FOREACH STEP ????
        enddo
        !
        !> Write dsigma to <dsigma_file_name> file path
        call sigma%write()
        !
        !> Flush due local variables
        deallocate( sigma, JmHat )
        !
        ! Verbose
        write( *, * ) "     - Finish jobInversionDCG"
        !
    end subroutine jobInversionDCG_FULL
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
        if( .NOT. allocated( all_predicted_data ) ) then
            !
            write( *, * ) "     - Create predicted data array"
            !
            !> Create an array of DataGroupTx to store the predicted data in the same format as the measured data
            !> Enabling the use of grouped predicted data in future jobs (all_predicted_data)
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
                call updateDataGroupTxArray( all_predicted_data, tx_data )
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
                call jobInversionDCG()
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
                    write( *, * ) "               ", i, air_layer%dz(i)
                enddo
                !
                write( *, * ) "          Air layers from the method: ", air_layer%method, "."
                !
                write( *, * ) "          Top of the air layers is at ", sum( air_layer%Dz ) / 1000, " km."
                !
                write( *, "(A, F10.2, A2, F10.2, A2, F10.2, A4, F10.2)" ) "           o(x_data,y,z) * rotDeg: (", main_grid%ox, ", ", main_grid%oy, ", ", main_grid%oz, ") * ", main_grid%rotDeg
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
             write( *, * ) "Number of Tx mismatched from Header :[", size( transmitters ), " and ", data_file_standard%n_tx, "]"
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
        JmHat_data_file_name     = "JmHat.dat"
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
        !> Flush memory used for global data arrays
        deallocate( measured_data, all_predicted_data )
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
        character(*), intent( in ) :: file_name
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
                            stop "Error: writeDataGroupArray: Unclassified data_group!"
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
    subroutine printUsage()
        implicit none
        !
        write( *, * ) "ModEM Minimal Usage:"
        write( *, * ) ""
        write( *, * ) "    Forward Modeling:"
        write( *, * ) "        <ModEM> -f -m <rFile_Model> -d <rFile_Data>"
        write( *, * ) "    Outputs:"
        write( *, * ) "        - predicted_data.dat or the path specified by [-pd]"
        write( *, * ) "        - esolution.bin or the path specified by      [-es]"
        write( *, * ) ""
        write( *, * ) "    JMult:"
        write( *, * ) "        <ModEM> -j -m <rFile_Model> -pm <rFile_pModel> -d <rFile_Data>"
        write( *, * ) "    Outputs:"
        write( *, * ) "        - JmHat.dat or the path specified by          [-jm]"
        write( *, * ) "        - predicted_data.dat or the path specified by [-pd]"
        write( *, * ) "        - esolution.bin or the path specified by      [-es]"
        write( *, * ) ""
        write( *, * ) "    JMult_T:"
        write( *, * ) "        <ModEM> -jt -m <rFile_Model> -d <rFile_Data>"
        write( *, * ) "    Output:"
        write( *, * ) "        - dsigma.mod or the path specified by         [-dm]"
        write( *, * ) "        - predicted_data.dat or the path specified by [-pd]"
        write( *, * ) "        - esolution.bin or the path specified by      [-es]"
        write( *, * ) ""
        write( *, * ) "    Inversion:"
        write( *, * ) "        <ModEM> -i -m <rFile_Model> -d <rFile_Data>"
        write( *, * ) "    Output:"
        write( *, * ) "        - dsigma.mod or the path specified by         [-dm]"
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
        write( *, * ) "        [-d],  [--data]      :  Flag for the input data file path."
        write( *, * ) "        [-m],  [--model]     :  Flag for the input model file path."
        write( *, * ) "        [-pm], [--pmodel]    :  Flag for the input perturbation model file path."
        write( *, * ) "        [-c],  [--control]   :  Flag for the user control file path."
        write( *, * ) "        [-dm], [--dmodel]    :  Flag for the output dsigma model file path."
        write( *, * ) "        [-pd], [--predicted] :  Flag for the output predicted data file path."
        write( *, * ) "        [-jm], [--jmhat]     :  Flag for the output JmHat data file path."
        write( *, * ) "        [-es], [--esolution] :  Flag for the binary output e-solution file path."
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