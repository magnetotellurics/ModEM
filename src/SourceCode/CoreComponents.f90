!
!> Module to group the main Global components (variables and routines) of the program.
!
module CoreComponents
    !
    use ForwardControlFile
    use InversionControlFile
    !
    use ModelOperator_MF_SG
    use ModelOperator_SP_V1
    use ModelOperator_SP_V2
    !
    use ModelParameterCell_SG
    use ModelParameterCell_MR
    !
    use ModelCovariance
    !
    use ForwardSolver_IT_DC
    !
    use SourceMT_1D
    use SourceMT_2D
    use SourceCSEM_EM1D
    use SourceCSEM_Dipole1D
    use SourceAdjoint
    !
    use ReceiverFullImpedance
    use ReceiverFullVerticalMagnetic
    use ReceiverOffDiagonalImpedance
    use ReceiverSingleField
    !
    use DataGroupTxArray
    !
    use ModelReader_Weerachai
    !
    use DataFileStandard
    !
    !> Program Control Variables
    character(8) :: str_date
    character(6) :: str_time
    character(50) :: outdir_name
    !
    character(:), allocatable :: modem_job
    !
    character(:), allocatable :: fwd_control_file_name, inv_control_file_name
    character(:), allocatable :: model_file_name, pmodel_file_name
    character(:), allocatable :: data_file_name
    character(:), allocatable :: dsigma_file_name
    character(:), allocatable :: cov_file_name
    character(:), allocatable :: e_solution_file_name
    !
    !> Program Control Flags
    logical :: has_outdir_name
    logical :: has_fwd_control_file, has_inv_control_file
    logical :: has_model_file, has_pmodel_file
    logical :: has_cov_file
    logical :: has_data_file
    logical :: has_e_solution_file
    logical :: verbosis
    !
    !> Public Module Routines
    public :: handleModelFile, handlePModelFile, handleDataFile
    public :: handleForwardControlFile, handleInversionControlFile
    public :: handleArguments
    public :: setupDefaultParameters
    public :: createOutputDirectory
    public :: getLiteralTime
    public :: garbageCollector
    public :: printUsage, printHelp
    public :: printForwardControlFileTemplate
    public :: printInversionControlFileTemplate
    !
contains
    !
    !> Read Model File and instantiate global variables: main_grid, model_operator and sigma0
    !
    subroutine handleModelFile( sigma0 )
        implicit none
        !
        class( ModelParameter_t ), allocatable, intent( out ) :: sigma0
        !
        type( ModelReader_Weerachai_t ) :: model_reader
        type( TAirLayers ) :: air_layer
        !
        integer :: i
        !
        ! Verbose
        write( *, * ) "     < Model File: [", model_file_name, "]"
        !
        !> Initialize main_grid and sigma0 with ModelReader
        !> Only ModelReader_Weerachai by now ????
        call model_reader%read( model_file_name, main_grid, sigma0 ) 
        !
        call main_grid%setupAirLayers( air_layer, model_method, model_n_air_layer, model_max_height )
        !
        call main_grid%updateAirLayers( air_layer%nz, air_layer%dz )
        !
        ! Verbose
        write( *, "( a39 )" ) "Model Air Layers [i, dz(i)]:"
        !
        do i = air_layer%nz, 1, -(1)
            write( *, "( i20, f20.3 )" ) i, air_layer%dz(i)
        enddo
        !
        write( *, * ) "          Air layers from the method [", trim( air_layer%method ), "]"
        !
        write( *, "( a39, f12.3, a4 )" ) "Top of the air layers is at ", ( sum( air_layer%Dz ) / 1000. ), " km."
        !
        write( *, "( a33, f12.3 )" ) "Air layers Max Height ", air_layer%maxHeight
        !
        write( *, "( a23, i6, a2, i6, a2, i6, a1 )" ) "dim(x,y,z):(", main_grid%nx, ", ", main_grid%ny, ", ", main_grid%nz, ")"
        !
        write( *, "( a30, f16.3, a2, f16.3, a2, f16.3, a4, f16.3 )" ) "o(x,y,z) * rotDeg:(", main_grid%ox, ", ", main_grid%oy, ", ", main_grid%oz, ") * ", main_grid%rotDeg
        !
        !> Instantiate model_operator
        !> Specific type can be chosen via fwd control file
        select case( model_operator_type )
            !
            case( MODELOP_MF )
                !
                allocate( model_operator, source = ModelOperator_MF_SG_t( main_grid ) )
                !
            case( MODELOP_SP )
                !
                allocate( model_operator, source = ModelOperator_SP_V1_t( main_grid ) )
                !
            case( MODELOP_SP2 )
                !
                allocate( model_operator, source = ModelOperator_SP_V2_t( main_grid ) )
                !
            case( "" )
                !
                call warning( "handleModelFile > model_operator_type not provided, using ModelOperator_MF_SG_t." )
                !
                allocate( model_operator, source = ModelOperator_MF_SG_t( main_grid ) )
                !
            case default
                !
                call errStop( "handleModelFile > Wrong Model Operator type: ["//model_operator_type//"]" )
                !
        end select
        !
        call model_operator%setEquations
        !
        call sigma0%setMetric( model_operator%metric )
        !
    end subroutine handleModelFile
    !
    !> Read Perturbation Model File: instantiate pmodel with ModelReader
    !> Only ModelReader_Weerachai by now ????
    !
    subroutine handlePModelFile( pmodel )
        implicit none
        !
        class( ModelParameter_t ), allocatable, intent( out ) :: pmodel
        !
        type( ModelReader_Weerachai_t ) :: model_reader
        class( Grid_t ), allocatable :: prior_grid
        !
        ! Verbose
        write( *, * ) "     < PModel File: [", pmodel_file_name, "]"
        !
        !> Read prior_grid and pmodel with ModelReader_Weerachai
        call model_reader%Read( pmodel_file_name, prior_grid, pmodel ) 
        !
        deallocate( prior_grid )
        !
    end subroutine handlePModelFile
    !
    !> Read Data File: instantiate Txs and Rxs and build the Data relation between them
    !
    subroutine handleDataFile()
        implicit none
        !
        !> Local object to dealt data, self-destructs at the end of the subroutine
        type( DataFileStandard_t ) :: data_file_standard
        !
        character( len=100 ) :: str_msg
        integer :: i_tx, i_rx, n_rx, n_tx
        !
        !>
        write( *, * ) "     < Data File: [", data_file_name, "]"
        !
        data_file_standard = DataFileStandard_t( ioStartup, data_file_name )
        !
        n_rx = size( receivers )
        !
        n_tx = size( transmitters )
        !
        ! Verbose
        if( n_rx == data_file_standard%n_rx ) then
            !
            write( *, * ) "          Checked ", n_rx, " Receivers."
            !
        else
            !
            write( str_msg, "( A55, I2, A5, I2, A1 )") "handleDataFile > Number of Rx mismatched from Header :[", n_rx, " and ", data_file_standard%n_rx, "]"
            call errStop( str_msg )
            !
        endif
            !
        if( n_tx == data_file_standard%n_tx ) then
            !
            write( *, * ) "          Checked ", n_tx, " Transmitters."
            !
            do i_tx = 1, n_tx
                call transmitters( i_tx )%Tx%print
            enddo
            !
        else
            !
            write( str_msg, "( A54, I2, A5, I2, A1 )") "handleDataFile > Number of Tx mismatched from Header :[", n_tx, " and ", data_file_standard%n_tx, "]"
            call errStop( str_msg )
            !
        endif
        !
        write( *, * ) "          Checked ", countData( all_measured_data ), " DataGroups."
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
    !> Read Forward Modeling Control File: 
    !> set the parameters as described in the file
    !
    subroutine handleForwardControlFile()
        implicit none
        !
        write( *, * ) "     < FWD control File: [", fwd_control_file_name, "]"
        !
        allocate( fwd_control_file, source = ForwardControlFile_t( ioStartup, fwd_control_file_name ) )
        !
    end subroutine handleForwardControlFile
    !
    !> Read Inversion Control File:
    !> set the parameters as described in the file
    !
    subroutine handleInversionControlFile()
        implicit none
        !
        write( *, * ) "     < INV control File: [", inv_control_file_name, "]"
        !
        allocate( inv_control_file, source = InversionControlFile_t( ioStartup, inv_control_file_name ) )
        !
    end subroutine handleInversionControlFile
    !
    !> Routine to run a full ForwardModeling job 
    !> and deliver the result(PredictedData) in a text file <all_predicted_data.dat>
    !
    subroutine jobSplitModel()
        implicit none
        !
        class( ModelParameter_t ), allocatable :: model, aux_model
        type( rScalar3D_SG_t ) :: cell_cond
        !
        ! Verbose
        write( *, * ) "     - Start jobSplitModel"
        !
        if( has_model_file ) then
            !
            call handleModelFile( model )
        !
        else
            call errStop( "jobSplitModel > Missing Model file!" )
        endif
        !
        !> Just VTI implemented yet
        if( model%anisotropic_level /= 2 ) then
            !
            call errStop( "jobSplitModel unsupported for this model" )
            !
        else
            !
            !> Create new isotropic model with target horizontal cond
            cell_cond = model%getCond(1)
            !
            select type( grid => model%metric%grid )
                !
                class is( Grid3D_SG_t )
                    allocate( aux_model, source = ModelParameterCell_SG_t( cell_cond, 1, model%param_type ) )
                class is( Grid3D_MR_t )
                    allocate( aux_model, source = ModelParameterCell_MR_t( cell_cond, 1, model%param_type, grid%cs ) )
                class default
                    call errStop( "jobSplitModel > Unknow grid for ModelParameter" )
                !
            end select
            !
            call aux_model%setMetric( model_operator%metric )
            !
            call aux_model%write( model_file_name//"_h" )
            !
            write( *, * ) "               < Created ["//model_file_name//"_h] file."
            !
            !> Set new model horizontal as the target vertical cond
            cell_cond = model%getCond(2)
            !
            call aux_model%setCond( cell_cond, 1 )
            !
            call aux_model%write( model_file_name//"_v" )
            !
            write( *, * ) "               < Created ["//model_file_name//"_v] file."
            !
            deallocate( model, aux_model )
            !
            ! Verbose
            write( *, * ) "     - Finish jobSplitModel"
            !
        endif
        !
    end subroutine jobSplitModel
    !
    !> Handle all the arguments passed in the execution command line
    !
    subroutine handleArguments()
        implicit none
        !
        character( len=200 ) :: argument
        integer :: argument_index
        !
        if( command_argument_count() == 0 ) then
            !
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
                select case( argument )
                    !
                    case( "-cf", "--ctrl_fwd" )
                        !
                        call get_command_argument( argument_index + 1, argument )
                        fwd_control_file_name = trim( argument )
                        !
                        if( len( fwd_control_file_name ) > 0 ) has_fwd_control_file = .TRUE.
                        !
                        argument_index = argument_index + 2
                        !
                    case( "-ci", "--ctrl_inv" )
                        !
                        call get_command_argument( argument_index + 1, argument )
                        inv_control_file_name = trim( argument )
                        !
                        if( len( inv_control_file_name ) > 0 ) has_inv_control_file = .TRUE.
                        !
                        argument_index = argument_index + 2
                        !
                    case( "-d", "--data" )
                        !
                        call get_command_argument( argument_index + 1, argument )
                        data_file_name = trim( argument )
                        !
                        if( len( data_file_name ) > 0 ) has_data_file = .TRUE.
                        !
                        argument_index = argument_index + 2
                        !
                    case( "-f", "--forward" )
                        !
                        modem_job = "JobForwardModeling"
                        !
                        argument_index = argument_index + 1
                        !
                    case( "-j", "--jmult" )
                        !
                        modem_job = "JobJMult"
                        !
                        argument_index = argument_index + 1
                        !
                    case( "-jt", "--jmult_t" )
                        !
                        modem_job = "JobJMult_T"
                        !
                        argument_index = argument_index + 1
                        !
                    case( "-i", "--inversion" )
                        !
                        modem_job = "JobInversion"
                        !
                        argument_index = argument_index + 1
                        !
                    case( "-s", "--split" )
                        !
                        modem_job = "JobSplitModel"
                        !
                        argument_index = argument_index + 1
                        !
                    case( "-m", "--model" )
                        !
                        call get_command_argument( argument_index + 1, argument )
                        model_file_name = trim( argument )
                        !
                        if( len( model_file_name ) > 0 ) has_model_file = .TRUE.
                        !
                        argument_index = argument_index + 2
                        !
                    case( "-c", "--cov" )
                        !
                        call get_command_argument( argument_index + 1, argument )
                        cov_file_name = trim( argument )
                        !
                        if( len( cov_file_name ) > 0 ) has_cov_file = .TRUE.
                        !
                        argument_index = argument_index + 2
                        !
                    case( "-o", "--outdir" )
                        !
                        call get_command_argument( argument_index + 1, argument )
                        outdir_name = trim( argument )
                        !
                        if( len( outdir_name ) > 0 ) has_outdir_name = .TRUE.
                        !
                        argument_index = argument_index + 2
                        !
                    case( "-pm", "--pmodel" )
                        !
                        call get_command_argument( argument_index + 1, argument )
                        pmodel_file_name = trim( argument )
                        !
                        if( len( pmodel_file_name ) > 0 ) has_pmodel_file = .TRUE.
                        !
                        argument_index = argument_index + 2
                        !
                    case( "-dm", "--dsigma" )
                        !
                        call get_command_argument( argument_index + 1, argument )
                        dsigma_file_name = trim( argument )
                        !
                        argument_index = argument_index + 2
                        !
                    case( "-pd", "--predicted" )
                        !
                        call get_command_argument( argument_index + 1, argument )
                        predicted_data_file_name = trim( argument )
                        !
                        argument_index = argument_index + 2
                        !
                    case( "-jm", "--jmhat" )
                        !
                        call get_command_argument( argument_index + 1, argument )
                        jmhat_data_file_name = trim( argument )
                        !
                        argument_index = argument_index + 2
                        !
                    case( "-es", "--esolution" )
                        !
                        call get_command_argument( argument_index + 1, argument )
                        e_solution_file_name = trim( argument )
                        !
                        if( len( e_solution_file_name ) > 0 ) has_e_solution_file = .TRUE.
                        !
                        argument_index = argument_index + 2
                        !
                    case( "-v", "--version" )
                        !
                        write( *, * ) "    + ModEM-OO version "//VERSION
                        stop
                        !
                    case( "-h", "--help" )
                        !
                        call printHelp()
                        stop
                        !
                    case( "-tmp", "--template" )
                        !
                        call printForwardControlFileTemplate
                        call printInversionControlFileTemplate
                        stop
                        !
                    case( "-vb", "--verbose" )
                        !
                        call errStop( "handleArguments > Verbose level not implemented yet!" )
                        !
                    case default
                        !
                        write( *, * ) "     "//achar(27)//"[31m# Error:"//achar(27)//"[0m Unknown Argument: [", trim( argument ), "]"
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
    !> Set default values for all the control variables
    !
    subroutine setupDefaultParameters()
        implicit none
        !
        ! I|O
        predicted_data_file_name = "all_predicted_data.dat"
        jmhat_data_file_name = "jmhat.dat"
        e_solution_file_name = "esolution.bin"
        dsigma_file_name = "dsigma.rho"
        !
        ! Control flags
        has_outdir_name = .FALSE.
        !
        has_fwd_control_file = .FALSE.
        has_inv_control_file = .FALSE.
        has_model_file = .FALSE.
        has_pmodel_file = .FALSE.
        has_data_file = .FALSE.
        has_cov_file = .FALSE.
        has_e_solution_file = .FALSE.
        verbosis = .FALSE.
        !
        ! Solver parameters
        max_solver_iters = 80
        max_solver_calls = 20
        max_divcor_iters = 100
        tolerance_divcor = 1E-5
        tolerance_solver = 1E-7
        !
        ! Inversion Parameters
        inversion_type = NLCG
        !
        ! Forward Modeling Parameters
        solver_type = ""
        forward_solver_type = FWD_IT_DC
        model_operator_type = MODELOP_MF
        !
        ! Source parameters
        source_type_mt = SRC_MT_1D
        get_1d_from = "Geometric_mean"
        !
        ! Model parameters
        model_method = MM_METHOD_FIXED_H
        model_n_air_layer = 10
        model_max_height = 200.0
        !
    end subroutine setupDefaultParameters
    !
    !> Create a directory to store output files
    !> Named by date and runtime or as specified by the input argument [-o]
    !
    subroutine createOutputDirectory()
        implicit none
        !
        if( .NOT. has_outdir_name ) then
            !
            select case( inversion_type )
                !
                case( DCG )
                    !
                    write( outdir_name, "(a11, a8, a1, a6)" ) "DCG_Output_", str_date, "_", str_time
                    !
                case( NLCG )
                    !
                    write( outdir_name, "(a12, a8, a1, a6)" ) "NLCG_Output_", str_date, "_", str_time
                    !
                case default
                    !
                    call errStop( "jobInversion > Undefined inversion_type" )
                    !
            end select
            !
            !write( *, * ) "outdir_name: [", trim( outdir_name ), "]"
            !
        endif
        !
        !> Create a directory with the specified path <outdir_name>
        !> With new f2008 system command
        call EXECUTE_COMMAND_LINE( "mkdir -p "//outdir_name )
        !
    end subroutine createOutputDirectory
    !
    !> Translate the input in seconds to a formated string 0d0h0m0s
    !
    function getLiteralTime( seconds ) result( str_time )
        implicit none
        !
        integer, intent( inout ) :: seconds
        !
        character(50) :: str_time
        !
        integer :: days, hours, minutes
        !
        days = seconds / 86400
        !
        seconds = mod( seconds, 86400 )
        !
        hours = seconds / 3600
        !
        seconds = mod( seconds, 3600 )
        !
        minutes = seconds / 60
        !
        seconds = mod( seconds, 60 )
        !
        if( days .EQ. 0 .AND. hours .EQ. 0 .AND. minutes .EQ. 0 ) then
            write( str_time, "(i2, a1)" ) seconds, "s"
        elseif( days .EQ. 0 .AND. hours .EQ. 0 ) then
            write( str_time, "(i2, a2, i2, a1)" ) minutes, "m,", seconds, "s"
        elseif( days .EQ. 0 .AND. minutes .EQ. 0 ) then
            write( str_time, "(i2, a2, i2, a1)" ) hours, "h,", seconds, "s"
        elseif( hours .EQ. 0 .AND. minutes .EQ. 0  ) then
            write( str_time, "(i2, a2, i2, a2, i2, a2, i2, a1)" ) days, "d,", seconds, "s"
        elseif( days .EQ. 0 ) then
            write( str_time, "(i2, a2, i2, a2, i2, a1)" ) hours, "h,", minutes, "m,", seconds, "s"
        elseif( hours .EQ. 0 ) then
            write( str_time, "(i2, a2, i2, a2, i2, a1)" ) days, "d,", minutes, "m,", seconds, "s"
        elseif( minutes .EQ. 0  ) then
            write( str_time, "(i2, a2, i2, a2, i2, a2, i2, a1)" ) days, "d,", hours, "h,", seconds, "s"
        else
            write( str_time, "(i2, a2, i2, a2, i2, a2, i2, a1)" ) days, "d,", hours, "h,", minutes, "m,", seconds, "s"
        endif
        !
    end function getLiteralTime
    !
    !> Deallocate remaining unallocated global memory
    !
    subroutine garbageCollector()
        implicit none
        !
        !> Deallocate global array of measured data
        if( allocated( all_measured_data ) ) call deallocateData( all_measured_data )
        !
        !> Deallocate global array of Receivers
        if( allocated( receivers ) ) call deallocateReceiverArray()
        !
        !> Deallocate global array of Transmitters
        if( allocated( transmitters ) ) call deallocateTransmitterArray()
        !
        !> Deallocate global components
        if( allocated( forward_solver ) ) deallocate( forward_solver )
        if( allocated( model_operator ) ) deallocate( model_operator )
        if( allocated( main_grid ) ) deallocate( main_grid )
        !
        !> Deallocate global model_cov, if its the case
        if( allocated( model_cov ) ) deallocate( model_cov )
        !
        !> Flush memory used by main program control variables and flags
        if( allocated( inversion_type ) ) deallocate( inversion_type )
        if( allocated( joint_type ) ) deallocate( joint_type )
        !
        if( allocated( grid_format ) ) deallocate( grid_format )
        if( allocated( model_operator_type ) ) deallocate( model_operator_type )
        if( allocated( solver_type ) ) deallocate( solver_type )
        if( allocated( forward_solver_type ) ) deallocate( forward_solver_type )
        if( allocated( source_type_mt ) ) deallocate( source_type_mt )
        if( allocated( source_type_csem ) ) deallocate( source_type_csem )
        if( allocated( model_method ) ) deallocate( model_method )
        if( allocated( get_1d_from ) ) deallocate( get_1d_from )
        !
        if( allocated( predicted_data_file_name ) ) deallocate( predicted_data_file_name )
        if( allocated( jmhat_data_file_name ) ) deallocate( jmhat_data_file_name )
        if( allocated( e_solution_file_name ) ) deallocate( e_solution_file_name )
        if( allocated( model_file_name ) ) deallocate( model_file_name )
        if( allocated( pmodel_file_name ) ) deallocate( pmodel_file_name )
        if( allocated( data_file_name ) ) deallocate( data_file_name )
        if( allocated( dsigma_file_name ) ) deallocate( dsigma_file_name )
        if( allocated( cov_file_name ) ) deallocate( cov_file_name )
        !
        if( allocated( fwd_control_file_name ) ) deallocate( fwd_control_file_name )
        if( allocated( fwd_control_file ) ) deallocate( fwd_control_file )
        if( allocated( inv_control_file_name ) ) deallocate( inv_control_file_name )
        if( allocated( inv_control_file ) ) deallocate( inv_control_file )
        !
        if( allocated( units_in_file ) ) deallocate( units_in_file )
        !
        if( allocated( modem_job ) ) deallocate( modem_job )
        !
    end subroutine garbageCollector
    !
    !> Print on the screen basic information on how to use the program.
    !
    subroutine printUsage()
        implicit none
        !
        write( *, * ) "ModEM_"//VERSION//" Minimal Usage:"
        write( *, * ) ""
        write( *, * ) "    Forward Modeling (FWD):"
        write( *, * ) "        <ModEM> -f -m <rFile_Model> -d <rFile_Data>"
        write( *, * ) "        Output:"
        write( *, * ) "        - 'all_predicted_data.dat' or the path specified by      [-pd]"
        write( *, * ) ""
        write( *, * ) "    Jacobian Multiplication (JMult):"
        write( *, * ) "        <ModEM> -j -m <rFile_Model> -pm <rFile_pModel> -d <rFile_Data>"
        write( *, * ) "        Output:"
        write( *, * ) "        - 'jmhat.dat' or the path specified by                   [-jm]"
        write( *, * ) ""
        write( *, * ) "    Transposed J Multiplication (JMult_T):"
        write( *, * ) "        <ModEM> -jt -m <rFile_Model> -d <rFile_Data>"
        write( *, * ) "        Output:"
        write( *, * ) "        - 'dsigma.rho' or the path specified by                  [-dm]"
        write( *, * ) ""
        write( *, * ) "    Inversion:"
        write( *, * ) "        <ModEM> -i -m <rFile_Model> -d <rFile_Data>"
        write( *, * ) "        Output:"
        write( *, * ) "        - directory named 'Output_<date>_<time>' or specified by [-o]"
        write( *, * ) ""
        write( *, * ) "    Other options:"
        write( *, * ) "        <ModEM> -h or <ModEM> --help"
        write( *, * ) ""
        !
    end subroutine printUsage
    !
    !> Print on the screen all the supported input parameters.
    !
    subroutine printHelp()
        implicit none
        !
        write( *, * ) "ModEM_"//VERSION//" Options:"
        write( *, * ) ""
        write( *, * ) "    Flags to define a job:"
        write( *, * ) "        [-f],  [--forward]   :  Forward Modeling."
        write( *, * ) "        [-j],  [--jmult]     :  Jacobian Multiplication."
        write( *, * ) "        [-jt], [--jmult_t]   :  Transposed Jacobian Multiplication."
        write( *, * ) "        [-i],  [--inversion] :  Inversion."
        write( *, * ) "        [-s],  [--split]     :  Split one input model into single axes."
        write( *, * )
        write( *, * ) "    Other arguments:"
        write( *, * ) "        [-d],  [--data]      :  Flag to precede data file path."
        write( *, * ) "        [-m],  [--model]     :  Flag to precede model file path."
        write( *, * ) "        [-pm], [--pmodel]    :  Flag to precede perturbation model file path."
        write( *, * ) "        [-c],  [--cov]       :  Flag to precede covariance file path."
        write( *, * ) "        [-cf], [--ctrl_fwd]  :  Flag to precede forward control file path."
        write( *, * ) "        [-ci], [--ctrl_inv]  :  Flag to precede inversion control file path."
        write( *, * ) "        [-o],  [--outdir]    :  Flag to precede output directory path."
        write( *, * ) "        [-dm], [--dmodel]    :  Flag to precede output dsigma model file path."
        write( *, * ) "        [-pd], [--predicted] :  Flag to precede output predicted data file path."
        write( *, * ) "        [-jm], [--jmhat]     :  Flag to precede output JmHat data file path."
        write( *, * ) "        [-es], [--esolution] :  Flag to precede binary output e-solution file path."
        write( *, * ) "        [-v],  [--version]   :  Print version."
        write( *, * ) "        [-h],  [--help]      :  Print this information."
        write( *, * ) "        [-tmp],[--template]  :  Create control file templates."
        !
    end subroutine printHelp
    !
    !> Create a template text file called [fwd_ctrl_template.txt] 
    !> with all supported parameters for the Forward Modeling Control File
    !
    subroutine printForwardControlFileTemplate()
        implicit none
        !
        integer :: ios
        !
        open( unit = ioFwdTmp, file = "fwd_ctrl_template.txt", status="unknown", iostat=ios )
        !
        if( ios == 0 ) then
            !
            write( ioFwdTmp, "(A46)" ) "##############################################"
            write( ioFwdTmp, "(A46)" ) "# ModEM Forward Modeling Control File Template"
            write( ioFwdTmp, "(A46)" ) "#     Here are all supported parameters       "
            write( ioFwdTmp, "(A46)" ) "#     Comment or remove to use default value  "
            write( ioFwdTmp, "(A46)" ) "##############################################"
            write( ioFwdTmp, "(A1)" )  "#"
            write( ioFwdTmp, "(A20)" ) "# <Field parameters>"
            write( ioFwdTmp, "(A1)" )  "#"
            write( ioFwdTmp, "(A36)" ) "model_operator_type [MF|SP|SP2] : MF"
            write( ioFwdTmp, "(A1)" )  "#"
            write( ioFwdTmp, "(A19)" ) "# <Grid parameters>"
            write( ioFwdTmp, "(A1)" )  "#"
            write( ioFwdTmp, "(A33)" ) "#grid_header [ModEM|HDF5] : ModEM"
            write( ioFwdTmp, "(A30)" ) "grid_format [SG|MR]       : SG"
            write( ioFwdTmp, "(A1)" )  "#"
            write( ioFwdTmp, "(A20)" ) "# <Model parameters>"
            write( ioFwdTmp, "(A1)" )  "#"
            write( ioFwdTmp, "(A49)" ) "model_method [mirror|fixed height] : fixed height"
            write( ioFwdTmp, "(A39)" ) "model_n_air_layer [10]             : 10"
            write( ioFwdTmp, "(A41)" ) "model_max_height [200.]            : 200."
            write( ioFwdTmp, "(A1)" )  "#"
            write( ioFwdTmp, "(A21)" ) "# <Source parameters>"
            write( ioFwdTmp, "(A1)" )  "#"
            write( ioFwdTmp, "(A66)" ) "source_type_mt [1D|2D]                                        : 1D"
            write( ioFwdTmp, "(A72)" ) "source_type_csem [EM1D|Dipole1D]                              : Dipole1D"
            write( ioFwdTmp, "(A78)" ) "get_1d_from [Fixed|Geometric_mean|Mean_around_Tx|Tx_Position] : Geometric_mean"
            write( ioFwdTmp, "(A1)" )  "#"
            write( ioFwdTmp, "(A21)" ) "# <Solver parameters>"
            write( ioFwdTmp, "(A1)" )  "#"
            write( ioFwdTmp, "(A36)" ) "solver_type [QMR|BICG]         : QMR"
            write( ioFwdTmp, "(A38)" ) "forward_solver_type [IT|IT_DC] : IT_DC"
            write( ioFwdTmp, "(A35)" ) "max_solver_iters [80]          : 80"
            write( ioFwdTmp, "(A35)" ) "max_solver_calls [20]          : 20"
            write( ioFwdTmp, "(A36)" ) "max_divcor_iters [100]         : 100"
            write( ioFwdTmp, "(A37)" ) "tolerance_solver [1E-7]        : 1E-7"
            write( ioFwdTmp, "(A37)" ) "tolerance_divcor [1E-5]        : 1E-5"
            write( ioFwdTmp, "(A1)" ) "#"
            !
            close( ioFwdTmp )
            !
        else
            !
            call errStop( "printInversionControlFileTemplate > opening [fwd_ctrl_template.txt]" )
            !
        endif
        !
    end subroutine printForwardControlFileTemplate
    !
    !> Create a template text file called [inv_ctrl_template.txt] 
    !> with all supported parameters for the Inversion Control File
    !
    subroutine printInversionControlFileTemplate()
        implicit none
        !
        integer :: ios
        !
        open( unit = ioInvTmp, file = "inv_ctrl_template.txt", status="unknown", iostat=ios )
        !
        if( ios == 0 ) then
            !
            write( ioInvTmp, "(A46)" ) "##############################################"
            write( ioInvTmp, "(A46)" ) "# ModEM Inversion Control File Template       "
            write( ioInvTmp, "(A46)" ) "#     Here are all supported parameters       "
            write( ioInvTmp, "(A46)" ) "#     Comment or remove to use default values "
            write( ioInvTmp, "(A46)" ) "##############################################"
            write( ioInvTmp, "(A1)" )  "#"
            write( ioInvTmp, "(A38)" ) "inversion_type [DCG|NLCG]       : NLCG"
            write( ioInvTmp, "(A44)" ) "joint_type [Unweighted|TxBased] : Unweighted"
            write( ioInvTmp, "(A37)" ) "max_inv_iters [100]             : 100"
            write( ioInvTmp, "(A36)" ) "max_grad_iters [20]             : 20"
            write( ioInvTmp, "(A38)" ) "error_tol [1E-3]                : 1E-3"
            write( ioInvTmp, "(A38)" ) "rms_tol [1.05]                  : 1.05"
            write( ioInvTmp, "(A36)" ) "lambda [1.]                     : 1."
            write( ioInvTmp, "(A40)" ) "lambda_tol [1.0e-4]             : 1.0e-4"
            write( ioInvTmp, "(A37)" ) "lambda_div [10.]                : 10."
            write( ioInvTmp, "(A37)" ) "startdm [20.]                   : 20."
            write( ioInvTmp, "(A1)" ) "#"
            !
            close( ioInvTmp )
            !
        else
            !
            call errStop( "printInversionControlFileTemplate > opening [inv_ctrl_template.txt]" )
            !
        endif
        !
    end subroutine printInversionControlFileTemplate
    !
end module CoreComponents
!