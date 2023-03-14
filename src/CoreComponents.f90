!
!> Module with the sensitivity routines serialJMult, JMult_Tx, serialJMult_T, JMult_T_Tx 
!
module CoreComponents
    !
    use Constants
    !
    use ForwardControlFile
    use InversionControlFile
    !
    use Grid3D_SG
    !
    use ModelOperator_MF
    use ModelOperator_SP
    !
    use ModelParameterCell_SG
    !
    use ModelCovarianceRec
    !
    use ForwardSolverIT_DC
    !
    use SourceMT_1D
    use SourceMT_2D
    use SourceCSEM_Dipole1D
    use SourceInteriorForce
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
    !> Global Variables
    type( ForwardControlFile_t ), allocatable :: fwd_control_file
    type( InversionControlFile_t ), allocatable :: inv_control_file
    !
    class( Grid_t ), allocatable, target :: main_grid
    !
    class( ModelOperator_t ), allocatable :: model_operator
    !
    class( ForwardSolver_t ), allocatable, target :: forward_solver
    !
    class( ModelCovarianceRec_t ), allocatable :: model_cov
    !
    !> Program Control Variables
    character(50) :: outdir_name
    !
    character(:), allocatable :: fwd_control_file_name
    character(:), allocatable :: inv_control_file_name
    character(:), allocatable :: model_file_name
    character(:), allocatable :: pmodel_file_name
    character(:), allocatable :: data_file_name
    character(:), allocatable :: dsigma_file_name
    character(:), allocatable :: e_solution_file_name
    character(:), allocatable :: modem_job
    !
    !> Program Control Flags
    logical :: has_outdir_name
    logical :: has_fwd_control_file
    logical :: has_inv_control_file
    logical :: has_model_file
    logical :: has_pmodel_file
    logical :: has_data_file
    logical :: has_e_solution_file
    logical :: verbosis
    !
    !> Public Module Routines
    public :: handleModelFile
    public :: handlePModelFile
    public :: handleDataFile
    public :: handleForwardControlFile
    public :: handleInversionControlFile
    public :: handleArguments
    public :: setupDefaultParameters
    public :: createOutputDirectory
    public :: garbageCollector
    public :: printUsage
    public :: printHelp
    public :: printForwardControlFileTemplate
    public :: printInversionControlFileTemplate
    !
contains
    !
    !> Read Model File and instantiate global variables: main_grid, model_operator and sigma0
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
        !> Initialize main_grid and sigma0 with ModelReader(Only ModelReader_Weerachai by now????)
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
                    write( *, "( i20, f20.5 )" ) i, air_layer%dz(i)
                enddo
                !
                write( *, * ) "          Air layers from the method [", trim( air_layer%method ), "]"
                !
                write( *, "( a39, f12.5, a4 )" ) "          Top of the air layers is at ", sum( air_layer%Dz ) / 1000, " km."
                !
                write( *, "( a31, f16.5, a2, f16.5, a2, f16.5, a4, f16.5 )" ) "          o(x,y,z) * rotDeg:(", main_grid%ox, ", ", main_grid%oy, ", ", main_grid%oz, ") * ", main_grid%rotDeg
                !
                allocate( model_operator, source = ModelOperator_MF_t( main_grid ) )
                !
                call model_operator%setEquations
                !
                call sigma0%setMetric( model_operator%metric )
                !
                call model_operator%setCond( sigma0 )
                !
            class default
                stop "Error: handleModelFile > Unclassified main_grid"
            !
        end select
        !
    end subroutine handleModelFile
    !
    !> Read Perturbation Model File: instantiate pmodel
    !
    subroutine handlePModelFile( pmodel )
        implicit none
        !
        class( ModelParameter_t ), allocatable, intent( out ) :: pmodel
        !
        type( ModelReader_Weerachai_t ) :: model_reader
        !
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
        integer :: i_tx, i_rx, n_rx, n_tx
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
        n_tx = size( transmitters )
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
             write( *, * ) "Error: Number of Tx mismatched from Header :[", n_tx, " and ", data_file_standard%n_tx, "]"
             stop
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
    !> Read Forward Modeling Control File: set the parameter values defined in this files
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
    !> Read Inversion Control File: set the parameter values defined in this files
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
    !> Treat all the arguments passed in the execution command line
    !
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
                         modem_job = "Forward"
                         !
                         argument_index = argument_index + 1
                         !
                      case( "-i", "--inversion" )
                         !
                         modem_job = "Inversion"
                         !
                         argument_index = argument_index + 1
                         !
                      case( "-j", "--jmult" )
                         !
                         modem_job = "serialJMult"
                         !
                         argument_index = argument_index + 1
                         !
                      case( "-jt", "--jmult_t" )
                         !
                         modem_job = "serialJMult_T"
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
                         write( *, * ) "    + ModEM-OO version 1.0.0"
                        stop
                         !
                      case( "-h", "--help" )
                         !
                         call printHelp()
                         call printUsage()
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
    !> Set default values for the control variables
    !> and Forward Modeling parameters
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
    !> Create a directory to store output files
    !> With a name standardized by date and runtime or specified by the input argument [-o]
    !
    subroutine createOutputDirectory()
        implicit none
        !
        character(8) str_date
        character(6) str_time
        !
        if( .NOT. has_outdir_name ) then
            !
            !>
            call date_and_time( str_date, str_time )
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
                    stop "Error: jobInversion > Undefined inversion_type"
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
        !
        if( allocated( forward_solver_type ) ) deallocate( forward_solver_type )
        if( allocated( source_type ) ) deallocate( source_type )
        if( allocated( model_method ) ) deallocate( model_method )
        if( allocated( get_1D_from ) ) deallocate( get_1D_from )
        if( allocated( predicted_data_file_name ) ) deallocate( predicted_data_file_name )
        if( allocated( jmhat_data_file_name ) ) deallocate( jmhat_data_file_name )
        if( allocated( e_solution_file_name ) ) deallocate( e_solution_file_name )
        !
        if( allocated( fwd_control_file_name ) ) deallocate( fwd_control_file_name )
        if( allocated( fwd_control_file ) ) deallocate( fwd_control_file )
        if( allocated( inv_control_file_name ) ) deallocate( inv_control_file_name )
        if( allocated( inv_control_file ) ) deallocate( inv_control_file )
        !
        if( allocated( model_file_name ) ) deallocate( model_file_name )
        if( allocated( pmodel_file_name ) ) deallocate( pmodel_file_name )
        if( allocated( data_file_name ) ) deallocate( data_file_name )
        if( allocated( dsigma_file_name ) ) deallocate( dsigma_file_name )
        if( allocated( modem_job ) ) deallocate( modem_job )
        !
    end subroutine garbageCollector
    !
    !> Print basic information on how to use the program on the screen
    !
    subroutine printUsage()
        implicit none
        !
        write( *, * ) "ModEM Minimal Usage:"
        write( *, * ) ""
        write( *, * ) "    Forward Modeling:"
        write( *, * ) "        <ModEM> -f -m <rFile_Model> -d <rFile_Data>"
        write( *, * ) "        Output:"
        write( *, * ) "        - 'all_predicted_data.dat' or the path specified by      [-pd]"
        write( *, * ) ""
        write( *, * ) "    JMult:"
        write( *, * ) "        <ModEM> -j -m <rFile_Model> -pm <rFile_pModel> -d <rFile_Data>"
        write( *, * ) "        Output:"
        write( *, * ) "        - 'jmhat.dat' or the path specified by                   [-jm]"
        write( *, * ) ""
        write( *, * ) "    JMult_T:"
        write( *, * ) "        <ModEM> -jt -m <rFile_Model> -d <rFile_Data>"
        write( *, * ) "        Output:"
        write( *, * ) "        - 'dsigma.rho' or the path specified by                  [-dm]"
        write( *, * ) ""
        write( *, * ) "    Inversion:"
        write( *, * ) "        <ModEM> -i -m <rFile_Model> -d <rFile_Data>"
        write( *, * ) "        Output:"
        write( *, * ) "        - directory named 'Output_<date>_<time>' or specified by [-o]"
        !
    end subroutine printUsage
    !
    !> Print all the supported input parameters  on the screen
    !
    subroutine printHelp()
        implicit none
        !
        write( *, * ) "ModEM Options:"
        write( *, * ) ""
        write( *, * ) "    Flags to define a job:"
        write( *, * ) "        [-f],  [--forward]   :  Forward Modeling."
        write( *, * ) "        [-j],  [--jmult]     :  serialJMult."
        write( *, * ) "        [-jt], [--jmult_t]   :  Transposed serialJMult."
        write( *, * ) "        [-i],  [--inversion] :  Inversion."
        write( *, * )
        write( *, * ) "    Other arguments:"
        write( *, * ) "        [-d],  [--data]      :  Flag to precede data file path."
        write( *, * ) "        [-m],  [--model]     :  Flag to precede model file path."
        write( *, * ) "        [-pm], [--pmodel]    :  Flag to precede perturbation model file path."
        write( *, * ) "        [-cf], [--ctrl_fwd]  :  Flag to precede forward control file path."
        write( *, * ) "        [-ci], [--ctrl_inv]  :  Flag to precede inversion control file path."
        write( *, * ) "        [-o],  [--outdir]    :  Flag to precede output directory path."
        write( *, * ) "        [-dm], [--dmodel]    :  Flag to precede output dsigma model file path."
        write( *, * ) "        [-pd], [--predicted] :  Flag to precede output predicted data file path."
        write( *, * ) "        [-jm], [--jmhat]     :  Flag to precede output JmHat data file path."
        write( *, * ) "        [-es], [--esolution] :  Flag to precede binary output e-solution file path."
        write( *, * ) "        [-v],  [--version]   :  Print version."
        write( *, * ) "        [-h],  [--help]      :  Print usage information."
        write( *, * ) "        [-tmp],[--template]  :  Create control file templates."
        !
        write( *, * ) ""
        write( *, * ) "Version 1.0.0"
        !
    end subroutine printHelp
    !
    !> Create a text file called [fwd_ctrl_template.txt] 
    !> with the supported content for set the Forward Modeling parameters
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
            write( ioFwdTmp, "(A46)" ) "#     Comment or remove to use default        "
            write( ioFwdTmp, "(A46)" ) "##############################################"
            write( ioFwdTmp, "(A1)" )  "#"
            write( ioFwdTmp, "(A18)" ) "# Grid parameters:"
            write( ioFwdTmp, "(A1)" )  "#"
            write( ioFwdTmp, "(A42)" ) "#grid_header [ModEM|HDF5]          : ModEM"
            write( ioFwdTmp, "(A39)" ) "#grid_type [SG|MR]                 : SG"
            write( ioFwdTmp, "(A1)" )  "#"
            write( ioFwdTmp, "(A19)" ) "# Model parameters:"
            write( ioFwdTmp, "(A1)" )  "#"
            write( ioFwdTmp, "(A49)" ) "model_method [mirror|fixed height] : fixed height"
            write( ioFwdTmp, "(A39)" ) "model_n_air_layer [10]             : 10"
            write( ioFwdTmp, "(A42)" ) "model_max_height [200.0]           : 200.0"
            write( ioFwdTmp, "(A1)" )  "#"
            write( ioFwdTmp, "(A20)" ) "# Source parameters:"
            write( ioFwdTmp, "(A1)" )  "#"
            write( ioFwdTmp, "(A39)" ) "source [1D|2D]                     : 1D"
            write( ioFwdTmp, "(A1)" )  "#"
            write( ioFwdTmp, "(A20)" ) "# Solver parameters:"
            write( ioFwdTmp, "(A1)" )  "#"
            write( ioFwdTmp, "(A39)" ) "QMR_iters [40]                     : 40"
            write( ioFwdTmp, "(A39)" ) "BCG_iters [80]                     : 80"
            write( ioFwdTmp, "(A39)" ) "max_divcor_calls [20]              : 20"
            write( ioFwdTmp, "(A40)" ) "max_divcor_iters [100]             : 100"
            write( ioFwdTmp, "(A41)" ) "tolerance_qmr [1E-7]               : 1E-7"
            write( ioFwdTmp, "(A41)" ) "tolerance_divcor [1E-5]            : 1E-5"
            write( ioFwdTmp, "(A42)" ) "forward_solver [IT|IT_DC]          : IT_DC"
            write( ioFwdTmp, "(A1)" )  "#"
            !
            close( ioFwdTmp )
            !
        else
            !
            stop "Error: printInversionControlFileTemplate > opening [inv_control_file_template.txt]"
            !
        endif
        !
    end subroutine printForwardControlFileTemplate
    !
    !> Create a text file called [inv_ctrl_template.txt] 
    !> with the supported content for set the Inversion parameters
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
            write( ioInvTmp, "(A39)" ) "#######################################"
            write( ioInvTmp, "(A39)" ) "# ModEM Inversion Control File Template"
            write( ioInvTmp, "(A39)" ) "#     Here are all supported parameters"
            write( ioInvTmp, "(A39)" ) "#     Comment or remove to use default "
            write( ioInvTmp, "(A39)" ) "#######################################"
            write( ioInvTmp, "(A1)" )  "#"
            write( ioInvTmp, "(A31)" ) "inversion_type [DCG|NLCG] : DCG"
            write( ioInvTmp, "(A29)" ) "max_inv_iters [5]         : 5"
            write( ioInvTmp, "(A30)" ) "max_grad_iters [20]       : 20"
            write( ioInvTmp, "(A32)" ) "error_tol [1E-3]    : 1E-3"
            write( ioInvTmp, "(A32)" ) "rms_tol [1.05]      : 1.05"
            write( ioInvTmp, "(A31)" ) "lambda [10.]              : 10."
            write( ioInvTmp, "(A1)" )  "#"
            !
            close( ioInvTmp )
            !
        else
            !
            stop "Error: printInversionControlFileTemplate > opening [inv_control_file_template.txt]"
            !
        endif
        !
    end subroutine printInversionControlFileTemplate
    !
end module CoreComponents
!