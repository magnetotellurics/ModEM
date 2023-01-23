!
!> Module with the sensitivity routines JMult, JMult_Tx, JMult_T, JMult_T_Tx 
!
module GlobalVariables
    !
    use Constants
    !
    use Grid3D_SG
    !
    use ModelParameterCell_SG
    !
    use ModelOperator_MF
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
	class( Source_t ), allocatable :: fwd_source
	!
    class( Grid_t ), allocatable, target :: main_grid
    !
    class( ModelOperator_t ), allocatable :: model_operator
    !
    class( ForwardSolver_t ), allocatable, target :: forward_solver
    !
    class( ModelCovarianceRec_t ), allocatable :: model_cov
    !
    !> Program control variables
    character(50) :: outdir_name
    !
    character(:), allocatable :: control_file_name
    character(:), allocatable :: model_file_name
    character(:), allocatable :: pmodel_file_name
    character(:), allocatable :: data_file_name
    character(:), allocatable :: dsigma_file_name
    character(:), allocatable :: e_solution_file_name
    character(:), allocatable :: modem_job
    !
    !> Program control flags
    logical :: has_outdir_name
    !
    logical :: has_control_file
    logical :: has_model_file
    logical :: has_pmodel_file
    logical :: has_data_file
    logical :: has_e_solution_file
    logical :: verbosis
    !
    !> Public module routines
    public :: handleModelFile
    public :: handlePModelFile
    public :: handleDataFile
    public :: createOutputDirectory
    public :: garbageCollector
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
                    write( *, "( i20, f20.5 )" ) i, air_layer%dz(i)
                enddo
                !
                write( *, * ) "          Air layers from the method [", trim( air_layer%method ), "]"
                !
                write( *, "( a39, f12.5, a4 )" ) "          Top of the air layers is at ", sum( air_layer%Dz ) / 1000, " km."
                !
                write( *, "( a31, f16.5, a2, f16.5, a2, f16.5, a4, f16.5 )" ) "          o(x,y,z) * rotDeg: (", main_grid%ox, ", ", main_grid%oy, ", ", main_grid%oz, ") * ", main_grid%rotDeg
                !
                allocate( model_operator, source = ModelOperator_MF_t( main_grid ) )
                !
                call model_operator%setEquations()
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
                call transmitters( i_tx )%Tx%print()
            enddo
            !
        else
             !
             write( *, * ) "Error: Number of Tx mismatched from Header :[", n_tx, " and ", data_file_standard%n_tx, "]"
             stop
             !
        endif
        !
        write( *, * ) "          Checked ", countDataGroupTxArray( all_measured_data ), " DataGroups."
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
            select case ( inversion_algorithm )
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
                    stop "Error: jobInversion > Undefined inversion_algorithm"
                    !
            end select
            !
            !write( *, * ) "outdir_name: [", trim( outdir_name ), "]"
            !
        endif
        !
        call system( "mkdir -p "//outdir_name )
        !
    end subroutine createOutputDirectory
    !
    !> Deallocate remaining unallocated global memory
    !
    subroutine garbageCollector()
        implicit none
        !
        !> Deallocate global array of measured data
        if( allocated( all_measured_data ) ) call deallocateDataGroupTxArray( all_measured_data )
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
        if( allocated( inversion_algorithm ) ) deallocate( inversion_algorithm )
        !
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
end module GlobalVariables
!