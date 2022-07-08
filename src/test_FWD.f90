program ModEM
    !
    use Constants
    use FileUnits
    !
    use ModEMControlFile
    !
    use Grid3D_SG
    !
    use ModelReader
    use ModelReader_Weerachai
    use ModelOperator_MF
    use ModelOperator_File
    use ModelParameterCell_SG
    !
    use DataFileStandard
    !
    use ForwardSolverIT_DC
    !
    use SourceMT_1D
    use SourceMT_2D
    use SourceCSEM_Dipole1D
    !
    use DataHandleFArray
    use DataHandleMT
    use DataHandleCSEM
    !
    ! 
    class( Grid_t ), allocatable           :: main_grid
    class( ModelParameter_t ), allocatable :: model_parameter
    class( ModelOperator_t ), allocatable  :: model_operator
    !
    character(:), allocatable :: control_file_name, model_file_name, data_file_name, modem_job
    logical                   :: has_control_file, has_model_file, has_data_file, verbosis
    !
    !
    modem_job = "unknow"
    !
    call setupDefaultParameters()
    !
    ! Validate arguments, set model_file_name, data_file_name, control_file_name, etc...
    call handleArguments()
    !
    write ( *, * )
    write ( *, * ) "Start ModEM-OO."
    write ( *, * )
    !
    ! Check parameters at the control file
    if( has_control_file ) call handleControlFile()
    !
    ! Execute the modem_job
    call handleJob()
    !
    ! Deallocate remaining main program memory
    call garbageCollector()
    !
    write ( *, * )
    write ( *, * ) "Finish ModEM-OO."
    write ( *, * )
    !
contains
    !
    subroutine ForwardModelling()
        implicit none
        !
        ! Use save ????
        class( ForwardSolver_t ), allocatable, target, save :: forward_solver
        !
        ! Temporary alias pointers
        class( Transmitter_t ), pointer :: Tx
        class( Receiver_t ), pointer    :: Rx
        !
        integer :: iTx, number_of_tx, iRx, iDh
        !
        ! Final array that contains the predicted data of all receivers of all transmitters
        type( Dh_t ), allocatable, dimension(:) :: all_data_handles
        !
        ! Verbosis ...
        write ( *, * ) "    > Start forward modelling."
        !
        ! Reads Model File: instantiates Grid, ModelOperator and ModelParameter
        if( .NOT. has_model_file ) then 
            stop " - Missing Model file!"
        else
            call handleModelFile()
        endif
        !
        ! Reads Data File: instantiates and builds the Data relation between Txs and Txs
        if( .NOT. has_data_file ) then 
            stop " - Missing Data file!"
        else
            call handleDataFile()
        endif
        !
        ! Instantiate the ForwardSolver - Specific type can be chosen via control file
        select case ( forward_solver_type )
            !
            case( FWD_IT_DC )
                allocate( forward_solver, source = ForwardSolverIT_DC_t( model_operator, QMR ) )
                !
            case default
                allocate( forward_solver, source = ForwardSolverIT_DC_t( model_operator, QMR ) )
            !
        end select
        !
        ! Forward Modelling
        !
        number_of_tx = size( transmitters )
        !
        ! Writes the first header of the ESolution binary file, according to the first transmitter
        Tx => getTransmitter(1)
        call writeEsolutionHeader( number_of_tx, Tx%n_pol )
        !
        ! Loop over all Transmitters
        do iTx = 1, number_of_tx
            !
            ! Points the Tx alias to the current loop transmitter
            Tx => getTransmitter( iTx )
            !
            ! Set Transmitter's ForwardSolver
            Tx%forward_solver => forward_solver
            !
            ! Set Transmitter's ForwardSolver Omega(Period) and Conductivity
            call Tx%forward_solver%setFrequency( model_parameter, Tx%period )
            !
            ! Instantiate Transmitter's Source - According to transmitter type or chosen via control file
            select type( Tx )
                !
                class is( TransmitterMT_t )
                    !
                    Tx%source = SourceMT_1D_t( model_operator, model_parameter, Tx%period )
                    !
                class is( TransmitterCSEM_t )
                    !
                    Tx%source = SourceCSEM_Dipole1D_t( model_operator, model_parameter, Tx%period, Tx%location, Tx%dip, Tx%azimuth, Tx%moment )
                    !
            end select
            !
            ! Solve Forward Modeling for this Transmitter
            call Tx%solveFWD()
            !
            ! Loop for each Receiver related to this Transmitter
            do iRx = 1, size( Tx%receiver_indexes )
                !
                ! Points the Rx alias to the current loop Receiver
                Rx => getReceiver( Tx%receiver_indexes( iRx ) )
                !
                ! Verbosis...
                !write( *, * ) "                        Rx Id:", Rx%id, "XYZ:", Rx%location
                !
                ! Calculates Rx predicted data and stores the result in the Receiver
                call Rx%predictedData( Tx )
                !
                ! For each predicted data stored in the Receiver
                do iDh = 1, size( Rx%predicted_data )
                    !
                    ! Store in the final data array
                    call updateDataHandleArray( all_data_handles, getDataHandle( Rx%predicted_data, iDh ) )
                    !
                end do
                !
            enddo
            !
        enddo
        !
        call deallocateTransmitterArray()
        !
        deallocate( forward_solver )
        !
        deallocate( model_operator )
        deallocate( model_parameter )
        deallocate( main_grid )
        !
        ! Verbosis...
        write( *, * ) "    -> Writing Predicted Data to file: [", predicted_data_file_name, "]"
        !
        ! Writes the final data array, with the proper Rx header, to the file <predicted_data_file_name>
        call writeDataHandleArray( all_data_handles )
        !
        call deallocateDataHandleArray( all_data_handles )
        !
        call deallocateReceiverArray()
        !
    end subroutine ForwardModelling
    !
    subroutine handleJob()
        implicit none
        !
        select case ( modem_job )
            !
        case ( "forward" )
            !
            call ForwardModelling()
            !
        case default
            !
            write( *, * ) "    - Unknow job: [", modem_job, "]"
            call printHelp()
            !
        end select
        !
    end subroutine handleJob
    !
    subroutine handleControlFile()
        implicit none
        !
        type( ModEMControlFile_t ) :: control_file
        !
        write( *, * ) "    -> Control File: [", control_file_name, "]"
        !
        ! Instantiate the ControlFile object
        ! Reads control file and sets the options in the Constants module
        control_file = ModEMControlFile_t( ioStartup, control_file_name )
        !
    end subroutine handleControlFile
    !
    !
    subroutine handleModelFile()
        implicit none
        !
        type( ModelReader_Weerachai_t ) :: model_reader
        type( TAirLayers )              :: air_layer
        !
        !
        write( *, * ) "    -> Model File: [", model_file_name, "]"
        !
        ! Read Grid and ModelParameter with ModelReader_Weerachai
        call model_reader%Read( model_file_name, main_grid, model_parameter ) 
        !
        ! Instantiate the ModelOperator object
        select type( main_grid )
            !
            class is( Grid3D_SG_t )
                !
                call main_grid%SetupAirLayers( air_layer, model_method, model_n_air_layer, model_max_height )
                !
                call main_grid%UpdateAirLayers( air_layer%nz, air_layer%dz )
                !
                allocate( model_operator, source = ModelOperator_MF_t( main_grid ) )
                !
                call model_parameter%setMetric( model_operator%metric )
                !
                call model_operator%SetEquations()
                !
                call model_operator%SetCond( model_parameter )
                !
            class default
                stop "Unclassified main_grid"
            !
        end select
        !
    end subroutine handleModelFile
    !
    !
    subroutine handleDataFile()
        implicit none
        !
        integer :: irx, nrx
        class( Receiver_t ), pointer :: Rx
        !
        ! Local object to dealt data, self-destructs at the end of the subroutine
        type( DataFileStandard_t ) :: data_file_standard
        !
        write( *, * ) "    -> Data File: [", data_file_name, "]"
        !
        data_file_standard = DataFileStandard_t( ioStartup, data_file_name )
        !
        nrx = size( receivers )
        !
        if( nrx == data_file_standard%nRx ) then
            !
            write( *, * ) nrx, " Receivers checked!"
            !
            write( *, * ) "    -> Creating Rx evaluation vectors"
            !
            do irx = 1, nrx
                Rx => getReceiver( irx )
                !
                call Rx%evaluationFunction( model_operator )
                !
            enddo
        else
             !
             write(*,*) "Number of Rx mismatched from Header :[", nrx, " and ", data_file_standard%nRx, "]"
             STOP "DataManager.f08: DataManager_ctor()"
             !
        endif
        !
        if( size( transmitters ) == data_file_standard%nTx ) then
             write( *, * ) size( transmitters ), " Transmitters checked!"
             !
             call printTransmitterArray()
        else
             !
             write(*,*) "Number of Tx mismatched from Header :[", size( transmitters ), " and ", data_file_standard%nTx, "]"
             STOP "DataManager.f08: DataManager_ctor()"
             !
        endif
        !
    end subroutine handleDataFile
    !
    subroutine handleArguments()
        implicit none
        !
        character( len=200 ) :: argument
        integer              :: argument_index
        !
        if ( command_argument_count() == 0 ) then
            !
            call printHelp()
            !
        else
            !
            argument_index = 1
            !
            do while(argument_index <= command_argument_count()) 
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
                         if ( len( control_file_name ) > 0 ) has_control_file = .TRUE.
                         !
                         argument_index = argument_index + 2
                         !
                      case ( "-d", "--data" )
                         !
                         call get_command_argument( argument_index + 1, argument )
                         data_file_name = trim( argument )
                         !
                         if ( len( data_file_name ) > 0 ) has_data_file = .TRUE.
                         !
                         argument_index = argument_index + 2
                         !
                      case ( "-f", "--forward" )
                         !
                         modem_job = "forward"
                         !
                         argument_index = argument_index + 1
                         !
                      case ( "-m", "--model" )
                         !
                         call get_command_argument( argument_index + 1, argument )
                         model_file_name = trim( argument )
                         !
                         if ( len( model_file_name ) > 0 ) has_model_file = .TRUE.
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
                         !
                      case ( "--verbosis" )
                         !
                         call printHelp()
                         !
                      case default
                         !
                         write( *, * ) "    - Unknow Argument: [", argument, "]"
                         call printHelp()
                         !
                 end select
                 !
                 argument = ""
                 !
            end do
            !
        end if
        !
    end subroutine handleArguments
    !
    subroutine setupDefaultParameters()
        implicit none
        !
        ! I|O
        predicted_data_file_name = "predicted_data.dat"
        e_solution_file_name     = "esolution.bin"
        has_control_file         = .FALSE.
        has_model_file           = .FALSE.
        has_data_file            = .FALSE.
        verbosis                 = .FALSE.
        !
        ! Solver
        !max_iter = 100
        !tolerance = TOL8
        !
        ! Source
        get_1D_from = "Geometric_mean"
        !
        ! Model
        model_method      = MM_METHOD_FIXED_H
        model_n_air_layer = 10
        model_max_height  = 200.0
        !
        source_type = SRC_MT_1D
        !
        forward_solver_type = FWD_IT_DC
        !
    end subroutine setupDefaultParameters
    !
    subroutine garbageCollector()
        implicit none
        !
        if( allocated( forward_solver_type ) ) deallocate( forward_solver_type )
        if( allocated( source_type ) ) deallocate( source_type )
        if( allocated( model_method ) ) deallocate( model_method )
          if( allocated( get_1D_from ) ) deallocate( get_1D_from )
        if( allocated( predicted_data_file_name ) ) deallocate( predicted_data_file_name )
        if( allocated( e_solution_file_name ) ) deallocate( e_solution_file_name )
        !
        if( allocated( control_file_name ) ) deallocate( control_file_name )
        if( allocated( model_file_name ) ) deallocate( model_file_name )
        if( allocated( data_file_name ) ) deallocate( data_file_name )
        if( allocated( modem_job ) ) deallocate( modem_job )
        !
    end subroutine garbageCollector
    subroutine writeEsolutionHeader( nTx, nMode )
        implicit none
        !
        ! implement separated routine
        integer, intent( in ) :: nTx, nMode
        integer               :: ios
        character(len=20)     :: version
        !
        version = "Modem-OO"
        !
        open( ioESolution, file = e_solution_file_name, action = "write", form = "unformatted", iostat = ios)
        !
        if( ios /= 0 ) then
            write( *, * ) "Error opening file in FileWriteInit: e_solution"
            stop
        else
            !
            write( ioESolution ) version, nTx, nMode, &
            main_grid%nx, main_grid%ny, main_grid%nz, main_grid%nzAir, &
            main_grid%ox, main_grid%oy, main_grid%oz, main_grid%rotdeg
            !
            write( ioESolution ) main_grid%dx
            write( ioESolution ) main_grid%dy
            write( ioESolution ) main_grid%dz
            !
            close( ioESolution )
            !
        endif
        !
        !
    end subroutine writeEsolutionHeader
    !
    recursive subroutine sortByReceiver( data_handle_array, first, last )
        implicit none
        !
        type( Dh_t ), allocatable, dimension(:), intent( inout ) :: data_handle_array
        type( Dh_t ), allocatable :: x_Dh, t_Dh
        class( DataHandle_t ), allocatable :: i_data_handle, j_data_handle, x_data_handle
        integer first, last
        integer i, j
        !
        x_Dh = data_handle_array( (first+last) / 2 )
        x_data_handle = x_Dh%Dh
        i = first
        j = last
        !
        do
            !
            i_data_handle = getDataHandle( data_handle_array, i )
            j_data_handle = getDataHandle( data_handle_array, j )
            !
            do while ( i_data_handle%code < x_data_handle%code )
                i=i+1
                i_data_handle = getDataHandle( data_handle_array, i )
            end do
            do while ( x_data_handle%code < j_data_handle%code )
                j=j-1
                j_data_handle = getDataHandle( data_handle_array, j )
            end do
            if (i >= j) exit
            t_Dh = data_handle_array(i)
            data_handle_array(i) = data_handle_array(j)
            data_handle_array(j) = t_Dh
            i=i+1
            i_data_handle = getDataHandle( data_handle_array, i )
            j=j-1
            j_data_handle = getDataHandle( data_handle_array, j )
        end do
        !
        if (first < i-1) call sortByReceiver( data_handle_array, first, i-1 )
        if (j+1 < last)  call sortByReceiver( data_handle_array, j+1, last )
        !
    end subroutine sortByReceiver
    !
    !
    subroutine writeDataHandleArray( data_handle_array )
        implicit none
        !
        type( Dh_t ), allocatable, dimension(:), intent( inout ) :: data_handle_array
        !
        class( DataHandle_t ), allocatable :: data_handle
        !
        integer :: receiver_type, i, array_size, ios
        !
        ! Order by receiver
        !call sortByReceiver( data_handle_array, 1, size( data_handle_array ) )
        !
        receiver_type = 0
        !
        open( unit = ioPredData, file = predicted_data_file_name, action = "write", form = "formatted", iostat = ios )
        !
        if( ios == 0 ) then
            !
            array_size = size( data_handle_array )
            !
            do i = 1, array_size
                !
                allocate( data_handle, source = getDataHandle( data_handle_array, i ) )
                !
                call writePredictedDataHeader( data_handle, receiver_type )
                !
                ! Instantiate the ModelOperator object
                select type( data_handle )
                    !
                    class is( DataHandleMT_t )
                        !
                        write( ioPredData, "(es12.6, 1X, A, 1X, f15.3, f15.3, f15.3, f15.3, f15.3, 1X, A, 1X, es16.6, es16.6, es16.6)" ) data_handle%period, data_handle%code, R_ZERO, R_ZERO, data_handle%rx_location(1), data_handle%rx_location(2), data_handle%rx_location(3), data_handle%component, data_handle%real, data_handle%imaginary, 1.0
                        !
                    class is( DataHandleCSEM_t )
                        !
                        write( ioPredData, "(A, 1X, es12.6, f15.3, f15.3, f15.3, f15.3, f15.3, f15.3, 1X, A, 1X, f15.3, f15.3, f15.3, 1X, A, 1X, es16.6, es16.6, es16.6)" ) data_handle%dipole, data_handle%period, data_handle%moment, data_handle%azimuth, data_handle%dip, data_handle%tx_location(1), data_handle%tx_location(2), data_handle%tx_location(3), data_handle%code, data_handle%rx_location(1), data_handle%rx_location(2), data_handle%rx_location(3), data_handle%component, data_handle%real, data_handle%imaginary, 1.0
                        !
                    class default
                        stop "Unclassified data_handle"
                    !
                end select
                !
                deallocate( data_handle )
                !
            enddo
            !
            close( ioPredData )
            !
        else
            write( *, * ) "Error opening [", predicted_data_file_name, "] in writeDataHandleArray"
            stop
        end if
        !
    end subroutine writeDataHandleArray
    !
    subroutine writePredictedDataHeader( data_handle, receiver_type )
        implicit none
        !
        !
        class( DataHandle_t ), intent( in ) :: data_handle
        !
        integer, intent( inout ) :: receiver_type
        !
        if( receiver_type /= data_handle%rx_type ) then
            !
            write( ioPredData, "(4A, 40A)" ) "#    ", DATA_FILE_TITLE
            !
            select case( data_handle%rx_type )
                !
                case( 1, 11, 12 )
                    write( ioPredData, "(74A)" ) "#    Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error"
                case( 2 )
                    write( ioPredData, "(60A)" ) "#    Missing title for Full_Interstation_TF"
                case( 3 )
                    write( ioPredData, "(60A)" ) "#    Missing title for Off_Diagonal_Rho_Phase"
                case( 4 )
                    write( ioPredData, "(60A)" ) "#    Missing title for Phase_Tensor"
                case( 5 )
                    write( ioPredData, "(60A)" ) "#    Missing title for Off_Diagonal_Impedance"
                case( 6, 7, 8, 9, 10 )
                    write( ioPredData, "(125A)" ) "#    Tx_Dipole Tx_Period(s) Tx_Moment(Am) Tx_Azi Tx_Dip Tx_X(m) Tx_Y(m) Tx_Z(m) Code X(m) Y(m) Z(m) Component Real Imag Error"
                case default
                    write( *, * ) "unknow receiver type :[", data_handle%rx_type, "]"
                    STOP "test_FWD.f90: writePredictedDataHeader()"
                !
            end select
            !
            write( ioPredData, "(4A, 100A)" ) ">    ", getStringReceiverType( data_handle%rx_type )
            write( ioPredData, "(4A, 100A)" ) ">    ", "exp(-i\omega t)"
            write( ioPredData, "(4A, 100A)" ) ">    ", "[V/m]/[T]"
            write( ioPredData, "(7A, 100A)" ) ">        ", "0.00"
            write( ioPredData, "(7A, 100A)" ) ">        ", "0.000    0.000"
            write( ioPredData, "(A3, i8, i8)" ) ">        ", size( transmitters ), size( receivers )
            !
            receiver_type = data_handle%rx_type
            !
        endif
        !
    end subroutine writePredictedDataHeader
    !
    subroutine printHelp()
        implicit none
        !
        write( *, * ) "ModEM-OO Usage:"
        write( *, * ) ""
        write( *, * ) "    Flags to define a job:"
        write( *, * ) "        [-f], [--forward]    :  Forward modelling."
        write( *, * ) "        [-i], [--inverse]    :  Inversion modelling."
        write( *, * )
        write( *, * ) "    Other arguments:"
        write( *, * ) "        [-d], [--data]       :  Flags for input data file path."
        write( *, * ) "        [-m], [--model]      :  Flags for input model file path."
        write( *, * ) "        [-c], [--control]    :  Flags for user control file path."
        write( *, * ) "        [-v], [--version]    :  Print version."
        write( *, * ) "        [-h], [--help]       :  Print usage information."
        write( *, * ) "        [-pd], [--predicted] :  Output data file path."
        write( *, * ) "        [-es], [--esolution] :  Output binary e-solution file path."
        write( *, * ) "        [--verbosis]         :  Print runtime information."
        !
        write( *, * ) ""
        write( *, * ) "Version 1.0.0"
        !
        stop
        !
    end subroutine printHelp
    !
end program ModEM
