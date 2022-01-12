program ModEM
   !
   use DataManager
   use ModEMControlFile
   !
   use ModelReader
   use ModelReader_Weerachai
   use ModelOperator_MF
   use ModelOperator_File
   use ModelParameterCell_SG
   !
   use Grid3D_SG
   !
   use ForwardSolverFromFile
   use ForwardSolverIT_DC
   !
   use SourceMT_1D
   use SourceMT_2D
   !
   ! 
   class( Grid_t ), allocatable           :: main_grid
   class( ModelParameter_t ), allocatable :: model_parameter
   class( ModelOperator_t ), allocatable  :: model_operator
   !
   character(:), allocatable :: control_file_name, model_file_name, data_file_name, modem_job
   logical                   :: has_control_file = .false., has_model_file = .false., has_data_file = .false.
   !
   !
   modem_job = "unknow"
   !
   ! Validate arguments, set model_file_name, data_file_name
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
   write ( *, * )
   write ( *, * ) "Finish ModEM-OO."
   write ( *, * )
   !
contains
   !
   subroutine ForwardModelling()
      implicit none
      !
      ! These objects are frequency dependent,
      ! must be instantiated on each Worker
      class( ForwardSolver_t ), allocatable, target, save :: fwd_solver
      !
      class( Source_t ), allocatable, target, save        :: fwd_source 
      !
      ! Temporary alias pointers
      class( Transmitter_t ), allocatable :: Tx
      class( Receiver_t ), allocatable    :: Rx
      !
      ! Local variables
      integer :: iTx, nTx, iRx, nRx
      character(:), allocatable :: transmitter_type
      !
      ! Verbosis
      write ( *, * ) "   > Start forward modelling."
      !
      ! Tx type for predicted_data header changes
      transmitter_type = "Unknow"
      !
      ! Reads Model File: instantiates Grid, ModelOperator and ModelParameter
      if( .not. has_model_file ) then 
         stop " - Missing Model file!"
      else
         call handleModelFile()
      endif
      !
      ! Reads Data File: instantiates and builds the Data relation between Txs and Txs
      if( .not. has_data_file ) then 
         stop " - Missing Data file!"
      else
         call handleDataFile()
      endif
      !
      ! High-level object instantiation
      ! Some types are chosen from the control file
      !
      ! ForwardSolver - Chosen from control file
      select case ( forward_solver_type )
         !
         case( FWD_FILE )
            fwd_solver = ForwardSolverFromFile_t( model_operator )
            !
         case( FWD_IT )
            fwd_solver = ForwardSolverIT_t( model_operator, "QMR" )
            !
         case( FWD_IT_DC )
            fwd_solver = ForwardSolverIT_DC_t( model_operator, "QMR" )
            !
         case default
            fwd_solver = ForwardSolverIT_DC_t( model_operator, "QMR" )
         !
      end select
      !
      ! Source - Chosen from control file
      select case ( source_type )
         !
         case( SRC_MT_1D )
            allocate( fwd_source, source = SourceMT_1D_t( model_operator, model_parameter ) )
            !
         case( SRC_MT_2D )
            allocate( fwd_source, source = SourceMT_2D_t( model_operator, model_parameter ) )
            !
         case default
            allocate( fwd_source, source = SourceMT_1D_t( model_operator, model_parameter ) )
            !
      end select
      !
      ! Forward Modelling
      !
      ! Loop over all Transmitters
      nTx = transmitters%size()
      !
      do iTx = 1, nTx
         !
         ! Temporary Transmitter alias
         Tx = transmitters%get( iTx )
         !
         ! Verbosis...
         write( *, * ) "   Tx Id:", Tx%id, "Period:", int( Tx%period )
         !
         ! According to its type,
         ! write the proper header in the 'predicted_data.dat' file
         call writePredictedDataHeader( Tx, transmitter_type )
         !
         ! Tx points to its due Source
         call Tx%setSource( fwd_source )
         !
         ! Tx points to its due ForwardSolver
         call Tx%setForwardSolver( fwd_solver )
         !
         ! Solve Tx Forward Modelling
         call Tx%solveFWD()
         !
         ! Loop over Receivers of each Transmitter
         nRx = Tx%getNRx()
         !
         do iRx = 1, nRx
            !
            ! Temporary Receiver alias
            Rx = receivers%get( Tx%get( iRx ) )
            !
            ! Verbosis...
            !write( *, * ) "                  Rx Id:", Rx%id, "XYZ:", Rx%location
            !
            ! Calculate Rx Predicted Data
            call Rx%predictedData( model_operator, Tx )
            !
            deallocate( Rx )
            !
         enddo
         !
         deallocate( Tx )
         !
      enddo
      !
      deallocate( data_groups )
      !
   end subroutine ForwardModelling
   !
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
         write( *, * ) "   - Unknow job: [", trim( modem_job ), "]"
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
      write( *, * ) "   -> Control File: [", control_file_name, "]"
      !
      ! Instantiate the ControlFile object
      ! Reads control file and sets the options in the Constants module
      control_file = ModEMControlFile_t( ioStartup, control_file_name )
      !
   end subroutine handleControlFile
   !
   !
   subroutine handleDataFile()
      implicit none
      !
      ! Local object to dealt data, self-destructs at the end of the subroutine
      type( DataManager_t ) :: data_manager
      !
      write( *, * ) "   -> Data File: [", data_file_name, "]"
      !
      data_manager = DataManager_t( data_file_name )
      !
   end subroutine handleDataFile
   !
   subroutine handleModelFile()
      implicit none
      !
      ! It remains to standardize ????
      type( ModelReader_Weerachai_t ) :: model_reader
      !
      character(:), allocatable :: fname
      !
      fname = "/mnt/c/Users/protew/Desktop/ON/GITLAB_PROJECTS/modem-oo/inputs/Full_A_Matrix_TinyModel"
      !fname = "/Users/garyegbert/Desktop/ModEM_ON/modem-oo/inputs/Full_A_Matrix_TinyModel"
      !
      write( *, * ) "   -> Model File: [", model_file_name, "]"
      !
      ! Read Grid and ModelParameter with ModelReader_Weerachai
      call model_reader%Read( model_file_name, main_grid, model_parameter ) 
      !
      ! Instantiate the ModelOperator object
      select type( main_grid )
         !
         class is( Grid3D_SG_t )
             ! 
             !model_operator = ModelOperator_File_t( main_grid, fname )
             !
             model_operator = ModelOperator_MF_t( main_grid )
             !
             ! complete model operator setup
             call model_operator%SetEquations()
             !
             call model_parameter%setMetric( model_operator%metric )
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
   subroutine handleArguments()
      implicit none
      !
      character(200) :: argument
      integer        :: argument_index
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
                   if ( len( control_file_name ) > 0 ) has_control_file = .true.
                   !
                   argument_index = argument_index + 2
                   !
                 case ( "-d", "--data" )
                   !
                   call get_command_argument( argument_index + 1, argument )
                   data_file_name = trim( argument )
                   !
                   if ( len( data_file_name ) > 0 ) has_data_file = .true.
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
                   model_file_name = trim(argument)
                   !
                   if ( len( model_file_name ) > 0 ) has_model_file = .true.
                   !
                   argument_index = argument_index + 2
                   !
                 case ( "-v", "--version" )
                   !
                   write( *, * ) "   + ModEM-OO version 1.0.0"
                   stop
                   !
                 case ( "-h", "--help" )
                   !
                   call printHelp()
                   !
                 case default
                   !
                   write( *, * ) "   - Unknow Argument: [", trim( argument ), "]"
                   call printHelp()
                   !
             end select
             !
         end do
         !
      end if
      !
   end subroutine handleArguments
   !
   !
   subroutine writePredictedDataHeader( Tx, transmitter_type )
      implicit none
      !
      class( Transmitter_t ), intent( in )       :: Tx
      character(:), allocatable, intent( inout ) :: transmitter_type
      !
      if( ( index( transmitter_type, "Unknow" ) /= 0 ) .OR. transmitter_type /= Tx%getType() ) then
         !
         transmitter_type = Tx%getType()
         !
         open( ioPredData, file = 'predicted_data.dat', action='write', position='append' )
         !
         write( ioPredData, * ) '#', DATA_FILE_TITLE
         write( ioPredData, * ) '#', Tx%DATA_TITLE
         !
         close( ioPredData )
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
      write( *, * ) "   Flags to define a job:"
      write( *, * ) "      [-f], [--forward] : Forward modelling."
      write( *, * ) "      [-i], [--inverse] : Inversion modelling."
      write( *, * )
      write( *, * ) "   Other arguments:"
      write( *, * ) "      [-d], [--data]    : Flags for data file path."
      write( *, * ) "      [-m], [--model]   : Flags for model file path."
      write( *, * ) "      [-c], [--control] : Flags for user control file path."
      write( *, * ) "      [-v], [--version] : Print version."
      write( *, * ) "      [-h], [--help]    : Print usage information."
      !
      write( *, * ) ""
      write( *, * ) "Version 1.0.0"
      !
      stop
      !
   end subroutine printHelp

end program ModEM
