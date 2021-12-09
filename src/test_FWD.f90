program ModEM
   !
   use DataManager
   use ModEMControlFile
   !
   use ModelReader
   use ModelReader_Weerachai
   use ModelOperator_MF
   use ModelParameterCell_SG
   
   use Grid3D_SG
   !
   use MetricElements_CSG
   !
   use Solver_QMR
   use Solver_PCG
   !
   use DivergenceCorrection
   !
   use ForwardSolverFromFile
   use ForwardSolverDC
   !
   use SourceMT_1D
   use SourceMT_2D
   !
   use PreConditioner_MF_CC
   use PreConditioner_MF_DC
   !
   class( ModEMControlFile_t ), allocatable :: control_file
   !
   class( Grid_t ), allocatable           :: main_grid
   class( MetricElements_t ), allocatable :: metric_elemements
   class( ModelParameter_t ), allocatable :: model_parameter
   class( ModelOperator_t ), allocatable  :: model_operator
   !
   character(200) :: argument
   integer        :: argument_index
   !
   character(:), allocatable :: control_file_name, model_file_name, data_file_name, modem_job
   logical                   :: has_control_file = .false., has_model_file = .false., has_data_file = .false.
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
   subroutine forward()
      implicit none
      !
      ! These objects must be instantiated only once in the Master
      type( DivergenceCorrection_t ), target, save :: divergence_correction
      !
      type( PreConditioner_MF_DC_t ), target, save :: preconditioner_dc
      !
      type( Solver_PCG_t ), target, save           :: solver_pcg
      !
      ! These objects are frequency dependent,
      ! must be instantiated on each Worker
      type( Solver_QMR_t ), target, save                  :: solver_qmr
      !
      class( ForwardSolver_t ), allocatable, target, save :: fwd_solver
      !
      class( Source_t ), allocatable, target, save        :: fwd_source 
      !
      type( PreConditioner_MF_CC_t ), target, save        :: preconditioner_cc
      !
      ! Temporary alias pointers
      class( Transmitter_t ), allocatable :: Tx
      class( Receiver_t ), allocatable    :: Rx
      !
      ! Local variables
      integer :: iTx, nTx, iRx, nRx
      real( kind=prec ) :: omega
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
      ! PreConditioners need to be instantiated within the selection case
      ! as they receive a specific ModelOperator
      select type( model_operator )
         class is( ModelOperator_MF_t )
           !
           ! PreConditioner CC
           preconditioner_cc = PreConditioner_MF_CC_t( model_operator )
           !
           ! PreConditioner DC
           preconditioner_dc = PreConditioner_MF_DC_t( model_operator )
           !
      end select
      !
      ! Specific Solver QMR
      solver_qmr = Solver_QMR_t( preconditioner_cc )
      !
      ! Specific Solver PCG
      solver_pcg = Solver_PCG_t( preconditioner_dc )
      !
      ! DivergenceCorrection has only one type for now
      divergence_correction = DivergenceCorrection_t( solver_pcg )
      !
      ! ForwardSolver - Chosen from control file
      select case ( forward_solver_type )
         !
         case( FWD_FILE )
            allocate( fwd_solver, source = ForwardSolverFromFile_t( model_operator ) )
            !
         case( FWD_DC )
            allocate( fwd_solver, source = ForwardSolverDC_t( solver_qmr, divergence_correction ) )
            !
         case default
            allocate( fwd_solver, source = ForwardSolverDC_t( solver_qmr, divergence_correction ) )
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
         ! Write the proper header in the 'predicted_data.dat' file
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
   end subroutine forward
   !
   !
   subroutine handleJob()
      implicit none
      !
      select case ( modem_job )
         !
      case ( "forward" )
         !
         call forward()
         !
      case default
         !
         write( *, * ) "   - Unknow job: [", trim( modem_job ), "]"
         call printHelp()
         stop
         !
      end select
      !
   end subroutine handleJob
   !
   subroutine handleControlFile()
      implicit none
      !
      write( *, * ) "   -> Control File: [", control_file_name, "]"
      !
      ! Instantiate the ControlFile object
      ! Read control file options
      allocate( control_file, source = ModEMControlFile_t( ioStartup, control_file_name ) )
      !
      ! SET TYPE OF SOLVER, PRECONDITIONER... AND POINT THEM
      !
      ! Dispose memory controls that will no longer be used.
      deallocate( control_file )
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
			 allocate( metric_elemements, source = MetricElements_CSG_t( main_grid ) )
             !
             allocate( model_operator, source = ModelOperator_MF_t( metric_elemements ) )
             !
             call model_operator%SetEquations()
             !
			 call model_parameter%setMetric( metric_elemements )
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
      if ( command_argument_count() == 0 ) then
         !
         call printHelp()
         stop
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
                   stop
                   !
                 case default
                   !
                   write( *, * ) "   - Unknow Argument: [", trim( argument ), "]"
                   call printHelp()
                   stop
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
         open( 666, file = 'predicted_data.dat', action='write', position='append' )
         !
         write( 666, * ) '#', DATA_FILE_TITLE
         write( 666, * ) '#', Tx%DATA_TITLE
         !
         close( 666 )
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
   end subroutine printHelp

end program ModEM
