program Test_FWD
   !
   use DataManager
   use ModEMControlFile
   !
   use ModelReader
   use ModelReader_Weerachai
   use ModelOperator_MF
   use ModelParameterCell_SG
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
      type( DivergenceCorrection_t ), target, save  :: divergence_correction
      !
      type( PreConditioner_MF_CC_t ), target, save  :: preconditioner_cc
      !
      type( PreConditioner_MF_DC_t ), target, save  :: preconditioner_dc
      !
      type( Solver_PCG_t ), target, save            :: solver_pcg
      !
      type( Solver_QMR_t ), target, save            :: solver_qmr
      !
      class( ForwardSolver_t ), allocatable, target, save :: fwd_solver
      !
      class( Source_t ), allocatable, target, save        :: fwd_source 
      !
      !
      class( Transmitter_t ), allocatable :: tx
      class( Receiver_t ), allocatable    :: rx
      !
      integer :: iTx, nTx, iRx, nRx
      real( kind=prec ) :: omega
      character(:), allocatable :: transmitter_type
      !
      transmitter_type = "Unknow"
      !
      if( .not. has_data_file ) stop " - Missing data file!"
      !
      if( .not. has_model_file ) stop " - Missing model file!"
      !
      write ( *, * ) "   > Start forward modelling."
      !
      ! Load model, define main_grid, model_operator, model_parameter
      call handleModelFile()
      !
      ! Load data, construct the relation between txs and rxs
      call handleDataFile()
      !
      ! Set right object types from control file:
      !
      select type( model_operator )
         class is( ModelOperator_MF_t )
            !
            ! 1) PreConditioner CC
            preconditioner_cc = PreConditioner_MF_CC_t( model_operator )
            !
            ! 2) PreConditioner DC
            preconditioner_dc = PreConditioner_MF_DC_t( model_operator )
            !
      end select
      !
      ! 3) Solver PCG
      solver_pcg = Solver_PCG_t( preconditioner_dc )
      !
      ! 4) DivergenceCorrection
      divergence_correction = DivergenceCorrection_t( solver_pcg )
      !
      ! 5) Solver QMR
      solver_qmr = Solver_QMR_t( preconditioner_cc )
      !
      ! 6) ForwardSolver
      !allocate( fwd_solver, source = ForwardSolverFromFile_t() )
      allocate( fwd_solver, source = ForwardSolverDC_t( solver_qmr, divergence_correction ) )
      !
      ! Allocate Forward Solver e_solution
      select type( grid => model_operator%grid )
          class is( Grid3D_SG_t )
              allocate( fwd_solver%e_solution, source=cVector3D_SG_t( grid, EDGE ) )
      end select
      !
      ! 7) Source
      allocate( fwd_source, source = SourceMT_1D_t( model_operator, model_parameter ) )
      !
      ! Loop over all Transmitters
      nTx = transmitters%size()
      write( *, * ) nTx, " Transmitters"
      !
      do iTx = 1, nTx
         !
         ! Temporary generic Tx
         tx = transmitters%get( iTx )
         !
         ! verbosis
         write( *, * ) "   Tx Id:", tx%id, "Period:", int( tx%period )
         !
         ! Write the proper header in the 'predicted_data.dat' file
         call writePredictedDataHeader( tx, transmitter_type )
         !
         ! Set forward solver on the Tx
       call fwd_solver%setSolver( solver_qmr )
         allocate( tx%forward_solver, source = fwd_solver )
         !
         ! Set source on the Tx
         allocate( tx%source, source = fwd_source )
         !
         ! Solve Tx Forward modelling
         call tx%solveFWD()
         !
         ! Loop over receivers of each transmitter
         nRx = tx%getNRx()
         write( *, * ) nRx, " Receivers"
         !
         do iRx = 1, nRx
            !
            rx = receivers%get( tx%get( iRx ) )
            !
            !write( *, * ) "                   Rx Id:", rx%id, "XYZ:", rx%location
            !
            call rx%predictedData( model_operator, tx )
            !
            deallocate( rx )
            !
         enddo
         !
         deallocate( tx )
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
      type( DataManager_t )    :: data_manager
     !
     data_manager = DataManager_t( data_file_name )
      !
      write( *, * ) "   -> Data File: [", data_file_name, "]"
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
              allocate( model_operator, source = ModelOperator_MF_t( main_grid ) )
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
      ! Instantiate the PreConditioner object
    !allocate( preconditioner, source = PreConditioner_MF_DC_t( model_operator ) )
     !
      write( *, * ) "   -> Readed Model File: [", model_file_name, "]"
      !
   end subroutine handleModelFile
   !
   !
   subroutine handleArguments()
      implicit none
      !
      if (command_argument_count() == 0 ) then
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
      !
   end subroutine handleArguments
   !
   !
   subroutine writePredictedDataHeader( tx, transmitter_type )
      implicit none
    class( Transmitter_t ), intent( in )          :: tx
    character(:), allocatable, intent( inout ) :: transmitter_type
      !
      if( ( index( transmitter_type, "Unknow" ) /= 0 ) .OR. transmitter_type /= tx%getType() ) then
         !
         transmitter_type = tx%getType()
         !
         open( 666, file = 'predicted_data.dat', action='write', position='append' )
         !
         write( 666, * ) '#', DATA_FILE_TITLE
         write( 666, * ) '#', tx%DATA_TITLE
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
      write( *, * ) "Arguments to define the ModEM task:"
      write( *, * ) "[-f], [--forward] : Forward modelling."
      write( *, * )
      write( *, * ) "Other arguments:"
      write( *, * ) "[-d], [--data]      : Preceding argument for data file path."
      write( *, * ) "[-m], [--model]    : Preceding argument for model file path."
      write( *, * ) "[-v], [--version] : Print version."
      write( *, * ) "[-h], [--help]      : Print usage information."
      !
   end subroutine printHelp

end program Test_FWD
