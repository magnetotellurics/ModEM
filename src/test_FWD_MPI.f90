program ModEM
   !
   use DeclarationMPI
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
   class( Grid_t ), allocatable           :: main_grid
   class( ModelParameter_t ), allocatable :: model_parameter
   class( ModelOperator_t ), allocatable, target  :: model_operator
   type( ModelOperator_MF_t ), pointer   :: model_operator_ptr => null()
   type( DataManager_t ) :: data_manager
   !
   integer, pointer, dimension(:), save :: test_mpi_ptr
   !
   !
   character(:), allocatable :: control_file_name, model_file_name, data_file_name, modem_job
   logical                   :: has_control_file = .false., has_model_file = .false., has_data_file = .false.
   !
   main_comm = MPI_COMM_WORLD
   !
   call MPI_Init( ierr )
   !
   ! SET mpi_size WITH n FROM MPIRUN FOR MPI_COMM_WORLD
   call MPI_Comm_size( main_comm, mpi_size, ierr )
   !
   ! SET mpi_rank WITH PROCESS ID FOR MPI_COMM_WORLD
   call MPI_Comm_rank( main_comm, mpi_rank, ierr )
   !
   ! SPLIT MPI_COMM_WORLD into shared subcommunicator: shared_comm
   call MPI_Comm_split_type( main_comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shared_comm, ierr )
   !
   call MPI_Get_processor_name( nodename, nodestringlen, ierr )
   !
   call MPI_Comm_size( shared_comm, node_size, ierr )
   !
   call MPI_Comm_rank( shared_comm, node_rank, ierr )
   !
   write( *, * ) "Rank ", mpi_rank," in COMM_WORLD is ", node_rank, &
             " in SHARED_COMM on Node: ", nodename(1:nodestringlen)
   !
   write( *, * ) "node_size: ", node_size
   !
   !
   ! MASTER
   !
   if ( mpi_rank == 0 ) then
       !
       modem_job = "unknow"
       !
       ! MASTER
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
   ! WORKER
   !
   else
      !
      call MPI_Win_allocate_shared( winsize, 0, MPI_INFO_NULL, shared_comm, baseptr, nodewin, ierr )
      !
      !
      if( ierr == MPI_SUCCESS ) then
         !
         write( *, * ) "!!!! WORKER MPI_Win_allocate_shared SUCCEDED: ", baseptr
         !
         do while ( job_master .ne. job_finish )
            !
            write ( *, * ) "WORKER: ", mpi_rank, " WAITING MASTER"
            !
            call receiveFrom( master_id )
            !
            write ( *, * ) "WORKER: ", mpi_rank, " RECEIVE JOB: ", job_info%name
            !
            job_master = trim( job_info%name )
            !
            select case ( job_master )
                !
                case ( "FORWARD" )
                    !
                    call workerForwardModelling()
                    !
                case ( "STOP_JOBS" )
                    !
                    job_info%name = job_ok
                    !
                    write ( *, * ) "WORKER: ", mpi_rank, " SEND JOB: ", job_info%name, " TO MASTER"
                    !
                    call sendTo( master_id )
                !
            end select
            !
         enddo
         !
         call MPI_Win_fence( 0, nodewin, ierr )
         !
         call MPI_FINALIZE( ierr )
         !
      else
         stop "WORKER MPI_Win_allocate_shared FAILS !!!"
      endif
      !
    endif
   !
contains
   !
   subroutine masterForwardModelling()
      implicit none
      !
      ! These objects are frequency dependent,
      ! must be instantiated on each Worker
      class( ForwardSolver_t ), allocatable, target, save :: fwd_solver
      !
      class( Source_t ), allocatable, target, save        :: fwd_source 
      !
      ! Temporary alias pointers
      class( Transmitter_t ), pointer  :: Tx
      class( Receiver_t ), allocatable :: Rx
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
           fwd_solver = ForwardSolverIT_t( model_operator, QMR )
           !
        case( FWD_IT_DC )
           fwd_solver = ForwardSolverIT_DC_t( model_operator, QMR )
           !
        case default
           fwd_solver = ForwardSolverIT_DC_t( model_operator, QMR )
           !
      end select
      !
      call fwd_solver%setCond( model_parameter )
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
      job_info%name = "FORWARD"
      !
      ! Loop over all Transmitters
      nTx = size( transmitters )
      !
      call writeEsolutionHeader( nTx, 2 )
      !
      do iTx = 1, nTx
         !
         ! Temporary Transmitter alias
         Tx => getTransmitter( iTx )
         !
         write( *, * ) "FOR PERIOD: ", Tx%period
         !
         ! According to Tx type,
         ! write the proper header in the "predicted_data.dat" file
         !call writePredictedDataHeader( Tx, transmitter_type )
         !
         ! Tx points to its due Source
         call Tx%setSource( fwd_source )
         !
         ! Set ForwardSolver Period
         call fwd_solver%setPeriod( Tx%period )
         !
         ! Tx points to its due ForwardSolver
         call Tx%setForwardSolver( fwd_solver )
         !
         write ( *, * ) "MASTER SEND JOB: ", job_info%name, " TO WORKER: ", iTx
         !
         ! SEND FORWARD JOB TO WORKER
         call sendTo( iTx )
         !
      enddo
      !
      ! RECEIVE FOREACH WORKER
      do iTx = 1, nTx
          !
          write ( *, * ) "MASTER WAITING WORKER: ", iTx
          !
          call receiveFrom( iTx )
          !
          write ( *, * ) "MASTER RECEIVED JOB: ", job_info%name, " FROM WORKER: ", iTx
          !
      enddo
      !
      ! DISTRIBUTE STOP_JOBS
      job_info%name = "STOP_JOBS"
      !
      ! SEND FOREACH WORK PROCESS
      do iTx = 1, nTx
          !
          write ( *, * ) "MASTER SEND JOB: ", job_info%name, " TO WORKER: ", iTx
          !
          call sendTo( iTx )
          !
      enddo
      !
      ! RECEIVE FOREACH WORK PROCESS
      do iTx = 1, nTx
          !
          write ( *, * ) "MASTER WAITING WORKER: ", iTx
          !
          call receiveFrom( iTx )
          !
          write ( *, * ) "MASTER RECEIVED JOB: ", job_info%name, " FROM WORKER: ", iTx
          !
      enddo
      !
      !
      call MPI_Win_fence( 0, nodewin, ierr )
      !
      do iTx = 1, 5
          !
          write ( *, * ) "MASTER test_mpi_ptr(iTx): ", test_mpi_ptr(iTx)
          !
      enddo
      !
      !call MPI_Win_fence( 0, nodewin, ierr )
      !
      call MPI_FINALIZE( ierr )
      !
   end subroutine masterForwardModelling
   !
   subroutine workerForwardModelling()
      implicit none
      !
      ! Temporary alias pointers
      class( Transmitter_t ), pointer  :: Tx
      class( Receiver_t ), allocatable :: Rx
      !
      ! Local variables
      integer :: iRx, nRx
      !
      !call sleep(  mpi_rank * 5 )
      !
      write( *, * ) "### WORKER: ", mpi_rank, " START JOB:", job_info%name
      !
      call MPI_Win_shared_query( nodewin, 0, winsize, disp_unit, baseptr, ierr )
      !
      if( ierr == MPI_SUCCESS ) then
         write( *, * ) "!!!! WORKER MPI_Win_shared_query SUCCEDED: ", winsize, baseptr
      else
         write( *, * ) "!!!! WORKER MPI_Win_shared_query FAILS: ", winsize, baseptr
      endif
      !
      !write( *, * ) "test_mpi_ptr:[", test_mpi_ptrl, "]"
      !
      call c_f_pointer( baseptr, test_mpi_ptr, (/winsize/) )
	  !call c_f_pointer( baseptr, model_operator_ptrl )
      !
      !write( *, * ) "model_operator_ptrl:[", model_operator_ptrl%grid%nx, model_operator_ptrl%grid%ny, model_operator_ptrl%grid%nz, "]"
      !
      test_mpi_ptr( node_rank+1 ) = node_rank * 10
      !
      ! Temporary Transmitter alias
      !Tx => getTransmitter( mpi_rank )
      !
      ! Verbosis...
      !write( *, * ) "WORKER Tx Id:", Tx%id, "Period:", int( Tx%period )
      !
      ! Solve Tx Forward Modelling
      !call Tx%solveFWD()
      !
      ! Loop over Receivers of each Transmitter
      !nRx = Tx%getNRx()
      !
      !do iRx = 1, nRx
         !
         ! Temporary Receiver alias
         !Rx = receivers%get( Tx%get( iRx ) )
         !
         ! Verbosis...
         !write( *, * ) "                Rx Id:", Rx%id, "XYZ:", Rx%location
         !
         ! Calculate Rx Predicted Data
         !call Rx%predictedData( model_operator, Tx )
         !
      !enddo
      !
      !deallocate( Tx )
      !
      job_info%name = job_ok
      !
      write ( *, * ) "WORKER: ", mpi_rank, " SEND JOB: ", job_info%name, " TO MASTER"
      !
      call sendTo( master_id )
      !
   end subroutine workerForwardModelling
   !
   subroutine handleJob()
      implicit none
      !
      select case ( modem_job )
        !
      case ( "forward" )
        !
        call masterForwardModelling()
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
      !type( DataManager_t ) :: data_manager
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
      character(:), allocatable      :: fname
      !
      type( TAirLayers )            :: air_layer
      !
      fname = "/mnt/c/Users/protew/Desktop/ON/GITLAB_PROJECTS/modem-oo/inputs/Full_A_Matrix_TinyModel"
      !fname = "/Users/garyegbert/Desktop/ModEM_ON/modem-oo/inputs/Full_A_Matrix_TinyModel"
      !
      write( *, * ) "   -> Model File: [", model_file_name, "]"
     !
     model_method = MM_METHOD_FIXED_H
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
           !   as coded have to use air_layer data structure to update grid
           call main_grid%UpdateAirLayers( air_layer%nz, air_layer%dz )
           !
           !model_operator = ModelOperator_File_t( main_grid, fname )
           !
           model_operator = ModelOperator_MF_t( main_grid )
           !
           call model_parameter%setMetric( model_operator%metric )
           !
           ! complete model operator setup
           call model_operator%SetEquations()
           !
           call model_operator%SetCond( model_parameter )
           !
           select type( model_operator )
              !
              class is( ModelOperator_MF_t )
              !
              allocate( test_mpi_ptr(5) )
              !
              !winsize = size( test_mpi_ptr )
			  winsize = sizeof( model_operator )
              !
              call MPI_Win_allocate_shared( winsize, disp_unit, MPI_INFO_NULL, shared_comm, baseptr, nodewin, ierr )
              !
              if( ierr == MPI_SUCCESS ) then
                 write( *, * ) "!!!! MASTER MPI_Win_allocate_shared SUCCEDED: ", winsize, baseptr
              else
                 write( *, * ) "!!!! MASTER MPI_Win_allocate_shared FAILS: ", winsize, baseptr
              endif
              !
              call c_f_pointer( baseptr, test_mpi_ptr, (/winsize/) )
			  !call c_f_pointer( baseptr, model_operator_ptr )
              !
			  test_mpi_ptr = (/-1, -1, -1, -1, -1/)
			  !allocate( model_operator_ptr, source = model_operator )
			  !
              !write( *, * ) "### MASTER model_operator_ptr:[", model_operator_ptr%grid%nx, model_operator_ptr%grid%ny, model_operator_ptr%grid%nz, "]"
              !
           end select
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
      integer       :: argument_index
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
   subroutine writeEsolutionHeader( nTx, nMode )
      implicit none
      !
      ! implement separated routine
      integer, intent( in ) :: nTx, nMode
      integer             :: ios
      character (len=20)    :: version
      !
      open( ioESolution, file = "e_solution", action="write", form ="unformatted", iostat=ios)
      !
      if( ios/=0) then
        write(0,*) "Error opening file in FileWriteInit: e_solution"
      else
        !
        version = ""
        ! write the header (contains the basic information for the forward
        ! modeling). the header is 4 lines
        write( ioESolution ) version, nTx, nMode, &
        main_grid%nx, main_grid%ny, main_grid%nz, main_grid%nzAir, &
        main_grid%ox, main_grid%oy, main_grid%oz, main_grid%rotdeg
        !
        write( ioESolution ) main_grid%dx
        write( ioESolution ) main_grid%dy
        write( ioESolution ) main_grid%dz
      endif
      !
      !
   end subroutine writeEsolutionHeader
   !
   subroutine writePredictedDataHeader( Tx, transmitter_type )
      implicit none
      !
      class( Transmitter_t ), intent( in )      :: Tx
      character(:), allocatable, intent( inout ) :: transmitter_type
      !
      integer :: nTx, nRx
      logical :: tx_changed = .false.
      !
      if( ( index( transmitter_type, "Unknow" ) /= 0 ) .OR. transmitter_type /= trim( Tx%type ) ) then
        !
        tx_changed = .true.
        !
      endif
      !
      if( ( index( transmitter_type, "Unknow" ) /= 0 ) ) then
        !
        open( ioPredData, file = "predicted_data.dat", action="write", form ="formatted" )
        !
      else if( transmitter_type /= trim( Tx%type ) ) then
        !
        open( ioPredData, file = "predicted_data.dat", action="write", form ="formatted", position="append" )
        !
      endif
      !
      if( tx_changed ) then
        !
        write( ioPredData, "(4A, 100A)" ) "#   ", DATA_FILE_TITLE
        write( ioPredData, "(4A, 100A)" ) "#   ", Tx%DATA_TITLE
        write( ioPredData, "(4A, 100A)" ) ">   ", trim( Tx%type )
        write( ioPredData, "(4A, 100A)" ) ">   ", "exp(-i\omega t)"
        write( ioPredData, "(4A, 100A)" ) ">   ", "[V/m]/[T]"
        write( ioPredData, "(7A, 100A)" ) ">      ", "0.00"
        write( ioPredData, "(7A, 100A)" ) ">      ", "0.000   0.000"
        write( ioPredData, "(A3, i8, i8)" ) ">      ", size( transmitters ), receivers%size()
        !
        close( ioPredData )
        !
        transmitter_type = trim( Tx%type )
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
