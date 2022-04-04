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
   character(:), allocatable :: process_name
   !
   class( Grid_t ), allocatable           :: main_grid
   class( ModelParameter_t ), allocatable :: model_parameter
   class( ModelOperator_t ), allocatable, target  :: model_operator
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
   ! SPLIT MPI_COMM_WORLD into shared subcommunicator: child_comm
   call MPI_Comm_split_type( main_comm, MPI_COMM_TYPE_SHARED, mpi_rank, MPI_INFO_NULL, child_comm, ierr )
   !
   call MPI_Get_processor_name( node_name, nodestringlen, ierr )
   !
   call MPI_Comm_size( child_comm, node_size, ierr )
   !
   call MPI_Comm_rank( child_comm, node_rank, ierr )
   !
   write( *, * ) "Rank ", mpi_rank," in COMM_WORLD (", mpi_size, ") is ", node_rank, &
             " in SHARED_COMM (", node_size, ") on Node: ", node_name(1:nodestringlen)
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
      call MPI_Win_allocate_shared( shared_window_size, disp_unit, MPI_INFO_NULL, child_comm, shared_c_ptr, shared_window, ierr )
      !
      if( ierr == MPI_SUCCESS ) then
         !
         write( *, * ) "!!!! WORKER MPI_Win_allocate_shared SUCCEDED: ", shared_window_size, disp_unit, shared_c_ptr
         !
         do while ( job_master .ne. job_finish )
            !
            write ( *, * ) "WORKER: ", mpi_rank, " WAITING MASTER"
            !
            call receiveFrom( master_id )
            !
            job_master = trim( job_info%name )
            !
            select case ( job_master )
                !
                case ( "SHARE_MEM" )
                    !
                    call workerQuerySharedMemory()
                    !
                case ( "FORWARD" )
                    !
                    call workerForwardModelling()
                    !
                    ! SEND FORWARD JOB TO WORKER
                    job_info%name = job_done
                    job_info%worker_rank = mpi_rank
                    !
                    call sendTo( master_id )
                    !
            end select
            !
         enddo
         !
		 write ( *, * )  "WORKER: ", mpi_rank, " DONE FOR THE DAY!"
		 !
         call MPI_Win_fence( 0, shared_window, ierr )
         !
         call MPI_Finalize( ierr )
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
      ! Local variables
      integer :: pid, worker_rank = 1, tx_received = 0, tx_index = 0
      character(:), allocatable :: transmitter_type
      !
      ! Verbosis
      write ( *, * ) "   > Start forward modelling."
      !
      ! Tx type for predicted_data header changes
      transmitter_type = "Unknow"
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
      ! SHARE MEM WITH ALL WORKERS
      call masterExposeSharedMemory()
      !
      do pid = 1, mpi_size - 1
          !
          job_info%name = "SHARE_MEM"
          !
          call sendTo( pid )
          !
      enddo
      !
      ! SEND 1 TRANSMITTER TO EACH WORKER
      do while ( worker_rank <= ( mpi_size - 1 ) )
         !
         tx_index = tx_index + 1
         !
         job_info%name = "FORWARD"
         job_info%worker_rank = worker_rank
         job_info%tx_index = tx_index
         !
         call sendTo( worker_rank )
         !
         worker_rank = worker_rank + 1
         !
      end do
      !
      ! SEND 1 TRANSMITTER TO FIRST AVAILABLE WORKER
      do while( tx_index < size( transmitters ) )
         !
         write ( *, * ) "THERES", ( size( transmitters ) - tx_index ), " TX LEFT!"
         !
         call receiveFromAny()
         !
         tx_received = tx_received + 1
         !
         tx_index = tx_index + 1
         !
         job_info%name = "FORWARD"
         job_info%tx_index = tx_index
         !
         call sendTo( job_info%worker_rank )
         !
      end do
      !
      ! RECEIVES job_done FROM EACH FINISHED WORKER
      do while ( tx_received < size( transmitters ) )
          !
          write ( *, * ) "MASTER WAITING ANY WORKER"
          !
          call receiveFromAny()
          !
          tx_received = tx_received + 1
          !
          job_info%name = job_finish
          !
          call sendTo( job_info%worker_rank )
          !
      enddo
      !
      call MPI_Win_fence( 0, shared_window, ierr )
      !
      process_name = "#### FINAL ####"
      call showProcessState( process_name )
      !
      ! Verbosis
      write ( *, * ) "   > Finish forward modelling."
      !
   end subroutine masterForwardModelling
   !
   subroutine showProcessState( process_name )
      implicit none
      !
      character(:), allocatable :: process_name
      !
      !
      write( *, * ) trim( process_name )
      !
      write( *, * ) "#### GRID:[", main_grid%allocated, main_grid%nx, main_grid%ny, main_grid%nz, "]"
      !
      write( *, * ) "#### MODEL_OPERATOR:"
      select type( model_operator )
        !
        class is( ModelOperator_MF_t )
            !
            call model_operator%Sigma_E%print()
            write( *, * ) "#### Sigma_E GRID TYPE:[", model_operator%Sigma_E%gridType, "]"
            call model_operator%c%print()
            !
      end select
      !
      write( *, * ) "#### MODEL PARAMETER:"
      select type( model_parameter )
        !
        class is( ModelParameterCell_SG_t )
            !
            write( *, * ) "#### PARAM GRID:[", model_parameter%paramGrid%nx, model_parameter%paramGrid%ny, model_parameter%paramGrid%nz, "]"
            !call model_parameter%cellCond%print()
            !
      end select
      !
   end subroutine showProcessState
   !
   !
   subroutine workerQuerySharedMemory()
      implicit none
      !
      call MPI_Win_shared_query( shared_window, master_id, shared_window_size, disp_unit, shared_c_ptr, ierr )
      !
      if( ierr == MPI_SUCCESS ) then
         write( *, * ) "!!!! WORKER [", mpi_rank, "] MPI_Win_shared_query SUCCEDED: ", shared_window_size, disp_unit, shared_c_ptr
      else
         write( *, * ) "!!!! WORKER MPI_Win_shared_query FAILS: ", shared_window_size, disp_unit, shared_c_ptr
         stop
      endif
      !
      call c_f_pointer( shared_c_ptr, shared_buffer, (/shared_window_size/) )
      !
      call unpackSharedBuffer( int( shared_window_size ), main_grid, model_operator, model_parameter )
      !
      process_name = "#### WORKER ####"
      call showProcessState( process_name )
      !
   end subroutine workerQuerySharedMemory
   !
   subroutine workerForwardModelling()
      implicit none
      !
      class( ForwardSolver_t ), allocatable :: fwd_solver
      !
      class( Source_t ), allocatable        :: fwd_source 
      !
      ! Temporary alias pointers
      class( Receiver_t ), allocatable     :: Rx
      class( Transmitter_t ), allocatable  :: Tx
      !
      ! Local variablesh
      integer :: iTx, iRx, nRx
      !
      ! High-level object instantiation
      ! Some types are chosen from the control file
      !
      write( *, * ) "INIT WORWER FWD"
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
      write( *, * ) "WORWER FWD SET SOLVER OK"
      !
      call fwd_solver%setCond( model_parameter )
      !
      write( *, * ) "WORWER FWD SET SOLVER COND OK"
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
      write( *, * ) "WORWER FWD SET SOURCE COND OK"
      !
      ! Forward Modelling
      !
      !call writeEsolutionHeader( size( transmitters ), 2 )
      !
      !
      ! Temporary Transmitter alias
      !Tx = getTransmitter( iTx )
      !
      ! Verbosis...
      !write( *, * ) "   Tx Id:", Tx%id, "Period:", int( Tx%period )
      !
      ! According to Tx type,
      ! write the proper header in the "predicted_data.dat" file
      !call writePredictedDataHeader( Tx, transmitter_type )
      !
      ! Tx points to its due Source
      !call Tx%setSource( fwd_source )
      !
      ! Set ForwardSolver Period
      !call fwd_solver%setPeriod( Tx%period )
      !
      ! Tx points to its due ForwardSolver
      !call Tx%setForwardSolver( fwd_solver )
      !
      ! Solve Tx Forward Modelling
      !call Tx%solveFWD()
      !
      ! Loop over Receivers of each Transmitter
      !nRx = size( Tx%receiver_indexes )
      !
      !do iRx = 1, nRx
      !
      ! Temporary Receiver alias
      !Rx = getReceiver( Tx%receiver_indexes( iRx ) )
      !
      ! Verbosis...
      !write( *, * ) "                  Rx Id:", Rx%id, "XYZ:", Rx%location
      !
      ! Calculate Rx Predicted Data
      !call Rx%predictedData( model_operator, Tx )
      !
      !enddo
      !
      !deallocate( Tx )
      !
      ! Loop over all Receivers
      !nRx = size( receivers )
      !
      !do iRx = 1, nRx
         !
         ! Temporary Receiver alias
         !Rx = getReceiver( iRx )
         !
         !call Rx%writePredictedData()
         !
         !deallocate( Rx )
         !
      !enddo
      !
      !deallocate( data_groups )
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
   subroutine handleDataFile()
      implicit none
      !
      type( DataManager_t ) :: data_manager
      !
      write( *, * ) "   -> Data File: [", data_file_name, "]"
      !
      data_manager = DataManager_t( data_file_name )
      !
      if( ( mpi_size - 1 ) > size( transmitters ) ) then
         write( *, * ) "There are more MPI worker processes than necessary!!!"
         write( *, * ) "    ", size( transmitters ), " Transmitters"
         write( *, * ) "    ", ( mpi_size - 1 ), " MPI processes"
         !
         call MPI_Abort( main_comm, ierr )
         !
         stop
      endif
      !
   end subroutine handleDataFile
   !
   subroutine handleModelFile()
      implicit none
      !
      ! It remains to standardize ????
      type( ModelReader_Weerachai_t ) :: model_reader
      !
      type( TAirLayers )            :: air_layer
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
           !
           call main_grid%UpdateAirLayers( air_layer%nz, air_layer%dz )
           !
           model_operator = ModelOperator_MF_t( main_grid )
           !
           call model_parameter%setMetric( model_operator%metric )
           !
           ! complete model operator setup
           call model_operator%SetEquations()
           !
           class default
              stop "Unclassified main_grid"
           !
      end select
      !
   end subroutine handleModelFile
   !
   subroutine masterExposeSharedMemory()
      implicit none
      !
      select type( model_operator )
        !
        class is( ModelOperator_MF_t )
            !
            call allocateSharedBuffer( main_grid, model_operator, model_parameter )
            !
            shared_window_size = shared_buffer_size
            disp_unit = 1
            !
            call MPI_Win_allocate_shared( shared_window_size, disp_unit, MPI_INFO_NULL, child_comm, shared_c_ptr, shared_window, ierr )
            !
            if( ierr == MPI_SUCCESS ) then
                write( *, * ) "!!!! MASTER MPI_Win_allocate_shared SUCCEDED: ", shared_window_size, disp_unit, shared_c_ptr
            else
                write( *, * ) "!!!! MASTER MPI_Win_allocate_shared FAILS: ", shared_window_size, disp_unit, shared_c_ptr
            endif
            !
            call c_f_pointer( shared_c_ptr, shared_buffer, (/shared_window_size/) )
            !
            call packSharedBuffer( main_grid, model_operator, model_parameter )
            !
            process_name = "#### MASTER ####"
            call showProcessState( process_name )
            !
      end select
      !
   end subroutine masterExposeSharedMemory
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
   subroutine writeEsolutionHeader( nMode )
      implicit none
      !
      ! implement separated routine
      integer, intent( in ) :: nMode
      integer               :: ios
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
        write( ioESolution ) version, size( transmitters ), nMode, &
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
      class( Transmitter_t ), intent( in )       :: Tx
      character(:), allocatable, intent( inout ) :: transmitter_type
      !
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
        write( ioPredData, "(A3, i8, i8)" ) ">      ", size( transmitters ), size( receivers )
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
