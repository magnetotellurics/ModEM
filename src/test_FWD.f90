program Test_FWD
  !
  use DataManager
  use ModEMControlFile
  use ModelReader
  use ModelReader_Weerachai
  use ModelOperator_MF
  !
  class( ModEMControlFile_t ), pointer  :: control_file
  !
  class( DataManager_t ), pointer       :: data_manager 
  class( Transmitter_t ), pointer       :: tx
  class( Receiver_t ), pointer          :: rx
  !
  class( Grid_t ), pointer              :: grid
  class( ModelParameter_t ), pointer    :: model_parameter
  class( ModelOperator_t ), allocatable :: model_operator
  type( ModelReader_Weerachai_t )       :: model_reader
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
    integer :: iTx, nTx, iRx, nRx
    real( kind=prec ) :: omega
    character(:), allocatable :: transmitter_type
    !
    transmitter_type = "Unknow"
    !
    if( .not. has_data_file ) stop " - Provide data file!"
    !
    if( .not. has_model_file ) stop " - Provide model file!"
    !
    write ( *, * ) "  > Start forward modelling."
    !
    ! Load model, define grid, model_operator, model_parameter
    call handleModelFile()
    !
    call model_operator%SetEquations()
    !
    call model_operator%SetCond( model_parameter )
    !
    ! Load data, construct the relation between txs and rxs
    call handleDataFile()
    !
    !Loop over all Transmitters
    nTx = transmitters%size()
    !
    do iTx = 1, nTx
      !
      tx => transmitters%get( iTx )
      !
      write( *, * ) "               Id:", tx%id, "Period:", int( tx%period )
      !
      call writePredictedDataHeader( tx, transmitter_type )
      !
      call tx%solveFWD( model_operator )
      !
      ! Loop over receivers of each transmitter
      nRx = tx%getNRx()
      !
      do iRx = 1, nRx
         !
         rx => receivers%get( tx%get( iRx ) )
         !
         call rx%predictedData( model_operator, tx )
         !
      enddo
      !
    enddo
    !
    deallocate( transmitters )
    deallocate( receivers )
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
      write( *, * ) "  - Unknow job: [", trim( modem_job ), "]"
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
    write( *, * ) "  -> Control File: [", control_file_name, "]"
    !
    ! Read control file
    control_file => ModEMControlFile_t( ioStartup, control_file_name ) 
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
    write( *, * ) "  -> Data File: [", data_file_name, "]"
    !
    ! Read data file and creates RXs <-> TXs
    data_manager => DataManager_t( data_file_name ) 
    !
    ! Dispose memory data that will no longer be used.
    deallocate( data_manager )
    !
  end subroutine handleDataFile
  !
  subroutine handleModelFile()
    implicit none
    !
    write( *, * ) "  -> Model File: [", model_file_name, "]"
    !
    call model_reader%Read( model_file_name, grid, model_parameter ) 
    !
    select type( grid )
    class is( Grid3D_SG_t )
       allocate( model_operator, source = ModelOperator_MF_t( grid ) )
    class default
       stop "Unclassified grid"
    end select
    !
    write( *, * ) "  -> Readed Model File: [", model_file_name, "]"
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
      write( *, * ) "  + ModEM-OO version 1.0.0"
      stop
      !
    case ( "-h", "--help" )
      !
      call printHelp()
      stop
      !
    case default
      !
      write( *, * ) "  - Unknow Argument: [", trim( argument ), "]"
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
   class( Transmitter_t ), intent( in )       :: tx
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
    write( *, * ) "[-d], [--data]    : Preceding argument for data file path."
    write( *, * ) "[-m], [--model]   : Preceding argument for model file path."
    write( *, * ) "[-v], [--version] : Print version."
    write( *, * ) "[-h], [--help]    : Print usage information."
    !
  end subroutine printHelp

end program Test_FWD
