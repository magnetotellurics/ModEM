!
!> Module with the Forward Modeling routines
!
module ForwardModeling
!
#ifdef MPI
    use MasterMPI
#else
    use CoreComponents
#endif
    !
    !> Public module routines
    public :: solveAll, solveTx
    public :: jobForwardModeling, serialForwardModeling
    public :: createDistributeForwardSolver
    public :: writeAllESolutionHeader
    !
contains
    !
    !> Calculate ESolution for all transmitters.
    !> ForwardSolver must be allocated!
    !
    subroutine solveAll( sigma )
        implicit none
        !
        class( ModelParameter_t ), intent( inout ) :: sigma
        !
        class( Transmitter_t ), pointer :: Tx
        integer :: i_tx
        !
        ! Verbose
        !write( *, * ) "          - Start EM Solve All"
        !
        !> Loop over all Transmitters
        do i_tx = 1, size( transmitters )
            !
            !> Point to the current Transmitter
            Tx => getTransmitter( i_tx )
            !
            call solveTx( sigma, Tx )
            !
        enddo
        !
        ! Verbose
        !write( *, * ) "          - Finish EM Solve All"
        !
    end subroutine solveAll
    !
    !> Calculate E_Solution(e_sol_0) for a single Transmitter
    !> Tx%forward_solver must be previously allocated !!!!
    !
    subroutine solveTx( sigma, Tx )
        implicit none
        !
        class( ModelParameter_t ), intent( inout ) :: sigma
        class( Transmitter_t ), pointer, intent( inout ) :: Tx
        !
        class( Source_t ), allocatable :: source
        !
        call Tx%forward_solver%setFrequency( sigma, Tx%period )
        !
        !> Instantiate Transmitter's Source - According to transmitter type and chosen via fwd control file
        select type( Tx )
            !
            class is( TransmitterMT_t )
                !
                !> Instantiate the MT Source - Specific type can be chosen via fwd control file
                select case( source_type_mt )
                    !
                    case( SRC_MT_1D )
                        !
                        allocate( source, source = SourceMT_1D_t( model_operator, sigma, Tx%period ) )
                        !
                    case( SRC_MT_2D )
                        !
                        !allocate( source, source = SourceMT_2D_t( model_operator, sigma, Tx%period ) )
                        !
                        call errStop( "solveTx > SourceMT_2D not implemented yet!" )
                        !
                    case( "" )
                        !
                        call warning( "solveTx > MT Source type not provided, using SourceMT_1D_t." )
                        !
                        allocate( source, source = SourceMT_1D_t( model_operator, sigma, Tx%period ) )
                        !
                    case default
                        !
                        call errStop( "solveTx > Wrong MT Source type: ["//source_type_mt//"]" )
                        !
                end select
                !
            class is( TransmitterCSEM_t )
                !
                !> Instantiate the CSEM Source - Specific type can be chosen via control file
                select case( source_type_csem )
                    !
                    case( SRC_CSEM_EM1D )
                        !
                        allocate( source, source = SourceCSEM_EM1D_t( model_operator, sigma, Tx%period, Tx%location, Tx%i_tx ) )
                        !
                    case( SRC_CSEM_DIPOLE1D )
                        !
                        allocate( source, source = SourceCSEM_Dipole1D_t( model_operator, sigma, Tx%period, Tx%location, Tx%dip, Tx%azimuth, Tx%moment ) )
                        !
                    case( "" )
                        !
                        call warning( "solveTx > CSEM Source type not provided, using EM1D." )
                        !
                        allocate( source, source = SourceCSEM_EM1D_t( model_operator, sigma, Tx%period, Tx%location, Tx%i_tx ) )
                        !
                    case default
                        !
                        call errStop( "solveTx > Wrong CSEM Source type: ["//source_type_csem//"]" )
                        !
                end select
                !
            class default
                call errStop( "solveTx > Unclassified Transmitter" )
            !
        end select
        !
        call Tx%setSource( source )
        !
        call Tx%source%createE
        !
        call Tx%solve
        !
        !>
        if( has_e_solution_file ) call Tx%writeESolution
        !
    end subroutine solveTx
    !
    !> Routine to execute a full ForwardModeling job 
    !> Outputting its result(PredictedData) into a text file
    !> Default <all_predicted_data.dat> or specified by argument -pd|--predicted
    !
    subroutine jobForwardModeling()
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: all_predicted_data
        !
        class( ModelParameter_t ), allocatable :: sigma
        !
        ! Verbose
        write( *, "( A30 )" ) "- Start Forward Modeling"
        write( *, * ) ""
        !
        if( has_model_file ) then
            !
            call handleModelFile( sigma )
        !
        else
            call errStop( "jobForwardModeling > Missing Model File!" )
        endif
        !
        if( has_data_file ) then
            !
            call handleDataFile
        !
        else
            call errStop( "jobForwardModeling > Missing Data File!" )
        endif
        !
        !> If path is specified by argument -es|-- esolution
        !> Write e_sol_0 to a binary file at this path
        if( has_e_solution_file ) call writeAllESolutionHeader( size( transmitters ), transmitters(1)%Tx%n_pol, e_solution_file_name )
        !
        all_predicted_data = all_measured_data
        !
#ifdef MPI
        !
        call broadcastBasicComponents
        !
        !> Deallocate global FWD Objects not used by the Master process
        deallocate( model_operator, main_grid )
        !
        call masterForwardModelling( sigma, all_predicted_data )
        !
        call broadcastFinish
        !
#else
        !
        call createDistributeForwardSolver
        !
        call serialForwardModeling( sigma, all_predicted_data )
        !
#endif
        !
        call writeData( all_predicted_data, predicted_data_file_name )
        !
        deallocate( sigma )
        !
        ! Verbose
        write( *, * ) ""
        write( *, "( A31 )" ) "- Finish Forward Modeling"
        !
    end subroutine jobForwardModeling
    !
    !> serialForwardModeling:
    !>     Calculate ESolution for all transmitters
    !>     Calculate the predicted data for each transmitter-receiver pair.
    !>     ForwardSolver must be allocated
    !
    subroutine serialForwardModeling( sigma, all_predicted_data, i_sol )
        implicit none
        !
        class( ModelParameter_t ), intent( inout ) :: sigma
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: all_predicted_data
        integer, intent( in ), optional :: i_sol
        !
        class( Transmitter_t ), pointer :: Tx
        class( Receiver_t ), pointer :: Rx
        type( DataGroup_t ) :: data_group
        integer :: i_tx, n_tx, i_rx, sol_index
        !
        sol_index = 0
        !
        !> Set i_sol if present
        if( present( i_sol ) ) sol_index = i_sol
        !
        !>
        n_tx = size( transmitters )
        !
        !> Loop over all Transmitters
        do i_tx = 1, n_tx
            !
            !> Pointer to the Transmitter
            Tx => getTransmitter( i_tx )
            !
            Tx%i_sol = sol_index
            !
            call solveTx( sigma, Tx )
            !
            ! Verbose
            write( *, "( A36 )" ) "- Calculate Predicted Data"
            !
            !> Loop for each Receiver related to this Transmitter
            do i_rx = 1, size( Tx%receiver_indexes )
                !
                !> Pointer to the Tx Receiver
                Rx => getReceiver( Tx%receiver_indexes( i_rx ) )
                !
                call Rx%predictedData( Tx, data_group )
                !
                call all_predicted_data( i_tx )%setValues( data_group )
                !
            enddo
            !
        enddo
        !
    end subroutine serialForwardModeling
    !
    !> Create the global ForwardSolver and point all transmitters to it.
    !
    subroutine createDistributeForwardSolver()
        implicit none
        !
        class( Transmitter_t ), pointer :: Tx
        integer :: i_tx
        !
        if( allocated( forward_solver ) ) deallocate( forward_solver )
        !
        !> Instantiate the ForwardSolver - Specific type can be set via control file
        select case( forward_solver_type )
            !
            case( FWD_IT )
                !
                allocate( forward_solver, source = ForwardSolver_IT_t( model_operator, solver_type ) )
                !
            case( FWD_IT_DC )
                !
                allocate( forward_solver, source = ForwardSolver_IT_DC_t( model_operator, solver_type ) )
                !
            case( "" )
                !
                !> Set the best forward_solver_type according to the model_operator_type.
                if( model_operator_type .EQ. MODELOP_MF .OR. model_operator_type .EQ. MODELOP_SP ) then
                    !
                    call warning( "Parameter forward_solver_type not provided, using IT_DC." )
                    !
                    allocate( forward_solver, source = ForwardSolver_IT_DC_t( model_operator, solver_type ) )
                    !
                elseif( model_operator_type .EQ. MODELOP_SP2 ) then
                    !
                    call warning( "Parameter forward_solver_type not provided, using IT." )
                    !
                    allocate( forward_solver, source = ForwardSolver_IT_t( model_operator, solver_type ) )
                    !
                endif
                !
            case default
                !
                call errStop( "createDistributeForwardSolver > Wrong forward_solver_type: ["//forward_solver_type//"]" )
                !
        end select
        !
        do i_tx = 1, size( transmitters )
            !
            Tx => getTransmitter( i_tx )
            !
            Tx%forward_solver => forward_solver
            !
        enddo
        !
    end subroutine createDistributeForwardSolver
    !
    !> Write a header into the ESolution binary file
    !
    subroutine writeAllESolutionHeader( n_tx, nMode, file_name )
        implicit none
        !
        integer, intent( in ) :: n_tx, nMode
        character(*), intent( in ) :: file_name
        integer :: ios
        character( len=20 ) :: version
        !
        version = "ModEM "//VERSION
        !
        open( ioESolution, file = file_name, action = "write", form = "unformatted", iostat = ios )
        !
        if( ios == 0 ) then
            !
            write( ioESolution ) version, n_tx, nMode, &
            main_grid%nx, main_grid%ny, main_grid%nz, main_grid%nzAir, &
            main_grid%ox, main_grid%oy, main_grid%oz, main_grid%rotdeg
            !
            write( ioESolution ) main_grid%dx
            write( ioESolution ) main_grid%dy
            write( ioESolution ) main_grid%dz
            !
            close( ioESolution )
            !
        else
            !
            call errStop( "writeAllESolutionHeader > Unable to open file ["//file_name//"]!" )
            !
        endif
        !
    end subroutine writeAllESolutionHeader
    !
end module ForwardModeling
!