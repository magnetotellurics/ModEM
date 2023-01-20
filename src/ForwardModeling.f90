!
!> Module with the ForwardModeling routines
!
module ForwardModeling
!
#ifdef MPI
    use MasterMPI
#else
    use GlobalVariables
#endif
    !
    !> Public module routines
    public :: solveAll, solveTx
    public :: jobForwardModeling, serialForwardModeling
    public :: createDistributeForwardSolver
    public :: writeAllESolution, writeAllESolutionHeader
    !
contains
    !
    !> Calculate ESolution for all transmitters.
    !> ForwardSolver must be allocated
    subroutine solveAll( sigma )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma
        !
        class( Transmitter_t ), pointer :: Tx
        integer :: i_tx, n_tx
        !
        ! Verbose
        write( *, * ) "          - Start EM Solve"
        !
        n_tx = size( transmitters )
        !
        !> Loop over all Transmitters
        do i_tx = 1, n_tx
            !
            !> Point to the current Transmitter
            Tx => getTransmitter( i_tx )
            !
            call solveTx( sigma, Tx )
            !
        enddo
        !
        ! Verbose
        write( *, * ) "          - Finish EM Solve"
        !
    end subroutine solveAll
    !
    !> Calculate E_Solution (e_sol) for a single Transmitter
    !> ForwardSolver must be allocated
    !
    subroutine solveTx( sigma, Tx )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma
        class( Transmitter_t ), pointer, intent( inout ) :: Tx
        !
        call Tx%forward_solver%setFrequency( sigma, Tx%period )
        !
        !> Instantiate Transmitter's Source - According to transmitter type and chosen via control file
        select type( Tx )
            !
            class is( TransmitterMT_t )
                !
                call Tx%setSource( SourceMT_1D_t( model_operator, sigma, Tx%period ) )
                !
            class is( TransmitterCSEM_t )
                !
                call Tx%setSource( SourceCSEM_Dipole1D_t( model_operator, sigma, Tx%period, Tx%location, Tx%dip, Tx%azimuth, Tx%moment ) )
                !
            class default
                stop "Error: solveTx > Unclassified Transmitter"
            !
        end select
        !
        call Tx%source%createE()
        !
        call Tx%solve()
        !
    end subroutine solveTx
    !
    !> Routine to run a full ForwardModeling job 
    !> and deliver the result (PredictedData) in a text file <all_predicted_data.dat>
    !
    subroutine jobForwardModeling()
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: all_predicted_data
        !
        class( ModelParameter_t ), allocatable :: sigma
        !
        ! Verbose
        write( *, * ) "     - Start jobForwardModeling"
        !
        if( has_model_file ) then
            !
            call handleModelFile( sigma )
            !
        else
            stop "Error: jobForwardModeling > Missing Model file!"
        endif
        !
        if( has_data_file ) then
            !
            call handleDataFile()
            !
        else
            stop "Error: jobForwardModeling > Missing Data file!"
        endif
        !
        all_predicted_data = all_measured_data
        !
#ifdef MPI
        !
        call broadcastBasicComponents()
        !
        !> Deallocate global FWD components not used by the Master process
        deallocate( model_operator, main_grid )
        !
        call masterForwardModelling( sigma, all_predicted_data )
        !
        call broadcastFinish
        !
#else
        !
        call createDistributeForwardSolver()
        !
        call serialForwardModeling( sigma, all_predicted_data )
        !
        if( has_e_solution_file ) call writeAllESolution( e_solution_file_name )
        !
#endif
        !
        call writeDataGroupTxArray( all_predicted_data, predicted_data_file_name )
        !
        deallocate( sigma )
        !
        ! Verbose
        write( *, * ) "     - Finish jobForwardModeling"
        !
    end subroutine jobForwardModeling
    !
    !> serialForwardModeling:
    !>     Calculate ESolution for all transmitters
    !>     Calculate the predicted data for each transmitter-receiver pair.
    !>     ForwardSolver must be allocated
    !
    subroutine serialForwardModeling( sigma, all_predicted_data, e_all )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: all_predicted_data
        type( EAllMTx_t ), optional, intent( inout ) :: e_all
        !
        class( Transmitter_t ), pointer :: Tx
        class( Receiver_t ), pointer :: Rx
        type( DataGroup_t ) :: data_group
        integer :: i_tx, n_tx, i_rx
        !
        ! Verbose
        write( *, * ) "          - Start Forward Modeling"
        !
        !>
        n_tx = size( transmitters )
        !
        !> Set e_all if present
        if( present( e_all ) ) then
            !
            !> SET e_all
            if( allocated( e_all%e ) ) deallocate( e_all%e )
            allocate( e_all%e( n_tx ) )
            !
        endif
        !
        !> Loop over all Transmitters
        do i_tx = 1, n_tx
            !
            !> Pointer to the Transmitter
            Tx => getTransmitter( i_tx )
            !
            call solveTx( sigma, Tx )
            !
            if( present( e_all ) ) then
                !
                e_all%e( i_tx )%pol = Tx%e_sol
                !
            endif
            !
            ! Verbose
            write( *, * ) "               - Calculate Predicted Data"
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
        ! Verbose
        write( *, * ) "          - Finish Forward Modeling"
        !
    end subroutine serialForwardModeling
    !
    !> Create the global ForwardSolver and point all transmitters to it.
    !
    subroutine createDistributeForwardSolver()
        implicit none
        !
        integer :: i_tx
        !
        if( allocated( forward_solver ) ) deallocate( forward_solver )
        !
        !> Instantiate the ForwardSolver - Specific type can be chosen via control file
        select case ( forward_solver_type )
            !
            case( FWD_IT_DC )
                allocate( forward_solver, source = ForwardSolverIT_DC_t( model_operator, QMR ) )
                !
            case default
                !
                write( *, * ) "Warning: createDistributeForwardSolver > Undefined forward_solver, using IT_DC"
                !
                allocate( forward_solver, source = ForwardSolverIT_DC_t( model_operator, QMR ) )
                !
        end select
        !
        do i_tx = 1, size( transmitters )
            !
            transmitters( i_tx )%Tx%forward_solver => forward_solver
            !
        enddo
        !
    end subroutine createDistributeForwardSolver
    !
    !> Write all e_solutions, with its proper Rx headers, 
    !> into to the binary file <predicted_data_file_name>
    subroutine writeAllESolution( file_name )
        implicit none
        !
        character(*), intent( in ) :: file_name
        !
        integer :: i_tx, n_tx
        !
        if( .NOT. allocated( transmitters ) .OR. size( transmitters ) .LT. 1 ) then
            stop "Error: writeAllESolution > Theres no transmitters to write!"
        endif
        !
        if( .NOT. allocated( transmitters(1)%Tx%e_sol ) .OR. size( transmitters(1)%Tx%e_sol ) .LT. 1 ) then
            stop "Error: writeAllESolution > Theres no e_sol to write!"
        endif
        !
        ! Verbose
        write( *, * ) "          > Write All ESolutions to file: [", file_name, "]"
        !
        !>
        n_tx = size( transmitters )
        !
        call writeAllESolutionHeader( n_tx, transmitters(1)%Tx%n_pol, file_name )
        !
        do i_tx = 1, n_tx
            !
            call transmitters( i_tx )%Tx%writeESolution()
            !
        enddo
        !
    end subroutine writeAllESolution
    !
    !> Write a header into the ESolution binary file
    !
    subroutine writeAllESolutionHeader( n_tx, nMode, file_name )
        implicit none
        !
        integer, intent( in ) :: n_tx, nMode
        character(*), intent( in ) :: file_name
        integer :: ios
        character(len=20) :: version
        !
        version = "Modem-OO"
        !
        open( ioESolution, file = file_name, action = "write", form = "unformatted", iostat = ios)
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
            write( *, * ) "Error opening file in writeAllESolutionHeader [", file_name, "]!"
            stop
            !
        endif
        !
        !
    end subroutine writeAllESolutionHeader
    !
end module ForwardModeling
!