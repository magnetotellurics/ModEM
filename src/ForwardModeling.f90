!
!> Module with the sensitivity routines JMult, JMult_Tx, JMult_T, JMult_T_Tx 
!
module ForwardModeling
    !
    use Constants
    !
    use FileUnits
    !
    use cVector3D_SG
    !
    use GlobalVariables
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
    use ReceiverArray
    !
    use TransmitterMT
    use TransmitterCSEM
    use TransmitterArray
    !
    use DataGroupArray
    use DataGroupTxArray
    !
    !> Global FWD Routines
    public :: runForwardModeling, writeEsolutionHeader
    !
contains
    !
    !> runForwardModeling: Calculate ESolution for all transmitters and...
    !>     Receive a flag indicating whether it will be used for Adjoint:
    !>         - False or None: Calculate only the predicted data for each transmitter-receiver pair.
    !>         - True: Calculate LRows in the receivers after calculating predicted data.
    !> Obs.: Require the previous definition of the global ForwardSolver (createForwardSolver())
    !
    subroutine runForwardModeling( sigma )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma
        !
        class( Transmitter_t ), pointer :: Tx
        class( Receiver_t ), pointer :: Rx
        integer :: i_tx, n_tx, i_rx
        !
        ! Verbose
        write( *, * ) "          - Start Forward Modeling"
        !
        !>
        n_tx = size( transmitters )
        !
        ! Verbose
        write( *, * ) "          > Write ESolution to file: [", e_solution_file_name, "]"
        !
        !> Write the first header of the ESolution binary file, according to the first transmitter
        call writeEsolutionHeader( n_tx, transmitters(1)%Tx%n_pol )
        !
        !> Loop over all Transmitters
        do i_tx = 1, n_tx
            !
            !> Pointer to the Transmitter
            Tx => getTransmitter( i_tx )
            !
            !> Set Transmitter's ForwardSolver Period) and Conductivity
            call Tx%forward_solver%setFrequency( sigma, Tx%period )
            !
            !> Instantiate Transmitter's Source - According to transmitter type or chosen via control file
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
                    stop "Error: ForwardModeling > Unclassified Transmitter"
                !
            end select
            !
            !> Build Source E according to source type
            call Tx%source%createE()
            !
            !> Solve e_all for this Transmitter
            call Tx%solve()
            !
            ! Verbose
            write( *, * ) "                   - Receivers calculation"
            !
            !> Loop for each Receiver related to this Transmitter
            do i_rx = 1, size( Tx%receiver_indexes )
                !
                !> Pointer to the Tx Receiver
                Rx => getReceiver( Tx%receiver_indexes( i_rx ) )
                !
                call Rx%predictedData( Tx )
                !
                call all_predicted_data( i_tx )%set( Rx%data_group )
                !
            enddo
            !
        enddo
        !
        ! Verbose
        write( *, * ) "          - Finish Forward Modeling"
        !
    end subroutine runForwardModeling
    !
    !> No subroutine briefing
    subroutine writeEsolutionHeader( n_tx, nMode )
        implicit none
        !
        integer, intent( in ) :: n_tx, nMode
        integer :: ios
        character(len=20) :: version
        !
        version = "Modem-OO"
        !
        open( ioESolution, file = e_solution_file_name, action = "write", form = "unformatted", iostat = ios)
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
            write( *, * ) "Error opening file in writeEsolutionHeader [", e_solution_file_name, "]!"
            stop
            !
        endif
        !
        !
    end subroutine writeEsolutionHeader
    !
end module ForwardModeling
!