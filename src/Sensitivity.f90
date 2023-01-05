!
!> Module with the sensitivity routines JMult, JMult_Tx, JMult_T, JMult_T_Tx 
!
module Sensitivity
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
    !> Global Sensitivity Routines
    public :: JMult, JMult_Tx, JMult_T, JMult_T_Tx
    !
contains
    !
    !> Get the JmHat for all transmitters represented in an array of DataGroupTx:
    !>     Call the setFrequency routine to the forward_solver pointer of the transmitter
    !>     Set the transmitter source by calling pMult
    !>     Call JMult_TX for the transmitter
    subroutine JMult( sigma, dsigma, JmHat )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma, dsigma
        type( DataGroupTx_t ), dimension(:), intent( inout ) :: JmHat
        !
        integer :: i_dtx
        class( Transmitter_t ), pointer :: Tx
        !
        ! Verbose
        !write( *, * ) "          - Start JMult"
        !
        !> Loop over All DataGroupTxs
        do i_dtx = 1, size( JmHat )
            !
            !> Pointer to the transmitter leading the current data
            Tx => getTransmitter( i_dtx )
            !
            !> Set Transmitter's ForwardSolver Omega(Period) and Conductivity
            call Tx%forward_solver%setFrequency( sigma, Tx%period )
            !
            !> Switch Transmitter's source to SourceInteriorForce
            call Tx%setSource( Tx%pMult( sigma, dsigma, model_operator ) )
            !
            !> Solve e_sens from with the new Source
            call Tx%solve()
			!
			call Tx%e_sol(1)%print( 3001, "JMult ESol(1)" )
			call Tx%e_sol(2)%print( 3002, "JMult ESol(2)" )
			!
			call Tx%e_sens(1)%print( 4001, "JMult ESens(1)" )
			call Tx%e_sens(2)%print( 4002, "JMult ESens(2)" )
			!
            !> Fill tx_data with JMult routine
            call JMult_Tx( JmHat( i_dtx ) )
            !
        enddo
        !
        ! Verbose
        !write( *, * ) "          - Finish JMult"
        !
    end subroutine JMult
    !
    !> Calculate JmHat for a single transmitter and store it in a DataGroupTx:
    !>     By the sum of all LRows * ESens
    subroutine JMult_Tx( JmHat_tx )
        implicit none
        !
        type( DataGroupTx_t ), intent( inout ) :: JmHat_tx
        !
        complex( kind=prec ) :: lrows_esens
        integer :: i_data, i_comp, i_pol
        class( Transmitter_t ), pointer :: Tx
        class( Receiver_t ), pointer :: Rx
        !
        !> Auxiliary pointer to modify a transmitter.
        Tx => getTransmitter( JmHat_tx%i_tx )
        !
        !> Loop over data
        do i_data = 1, size( JmHat_tx%data )
            !
            !> Pointer to the Data Receiver
            Rx => getReceiver( JmHat_tx%data( i_data )%i_rx )
            !
            call Rx%setLRows( Tx )
            !
            !> Loop over components
            do i_comp = 1, JmHat_tx%data( i_data )%n_comp
                !
                lrows_esens = C_ZERO
                !
                !> Loop over polarizations
                do i_pol = 1, Tx%n_pol
                    !
                    !> LRows .dot. ESens
                    lrows_esens = lrows_esens + Rx%lrows( i_pol, i_comp )%dotProd( Tx%e_sens( i_pol ) )
                    !
                enddo
                !
                !> Set the sum into the current data component, according to type
                !if( Rx%is_complex ) then
                    call JmHat_tx%data( i_data )%set( i_comp, real( lrows_esens, kind=prec ), real( aimag( lrows_esens ), kind=prec ) )
                !else
                    !call JmHat_tx%data( i_data )%set( i_comp, real( lrows_esens, kind=prec ), R_ZERO )
                !endif
                !
            enddo
            !
            JmHat_tx%data( i_data )%is_complex = Rx%is_complex 
            !
            JmHat_tx%data( i_data )%error_bar = .FALSE.
            !
        enddo
        !
    end subroutine JMult_Tx
    !
    !> Call JMult_T_Tx for for all transmitters:
    !>     Calculate residual data with predicted data for each transmitter.
    !>     Add the result obtained for each transmitter into a resulting ModelOperator DSigma.
    subroutine JMult_T( sigma, all_data, dsigma )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma
        type( DataGroupTx_t ), dimension(:), intent( in ) :: all_data
        class( ModelParameter_t ), allocatable, intent( inout ) :: dsigma
        !
        class( ModelParameter_t ), allocatable :: dsigma_tx
        integer :: i_tx
        !
        ! Verbose
        !write( *, * ) "          - Start JMult_T"
        !
        !> Initialize dsigma with Zeros
        if( sigma%is_allocated ) then
            !
            if( allocated( dsigma ) ) deallocate( dsigma )
            allocate( dsigma, source = sigma )
            !
            call dsigma%zeros()
            !
        else
            stop "Error: JMult_T > sigma not allocated"
        endif
        !
        !> Loop over all transmitters
        do i_tx = 1, size( transmitters )
            !
            !> Set current tx_dsigma from JMult_T_Tx
            call JMult_T_Tx( sigma, all_data( i_tx ), dsigma_tx )
            !
			!> Add dsigma_tx to dsigma
            call dsigma%linComb( ONE, ONE, dsigma_tx )
            !
            deallocate( dsigma_tx )
            !
        enddo
        !
        ! Verbose
        !write( *, * ) "          - Finish JMult_T"
        !
    end subroutine JMult_T
    !
    !> Get DSigma for a single transmitter:
    !>     Create a Rhs from LRows * residual data for all receivers related to the transmitter.
    !>     Solve ESens on the transmitter with SourceInteriorForce and the new Rhs.
    !>     Call Tx%PMult to get a new ModelParameter DSigma for the transmitter.
    subroutine JMult_T_Tx( sigma, tx_data, dsigma )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma
        type( DataGroupTx_t ), intent( in ) :: tx_data
        class( ModelParameter_t ), allocatable, intent( inout ) :: dsigma
        !
        class( Vector_t ), allocatable, dimension(:) :: bSrc
        class( Transmitter_t ), pointer :: Tx
        class( Receiver_t ), pointer :: Rx
        type( DataGroup_t ) :: data_group
        complex( kind=prec ) :: tx_data_cvalue
        integer :: i_data, i_comp, i_pol
        complex( kind=prec ) :: comega
        !
        !> Initialize dsigma with zeros
        allocate( dsigma, source = sigma )
        !
        call dsigma%zeros()
        !
        !> Pointer to the Data Transmitter
        Tx => getTransmitter( tx_data%i_tx )
        !
        comega = cmplx( 0.0, 1./ ( 2.0 * PI / Tx%period ), kind=prec )
        !
        !> Initialize bSrc( n_pol ) with zeros
        allocate( cVector3D_SG_t :: bSrc( Tx%n_pol ) )
        !
        do i_pol = 1, Tx%n_pol
            !
            bSrc( i_pol ) = cVector3D_SG_t( sigma%metric%grid, EDGE )
            call bSrc( i_pol )%zeros()
            !
        enddo
        !
        !> Loop over DataGroups of the transmitter
        do i_data = 1, size( tx_data%data )
            !
            data_group = tx_data%data( i_data )
            !
            !> Pointer to the Data Receiver
            Rx => getReceiver( tx_data%data( i_data )%i_rx )
            !
            !> Calculate lrows for this Receiver
            call Rx%setLRows( Tx )
            !
            !> Loop over the Data components
            do i_comp = 1, data_group%n_comp
                !
                if( Rx%is_complex ) then
                    tx_data_cvalue = cmplx( data_group%reals( i_comp ), -data_group%imaginaries( i_comp ), kind=prec )
                else
                    tx_data_cvalue = cmplx( data_group%reals( i_comp ), R_ZERO, kind=prec )
                endif
				!
				!
				write( *, * ) tx_data_cvalue
				!
                !> Loop over polarizations
                do i_pol = 1, Tx%n_pol
                    !
                    call Rx%lrows( i_pol, i_comp )%mult( tx_data_cvalue )
                    !
                    call bSrc( i_pol )%add( Rx%lrows( i_pol, i_comp ) )
                    !
                enddo
                !
            enddo
            !
        enddo
		!
		call bSrc( 1 )%print( 2001, "bSrc( 1 )" )
		call bSrc( 2 )%print( 2002, "bSrc( 2 )" )
		!
		!stop
        !
        !> Set Transmitter's ForwardSolver Omega(Period) and Conductivity
        call Tx%forward_solver%setFrequency( sigma, Tx%period )
        !
        !> Switch Transmitter's source to SourceInteriorForce
        call Tx%setSource( SourceInteriorForce_t( model_operator, sigma, Tx%period, .TRUE. ) )
        !
        !> Set E of the transmitter source and create Rhs from it
        call Tx%source%setE( bSrc )
        !
        deallocate( bSrc )
        !
        !> Solve Transmitter's e_sens with the new SourceInteriorForce
        call Tx%solve()
		!
		call Tx%e_sol(1)%print( 3001, "JMult_T ESol(1)" )
		call Tx%e_sol(2)%print( 3002, "JMult_T ESol(2)" )
        !
		call Tx%e_sens(1)%print( 4001, "JMult_T ESens(1)" )
		call Tx%e_sens(2)%print( 4002, "JMult_T ESens(2)" )
        !
        !> Get dsigma from pMult_t
        call Tx%pMult_t( sigma, dsigma )
        !
    end subroutine JMult_T_Tx
    !
end module Sensitivity
!