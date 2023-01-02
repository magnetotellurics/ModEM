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
    public :: JMult, JMult_Tx, JMult_T, JMult_T_Tx, getResidualRMS
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
                !> Set the sum into the current data component
                call JmHat_tx%data( i_data )%set( i_comp, real( lrows_esens, kind=prec ), real( aimag( lrows_esens ), kind=prec ) )
                !
            enddo
            !
        enddo
        !
    end subroutine JMult_Tx
    !
    !> Call JMult_T_Tx for for all transmitters:
    !>     Calculate residual data with predicted data for each transmitter.
    !>     Add the result obtained for each transmitter into a resulting ModelOperator DSigma.
    subroutine JMult_T( sigma, all_data_tx, dsigma )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma
        type( DataGroupTx_t ), dimension(:), intent( in ) :: all_data_tx
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
            dsigma = sigma
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
            call JMult_T_Tx( sigma, all_data_tx( i_tx ), dsigma_tx )
            !
            call dsigma%add( dsigma_tx )
            !
        enddo
        !
        deallocate( dsigma_tx )
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
        !
        !> Initialize dsigma with zeros
        dsigma = sigma
        !
        call dsigma%zeros()
        !
        !> Pointer to the Data Transmitter
        Tx => getTransmitter( tx_data%i_tx )
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
                    tx_data_cvalue = cmplx( data_group%reals( i_comp ), data_group%imaginaries( i_comp ), kind=prec )
                else
                    tx_data_cvalue = cmplx( data_group%reals( i_comp ), R_ZERO, kind=prec )
                endif
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
        !> Set Transmitter's ForwardSolver Omega(Period) and Conductivity
        call Tx%forward_solver%setFrequency( sigma, Tx%period )
        !
        !> Switch Transmitter's source to SourceInteriorForce
        call Tx%setSource( SourceInteriorForce_t( model_operator, sigma, Tx%period ) )
        !
        !> Set E of the transmitter source and create Rhs from it
        call Tx%source%setE( bSrc )
        !
        deallocate( bSrc )
        !
        !> Solve Transmitter's e_sens with the new SourceInteriorForce
        call Tx%solve()
        !
        !> Get dsigma from pMult_t
        call Tx%pMult_t( sigma, dsigma )
        !
    end subroutine JMult_T_Tx
    !
    !> Get DSigma for a single transmitter:
    !>     Create a Rhs from LRows * residual data for all receivers related to the transmitter.
    !>     Solve ESens on the transmitter with SourceInteriorForce and the new Rhs.
    !>     Call Tx%PMult to get a new ModelParameter DSigma for the transmitter.
    function getResidualRMS( predicted, residual ) result( rmsd )
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:), intent( in ) :: predicted
        type( DataGroupTx_t ), allocatable, dimension(:), intent( inout ) :: residual
        !
        real( kind=prec ) :: rmsd
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: n_residual
        !
        ! initialize res
        residual = all_measured_data
        !
        call linCombDataGroupTxArray( ONE, all_measured_data, MinusONE, predicted, residual )
        !
        n_residual = residual
        !
        call normalizeDataGroupTxArray( n_residual, 2 )
        !
        rmsd = sqrt( dotProdDataGroupTxArray( residual, n_residual ) / countDataGroupTxArray( residual ) )
        !
    end function getResidualRMS
    !
end module Sensitivity
!