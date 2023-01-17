!
!> Module with the sensitivity routines JMult, JMult_Tx, JMult_T, JMult_T_Tx 
!
module Sensitivity
    !
    use ForwardModeling
    !
    !> Global Sensitivity Routines
    public :: jobJMult, JMult, JMult_Tx
    public :: jobJMult_T, JMult_T, JMult_T_Tx
    !
contains
    !
    !> Routine to run a full JMult job and deliver the result (JmHat Data) in a text file
    subroutine jobJMult()
        implicit none
        !
        !> Local Data Array to store JMult output
        type( DataGroupTx_t ), allocatable, dimension(:) :: JmHat
        !
        class( ModelParameter_t ), allocatable :: sigma, dsigma
        !
        ! Verbose
        write( *, * ) "     - Start jobJMult"
        !
        !> Read Prior Model File: instantiate dsigma
        if( has_pmodel_file ) then
            !
            call handlePModelFile( dsigma )
            !
            !> Read Model File and instantiate global variables: main_grid, model_operator and Sigma0
            if( has_model_file ) then
                !
                call handleModelFile( sigma )
                !
                call dsigma%setMetric( model_operator%metric )
                !
            else
                stop "Error: jobJMult > Missing Model file!"
            endif
            !
        else
            stop "Error: jobJMult > Missing Prior Model file!"
        endif
        !
        !> Read Data File: instantiate Txs and Rxs and build the Data relation between them
        if( has_data_file ) then
            !
            call handleDataFile()
            !
        else
            stop "Error: jobJMult > Missing Data file!"
        endif
        !
        JmHat = all_measured_data
        !
#ifdef MPI
        !
        call broadcastBasicComponents()
        !
        call masterJMult( sigma, dsigma, JmHat )
        !
        call broadcastFinish
        !
#else
        !
        call createDistributeForwardSolver()
        !
        call JMult( sigma, dsigma, JmHat )
        !
#endif
        !
        !> Write JmHat to the file <jmhat_data_file_name>
        call writeDataGroupTxArray( JmHat, jmhat_data_file_name )
        !
        !> Flush local variables
        call deallocateDataGroupTxArray( JmHat )
        !
        deallocate( sigma, dsigma )
        !
        ! Verbose
        write( *, * ) "     - Finish jobJMult"
        !
    end subroutine jobJMult
    !
    !> Get the JmHat for all transmitters represented in an array of DataGroupTx:
    !>     Call the setFrequency routine to the forward_solver pointer of the transmitter
    !>     Set the transmitter source by calling PMult
    !>     Call JMult_TX for the transmitter
    subroutine JMult( sigma, dsigma, JmHat )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma, dsigma
        type( DataGroupTx_t ), dimension(:), intent( inout ) :: JmHat
        !
        integer :: i_data_tx
        class( Transmitter_t ), pointer :: Tx
        !
        ! Verbose
        !write( *, * ) "          - Start JMult"
        !
        !> Loop over All DataGroupTxs
        do i_data_tx = 1, size( JmHat )
            !
            !> Pointer to the transmitter leading the current data
            Tx => getTransmitter( i_data_tx )
            !
            call solveTx( sigma, Tx )
            !
            !> Switch Transmitter's source to SourceInteriorForce
            call Tx%setSource( Tx%PMult( sigma, dsigma, model_operator ) )
            !
            !> Solve e_sens from with the new Source
            call Tx%solve()
            !
            !> Fill tx_data with JMult routine
            call JMult_Tx( JmHat( i_data_tx ) )
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
    !
    subroutine JMult_Tx( JmHat_tx )
        implicit none
        !
        type( DataGroupTx_t ), intent( inout ) :: JmHat_tx
        !
        class( Vector_t ), allocatable, dimension(:,:) :: lrows
        complex( kind=prec ) :: lrows_x_esens
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
            call Rx%setLRows( Tx, lrows )
            !
            !> Loop over components
            do i_comp = 1, JmHat_tx%data( i_data )%n_comp
                !
                lrows_x_esens = C_ZERO
                !
                !> Loop over polarizations
                do i_pol = 1, Tx%n_pol
                    !
                    !> NECESSARY FOR FULL VECTOR LROWS ????
                    call lrows( i_pol, i_comp )%conjugate()
                    !
                    lrows_x_esens = lrows_x_esens + Tx%e_sens( i_pol )%dotProd( lrows( i_pol, i_comp ) )
                    !
                enddo
                !
                call JmHat_tx%data( i_data )%set( i_comp, -real( lrows_x_esens, kind=prec ), real( aimag( lrows_x_esens ), kind=prec ) )
                !
                !write( *, * ) "JMult Z: ", JmHat_tx%data( i_data )%reals( i_comp ), JmHat_tx%data( i_data )%imaginaries( i_comp )
                !
            enddo
            !
            deallocate( lrows )
            !
            JmHat_tx%data( i_data )%error_bar = .FALSE.
            !
        enddo
        !
    end subroutine JMult_Tx
    !
    !> Routine to run a full JMult_T job and deliver the result (DSigma model) in a text file
    subroutine jobJMult_T()
        implicit none
        !
        class( ModelParameter_t ), allocatable :: sigma, dsigma
        !
        ! Verbose
        write( *, * ) "     - Start jobJMult_T"
        !
        !> Read Model File and instantiate global variables: main_grid, model_operator and Sigma0
        if( has_model_file ) then 
            !
            call handleModelFile( sigma )
            !
        else
            stop "Error: jobJMult_T > Missing Model file!"
        endif
        !
        !> Read Data File: instantiate Txs and Rxs and build the Data relation between them
        !> Initialize a DataGroupTxArray to hold the measured data
        if( has_data_file ) then 
            !
            call handleDataFile()
            !
        else
            stop "Error: jobJMult_T > Missing Data file!"
        endif
        !
#ifdef MPI
        !
        call broadcastBasicComponents()
        !
        call masterJMult_T( sigma, all_measured_data, dsigma )
        !
        call broadcastFinish
        !
#else
        !
        call createDistributeForwardSolver()
        !
        call JMult_T( sigma, all_measured_data, dsigma )
        !
#endif
        !
        !> Write dsigma to <dsigma_file_name> file path
        call dsigma%write( dsigma_file_name )
        !
        !> Flush local variables
        deallocate( sigma, dsigma )
        !
        ! Verbose
        write( *, * ) "     - Finish jobJMult_T"
        !
    end subroutine jobJMult_T
    !
    !> Call JMult_T_Tx for for all transmitters:
    !>     Calculate residual data with predicted data for each transmitter.
    !>     Add the result obtained for each transmitter into a resulting ModelOperator DSigma.
    !
    subroutine JMult_T( sigma, all_data, dsigma )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma
        type( DataGroupTx_t ), dimension(:), intent( in ) :: all_data
        class( ModelParameter_t ), allocatable, intent( out ) :: dsigma
        !
        class( Transmitter_t ), pointer :: Tx
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
            !> Pointer to the transmitter leading the current data
            Tx => getTransmitter( i_tx )
            !
            call solveTx( sigma, Tx )
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
    !
    subroutine JMult_T_Tx( sigma, tx_data, dsigma )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma
        type( DataGroupTx_t ), intent( in ) :: tx_data
        class( ModelParameter_t ), allocatable, intent( inout ) :: dsigma
        !
        class( Vector_t ), allocatable, dimension(:,:) :: lrows
        class( Vector_t ), allocatable, dimension(:) :: bSrc
        class( Transmitter_t ), pointer :: Tx
        class( Receiver_t ), pointer :: Rx
        type( DataGroup_t ) :: data_group
        complex( kind=prec ) :: tx_data_cvalue
        integer :: i_data, i_comp, i_pol
        !
        !> Initialize dsigma with zeros
        allocate( dsigma, source = sigma )
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
            call Rx%setLRows( Tx, lrows )
            !
            !> Loop over the Data components
            do i_comp = 1, data_group%n_comp
                !
                if( Rx%is_complex ) then
                    !
                    tx_data_cvalue = cmplx( data_group%reals( i_comp ), -data_group%imaginaries( i_comp ), kind=prec )
                else
                    tx_data_cvalue = cmplx( data_group%reals( i_comp ), R_ZERO, kind=prec )
                endif
                !
                !write( *, * ) "JMult_T Z: ", tx_data_cvalue
                !
                !> Loop over polarizations
                do i_pol = 1, Tx%n_pol
                    !
                    call lrows( i_pol, i_comp )%mult( tx_data_cvalue )
                    !
                    call bSrc( i_pol )%add( lrows( i_pol, i_comp ) )
                    !
                enddo
                !
            enddo
            !
            deallocate( lrows )
            !
        enddo
        !
        !> NECESSARY FOR FULL VECTOR LROWS ????
        do i_pol = 1, Tx%n_pol
            !
            call bSrc( i_pol )%mult( C_MinusOne )
            !
        enddo
        !
        !> Set Transmitter's ForwardSolver Omega(Period) and Conductivity
        call Tx%forward_solver%setFrequency( sigma, Tx%period )
        !
        !> Switch Transmitter's source to SourceInteriorForce, with trans = .TRUE.
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
        !> Get dsigma from PMult_t
        call Tx%PMult_t( sigma, dsigma )
        !
    end subroutine JMult_T_Tx
    !
end module Sensitivity
!