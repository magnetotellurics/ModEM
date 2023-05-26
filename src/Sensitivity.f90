!
!> Module with the Sensitivity routines
!
module Sensitivity
    !
    use ForwardModeling
    !
    !> Public module routines
    public :: jobJMult, serialJMult, JMult_Tx
    public :: jobJMult_T, serialJMult_T, JMult_T_Tx
    !
contains
    !
    !> Routine to run a full serialJMult job 
    !> and deliver the result(JmHat) in a text file <jmhat.dat>
    !
    subroutine jobJMult()
        implicit none
        !
        type( DataGroupTx_t ), allocatable, dimension(:) :: JmHat
        !
        class( ModelParameter_t ), allocatable :: sigma, dsigma
        !
        ! Verbose
        !
        write( *, * ) "     - Start jobJMult"
        !
        if( has_pmodel_file ) then
            !
            call handlePModelFile( dsigma )
            !
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
            stop "Error: jobJMult > Missing Perturbation Model file!"
        endif
        !
        if( has_data_file ) then
            !
            call handleDataFile
            !
        else
            stop "Error: jobJMult > Missing Data file!"
        endif
        !
#ifdef MPI
        !
        call broadcastBasicComponents
        !
        call masterSolveAll( sigma )
        !
        call masterJMult( sigma, dsigma, JmHat )
        !
        call broadcastFinish
        !
#else
        !
        call createDistributeForwardSolver
        !
        call solveAll( sigma )
        !
        call serialJMult( sigma, dsigma, JmHat )
        !
#endif
        !
        call writeData( JmHat, jmhat_data_file_name )
        !
        deallocate( sigma, dsigma )
        !
        ! Verbose
        write( *, * ) "     - Finish jobJMult"
        !
    end subroutine jobJMult
    !
    !> Calculate JmHat for all transmitters represented in an array of DataGroupTx:
    !>     Call the setFrequency routine to the forward_solver pointer of the transmitter
    !>     Set the transmitter source by calling PMult
    !>     Call JMult_TX for the transmitter
    !
    subroutine serialJMult( sigma, dsigma, JmHat )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma, dsigma
        type( DataGroupTx_t ), allocatable, dimension(:), intent( out ) :: JmHat
        !
        integer :: i_data_tx
        class( Transmitter_t ), pointer :: Tx
        !
        if( allocated( JmHat ) ) deallocate( JmHat )
        !
        JmHat = all_measured_data
        !
        !> Loop over All DataGroupTxs
        do i_data_tx = 1, size( JmHat )
            !
            !> Pointer to the transmitter leading the current data
            Tx => getTransmitter( i_data_tx )
            !
            call Tx%forward_solver%setFrequency( sigma, Tx%period )
            !
            !> Switch Transmitter's source to SourceInteriorForce
            call Tx%setSource( Tx%PMult( sigma, dsigma, model_operator ) )
            !
            call Tx%solve
            !
            call JMult_Tx( JmHat( i_data_tx ) )
            !
        enddo
        !
    end subroutine serialJMult
    !
    !> Calculate JmHat for a single transmitter and store it in a DataGroupTx:
    !>     By the sum of all LRows * ESens
    !
    subroutine JMult_Tx( JmHat_tx )
        implicit none
        !
        type( DataGroupTx_t ), intent( inout ) :: JmHat_tx
        !
        class( Vector_t ), allocatable :: lrows
        complex( kind=prec ) :: lrows_x_esens
        integer :: i_data, i_comp, i_pol
        class( Transmitter_t ), pointer :: Tx
        class( Receiver_t ), pointer :: Rx
        !
        !> Point to the JmHat's Transmitter
        Tx => getTransmitter( JmHat_tx%i_tx )
        !
        !> Loop over data
        do i_data = 1, size( JmHat_tx%data )
            !
            !> Pointer to the data's Receiver
            Rx => getReceiver( JmHat_tx%data( i_data )%i_rx )
            !
            call Rx%setLRows( Tx )
            !
            !> Loop over components
            do i_comp = 1, JmHat_tx%data( i_data )%n_comp
                !
                lrows_x_esens = C_ZERO
                !
                !> Loop over polarizations
                do i_pol = 1, Tx%n_pol
                    !
                    allocate( lrows, source = Rx%lrows( i_pol, i_comp ) )
                    !
                    !> NECESSARY FOR FULL VECTOR LROWS ????
                    call lrows%conjugate
                    !
                    lrows_x_esens = lrows_x_esens + Tx%e_sens( i_pol )%dotProd( lrows )
                    !
                    deallocate( lrows )
                    !
                enddo
                !
                call JmHat_tx%data( i_data )%set( i_comp, -real( lrows_x_esens, kind=prec ), real( aimag( lrows_x_esens ), kind=prec ) )
                !
                !write( *, * ) "serialJMult Z: ", JmHat_tx%data( i_data )%reals( i_comp ), JmHat_tx%data( i_data )%imaginaries( i_comp )
                !
            enddo
            !
            JmHat_tx%data( i_data )%error_bar = .FALSE.
            !
        enddo
        !
    end subroutine JMult_Tx
    !
    !> Routine to run a full JMult_T job 
    !> and deliver the result(dsigma) in a text file <dsigma.rho>
    !
    subroutine jobJMult_T()
        implicit none
        !
        class( ModelParameter_t ), allocatable :: sigma, dsigma
        !
        ! Verbose
        write( *, * ) "     - Start jobJMult_T"
        !
        if( has_model_file ) then 
            !
            call handleModelFile( sigma )
            !
        else
            stop "Error: jobJMult_T > Missing Model file!"
        endif
        !
        if( has_data_file ) then 
            !
            call handleDataFile
            !
        else
            stop "Error: jobJMult_T > Missing Data file!"
        endif
        !
#ifdef MPI
        !
        call broadcastBasicComponents
        !
        call masterSolveAll( sigma )
        !
        call masterJMult_T( sigma, all_measured_data, dsigma )
        !
        call broadcastFinish
        !
#else
        !
        call createDistributeForwardSolver
        !
        call solveAll( sigma )
        !
        call serialJMult_T( sigma, all_measured_data, dsigma )
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
    !> Initialize dsigma with zeros
    !> Call JMult_T_Tx with measured data for for all transmitters
    !> Add the result obtained for each transmitter into dsigma
    !
    subroutine serialJMult_T( sigma, all_data, dsigma, i_sol, s_hat )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma
        type( DataGroupTx_t ), dimension(:), intent( in ) :: all_data
        class( ModelParameter_t ), allocatable, intent( out ) :: dsigma
        integer, intent( in ), optional :: i_sol
        class( Scalar_t ), allocatable, dimension(:), intent( inout ), optional :: s_hat
        !
        class( Transmitter_t ), pointer :: Tx
        class( ModelParameter_t ), allocatable :: dsigma_tx
        class( Scalar_t ), allocatable :: temp_scalar
        integer :: i_tx, sol_index
        !
        ! Verbose
        !write( *, * ) "          - Start serialJMult_T"
        !
        sol_index = 0
        !
        !> Set i_sol if present
        if( present( i_sol ) ) sol_index = i_sol
        !
        !> Initialize dsigma with zeros
        if( sigma%is_allocated ) then
            !
            if( allocated( dsigma ) ) deallocate( dsigma )
            allocate( dsigma, source = sigma )
            !
            call dsigma%zeros
            !
        else
            stop "Error: serialJMult_T > sigma not allocated"
        endif
        !
        !> Loop over all transmitters
        do i_tx = 1, size( transmitters )
            !
            call JMult_T_Tx( sigma, all_data( i_tx ), dsigma_tx, sol_index )
            !
            if( present( s_hat ) ) then
                !
                call dsigma_tx%getCond( temp_scalar )
                !
                s_hat( i_tx ) = temp_scalar
                !
            endif
            !
            !> Add dsigma_tx to dsigma
            call dsigma%linComb( ONE, ONE, dsigma_tx )
            !
            deallocate( dsigma_tx )
            !
        enddo
        !
        ! Verbose
        !write( *, * ) "          - Finish serialJMult_T"
        !
    end subroutine serialJMult_T
    !
    !> Calculate dsigma for the data_tx's transmitter:
    !>     Create a rhs from LRows * residual data for all receivers related to the transmitter.
    !>     Solve ESens on the transmitter using a transpose SourceInteriorForce, with the new rhs.
    !>     Call Tx%PMult to get a new ModelParameter dsigma.
    !
    subroutine JMult_T_Tx( sigma, tx_data, dsigma, i_sol )
        implicit none
        !
        class( ModelParameter_t ), intent( in ) :: sigma
        type( DataGroupTx_t ), intent( in ) :: tx_data
        class( ModelParameter_t ), allocatable, intent( inout ) :: dsigma
        integer, intent( in ), optional :: i_sol
        !
        class( Vector_t ), allocatable :: lrows
        class( Vector_t ), allocatable, dimension(:) :: bSrc
        class( Transmitter_t ), pointer :: Tx
        class( Receiver_t ), pointer :: Rx
        type( DataGroup_t ) :: data_group
        complex( kind=prec ) :: tx_data_cvalue
        integer :: i_data, i_comp, i_pol, sol_index
        !
        sol_index = 0
        !
        !> Set i_sol if present
        if( present( i_sol ) ) sol_index = i_sol
        !
        !> Initialize dsigma with zeros
        allocate( dsigma, source = sigma )
        !
        call dsigma%zeros
        !
        !> Pointer to the tx_data's Transmitter
        Tx => getTransmitter( tx_data%i_tx )
        !
        Tx%i_sol = sol_index
        !
        !> Initialize bSrc( n_pol ) with zeros
        allocate( cVector3D_SG_t :: bSrc( Tx%n_pol ) )
        !
        do i_pol = 1, Tx%n_pol
            !
            bSrc( i_pol ) = cVector3D_SG_t( sigma%metric%grid, EDGE )
            call bSrc( i_pol )%zeros
            !
        enddo
        !
        !> Loop over all data of the transmitter
        do i_data = 1, size( tx_data%data )
            !
            data_group = tx_data%data( i_data )
            !
            !> Pointer to the data's Receiver
            Rx => getReceiver( tx_data%data( i_data )%i_rx )
            !
            call Rx%setLRows( Tx )
            !
            !> Loop over the data components
            do i_comp = 1, data_group%n_comp
                !
                if( Rx%is_complex ) then
                    !
                    tx_data_cvalue = cmplx( data_group%reals( i_comp ), -data_group%imaginaries( i_comp ), kind=prec )
                else
                    tx_data_cvalue = cmplx( data_group%reals( i_comp ), R_ZERO, kind=prec )
                endif
                !
                !write( *, * ) "serialJMult_T Z: ", tx_data_cvalue
                !
                !> Loop over polarizations
                do i_pol = 1, Tx%n_pol
                    !
                    allocate( lrows, source = Rx%lrows( i_pol, i_comp ) )
                    !
                    call lrows%mult( tx_data_cvalue )
                    !
                    call bSrc( i_pol )%add( lrows )
                    !
                    deallocate( lrows )
                    !
                enddo
                !
            enddo
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
        call Tx%forward_solver%setFrequency( sigma, Tx%period )
        !
        !> Switch Transmitter's source to SourceInteriorForce, with transpose = .TRUE.
        call Tx%setSource( SourceInteriorForce_t( model_operator, sigma, Tx%period, .TRUE. ) )
        !
        call Tx%source%setE( bSrc )
        !
        deallocate( bSrc )
        !
        !> Solve Transmitter's e_sens with the new SourceInteriorForce
        call Tx%solve
        !
        call Tx%PMult_t( sigma, dsigma )
        !
    end subroutine JMult_T_Tx
    !
end module Sensitivity
!