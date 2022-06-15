! *************
! 
! Derived class to define a Full_Impedance Receiver
! 
! Last modified at 08/06/2021 by Paulo Werdt
! 
! *************
! 
module ReceiverOffDiagonalImpedance
    !
    use Receiver
    use DataHandleMT
    !
    type, extends( Receiver_t ), public :: ReceiverOffDiagonalImpedance_t
        !
        ! PROPERTIES HERE
        !
        contains
            !
            final :: ReceiverOffDiagonalImpedance_dtor
            !
            procedure, public :: isEqualRx => isEqualOffDiagonalImpedance
            !
            procedure, public :: predictedData => predictedDataOffDiagonalImpedance
            !
            procedure, public :: savePredictedData => savePredictedDataOffDiagonalImpedance
            !
            procedure, public :: write => writeReceiverOffDiagonalImpedance
            !
    end type ReceiverOffDiagonalImpedance_t
    !
    interface ReceiverOffDiagonalImpedance_t
        module procedure ReceiverOffDiagonalImpedance_ctor
    end interface ReceiverOffDiagonalImpedance_t
    !
contains
    !
    function ReceiverOffDiagonalImpedance_ctor( location, rx_type ) result( self )
        implicit none
        !
        real( kind=prec ), intent( in ) :: location(3)
        integer, optional, intent( in ) :: rx_type
        !
        type( ReceiverOffDiagonalImpedance_t ) :: self
        !
        character(:), allocatable :: aux_str
        !
        ! write(*,*) "Constructor ReceiverOffDiagonalImpedance_t"
        !
        call self%init()
        !
        self%location = location
        !
        self%rx_type = rx_type
        !
        self%n_comp = 4
        self%is_complex = .TRUE.
        !
        allocate( self%EHxy( 4 ) )
        !
        ! components required to get the full impdence tensor self%response [Zxx, Zxy, Zyx, Zyy]
        self%EHxy(1)%str = "Ex"
        self%EHxy(2)%str = "Ey"
        self%EHxy(3)%str = "Bx"
        self%EHxy(4)%str = "By"
        !
    end function ReceiverOffDiagonalImpedance_ctor
    !
    subroutine ReceiverOffDiagonalImpedance_dtor( self )
        implicit none
        !
        type( ReceiverOffDiagonalImpedance_t ), intent( inout ) :: self
        !
        ! write(*,*) "Destructor ReceiverOffDiagonalImpedance_t"
        !
        call self%dealloc()
        !
    end subroutine ReceiverOffDiagonalImpedance_dtor
    !
    function isEqualOffDiagonalImpedance( self, other ) result( equal )
        implicit none
        !
        class( ReceiverOffDiagonalImpedance_t ), intent( in ) :: self
        class( Receiver_t ), intent( in ) :: other
        !
        logical :: equal
        !
        equal = .FALSE.
        !
        select type( other )
            !
            class is( ReceiverOffDiagonalImpedance_t )
                !
                if( self%code == other%code .AND.   &
                    self%location(1) == other%location(1) .AND.    &
                    self%location(2) == other%location(2) .AND.    &
                    self%location(3) == other%location(3) ) then
                    equal = .TRUE.
                endif
                !
            class default
                equal = .FALSE.
            !
        end select
        !
    end function isEqualOffDiagonalImpedance
    !
    subroutine predictedDataOffDiagonalImpedance( self, transmitter )
        implicit none
        !
        class( ReceiverOffDiagonalImpedance_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in )                     :: transmitter
        !
        complex( kind=prec ) :: comega, det
        !
        complex( kind=prec ), allocatable :: BB(:,:), I_BB(:,:), EE(:,:)
        !
        integer :: i, j, ij
        !
        comega = cmplx( 0.0, 1./ ( 2.0 * PI / transmitter%period ), kind=prec )
        !
        !
        allocate( EE( 2, 2 ) )
        !
        select type( tx_e_1 => transmitter%e_all( 1 ) )
            class is( cVector3D_SG_t )
                !
                select type( tx_e_2 => transmitter%e_all( 2 ) )
                    class is( cVector3D_SG_t )
                        !
                        EE(1,1) = dotProdSparse( self%Lex, tx_e_1 )
                        EE(2,1) = dotProdSparse( self%Ley, tx_e_1 )
                        EE(1,2) = dotProdSparse( self%Lex, tx_e_2 )
                        EE(2,2) = dotProdSparse( self%Ley, tx_e_2 )
                        !
                        allocate( BB( 2, 2 ) )
                        !
                        BB(1,1) = dotProdSparse( self%Lbx, tx_e_1 )
                        BB(2,1) = dotProdSparse( self%Lby, tx_e_1 )
                        BB(1,2) = dotProdSparse( self%Lbx, tx_e_2 )
                        BB(2,2) = dotProdSparse( self%Lby, tx_e_2 )
                        BB = isign * BB * comega
                        !
                        !invert horizontal B matrix using Kramer's rule.
                        det = BB(1,1) * BB(2,2) - BB(1,2) * BB(2,1)
                        !
                        allocate( I_BB( 2, 2 ) )
                        !
                        if( det /= 0 ) then
                            I_BB( 1, 1 ) = BB( 2, 2 ) / det
                            I_BB( 2, 2 ) = BB( 1, 1 ) / det
                            I_BB( 1, 2 ) = -BB( 1, 2 ) / det
                            I_BB( 2, 1 ) = -BB( 2, 1 ) / det
                        else
                            STOP "ReceiverOffDiagonalImpedance.f90: Determinant is Zero!"
                        endif
                        !
                        deallocate( BB )
                        !
                        allocate( self%response( 2 ) )
                        !
                        self%response(1) = EE(1,1) * I_BB(1,2) + EE(1,2) * I_BB(2,2)
                        self%response(2) = EE(2,1) * I_BB(1,1) + EE(2,2) * I_BB(2,1)
                        !
                        deallocate( EE )
                        deallocate( I_BB )
                        !
                        ! WRITE ON PredictedFile.dat
                        call self%savePredictedData( transmitter )
                        !
                        deallocate( self%response )
                        !
                    class default
                        stop "evaluationFunctionRx: Unclassified transmitter%e_all_2"
                end select
                !
            class default
                stop "evaluationFunctionRx: Unclassified transmitter%e_all_1"
        end select
        !
    end subroutine predictedDataOffDiagonalImpedance
    !
    subroutine savePredictedDataOffDiagonalImpedance( self, tx )
        implicit none
        !
        class( ReceiverOffDiagonalImpedance_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in ) :: tx
        !
        character(:), allocatable :: code, component
        real( kind=prec )         :: period, real_part, imaginary, rx_location(3)
        integer                   :: i, rx_type
        !
        !#Period(s) Code GG_Lat GG_Lon X(m) Y(m) self%response(m) Component Real Imag Error
        !
        if( allocated( self%predicted_data ) ) call deallocateDataHandleArray( self%predicted_data )
        !
        do i = 1, self%n_comp
            !
            rx_type = int( self%rx_type )
            period = real( tx%period, kind=prec )
            code = trim( self%code )
            rx_location = (/real( self%location( 1 ), kind=prec ), real( self%location( 2 ), kind=prec ), real( self%location( 3 ), kind=prec )/)
            component = trim( self%comp_names( i )%str )
            real_part = real( self%response( i ), kind=prec )
            imaginary = real( imag( self%response( i ) ), kind=prec )
            !
            call updateDataHandleArray( self%predicted_data, DataHandleMT_t( rx_type, code, component, period, rx_location, real_part, imaginary ) )
            !
        enddo
        !
    end subroutine savePredictedDataOffDiagonalImpedance
    !
    subroutine writeReceiverOffDiagonalImpedance( self )
        implicit none
        !
        class( ReceiverOffDiagonalImpedance_t ), intent( in ) :: self
        !
        write(*,*) "Write ReceiverOffDiagonalImpedance_t: ", self%id
        !
    end subroutine writeReceiverOffDiagonalImpedance
    !
end module ReceiverOffDiagonalImpedance
