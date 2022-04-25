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
    function ReceiverOffDiagonalImpedance_ctor( location, type_name ) result( self )
        implicit none
        !
        real( kind=prec ), intent( in )                   :: location(3)
		character(:), allocatable, optional, intent( in ) :: type_name
        !
        type( ReceiverOffDiagonalImpedance_t ) :: self
        !
        ! write(*,*) "Constructor ReceiverOffDiagonalImpedance_t"
        !
        call self%init()
        !
        self%location = location
        !
        if( present( type_name ) ) then
            self%type_name = type_name
        else
            self%type_name = "ReceiverOffDiagonalImpedance"
        endif
        !
        self%DATA_TITLE = "Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error"
        !
        self%n_comp = 4
        self%is_complex = .TRUE.
        !
        allocate( self%EHxy( self%n_comp ) )
        !
        ! components required to get the full impdence tensor response [Zxx, Zxy, Zyx, Zyy]
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
    subroutine predictedDataOffDiagonalImpedance( self, model_operator, transmitter )
        implicit none
        !
        class( ReceiverOffDiagonalImpedance_t ), intent( inout ) :: self
        class( ModelOperator_t ), intent( in )                         :: model_operator
        class( Transmitter_t ), intent( in )                            :: transmitter
        !
        class( cVector_t ), allocatable    :: e_tx_pol_1, e_tx_pol_2
        complex( kind=prec ), allocatable :: BB(:,:), det
        real( kind=prec )                      :: omega
        integer                                    :: i, j, ij
        !
        omega = ( 2.0 * PI / transmitter%period )
        !
        ! Set Vectors Lex, Ley, Lbx, Lby
        call self%evaluationFunction( model_operator, omega )
        !
        ! get e_all from the Tx 1st polarization
        allocate( e_tx_pol_1, source = transmitter%e_all( 1 ) )
        !
        ! get e_all from the Tx 2nd polarization
        allocate( e_tx_pol_2, source = transmitter%e_all( 2 ) )
        !
        !
        allocate( complex( kind=prec ) :: self%EE( 2, 2 ) )
        !
        self%EE(1,1) = self%Lex .dot. e_tx_pol_1
        self%EE(2,1) = self%Ley .dot. e_tx_pol_1
        self%EE(1,2) = self%Lex .dot. e_tx_pol_2
        self%EE(2,2) = self%Ley .dot. e_tx_pol_2
        !
        allocate( complex( kind=prec ) :: BB( 2, 2 ) )
        !
        BB(1,1) = self%Lbx .dot. e_tx_pol_1
        BB(2,1) = self%Lby .dot. e_tx_pol_1
        BB(1,2) = self%Lbx .dot. e_tx_pol_2
        BB(2,2) = self%Lby .dot. e_tx_pol_2
        !
        deallocate( e_tx_pol_1 )
        deallocate( e_tx_pol_2 )
        !
        !invert horizontal B matrix using Kramer's rule.
        det = BB(1,1) * BB(2,2) - BB(1,2) * BB(2,1)
        !
        allocate( complex( kind=prec ) :: self%I_BB( 2, 2 ) )
        !
        if( det /= 0 ) then
            self%I_BB( 1, 1 ) = BB( 2, 2 ) / det
            self%I_BB( 2, 2 ) = BB( 1, 1 ) / det
            self%I_BB( 1, 2 ) = -BB( 1, 2 ) / det
            self%I_BB( 2, 1 ) = -BB( 2, 1 ) / det
        else
            STOP "ReceiverOffDiagonalImpedance.f90: Determinant is Zero!"
        endif
        !
        allocate( complex( kind=prec ) :: self%response( 2 ) )
        !
        self%response(1) = self%EE(1,1) * self%I_BB(1,2) + self%EE(1,2) * self%I_BB(2,2)
        self%response(2) = self%EE(2,1) * self%I_BB(1,1) + self%EE(2,2) * self%I_BB(2,1)
        !
        ! WRITE ON PredictedFile.dat
        call self%savePredictedData( transmitter )
        !
        deallocate( self%EE )
        deallocate( BB )
        deallocate( self%I_BB )
        deallocate( self%response )
        !
    end subroutine predictedDataOffDiagonalImpedance
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
