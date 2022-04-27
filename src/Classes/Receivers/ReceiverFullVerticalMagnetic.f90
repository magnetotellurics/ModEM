! *************
! 
! Derived class to define a Full_Impedance Receiver
! 
! Last modified at 08/06/2021 by Paulo Werdt
! 
! *************
! 
module ReceiverFullVerticalMagnetic
    !
    use Receiver
    !
    type, extends( Receiver_t ), public :: ReceiverFullVerticalMagnetic_t
        !
        ! PROPERTIES HERE
        !
        contains
            !
            final :: ReceiverFullVerticalMagnetic_dtor
            !
            procedure, public :: isEqualRx => isEqualFullVerticalMagnetic
            !
            procedure, public :: predictedData => predictedDataFullVerticalMagnetic
            !
            procedure, public :: write => writeReceiverFullVerticalMagnetic
            !
    end type ReceiverFullVerticalMagnetic_t
    !
    interface ReceiverFullVerticalMagnetic_t
        module procedure ReceiverFullVerticalMagnetic_ctor
    end interface ReceiverFullVerticalMagnetic_t
    !
contains
    !
    function ReceiverFullVerticalMagnetic_ctor( location, rx_type ) result( self )
        implicit none
        !
        real( kind=prec ), intent( in ) :: location(3)
        integer, intent( in ) 			:: rx_type
        !
        type( ReceiverFullVerticalMagnetic_t ) :: self
        !
        character(:), allocatable :: aux_str
        !
        !write(*,*) "Constructor ReceiverFullVerticalMagnetic_t"
        !
        call self%init()
        !
        self%location = location
        !
        self%rx_type = rx_type
        !
        self%n_comp = 2
        self%is_complex = .TRUE.
        !
        allocate( self%EHxy( 3 ) )
        !
        ! components required to get the full impdence tensor response [Zxx, Zxy, Zyx, Zyy]
        self%EHxy(1)%str = "Bx"
        self%EHxy(2)%str = "By"
        self%EHxy(3)%str = "Bz"
        !
        ! components required to get the full impdedance tensor response [Zxx, Zxy, Zyx, Zyy]
        allocate( self%comp_names( 2 ) )
        !
        self%comp_names(1)%str = "TX"
        self%comp_names(2)%str = "TY"
        !
    end function ReceiverFullVerticalMagnetic_ctor
    !
    subroutine ReceiverFullVerticalMagnetic_dtor( self )
        implicit none
        !
        type( ReceiverFullVerticalMagnetic_t ), intent( inout ) :: self
        !
        !write(*,*) "Destructor ReceiverFullVerticalMagnetic_t"
        !
        call self%dealloc()
        !
    end subroutine ReceiverFullVerticalMagnetic_dtor
    !
    function isEqualFullVerticalMagnetic( self, other ) result( equal )
        implicit none
        !
        class( ReceiverFullVerticalMagnetic_t ), intent( in ) :: self
        class( Receiver_t ), intent( in ) :: other
        !
        logical :: equal
        !
        equal = .FALSE.
        !
        select type( other )
            !
            class is( ReceiverFullVerticalMagnetic_t )
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
    end function isEqualFullVerticalMagnetic
    !
    subroutine writeReceiverFullVerticalMagnetic( self )
        implicit none
        !
        class( ReceiverFullVerticalMagnetic_t ), intent( in ) :: self
        !
        write( *, * ) "Write ReceiverFullVerticalMagnetic_t: ", self%id
        !
    end subroutine writeReceiverFullVerticalMagnetic
    !
    subroutine predictedDataFullVerticalMagnetic( self, model_operator, transmitter )
        implicit none
        !
        class( ReceiverFullVerticalMagnetic_t ), intent( inout ) :: self
        class( ModelOperator_t ), intent( in )                   :: model_operator
        class( Transmitter_t ), intent( in )                     :: transmitter
        !
        complex( kind=prec ), allocatable :: BB(:,:), det
        real( kind=prec )                 :: omega
        !
        omega = ( 2.0 * PI / transmitter%period )
        !
        ! Set Vectors Lex, Ley, Lbx, Lby
        call self%evaluationFunction( model_operator, omega )
        !
        allocate( BB( 3, 2 ) )
        !
        BB(1,1) = self%Lbx .dot. transmitter%e_all( 1 )
        BB(2,1) = self%Lby .dot. transmitter%e_all( 1 )
        BB(1,2) = self%Lbx .dot. transmitter%e_all( 2 )
        BB(2,2) = self%Lby .dot. transmitter%e_all( 2 )
        BB(3,1) = self%Lbz .dot. transmitter%e_all( 1 )
        BB(3,2) = self%Lbz .dot. transmitter%e_all( 2 )
        !
        deallocate( self%Lbx )
        deallocate( self%Lby )
        deallocate( self%Lbz )
        !
        !invert horizontal B matrix using Kramer"s rule.
        det = BB( 1, 1 ) * BB( 2, 2 ) - BB( 1, 2 ) * BB( 2, 1 )
        !
        !write(*,*) "det:", det
        !
        allocate( self%I_BB( 2, 2 ) )
        !
        if( det /= 0 ) then
            self%I_BB( 1, 1 ) =  BB( 2, 2 ) / det
            self%I_BB( 2, 2 ) =  BB( 1, 1 ) / det
            self%I_BB( 1, 2 ) = -BB( 1, 2 ) / det
            self%I_BB( 2, 1 ) = -BB( 2, 1 ) / det
        else
            STOP "ReceiverFullImpedance.f90: Determinant is Zero!"
        endif
        !
        allocate( self%response( 2 ) )
        !
        self%response(1) = BB(3,1) * self%I_BB(1,1) + BB(3,2) * self%I_BB(2,1)
        self%response(2) = BB(3,1) * self%I_BB(1,2) + BB(3,2) * self%I_BB(2,2)
        !
        deallocate( BB )
        deallocate( self%I_BB )
        !
        ! WRITE ON PredictedFile.dat
        call self%savePredictedData( transmitter )
        !
        deallocate( self%response )
        !
    end subroutine predictedDataFullVerticalMagnetic
    !
end module ReceiverFullVerticalMagnetic
