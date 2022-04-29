! *************
! 
! Derived class to define a Single E or B Field Receiver
! 
! Last modified at 30/06/2021 by Paulo Werdt
! 
! *************
! 
module ReceiverSingleField
    !
    use Receiver
    !
    type, extends( Receiver_t ), public :: ReceiverSingleField_t
        !
        real( kind=prec ) :: azimuth
        !
        contains
            !
            final :: ReceiverSingleField_dtor
            !
            procedure, public :: isEqualRx => isEqualSingleField
            !
            procedure, public :: predictedData => predictedDataSingleField
            !
            procedure, public :: write => writeReceiverSingleField
            !
    end type ReceiverSingleField_t
    !
    interface ReceiverSingleField_t
        module procedure ReceiverSingleField_ctor
    end interface ReceiverSingleField_t
    !
contains
    !
    function ReceiverSingleField_ctor( location, azimuth, rx_type ) result( self )
        implicit none
        !
        real( kind=prec ), intent( in ) :: location(3)
        real( kind=prec ), intent( in ) :: azimuth
        integer, intent( in )           :: rx_type
        !
        type( ReceiverSingleField_t )    :: self
        !
        ! write(*,*) "Constructor ReceiverSingleField_t"
        !
        call self%init()
        !
        self%location = location
        self%azimuth = azimuth
        !
        self%rx_type = rx_type
        !
        self%n_comp = 1
        self%is_complex = .TRUE.
        !
        !
        allocate( self%EHxy( self%n_comp ) )
        !
        ! components required to get the full impdence tensor response [Zxx, Zxy, Zyx, Zyy]
        !
        if( azimuth == 1.0 ) self%EHxy(1)%str = "Ex"
        if( azimuth == 2.0 ) self%EHxy(1)%str = "Ey"
        if( azimuth == 3.0 ) self%EHxy(1)%str = "Bx"
        if( azimuth == 4.0 ) self%EHxy(1)%str = "By"
        if( azimuth == 5.0 ) self%EHxy(1)%str = "Bz"
        !
    end function ReceiverSingleField_ctor
    !
    subroutine ReceiverSingleField_dtor( self )
        implicit none
        !
        type( ReceiverSingleField_t ), intent( inout ) :: self
        !
        ! write(*,*) "Destructor ReceiverSingleField_t"
        !
        call self%dealloc()
        !
    end subroutine ReceiverSingleField_dtor
    !
    function isEqualSingleField( self, other ) result( equal )
        implicit none
        !
        class( ReceiverSingleField_t ), intent( in ) :: self
        class( Receiver_t ), intent( in ) :: other
        !
        logical :: equal
        !
        equal = .FALSE.
        !
        select type( other )
            !
            class is( ReceiverSingleField_t )
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
    end function isEqualSingleField
    !
    subroutine predictedDataSingleField( self, transmitter )
        implicit none
        !
        class( ReceiverSingleField_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in )            :: transmitter
        !
        complex( kind=prec ) :: det, ctemp
        !
        write(*,*) "Implement predictedData ReceiverSingleField_t: ", self%id
        !
        allocate( complex( kind=prec ) :: self%response( 1 ) )
        !
        ! FIND BEST WAY TO IMPLEMENT
        !response(1) = dotProd_noConj_scvector_f( Lex,ef%pol(1) )
        !
    end subroutine predictedDataSingleField
    !
    subroutine writeReceiverSingleField( self )
        implicit none
        !
        class( ReceiverSingleField_t ), intent( in ) :: self
        !
        write(*,*) "Write ReceiverSingleField_t: ", self%id
        !
    end subroutine writeReceiverSingleField
    !
end module ReceiverSingleField
