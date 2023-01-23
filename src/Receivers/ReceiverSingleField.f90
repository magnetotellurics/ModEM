!
!> Derived class to define a Single E or B Field Receiver
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
            procedure, public :: setLRows => setLRowsSingleField
            !
            procedure, public :: predictedData => predictedDataSingleField
            !
            procedure, public :: isEqualRx => isEqualSingleField
            !
            procedure, public :: print => printReceiverSingleField
            !
    end type ReceiverSingleField_t
    !
    interface ReceiverSingleField_t
        module procedure ReceiverSingleField_ctor
    end interface ReceiverSingleField_t
    !
contains
    !
    !> No function briefing
    function ReceiverSingleField_ctor( location, azimuth, rx_type ) result( self )
        implicit none
        !
        real( kind=prec ), intent( in ) :: location(3)
        real( kind=prec ), intent( in ) :: azimuth
        integer, intent( in ) :: rx_type
        !
        type( ReceiverSingleField_t ) :: self
        !
        integer :: i, asize
        !
        !> write( *, * ) "Constructor ReceiverSingleField_t"
        !
        call self%init()
        !
        self%location = location
        self%azimuth = azimuth
        !
        self%rx_type = rx_type
        !
        self%n_comp = 1
        !
        self%is_complex = .TRUE.
        !
        !> components required to get the full impedance evaluation vectors [Ex, Ey, Bx, By]
        if( allocated( self%EHxy ) ) then
            !
            asize = size( self%EHxy )
            do i = asize, 1, -(1)
                deallocate( self%EHxy(i)%str )
            enddo
            deallocate( self%EHxy )
            !
        endif
        !
        allocate( self%EHxy(1) )
        !
        if( azimuth == 1.0 ) then
            self%EHxy(1)%str = "Ex"
        endif
        !
        if( azimuth == 2.0 ) then
            self%EHxy(1)%str = "Ey"
        endif
        !
        if( azimuth == 3.0 ) then
            self%EHxy(1)%str = "Bx"
        endif
        !
        if( azimuth == 4.0 ) then
            self%EHxy(1)%str = "By"
        endif
        !
        if( azimuth == 5.0 ) then
            self%EHxy(1)%str = "Bz"
        endif
        !
        !> components required to get the full impedance tensor self%response [Zxx, Zxy, Zyx, Zyy]
        if( allocated( self%comp_names ) ) then
            !
            asize = size( self%comp_names )
            do i = asize, 1, -(1)
                deallocate( self%comp_names(i)%str )
            enddo
            deallocate( self%comp_names )
            !
        endif
        allocate( self%comp_names(1) )
        !
        if( azimuth == 1.0 ) self%comp_names(1)%str = "Ex_Field"
        if( azimuth == 2.0 ) self%comp_names(1)%str = "Ey_Field"
        if( azimuth == 3.0 ) self%comp_names(1)%str = "Bx_Field"
        if( azimuth == 4.0 ) self%comp_names(1)%str = "By_Field"
        if( azimuth == 5.0 ) self%comp_names(1)%str = "Bz_Field"
        !
    end function ReceiverSingleField_ctor
    !
    !> No subroutine briefing
    subroutine ReceiverSingleField_dtor( self )
        implicit none
        !
        type( ReceiverSingleField_t ), intent( inout ) :: self
        !
        !> write( *, * ) "Destructor ReceiverSingleField_t"
        !
        call self%dealloc()
        !
    end subroutine ReceiverSingleField_dtor
    !
    !> No subroutine briefing
    subroutine setLRowsSingleField( self, transmitter )
        implicit none
        !
        class( ReceiverSingleField_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in ) :: transmitter
        !
        complex( kind=prec ) :: comega
        !
        comega = cmplx( 0.0, 1./ ( 2.0 * PI / transmitter%period ), kind=prec )
        !
        !> It's not needed for the LRows calculation itself
        !> but it's called here to maintain the general encapsulation of the program
        !> And also to debug the predicted data
        call self%predictedData( transmitter )
        !
        !> Allocate LRows matrix [ npol = 1, n_comp = 1 ]
        if( allocated( self%lrows ) ) deallocate( self%lrows )
        allocate( cVector3D_SG_t :: self%lrows( transmitter%n_pol, self%n_comp ) )
        !
        if( self%azimuth == 1.0 ) self%lrows( 1, 1 ) = self%Lex%getFullVector()
        if( self%azimuth == 2.0 ) self%lrows( 1, 1 ) = self%Ley%getFullVector()
        !
        if( self%azimuth == 3.0 ) self%lrows( 1, 1 ) = self%Lbx%getFullVector()
        if( self%azimuth == 4.0 ) self%lrows( 1, 1 ) = self%Lby%getFullVector()
        if( self%azimuth == 5.0 ) self%lrows( 1, 1 ) = self%Lbz%getFullVector()
        !
        if( self%azimuth == 3.0 .OR. self%azimuth == 4.0 .OR. self%azimuth == 5.0 ) then
            call self%lrows( 1, 1 )%mult( isign * comega )
        endif
        !
    end subroutine setLRowsSingleField
    !
    !> No subroutine briefing
    subroutine predictedDataSingleField( self, transmitter, data_group )
        implicit none
        !
        class( ReceiverSingleField_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in ) :: transmitter
        type( DataGroup_t ), intent( out ), optional :: data_group
        !
        complex( kind=prec ) :: comega
        !
        comega = cmplx( 0.0, 1./ ( 2.0 * PI / transmitter%period ), kind=prec )
        !
        select type( tx_e_1 => transmitter%e_sol(1) )
            !
            class is( cVector3D_SG_t )
                !
                if( allocated( self%response ) ) deallocate( self%response )
                allocate( self%response(1) )
                !
                select case ( self%EHxy(1)%str )
                    case( "Ex" )
                        self%response(1) = self%Lex%dotProd( tx_e_1 )
                    case( "Ey" )
                        self%response(1) = self%Ley%dotProd( tx_e_1 )
                    case( "Bx" )
                        self%response(1) = isign * self%Lbx%dotProd( tx_e_1 ) * comega
                    case( "By" )
                        self%response(1) = isign * self%Lbx%dotProd( tx_e_1 ) * comega
                    case( "Bz" )
                        self%response(1) = isign * self%Lbx%dotProd( tx_e_1 ) * comega
                end select
                !
                if( present( data_group ) ) then
                    !
                    call self%savePredictedData( transmitter, data_group )
                    !
                endif
                !
            class default
                stop "evaluationFunctionRx: Unclassified temp_full_vec_ey"
            !
        end select
        !
    end subroutine predictedDataSingleField
    !
    !> No function briefing
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
    !> No subroutine briefing
    subroutine printReceiverSingleField( self )
        implicit none
        !
        class( ReceiverSingleField_t ), intent( in ) :: self
        !
        write( *, * ) "Print ReceiverSingleField_t: ", self%i_rx
        !
    end subroutine printReceiverSingleField
    !
end module ReceiverSingleField
