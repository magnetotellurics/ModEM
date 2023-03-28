!
!> Derived class to define a Single E or B Field Receiver
!
module ReceiverExyAmpliPhase
    !
    use Receiver
    !
    type, extends( Receiver_t ), public :: ReceiverExyAmpliPhase_t
        !
        real( kind=prec ) :: azimuth
        !
        contains
            !
            final :: ReceiverExyAmpliPhase_dtor
            !
            procedure, public :: setLRows => setLRowsExyAmpliPhase
            !
            procedure, public :: predictedData => predictedDataExyAmpliPhase
            !
            procedure, public :: isEqualRx => isEqualExyAmpliPhase
            !
            procedure, public :: print => printReceiverExyAmpliPhase
            !
    end type ReceiverExyAmpliPhase_t
    !
    interface ReceiverExyAmpliPhase_t
        module procedure ReceiverExyAmpliPhase_ctor
    end interface ReceiverExyAmpliPhase_t
    !
contains
    !
    !> No subroutine briefing
    !
    function ReceiverExyAmpliPhase_ctor( location, azimuth, rx_type ) result( self )
        implicit none
        !
        real( kind=prec ), intent( in ) :: location(3)
        real( kind=prec ) :: azimuth
        integer, intent( in ) :: rx_type
        !
        type( ReceiverExyAmpliPhase_t ) :: self
        !
        integer :: i, asize
        !
        !> write( *, * ) "Constructor ReceiverExyAmpliPhase_t"
        !
        call self%init
        !
        self%location = location
        !
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
        allocate( self%EHxy(2) )
        !
        self%EHxy(1)%str = "Ex"
        self%EHxy(2)%str = "Ey"
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
        !
        allocate( self%comp_names(1) )
        !
        self%comp_names(1)%str = "Exy_Ampli"
        !
    end function ReceiverExyAmpliPhase_ctor
    !
    !> No subroutine briefing
    subroutine ReceiverExyAmpliPhase_dtor( self )
        implicit none
        !
        type( ReceiverExyAmpliPhase_t ), intent( inout ) :: self
        !
        !> write( *, * ) "Destructor ReceiverExyAmpliPhase_t"
        !
        call self%dealloc
        !
    end subroutine ReceiverExyAmpliPhase_dtor
    !
    !> No subroutine briefing
    subroutine setLRowsExyAmpliPhase( self, transmitter )
        implicit none
        !
        class( ReceiverExyAmpliPhase_t ), intent( inout ) :: self
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
        self%lrows( 1, 1 ) = self%Lex%getFullVector()
        self%lrows( 1, 2 ) = self%Ley%getFullVector()
        !
    end subroutine setLRowsExyAmpliPhase
    !
    !> No subroutine briefing
    subroutine predictedDataExyAmpliPhase( self, transmitter, data_group )
        implicit none
        !
        class( ReceiverExyAmpliPhase_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in ) :: transmitter
        type( DataGroup_t ), intent( out ), optional :: data_group
        !
        complex( kind=prec ) :: comega
        complex( kind=prec )    :: Ex_res, Ey_res, tempZ
        class( Vector_t ), pointer :: tx_e_1
        !
        comega = cmplx( 0.0, 1./ ( 2.0 * PI / transmitter%period ), kind=prec )
        !
        call transmitter%getSolutionVector( 1, tx_e_1 )
        !
        select type( tx_e_1 )
            !
            class is( cVector3D_SG_t )
                !
                if( allocated( self%response ) ) deallocate( self%response )
                allocate( self%response(2) )
                !
                Ex_res = self%Lex%dotProd( tx_e_1 )
                Ey_res = self%Ley%dotProd( tx_e_1 )
                !
                tempZ = ( ( ( cos( D2R * self%azimuth ) ) * Ex_res ) + &
                          ( ( sin( D2R * self%azimuth ) ) * Ey_res ) )
                !
                self%response(1) = log10( abs( tempZ ) )
                self%response(2) = atan2( isign * aimag( tempZ ), real( tempZ, kind=prec ) ) * R2D
                !
                if( present( data_group ) ) then
                    !
                    call self%savePredictedData( transmitter, data_group )
                    !
                endif
                !
            class default
                stop "predictedDataExyAmpliPhase: Unclassified tx_e_1"
            !
        end select
        !
    end subroutine predictedDataExyAmpliPhase
    !
    !> No subroutine briefing
    !
    function isEqualExyAmpliPhase( self, other ) result( equal )
        implicit none
        !
        class( ReceiverExyAmpliPhase_t ), intent( in ) :: self
        class( Receiver_t ), intent( in ) :: other
        !
        logical :: equal
        !
        equal = .FALSE.
        !
        select type( other )
            !
            class is( ReceiverExyAmpliPhase_t )
                !
                if( self%code == other%code .AND.   &
                    self%azimuth == other%azimuth .AND.    &
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
    end function isEqualExyAmpliPhase
    !
    !> No subroutine briefing
    subroutine printReceiverExyAmpliPhase( self )
        implicit none
        !
        class( ReceiverExyAmpliPhase_t ), intent( in ) :: self
        !
        write( *, * ) "Print ReceiverExyAmpliPhase_t: ", self%i_rx
        !
    end subroutine printReceiverExyAmpliPhase
    !
end module ReceiverExyAmpliPhase
