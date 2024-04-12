!
!> Derived class to define a Full_Impedance Receiver
!
module ReceiverFullVerticalMagnetic
    !
    use Receiver
    !
    type, extends( Receiver_t ), public :: ReceiverFullVerticalMagnetic_t
        !
        !> No derived properties
        !
        contains
            !
            final :: ReceiverFullVerticalMagnetic_dtor
            !
            procedure, public :: predictedData => predictedData_FullVerticalMagnetic
            !
            procedure, public :: setLRows => setLRows_FullVerticalMagnetic
            !
            procedure, public :: isEqualRx => isEqual_FullVerticalMagnetic
            !
            procedure, public :: print => print_FullVerticalMagnetic
            !
    end type ReceiverFullVerticalMagnetic_t
    !
    interface ReceiverFullVerticalMagnetic_t
        module procedure ReceiverFullVerticalMagnetic_ctor
    end interface ReceiverFullVerticalMagnetic_t
    !
contains
    !
    !> No subroutine briefing
    !
    function ReceiverFullVerticalMagnetic_ctor( location, rx_type ) result( self )
        implicit none
        !
        real( kind=prec ), intent( in ) :: location(3)
        integer, intent( in ) :: rx_type
        !
        type( ReceiverFullVerticalMagnetic_t ) :: self
        !
        integer :: i, asize
        !
        !write( *, * ) "Constructor ReceiverFullVerticalMagnetic_t"
        !
        call self%baseInit
        !
        self%location = location
        !
        self%rx_type = rx_type
        !
        self%n_comp = 2
        self%is_complex = .TRUE.
        !
        !> components required to get the full impedance evaluation vectors [Bx, By, Bz]
        if( allocated( self%EHxy ) ) then
            !
            asize = size( self%EHxy )
            do i = asize, 1, -(1)
                deallocate( self%EHxy(i)%str )
            enddo
            deallocate( self%EHxy )
            !
        endif
        allocate( self%EHxy( 3 ) )
        !
        self%EHxy(1)%str = "Bx"
        self%EHxy(2)%str = "By"
        self%EHxy(3)%str = "Bz"
        !
        !> components required to get the full impedance tensor self%response [Tx, Ty]
        if( allocated( self%comp_names ) ) then
            !
            asize = size( self%comp_names )
            do i = asize, 1, -(1)
                deallocate( self%comp_names(i)%str )
            enddo
            deallocate( self%comp_names )
            !
        endif
        allocate( self%comp_names( 2 ) )
        !
        self%comp_names(1)%str = "TX"
        self%comp_names(2)%str = "TY"
        !
    end function ReceiverFullVerticalMagnetic_ctor
    !
    !> No subroutine briefing
    subroutine ReceiverFullVerticalMagnetic_dtor( self )
        implicit none
        !
        type( ReceiverFullVerticalMagnetic_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor ReceiverFullVerticalMagnetic_t"
        !
        call self%baseDealloc
        !
    end subroutine ReceiverFullVerticalMagnetic_dtor
    !
    !> No subroutine briefing
    subroutine setLRows_FullVerticalMagnetic( self, transmitter )
        implicit none
        !
        class( ReceiverFullVerticalMagnetic_t ), intent( inout ) :: self
        type( TransmitterMT_t ), intent( in ) :: transmitter
        !
        call errStop( "setLRows_FullVerticalMagnetic to be implemented" )
        !
        !if( allocated( self%lrows ) ) deallocate( self%lrows )
        !allocate( self%lrows( transmitter%n_pol, self%n_comp ) )
        !
    end subroutine setLRows_FullVerticalMagnetic
    !
    !> No subroutine briefing
    subroutine predictedData_FullVerticalMagnetic( self, transmitter, data_group )
        implicit none
        !
        class( ReceiverFullVerticalMagnetic_t ), intent( inout ) :: self
        type( TransmitterMT_t ), intent( in ) :: transmitter
        type( DataGroup_t ), intent( out ), optional :: data_group
        !
        complex( kind=prec ) :: comega, det
        complex( kind=prec ), allocatable :: BB(:,:), I_BB(:,:)
        class( Vector_t ), allocatable :: tx_e_1, tx_e_2
        !
        call transmitter%getSolution( 1, tx_e_1 )
        call transmitter%getSolution( 2, tx_e_2 )
        !
        comega = cmplx( 0.0, 1./ ( 2.0 * PI / transmitter%period ), kind=prec )
        !
        allocate( BB( 3, 2 ) )
        !
        BB(1,1) = self%Lbx%dotProd( tx_e_1 )
        BB(2,1) = self%Lby%dotProd( tx_e_1 )
        BB(1,2) = self%Lbx%dotProd( tx_e_2 )
        BB(2,2) = self%Lby%dotProd( tx_e_2 )
        BB(3,1) = self%Lbz%dotProd( tx_e_1 )
        BB(3,2) = self%Lbz%dotProd( tx_e_2 )
        !
        deallocate( tx_e_1, tx_e_2 )
        !
        BB = isign * BB * comega
        !
        det = BB(1,1) * BB(2,2) - BB(1,2) * BB(2,1)
        !
        allocate( I_BB(2,2) )
        !
        if( det /= 0 ) then
            I_BB(1,1) =  BB(2,2) / det
            I_BB(2,2) =  BB(1,1) / det
            I_BB(1,2) = -BB(1,2) / det
            I_BB(2,1) = -BB(2,1) / det
        else
            call errStop( "predictedData_FullVerticalMagnetic > Determinant is Zero!" )
        endif
        !
        allocate( self%response(2) )
        !
        self%response(1) = BB(3,1) * I_BB(1,1) + BB(3,2) * I_BB(2,1)
        self%response(2) = BB(3,1) * I_BB(1,2) + BB(3,2) * I_BB(2,2)
        !
        deallocate( BB )
        deallocate( I_BB )
        !
        if( present( data_group ) ) then
            !
            call self%savePredictedData( transmitter, data_group )
            !
        endif
        !
    end subroutine predictedData_FullVerticalMagnetic
    !
    !> No subroutine briefing
    !
    function isEqual_FullVerticalMagnetic( self, other ) result( equal )
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
    end function isEqual_FullVerticalMagnetic
    !
    !> No subroutine briefing
    subroutine print_FullVerticalMagnetic( self )
        implicit none
        !
        class( ReceiverFullVerticalMagnetic_t ), intent( in ) :: self
        !
        write( *, * ) "Print ReceiverFullVerticalMagnetic_t: ", self%i_rx
        !
    end subroutine print_FullVerticalMagnetic
    !
end module ReceiverFullVerticalMagnetic
