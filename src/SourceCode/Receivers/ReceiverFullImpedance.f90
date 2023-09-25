!
!> Derived class to define a Full_Impedance Receiver
!
module ReceiverFullImpedance
    !
    use Receiver
    !
    type, extends( Receiver_t ), public :: ReceiverFullImpedance_t
        !
        !> No derived properties
        !
        contains
            !
            final :: ReceiverFullImpedance_dtor
            !
            procedure, public :: predictedData => predictedData_FullImpedance
            !
            procedure, public :: setLRows => setLRows_FullImpedance
            !
            procedure, public :: isEqualRx => isEqual_FullImpedance
            !
            procedure, public :: print => print_FullImpedance
            !
    end type ReceiverFullImpedance_t
    !
    interface ReceiverFullImpedance_t
        module procedure ReceiverFullImpedance_ctor
    end interface ReceiverFullImpedance_t
    !
contains
    !
    !> No subroutine briefing
    !
    function ReceiverFullImpedance_ctor( location, rx_type ) result( self )
        implicit none
        !
        real( kind=prec ), intent( in ) :: location(3)
        integer, intent( in ) :: rx_type
        !
        type( ReceiverFullImpedance_t ) :: self
        !
        !write( *, * ) "Constructor ReceiverFullImpedance_t"
        !
        call self%baseInit
        !
        self%location = location
        !
        self%rx_type = rx_type
        !
        self%n_comp = 4
        !
        self%is_complex = .TRUE.
        !
        !> components required to get the full impedance evaluation vectors [Ex, Ey, Bx, By]
        allocate( self%EHxy( 4 ) )
        !
        self%EHxy(1)%str = "Ex"
        self%EHxy(2)%str = "Ey"
        self%EHxy(3)%str = "Bx"
        self%EHxy(4)%str = "By"
        !
        allocate( self%comp_names( 4 ) )
        !
        self%comp_names(1)%str = "ZXX"
        self%comp_names(2)%str = "ZXY"
        self%comp_names(3)%str = "ZYX"
        self%comp_names(4)%str = "ZYY"
        !
    end function ReceiverFullImpedance_ctor
    !
    !> No subroutine briefing
    subroutine ReceiverFullImpedance_dtor( self )
        implicit none
        !
        type( ReceiverFullImpedance_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor ReceiverFullImpedance_t:", self%id
        !
        call self%baseDealloc
        !
    end subroutine ReceiverFullImpedance_dtor
    !
    !> No subroutine briefing
     !
    subroutine predictedData_FullImpedance( self, transmitter, data_group )
        implicit none
        !
        class( ReceiverFullImpedance_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in ) :: transmitter
        type( DataGroup_t ), intent( out ), optional :: data_group
        !
        integer :: i, j, ij
        complex( kind=prec ) :: comega, det
        complex( kind=prec ), allocatable :: BB(:,:), EE(:,:)
        class( Vector_t ), allocatable :: tx_e_1, tx_e_2
        !
        call transmitter%getSolutionVector( 1, tx_e_1 )
        call transmitter%getSolutionVector( 2, tx_e_2 )
        !
        !call tx_e_1%print( 6666 )
        !call tx_e_2%print( 6667 )
        !
        comega = cmplx( 0.0, 1. / ( 2.0 * PI / transmitter%period ), kind=prec )
        !
        allocate( EE(2,2) )
        !
        EE(1,1) = self%Lex%dotProd( tx_e_1 )
        !
        EE(2,1) = self%Ley%dotProd( tx_e_1 )
        !
        EE(1,2) = self%Lex%dotProd( tx_e_2 )
        !
        EE(2,2) = self%Ley%dotProd( tx_e_2 )
        !
        !write( *, * ) "EE"
        !write( *, * ) EE(1,1), EE(1,2)
        !write( *, * ) EE(2,1), EE(2,2)
        !
        allocate( BB( 2, 2 ) )
        BB(1,1) = self%Lbx%dotProd( tx_e_1 )
        BB(2,1) = self%Lby%dotProd( tx_e_1 )
        BB(1,2) = self%Lbx%dotProd( tx_e_2 )
        BB(2,2) = self%Lby%dotProd( tx_e_2 )
        !
        deallocate( tx_e_1, tx_e_2 )
        !
        !write( *, * ) "BB"
        !write( *, * ) BB(1,1), BB(1,2)
        !write( *, * ) BB(2,1), BB(2,2)
        !
        BB = isign * BB * comega
        !
        det = BB(1,1) * BB(2,2) - BB(1,2) * BB(2,1)
        !
        if( allocated( self%I_BB ) ) deallocate( self%I_BB )
        allocate( self%I_BB(2,2) )
        !
        if( det /= 0 ) then
            self%I_BB(1,1) =  BB(2,2) / det
            self%I_BB(2,2) =  BB(1,1) / det
            self%I_BB(1,2) = -BB(1,2) / det
            self%I_BB(2,1) = -BB(2,1) / det
        else
            call errStop( "predictedData_FullImpedance > Determinant is Zero!" )
        endif
        !
        deallocate( BB )
        !
        if( allocated( self%response ) ) deallocate( self%response )
        allocate( self%response(4) )
        !
        do j = 1, 2
             do i = 1, 2
                 ij = 2 * ( i-1 ) + j
                 self%response(ij) = EE(i,1) * self%I_BB(1,j) + EE(i,2) * self%I_BB(2,j)
             enddo
        enddo
        !
        deallocate( EE )
        !
        if( present( data_group ) ) then
            !
            call self%savePredictedData( transmitter, data_group )
            !
        endif
        !
    end subroutine predictedData_FullImpedance
    !
    !> No subroutine briefing
    subroutine setLRows_FullImpedance( self, transmitter )
        implicit none
        !
        class( ReceiverFullImpedance_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in ) :: transmitter
        !
        class( Vector_t ), allocatable :: Le, full_lex, full_ley, full_lbx, full_lby
        integer :: Ei, row, pol, comp
        complex( kind=prec ) :: comega
        !
        comega = isign * cmplx( 0.0, 1. / ( 2.0 * PI / transmitter%period ), kind=prec )
        !
        !> Call the predicted data routine to calculate responses
        call self%predictedData( transmitter )
        !
        !> Convert Le and Lb to Full Vectors (In the future they will be Sparse)
        allocate( full_lex, source = self%Lex%getFullVector() )
        allocate( full_ley, source = self%Ley%getFullVector() )
        !
        allocate( full_lbx, source = self%Lbx%getFullVector() )
        allocate( full_lby, source = self%Lby%getFullVector() )
        !
        !> Lrows{j,ki} = Hinv(j,i) * ( lE - Z(k,1) * 1/omega * Rx.Lhx - Z(k,2) * 1/omega * Rx.Lhy )
        !
        !> Loop over Ex and Ey
        do Ei = 1, 2
            !
            if( Ei == 1 ) then
                allocate( Le, source = full_lex )
            else
                allocate( Le, source = full_ley )
            endif
            !
            !> ????
            call Le%multAdd( -self%response( 2 * (Ei-1) + 1 ) * comega, full_lbx ) ! 1 & 3
            !
            call Le%multAdd( -self%response( 2 * Ei ) * comega, full_lby )         ! 2 & 4
            !
            !> Loop over two impedance rows
            do row = 1, 2
                !
                !> comp = 1, 3, 2, 4
                comp = 2 * (Ei-1) + row
                !
                !> Loop over two polarizations
                do pol = 1, 2
                    !
                    self%lrows( pol, comp )%v = Le
                    !
                    !> ????
                    call self%lrows( pol, comp )%v%mult( -self%I_BB( pol, row ) )
                    !
                enddo
                !
            enddo
            !
            deallocate( Le )
            !
        enddo
        !
        deallocate( full_lex, full_ley )
        !
        deallocate( full_lbx, full_lby )
        !
        !deallocate( self%I_BB, self%response )
        !
    end subroutine setLRows_FullImpedance
     !
    !> No subroutine briefing
    !
    function isEqual_FullImpedance( self, other ) result( equal )
        implicit none
        !
        class( ReceiverFullImpedance_t ), intent( in ) :: self
        class( Receiver_t ), intent( in ) :: other
        !
        logical :: equal
        !
        equal = .FALSE.
        !
        select type( other )
            !
            class is( ReceiverFullImpedance_t )
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
    end function isEqual_FullImpedance
    !
    !> No subroutine briefing
    subroutine print_FullImpedance( self )
        implicit none
        !
        class( ReceiverFullImpedance_t ), intent( in ) :: self
        !
        write( *, * ) "ReceiverFullImpedance_t: ", self%i_rx, self%rx_type, self%n_comp
        !
    end subroutine print_FullImpedance
    !
end module ReceiverFullImpedance
