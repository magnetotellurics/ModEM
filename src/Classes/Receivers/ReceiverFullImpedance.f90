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
            procedure, public :: setLRows => setLRowsFullImpedance
            !
            procedure, public :: predictedData => predictedDataFullImpedance
            !
            procedure, public :: isEqualRx => isEqualFullImpedance
            !
            procedure, public :: print => printReceiverFullImpedance
            !
    end type ReceiverFullImpedance_t
    !
    interface ReceiverFullImpedance_t
        module procedure ReceiverFullImpedance_ctor
    end interface ReceiverFullImpedance_t
    !
contains
    !
    !> No function briefing
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
        call self%init()
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
        call self%dealloc()
        !
    end subroutine ReceiverFullImpedance_dtor
    !
    !> No subroutine briefing
    subroutine setLRowsFullImpedance( self, transmitter )
        implicit none
        !
        class( ReceiverFullImpedance_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in ) :: transmitter
        !
        type( cVector3D_SG_t ) :: Le, full_lex, full_ley, full_lbx, full_lby
        integer :: Ei, row, pol, comp
        complex( kind=prec ) :: comega
        !
        comega = cmplx( 0.0, 1. / transmitter%omega, kind=prec )
        !
        !> Call the predicted data routine to calculate responses
        call self%predictedData( transmitter )
        !
        !> Allocate LRows matrix [ n_pol = 2, n_comp = 4 ]
        if( allocated( self%lrows ) ) deallocate( self%lrows )
        allocate( cVector3D_SG_t :: self%lrows( transmitter%n_pol, self%n_comp ) )
        !
        !> Convert Le and Lb to Full Vectors (In the future they will be Sparse)
        full_lex = self%Lex%getFullVector()
        full_ley = self%Ley%getFullVector()
        !
        full_lbx = self%Lbx%getFullVector()
        full_lby = self%Lby%getFullVector()
        !
        !> 
        !> Lrows{j,ki} = Hinv(j,i) * ( lE - Z(k,1) * 1/omega * Rx.Lhx - Z(k,2) * 1/omega * Rx.Lhy )
        !>
        !
        !> Loop over Ex and Ey
        do Ei = 1, 2
            !
            if( Ei == 1 ) then
                Le = full_lex
            else
                Le = full_ley
            endif
            !
            !> ????
            call Le%multAdd( C_MinusOne * self%response( 2 * (Ei-1) + 1 ) * comega, full_lbx )    ! 1 & 3
            !
            call Le%multAdd( C_MinusOne * self%response( 2 * Ei ) * comega, full_lby )            ! 2 & 4
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
                    self%lrows( pol, comp ) = Le
                    !
                    !> ????
                    call self%lrows( pol, comp )%mult( -self%I_BB( pol, row ) )
                    !
                    !call self%lrows( pol, comp )%print( 1000, "OO LRows" )
                    !
                enddo
                !
            enddo
            !
        enddo
        !
        deallocate( self%I_BB, self%response )
        !
    end subroutine setLRowsFullImpedance
    !
    !> No subroutine briefing
    subroutine predictedDataFullImpedance( self, transmitter )
        implicit none
        !
        class( ReceiverFullImpedance_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in ) :: transmitter
        !
        integer :: i, j, ij
        complex( kind=prec ) :: comega, det
        complex( kind=prec ), allocatable :: BB(:,:), EE(:,:)
        !
        comega = cmplx( 0.0, 1./ transmitter%omega, kind=prec )
        !
        allocate( EE(2,2) )
        !
        select type( tx_e_1 => transmitter%e_sol(1) )
            class is( cVector3D_SG_t )
                !
                select type( tx_e_2 => transmitter%e_sol(2) )
                    class is( cVector3D_SG_t )
                        !
                        EE(1,1) = self%Lex%dotProd( tx_e_1 )
                        EE(2,1) = self%Ley%dotProd( tx_e_1 )
                        EE(1,2) = self%Lex%dotProd( tx_e_2 )
                        EE(2,2) = self%Ley%dotProd( tx_e_2 )
                        !
                        allocate( BB( 2, 2 ) )
                        !
                        BB(1,1) = self%Lbx%dotProd( tx_e_1 )
                        BB(2,1) = self%Lby%dotProd( tx_e_1 )
                        BB(1,2) = self%Lbx%dotProd( tx_e_2 )
                        BB(2,2) = self%Lby%dotProd( tx_e_2 )
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
                            stop "Error: predictedDataFullImpedance > Determinant is Zero!"
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
                        call self%savePredictedData( transmitter )
                        !
                    class default
                        stop "Error: predictedDataFullImpedance: Unclassified tx_e_2"
                end select
                !
            class default
                stop "Error: predictedDataFullImpedance: Unclassified tx_e_1"
        end select
        !
    end subroutine predictedDataFullImpedance
    !
    !> No function briefing
    function isEqualFullImpedance( self, other ) result( equal )
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
    end function isEqualFullImpedance
    !
    !> No subroutine briefing
    subroutine printReceiverFullImpedance( self )
        implicit none
        !
        class( ReceiverFullImpedance_t ), intent( in ) :: self
        !
        write( *, * ) "ReceiverFullImpedance_t: ", self%id, self%rx_type, self%n_comp
        !
    end subroutine printReceiverFullImpedance
    !
end module ReceiverFullImpedance
