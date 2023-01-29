!
!> Derived class to define a Full_Impedance Receiver
!
module ReceiverOffDiagonalImpedance
    !
    use Receiver
    !
    type, extends( Receiver_t ), public :: ReceiverOffDiagonalImpedance_t
        !
        !> No derived properties
        !
        contains
            !
            final :: ReceiverOffDiagonalImpedance_dtor
            !
            procedure, public :: setLRows => setLRowsOffDiagonalImpedance
            !
            procedure, public :: predictedData => predictedDataOffDiagonalImpedance
            !
            procedure, public :: isEqualRx => isEqualOffDiagonalImpedance
            !
            procedure, public :: print => printReceiverOffDiagonalImpedance
            !
    end type ReceiverOffDiagonalImpedance_t
    !
    interface ReceiverOffDiagonalImpedance_t
        module procedure ReceiverOffDiagonalImpedance_ctor
    end interface ReceiverOffDiagonalImpedance_t
    !
contains
    !
    !> No subroutine briefing
	!
    function ReceiverOffDiagonalImpedance_ctor( location, rx_type ) result( self )
        implicit none
        !
        real( kind=prec ), intent( in ) :: location(3)
        integer, optional, intent( in ) :: rx_type
        !
        type( ReceiverOffDiagonalImpedance_t ) :: self
        !
        integer :: i, asize
        !
        !> write( *, * ) "Constructor ReceiverOffDiagonalImpedance_t"
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
        allocate( self%EHxy( 4 ) )
        !
        self%EHxy(1)%str = "Ex"
        self%EHxy(2)%str = "Ey"
        self%EHxy(3)%str = "Bx"
        self%EHxy(4)%str = "By"
        !
    end function ReceiverOffDiagonalImpedance_ctor
    !
    !> No subroutine briefing
    subroutine ReceiverOffDiagonalImpedance_dtor( self )
        implicit none
        !
        type( ReceiverOffDiagonalImpedance_t ), intent( inout ) :: self
        !
        !> write( *, * ) "Destructor ReceiverOffDiagonalImpedance_t"
        !
        call self%dealloc()
        !
    end subroutine ReceiverOffDiagonalImpedance_dtor
    !
    !> No subroutine briefing
    subroutine setLRowsOffDiagonalImpedance( self, transmitter )
        implicit none
        !
        class( ReceiverOffDiagonalImpedance_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in ) :: transmitter
        !
        stop "setLRowsOffDiagonalImpedance to be implemented"
        !
        if( allocated( self%lrows ) ) deallocate( self%lrows )
        allocate( cVector3D_SG_t :: self%lrows( transmitter%n_pol, self%n_comp ) )
        !
    end subroutine setLRowsOffDiagonalImpedance
    !
    !> No subroutine briefing
    subroutine predictedDataOffDiagonalImpedance( self, transmitter, data_group )
        implicit none
        !
        class( ReceiverOffDiagonalImpedance_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in ) :: transmitter
        type( DataGroup_t ), intent( out ), optional :: data_group
        !
        integer :: i, j, ij
        complex( kind=prec ) :: comega, det
        complex( kind=prec ), allocatable :: BB(:,:), I_BB(:,:), EE(:,:)
        !
        comega = cmplx( 0.0, 1./ ( 2.0 * PI / transmitter%period ), kind=prec )
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
                        allocate( BB(2,2) )
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
                        allocate( I_BB(2,2) )
                        !
                        if( det /= 0 ) then
                            I_BB(1,1) =  BB(2,2) / det
                            I_BB(2,2) =  BB(1,1) / det
                            I_BB(1,2) = -BB(1,2) / det
                            I_BB(2,1) = -BB(2,1) / det
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
                        if( present( data_group ) ) then
                            !
                            call self%savePredictedData( transmitter, data_group )
                            !
                        endif
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
    !> No subroutine briefing
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
    !> No subroutine briefing
    subroutine printReceiverOffDiagonalImpedance( self )
        implicit none
        !
        class( ReceiverOffDiagonalImpedance_t ), intent( in ) :: self
        !
        write( *, * ) "Print ReceiverOffDiagonalImpedance_t: ", self%i_rx
        !
    end subroutine printReceiverOffDiagonalImpedance
    !
end module ReceiverOffDiagonalImpedance
