! *************
! 
! Derived class to define a Full_Impedance Receiver
! 
! *************
! 
module ReceiverFullImpedance
    !
    use Receiver
    !
    type, extends( Receiver_t ), public :: ReceiverFullImpedance_t
        !
        ! Specific properties here
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
    function ReceiverFullImpedance_ctor( location, rx_type ) result( self )
        implicit none
        !
        real( kind=prec ), intent( in ) :: location(3)
        integer, intent( in )           :: rx_type
        !
        type( ReceiverFullImpedance_t ) :: self
        !
        integer :: i, asize
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
        self%is_complex = .TRUE.
        !
        ! components required to get the full impedance evaluation vectors [Ex, Ey, Bx, By]
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
        allocate( self%Lex, source = cSparsevector3D_SG_t() )
        allocate( self%Ley, source = cSparsevector3D_SG_t() )
        allocate( self%Lbx, source = cSparsevector3D_SG_t() )
        allocate( self%Lby, source = cSparsevector3D_SG_t() )
        !
        ! components required to get the full impedance tensor self%response [Zxx, Zxy, Zyx, Zyy]
        if( allocated( self%comp_names ) ) then
            !
            asize = size( self%comp_names )
            do i = asize, 1, -(1)
                deallocate( self%comp_names(i)%str )
            enddo
            deallocate( self%comp_names )
            !
        endif
        allocate( self%comp_names( 4 ) )
        !
        self%comp_names(1)%str = "ZXX"
        self%comp_names(2)%str = "ZXY"
        self%comp_names(3)%str = "ZYX"
        self%comp_names(4)%str = "ZYY"
        !
    end function ReceiverFullImpedance_ctor
    !
    subroutine ReceiverFullImpedance_dtor( self )
        implicit none
        !
        type( ReceiverFullImpedance_t ), intent( in out ) :: self
        !
        !write( *, * ) "Destructor ReceiverFullImpedance_t:", self%id
        !
        call self%dealloc()
        !
    end subroutine ReceiverFullImpedance_dtor
    !
    subroutine setLRowsFullImpedance( self, transmitter )
        implicit none
        !
        class( ReceiverFullImpedance_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in )              :: transmitter
        !
        type( cVector3D_SG_t ) :: Le, full_lex, full_ley, full_lbx, full_lby, aux_full_vec
        type( cSparsevector3D_SG_t ) :: aux_sparse_vec
        integer :: i, j, k, ki, kj
        !
        write( *, * ) "implementing setLRowsFullImpedance:"
        !
        allocate( self%lrows( transmitter%n_pol, self%n_comp ) )
        !
        ki = 0
        !
        full_lex = self%Lex%getFullVector()
        full_ley = self%Ley%getFullVector()
        !
        full_lbx = self%Lbx%getFullVector()
        full_lby = self%Lby%getFullVector()
        !
        do k = 1, 2
            !
            if( k == 1 ) then
                Le = self%Lex%getFullVector()
            else
                Le = self%Ley%getFullVector()
            endif
            !
            do i = 1, 2
                ki = ki + 1
                do j = 1, 2
                    !
                    kj = 2 * ( k-1 ) + j
                    !
                    call full_lbx%mult( self%response( kj ) )
                    call full_lby%mult( self%response( kj ) )
                    !
                    aux_full_vec = ( Le - full_lbx - full_lby )
                    !
                    call aux_full_vec%mult( self%I_BB( j, i ) )
                    !
                    aux_sparse_vec = cSparsevector3D_SG_t()
                    call aux_sparse_vec%fromFullVector( aux_full_vec )
                    !
                    self%lrows( j, ki ) = aux_sparse_vec
                enddo
            enddo
            !
        enddo
        !
    end subroutine setLRowsFullImpedance
    !
    subroutine predictedDataFullImpedance( self, transmitter )
        implicit none
        !
        class( ReceiverFullImpedance_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in )              :: transmitter
        !
        integer :: i, j, ij
        complex( kind=prec ) :: comega, det
        complex( kind=prec ), allocatable :: BB(:,:), EE(:,:)
        !
        !
        comega = cmplx( 0.0, 1./ ( 2.0 * PI / transmitter%period ), kind=prec )
        !
        allocate( EE(2,2) )
        !
        select type( tx_e_1 => transmitter%e_all(1) )
            class is( cVector3D_SG_t )
                !
                select type( tx_e_2 => transmitter%e_all(2) )
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
                            STOP "ReceiverFullImpedance.f90: Determinant is Zero!"
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
                        deallocate( self%response )
                        !
                    class default
                        stop "evaluationFunctionRx: Unclassified temp_full_vec_ey"
                end select
                !
            class default
                stop "evaluationFunctionRx: Unclassified temp_full_vec_ey"
        end select
        !
    end subroutine predictedDataFullImpedance
    !
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
    subroutine printReceiverFullImpedance( self )
        implicit none
        !
        class( ReceiverFullImpedance_t ), intent( in ) :: self
        !
        write( *, * ) "Print ReceiverFullImpedance_t: ", self%id
        !
    end subroutine printReceiverFullImpedance
    !
end module ReceiverFullImpedance
