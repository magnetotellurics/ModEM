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
    function ReceiverOffDiagonalImpedance_ctor( location, rx_type ) result( self )
        implicit none
        !
        real( kind=prec ), intent( in ) :: location(3)
        integer, optional, intent( in ) :: rx_type
        !
        type( ReceiverOffDiagonalImpedance_t ) :: self
        !
        character(:), allocatable :: aux_str
        !
        ! write(*,*) "Constructor ReceiverOffDiagonalImpedance_t"
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
    subroutine predictedDataOffDiagonalImpedance( self, transmitter )
        implicit none
        !
        class( ReceiverOffDiagonalImpedance_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in )                     :: transmitter
        !
        complex( kind=prec ) :: comega
        !
        complex( kind=prec ), allocatable :: BB(:,:), det
        real( kind=prec ) :: omega
        integer           :: i, j, ij
        !
        comega = cmplx( 0.0, 1./ ( 2.0 * PI / transmitter%period ), kind=prec )
        !
        !
        allocate( complex( kind=prec ) :: self%EE( 2, 2 ) )
        !
        select type( tx_e_1 => transmitter%e_all( 1 ) )
            class is( cVector3D_SG_t )
                !
                select type( tx_e_2 => transmitter%e_all( 2 ) )
                    class is( cVector3D_SG_t )
                        !
                        self%EE(1,1) = dotProdSparse( self%Lex, tx_e_1 )
                        self%EE(2,1) = dotProdSparse( self%Ley, tx_e_1 )
                        self%EE(1,2) = dotProdSparse( self%Lex, tx_e_2 )
                        self%EE(2,2) = dotProdSparse( self%Ley, tx_e_2 )
                        !
                        allocate( complex( kind=prec ) :: BB( 2, 2 ) )
                        !
                        BB(1,1) = dotProdSparse( self%Lbx, tx_e_1 )
                        BB(2,1) = dotProdSparse( self%Lby, tx_e_1 )
                        BB(1,2) = dotProdSparse( self%Lbx, tx_e_2 )
                        BB(2,2) = dotProdSparse( self%Lby, tx_e_2 )
                        BB = isign * BB * comega
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
