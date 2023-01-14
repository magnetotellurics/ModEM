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
            procedure, public :: setLRows => setLRowsFullVerticalMagnetic
            !
            procedure, public :: predictedData => predictedDataFullVerticalMagnetic
            !
            procedure, public :: isEqualRx => isEqualFullVerticalMagnetic
            !
            procedure, public :: print => printReceiverFullVerticalMagnetic
            !
    end type ReceiverFullVerticalMagnetic_t
    !
    interface ReceiverFullVerticalMagnetic_t
        module procedure ReceiverFullVerticalMagnetic_ctor
    end interface ReceiverFullVerticalMagnetic_t
    !
contains
    !
    !> No function briefing
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
        call self%init()
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
        call self%dealloc()
        !
    end subroutine ReceiverFullVerticalMagnetic_dtor
    !
    !> No subroutine briefing
    subroutine setLRowsFullVerticalMagnetic( self, transmitter, lrows )
        implicit none
        !
        class( ReceiverFullVerticalMagnetic_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in ) :: transmitter
        class( Vector_t ), allocatable, dimension(:,:), intent( out ) :: lrows
        !
        stop "setLRowsFullVerticalMagnetic to be implemented"
        !
    end subroutine setLRowsFullVerticalMagnetic
    !
    !> No subroutine briefing
    subroutine predictedDataFullVerticalMagnetic( self, transmitter, data_group )
        implicit none
        !
        class( ReceiverFullVerticalMagnetic_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in ) :: transmitter
        type( DataGroup_t ), intent( out ), optional :: data_group
        !
        complex( kind=prec ) :: comega, det
        complex( kind=prec ), allocatable :: BB(:,:), I_BB(:,:)
        !
        comega = cmplx( 0.0, 1./ ( 2.0 * PI / transmitter%period ), kind=prec )
        !
        allocate( BB( 3, 2 ) )
        select type( tx_e_1 => transmitter%e_sol(1) )
            class is( cVector3D_SG_t )
                !
                select type( tx_e_2 => transmitter%e_sol(2) )
                    class is( cVector3D_SG_t )
                        !
                        BB(1,1) = self%Lbx%dotProd( tx_e_1 )
                        BB(2,1) = self%Lby%dotProd( tx_e_1 )
                        BB(1,2) = self%Lbx%dotProd( tx_e_2 )
                        BB(2,2) = self%Lby%dotProd( tx_e_2 )
                        BB(3,1) = self%Lbz%dotProd( tx_e_1 )
                        BB(3,2) = self%Lbz%dotProd( tx_e_2 )
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
                            stop "Error: ReceiverFullVerticalMagnetic.f90: Determinant is Zero!"
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
                    class default
                        stop "Error: evaluationFunctionRx: Unclassified transmitter%e_all_2"
                end select
                !
            class default
                stop "Error: evaluationFunctionRx: Unclassified transmitter%e_all_1"
        end select
        !
    end subroutine predictedDataFullVerticalMagnetic
    !
    !> No function briefing
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
    !> No subroutine briefing
    subroutine printReceiverFullVerticalMagnetic( self )
        implicit none
        !
        class( ReceiverFullVerticalMagnetic_t ), intent( in ) :: self
        !
        write( *, * ) "Print ReceiverFullVerticalMagnetic_t: ", self%i_rx
        !
    end subroutine printReceiverFullVerticalMagnetic
    !
end module ReceiverFullVerticalMagnetic
