! *************
! 
! Derived class to define a Full_Impedance Receiver
! 
! Last modified at 08/06/2021 by Paulo Werdt
! 
! *************
! 
module ReceiverFullVerticalMagnetic
    !
    use Receiver
    use DataHandleMT
    !
    type, extends( Receiver_t ), public :: ReceiverFullVerticalMagnetic_t
        !
        ! PROPERTIES HERE
        !
        contains
            !
            final :: ReceiverFullVerticalMagnetic_dtor
            !
            procedure, public :: isEqualRx => isEqualFullVerticalMagnetic
            !
            procedure, public :: predictedData => predictedDataFullVerticalMagnetic
            !
            procedure, public :: savePredictedData => savePredictedDataFullVerticalMagnetic
            !
            procedure, public :: write => writeReceiverFullVerticalMagnetic
            !
    end type ReceiverFullVerticalMagnetic_t
    !
    interface ReceiverFullVerticalMagnetic_t
        module procedure ReceiverFullVerticalMagnetic_ctor
    end interface ReceiverFullVerticalMagnetic_t
    !
contains
    !
    function ReceiverFullVerticalMagnetic_ctor( location, rx_type ) result( self )
        implicit none
        !
        real( kind=prec ), intent( in ) :: location(3)
        integer, intent( in )           :: rx_type
        !
        type( ReceiverFullVerticalMagnetic_t ) :: self
        !
        character(:), allocatable :: aux_str
        !
        !write(*,*) "Constructor ReceiverFullVerticalMagnetic_t"
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
        allocate( self%EHxy( 3 ) )
        !
        ! components required to get the full impdence tensor self%response [Zxx, Zxy, Zyx, Zyy]
        self%EHxy(1)%str = "Bx"
        self%EHxy(2)%str = "By"
        self%EHxy(3)%str = "Bz"
        !
        ! components required to get the full impdedance tensor self%response [Zxx, Zxy, Zyx, Zyy]
        allocate( self%comp_names( 2 ) )
        !
        self%comp_names(1)%str = "TX"
        self%comp_names(2)%str = "TY"
        !
    end function ReceiverFullVerticalMagnetic_ctor
    !
    subroutine ReceiverFullVerticalMagnetic_dtor( self )
        implicit none
        !
        type( ReceiverFullVerticalMagnetic_t ), intent( inout ) :: self
        !
        !write(*,*) "Destructor ReceiverFullVerticalMagnetic_t"
        !
        call self%dealloc()
        !
    end subroutine ReceiverFullVerticalMagnetic_dtor
    !
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
    subroutine writeReceiverFullVerticalMagnetic( self )
        implicit none
        !
        class( ReceiverFullVerticalMagnetic_t ), intent( in ) :: self
        !
        write( *, * ) "Write ReceiverFullVerticalMagnetic_t: ", self%id
        !
    end subroutine writeReceiverFullVerticalMagnetic
    !
    subroutine predictedDataFullVerticalMagnetic( self, transmitter )
        implicit none
        !
        class( ReceiverFullVerticalMagnetic_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in )                     :: transmitter
        !
        complex( kind=prec ) :: comega, det
        !
        complex( kind=prec ), allocatable :: BB(:,:), I_BB(:,:)
        !
        comega = cmplx( 0.0, 1./ ( 2.0 * PI / transmitter%period ), kind=prec )
        !
        !
        allocate( BB( 3, 2 ) )
        select type( tx_e_1 => transmitter%e_all( 1 ) )
            class is( cVector3D_SG_t )
                !
                select type( tx_e_2 => transmitter%e_all( 2 ) )
                    class is( cVector3D_SG_t )
                        !
                        BB(1,1) = dotProdSparse( self%Lbx, tx_e_1 )
                        BB(2,1) = dotProdSparse( self%Lby, tx_e_1 )
                        BB(1,2) = dotProdSparse( self%Lbx, tx_e_2 )
                        BB(2,2) = dotProdSparse( self%Lby, tx_e_2 )
                        BB(3,1) = dotProdSparse( self%Lbz, tx_e_1 )
                        BB(3,2) = dotProdSparse( self%Lbz, tx_e_2 )
                        BB = isign * BB * comega
                        !
                        !invert horizontal B matrix using Kramer"s rule.
                        det = BB( 1, 1 ) * BB( 2, 2 ) - BB( 1, 2 ) * BB( 2, 1 )
                        !
                        !write(*,*) "det:", det
                        !
                        allocate( I_BB( 2, 2 ) )
                        !
                        if( det /= 0 ) then
                            I_BB( 1, 1 ) =  BB( 2, 2 ) / det
                            I_BB( 2, 2 ) =  BB( 1, 1 ) / det
                            I_BB( 1, 2 ) = -BB( 1, 2 ) / det
                            I_BB( 2, 1 ) = -BB( 2, 1 ) / det
                        else
                            stop "ReceiverFullImpedance.f90: Determinant is Zero!"
                        endif
                        !
                        allocate( self%response( 2 ) )
                        !
                        self%response(1) = BB(3,1) * I_BB(1,1) + BB(3,2) * I_BB(2,1)
                        self%response(2) = BB(3,1) * I_BB(1,2) + BB(3,2) * I_BB(2,2)
                        !
                        deallocate( BB )
                        deallocate( I_BB )
                        !
                        ! WRITE ON PredictedFile.dat
                        call self%savePredictedData( transmitter )
                        !
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
    end subroutine predictedDataFullVerticalMagnetic
    !
    subroutine savePredictedDataFullVerticalMagnetic( self, tx )
        implicit none
        !
        class( ReceiverFullVerticalMagnetic_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in ) :: tx
        !
        character(:), allocatable :: code, component
        real( kind=prec )         :: period, real_part, imaginary, rx_location(3)
        integer                   :: i, rx_type
        !
        !#Period(s) Code GG_Lat GG_Lon X(m) Y(m) self%response(m) Component Real Imag Error
        !
        if( associated( self%predicted_data ) ) call deallocateDataHandleArray( self%predicted_data )
        !
        do i = 1, self%n_comp
            !
            rx_type = int( self%rx_type )
            period = real( tx%period, kind=prec )
            code = trim( self%code )
            rx_location = (/real( self%location( 1 ), kind=prec ), real( self%location( 2 ), kind=prec ), real( self%location( 3 ), kind=prec )/)
            component = trim( self%comp_names( i )%str )
            real_part = real( self%response( i ), kind=prec )
            imaginary = real( imag( self%response( i ) ), kind=prec )
            !
            call updateDataHandleArray( self%predicted_data, DataHandleMT_t( rx_type, code, component, period, rx_location, real_part, imaginary ) )
            !
        enddo
        !
    end subroutine savePredictedDataFullVerticalMagnetic
    !
end module ReceiverFullVerticalMagnetic
