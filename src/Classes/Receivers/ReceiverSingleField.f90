! *************
! 
! Derived class to define a Single E or B Field Receiver
! 
! Last modified at 30/06/2021 by Paulo Werdt
! 
! *************
! 
module ReceiverSingleField
    !
    use Receiver
    use DataHandleCSEM
    use TransmitterCSEM
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
            procedure, public :: savePredictedData => savePredictedDataSingleField
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
    function ReceiverSingleField_ctor( location, azimuth, rx_type ) result( self )
        implicit none
        !
        real( kind=prec ), intent( in ) :: location(3)
        real( kind=prec ), intent( in ) :: azimuth
        integer, intent( in )           :: rx_type
        !
        type( ReceiverSingleField_t )    :: self
        !
        integer :: i, asize
        !
        ! write( *, * ) "Constructor ReceiverSingleField_t"
        !
        call self%init()
        !
        self%location = location
        self%azimuth = azimuth
        !
        self%rx_type = rx_type
        !
        self%n_comp = 1
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
        allocate( self%EHxy(1) )
        !
        if( azimuth == 1.0 ) then
            self%EHxy(1)%str = "Ex"
            allocate( self%Lex, source = cSparsevector3D_SG_t() )
        endif
        !
        if( azimuth == 2.0 ) then
            self%EHxy(1)%str = "Ey"
            allocate( self%Ley, source = cSparsevector3D_SG_t() )
        endif
        !
        if( azimuth == 3.0 ) then
            self%EHxy(1)%str = "Bx"
            allocate( self%Lbx, source = cSparsevector3D_SG_t() )
        endif
        !
        if( azimuth == 4.0 ) then
            self%EHxy(1)%str = "By"
            allocate( self%Lby, source = cSparsevector3D_SG_t() )
        endif
        !
        if( azimuth == 5.0 ) then
            self%EHxy(1)%str = "Bz"
            allocate( self%Lbz, source = cSparsevector3D_SG_t() )
        endif
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
    subroutine ReceiverSingleField_dtor( self )
        implicit none
        !
        type( ReceiverSingleField_t ), intent( inout ) :: self
        !
        ! write( *, * ) "Destructor ReceiverSingleField_t"
        !
        call self%dealloc()
        !
    end subroutine ReceiverSingleField_dtor
    !
    subroutine setLRowsSingleField( self, transmitter )
        implicit none
        !
        class( ReceiverSingleField_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in )            :: transmitter
        !
        write( *, * ) "setLRowsSingleField to be implemented"
        !
    end subroutine setLRowsSingleField
    !
    subroutine predictedDataSingleField( self, transmitter )
        implicit none
        !
        class( ReceiverSingleField_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in )            :: transmitter
        !
        integer :: i, j, ij
        complex( kind=prec ) :: comega, det
        !
        !
        comega = cmplx( 0.0, 1./ ( 2.0 * PI / transmitter%period ), kind=prec )
        !
        select type( tx_e_1 => transmitter%e_all(1) )
            class is( cVector3D_SG_t )
                !
                allocate( self%response(1) )
                !
                select case ( self%EHxy(1)%str )
                    case( "Ex" )
                        self%response(1) = self%Lex%dotProd( tx_e_1 )
                    case( "Ey" )
                        self%response(1) = self%Ley%dotProd( tx_e_1 )
                    case( "Bx" )
                        self%response(1) = self%Lbx%dotProd( tx_e_1 )
                        self%response(1) = isign * self%response(1) * comega
                    case( "By" )
                        self%response(1) = self%Lby%dotProd( tx_e_1 )
                        self%response(1) = isign * self%response(1) * comega
                    case( "Bz" )
                        self%response(1) = self%Lbz%dotProd( tx_e_1 )
                        self%response(1) = isign * self%response(1) * comega
                end select
                !
                ! WRITE ON PredictedFile.dat
                call self%savePredictedData( transmitter )
                !
                deallocate( self%response )
                !
            class default
                stop "evaluationFunctionRx: Unclassified temp_full_vec_ey"
        end select
        !
    end subroutine predictedDataSingleField
    !
    subroutine savePredictedDataSingleField( self, tx )
        implicit none
        !
        class( ReceiverSingleField_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in ) :: tx
        !
        integer :: rx_type
        character(:), allocatable :: code, component, dipole
        real( kind=prec )         :: period, tx_location(3), azimuth, dip, moment, real_part, imaginary, rx_location(3)
        !
        select type( tx )
            !
            class is( TransmitterCSEM_t )
                !
                rx_type = int( self%rx_type )
                period = real( tx%period, kind=prec )
                azimuth = real( tx%azimuth, kind=prec )
                dip = real( tx%dip, kind=prec )
                moment = real( tx%moment, kind=prec )
                code = trim( self%code )
                tx_location = (/real( tx%location(1), kind=prec ), real( tx%location( 2 ), kind=prec ), real( tx%location( 3 ), kind=prec )/)
                rx_location = (/real( self%location(1), kind=prec ), real( self%location( 2 ), kind=prec ), real( self%location( 3 ), kind=prec )/)
                dipole = trim( tx%dipole )
                component = trim( self%comp_names(1)%str )
                real_part = real( self%response(1), kind=prec )
                imaginary = real( imag( self%response(1) ), kind=prec )
                !
                if( allocated( self%predicted_data ) ) call deallocateDataHandleArray( self%predicted_data )
                !
                call updateDataHandleArray( self%predicted_data, DataHandleCSEM_t( rx_type, code, component, period, tx_location, azimuth, dip, moment, dipole, rx_location, real_part, imaginary ) )
                !
            class default
                stop "Wrong transmitter for ReceiverSingleField"
            !
        end select
        !
    end subroutine savePredictedDataSingleField
    !
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
    subroutine printReceiverSingleField( self )
        implicit none
        !
        class( ReceiverSingleField_t ), intent( in ) :: self
        !
        write( *, * ) "Print ReceiverSingleField_t: ", self%id
        !
    end subroutine printReceiverSingleField
    !
end module ReceiverSingleField
