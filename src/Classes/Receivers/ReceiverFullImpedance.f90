! *************
! 
! Derived class to define a Full_Impedance Receiver
! 
! *************
! 
module ReceiverFullImpedance
    !
    use Receiver
	use DataHandleMT
    !
    type, extends( Receiver_t ), public :: ReceiverFullImpedance_t
        !
        ! SPECIFIC PROPERTIES HERE
        !
        contains
            !
            final :: ReceiverFullImpedance_dtor
            !
            procedure, public :: isEqualRx => isEqualFullImpedance
            !
            procedure, public :: predictedData => predictedDataFullImpedance
			!
			procedure, public :: savePredictedData => savePredictedDataFullImpedance
            !
            procedure, public :: write => writeReceiverFullImpedance
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
        character(:), allocatable :: aux_str
        !
        !write(*,*) "Constructor ReceiverFullImpedance_t"
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
        allocate( self%EHxy( 4 ) )
        !
        self%EHxy(1)%str = "Ex"
        self%EHxy(2)%str = "Ey"
        self%EHxy(3)%str = "Bx"
        self%EHxy(4)%str = "By"
        !
        ! components required to get the full impdedance tensor response [Zxx, Zxy, Zyx, Zyy]
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
        !write(*,*) "Destructor ReceiverFullImpedance_t:", self%id
        !
        call self%dealloc()
        !
    end subroutine ReceiverFullImpedance_dtor
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
    subroutine predictedDataFullImpedance( self, transmitter )
        implicit none
        !
        class( ReceiverFullImpedance_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in )              :: transmitter
        !
        complex( kind=prec ) :: comega
        !
        complex( kind=prec ), allocatable :: BB(:,:), det
        integer                           :: i, j, ij
        !
        comega = cmplx( 0.0, 1./ ( 2.0 * PI / transmitter%period ), kind=prec )
        !
        !
        allocate( self%EE( 2, 2 ) )
        !
        select type( tx_e_1 => transmitter%e_all( 1 ) )
            class is( cVector3D_SG_t )
                !
                select type( tx_e_2 => transmitter%e_all( 2 ) )
                    class is( cVector3D_SG_t )
                        !
                        self%EE( 1, 1 ) = dotProdSparse( self%Lex, tx_e_1 )
                        self%EE( 2, 1 ) = dotProdSparse( self%Ley, tx_e_1 )
                        self%EE( 1, 2 ) = dotProdSparse( self%Lex, tx_e_2 )
                        self%EE( 2, 2 ) = dotProdSparse( self%Ley, tx_e_2 )
                        !
                        !write(*,*) "EE:"
                        !write(*,*) self%EE( 1, 1 ), self%EE( 1, 2 )
                        !write(*,*) self%EE( 2, 1 ), self%EE( 2, 2 )
                        !
                        allocate( BB( 2, 2 ) )
                        !
                        BB( 1, 1 ) = dotProdSparse( self%Lbx, tx_e_1 )
                        BB( 2, 1 ) = dotProdSparse( self%Lby, tx_e_1 )
                        BB( 1, 2 ) = dotProdSparse( self%Lbx, tx_e_2 )
                        BB( 2, 2 ) = dotProdSparse( self%Lby, tx_e_2 )
                        BB = isign * BB * comega
                        !
                        !write(*,*) "BB:"
                        !write(*,*) BB( 1, 1 ), BB( 1, 2 )
                        !write(*,*) BB( 2, 1 ), BB( 2, 2 )
                        !
                        !invert horizontal B matrix using Kramer"s rule.
                        det = BB( 1, 1 ) * BB( 2, 2 ) - BB( 1, 2 ) * BB( 2, 1 )
                        !
                        !write(*,*) "det:", det
                        !
                        allocate( self%I_BB( 2, 2 ) )
                        !
                        if( det /= 0 ) then
                            self%I_BB( 1, 1 ) =  BB( 2, 2 ) / det
                            self%I_BB( 2, 2 ) =  BB( 1, 1 ) / det
                            self%I_BB( 1, 2 ) = -BB( 1, 2 ) / det
                            self%I_BB( 2, 1 ) = -BB( 2, 1 ) / det
                        else
                            STOP "ReceiverFullImpedance.f90: Determinant is Zero!"
                        endif
                        !
                        !write(*,*) "Inverse BB:"
                        !write(*,*) self%I_BB( 1, 1 ), self%I_BB( 1, 2 )
                        !write(*,*) self%I_BB( 2, 1 ), self%I_BB( 2, 2 )
                        !
                        deallocate( BB )
                        !
                        allocate( self%response( 4 ) )
                        !
                        do j = 1, 2
                             do i = 1, 2
                                 ij = 2 * ( i-1 ) + j
                                 self%response( ij ) = self%EE( i, 1 ) * self%I_BB( 1, j ) + self%EE( i, 2 ) * self%I_BB( 2, j )
                             enddo
                        enddo
                        !
                        deallocate( self%EE )
                        deallocate( self%I_BB )
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
            class default
                stop "evaluationFunctionRx: Unclassified temp_full_vec_ey"
        end select
        !
    end subroutine predictedDataFullImpedance
    !
    subroutine savePredictedDataFullImpedance( self, tx )
        implicit none
        !
        class( ReceiverFullImpedance_t ), intent( inout ) :: self
        class( Transmitter_t ), intent( in )              :: tx
        !
        character(:), allocatable :: code, component
        real( kind=prec )         :: period, real_part, imaginary, rx_location(3)
        integer                   :: i, rx_type
        !
        !#Period(s) Code GG_Lat GG_Lon X(m) Y(m) response(m) Component Real Imag Error
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
    end subroutine savePredictedDataFullImpedance
    !
    subroutine writeReceiverFullImpedance( self )
        implicit none
        !
        class( ReceiverFullImpedance_t ), intent( in ) :: self
        !
        write(*,*) "Write ReceiverFullImpedance_t: ", self%id
        !
    end subroutine writeReceiverFullImpedance
    !
end module ReceiverFullImpedance
