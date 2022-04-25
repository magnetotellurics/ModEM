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
    function ReceiverFullImpedance_ctor( location, type_name ) result( self )
        implicit none
        !
        real( kind=prec ), intent( in ) :: location(3)
		character(:), allocatable, optional, intent( in ) :: type_name
        !
		type( ReceiverFullImpedance_t ) :: self
		!
        !write(*,*) "Constructor ReceiverFullImpedance_t"
        !
        call self%init()
        !
        self%location = location
        !
        if( present( type_name ) ) then
            self%type_name = type_name
        else
            self%type_name = "ReceiverFullImpedance"
        endif
        !
        self%DATA_TITLE = "Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error"
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
    subroutine predictedDataFullImpedance( self, model_operator, transmitter )
        implicit none
        !
        class( ReceiverFullImpedance_t ), intent( inout ) :: self
        class( ModelOperator_t ), intent( in )            :: model_operator
        class( Transmitter_t ), intent( in )              :: transmitter
        !
        complex( kind=prec ), allocatable :: BB(:,:), det
        real( kind=prec )                 :: omega
        integer                           :: i, j, ij
        !
        omega = ( 2.0 * PI / transmitter%period )
        !
        ! Set Vectors Lex, Ley, Lbx, Lby
        call self%evaluationFunction( model_operator, omega )
        !
        allocate( self%EE( 2, 2 ) )
        !
        self%EE( 1, 1 ) = self%Lex .dot. transmitter%e_all( 1 )
        self%EE( 2, 1 ) = self%Ley .dot. transmitter%e_all( 1 )
        self%EE( 1, 2 ) = self%Lex .dot. transmitter%e_all( 2 )
        self%EE( 2, 2 ) = self%Ley .dot. transmitter%e_all( 2 )
        !
        deallocate( self%Lex )
        deallocate( self%LeY )
        !
        !write(*,*) "EE:"
        !write(*,*) self%EE( 1, 1 ), self%EE( 1, 2 )
        !write(*,*) self%EE( 2, 1 ), self%EE( 2, 2 )
        !
        allocate( BB( 2, 2 ) )
        !
        BB( 1, 1 ) = self%Lbx .dot. transmitter%e_all( 1 )
        BB( 2, 1 ) = self%Lby .dot. transmitter%e_all( 1 )
        BB( 1, 2 ) = self%Lbx .dot. transmitter%e_all( 2 )
        BB( 2, 2 ) = self%Lby .dot. transmitter%e_all( 2 )
        !
        deallocate( self%Lbx )
        deallocate( self%Lby )
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
    end subroutine predictedDataFullImpedance
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
