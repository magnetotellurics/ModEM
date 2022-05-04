! *************
! 
! Derived class to define a CSEM Transmitter
! 
! Last modified at 08/06/2021 by Paulo Werdt
! 
! *************
! 
module TransmitterCSEM
    ! 
    use FileUnits
    use Transmitter 
    use cVector3D_SG
    !
    type, extends( Transmitter_t ), public :: TransmitterCSEM_t
        !
        real( kind=prec )  :: location(3), azimuth, dip, moment
        character(:), allocatable :: dipole
        !
        contains
            !
            final    :: TransmitterCSEM_dtor
            !
            procedure, public :: solveFWD => solveFWDTransmitterCSEM
            !
            procedure, public :: isEqual => isEqualTransmitterCSEM
            !
            procedure, public :: write => writeTransmitterCSEM
            !
    end type TransmitterCSEM_t
    !
    interface TransmitterCSEM_t
        module procedure TransmitterCSEM_ctor
    end interface TransmitterCSEM_t
    !
contains
    !
    ! Parametrized constructor
    function TransmitterCSEM_ctor( period, location, azimuth, dip, moment, dipole ) result ( self )
        !
        type( TransmitterCSEM_t ) :: self
        !
        real( kind=prec ), intent( in )  :: period, azimuth, dip, moment
        real( kind=prec ), intent( in )  :: location(3)
        character(:), allocatable, intent( in ) :: dipole
        !
        ! write(*,*) "Constructor TransmitterCSEM_t"
        !
        call self%init()
        !
        self%n_pol = 1
        self%period = period
        self%location = location
        self%azimuth = azimuth
        self%dip = dip
        self%moment = moment
        self%dipole = dipole
        !
    end function TransmitterCSEM_ctor
    !
    ! Destructor
    subroutine TransmitterCSEM_dtor( self )
        implicit none
        !
        type( TransmitterCSEM_t )    :: self
        !
        ! write(*,*) "Destructor TransmitterCSEM_t"
        !
        call self%dealloc()
        !
    end subroutine TransmitterCSEM_dtor
    !
    subroutine solveFWDTransmitterCSEM( self )
        implicit none
        !
        class( TransmitterCSEM_t ), intent( inout ) :: self
        !
        integer           :: ios
        real( kind=prec ) :: omega
        !
        character( len=20 ) :: ModeName
        !
        !
        omega = 2.0 * PI / self%period
        !
        allocate( cVector3D_SG_t :: self%e_all( 1 ) )
        !
        ! Verbosis...
        write( *, "(A20, I8, A20, es20.6, A20, I8)" ) "SolveFWD for Tx:", self%id, " -> Period:", self%period, " - Polarization:", 1
        !
        call self%source%setE( 1 )
        !
        select type( mgrid => self%source%model_operator%metric%grid )
            class is( Grid3D_SG_t )
                !
                self%e_all( 1 ) = cVector3D_SG_t( mgrid, EDGE )
                !
        end select
        !
        call self%forward_solver%getESolution( self%source, self%e_all( 1 ) )
        !
        self%e_all( 1 ) = self%e_all( 1 ) + self%source%E
        !
        ModeName = "Ey"
        !
        open( ioESolution, file = e_solution_file_name, action = "write", position = "append", form = "unformatted", iostat = ios )
        !
        if( ios /= 0 ) then
            stop "Error opening file in solveFWDTransmitterCSEM: e_solution"
        else
            !
            ! write the frequency header - 1 record
            write( ioESolution ) omega, self%id, 1, ModeName
            !
            call self%e_all( 1 )%write( ioESolution )
            !
            close( ioESolution )
            !
        endif
        !
    end subroutine solveFWDTransmitterCSEM
    !
    ! Compare two transmitters
    function isEqualTransmitterCSEM( self, other ) result( equal )
        implicit none
        !
        class( TransmitterCSEM_t ), intent( in ) :: self
        class( Transmitter_t ), intent( in )     :: other
        logical                                  :: equal
        !
        equal = .FALSE.
        !
        select type( other )
            !
            class is( TransmitterCSEM_t )
                !
                if( self%period == other%period .AND.   &
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
    end function isEqualTransmitterCSEM
    !
    subroutine writeTransmitterCSEM( self )
        !
        class( TransmitterCSEM_t ), intent(in)    :: self
        integer                                    :: iRx
        !
        write( *, "(A20, I8, A10, es12.6, A20, I8)") "TransmitterCSEM: ", self%id,    &
        " Period: ",    self%period,    &
        " N Receivers: ", size( self%receiver_indexes )
        !
    end subroutine writeTransmitterCSEM
    !
end module TransmitterCSEM
