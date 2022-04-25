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
    use Transmitter 
    !
    type, extends( Transmitter_t ), public :: TransmitterCSEM_t
        !
        real( kind=prec ) :: location(3), azimuth
        !
        contains
            !
            final    :: TransmitterCSEM_dtor
            !
            procedure, public    :: solveFWD => solveFWDTransmitterCSEM
            procedure, public    :: getSource => getSourceTransmitterCSEM
            !
            procedure, public    :: write => writeTransmitterCSEM
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
    function TransmitterCSEM_ctor( period, location ) result ( self )
        !
        type( TransmitterCSEM_t ) :: self
        !
        real( kind=prec ), intent( in )                         :: period
        real( kind=prec ), intent( in )                         :: location(3)
        !
        ! write(*,*) "Constructor TransmitterCSEM_t"
        !
        call self%init()
        !
        self%n_pol = 1
        self%period = period
        self%location = location
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
        !
        class( TransmitterCSEM_t ), intent( inout ) :: self
        !
        write(*,*) "implementing solveFWD TransmitterCSEM_t: ", self%id
        !
    end subroutine solveFWDTransmitterCSEM
    !
    !
    subroutine getSourceTransmitterCSEM( self )
        !
        class( TransmitterCSEM_t ), intent(in)    :: self
        !
        write(*,*) "getSource TransmitterCSEM_t: ", self%location
        !
    end subroutine getSourceTransmitterCSEM
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
    !function sizeOfTransmitterCSEM( self ) result( size )
        !
        !class( TransmitterCSEM_t ), intent( in ) :: self
        !integer                                          :: size
        !
        !size = sizeof( self%id ) + &
                 !sizeof( self%n_pol ) + &
                 !sizeof( self%fwd_key ) + &
                 !sizeof( self%type ) + &
                 !sizeof( self%period ) + &
                 !sizeof( self%forward_solver ) + &
                 !sizeof( self%e_all ) + &
                 !sizeof( self%receiver_indexes ) + &
                 !sizeof( self%DATA_TITLE )
                 !sizeof( self%location ) + &
                 !sizeof( self%azimuth )
                 !
        !write( *, * ) "Size: ", size
        !
    !end function sizeOfTransmitterCSEM
    !
end module TransmitterCSEM
