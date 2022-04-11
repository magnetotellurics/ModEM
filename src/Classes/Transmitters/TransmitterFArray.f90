!*************
!
! Class to provide a dynamic and polymorphic transmitters of Transmitter_t objects
!
! Last modified at 04/2022 by Paulo Werdt
!
!*************
!
module TransmitterFArray
    !
    use Constants
    !
    use TransmitterMT
    use TransmitterCSEM
    !
    implicit none
    !
    ! Allocatable Transmitter element of the array, for Old Fortran polymorphism !!!
    type, public :: Tx_t
        !
        class( Transmitter_t ), allocatable :: Tx
        !
    end type Tx_t
    !
    ! Global Array of Transmitters
    type( Tx_t ), pointer, dimension(:), save, public :: transmitters => null()
    !
    public :: getTransmitter, printTransmitterArray
    public :: updateTransmitterArray
    !
contains
    !
    ! Add a new Transmitter_t and initialize it if necessary
    subroutine updateTransmitterArray( new_tx )
        implicit none
        !
        class( Transmitter_t ), intent( in )    :: new_tx
        !
        integer                                 :: iTx, istat
        type( Tx_t ), allocatable, dimension(:)    :: temp_array
        type( Tx_t ), allocatable                  :: temp_tx
        !
        if( .NOT. associated( transmitters ) ) then
            allocate( transmitters( 1 ) )
            allocate( Tx_t :: temp_tx )
            temp_tx%Tx = new_tx
            temp_tx%Tx%id = 1
            transmitters( 1 ) = temp_tx
        else
            !
            do iTx = 1, size( transmitters )
                if ( new_tx%isEqual( transmitters( iTx )%Tx ) ) then
                    return
                end if
            end do
            !
            allocate( temp_array( size( transmitters ) + 1 ), STAT=istat )
            temp_array( 1 : size( transmitters ) ) = transmitters
            allocate( Tx_t :: temp_tx )
            temp_tx%Tx = new_tx
            temp_tx%Tx%id = size( transmitters ) + 1
            !
            temp_array( size( transmitters ) + 1 ) = temp_tx
            !
            allocate( transmitters, source = temp_array, STAT=istat )
            !
        endif
        !
    end subroutine updateTransmitterArray
    !
    function getTransmitter( iTx ) result( tx )
        !
        integer                                        :: iTx
        !
        class( Transmitter_t ), pointer            :: tx
        !
        tx => transmitters( iTx )%Tx
        !
    end function getTransmitter
    !
    ! Prints the content of the transmitters on screen
    subroutine printTransmitterArray()
        integer                    :: itx
        class( Tx_t ), allocatable :: alloc_tx
        !
        print *, size( transmitters ), " TransmitterFArray_t:"
        !
        do itx = 1, size( transmitters )
            alloc_tx = transmitters( itx )
            call alloc_tx%Tx%write()
        end do
        !
    end subroutine printTransmitterArray
    !
end module TransmitterFArray
