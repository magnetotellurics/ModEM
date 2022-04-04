!*************
!
! Class to provide a dynamic and polymorphic receivers of Receiver_t objects
!
! Last modified at 04/2022 by Paulo Werdt
!
!*************
!
module ReceiverFArray
    !
    use Constants
    !
    use Receiver
    !
    implicit none
    !
    ! Allocatable Receiver element of the array, for Old Fortran polymorphism !!!
    type, public :: Rx_t
        !
        class( Receiver_t ), allocatable :: Rx
        !
    end type Rx_t
    !
    ! Global Array of Receivers
    type( Rx_t ), pointer, dimension(:), save, public :: receivers => null()
    !
    public :: getReceiver, printReceiverArray
    public :: updateReceiverArray
    !
contains
    !
    ! Add a new Receiver_t and initialize it if necessary
    subroutine updateReceiverArray( new_rx )
        implicit none
        !
        class( Receiver_t ), intent( in )    :: new_rx
        !
        integer                                 :: iRx, istat
        type( Rx_t ), allocatable, dimension(:)    :: temp_array
        type( Rx_t ), allocatable                  :: temp_rx
        !
        if( .NOT. associated( receivers ) ) then
            allocate( receivers( 1 ) )
            allocate( Rx_t :: temp_rx )
            temp_rx%Rx = new_rx
            temp_rx%Rx%id = 1
            receivers( 1 ) = temp_rx
        else
            !
            do iRx = 1, size( receivers )
                if( new_rx%isEqual( receivers( iRx )%Rx ) ) then
                    return
                end if
            end do
            !
            allocate( temp_array( size( receivers ) + 1 ), STAT=istat )
            temp_array( 1 : size( receivers ) ) = receivers
            allocate( Rx_t :: temp_rx )
            temp_rx%Rx = new_rx
            temp_rx%Rx%id = size( receivers ) + 1
            !
            temp_array( size( receivers ) + 1 ) = temp_rx
            !
            allocate( receivers, source = temp_array, STAT=istat )
            !
        endif
        !
    end subroutine updateReceiverArray
    !
    function getReceiver( iRx ) result( rx )
        !
        integer                                        :: iRx
        !
        class( Receiver_t ), pointer            :: rx
        !
        rx => receivers( iRx )%Rx
        !
    end function getReceiver
    !
    ! Prints the content of the receivers on screen
    subroutine printReceiverArray()
        integer                    :: irx
        class( Rx_t ), allocatable :: alloc_rx
        !
        print *, size( receivers ), " ReceiverFArray_t:"
        !
        do irx = 1, size( receivers )
            alloc_rx = receivers( irx )
            call alloc_rx%Rx%write()
        end do
        !
    end subroutine printReceiverArray
    !
end module ReceiverFArray
