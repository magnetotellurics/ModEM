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
    public :: updateReceiverArray, deallocateReceiverArray
    !
contains
    !
    ! Add a new Receiver_t and initialize it if necessary
    function updateReceiverArray( new_rx ) result( id )
        implicit none
        !
        class( Receiver_t ), intent( in ) :: new_rx
        integer                           :: id
        !
        integer                             :: iRx, nRx
        type( Rx_t ), pointer, dimension(:) :: temp_array
        type( Rx_t ) :: temp_rx
        !
        id = 0
        !
        if( .NOT. associated( receivers ) ) then
            allocate( receivers( 1 ) )
            !
            temp_rx%Rx = new_rx
            !
            id = 1
            !
            temp_rx%Rx%id = 1
            receivers( 1 ) = temp_rx
            !
        else
            ! 
            nRx = size( receivers )
            !
            do iRx = 1, nRx
                if( new_rx%isEqualRx( receivers( iRx )%Rx ) ) then
                    id = iRx
                    return
                end if
            end do
            !
            allocate( temp_array( nRx + 1 ) )
			!
            temp_array( 1 : nRx ) => receivers(:)
            !
            temp_rx%Rx = new_rx
            temp_rx%Rx%id = nRx + 1
            id = nRx + 1
            !
            allocate( temp_array( nRx + 1 ), source = temp_rx )
            !
            receivers => temp_array
            !
            nullify( temp_array )
            !
        endif
        !
    end function updateReceiverArray
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
    !
    subroutine deallocateReceiverArray()
        integer                    :: nrx, irx
        class( Rx_t ), allocatable :: alloc_rx
        !
        !write( *, * ) "deallocateReceiverArray:", size( receivers )
        !
        nrx = size( receivers )
        do irx = 1, nrx
            alloc_rx = receivers( irx )
            deallocate( alloc_rx )
        end do
        !
        nullify( receivers )
        !
    end subroutine deallocateReceiverArray
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
