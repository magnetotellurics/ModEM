!
!> Class to provide a dynamic and polymorphic receivers of Receiver_t objects
!
module ReceiverArray
    !
    use Constants
    use Receiver
    !
    !> Allocatable Receiver element of the array, for Old Fortran polymorphism !!!
    type, public :: Rx_t
        !
        class( Receiver_t ), allocatable :: Rx
        !
    end type Rx_t
    !
    !> Global generic array of Receivers
    type( Rx_t ), allocatable, target, dimension(:) :: receivers
    !
    public :: getReceiver, printReceiverArray
    public :: updateReceiverArray, deallocateReceiverArray
    !
contains
    !
    !> Add a new Receiver_t and initialize it if necessary
    function updateReceiverArray( new_rx ) result( i_rx )
        implicit none
        !
        class( Receiver_t ), intent( in ) :: new_rx
        integer :: i_rx
        !
        integer :: iRx, n_rx
        type( Rx_t ), allocatable, dimension(:) :: temp_array
        type( Rx_t ), allocatable :: temp_rx
        !
        if( .NOT. allocated( receivers ) ) then
            !
            allocate( receivers(1) )
            allocate( Rx_t :: temp_rx )
            temp_rx%Rx = new_rx
            i_rx = 1
            temp_rx%Rx%i_rx = 1
            receivers(1) = temp_rx
            !
            deallocate( temp_rx )
            !
        else
            !
            n_rx = size( receivers )
            !
            do iRx = 1, n_rx
                if( new_rx%isEqualRx( receivers( iRx )%Rx ) ) then
                    i_rx = iRx
                    return
                endif
            enddo
            !
            allocate( temp_array( n_rx + 1 ) )
            temp_array( 1 : n_rx ) = receivers
            allocate( Rx_t :: temp_rx )
            temp_rx%Rx = new_rx
            temp_rx%Rx%i_rx = n_rx + 1
            i_rx = n_rx + 1
            !
            temp_array( n_rx + 1 ) = temp_rx
            !
            call deallocateReceiverArray()
            allocate( receivers, source = temp_array )
            !
            deallocate( temp_rx, temp_array )
            !
        endif
        !
    end function updateReceiverArray
    !
    !> No function briefing
    function getReceiver( iRx ) result( rx )
        implicit none
        !
        integer :: iRx
        !
        class( Receiver_t ), pointer :: rx
        !
        rx => receivers( iRx )%Rx
        !
    end function getReceiver
    !
    !> No subroutine briefing
    subroutine deallocateReceiverArray()
        implicit none
        !
        integer :: nrx, irx
        !
        !write( *, * ) "deallocateReceiverArray:", size( receivers )
        !
        if( allocated( receivers ) ) then
            !
            nrx = size( receivers )
            !
            if( nrx == 1 ) then
                if( allocated( receivers(1)%Rx ) ) deallocate( receivers(1)%Rx )
            else
                do irx = nrx, 1, -(1)
                    if( allocated( receivers( irx )%Rx ) ) deallocate( receivers( irx )%Rx )
                enddo
            endif
            !
            deallocate( receivers )
            !
        endif
        !
    end subroutine deallocateReceiverArray
    !
    !> Prints the content of the receivers on screen
    subroutine printReceiverArray()
        implicit none
        !
        integer :: irx
        !
        write( *, * ) size( receivers ), " ReceiverArray_t:"
        !
        do irx = 1, size( receivers )
            call receivers( irx )%Rx%print()
        enddo
        !
    end subroutine printReceiverArray
    !
end module ReceiverArray
