!*************
!
! Class to provide a dynamic and polymorphic receivers of Receiver_t objects
!
!*************
!
module ReceiverFArray
    !
    use Constants
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
    type( Rx_t ), allocatable, target, dimension(:), save, public :: receivers
    !
    public :: getReceiver, printReceiverArray
    public :: updateReceiverArray, deallocateReceiverArray
    !
contains
    !
    ! Add a new Receiver_t and initialize it if necessary
    function updateReceiverArray( new_rx ) result( id_rx )
        implicit none
        !
        class( Receiver_t ), intent( in ) :: new_rx
        integer                           :: id_rx
        !
        integer                                 :: iRx, nRx
        type( Rx_t ), allocatable, dimension(:) :: temp_array
        type( Rx_t ), allocatable               :: temp_rx
        !
        if( .NOT. allocated( receivers ) ) then
            allocate( receivers( 1 ) )
            allocate( Rx_t :: temp_rx )
            temp_rx%Rx = new_rx
            id_rx = 1
            temp_rx%Rx%id = 1
            receivers( 1 ) = temp_rx
            deallocate( temp_rx )
        else
            ! 
            nRx = size( receivers )
            !
            do iRx = 1, nRx
                if( new_rx%isEqualRx( receivers( iRx )%Rx ) ) then
                    id_rx = iRx
                    return
                end if
            end do
            !
            allocate( temp_array( nRx + 1 ) )
            temp_array( 1 : nRx ) = receivers
            allocate( Rx_t :: temp_rx )
            temp_rx%Rx = new_rx
            temp_rx%Rx%id = nRx + 1
            id_rx = nRx + 1
            !
            temp_array( nRx + 1 ) = temp_rx
            !
            if( allocated( receivers ) ) call deallocateReceiverArray()
            allocate( receivers, source = temp_array )
            !
            deallocate( temp_rx )
            deallocate( temp_array )
            !
        endif
        !
    end function updateReceiverArray
    !
    function getReceiver( iRx ) result( rx )
        implicit none
        !
        integer                      :: iRx
        !
        class( Receiver_t ), pointer :: rx
        !
        rx => receivers( iRx )%Rx
        !
    end function getReceiver
    !
    !
    subroutine deallocateReceiverArray()
        implicit none
        !
        integer :: nrx, irx
        !
        !write( *, * ) "deallocateReceiverArray:", size( receivers )
        !
        nrx = size( receivers )
        !
        if( nrx == 1 ) then
            deallocate( receivers(1)%Rx )
        else
            do irx = nrx, 1, -(1)
                deallocate( receivers(irx)%Rx )
            end do
        endif
        !
        deallocate( receivers )
        !
    end subroutine deallocateReceiverArray
    !
    ! Prints the content of the receivers on screen
    subroutine printReceiverArray()
        implicit none
        !
        integer :: irx
        !
        write( *, * ) size( receivers ), " ReceiverFArray_t:"
        !
        do irx = 1, size( receivers )
            call receivers(irx)%Rx%print()
        end do
        !
    end subroutine printReceiverArray
    !
end module ReceiverFArray
