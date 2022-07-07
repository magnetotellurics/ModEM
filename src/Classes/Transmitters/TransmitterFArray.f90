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
    type( Tx_t ), allocatable, target, dimension(:), save, public :: transmitters
    !
    public :: getTransmitter, printTransmitterArray
    public :: updateTransmitterArray, deallocateTransmitterArray
    !
contains
    !
    ! Add a new Transmitter_t and initialize it if necessary
    subroutine updateTransmitterArray( new_tx )
        implicit none
        !
        class( Transmitter_t ), intent( in )    :: new_tx
        !
        integer                                 :: iTx, nTx
        type( Tx_t ), allocatable, dimension(:) :: temp_array
        type( Tx_t ), allocatable               :: temp_tx
        !
        if( .NOT. allocated( transmitters ) ) then
            allocate( transmitters(1) )
            allocate( Tx_t :: temp_tx )
            temp_tx%Tx = new_tx
            temp_tx%Tx%id = 1
            transmitters( 1 ) = temp_tx
            deallocate( temp_tx )
        else
            !
            nTx = size( transmitters )
            !
            do iTx = 1, size( transmitters )
                if( new_tx%isEqualTx( transmitters( iTx )%Tx ) ) then
                    return
                end if
            end do
            !
            allocate( temp_array( nTx + 1 ) )
            temp_array( 1 : nTx ) = transmitters
            allocate( Tx_t :: temp_tx )
            temp_tx%Tx = new_tx
            temp_tx%Tx%id = nTx + 1
            !
            temp_array( nTx + 1 ) = temp_tx
            !
            if( allocated( transmitters ) ) deallocate( transmitters )
            allocate( transmitters, source = temp_array )
            !
            deallocate( temp_tx )
            deallocate( temp_array )
            !
        endif
        !
    end subroutine updateTransmitterArray
    !
    function getTransmitter( iTx ) result( tx )
        implicit none
        !
        integer                         :: iTx
        !
        class( Transmitter_t ), pointer :: tx
        !
        tx => transmitters( iTx )%Tx
        !
    end function getTransmitter
    !
    !
    subroutine deallocateTransmitterArray()
      implicit none
        !
        integer                :: ntx, itx
        class( Tx_t ), pointer :: alloc_tx
        !
        !write( *, * ) "deallocateTransmitterArray:", size( transmitters )
        !
        ntx = size( transmitters )
        do itx = ntx, 1, -(1)
            alloc_tx => transmitters( itx )
            if( associated( alloc_tx ) ) nullify( alloc_tx )
        end do
        !
        !if( allocated( transmitters ) ) deallocate( transmitters )
        !
    end subroutine deallocateTransmitterArray
    !
    ! Prints the content of the transmitters on screen
    subroutine printTransmitterArray()
      implicit none
        !
        integer                :: itx
        class( Tx_t ), pointer :: alloc_tx
        !
        print *, size( transmitters ), " TransmitterFArray_t:"
        !
        do itx = 1, size( transmitters )
            alloc_tx => transmitters( itx )
            call alloc_tx%Tx%write()
        end do
        !
    end subroutine printTransmitterArray
    !
end module TransmitterFArray
