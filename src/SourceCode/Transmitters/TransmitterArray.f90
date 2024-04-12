!
!> Class to provide a dynamic and polymorphic transmitters of Transmitter_t objects
!
module TransmitterArray
    !
    use TransmitterMT
    !
    !> Global generic array of Transmitters
    type( TransmitterMT_t ), allocatable, dimension(:) :: transmitters
    !
    public :: printTransmitterArray, countTransmitterTypes
    public :: updateTransmitterArray!, deallocateTransmitterArray
    !
contains
    !
    !> Add a new Transmitter_t and initialize it if necessary
    function updateTransmitterArray( new_tx ) result( i_tx )
        implicit none
        !
        type( TransmitterMT_t ), intent( in ) :: new_tx
        integer :: i_tx
        !
        integer :: iTx, n_tx
        type( TransmitterMT_t ), allocatable, dimension(:) :: temp_array
        !
        if( .NOT. allocated( transmitters ) ) then
            !
            allocate( transmitters(1) )
            !
            transmitters(1) = new_tx
            !
            i_tx = 1
            !
            transmitters(1)%i_tx = 1
            !
        else
            !
            n_tx = size( transmitters )
            !
            do iTx = 1, n_tx
                if( new_tx%isEqual( transmitters( iTx ) ) ) then
                    i_tx = iTx
                    return
                endif
            enddo
            !
            allocate( temp_array( n_tx + 1 ) )
            !
            temp_array( 1 : n_tx ) = transmitters
            !
            temp_array( n_tx + 1 ) = new_tx
            !
            temp_array( n_tx + 1 )%i_tx = n_tx + 1
            !
            i_tx = n_tx + 1
            !
            !call deallocateTransmitterArray()
            !
            deallocate( transmitters )
            !
            allocate( transmitters, source = temp_array )
            !
            deallocate( temp_array )
            !
        endif
        !
    end function updateTransmitterArray
    !
    ! !> No subroutine briefing
    ! subroutine deallocateTransmitterArray()
        ! implicit none
        ! !
        ! integer :: ntx, itx
        ! !
        ! if( allocated( transmitters ) ) then
            ! !
            ! ntx = size( transmitters )
            ! !
            ! if( ntx == 1 ) then
                ! if( allocated( transmitters(1) ) ) deallocate( transmitters(1) )
            ! else
                ! do itx = ntx, 1, -(1)
                    ! !
                    ! if( allocated( transmitters( itx ) ) ) deallocate( transmitters( itx ) )
                    ! !
                ! enddo
            ! endif
            ! !
            ! deallocate( transmitters )
            ! !
        ! endif
        ! !
    ! end subroutine deallocateTransmitterArray
    ! !
    !> Prints the content of the transmitters on screen
    subroutine printTransmitterArray()
        implicit none
        !
        integer :: itx
        !
        write( *, * ) size( transmitters ), " TransmitterArray_t:"
        !
        do itx = 1, size( transmitters )
            call transmitters( itx )%print
        enddo
        !
    end subroutine printTransmitterArray
    !
end module TransmitterArray
!