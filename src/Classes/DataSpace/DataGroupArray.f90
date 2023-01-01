!
!> Module with routines to handle DataGroup arrays
!
module DataGroupArray
    !
    use DataGroup
    !
    !> Global DataGroup array for all data measured from a file
    type( DataGroup_t ), allocatable, dimension(:) :: measured_data
    !
    public :: updateDataGroupArray
    public :: setDataGroup
    public :: getDataGroupByIndex
    public :: getDataGroupByRxTx
    public :: printDataGroupArray
    !
contains
    !
    !> Dynamically add a new DataGroup to the array, always via reallocation.
    subroutine updateDataGroupArray( data_array, data_group )
        implicit none
        !
        type( DataGroup_t ), allocatable, dimension(:), intent( inout ) :: data_array
        type( DataGroup_t ), intent( in ) :: data_group
        !
        integer :: n_data
        !
        type( DataGroup_t ), allocatable, dimension(:) :: temp_array
        !
        if( .NOT. allocated( data_array ) ) then
            !
            allocate( data_array(1) )
            !
            data_array(1) = data_group
            !
        else
            !
            n_data = size( data_array )
            !
            allocate( temp_array( n_data + 1 ) )
            !
            temp_array( 1 : n_data ) = data_array(:)
            !
            temp_array( n_data + 1 ) = data_group
            !
            deallocate( data_array )
            !
            allocate( data_array, source = temp_array )
            !
            deallocate( temp_array )
            !
        endif
        !
    end subroutine updateDataGroupArray
    !
    !> Replaces a specific DataGroup in the array with a new one,
    !> if it has the same transmitter-receiver indexes.
    subroutine setDataGroup( data_array, data_group )
        implicit none
        !
        type( DataGroup_t ), allocatable, dimension(:), intent( inout ) :: data_array
        type( DataGroup_t ), intent( in ) :: data_group
        !
        integer :: i_data, n_data
        !
        n_data = size( data_array )
        !
        do i_data = 1, n_data
            !
            if( data_array( i_data )%i_rx == data_group%i_rx .AND. &
                data_array( i_data )%i_tx == data_group%i_tx ) then
                !
                data_array( i_data ) = data_group
                !
                return
                !
            endif
            !
        enddo
        !
    end subroutine setDataGroup
    !
    !> Returns a pointer, allowing modifications directly to a DataGroup at a given index
    function getDataGroupByIndex( data_array, i_data ) result( data_group )
        implicit none
        !
        type( DataGroup_t ), dimension(:), target, intent( in ) :: data_array
        integer, intent( in ) :: i_data
        !
        type( DataGroup_t ), pointer :: data_group
        !
        data_group => data_array( i_data )
        !
    end function getDataGroupByIndex
    !
    !> Returns a pointer, allowing modifications directly to a DataGroup at a given transmitter-receiver pair indexes
    function getDataGroupByRxTx( data_array, i_rx, i_tx ) result( data_group )
        implicit none
        !
        type( DataGroup_t ), target, dimension(:), intent( in ) :: data_array
        integer, intent( in ) :: i_rx, i_tx
        !
        type( DataGroup_t ), pointer :: data_group
        !
        integer :: i_data, n_data
        !
        n_data = size( data_array )
        !
        do i_data = 1, n_data
            if( data_array( i_data )%i_rx == i_rx .AND. &
                data_array( i_data )%i_tx == i_tx ) then
                !
                data_group => data_array( i_data )
                !
                return
                !
            endif
            !
        enddo
        !
    end function getDataGroupByRxTx
    !
    !> Call the print routine of each DataGroup in the array
    subroutine printDataGroupArray( data_array )
        implicit none
        !
        type( DataGroup_t ), dimension(:), intent( in ) :: data_array
        !
        integer :: i_data
        !
        do i_data = 1, size( data_array )
            !
            call data_array( i_data )%print()
            !
        enddo
        !
    end subroutine printDataGroupArray
    !
end module DataGroupArray
!