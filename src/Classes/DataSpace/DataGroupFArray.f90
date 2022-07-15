!*************
!
! Class to provide a dynamic and polymorphic data_group_array of DataGroup_t objects
!
!*************
!
module DataGroupFArray
    !
    use DataGroup
    !
    ! Global Array of DataGroups
    type( DataGroup_t ), allocatable, target, dimension(:), save, public :: original_data
    type( DataGroup_t ), allocatable, target, dimension(:), save, public :: predicted_data
    !
    public :: getDataGroupByIndex, getDataGroupByRxTx, setDataGroup, updateDataGroupArray
    !
contains
    !
    ! Add a new DataGroup_t and initialize it if necessary
    subroutine updateDataGroupArray( data_group_array, new_dg )
        implicit none
        !
        type( DataGroup_t ), allocatable, dimension(:), intent( inout ) :: data_group_array
        type( DataGroup_t ), intent( in ) :: new_dg
        !
        integer :: iDg, nDg
        type( DataGroup_t ), allocatable, dimension(:) :: temp_array
        !
        if( .NOT. allocated( data_group_array ) ) then
            !
            allocate( data_group_array( 1 ) )
            !
            data_group_array( 1 ) = new_dg
            !
        else
            !
            nDg = size( data_group_array )
            !
            allocate( temp_array( nDg + 1 ) )
            !
            temp_array( 1 : nDg ) = data_group_array(:)
            !
            temp_array( nDg + 1 ) = new_dg
            !
            if( allocated( data_group_array ) ) deallocate( data_group_array )
            allocate( data_group_array, source = temp_array )
            !
            deallocate( temp_array )
            !
        endif
        !
    end subroutine updateDataGroupArray
!
    ! Add a new DataGroup_t and initialize it if necessary
    subroutine setDataGroup( data_group_array, new_dg )
        implicit none
        !
        type( DataGroup_t ), allocatable, dimension(:), intent( inout ) :: data_group_array
        type( DataGroup_t ), intent( in ) :: new_dg
        !
        integer :: idg, ndg
        !
        ndg = size( data_group_array )
        !
        do idg = 1, ndg
            if( data_group_array( idg )%id_rx == new_dg%id_rx .AND. &
                data_group_array( idg )%id_tx == new_dg%id_tx ) then
                !
                data_group_array( idg )%reals = new_dg%reals
                data_group_array( idg )%imaginaries = new_dg%imaginaries
                data_group_array( idg )%errors = new_dg%errors
                !
                return
                !
            endif
            !
        enddo
        !
    end subroutine setDataGroup
    !
    function getDataGroupByIndex( data_group_array, iDg ) result( data_group )
        implicit none
        !
        type( DataGroup_t ), target, dimension(:), intent( in ) :: data_group_array
        integer                                                 :: iDg
        !
        type( DataGroup_t ), pointer :: data_group
        !
        data_group => data_group_array( iDg )
        !
    end function getDataGroupByIndex
    !
    function getDataGroupByRxTx( data_group_array, id_rx, id_tx ) result( data_group )
        implicit none
        !
        type( DataGroup_t ), target, dimension(:), intent( in ) :: data_group_array
        integer                                                 :: id_rx, id_tx
        !
        type( DataGroup_t ), pointer :: data_group
        !
        integer :: idg, ndg
        !
        ndg = size( data_group_array )
        !
        do idg = 1, ndg
            if( data_group_array( idg )%id_rx == id_rx .AND. &
                data_group_array( idg )%id_tx == id_tx ) then
                !
                data_group => data_group_array( iDg )
                !
                return
                !
            endif
            !
        enddo
        !
    end function getDataGroupByRxTx
    !
end module DataGroupFArray
