!*************
!
! Class to provide a dynamic and polymorphic data_groups of DataGroup_t objects
!
!*************
!
module DataGroupFArray
    !
    use Constants
    !
    use DataGroup
    !
    implicit none
    !
    ! Global Array of DataGroups
    type( DataGroup_t ), allocatable, target, dimension(:), save, public :: data_groups
    !
    public :: getDataGroup, printDataGroupArray
    public :: updateDataGroupArray
    !
contains
    !
    ! Add a new DataGroup_t and initialize it if necessary
    subroutine updateDataGroupArray( new_data_group )
        implicit none
        !
        class( DataGroup_t ), intent( in ) :: new_data_group
        !
        integer :: iDg, nDg
        !
        type( DataGroup_t ), allocatable, dimension(:) :: temp_array
        class( DataGroup_t ), allocatable              :: temp_data_group
        !
        if( .NOT. allocated( data_groups ) ) then
            allocate( data_groups(1) )
            allocate( DataGroup_t :: temp_data_group )
            temp_data_group = new_data_group
            temp_data_group%id = 1
            data_groups(1) = temp_data_group
            deallocate( temp_data_group )
        else
            !
            nDg = size( data_groups )
            !
            do iDg = 1, size( data_groups )
                if( new_data_group%isEqual( data_groups( iDg ) ) ) then
                    return
                end if
            end do
            !
            allocate( temp_array( nDg + 1 ) )
            temp_array( 1 : nDg ) = data_groups
            allocate( DataGroup_t :: temp_data_group )
            temp_data_group = new_data_group
            temp_data_group%id = nDg + 1
            !
            temp_array( nDg + 1 ) = temp_data_group
            !
            if( allocated( data_groups ) ) deallocate( data_groups )
            allocate( data_groups, source = temp_array )
            !
            deallocate( temp_data_group )
            deallocate( temp_array )
            !
        endif
        !
    end subroutine updateDataGroupArray
    !
    function getDataGroup( iDg ) result( data_group )
        implicit none
        !
        integer :: iDg
        !
        class( DataGroup_t ), pointer :: data_group
        !
        data_group => data_groups( iDg )
        !
    end function getDataGroup
    !
    ! Prints the content of the data_groups on screen
    subroutine printDataGroupArray()
        implicit none
        !
        integer :: idata_group
        !
        write( *, * ) "          Checked ", size( data_groups ), " DataGroups:"
        !
        do idata_group = 1, size( data_groups )
            call data_groups( idata_group )%print()
        end do
        !
    end subroutine printDataGroupArray
    !
end module DataGroupFArray
