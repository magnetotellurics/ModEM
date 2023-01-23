!
!> Class to provide a dynamic and polymorphic data_entry_array of DataEntry_t objects
!
module DataEntryArray
    !
    use DataEntry
    !
    !> Allocatable DataEntry element of the array, for Old Fortran polymorphism !!!
    type, public :: De_t
        !
        class( DataEntry_t ), allocatable :: De
        !
    end type De_t
    !
    public :: getDataEntry, updateDataEntryArray, hasDataEntry, deallocateDataEntryArray
    !
contains
    !
    !> Add a new DataEntry_t and initialize it if necessary
    subroutine updateDataEntryArray( data_entry_array, new_De )
        implicit none
        !
        type( De_t ), allocatable, dimension(:), intent( inout ) :: data_entry_array
        class( DataEntry_t ), intent( in ) :: new_De
        !
        integer :: iDe, nDe
        type( De_t ), allocatable, dimension(:) :: temp_array
        type( De_t ) :: temp_De
        !
        if( .NOT. allocated( data_entry_array ) ) then
            !
            allocate( data_entry_array(1) )
            !
            temp_De%De = new_De
            !
            data_entry_array(1) = temp_De
            !
        else
            !
            nDe = size( data_entry_array )
            !
            allocate( temp_array( nDe + 1 ) )
            !
            temp_array( 1 : nDe ) = data_entry_array(:)
            !
            temp_De%De = new_De
            !
            temp_array( nDe + 1 ) = temp_De
            !
            if( allocated( data_entry_array ) ) deallocate( data_entry_array )
            allocate( data_entry_array, source = temp_array )
            !
            deallocate( temp_array )
            !
        endif
        !
    end subroutine updateDataEntryArray
    !
    !> No function briefing
    function hasDataEntry( data_entry_array, data_entry ) result( exist )
        implicit none
        !
        type( De_t ), dimension(:), intent( in ) :: data_entry_array
        class( DataEntry_t ), intent( in ) :: data_entry
        !
        logical :: exist
        integer :: iDe, nDe
        !
        exist = .FALSE.
        !
        nDe = size( data_entry_array )
        do iDe = 1, nDe
            if( data_entry%isEqual( data_entry_array( iDe )%De ) ) then
                exist = .TRUE.
                exit
            endif
        enddo
        !
    end function hasDataEntry
    !
    !> No subroutine briefing
    subroutine setDataEntry( data_entry_array, iDe, De )
        implicit none
        !
        type( De_t ), dimension(:), intent( inout ) :: data_entry_array
        integer, intent( in ) :: iDe
        !
        class( DataEntry_t ), allocatable, intent( in ) :: De
        !
        data_entry_array( iDe )%De = De
        !
    end subroutine setDataEntry
    !
    !> No function briefing
    function getDataEntry( data_entry_array, iDe ) result( De )
        implicit none
        !
        type( De_t ), target, dimension(:), intent( in ) :: data_entry_array
        integer :: iDe
        !
        class( DataEntry_t ), pointer :: De
        !
        De => data_entry_array( iDe )%De
        !
    end function getDataEntry
    !
    !> No subroutine briefing
    subroutine deallocateDataEntryArray( data_entry_array )
        implicit none
        !
        type( De_t ), allocatable, dimension(:), intent( inout ) :: data_entry_array
        !
        integer :: nDe, iDe
        !
        !write( *, * ) "deallocateDataEntryArray:", size( data_entry_array )
        !
        nDe = size( data_entry_array )
        !
        if( nDe == 1 ) then
            deallocate( data_entry_array(1)%De )
        else
            do iDe = nDe, 1, -(1)
                deallocate( data_entry_array( iDe )%De )
            enddo
        endif
        !
        deallocate( data_entry_array )
        !
    end subroutine deallocateDataEntryArray
    !
end module DataEntryArray
