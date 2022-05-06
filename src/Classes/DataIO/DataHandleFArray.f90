!*************
!
! Class to provide a dynamic and polymorphic data_handle_array of DataHandle_t objects
!
! Last modified at 04/2022 by Paulo Werdt
!
!*************
!
module DataHandleFArray
    !
    use Constants
    !
    use DataHandleMT
    use DataHandleCSEM
    !
    ! Allocatable DataHandle element of the array, for Old Fortran polymorphism !!!
    type, public :: Dh_t
        !
        class( DataHandle_t ), allocatable :: Dh
        !
    end type Dh_t
    !
    public :: getDataHandle, updateDataHandleArray, deallocateDataHandleArray
    !
contains
    !
    ! Add a new DataHandle_t and initialize it if necessary
    subroutine updateDataHandleArray( data_handle_array, new_dh )
        implicit none
        !
        type( Dh_t ), pointer, dimension(:), intent( inout ) :: data_handle_array
        class( DataHandle_t ), intent( in ) :: new_dh
        !
        integer                             :: iDh, nDh
        type( Dh_t ), pointer, dimension(:) :: temp_array
        type( Dh_t ) :: temp_dh
        !
        if( .NOT. associated( data_handle_array ) ) then
            !
            allocate( data_handle_array( 1 ) )
            !
            temp_dh%Dh = new_dh
            !
            data_handle_array( 1 ) = temp_dh
            !
        else
            !
            nDh = size( data_handle_array )
            !
            allocate( temp_array( nDh + 1 ) )
            !
            temp_array( 1 : nDh ) => data_handle_array(:)
            !
            temp_dh%Dh = new_dh
            !
            allocate( temp_array( nDh + 1 ), source = temp_dh )
            !
            data_handle_array => temp_array
            !
            nullify( temp_array )
            !
        endif
        !
    end subroutine updateDataHandleArray
    !
    subroutine setDataHandle( data_handle_array, iDh, Dh )
        implicit none
        !
        !
        type( Dh_t ), allocatable, dimension(:), intent( inout ) :: data_handle_array
        integer, intent( in )                                    :: iDh
        !
        class( DataHandle_t ), allocatable, intent( in )         :: Dh
        !
        data_handle_array( iDh )%Dh = Dh
        !
    end subroutine setDataHandle
    !
    function getDataHandle( data_handle_array, iDh ) result( dh )
        implicit none
        !
        type( Dh_t ), pointer, dimension(:), intent( in ) :: data_handle_array
        integer                                 :: iDh
        !
        class( DataHandle_t ), allocatable      :: dh
        !
        dh = data_handle_array( iDh )%Dh
        !
    end function getDataHandle
    !
    !
    subroutine deallocateDataHandleArray( data_handle_array )
        implicit none
        !
        type( Dh_t ), pointer, dimension(:), intent( inout ) :: data_handle_array
        !
        integer                    :: ndh, idh
        class( Dh_t ), allocatable :: alloc_dh
        !
        !write( *, * ) "deallocateDataHandleArray:", size( data_handle_array )
        !
        ndh = size( data_handle_array )
        do idh = 1, ndh
            alloc_dh = data_handle_array( idh )
            deallocate( alloc_dh )
        end do
        !
        nullify( data_handle_array )
        !
    end subroutine deallocateDataHandleArray
    !
end module DataHandleFArray
