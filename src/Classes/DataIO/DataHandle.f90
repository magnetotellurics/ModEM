!
!
!
module DataHandle
    !
    use Constants
    use FileUnits
    use DataEntryMT
    use DataEntryMT_REF
    use DataEntryCSEM
    !
	! Global file name for predicted data file
	character(:), allocatable :: predicted_data_file_name
	!
    type :: DataHandle_t
        !
        integer                   :: rx_type
        character(:), allocatable :: code, component
        real( kind=prec )         :: period, xyz(3)
        real( kind=prec )         :: real, imaginary
        !
    end type DataHandle_t
    !
    public :: buildDataHandle, updateDataHandleArray
    !
contains
    !
    function buildDataHandle( rx_type, code, component, period, xyz, real, imaginary ) result( data_handle )
        implicit none
        !
        integer, intent( in )                   :: rx_type
        character(:), allocatable, intent( in ) :: code, component
        real( kind=prec ), intent( in )         :: period, xyz(3), real, imaginary
        !
        type( DataHandle_t ) :: data_handle
        !
        data_handle%rx_type     = rx_type
        data_handle%code        = code
        data_handle%component   = component
        data_handle%period      = period
        data_handle%xyz         = xyz
        data_handle%real        = real
        data_handle%imaginary   = imaginary
        !
    end function buildDataHandle
    !
    !
    subroutine updateDataHandleArray( data_handle_array, new_data )
        implicit none
        !
        type( DataHandle_t ), allocatable, intent( inout ) :: data_handle_array(:)
        type( DataHandle_t ), intent( in ) :: new_data
        !
        !
        type( DataHandle_t ), allocatable :: temp_array(:)
        integer :: istat, narray
        !
        if( .NOT. allocated( data_handle_array ) ) then
            allocate( data_handle_array(1) )
            data_handle_array(1) = new_data
        else
            !
            narray = size( data_handle_array )
            allocate( temp_array( narray + 1 ), STAT = istat )
            !
            temp_array( 1 : narray ) = data_handle_array
            temp_array( narray + 1 ) = new_data
            data_handle_array = temp_array
            !
            deallocate( temp_array )
            !
        endif
        !
    end subroutine updateDataHandleArray
    !
end module DataHandle
