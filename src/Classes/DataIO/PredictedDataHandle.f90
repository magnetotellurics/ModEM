!
!
!
module PredictedDataHandle
    !
    use Constants
    use FileUnits
    use DataEntryMT
    use DataEntryMT_REF
    use DataEntryCSEM
    !
    type :: PredictedDataHandle_t
        !
        integer                   :: rx_id
        character(:), allocatable :: code, component
        real( kind=prec )         :: period, xyz(3)
        real( kind=prec )         :: real, imaginary
        !
    end type PredictedDataHandle_t
    !
    public :: buildDataHandle, updateDataHandleArray, writeDataHandleArray
    !
contains
    !
    function buildDataHandle( rx_id, code, component, period, xyz, real, imaginary ) result( data_entry_handle )
        implicit none
        !
        integer, intent( in )                   :: rx_id
        character(:), allocatable, intent( in ) :: code, component
        real( kind=prec ), intent( in )         :: period, xyz(3), real, imaginary
        !
        type( PredictedDataHandle_t ) :: data_entry_handle
        !
        data_entry_handle%rx_id       = rx_id
        data_entry_handle%code        = code
        data_entry_handle%component   = component
        data_entry_handle%period      = period
        data_entry_handle%xyz         = xyz
        data_entry_handle%real        = real
        data_entry_handle%imaginary   = imaginary
        !
    end function buildDataHandle
    !
    !
    subroutine updateDataHandleArray( data_entries, new_data )
          implicit none
          !
          type( PredictedDataHandle_t ), allocatable, intent( inout ) :: data_entries(:)
          type( PredictedDataHandle_t ), intent( in ) :: new_data
          !
          !
          type( PredictedDataHandle_t ), allocatable, dimension(:) :: temp_array
          integer :: istat
          !
          if( .NOT. allocated( data_entries )  ) then
                allocate( data_entries(1) )
                data_entries(1) = new_data
          else
                !
                allocate( temp_array( size( data_entries ) + 1 ), STAT = istat )
                temp_array( 1 : size( data_entries ) ) = data_entries
                temp_array( size( data_entries ) + 1 ) = new_data
                data_entries = temp_array
                !
          endif
          !
     end subroutine updateDataHandleArray
     !
     !
    subroutine writeDataHandleArray( data_handles )
        implicit none
        !
        type( PredictedDataHandle_t ), allocatable, intent( inout ) :: data_handles(:)
        type( PredictedDataHandle_t ) :: aux_data_entry
        !
        integer :: i, j, ios
        !
      ! Order by receiver
        do i = 1, size( data_handles ) - 1
            !
            do j = i + 1, size( data_handles )
                !
                if( data_handles(i)%rx_id > data_handles(j)%rx_id ) then
                    aux_data_entry = data_handles(i)
                    data_handles(i) = data_handles(j)
                    data_handles(j) = aux_data_entry
                endif
                !
            enddo
        enddo
        !
        open( ioPredData, file = "predicted_data.dat", action = "write", form = "formatted", position = "append", iostat = ios )
        !
        if( ios == 0 ) then
            !
            do i = 1, size( data_handles )
                !
                write( ioPredData, "(es12.6, A20, f15.3, f15.3, f15.3, f15.3, f15.3, A20, es16.6, es16.6, es16.6)" ) data_handles(i)%period, data_handles(i)%code, R_ZERO, R_ZERO, data_handles(i)%xyz(1), data_handles(i)%xyz(2), data_handles(i)%xyz(3), data_handles(i)%component, data_handles(i)%real, data_handles(i)%imaginary, 1.0
                !
            enddo
            !
            close( ioPredData )
            !
        else
            stop "Error opening predicted_data.dat in writeDataHandleArray"
        end if
        !
    end subroutine writeDataHandleArray
    !
end module PredictedDataHandle
