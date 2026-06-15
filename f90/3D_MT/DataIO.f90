! *****************************************************************************
module DataIO

    use utilities
    use dataspace

    implicit none

    public :: read_dataVectorMTX, write_dataVectorMTX, deall_dataFileInfo

    character(len=512), private :: DATA_INPUT_FTYPE = DATA_FILE_TYPE_ASCII
    character(len=512), private :: DATA_OUTPUT_FTYPE = DATA_FILE_TYPE_ASCII

    character(len=*), private, parameter :: INPUT = 'input'
    character(len=*), private, parameter :: OUTPUT = 'output'

   type :: data_file_block

       ! this block of information constitutes user preferences about the data format;
       ! there is one entry per each transmitter type and data type... (iTxt,iDt)
       ! if there are multiple data blocks of the same transmitter & data types,
       ! the last value is used.
       character(200) :: info_in_file
       character(20)  :: sign_info_in_file
       integer        :: sign_in_file
       character(20)  :: units_in_file
       real           :: origin_in_file(2)
       real           :: geographic_orientation

      ! these lists contain the indices into the data vector for each data type;
      ! they make it possible to sort the data by receiver for output.
      ! no data denoted by zero index; dimensions (nTx) and (nTx,nRx).
      ! these indices are typically allocated as we read the data file
      integer, pointer, dimension(:)   :: tx_index
      integer, pointer, dimension(:)   :: dt_index
      integer, pointer, dimension(:,:) :: rx_index

      ! some transmitter types and data types don't go together
      logical         :: defined

   end type data_file_block

  type (data_file_block), pointer, save, private, dimension(:,:) :: fileInfo

    ! Reading dataVectorMTX
    interface
        module subroutine read_Z_list(allData, cfile)
            implicit none
            character(*), intent(in) :: cfile
            type (dataVectorMTX_t), intent(inout) :: allData
        end subroutine read_Z_list

#ifdef HDF5
        module subroutine read_hdf5_data(allData, cfile)
            implicit none
            character(*), intent(in) :: cfile
            type (dataVectorMTX_t), intent(inout) :: allData
        end subroutine read_hdf5_data
#endif
    end interface

    ! Writing dataVectorMTX
    interface
        module subroutine write_z_list(allData, cfile)
            implicit none
            character(*), intent(in) :: cfile
            type (dataVectorMTX_t), intent(in) :: allData
        end subroutine write_z_list

#ifdef HDF5
        module subroutine write_hdf5_data(allData, cfile)
            implicit none
            character(*), intent(in) :: cfile
            type (dataVectorMTX_t), intent(in) :: allData
        end subroutine write_hdf5_data
#endif
    end interface

    ! Cleanup - Dallocating
    interface
        module subroutine deall_fileInfo_ascii()
            implicit none
        end subroutine deall_fileInfo_ascii
    end interface

contains

subroutine unsupported_DataIO_error(ftype_choice, input_or_output)

    implicit none

    character(len=*), intent(in) :: ftype_choice
    character(len=*), intent(in) :: input_or_output

    ! TODO: Finish this error message
    write(0,*) "ERROR: DataIO - Unsupported data filetype choice '", trim(ftype_choice), "' for data ", trim(input_or_output)
    write(0,*) "ERRRO: Supported types are: '", trim(DATA_FILE_TYPE_ASCII), "' or '", trim(DATA_FILE_TYPE_HDF5), "'"
    call ModEM_Abort()

end subroutine unsupported_DataIO_error

subroutine setup_DataIO(input_ftype, output_ftype)

    implicit none

    character(len=*), intent(in) :: input_ftype
    character(len=*), intent(in) :: output_ftype
    integer :: i

    select case(trim(input_ftype))
        case (DATA_FILE_TYPE_ASCII)
            DATA_INPUT_FTYPE = DATA_FILE_TYPE_ASCII
        case (DATA_FILE_TYPE_HDF5)
            call compiled_with_HDF5_check()
            DATA_INPUT_FTYPE = DATA_FILE_TYPE_HDF5
        case default
            call unsupported_DataIO_error(input_ftype, input)
    end select

    select case(trim(output_ftype))
        case (DATA_FILE_TYPE_ASCII)
            DATA_OUTPUT_FTYPE = DATA_FILE_TYPE_ASCII
        case (DATA_FILE_TYPE_HDF5)
            call compiled_with_HDF5_check()
            DATA_OUTPUT_FTYPE = DATA_FILE_TYPE_HDF5
        case default
            call unsupported_DataIO_error(output, output)
    end select

end subroutine setup_DataIO

subroutine write_dataVectorMTX(allData,cfile)

    implicit none

    type(dataVectorMTX_t), intent(in)         :: allData
    character(*), intent(in)                  :: cfile

    select case(DATA_OUTPUT_FTYPE)
        case (DATA_FILE_TYPE_ASCII)
            call write_z_list(allData, cfile)
        case (DATA_FILE_TYPE_HDF5)
            call compiled_with_HDF5_check()
#ifdef HDF5
            call write_hdf5_data(allData, cfile)
#endif
        case default
    end select


end subroutine write_dataVectorMTX


subroutine read_dataVectorMTX(allData,cfile)

    implicit none

    type(dataVectorMTX_t), intent(inout)   :: allData
    character(*), intent(in)               :: cfile

    select case(DATA_INPUT_FTYPE)
        case (DATA_FILE_TYPE_ASCII)
            call read_z_list(allData, cfile)
        case (DATA_FILE_TYPE_HDF5)
            call compiled_with_HDF5_check()
#ifdef HDF5
            call read_hdf5_data(allData, cfile)
#endif
        case default
    end select

end subroutine 

subroutine deall_dataFileInfo()

    implicit none

    select case(DATA_OUTPUT_FTYPE)
        case (DATA_FILE_TYPE_ASCII)
            call deall_fileInfo_ascii()
        case (DATA_FILE_TYPE_HDF5)
        case default
    end select

end subroutine deall_dataFileInfo


end module DataIO
