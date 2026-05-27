! *****************************************************************************
module DataIO

    use dataspace
    use DataIO_ASCII
    use DataIO_HDF5

    implicit none

    public :: read_dataVectorMTX, write_dataVectorMTX, deall_dataFileInfo

    character(len=*), parameter :: DATA_FILE_TYPE_ASCII = 'DATA_FILE_TYPE_ASCII'
    character(len=*), parameter :: DATA_FILE_TYPE_HDF5 = 'DATA_FILE_TYPE_HDF5'

    character(len=512), private :: DATA_INPUT_FTYPE = DATA_FILE_TYPE_ASCII
    character(len=512), private :: DATA_OUTPUT_FTYPE = DATA_FILE_TYPE_HDF5

contains


subroutine write_dataVectorMTX(allData,cfile)

    implicit none

    type(dataVectorMTX_t), intent(in)         :: allData
    character(*), intent(in)                  :: cfile

    write(0,*) 'DATA_OUTPUT_FTYPE: ', trim(DATA_OUTPUT_FTYPE)

    select case(DATA_OUTPUT_FTYPE)
        case (DATA_FILE_TYPE_ASCII)
            call write_z_list(allData, cfile)
        case (DATA_FILE_TYPE_HDF5)
            call write_hdf5_data(allData, cfile)
        case default
    end select


end subroutine write_dataVectorMTX


subroutine read_dataVectorMTX(allData,cfile)

    implicit none

    type(dataVectorMTX_t), intent(inout)   :: allData
    character(*), intent(in)               :: cfile

    write(0,*) 'DATA_INPUT_FTYPE: ', trim(DATA_INPUT_FTYPE)

    select case(DATA_INPUT_FTYPE)
        case (DATA_FILE_TYPE_ASCII)
            call read_z_list(allData, cfile)
        case (DATA_FILE_TYPE_HDF5)
            call read_hdf5_data(allData, cfile)
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
