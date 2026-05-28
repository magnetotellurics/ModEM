! *****************************************************************************
module DataIO

    use dataspace
    use DataIO_ASCII
    use DataIO_HDF5

    implicit none

    public :: read_dataVectorMTX, write_dataVectorMTX, deall_dataFileInfo

    character(len=*), parameter :: DATA_FILE_TYPE_ASCII = 'ASCII_LIST_FORMAT'
    character(len=*), parameter :: DATA_FILE_TYPE_HDF5 = 'HDF5'

    character(len=512), private :: DATA_INPUT_FTYPE = DATA_FILE_TYPE_ASCII
    character(len=512), private :: DATA_OUTPUT_FTYPE = DATA_FILE_TYPE_ASCII

    character(len=*), private, parameter :: INPUT = 'input'
    character(len=*), private, parameter :: OUTPUT = 'output'

contains

subroutine unsupported_DataIO_error(ftype_choice, input_or_output)

    implicit none

    character(len=*), intent(in) :: ftype_choice
    character(len=*), intent(in) :: input_or_output

    ! TODO: Finish this error message
    write(0,*) 'ERROR: DataIO'
    call errStop('Modem not compiled with HDF5')

end subroutine unsupported_DataIO_error

subroutine setup_DataIO(input_ftype, output_ftype)

    implicit none

    character(len=*), intent(in) :: input_ftype
    character(len=*), intent(in) :: output_ftype

    write(0,*) 'setup_DataIO - ', trim(input_ftype), ' ', trim(output_ftype)

    select case(input_ftype)
        case (DATA_FILE_TYPE_ASCII)
            DATA_INPUT_FTYPE = DATA_FILE_TYPE_ASCII
        case (DATA_FILE_TYPE_HDF5)
            call compiled_with_HDF5_check()
            DATA_INPUT_FTYPE = DATA_FILE_TYPE_HDF5
        case default
            call unsupported_DataIO_error(input_ftype, input)
    end select

    select case(output_ftype)
        case (DATA_FILE_TYPE_ASCII)
            DATA_OUTPUT_FTYPE = DATA_FILE_TYPE_ASCII
        case (DATA_FILE_TYPE_HDF5)
            call compiled_with_HDF5_check()
            DATA_OUTPUT_FTYPE = DATA_FILE_TYPE_HDF5
        case default
            call unsupported_DataIO_error(output, output)
    end select

end subroutine setup_DataIO


subroutine compiled_with_HDF5_check()

    implicit none
#ifdef HDF5
    ! ModEM is built with HDF5
#else
    call HDF5_not_built_error()
#endif

end subroutine compiled_with_HDF5_check

subroutine HDF5_not_built_error()

    implicit none

    ! TODO: Finish this error message
    write(0,*) 'ERROR: ModEM not compiled with HDF5'
    call errStop('Modem not compiled with HDF5')

end subroutine HDF5_not_built_error


subroutine write_dataVectorMTX(allData,cfile)

    implicit none

    type(dataVectorMTX_t), intent(in)         :: allData
    character(*), intent(in)                  :: cfile

    write(0,*) 'DATA_OUTPUT_FTYPE: ', trim(DATA_OUTPUT_FTYPE)

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

    write(0,*) 'DATA_INPUT_FTYPE: ', trim(DATA_INPUT_FTYPE)

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
