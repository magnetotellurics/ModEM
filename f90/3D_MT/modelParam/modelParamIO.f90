module ModelSpaceIO

    use Declaration_MPI
    use utilities
    use ModelSpace
    use ModelParam_IO_WS
    use ModelParam_IO_HDF5
    use ModelParam_IO_Mackie
    use ModelParam_IO_Binary

    implicit none

    character(len=*), parameter :: WS_FILE_TYPE = 'WS_FILE_TYPE'
    character(len=*), parameter :: HDF5_FILE_TYPE = 'HDF5_FILE_TYPE'
    character(len=*), parameter :: FTRAN_BINARY_FILE_TYPE = 'FTRAN_BINARY_FILE_TYPE'
    character(len=*), parameter :: MACKIE_FILE_TYPE = 'MACKIE_FILE_TYPE'

    character(len=*), private, parameter :: INPUT = 'input'
    character(len=*), private, parameter :: OUTPUT = 'output'

    character(len=512), private :: INPUT_FILE_TYPE
    character(len=512), private :: OUTPUT_FILE_TYPE

contains

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

    subroutine unsupported_modelParamIO_error(ftype_choice, input_or_output)

        implicit none

        character(len=*), intent(in) :: ftype_choice
        character(len=*), intent(in) :: input_or_output

        ! TODO: Finish this error message
        write(0,*) 'ERROR: '
        call errStop('Modem not compiled with HDF5')

    end subroutine unsupported_modelParamIO_error

    subroutine setup_modelParamIO(input_ftype, output_ftype)

        implicit none

        character(len=*), intent(in) :: input_ftype
        character(len=*), intent(in) :: output_ftype

        write(0,*) 'Setting modelParamIO input_type to: ', trim(input_ftype)
        write(0,*) 'Setting modelParamIO output_type to: ', trim(output_ftype)


        select case(input_ftype)
            case (WS_FILE_TYPE)
                INPUT_FILE_TYPE = WS_FILE_TYPE 
            case (HDF5_FILE_TYPE)
                call compiled_with_HDF5_check()
                INPUT_FILE_TYPE = HDF5_FILE_TYPE
            case (FTRAN_BINARY_FILE_TYPE)
                INPUT_FILE_TYPE = FTRAN_BINARY_FILE_TYPE
            case (MACKIE_FILE_TYPE)
                INPUT_FILE_TYPE = MACKIE_FILE_TYPE
            case default
                call unsupported_modelParamIO_error(input_ftype, INPUT)
        end select

        select case(output_ftype)
            case (WS_FILE_TYPE)
                OUTPUT_FILE_TYPE = WS_FILE_TYPE
            case (HDF5_FILE_TYPE)
                call compiled_with_HDF5_check()
                OUTPUT_FILE_TYPE = HDF5_FILE_TYPE
            case (FTRAN_BINARY_FILE_TYPE)
                OUTPUT_FILE_TYPE = FTRAN_BINARY_FILE_TYPE
            case (MACKIE_FILE_TYPE)
                OUTPUT_FILE_TYPE = MACKIE_FILE_TYPE
            case default
                call unsupported_modelParamIO_error(output_ftype, OUTPUT)
        end select

    end subroutine setup_modelParamIO

	subroutine write_modelParam(m,cfile,comment)

        implicit none

        type(modelParam_t), intent(in)	   :: m
        character(*), intent(in)             :: cfile
        character(*), intent(in), optional   :: comment

        if (taskid /= 0) then
            return
        end if

        write(0,*) taskid, 'Writing file: ', trim(cfile)

        select case(OUTPUT_FILE_TYPE)
            case (WS_FILE_TYPE)
                call write_modelParam_WS(m, cfile, comment)
            case (HDF5_FILE_TYPE)
                call compiled_with_HDF5_check()
                call write_modelParam_hdf5(m, cfile, comment)
            case (FTRAN_BINARY_FILE_TYPE)
                call write_modelParam_binary(m, cfile, comment)
            case (MACKIE_FILE_TYPE)
                call write_modelParam_mackie(m, cfile)
            case default
                call unsupported_modelParamIO_error(OUTPUT_FILE_TYPE, OUTPUT)
        end select

    end subroutine write_modelParam

    subroutine read_modelParam(grid,airLayers,m,cfile)

        type(grid_t), target, intent(inout)  :: grid
	    type(airLayers_t), intent(inout)	   :: airLayers
        type(modelParam_t), intent(out)	   :: m
        character(*), intent(in)             :: cfile

        write(0,*) "Opening modelParam for file: ", trim(cfile)

        select case(INPUT_FILE_TYPE)
            case (WS_FILE_TYPE)
                call read_modelParam_ws(grid, airLayers, m, cfile)
            case (HDF5_FILE_TYPE)
                call compiled_with_HDF5_check()
                call read_modelParam_hdf5(grid, airLayers, m, cfile)
            case (FTRAN_BINARY_FILE_TYPE)
                call read_modelParam_binary(grid, m, cfile)
            case (MACKIE_FILE_TYPE)
                call read_modelParam_mackie(grid, m, cfile)
            case default
                call unsupported_modelParamIO_error(INPUT_FILE_TYPE, INPUT)
        end select

    end subroutine read_modelParam

    subroutine readVec_modelParam(grid,nSigma,sigma,header,cfile)

      type(grid_t), target, intent(inout)   :: grid
      integer, intent(in)		:: nSigma
      type(modelParam_t), intent(inout) 	:: sigma(nSigma)
      character(*), intent(out)		:: header
      character(*), intent(in)		:: cfile

    end subroutine readVec_modelParam

    subroutine writeVec_modelParam(nSigma,sigma,header,cfile)

        integer, intent(in)		:: nSigma
        character(*), intent(in)		:: header, cfile
        type(modelParam_t), intent(in)	:: sigma(nSigma)

    end subroutine writeVec_modelParam

end module ModelSpaceIO
