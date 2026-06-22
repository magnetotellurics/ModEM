submodule (ModelSpace) ModelSpaceIO

    use utilities
    use Declaration_MPI

    implicit none

    character(len=*), parameter :: INPUT = 'input'
    character(len=*), parameter :: OUTPUT = 'output'

    character(len=512) :: INPUT_FILE_TYPE = WS_FILE_TYPE
    character(len=512) :: OUTPUT_FILE_TYPE = WS_FILE_TYPE

    ! Reading interfaces
    interface
        module subroutine write_modelParam_binary(m, cfile, comment)
            implicit none
            character(*), intent(in) :: cfile
            type(modelParam_t), intent(in) :: m
            character(len=*), intent(in), optional :: comment
        end subroutine write_modelParam_binary

        module subroutine write_modelParam_ws(m,cfile,comment)
            implicit none
            type(modelParam_t), intent(in) :: m
            character(*), intent(in)             :: cfile
            character(*), intent(in), optional   :: comment
        end subroutine write_modelParam_ws

        module subroutine write_modelParam_mackie(m,cfile)
            implicit none
            type(modelParam_t), intent(in) :: m
            character(*), intent(in)             :: cfile
        end subroutine write_modelParam_mackie

#ifdef HDF5
        module subroutine write_modelParam_hdf5(m,cfile,comment)
            implicit none
            type(modelParam_t), intent(in) :: m
            character(*), intent(in)             :: cfile
            character(*), intent(in), optional   :: comment
        end subroutine write_modelParam_hdf5
#endif
    end interface

    ! Reading Intefaces
    interface
        module subroutine read_modelParam_binary(grid, m, cfile)
            use GridDef, only : grid_t
            implicit none
            type(grid_t), target, intent(inout)  :: grid
            type(modelParam_t), intent(out)      :: m
            character(*), intent(in)  :: cfile
        end subroutine read_modelParam_binary

        module subroutine read_modelParam_ws(grid, airLayers, m, cfile)
            use GridDef, only : grid_t, airLayers_t
            implicit none
            type(grid_t), target, intent(inout) :: grid
            type(airLayers_t), intent(inout) :: airLayers
            type(modelParam_t), intent(out) :: m
            character(*), intent(in) :: cfile
        end subroutine read_modelParam_WS

        module subroutine read_modelParam_mackie(grid, m, cfile)
            implicit none
            type(grid_t), target, intent(inout) :: grid
            type(modelParam_t), intent(out) :: m
            character(*), intent(in) :: cfile
        end subroutine read_modelParam_mackie

#ifdef HDF5
        module subroutine read_modelParam_hdf5(grid,airLayers,m,cfile)
            implicit none
            type(grid_t), target, intent(inout) :: grid
            type(airLayers_t), intent(inout) :: airLayers
            type(modelParam_t), intent(out) :: m
            character(*), intent(in) :: cfile
        end subroutine read_modelParam_HDF5
#endif
    end interface

contains

    subroutine unsupported_modelParamIO_error(ftype_choice, input_or_output)

        implicit none

        character(len=*), intent(in) :: ftype_choice
        character(len=*), intent(in) :: input_or_output

        ! TODO: Finish this error message
        write(0,*) "ERROR: Unsupported ftype choice: '", trim(ftype_choice), "' for ", trim(input_or_output)
        write(0,*) "ERROR: Supported choiced are: "
        write(0,*) "ERROR:   - '", trim(WS_FILE_TYPE), "'"
        write(0,*) "ERROR:   - '", trim(FTRAN_BINARY_FILE_TYPE), "'"
        write(0,*) "ERROR:   - '", trim(MACKIE_FILE_TYPE), "'"
        write(0,*) "ERROR:   - '", trim(HDF5_FILE_TYPE), "' - Note: If compiled with HDF5"
        write(0,*) "ERROR: See documentation for more information"
        call errStop('Unsupported file type for model '//trim(input_or_output))

    end subroutine unsupported_modelParamIO_error

    module subroutine setup_modelParamIO(input_ftype, output_ftype)

        implicit none

        character(len=*), intent(in) :: input_ftype
        character(len=*), intent(in) :: output_ftype

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

    subroutine write_modelParam(m, cfile, comment)

        implicit none

        type(modelParam_t), intent(in):: m
        character(len=*), intent(in) :: cfile
        character(len=*), intent(in), optional :: comment

        if (taskid /= 0) then
            return
        end if

        select case(OUTPUT_FILE_TYPE)
            case (WS_FILE_TYPE)
                call write_modelParam_WS(m, cfile, comment)
            case (HDF5_FILE_TYPE)
                call compiled_with_HDF5_check()
#ifdef HDF5
                call write_modelParam_hdf5(m, cfile, comment)
#endif
            case (FTRAN_BINARY_FILE_TYPE)
                call write_modelParam_binary(m, cfile, comment)
            case (MACKIE_FILE_TYPE)
                call write_modelParam_mackie(m, cfile)
            case default
                call unsupported_modelParamIO_error(OUTPUT_FILE_TYPE, OUTPUT)
        end select

    end subroutine write_modelParam

    subroutine read_modelParam(grid,airLayers,m,cfile)

        implicit none

        type(grid_t), target, intent(inout)  :: grid
        type(airLayers_t), intent(inout) :: airLayers
        type(modelParam_t), intent(out) :: m
        character(*), intent(in)             :: cfile

        select case(INPUT_FILE_TYPE)
            case (WS_FILE_TYPE)
                call read_modelParam_ws(grid, airLayers, m, cfile)
            case (HDF5_FILE_TYPE)
                call compiled_with_HDF5_check()
#ifdef HDF5
                call read_modelParam_hdf5(grid, airLayers, m, cfile)
#endif
            case (FTRAN_BINARY_FILE_TYPE)
                call read_modelParam_binary(grid, m, cfile)
            case (MACKIE_FILE_TYPE)
                call read_modelParam_mackie(grid, m, cfile)
            case default
                call unsupported_modelParamIO_error(INPUT_FILE_TYPE, INPUT)
        end select

    end subroutine read_modelParam

    subroutine readVec_modelParam(grid,nSigma,sigma,header,cfile)

      implicit none

      type(grid_t), target, intent(inout)   :: grid
      integer, intent(in)		:: nSigma
      type(modelParam_t), intent(inout) 	:: sigma(nSigma)
      character(*), intent(out)		:: header
      character(*), intent(in)		:: cfile

    end subroutine readVec_modelParam

    subroutine writeVec_modelParam(nSigma,sigma,header,cfile)

        implicit none

        integer, intent(in)		:: nSigma
        character(*), intent(in)		:: header, cfile
        type(modelParam_t), intent(in)	:: sigma(nSigma)

    end subroutine writeVec_modelParam

end submodule ModelSpaceIO
