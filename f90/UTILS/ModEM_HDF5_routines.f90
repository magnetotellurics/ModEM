module ModEM_HDF5

    use hdf5
    use utilities
    use math_constants

    implicit none

    interface ModEM_HDF5_write_dataset
       MODULE PROCEDURE ModEM_HDF5_write_dataset_string
       MODULE PROCEDURE ModEM_HDF5_write_dataset_int_1D
       MODULE PROCEDURE ModEM_HDF5_write_dataset_real_double_1D
       MODULE PROCEDURE ModEM_HDF5_write_dataset_real_double_2D
    end interface

    interface ModEM_HDF5_read_dataset
       MODULE PROCEDURE ModEM_HDF5_read_dataset_real_double_1D
       MODULE PROCEDURE ModEM_HDF5_read_dataset_real_double_2D
    end interface

    interface ModEM_HDF5_add_attr
        MODULE PROCEDURE ModEM_HDF5_add_attr_string
        MODULE PROCEDURE ModEM_HDF5_add_attr_int
        MODULE PROCEDURE ModEM_HDF5_add_attr_real_double
        MODULE PROCEDURE ModEM_HDF5_add_attr_real_1D
        MODULE PROCEDURE ModEM_HDF5_add_attr_real_double_2D
    end interface

    abstract interface
       integer function ModEM_h5_itr_cb_interf(loc_id, name, info, parent_id) bind(c)
           use iso_c_binding
           use hdf5
           implicit none
           integer(hid_t), value :: loc_id
           character(len=1), dimension(1:10) :: name ! must have LEN=1 for bind(C) strings
           type(c_ptr) :: info
           integer (kind=HID_T) :: parent_id
        end function ModEM_h5_itr_cb_interf
    end interface

contains

subroutine ModEM_HDF5_init()

    implicit none

    integer :: hdferr

    write(0,*) 'Calling h5open_f'
    ! Initalize HDF5 Interface
    call h5open_f(hdferr)

    if (hdferr < 0) then
        write(0,*) "ERROR: HDF5 Error on h5open_f (initalization)"
        call h5eprint_f(h5e_default_f, hdferr)
        call ModEM_Abort()
    end if

end subroutine ModEM_HDF5_init


subroutine ModEM_HDF5_finalize()

    implicit none

    integer :: hdferr

#ifdef HDF5
    call h5close_f(hdferr)
    if (hdferr < 0) then
       write(0,*) "ERROR: HDF5 Error on h5close_f (shutdown)"
       call h5eprint_f(h5e_default_f, hdferr)
    end if
#endif

end subroutine ModEM_HDF5_finalize

subroutine ModEM_HDF5_open(fname, file_id, mode, hdferr)

    implicit none

    character (len=*), intent(in) :: fname
    integer (kind=HID_T), intent(out) :: file_id
    integer, intent(in) :: mode 
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    call h5fopen_f(fname, mode, file_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when opening file in ModEM_HDF5_open_file"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_open

subroutine ModEM_HDF5_create_file(fname, mode, file_id, hdferr)

    implicit none

    character (len=*), intent(in) :: fname
    integer, intent(in) :: mode
    integer (kind=HID_T), intent(out) :: file_id
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

    ! call create
    call h5fcreate_f(fname, mode, file_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when creating file in ModEM_HDF5_create_file"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_create_file

subroutine ModEM_HDF5_close_file(file_id, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: file_id
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

    call h5fclose_f(file_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when close file in ModEM_HDF5_close_file"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_close_file

subroutine ModEM_HDF5_open_group(file_id, group_name, group_id, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: file_id
    character (len=*), intent(in) :: group_name
    integer (kind=HID_T), intent(out) :: group_id
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

    call h5gopen_f(file_id, group_name, group_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when opening group in ModEM_HDF5_open_group"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_open_group

subroutine ModEM_HDF5_create_group(file_id, group_name, group_id, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: file_id
    character (len=*), intent(in) :: group_name
    integer (kind=HID_T), intent(out) :: group_id
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

    call h5gcreate_f(file_id, trim(group_name), group_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when creating group in ModEM_HDF5_create_group"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_create_group

subroutine ModEM_HDF5_close_group(group_id, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: group_id
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

    call h5gclose_f(group_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when closing group in ModEM_HDF5_close_group"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_close_group

subroutine ModEM_HDF5_does_group_exist(file_id, group_name, exists, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: file_id
    character (len=*), intent(in) :: group_name
    logical, intent(out) :: exists
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

end subroutine ModEM_HDF5_does_group_exist

subroutine ModEM_HDF5_iterate_group(loc_id, callback_function, hdferr)

    use iso_c_binding, only : c_funloc, C_NULL_PTR

    implicit none

    integer (kind=HID_T), intent(in) :: loc_id
    procedure(ModEM_h5_itr_cb_interf) :: callback_function
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    integer (kind=HSIZE_T) :: idx = 1
    integer :: return_value

    type (c_ptr) :: ptr
    type (c_funptr) :: funptr

    integer(kind=HID_T), target :: pid = 0 ! This should be an input variable probably...

    raise_error = present(hdferr)

    pid = loc_id
    ptr = c_loc(pid)
    idx = 0
    funptr = c_funloc(callback_function)

    call h5literate_f(loc_id, H5_INDEX_NAME_F, H5_ITER_NATIVE_F, idx, funptr, ptr, &
            return_value, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5literate_f in ModEM_HDF5_iterate_group"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_iterate_group
    
subroutine ModEM_HDF5_get_dataspace(dset_id, dspace_id, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: dset_id
    integer (kind=HID_T), intent(out) :: dspace_id
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

    call h5dget_space_f(dset_id, dspace_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when getting dataspace in ModEM_HDF5_get_dataspace"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_get_dataspace

subroutine ModEM_HDF5_get_dataspace_size(dspace_id, npoints, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: dspace_id
    integer (kind=HSIZE_T), intent(out) :: npoints
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

    call h5sget_simple_extent_npoints_f(dspace_id, npoints, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when getting data space size in ModEM_HDF5_get_dataspace_size"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_get_dataspace_size

subroutine ModEM_HDF5_get_dataspace_dims(dspace_id, dims, maxdims, rank, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: dspace_id
    integer(HSIZE_T), allocatable, intent(out) :: dims(:)
    integer(HSIZE_T), allocatable, intent(out) :: maxdims(:)
    integer, intent(out) :: rank
    integer, optional, intent(out) ::	hdferr 
    
    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

    ! Discover dataspace rank first
    call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferr_lcl)

    allocate(dims(rank))
    allocate(maxdims(rank))

    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, hdferr_lcl)
    if (hdferr_lcl < 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when getting data space size in ModEM_HDF5_get_dataspace_size"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_get_dataspace_dims

subroutine ModEM_HDF5_create_dataspace(rank, dims, dspace_id, hdferr)

    implicit none

    integer, intent(in) :: rank
    integer (kind=HSIZE_T), dimension(:), intent(in) :: dims
    integer (kind=HID_T), intent(out) :: dspace_id
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

    call h5screate_simple_f(rank, dims, dspace_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when creating data space in ModEM_HDF5_create_dataspace"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_create_dataspace

subroutine ModEM_HDF5_close_dataspace(dspace_id, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: dspace_id
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

    call h5sclose_f(dspace_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when creating data space in ModEM_HDF5_create_dataspace"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_close_dataspace

subroutine ModEM_HDF5_create_string_type(type_id, str_len, hdferr)

    implicit none

    integer (kind=HID_T), intent(out) :: type_id
    integer (kind=HSIZE_T), intent(in) :: str_len
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

    call h5tcopy_f(H5T_NATIVE_CHARACTER, type_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5tcopy_f in ModEM_HDF5_create_string_type"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5tset_size_f(type_id, str_len, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5tset_size_f in ModEM_HDF5_create_string_type"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if


end subroutine ModEM_HDF5_create_string_type

subroutine ModEM_HDF5_open_dataset(loc_id, dataset_name, dset_id, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: dataset_name
    integer (kind=HID_T), intent(out) :: dset_id
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

    call h5dopen_f(loc_id, dataset_name, dset_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when opening data set in ModEM_HDF5_open_dataset"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_open_dataset

subroutine ModEM_HDF5_create_dataset(loc_id, dataset_name, type, dspace_id, dset_id, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: dataset_name
    integer (kind=HID_T), intent(in) :: type
    integer (kind=HID_T), intent(in) :: dspace_id
    integer (kind=HID_T), intent(out) :: dset_id
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

    call h5dcreate_f(loc_id, dataset_name, type, dspace_id, dset_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when creating data set in ModEM_HDF5_create_dataset"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_create_dataset

subroutine ModEM_HDF5_close_dataset(dset_id, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: dset_id
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)
    call h5dclose_f(dset_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        write(0,*) 'ERROR!'
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when closing data set in ModEM_HDF5_close_dataset"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_close_dataset

subroutine ModEM_HDF5_write_dataset_cptr(dset_id, type, buf, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: dset_id
    integer (kind=HID_T), intent(in) :: type
    type (c_ptr), intent(in) :: buf
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

    call h5dwrite_f(dset_id, type, buf, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when closing data set in ModEM_HDF5_close_dataset"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    hdferr = hdferr_lcl

end subroutine ModEM_HDF5_write_dataset_cptr

subroutine ModEM_HDF5_write_dataset_string(dset_id, type, buf, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    integer (kind=HID_T), intent(in) :: dset_id
    integer (kind=HID_T), intent(in) :: type
    character(len=*), pointer, dimension(:), intent(in) :: buf
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl
    integer (kind=HID_T) :: str_type_id

    type (c_ptr) :: buf_ptr
    
    buf_ptr = c_loc(buf(1))

    call ModEM_HDF5_write_dataset_cptr(dset_id, type, buf_ptr, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error calling ModEM_HDF5_write_dataset_cptr in ModEM_HDF5_write_dataset_string"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_write_dataset_string

subroutine ModEM_HDF5_write_dataset_real_double_1D(dset_id, type, buf, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    integer (kind=HID_T), intent(in) :: dset_id
    integer (kind=HID_T), intent(in) :: type
    real (kind=prec), pointer, dimension(:), intent(in) :: buf
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: buf_ptr

    raise_error = present(hdferr)

    buf_ptr = c_loc(buf(1))
    call ModEM_HDF5_write_dataset_cptr(dset_id, type, buf_ptr, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        write(0,*) "ERROR!"
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when writing dataset set in ModEM_HDF5_write_dataset_real_double_1D"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if


end subroutine ModEM_HDF5_write_dataset_real_double_1D

subroutine ModEM_HDF5_write_dataset_real_double_2D(dset_id, type, buf, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    integer (kind=HID_T), intent(in) :: dset_id
    integer (kind=HID_T), intent(in) :: type
    real (kind=prec), pointer, dimension(:,:), intent(in) :: buf
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: buf_ptr

    raise_error = present(hdferr)

    buf_ptr = c_loc(buf(1, 1))
    call ModEM_HDF5_write_dataset_cptr(dset_id, type, buf_ptr, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when writing dataset set in ModEM_HDF5_close_dataset"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if


end subroutine ModEM_HDF5_write_dataset_real_double_2D

subroutine ModEM_HDF5_write_dataset_int_1D(dset_id, type, buf, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    integer (kind=HID_T), intent(in) :: dset_id
    integer (kind=HID_T), intent(in) :: type
    integer, pointer, dimension(:), intent(in) :: buf
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: buf_ptr

    raise_error = present(hdferr)

    buf_ptr = c_loc(buf(1))
    call ModEM_HDF5_write_dataset_cptr(dset_id, type, buf_ptr, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when writing dataset set in ModEM_HDF5_close_dataset"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_write_dataset_int_1D

subroutine ModEM_HDF5_read_dataset_cptr(dset_id, type, buf, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: dset_id
    integer (kind=HID_T), intent(in) :: type
    type (c_ptr), intent(out) :: buf
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

    call h5dread_f(dset_id, type, buf, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when reading dataset in ModEM_HDF5_read_dataset_cptr"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    hdferr = 0

end subroutine ModEM_HDF5_read_dataset_cptr

subroutine ModEM_HDF5_read_dataset_real_double_1D(dset_id, type, buf, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: dset_id
    integer (kind=HID_T), intent(in) :: type
    real (kind=prec), dimension(:), target :: buf
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: buf_ptr

    raise_error = present(hdferr)

    buf_ptr = c_loc(buf(1))
    call ModEM_HDF5_read_dataset_cptr(dset_id, type, buf_ptr, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when reading dataset set in ModEM_HDF5_read_dataset_1D_real_double"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_read_dataset_real_double_1D

subroutine ModEM_HDF5_read_dataset_real_double_2D(dset_id, type, buf, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: dset_id
    integer (kind=HID_T), intent(in) :: type
    real (kind=prec), dimension(:,:), target :: buf
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: buf_ptr

    raise_error = present(hdferr)

    buf_ptr = c_loc(buf(1,1))
    call ModEM_HDF5_read_dataset_cptr(dset_id, type, buf_ptr, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when reading dataset set in ModEM_HDF5_read_dataset_1D_real_double"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_read_dataset_real_double_2D

subroutine ModEM_HDF5_add_attr_string(loc_id, attr_name, attr_value, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: attr_name
    character (len=*), intent(in) :: attr_value
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    integer (kind=HID_T) :: aspace_id, atype_id, attr_id
    integer (HSIZE_T) :: admins(1)

    raise_error = present(hdferr)

    ! Create Fortran string type that matches the size of our string
    admins = [1_hsize_t]
    call h5screate_simple_f(1, admins, aspace_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5screate_simple_f in ModEM_HDF5_add_attr_string"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5tcopy_f in ModEM_HDF5_add_attr_string", hdferr_lcl
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5tset_size_f(atype_id, len_trim(attr_value, kind=HSIZE_T), hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5tset_size in ModEM_HDF5_add_attr_string"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    ! Create the attribute
    call h5acreate_f(loc_id, attr_name, atype_id, aspace_id, attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5acreate_f in ModEM_HDF5_add_attr_string"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    ! Write the attribute
    call h5awrite_f(attr_id, atype_id, attr_value, admins, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5awrite_f in ModEM_HDF5_add_attr_string"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if


    ! Close the attribute and the attribute space
    call h5aclose_f(attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5aclose_f in ModEM_HDF5_add_attr_string"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5sclose_f(aspace_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5sclose_f in ModEM_HDF5_add_attr_string"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine MOdEM_HDF5_add_attr_string

subroutine ModEM_HDF5_add_attr_int(loc_id, attr_name, attr_value, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    implicit none

    integer (kind=HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: attr_name
    integer, target :: attr_value 
    integer, optional, intent(out) :: hdferr

    integer (kind=HID_T) :: aspace_id, attr_id

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: attr_value_ptr

    raise_error = present(hdferr)

    call h5screate_simple_f(rank(attr_value), shape(attr_value, kind=HSIZE_T), aspace_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5screate_simple_f in ModEM_HDF5_add_attr_int"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5acreate_f(loc_id, attr_name, H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5acreate_f in ModEM_HDF5_add_attr_int"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    attr_value_ptr = c_loc(attr_value)

    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_value_ptr, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5awrite_f in ModEM_HDF5_add_attr_int"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5aclose_f(attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5aclose_f in ModEM_HDF5_add_attr_int"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5sclose_f(aspace_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5sclose_f in ModEM_HDF5_add_attr_int"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_add_attr_int


subroutine ModEM_HDF5_add_attr_real_double(loc_id, attr_name, attr_value, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    implicit none

    integer (kind=HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: attr_name
    real(kind=prec), target :: attr_value 
    integer, optional, intent(out) :: hdferr

    integer (kind=HID_T) :: aspace_id, attr_id

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: attr_value_ptr

    raise_error = present(hdferr)

    call h5screate_simple_f(rank(attr_value), shape(attr_value, kind=HSIZE_T), aspace_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5screate_simple_f in ModEM_HDF5_add_attr_real_double"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5acreate_f(loc_id, attr_name, H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5acreate_f in ModEM_HDF5_add_attr_real_double"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    attr_value_ptr = c_loc(attr_value)

    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_value_ptr, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5awrite_f in ModEM_HDF5_add_attr_real_double"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5aclose_f(attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5aclose_f in ModEM_HDF5_add_attr_real_double"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5sclose_f(aspace_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5sclose_f in ModEM_HDF5_add_attr_real_double"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_add_attr_real_double

subroutine ModEM_HDF5_add_attr_real_1D(loc_id, attr_name, attr_value, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    implicit none

    integer (kind=HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: attr_name
    real, dimension(:), target :: attr_value 
    integer, optional, intent(out) :: hdferr

    integer (kind=HID_T) :: aspace_id, attr_id

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: attr_value_ptr

    raise_error = present(hdferr)

    call h5screate_simple_f(rank(attr_value), shape(attr_value, kind=HSIZE_T), aspace_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5screate_simple_f in ModEM_HDF5_add_attr_real_1D"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5acreate_f(loc_id, attr_name, H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5acreate_f in ModEM_HDF5_add_attr_real_1D"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    attr_value_ptr = c_loc(attr_value)

    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_value_ptr, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5awrite_f in ModEM_HDF5_add_attr_real_1D"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5aclose_f(attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5aclose_f in ModEM_HDF5_add_attr_real_1D"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5sclose_f(aspace_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5sclose_f in ModEM_HDF5_add_attr_real_1D"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_add_attr_real_1D

subroutine ModEM_HDF5_add_attr_real_double_2D(loc_id, attr_name, attr_value, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    implicit none

    integer (kind=HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: attr_name
    real (kind=prec), dimension(:), pointer :: attr_value 
    integer, optional, intent(out) :: hdferr

    integer (kind=HID_T) :: aspace_id, attr_id

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: attr_value_ptr

    raise_error = present(hdferr)

    call h5screate_simple_f(rank(attr_value), shape(attr_value, kind=HSIZE_T), aspace_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5screate_simple_f in ModEM_HDF5_add_attr_real_double_2D"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5acreate_f(loc_id, attr_name, H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5acreate_f in ModEM_HDF5_add_attr_real_double_2D"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    attr_value_ptr = c_loc(attr_value)

    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_value_ptr, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5awrite_f in ModEM_HDF5_add_attr_real_double_2D"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5aclose_f(attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5aclose_f in ModEM_HDF5_add_attr_real_double_2D"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5sclose_f(aspace_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5sclose_f in ModEM_HDF5_add_attr_real_double_2D"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_add_attr_real_double_2D

end module ModEM_HDF5
