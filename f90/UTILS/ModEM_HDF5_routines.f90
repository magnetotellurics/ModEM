module ModEM_HDF5

    use hdf5
    use h5ds
    use utilities
    use math_constants

    implicit none

    interface ModEM_HDF5_write_dataset
        module procedure ModEM_HDF5_write_dataset_string
        module procedure ModEM_HDF5_write_dataset_int
        module procedure ModEM_HDF5_write_dataset_int_1D
        module procedure ModEM_HDF5_write_dataset_real_double
        module procedure ModEM_HDF5_write_dataset_real_double_1D
        module procedure ModEM_HDF5_write_dataset_real_double_2D
        module procedure ModEM_HDF5_write_dataset_real_double_3D
        module procedure ModEM_HDF5_write_dataset_complex_real
    end interface

    interface ModEM_HDF5_read_dataset
        module procedure ModEM_HDF5_read_dataset_string
        module procedure ModEM_HDF5_read_dataset_int_1D
        module procedure ModEM_HDF5_read_dataset_real_double_1D
        module procedure ModEM_HDF5_read_dataset_real_double_2D
        module procedure ModEM_HDF5_read_dataset_real_double_3D
        module procedure ModEM_HDF5_read_dataset_complex_real
    end interface

    interface ModEM_HDF5_add_attr
        module procedure ModEM_HDF5_add_attr_string
        module procedure ModEM_HDF5_add_attr_int
        module procedure ModEM_HDF5_add_attr_int_1D_HID
        module procedure ModEM_HDF5_add_attr_real
        module procedure ModEM_HDF5_add_attr_real_double
        module procedure ModEM_HDF5_add_attr_real_1D
        module procedure ModEM_HDF5_add_attr_real_double_2D
    end interface

    interface ModEM_HDF5_read_attr
        module procedure ModEM_HDF5_read_attr_string
        module procedure ModEM_HDF5_read_attr_int
        module procedure ModEM_HDF5_read_attr_real
        module procedure ModEM_HDF5_read_attr_real_double
        module procedure ModEM_HDF5_read_attr_real_1D
        module procedure ModEM_HDF5_read_attr_real_double_2D
    end interface

    interface ModEM_HDF5_attach_dim
        module procedure ModEM_HDF5_attach_dim_by_name
        module procedure ModEM_HDF5_attach_dim_by_id
        module procedure ModEM_HDF5_attach_dim_by_names
    end interface

    interface ModEM_HDF5_create_scalar_variable
        module procedure ModEM_HDF5_create_scalar_variable_real_double
        module procedure ModEM_HDF5_create_scalar_variable_logical
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

    character (len=*), parameter :: DIMENSION_SCALE_CLASS = 'DIMENSION_SCALE'
    character (len=*), parameter :: NETCDF_PURE_DIMENSION_ATTR_VALUE = 'This is a netCDF dimension but not a netCDF variable.'

    integer :: FILE_CREATE_TRUNC_MODE
    integer :: FILE_READ_ONLY_MODE 

    integer :: MODEM_HID_T = HID_T

contains

subroutine ModEM_HDF5_init()

    implicit none

    integer :: hdferr

    ! Initalize HDF5 Interface
    call h5open_f(hdferr)
    if (hdferr < 0) then
        write(0,*) "ERROR: HDF5 Error on h5open_f (initalization)"
        call h5eprint_f(h5e_default_f, hdferr)
        call ModEM_Abort()
    end if

    FILE_CREATE_TRUNC_MODE = H5F_ACC_TRUNC_F
    FILE_READ_ONLY_MODE = H5F_ACC_RDONLY_F

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

subroutine ModEM_HDF5_get_version(version_str)

    implicit none

    character(*), intent(out) :: version_str

    integer :: majnum, minnum, relnum, error

    call h5get_libversion_f(majnum, minnum, relnum, error)
    write(version_str, '(I0,".",I0,".",I0)') majnum, minnum, relnum

end subroutine ModEM_HDF5_get_version

subroutine ModEM_HDF5_open(fname, file_id, mode, hdferr)

    implicit none

    character (len=*), intent(in) :: fname
    integer (kind=HID_T), intent(out) :: file_id
    integer, intent(in) :: mode 
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

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

subroutine ModEM_HDF5_open_file(fname, mode, file_id, hdferr)

    implicit none

    character (len=*), intent(in) :: fname
    integer, intent(in) :: mode
    integer (kind=HID_T), intent(out) :: file_id
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)
    write(0,*) "I am inside open!"

    ! call create
    call h5fopen_f(trim(fname), mode, file_id, hdferr_lcl)
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

    write(0,*) "done with open.."


end subroutine ModEM_HDF5_open_file

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
    call h5lexists_f(file_id, group_name, exists, hdferr)

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

! Convenience function to help us create NETCDF4-CLASSIC pure dimensions.
! A pure dimension is a dimension that only appears in the dimension of a 
! NetCDF file and not inside of the 'variables' section. 
subroutine ModEM_HDF5_create_pure_dimension(loc_id, dim_name, size, dim_id, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: dim_name
    integer (kind=HSIZE_T), intent(in) :: size
    integer (kind=HID_T), optional, intent(out) :: dim_id
    integer, optional, intent(out) :: hdferr

    integer (kind=HID_T) :: dim_dspace_id, dim_dset_id

    call ModEM_HDF5_create_dataspace(1, (/int(size, kind=HSIZE_T)/), dim_dspace_id)
    ! TODO: Error checking...

    ! It doesn't actually matter what 'type' we use here (I think)
    call ModEM_HDF5_create_dataset(loc_id, trim(dim_name), H5T_NATIVE_INTEGER, dim_dspace_id, dim_dset_id)
    ! TODO: Error checking..

    ! Now make the dimension a NETCDF4-CLASSIC dimension by calling our helper functions
    ! note that the below call also calls `ModEM_HDF5_make_dataset_dimension`. These 
    ! simply set HDF5 attributes that tell NetCDF readers that its a dimensinon
    call ModEM_HDF5_make_dataset_pure_dimension(dim_dset_id)
    ! TODO: Error checking
    call ModEM_HDF5_close_dataspace(dim_dspace_id)
    ! TODO: Error checking
    call ModEM_HDF5_close_dataset(dim_dset_id)
    ! TODO: Error checking

end subroutine ModEM_HDF5_create_pure_dimension


! Make an HDF5 dataset a dimemsion in the eyes of NetCDF4. 
subroutine ModEM_HDF5_make_dataset_dimension(dset_id, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: dset_id
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    integer (kind=HID_T) :: type_id, class_attr_id, class_dspace_id
    integer (kind=size_t) :: str_len

    integer(hsize_t), dimension(1) :: dim = (/ 1 /)

    raise_error = present(hdferr)

    ! CLASS must be a null-terminated C string. A space-padded Fortran string
    ! makes h5dsattach_scale_f leak an attribute handle, which keeps the file
    ! open after h5fclose_f.
    call h5tcopy_f(H5T_C_S1, type_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when copying type in ModEM_HDF5_make_dataset_dimension"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    ! Length of the string plus one for the null terminator
    str_len = len(DIMENSION_SCALE_CLASS) + 1
    call h5tset_size_f(type_id, str_len, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else
            write(0,*) "ERROR: HDF5 Error when setting type size in ModEM_HDF5_make_dataset_dimension"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5tset_strpad_f(type_id, H5T_STR_NULLTERM_F, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else
            write(0,*) "ERROR: HDF5 Error when setting strpad in ModEM_HDF5_make_dataset_dimension"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5screate_f(H5S_SCALAR_F, class_dspace_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else
            write(0,*) "ERROR: HDF5 Error when creating scalar dataspace in ModEM_HDF5_make_dataset_dimension"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5acreate_f(dset_id, "CLASS", type_id, class_dspace_id, class_attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else
            write(0,*) "ERROR: HDF5 Error when creating CLASS attribute in ModEM_HDF5_make_dataset_dimension"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5awrite_f(class_attr_id, type_id, DIMENSION_SCALE_CLASS, dim, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else
            write(0,*) "ERROR: HDF5 Error when writing CLASS attribute in ModEM_HDF5_make_dataset_dimension"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5aclose_f(class_attr_id, hdferr_lcl)
    call h5sclose_f(class_dspace_id, hdferr_lcl)
    call h5tclose_f(type_id, hdferr_lcl)

end subroutine ModEM_HDF5_make_dataset_dimension

subroutine ModEM_HDF5_make_dataset_pure_dimension(dset_id, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: dset_id
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    integer (kind=HID_T) :: type_id, name_dspace_id, name_attr_id
    integer (kind=size_t) :: str_len

    integer(hsize_t), dimension(1) :: dim = (/ 1 /)

    raise_error = present(hdferr)

    ! Assuming this correctly adds CLASS="DIMENSION_SCALE"
    call ModEM_HDF5_make_dataset_dimension(dset_id)

    ! Use H5T_C_S1 instead of NATIVE_CHARACTER for strict NetCDF string compatibility
    call h5tcopy_f(H5T_C_S1, type_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when copying type in ModEM_HDF5_make_dataset_pure_dimension"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    ! Length of the magic string plus one for the null terminator
    str_len = len_trim(NETCDF_PURE_DIMENSION_ATTR_VALUE) + 1
    call h5tset_size_f(type_id, str_len, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when setting type size in ModEM_HDF5_make_dataset_pure_dimension"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    ! Explicitly tell HDF5/NetCDF that the string pad is NULLTERM
    call h5tset_strpad_f(type_id, H5T_STR_NULLTERM_F, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when setting strpad in ModEM_HDF5_make_dataset_pure_dimension"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    ! Create SCALAR dataspace (Crucial for NetCDF3/Classic)
    call h5screate_f(H5S_SCALAR_F, name_dspace_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when creating scalar dataspace in ModEM_HDF5_make_dataset_pure_dimension"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    ! Create the "NAME" attribute
    call h5acreate_f(dset_id, "NAME", type_id, name_dspace_id, name_attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when creating attribute in ModEM_HDF5_make_dataset_pure_dimension"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    ! Write the raw string (do NOT append c_null_char manually)
    call h5awrite_f(name_attr_id, type_id, trim(NETCDF_PURE_DIMENSION_ATTR_VALUE), dim, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when writing attribute in ModEM_HDF5_make_dataset_pure_dimension"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    ! Cleanup
    call h5aclose_f(name_attr_id, hdferr_lcl)
    call h5sclose_f(name_dspace_id, hdferr_lcl)
    call h5tclose_f(type_id, hdferr_lcl)

end subroutine ModEM_HDF5_make_dataset_pure_dimension

subroutine ModEM_HDF5_attach_dim_by_name(loc_id, dim_name, dset_id, dim, hdferr)
    
    integer (kind=HID_T), intent(in) :: loc_id
    character(len=*), intent(in) :: dim_name
    integer (kind=HID_T), intent(in) :: dset_id
    integer, intent(in) :: dim
    integer, optional, intent(out) :: hdferr

    integer (kind=HID_T) :: dim_dset_id
    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

    call h5dopen_f(loc_id, dim_name, dim_dset_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when opening dimension dataset in ModEM_HDF5_attach_dim"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5dsattach_scale_f(dset_id, dim_dset_id, dim, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when attaching dimension in ModEM_HDF5_attach_dim"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5dclose_f(dim_dset_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        write(0,*) 'ERROR WAS NOT CLOSED!'
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when closing dimension dataset in ModEM_HDF5_attach_dim_by_name"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_attach_dim_by_name

subroutine ModEM_HDF5_attach_dim_by_names(dim_gid, dim_name, var_gid, var_name, dim, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: dim_gid
    character (len=*), intent(in) :: dim_name
    integer (kind=HID_T), intent(in) :: var_gid
    character (len=*), intent(in) :: var_name
    integer, intent(in) :: dim
    integer, optional, intent(out) :: hdferr

    integer (kind=HID_T) :: dim_dset_id
    integer (kind=HID_T) :: var_dset_id

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

    ! First, get the dimension dset id
    call h5dopen_f(dim_gid, dim_name, dim_dset_id, hdferr_lcl)
    ! TODO: Error checking..

    ! Next, get the dimension dset id
    call h5dopen_f(var_gid, var_name, var_dset_id, hdferr_lcl)
    ! TODO: Error checking

    ! Finally, attach the dimensions
    call h5dsattach_scale_f(var_dset_id, dim_dset_id, dim, hdferr_lcl)
    ! TODO: Error checking

    call ModEM_HDF5_close_dataset(dim_dset_id)
    call ModEM_HDF5_close_dataset(var_dset_id)

end subroutine ModEM_HDF5_attach_dim_by_names

subroutine ModEM_HDF5_attach_dim_by_id(dim_id, dset_id, dim, hdferr)
    
    integer (kind=HID_T), intent(in) :: dim_id 
    integer (kind=HID_T), intent(in) :: dset_id
    integer, intent(in) :: dim
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

    call h5dsattach_scale_f(dset_id, dim_id, dim, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when attaching dimension in ModEM_HDF5_attach_dim_by_id"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_attach_dim_by_id


subroutine ModEM_HDF5_get_dataset_type(dset_id, type_id, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: dset_id
    integer (kind=HID_T), intent(out) :: type_id
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

    call h5dget_type_f(dset_id, type_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when getting dataset type in ModEM_HDF5_get_dataset_type"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_get_dataset_type

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

subroutine ModEM_HDF5_create_scalar_dataspace(dspace_id, hdferr)

    implicit none

    integer (kind=HID_T), intent(out) :: dspace_id
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

    call h5screate_f(H5S_SCALAR_F, dspace_id, hdferr_lcl)
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

end subroutine ModEM_HDF5_create_scalar_dataspace

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
            write(0,*) "ERROR: HDF5 Error when closing data space in ModEM_HDF5_close_dataspace"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_close_dataspace

subroutine ModEM_HDF5_close_dataspaces(dspace_ids, hdferr)

    integer (kind=HID_T), dimension(:), intent(in) :: dspace_ids
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    integer :: i

    raise_error = present(hdferr)

    do i = 1, size(dspace_ids)
        call ModEM_HDF5_close_dataspace(dspace_ids(i), hdferr_lcl)
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when closing data space in ModEM_HDF5_close_dataspaces for the ", i, " dataspace"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end do

end subroutine ModEM_HDF5_close_dataspaces

subroutine ModEM_HDF5_create_string_type(type_id, str_len, hdferr)

    implicit none

    integer (kind=HID_T), intent(out) :: type_id
    integer (kind=HSIZE_T), intent(in) :: str_len
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

    call h5tcopy_f(H5T_FORTRAN_S1, type_id, hdferr_lcl)
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

subroutine ModEM_HDF5_close_type(type_id, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: type_id
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

    call h5tclose_f(type_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else
            write(0,*) "ERROR: HDF5 Error when closing type in MOdME_HDF5_close_type"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_close_type

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

    call h5dcreate_f(loc_id, trim(dataset_name), type, dspace_id, dset_id, hdferr_lcl)
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

subroutine ModEM_HDF5_create_complex_dataset(loc_id, real_dataset_name, imag_dataset_name, dspace_id, &
        real_dset_id, imag_dset_id, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: real_dataset_name
    character (len=*), intent(in) :: imag_dataset_name
    integer (kind=HID_T), intent(in) :: dspace_id
    integer (kind=HID_T), intent(out) :: real_dset_id
    integer (kind=HID_T), intent(out) :: imag_dset_id
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    raise_error = present(hdferr)

    hdferr_lcl = 0

    write(0,*) "Creating real part of complex dataset: ", trim(real_dataset_name)
    call ModEM_HDF5_create_dataset(loc_id, real_dataset_name, H5T_NATIVE_DOUBLE, dspace_id, real_dset_id)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call ModEM_HDF5_create_dataset(loc_id, imag_dataset_name, H5T_NATIVE_DOUBLE, dspace_id, imag_dset_id)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when creating imaginary part of complex data set in ModEM_HDF5_create_complex_dataset"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_create_complex_dataset

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

subroutine ModEM_HDF5_add_string_variable(loc_id, name, buf, hdferr)

    use iso_c_binding, only : c_size_t

    implicit none

    integer (kind=HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: name
    character (len=*), intent(in) :: buf
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    integer (kind=HID_T) :: vspace_id, vtype_id, vdset_id
    integer (HSIZE_T) :: admins(1)

    character(kind=c_char), dimension(:), allocatable :: str_value_c

    raise_error = present(hdferr)

    admins = [1_hsize_t]
    call h5screate_f(H5S_SCALAR_F, vspace_id, hdferr_lcl)
    call h5tcopy_f(H5T_C_S1, vtype_id, hdferr_lcl)
    allocate(str_value_c(len_trim(buf)+1))
    str_value_c = trim(buf) // c_null_char

    call h5tset_size_f(vtype_id, int(1, kind=c_size_t), hdferr_lcl)
    call h5tset_strpad_f(vtype_id, H5T_STR_NULLTERM_F, hdferr_lcl)
    write(0,*) 'str_value_c size: ', len(str_value_c), str_value_c
    call h5screate_simple_f(1, (/len(buf, kind=HSIZE_T)+1/), vspace_id, hdferr_lcl)
    call h5dcreate_f(loc_id, name, vtype_id, vspace_id, vdset_id, hdferr_lcl)
    call h5dwrite_f(vdset_id, vtype_id, buf, (/1_HSIZE_T/), hdferr_lcl)

    deallocate(str_value_c)

    call ModEM_HDF5_close_dataset(vdset_id)
    call ModEM_HDF5_close_dataspace(vspace_id)
    call h5tclose_f(vtype_id, hdferr_lcl)

end subroutine ModEM_HDF5_add_string_variable


subroutine ModEM_HDF5_create_scalar_variable_real_double(loc_id, var_name, var_value, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: var_name
    real (kind=prec), pointer, intent(in) :: var_value
    integer, optional, intent(out) :: hdferr

    integer (kind=HID_T) :: var_dspace_id, var_dset_id

    call ModEM_HDF5_create_scalar_dataspace(var_dspace_id)
    ! TODO: Error checking

    call ModEM_HDF5_create_dataset(loc_id, var_name, H5T_NATIVE_DOUBLE, var_dspace_id, var_dset_id)
    ! TODO: Error checking
    
    call ModEM_HDF5_write_dataset(var_dset_id, var_value)
    ! TODO: Error checking

    call ModEM_HDF5_close_dataset(var_dset_id)
    ! TODO: Error checking

    call ModEM_HDf5_close_dataspace(var_dspace_id)
    ! TODO: Error checking

end subroutine ModEM_HDF5_create_scalar_variable_real_double

subroutine ModEM_HDf5_create_scalar_variable_logical(loc_id, var_name, var_value, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: var_name
    logical, intent(in) :: var_value
    integer, optional, intent(out) :: hdferr

    integer (kind=HID_T) :: var_dspace_id, var_dset_id

    call ModEM_HDF5_create_scalar_dataspace(var_dspace_id)
    ! TODO: Error checking

    call ModEM_HDF5_create_dataset(loc_id, var_name, H5T_NATIVE_INTEGER, var_dspace_id, var_dset_id)
    ! TODO: Error checking
    
    if (var_value) then
        call ModEM_HDF5_write_dataset(var_dset_id, 1)
    else 
        call ModEM_HDF5_write_dataset(var_dset_id, 0)
    end if
    ! TODO: Error checking

    call ModEM_HDF5_close_dataset(var_dset_id)
    ! TODO: Error checking

    call ModEM_HDf5_close_dataspace(var_dspace_id)
    ! TODO: Error checking

end subroutine ModEM_HDF5_create_scalar_variable_logical

subroutine ModEM_HDF5_write_dataset_cptr(dset_id, type, buf, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: dset_id
    integer (kind=HID_T), intent(in) :: type
    type (c_ptr), value :: buf
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

subroutine ModEM_HDF5_write_dataset_real_double(dset_id, buf, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    integer (kind=HID_T), intent(in) :: dset_id
    real (kind=prec), pointer, intent(in) :: buf
    integer, optional, intent(out) :: hdferr

    integer (kind=HID_T) :: type

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: buf_ptr

    raise_error = present(hdferr)

    buf_ptr = c_loc(buf)
    call ModEM_HDF5_write_dataset_cptr(dset_id, H5T_NATIVE_DOUBLE, buf_ptr, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        write(0,*) "ERROR!"
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when writing dataset set in ModEM_HDF5_write_dataset_real_double"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if


end subroutine ModEM_HDF5_write_dataset_real_double

subroutine ModEM_HDF5_write_dataset_real_double_1D(dset_id, buf, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    integer (kind=HID_T), intent(in) :: dset_id
    real (kind=prec), pointer, dimension(:), intent(in) :: buf
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: buf_ptr

    raise_error = present(hdferr)

    buf_ptr = c_loc(buf(1))
    call ModEM_HDF5_write_dataset_cptr(dset_id, H5T_NATIVE_DOUBLE, buf_ptr, hdferr_lcl)
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

subroutine ModEM_HDF5_write_dataset_real_double_2D(dset_id, buf, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    integer (kind=HID_T), intent(in) :: dset_id
    real (kind=prec), pointer, dimension(:,:), intent(in) :: buf
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: buf_ptr

    raise_error = present(hdferr)

    buf_ptr = c_loc(buf(1, 1))
    call ModEM_HDF5_write_dataset_cptr(dset_id, H5T_NATIVE_DOUBLE, buf_ptr, hdferr_lcl)
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

subroutine ModEM_HDF5_write_dataset_int(dset_id, buf, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    integer (kind=HID_T), intent(in) :: dset_id
    integer, target, intent(in) :: buf
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: buf_ptr

    raise_error = present(hdferr)

    buf_ptr = c_loc(buf)
    call ModEM_HDF5_write_dataset_cptr(dset_id, H5T_NATIVE_INTEGER, buf_ptr, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when writing dataset set in ModEM_HDF5_write_dataset_int"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_write_dataset_int


subroutine ModEM_HDF5_write_dataset_int_1D(dset_id, buf, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    integer (kind=HID_T), intent(in) :: dset_id
    integer, pointer, dimension(:), intent(in) :: buf
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: buf_ptr

    raise_error = present(hdferr)

    buf_ptr = c_loc(buf(1))
    call ModEM_HDF5_write_dataset_cptr(dset_id, H5T_NATIVE_INTEGER, buf_ptr, hdferr_lcl)
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

subroutine ModEM_HDF5_write_dataset_real_double_3D(dset_id, buf, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    integer (kind=HID_T), intent(in) :: dset_id
    real (kind=prec), pointer, dimension(:,:,:), intent(in) :: buf
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: buf_ptr

    raise_error = present(hdferr)

    buf_ptr = c_loc(buf(1, 1, 1))
    call ModEM_HDF5_write_dataset_cptr(dset_id, H5T_NATIVE_DOUBLE, buf_ptr, hdferr_lcl)
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

end subroutine ModEM_HDF5_write_dataset_real_double_3D

subroutine ModEM_HDF5_write_dataset_complex_real(real_dset_id, imag_dset_id, buf, memspace_id, fspace_id, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    integer (kind=HID_T), intent(in) :: real_dset_id
    integer (kind=HID_T), intent(in) :: imag_dset_id
    complex (kind=prec), pointer, dimension(:,:,:), intent(in) :: buf
    integer (kind=HID_T), intent(in) :: memspace_id
    integer (kind=HID_T), intent(in) :: fspace_id
    integer, optional, intent(out) :: hdferr

    integer (kind=HSIZE_T) :: start(4), count(4), block(1), stride(1)

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: buf_ptr

    raise_error = present(hdferr)

    start = [0_HSIZE_T, 0_HSIZE_T, 0_HSIZE_T, 0_HSIZE_T]
    count = [1_HSIZE_T, size(buf, 1, kind=HSIZE_T), size(buf, 2, kind=HSIZE_T), size(buf, 3, kind=HSIZE_T)]
    buf_ptr = c_loc(buf(1, 1, 1))

    call h5sselect_hyperslab_f(memspace_id, H5S_SELECT_SET_F, start, count, hdferr_lcl)
    call h5dwrite_f(real_dset_id, H5T_NATIVE_DOUBLE, buf_ptr, hdferr_lcl, memspace_id, fspace_id)

    start = [1_HSIZE_T, 0_HSIZE_T, 0_HSIZE_T, 0_HSIZE_T]

    call h5sselect_hyperslab_f(memspace_id, H5S_SELECT_SET_F, start, count, hdferr_lcl)
    call h5dwrite_f(imag_dset_id, H5T_NATIVE_DOUBLE, buf_ptr, hdferr_lcl, memspace_id, fspace_id)

end subroutine ModEM_HDF5_write_dataset_complex_real

subroutine ModEM_HDF5_read_dataset_cptr(dset_id, type, buf, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: dset_id
    integer (kind=HID_T), intent(in) :: type
    type (c_ptr), value  :: buf
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

subroutine ModEM_HDF5_read_dataset_string(dset_id, type, buf, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    implicit none

    integer (kind=HID_T), intent(in) :: dset_id
    integer (kind=HID_T), intent(in) :: type
    character(len=*), pointer, dimension(:), intent(out) :: buf
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
            write(0,*) "ERROR: HDF5 Error when reading dataset set in ModEM_HDF5_read_dataset_string"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_read_dataset_string

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

subroutine ModEM_HDF5_read_dataset_int(dset_id, type, buf, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: dset_id
    integer (kind=HID_T), intent(in) :: type
    integer, target :: buf
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: buf_ptr

    raise_error = present(hdferr)

    buf_ptr = c_loc(buf)
    call ModEM_HDF5_read_dataset_cptr(dset_id, type, buf_ptr, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else
            write(0,*) "ERROR: HDF5 Error when reading dataset set in ModEM_HDF5_read_dataset_int"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_read_dataset_int

subroutine ModEM_HDF5_read_dataset_int_1D(dset_id, type, buf, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: dset_id
    integer (kind=HID_T), intent(in) :: type
    integer, dimension(:), target :: buf
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
            write(0,*) "ERROR: HDF5 Error when reading dataset set in ModEM_HDF5_read_dataset_int_1D"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_read_dataset_int_1D

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

subroutine ModEM_HDF5_read_dataset_real_double_3D(dset_id, type, buf, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: dset_id
    integer (kind=HID_T), intent(in) :: type
    real (kind=prec), dimension(:,:,:), target :: buf
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: buf_ptr

    raise_error = present(hdferr)

    buf_ptr = c_loc(buf(1,1,1))
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

end subroutine ModEM_HDF5_read_dataset_real_double_3D

subroutine ModEM_HDF5_read_dataset_complex_real(real_dset_id, imag_dset_id, buf, memspace_id, fspace_id, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    implicit none

    integer (kind=HID_T), intent(in) :: real_dset_id
    integer (kind=HID_T), intent(in) :: imag_dset_id
    complex (kind=prec), pointer, dimension(:,:,:), intent(in) :: buf
    integer (kind=HID_T), intent(in) :: memspace_id
    integer (kind=HID_T), intent(in) :: fspace_id
    integer, optional, intent(out) :: hdferr

    integer (kind=HSIZE_T) :: start(4), count(4), block(1), stride(1)

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: buf_ptr

    start = [0_HSIZE_T, 0_HSIZE_T, 0_HSIZE_T, 0_HSIZE_T]
    count = [1_HSIZE_T, size(buf, 1, kind=HSIZE_T), size(buf, 2, kind=HSIZE_T), size(buf, 3, kind=HSIZE_T)]
    buf_ptr = c_loc(buf(1, 1, 1))

    call h5sselect_hyperslab_f(memspace_id, H5S_SELECT_SET_F, start, count, hdferr_lcl)
    call h5dread_f(real_dset_id, H5T_NATIVE_DOUBLE, buf_ptr, hdferr_lcl, memspace_id, fspace_id)

    start = [1_HSIZE_T, 0_HSIZE_T, 0_HSIZE_T, 0_HSIZE_T]

    call h5sselect_hyperslab_f(memspace_id, H5S_SELECT_SET_F, start, count, hdferr_lcl)
    call h5dread_f(imag_dset_id, H5T_NATIVE_DOUBLE, buf_ptr, hdferr_lcl, memspace_id, fspace_id)

end subroutine ModEM_HDF5_read_dataset_complex_real

!subroutine ModEM_HDF5_write_dataset_complex_real(real_dset_id, imag_dset_id, buf, memspace_id, fspace_id, hdferr)


subroutine ModEM_HDF5_add_attr_string(loc_id, attr_name, attr_value, hdferr)

    use iso_c_binding, only : c_null_char

    implicit none

    integer (kind=HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: attr_name
    character (len=*), intent(in) :: attr_value
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    integer (kind=HID_T) :: aspace_id, atype_id, attr_id
    integer (HSIZE_T) :: admins(1)
    character(len=:), allocatable :: attr_value_c

    raise_error = present(hdferr)

    ! Create Fortran string type that matches the size of our string
    admins = [1_hsize_t]
    call h5screate_f(H5S_SCALAR_F, aspace_id, hdferr_lcl)
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

    call h5tcopy_f(H5T_C_S1, atype_id, hdferr_lcl)
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

    attr_value_c = trim(attr_value) // c_null_char

    call h5tset_size_f(atype_id, len(attr_value_c, kind=HSIZE_T), hdferr_lcl)
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

    call h5tset_strpad_f(atype_id, H5T_STR_NULLTERM_F, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else
            write(0,*) "ERROR: HDF5 Error when calling h5tset_strpad_f in ModEM_HDF5_add_attr_string"
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
    call h5awrite_f(attr_id, atype_id, attr_value_c, admins, hdferr_lcl)
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

    call h5tclose_f(atype_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else
            write(0,*) "ERROR: HDF5 Error when calling h5tclose_f in ModEM_HDF5_add_attr_string"
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

subroutine ModEM_HDF5_add_attr_int_1D_HID(loc_id, attr_name, attr_value, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    implicit none

    integer (kind=HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: attr_name
    integer (kind=HID_T), dimension(:), target :: attr_value 
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

end subroutine ModEM_HDF5_add_attr_int_1D_HID

subroutine ModEM_HDF5_add_attr_real(loc_id, attr_name, attr_value, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    implicit none

    integer (kind=HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: attr_name
    real, target :: attr_value
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
            write(0,*) "ERROR: HDF5 Error when calling h5screate_simple_f in ModEM_HDF5_add_attr_real"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    call h5acreate_f(loc_id, attr_name, H5T_NATIVE_REAL, aspace_id, attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else
            write(0,*) "ERROR: HDF5 Error when calling h5acreate_f in ModEM_HDF5_add_attr_real"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    attr_value_ptr = c_loc(attr_value)

    call h5awrite_f(attr_id, H5T_NATIVE_REAL, attr_value_ptr, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else
            write(0,*) "ERROR: HDF5 Error when calling h5awrite_f in ModEM_HDF5_add_attr_real"
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
            write(0,*) "ERROR: HDF5 Error when calling h5aclose_f in ModEM_HDF5_add_attr_real"
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
            write(0,*) "ERROR: HDF5 Error when calling h5sclose_f in ModEM_HDF5_add_attr_real"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_add_attr_real

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

function ModEM_HDF5_attr_exists(loc_id, attr_name, hdferr) result(is_present)

    implicit none

    integer (kind=HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: attr_name
    integer, optional, intent(out) :: hdferr

    integer :: hdferr_lcl

    logical :: is_present
    logical :: raise_error

    raise_error = present(hdferr)

    is_present = .false.

    ! Check if the attribute exists
    call h5aexists_f(loc_id, attr_name, is_present, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
        else
            write(0,*) "ERROR: HDF5 Error when calling h5aexists_f in ModEM_HDF5_attr_exists"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end function ModEM_HDF5_attr_exists

subroutine ModEM_HDF5_read_attr_string(loc_id, attr_name, attr_value, hdferr)

    implicit none

    integer (kind=HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: attr_name
    character (len=*), intent(out) :: attr_value
    integer, optional, intent(out) :: hdferr

    logical :: raise_error
    integer :: hdferr_lcl

    integer (kind=HID_T) :: attr_id, atype_id, aspace_id

    raise_error = present(hdferr)

    attr_value = ""

    ! Open the attribute
    call h5aopen_f(loc_id, attr_name, attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5aopen_f in ModEM_HDF5_read_attr_string"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    ! Get the attribute type
    call h5aget_type_f(attr_id, atype_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5aget_type_f in ModEM_HDF5_read_attr_string"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    ! Read the attribute
    call h5aread_f(attr_id, atype_id, attr_value, [1_hsize_t], hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5aread_f in ModEM_HDF5_read_attr_string"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    ! Close the attribute
    call h5aclose_f(attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5aclose_f in ModEM_HDF5_read_attr_string"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_read_attr_string

subroutine ModEM_HDF5_read_attr_int(loc_id, attr_name, attr_value, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    implicit none

    integer (kind=HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: attr_name
    integer, target, intent(out) :: attr_value 
    integer, optional, intent(out) :: hdferr

    integer (kind=HID_T) :: attr_id

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: attr_value_ptr

    raise_error = present(hdferr)

    ! Open the attribute
    call h5aopen_f(loc_id, attr_name, attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5aopen_f in ModEM_HDF5_read_attr_int"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    attr_value_ptr = c_loc(attr_value)

    ! Read the attribute
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, attr_value_ptr, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5aread_f in ModEM_HDF5_read_attr_int"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    ! Close the attribute
    call h5aclose_f(attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5aclose_f in ModEM_HDF5_read_attr_int"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_read_attr_int

subroutine ModEM_HDF5_read_attr_real(loc_id, attr_name, attr_value, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    implicit none

    integer (kind=HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: attr_name
    real, target, intent(out) :: attr_value
    integer, optional, intent(out) :: hdferr

    integer (kind=HID_T) :: attr_id

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: attr_value_ptr

    raise_error = present(hdferr)

    ! Open the attribute
    call h5aopen_f(loc_id, attr_name, attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else
            write(0,*) "ERROR: HDF5 Error when calling h5aopen_f in ModEM_HDF5_read_attr_real"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    attr_value_ptr = c_loc(attr_value)

    ! Read the attribute
    call h5aread_f(attr_id, H5T_NATIVE_REAL, attr_value_ptr, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else
            write(0,*) "ERROR: HDF5 Error when calling h5aread_f in ModEM_HDF5_read_attr_real"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    ! Close the attribute
    call h5aclose_f(attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else
            write(0,*) "ERROR: HDF5 Error when calling h5aclose_f in ModEM_HDF5_read_attr_real"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_read_attr_real

subroutine ModEM_HDF5_read_attr_real_double(loc_id, attr_name, attr_value, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    implicit none

    integer (kind=HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: attr_name
    real(kind=prec), target, intent(out) :: attr_value 
    integer, optional, intent(out) :: hdferr

    integer (kind=HID_T) :: attr_id

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: attr_value_ptr

    raise_error = present(hdferr)

    ! Open the attribute
    call h5aopen_f(loc_id, attr_name, attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5aopen_f in ModEM_HDF5_read_attr_real_double"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    attr_value_ptr = c_loc(attr_value)

    ! Read the attribute
    call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, attr_value_ptr, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5aread_f in ModEM_HDF5_read_attr_real_double"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    ! Close the attribute
    call h5aclose_f(attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5aclose_f in ModEM_HDF5_read_attr_real_double"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_read_attr_real_double

subroutine ModEM_HDF5_read_attr_real_1D(loc_id, attr_name, attr_value, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    implicit none

    integer (kind=HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: attr_name
    real, dimension(:), target, intent(out) :: attr_value 
    integer, optional, intent(out) :: hdferr

    integer (kind=HID_T) :: attr_id

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: attr_value_ptr

    raise_error = present(hdferr)

    ! Open the attribute
    call h5aopen_f(loc_id, attr_name, attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5aopen_f in ModEM_HDF5_read_attr_real_1D"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    attr_value_ptr = c_loc(attr_value(1))

    ! Read the attribute
    call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, attr_value_ptr, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5aread_f in ModEM_HDF5_read_attr_real_1D"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    ! Close the attribute
    call h5aclose_f(attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5aclose_f in ModEM_HDF5_read_attr_real_1D"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_read_attr_real_1D

subroutine ModEM_HDF5_read_attr_real_double_2D(loc_id, attr_name, attr_value, hdferr)

    use iso_c_binding, only : c_loc, c_ptr

    implicit none

    integer (kind=HID_T), intent(in) :: loc_id
    character (len=*), intent(in) :: attr_name
    real (kind=prec), dimension(:,:), pointer, intent(out) :: attr_value 
    integer, optional, intent(out) :: hdferr

    integer (kind=HID_T) :: attr_id

    logical :: raise_error
    integer :: hdferr_lcl

    type (c_ptr) :: attr_value_ptr

    raise_error = present(hdferr)

    ! Open the attribute
    call h5aopen_f(loc_id, attr_name, attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5aopen_f in ModEM_HDF5_read_attr_real_double_2D"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    attr_value_ptr = c_loc(attr_value(1,1))

    ! Read the attribute
    call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, attr_value_ptr, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5aread_f in ModEM_HDF5_read_attr_real_double_2D"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

    ! Close the attribute
    call h5aclose_f(attr_id, hdferr_lcl)
    if (hdferr_lcl /= 0) then
        if (raise_error) then
            hdferr = hdferr_lcl
            return
        else 
            write(0,*) "ERROR: HDF5 Error when calling h5aclose_f in ModEM_HDF5_read_attr_real_double_2D"
            call h5eprint_f(h5e_default_f, hdferr_lcl)
            call ModEM_abort()
        end if
    end if

end subroutine ModEM_HDF5_read_attr_real_double_2D

end module ModEM_HDF5
