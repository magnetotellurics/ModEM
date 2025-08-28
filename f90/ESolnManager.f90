module ESolnManager

    use mpi
    use Declaration_MPI
    use GridDef
    use SolnSpace
    use utilities

    implicit None

    private

    type (solnVectorMTX_t) :: EsMgr_eAll
    type (ModEM_mpi_context_t), pointer :: EsMgr_ctx

    character(len=*), parameter :: FTYPE_ASCII = "ASCII"
    character(len=*), parameter :: FTYPE_BINARY = "BINARY"
    character(len=*), parameter :: FTYPE_NETCDF = "NETCDF"
    character(len=*), parameter :: FTYPE_HDF5 = "HDF5"
    
    type (grid_t), pointer :: EsMgr_grid => null()
    character(len=25) :: EsMgr_ftype
    character(len=256) :: EsMgr_prefix
    logical :: EsMgr_save_in_file

    public :: FTYPE_ASCII, FTYPE_BINARY, FTYPE_NETCDF, FTYPE_HDF5
    public :: EsMgr_init
    public :: EsMgr_create_eAll, EsMgr_create_e 
    public :: EsMgr_get
    public :: EsMgr_save
    public :: EsMgr_save_in_file

contains

    ! EsMgr_init - Initalizes the manager sets the default filetype
    ! and determines if the Esoln files should be saved or not
    subroutine EsMgr_init(grid, context, save_in_file, prefix, ftype)

        implicit none

        type (grid_t), target, intent(in) :: grid
        type (ModEM_mpi_context_t), target, intent(in) :: context
        logical, intent(in), optional :: save_in_file
        character(len=*), intent(in), optional :: prefix
        character(len=*), intent(in), optional :: ftype

        logical :: save_in_file_lcl
        character(len=256) :: prefix_lcl
        character(len=256) :: ftype_lcl

        EsMgr_ctx => context

        if (present(save_in_file)) then
            save_in_file_lcl = save_in_file
        else
            save_in_file_lcl = .false.
        end if

        if (.not. save_in_file_lcl .and. present(prefix)) then
            write(0,*) "Warning: Argument 'prefix' was passed, but 'save_in_file' was not present"
            write(0,*) "Warning: 'prefix' will not have an effect. Set 'save_in_file' to true to save"
            write(0,*) "Warning: esolns in files"
        end if

        if (.not. save_in_file_lcl .and. present(ftype)) then
            write(0,*) "Warning: Argument 'ftype' was passed, but 'save_in_file' was not present"
            write(0,*) "Warning: 'ftype' will not have an effect. Set 'save_in_file' to true to save"
            write(0,*) "Warning: esolns in files"
        end if

        if (present(prefix)) then
            prefix_lcl = prefix
        else
            prefix_lcl = "esoln"
        end if

        if (present(ftype)) then
            ftype_lcl = ftype
        else
            ftype_lcl = FTYPE_ASCII
        end if

        select case (ftype_lcl)
            case (FTYPE_ASCII : FTYPE_BINARY)
            case (FTYPE_NETCDF)
                write(0,*) "ERROR: NetCDF not implemented for storing esoln files"
                write(0,*) "ERROR: Valid options are: [", trim(FTYPE_ASCII), " | ", trim(FTYPE_BINARY), "]"
                call ModEM_abort()
            case (FTYPE_HDF5)
                write(0,*) "ERROR: HDF5 not implemented for storing esoln files"
                write(0,*) "ERROR: Valid options are: [", trim(FTYPE_ASCII), " | ", trim(FTYPE_BINARY), "]"
                call ModEM_abort()
            case DEFAULT
                write(0,*) "ERROR: ", trim(ftype_lcl), " is not a valid file type for Esoln"
                write(0,*) "ERROR: Valid options are: [", trim(FTYPE_ASCII), " | ", trim(FTYPE_BINARY), "]"
                call ModEM_abort()
        end select

        EsMgr_grid => grid
        EsMgr_save_in_file = save_in_file_lcl
        EsMgr_prefix = prefix_lcl
        EsMgr_ftype = ftype_lcl

        call MPI_Barrier(MPI_COMM_WORLD, ierr)

    end subroutine EsMgr_init

    subroutine EsMgr_create_e(e, iTx, nPol, place_holder)

        implicit none

        type (solnVector_t), intent(inout) :: e
        integer, intent(in) :: iTx
        integer, intent(in), optional :: nPol
        logical, intent(in), optional :: place_holder

        logical :: place_holder_lcl

        if (present(place_holder)) then
            place_holder_lcl = place_holder
        else
            place_holder_lcl = EsMgr_save_in_file
        end if

        call create_solnVector(EsMgr_grid, iTx, e, place_holder=EsMgr_save_in_file)

    end subroutine EsMgr_create_e

    subroutine EsMgr_create_eAll(eAll, nTx, place_holder)

        implicit none

        type (solnVectorMTX_t), intent(inout) :: eAll
        integer, intent(in) :: nTx
        logical, intent(in), optional :: place_holder

        integer :: iTx
        type(solnVector_t) :: e
        logical :: place_holder_lcl

        if (present(place_holder)) then
            place_holder_lcl = place_holder
        else
            place_holder_lcl = EsMgr_save_in_file
        end if

        call create_solnVectorMTX(nTx, eAll)

        do iTx = 1, nTx
            call create_solnVector(EsMgr_grid, iTx, e, place_holder=place_holder_lcl)
            call copy_solnVector(eAll % solns(iTx), e)
        end do

    end subroutine EsMgr_create_eAll

    subroutine EsMgr_get(e, iTx, pol_index, from, prefix)

        implicit none

        type (solnVector_t), intent(inout) :: e
        integer, intent(in) :: iTx
        integer, intent(in), optional :: pol_index
        integer, intent(in), optional :: from
        character(len=*), intent(in), optional :: prefix

        ! If we are reading and writing files, do nothing
        if (EsMgr_save_in_file .and. EsMgr_ctx % rank_world == 0) then
            return
        end if

        if (EsMgr_save_in_file .and. EsMgr_ctx % rank_world /= 0) then
            call read_esoln_from_file(e, iTx, pol_index, prefix=prefix)
            return
        end if

        call EsMgr_recv_e(e, from)

    end subroutine EsMgr_get

    subroutine EsMgr_save(e, to, prefix)

        implicit none

        type (solnVector_t), intent(inout) :: e
        integer, intent(in), optional :: to
        character(len=*), intent(in), optional :: prefix

        if (EsMgr_save_in_file .and. EsMgr_ctx % rank_world == 0) then
            return
        end if


        if (EsMgr_save_in_file .and. .not. EsMgr_ctx % rank_world == 0) then
            call write_soln_to_file(e, prefix)
            return
        end if

        call EsMgr_send_e(e, to)

    end subroutine EsMgr_save

    subroutine write_soln_to_file(e, prefix, iPol)

        implicit none

        type (solnVector_t), intent(in) :: e
        character(len=*), intent(in), optional :: prefix
        integer, intent(in), optional :: iPol

        integer :: iPol_lcl
        character(len=256) :: prefix_lcl

        if (present(iPol)) then
            iPol_lcl = iPol
        else
            iPol_lcl = 1
        end if

        if (present(prefix)) then
            prefix_lcl = prefix 
        else
            prefix_lcl = ""
        end if

        call write_solnVector(e, trim(EsMgr_prefix)//trim(prefix_lcl), ftype=EsMgr_ftype, pol_index=iPol_lcl)

    end subroutine write_soln_to_file

    subroutine read_esoln_from_file(e, iTx, iPol, prefix)

        implicit none

        type (solnVector_t), intent(inout) :: e
        integer, intent(in) :: iTx
        integer, intent(in), optional :: iPol
        character(len=*), intent(in), optional :: prefix
        
        integer :: iPol_lcl
        character(len=256) :: prefix_lcl

        if (present(iPol)) then
            iPol_lcl = iPol
        else
            iPol_lcl = 1
        end if

        if (present(prefix)) then
            prefix_lcl = prefix 
        else
            prefix_lcl = ""
        end if

        call read_solnVector(e, trim(EsMgr_prefix)//trim(prefix_lcl), ftype=EsMgr_ftype, pol_index=iPol_lcl)

    end subroutine read_esoln_from_file

    subroutine EsMgr_send_e(e, to)

        implicit none

        type (solnVector_t), intent(inout) :: e
        integer, intent(in) :: to
        integer :: tag
        character, pointer, dimension(:) :: e_packed => null()

        if (to == 0) then
            tag = FROM_WORKER
        else
            tag = FROM_MASTER
        end if

        call int_create_e_param_place_holder(e, e_packed) 
        call int_Pack_packed_vec(e, e_packed)
        call MPI_Send(e_packed, Nbytes, MPI_PACKED, to, tag, EsMgr_ctx % comm_current, ierr) 

    end subroutine EsMgr_send_e

    subroutine EsMgr_recv_e(e, from)

        implicit none

        type (solnVector_t), intent(inout) :: e
        integer, intent(in) :: from
        integer :: tag
        character, pointer, dimension(:) :: e_packed => null()

        if (from == 0) then
            tag = FROM_MASTER
        else
            tag = FROM_WORKER
        end if

        call int_create_e_param_place_holder(e, e_packed)
        call MPI_Recv(e_packed, Nbytes, MPI_PACKED, from, tag, EsMgr_ctx % comm_current, STATUS, ierr)
        call int_Unpack_packed_vec(e, e_packed)

    end subroutine EsMgr_recv_e



    subroutine int_create_e_param_place_holder(e, holder)

         implicit none
         type(solnVector_t), intent(in)	:: e
         character, pointer, dimension(:), intent(inout) :: holder 
         integer                        :: Ex_size,Ey_size,Ez_size,Nbytes1,Nbytes2,Nbytes3,Nbytes4

         Ex_size=size(e%pol(1)%x)
         CALL MPI_PACK_SIZE(Ex_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes1,  ierr)
         Ey_size=size(e%pol(1)%y)
         CALL MPI_PACK_SIZE(Ey_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes2,  ierr)
         Ez_size=size(e%pol(1)%z)
         CALL MPI_PACK_SIZE(Ez_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes3,  ierr)
         CALL MPI_PACK_SIZE(1, MPI_INTEGER, MPI_COMM_WORLD, Nbytes4,  ierr)
         Nbytes=((Nbytes1+Nbytes2+Nbytes3+Nbytes4))+1

         if(associated(holder)) then
            deallocate(holder)
         end if
         allocate(holder(Nbytes))       

     end subroutine int_create_e_param_place_holder

     subroutine int_pack_packed_vec(e, unpacked_vec)

        implicit none

        type(solnVector_t), intent(in)	:: e
        character, pointer, dimension(:), intent(inout) :: unpacked_vec
        integer index,Ex_size,Ey_size,Ez_size

        Ex_size=size(e%pol(1)%x)
        Ey_size=size(e%pol(1)%y)
        Ez_size=size(e%pol(1)%z)
        index=1

        call MPI_Pack(e%pol(which_pol)%x(1,1,1),Ex_size, MPI_DOUBLE_COMPLEX, unpacked_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(e%pol(which_pol)%y(1,1,1),Ey_size, MPI_DOUBLE_COMPLEX, unpacked_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(e%pol(which_pol)%z(1,1,1),Ez_size, MPI_DOUBLE_COMPLEX, unpacked_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(e%tx, 1, MPI_INTEGER, unpacked_vec, Nbytes, index, MPI_COMM_WORLD, ierr)

    end subroutine int_pack_packed_vec


    subroutine int_unpack_packed_vec(e, packed_vec)

        implicit none

        type(solnVector_t), intent(inout)	:: e
        character, pointer, dimension(:), intent(inout) :: packed_vec

        integer index,Ex_size,Ey_size,Ez_size

        Ex_size=size(e%pol(1)%x)
        Ey_size=size(e%pol(1)%y)
        Ez_size=size(e%pol(1)%z)
        index=1

        call MPI_Unpack(packed_vec, Nbytes, index, e%pol(which_pol)%x(1,1,1),Ex_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(packed_vec, Nbytes, index, e%pol(which_pol)%y(1,1,1),Ey_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(packed_vec, Nbytes, index, e%pol(which_pol)%z(1,1,1),Ez_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(packed_vec, Nbytes, index, e%tx,1, MPI_INTEGER,MPI_COMM_WORLD, ierr)

    end subroutine int_unpack_packed_vec
    

end module ESolnManager
