! EsolnManager Module
! ====================
!
! The EsolnManager module is a module that helps with transferring electric 
! and adjoint fields to and from worker tasks. It does so in two ways:
! 
! 1. MPI Transferring - Original ModEM MPI Behavior
!   
!   If specified, main tasks and worker tasks will communicate their 
!   electric fields/adjoint solutions to each other using MPI send and
!   receives.
!
!   This will cause the main task to allocate space for eAll.
!
! 2. Saving Electric fields/Adjoin Solutions to File
!
!   If specified, the worker tasks will read and write out electric fields
!   when they are told to access them by the main task rather than communicating
!   the electric fields back to the main task.
!
!   
! How to use this Module:
!
! To use this module, you will want to first call EsMgr_init with either:
!
! save_in_file = .true. or .false. depending on how you wish electric fields to
! be saved. 
!
! Regardless of value of save_in_file, you will need to call EsMgr_get and EsMgr_save
! in order to complete a 'transaction'.
!
! With save_in_file = .false., EsMgr_get will call a MPI_Recv and EsMgr_save will call
! a corresponding MPI_Send
!
! With save_in_file = .ture., the worker tasks will write out their e fields to disk via 
! EsMgr_save and will read e fields using EsMgr_get. The main task will still need to call
! EsMgr_get/EsMgr_save as the worker task will communicate to the main task when it is done 
! reading/writing.

module ESolnManager

    use mpi
    use Declaration_MPI
    use GridDef
    use SolnSpace
    use utilities

    implicit None

    private

    type (ModEM_mpi_context_t), pointer :: EsMgr_ctx

    character(len=*), parameter :: FTYPE_ASCII = "ascii"
    character(len=*), parameter :: FTYPE_BINARY = "binary"

    character(len=*), parameter :: E_FIELD_TYPE_FWD = "FWD"
    character(len=*), parameter :: E_FIELD_TYPE_JMULTT = "JmultT"

    type (grid_t), pointer :: EsMgr_grid => null()
    character(len=25) :: EsMgr_ftype
    character(len=256) :: EsMgr_prefix
    logical :: EsMgr_save_in_file

    character, pointer, dimension(:) :: esmgr_holder => null()
    integer :: esmgr_holder_bytes
    logical :: esmgr_holder_allocated

    public :: FTYPE_ASCII, FTYPE_BINARY
    public :: E_FIELD_TYPE_FWD, E_FIELD_TYPE_JMULTT
    public :: EsMgr_init
    public :: EsMgr_create_solnVectorMTX, EsMgr_create_e
    public :: EsMgr_get
    public :: EsMgr_save
    public :: EsMgr_save_in_file

contains

    ! EsMgr_init - Initializes the ESolnManager sets several module variables
    ! and determines if the Esoln files should be saved or not. The main task and all
    ! and worker tasks must call this before calling any other subroutines in this 
    ! module.
    !
    ! Arguments:
    !  * grid - grid_t target
    !       The grid
    !  * context - ModEM_mpi_context_t target
    !       An initalized ModEM_mpi_context_t which the ESolnManager should use
    !       when doing MPI communication routines.
    !  * save_in_file - logical - optional - Default: .false.
    !       If save_in_file is present and .true., the electric fields will be saved in a file
    !       rather than passed from the main task and worker tasks via MPI.
    !
    !       If .false. or not present, the electric fields will be passed from the main and
    !       worker task via MPI.
    !  * prefix - character - optional - Default: 'esoln'
    !       A character string prefix which to label the electric field/adjoint solution 
    !       filenames if they are written out to disk.
    !       
    !       If save_in_file is .false. or not present, prefix has no effect.
    !  * ftype - character - optional - Default: 'FTYPE_BINARY'
    !       Determines the type of file the cvectors should be written as if `save_in_file`
    !       is true.
    !
    !       Options are above in this module and are: 
    !           * FTYPE_ASCII - ASCII representation
    !           * FTYPE_BINARY - Fortran unformatted binary
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

        if ( .not. (save_in_file_lcl .and. present(prefix)) ) then
            if (EsMgr_ctx % rank_world == 0) then
                write(0,*) "Warning: Argument 'prefix' was passed, but 'save_in_file' was not present"
                write(0,*) "Warning: 'prefix' will not have an effect. Set 'save_in_file' to true to save"
                write(0,*) "Warning: esolns in files"
            end if
        end if

        if ( .not. (save_in_file_lcl .and. present(ftype)) ) then
            if (EsMgr_ctx % rank_world == 0) then
                write(0,*) "Warning: Argument 'ftype' was passed, but 'save_in_file' was not present"
                write(0,*) "Warning: 'ftype' will not have an effect. Set 'save_in_file' to true when calling"
                write(0,*) "WARNING: EsMgr_init to save esolns in files"
            end if
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

    ! EsMgr_create_eAll
    !
    !  Convenience function to initialize an solnVectorMTX_t with solnVectors. 
    !
    !  If EsMgr_init was called with save_in_file=.true., then the solnVectors_t
    !  created in this subroutine will be created as 'place_holders' and their 
    !  pol cvectors variables will also be created as 'place_holders'. These 
    !  cvectors will not have their x, y, z arrays will not be allocated.
    !
    !  However, if save_in_file=.false. then these arrays will be allocated.
    !
    ! Note: EsMgr_init should be called before calling this function.
    !   
    ! Arguments:
    !   eAll - solnVectorMTX_T
    !       The solnVectorMTX_t to hold eAll. See description above.
    !   nTx - integer
    !       The number of transmitters to create for the SolnVectorMTX_t
    !   grid - grid_t target
    !       The grid that should be associated with this SolnVectorMTX_t
    !   place_holder - logical  - Defaults to save_in_file passed in during init
    !       Overrides the default action described above.
    !       
    subroutine EsMgr_create_solnVectorMTX(eAll, nTx, grid, place_holder)

        implicit none

        type (solnVectorMTX_t), intent(inout) :: eAll
        integer, intent(in) :: nTx
        type(grid_t), intent(in), target :: grid
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
            call create_solnVector(grid, iTx, eall % solns(iTx), place_holder=place_holder_lcl)
        end do

    end subroutine EsMgr_create_solnVectorMTX

    ! EsMgr_create_e
    !
    ! Convenience function to create a single solnVector_t. 
    !
    !  If EsMgr_init was called with save_in_file=.true., then the solnVector_t
    !  created in this subroutine by the main task will be created as 'place_holders' 
    !  and their pol cvectors variables will also be created as 'place_holders'. 
    !  These cvectors will not have their x, y, z arrays will not be allocated.
    !
    !  If save_in_file=.true. and this routine is called by a worker task, then
    !  the this will routine will completly allocate the cvectors associated with
    !  the solnVector_t.
    !
    !  However, if save_in_file=.false. then these arrays will be allocated.
    !
    !  Note: EsMgr_init should be called before calling this function.
    !
    !  Arguments:
    !   e - solnVector_t 
    !       The solnVector_t that will be created
    !   iTx - integer
    !       The transmitter number to associate with this solnVector_t
    !   place_holder - logical - Defaults to save_in_file passed in during init
    !       Overides the default action described above.
    !
    subroutine EsMgr_create_e(e, iTx, place_holder)

        implicit none

        type (solnVector_t), intent(inout) :: e
        integer, intent(in) :: iTx
        logical, intent(in), optional :: place_holder

        logical :: place_holder_lcl

        if (present(place_holder)) then
            place_holder_lcl = place_holder
        else
            if (EsMgr_ctx % rank_current == 0) then
                place_holder_lcl = EsMgr_save_in_file
            else
                ! This is a worker task - Don't make this a place holder
                place_holder_lcl = .false.
            end if
        end if

        call create_solnVector(EsMgr_grid, iTx, e, place_holder=place_holder_lcl)

    end subroutine EsMgr_create_e

    ! create_prefix - Generate a prefix for a solnVector file name
    !
    ! This routine will create a prefix that can be used to name the
    ! solnVector file.
    !
    ! The prefix will always start with EsMgr_prefix, which was
    ! passed in during the init call, the rest will be in one of
    ! depending on which optional arguments are present or not:
    !
    ! - $EsMgr_prefix.$E_field_type.$label
    ! - $EsMgr_prefix.$E_field_type
    ! - $EsMgr_prefix.$label
    ! - $EsMgr_prefix
    !
    subroutine create_prefix(prefix, E_field_type, label)

        implicit none

        character(len=*), intent(out) :: prefix
        character(len=*), intent(in), optional :: E_field_type
        character(len=*), intent(in), optional :: label

        character(len=512) :: PREFIX_FNAME_FORMAT
        logical :: use_label

        use_label = .false.

        if (present(label)) then
            use_label = trim(label) /= "NULL" .and. len_trim(label) /= 0
        end if

        write(0,*) 'LABEL: ', use_label, trim(label)

        if (present(E_field_type) .and. use_label) then
            PREFIX_FNAME_FORMAT = '(A, A, A, A, A)'
            write(prefix, PREFIX_FNAME_FORMAT) trim(EsMgr_prefix), ".", trim(E_field_type), ".", trim(label)
        else if (present(E_field_type) .and. .not. use_label) then
            PREFIX_FNAME_FORMAT = '(A, A, A)'
            write(prefix, PREFIX_FNAME_FORMAT) trim(EsMgr_prefix), ".", trim(E_field_type)
        else if (.not. present(E_field_type) .and. use_label) then
            PREFIX_FNAME_FORMAT = '(A, A, A)'
            write(prefix, PREFIX_FNAME_FORMAT) trim(EsMgr_prefix), ".", trim(label)
        else
            PREFIX_FNAME_FORMAT = '(A, A)'
            write(prefix, PREFIX_FNAME_FORMAT) trim(EsMgr_prefix)
        end if

    end subroutine create_prefix

    ! EsMgr_get
    !
    !  Populate e with the transmitter iTx and polarization pol_index
    !  via MPI or by reading a file.
    !
    !  EsMgr_save_in_file == .false.
    !  ======================
    !
    !  If `EsMgr_save_in_file == .false.`, then a task will receive the
    !  solnVector via MPI from the MPI tasks specified in the `from`,
    !  argument.
    !
    !  Note: The task in the `from` argument must call EsMgr_save
    !  with a corosponding task number in the `to` argument.
    !
    !  EsMgr_save_in_file == .true.
    !  ======================
    !
    !  If the EsolnManager was initialized with `EsMgr_save_in_file == .true.`,
    !  then worker_tasks that call this function will read the corresponding
    !  file that matches the requested arguments.
    !
    !  After reading, the worker task will communicate with the main task
    !  to indicate that it is done reading.
    !
    !  Main Tasks:
    !
    !  If EsMgr_save_in_file is .true., then the main task will wait for
    !  the task speciried in from to finish writing it's solnVector_t
    !  that it made in ESoln_save.
    !
    ! Note: EsMgr_init should be called before calling this function.
    !
    ! Arguments:
    !  e - solnVector_t
    !    The solnVector_t to read the polarization pol_index of iTx into
    !  iTx - integer
    !    The transmitter number to read
    !  pol_index - integer - optional
    !    The polarization to read
    !  from - integer - optional
    !    The MPI task to receive the solnVector_t from if EsMgr_save_in_file is .false.
    !  E_field_type - character - optional
    !    The type of the electric field that the SolnVector was saved with,
    !    normally this is either E_FIELD_TYPE_FWD or E_FIELD_TYPE_JMULTT.
    !  label - character - optional
    !    The additonal label that was used to save this SolnVector. Normally
    !    this is used to label the specific line search calulation (e.g.
    !    trial, fwd01, fwd02 etc.).
    subroutine EsMgr_get(e, iTx, pol_index, from, E_field_type, label)

        implicit none

        type (solnVector_t), intent(inout) :: e
        integer, intent(in) :: iTx
        integer, intent(in), optional :: pol_index
        integer, intent(in), optional :: from
        character(len=*), intent(in), optional :: E_field_type
        character(len=*), intent(in), optional :: label
        character(len=512) :: prefix

        e % tx = iTx

        ! If we are reading and writing files, do nothing
        if (EsMgr_save_in_file .and. EsMgr_ctx % rank_world == 0) then
            call wait_on_task(from)
            return
        end if

        if (EsMgr_save_in_file .and. EsMgr_ctx % rank_world /= 0) then
            call create_prefix(prefix, E_field_type, label)
            call read_esoln_from_file(e, iTx, pol_index, prefix=prefix)
            return
        end if

        call EsMgr_recv_e(e, from)

    end subroutine EsMgr_get

    ! EsMgr_save
    !
    ! Save `e` by either writing it to disk or sending it to the main task via MPI.
    !
    ! EsMgr_save_in_file = .false.
    ! =============================
    !
    ! If `EsMgr_save_in_file == .false.`, then a calling task will send `e` to the task
    ! specified in the `to` argument.
    !
    ! To receive this `e`, a task will need to call EsMgr_get with the correspond from.
    !
    ! EsMgr_save_in_file = .true.
    ! =============================
    !
    ! If EsMgr_save_in_file is .true., then when a worker task that calls this 
    ! subroutine will write out e to a corresponding file. The worker task will then
    ! communicate to the main task to declare that it is done writing.
    !
    ! If the main task calls this routine with EsMgr_save_in_file == .true., then
    ! it will return immediately as no work is needed.
    !
    ! Arguments
    !  e - solnVector_t
    !    The solnVector_t to save
    !  to - integer - optional 
    !    The task to send e to
    !  E_field_type - character - optional
    !    The type of the electric field, normally this is either E_FIELD_TYPE_FWD
    !    or E_FIELD_TYPE_JMULTT, but could be anything. However, the same E_field_type
    !    will need to be used in EsMgr_get to read the same solnVector.
    !  label - character - optional
    !    An additional label to use to label the solnVector file (if
    !    EsMgr_save_in_file == .true.). Normally, this is used to label the specific
    !    line search calculations (e.g. trial, fwd01, fwd02, etc.). Again, the same
    !    label will need to be specified when calling the coorosponding EsMgr_get.
    !
    subroutine EsMgr_save(e, to, E_field_type, label)

        implicit none

        type (solnVector_t), intent(inout) :: e
        integer, intent(in), optional :: to
        character(len=*), intent(in), optional :: E_field_type
        character(len=*), intent(in), optional :: label

        character(len=512) :: prefix

        if (EsMgr_save_in_file .and. EsMgr_ctx % rank_world == 0) then
            return
        end if

        if (EsMgr_save_in_file .and. .not. EsMgr_ctx % rank_world == 0) then
            call create_prefix(prefix, E_field_type, label)
            call EsMgr_write_to_file(e, prefix)
            call communicate_file_done_writing()
            return
        end if

        call EsMgr_send_e(e, to)

    end subroutine EsMgr_save

    subroutine communicate_file_done_writing()

        implicit none

        integer :: dummy_buffer
        integer :: tag
        integer :: count
        integer :: to

        count = 0
        tag = 0
        to = 0

        call MPI_Send(dummy_buffer, count, MPI_INTEGER, to, tag, EsMgr_ctx % comm_current, ierr)

    end subroutine communicate_file_done_writing

    subroutine wait_on_task(task)

        implicit none

        integer, intent(in) :: task

        integer :: dummy_buff
        integer :: count
        integer :: tag

        count = 0
        tag = 0

        call MPI_Recv(dummy_buff, count, MPI_INTEGER, task, tag, EsMgr_ctx % comm_current, MPI_STATUS_IGNORE, ierr)

    end subroutine wait_on_task

    subroutine EsMgr_write_to_file(e, prefix, iPol, ftype)

        implicit none

        type (solnVector_t), intent(in) :: e
        character(len=*), intent(in), optional :: prefix
        integer, intent(in), optional :: iPol
        character(len=*), optional, intent(in) :: ftype

        integer :: iPol_lcl
        character(len=256) :: prefix_lcl
        character(len=256) :: ftype_lcl

        if (present(iPol)) then
            iPol_lcl = iPol
        else
            iPol_lcl = 1
        end if

        if (present(prefix)) then
            prefix_lcl = prefix 
        else
            prefix_lcl = "esoln"
        end if

        if (present(ftype)) then
            ftype_lcl = ftype
        else
            ftype_lcl = EsMgr_ftype
        endif

        call write_solnVector(e, trim(prefix_lcl), ftype=ftype_lcl, pol_index=iPol_lcl)

    end subroutine EsMgr_write_to_file

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
            prefix_lcl = "esoln"
        end if

        e % tx = iTx
        call read_solnVector(e, trim(prefix_lcl), ftype=EsMgr_ftype, pol_index=iPol_lcl)

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

        call int_create_e_param_place_holder(e) 
        call int_Pack_packed_vec(e)
        call MPI_Send(esmgr_holder, esmgr_holder_bytes, MPI_PACKED, to, tag, EsMgr_ctx % comm_current, ierr) 

    end subroutine EsMgr_send_e

    subroutine EsMgr_recv_e(e, from)

        implicit none

        type (solnVector_t), intent(inout) :: e
        integer, intent(in) :: from
        integer :: tag

        if (from == 0) then
            tag = FROM_MASTER
        else
            tag = FROM_WORKER
        end if

        call int_create_e_param_place_holder(e)
        call MPI_Recv(esmgr_holder, esmgr_holder_bytes, MPI_PACKED, from, tag, EsMgr_ctx % comm_current, STATUS, ierr)
        call int_Unpack_packed_vec(e)

    end subroutine EsMgr_recv_e

    subroutine int_create_e_param_place_holder(e)

         implicit none
         type(solnVector_t), intent(in)	:: e
         integer                        :: Ex_size,Ey_size,Ez_size,Nbytes1,Nbytes2,Nbytes3,Nbytes4

         if (esmgr_holder_allocated) then
             return
         end if

         Ex_size=size(e%pol(1)%x)
         CALL MPI_PACK_SIZE(Ex_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes1,  ierr)
         Ey_size=size(e%pol(1)%y)
         CALL MPI_PACK_SIZE(Ey_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes2,  ierr)
         Ez_size=size(e%pol(1)%z)
         CALL MPI_PACK_SIZE(Ez_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes3,  ierr)
         CALL MPI_PACK_SIZE(1, MPI_INTEGER, MPI_COMM_WORLD, Nbytes4,  ierr)
         esmgr_holder_bytes=((Nbytes1+Nbytes2+Nbytes3+Nbytes4))+1

         allocate(esmgr_holder(esmgr_holder_bytes))       
         esmgr_holder_allocated = .true.

     end subroutine int_create_e_param_place_holder

     subroutine int_pack_packed_vec(e)

        implicit none

        type(solnVector_t), intent(in)	:: e
        integer index,Ex_size,Ey_size,Ez_size

        Ex_size=size(e%pol(1)%x)
        Ey_size=size(e%pol(1)%y)
        Ez_size=size(e%pol(1)%z)
        index=1

        call MPI_Pack(e%pol(which_pol)%x(1,1,1),Ex_size, MPI_DOUBLE_COMPLEX, esmgr_holder, esmgr_holder_bytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(e%pol(which_pol)%y(1,1,1),Ey_size, MPI_DOUBLE_COMPLEX, esmgr_holder, esmgr_holder_bytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(e%pol(which_pol)%z(1,1,1),Ez_size, MPI_DOUBLE_COMPLEX, esmgr_holder, esmgr_holder_bytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(e%tx, 1, MPI_INTEGER, esmgr_holder, esmgr_holder_bytes, index, MPI_COMM_WORLD, ierr)

    end subroutine int_pack_packed_vec


    subroutine int_unpack_packed_vec(e)

        implicit none

        type(solnVector_t), intent(inout)	:: e

        integer index,Ex_size,Ey_size,Ez_size

        Ex_size=size(e%pol(1)%x)
        Ey_size=size(e%pol(1)%y)
        Ez_size=size(e%pol(1)%z)
        index=1

        call MPI_Unpack(esmgr_holder, esmgr_holder_bytes, index, e%pol(which_pol)%x(1,1,1),Ex_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(esmgr_holder, esmgr_holder_bytes, index, e%pol(which_pol)%y(1,1,1),Ey_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(esmgr_holder, esmgr_holder_bytes, index, e%pol(which_pol)%z(1,1,1),Ez_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(esmgr_holder, esmgr_holder_bytes, index, e%tx,1, MPI_INTEGER,MPI_COMM_WORLD, ierr)
    end subroutine int_unpack_packed_vec
    

end module ESolnManager
