
Module Declaration_MPI
     use file_units
     use math_constants
#ifdef MPI
     use mpi
#if defined(CUDA)
     use cudaFortMap
#elif defined(HIP)
     use hipFortMap
#endif
     implicit none

! Declaration of general MPI stuff
!********************************************************************
Integer        :: taskid,total_number_of_Proc,number_of_workers
Integer        :: MASTER, FROM_MASTER, FROM_WORKER,TAG,ierr,dest
Integer        :: STATUS(MPI_STATUS_SIZE)
parameter         (MASTER=0,FROM_MASTER=1,FROM_WORKER=2,Tag=1)
!********************************************************************
! additional parameters needed by two-layered parallelization
!********************************************************************
integer        :: comm_world, comm_leader, comm_local
integer        :: rank_world, rank_leader, rank_local, rank_node
integer        :: size_world, size_leader, size_local, size_node
integer        :: info_world, info_leader, info_local
integer        :: group_world, group_leader
! this is used to store the name of proc/cpu for current rank and identify 
! ranks from different hosts, useful when grouping cpus according to 
! topology - note hostname_len cannot exceed 40 (hard coded here)
character*(40) :: hostname_MPI
integer        :: hostname_len, ngroup
! one master switch for parallel paradigm - currently not supported for 
! on-the-fly modification, 
! TODO: need to think about it - if we really need an on-the-fly config
! 0 : simple grouping, 1 proc per each fwd/trn task (default)
! 1 : topology grouping, 1 node per each fwd/trn task, for debug purpose
! 2 : equal grouping, n procs per each fwd/trn task
! 3 : dynamic grouping, variable number of procs per each fwd/trn task, 
!     for load-balancing
#if defined(FG) && (defined(CUDA) || defined(HIP))
integer        :: para_method = 2 ! use 1 for debug only
#elif defined(FG)
integer        :: para_method = 2
#elif defined(PETSC) 
integer        :: para_method = 2
#else
integer        :: para_method = 0
#endif
!********************************************************************
! additional parameters needed by CUDA acceleration
!********************************************************************
   integer,target         :: size_gpu = 0 
   integer,target         :: size_gpu_total = 0 
#if defined(CUDA) || defined(HIP)
   type(c_ptr)            :: size_gpuPtr
   type(c_ptr)            :: size_gpu_totalPtr
#endif
! change the cpus_per_gpu if you want to use more than one cpus to 
! feed one gpu, modify at your own risk! 
#if defined(FG) && (defined(CUDA) || defined(HIP))
   type(ncclComm)         :: comm_nccl        ! nccl/rccl communicator
   integer                :: rank_nccl        ! nccl/rccl rank
   integer                :: size_nccl        ! nccl/rccl size
   integer                :: ncclIsInit=0     ! flag
   type(ncclUniqueId)     :: uid              ! nccl id
   integer                :: cpus_per_gpu = 3 ! override for GPU+FG
#else
   integer                :: cpus_per_gpu = 3 ! hard coded here 
#endif
integer                   :: device_id = -1
! this is used to store the timer of each mpi sub-process
DOUBLE PRECISION          :: previous_time
integer, pointer, dimension(:) :: group_sizes

!********************************************************************
! Parameters required to create an MPI derived data types.
!********************************************************************
Integer        :: cvector_mpi_3D,gridDef3D_mpi,eAll_mpi,dvecMTX_mpi
Integer        :: dvec_mpi,grid_t_mpi,worker_job_task_mpi
Integer        :: modelParam_t_mpi_sing,userdef_control_MPI
Integer        :: modelParam_t_mpi,extent
integer        :: oldtypes(0:20), blockcounts(0:20),offsets(0:20)
integer        :: block_lengths(0:20)
integer        :: displacements(0:20)
integer        :: address(0:21)
integer        :: typelist(0:21)
!********************************************************************



! Parameters used in communication
!********************************************************************
Integer        :: answers_to_receive,received_answers,recv_loop
Integer        :: who, which_stn,which_per,which_dt,which_pol,orginal_nPol
Integer , pointer, dimension(:)  :: eAll_location
logical                          :: eAll_exist=.false.
real*8,   pointer, dimension(:)  :: model_para_vec
character, pointer, dimension(:) :: eAll_para_vec       !! needed for MPI_pack/MPI_unpack; counted in bytes
character, pointer, dimension(:) :: e_para_vec          !! needed for MPI_pack/MPI_unpack; counted in bytes
character, pointer, dimension(:) :: sigma_para_vec      !! needed for MPI_pack/MPI_unpack; counted in bytes
character, pointer, dimension(:) :: data_para_vec       !! needed for MPI_pack/MPI_unpack; counted in bytes
character, pointer, dimension(:) :: worker_job_package  !! needed for MPI_pack/MPI_unpack; counted in bytes
character, pointer, dimension(:) :: userdef_control_package !! needed for MPI_pack/MPI_unpack; counted in bytes

Integer                          :: Nbytes              !! used in all MPI_pack/MPI_unpack
!********************************************************************


Integer                          :: nPol_MPI

! Time measuring
!********************************************************************
DOUBLE PRECISION    :: starttime,endtime,time_used
DOUBLE PRECISION    :: starttime_total,endtime_total
!********************************************************************


! A drived data type to distribute Info. between Processors.
!********************************************************************
type :: define_worker_job
     SEQUENCE
     character*80  :: what_to_do='NOTHING'
     Integer       :: per_index,Stn_index,pol_index,data_type_index,data_type
     Integer       :: taskid
     logical       :: keep_E_soln=.false.
     logical       :: several_Tx=.false.
     logical       :: create_your_own_e0=.false.
     ! 2022.10.06, Liu Zhongyin, add iSite storing the site index in rx of dataBlock_t
     Integer       :: iSite
 end type define_worker_job
type(define_worker_job), save :: worker_job_task
!********************************************************************

type :: ModEM_mpi_context_t
    integer :: comm_current
    integer :: comm_world
    integer :: comm_local
    integer :: rank_world
    integer :: rank_local
    integer :: rank_current
    integer :: size_local
    integer :: number_of_workers
end type ModEM_mpi_context_t

type (ModEM_mpi_context_t), target, save :: modem_ctx

Contains

!##########################################################################
subroutine create_worker_job_task_place_holder

     implicit none
     integer index,Nbytes1,Nbytes2,Nbytes3

       CALL MPI_PACK_SIZE(80, MPI_CHARACTER, MPI_COMM_WORLD, Nbytes1,  ierr)
       ! 2019.05.10, Liu Zhongyin, replace 6 with 7, for the iSite
       ! CALL MPI_PACK_SIZE(6, MPI_INTEGER, MPI_COMM_WORLD, Nbytes2,  ierr)
       CALL MPI_PACK_SIZE(7, MPI_INTEGER, MPI_COMM_WORLD, Nbytes2,  ierr)
       CALL MPI_PACK_SIZE(3, MPI_LOGICAL, MPI_COMM_WORLD, Nbytes3,  ierr)

         Nbytes=(Nbytes1+Nbytes2+Nbytes3)+1

         if(.not. associated(worker_job_package)) then
            allocate(worker_job_package(Nbytes))
         end if

end subroutine create_worker_job_task_place_holder
!*******************************************************************************

subroutine Pack_worker_job_task
implicit none
integer index

index=1

        call MPI_Pack(worker_job_task%what_to_do,80, MPI_CHARACTER, worker_job_package, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(worker_job_task%per_index ,1 , 		MPI_INTEGER  , worker_job_package, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(worker_job_task%Stn_index ,1 ,	 	MPI_INTEGER  , worker_job_package, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(worker_job_task%pol_index ,1 , 		MPI_INTEGER  , worker_job_package, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(worker_job_task%data_type_index ,1 , 	MPI_INTEGER  , worker_job_package, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(worker_job_task%data_type ,1 , 		MPI_INTEGER  , worker_job_package, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(worker_job_task%taskid ,1 , 			MPI_INTEGER  , worker_job_package, Nbytes, index, MPI_COMM_WORLD, ierr)

        call MPI_Pack(worker_job_task%keep_E_soln,1, 		MPI_LOGICAL, worker_job_package, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(worker_job_task%several_Tx,1, 		MPI_LOGICAL, worker_job_package, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(worker_job_task%create_your_own_e0,1, MPI_LOGICAL, worker_job_package, Nbytes, index, MPI_COMM_WORLD, ierr)

        ! 2019.05.08, Liu Zhongyin, add iSite for rx in dataBlock_t
        call MPI_Pack(worker_job_task%iSite, 1,             MPI_INTEGER, worker_job_package, Nbytes, index, MPI_COMM_WORLD, ierr)

end subroutine Pack_worker_job_task

subroutine Unpack_worker_job_task
implicit none
integer index
index=1
        call MPI_Unpack(worker_job_package, Nbytes, index, worker_job_task%what_to_do,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)

        call MPI_Unpack(worker_job_package, Nbytes, index, worker_job_task%per_index ,1 , MPI_INTEGER,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(worker_job_package, Nbytes, index, worker_job_task%Stn_index ,1 , MPI_INTEGER,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(worker_job_package, Nbytes, index, worker_job_task%pol_index ,1 , MPI_INTEGER,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(worker_job_package, Nbytes, index, worker_job_task%data_type_index ,1 , MPI_INTEGER,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(worker_job_package, Nbytes, index, worker_job_task%data_type ,1 , MPI_INTEGER,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(worker_job_package, Nbytes, index, worker_job_task%taskid ,1 , MPI_INTEGER,MPI_COMM_WORLD, ierr)
        
        call MPI_Unpack(worker_job_package, Nbytes, index, worker_job_task%keep_E_soln,1, MPI_LOGICAL,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(worker_job_package, Nbytes, index, worker_job_task%several_Tx,1, MPI_LOGICAL,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(worker_job_package, Nbytes, index, worker_job_task%create_your_own_e0,1, MPI_LOGICAL,MPI_COMM_WORLD, ierr)

        ! 2019.05.08, Liu Zhongyin, add iSite for rx in dataBlock_t
        call MPI_Unpack(worker_job_package, Nbytes, index, worker_job_task%iSite, 1, MPI_INTEGER,MPI_COMM_WORLD, ierr)

end subroutine Unpack_worker_job_task

subroutine gather_runtime(comm_current,time_passed,time_buff)
! simple subroutine to get the runtime of each sub tasks to access the
! parallel efficiency 
! collective on comm_current
      implicit none
      real(kind=prec), intent(in)                :: time_passed
      integer,intent(in)                         :: comm_current
      real(kind=prec), intent(out),pointer,dimension(:)   :: time_buff
      integer                                    :: current_rank
      integer                                    :: current_size,root=0
      call MPI_COMM_RANK(comm_current,current_rank,ierr)
      call MPI_COMM_SIZE(comm_current,current_size,ierr)
      allocate(time_buff(current_size))
      call MPI_Gather(time_passed, 1, MPI_DOUBLE_PRECISION, time_buff, 1,     &
     &     MPI_DOUBLE_PRECISION, root,comm_current,ierr) 
      return
end subroutine gather_runtime

subroutine get_host_topology(nTx, nPol, host_sizes)
    ! a silly subroutine to set get topology according to the hostnames
    ! the worker procs with the same hostname, are grouped togather. 
    ! e.g. if we have a set-up like this:
    ! machine1 0 1 2 3 
    ! machine2 4 5 6 7  
    ! machine3 8 9   
    ! here we consider the master should be in a seperate host, so
    ! the host topology size will be determined as 1 3 4 2 (1+3 hosts)
    ! collective on comm_world
    ! note we didn't use the MPI-3.0 features (MPI_GET/PUT)
    ! for compatiblity considerations
    ! NOTE: one should only use this with comm_world
     implicit none
     integer, intent(in)                        :: nTx,nPol
     integer, intent(out),pointer,dimension(:)  :: host_sizes
    ! local variables 
     integer                                 :: temp_sizes(nTx*nPol+1)
     integer                                 :: iProc, ierr, nHost, iHost
     integer                                 :: current_host, procs_in_host
     character *(40)                         :: crnt_hostname, prev_hostname

     ! of course we cannot have more hosts than the number of workers
     if (rank_world.eq.0) then ! the root process determines the grouping
         ! the first group always has a size of 1
         temp_sizes(1) = 1 
         ! start from the second group (should have least 1 proc)
         current_host = 2
         procs_in_host = 1
         ! receive the host name from all workers
         call MPI_RECV(prev_hostname, 40, MPI_CHARACTER, 1,  &
    &        FROM_WORKER, comm_world, STATUS, ierr)
         ! loop through all workers
         do iProc = 2, number_of_workers
             call MPI_RECV(crnt_hostname, 40, MPI_CHARACTER, iProc,  &
    &            FROM_WORKER, comm_world, STATUS, ierr)
             if (trim(crnt_hostname).eq. trim(prev_hostname)) then
                 procs_in_host = procs_in_host + 1
             else
                 prev_hostname = crnt_hostname
                 temp_sizes(current_host) = procs_in_host
                 current_host = current_host + 1
                 procs_in_host = 1
             end if
         end do
         ! last host
         temp_sizes(current_host) = procs_in_host
         ! now allocate the new 
         nHost = current_host
         call MPI_BCAST(nHost, 1, MPI_INTEGER, 0, comm_world, ierr)
         allocate(host_sizes(nHost))
         host_sizes = temp_sizes(1:nHost)
     else
         ! workers - only send the host name to MASTER
         call MPI_SEND(hostname_MPI, 40, MPI_CHARACTER, 0,     &
    &        FROM_WORKER, comm_world, ierr)
         call MPI_BCAST(nHost, 1, MPI_INTEGER, 0, comm_world, ierr)
         allocate(host_sizes(nHost))
     endif
     call MPI_BCAST(host_sizes,nHost,MPI_INTEGER, 0,comm_world,ierr)
     ! now try to determine the rank of current proc in
     ! current machine node
     ! rank_node is the rank of the process in the entire node
     ! size_node is the number of processes in that node
     size_node = 1
     do iHost = 1, nHost
         if (rank_world .lt. SUM(host_sizes(1:iHost))) then
             size_node = host_sizes(iHost)
             exit
         endif
     end do
     if (iHost == 1) then
         rank_node = rank_world
     else
         rank_node = rank_world - SUM(host_sizes(1:iHost-1))
     end if
     if (output_level .gt. 3) then
         write(6,*) 'rank_world = ', rank_world, ' size_node = ', size_node, &
            ' rank_node = ', rank_node, ' nHost = ', nHost
     endif
end subroutine get_host_topology

! stub, reserved for multi-gpu
! I am not yet sure how to convert the current, ad-hoc setup to a 
! more reasonale configuration
! subroutine get_gpu_affinity(host_sizes, gpu_sizes, comm_current)
!      ! a silly subroutine to get to know the number of gpus in 
!      ! each group/host - useful for multi-machine, multi-card setup
!      integer, intent(in),pointer,dimension(:)     :: host_sizes
!      integer, intent(inout),pointer,dimension(:)  :: gpu_sizes
!      integer, intent(in)                          :: comm_current
!      ! local variables
!      integer                                      :: size_gpu
!      integer                                      :: nHost
!      nHost = size(host_sizes)
!      size_gpuPtr = c_loc(size_gpu) ! kind of crude here
!      ierr = hipGetDeviceCount(size_gpuPtr)
!      if ((ctrl%output_level .gt. 3).and. (taskid .eq. 0)) then
!          write(6,*) 'number of GPU devices = ', size_gpu, 'err = ', ierr
!      endif
! end subroutine get_gpu_affinity

#endif

end module Declaration_MPI
