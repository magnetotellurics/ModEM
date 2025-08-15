
! the develop branch, which worked a testing bed for the PETSc version
! the main purpose is to utilize more cpus for the FWD/ADJ calculations
!
! - this is now supposed to replace the old Main_MPI version
! for both the MF and SP versions
! (the number was restricted to one per TX in the old version)
! this also includes a more flexible routine to enable a heteogeneous
! distribution of CPU power.
module Main_MPI
#ifdef MPI
  use UserCtrl, only : UserCtrl_ctrl
  use file_units
  use utilities
  use datasens
  use SolverSens
  use ForwardSolver
  use SensComp

  use Declaration_MPI
  use Sub_MPI
  ! use ioascii

  implicit none

  ! temporary EM fields, that are saved for efficiency - to avoid
  !  memory allocation & deallocation for each transmitter
  type(solnVector_t), save, private    :: e,e0
  type(rhsVector_t) , save, private    :: b0,comb
  type(grid_t), target, save, private  :: grid


Contains

!##########################  MPI_initialization  ############################

Subroutine Constructor_MPI

      implicit none
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, taskid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, total_number_of_Proc, ierr )
      call MPI_COMM_GROUP( MPI_COMM_WORLD, group_world, ierr )
      number_of_workers = total_number_of_Proc-1
      if (taskid .eq. 0 )then
          ! for debug
          ! write(6,*)'Reporting from Node #', taskid
          if (number_of_workers .lt. 1) then 
              write(6,*)'WARNING: not enough processes for MPI version'
              write(6,*)'please use at least TWO parallel processes'
          end if
      end if
      ! some translation works
      comm_world = MPI_COMM_WORLD
      rank_world = taskid
      size_world = total_number_of_Proc
      ! intra-group communicator including only the group members
      comm_local = MPI_COMM_NULL
      rank_local = -1
      size_local = -1
      ! branch communicator including only the master and leaders of each group
      ! by default every thread is in its own group, so everyone is in this
      ! communicator, at the beginning of the MPI program
      comm_leader= comm_world
      rank_leader= rank_world
      size_leader= size_world
      ! group identifier for all leaders
      group_leader = MPI_GROUP_NULL
      ! get to know your host name (useful for grouping)
      call MPI_GET_PROCESSOR_NAME( hostname_MPI,hostname_len,ierr )
      if (output_level .gt. 3) then 
      ! for debug
          write(6,*) 'Node #', taskid, ' hostname= ',                  &
    &         trim(hostname_MPI)
          write(6,*) ' hostname_len= ', hostname_len
      end if
      ! number of GPUs, for now we assume there are no GPUs
      ! the real number is determined after the split_MPI_group call
      ! for *master* we also always assume there is no gpu
      size_gpu = 0 
      ! tic-toc
      previous_time = MPI_Wtime()
End Subroutine Constructor_MPI
!----------------------------------------------------------------------------

!#######################  split MPI groups   ######################
Subroutine split_MPI_groups(nTx,nPol,group_sizes)
! split global MPI_COMM to a series of LOCAL communicators to implement
! two-layered (hierarchical) parallelization
! comm_world is the global communicator 
! comm_local is a local communicator within a group/node
! comm_leader is the exlusive branch communicator for the master and leaders 
! of each group 
!
! the default numbering of processes in this subroutine is linear - which means
! the Node groups are arranged as: 
! [1 2 3 4] [5 6 7 8] [9 10 11 12]
! this is different from the default "round-robin" approach like:
! [1 4 7 10] [2 5 8 11] [3 6 9 12]
! Note that any previous local and branch communicators will be reset after 
! calling this subroutine
! also note the group_sizes is not checked for compatability 
! so be careful when you use your own setup of groups
     implicit none
     integer, intent(in)    :: nTx, nPol
     ! cannot use numbers > nTx*nPol+1
     integer, intent(inout), pointer, dimension(:)  :: group_sizes 
     ! Local
     integer ::  Ngroup, Ntask, igroup
     integer ::  total_nTx_nPol, resultlen
     integer, allocatable, dimension(:) :: leaders

! comment from Naser
! Following the same prosedure as the parallelization over periods and 
! polarizations, 
! we would like to have the number of groups smaller than nTx*nPol + 1
! (this includes the master itself)
! The total number of worker is then distributed among these groups.
! this subroutine is thread-safe
! 
     ! check the basic input parameters
     if (size(group_sizes) .gt. (nTx*nPol+1)) then
         write(0,*) 'number of group is ', size(group_sizes) 
         write(0,*) 'number of tasks is ', (nTx*nPol+1)
         call errStop('ERROR: number of groups incompatible with the &
    &            number of tasks')
     else ! just group according to the input sizes
         ! start splitting here
         total_nTx_nPol = nTx*nPol ! total number of independent tasks
         ! firstly see how many groups we have here
         Ngroup = size(group_sizes)
     endif
     ! now reset previous branch and local communicators
     if (comm_leader .ne. MPI_COMM_NULL) then
         if (comm_leader .ne. MPI_COMM_WORLD) then
             call MPI_COMM_FREE(comm_leader,ierr)
         endif
     endif 
     ! free the local communicators 
     if (comm_local .ne. MPI_COMM_NULL) then
         call MPI_COMM_FREE(comm_local,ierr)
     endif 
     ! now start to see if we have enough workers 
     if (number_of_workers .lt. (Ngroup-1)) then 
         ! we do not have enough workers, each one is in its own group now
         Ngroup = number_of_workers+1
         ! re-allocate the groups
         deallocate(group_sizes)
         allocate(group_sizes(Ngroup))
         group_sizes = 1
         allocate(leaders(Ngroup))
     else
         if (sum(group_sizes) .gt. (number_of_workers+1)) then
         ! see if the sizes are compatiable
             call errStop('ERROR: grouping incompatible with the number &
    &            of workers')
         else if (group_sizes(1) .ne. 1) then 
         ! the master should always be by itself 
             call errStop('ERROR: grouping incompatible with the number &
    &            of master')
         else
         ! try allocate the leader group
             allocate(leaders(Ngroup))
         endif
     endif 
     ! start dividing groups accordingly...
     ! master is also treated a leader, of course
     leaders(1) = 0
     do igroup = 2, Ngroup      ! store the rank_world of all the leaders
         ! f95 feature here! 
         leaders(igroup) = SUM(group_sizes(1:igroup-1))
     end do 
     if (rank_world .eq. 0) then ! the master always has to be group 0
         igroup = 0
     else ! the below expression is a little hard to understand...
         do igroup = 2, Ngroup
             if (rank_world+1 .le. SUM(group_sizes(1:igroup))) then
                 exit 
             end if
         end do
         igroup = igroup -1
     end if
!    start spliting here
!    for generality, I have to refrain to use the Open-MPI only type
!    like OMPI_COMM_TYPE_NUMA 
     call MPI_COMM_SPLIT(comm_world,igroup,rank_world,comm_local,ierr)
     call MPI_COMM_RANK(comm_local, rank_local, ierr)
     call MPI_COMM_SIZE(comm_local, size_local, ierr)
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!    for debug
     if (output_level .ge. 3) then 
         write(6,100) rank_world, rank_local, igroup, size_local
100      format("Global rank: ",i4," local rank: ",i4,                      &
    &          " size of group #"  i4, " is ",i4)
     endif
     call reset_MPI_leader_group(comm_world, group_world, Ngroup,   &
    &       leaders, comm_leader, group_leader)
     if (rank_local .eq. 0) then
         call MPI_COMM_RANK(comm_leader, rank_leader, ierr)
         call MPI_COMM_SIZE(comm_leader, size_leader, ierr)
     end if 
end subroutine split_MPI_groups
!----------------------------------------------------------------------------

!#########################  reset MPI leader group  ###########################

subroutine reset_MPI_leader_group(comm_world,group_world, Ngroup, leaders &
     &                   ,comm_leader,group_leader)
!-----------------------------------------------------------------------------
     implicit none
     integer, intent(in) :: comm_world,group_world,Ngroup
     integer, intent(in) :: leaders(Ngroup)
     integer, intent(inout) :: comm_leader,group_leader
!    local variables
     integer ierr
!----------------------------------------------------------------------------
!    create a seperate group for the local group leaders
!    release the group if already existed
     if (group_leader .ne. MPI_GROUP_NULL) then
         call MPI_GROUP_FREE(group_leader,ierr)
     end if
     if (comm_leader .ne. MPI_COMM_NULL) then
         if (comm_leader .ne. MPI_COMM_WORLD) then
             call MPI_COMM_FREE(comm_leader,ierr)
         end if 
     end if
     call MPI_GROUP_INCL(group_world,  Ngroup, leaders ,                &
    &     group_leader,ierr)
!     create a seperate communicator for group_leader
     call MPI_COMM_CREATE(comm_world,group_leader,comm_leader,ierr)

end subroutine reset_MPI_leader_group
!----------------------------------------------------------------------------

!########################   set_group_sizes  ###############################
Subroutine set_group_sizes(nTx,nPol,comm_current,group_sizes,walltime)
! this is a quick and dirty way to divide workers into nG groups
! we will try to (dynamicaly) distribute more processing power to 
! groups that consumpt more time in the previous task 
! the ultimate target of load balancing is minimize the wall time...
! so we don't really care about how fast is the fastest task - but how slow is
! the slowest task
! this subroutine is thread safe - (not really)
! collective on comm_current
! def = 0 - simple settings, each proc has its own group (default setting)
!           use only one process for each FWD problem - this is very much
!           like your old 1-layer MPI. Each proc directly reports to 
!           the master. This could still benefit from the memory affinity
!           (for better memory bandwidth) on a same physical node 
!           you should use this (or #2) when GPUs are attached.
! def = 1 - topology-based settings, maps procs from a certain physical node 
!           in a group - efficient when you have a lot of nodes at service. 
!           use one entire node (many procs) to calculate one FWD problem
!           currently through PETSc, or FG only. In this case only the leaders
!           communicate with the master.
!           Mainly for debug purpose.
! def = 2 - equal group settings - the idea is to use a same number of
!           Lprocs per task - let's say we have 3 processes per Tx/pol, 
!           we will group Nworkers into  Nworkers/3 groups
!           currently through PETSc or FG only. 
! def = 3 - dynamic group settings - dynamicly distribute the
!           procs to different tasks (consider removing
!           this) - currently through PETSc or FG only. 
     implicit none
     integer, intent(in)                     :: nTx,nPol,comm_current
     integer, intent(inout),pointer,dimension(:)                 :: group_sizes
     real(kind=prec), optional, pointer, dimension(:), intent(in):: walltime
     ! local variables
     integer                                 :: i,j,k,a_size,nG,nMore,midx
     integer                                 :: nBand,nS,nM,nL,root=0
     real(kind=prec)                         :: time_sort(nTx*nPol)
     integer                                 :: def, nTask, total, diff
     integer, allocatable, dimension(:)      :: idx
     real(kind=prec), pointer, dimension(:)  :: cputime
     real(kind=prec)                         :: workload
     real(kind=prec)                         :: ratio
     integer, pointer, dimension(:)          :: host_sizes
     ! currently just read from the hard-coded preference
     def = para_method
     ! get the topology array
     if ((rank_world .eq. 0) .and. (output_level .gt. 2)) then
         write(6,*) '# determine the system topology...'
         write(6,*) '# and group the processes by physical nodes...'
     end if
     ! this is useful for "production" setups with
     ! hybrid CPU-GPU calculations
     call get_host_topology(nTx,nPol,host_sizes)
     ! number of total tasks
     nTask = nTx * nPol
     ! need to firstly decide how many groups do we have 
     if (rank_world.eq.0) then ! the root process determine the grouping
         ! try to find the slowest groups
         if (present(walltime).and.(size(walltime).eq.size(group_sizes))) then 
             ! overide the default settings
             ! calculate the appoximate workload for each group
             ! NOTE: 
             ! here we try to deduce the time supposed to be consumed 
             ! with one CPU process - stored in cputime array
             if ((def .eq. 3) .and. (walltime(2).eq.-1.0)) then 
                 ! first run
                 ! procs are divided into equal-sized groups
             elseif ((def .eq. 3) .and.(size(walltime) .lt. nTask)) then
                 def = 2 ! equal group size
                 ! no need to do dynamic, not enough groups
             elseif (def .eq. 3) then 
                 nG = size(walltime)
                 allocate(idx(nG))
                 idx = (/(i,i=1,nG,1)/)
                 allocate(cputime(nG))
                 ! cputime - unit is core second
                 cputime = 0
                 do i = 2,nG
                      cputime(i) = log(2.0+group_sizes(i))/log(3.0)
                 end do
                 cputime = walltime*cputime
                 ! estimate the unit workload - 
                 call QSort(cputime,idx) ! in ascend order 
                 ! fast <-----------------------> slow
                 ! write(6,*) 'sorted cputime =', cputime
                 ! total workload
                 workload=sum(cputime(2:nG))
                 ! average workload
                 workload=workload/sum(log(2.0+group_sizes(2:nG))/log(3.0))
             endif
         endif
         if (def .eq. 0) then
             ! just put each one in its own group, default
             nG = number_of_workers+1
         elseif (def .eq. 1) then
             ! group according to topography - put all procs from one host in 
             ! a group, useful for GPU 
             nG = size(host_sizes)
         elseif (def .eq. 2) then  
             ! equal group, but also consider the topology, i.e. one group 
             ! should not be distributed to more than one host
             ! for simplicity, we use a groups of 3 for each per/pol
             if (number_of_workers .gt. nTask) then
                 nG = nTask+1
             else
                 ! we don't have that many workers
#if defined(CUDA)||defined(HIP)
                 ! k = 2, hard coded here
                 nG = number_of_workers / 2 + 1
                 if (nG .le. 1) then
                     nG = 2
                 endif
#else
                 nG = number_of_workers+1
                 ! fall back to para_method = 0
                 def = 0
#endif
             endif
         elseif (def .eq. 3) then
             ! dynamic - give more procs to harder jobs
             ! note this is only possible when nG > number_of_workers+1
             nG = size(walltime)
         endif
     endif
     ! now distribute the number of groups to everybody
     call MPI_BCAST(def, 1, MPI_INTEGER, root, comm_current, ierr)
     if ((rank_world .eq. 0) .and. (output_level .gt. 2)) then 
         ! master 
         write(6,*) '# current grouping strategy is: ', def
     endif
     call MPI_BCAST(nG, 1, MPI_INTEGER, root, comm_current, ierr)
     if ((rank_world .eq. 0) .and. (output_level .gt. 2)) then 
         ! master 
         write(6,*) '# current number of groups is : ', nG
     endif
     allocate(group_sizes(nG))
     ! now determine the detailed grouping
     if (rank_world.eq.0) then ! the root process determine the grouping
         if (def.eq.0) then ! 1 group for everyone
             group_sizes = 1
         elseif (def .eq. 1) then ! 1 group for each (physical) machine
             group_sizes = 1 
             group_sizes(2:nG) = host_sizes(2:nG)
         elseif (def .eq. 2) then ! equally distribute the workers
             ! note that a_size is an integer
             a_size = number_of_workers/(nG-1)
             if (a_size.ge.1) then
                 ! we have enough workers for each task 
                 nMore = number_of_workers-((nG-1)*a_size) ! number of workers 
                 group_sizes = a_size ! every group gets a few workers...
                 ! some periods (hopefully the longer ones) get more (+1) procs
                 group_sizes(1) = 1
                 do j = nG-nMore+1,nG 
                     group_sizes(j) = group_sizes(j) + 1 
                 end do
             else ! we have fewer workers than tasks
                 ! each task gets one worker 
                 group_sizes = 1
             endif !asize
         elseif (def .eq. 3) then ! dynamic balancing the workers
             group_sizes = 0
             ! root always has its own group
             group_sizes(1) = 1
             ! for debug
             ! write(6,*) 'walltime = ', walltime(2:nG+1)
             ! write(6,*) 'cputime = ', cputime
             ! write(6,*) 'workload= ', workload
             do i = 1,nG
                 ! ratio should always be larger (or equal) than 1 
                 ratio = cputime(i)/workload
                 ! for debug
                 ! write(6,*) 'group ', idx(i), 'ratio =', ratio
                 ratio = 3.0**ratio - 2.0
                 ! write(6,*) 'group ', idx(i), 'ncpu =', ratio
                 group_sizes(idx(i)+1) = nint(ratio)
                 if (group_sizes(idx(i)+1).le.0) then
                     ! assign at least one cpu
                     group_sizes(idx(i)+1) = 1
                 endif
             enddo
             ! now start a bang-bang control procedure
             ! the idea is always giving more cpus to the slowest tasks
             nMore = number_of_workers - sum(group_sizes(2:nG+1))
             ! if we happend to come across no extra workers, then go on
             if (nMore.lt.0) then ! we don't have that many workers
                 do i = 1,nG
                     ! subtract workers from easier tasks
                     if (group_sizes(idx(i)+1).gt.1) then
                         group_sizes(idx(i)+1) = group_sizes(idx(i)+1)-1
                         nMore = nMore + 1
                     endif
                     if (nMore.ge.0) then
                         exit
                     endif
                 enddo
             elseif (nMore.gt.0) then ! too many workers
                 do i = nG,1,-1
                     ! send extra workers to most difficult tasks
                     group_sizes(idx(i)+1) = group_sizes(idx(i)+1)+1
                     nMore = nMore - 1
                     if (nMore.le.0) then
                         exit
                     endif
                 enddo
             endif
             ! now do a final check of compatibility
             if (sum(group_sizes(2:nG)).ne.number_of_workers) then
                 ! throw an error 
                 call errStop('ERROR: incompatible grouping settings &
    &                for workers')
             endif
             deallocate(cputime)
             deallocate(idx)
         else
             call errStop('ERROR: unknown grouping strategy')
         endif !def
         ! here we also check the host topology
         ! should avoid splitting a group onto different physical machines
         if (def.gt.1) then
             k = 2 ! this includes the first host (always 1)
             ! for debug
             ! write(6,*) 'host_sizes =', host_sizes
             ! write(6,*) 'old group_sizes =', group_sizes
             do j = 2,nG
             ! FIXME: there is a possiblity that one group is larger 
             ! than the number of procs of a host 
             ! the grouping could be problematic in this case
                 total = sum(group_sizes(1:j)) 
                 if (total .lt. sum(host_sizes(1:k))) then
                     ! do nothing and continue  
                 elseif (total .eq. sum(host_sizes(1:k))) then
                     ! go to the next host
                     k = k + 1
                 else ! we have one group spliting into two hosts 
                     i = 1
                     diff =  total - sum(host_sizes(1:k))
                     ! we have two choices -  
                     if ((group_sizes(j) - diff) .lt. diff) then 
                         ! better to add more workers to previous groups
                         diff = sum(host_sizes(1:k)) - sum(group_sizes(1:j-1))
                         do while (diff .gt. 0)
                             group_sizes(j-i) = group_sizes(j-i) + 1
                             diff = diff - 1
                             i = i + 1
                         end do
                     else
                         ! better to remove workers from previous groups
                         do while (diff .gt. 0)
                             group_sizes(j-i+1) = group_sizes(j-i+1) - 1
                             if (group_sizes(j-i+1) .lt. 1) then
                                 group_sizes(j-i+1) = 1
                             else
                                 diff = diff - 1
                             endif
                             i = i + 1
                         end do
                     endif
                     ! go to the next host
                     k = k + 1
                 endif
                 if (k .gt. size(host_sizes)) then
                     ! now just let it go 
                     exit
                 endif
             end do
             ! now check the last node
             diff =  sum(host_sizes) - sum(group_sizes)
             ! and last group
             j = nG
             ! start here
             i = 1
             if (diff .gt. 0) then
                 ! need to add more workers 
                 do while (diff .gt. 0)
                     group_sizes(j-i+1) = group_sizes(j-i+1) + 1
                     diff = diff - 1
                     i = i + 1
                 end do
             elseif (diff .lt. 0) then
                 ! need to reduce workers 
                 diff =  - diff
                 do while (diff .gt. 0)
                     group_sizes(j-i+1) = group_sizes(j-i+1) - 1
                     if (group_sizes(j-i+1) .lt. 1) then
                         group_sizes(j-i+1) = 1
                     else
                         diff = diff - 1
                     endif
                     i = i + 1
                 end do
             endif
             ! for debug
             ! write(6,*) 'new group_sizes =', group_sizes
             if (sum(host_sizes) .ne. sum(group_sizes)) then
                 ! throw an error
                 write(0,*) 'host sizes =', sum(host_sizes)
                 write(0,*) 'group sizes =', sum(group_sizes)
                 call errStop('ERROR: number of groups incompatible with &
    &                 size of hosts')
             endif
         endif
     endif ! rank = 0
     ! now spread group sizes (of size nG) to all
     call MPI_BCAST(group_sizes,nG,MPI_INTEGER,root,comm_current,ierr)
     ! for debug
     ! write(6,*) 'current group sizes are as follows:',group_sizes 
     deallocate(host_sizes)
     return
end subroutine set_group_sizes

!########################   Master_job_Regroup   ############################
Subroutine Master_job_Regroup(nTx, nPol, comm)
     implicit none
     integer, intent(in)                :: nTx,nPol
     integer, intent(in), optional      :: comm
     ! local
     integer                            :: idest,i,j
     DOUBLE PRECISION                   :: time_passed, now 
     DOUBLE PRECISION, pointer,dimension(:) :: time_buff ! nGroup by 1 array
     integer                            :: size_current, comm_current
     
     ! over-ride the default communicator, if needed
     if (present(comm)) then ! given communicator
         if (comm .eq. MPI_COMM_NULL) then
             comm_current = comm_world
         else
             comm_current = comm
         endif
     else
         comm_current = comm_world
     end if
     call MPI_COMM_SIZE( comm_current, size_current, ierr )

     ! if (.not.associated(group_sizes)) then
     !     allocate(group_sizes(nTx*nPol+1))
     ! elseif (size(group_sizes).ne.nTx*nPol+1) then
     !     deallocate(group_sizes)
     !     allocate(group_sizes(nTx*nPol+1))
     ! endif ! group_sizes(1) = 0
     now = MPI_Wtime()
     time_passed = now - previous_time
     previous_time =  now
     worker_job_task%what_to_do='REGROUP'
     worker_job_task%per_index=nTx
     worker_job_task%pol_index=nPol
     call create_worker_job_task_place_holder
     call Pack_worker_job_task
     do idest = 1, size_current - 1
         call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED, idest,     &
    &         FROM_MASTER, comm_current, ierr)
     enddo
     call gather_runtime(comm_current,time_passed,time_buff)
     !===== add the time_buff variable to test the dynamic load balancer=====!
     call set_group_sizes(nTx,nPol,comm_world,group_sizes,time_buff)
     call split_MPI_groups(nTx,nPol,group_sizes)
     ! for debug
     ! do i = 1,nTx
     !    do j = 1,nPol
     !        if (size(time_buff).lt.(i-1)*nPol+j+1) then
     !            ! only display time from available workers
     !            exit
     !        endif
     !        write(6,*) 'Tx:', i, 'Pol:', j,  'time elapsed:',            &
     !            time_buff((i-1)*nPol+j+1), 's'
     !        write(6,*)  group_sizes((i-1)*nPol+j+1), 'process(es) ',&
     !            'will be used in the next stage.'
     !    end do
     ! end do
     deallocate(time_buff)

end subroutine Master_job_Regroup

!----------------------------------------------------------------------------
!##########################  Master_job_fwdPred ############################

Subroutine Master_job_fwdPred(sigma,d1,eAll,comm)

     implicit none
     type(modelParam_t), intent(in)       :: sigma
     type(dataVectorMTX_t), intent(inout) :: d1
     type(solnVectorMTX_t), intent(inout) :: eAll
     integer, intent(inout),optional      :: comm
     integer nTx


     ! local variables
     Integer        :: iper, comm_current
     Integer        :: per_index,pol_index,stn_index,iTx,i,iDt,j
     character(80)  :: job_name

     ! nTX is number of transmitters;
     nTx = d1%nTx
     starttime = MPI_Wtime()
     ! over-ride the default communicator, if needed
     if (present(comm)) then ! given communicator
         if (comm .eq. MPI_COMM_NULL) then
             comm_current = comm_world
         else
             comm_current = comm
         endif
     else
         if (para_method.eq.0) then
             comm_current = comm_world
         else 
             comm_current = comm_leader
         end if
     end if

     ! First, distribute the current model to all workers
     call Master_job_Distribute_Model(sigma)
     ! call Master_job_Distribute_Data(d1)
     if(.not. eAll%allocated) then
     ! call deall(eAll)
     ! end if
         call create_solnVectorMTX(d1%nTx,eAll)
            do iTx=1,nTx
                call create_solnVector(grid,iTx,e0)
                call copy_solnVector(eAll%solns(iTx),e0)
            end do
     end if
     job_name= 'FORWARD'
     call Master_job_Distribute_Taskes(job_name,nTx,sigma,eAll,comm_current)

     ! Initialize only those grid elements on the master that are used in
     ! EMfieldInterp
     ! (obviously a quick patch, needs to be fixed in a major way)
     ! A.Kelbert 2018-01-28

     Call EdgeLength(grid, l_E)
     Call FaceArea(grid, S_F)

     ! Compute the model Responces
     do iTx=1,nTx
         do i = 1,d1%d(iTx)%nDt
             d1%d(iTx)%data(i)%errorBar = .false.
             iDt = d1%d(iTx)%data(i)%dataType
             do j = 1,d1%d(iTx)%data(i)%nSite
                 call dataResp(eAll%solns(iTx),sigma,iDt,d1%d(iTx)%data(i)%rx(j),d1%d(iTx)%data(i)%value(:,j), &
                           d1%d(iTx)%data(i)%orient(j))
             end do
         end do
     end do
     ! clean up the grid elements stored in GridCalc on the master node
     call deall_rvector(l_E)
     call deall_rvector(S_F)
     write(ioMPI,*)'FWD: Finished calculating for (', nTx , ') Transmitters '
     endtime=MPI_Wtime()
     time_used = endtime-starttime
     write(ioMPI,*)'FWD: TIME REQUIERED: ',time_used ,'s'
     call deall (e0)

end subroutine Master_job_fwdPred


!#########################   Master_job_Compute_J ##########################

Subroutine Master_job_calcJ(d,sigma,sens,eAll,comm)

     implicit none
     type(modelParam_t), intent(in)               :: sigma
     type(dataVectorMTX_t), intent(in)            :: d
     type(sensMatrix_t), pointer                  :: sens(:)
     type(solnVectorMTX_t), intent(in), optional  :: eAll
     integer, intent(in), optional                :: comm
     !Local
     logical        :: savedSolns
     Integer        :: iper,idt,istn, comm_current
     Integer        :: per_index,dt_index,stn_index,ipol1
     integer        :: nTx,nDt,nStn,dest,answers_to_receive
     integer        :: received_answers,nComp,nFunc,ii,iFunc,istat
     logical        :: isComplex  
     type(modelParam_t), pointer   :: Jreal(:),Jimag(:)           

     ! over-ride the default communicator, if needed
     if (present(comm)) then ! given communicator
         if (comm .eq. MPI_COMM_NULL) then
             comm_current = comm_world
         else
             comm_current = comm
         endif
     else
         if (para_method.eq.0) then
             comm_current = comm_world
         else 
             comm_current = comm_leader
         end if
     end if
     starttime = MPI_Wtime()
     ! now, allocate for sensitivity values, if necessary
     if (.not. associated(sens)) then
         call create_sensMatrixMTX(d, sigma, sens)
     endif
     ! Check if an Esoln is passed    
     savedSolns = present(eAll)  
     if (.not. savedSolns )then
         ! call Master_job_fwdPred(sigma,d,eAll)
     end if
      
     dest=0
     worker_job_task%what_to_do='COMPUTE_J' 
     ! now loop over Periods
     nTx = d%nTx  
     per_index=0 
     do iper=1,nTx
         per_index=per_index+1
         worker_job_task%per_index= per_index
         call get_nPol_MPI(eAll%solns(per_index)) 
         ! now loop over data types
         dt_index=0
         nDt=d%d(iper)%nDt
         do idt = 1,nDt
             dt_index=dt_index+1
             worker_job_task%data_type= d%d(iper)%data(idt)%dataType
             worker_job_task%data_type_index= dt_index
             
             ! now loop over stations
             stn_index=0
             nStn= d%d(iper)%data(idt)%nSite
             do istn=1,nStn
                 stn_index=stn_index+1
                 worker_job_task%Stn_index= stn_index  
                 dest=dest+1
                 call create_worker_job_task_place_holder
                 call Pack_worker_job_task
                 call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,dest,& 
    &                 FROM_MASTER, comm_current, ierr)
                 which_per=per_index
    
                 do ipol1=1,nPol_MPI
                     which_pol=ipol1
                     call create_e_param_place_holder(eAll%solns(which_per))
                     call Pack_e_para_vec(eAll%solns(which_per))
                     call MPI_SEND(e_para_vec, Nbytes, MPI_PACKED,     &
    &                     dest,FROM_MASTER, comm_current, ierr) 
                 end do
                 write(ioMPI,8550) per_index ,dt_index,stn_index,dest  
                 if (dest .ge. number_of_workers) then 
                     goto 10
                 end if
             end do
         end do
     end do
      
     10    continue   
                 
     ! Count hom many rows of J we are going to receive
     answers_to_receive=0
     do iper=1,nTx 
         nDt=d%d(iper)%nDt          
         do idt = 1,nDt 
             nStn= d%d(iper)%data(idt)%nSite   
             do istn=1,nStn
                 answers_to_receive=answers_to_receive+1
             end do
         end do
     end do
! Start the PING PONG procedure:
! 1- Receive an answer.
! 2 -Check if all stations, all data types and all periods haven been sent.
! 3- Send indicies if required.  

     received_answers = 0
     write(6,*)'answers_to_receive=', answers_to_receive
     do while (received_answers .lt. answers_to_receive) 

!    Receive worker INFO:                       
     call create_worker_job_task_place_holder
     call MPI_RECV(worker_job_package, Nbytes, MPI_PACKED ,MPI_ANY_SOURCE,&
    &              FROM_WORKER,comm_current,STATUS, ierr)
     call Unpack_worker_job_task

     who=worker_job_task%taskid
     which_per=worker_job_task%per_index
     which_dt=worker_job_task%data_type_index                      
     which_stn=worker_job_task%Stn_index
     write(ioMPI,8552) which_per ,which_dt,which_stn,who             
!    Receive results from a worker:
!    Store the result in sens
     nComp = d%d(which_per)%data(which_dt)%nComp           
     isComplex = d%d(which_per)%data(which_dt)%isComplex

     if (isComplex) then
         !  data are complex; one sensitivity calculation can be
         !   used for both real and imaginary parts
         if (mod(nComp,2).ne.0) then
              call errStop('for complex data # of components must be even &
    &            in calcJ')
         endif
         nFunc = nComp/2
     else
         !  data are treated as real: full sensitivity computation is required
         !   for each component
         nFunc = nComp
     endif
     allocate(Jreal(nFunc),STAT=istat)
     allocate(Jimag(nFunc),STAT=istat)    
! allocate and initialize sensitivity values
     do iFunc = 1,nFunc
         ! this makes a copy of modelParam, then zeroes it
         Jreal(iFunc) = sigma
         call zero(Jreal(iFunc))
         Jimag(iFunc) = sigma
         call zero(Jimag(iFunc))
     enddo

     do iFunc = 1,nFunc    
         !real part
         call create_model_param_place_holder(Jreal(iFunc))
         call MPI_RECV(sigma_para_vec, Nbytes, MPI_PACKED, who,          &
    &          FROM_WORKER, comm_current, STATUS, ierr)
         call unpack_model_para_values(Jreal(iFunc))
         !image part
         call create_model_param_place_holder(Jimag(iFunc))
         call MPI_RECV(sigma_para_vec, Nbytes, MPI_PACKED, who,          &
    &          FROM_WORKER, comm_current, STATUS, ierr)
         call unpack_model_para_values(Jimag(iFunc))
     end do    
     ! store in the full sensitivity matrix
     ii = 1
     do iFunc = 1,nFunc
         if(isComplex) then
             sens(which_per)%v(which_dt)%dm(ii,which_stn)   = Jreal(iFunc)
             sens(which_per)%v(which_dt)%dm(ii+1,which_stn) = Jimag(iFunc)
             ii = ii + 2
         else
             ! for real data, throw away the imaginary part
             sens(which_per)%v(which_dt)%dm(ii,which_stn)   = Jreal(iFunc)
             ii = ii + 1
         endif
     enddo

     ! deallocate temporary vectors
     do iFunc = 1,nFunc
         call deall_modelParam(Jreal(iFunc))
         call deall_modelParam(Jimag(iFunc))
     enddo
     deallocate(Jreal, STAT=istat)
     deallocate(Jimag, STAT=istat)
            
     received_answers=received_answers+1
      
     if (Per_index ==  nTx .and. dt_index == d%d(Per_index)%nDt .and. &
    &      stn_index == d%d(Per_index)%data(dt_index)%nSite ) goto 300 
      
     stn_index=stn_index+1 
     ! Check if we sent everything
     if (stn_index .gt. d%d(Per_index)%data(dt_index)%nSite ) then
         dt_index=dt_index+1
         stn_index=1   
     elseif ( stn_index .le. d%d(Per_index)%data(dt_index)%nSite) then
         dt_index=dt_index
     end if
     
     if (dt_index .gt. d%d(Per_index)%nDt ) then
         per_index=per_index+1
         dt_index=1   
     elseif (dt_index .le. d%d(Per_index)%nDt) then
         per_index=per_index
     end if
     
     if (Per_index .gt. nTx ) goto 300
         worker_job_task%Stn_index= stn_index  
         worker_job_task%data_type_index=dt_index
         worker_job_task%data_type=d%d(per_index)%data(dt_index)%dataType
         worker_job_task%per_index= per_index
         write(ioMPI,8551) per_index ,worker_job_task%data_type_index,  &
    &         stn_index,who    
     !   Send Indices to who (the worker who just send back an answer)
         call create_worker_job_task_place_holder
         call Pack_worker_job_task
         call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,who,       &
    &    FROM_MASTER, comm_current, ierr)
     
         do ipol1=1,nPol_MPI
             which_pol=ipol1
             call create_e_param_place_holder(eAll%solns(per_index))
             call Pack_e_para_vec(eAll%solns(per_index))
             call MPI_SEND(e_para_vec, Nbytes, MPI_PACKED, who,         &
    &    FROM_MASTER, comm_current, ierr) 
         end do          
 300     continue                                   
     end do

8550 FORMAT('COMPUTE_J: Send Per. # ',i5, ' , data type # ',i5,        &
    &       ' and Site #',i5,' to node # ',i5  )
8551 FORMAT('COUNT--> COMPUTE_J: Send Per. # ',i5, ' , data type # ',  &
    &       i5, ' and Site #',i5,' to node # ',i5  )
8552 FORMAT('COMPUTE_J: RECV Per. # ',i5, ' , data type # ',i5,        &
    &       ' and Site #',i5,' from node # ',i5  )

     endtime=MPI_Wtime()
     time_used = endtime-starttime
     !DONE: Received soln for all transmitter from all nodes
     write(ioMPI,*)'COMPUTE_J: Finished computing for (',             &
    &     answers_to_receive , ';[Transmitters*data type*site]) '
     endtime=MPI_Wtime()
     time_used = endtime-starttime
     write(ioMPI,*)'COMPUTE_J: TIME REQUIERED: ',time_used ,'s'

end subroutine Master_job_calcJ


!#########################    Master_job_JmultT ############################
Subroutine Master_job_JmultT(sigma,d,dsigma,eAll,s_hat,comm)

     implicit none
     type(modelParam_t), intent(in)       :: sigma
     type(dataVectorMTX_t), intent(in)    :: d
     type(modelParam_t), intent(Out)      :: dsigma
     type(solnVectorMTX_t), intent(in), optional                     :: eAll
     type(modelParam_t),intent(inout),pointer,dimension(:), optional :: s_hat
     integer, intent(in),optional         :: comm
  
     ! Local
     type(modelParam_t)           :: dsigma_temp
     type(modelParam_t)           :: Qcomb
     type(solnVectorMTX_t)        :: eAll_out 
     type(solnVectorMTX_t)        :: eAll_temp
     type(dataVectorMTX_t)        :: d_temp
     
     logical        :: savedSolns,returne_m_vectors
     Integer        :: iper,ipol,nTx,iTx
     Integer        :: per_index,pol_index,stn_index
     character(80)  :: job_name,file_name
     Integer        :: comm_current

     savedSolns = present(eAll)
     returne_m_vectors= present(s_hat)
     ! over-ride the default communicator, if needed
     if (present(comm)) then ! given communicator
         if (comm .eq. MPI_COMM_NULL) then
             comm_current = comm_world
         else
             comm_current = comm
         endif
     else
         if (para_method.eq.0) then
             comm_current = comm_world
         else 
             comm_current = comm_leader
         endif
     end if
     ! nTX is number of transmitters;
     nTx = d%nTx
     if(.not. eAll_temp%allocated) then
         call create_solnVectorMTX(d%nTx,eAll_temp)
     end if 
     if (.not. savedSolns )then
         d_temp=d
         call Master_job_fwdPred(sigma,d_temp,eAll_temp)
     else
         eAll_temp=eAll 
     end if
     if(.not. eAll_out%allocated) then
         call create_solnVectorMTX(d%nTx,eAll_out)
         do iTx=1,nTx
             call create_solnVector(grid,iTx,e0)
             call copy_solnVector(eAll_out%solns(iTx),e0) 
             call deall (e0)  
         end do 
     end if 


     if (returne_m_vectors) then
         if (.not. associated(s_hat)) then
             allocate(s_hat(nTx))
         end if 
         do iper=1,nTx
             s_hat(iper)=sigma
             call zero(s_hat(iper))
         end do
     end if
     starttime = MPI_Wtime()
     ! First ditribute both model parameters and data
     call Master_job_Distribute_Model(sigma)
     call Master_job_Distribute_Data(d)

     dsigma_temp = sigma
     dsigma      = sigma
     call zero(dsigma_temp)
     call zero(dsigma)
     Qcomb = sigma
     call zero(Qcomb)
     job_name= 'JmultT'
     call Master_job_Distribute_Taskes(job_name,nTx,sigma,eAll_out,      &
    &     comm_current, eAll_temp)
       
     file_name='e0.soln'
     !call write_solnVectorMTX(10,file_name,eAll_temp)
     file_name='e.soln'
     !call write_solnVectorMTX(20,file_name,eAll_out)

     do iper=1,nTx
         !e0=eAll%solns(iper)  
         !e =eAll_out%solns(iper)
         call PmultT(eAll_temp%solns(iper),sigma,eAll_out%solns(iper)    &
    &         ,dsigma_temp)
         call QmultT(eAll_temp%solns(iper),sigma,d%d(iper),Qcomb)
         call scMultAdd(ONE,Qcomb,dsigma_temp)
         if (returne_m_vectors) then
             s_hat(iper)=dsigma_temp
         end if
         call linComb_modelParam(ONE,dsigma,ONE,dsigma_temp,dsigma)
     end do

     endtime=MPI_Wtime()
     time_used = endtime-starttime
     !DONE: Received soln for all transmitter from all nodes
     write(ioMPI,*)'JmultT: Finished calculating for (', d%nTx ,        &
    &     ') Transmitters '
     endtime=MPI_Wtime()
     time_used = endtime-starttime
     write(ioMPI,*)'JmultT: TIME REQUIERED: ',time_used ,'s'

     ! clean up
     call deall_modelParam(dsigma_temp)
     call deall_modelParam(Qcomb)
     call deall (eAll_out)
     call deall (eAll_temp) 
     ! call deall (e0)   
     call deall_dataVectorMTX(d_temp)
   
end Subroutine Master_job_JmultT


!########################    Master_job_Jmult ##############################
Subroutine Master_job_Jmult(mHat,m,d,eAll,comm)

     implicit none
     type(dataVectorMTX_t), intent(inout)              :: d
     type(modelParam_t), intent(in)                    :: mHat,m
     type(solnVectorMTX_t), intent(in), optional       :: eAll
     integer, intent(in), optional                     :: comm
     !Local
     integer        :: nTx,nTot,m_dimension,iDT,iTx,ndata,ndt
     logical        :: savedSolns
     Integer        :: iper, comm_current
     Integer        :: per_index,pol_index,stn_index
     type(dataVector_t)                                :: d1,d2
     type(dataVectorMTX_t)                             :: d_temp
     type(solnVectorMTX_t)                             :: eAll_out 
     type(solnVectorMTX_t)                             :: eAll_temp 
     character(80)                                     :: job_name
  
   
     savedSolns = present(eAll)
     ! over-ride the default communicator, if needed
     if (present(comm)) then ! given communicator
         if (comm .eq. MPI_COMM_NULL) then
             comm_current = comm_world
         else
             comm_current = comm
         endif
     else
         if (para_method.eq.0) then
             comm_current = comm_world
         else 
             comm_current = comm_leader
         end if
     end if
     ! nTX is number of transmitters;
     nTx = d%nTx
     ! nTot total number of data points
     nTot = countData(d) !d%Ndata
     starttime = MPI_Wtime()
     !  initialize the temporary data vectors
     d1 = d%d(1)
     d2 = d%d(1) 
     if (.not. eAll_out%allocated) then
         call create_solnVectorMTX(d%nTx,eAll_out)
         do iTx=1,nTx
             call create_solnVector(grid,iTx,e0)
             call copy_solnVector(eAll_out%solns(iTx),e0) 
             call deall (e0)  
         end do 
     end if     
     if (.not. eAll_temp%allocated) then
         call create_solnVectorMTX(d%nTx,eAll_temp)
     end if
     if (.not. savedSolns )then
         d_temp=d
         call Master_job_fwdPred(m,d_temp,eAll_temp,comm_current)
     else
         eAll_temp=eAll 
     end if
   ! First distribute m, mHat and d
     call Master_job_Distribute_Model(m,mHat,comm_current)
     call Master_job_Distribute_Data(d,comm_current)

     job_name= 'Jmult'
     call Master_job_Distribute_Taskes(job_name,nTx,m,eAll_out, &
    &    comm_current, eAll_temp)

     do iper=1,nTx
         !e0=eAll%solns(iper)  
         !e =eAll_out%solns(iper)
         d1 = d%d(iper)
         d2 = d%d(iper)
         call Lmult(eAll%solns(iper)  ,m,eAll_out%solns(iper),d1)
         call Qmult(eAll%solns(iper)  ,m,mHat,d2)
         call linComb_dataVector(ONE,d1,ONE,d2,d%d(iper))
     end do

     call deall_dataVector(d1)
     call deall_dataVector(d2)       
     call deall_dataVectorMTX(d_temp)
     call deall (eAll_out)
     call deall (eAll_temp)
     ! call deall (e0)  
     !DONE: Received soln for all transmitter from all nodes
     write(ioMPI,*)'Jmult: Finished calculating for (', d%nTx ,         &
    &         ') Transmitters '
     endtime=MPI_Wtime()
     time_used = endtime-starttime
     write(ioMPI,*)'Jmult: TIME REQUIERED: ',time_used ,'s'
 
end Subroutine Master_job_Jmult

!###################### Master_job_Distribute_Data #########################
Subroutine Master_job_Distribute_Data(d, comm)
     implicit none
     type(dataVectorMTX_t), intent(in)       :: d
     integer, intent(in), optional           :: comm
     ! local 
     integer                                 :: nTx
     integer                                 :: size_current, comm_current

     ! over-ride the default communicator, if needed
     if (present(comm)) then ! given communicator
         if (comm .eq. MPI_COMM_NULL) then
             comm_current = comm_world
         else
             comm_current = comm
         endif
     else
         if (para_method.eq.0) then
             comm_current = comm_world
         else 
             comm_current = comm_leader
         end if
     end if
     call MPI_COMM_SIZE( comm_current, size_current, ierr )
     nTx=d%nTx
     do dest=1,size_current-1
         worker_job_task%what_to_do='Distribute nTx'
         call create_worker_job_task_place_holder
         call Pack_worker_job_task
         call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,dest,      &
    &         FROM_MASTER, comm_current, ierr)
     end do
     ! distribute the number of Txs
     call MPI_BCAST(nTx,1, MPI_INTEGER,0, comm_current,ierr)

     do dest=1,size_current-1
         worker_job_task%what_to_do='Distribute Data'
         call create_worker_job_task_place_holder
         call Pack_worker_job_task
         call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,dest,      &
    &         FROM_MASTER, comm_current, ierr)
     end do
     call create_data_vec_place_holder(d)
     call Pack_data_para_vec(d)
     ! distribute the packed data vector
     ! note that the *_para_vec are all public in declarition_MPI module
     call MPI_BCAST(data_para_vec,Nbytes, MPI_PACKED,0, comm_current,ierr)
     if (associated(sigma_para_vec)) then
        deallocate(sigma_para_vec)
     end if

end Subroutine Master_job_Distribute_Data


!######################## Master_job_Distribute_Model #######################
Subroutine Master_job_Distribute_Model(sigma,delSigma,comm)
     implicit none
     type(modelParam_t), intent(in)           :: sigma
     type(modelParam_t), intent(in), optional :: delSigma
     integer, intent(in), optional            :: comm
     !local
     type(modelParam_t)                       :: sigma_temp
     Integer,pointer:: buffer(:)
     integer        :: buffer_size,m_dimension
     logical        :: send_delSigma
     Integer        :: iper, comm_current, size_current

     send_delSigma = present(delSigma)

     ! over-ride the default communicator, if needed
     if (present(comm)) then ! given communicator
         if (comm .eq. MPI_COMM_NULL) then
             comm_current = comm_world
         else
             comm_current = comm
         endif
     else
         if (para_method.eq.0) then
             comm_current = comm_world
         else 
             comm_current = comm_leader
         end if
     end if
     call MPI_COMM_SIZE( comm_current, size_current, ierr )

     if (send_delSigma) then !send pertubation model 
         do dest=1,size_current-1
             worker_job_task%what_to_do='Distribute delSigma'
             call create_worker_job_task_place_holder
             call Pack_worker_job_task
             call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,dest,    &
    &             FROM_MASTER, comm_current, ierr)
         end do
         call create_model_param_place_holder(delSigma)
         call pack_model_para_values(delSigma)
         call MPI_BCAST(sigma_para_vec,Nbytes, MPI_PACKED,0,              &
    &         comm_current,ierr)

         do dest=1,size_current-1
             worker_job_task%what_to_do='Distribute Model'
             call create_worker_job_task_place_holder
             call Pack_worker_job_task
             call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,dest,    &
    &             FROM_MASTER, comm_current, ierr)
         end do

         call create_model_param_place_holder(sigma)
         call pack_model_para_values(sigma)
         call MPI_BCAST(sigma_para_vec,Nbytes, MPI_PACKED,0,             &
    &         comm_current,ierr)

     else! do not send purtubation model
         do dest=1,size_current-1
             worker_job_task%what_to_do='Distribute Model'
             call create_worker_job_task_place_holder
             call Pack_worker_job_task
             call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,dest,   &
    &             FROM_MASTER, comm_current, ierr)
         end do
         call create_model_param_place_holder(sigma)
         call pack_model_para_values(sigma)
         call MPI_BCAST(sigma_para_vec,Nbytes, MPI_PACKED,0,             &
              comm_current,ierr)
     end if
     if (associated(sigma_para_vec)) then
         deallocate(sigma_para_vec)
     end if

end Subroutine Master_job_Distribute_Model

!######################## Master_job_Distribute_eAll #########################
Subroutine Master_job_Distribute_eAll(d,eAll, comm)
     implicit none
     type(dataVectorMTX_t), intent(in)    :: d
     type(solnVectorMTX_t), intent(in)    :: eAll
     integer, intent(in), optional        :: comm
     !Local
     integer                              :: nTx,nTot
     Integer                              :: iper
     Integer                              :: comm_current, size_current

     ! over-ride the default communicator, if needed
     if (present(comm)) then ! given communicator
         if (comm .eq. MPI_COMM_NULL) then
             comm_current = comm_world
         else
             comm_current = comm
         endif
     else
         if (para_method.eq.0) then
             comm_current = comm_world
         else 
             comm_current = comm_leader
         end if
     end if
     call MPI_COMM_SIZE( comm_current, size_current, ierr)
     nTx = d%nTx
     do dest=1,size_current-1
         worker_job_task%what_to_do='Distribute eAll'
         call create_worker_job_task_place_holder
         call Pack_worker_job_task
         call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,dest,       &
    &         FROM_MASTER, comm_current, ierr)
     end do
     call MPI_COMM_SIZE( comm_current, size_current, ierr )
     do iper=1,d%nTx
         which_per=iper
         do dest=1,size_current-1
             call create_eAll_param_place_holder(e0)
             call Pack_eAll_para_vec(e0)
             call MPI_SEND(eAll_para_vec, Nbytes, MPI_PACKED,dest,        &
    &             FROM_MASTER,comm_current, ierr)
         end do
     end do
     if (associated(eAll_para_vec)) then
         deallocate(eAll_para_vec)
     end if
end Subroutine Master_job_Distribute_eAll

!######################## Master_job_Collect_eAll ############################

Subroutine Master_job_Collect_eAll(d,eAll, comm)
     implicit none
     type(dataVectorMTX_t), intent(in)    :: d
     type(solnVectorMTX_t), intent(inout) :: eAll
     integer, intent(inout),optional      :: comm
     ! Local
     integer nTx,nTot,iTx
     Integer                              :: iper
     Integer                              :: comm_current, size_current
   
     ! over-ride the default communicator, if needed
     if (present(comm)) then ! given communicator
         if (comm .eq. MPI_COMM_NULL) then
             comm_current = comm_world
         else
             comm_current = comm
         endif
     else
         if (para_method.eq.0) then
             comm_current = comm_world
         else 
             comm_current = comm_leader
         end if
     end if
     nTx = d%nTx
     call MPI_COMM_SIZE( comm_current, size_current, ierr )

     if(.not. eAll%allocated) then
         call create_solnVectorMTX(d%nTx,eAll)
     else if(d%nTx .ne. eAll%nTx) then
         call errStop('dimensions of eAll and d do not agree in          &
    &        Master_job_Collect_eAll')
     endif

     do iTx=1,nTx
         call create_solnVector(grid,iTx,e0)
         call copy_solnVector(eAll%solns(iTx),e0)
     end do

     do iper=1,d%nTx
         worker_job_task%what_to_do='Send eAll to Master'
         worker_job_task%per_index=iper
         who=1
         which_per=iper
         call create_worker_job_task_place_holder
         call Pack_worker_job_task
         call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,who,        &
    &         FROM_MASTER, comm_current, ierr)
         call create_eAll_param_place_holder(e0)
         call MPI_RECV(eAll_para_vec, Nbytes, MPI_PACKED, who,           &
    &         FROM_WORKER, comm_current, STATUS, ierr)
         call Unpack_eAll_para_vec(e0)
     end do

     if (associated(eAll_para_vec)) then
         deallocate(eAll_para_vec)
     end if

end Subroutine Master_job_Collect_eAll

!######################## Master_job_keep_prev_eAll #########################
subroutine Master_job_keep_prev_eAll(comm)
     implicit none
     integer, intent(in), optional :: comm
     !local
     Integer                       :: comm_current, size_current
     ! over-ride the default communicator, if needed
     if (present(comm)) then ! given communicator
         if (comm .eq. MPI_COMM_NULL) then
             comm_current = comm_world
         else
             comm_current = comm
         endif
     else
         if (para_method.eq.0) then
             comm_current = comm_world
         else 
             comm_current = comm_leader
         end if
     end if
     call MPI_COMM_SIZE( comm_current, size_current, ierr )
     do dest=1,size_current-1
         worker_job_task%what_to_do='keep_prev_eAll'
         call create_worker_job_task_place_holder
         call Pack_worker_job_task
         call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,dest,       &
    &         FROM_MASTER, comm_current, ierr)
     end do

end subroutine Master_job_keep_prev_eAll


!################### Master_job_Distribute_userdef_control###################
Subroutine Master_job_Distribute_userdef_control(ctrl, comm)
     implicit none
     type(userdef_control), intent(in) :: ctrl
     integer, intent(in), optional :: comm
     !local
     character(20)                 :: which_proc
     integer                       :: comm_current, size_current
     ! over-ride the default communicator, if needed
     if (present(comm)) then ! given communicator
         if (comm .eq. MPI_COMM_NULL) then
             comm_current = comm_world
         else
             comm_current = comm
         endif
     else
         comm_current = comm_world
     end if
     call MPI_COMM_SIZE( comm_current, size_current, ierr )
     which_proc='Master'
     call check_userdef_control_MPI (which_proc,ctrl)
     call create_userdef_control_place_holder
     call pack_userdef_control(ctrl)
     do dest=1,size_current -1 
         call MPI_SEND(userdef_control_package,Nbytes, MPI_PACKED,dest,  &
              FROM_MASTER, comm_current, ierr)
     end do
     if (associated(userdef_control_package)) then
         deallocate(userdef_control_package)
     end if

end Subroutine Master_job_Distribute_userdef_control

!#########################  Master_job_Clean Memory ##########################
Subroutine Master_job_Clean_Memory(comm)

     implicit none
     integer, intent(in), optional :: comm
     !local
     integer                       :: comm_current, size_current
     ! over-ride the default communicator, if needed
     if (present(comm)) then ! given communicator
         if (comm .eq. MPI_COMM_NULL) then
             comm_current = comm_world
         else
             comm_current = comm
         endif
     else
         if (para_method.eq.0) then
             comm_current = comm_world
         else 
             comm_current = comm_leader
         end if
     end if
     call MPI_COMM_SIZE( comm_current, size_current, ierr )
     write(ioMPI,*)'Sending Clean memory message to all nodes'
     do dest=1,size_current-1
         worker_job_task%what_to_do='Clean memory'
         call create_worker_job_task_place_holder
         call Pack_worker_job_task
         call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,dest,       &
    &         FROM_MASTER, comm_current, ierr)
         call MPI_RECV(worker_job_package, Nbytes, MPI_PACKED ,dest,     &
    &         FROM_WORKER, comm_current,STATUS, ierr)
         call Unpack_worker_job_task
         write(ioMPI,*)'Node  :',worker_job_task%taskid, ' status=  ',   &
    &         trim(worker_job_task%what_to_do)
     end do

end Subroutine Master_job_Clean_Memory

!#######################    Master_job_Stop_MESSAGE #########################
Subroutine Master_job_Stop_MESSAGE(comm)

     implicit none
     integer, intent(in), optional :: comm
     !local
     integer                       :: comm_current,size_current
     ! over-ride the default communicator, if needed
     if (present(comm)) then ! given communicator
         if (comm .eq. MPI_COMM_NULL) then
             comm_current = comm_world
         else
             comm_current = comm
         endif
     else
         if (para_method.eq.0) then
             comm_current = comm_world
         else 
             comm_current = comm_leader
         end if
     end if
     call MPI_COMM_SIZE( comm_current, size_current, ierr )

     write(ioMPI,*)'FWD: Sending stop message to all nodes'
     do dest=1,size_current-1 
         worker_job_task%what_to_do='STOP'
         call create_worker_job_task_place_holder
         call Pack_worker_job_task
         call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,dest,       &
    &         FROM_MASTER, comm_current, ierr)
         call MPI_RECV(worker_job_package, Nbytes, MPI_PACKED ,dest,     &
    &         FROM_WORKER, comm_current,STATUS, ierr)
         call Unpack_worker_job_task
         write(ioMPI,*)'Node/Group :',worker_job_task%taskid,            &
    &         ' status=  ', trim(worker_job_task%what_to_do)
     end do
     ! all leaders reported 
end Subroutine Master_job_Stop_MESSAGE

!********************** Master Distribute Tasks *****************************

subroutine Master_job_Distribute_Taskes(job_name,nTx,sigma,eAll_out, &
    &   comm, eAll_in)

     implicit none
     character(80) , intent(in)                          :: job_name
     Integer       , intent(in)                          :: nTx
     type(modelParam_t),intent(in)                       :: sigma
     type(solnVectorMTX_t), intent(in), optional         :: eAll_in
     integer, intent(in), optional                       :: comm
     type(solnVectorMTX_t), intent(inout), optional      :: eAll_out     
     !Local
     Integer        :: iper,ipol,ipol1,ijob,total_jobs
     Integer        :: bcounter, tcounter, robin
     Integer        :: per_index,pol_index,des_index
     logical        :: keep_soln,savedSolns,ascend
     integer        :: comm_current, size_current
     double precision :: now, just, time_passed

     savedSolns = present(eAll_in)
     ! over-ride the default communicator, if needed
     if (present(comm)) then ! given communicator
         if (comm .eq. MPI_COMM_NULL) then
             comm_current = comm_world
         else
             comm_current = comm
         endif
     else
         comm_current = comm_world
     end if
     call get_nPol_MPI(eAll_out%solns(1)) 
     if (rank_local.eq.-1) then ! first run!
     ! run initial regroup -- note this requires the comm to be
     ! comm_world as it is the only comm that works in this stage
         call Master_job_Regroup(nTx,nPol_MPI,comm)
         if (para_method.gt.0) then
             comm_current = comm_leader
         else
             comm_current = comm_world
         end if
     endif
     call MPI_COMM_SIZE( comm_current, size_current, ierr )
     who = 0
     worker_job_task%what_to_do=trim(job_name) 
     call count_number_of_messages_to_RECV(eAll_out)
     total_jobs = answers_to_receive
     ! counter to locate the task 
     bcounter = 1
     tcounter = 1
     ! send easier tasks first
     ascend = .true.
     ! here we use an mechanism to determine the number of physical nodes
     ! and to distribute the tasks through a round-robin way
     ! to prevent the longer tasks from clustering on a same node
     robin = total_jobs/(size_leader-1)
     if (robin.lt.1) then
         robin = 1
     end if
     ! for debug
     ! if (trim(job_name).eq. 'JmultT') then
     !     write(6,*)'robin = ', robin
     !    write(6,*)'total jobs = ', total_jobs
     ! endif
     do ijob=1,total_jobs !loop through all jobs
         ! loop through all jobs, until we run out of workers
         who=who+1
         ! count from head
         call find_next_job(nTx,total_jobs,tcounter,ascend,eAll_out,&
    &        per_index, pol_index)
         tcounter = tcounter + robin
         if (tcounter.gt.total_jobs) then 
             bcounter = bcounter + 1
             tcounter = bcounter
         end if
         worker_job_task%per_index = per_index
         worker_job_task%pol_index = pol_index
         call create_worker_job_task_place_holder
         call Pack_worker_job_task
         call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED, who, &
    &        FROM_MASTER, comm_current, ierr)
         if (trim(job_name).eq. 'JmultT' .or. trim(job_name).eq.     &
    &        'Jmult') then
             ! In case of JmultT and Jmult the compleate eAll solution 
             ! (the background solution) must be sent to each node.
             ! In the 3D MT case there are two polarisations.
             ! In the 3D CSEM case there are one polarisation (for now).
             ! In the 2D MT case there are one polarisation.
             which_per=per_index
             do ipol1=1,nPol_MPI
                 which_pol=ipol1
                 call create_e_param_place_holder(eAll_in%solns(which_per))
                 call Pack_e_para_vec(eAll_in%solns(which_per))
                 call MPI_SEND(e_para_vec, Nbytes, MPI_PACKED, who,  &
    &             FROM_MASTER, comm_current, ierr) 
             end do   
         end if  
         write(ioMPI,'(a10,a16,i5,a8,i5,a11,i5)')trim(job_name),     &
    &         ': Send Per. # ',per_index , ' : Pol #', pol_index,    &
    &         ' to node # ', who
         if (who .ge. (size_current-1)) then
             ! break, we are done distributing the first batch
             ! now wait for the first group to report back
                 goto 10
         end if
     end do
                   
10   continue

     !answers_to_receive = nTx*nPol_MPI
     received_answers = 0
     do while (received_answers .lt. answers_to_receive)
     ! now loop until all tasks are finished
         call create_worker_job_task_place_holder
         call MPI_RECV(worker_job_package, Nbytes, MPI_PACKED ,          &
    &         MPI_ANY_SOURCE, FROM_WORKER,comm_current,STATUS, ierr)
         call Unpack_worker_job_task
         who=worker_job_task%taskid
         which_per=worker_job_task%per_index
         which_pol=worker_job_task%pol_index
                  
         call create_e_param_place_holder(eAll_out%solns(which_per))
         call MPI_RECV(e_para_vec, Nbytes, MPI_PACKED, who,FROM_WORKER,  &
    &         comm_current, STATUS, ierr)
         ! call get_nPol_MPI(eAll_out%solns(which_per)) 
         ! if (nPol_MPI==1)  which_pol=1

         call Unpack_e_para_vec(eAll_out%solns(which_per))

         write(ioMPI,'(a10,a16,i5,a8,i5,a11,i5)')trim(job_name) ,        &
    &   ': Receive Per # ',which_per ,' and Pol # ', which_pol ,' from ',&
    &    who 
         write(6,'(a10,a16,i5,a8,i5,a11,i5)')trim(job_name) ,            &
    &   ': Receive Per # ',which_per ,' and Pol # ', which_pol ,' from ',&
    &    who 

         received_answers=received_answers+1
                  
         write(ioMPI,'(a10,a16,i5,a8,i5)') trim(job_name),               &
    &   ': TX finished # ', received_answers, ' out of ',                &
    &    answers_to_receive 
         write(6,'(a10,a16,i5,a8,i5)') trim(job_name),                   &
    &   ': TX finished # ', received_answers, ' out of ',                &
    &    answers_to_receive 

         ! Check if we send all transmitters and polarizations, if not 
         ! then send the next transmitter to the worker who is free now....
         ! This part is very important if we have less worker groups than
         ! transmitters. 

         if (received_answers .ge. answers_to_receive) then 
             goto 1500 ! continue, just wait all the tasks to end 
         else
             if (bcounter.gt.robin) then
                 ! we have sent enough jobs, now just wait
                 goto 1500 
             end if
             ! send new jobs to folks 
             call find_next_job(nTx,total_jobs,tcounter,ascend,eAll_out,&
    &            per_index, pol_index)
             tcounter = tcounter + robin
             if (tcounter.gt.total_jobs) then 
                 bcounter = bcounter + 1
                 tcounter = bcounter
             end if
         end if
         write(6,*) 'Per_index ',Per_index,'pol_index ',pol_index,      &
    &         'Going to send'
          
         worker_job_task%per_index= per_index
         worker_job_task%pol_index= pol_index
         worker_job_task%what_to_do= trim(job_name) 

         call create_worker_job_task_place_holder
         call Pack_worker_job_task
         call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,who,       &
    &         FROM_MASTER, comm_current, ierr)

         if (trim(job_name).eq. 'JmultT' .or. trim(job_name).eq.'Jmult')&
    &        then
             which_per=per_index
             call get_nPol_MPI(eAll_out%solns(per_index)) 
             do ipol1=1,nPol_MPI
                 which_pol=ipol1
                 call create_e_param_place_holder(eAll_in%solns(      &
    &                 which_per))
                 call Pack_e_para_vec(eAll_in%solns(which_per))
                 call MPI_SEND(e_para_vec, Nbytes, MPI_PACKED,        &
    &                 who, FROM_MASTER, comm_current, ierr) 
             end do  
         end if 
         write(ioMPI,'(a10,a16,i5,a8,i5,a11,i5)')trim(job_name),      &
    &         ': Sent Per. # ',per_index , ' : Pol #', pol_index,     &
    &         ' to node # ',who
         write(6,'(a10,a16,i5,a8,i5,a11,i5)')trim(job_name),          &
    &         ': Sent Per. # ',per_index , ' : Pol #', pol_index,     &
    &         ' to node # ',who
1500 continue                                   
     end do

     if (trim(job_name).eq. 'JmultT') then
         ! only regroup after the ADJ calculation
         ! avoid the extra overheads
         call Master_job_Regroup(nTx,nPol_MPI,comm)
     endif
     if (associated(eAll_para_vec)) then
         deallocate(eAll_para_vec)
     end if

     if (associated(e_para_vec)) then
         deallocate(e_para_vec)
     end if
     !call deall(e0)

end subroutine Master_job_Distribute_Taskes

!##################   find next job -- from back or front   ##################
    Subroutine find_next_job(nTx,total_jobs,counter,fromhead, eAll_out, &
     &    per_index,pol_index)
     implicit none
     Integer, intent(in)                      :: nTx
     Integer, intent(in)                      :: total_jobs 
     Integer, intent(in)                      :: counter
     Logical, intent(in)                      :: fromhead
     type(solnVectorMTX_t), intent(in)        :: eAll_out     
     Integer, intent(out)                     :: per_index,pol_index
     !Local
     Integer                                  :: iper, ipol, jobs
     ! a very silly (literally) subroutine to find the next job 
     ! normally there should be nTx*nPol jobs
     ! but it is possible (in principle) that we don't have 2 pols at one per
     ! note this is not thread-safe - should only be called by the master 
     jobs = 0
     if (fromhead) then
         ! count from start
         head: do iper=1,nTx !loop through all transmitters 
             call get_nPol_MPI(eAll_out%solns(iper)) 
             do ipol=1,nPol_MPI ! loop through all pols
                 ! count the total tasks
                 jobs=jobs+1
                 if (jobs.ge.counter) then
                     per_index=iper
                     pol_index=ipol
                     exit head
                 end if
             end do
         end do head
     else
         ! count from back
         tail: do iper=nTx,1,-1 !loop through all transmitters 
             call get_nPol_MPI(eAll_out%solns(iper)) 
             do ipol=nPol_MPI,1,-1 ! loop through all pols
                 ! count the total tasks
                 jobs=jobs+1
                 if (jobs.ge.counter) then
                     per_index=iper
                     pol_index=ipol
                     exit tail
                 end if
             end do
         end do tail
     end if
    end subroutine find_next_job

!###################   Worker_job: High Level Subroutine   ###################
Subroutine Worker_job(sigma,d)
     ! subroutine for *all* worker jobs -
     ! the general idea (from Naser, I believe) is:  
     ! 
     ! 1. a worker firstly prepare a worker_job info (structure datatype)
     ! 2. listen to the boss and find out what to do (receive worker_job) 
     ! 3. receive the data package from the boss
     ! 4. do the job and return the data to the boss
     ! 5. return to 1 and loop, untill the STOP instruction
     ! 
     ! but this is already growing too large... 
     ! 
     ! an obvious problem for this, is each proc must maintain its own set 
     ! of model/data structures, which takes a lot of memory
     ! TODO: test shared memory model - it is probably better to just declare
     ! everything private (like in the OMP paradigm) and only share the 
     ! grid, sigma, etc. 
     implicit none

     type(modelParam_t),intent(inout)       :: sigma
     type(dataVectorMTX_t) ,intent(inout)   :: d
   
     ! Local 
     type(modelParam_t)                     :: delSigma
     type(userdef_control)                  :: ctrl
     Integer                                :: nTx,m_dimension,ndata
     Integer                                :: itx, ndt, dt, dt_index
   
     Integer                                :: iper,ipol,i,des_index
     Integer                                :: per_index,per_index_pre 
     Integer                                :: pol_index, stn_index
     Integer                                :: eAll_vec_size
     Integer                                :: comm_current, rank_current
     Integer                                :: cpu_only_ranks
     Integer,pointer,dimension(:)           :: group_sizes
     character(20)                          :: which_proc
     character(80)                          :: paramType,previous_message
   
     ! sensitivity
     type(modelParam_t), pointer            :: Jreal(:),Jimag(:)
     integer                                :: nComp,nFunc,iFunc,istat
     logical                                :: isComplex  
     type(sparseVector_t), pointer          :: L(:)
     type(modelParam_t), pointer            :: Qreal(:),Qimag(:)
     logical                                :: Qzero
     type(orient_t)               :: orient
 
     ! 2019.05.08, Liu Zhongyin, add isite for rx in dataBlock_t
     integer                       :: isite

      
     ! time
     DOUBLE PRECISION                       :: time_passed, now, just
     DOUBLE PRECISION, pointer,dimension(:) :: time_buff
       
     nTx=d%nTx
     recv_loop=0
     previous_message=''
     write(node_info,'(a5,i3.3,a4)') 'node[',taskid,']:  '

     do  ! the major loop
         recv_loop=recv_loop+1
         ! prepare the job info structure 
         call create_worker_job_task_place_holder
         ! firstly need to know who pays the check for you... 
         if (para_method .eq. 0) then
             ! override - everyone reports back to master
             comm_current = comm_world
             rank_current = rank_world
             write(6,'(a12,a35)') node_info,                              &
    &            ' Waiting for a message from Master'
             write(6,'(a12, a28, f12.6)') node_info,        &
                 ' total time elapsed (sec): ', time_passed
         else
             if (rank_local .eq. 0) then
                 ! group leader, reports to master
                 comm_current = comm_leader
                 rank_current = rank_leader
                 write(6,'(a12,a35)') node_info,                          &
    &                ' Waiting for a message from Master'
                 write(6,'(a12, a22, f12.6)') node_info,        &
    &                ' time elapsed (sec): ', time_passed
             elseif (rank_local .gt. 0) then 
                 ! group worker, reports to group leader
                 comm_current = comm_local
                 rank_current = rank_local
                 write(6,'(a12,a35)') node_info,                          &
    &        ' Waiting for a message from Leader'
             else ! -1
                 ! uninitialized, reports to master first 
                 comm_current = comm_world
                 rank_current = rank_world
                 write(6,'(a12,a35)') node_info,                          &
    &        ' Waiting for a message from Someone'
             end if
         end if
         ! reset the timer
         previous_time = now
         now = MPI_Wtime()
         ! receive the job info from someone
         call MPI_RECV(worker_job_package, Nbytes, MPI_PACKED ,0,     &
    &        FROM_MASTER,comm_current,STATUS, ierr)
         ! unpack message including what to do and other info. required 
         ! (i.e per_index_stn_index, ...,etc)
         call Unpack_worker_job_task
         write(6,'(a12,a12,a30,a16,i5)') node_info,' MPI TASK [',         &
    &        trim(worker_job_task%what_to_do),'] received from ',        &
    &        STATUS(MPI_SOURCE)
         ! for debug
         ! write(6,*) 'source = ', MPI_SOURCE
         ! write(6,*) 'tag = ', MPI_TAG
         ! write(6,*) 'err = ', MPI_ERROR
         ! write(6,*) node_info,' MPI INFO [keep soln = ',               &
         !            (worker_job_task%keep_E_soln), &
         !            '; several TX = ',worker_job_task%several_Tx,']'

         if (trim(worker_job_task%what_to_do) .eq. 'FORWARD') then
             ! forward modelling
             per_index=worker_job_task%per_index
             pol_index=worker_job_task%pol_index
             if ((size_local.gt.1).and.(para_method.gt.0).and.          &
    &            (rank_local.eq.0)) then 
                 ! group leader passing the command to workers
                 do des_index=1, size_local-1
                     call MPI_SEND(worker_job_package,Nbytes,           &
    &                    MPI_PACKED, des_index,  FROM_MASTER, comm_local, &
    &                    ierr)
                 end do
             end if
             if ((rank_local.eq.0) .or. (para_method.eq.0)) then
                 ! leader prepares the basic data structure
                 call initSolver(per_index,sigma,grid,e0,size_local)
                 call set_e_soln(pol_index,e0)
                 call fwdSetup(per_index,e0,b0)
             else
                 ! worker just fills in some dummy parameters
                 iTx = 1
                 call create_solnVector(grid,iTx,e0)
                 call set_e_soln(pol_index,e0)
                 call create_rhsVector(grid,iTx,b0)
             end if
             if ((para_method.eq.0).or.(size_local.eq.1)) then
                 ! you are on your own, bro!
                 call fwdSolve(per_index,e0,b0,device_id) 
             else
#ifdef PETSC
                 call fwdSolve(per_index,e0,b0,device_id,comm_local) 
#elif defined(FG)
                 call fwdSolve(per_index,e0,b0,device_id,comm_local) 
#else
                 if (rank_local.eq.0) then
                     call fwdSolve(per_index,e0,b0,device_id) 
                 else
                     write(6,'(a12,a18)') node_info, ' hanging around...'
                     write(6,*) ' WARNING: more than enough CPU(s) detected'
                     write(6,*) ' Please consider recompiling with PETSC   '
                     write(6,*) ' or FG configurations '
                 endif
#endif
             end if

             if ((rank_local.eq.0) .or. (para_method.eq.0)) then


                 ! leader reports back to master
                 call create_worker_job_task_place_holder 
                 worker_job_task%taskid=rank_current
                 call Pack_worker_job_task
                 call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,0,   &
    &                FROM_WORKER, comm_current, ierr)

                 ! Create e0_temp package (one Period and one Polarization) 
                 ! and send it to the master
                 write(0,*) "User Ctrl: ", UserCtrl_ctrl % storeSolnsInFile, UserCtrl_ctrl % prefix
                 if (UserCtrl_ctrl % storeSolnsInFile) then
                     call write_solnVector(e0, UserCtrl_ctrl % prefix)
                 else
                     which_pol=1
                     call create_e_param_place_holder(e0)
                     call Pack_e_para_vec(e0)
                     call MPI_SEND(e_para_vec, Nbytes, MPI_PACKED, 0,         &
        &                FROM_WORKER, comm_current, ierr)
                     call reset_e_soln(e0)
                 end if
             end if

             ! so long!
             now = MPI_Wtime()
             time_passed =  now - previous_time
             previous_time = now

         elseif (trim(worker_job_task%what_to_do) .eq. 'COMPUTE_J') then
             ! compute (explicit) J
             if ((size_local.gt.1).and.(para_method.gt.0).and.          &
    &            (rank_local.eq.0)) then 
                 ! group leader passing the command to workers
                 do des_index=1, size_local-1
                     call MPI_SEND(worker_job_package,Nbytes,           &
    &                    MPI_PACKED, des_index,  FROM_MASTER, comm_local, &
    &                    ierr)
                 end do
             end if
             per_index=worker_job_task%per_index
             stn_index=worker_job_task%stn_index
             dt_index=worker_job_task%data_type_index
             dt=worker_job_task%data_type
             nComp = d%d(per_index)%data(dt_index)%nComp           
             isComplex = d%d(per_index)%data(dt_index)%isComplex
   
             if (isComplex) then
             !   data are complex; one sensitivity calculation can be
             !   used for both real and imaginary parts
                 if (mod(nComp,2).ne.0) then
                     call errStop('for complex data # of components must  &
    &                    be even in calcJ')
                 endif
                 nFunc = nComp/2
             else
             !   data are treated as real: full sensitivity computation is
             !   required  for each component
                 nFunc = nComp
             endif
             if ((rank_local.eq.0) .or. (para_method.eq.0)) then
                 ! only leader allocates basic sensitivity structure
                 allocate(Jreal(nFunc),STAT=istat)
                 allocate(Jimag(nFunc),STAT=istat)
                 ! allocate and initialize sensitivity values
                 do iFunc = 1,nFunc
                     ! this makes a copy of modelParam, then zeroes it
                     Jreal(iFunc) = sigma
                     call zero(Jreal(iFunc))
                     Jimag(iFunc) = sigma
                     call zero(Jimag(iFunc))
                 enddo
             end if
             ! initialize solver  
             if ((rank_local.eq.0) .or. (para_method.eq.0)) then
                 ! leader prepares the basic data structure
                 call initSolver(per_index,sigma,grid,e0,size_local,e,comb)
                 call get_nPol_MPI(e0)
             else
                 ! worker just fills in some dummy parameters
                 iTx = 1
                 call create_solnVector(grid,iTx,e0)
                 call set_e_soln(pol_index,e0)
                 call create_rhsVector(grid,iTx,comb)
             end if
             if ((rank_local.eq.0) .or. (para_method.eq.0)) then !leader
                 write(6,'(a12,a18,i5,a12)') node_info,' Start Receiving ',&
    &                orginal_nPol, ' from Master'
                 do ipol=1,nPol_MPI 
                     which_pol=ipol
                     call create_e_param_place_holder(e0)
                     call MPI_RECV(e_para_vec, Nbytes, MPI_PACKED, 0,      &
    &                    FROM_MASTER,comm_current, STATUS, ierr)
                     call Unpack_e_para_vec(e0)
                 end do
                 write(6,'(a12,a18,i5,a12)') node_info,' Finished Receiving '&
    &                , orginal_nPol, ' from Master'
                 allocate(L(nFunc),STAT=istat)
                 allocate(Qreal(nFunc),STAT=istat)
                 allocate(Qimag(nFunc),STAT=istat)
                 do iFunc=1,nFunc
                     call create_sparseVector(e0%grid,per_index,L(iFunc))
                 end do
                 ! Initialize only those grid elements on the leader that 
                 ! are used in EMfieldInterp
                 ! (obviously a quick patch, needs to be fixed in a major way)
                 ! A.Kelbert 2018-01-28
                 Call EdgeLength(e0%grid, l_E)
                 Call FaceArea(e0%grid, S_F)
                 ! compute linearized data functional(s) : L
                 !call Lrows(e0,sigma,dt,stn_index,L)
                 ! 2022.10.06, Liu Zhongyin, Add Azimuth
                 call Lrows(e0,sigma,dt,stn_index,orient,L)
                 ! compute linearized data functional(s) : Q
                 call Qrows(e0,sigma,dt,stn_index,Qzero,Qreal,Qimag)
                 ! clean up the grid elements stored in GridCalc on the 
                 ! leader node
                 call deall_rvector(l_E)
                 call deall_rvector(S_F) 
             end if
             ! loop over functionals  (e.g., for 2D TE/TM impedances 
             ! nFunc = 1)
             do iFunc = 1,nFunc
                 ! solve transpose problem for each of nFunc functionals
                 call zero_rhsVector(comb)
                 if ((rank_local.eq.0) .or. (para_method.eq.0)) then
                     call add_sparseVrhsV(C_ONE,L(iFunc),comb)
                 end if
                 if ((para_method.eq.0).or.(size_local.eq.1)) then
                     ! you are on your own, tovarishch
                     call sensSolve(per_index,TRN,e,comb,device_id)
                 else
#ifdef PETSC
                   call sensSolve(per_index,TRN,e,comb,device_id,comm_local)
#elif defined(FG)
                   call sensSolve(per_index,TRN,e,comb,device_id,comm_local)
#else
                   if (rank_local.eq.0) then
                     call sensSolve(per_index,TRN,e,comb,device_id)
                   else
                     write(6,'(a12,a18)') node_info, ' hanging around...'
                     write(6,*) ' WARNING: more than enough CPU(s) detected'
                     write(6,*) ' Please consider recompiling with PETSC   '
                     write(6,*) ' or FG configurations '
                   endif
#endif
                 end if
                 if ((rank_local.eq.0) .or. (para_method.eq.0)) then
                     ! multiply by P^T and add the rows of Q
                     call PmultT(e0,sigma,e,Jreal(iFunc),Jimag(iFunc))
                     if (.not. Qzero) then
                         call scMultAdd(ONE,Qreal(iFunc),Jreal(iFunc))
                         call scMultAdd(ONE,Qimag(iFunc),Jimag(iFunc))
                     endif
                     ! deallocate temporary vectors
                     call deall_sparseVector(L(iFunc))
                     call deall_modelParam(Qreal(iFunc))
                     call deall_modelParam(Qimag(iFunc))
                 end if
             enddo  ! iFunc
             
             if ((rank_local.eq.0) .or. (para_method.eq.0)) then
                 ! leader do all other calculations
                 ! deallocate local arrays
                 deallocate(L,STAT=istat)
                 deallocate(Qreal,STAT=istat)
                 deallocate(Qimag,STAT=istat)
                 ! call Jrows(per_index,dt_index,stn_index,sigma,e0,
                 ! Jreal,Jimag)        
                 ! Create worker job package and send it to the master
                 ! now send the calculated J back to master
                 call create_worker_job_task_place_holder
                 worker_job_task%taskid=rank_current
                 call Pack_worker_job_task
                 call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,0,   &
    &                FROM_WORKER, comm_current, ierr)
                 do iFunc = 1,nFunc 
                 ! Create worker model package for Jreal and send it to 
                 ! the master       
                     call create_model_param_place_holder(Jreal(iFunc))
                     call pack_model_para_values(Jreal(iFunc))
                     call MPI_SEND(sigma_para_vec, Nbytes, MPI_PACKED, 0, &
    &                    FROM_WORKER, comm_current, ierr)
                 ! Create worker model  package for Jimag and send it to
                 ! the master       
                     call create_model_param_place_holder(Jimag(iFunc))
                     call pack_model_para_values(Jimag(iFunc))
                     call MPI_SEND(sigma_para_vec, Nbytes, MPI_PACKED, 0, &
    &                    FROM_WORKER, comm_current, ierr)
                 end do    
             end if
             ! Das vidania
             now = MPI_Wtime()
             time_passed = now - previous_time
             previous_time = now
             
         elseif (trim(worker_job_task%what_to_do) .eq. 'JmultT') then

             ! calculate JmultT (a.k.a. adjoint)
             per_index=worker_job_task%per_index
             pol_index=worker_job_task%pol_index
             if ((size_local.gt.1).and.(para_method.gt.0).and.          &
    &            (rank_local.eq.0)) then 
                 ! group leader passing the command to workers
                 do des_index=1, size_local-1
                     call MPI_SEND(worker_job_package,Nbytes,           &
    &                    MPI_PACKED, des_index,  FROM_MASTER, comm_local, &
    &                    ierr)
                 end do
             end if
             if ((rank_local.eq.0) .or. (para_method.eq.0)) then
                 ! leader prepares the basic data structure
                 ! firstly fill some dummy values
                 iTx = 1
                 call create_solnVector(grid,iTx,e0)
                 call get_nPol_MPI(e0)
                 write(6,'(a12,a18,i5,a12)') node_info,                   &
    &                ' Start Receiving ' , orginal_nPol, ' from Master'
                 do ipol=1,nPol_MPI 
                     which_pol=ipol
                     call create_e_param_place_holder(e0)
                     call MPI_RECV(e_para_vec, Nbytes, MPI_PACKED, 0,     &
    &                    FROM_MASTER,comm_current, STATUS, ierr)
                     call Unpack_e_para_vec(e0)
                 end do
                 call initSolverWithOutE0(per_index,sigma,grid,size_local,&
                     e,comb)
                 write(6,'(a12,a18,i5,a12)') node_info,                   &
    &                ' Finished Receiving ', orginal_nPol, ' from Master'
                 call LmultT(e0,sigma,d%d(per_index),comb)
                 call set_e_soln(pol_index,e)
             else
                 ! worker just fills in some dummy parameters
                 iTx = 1
                 call create_solnVector(grid,iTx,e)
                 call set_e_soln(pol_index,e)
                 call create_rhsVector(grid,iTx,comb)
             end if

             if ((para_method.eq.0).or.(size_local.eq.1)) then
                 ! you are on your own, amigo!
                 call sensSolve(per_index,TRN,e,comb,device_id)
             else
#ifdef PETSC
                 call sensSolve(per_index,TRN,e,comb,device_id,comm_local) 
#elif defined(FG)
                 call sensSolve(per_index,TRN,e,comb,device_id,comm_local) 
#else
                 if (rank_local.eq.0) then 
                     call sensSolve(per_index,TRN,e,comb,device_id)
                 else
                     write(6,'(a12,a18)') node_info, ' hanging around...'
                     write(6,*) ' WARNING: more than enough CPU(s) detected'
                     write(6,*) ' Please consider recompiling with PETSC   '
                     write(6,*) ' or FG configurations '
                 end if
#endif
             end if
             call reset_e_soln(e)
             if ((rank_local.eq.0) .or. (para_method.eq.0)) then !leader


                 ! leader reports back to master
                 call create_worker_job_task_place_holder
                 worker_job_task%taskid=rank_current
                 call Pack_worker_job_task
                 call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,0,   &
    &                FROM_WORKER, comm_current, ierr)
                 which_pol=1

                 if (UserCtrl_ctrl % storeSolnsInFile) then
                    call write_solnVector(e, trim(UserCtrl_ctrl % prefix)//".JmultT")
                 else
                     call create_e_param_place_holder(e)
                     call Pack_e_para_vec(e)
                     call MPI_SEND(e_para_vec, Nbytes, MPI_PACKED, 0,         &
        &                FROM_WORKER, comm_current, ierr)
                 end if
                 !deallocate(e_para_vec,worker_job_package)
             end if
             ! hasta la vista!
             now = MPI_Wtime()
             time_passed = now - previous_time
             previous_time = now
                    
         elseif (trim(worker_job_task%what_to_do) .eq. 'Jmult') then

             ! calculate Jmult (probably only used in DCG)
             ! need to test it thoroughly before merging back to stable
             per_index=worker_job_task%per_index
             pol_index=worker_job_task%pol_index
             worker_job_task%taskid=rank_current
             if ((size_local.gt.1).and.(para_method.gt.0).and.          &
    &            (rank_local.eq.0)) then 
                 ! group leader passing the command to workers
                 do des_index=1, size_local-1
                     call MPI_SEND(worker_job_package,Nbytes,           &
    &                    MPI_PACKED, des_index,  FROM_MASTER, comm_local, &
    &                    ierr)
                 end do
             end if

             if ((rank_local.eq.0) .or. (para_method.eq.0)) then
                 ! leader prepares the basic data structure
                 call initSolver(per_index,sigma,grid,e0,size_local,e,comb)
                 write(6,'(a12,a18,i5,a12)') node_info,                  &
    &                ' Start Receiving    ' , orginal_nPol, ' from Master'
                 do ipol=1,orginal_nPol 
                     which_pol=ipol
                     call create_e_param_place_holder(e0)
                     call MPI_RECV(e_para_vec, Nbytes, MPI_PACKED, 0,         &
    &                    FROM_MASTER,comm_current, STATUS, ierr)
                     call Unpack_e_para_vec(e0)
                 end do
                 write(6,'(a12,a18,i5,a12)') node_info,                  &
    &            ' Finished Receiving ' , orginal_nPol, ' from Master'
                 call Pmult(e0,sigma,delSigma,comb)
             else
                 ! worker just fills in some dummy parameters
                 iTx = 1
                 call create_solnVector(grid,iTx,e)
                 call create_rhsVector(grid,iTx,comb)
             end if
             call set_e_soln(pol_index,e)
             if ((para_method.eq.0).or.(size_local.eq.1)) then
                 ! you are on your own, aibo!
                 call sensSolve(per_index,FWD,e,comb,device_id)
             else
#ifdef PETSC
                 call sensSolve(per_index,FWD,e,comb,device_id,comm_local) 
#elif defined(FG)
                 call sensSolve(per_index,TRN,e,comb,device_id,comm_local) 
#else
                 if (rank_local.eq.0) then
                     call sensSolve(per_index,FWD,e,comb,device_id)
                 else 
                     write(6,'(a12,a18)') node_info, ' hanging around...'
                     write(6,*) ' WARNING: more than enough CPU(s) detected'
                     write(6,*) ' Please consider recompiling with PETSC   '
                     write(6,*) ' or FG configurations '
                 end if
#endif
             end if
             call reset_e_soln(e)

             if ((rank_local.eq.0) .or. (para_method.eq.0)) then !leader

                 ! leader reports back to master
                 call create_worker_job_task_place_holder
                 call Pack_worker_job_task
                 call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,0,   &
    &                FROM_WORKER, comm_current, ierr)
                 ! send the results back
                 if (UserCtrl_ctrl % storeSolnsInFile) then
                    call write_solnVector(e, trim(UserCtrl_ctrl % prefix)//".Jmult")
                 else
                     which_pol=1
                     call create_e_param_place_holder(e)
                     call Pack_e_para_vec(e)
                     call MPI_SEND(e_para_vec, Nbytes, MPI_PACKED, 0,         &
        &                FROM_WORKER, comm_current, ierr)
                end if
             end if
             ! Aba yo!
             now = MPI_Wtime()
             time_passed = now - previous_time
             previous_time = now

         elseif (trim(worker_job_task%what_to_do) .eq. 'Distribute nTx')&
    &            then
             if ((size_local.gt.1).and.(para_method.gt.0).and.          &
    &            (rank_local.eq.0)) then 
                 ! group leader passing the command to workers
                 do des_index=1, size_local-1
                     call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,&
    &                    des_index, FROM_MASTER, comm_local, ierr)
                 end do
             end if
             ! get the nTx from master
             call MPI_BCAST(nTx, 1, MPI_INTEGER, 0, comm_current, ierr)
             if ((size_local.gt.1).and.(para_method.gt.0).and.          &
    &            (rank_local.eq.0)) then 
                 ! now broadcast nTx to the fellow workers
                 call MPI_BCAST(nTx, 1, MPI_INTEGER, 0, comm_local, ierr)
             endif
             d%nTx=nTx
         elseif (trim(worker_job_task%what_to_do) .eq. 'Distribute Data')&
    &            then
             if ((size_local.gt.1).and.(para_method.gt.0).and.          &
    &            (rank_local.eq.0)) then 
                 ! group leader passing the command to workers
                 do des_index=1, size_local-1
                     call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,&
    &                    des_index, FROM_MASTER, comm_local, ierr)
                 end do
             endif
             call create_data_vec_place_holder(d)
             ! get the vector from master
             call MPI_BCAST(data_para_vec, Nbytes, MPI_PACKED, 0,     &
    &             comm_current,ierr)
             if ((size_local.gt.1).and.(para_method.gt.0).and.          &
    &            (rank_local.eq.0)) then 
                 ! now broadcast data to the fellow workers
                 call MPI_BCAST(data_para_vec, Nbytes, MPI_PACKED, 0,     &
    &                 comm_local,ierr)
             endif
             call UnPack_data_para_vec(d)
             ! clear the packed data vector
             if (associated(data_para_vec)) then
                 deallocate(data_para_vec)
             endif
              
         elseif (trim(worker_job_task%what_to_do) .eq. 'Distribute eAll') &
    &            then
             ! note that a worker (non-leader in a group) should never get this
             do iper=1,d%nTx
                 which_per=iper
                 call create_eAll_param_place_holder(e0)
                 call MPI_RECV(eAll_para_vec, Nbytes, MPI_PACKED ,0,  &
    &                 FROM_MASTER,comm_current, STATUS, ierr)
                 call Unpack_eAll_para_vec(e0)
             end do
             eAll_exist=.true.
             if (associated(eAll_para_vec)) then
                 deallocate(eAll_para_vec)
             endif

         elseif (trim(worker_job_task%what_to_do) .eq. 'Distribute Model')&
    &            then
             if ((size_local.gt.1).and.(para_method.gt.0).and.          &
    &            (rank_local.eq.0)) then 
                 ! group leader passing the command to workers
                 do des_index=1, size_local-1
                     call MPI_SEND(worker_job_package,Nbytes,             &
    &                   MPI_PACKED, des_index,  FROM_MASTER, comm_local,  &
    &                   ierr)
                 end do
             endif
             call create_model_param_place_holder(sigma)
             ! get the vector from master
             call MPI_BCAST(sigma_para_vec, Nbytes, MPI_PACKED, 0,     &
    &            comm_current,ierr)
             if ((size_local.gt.1).and.(para_method.gt.0).and.          &
    &            (rank_local.eq.0)) then 
                 ! now broadcast data to the fellow workers
                 call MPI_BCAST(sigma_para_vec, Nbytes, MPI_PACKED, 0,     &
    &                 comm_local,ierr)
             end if
             call unpack_model_para_values(sigma)
             if (associated(sigma_para_vec)) then
                 deallocate(sigma_para_vec)
             endif

         elseif (trim(worker_job_task%what_to_do).eq.'Distribute delSigma')&
    &            then
             ! distribute the conductivity pertubation to all procs
             if ((size_local.gt.1).and.(para_method.gt.0).and.          &
    &            (rank_local.eq.0)) then 
                 !passing the command to workers
                 do des_index=1, size_local-1
                     call MPI_SEND(worker_job_package,Nbytes,             &
    &                   MPI_PACKED, des_index,  FROM_MASTER, comm_local,  &
    &                   ierr)
                 end do
             end if
             call create_model_param_place_holder(sigma)
             ! get the vector from master
             call MPI_BCAST(sigma_para_vec, Nbytes, MPI_PACKED, 0,     &
    &            comm_current,ierr)
             if ((size_local.gt.1).and.(para_method.gt.0).and.          &
    &            (rank_local.eq.0)) then 
                 ! now broadcast data to the fellow workers
                 call MPI_BCAST(sigma_para_vec, Nbytes, MPI_PACKED, 0, &
    &                comm_local,ierr)
             endif
             call copy_ModelParam(delSigma,sigma)
             call unpack_model_para_values(delSigma)
             if (associated(sigma_para_vec)) then
                 deallocate(sigma_para_vec)
             endif
         elseif (trim(worker_job_task%what_to_do) .eq.                  &
    &            'Send eAll to Master' ) then
             !note that a worker (non-leader in a group) should never get this
             per_index=worker_job_task%per_index
             worker_job_task%taskid=taskid
             which_per=per_index
             call create_eAll_param_place_holder(e0)
             call Pack_eAll_para_vec(e0)
             call MPI_SEND(eAll_para_vec, Nbytes, MPI_PACKED, 0,    &
    &             FROM_WORKER, comm_leader, ierr)
             deallocate(eAll_para_vec)

         elseif (trim(worker_job_task%what_to_do) .eq. 'Clean memory' ) then
             ! clean memory before exiting, but wait - it didn't do anything
             ! useful (except telling the master so)
             if ((size_local.gt.1).and.(para_method.gt.0).and.          &
    &            (rank_local.eq.0)) then 
                 ! group leader passing the command to workers
                 if (size_local.gt.1) then !passing the command to workers
                     do des_index=1, size_local-1
                         call MPI_SEND(worker_job_package,Nbytes,         &
    &                       MPI_PACKED, des_index,  FROM_MASTER,          &    
    &                       comm_local, ierr)
                         call MPI_RECV(worker_job_package, Nbytes,        &
    &                        MPI_PACKED,des_index, FROM_WORKER,           &
    &                        comm_local, STATUS, ierr)
                         call Unpack_worker_job_task
                     end do
                 end if 
             end if
             worker_job_task%what_to_do='Cleaned Memory and Waiting'
             worker_job_task%taskid=taskid
             call Pack_worker_job_task
             call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,0,   &
    &             FROM_WORKER, comm_current, ierr)
         elseif (trim(worker_job_task%what_to_do) .eq. 'REGROUP') then
             ! calculate the time between two regroup events
             if ((size_local.gt.1).and.(para_method.gt.0).and.          &
    &            (rank_local.eq.0)) then 
                 ! group leader passing the command to workers
                 if (size_local.gt.1) then !passing the command to workers
                     do des_index=1, size_local-1
                         call MPI_SEND(worker_job_package,Nbytes,         &
    &                       MPI_PACKED, des_index,  FROM_MASTER,          &    
    &                       comm_local, ierr)
                     end do
                 end if 
             end if
             if (rank_local.eq.-1) then ! initial run
                 time_passed = -1.0
             end if
             call gather_runtime(comm_current,time_passed,time_buff)
             ! strangely, the which_per and which_pol stores the maximum
             ! number of periods and polarizations
             which_per=worker_job_task%per_index
             which_pol=worker_job_task%pol_index
             ! the following commands are always called with comm_world
             call set_group_sizes(which_per,which_pol,comm_world,group_sizes)
             call split_MPI_groups(which_per,which_pol,group_sizes)
             ! for debug
             ! write(6,*) 'rank_world= ', rank_world, 'rank_local= ', &
             !     rank_local, 'rank_leader = ', rank_leader
#if defined(CUDA) || defined(HIP)
             ! override the gpu setting after each regroup
             ! this is kind of awkward as even for a "gpu-aware" mpi, 
             ! the program will NOT be able to tell the gpu number
             ! here we have to call the cuda libs to get 
             ! the number of devices
             ! it should be noted that this only tells you the number of GPUs
             ! on current machine, i.e. could be different for different 
             ! physical nodes 
             size_gpuPtr = c_loc(size_gpu) ! kind of crude here
             ierr = cudaGetDeviceCount(size_gpuPtr)
             if ((output_level .gt. 3).and. (rank_local .eq. 0)) then
                 write(6,*) 'number of available GPU devices = ', size_gpu
             end if
             ! see if we have at least one GPU to spare in this group
             ! note this could be problematic, if the grouping is 
             ! not based on topology
             if (size_gpu*cpus_per_gpu .ge. size_node) then
                 ! everyone can get GPU acceleration
                 cpu_only_ranks = 0
             else
                 cpu_only_ranks = size_node-size_gpu*cpus_per_gpu
             end if
             ! currently use cpus_per_gpu cpus for one GPU 
             ! (hard coded in Declaritiion_MPI.f90)
             ! i.e. if we have n cpus and m gpus on a physical node
             ! the last m*cpus_per_gpu cpus will get accelerated by GPU
#else
             size_gpu = 0
             cpu_only_ranks = size_node
#endif
             if (rank_node.ge.cpu_only_ranks) then
                 ! assign the GPU(s) to the last leaders
                 device_id = mod((rank_node - cpu_only_ranks),&
     &               size_gpu)
             else
                 ! no gpu left, use CPU to calculate
                 device_id = -1
             endif
#if defined(CUDA) || defined(HIP)
             ! if we have a device to use...
             if (device_id .ge. 0) then
                 if ((output_level .gt. 3).and. (rank_local .eq. 0)) then
                     write(6,*) ' hooking up the device #', device_id
                 end if
                 ierr = cudaSetDevice(device_id);
             endif
#endif
#if defined(FG) && (defined(CUDA) || defined(HIP))
             if ((size_local .gt. 1) .and. (ncclIsInit .eq. 0)) then
               ! we have enough workers 
               ! now init the NCCL communicator
               if (rank_local .eq. 0) then
                ! leader generating the uniqueId
                 ierr = ncclGetUniqueId(uid)
                 if (ierr .ne. 0) then
                     write(6,*) 'error getting NCCL unique id on device:', ierr
                     stop
                 endif
               endif
               ! distribute the ID to all workers, using MPI
               call MPI_BCAST(uid%internal, NCCL_UNIQUE_ID_BYTES, MPI_CHAR, & 
   &                0, comm_local, ierr)
               rank_nccl = rank_local
               size_nccl = size_local
               ! for debug
               ! write(6,*) 'uid = ', uid%internal, ' @rank: ', rank_nccl
               ierr = ncclCommInitRank(comm_nccl, size_nccl, uid, rank_nccl)
               ! for debug
               ! write(6,*) 'initializing NCCL communicator... @ ', rank_nccl
               if (ierr.ne.0) then
                   write(0,'(A, I4)') 'Error initializing nccl ',ierr
                   stop
               end if 
               ncclIsInit = 1 ! set the flag
             end if
#endif
             ! now reset the timer
             time_passed = 0.0
             if (associated(time_buff)) then 
                 deallocate(time_buff) 
             endif
             if (associated(group_sizes)) then 
                 deallocate(group_sizes) 
             endif
         elseif (trim(worker_job_task%what_to_do) .eq. 'STOP' ) then
             ! clear all the temp packages and stop
             if (associated(sigma_para_vec)) then
                 deallocate(sigma_para_vec)
             end if
             if (associated(e_para_vec)) then
                 deallocate(e_para_vec)
             end if
             if (associated(eAll_para_vec)) then
                 deallocate(eAll_para_vec)
             end if
             if ((size_local.gt.1).and.(para_method.gt.0).and.          &
    &            (rank_local.eq.0)) then 
                 ! group leader passing the command to workers
                 do des_index=1, size_local-1
                     worker_job_task%what_to_do='STOP'
                     worker_job_task%taskid=taskid
                     call Pack_worker_job_task
                     call MPI_SEND(worker_job_package,Nbytes,             &
    &                   MPI_PACKED, des_index,  FROM_MASTER, comm_local, &
    &                   ierr)
                     call MPI_RECV(worker_job_package, Nbytes, MPI_PACKED &
    &                  ,des_index, FROM_WORKER, comm_local, STATUS, ierr)
                     call Unpack_worker_job_task
                 end do
                   ! all workers report finished
             end if
#if defined(FG) && (defined(CUDA) || defined(HIP))
             if ((size_local .gt. 1) .and. (ncclIsInit .ne. 0)) then
             ! for debug
             ! write(6,*) 'finalizing NCCL communicator... @ ', rank_nccl
                 ierr = ncclCommFinalize(comm_nccl) 
                 if (ierr .ne. 0 ) then
                     write(6, '(A, I2)') " nccl comm finalize error: ", ierr
                     stop
                 end if
                 ierr = ncclCommDestroy(comm_nccl) 
                 if (ierr .ne. 0 ) then
                     write(6, '(A, I2)') " nccl comm destroy error: ", ierr
                     stop
                 end if
                 ncclIsInit = 0 ! set the flag
             endif
#endif
             worker_job_task%what_to_do='Job Completed'
             worker_job_task%taskid=taskid
             call Pack_worker_job_task
             call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED, 0,  &
    &            FROM_WORKER, comm_current, ierr)
             exit
         endif
         !previous_message=trim(worker_job_task%what_to_do)
         !write(6,'(a12,a12,a30,a12)') node_info,' MPI TASK [',
         !trim(worker_job_task%what_to_do),'] successful'
         worker_job_task%what_to_do='Waiting for new message'
     end do

End Subroutine Worker_job

!******************************************************************************

subroutine create_data_vec_place_holder(d)

     implicit none
     integer Nbytes1,Nbytes2,ndata,iper,ndt,sum1,sum2
     type(dataVectorMTX_t), intent(in)           :: d
     sum1=0
     sum2=0
     do iper=1,d%nTx
         do ndt=1,d%d(iper)%ndt
             ndata=size(d%d(iper)%data(ndt)%value)
             CALL MPI_PACK_SIZE(ndata, MPI_DOUBLE_PRECISION,             &
    &             MPI_COMM_WORLD, Nbytes1,  ierr)
             CALL MPI_PACK_SIZE(1, MPI_LOGICAL, MPI_COMM_WORLD,          &
    &              Nbytes2,  ierr)
             sum1=sum1+Nbytes1
             sum2=sum2+Nbytes2
         end do
     end do
     Nbytes=((2*sum1)+(2*sum2))+1
     if (.not. associated(data_para_vec)) then
         allocate(data_para_vec(Nbytes))
     end if
         
end subroutine create_data_vec_place_holder
!****************************************************************************** 
 subroutine Pack_data_para_vec(d)
     implicit none

     type(dataVectorMTX_t), intent(in) :: d
     integer index
     integer ndata,iper,ndt
     index=1
     do iper=1,d%nTx
         do ndt=1,d%d(iper)%ndt
             ndata=size(d%d(iper)%data(ndt)%value)            
             call MPI_Pack(d%d(iper)%data(ndt)%value(1,1),ndata,          &
    &             MPI_DOUBLE_PRECISION, data_para_vec, Nbytes, index,     &
    &             MPI_COMM_WORLD, ierr)
             call MPI_Pack(d%d(iper)%data(ndt)%error(1,1),ndata,          &
    &             MPI_DOUBLE_PRECISION, data_para_vec, Nbytes, index,     &
    &             MPI_COMM_WORLD, ierr)
             call MPI_Pack(d%d(iper)%data(ndt)%errorBar,1, MPI_LOGICAL,   &
    &             data_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
             call MPI_Pack(d%d(iper)%data(ndt)%allocated,1,MPI_LOGICAL,   &
    &             data_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
       end do
     end do  

end subroutine Pack_data_para_vec
!****************************************************************************** 
subroutine UnPack_data_para_vec(d)
     implicit none
     type(dataVectorMTX_t), intent(inout)  :: d
     ! Local
     integer index
     integer ndata,iper,ndt

     index=1
     do iper=1,d%nTx
         do ndt=1,d%d(iper)%ndt
             ndata=size(d%d(iper)%data(ndt)%value)            
             call MPI_Unpack(data_para_vec, Nbytes, index,                &
    &            d%d(iper)%data(ndt)%value(1,1),ndata,                    &
    &            MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
             call MPI_Unpack(data_para_vec, Nbytes, index,                &
    &            d%d(iper)%data(ndt)%error(1,1),ndata,                    &
    &            MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
             call MPI_Unpack(data_para_vec, Nbytes, index,                &
    &            d%d(iper)%data(ndt)%errorBar  ,1, MPI_LOGICAL,           &
    &            MPI_COMM_WORLD, ierr)
             call MPI_Unpack(data_para_vec, Nbytes, index,                &
    &            d%d(iper)%data(ndt)%allocated ,1, MPI_LOGICAL,           &
    &            MPI_COMM_WORLD, ierr)
         end do
     end do  

     if (associated(data_para_vec)) then
         deallocate (data_para_vec)
     end if

end subroutine UnPack_data_para_vec

!***************************************************************************

subroutine RECV_cUserDef(cUserDef)
     implicit none
     type (userdef_control),intent(inout)   :: cUserDef
     character(20)                          :: which_proc

     call create_userdef_control_place_holder
     call MPI_RECV(userdef_control_package, Nbytes, MPI_PACKED, MASTER,   &
    &     FROM_MASTER,comm_world,STATUS, ierr)
     
     call unpack_userdef_control (cUserDef)
     if (taskid==1 ) then !first worker
         which_proc='Worker'
         call check_userdef_control_MPI (which_proc,cUserDef)
     end if
     if (associated(userdef_control_package)) then
         deallocate (userdef_control_package)
     end if
   
end subroutine RECV_cUserDef

!******************************************************************************

subroutine setGrid_MPI(newgrid)

   !  Use to set and/or update the numerical grid, that is then used
   !   all computations in this module;
   !   This is not a pointer target.
   !  Might also have to run exitSolver at this point, if we are updating
   !   the grid during an inversion; that restarts the ForwardSolver module.

   type(grid_t), intent(in)     :: newgrid

   grid = newgrid

   if (.not. grid%allocated) then
    call errStop('grid is not allocated in setGrid_MPI; exiting')
   else if ((grid%Nx <= 0) .or. (grid%Ny <= 0) .or. (grid%Nz <= 0)) then
    write(0,*) 'Grid information: Nx=',grid%Nx,' Ny=',grid%Ny,' Nz=',grid%Nz
    call errStop('grid is not set up properly in setGrid_MPI; exiting')
   end if

end subroutine setGrid_MPI

!*****************************************************************************

subroutine cleanUp_MPI()

     ! Subroutine to deallocate all memory stored in this module

     call exitSolver(e0,e,comb)
     call deall_grid(grid)

end subroutine cleanUp_MPI
  
subroutine destructor_MPI

     ! call MPI_FINALIZE (or the mpiexec might complain)
     call MPI_BARRIER(comm_world, ierr)
     call MPI_FINALIZE(ierr)

end subroutine destructor_MPI

#endif

end module Main_MPI
