program rmatest

  use iso_c_binding, only: c_ptr, c_f_pointer

  use mpi

  implicit none

! Set the size of the road

  integer, parameter :: nlocal = 5
  integer :: i, n
  integer, dimension(MPI_STATUS_SIZE) :: status

  integer, pointer, dimension(:) :: rma

  integer :: comm, nodecomm, nodewin
  integer :: ierr, size, rank, nodesize, noderank, nodestringlen
  integer(MPI_ADDRESS_KIND) :: winsize
  integer :: intsize, disp_unit
  character*(MPI_MAX_PROCESSOR_NAME) :: nodename
  type(c_ptr) :: baseptr

  comm = MPI_COMM_WORLD

  call MPI_Init(ierr)

  call MPI_Comm_size(comm, size, ierr)
  call MPI_Comm_rank(comm, rank, ierr)

  ! Create node-local communicator

  call MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, rank, &
                           MPI_INFO_NULL, nodecomm, ierr)

  ! Check it all went as expected

  call MPI_Get_processor_name(nodename, nodestringlen, ierr)
  call MPI_Comm_size(nodecomm, nodesize, ierr)
  call MPI_Comm_rank(nodecomm, noderank, ierr)

  n = nlocal*nodesize

  if (rank == 0) then

     write(*,*) "Running on ", size, " processes with n = ", n

  end if

  write(*,*) "Rank ", rank," in COMM_WORLD is rank ", noderank, &
             " in nodecomm on node ", nodename(1:nodestringlen)

  call MPI_Type_size(MPI_INTEGER, intsize, ierr)

  winsize = nlocal*intsize

  ! displacements counted in units of integers

  disp_unit = intsize

  call MPI_Win_allocate_shared(winsize, disp_unit, &
       MPI_INFO_NULL, nodecomm, baseptr, nodewin, ierr)

  ! coerce baseptr to a Fortran array: global on rank 0, local on others

  if (noderank == 0) then

     call c_f_pointer(baseptr, rma, [n])

  else

     call c_f_pointer(baseptr, rma, [nlocal])

  end if

  ! Set the local arrays

  rma(1:nlocal) = 0

  ! Set values on noderank 0

  call MPI_Win_fence(0, nodewin, ierr)

  if (rank == 0) then
     do i = 1, n
        rma(i) = i
     end do
  end if

  call MPI_Win_fence(0, nodewin, ierr)

  ! Print the values  

  write(*,*) "rank, noderank, arr: ", rank, noderank, (rma(i), i=1,nlocal)

  call MPI_Finalize(ierr)

end program rmatest