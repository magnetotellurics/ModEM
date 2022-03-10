!*************
!
! MPI declarations
!
! Last modified at 03/2022 by Paulo Werdt
!
!*************
!
module DeclarationMPI
	!
	use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
    !
	!use mpi
	use mpi_f08
	!
	type( MPI_Win )  :: nodewin
	type( MPI_Comm ) :: main_comm, shared_comm
	!
	character, pointer, dimension(:) :: job_package
	!
	integer							 ::	nbytes, &
										disp_unit = 0, node_size, node_rank, nodestringlen, &
										mpi_rank, mpi_size, ierr
	!
    integer( MPI_ADDRESS_KIND ) :: winsize
    !
    character*( MPI_MAX_PROCESSOR_NAME ) :: nodename
    type( c_ptr ) :: baseptr
	!
	! PROGRAM GLOBAL VARIABLES
	integer			:: tag = 666, master_id = 0
	!
	character*15	:: job_master = "OKAY_BOSS", job_worker, job_finish = "STOP_JOBS", job_ok = "OKAY_BOSS"
	!
	! STRUCT job_info
	type :: struct_job_info
		SEQUENCE
		character*15	:: name = "FORWARD"
		integer       	:: id_rank
	end type struct_job_info
	!
	type( struct_job_info ), save :: job_info
	!
	contains
	!
	! ALLOCATE job_package
	subroutine allocateJobPackage
		!
		integer nbytes1, nbytes2
		!
		call MPI_PACK_SIZE( 15, MPI_CHARACTER, shared_comm, nbytes1, ierr )
		call MPI_PACK_SIZE( 1, MPI_INTEGER, shared_comm, nbytes2, ierr )
		!
		nbytes = ( nbytes1 + nbytes2 ) + 1
		!
		if( .not. associated( job_package ) ) then
			allocate( job_package( nbytes ) )
		endif
		!
	end subroutine allocateJobPackage
	!
	! PACK job_info STRUCT TO job_package
	subroutine packJobTask
		!
		integer :: index
		!
		index = 1
		!
		call MPI_PACK( job_info%name, 15, MPI_CHARACTER, job_package, nbytes, index, shared_comm, ierr )
		call MPI_PACK( job_info%id_rank, 1, MPI_INTEGER, job_package, nbytes, index, shared_comm, ierr )
		!
	end subroutine packJobTask
	!
	! UNPACK job_package TO job_info STRUCT
	subroutine unpackJobTask
		!
		integer :: index
		!
		index = 1
		!
		call MPI_UNPACK( job_package, nbytes, index, job_info%name, 15, MPI_CHARACTER, shared_comm, ierr )
		call MPI_UNPACK( job_package, nbytes, index, job_info%id_rank , 1, MPI_INTEGER, shared_comm, ierr )
		!
	end subroutine unpackJobTask
	!
	! RECEIVE job_info FROM target_id
	subroutine receiveFrom( target_id )
		!
		integer, intent(in)	:: target_id
		!
		!write( *, * ) "<<<< ", mpi_rank, " RECV: ", job_info%name, " FROM: ", target_id
		!
		call allocateJobPackage
		call MPI_RECV( job_package, nbytes, MPI_PACKED, target_id, tag, shared_comm, MPI_STATUS_IGNORE, ierr)
		call unpackJobTask
		!
	end subroutine receiveFrom
	!
	! SEND job_info FROM target_id
	subroutine sendTo( target_id )
		!
		integer, intent(in)	:: target_id
		!
		!write( *, * ) ">>>> ", mpi_rank, " SEND: ", job_info%name, " TO: ", target_id
		!
		call allocateJobPackage
		call packJobTask
		call MPI_SEND( job_package, nbytes, MPI_PACKED, target_id, tag, shared_comm, ierr )
		!
	end subroutine sendTo
	!
end module DeclarationMPI