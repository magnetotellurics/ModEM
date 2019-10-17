
#ifdef MPI

!#####################################################################################
subroutine create_model_param_place_holder(sigma_mpi)
#ifdef MPIMOD
    use mpi
    implicit none
#else
    implicit none
    include 'mpif.h'
#endif

     type (modelParam_t), intent(in)   	        :: sigma_mpi
     integer                                    ::v_size,Nbytes1                  
         v_size=size(sigma_mpi%v)
		 
		 Nbytes=0
		 CALL MPI_PACK_SIZE(2, MPI_INTEGER, MPI_COMM_WORLD, Nbytes1,  ierr)
		 Nbytes=Nbytes+Nbytes1
		 
         CALL MPI_PACK_SIZE(v_size, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, Nbytes1,  ierr)
		 Nbytes=Nbytes+Nbytes1
		 
		 CALL MPI_PACK_SIZE(1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, Nbytes1,  ierr)
		 Nbytes=Nbytes+Nbytes1
		 
		 CALL MPI_PACK_SIZE(3, MPI_LOGICAL, MPI_COMM_WORLD, Nbytes1,  ierr)
         Nbytes=Nbytes+Nbytes1
		 
		 CALL MPI_PACK_SIZE(80, MPI_CHARACTER, MPI_COMM_WORLD, Nbytes1,  ierr)
         Nbytes=Nbytes+Nbytes1		 
		 
		 
         Nbytes=Nbytes+1  

         if(associated(sigma_para_vec)) then
             deallocate(sigma_para_vec)
         end if
             allocate(sigma_para_vec(Nbytes))

end subroutine create_model_param_place_holder
!#####################################################################################

subroutine Pack_model_para_values(sigma_mpi)
#ifdef MPIMOD
    use mpi
    implicit none
#else
    implicit none
    include 'mpif.h'
#endif
    
    type (modelParam_t), intent(in)   	        :: sigma_mpi
    integer                                     :: index,v_size

       v_size=size(sigma_mpi%v)
       index=1
	   call MPI_Pack(sigma_mpi%ny,1, MPI_INTEGER, sigma_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
	   call MPI_Pack(sigma_mpi%NzEarth,1, MPI_INTEGER, sigma_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
       call MPI_Pack(sigma_mpi%v(1,1),v_size, MPI_DOUBLE_PRECISION, sigma_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
	   call MPI_Pack(sigma_mpi%AirCond,1, MPI_DOUBLE_PRECISION, sigma_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
	   call MPI_Pack(sigma_mpi%allocated,1, MPI_LOGICAL, sigma_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
	   call MPI_Pack(sigma_mpi%zeroValued,1, MPI_LOGICAL, sigma_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
	   call MPI_Pack(sigma_mpi%temporary,1, MPI_LOGICAL, sigma_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
	   call MPI_Pack(sigma_mpi%paramType,80, MPI_CHARACTER, sigma_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
        
end subroutine Pack_model_para_values
!#####################################################################################

subroutine Unpack_model_para_values(sigma_mpi)
#ifdef MPIMOD
    use mpi
    implicit none
#else
    implicit none
    include 'mpif.h'
#endif
    
    type (modelParam_t), intent(inout)   	    :: sigma_mpi
    integer                                     :: index,v_size

       v_size=size(sigma_mpi%v)
       index=1
	   call MPI_Unpack(sigma_para_vec, Nbytes, index,sigma_mpi%ny,1, MPI_INTEGER,MPI_COMM_WORLD, ierr)
	   call MPI_Unpack(sigma_para_vec, Nbytes, index,sigma_mpi%NzEarth,1, MPI_INTEGER,MPI_COMM_WORLD, ierr)
       call MPI_Unpack(sigma_para_vec, Nbytes, index,sigma_mpi%v(1,1),v_size, MPI_DOUBLE_PRECISION,MPI_COMM_WORLD, ierr)
	   call MPI_Unpack(sigma_para_vec, Nbytes, index,sigma_mpi%AirCond,1, MPI_DOUBLE_PRECISION,MPI_COMM_WORLD, ierr)
	   call MPI_Unpack(sigma_para_vec, Nbytes, index,sigma_mpi%allocated,1, MPI_LOGICAL,MPI_COMM_WORLD, ierr)
	   call MPI_Unpack(sigma_para_vec, Nbytes, index,sigma_mpi%zeroValued,1, MPI_LOGICAL,MPI_COMM_WORLD, ierr)
	   call MPI_Unpack(sigma_para_vec, Nbytes, index,sigma_mpi%temporary,1, MPI_LOGICAL,MPI_COMM_WORLD, ierr)
	   call MPI_Unpack(sigma_para_vec, Nbytes, index,sigma_mpi%paramType,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
        
end subroutine unPack_model_para_values
!#####################################################################################

#endif