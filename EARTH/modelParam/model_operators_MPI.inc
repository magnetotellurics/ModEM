
#ifdef MPI
  !Create a vector suitable for MPI messages from modelParam_t

!#####################################################################################
subroutine create_model_param_place_holder(sigma_mpi)
    implicit none
    include 'mpif.h'

     type (modelParam_t), intent(in)            :: sigma_mpi
     integer                                    ::v_size,Nbytes1

         v_size = count(sigma_mpi%c%exists)
         CALL MPI_PACK_SIZE(v_size, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, Nbytes1,  ierr)

         Nbytes=((Nbytes1))+1

         if(associated(sigma_para_vec)) then
             deallocate(sigma_para_vec)
         end if
         allocate(sigma_para_vec(Nbytes))

end subroutine create_model_param_place_holder
!#####################################################################################

subroutine Pack_model_para_values(sigma_mpi)
    implicit none
    include 'mpif.h'

    type (modelParam_t), intent(in)             :: sigma_mpi
    real(8), allocatable, dimension(:)          :: sigma_values
    integer                                     :: index,istat,v_size

       v_size=count(sigma_mpi%c%exists)
       allocate(sigma_values(v_size),STAT=istat)
       call getParamValues_modelParam(sigma_mpi, sigma_values)
       index=1
       call MPI_Pack(sigma_values(1),v_size, MPI_DOUBLE_PRECISION, sigma_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
       deallocate(sigma_values,STAT=istat)

end subroutine Pack_model_para_values
!#####################################################################################

subroutine Unpack_model_para_values(sigma_mpi)
    implicit none
    include 'mpif.h'

    type (modelParam_t), intent(inout)          :: sigma_mpi
    real(8), allocatable, dimension(:)          :: sigma_values
    integer                                     :: index,istat,v_size

       v_size=count(sigma_mpi%c%exists)
       allocate(sigma_values(v_size),STAT=istat)
       index=1
       call MPI_Unpack(sigma_para_vec, Nbytes, index,sigma_values(1),v_size, MPI_DOUBLE_PRECISION,MPI_COMM_WORLD, ierr)
       call fillParamValues_modelParam(sigma_mpi, sigma_values)
       deallocate(sigma_values,STAT=istat)

end subroutine unPack_model_para_values
!#####################################################################################


#endif
