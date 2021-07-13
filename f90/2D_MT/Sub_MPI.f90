
module Sub_MPI
#ifdef MPI

use math_constants
use utilities
use SolnSpace
use UserCtrl




use Declaration_MPI

Contains



!*****************************************************************************************
subroutine set_e_soln(pol_index,e)
    Integer, intent(in)                :: pol_index
    type(solnVector_t), intent(inout)  :: e

                  orginal_nPol=1


end subroutine set_e_soln
!*****************************************************************************************
subroutine count_number_of_meaasges_to_RECV(eAll1)

  type(solnVectorMTX_t), intent(in)  :: eAll1
  integer                            :: itx
            answers_to_receive=0
            do itx=1,eAll1%nTx
              answers_to_receive=answers_to_receive+(1)
			end do
			

end subroutine count_number_of_meaasges_to_RECV
!*****************************************************************************************
subroutine reset_e_soln(emsoln)
    type(solnVector_t), intent(inout)  :: emsoln



end subroutine reset_e_soln
!*****************************************************************************************
subroutine get_nPol_MPI(emsoln)

    type(solnVector_t), intent(in)  :: emsoln
             
            nPol_MPI= 1

end subroutine get_nPol_MPI





!*************************************************************************
!Packing userdef_control in a Package
  !1- Allocate a place holder
 subroutine create_userdef_control_place_holder

     implicit none
     integer Nbytes1,Nbytes2,Nbytes3,Nbytes4




       CALL MPI_PACK_SIZE(80*21, MPI_CHARACTER,        MPI_COMM_WORLD, Nbytes1,  ierr)
       CALL MPI_PACK_SIZE(3,     MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, Nbytes2,  ierr)
       CALL MPI_PACK_SIZE(1,     MPI_INTEGER,          MPI_COMM_WORLD, Nbytes3,  ierr)
        Nbytes=(Nbytes1+Nbytes2+Nbytes3)+1

         if(associated(userdef_control_package)) then
             deallocate(userdef_control_package)
         end if
             allocate(userdef_control_package(Nbytes))

 end subroutine create_userdef_control_place_holder

    !*************************************************************************
    !2- Pack ctrl into userdef_control_package
 subroutine pack_userdef_control(ctrl)
    implicit none

     	type(userdef_control), intent(in)   :: ctrl
        integer index

       index=1

        call MPI_Pack(ctrl%job,80*21, MPI_CHARACTER, userdef_control_package, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(ctrl%lambda,3, MPI_DOUBLE_PRECISION, userdef_control_package, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(ctrl%output_level,1, MPI_INTEGER, userdef_control_package, Nbytes, index, MPI_COMM_WORLD, ierr)

end subroutine pack_userdef_control

    !*************************************************************************
    !3- Unpack userdef_control_package into ctrl
 subroutine unpack_userdef_control(ctrl)
    implicit none

     	type(userdef_control), intent(in)   :: ctrl
        integer index

       index=1

   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%job,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%rFile_invCtrl,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%rFile_fwdCtrl,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%wFile_MPI,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%rFile_Grid,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%rFile_Model,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%rFile_Data,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%rFile_dModel,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%rFile_EMsoln,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%rFile_EMrhs,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%rFile_Prior,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%wFile_Grid,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%wFile_Model,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%wFile_Data,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%wFile_dModel,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%wFile_EMsoln,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%wFile_EMrhs,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%wFile_Sens,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%rFile_Cov,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%search,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%option,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)

   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%lambda,1, MPI_DOUBLE_PRECISION,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%eps,1, MPI_DOUBLE_PRECISION,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%delta,1, MPI_DOUBLE_PRECISION,MPI_COMM_WORLD, ierr)

   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%output_level,1, MPI_INTEGER,MPI_COMM_WORLD, ierr)

end subroutine unpack_userdef_control

!********************************************************************
subroutine check_userdef_control_MPI (which_proc,ctrl)

	type(userdef_control), intent(in)   :: ctrl
	character(20), intent(in)           :: which_proc

       write(6,*)trim(which_proc),' : ctrl%wFile_Sens ',trim(ctrl%wFile_Sens)
       write(6,*)trim(which_proc),' : ctrl%lambda ',(ctrl%lambda)
       write(6,*)trim(which_proc),' : ctrl%eps ',(ctrl%eps)
       write(6,*)trim(which_proc),' : ctrl%rFile_Cov ',trim(ctrl%rFile_Cov)
       write(6,*)trim(which_proc),' : ctrl%search ',trim(ctrl%search)
       write(6,*)trim(which_proc),' : ctrl%output_level ',ctrl%output_level
       write(6,*)trim(which_proc),' : ctrl%rFile_fwdCtrl ',trim(ctrl%rFile_fwdCtrl)
       write(6,*)trim(which_proc),' : ctrl%rFile_invCtrl ',trim(ctrl%rFile_invCtrl)


end subroutine check_userdef_control_MPI
!********************************************************************
subroutine create_e_param_place_holder(e)

     implicit none
     type(solnVector_t), intent(in)	:: e
     integer                        ::v_size,Nbytes1,Nbytes2






         v_size=size(e%vec%v)
         CALL MPI_PACK_SIZE(v_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes1,  ierr)
		 CALL MPI_PACK_SIZE(1, MPI_INTEGER, MPI_COMM_WORLD, Nbytes2,  ierr)

         Nbytes=((Nbytes1+Nbytes2))+1

         if(associated(e_para_vec)) then
             deallocate(e_para_vec)
         end if
             allocate(e_para_vec(Nbytes))


 end subroutine create_e_param_place_holder
 !********************************************************************
 subroutine Pack_e_para_vec(e)
    implicit none

     type(solnVector_t), intent(in)	:: e
     integer index,v_size



       v_size=size(e%vec%v)

       index=1

        call MPI_Pack(e%vec%v(1,1),v_size, MPI_DOUBLE_COMPLEX, e_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
		call MPI_Pack(e%tx,1,             MPI_INTEGER,        e_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)


end subroutine Pack_e_para_vec
!********************************************************************
subroutine Unpack_e_para_vec(e)
    implicit none

     type(solnVector_t), intent(inout)	:: e

     integer index,v_size


       v_size=size(e%vec%v)
       index=1

        call MPI_Unpack(e_para_vec, Nbytes, index, e%vec%v(1,1),v_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
		call MPI_Unpack(e_para_vec, Nbytes, index, e%tx,1, MPI_INTEGER,MPI_COMM_WORLD, ierr)

end subroutine Unpack_e_para_vec

!********************************************************************
subroutine create_eAll_param_place_holder(e)

     implicit none
     type(solnVector_t), intent(in)	:: e
     integer v_size,Nbytes1,Nbytes2


         v_size=size(e%vec%v)
         CALL MPI_PACK_SIZE(v_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes1,  ierr)
		 CALL MPI_PACK_SIZE(1, MPI_INTEGER, MPI_COMM_WORLD, Nbytes2,  ierr)

         Nbytes=((Nbytes1+Nbytes2))+1



         if(associated(eAll_para_vec)) then
             deallocate(eAll_para_vec)
         end if
             allocate(eAll_para_vec(Nbytes))



 end subroutine create_eAll_param_place_holder



!********************************************************************
subroutine pack_eAll_para_vec(e)
    implicit none

     type(solnVector_t), intent(in)	:: e
     integer index,v_size



       v_size=size(e%vec%v)

       index=1

        call MPI_Pack(e%vec%v(1,1),v_size, MPI_DOUBLE_COMPLEX, eAll_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(e%tx,1,             MPI_INTEGER,        eAll_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)

end subroutine pack_eAll_para_vec


!********************************************************************
subroutine Unpack_eAll_para_vec(e)
    implicit none

     type(solnVector_t), intent(in)	:: e

     integer index,v_size


       v_size=size(e%vec%v)
       index=1

        call MPI_Unpack(eAll_para_vec, Nbytes, index, e%vec%v(1,1),v_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
		call MPI_Unpack(eAll_para_vec, Nbytes, index, e%tx,1, MPI_INTEGER,MPI_COMM_WORLD, ierr)



end subroutine Unpack_eAll_para_vec


#endif

end module Sub_MPI
