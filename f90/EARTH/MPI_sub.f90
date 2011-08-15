
module MPI_sub
#ifdef MPI

use math_constants
use utilities
use SolnSpace
use ioTypes
use MPI_declaration

implicit none




Contains

!*************************************************************************
!Packing userdef_control in a Package
  !1- Allocate a place holder
 subroutine create_userdef_control_place_holder

     implicit none
     integer Nbytes1,Nbytes2,Nbytes3,Nbytes4

       CALL MPI_PACK_SIZE(80*(21+5), MPI_CHARACTER,        MPI_COMM_WORLD, Nbytes1,  ierr)
       CALL MPI_PACK_SIZE(2,     MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, Nbytes2,  ierr)
        Nbytes=(Nbytes1+Nbytes2)+1

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

        call MPI_Pack(ctrl%paramname,80*(21+5), MPI_CHARACTER, userdef_control_package, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(ctrl%damping,2, MPI_DOUBLE_PRECISION, userdef_control_package, Nbytes, index, MPI_COMM_WORLD, ierr)

end subroutine pack_userdef_control

    !*************************************************************************
    !3- Unpack userdef_control_package into ctrl
 subroutine unpack_userdef_control(ctrl)
    implicit none

        type(userdef_control), intent(in)   :: ctrl
        integer index

       index=1

   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%paramname,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%modelname,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%verbose,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%calculate,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%secondary_field,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%fn_thinsheet,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%fn_grid,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%fn_param0,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%fn_param,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%fn_period,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%fn_extsource,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%fn_intsource,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%fn_fwdctrl,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%fn_invctrl,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%fn_slices,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%fn_coords,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%fn_func,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%fn_cdata,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%fn_ddata,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%fn_hdata,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%fn_href,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)

   ! obsolete
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%fn_field,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%fn_precond,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%fn_misfit,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%fn_gradient,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%fn_point,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%damping,1, MPI_DOUBLE_PRECISION,MPI_COMM_WORLD, ierr)
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%step_size,1, MPI_DOUBLE_PRECISION,MPI_COMM_WORLD, ierr)

end subroutine unpack_userdef_control

!********************************************************************
subroutine check_userdef_control_MPI (which_proc,ctrl)

    type(userdef_control), intent(in)   :: ctrl
    character(20), intent(in)           :: which_proc

       write(6,*) '[',trim(which_proc),'] Parametrization: ', trim(ctrl%paramname)
       write(6,*) '[',trim(which_proc),'] Output filename: ', trim(ctrl%modelname)
       write(6,*) '[',trim(which_proc),'] Level of output: ', trim(ctrl%verbose)
       write(6,*) '[',trim(which_proc),'] What to compute: ', trim(ctrl%calculate)
       write(6,*) '[',trim(which_proc),'] Secondary field: ', trim(ctrl%secondary_field)
       write(6,*) '[',trim(which_proc),'] Thinsheet file : ', trim(ctrl%fn_thinsheet)
       write(6,*) '[',trim(which_proc),'] Numerical grid : ', trim(ctrl%fn_grid)
       write(6,*) '[',trim(which_proc),'] Base model info: ', trim(ctrl%fn_param0)
       write(6,*) '[',trim(which_proc),'] Parameters info: ', trim(ctrl%fn_param)
       write(6,*) '[',trim(which_proc),'] Period filename: ', trim(ctrl%fn_period)
       write(6,*) '[',trim(which_proc),'] External source: ', trim(ctrl%fn_extsource)
       write(6,*) '[',trim(which_proc),'] Interior source: ', trim(ctrl%fn_intsource)
       write(6,*) '[',trim(which_proc),'] Forward control: ', trim(ctrl%fn_fwdctrl)
       write(6,*) '[',trim(which_proc),'] Inverse control: ', trim(ctrl%fn_invctrl)
       write(6,*) '[',trim(which_proc),'] Radii to output: ', trim(ctrl%fn_slices)
       write(6,*) '[',trim(which_proc),'] Obs coordinates: ', trim(ctrl%fn_coords)
       write(6,*) '[',trim(which_proc),'] Functional info: ', trim(ctrl%fn_func)
       write(6,*) '[',trim(which_proc),'] C response data: ', trim(ctrl%fn_cdata)
       write(6,*) '[',trim(which_proc),'] D response data: ', trim(ctrl%fn_ddata)
       write(6,*) '[',trim(which_proc),'] Field data file: ', trim(ctrl%fn_hdata)
       write(6,*) '[',trim(which_proc),'] Field reference: ', trim(ctrl%fn_href)

end subroutine check_userdef_control_MPI

!*****************************************************************************************
subroutine set_e_soln(pol_index,e)
    Integer, intent(in)                :: pol_index
    type(solnVector_t), intent(inout)  :: e

            ! Empty stub for global code
            ! e%nPol=1
            ! e%Pol_index(1)=pol_index

end subroutine set_e_soln
!*****************************************************************************************
subroutine reset_e_soln(emsoln)
    type(solnVector_t), intent(inout)  :: emsoln

 ! Empty stub for global code

end subroutine reset_e_soln
!*****************************************************************************************

subroutine get_nPol_MPI(eAll)

   type(solnVectorMTX_t), intent(in)    :: eAll


            nPol_MPI= 1

end subroutine get_nPol_MPI

!********************************************************************
subroutine create_e_param_place_holder(e)

     implicit none
     type(solnVector_t), intent(in) :: e
     integer                        :: Ex_size,Ey_size,Ez_size,Nbytes1,Nbytes2,Nbytes3

       Ex_size=size(e%vec%x)
       CALL MPI_PACK_SIZE(Ex_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes1,  ierr)
       Ey_size=size(e%vec%y)
       CALL MPI_PACK_SIZE(Ey_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes2,  ierr)
       Ez_size=size(e%vec%z)
       CALL MPI_PACK_SIZE(Ez_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes3,  ierr)

         Nbytes=((Nbytes1+Nbytes2+Nbytes3))+1

         if(associated(e_para_vec)) then
             deallocate(e_para_vec)
         end if
         allocate(e_para_vec(Nbytes))

 end subroutine create_e_param_place_holder
!********************************************************************
 subroutine Pack_e_para_vec(e)
    implicit none

     type(solnVector_t), intent(in) :: e
     integer index,Ex_size,Ey_size,Ez_size



       Ex_size=size(e%vec%x)
       Ey_size=size(e%vec%y)
       Ez_size=size(e%vec%z)
       index=1

        call MPI_Pack(e%vec%x(1,1,1),Ex_size, MPI_DOUBLE_COMPLEX, e_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(e%vec%y(1,1,1),Ey_size, MPI_DOUBLE_COMPLEX, e_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(e%vec%z(1,1,1),Ez_size, MPI_DOUBLE_COMPLEX, e_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)



end subroutine Pack_e_para_vec
!********************************************************************
subroutine Unpack_e_para_vec(e)
    implicit none

     type(solnVector_t), intent(inout)  :: e

     integer index,Ex_size,Ey_size,Ez_size


       Ex_size=size(e%vec%x)
       Ey_size=size(e%vec%y)
       Ez_size=size(e%vec%z)
       index=1

        call MPI_Unpack(e_para_vec, Nbytes, index, e%vec%x(1,1,1),Ex_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(e_para_vec, Nbytes, index, e%vec%y(1,1,1),Ey_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(e_para_vec, Nbytes, index, e%vec%z(1,1,1),Ez_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)





end subroutine Unpack_e_para_vec

!********************************************************************
subroutine create_eAll_param_place_holder(e)

     implicit none
     type(solnVector_t), intent(in)  :: e
     integer Ex_size,Ey_size,Ez_size,Nbytes1,Nbytes2,Nbytes3



       Ex_size=size(e%vec%x)
       CALL MPI_PACK_SIZE(Ex_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes1,  ierr)
       Ey_size=size(e%vec%y)
       CALL MPI_PACK_SIZE(Ey_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes2,  ierr)
       Ez_size=size(e%vec%z)
       CALL MPI_PACK_SIZE(Ez_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes3,  ierr)

         Nbytes=(2*(Nbytes1+Nbytes2+Nbytes3))+1  ! Multiple by 2 for both polarizations

         if(associated(eAll_para_vec)) then
             deallocate(eAll_para_vec)
         end if
         allocate(eAll_para_vec(Nbytes))



 end subroutine create_eAll_param_place_holder



!********************************************************************
subroutine pack_eAll_para_vec(e)
    implicit none

     type(solnVector_t), intent(in)  :: e
     integer index,Ex_size,Ey_size,Ez_size



       Ex_size=size(e%vec%x)
       Ey_size=size(e%vec%y)
       Ez_size=size(e%vec%z)
       index=1


        call MPI_Pack(e%vec%x(1,1,1),Ex_size, MPI_DOUBLE_COMPLEX, eAll_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(e%vec%y(1,1,1),Ey_size, MPI_DOUBLE_COMPLEX, eAll_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(e%vec%z(1,1,1),Ez_size, MPI_DOUBLE_COMPLEX, eAll_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)

end subroutine pack_eAll_para_vec


!********************************************************************
subroutine Unpack_eAll_para_vec(e)
    implicit none

     type(solnVector_t), intent(inout)   :: e

     integer index,Ex_size,Ey_size,Ez_size


       Ex_size=size(e%vec%x)
       Ey_size=size(e%vec%y)
       Ez_size=size(e%vec%z)
       index=1



        call MPI_Unpack(eAll_para_vec, Nbytes, index, e%vec%x(1,1,1),Ex_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(eAll_para_vec, Nbytes, index, e%vec%y(1,1,1),Ey_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(eAll_para_vec, Nbytes, index, e%vec%z(1,1,1),Ez_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)



end subroutine Unpack_eAll_para_vec


#endif

end module MPI_sub
