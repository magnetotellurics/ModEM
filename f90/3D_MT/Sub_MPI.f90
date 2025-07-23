
module Sub_MPI
#ifdef MPI

use math_constants
use utilities
use SolnSpace
use UserCtrl
use ForwardSolver

use Declaration_MPI

Contains

!##############################################   Start Nested_parameters #########################################################
subroutine Master_job_Distribute_nTx_nPol(nTx_nPol)
   integer, intent(in) :: nTx_nPol
   integer             :: j

           call MPI_BCAST(nTx_nPol,1, MPI_INTEGER,0, MPI_COMM_WORLD,ierr) 

       do j=1, nTx_nPol    
            call create_boundary_param_place_holder(BC_from_file(j))
            call Pack_boundary_para_vec(BC_from_file(j))
            call MPI_BCAST(e_para_vec,Nbytes, MPI_PACKED,0, MPI_COMM_WORLD,ierr)
			deallocate( e_para_vec )
       end do     
           
   
end subroutine Master_job_Distribute_nTx_nPol   
!*****************************************************************************************
subroutine RECV_nTx_nPol
        
       call MPI_BCAST(nTx_nPol,1, MPI_INTEGER,0, MPI_COMM_WORLD,ierr)           
            
end subroutine RECV_nTx_nPol
!***************************************************************************************** 
subroutine RECV_BC_form_Master
integer          :: j

        do j=1, nTx_nPol    
            call create_boundary_param_place_holder(BC_from_file(j))
            call MPI_BCAST(e_para_vec,Nbytes, MPI_PACKED,0, MPI_COMM_WORLD,ierr) 
            call unPack_boundary_para_vec(BC_from_file(j))   
        end do  
        
end subroutine RECV_BC_form_Master      
!***************************************************************************************** 
!##############################################   End Nested_parameters #########################################################
subroutine set_e_soln(pol_index,emsoln)
    Integer, intent(in)                :: pol_index
    type(solnVector_t), intent(inout)  :: emsoln
	
            orginal_nPol=emsoln%nPol
            emsoln%nPol=1
		    emsoln%Pol_index(1)=pol_index

end subroutine set_e_soln
!*****************************************************************************************

!*****************************************************************************************
subroutine reset_e_soln(emsoln)
    type(solnVector_t), intent(inout)  :: emsoln

            emsoln%nPol=orginal_nPol

end subroutine reset_e_soln
!*****************************************************************************************

subroutine get_nPol_MPI(emsoln)

  type(solnVector_t), intent(in)  :: emsoln


            nPol_MPI= emsoln%nPol

end subroutine get_nPol_MPI
!*****************************************************************************************

subroutine count_number_of_meaasges_to_RECV(eAll1)

  type(solnVectorMTX_t), intent(in)  :: eAll1
  integer                            :: itx
            answers_to_receive=0
            do itx=1,eAll1%nTx
              answers_to_receive=answers_to_receive+(eAll1%solns(itx)%nPol)
			end do
			

end subroutine count_number_of_meaasges_to_RECV




!*************************************************************************
!Packing userdef_control in a Package
  !1- Allocate a place holder
 subroutine create_userdef_control_place_holder

	implicit none
	integer Nbytes1,Nbytes2,Nbytes3,Nbytes4
	!
	CALL MPI_PACK_SIZE(80*22, MPI_CHARACTER,        MPI_COMM_WORLD, Nbytes1,  ierr)
	CALL MPI_PACK_SIZE(3,     MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, Nbytes2,  ierr)
	CALL MPI_PACK_SIZE(1,     MPI_INTEGER,          MPI_COMM_WORLD, Nbytes3,  ierr)
	Nbytes=(Nbytes1+Nbytes2+Nbytes3)+1

	if( associated( userdef_control_package ) ) deallocate( userdef_control_package )
	allocate( userdef_control_package( Nbytes ) )
	!
 end subroutine create_userdef_control_place_holder

    !*************************************************************************
    !2- Pack ctrl into userdef_control_package
 subroutine pack_userdef_control(ctrl)
    implicit none

     	type(userdef_control), intent(in)   :: ctrl
        integer index

		index=1

		call MPI_Pack(ctrl%job,80*22, MPI_CHARACTER, userdef_control_package, Nbytes, index, MPI_COMM_WORLD, ierr)
		!
        call MPI_Pack(ctrl%lambda,3, MPI_DOUBLE_PRECISION, userdef_control_package, Nbytes, index, MPI_COMM_WORLD, ierr)
		!
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
   call MPI_Unpack(userdef_control_package, Nbytes, index, ctrl%rFile_Config,80, MPI_CHARACTER,MPI_COMM_WORLD, ierr)
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
   
   deallocate( userdef_control_package )

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
	   write(6,*)trim(which_proc),' : ctrl%rFile_Config ',trim(ctrl%rFile_Config)


end subroutine check_userdef_control_MPI
!********************************************************************
subroutine create_boundary_param_place_holder(cb)
  type (cboundary), intent(in)     :: cb
  integer                        :: size1,Nbytes1,sum
  
     !complex (kind=prec), pointer, dimension(:,:)    :: xYMax, zYMax
     !complex (kind=prec), pointer, dimension(:,:)    :: xYMin, zYMin
     !complex (kind=prec), pointer, dimension(:,:)    :: yXMax, zXMax
     !complex (kind=prec), pointer, dimension(:,:)    :: yXMin, zXMin
     !complex (kind=prec), pointer, dimension(:,:)    :: xZMin, yZMin
     !complex (kind=prec), pointer, dimension(:,:)    :: xZMax, yZMax 
     
  
     Nbytes=0
     size1=size(cb%xYMax)
     CALL MPI_PACK_SIZE(size1, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes1,  ierr)
     Nbytes=Nbytes+Nbytes1
     
     size1=size(cb%zYMax)
     CALL MPI_PACK_SIZE(size1, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes1,  ierr)
     Nbytes=Nbytes+Nbytes1
     
     size1=size(cb%xYMin)
     CALL MPI_PACK_SIZE(size1, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes1,  ierr)     
     Nbytes=Nbytes+Nbytes1
     
     size1=size(cb%zYMin)
     CALL MPI_PACK_SIZE(size1, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes1,  ierr)    
     Nbytes=Nbytes+Nbytes1
     
     size1=size(cb%yXMax)
     CALL MPI_PACK_SIZE(size1, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes1,  ierr) 
     Nbytes=Nbytes+Nbytes1
     
     size1=size(cb%zXMax)
     CALL MPI_PACK_SIZE(size1, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes1,  ierr) 
     Nbytes=Nbytes+Nbytes1
          
     size1=size(cb%yXMin)
     CALL MPI_PACK_SIZE(size1, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes1,  ierr) 
     Nbytes=Nbytes+Nbytes1
     
     size1=size(cb%zXMin)
     CALL MPI_PACK_SIZE(size1, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes1,  ierr) 
     Nbytes=Nbytes+Nbytes1
     
     size1=size(cb%xZMin)
     CALL MPI_PACK_SIZE(size1, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes1,  ierr)      
     Nbytes=Nbytes+Nbytes1
     
     size1=size(cb%yZMin)
     CALL MPI_PACK_SIZE(size1, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes1,  ierr)
     Nbytes=Nbytes+Nbytes1     
     
     size1=size(cb%xZMax)
     CALL MPI_PACK_SIZE(size1, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes1,  ierr)
     Nbytes=Nbytes+Nbytes1
      
     size1=size(cb%yZMax)
     CALL MPI_PACK_SIZE(size1, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes1,  ierr)    
     Nbytes=Nbytes+Nbytes1+1

	 if( associated( e_para_vec ) ) deallocate( e_para_vec )
	 allocate( e_para_vec( Nbytes ) )
	 !
    
end subroutine create_boundary_param_place_holder
 !********************************************************************
 subroutine Pack_boundary_para_vec(cb)
   type (cboundary), intent(in)     :: cb
   integer                          :: size1,index
   
     !complex (kind=prec), pointer, dimension(:,:)    :: xYMax, zYMax
     !complex (kind=prec), pointer, dimension(:,:)    :: xYMin, zYMin
     !complex (kind=prec), pointer, dimension(:,:)    :: yXMax, zXMax
     !complex (kind=prec), pointer, dimension(:,:)    :: yXMin, zXMin
     !complex (kind=prec), pointer, dimension(:,:)    :: xZMin, yZMin
     !complex (kind=prec), pointer, dimension(:,:)    :: xZMax, yZMax 
   
    index=1
    size1=size(cb%xYMax) ! xYMax
    call MPI_Pack(cb%xYMax(1,1),size1, MPI_DOUBLE_COMPLEX, e_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
    
    size1=size(cb%zYMax) ! zYMax
    call MPI_Pack(cb%zYMax(1,1),size1, MPI_DOUBLE_COMPLEX, e_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
    
    size1=size(cb%xYMin) ! xYMin
    call MPI_Pack(cb%xYMin(1,1),size1, MPI_DOUBLE_COMPLEX, e_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
    
    size1=size(cb%zYMin) ! zYMin
    call MPI_Pack(cb%zYMin(1,1),size1, MPI_DOUBLE_COMPLEX, e_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
    
    size1=size(cb%yXMax) ! yXMax
    call MPI_Pack(cb%yXMax(1,1),size1, MPI_DOUBLE_COMPLEX, e_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)    
  
    size1=size(cb%zXMax) ! zXMax
    call MPI_Pack(cb%zXMax(1,1),size1, MPI_DOUBLE_COMPLEX, e_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
    
    size1=size(cb%yXMin) ! yXMin
    call MPI_Pack(cb%yXMin(1,1),size1, MPI_DOUBLE_COMPLEX, e_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
    
    size1=size(cb%zXMin) ! zXMin
    call MPI_Pack(cb%zXMin(1,1),size1, MPI_DOUBLE_COMPLEX, e_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
    
    size1=size(cb%xZMin) ! xZMin
    call MPI_Pack(cb%xZMin(1,1),size1, MPI_DOUBLE_COMPLEX, e_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
    
    size1=size(cb%yZMin) ! yZMin
    call MPI_Pack(cb%yZMin(1,1),size1, MPI_DOUBLE_COMPLEX, e_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)      
    
    size1=size(cb%xZMax) ! xZMax
    call MPI_Pack(cb%xZMax(1,1),size1, MPI_DOUBLE_COMPLEX, e_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
    
    size1=size(cb%yZMax) !yZMax
    call MPI_Pack(cb%yZMax(1,1),size1, MPI_DOUBLE_COMPLEX, e_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)     
    
    
 
 end subroutine Pack_boundary_para_vec
!********************************************************************
 subroutine UnPack_boundary_para_vec(cb)
   type (cboundary), intent(inout)     :: cb
   integer                          :: size1,index
 
     index=1
     size1=size(cb%xYMax) ! xYMax
     call MPI_Unpack(e_para_vec, Nbytes, index, cb%xYMax(1,1),size1, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
 
     size1=size(cb%zYMax) ! zYMax
     call MPI_Unpack(e_para_vec, Nbytes, index, cb%zYMax(1,1),size1, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
    
     size1=size(cb%xYMin) ! xYMin
     call MPI_Unpack(e_para_vec, Nbytes, index, cb%xYMin(1,1),size1, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
    
     size1=size(cb%zYMin) ! zYMin
     call MPI_Unpack(e_para_vec, Nbytes, index, cb%zYMin(1,1),size1, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
    
     size1=size(cb%yXMax) ! yXMax
     call MPI_Unpack(e_para_vec, Nbytes, index, cb%yXMax(1,1),size1, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)  
  
     size1=size(cb%zXMax) ! zXMax
     call MPI_Unpack(e_para_vec, Nbytes, index, cb%zXMax(1,1),size1, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
    
     size1=size(cb%yXMin) ! yXMin
     call MPI_Unpack(e_para_vec, Nbytes, index, cb%yXMin(1,1),size1, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
    
     size1=size(cb%zXMin) ! zXMin
     call MPI_Unpack(e_para_vec, Nbytes, index, cb%zXMin(1,1),size1, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
    
     size1=size(cb%xZMin) ! xZMin
     call MPI_Unpack(e_para_vec, Nbytes, index, cb%xZMin(1,1),size1, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
    
     size1=size(cb%yZMin) ! yZMin
     call MPI_Unpack(e_para_vec, Nbytes, index, cb%yZMin(1,1),size1, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)     
    
     size1=size(cb%xZMax) ! xZMax
     call MPI_Unpack(e_para_vec, Nbytes, index, cb%xZMax(1,1),size1, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr) 
    
     size1=size(cb%yZMax) !yZMax
     call MPI_Unpack(e_para_vec, Nbytes, index, cb%xZMax(1,1),size1, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
	 
	 
   deallocate( e_para_vec )
     
     
     
 end subroutine UnPack_boundary_para_vec  
!********************************************************************  
subroutine create_e_param_place_holder(e)

     implicit none
     type(solnVector_t), intent(in)	:: e
     integer                        :: Ex_size,Ey_size,Ez_size,Nbytes1,Nbytes2,Nbytes3,Nbytes4






       Ex_size=size(e%pol(1)%x)
       CALL MPI_PACK_SIZE(Ex_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes1,  ierr)
       Ey_size=size(e%pol(1)%y)
       CALL MPI_PACK_SIZE(Ey_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes2,  ierr)
       Ez_size=size(e%pol(1)%z)
       CALL MPI_PACK_SIZE(Ez_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes3,  ierr)
       CALL MPI_PACK_SIZE(1, MPI_INTEGER, MPI_COMM_WORLD, Nbytes4,  ierr)
         Nbytes=((Nbytes1+Nbytes2+Nbytes3+Nbytes4))+1

         if( associated( e_para_vec ) ) deallocate( e_para_vec )
         allocate( e_para_vec( Nbytes ) )
		 !

 end subroutine create_e_param_place_holder
 !********************************************************************
 subroutine Pack_e_para_vec(e)
    implicit none

     type(solnVector_t), intent(in)	:: e
     integer index,Ex_size,Ey_size,Ez_size



       Ex_size=size(e%pol(1)%x)
       Ey_size=size(e%pol(1)%y)
       Ez_size=size(e%pol(1)%z)
       index=1

        call MPI_Pack(e%pol(which_pol)%x(1,1,1),Ex_size, MPI_DOUBLE_COMPLEX, e_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(e%pol(which_pol)%y(1,1,1),Ey_size, MPI_DOUBLE_COMPLEX, e_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(e%pol(which_pol)%z(1,1,1),Ez_size, MPI_DOUBLE_COMPLEX, e_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(e%tx,1,             MPI_INTEGER,        e_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
  

end subroutine Pack_e_para_vec
!********************************************************************
subroutine Unpack_e_para_vec(e)
    implicit none

     type(solnVector_t), intent(inout)	:: e

     integer index,Ex_size,Ey_size,Ez_size


       Ex_size=size(e%pol(1)%x)
       Ey_size=size(e%pol(1)%y)
       Ez_size=size(e%pol(1)%z)
       index=1

        call MPI_Unpack(e_para_vec, Nbytes, index, e%pol(which_pol)%x(1,1,1),Ex_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(e_para_vec, Nbytes, index, e%pol(which_pol)%y(1,1,1),Ey_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(e_para_vec, Nbytes, index, e%pol(which_pol)%z(1,1,1),Ez_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(e_para_vec, Nbytes, index, e%tx,1, MPI_INTEGER,MPI_COMM_WORLD, ierr)

        deallocate( e_para_vec )


end subroutine Unpack_e_para_vec

!********************************************************************
subroutine create_eAll_param_place_holder(e)

     implicit none
     type(solnVector_t), intent(in)	:: e
     integer Ex_size,Ey_size,Ez_size,Nbytes1,Nbytes2,Nbytes3



       Ex_size=size(e%pol(1)%x)
       CALL MPI_PACK_SIZE(Ex_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes1,  ierr)
       Ey_size=size(e%pol(1)%y)
       CALL MPI_PACK_SIZE(Ey_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes2,  ierr)
       Ez_size=size(e%pol(1)%z)
       CALL MPI_PACK_SIZE(Ez_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, Nbytes3,  ierr)

         Nbytes=(2*(Nbytes1+Nbytes2+Nbytes3))+1  ! Multiple by 2 for both polarizations

         if( associated( eAll_para_vec ) ) deallocate( eAll_para_vec )
         allocate( eAll_para_vec( Nbytes ) )
		 !



 end subroutine create_eAll_param_place_holder


! for now keep eAll_para_vec subroutines. Later we will not use them.
!********************************************************************
subroutine pack_eAll_para_vec(e)
    implicit none

     type(solnVector_t), intent(in)	:: e
     integer index,Ex_size,Ey_size,Ez_size



       Ex_size=size(e%pol(1)%x)
       Ey_size=size(e%pol(1)%y)
       Ez_size=size(e%pol(1)%z)
       index=1


        call MPI_Pack(e%pol(1)%x(1,1,1),Ex_size, MPI_DOUBLE_COMPLEX, eAll_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(e%pol(1)%y(1,1,1),Ey_size, MPI_DOUBLE_COMPLEX, eAll_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
        call MPI_Pack(e%pol(1)%z(1,1,1),Ez_size, MPI_DOUBLE_COMPLEX, eAll_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)

end subroutine pack_eAll_para_vec


!********************************************************************
subroutine Unpack_eAll_para_vec(e)
    implicit none

     type(solnVector_t), intent(inout)	:: e

     integer index,Ex_size,Ey_size,Ez_size


       Ex_size=size(e%pol(1)%x)
       Ey_size=size(e%pol(1)%y)
       Ez_size=size(e%pol(1)%z)
       index=1



        call MPI_Unpack(eAll_para_vec, Nbytes, index, e%pol(1)%x(1,1,1),Ex_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(eAll_para_vec, Nbytes, index, e%pol(1)%y(1,1,1),Ey_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)
        call MPI_Unpack(eAll_para_vec, Nbytes, index, e%pol(1)%z(1,1,1),Ez_size, MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD, ierr)

        deallocate( eAll_para_vec )
		
end subroutine Unpack_eAll_para_vec


#endif

end module Sub_MPI
