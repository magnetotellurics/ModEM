module Main_MPI_files
#ifdef MPI

  use math_constants
  use file_units
  use utilities
  use datasens	 
  use SolverSens  
  use ForwardSolver
  use SensComp
  use ioAscii
  use Declaration_MPI
  use Sub_MPI
      !use ioascii

  implicit none



  ! temporary EM fields, that are saved for efficiency - to avoid
  !  memory allocation & deallocation for each transmitter
  type(solnVector_t), save, private		    :: e,e0
  type(rhsVector_t) , save, private		    :: b0,comb
  type (grid_t), target, save, private     :: grid
  
  
Contains








!###########################################  MPI_initialization   ############################################################

Subroutine constructor_MPI

    implicit none
    include 'mpif.h'

          call MPI_INIT( ierr )
          call MPI_COMM_RANK( MPI_COMM_WORLD, taskid, ierr )
          call MPI_COMM_SIZE( MPI_COMM_WORLD, total_number_of_Proc, ierr )
          number_of_workers = total_number_of_Proc-1

End Subroutine constructor_MPI




!##############################################   Master_Job_FORWARD #########################################################



Subroutine Master_Job_fwdPred(sigma,d1,eAll)

    implicit none
    include 'mpif.h'
   type(modelParam_t), intent(in)	    :: sigma
   type(dataVectorMTX_t), intent(inout)	:: d1
   type(solnVectorMTX_t), intent(inout)	:: eAll
   integer nTx



   Integer        :: iper,ipol
   Integer        :: per_index,pol_index,stn_index,iTx,i,iDt,j,ios
   character(80)                        :: job_name
    character(20)                               :: which_proc,sys_command_line
  type(solnVectorMTX_t) eAll_temp
   Character(len=50) :: filename 



   ! nTX is number of transmitters;
   nTx = d1%nTx
   

     starttime = MPI_Wtime()




     write(115,*)"eAll%SolnIndex: ",eAll%SolnIndex,eAll%allocated


     ! First, distribute the current model to all workers
       call Master_job_Distribute_Model(sigma)
       !call Master_job_Distribute_Data(d1)
	   !call deall (eAll)
	   
      if(.not. eAll%allocated) then      
         call create_solnVectorMTX(1,eAll)
            do iTx=1,1
         		call create_solnVector(grid,iTx,e0)
        		call copy_solnVector(eAll%solns(iTx),e0) 
        	 end do 
      end if 	 

        
    
        job_name= 'FORWARD'      
        call Master_job_Distribute_Taskes(job_name,nTx,sigma,eAll)   

 ! Initialize only those grid elements on the master that are used in EMfieldInterp
 ! (obviously a quick patch, needs to be fixed in a major way)
 ! A.Kelbert 2018-01-28
    Call EdgeLength(grid, l_E)
    Call FaceArea(grid, S_F)
         
		 
 ! Compute the model Responces           
   do iTx=1,nTx
        !Read results in a file. Naser changes: 14.06.2022.
		 nPol_MPI=txDict(iTx)%nPol
		 eAll%solns(1)%tx=iTx
				  do ipol=1,nPol_MPI
				 	WRITE (filename, '(a,I3.3,a,I1,a,I1,a)') 'FWD_E_solution_Per',iTx, '_Pol',ipol,'_Index_',eAll%SolnIndex,'.soln'
					open (unit=ioE,file=filename,status='unknown', form ='unformatted',iostat=ios)
					call read_cvector(ioE, eAll%solns(1)%pol(ipol), 'binary')
					close(ioE)
				  end do	
	 
				  do i = 1,d1%d(iTx)%nDt
					 d1%d(iTx)%data(i)%errorBar = .false.
					 iDt = d1%d(iTx)%data(i)%dataType
						 do j = 1,d1%d(iTx)%data(i)%nSite
							call dataResp(eAll%solns(1),sigma,iDt,d1%d(iTx)%data(i)%rx(j),d1%d(iTx)%data(i)%value(:,j))
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
        !call deall (eAll_temp)
		!call deall (eAll)
end subroutine Master_Job_fwdPred


!##############################################   Master_Job_Compute_J #########################################################

Subroutine Master_job_calcJ(d,sigma,sens,eAll)

   implicit none
   type(modelParam_t), intent(in)	:: sigma
   type(dataVectorMTX_t), intent(in)		:: d
   type(solnVectorMTX_t), intent(in), optional	:: eAll
   type(sensMatrix_t), pointer 		      :: sens(:)

   !Local
    logical        :: savedSolns
    Integer        :: iper,idt,istn
    Integer        :: per_index,dt_index,stn_index,ipol1
    integer        :: nTx,nDt,nStn,dest,answers_to_receive,received_answers,nComp,nFunc,ii,iFunc,istat
    logical      :: isComplex  
    type(modelParam_t), pointer   :: Jreal(:),Jimag(:)           
    type(solnVectorMTX_t) eAll_temp
    type(dataVectorMTX_t)		:: d_Temp
   starttime = MPI_Wtime()
   
   
! now, allocate for sensitivity values, if necessary
if(.not. associated(sens)) then
     call create_sensMatrixMTX(d, sigma, sens)
endif
 

	         
! Check if an Esoln is passed    
savedSolns = present(eAll)  
if (.not. savedSolns )then
    d_temp=d
    call Master_Job_fwdPred(sigma,d_temp,eAll_temp)
else
      eAll_temp=eAll 
end if

dest=0
worker_job_task%what_to_do='COMPUTE_J' 
! now loop over Periods
nTx = d%nTx  
per_index=0 
do iper=1,nTx
     per_index=per_index+1
     worker_job_task%per_index= per_index
      call get_nPol_MPI(eAll_temp%solns(per_index)) 
     
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
	            call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)
				deallocate( worker_job_package )					
	            	 which_per=per_index

								        do ipol1=1,nPol_MPI
                                           which_pol=ipol1
			        				       call create_e_param_place_holder(eAll_temp%solns(which_per))
				    				       call Pack_e_para_vec(eAll_temp%solns(which_per))
				    				       call MPI_SEND(e_para_vec, Nbytes, MPI_PACKED, dest,FROM_MASTER, MPI_COMM_WORLD, ierr) 
										end do
	             write(ioMPI,8550) per_index ,dt_index,stn_index,dest  
	            if (dest .ge. number_of_workers) then 
	               goto 10
	            end if			                         
           end do
      end do
end do

10    continue   
           
! Count hom many rows of J we are going to recieve
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
! 1- Recieve an answer.
! 2 -Check if all stations, all data types and all periods haven been sent.
! 3- Send indicies if required.  

received_answers = 0
write(6,*)'answers_to_receive=', answers_to_receive
do while (received_answers .lt. answers_to_receive) 

!    Recieve worker INFO:                       
            call create_worker_job_task_place_holder
            call MPI_RECV(worker_job_package, Nbytes, MPI_PACKED ,MPI_ANY_SOURCE, FROM_WORKER,MPI_COMM_WORLD,STATUS, ierr)
            call Unpack_worker_job_task

                       who=worker_job_task%taskid
                       which_per=worker_job_task%per_index
                       which_dt=worker_job_task%data_type_index                      
                       which_stn=worker_job_task%Stn_index
                       write(ioMPI,8552) which_per ,which_dt,which_stn,who             
         
!    Recieve results from a worker:

! Store the result in sens

nComp = d%d(which_per)%data(which_dt)%nComp           
isComplex = d%d(which_per)%data(which_dt)%isComplex

		    if(isComplex) then
		       !  data are complex; one sensitivity calculation can be
		       !   used for both real and imaginary parts
		       if(mod(nComp,2).ne.0) then
		         call errStop('for complex data # of components must be even in calcJ')
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
            call MPI_RECV(sigma_para_vec, Nbytes, MPI_PACKED, who, FROM_WORKER,MPI_COMM_WORLD, STATUS, ierr)
            call unpack_model_para_values(Jreal(iFunc))
            
            !image part
            call create_model_param_place_holder(Jimag(iFunc))
            call MPI_RECV(sigma_para_vec, Nbytes, MPI_PACKED, who, FROM_WORKER,MPI_COMM_WORLD, STATUS, ierr)
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

if (Per_index ==  nTx .and. dt_index == d%d(Per_index)%nDt .and. stn_index == d%d(Per_index)%data(dt_index)%nSite ) goto 300 


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
 write(ioMPI,8551) per_index ,worker_job_task%data_type_index,stn_index,who    
             
! Send Indices to who (the worker who just send back an answer)
            call create_worker_job_task_place_holder
            call Pack_worker_job_task
            call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,who, FROM_MASTER, MPI_COMM_WORLD, ierr)
			deallocate( worker_job_package )					   
            

								        do ipol1=1,nPol_MPI
                                           which_pol=ipol1
			        				       call create_e_param_place_holder(eAll_temp%solns(per_index))
				    				       call Pack_e_para_vec(eAll_temp%solns(per_index))
				    				       call MPI_SEND(e_para_vec, Nbytes, MPI_PACKED, who,FROM_MASTER, MPI_COMM_WORLD, ierr) 
										end do          
                   
    
 

300    continue                                   
                   
end do

 8550  FORMAT('COMPUTE_J: Send Per. # ',i5, ' , data type # ',i5, ' and Site #',i5,' to node # ',i5  )
 8551  FORMAT('COUNT--> COMPUTE_J: Send Per. # ',i5, ' , data type # ',i5, ' and Site #',i5,' to node # ',i5  )
 8552  FORMAT('COMPUTE_J: RECV Per. # ',i5, ' , data type # ',i5, ' and Site #',i5,' from node # ',i5  )

        endtime=MPI_Wtime()
        time_used = endtime-starttime
        !DONE: Received soln for all transmitter from all nodes
        write(ioMPI,*)'COMPUTE_J: Finished computing for (',answers_to_receive , ';[Transmitters*data type*site]) '
        endtime=MPI_Wtime()
        time_used = endtime-starttime
        write(ioMPI,*)'COMPUTE_J: TIME REQUIERED: ',time_used ,'s'





end subroutine Master_job_calcJ



!##############################################    Master_job_JmultT #########################################################
Subroutine Master_job_JmultT(sigma,d,dsigma,eAll,s_hat)

   implicit none
    include 'mpif.h'

   type(modelParam_t), intent(in)	:: sigma
   type(dataVectorMTX_t), intent(in)		:: d
   type(modelParam_t), intent(Out)  	:: dsigma
   type(solnVectorMTX_t), intent(in), optional	:: eAll
   type(modelParam_t),intent(inout),pointer,dimension(:), optional :: s_hat

   ! Local
   type(modelParam_t)           :: dsigma_temp
   type(modelParam_t)           :: Qcomb
   type(solnVectorMTX_t)      	:: eAll_out 
   type(solnVectorMTX_t)      	:: eAll_temp
   type(dataVectorMTX_t)    	  :: d_temp
   
   logical        :: savedSolns,returne_m_vectors
   Integer        :: iper,ipol,nTx,iTx
   Integer        :: per_index,pol_index,stn_index,ios
   character(80)  :: job_name 
   Character(len=50) :: filename 
 

   
   savedSolns = present(eAll)
   returne_m_vectors= present(s_hat)
  ! nTX is number of transmitters;
   nTx = d%nTx
      !if(.not. eAll_temp%allocated) then
         call create_solnVectorMTX(1,eAll_temp)
            do iTx=1,1
         		call create_solnVector(grid,iTx,e0)
        		call copy_solnVector(eAll_temp%solns(iTx),e0) 
        	 end do 
      !end if 
	  
	  if (.not. savedSolns )then
			d_temp=d
			eAll_temp%SolnIndex=0
			call Master_Job_fwdPred(sigma,d_temp,eAll_temp)
		else
			  eAll_temp=eAll 
			  eAll_temp%SolnIndex=eAll%SolnIndex
	  end if
		
		
      !if(.not. eAll_out%allocated) then
         call create_solnVectorMTX(1,eAll_out)
            do iTx=1,1
         		call create_solnVector(grid,iTx,e0)
        		call copy_solnVector(eAll_out%solns(iTx),e0) 
                call deall (e0)  
        	 end do 
     ! end if 


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
    
   !First ditribute both model parameters and data
        call Master_job_Distribute_Model(sigma)
        call Master_job_Distribute_Data(d)

   dsigma_temp = sigma
   dsigma 	   = sigma
   call zero(dsigma_temp)
   call zero(dsigma)
   Qcomb = sigma
   call zero(Qcomb)
 
  
  
  
  
       job_name= 'JmultT'
       call Master_job_Distribute_Taskes(job_name,nTx,sigma,eAll_out,eAll_temp)
       


			 
          do iper=1,nTx
            !e0=eAll_temp%solns(iper)  
            !e =eAll_out%solns(iper)
			
			
				 nPol_MPI=txDict(iper)%nPol
				 eAll_temp%solns(1)%tx=iper
				 
				  do ipol=1,nPol_MPI
				 	WRITE (filename, '(a,I3.3,a,I1,a,I1,a)') 'FWD_E_solution_Per',iper, '_Pol',ipol,'_Index_',eAll_temp%SolnIndex,'.soln'
					open (unit=ioE,file=filename,status='unknown', form ='unformatted',iostat=ios)
					call read_cvector(ioE, eAll_temp%solns(1)%pol(ipol), 'binary')
					close(ioE)
				  end do	
					
				 eAll_out%solns(1)%tx=iper
				 do ipol=1,nPol_MPI
					WRITE (filename, '(a,I3.3,a,I1,a)') 'JmultT_E_solution_Per',iper, '_Pol',ipol,'.soln'
					open (unit=ioE,file=filename,status='unknown', form ='unformatted',iostat=ios)
					call read_cvector(ioE, eAll_out%solns(1)%pol(ipol), 'binary')
					close(ioE)
				 end do
	 
			
			
			
            call PmultT(eAll_temp%solns(1)  ,sigma,eAll_out%solns(1),dsigma_temp)
            call QmultT(eAll_temp%solns(1)  ,sigma,d%d(iper),Qcomb)
    		    call scMultAdd(ONE,Qcomb,dsigma_temp)
    		         if (returne_m_vectors) then
                       s_hat(iper)=dsigma_temp
                 end if
    		   call linComb_modelParam(ONE,dsigma,ONE,dsigma_temp,dsigma)
         end do		

                endtime=MPI_Wtime()
                time_used = endtime-starttime
        !DONE: Received soln for all transmitter from all nodes
        write(ioMPI,*)'JmultT: Finished calculating for (', d%nTx , ') Transmitters '
        endtime=MPI_Wtime()
        time_used = endtime-starttime
        write(ioMPI,*)'JmultT: TIME REQUIERED: ',time_used ,'s'



   !  clean up
   call deall_modelParam(dsigma_temp)
   call deall_modelParam(Qcomb)
   call deall (eAll_out)
   call deall (eAll_temp) 
   !call deall (e0)   
   call deall_dataVectorMTX(d_temp)
   

end Subroutine Master_job_JmultT

!##############################################    Master_job_Jmult #########################################################
Subroutine Master_job_Jmult(mHat,m,d,eAll)

    implicit none
    include 'mpif.h'

   type(dataVectorMTX_t), intent(inout)		:: d
   type(modelParam_t), intent(in)			:: mHat,m
   type(solnVectorMTX_t), intent(in), optional	:: eAll

   ! Local
   type(modelParam_t)           :: dsigma_temp
   type(modelParam_t)           :: Qcomb
   type(solnVectorMTX_t)        :: eAll_out
   type(solnVectorMTX_t)        :: eAll_temp
   type(dataVectorMTX_t)        :: d_temp

   integer nTx,nTot,m_dimension,iDT,iTx,ndata,ndt
   logical savedSolns
   Integer        :: iper
   Integer        :: per_index,pol_index,stn_index
   type(dataVector_t) :: d1,d2
   character(80)  :: job_name

   savedSolns = present(eAll)
   ! nTot total number of data points
   nTot = countData(d) !d%Ndata
   starttime = MPI_Wtime()
	  !  initialize the temporary data vectors
	  d1 = d%d(1)
	  d2 = d%d(1) 

  ! nTX is number of transmitters;
   nTx = d%nTx
      if(.not. eAll_temp%allocated) then
         call create_solnVectorMTX(d%nTx,eAll_temp)
            do iTx=1,nTx
                call create_solnVector(grid,iTx,e0)
                call copy_solnVector(eAll_temp%solns(iTx),e0)
             end do
      end if
if (.not. savedSolns )then
    d_temp=d
    call Master_Job_fwdPred(m,d_temp,eAll_temp)
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
   ! First distribute m, mHat and d
    	    call Master_job_Distribute_Model(m,mHat)
	        call Master_job_Distribute_Data(d)


       job_name= 'Jmult'
       call Master_job_Distribute_Taskes(job_name,nTx,m,eAll_out,eAll_temp)

  



          do iper=1,nTx
            !e0=eAll%solns(iper)  
            !e =eAll_out%solns(iper)
            d1 = d%d(iper)
	        d2 = d%d(iper)
	        call Lmult(eAll_temp%solns(iper)  ,m,eAll_out%solns(iper),d1)
	        call Qmult(eAll_temp%solns(iper)  ,m,mHat,d2)
	        call linComb_dataVector(ONE,d1,ONE,d2,d%d(iper))
         end do	
         

  	  call deall_dataVector(d1)
	  call deall_dataVector(d2)       
   call deall (eAll_out)
   call deall (eAll_temp)
   call deall (e0)  
        !DONE: Received soln for all transmitter from all nodes
        write(ioMPI,*)'Jmult: Finished calculating for (', d%nTx , ') Transmitters '

                endtime=MPI_Wtime()
                time_used = endtime-starttime
        write(ioMPI,*)'Jmult: TIME REQUIERED: ',time_used ,'s'
end Subroutine Master_job_Jmult

!############################################## Master_job_Distribute_Data #########################################################
Subroutine Master_job_Distribute_Data(d)
    implicit none
    include 'mpif.h'
    type(dataVectorMTX_t), intent(in)		:: d
    integer                                 :: nTx

       nTx=d%nTx
       
        do dest=1,number_of_workers
            worker_job_task%what_to_do='Distribute nTx'
            call create_worker_job_task_place_holder
            call Pack_worker_job_task
            call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)
			deallocate( worker_job_package )					   
        end do
        
          call MPI_BCAST(nTx,1, MPI_INTEGER,0, MPI_COMM_WORLD,ierr)


           

        do dest=1,number_of_workers
            worker_job_task%what_to_do='Distribute Data'
            call create_worker_job_task_place_holder
            call Pack_worker_job_task
            call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)
			deallocate( worker_job_package )					   
        end do






            call create_data_vec_place_holder(d)
            call Pack_data_para_vec(d)
            call MPI_BCAST(data_para_vec,Nbytes, MPI_PACKED,0, MPI_COMM_WORLD,ierr)
            deallocate (data_para_vec)



end Subroutine Master_job_Distribute_Data


!############################################## Master_job_Distribute_Model #########################################################
Subroutine Master_job_Distribute_Model(sigma,delSigma)
    implicit none
    include 'mpif.h'
    type(modelParam_t), intent(in) 	:: sigma
    type(modelParam_t), intent(in), optional :: delSigma
    !local
    type(modelParam_t)          	:: sigma_temp
    Integer,pointer:: buffer(:)
    integer buffer_size,m_dimension
    logical send_delSigma
    Integer        :: iper

 send_delSigma = present(delSigma)

 if (send_delSigma) then
         do dest=1,number_of_workers
            worker_job_task%what_to_do='Distribute delSigma'
            call create_worker_job_task_place_holder
            call Pack_worker_job_task
            call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)
			deallocate( worker_job_package )					   
        end do
         call create_model_param_place_holder(delSigma)
         call pack_model_para_values(delSigma)
         call MPI_BCAST(sigma_para_vec,Nbytes, MPI_PACKED,0, MPI_COMM_WORLD,ierr)


        do dest=1,number_of_workers
            worker_job_task%what_to_do='Distribute Model'
            call create_worker_job_task_place_holder
            call Pack_worker_job_task
            call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)
			deallocate( worker_job_package )					   
        end do

         call create_model_param_place_holder(sigma)
         call pack_model_para_values(sigma)
         call MPI_BCAST(sigma_para_vec,Nbytes, MPI_PACKED,0, MPI_COMM_WORLD,ierr)
		 !deallocate(sigma_para_vec)

else
        do dest=1,number_of_workers
            worker_job_task%what_to_do='Distribute Model'
            call create_worker_job_task_place_holder
            call Pack_worker_job_task
            call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)
        end do

         call create_model_param_place_holder(sigma)
         call pack_model_para_values(sigma)
         call MPI_BCAST(sigma_para_vec,Nbytes, MPI_PACKED,0, MPI_COMM_WORLD,ierr)
		 !deallocate(sigma_para_vec)
end if

end Subroutine Master_job_Distribute_Model


!############################################## Master_job_Distribute_eAll #########################################################
Subroutine Master_job_Distribute_eAll(d,eAll)
    implicit none
    include 'mpif.h'
   type(dataVectorMTX_t), intent(in)		:: d
   type(solnVectorMTX_t), intent(in)  	:: eAll
       integer nTx,nTot
       Integer        :: iper



    nTx = d%nTx
        do dest=1,number_of_workers
            worker_job_task%what_to_do='Distribute eAll'
            call create_worker_job_task_place_holder
            call Pack_worker_job_task
            call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)
			deallocate( worker_job_package )					   
        end do

  do iper=1,d%nTx
       which_per=iper
       do dest=1,number_of_workers
          call create_eAll_param_place_holder(e0)
          call Pack_eAll_para_vec(e0)
          call MPI_SEND(eAll_para_vec, Nbytes, MPI_PACKED,dest, FROM_MASTER,MPI_COMM_WORLD, ierr)
		 deallocate( eAll_para_vec )					   
        end do
  end do
end Subroutine Master_job_Distribute_eAll

!############################################## Master_job_Collect_eAll #########################################################
Subroutine Master_job_Collect_eAll(d,eAll)
    implicit none
    include 'mpif.h'
   type(dataVectorMTX_t), intent(in)		:: d
   type(solnVectorMTX_t), intent(inout)	:: eAll
   integer nTx,nTot,iTx
   Integer        :: iper

    nTx = d%nTx





      if(.not. eAll%allocated) then
         call create_solnVectorMTX(d%nTx,eAll)
      else if(d%nTx .ne. eAll%nTx) then
         call errStop('dimensions of eAll and d do not agree in Master_job_Collect_eAll')
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
            call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,who, FROM_MASTER, MPI_COMM_WORLD, ierr)
			deallocate( worker_job_package )					   
            call create_eAll_param_place_holder(e0)
            call MPI_RECV(eAll_para_vec, Nbytes, MPI_PACKED, who, FROM_WORKER,MPI_COMM_WORLD, STATUS, ierr)
            call Unpack_eAll_para_vec(e0)
      end do


end Subroutine Master_job_Collect_eAll
!############################################## Master_job_keep_prev_eAll #########################################################
subroutine Master_job_keep_prev_eAll
    implicit none
    include 'mpif.h'

        do dest=1,number_of_workers
            worker_job_task%what_to_do='keep_prev_eAll'
            call create_worker_job_task_place_holder
            call Pack_worker_job_task
            call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)
			deallocate( worker_job_package )					   
        end do

end  subroutine Master_job_keep_prev_eAll



!############################################## Master_job_Distribute_userdef_control#########################################################
Subroutine Master_job_Distribute_userdef_control(ctrl)
	!
	implicit none
	include 'mpif.h'
	!
	type(userdef_control), intent(in)		:: ctrl
	character(20)                           :: which_proc
	!
	which_proc='Master'
	!
	call check_userdef_control_MPI (which_proc,ctrl)
	!
	call create_userdef_control_place_holder
	!
	call pack_userdef_control(ctrl)
	!
	do dest=1,number_of_workers
		call MPI_SEND(userdef_control_package,Nbytes, MPI_PACKED,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)
	end do
	!
	deallocate( userdef_control_package )								  
end Subroutine Master_job_Distribute_userdef_control

!##############################################    Master_job_Clean Memory ########################################################
Subroutine Master_job_Clean_Memory

    implicit none
    include 'mpif.h'

          write(ioMPI,*)'Sending Clean memory message to all nodes'

       do dest=1,number_of_workers
           worker_job_task%what_to_do='Clean memory'
            call create_worker_job_task_place_holder
            call Pack_worker_job_task
            call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)
            call MPI_RECV(worker_job_package, Nbytes, MPI_PACKED ,dest, FROM_WORKER,MPI_COMM_WORLD,STATUS, ierr)
            call Unpack_worker_job_task
           write(ioMPI,*)'Node       :',worker_job_task%taskid, ' status=  ', trim(worker_job_task%what_to_do)

        end do

end Subroutine Master_job_Clean_Memory

!##############################################    Master_job_Stop_MESSAGE ########################################################
Subroutine Master_job_Stop_MESSAGE

    implicit none
    include 'mpif.h'

          write(ioMPI,*)'FWD: Sending stop message to all nodes'

       do dest=1,number_of_workers
           worker_job_task%what_to_do='STOP'
            call create_worker_job_task_place_holder
            call Pack_worker_job_task
            call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)
            call MPI_RECV(worker_job_package, Nbytes, MPI_PACKED ,dest, FROM_WORKER,MPI_COMM_WORLD,STATUS, ierr)
            call Unpack_worker_job_task
           write(ioMPI,*)'Node       :',worker_job_task%taskid, ' status=  ', trim(worker_job_task%what_to_do)

        end do

end Subroutine Master_job_Stop_MESSAGE













!############################################################   Worker_job :High Level Subroutine   #####################################################################
Subroutine Worker_job (sigma,d)
    implicit none
    include 'mpif.h'



   type(modelParam_t),intent(inout)	            :: sigma
   type(dataVectorMTX_t) ,intent(inout)    	    :: d
   
   
   
   !Local 
   type(modelParam_t)           	            :: delSigma
   type(userdef_control)                        :: ctrl

      
   Integer nTx,m_dimension,ndata,itx,ndt,dt_index,per_index_pre,dt
   character(80) 		  :: paramType,previous_message


   Integer        :: iper,ipol,i
   Integer        :: per_index,pol_index,stn_index,eAll_vec_size,eAll_SolnIndex
   character(20)                               :: which_proc
 
   type(modelParam_t), pointer   :: Jreal(:),Jimag(:)
   integer                       ::nComp,nFunc,iFunc,istat,ios
   logical                       :: isComplex  
   type(sparseVector_t), pointer	:: L(:)
   type(modelParam_t), pointer    :: Qreal(:),Qimag(:)
   logical      :: Qzero
   Character(len=50) :: filename 
   type(solnVectorMTX_t)	:: eAll_temp
   character(20)                               :: sys_command_line
    
nTx=d%nTx
recv_loop=0
previous_message=''

write(node_info,'(a5,i3.3,a4)') 'node[',taskid,']:  '

 do
          recv_loop=recv_loop+1
          write(6,'(a12,a35)') node_info,' Waiting for a message from Master'
          call create_worker_job_task_place_holder
          call MPI_RECV(worker_job_package, Nbytes, MPI_PACKED ,0, FROM_MASTER,MPI_COMM_WORLD,STATUS, ierr)
          call Unpack_worker_job_task

          !Receive message including what to do and another info. requiered (i.e per_index_stn_index, ...,etc)
          !call MPI_RECV(worker_job_task,1,worker_job_task_mpi,0, FROM_MASTER,MPI_COMM_WORLD, STATUS, ierr)


			write(6,'(a12,a12,a30,a16,i5)') node_info,' MPI TASK [',trim(worker_job_task%what_to_do),'] received from ',STATUS(MPI_SOURCE)
			!write(6,*) node_info,' MPI INFO [keep soln = ',(worker_job_task%keep_E_soln), &
			! '; several TX = ',worker_job_task%several_Tx,']'
            

if (trim(worker_job_task%what_to_do) .eq. 'FORWARD') then



          per_index=worker_job_task%per_index
          pol_index=worker_job_task%pol_index
		  eAll_SolnIndex=worker_job_task%eAll_SolnIndex
          worker_job_task%taskid=taskid

		       call initSolver(per_index,sigma,grid,e0)
		       call set_e_soln(pol_index,e0)
		       

		       call fwdSetup(per_index,e0,b0)
		       call fwdSolve(per_index,e0,b0)
               call reset_e_soln(e0)

              
                    !Write results in a file. Naser changes: 14.06.2021.  
					
					WRITE (filename, '(a,I3.3,a,I1,a,I1,a)') 'FWD_E_solution_Per',per_index, '_Pol',pol_index,'_Index_',eAll_SolnIndex,'.soln'
					open (unit=ioE,file=filename,status='unknown', form ='unformatted',iostat=ios)
					call write_cvector(ioE, e0%pol(1), 'binary')
					close(ioE)
					
					! call create_solnVectorMTX(1,eAll_temp)
					! call create_solnVector(grid,1,e0)
					! call copy_solnVector(eAll_temp%solns(1),e0) 
					! call write_solnVectorMTX(eAll_temp,filename)
					! call deall (eAll_temp) 
				  	
					
                    					
        
 		      ! Create worker job package and send it to the master
		            call create_worker_job_task_place_holder
		            call Pack_worker_job_task
		            call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,0,FROM_WORKER, MPI_COMM_WORLD, ierr)
		            deallocate( worker_job_package )

		   
		      ! Create e0_temp package (one Period and one Polarization) and send it to the master
              which_pol=1
			  !STOP Sending e solution from workers. Naser changes: 14.06.2021. The following 3 lines are deactivated
		      !call create_e_param_place_holder(e0) 
		      !call Pack_e_para_vec(e0)
		      !call MPI_SEND(e_para_vec, Nbytes, MPI_PACKED, 0,FROM_WORKER, MPI_COMM_WORLD, ierr) 

              !deallocate(e_para_vec,worker_job_package)
              call exitSolver(e0)

elseif (trim(worker_job_task%what_to_do) .eq. 'COMPUTE_J') then

          per_index=worker_job_task%per_index
          stn_index=worker_job_task%stn_index
          dt_index=worker_job_task%data_type_index
          dt=worker_job_task%data_type
          worker_job_task%taskid=taskid
          
nComp = d%d(per_index)%data(dt_index)%nComp           
isComplex = d%d(per_index)%data(dt_index)%isComplex

		    if(isComplex) then
		       !  data are complex; one sensitivity calculation can be
		       !   used for both real and imaginary parts
		       if(mod(nComp,2).ne.0) then
		         call errStop('for complex data # of components must be even in calcJ')
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
              	              
! Do some computation
                    call initSolver(per_index,sigma,grid,e0,e,comb) 
                   call get_nPol_MPI(e0)
                    write(6,'(a12,a18,i5,a12)') node_info, ' Start Receiving ' , orginal_nPol, ' from Master'
		          do ipol=1,nPol_MPI 
                    which_pol=ipol				  
		            call create_e_param_place_holder(e0)
		            call MPI_RECV(e_para_vec, Nbytes, MPI_PACKED, 0, FROM_MASTER,MPI_COMM_WORLD, STATUS, ierr)
                    call Unpack_e_para_vec(e0)
                  end do
                      write(6,'(a12,a18,i5,a12)') node_info, ' Finished Receiving ' , orginal_nPol, ' from Master'

   allocate(L(nFunc),STAT=istat)
   allocate(Qreal(nFunc),STAT=istat)
   allocate(Qimag(nFunc),STAT=istat)
   
	  do iFunc=1,nFunc
		  call create_sparseVector(e0%grid,per_index,L(iFunc))
	  end do
	  
 ! Initialize only those grid elements on the master that are used in EMfieldInterp
 ! (obviously a quick patch, needs to be fixed in a major way)
 ! A.Kelbert 2018-01-28
    Call EdgeLength(e0%grid, l_E)
    Call FaceArea(e0%grid, S_F)

   ! compute linearized data functional(s) : L
   call Lrows(e0,sigma,dt,stn_index,L)
   ! compute linearized data functional(s) : Q
   call Qrows(e0,sigma,dt,stn_index,Qzero,Qreal,Qimag)

 ! clean up the grid elements stored in GridCalc on the master node
    call deall_rvector(l_E)
    call deall_rvector(S_F)

   ! loop over functionals  (e.g., for 2D TE/TM impedances nFunc = 1)
   do iFunc = 1,nFunc

      ! solve transpose problem for each of nFunc functionals
      call zero_rhsVector(comb)
      call add_sparseVrhsV(C_ONE,L(iFunc),comb)

      call sensSolve(per_index,TRN,e,comb)

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

   enddo  ! iFunc

   !  deallocate local arrays
   deallocate(L,STAT=istat)
   deallocate(Qreal,STAT=istat)
   deallocate(Qimag,STAT=istat)
   		                            
                ! call Jrows(per_index,dt_index,stn_index,sigma,e0,Jreal,Jimag)        

                
                
 		      ! Create worker job package and send it to the master
		            call create_worker_job_task_place_holder
		            call Pack_worker_job_task
		            call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,0,FROM_WORKER, MPI_COMM_WORLD, ierr)
					deallocate( worker_job_package )				 
		            
			do iFunc = 1,nFunc 	            
		      ! Create worker model  package for Jreal and send it to the master       
                   call create_model_param_place_holder(Jreal(iFunc))
                   call pack_model_para_values(Jreal(iFunc))
                   call MPI_SEND(sigma_para_vec, Nbytes, MPI_PACKED, 0,FROM_WORKER, MPI_COMM_WORLD, ierr)
                   
		      ! Create worker model  package for Jimag and send it to the master       
                   call create_model_param_place_holder(Jimag(iFunc))
                   call pack_model_para_values(Jimag(iFunc))
                   call MPI_SEND(sigma_para_vec, Nbytes, MPI_PACKED, 0,FROM_WORKER, MPI_COMM_WORLD, ierr)
            end do    
            
            call exitSolver(e0,e,comb)                      		                      
elseif (trim(worker_job_task%what_to_do) .eq. 'JmultT') then


                       per_index=worker_job_task%per_index
                       pol_index=worker_job_task%pol_index
					   eAll_SolnIndex=worker_job_task%eAll_SolnIndex
                       worker_job_task%taskid=taskid
            
                    call initSolver(per_index,sigma,grid,e0,e,comb)   
                    call get_nPol_MPI(e0)
                    write(6,'(a12,a18,i5,a12)') node_info, ' Start Receiving ' , orginal_nPol, ' from Master'
		          do ipol=1,nPol_MPI 
                    which_pol=ipol				  
		            ! call create_e_param_place_holder(e0)
		            ! call MPI_RECV(e_para_vec, Nbytes, MPI_PACKED, 0, FROM_MASTER,MPI_COMM_WORLD, STATUS, ierr)
                    ! call Unpack_e_para_vec(e0)
					
						WRITE (filename, '(a,I3.3,a,I1,a,I1,a)') 'FWD_E_solution_Per',per_index, '_Pol',ipol,'_Index_',eAll_SolnIndex,'.soln'
						open (unit=ioE,file=filename,status='unknown', form ='unformatted',iostat=ios)
						call read_cvector(ioE, e0%pol(ipol), 'binary')
						close(ioE)
						
						
						! WRITE (filename, '(a,I3.3,a,I1,a,I1,a)') 'FWD_E_solution_Per',per_index, '_Pol',2,'_Index_',eAll_SolnIndex,'.soln'
						! open (unit=ioE,file=filename,status='unknown', form ='unformatted',iostat=ios)
						! call read_cvector(ioE, eAll%solns(1)%pol(2), 'binary')
						! close(ioE)						
					
					
					
                  end do
				  
                      write(6,'(a12,a18,i5,a12)') node_info, ' Finished Receiving ' , orginal_nPol, ' from Master'

            call LmultT(e0,sigma,d%d(per_index),comb)
            call set_e_soln(pol_index,e)
            call sensSolve(per_index,TRN,e,comb)
            call reset_e_soln(e)
          
					WRITE (filename, '(a,I3.3,a,I1,a)') 'JmultT_E_solution_Per',per_index, '_Pol',pol_index,'.soln'
					open (unit=ioE,file=filename,status='unknown', form ='unformatted',iostat=ios)
					call write_cvector(ioE, e%pol(1), 'binary')
					close(ioE)
  

             ! Send Info. about the current worker.
            call create_worker_job_task_place_holder
            call Pack_worker_job_task
            call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,0,FROM_WORKER, MPI_COMM_WORLD, ierr)
			deallocate( worker_job_package )					   


                   which_pol=1
                   ! call create_e_param_place_holder(e)
                   ! call Pack_e_para_vec(e)
                   ! call MPI_SEND(e_para_vec, Nbytes, MPI_PACKED, 0,FROM_WORKER, MPI_COMM_WORLD, ierr)
                   
                 !deallocate(e_para_vec,worker_job_package)
                 call exitSolver(e0,e,comb)

                   
elseif (trim(worker_job_task%what_to_do) .eq. 'Jmult') then

                       per_index=worker_job_task%per_index
                       pol_index=worker_job_task%pol_index
                       worker_job_task%taskid=taskid
	 
                    call initSolver(per_index,sigma,grid,e0,e,comb) 
		           
                   
                    write(6,'(a12,a18,i5,a12)') node_info, ' Start Receiving ' , orginal_nPol, ' from Master'
		          do ipol=1,orginal_nPol 
                    which_pol=ipol				  
		            call create_e_param_place_holder(e0)
		            call MPI_RECV(e_para_vec, Nbytes, MPI_PACKED, 0, FROM_MASTER,MPI_COMM_WORLD, STATUS, ierr)
		            call Unpack_e_para_vec(e0)
                  end do
                      write(6,'(a12,a18,i5,a12)') node_info, ' Finished Receiving ' , orginal_nPol, ' from Master'
		            
               
            
           
            call Pmult(e0,sigma,delSigma,comb)
            call set_e_soln(pol_index,e)
	        call sensSolve(per_index,FWD,e,comb)
            call reset_e_soln(e)

	                                                    
          
  

             ! Send Info. about the current worker.
            call create_worker_job_task_place_holder
            call Pack_worker_job_task
            call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,0,FROM_WORKER, MPI_COMM_WORLD, ierr)
			deallocate( worker_job_package )					   


                   which_pol=1
                   call create_e_param_place_holder(e)
                   call Pack_e_para_vec(e)
                   call MPI_SEND(e_para_vec, Nbytes, MPI_PACKED, 0,FROM_WORKER, MPI_COMM_WORLD, ierr)
				   deallocate( e_para_vec )		   

				   call exitSolver(e0,e,comb)

elseif (trim(worker_job_task%what_to_do) .eq. 'Distribute nTx') then
     call MPI_BCAST(nTx,1, MPI_INTEGER,0, MPI_COMM_WORLD,ierr)
     d%nTx=nTx
elseif (trim(worker_job_task%what_to_do) .eq. 'Distribute Data') then


            call create_data_vec_place_holder(d)
            call MPI_BCAST(data_para_vec,Nbytes, MPI_PACKED,0, MPI_COMM_WORLD,ierr)
            call UnPack_data_para_vec(d)
            !deallocate (data_para_vec)           
              
elseif (trim(worker_job_task%what_to_do) .eq. 'Distribute eAll') then

 
              do iper=1,d%nTx
                  which_per=iper
                  call create_eAll_param_place_holder(e0)
                  call MPI_RECV(eAll_para_vec, Nbytes, MPI_PACKED ,0, FROM_MASTER,MPI_COMM_WORLD,STATUS, ierr)
                  call Unpack_eAll_para_vec(e0)
              end do
         eAll_exist=.true.

elseif (trim(worker_job_task%what_to_do) .eq. 'Distribute Model') then

            call create_model_param_place_holder(sigma)
            call MPI_BCAST(sigma_para_vec,Nbytes, MPI_PACKED,0, MPI_COMM_WORLD,ierr)
            call unpack_model_para_values(sigma)
            !deallocate (sigma_para_vec)

elseif (trim(worker_job_task%what_to_do) .eq. 'Distribute delSigma') then

            call create_model_param_place_holder(sigma)
            call MPI_BCAST(sigma_para_vec,Nbytes, MPI_PACKED,0, MPI_COMM_WORLD,ierr)
            call copy_ModelParam(delSigma,sigma)
            call unpack_model_para_values(delSigma)
            !deallocate (sigma_para_vec)

elseif (trim(worker_job_task%what_to_do) .eq. 'Send eAll to Master' ) then



                   per_index=worker_job_task%per_index
                   worker_job_task%taskid=taskid

                   which_per=per_index
                   call create_eAll_param_place_holder(e0)
                   call Pack_eAll_para_vec(e0)
                   call MPI_SEND(eAll_para_vec, Nbytes, MPI_PACKED, 0,FROM_WORKER, MPI_COMM_WORLD, ierr)
					 deallocate( eAll_para_vec )			  


elseif (trim(worker_job_task%what_to_do) .eq. 'Clean memory' ) then




         worker_job_task%what_to_do='Cleaned Memory and Waiting'
         worker_job_task%taskid=taskid
            call create_worker_job_task_place_holder
            call Pack_worker_job_task
            call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,0,FROM_WORKER, MPI_COMM_WORLD, ierr)
			deallocate( worker_job_package )					   


elseif (trim(worker_job_task%what_to_do) .eq. 'STOP' ) then

                             worker_job_task%what_to_do='Job Completed'
                             worker_job_task%taskid=taskid
            call create_worker_job_task_place_holder
            call Pack_worker_job_task
            call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,0,FROM_WORKER, MPI_COMM_WORLD, ierr)
			deallocate( worker_job_package )					   
                             exit

end if
!previous_message=trim(worker_job_task%what_to_do)
!write(6,'(a12,a12,a30,a12)') node_info,' MPI TASK [',trim(worker_job_task%what_to_do),'] successful'
worker_job_task%what_to_do='Waiting for new message'

end do

call deall_solnVector(e0)
call deall_rhsVector(b0)						

End Subroutine Worker_job

subroutine Master_job_Distribute_Taskes(job_name,nTx,sigma,eAll_out,eAll_in)
     implicit none
        character(80) , intent(in)                          :: job_name
        Integer    , intent(in)                            :: nTx
        type(modelParam_t),intent(in)	            :: sigma
        type(solnVectorMTX_t), intent(in), optional	        :: eAll_in
        type(solnVectorMTX_t), intent(inout), optional	    :: eAll_out     
        !Local
        Integer        :: iper,ipol,ipol1,itx
        Integer        :: per_index,pol_index,solver_residual_iter
        logical keep_soln,savedSolns

        
        

 savedSolns = present(eAll_in)
  
                dest=0
                per_index=0
                worker_job_task%what_to_do=trim(job_name) 
                    do iper=1,nTx
                          per_index=per_index+1
                          worker_job_task%per_index= per_index
						   worker_job_task%eAll_SolnIndex=eAll_out%SolnIndex
                          pol_index=0
                          !call get_nPol_MPI(eAll_out%solns(per_index)) 
						  nPol_MPI=txDict(per_index)%nPol
                          do ipol=1,nPol_MPI
                              pol_index=pol_index+1
                              worker_job_task%pol_index= pol_index
	                          dest=dest+1

                  	             
								 						 
	            			     if (trim(job_name).eq. 'JmultT' .or. trim(job_name).eq. 'Jmult')then
								 ! In case of JmultT and Jmult the compleate eAll solution(the background solution) must be send to each node.
								 ! In the 3D MT case there are two polarisation.
								 ! In the 3D CSEM case there are one polarisation (for now).
								 ! In the 2D MT case there are one polarisation.
								        worker_job_task%eAll_SolnIndex=eAll_in%SolnIndex
								 	    which_per=per_index
								        do ipol1=1,nPol_MPI
                                           which_pol=ipol1
			        				       !call create_e_param_place_holder(eAll_in%solns(which_per))
				    				       !call Pack_e_para_vec(eAll_in%solns(which_per))
				    				       !call MPI_SEND(e_para_vec, Nbytes, MPI_PACKED, dest,FROM_MASTER, MPI_COMM_WORLD, ierr) 
										end do   
                                 end if	  
								 
	            			  call create_worker_job_task_place_holder
	            			  call Pack_worker_job_task
	            			  call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,dest, FROM_MASTER, MPI_COMM_WORLD, ierr)								 
							  deallocate( worker_job_package )		 
								 
	                          write(ioMPI,'(a10,a16,i5,a8,i5,a11,i5)')trim(job_name),': Send Per. # ',per_index , ' : Pol #', pol_index,' to node # ',dest
	                          if (dest .ge. number_of_workers) then
	                            goto 10
	                          end if
	                      end do    
                    end do
                    
10         continue


          !call count_number_of_meaasges_to_RECV(eAll_out) ! Because we are NOT passing FULL eAll_out
		  !answers_to_receive = nTx*nPol_MPI
		   
            answers_to_receive=0
            do itx=1,nTx
              answers_to_receive=answers_to_receive+txDict(itx)%nPol
			end do
			
       
        received_answers = 0
        do while (received_answers .lt. answers_to_receive)

            call create_worker_job_task_place_holder
            call MPI_RECV(worker_job_package, Nbytes, MPI_PACKED ,MPI_ANY_SOURCE, FROM_WORKER,MPI_COMM_WORLD,STATUS, ierr)
            call Unpack_worker_job_task

                       who=worker_job_task%taskid
                       which_per=worker_job_task%per_index
                       which_pol=worker_job_task%pol_index

            !STOP Receiving e solution from workers. Naser changes: 14.06.2021. The following 3 lines are deactivated    
            !call create_e_param_place_holder(eAll_out%solns(which_per))
            !call MPI_RECV(e_para_vec, Nbytes, MPI_PACKED, who, FROM_WORKER,MPI_COMM_WORLD, STATUS, ierr)
            !call Unpack_e_para_vec(eAll_out%solns(which_per))
			
                   write(ioMPI,'(a10,a16,i5,a8,i5,a11,i5)')trim(job_name) ,': Recieve Per # ',which_per ,' and Pol # ', which_pol ,' from ', who 
                   ! Writting into the solver's diagonestic file.
				   ! Naser and Paulo 02.10.2019
				   
				   do solver_residual_iter = 1,size(worker_job_task%solver_residual_vec)
					   if (worker_job_task%solver_residual_vec(solver_residual_iter) .gt. R_ZERO) then
						   write(ioSolverStat,'(a20,g20.7,i20,i20,i20,g20.7)')trim(job_name), worker_job_task%period, which_per, which_pol, solver_residual_iter, worker_job_task%solver_residual_vec(solver_residual_iter)
					   end if
				   end do
				   received_answers=received_answers+1
                    
                   
        ! Check if we send all transmitters and polarizations, if not then send the next transmitter to the worker who is free now....
        ! This part is very important if we have less workers than transmitters.           
!write(6,*)'Per_index ',Per_index,'pol_index ',pol_index,'received_answers',received_answers,' out of ',answers_to_receive

if (Per_index ==  nTx .and. pol_index ==nPol_MPI) goto 1500 
                      
write(6,*)'Per_index ',Per_index,'pol_index ',pol_index,'Going to send'


pol_index=pol_index+1 

 if ( pol_index .gt. nPol_MPI ) then
          Per_index=Per_index+1
          pol_index=1   
 elseif ( pol_index .le. nPol_MPI) then
          Per_index=Per_index
 end if
 
 if (Per_index .gt. nTx ) goto 1500     
        
                    !call get_nPol_MPI(eAll_out%solns(per_index)) 
					 nPol_MPI=txDict(per_index)%nPol
                    if (nPol_MPI==1)  pol_index=1
                   
                           worker_job_task%per_index= per_index
                           worker_job_task%pol_index= pol_index
						   worker_job_task%eAll_SolnIndex=eAll_out%SolnIndex
                           worker_job_task%what_to_do=trim(job_name) 
            call create_worker_job_task_place_holder
            call Pack_worker_job_task
            call MPI_SEND(worker_job_package,Nbytes, MPI_PACKED,who, FROM_MASTER, MPI_COMM_WORLD, ierr)
			deallocate( worker_job_package )					   
            if (trim(job_name).eq. 'JmultT' .or. trim(job_name).eq. 'Jmult')then
	            which_per=per_index
                !call get_nPol_MPI(eAll_out%solns(per_index))
                nPol_MPI=txDict(per_index)%nPol				
				do ipol1=1,nPol_MPI
                    which_pol=ipol1
			        !call create_e_param_place_holder(eAll_in%solns(which_per))
				    !call Pack_e_para_vec(eAll_in%solns(which_per))
				    !call MPI_SEND(e_para_vec, Nbytes, MPI_PACKED, who,FROM_MASTER, MPI_COMM_WORLD, ierr) 
				end do  
		    end if 
		     	                      
           write(ioMPI,'(a10,a16,i5,a8,i5,a11,i5)')trim(job_name),': Send Per. # ',per_index , ' : Pol #', pol_index,' to node # ',who
                           
                           

    1500         continue                                   
                   
        end do

if (associated(eAll_para_vec)) then
!deallocate(eAll_para_vec)
end if

if (associated(e_para_vec)) then
!deallocate(e_para_vec)
end if
  !call deall(e0)

end subroutine Master_job_Distribute_Taskes
!*****************************************************************************************

subroutine create_data_vec_place_holder(d)


     implicit none
     integer Nbytes1,Nbytes2,ndata,iper,ndt,sum1,sum2
     type(dataVectorMTX_t), intent(in)		:: d




              sum1=0
              sum2=0
              do iper=1,d%nTx
                    do ndt=1,d%d(iper)%ndt
                      ndata=size(d%d(iper)%data(ndt)%value)
                      CALL MPI_PACK_SIZE(ndata, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, Nbytes1,  ierr)
                      CALL MPI_PACK_SIZE(1,     MPI_LOGICAL,          MPI_COMM_WORLD, Nbytes2,  ierr)
                      sum1=sum1+Nbytes1
                      sum2=sum2+Nbytes2
                    end do
              end do
                    
        Nbytes=((2*sum1)+(2*sum2))+1

         if( associated( data_para_vec ) ) deallocate( data_para_vec )
         allocate( data_para_vec( Nbytes ) )
         
end subroutine create_data_vec_place_holder
!***************************************************************************************** 
 subroutine Pack_data_para_vec(d)
    implicit none

     type(dataVectorMTX_t), intent(in)		:: d
     integer index
     integer ndata,iper,ndt



       index=1
      do iper=1,d%nTx
        do ndt=1,d%d(iper)%ndt
             ndata=size(d%d(iper)%data(ndt)%value)            
             call MPI_Pack(d%d(iper)%data(ndt)%value(1,1),ndata, MPI_DOUBLE_PRECISION, data_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
             call MPI_Pack(d%d(iper)%data(ndt)%error(1,1),ndata, MPI_DOUBLE_PRECISION, data_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
             call MPI_Pack(d%d(iper)%data(ndt)%errorBar,1,       MPI_LOGICAL,          data_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
             call MPI_Pack(d%d(iper)%data(ndt)%allocated,1,       MPI_LOGICAL,          data_para_vec, Nbytes, index, MPI_COMM_WORLD, ierr)
       end do
     end do  


end subroutine Pack_data_para_vec
!***************************************************************************************** 
subroutine UnPack_data_para_vec(d)
    implicit none

     type(dataVectorMTX_t), intent(inout)		:: d
     integer index
     integer ndata,iper,ndt


       index=1
      do iper=1,d%nTx
        do ndt=1,d%d(iper)%ndt
             ndata=size(d%d(iper)%data(ndt)%value)            
             call MPI_Unpack(data_para_vec, Nbytes, index,d%d(iper)%data(ndt)%value(1,1),ndata, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
             call MPI_Unpack(data_para_vec, Nbytes, index,d%d(iper)%data(ndt)%error(1,1),ndata, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
             call MPI_Unpack(data_para_vec, Nbytes, index,d%d(iper)%data(ndt)%errorBar  ,1,     MPI_LOGICAL,          MPI_COMM_WORLD, ierr)
             call MPI_Unpack(data_para_vec, Nbytes, index,d%d(iper)%data(ndt)%allocated  ,1,     MPI_LOGICAL,          MPI_COMM_WORLD, ierr)
       end do
     end do  

	 deallocate( data_para_vec )							





end subroutine UnPack_data_para_vec
!***************************************************************************************** 
subroutine RECV_cUserDef(cUserDef)
implicit none
type (userdef_control),intent(inout)	:: cUserDef
 character(20)                               :: which_proc


             
 
          call create_userdef_control_place_holder
          call MPI_RECV(userdef_control_package, Nbytes, MPI_PACKED ,0, FROM_MASTER,MPI_COMM_WORLD,STATUS, ierr)
          call unpack_userdef_control (cUserDef)

	   if (taskid==1 ) then
        which_proc='Worker'
        call check_userdef_control_MPI (which_proc,cUserDef)
	   end if
	   !deallocate (userdef_control_package)
	   

end subroutine RECV_cUserDef
!*****************************************************************************************
subroutine setGrid_MPI(newgrid)

   !  Use to set and/or update the numerical grid, that is then used
   !   all computations in this module;
   !   This is not a pointer target.
   !  Might also have to run exitSolver at this point, if we are updating
   !   the grid during an inversion; that restarts the ForwardSolver module.

   type(grid_t), intent(in)     :: newgrid

   grid = newgrid

end subroutine setGrid_MPI
!*****************************************************************************************
  subroutine cleanUp_MPI()

   ! Subroutine to deallocate all memory stored in this module

   !call exitSolver(e0,e,comb)
   call deall_grid(grid)

  end subroutine cleanUp_MPI
  
subroutine destructor_MPI
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_FINALIZE(ierr)

end subroutine destructor_MPI
#endif
end module Main_MPI_files
