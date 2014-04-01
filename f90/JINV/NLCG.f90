module NLCG

!use math_constants
!use utilities
use senscomp
use dataio

#ifdef MPI
  use MPI_main
  use MPI_sub
#endif


   ! inherits datasens,  dataspace, dataFunc, SolnSpace,
   !            modelspace, soln2d

implicit none

public  :: NLCGsolver

! iteration control for the NLCG solver is initialized once
! and saved in the module to be used by most subroutines

  type  :: NLCGiterControl_t
     ! maximum number of iterations in one call to iterative solver
     integer					:: maxIter
     ! convergence criteria: return from solver if rms < rmsTol
     real (kind=prec)	:: rmsTol
     ! the condition to identify when the inversion stalls
     real (kind=prec)   :: fdiffTol
     ! initial value of lambda (will not override the NLCG input argument)
     real (kind=prec)   :: lambda
     ! exit if lambda < lambdaTol approx. 1e-4
     real (kind=prec)   :: lambdaTol
     ! set lambda_i = lambda_{i-1}/k when the inversion stalls
     real (kind=prec)   :: k
     ! the factor that ensures sufficient decrease in the line search
     real (kind=prec)   :: c
     ! restart CG every nCGmax iterations to ensure conjugacy
     integer                    :: nCGmax
     ! restart CG if orthogonality is lost (not necessarily needed)
     ! real (kind=prec)   :: delta ! 0.5
     ! the starting step for the line search
     real (kind=prec)   :: alpha_1
     ! if alpha_{i+1} < alpha_i * k_{alpha}, set alpha_{i+1} = alpha_i/2
     ! real (kind=prec)   :: alpha_k ! 0.1
     ! if alpha_{i+1} - alpha_i < tol_{alpha}, set alpha_{i+1} = alpha_i/2
     ! real (kind=prec)   :: alpha_tol ! 1.0e-2
     ! maximum initial delta mHat (overrides alpha_1)
     real (kind=prec)   :: startdm
     ! optional relaxation parameter (Renormalised Steepest Descent algorithm)
     real (kind=prec)   :: gamma
     ! model and data output file name
     character(80)              :: fname
  end type NLCGiterControl_t

  type(NLCGiterControl_t), private, save :: iterControl
  real (kind=prec), private, save        :: scaling_factor_MT,scaling_factor_CSEM,scaling_factor_DC,scaling_factor
  integer, private, save                 :: Ndata_MT,Ndata_CSEM,Ndata_DC,which_iter
  type(dataVectorMTX_t), save                   :: MT_d,CSEM_d,DC_d,data_temp
  real(kind=prec)	,pointer, dimension(:)          :: approx_res  

Contains

!**********************************************************************
   subroutine set_NLCGiterControl(iterControl)
   !  at present this is a simple routine to setup for
   !  iteration control; allowing parameters to be read from
   !   a file would be desireable.

   type(NLCGiterControl_t), intent(inout)	:: iterControl

     ! maximum number of iterations in one call to iterative solver
     iterControl%maxIter = 200
     ! convergence criteria: return from solver if rms < rmsTol
     iterControl%rmsTol  = 1.05
     ! inversion stalls when abs(rms - rmsPrev) < fdiffTol (2e-3 works well)
     iterControl%fdiffTol = 2.0e-3
     ! initial value of lambda (will not override the NLCG input argument)
     iterControl%lambda = 1.
     ! exit if lambda < lambdaTol approx. 1e-4
     iterControl%lambdaTol = 1.0e-8
     ! set lambda_i = lambda_{i-1}/k when the inversion stalls
     iterControl%k = 10.
     ! the factor that ensures sufficient decrease in the line search >=1e-4
     iterControl%c = 1.0e-4
     ! restart CG every nCGmax iterations to ensure conjugacy
     iterControl%nCGmax = 8
     ! the starting step for the line search
     iterControl%alpha_1 = 20.
     ! maximum initial delta mHat (overrides alpha_1)
     iterControl%startdm = 20.
     ! optional relaxation parameter (Renormalised Steepest Descent algorithm)
     iterControl%gamma = 0.99
     ! model and data output file name
     iterControl%fname = 'Modular'

   end subroutine set_NLCGiterControl


   ! ***************************************************************************
   ! * read_NLCGiterControl reads the inverse solver configuration from file

   subroutine read_NLCGiterControl(iterControl,rFile,fileExists)

	type(NLCGiterControl_t), intent(inout)	:: iterControl
    character(*), intent(in)		        :: rFile
	logical, intent(out), optional          :: fileExists
    integer									:: ios
	logical                             	:: exists
	character(80)							:: string

    ! Initialize inverse solver configuration

    call set_NLCGiterControl(iterControl)

    inquire(FILE=rFile,EXIST=exists)
    if (present(fileExists)) then
       fileExists = exists
    end if

    if (.not. exists) then
       return
    else
       write(*,*) 'Reading inverse configuration from file ',trim(rFile)
    end if

    open (unit=ioInvCtrl,file=rFile,status='old',iostat=ios)

    if(ios/=0) then
       write(0,*) 'Error opening file: ', rFile
    end if

    ! This is the list of options specified in the startup file

    read (ioInvCtrl,'(a36,a80)') string,iterControl%fname
    if (output_level > 2) then
       write (*,*)
       write (*,'(a36,a80)') string,iterControl%fname
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%lambda
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%lambda
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%k
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%k
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%startdm
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%startdm
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%fdiffTol
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%fdiffTol
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%rmsTol
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%rmsTol
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%lambdaTol
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%lambdaTol
    end if
    read (ioInvCtrl,'(a36,i4)') string,iterControl%maxIter
    if (output_level > 2) then
       write (*,'(a36,i4)') string,iterControl%maxIter
       write (*,*)
    end if
    close(ioInvCtrl)

   end subroutine read_NLCGiterControl


!**********************************************************************
   subroutine printf(comment,lambda,alpha,f,mNorm,rms,logfile)

   ! Compute the full penalty functional F
   ! Also output the predicted data and the EM solution
   ! that can be used for evaluating the gradient
   ! Assuming that the model norm is already scaled by Nmodel

    character(*), intent(in)               :: comment
    real(kind=prec), intent(in)  :: lambda, alpha, f, mNorm, rms
    character(*), intent(in), optional	:: logfile
    integer  :: io_unit, ios
    logical  :: opened

    if (present(logfile)) then
    	io_unit = ioLog
    	inquire(file=logfile,opened=opened)
    	if (.not. opened) then
    		open (unit=ioLog,file=logfile,status='unknown',position='append',iostat=ios)
    	end if
    else
    	io_unit = 6
    end if

	write(io_unit,'(a10)',advance='no') trim(comment)//':'
	write(io_unit,'(a3,es12.6)',advance='no') ' f=',f
	write(io_unit,'(a4,es12.6)',advance='no') ' m2=',mNorm
	write(io_unit,'(a5,f11.6)',advance='no') ' rms=',rms
	write(io_unit,'(a8,es12.6)',advance='no') ' lambda=',lambda
	write(io_unit,'(a7,es12.6)') ' alpha=',alpha

	! flush(io_unit): this has the effect of flushing the buffer
	if (present(logfile)) then
		close(io_unit)
		open (unit=ioLog,file=logfile,status='old',position='append',iostat=ios)
	end if

   end subroutine printf


!**********************************************************************
   subroutine func(lambda,d,m0,mHat,F,mNorm,dHat,eAll,RMS)

   ! Compute the full penalty functional F
   ! Also output the predicted data and the EM solution
   ! that can be used for evaluating the gradient

   real(kind=prec), intent(in)  :: lambda
   type(dataVectorMTX_t), intent(in)              :: d
   type(modelParam_t), intent(in)           :: m0
   type(modelParam_t), intent(in)           :: mHat
   real(kind=prec), intent(out) :: F, mNorm
   type(dataVectorMTX_t), optional, intent(inout)   :: dHat
   type(solnVectorMTX_t), optional, intent(inout) :: eAll
   real(kind=prec), optional, intent(out) :: RMS

   !  local variables
   type(dataVectorMTX_t)    :: res,Nres
   type(modelParam_t) :: m,JTd
   real(kind=prec) :: SS
   integer :: Ndata, Nmodel

   ! compute the smoothed model parameter vector
   call CmSqrtMult(mHat,m)

   ! overwriting input with output
   call linComb(ONE,m,ONE,m0,m)

   ! initialize dHat
   dHat = d

   !  compute predicted data for current model parameter m
   !   also sets up forward solutions for all transmitters in eAll
   !   (which is created on the fly if it doesn't exist)
#ifdef MPI
      call Master_Job_fwdPred(m,dHat,eAll)
#else
      call fwdPred(m,dHat,eAll)
#endif

!	call write_Z_ascii(fidWrite,cfile,nPer,periods,modes, &
!			nSites,sites,allData)

   ! initialize res
   res = d

   ! compute residual: res = d-dHat
   call linComb(ONE,d,MinusONE,dHat,res)
    

   ! normalize residuals, compute sum of squares
   call CdInvMult(res,Nres) 
	 
   SS = dotProd(res,Nres)
   Ndata = countData(res)

   ! compute the model norm
   mNorm = dotProd(mHat,mHat)
   Nmodel = countModelParam(mHat)

   ! penalty functional = sum of squares + scaled model norm
   F = SS/Ndata + (lambda * mNorm/Nmodel)

   ! scale mNorm for output
   mNorm = mNorm/Nmodel

   ! if required, compute the Root Mean Squared misfit
   if (present(RMS)) then
   	RMS = sqrt(SS/Ndata)
   end if

   call deall_dataVectorMTX(res)
   call deall_dataVectorMTX(Nres)
   call deall_modelParam(m)
   call deall_modelParam(JTd)

   end subroutine func
   !***********************************************************************
  subroutine normalize_with_dataVecMTX(d_in_out,d,N)
  
     type(dataVectorMTX_t), intent(in)              :: d
     type(dataVectorMTX_t), intent(inout)           :: d_in_out
     integer, optional, intent(in)               :: N
     integer                                     :: i,j,k,iDt
     

 	            do i=1,d%nTx
	             do iDt=1,d%d(i)%nDt
	              do j=1,d%d(i)%data(iDt)%nSite
	               do k=1,d%d(i)%data(iDt)%nComp
	                      d_in_out%d(i)%data(iDt)%value(k,j)=  (d_in_out%d(i)%data(iDt)%value(k,j)/d%d(i)%data(iDt)%error(k,j)**N)
	                      d_in_out%d(i)%data(iDt)%errorBar=.true.
	                end do      
	              end do
	            end do                                       
	        end do    
     

     
  
  end subroutine normalize_with_dataVecMTX
!**********************************************************************
   subroutine gradient(lambda,d,m0,mHat,grad,dHat,eAll)

   !  Computes the gradient of the penalty functional,
   !  using EM solution (eAll) and the predicted data (dHat)
   !  Here, mHat denotes the non-regularized model parameter that
   !  is normally referred to as \tilde{m} = C_m^{-1/2}(m - m_0),
   !  and the gradient is computed with respect to \tilde{m}.
   !  Before calling this routine, the forward solver must be run:
   !  call CmSqrtMult(mHat,m)
   !  call linComb(ONE,m,ONE,m0,m)
   !  call fwdPred(m,dHat,eAll)

   real(kind=prec), intent(in)  :: lambda
   type(dataVectorMTX_t), intent(inout)              :: d
   type(modelParam_t), intent(in)           :: m0
   type(modelParam_t), intent(inout)           :: mHat
   type(modelParam_t), intent(inout)          :: grad
   type(dataVectorMTX_t), intent(inout)              :: dHat
   type(solnVectorMTX_t), intent(inout)            :: eAll

   !  local variables
   real(kind=prec)       :: Ndata,Nmodel,mNorm,mHatNorm,resNorm,CmJTd_norm_before,CmJTd_norm_after,CmJTd_DC_norm,CmJTd_CSEM_norm,CmJTd_MT_norm,CmJTd_DC_CSEM_norm,CmJTd_MT_CSEM_norm,CmJTd_DC_CSEM_MT_norm,vAir
   real (kind=prec)      :: RMS_DC,RMS_MT,RMS_CSEM,RMS_total,lambda1,sum1,sum0,RMS,F,temp_norm,term,b1,b2,beta_MT,beta_CSEM 
   type(dataVectorMTX_t)    :: res,Nres,DC_res,MT_res,CSEM_res,unite_d
   type(modelParam_t) :: m,m_sub,mHat_sub,CmHat_sub,CmJTd_MT_N,CmJTd_CSEM_N,JTd,CmJTd,JTd_DC,JTd_MT,JTd_CSEM,CmJTd_DC,CmJTd_MT,CmJTd_CSEM,CmJTd_MT_CSEM,CmJTd_DC_CSEM,CmJTd_DC_CSEM_MT,CmHat,Cm_s_hat,q_temp
   integer            :: iTx,Tx_MT_counter,Tx_DC_counter,Tx_CSEM_counter,counter_DC,counter_MT,counter_CSEM,idt,nTx,k,j,i,ios
   character(100)       :: mFile,  dataFile
   type(modelParam_t),pointer, dimension(:) :: s_hat
   character(3)         :: iterChar
   character(10)        :: set_reset
   type(rscalar)        ::       model_MT,model_CSEM,model_MIX
   character(len=80)    ::       paramtype
   ! integer :: j, Ny, NzEarth
      type(dataVectorMTX_t)    :: MT_d_temp,CSEM_d_temp,DC_d_temp,MT_d_N,CSEM_d_N,DC_d_N
    
		real(kind=prec)	,pointer, dimension(:,:)        :: JJT,JJT_temp
	    real(kind=prec)	,pointer, dimension(:)          :: b,AP,b_Ax
	    Integer                                         ::INFO,datatype
        real(kind=prec)                                 :: a_sub(2,2),b_sub(2),x_sub(2),a_sub_INV(2,2),det,w1,w2
	

	
       write(iterChar,'(i3.3)') which_iter
   
     if (which_iter == 0) then
       open (unit=ioWeighting,file=trim(iterControl%fname)//'_Weighting.log',status='unknown',iostat=ios)
       write(ioWeighting,'(a14,12a30)') 'Iteration No. ', 'CmJTd_norm_before','CmJTd_norm_after','CmJTd_MT_norm','CmJTd_CSEM_norm', 'scaling_factor_MT','scaling_factor_CSEM','RMS_MT','RMS_CSEM','RMS_total','CmJTd_MT_norm/CmJTd_CSEM_norm','CmJTd_CSEM_norm/CmJTd_MT_norm','scaling_factor'
	   data_temp=d
	 end if
	 

        
		
        nTx=d%nTx
	    JTd_MT 	       = m0
	    JTd_CSEM 	   = m0
	    JTd_DC  	   = m0
        call zero(JTd_MT)
        call zero(JTd_CSEM)	
        call zero(JTd_DC)	
		
	    CmJTd_MT 	   = m0
	    CmJTd_CSEM 	   = m0
	    CmJTd_DC 	   = m0
        CmJTd_MT_N     = m0
        CmJTd_CSEM_N     = m0
		CmHat          = m0
        mHat_sub      =m0
        CmHat_sub      =m0
        m_sub      =m0
		q_temp         = m0

        call zero(CmJTd_MT)
        call zero(CmJTd_CSEM)
        call zero(CmJTd_DC)
        call zero(CmJTd_MT_N)
        call zero(CmJTd_CSEM_N)
        call zero(m_sub)
        call zero(CmHat_sub)
          call zero(mHat_sub)      
		CmJTd_MT_CSEM    = m0
		CmJTd_DC_CSEM    = m0
		CmJTd_DC_CSEM_MT  =m0
		
		call zero(CmJTd_MT_CSEM)
		call zero(CmJTd_DC_CSEM)
		call zero(CmJTd_DC_CSEM_MT)
		
   ! compute the smoothed model parameter vector
   call CmSqrtMult(mHat,m)
    call CmSqrtMult(mHat,CmHat)

   ! overwriting the input with output
   call linComb(ONE,m,ONE,m0,m)

                          

						
   
   ! initialize res
   res = d

   ! compute residual: res = (d-dHat)/Ndata
   call linComb(ONE,d,MinusONE,dHat,res)
   
	!set_reset='reset'
	!call set_reset_data_error(res,set_reset) 
   
! Spilt res before and after normalizing to compute RMS
       call split_dataVectorMTX(res)
if (Ndata_CSEM .ne. 0 .and. Ndata_DC .ne. 0 .and. Ndata_MT .ne. 0) then
	                   MT_d_temp= MT_d
					   DC_d_temp= DC_d
					   CSEM_d_temp= CSEM_d
					   call deall_dataVectorMTX(MT_d)
					   call deall_dataVectorMTX(DC_d)
					   call deall_dataVectorMTX(CSEM_d)
				   !Normalize res	   
					  call CdInvMult(res,Nres)
				   
					  call split_dataVectorMTX(Nres)  
					  DC_d_N= DC_d
					  MT_d_N= MT_d
					  CSEM_d_N=CSEM_d
						 call deall_dataVectorMTX(DC_d)
						 call deall_dataVectorMTX(MT_d)
						 call deall_dataVectorMTX(CSEM_d) 
						 
				!Compute RMS value for each data set
						 RMS_DC   = sqrt(dotProd(DC_d_temp,DC_d_N)/Ndata_DC)
                         RMS_MT   = sqrt(dotProd(MT_d_temp,MT_d_N)/Ndata_MT)						 
						 RMS_CSEM = sqrt(dotProd(CSEM_d_temp,CSEM_d_N)/Ndata_CSEM)
						 RMS_total= sqrt(dotProd(res,Nres)/(Ndata_CSEM+Ndata_DC+Ndata_MT)) 
						
						 
						 res=Nres


					 

                        !set_reset='set'					   
                        !call set_reset_data_error(res,res,set_reset) 						 
						!call set_reset_data_error(d,data_temp,set_reset)  
						 
						 ! call CdInvMult(res)
						 
				   !initialize s_hat to used later in spiltting JTd for each Tx
					 allocate(s_hat(nTx))
					 do iTx=1,nTx
						s_hat(iTx)=m
						call zero(s_hat(iTx))
					  end do
					  
				 ! multiply the normalized residual (res) by J^T   
#ifdef MPI
						call Master_job_JmultT(m,res,JTd,eAll,s_hat) ! s_hat contains the JTd for each Tx. Thus, s_hat is a vector of length nTx of model parameters
#else
						call JmultT(m,res,JTd,eAll)
#endif

                ! Compute and write the norm of each  JTd correspond to each Tx
				open(550,file='Norms_of_s_hat_'//iterChar//'.dat',status='unknown',iostat=ios)
				 do iTx=1,nTx
				 temp_norm = sqrt(dotProd(s_hat(iTx),s_hat(iTx)))
				  if (txDict(iTx)%Tx_type=='DC') then
				       write(550,*)iTx,temp_norm,txDict(iTx)%period,'DC'
				  elseif (txDict(iTx)%Tx_type=='CSEM') then
				       write(550,*)iTx,temp_norm,'CSEM'
				  elseif (txDict(iTx)%Tx_type=='MT') then
				       write(550,*)iTx,temp_norm,'MT'
				   end if					   
				 
				 		write(iterChar,'(i3.3)') iTx
					  if (txDict(iTx)%Tx_type=='DC') then
					  	  mFile = 'JTd_DC_'//iterChar//'.rho'
						  call write_modelParam(s_hat(iTx),trim(mFile))
						  call linComb_modelParam(ONE,JTd_DC,ONE,s_hat(iTx),JTd_DC) ! get JTd_DC
					  elseif (txDict(iTx)%Tx_type=='CSEM') then
					  	  mFile = 'JTd_CSEM_'//iterChar//'.rho'
						  call write_modelParam(s_hat(iTx),trim(mFile))
						  call linComb_modelParam(ONE,JTd_CSEM,ONE,s_hat(iTx),JTd_CSEM) ! get JTd_CSEM
					  elseif (txDict(iTx)%Tx_type=='MT') then
					  	  mFile = 'JTd_MT_'//iterChar//'.rho'
						  call write_modelParam(s_hat(iTx),trim(mFile))
						  call linComb_modelParam(ONE,JTd_MT,ONE,s_hat(iTx),JTd_MT) ! get JTd_MT
					  end if
				 end do	  
				 close(550)
				 ! Up here we have three data gradient: 
				 ! 1) JTd for ALL nTx
				 ! 2) CmJTd_MT for ONLY MT transmitters (Periods)
				 ! 3) CmJTd_CSEM for ONLY CSEM transmitters and their periods

				 ! Smooth each of the data gradient
				   call CmSqrtMult(JTd,CmJTd)
				   !datatype=2
				   !call create_CmSqrt(JTd_DC,datatype)
				   call CmSqrtMult(JTd_DC,CmJTd_DC)
				   !datatype=2
				  ! call create_CmSqrt(JTd_DC,datatype)
				   call CmSqrtMult(JTd_MT,CmJTd_MT)
				   !call create_CmSqrt(JTd_DC,datatype)
				   call CmSqrtMult(JTd_CSEM,CmJTd_CSEM)



					
				  
				   !call split_dataVectorMTX(res)

				   CmJTd_norm_before = sqrt(dotProd(CmJTd,CmJTd))
				   CmJTd_DC_norm     = sqrt(dotProd(CmJTd_DC,CmJTd_DC))
				   CmJTd_MT_norm     = sqrt(dotProd(CmJTd_MT,CmJTd_MT))
				   CmJTd_CSEM_norm   = sqrt(dotProd(CmJTd_CSEM,CmJTd_CSEM))

				   



					 call linComb_modelParam(ONE,CmJTd_DC,ONE,CmJTd_CSEM,CmJTd_DC_CSEM)
					 call linComb_modelParam(ONE,CmJTd_DC_CSEM,ONE,CmJTd_MT,CmJTd_DC_CSEM_MT)
					 
					 CmJTd_DC_CSEM_MT_norm = sqrt(dotProd(CmJTd_DC_CSEM_MT,CmJTd_DC_CSEM_MT))
						 

						 scaling_factor=1.0!/3.0
						 

						 
						! if (RMS_DC .le. 2.0) then
						!  scaling_factor= 0.5
						! end if		
						! 
						!
						!if (RMS_DC .le. 1.2) then
						!  scaling_factor= 0.1
						!end if
						!
						!if (RMS_CSEM .le. 1.05) then
						!  scaling_factor= 0.75
						!end if	
						
			  
					   
                        !set_reset='set'					   
                        !call set_reset_data_error(d,data_temp,set_reset) 
			   


					 
					  scaling_factor_DC =1.0 !(scaling_factor*CmJTd_norm_before)/CmJTd_DC_norm    
                      scaling_factor_CSEM = 1.0 !((scaling_factor)*CmJTd_norm_before)/CmJTd_CSEM_norm    				  
                      scaling_factor_MT =1.0 !(scaling_factor*CmJTd_norm_before)/CmJTd_MT_norm 
					  
					  ! Scale each JTd (CmJTd_MT and CmJTd_CSEM) and sum to get CmJTd_MT_CSEM
					  call linComb_modelParam(scaling_factor_DC,CmJTd_DC,scaling_factor_CSEM,CmJTd_CSEM,CmJTd_DC_CSEM)
                      call linComb_modelParam(ONE,CmJTd_DC_CSEM,scaling_factor_MT,CmJTd_MT,CmJTd_DC_CSEM_MT)
					  
					  CmJTd_DC_CSEM_MT_norm = sqrt(dotProd(CmJTd_DC_CSEM_MT,CmJTd_DC_CSEM_MT))
					  
					  !Set CmJTd = CmJTd_DC_CSEM
					   CmJTd=CmJTd_DC_CSEM_MT
					   CmJTd_norm_after = sqrt(dotProd(CmJTd,CmJTd))  
					   


					write(ioWeighting,'(i5,13(4x,es12.6))') which_iter, CmJTd_norm_before,CmJTd_norm_after,CmJTd_DC_norm,CmJTd_CSEM_norm,CmJTd_MT_norm, scaling_factor_DC,scaling_factor_CSEM,scaling_factor_MT,RMS_DC,RMS_CSEM,RMS_MT,RMS_total,scaling_factor


				! initialize grad
					grad = m	
				! compute the number of data and model parameters for scaling	
					 Ndata = countData(res)
					 Nmodel = countModelParam(mHat)
					 
					 call linComb(MinusTWO/Ndata,CmJTd,TWO*lambda/Nmodel,mHat,grad)
					
					
						set_reset='set'
						!call set_reset_data_error(d,set_reset) 

					  do iTx=1,nTx
						call deall_modelParam(s_hat(iTx))
					  end do
					  deallocate(s_hat)   	   
elseif (Ndata_CSEM .ne. 0 .and. Ndata_DC .ne. 0) then
	
					   DC_d_temp= DC_d
					   CSEM_d_temp= CSEM_d
					   call deall_dataVectorMTX(DC_d)
					   call deall_dataVectorMTX(CSEM_d)
				   !Normalize res	   
					  call CdInvMult(res,Nres)
				   
					  call split_dataVectorMTX(Nres)  
					  DC_d_N= DC_d
					  CSEM_d_N=CSEM_d
						 call deall_dataVectorMTX(DC_d)
						 call deall_dataVectorMTX(CSEM_d) 
						 
				!Compute RMS value for each data set
						 RMS_DC   = sqrt(dotProd(DC_d_temp,DC_d_N)/Ndata_DC)		 
						 RMS_CSEM = sqrt(dotProd(CSEM_d_temp,CSEM_d_N)/Ndata_CSEM)
						 RMS_total= sqrt(dotProd(res,Nres)/(Ndata_CSEM+Ndata_DC)) 
						
						 
						 res=Nres


					 

                        !set_reset='set'					   
                        !call set_reset_data_error(res,res,set_reset) 						 
						!call set_reset_data_error(d,data_temp,set_reset)  
						 
						 ! call CdInvMult(res)
						 
				   !initialize s_hat to used later in spiltting JTd for each Tx
					 allocate(s_hat(nTx))
					 do iTx=1,nTx
						s_hat(iTx)=m
						call zero(s_hat(iTx))
					  end do
					  
				 ! multiply the normalized residual (res) by J^T   
#ifdef MPI
						call Master_job_JmultT(m,res,JTd,eAll,s_hat) ! s_hat contains the JTd for each Tx. Thus, s_hat is a vector of length nTx of model parameters
#else
						call JmultT(m,res,JTd,eAll)
#endif

                ! Compute and write the norm of each  JTd correspond to each Tx
				open(550,file='Norms_of_s_hat_'//iterChar//'.dat',status='unknown',iostat=ios)
				 do iTx=1,nTx
				 temp_norm = sqrt(dotProd(s_hat(iTx),s_hat(iTx)))
				  if (txDict(iTx)%Tx_type=='DC') then
				       write(550,*)iTx,temp_norm,txDict(iTx)%period,'DC'
				  elseif (txDict(iTx)%Tx_type=='CSEM') then
				       write(550,*)iTx,temp_norm,'CSEM'
                  end if					   
				 
				 		write(iterChar,'(i3.3)') iTx
					  if (txDict(iTx)%Tx_type=='DC') then
					  	  mFile = 'JTd_DC_'//iterChar//'.rho'
						  call write_modelParam(s_hat(iTx),trim(mFile))
						  call linComb_modelParam(ONE,JTd_DC,ONE,s_hat(iTx),JTd_DC) ! get JTd_MT
					  elseif (txDict(iTx)%Tx_type=='CSEM') then
					  	  mFile = 'JTd_CSEM_'//iterChar//'.rho'
						  call write_modelParam(s_hat(iTx),trim(mFile))
						 call linComb_modelParam(ONE,JTd_CSEM,ONE,s_hat(iTx),JTd_CSEM) ! get JTd_CSEM
					  end if
				 end do	  
				 close(550)
				 ! Up here we have three data gradient: 
				 ! 1) JTd for ALL nTx
				 ! 2) CmJTd_MT for ONLY MT transmitters (Periods)
				 ! 3) CmJTd_CSEM for ONLY CSEM transmitters and their periods

				 ! Smooth each of the data gradient
				   call CmSqrtMult(JTd,CmJTd)
				   call CmSqrtMult(JTd_DC,CmJTd_DC)
				   call CmSqrtMult(JTd_CSEM,CmJTd_CSEM)



					
				  
				   !call split_dataVectorMTX(res)

				   CmJTd_norm_before = sqrt(dotProd(CmJTd,CmJTd))
				   CmJTd_DC_norm     = sqrt(dotProd(CmJTd_DC,CmJTd_DC))
				   CmJTd_CSEM_norm   = sqrt(dotProd(CmJTd_CSEM,CmJTd_CSEM))

				   



					 call linComb_modelParam(ONE,CmJTd_DC,ONE,CmJTd_CSEM,CmJTd_DC_CSEM)
					 CmJTd_DC_CSEM_norm = sqrt(dotProd(CmJTd_DC_CSEM,CmJTd_DC_CSEM))
						 

						 scaling_factor=0.25
						 

						 
						 if (RMS_DC .le. 2.0) then
						  scaling_factor= 0.5
						 end if		
						 
						
						if (RMS_DC .le. 1.2) then
						  scaling_factor= 0.1
						end if
						
						if (RMS_CSEM .le. 1.05) then
						  scaling_factor= 0.75
						end if	
						
			  
					   
                        !set_reset='set'					   
                        !call set_reset_data_error(d,data_temp,set_reset) 
			   


					 
					  scaling_factor_DC =(scaling_factor*CmJTd_norm_before)/CmJTd_DC_norm    
                      scaling_factor_CSEM = ((1-scaling_factor)*CmJTd_norm_before)/CmJTd_CSEM_norm    				  

					  ! Scale each JTd (CmJTd_MT and CmJTd_CSEM) and sum to get CmJTd_MT_CSEM
					  call linComb_modelParam(scaling_factor_DC,CmJTd_DC,scaling_factor_CSEM,CmJTd_CSEM,CmJTd_DC_CSEM)

					  
					  CmJTd_DC_CSEM_norm = sqrt(dotProd(CmJTd_DC_CSEM,CmJTd_DC_CSEM))
					  
					  !Set CmJTd = CmJTd_DC_CSEM
					   CmJTd=CmJTd_DC_CSEM
					   CmJTd_norm_after = sqrt(dotProd(CmJTd,CmJTd))  
					   


						write(ioWeighting,'(i5,12(4x,es12.6))') which_iter, CmJTd_norm_before,CmJTd_norm_after,CmJTd_DC_norm,CmJTd_CSEM_norm, scaling_factor_DC,scaling_factor_CSEM,RMS_DC,RMS_CSEM,RMS_total,CmJTd_DC_norm/CmJTd_CSEM_norm,CmJTd_CSEM_norm/CmJTd_DC_norm,scaling_factor


				! initialize grad
					grad = m	
				! compute the number of data and model parameters for scaling	
					 Ndata = countData(res)
					 Nmodel = countModelParam(mHat)
					 
					 call linComb(MinusTWO/Ndata,CmJTd,TWO*lambda/Nmodel,mHat,grad)
					
					
						set_reset='set'
						!call set_reset_data_error(d,set_reset) 

					  do iTx=1,nTx
						call deall_modelParam(s_hat(iTx))
					  end do
					  deallocate(s_hat)   
elseif (Ndata_CSEM .ne. 0 .and. Ndata_MT .ne. 0) then
	
					   MT_d_temp= MT_d
					   CSEM_d_temp= CSEM_d
					   call deall_dataVectorMTX(MT_d)
					   call deall_dataVectorMTX(CSEM_d)
				   !Normalize res	   
					  call CdInvMult(res,Nres)
				   
					  call split_dataVectorMTX(Nres)  
					  MT_d_N= MT_d
					  CSEM_d_N=CSEM_d
						 call deall_dataVectorMTX(MT_d)
						 call deall_dataVectorMTX(CSEM_d) 
						 
				!Compute RMS value for each data set
						 RMS_MT   = sqrt(dotProd(MT_d_temp,MT_d_N)/Ndata_MT)		 
						 RMS_CSEM = sqrt(dotProd(CSEM_d_temp,CSEM_d_N)/Ndata_CSEM)
						 RMS_total= sqrt(dotProd(res,Nres)/(Ndata_CSEM+Ndata_MT)) 
						
						 
						 res=Nres


					  

                        !set_reset='set'					   
                        !call set_reset_data_error(res,res,set_reset) 						 
						!call set_reset_data_error(d,data_temp,set_reset)  
						 
						 ! call CdInvMult(res)
						 
				   !initialize s_hat to used later in spiltting JTd for each Tx
					 allocate(s_hat(nTx))
					 do iTx=1,nTx
						s_hat(iTx)=m
						call zero(s_hat(iTx))
					  end do
					  
				 ! multiply the normalized residual (res) by J^T   
#ifdef MPI
						call Master_job_JmultT(m,res,JTd,eAll,s_hat) ! s_hat contains the JTd for each Tx. Thus, s_hat is a vector of length nTx of model parameters
#else
						call JmultT(m,res,JTd,eAll,s_hat)
#endif

                ! Compute and write the norm of each  JTd correspond to each Tx
				open(550,file='Norms_of_s_hat_'//iterChar//'.dat',status='unknown',iostat=ios)
				 do iTx=1,nTx
				 temp_norm = sqrt(dotProd(s_hat(iTx),s_hat(iTx)))
				  if (txDict(iTx)%Tx_type=='MT') then
				       write(550,*)iTx,temp_norm,txDict(iTx)%period,'MT'
				  elseif (txDict(iTx)%Tx_type=='CSEM') then
				       write(550,*)iTx,temp_norm,'CSEM'
                  end if					   
				 
				 		write(iterChar,'(i3.3)') iTx
					  if (txDict(iTx)%Tx_type=='MT') then
					  	  !mFile = 'JTd_MT_'//iterChar//'.rho'
						  !call write_modelParam(s_hat(iTx),trim(mFile))
						 call linComb_modelParam(ONE,JTd_MT,ONE,s_hat(iTx),JTd_MT) ! get JTd_MT
					  elseif (txDict(iTx)%Tx_type=='CSEM') then
					  	  !mFile = 'JTd_CSEM_'//iterChar//'.rho'
						  !call write_modelParam(s_hat(iTx),trim(mFile))
						 call linComb_modelParam(ONE,JTd_CSEM,ONE,s_hat(iTx),JTd_CSEM) ! get JTd_CSEM
					  end if
				 end do	  
				 close(550)
				 ! Up here we have three data gradient: 
				 ! 1) JTd for ALL nTx
				 ! 2) CmJTd_MT for ONLY MT transmitters (Periods)
				 ! 3) CmJTd_CSEM for ONLY CSEM transmitters and their periods

				 ! Smooth each of the data gradient
				   call CmSqrtMult(JTd,CmJTd)
				   call CmSqrtMult(JTd_MT,CmJTd_MT)
				   call CmSqrtMult(JTd_CSEM,CmJTd_CSEM)



					
				  
				   !call split_dataVectorMTX(res)

				   CmJTd_norm_before = sqrt(dotProd(CmJTd,CmJTd))
				   CmJTd_MT_norm     = sqrt(dotProd(CmJTd_MT,CmJTd_MT))
				   CmJTd_CSEM_norm   = sqrt(dotProd(CmJTd_CSEM,CmJTd_CSEM))

				   



					 call linComb_modelParam(ONE,CmJTd_MT,ONE,CmJTd_CSEM,CmJTd_MT_CSEM)
					 CmJTd_MT_CSEM_norm = sqrt(dotProd(CmJTd_MT_CSEM,CmJTd_MT_CSEM))
						 

						 scaling_factor=0.5
						 

						 
						! if (RMS_MT .le. 2.0) then
						!  scaling_factor= 0.5
						! end if		
						! 
						!
						!if (RMS_MT .le. 1.2) then
						!  scaling_factor= 0.1
						!end if
						!
						if (RMS_MT .le. 1.05) then
						  scaling_factor= 0.25
						end if	
						
			  
					   
                        !set_reset='set'					   
                        !call set_reset_data_error(d,data_temp,set_reset) 
			   
                     res=d
                     call linComb(ONE,d,MinusONE,dHat,res)
                     call normalize_dataVectorMTX(res,1)
                     call split_dataVectorMTX(res)  
                     MT_d_N= MT_d
                     CSEM_d_N=CSEM_d

                     beta_MT=sqrt(dotProd(MT_d_N,MT_d_N))   
                     beta_CSEM=sqrt(dotProd(CSEM_d_N,CSEM_d_N))
                     
                     !call scMult(1.0/beta_MT,CmJTd_MT,CmJTd_MT_N) 
                     CmJTd_MT_N=CmJTd_MT
                     call scMult(1.0/beta_MT,MT_d_N,MT_d_temp)
 
                     !call scMult(1.0/beta_CSEM,CmJTd_CSEM,CmJTd_CSEM_N) 
                     CmJTd_CSEM_N=CmJTd_CSEM
                     call scMult(1.0/beta_CSEM,CSEM_d_N,CSEM_d_temp)
                     
                     a_sub(1,1)=dotProd(CmJTd_MT_N,CmJTd_MT_N)
                     a_sub(1,2)=dotProd(CmJTd_CSEM_N,CmJTd_MT_N)  
                     a_sub(2,1)=dotProd(CmJTd_MT_N,CmJTd_CSEM_N) 
                     a_sub(2,2)=dotProd(CmJTd_CSEM_N,CmJTd_CSEM_N)
                     
                     b_sub(1)=beta_MT**2 
                     b_sub(2)=beta_CSEM**2 
                     
                      det = (a_sub(1,1)*a_sub(2,2))-(a_sub(1,2)*a_sub(2,1))
			          a_sub_INV(1,1) =  a_sub(2,2)/det
			          a_sub_INV(2,2) =  a_sub(1,1)/det
			          a_sub_INV(1,2) = -a_sub(1,2)/det
			          a_sub_INV(2,1) = -a_sub(2,1)/det
                     
                      x_sub(1)=a_sub_INV(1,1)*b_sub(1)+a_sub_INV(1,2)*b_sub(2)
                      x_sub(2)=a_sub_INV(2,1)*b_sub(1)+a_sub_INV(2,2)*b_sub(2)
                    
                      call linComb_modelParam(x_sub(1),CmJTd_MT_N,x_sub(2),CmJTd_CSEM_N,mHat_sub)       
  	      call CmSqrtMult(mHat_sub,CmHat) 
	      mHat_sub=CmHat
	        	     
	     call linComb_modelParam(ONE,m0,ONE,mHat_sub,m_sub)
                          write(iterChar,'(i3.3)') which_iter
  					  	  mFile = 'Model_Project_'//iterChar//'.rho'
						  call write_modelParam(m_sub,trim(mFile))       
         if (RMS_MT .le. 1.05) then
					  scaling_factor = x_sub(1)/(x_sub(1)+x_sub(2))  
                      w2=scaling_factor
                      w1=1.0-scaling_factor
         else
                      scaling_factor = x_sub(1)/(x_sub(1)+x_sub(2))
                      w1=scaling_factor
                      w2=1.0-scaling_factor
         end if 
                      scaling_factor = x_sub(2)/(x_sub(1)+x_sub(2))
                      w1=scaling_factor
                      w2=1.0-scaling_factor
                      
                       !scaling_factor = x_sub(1)/(x_sub(1)+x_sub(2))
                      ! (scaling_factor*CmJTd_norm_before)/CmJTd_MT_norm    
                      !scaling_factor_CSEM =x_sub(2)/(x_sub(1)+x_sub(2)) !x_sub(2)/beta_CSEM ! ((1-scaling_factor)*CmJTd_norm_before)/CmJTd_CSEM_norm    				  

                      
                      scaling_factor_MT =  (w1)*CmJTd_norm_before/CmJTd_MT_norm    
                      scaling_factor_CSEM =(w2)*CmJTd_norm_before/CmJTd_CSEM_norm    				  

                      
					  ! Scale each JTd (CmJTd_MT and CmJTd_CSEM) and sum to get CmJTd_MT_CSEM
					  call linComb_modelParam(scaling_factor_MT,CmJTd_MT,scaling_factor_CSEM,CmJTd_CSEM,CmJTd_MT_CSEM)

					  
					  CmJTd_MT_CSEM_norm = sqrt(dotProd(CmJTd_MT_CSEM,CmJTd_MT_CSEM))
					  
					  !Set CmJTd = CmJTd_MT_CSEM
					   CmJTd=CmJTd_MT_CSEM
					   CmJTd_norm_after = sqrt(dotProd(CmJTd,CmJTd))  
					   


						write(ioWeighting,'(i5,12(4x,es12.6))') which_iter, CmJTd_norm_before,CmJTd_norm_after,CmJTd_MT_norm,CmJTd_CSEM_norm, scaling_factor_MT,scaling_factor_CSEM,RMS_MT,RMS_CSEM,RMS_total,CmJTd_MT_norm/CmJTd_CSEM_norm,CmJTd_CSEM_norm/CmJTd_MT_norm,scaling_factor
                        !write(ioWeighting,*)a_sub(1,1),a_sub(1,2)
                        !write(ioWeighting,*)a_sub(2,1),a_sub(2,2)
                        !write(ioWeighting,*)b_sub(1),b_sub(2)
                        !write(ioWeighting,*)b_sub(1)/a_sub(1,1)

                        
           
                     
                        
                        
                        
                        
                        
				! initialize grad
					grad = m	
				! compute the number of data and model parameters for scaling	
					 Ndata = countData(res)
					 Nmodel = countModelParam(mHat)
					 
					 call linComb(MinusTWO/Ndata,CmJTd,TWO*lambda/Nmodel,mHat,grad)
					
					
						set_reset='set'
						!call set_reset_data_error(d,set_reset) 

					  do iTx=1,nTx
						call deall_modelParam(s_hat(iTx))
					  end do
					  deallocate(s_hat)
else

                 call CdInvMult(res)
#ifdef MPI
						call Master_job_JmultT(m,res,JTd,eAll) 
#else
						call JmultT(m,res,JTd,eAll)
#endif
                         
					     mFile = 'JTd_before_smoothing'//iterChar//'.rho'
						 call write_modelParam(JTd,trim(mFile))
                        
    					 call CmSqrtMult(JTd,CmJTd)
				         
						 mFile = 'JTd_after_smoothing'//iterChar//'.rho'
						 call write_modelParam(CmJTd,trim(mFile))
						 
				 ! initialize grad
					grad = m	
				! compute the number of data and model parameters for scaling	
					 Ndata = countData(res)
					 Nmodel = countModelParam(mHat)
					 call linComb(MinusTWO/Ndata,CmJTd,TWO*lambda/Nmodel,mHat,grad)	  
	  
end if
	  
   call deall_dataVectorMTX(res)
   call deall_modelParam(m)
   call deall_modelParam(JTd)
   call deall_modelParam(CmJTd)
   !call deall(eAll)
   end subroutine gradient
  !**********************************************************************
subroutine set_unite_data(d_in) 

      type(dataVectorMTX_t), intent(inout)           :: d_in
	  
	  !LOCAL
	  Integer                     :: iTx,iDt,j,k,counter
	  real(kind=prec)             :: term_MT,term_CSEM,sum1
	  
		allocate(approx_res(d_in%nTx))
		
   	       do iTx=1,d_in%nTx
			 if (txDict(iTx)%Tx_type=='MT') then
			  sum1=0.0
			  counter=0
					do iDt=1,d_in%d(iTx)%nDt
					 do j=1,d_in%d(iTx)%data(iDt)%nSite
					   do k=1,d_in%d(iTx)%data(iDt)%nComp
					         sum1=sum1+d_in%d(iTx)%data(iDt)%value(k,j)
							 counter=counter+1
							  d_in%d(iTx)%data(iDt)%value(k,j)= ONE
					  end do
					end do
				   end do 
			 elseif (txDict(iTx)%Tx_type=='CSEM') then
			  sum1=0.0
			  counter=0			 
					do iDt=1,d_in%d(iTx)%nDt
					 do j=1,d_in%d(iTx)%data(iDt)%nSite
					   do k=1,d_in%d(iTx)%data(iDt)%nComp
					         sum1=sum1+d_in%d(iTx)%data(iDt)%value(k,j)
							 counter=counter+1					   
							  d_in%d(iTx)%data(iDt)%value(k,j)= ONE
					  end do
					end do
				   end do 
		     end if		 
          approx_res(iTx)=sum1			 
         end do 
		
end subroutine set_unite_data		
   !**********************************************************************
subroutine set_reset_data_error(d_in,d_orgi,set_reset) 

      type(dataVectorMTX_t), intent(inout)           :: d_in
	  type(dataVectorMTX_t), intent(in)              :: d_orgi
	  character(10),intent(in)                    ::  set_reset
	  
	  !LOCAL
	  Integer                     :: iTx,iDt,j,k
	  real(kind=prec)             :: term_MT,term_CSEM
	  
	  if (trim(set_reset)=='set') then
	    !term_MT   = sqrt(scaling_factor_MT)
	    !term_CSEM = sqrt(scaling_factor_CSEM)
	    term_MT   = (scaling_factor_MT)
	    term_CSEM = (scaling_factor_CSEM)		
	 elseif (trim(set_reset)=='reset') then
	    !term_MT   = ONE/sqrt(scaling_factor_MT)
	    !term_CSEM = ONE/sqrt(scaling_factor_CSEM)	
		term_MT   = ONE
	    term_CSEM = ONE
	  end if	
	    term_MT   = (scaling_factor_MT)
	    term_CSEM = (scaling_factor_CSEM)
		
   	       do iTx=1,d_in%nTx
			 if (txDict(iTx)%Tx_type=='MT') then
					do iDt=1,d_in%d(iTx)%nDt
					 do j=1,d_in%d(iTx)%data(iDt)%nSite
					   do k=1,d_in%d(iTx)%data(iDt)%nComp
							  d_in%d(iTx)%data(iDt)%error(k,j)= d_orgi%d(iTx)%data(iDt)%error(k,j)*term_MT
					  end do
					end do
				   end do 
			 elseif (txDict(iTx)%Tx_type=='CSEM') then
					do iDt=1,d_in%d(iTx)%nDt
					 do j=1,d_in%d(iTx)%data(iDt)%nSite
					   do k=1,d_in%d(iTx)%data(iDt)%nComp
							  d_in%d(iTx)%data(iDt)%error(k,j)= d_orgi%d(iTx)%data(iDt)%error(k,j)*term_CSEM
					  end do
					end do
				   end do 
		     end if		 	   
        end do 
		
end subroutine set_reset_data_error		
		
!**********************************************************************
subroutine split_dataVectorMTX(d_in)
   type(dataVectorMTX_t), intent(inout)           :: d_in
   
   !LOCAL
   integer                  :: iTx,Tx_MT_counter,Tx_DC_counter,Tx_CSEM_counter,counter_MT,counter_CSEM,counter_DC
   character(100)           ::  dataFile
	
	

   
 !Split d into two data: CSEM and MT 
   Tx_MT_counter=0
   Tx_CSEM_counter=0  
   Tx_DC_counter=0   
   
   do iTx=1,d_in%nTx
	 if (txDict(iTx)%Tx_type=='MT') then
	   Tx_MT_counter=Tx_MT_counter+1
	 elseif (txDict(iTx)%Tx_type=='CSEM') then
	   Tx_CSEM_counter=Tx_CSEM_counter+1  
	 elseif (txDict(iTx)%Tx_type=='DC') then
	   Tx_DC_counter=Tx_DC_counter+1
     end if	 
  enddo
if (Tx_MT_counter .ne. 0  ) then   

           if (MT_d%allocated) then
		     call deall_dataVectorMTX(MT_d)
		   end if
		   
		   call create_dataVectorMTX(Tx_MT_counter,MT_d)

		   
		   counter_MT=0
		   
		 do iTx=1,d_in%nTx
			if (txDict(iTx)%Tx_type=='MT') then
				 counter_MT=counter_MT+1
				  call create_dataVector(d_in%d(iTx)%ndt,MT_d%d(counter_MT))
				  call copy_dataVector(MT_d%d(counter_MT), d_in%d(iTx))		 
			 end if
		end do	
				  MT_d%allocated = .true.
				  
				 Ndata_MT = countData(MT_d)
end if
if (Tx_CSEM_counter .ne. 0 ) then   


           if (CSEM_d%allocated) then
		     call deall_dataVectorMTX(CSEM_d)
		   end if
		   
		   call create_dataVectorMTX(Tx_CSEM_counter,CSEM_d)
		   counter_CSEM=0

		   
		 do iTx=1,d_in%nTx
           if (txDict(iTx)%Tx_type=='CSEM') then
				  counter_CSEM=counter_CSEM+1
				  call create_dataVector(d_in%d(iTx)%ndt,CSEM_d%d(counter_CSEM))
				  call copy_dataVector(CSEM_d%d(counter_CSEM), d_in%d(iTx))
			 end if
		end do	
				  CSEM_d%allocated = .true.
				 Ndata_CSEM = countData(CSEM_d)
end if

if (Tx_DC_counter .ne. 0 ) then   

           if (DC_d%allocated) then
		     call deall_dataVectorMTX(DC_d)
		   end if
		   
		   call create_dataVectorMTX(Tx_DC_counter,DC_d)
		   

		   counter_DC=0
		   
		 do iTx=1,d_in%nTx
            if (txDict(iTx)%Tx_type=='DC') then
				  counter_DC=counter_DC+1
				  call create_dataVector(d_in%d(iTx)%ndt,DC_d%d(counter_DC))
				  call copy_dataVector(DC_d%d(counter_DC), d_in%d(iTx))
			 end if
		end do	

				  DC_d%allocated = .true.
				  
				 Ndata_DC = countData(DC_d)
end if
		 
         !write(40,*) 'RMS_MT= ',sqrt(dotProd(MT_d,MT_d)/Ndata_MT),'Ndata_CSEM= ',Ndata_CSEM
		 !Splitting the data was OK. But when writting the data into a file there was a problem in writting the second data type!!!
		 
		 
		 

end subroutine split_dataVectorMTX
!**********************************************************************
subroutine scale_dataVectorMTX(d_in,d_out,m0,eAll)

   type(dataVectorMTX_t), intent(inout)           :: d_in
   type(dataVectorMTX_t), optional, intent(out)   :: d_out
   type(modelParam_t),optional, intent(in)           :: m0
   type(solnVectorMTX_t), intent(inout)            :: eAll
   type(dataVectorMTX_t)                          :: d
   
   !  local variables
   real(kind=prec)       :: Ndata,Nmodel,mNorm,resNorm,CSEM_data_norm,MT_data_norm,scale_norm,scaling_factor,Ndata_CSEM,Ndata_MT,grad_CSEM,grad_MT
   type(dataVectorMTX_t)    :: res,Nres,MT_res,CSEM_res,Jm_MT,Jm_CSEM,Jm
   integer            :: iTx,Tx_MT_counter,Tx_CSEM_counter,counter_MT,counter_CSEM,iDt,k,j,methode
   character(100)       :: mFile,  dataFile
   type(dataVectorMTX_t)    :: MT_d,CSEM_d
  type(modelParam_t)        ::  mHat
   
   
    d = d_in



!################################################################## 

		 
		 
! Compute data gradient of MT and CSEM	 J m
        Jm =d	 
		mHat=m0
		!call zero(mHat)
		
#ifdef MPI
            !call zero_dataVectorMTX(Jm)
	        !call Master_job_Jmult(mHat,m0,Jm,eAll) 
#else
	        !call Jmult(mHat,m0,Jm,eAll)
#endif


!Split Jm into two data: CSEM and MT 
   Tx_MT_counter=0
   Tx_CSEM_counter=0   
   
   do iTx=1,d%nTx
	 if (txDict(iTx)%Tx_type=='MT') then
	   Tx_MT_counter=Tx_MT_counter+1
	 elseif (txDict(iTx)%Tx_type=='CSEM') then
	   Tx_CSEM_counter=Tx_CSEM_counter+1
     end if	   
   enddo
   
   call create_dataVectorMTX(Tx_MT_counter,MT_d)
   call create_dataVectorMTX(Tx_CSEM_counter,CSEM_d)
   
   counter_MT=0
   counter_CSEM=0
   
 do iTx=1,Jm%nTx
	  if (txDict(iTx)%Tx_type=='MT') then
	     counter_MT=counter_MT+1
	      call create_dataVector(Jm%d(iTx)%ndt,MT_d%d(counter_MT))
	 	 do idt=1,Jm%d(iTx)%ndt
		   call copy_dataBlock(MT_d%d(counter_MT)%data(idt), Jm%d(iTx)%data(idt))
		 end do
		 MT_d%d(counter_MT)%allocated = .true.
		 
	elseif (txDict(iTx)%Tx_type=='CSEM') then
	 	  counter_CSEM=counter_CSEM+1
	      call create_dataVector(Jm%d(iTx)%ndt,CSEM_d%d(counter_CSEM))
		  do idt=1,Jm%d(iTx)%ndt
			call copy_dataBlock(CSEM_d%d(counter_CSEM)%data(idt), Jm%d(iTx)%data(idt))
		  end do
		  CSEM_d%d(counter_CSEM)%allocated = .true.
	end if
end do	
          MT_d%allocated = .true.
		  CSEM_d%allocated = .true.
		  
         Ndata_MT = countData(MT_d)
		 Ndata_CSEM = countData(CSEM_d)
		 
         write(20,*) 'Ndata_MT= ',Ndata_MT,'Ndata_CSEM= ',Ndata_CSEM
		 
         !MT_data_norm = sqrt(dotProd(MT_d,MT_d))
         !CSEM_data_norm = sqrt(dotProd(CSEM_d,CSEM_d))


	grad_MT= MT_data_norm
    grad_CSEM= CSEM_data_norm
	         write(20,*) 'grad_MT= ',grad_MT,'grad_CSEM= ',grad_CSEM
			 
			 
methode=3	

	!Applying the weighting scheme number 1 used in Commer and Newmann: This method based on the number of data of each data type. The data set with small data point will be
    ! up-weighted by multiplying its error by the factor sqrt(N1/N2) where N1 > N2.
if (methode==1)then	
			if (Ndata_MT .gt. Ndata_CSEM) then
			  scaling_factor_CSEM=sqrt(Ndata_MT/Ndata_CSEM)
			  scaling_factor_MT=1.0
			else
			 scaling_factor_MT=sqrt(Ndata_CSEM/Ndata_MT)
			 scaling_factor_CSEM=1.0
			end if
elseif (methode==2) then			
			!Applying the weighting scheme number 2 used in Commer and Newmann: This method based on the data gradient of data of each data type. 
			
			if (Ndata_MT .gt. Ndata_CSEM) then
			  scaling_factor_CSEM=sqrt(grad_MT/grad_CSEM)
			  scaling_factor_MT=1.0
			else
			 scaling_factor_MT=sqrt(grad_CSEM/grad_MT)
			 scaling_factor_CSEM=1.0
			end if
elseif (methode==3) then			
			!Applying the weighting scheme number 2 used in Commer and Newmann: This method based on the data gradient of data of each data type. 
			
			if (Ndata_MT .gt. Ndata_CSEM) then
			  scaling_factor_CSEM=1.0
			  scaling_factor_MT=1.0
			else
			 scaling_factor_MT=1.0
			 scaling_factor_CSEM=1.0
			end if
end if
	
	
	write(20,*) 'scalling_MT=', scaling_factor_MT,'scalling_CSEM=', scaling_factor_CSEM
	
	       do iTx=1,d%nTx
			 if (txDict(iTx)%Tx_type=='MT') then
					do iDt=1,d%d(iTx)%nDt
					 do j=1,d%d(iTx)%data(iDt)%nSite
					   do k=1,d%d(iTx)%data(iDt)%nComp
							  d%d(iTx)%data(iDt)%error(k,j)= d%d(iTx)%data(iDt)%error(k,j)/scaling_factor_MT
					  end do
					end do
				   end do 
			 elseif (txDict(iTx)%Tx_type=='CSEM') then
					do iDt=1,d%d(iTx)%nDt
					 do j=1,d%d(iTx)%data(iDt)%nSite
					   do k=1,d%d(iTx)%data(iDt)%nComp
							  d%d(iTx)%data(iDt)%error(k,j)= d%d(iTx)%data(iDt)%error(k,j)/scaling_factor_CSEM
					  end do
					end do
				   end do 
		     end if		 	   
        end do 
		
		
     dataFile = 'weighted_error.dat'
     call write_dataVectorMTX(d,trim(dataFile))	
	
   	 
		 
		 

	 
scaling_factor_MT=scaling_factor
scaling_factor_CSEM=1.0-scaling_factor_MT
scaling_factor_MT=scaling_factor_MT*(CSEM_data_norm/MT_data_norm)

 do iTx=1,d%nTx
	  if (txDict(iTx)%Tx_type=='MT') then
		 do idt=1,d%d(iTx)%ndt
		  ! call scMult_dataBlock(scaling_factor_MT,d%d(iTx)%data(idt),d%d(iTx)%data(idt))
		 end do
	elseif (txDict(iTx)%Tx_type=='CSEM') then
		  do idt=1,d%d(iTx)%ndt
          ! call scMult_dataBlock(scaling_factor_CSEM,d%d(iTx)%data(idt),d%d(iTx)%data(idt))
		  end do
	end if
end do	

      !mNorm = sqrt(dotProd(mHat,mHat))
   !resNorm = sqrt(dotProd(res,Nres))
   !write(6,*) resNorm,sqrt(dotProd(res,Nres)/Ndata),mNorm,sqrt(dotProd(m,m))
   ! multiply by 2 (to be consistent with the formula)
   ! and add the gradient of the model norm
   !if (mNorm .eq. 0.0) then
   !  mNorm=Nmodel
   !end if
   !if (resNorm .eq. 0.0) then
   !  resNorm=1.0
   !end if  
    !scaling_factor=resNorm/mNorm 

	 
	 !##################################################################  
	 
	 
	    	
	if (present(d_out)) then
   		d_out = d
   	else
   	    d_in = d
   	end if

   	call deall(d)
   	call deall(MT_d)
   	call deall(CSEM_d)	
end subroutine scale_dataVectorMTX
   
!**********************************************************************
   subroutine update_damping_parameter(lambda,mHat,F,grad)

   real(kind=prec), intent(inout)  :: lambda
   type(modelParam_t), intent(in)              :: mHat
   real(kind=prec), intent(inout)  :: F
   type(modelParam_t), intent(inout)             :: grad

   real(kind=prec) :: SS, mNorm, Nmodel
	 type(modelParam_t)          :: dSS

   ! compute the model norm
   mNorm = dotProd(mHat,mHat)
   Nmodel = countModelParam(mHat)

   ! (scaled) sum of squares = penalty functional - scaled model norm
   SS = F - (lambda * mNorm/Nmodel)

   ! initialize
   dSS = mHat

   ! subtract the model norm derivative from the gradient of the penalty functional
   call linComb(ONE,grad,MinusTWO*lambda/Nmodel,mHat,dSS)

   ! update the damping parameter lambda
   lambda = lambda/iterControl%k

   ! penalty functional = (scaled) sum of squares + scaled model norm
   F = SS + (lambda * mNorm/Nmodel)

	 ! add the model norm derivative to the gradient of the penalty functional
   call linComb(ONE,dSS,TWO*lambda/Nmodel,mHat,grad)

   call deall_modelParam(dSS)

   end subroutine update_damping_parameter
!**********************************************************************
  subroutine un_normalize_with_dataVecMTX(d_in_out,d,N)
  
     type(dataVectorMTX_t), intent(in)              :: d
     type(dataVectorMTX_t), intent(inout)           :: d_in_out
     integer, optional, intent(in)               :: N
     integer                                     :: i,j,k,iDt
     

 	            do i=1,d%nTx
	             do iDt=1,d%d(i)%nDt
	              do j=1,d%d(i)%data(iDt)%nSite
	               do k=1,d%d(i)%data(iDt)%nComp	 
						   d_in_out%d(i)%data(iDt)%value(k,j)=  (d_in_out%d(i)%data(iDt)%value(k,j)*d%d(i)%data(iDt)%error(k,j)**N)
	                end do      
	              end do
				           d_in_out%d(i)%data(iDt)%errorBar=.true.
	            end do                                       
	        end do    
     

     
  
  end subroutine un_normalize_with_dataVecMTX
!**********************************************************************
!**********************************************************************
   subroutine CdInvMult(d_in,d_out)

   ! Divides by the data covariance C_d, which is a diagonal
   ! operator. Divides by the variances (squared error bars)
   ! and scales by the number of data (degrees of freedom).

   type(dataVectorMTX_t), intent(inout)           :: d_in
   type(dataVectorMTX_t), optional, intent(out)   :: d_out
   type(dataVectorMTX_t)                          :: d
   !integer                                :: Ndata

    d = d_in

    ! divide each data component by its variance
    call normalize_dataVectorMTX(d,2)

    ! divide by the number of data
    !Ndata = countData(d)
    !d = scMult(ONE/Ndata,d)

   	if (present(d_out)) then
   		d_out = d
   	else
   	    d_in = d
   	end if

   	call deall(d)

   end subroutine CdInvMult


!**********************************************************************
   subroutine CmSqrtMult(m_in,m_out)

   ! Multiplies by the square root of the model covariance,
   ! which is viewed as a smoothing operator. Intended
   ! to be used to compute m = C_m^{1/2} \tilde{m} + m_0.
   ! For efficiency, CmSqrt is a saved, private variable inside
   ! the modelParam module. Before this routine can be called,
   ! it has to be initialized by calling create_CmSqrt(m).
   ! Now that multBy_CmSqrt routine exists in modelParam,
   ! this routine is no longer needed. Leaving it here for now,
   ! to minimize changes.

   type(modelParam_t), intent(in)              :: m_in
   type(modelParam_t), intent(out)             :: m_out

	! apply the operator Cm^(1/2) here
	! m_out = m_in
	m_out = multBy_CmSqrt(m_in)

   end subroutine CmSqrtMult

!**********************************************************************
   subroutine NLCGsolver(d,lambda,m0,m,fname)

   ! computes inverse solution minimizing penalty functional
   !   for fixed value of regularization parameter, using
   !   a variant of non-linear conjugate gradient search.
   !   Various flavours of the algorithm and of the line search
   !   can be called from this routine
   !
   !  Note about the starting model:
   !  The starting model has to be in the smoothed model space,
   !  i.e. of the form m = C_m^{1/2} \tilde{m} + m_0.
   !  In order to compute \tilde{m} from the starting model,
   !  C_m^{-1/2} has to be implemented. To avoid this issue
   !  altogether, we are always starting from the prior,
   !  with \tilde{m} = 0. However, in general we could also
   !  start with the result of a previous search.

   !  d is data; on output it contains the responses for the inverse model
   type(dataVectorMTX_t), intent(inout)		   :: d
   !  lambda is regularization parameter
   real(kind=prec), intent(inout)  :: lambda
   !  m0 is prior model parameter
   type(modelParam_t), intent(in)		       :: m0
   !  m is solution parameter ... on input m contains starting guess
   type(modelParam_t), intent(inout)	       :: m
   !  flavor is a string that specifies the algorithm to use
   character(*), intent(in), optional        :: fname
   !  initial step size in the line search direction in model units
   real(kind=prec)                           :: startdm
   !  flavor is a string that specifies the algorithm to use
   character(80)                           :: flavor = 'Cubic'

   !  local variables
   type(dataVectorMTX_t)			:: dHat, res
   type(modelParam_t)			:: mHat, m_minus_m0, grad, g, h, gPrev
   !type(NLCGiterControl_t)			:: iterControl
   real(kind=prec)		:: value, valuePrev, rms, rmsPrev, alpha, beta, gnorm, mNorm, Nmodel
   real(kind=prec)      :: grad_dot_h, g_dot_g, g_dot_gPrev, gPrev_dot_gPrev, g_dot_h
   integer				:: iter, nCG, nLS, nfunc, ios
   logical              :: ok
   character(3)         :: iterChar
   character(100)       :: mFile, mHatFile, gradFile, dataFile, resFile, logFile
   type(solnVectorMTX_t)      :: eAll

   
   
         scaling_factor_CSEM = ONE 
         scaling_factor_MT = ONE
	  
	  
   if (present(fname)) then
      call read_NLCGiterControl(iterControl,fname,ok)
      if (ok) then
         lambda = iterControl%lambda
      end if
   else
      call set_NLCGiterControl(iterControl)
   end if
   !call scale_dataVectorMTX(d)
    !call scale_dataVectorMTX(d,d,m0,eAll)

   ! initialize the output to log file
   logFile = trim(iterControl%fname)//'_NLCG.log'
   open (unit=ioLog,file=logFile,status='unknown',position='append',iostat=ios)
   
   ! initialize the line search
   alpha = iterControl%alpha_1
   startdm = iterControl%startdm

   write(*,'(a41,es8.1)') 'The initial damping parameter lambda is ',lambda
   write(*,'(a55,f12.6)') 'The initial line search step size (in model units) is ',startdm

   write(ioLog,'(a41,es8.1)') 'The initial damping parameter lambda is ',lambda
   write(ioLog,'(a55,f12.6)') 'The initial line search step size (in model units) is ',startdm


   ! starting from the prior hardcoded by setting mHat = 0 and m = m0
   ! m = m0
   ! mHat = m0
   ! call zero(mHat)

   ! starting model contains the rough deviations from the prior
   mHat = m

   !  compute the penalty functional and predicted data
   call func(lambda,d,m0,mHat,value,mNorm,dHat,eAll,rms)
   call printf('START',lambda,alpha,value,mNorm,rms)
   call printf('START',lambda,alpha,value,mNorm,rms,logFile)
	 nfunc = 1
   write(iterChar,'(i3.3)') 0
   which_iter=0

   ! output initial model and responses for later reference
   if (output_level > 1) then
     mFile = trim(iterControl%fname)//'_NLCG_'//iterChar//'.rho'
     call write_modelParam(m,trim(mFile))
   end if
   if (output_level > 2) then
     dataFile = trim(iterControl%fname)//'_NLCG_'//iterChar//'.dat'
     call write_dataVectorMTX(dHat,trim(dataFile))
   end if

   ! compute gradient of the full penalty functional
   call gradient(lambda,d,m0,mHat,grad,dHat,eAll)
   if (output_level > 4) then
     gradFile = trim(iterControl%fname)//'_NLCG_'//iterChar//'.grt'
     call write_modelParam(grad,trim(gradFile))
   end if

   ! update the initial value of alpha if necessary
   gnorm = sqrt(dotProd(grad,grad))
   write(*,'(a37,es12.6)') 'The initial norm of the gradient is ',gnorm
   write(ioLog,'(a37,es12.6)') 'The initial norm of the gradient is ',gnorm
   if (gnorm < TOL6) then
      call errStop('Problem with your gradient computations: first gradient is zero')
   else !if (alpha * gnorm > startdm) then
      alpha = startdm / gnorm
      write(*,'(a39,es12.6)') 'The initial value of alpha updated to ',alpha
      write(ioLog,'(a39,es12.6)') 'The initial value of alpha updated to ',alpha
   end if

   ! initialize CG: g = - grad; h = g
   nCG = 0
   iter = 0
   g = grad
   call linComb(MinusONE,grad,R_ZERO,grad,g)
   h = g

   do
      !  test for convergence ...
      if((rms.lt.iterControl%rmsTol).or.(iter.ge.iterControl%maxIter)) then
         exit
      end if

	  iter = iter + 1
      which_iter=iter
	  ! save the values of the functional and the directional derivative
		rmsPrev = rms
	  valuePrev = value
	  grad_dot_h = dotProd(grad,h)

	  ! at the end of line search, set mHat to the new value
	  ! mHat = mHat + alpha*h  and evaluate gradient at new mHat
	  ! data and solnVector only needed for output
      write(*,'(a23)') 'Starting line search...'
      write(ioLog,'(a23)') 'Starting line search...'
	  select case (flavor)
	  case ('Cubic')
	  	call lineSearchCubic(lambda,d,m0,h,alpha,mHat,value,grad,rms,nLS,dHat,eAll)
	  	!call deall(eAll)
	  case ('Quadratic')
	  	call lineSearchQuadratic(lambda,d,m0,h,alpha,mHat,value,grad,rms,nLS,dHat,eAll)
	  	!call deall(eAll)
	  case default
        call errStop('Unknown line search requested in NLCG')
	  end select
		nfunc = nfunc + nLS
	  gPrev = g
	  call linComb(MinusONE,grad,R_ZERO,grad,g)

	  ! compute the starting step for the next line search
	  alpha = 2*(value - valuePrev)/grad_dot_h

	  ! adjust the starting step to ensure superlinear convergence properties
	  alpha = (ONE+0.01)*alpha
	  write(*,'(a25,i5)') 'Completed NLCG iteration ',iter
	  write(ioLog,'(a25,i5)') 'Completed NLCG iteration ',iter
	  Nmodel = countModelParam(mHat)
	  mNorm = dotProd(mHat,mHat)/Nmodel
      call printf('with',lambda,alpha,value,mNorm,rms)
      call printf('with',lambda,alpha,value,mNorm,rms,logFile)

      ! write out the intermediate model solution and responses
      call CmSqrtMult(mHat,m_minus_m0)
   	  call linComb(ONE,m_minus_m0,ONE,m0,m)
   	  write(iterChar,'(i3.3)') iter
   	  if (output_level > 1) then
   	    mFile = trim(iterControl%fname)//'_NLCG_'//iterChar//'.rho'
        call write_modelParam(m,trim(mFile))
      end if
   	  if (output_level > 2) then
   	    mHatFile = trim(iterControl%fname)//'_NLCG_'//iterChar//'.prm'
        call write_modelParam(mHat,trim(mHatFile))
      end if
   	  if (output_level > 2) then
   	    dataFile = trim(iterControl%fname)//'_NLCG_'//iterChar//'.dat'
        call write_dataVectorMTX(dHat,trim(dataFile))
      end if
      ! compute residual for output: res = d-dHat; do not normalize by errors
   	  if (output_level > 2) then
        res = d
        call linComb(ONE,d,MinusONE,dHat,res)
   	    resFile = trim(iterControl%fname)//'_NLCG_'//iterChar//'.res'
        call write_dataVectorMTX(res,trim(resFile))
      end if

	  ! if alpha is too small, we are not making progress: update lambda
      if (abs(rmsPrev - rms) < iterControl%fdiffTol) then
      		! update lambda, penalty functional and gradient
      		call update_damping_parameter(lambda,mHat,value,grad)
      		! update alpha
      		gnorm = sqrt(dotProd(grad,grad))
            write(*,'(a34,es12.6)') 'The norm of the last gradient is ',gnorm
            write(ioLog,'(a34,es12.6)') 'The norm of the last gradient is ',gnorm
            !alpha = min(iterControl%alpha_1,startdm/gnorm)
      		alpha = min(ONE,startdm)/gnorm
      		write(*,'(a48,es12.6)') 'The value of line search step alpha updated to ',alpha
            write(ioLog,'(a48,es12.6)') 'The value of line search step alpha updated to ',alpha
      		! g = - grad
			call linComb(MinusONE,grad,R_ZERO,grad,g)
			! check that lambda is still at a reasonable value
			if (lambda < iterControl%lambdaTol) then
				write(*,'(a55)') 'Unable to get out of a local minimum. Exiting...'
                write(ioLog,'(a55)') 'Unable to get out of a local minimum. Exiting...'
				! multiply by C^{1/2} and add m_0
                call CmSqrtMult(mHat,m_minus_m0)
                call linComb(ONE,m_minus_m0,ONE,m0,m)
                d = dHat
				return
			end if
	  	! restart
			write(*,'(a55)') 'Restarting NLCG with the damping parameter updated'
			call printf('to',lambda,alpha,value,mNorm,rms)
			write(ioLog,'(a55)') 'Restarting NLCG with the damping parameter updated'
			call printf('to',lambda,alpha,value,mNorm,rms,logFile)
	  	h = g
	  	nCG = 0
	  	cycle
	  end if

	  g_dot_g = dotProd(g,g)
	  g_dot_gPrev = dotProd(g,gPrev)
	  gPrev_dot_gPrev = dotProd(gPrev,gPrev)
	  g_dot_h = dotProd(g,h)

	  ! Polak-Ribiere variant
	  beta = ( g_dot_g - g_dot_gPrev )/gPrev_dot_gPrev

	  ! restart CG if the orthogonality conditions fail. Using the fact that
		! h_{i+1} = g_{i+1} + beta * h_i. In order for the next directional
		! derivative = -g_{i+1}.dot.h_{i+1} to be negative, the condition
		! g_{i+1}.dot.(g_{i+1}+beta*h_i) > 0 must hold. Alternatively, books
		! say we can take beta > 0 (didn't work as well)
	  if ((g_dot_g + beta*g_dot_h > 0).and.(nCG < iterControl%nCGmax)) then
      	call linComb(ONE,g,beta,h,h)
      	nCG = nCG + 1
	  else
   	    ! restart
		write(*,'(a45)') 'Restarting NLCG to restore orthogonality'
		write(ioLog,'(a45)') 'Restarting NLCG to restore orthogonality'
        h = g
        nCG = 0
   	  end if

   end do

   ! multiply by C^{1/2} and add m_0
   call CmSqrtMult(mHat,m_minus_m0)
   call linComb(ONE,m_minus_m0,ONE,m0,m)
   d = dHat
   write(*,'(a25,i5,a25,i5)') 'NLCG iterations:',iter,' function evaluations:',nfunc
   write(ioLog,'(a25,i5,a25,i5)') 'NLCG iterations:',iter,' function evaluations:',nfunc
   close(ioLog,iostat=ios)
   close(ioWeighting,iostat=ios)

   ! cleaning up
   call deall_dataVectorMTX(dHat)
   call deall_dataVectorMTX(res)
   call deall_modelParam(mHat)
   call deall_modelParam(m_minus_m0)
   call deall_modelParam(grad)
   call deall_modelParam(g)
   call deall_modelParam(h)
   call deall_modelParam(gPrev)
   call deall_solnVectorMTX(eAll)

   end subroutine NLCGsolver

!**********************************************************************
  subroutine lineSearchQuadratic(lambda,d,m0,h,alpha,mHat,f,grad, &
  								rms,niter,dHat,eAll,gamma)

   ! Line search that imitates the strategy of Newman & Alumbaugh (2000),
   ! except without the errors. In particular, we only test the sufficient
   ! decrease (Armijo) condition (ignoring the curvature condition) and
   ! we use quadratic (not cubic) interpolation for backtracking.
   ! This strategy only requires one gradient evaluation (but so does
   ! the cubic interpolation method described in the Numerical Recipes).
   ! This is likely to be less efficient than the cubic interpolation,
   ! but it is simple to implement, and in some cases will work just as
   ! well (assuming an adequate initial step size has been chosen).
   !
   ! The initial step size is set outside of this routine (in the NLCG)
   ! but these are the good choices (ref. Michael Ferris, Chapter 3, p 59):
   ! alpha_1 = alpha_{k-1} dotProd(grad_{k-1},h_{k-1})/dotProd(grad_k,h_k})
   ! or interpolate the quadratic to f(m_{k-1}), f(m_k) and
   ! dotProd(grad_{k-1},h_{k-1}) and find the minimizer
   ! alpha_1 = 2(f_k - f_{k-1})/dotProd(grad_{k-1},h_{k-1}),
   ! the update alpha_1 <- min(1.00,1.01 * alpha_1).
   !
   ! Set f(alpha) = func(mHat + alpha*h). Fit the quadratic
   !     f_q(alpha) = a alpha^2 + b alpha + f(0)
   ! using the information f(0), f'(0) and f(alpha_1) to obtain
   ! a = (f(alpha_1) - f(0) - f'(0) alpha_1)/(alpha_1 * alpha_1),
   ! b = f'(0).
   ! Then, the minimum point of the quadratic is alpha_q = -b/(2a),
   ! assuming that a > 0. If this try is not successful, fit another
   ! quadratic using f(0), f'(0) and f(alpha_q). The new quadratic
   ! is not identical to the previous quadratic since f_q is only
   ! an approximation to f: in general, f(alpha_q) /= f_q(alpha_q),
   ! hence the new point does not lie on the same quadratic curve.
   !
   ! Our solution has to satisfy the sufficient decrease condition
   !     f(alpha) < f(0) + c alpha f'(0).
   !
   ! The optional relaxation parameter gamma is needed for algorithms
   ! like the Renormalised Steepest Descent (RSD). See the dynamical
   ! systems in optimisation research (Pronzato et al [2000, 2001]).
   ! To the best of my knowledge, it is not useful for NLCG.

   real(kind=prec), intent(in)     :: lambda
   type(dataVectorMTX_t), intent(inout)		       :: d
   type(modelParam_t), intent(in)		       :: m0
   type(modelParam_t), intent(in)            :: h  ! search direction
   real(kind=prec), intent(inout)  :: alpha ! step size
   type(modelParam_t), intent(inout)         :: mHat
   real(kind=prec), intent(inout)  :: f
   type(modelParam_t), intent(inout)         :: grad
   real(kind=prec), intent(out)    :: rms
   integer,intent(out)                     :: niter
   type(dataVectorMTX_t), intent(out)         :: dHat
   type(solnVectorMTX_t), intent(inout)          :: eAll

   ! optionally add relaxation (e.g. for Renormalised Steepest Descent)
   real(kind=prec), intent(in), optional :: gamma

   ! local variables
   real(kind=prec)                 :: alpha_1,alpha_i,mNorm
   logical                                 :: starting_guess
   logical                                 :: relaxation
   real(kind=prec)                 :: eps,k,c,a,b
   real(kind=prec)                 :: g_0,f_0,f_1,f_i,rms_1,mNorm_1
   type(modelParam_t)                        :: mHat_0,mHat_1
   type(dataVectorMTX_t)                           :: dHat_1
   type(solnVectorMTX_t)                         :: eAll_1
   character(100)							:: logFile

   ! parameters
   c = iterControl%c
   !k = iterControl%alpha_k
   !eps = iterControl%alpha_tol
   logFile = trim(iterControl%fname)//'_NLCG.log'

   ! initialize the line search
   niter = 0
   mHat_0 = mHat
   f_0 = f
   starting_guess = .false.

   ! rescale the search direction
   !h_dot_h = dotProd(h,h)
   !h = scMult_modelParam(ONE/sqrt(h_dot_h),h)

   ! g_0 is the directional derivative f'(0) = (df/dm).dot.h
   g_0 = dotProd(grad,h)

   ! alpha_1 is the initial step size, which is set in NLCG
   !alpha_1 = ONE/maxNorm_modelParam(h)
   alpha_1 = alpha

   ! with relaxation, we specify gamma = 1 - eps, eps > 0 small; then the final
   ! solution is f(gamma*alpha) = func(mHat + gamma*alpha*h)
   if (present(gamma)) then
      relaxation = .true.
   else
      relaxation = .false.
   end if

   ! initialize
   mHat_1 = mHat_0
   !  compute the trial parameter mHat_1
   call linComb(ONE,mHat_0,alpha_1,h,mHat_1)

   !  compute the penalty functional and predicted data at mHat_1
   call func(lambda,d,m0,mHat_1,f_1,mNorm_1,dHat_1,eAll_1,rms_1)
   call printf('STARTLS',lambda,alpha,f_1,mNorm_1,rms_1)
   call printf('STARTLS',lambda,alpha,f_1,mNorm_1,rms_1,logFile)
   niter = niter + 1

   if (f_1 - f_0 >= LARGE_REAL) then
      print *, 'Try a smaller starting value of alpha.'
      print *, 'Exiting...'
      stop
   end if

   f_i = f_1
   alpha_i = alpha_1

   fit_quadratic: do
    a = (f_i - f_0 - g_0*alpha_i)/(alpha_i**2)
    b = g_0
    alpha = - b/(TWO*a) ! the minimizer of the quadratic
    ! if the quadratic has negative curvature & no minimum, exit
    if (a < 0) then
    	starting_guess = .true.
    	exit
    end if
	!	The step size alpha should not be adjusted manually at all!
	! Even when it is too small or too close to the previous try,
	! adjusting it won't result in an improvement (it's better to exit
	! the line search in that case, if anything)...
  !  if ((alpha_i - alpha < eps).or.(alpha < k*alpha_i)) then
  !  	alpha = alpha_i/TWO ! reset alpha to ensure progress
  !  end if
    call linComb(ONE,mHat_0,alpha,h,mHat)
    call func(lambda,d,m0,mHat,f,mNorm,dHat,eAll,rms)
    call printf('QUADLS',lambda,alpha,f,mNorm,rms)
    call printf('QUADLS',lambda,alpha,f,mNorm,rms,logFile)
    niter = niter + 1
    ! check whether the solution satisfies the sufficient decrease condition
    if (f < f_0 + c * alpha * g_0) then
        write(*,'(a60)') 'Good enough value found, exiting line search'
        write(ioLog,'(a60)') 'Good enough value found, exiting line search'
    	exit
    end if
    ! this should not happen, but in practice it is possible to end up with
    ! a function increase at this point (e.g. in the current global code).
    ! Most likely, this is due to an inaccuracy in the gradient computations.
    ! In this case, we avoid an infinite loop by exiting the line search.
    if (f > f_0) then
        write(*,'(a75)') 'Unable to fit a quadratic due to bad gradient estimate, exiting line search'
        write(ioLog,'(a75)') 'Unable to fit a quadratic due to bad gradient estimate, exiting line search'
   		exit
    end if
    ! otherwise, iterate, using the most recent value of f & alpha
    alpha_i = alpha
    f_i = f
   end do fit_quadratic

   ! if the initial guess was better than what we found, take it
   if (f_1 < f) then
   	starting_guess = .true.
   end if

   if (starting_guess) then
   	alpha = alpha_1
   	dHat = dHat_1
   	eAll = eAll_1
   	mHat = mHat_1
   	rms = rms_1
   	f = f_1
   end if

   ! compute gradient of the full penalty functional and exit
    if (relaxation) then
   		call linComb(ONE,mHat_0,gamma*alpha,h,mHat)
    	call func(lambda,d,m0,mHat,f,mNorm,dHat,eAll,rms)
   		call printf('RELAX',lambda,gamma*alpha,f,mNorm,rms)
   		call printf('RELAX',lambda,gamma*alpha,f,mNorm,rms,logFile)
   	end if
    call gradient(lambda,d,m0,mHat,grad,dHat,eAll)
    write(*,'(a39)') 'Gradient computed, line search finished'
    write(ioLog,'(a39)') 'Gradient computed, line search finished'

   call deall_dataVectorMTX(dHat_1)
   call deall_modelParam(mHat_0)
   call deall_modelParam(mHat_1)
   call deall_solnVectorMTX(eAll_1)

  end subroutine lineSearchQuadratic


  !**********************************************************************
  subroutine lineSearchCubic(lambda,d,m0,h,alpha,mHat,f,grad, &
  							rms,niter,dHat,eAll,gamma)

   ! Line search that is based on the Numerical Recipes and on the
   ! text by Michael Ferris, Chapter 3, p 59. We only test the sufficient
   ! decrease (Armijo) condition (ignoring the curvature condition).
   ! We first interpolate using a quadratic approximation; if the
   ! solution does not satisfy the condition, we backtrack using
   ! cubic interpolation. This strategy only requires one gradient
   ! evaluation and is very efficient when computing gradients is
   ! expensive.
   !
   ! The initial step size is set outside of this routine (in the NLCG)
   ! but these are the good choices (ref. Michael Ferris, Chapter 3, p 59):
   ! alpha_1 = alpha_{k-1} dotProd(grad_{k-1},h_{k-1})/dotProd(grad_k,h_k})
   ! or interpolate the quadratic to f(m_{k-1}), f(m_k) and
   ! dotProd(grad_{k-1},h_{k-1}) and find the minimizer
   ! alpha_1 = 2(f_k - f_{k-1})/dotProd(grad_{k-1},h_{k-1}),
   ! the update alpha_1 <- min(1.00,1.01 * alpha_1).
   !
   ! Set f(alpha) = func(mHat + alpha*h). Fit the quadratic
   !     f_q(alpha) = a alpha^2 + b alpha + f(0)
   ! using the information f(0), f'(0) and f(alpha_1) to obtain
   ! a = (f(alpha_1) - f(0) - f'(0) alpha_1)/(alpha_1 * alpha_1),
   ! b = f'(0).
   ! Then, the minimum point of the quadratic is alpha_q = -b/(2a),
   ! assuming that a > 0. If this try is not successful, fit a cubic
   !     f_c(alpha) = a alpha^3 + b alpha^2 + f'(0) alpha + f(0)
   ! using f(0), f'(0), f(alpha_1) and f(alpha_q). Repeat as necessary.
   ! Here, a and b are as described in the code.
   ! A new cubic is not identical to a previous curve since f_c is only
   ! an approximation to f: in general, f(alpha_c) /= f_c(alpha_c),
   ! hence the new point does not lie on the approximating curve.
   !
   ! Our solution has to satisfy the sufficient decrease condition
   !     f(alpha) < f(0) + c alpha f'(0).
   !
   ! The optional relaxation parameter gamma is needed for algorithms
   ! like the Renormalised Steepest Descent (RSD). See the dynamical
   ! systems in optimisation research (Pronzato et al [2000, 2001]).
   ! To the best of my knowledge, it is not useful for NLCG.

   real(kind=prec), intent(in)     :: lambda
   type(dataVectorMTX_t), intent(inout)		       :: d
   type(modelParam_t), intent(in)		       :: m0
   type(modelParam_t), intent(in)            :: h  ! search direction
   real(kind=prec), intent(inout)  :: alpha ! step size
   type(modelParam_t), intent(inout)         :: mHat
   real(kind=prec), intent(inout)  :: f
   type(modelParam_t), intent(inout)         :: grad
   real(kind=prec), intent(out)    :: rms
   integer, intent(out)                    :: niter
   type(dataVectorMTX_t), intent(out)         :: dHat
   type(solnVectorMTX_t), intent(inout)          :: eAll

   ! optionally add relaxation (e.g. for Renormalised Steepest Descent)
   real(kind=prec), intent(in), optional :: gamma

    ! local variables
   real(kind=prec)                 :: alpha_1,alpha_i,alpha_j,mNorm
   logical                                 :: starting_guess
   logical                                 :: relaxation
   real(kind=prec)                 :: eps,k,c,a,b,q1,q2,q3
   real(kind=prec)                 :: g_0,f_0,f_1,f_i,f_j,rms_1,mNorm_1
   type(modelParam_t)                        :: mHat_0,mHat_1
   type(dataVectorMTX_t)                           :: dHat_1
   type(solnVectorMTX_t)                         :: eAll_1
   character(100)							:: logFile

   ! parameters
   c = iterControl%c
   !k = iterControl%alpha_k
   !eps = iterControl%alpha_tol
   logFile = trim(iterControl%fname)//'_NLCG.log'

   ! initialize the line search
   niter = 0
   mHat_0 = mHat
   f_0 = f
   starting_guess = .false.

   ! rescale the search direction
   !h_dot_h = dotProd(h,h)
   !h = scMult_modelParam(ONE/sqrt(h_dot_h),h)

   ! g_0 is the directional derivative f'(0) = (df/dm).dot.h
   g_0 = dotProd(grad,h)

   ! alpha_1 is the initial step size, which is set in NLCG
   alpha_1 = alpha

   ! with relaxation, we specify gamma = 1 - eps, eps > 0 small; then the final
   ! solution is f(gamma*alpha) = func(mHat + gamma*alpha*h)
   if (present(gamma)) then
      relaxation = .true.
   else
      relaxation = .false.
   end if

   ! compute the trial mHat, f, dHat, eAll, rms
   mHat_1 = mHat_0
   call linComb(ONE,mHat_0,alpha_1,h,mHat_1)
   call func(lambda,d,m0,mHat_1,f_1,mNorm_1,dHat_1,eAll_1,rms_1)
   call printf('STARTLS',lambda,alpha,f_1,mNorm_1,rms_1)
   call printf('STARTLS',lambda,alpha,f_1,mNorm_1,rms_1,logFile)
   niter = niter + 1

	 if (f_1 - f_0 >= LARGE_REAL) then
		print *, 'Try a smaller starting value of alpha.'
		print *, 'Exiting...'
		stop
	 end if

   ! try fitting a quadratic
   a = (f_1 - f_0 - g_0*alpha_1)/(alpha_1**2)
   b = g_0
   ! if the curvature is -ve, there is no minimum; take the initial guess
   if (a < 0) then
	starting_guess = .true.
  	alpha = alpha_1
   	dHat = dHat_1
    eAll = eAll_1
   	mHat = mHat_1
   	rms = rms_1
   	f = f_1
    ! compute the gradient and exit
    if (relaxation) then
   		call linComb(ONE,mHat_0,gamma*alpha,h,mHat)
    	call func(lambda,d,m0,mHat,f,mNorm,dHat,eAll,rms)
   		call printf('RELAX',lambda,gamma*alpha,f,mNorm,rms)
   		call printf('RELAX',lambda,gamma*alpha,f,mNorm,rms,logFile)
   	end if
    call gradient(lambda,d,m0,mHat,grad,dHat,eAll)
    write(*,'(a45)') 'Quadratic has no minimum, exiting line search'
    write(ioLog,'(a45)') 'Quadratic has no minimum, exiting line search'
	call deall_dataVectorMTX(dHat_1)
	call deall_modelParam(mHat_0)
	call deall_modelParam(mHat_1)
	call deall_solnVectorMTX(eAll_1)
   	return
   end if

   ! otherwise compute the functional at the minimizer of the quadratic
   alpha = - b/(TWO*a)
   call linComb(ONE,mHat_0,alpha,h,mHat)
   call func(lambda,d,m0,mHat,f,mNorm,dHat,eAll,rms)
   call printf('QUADLS',lambda,alpha,f,mNorm,rms)
   call printf('QUADLS',lambda,alpha,f,mNorm,rms,logFile)
   niter = niter + 1
   ! check whether the solution satisfies the sufficient decrease condition
   if (f < f_0 + c * alpha * g_0) then
    ! if the initial guess was better than what we found, take it
   	if (f_1 < f) then
   		starting_guess = .true.
   		alpha = alpha_1
   		dHat = dHat_1
     	eAll = eAll_1
   		mHat = mHat_1
   		rms = rms_1
   		f = f_1
    end if
    ! compute the gradient and exit
    if (relaxation) then
   		call linComb(ONE,mHat_0,gamma*alpha,h,mHat)
    	call func(lambda,d,m0,mHat,f,mNorm,dHat,eAll,rms)
   		call printf('RELAX',lambda,gamma*alpha,f,mNorm,rms)
   		call printf('RELAX',lambda,gamma*alpha,f,mNorm,rms,logFile)
   	end if

    call gradient(lambda,d,m0,mHat,grad,dHat,eAll)
    write(*,'(a60)') 'Sufficient decrease condition satisfied, exiting line search'
    write(ioLog,'(a60)') 'Sufficient decrease condition satisfied, exiting line search'
	call deall_dataVectorMTX(dHat_1)
	call deall_modelParam(mHat_0)
	call deall_modelParam(mHat_1)
	call deall_solnVectorMTX(eAll_1)
   	return
   end if

   ! this should not happen, but in practice it is possible to end up with
   ! a function increase at this point (e.g. in the current global code).
   ! Most likely, this is due to an inaccuracy in the gradient computations.
   ! In this case, we avoid an infinite loop by exiting line search.
   ! It is also possible that both f_1 and f are worse than the starting value!
   ! Then, take whichever is smaller. Ideally, want to decrease the tolerance
   ! for gradient computations if this happens.
   if (f > f_0) then

    write(*,'(a75)') 'Unable to fit a quadratic due to bad gradient estimate, exiting line search'
    write(ioLog,'(a75)') 'Unable to fit a quadratic due to bad gradient estimate, exiting line search'

   else
    ! fit a cubic and backtrack (initialize)
    alpha_i = alpha_1
    f_i = f_1
    alpha_j = alpha
    f_j = f
    fit_cubic: do
        ! compute the minimizer
   	    q1 = f_i - f_0 - g_0 * alpha_i
   	    q2 = f_j - f_0 - g_0 * alpha_j
   	    q3 = alpha_i**2 * alpha_j**2 * (alpha_j - alpha_i)
   	    a = (alpha_i**2 * q2 - alpha_j**2 * q1)/q3
   	    b = (alpha_j**3 * q1 - alpha_i**3 * q2)/q3
   	    alpha = (- b + sqrt(b*b - 3*a*g_0))/(3*a)
        ! if alpha is too close or too much smaller than its predecessor
        !  if ((alpha_j - alpha < eps).or.(alpha < k*alpha_j)) then
        !  	alpha = alpha_j/TWO ! reset alpha to ensure progress
        !  end if
        ! compute the penalty functional
        call linComb(ONE,mHat_0,alpha,h,mHat)
        call func(lambda,d,m0,mHat,f,mNorm,dHat,eAll,rms)
        call printf('CUBICLS',lambda,alpha,f,mNorm,rms)
        call printf('CUBICLS',lambda,alpha,f,mNorm,rms,logFile)
        niter = niter + 1
        ! check whether the solution satisfies the sufficient decrease condition
        if (f < f_0 + c * alpha * g_0) then
    	   exit
        end if
        ! if not, iterate, using the two most recent values of f & alpha
        alpha_i = alpha_j
        f_i = f_j
        alpha_j = alpha
        f_j = f
        ! check that the function still decreases to avoid infinite loops in case of a bug
        if (abs(f_j - f_i) < TOL8) then
           write(*,'(a69)') 'Warning: exiting cubic search since the function no longer decreases!'
           write(ioLog,'(a69)') 'Warning: exiting cubic search since the function no longer decreases!'
    	   exit
        end if
    end do fit_cubic
   end if

   if (f_1 < f) then
   	starting_guess = .true.
   end if

   ! if the initial guess was better than what we found, take it
   if (starting_guess) then
   	alpha = alpha_1
   	dHat = dHat_1
   	eAll = eAll_1
   	mHat = mHat_1
   	rms = rms_1
   	f = f_1
   end if

   ! compute gradient of the full penalty functional and exit
    if (relaxation) then
   		call linComb(ONE,mHat_0,gamma*alpha,h,mHat)
    	call func(lambda,d,m0,mHat,f,mNorm,dHat,eAll,rms)
   		call printf('RELAX',lambda,gamma*alpha,f,mNorm,rms)
   		call printf('RELAX',lambda,gamma*alpha,f,mNorm,rms,logFile)
   	end if
    call gradient(lambda,d,m0,mHat,grad,dHat,eAll)
	write(*,'(a39)') 'Gradient computed, line search finished'
    write(ioLog,'(a39)') 'Gradient computed, line search finished'

   call deall_dataVectorMTX(dHat_1)
   call deall_modelParam(mHat_0)
   call deall_modelParam(mHat_1)
   call deall_solnVectorMTX(eAll_1)

  end subroutine lineSearchCubic

end module NLCG
