module DCG

	use math_constants
 	use utilities
    use senscomp
    use main
#ifdef MPI
	Use Main_MPI
	use Sub_MPI
#endif
implicit none
  type  :: DCGiterControl_t
     !NOTE: For the standard DCG algorethem only the first 4 attributes of the DCGiterControl_t are used. 
     ! maximum number of iterations in one call to iterative solver
     integer            :: maxIter
     ! convergence criteria: return from solver if rms < rmsTol
     real (kind=prec)	:: rmsTol
     ! initial value of lambda (will not override the NLCG input argument)
     real (kind=prec)   :: lambda
     ! model and data output file name
     character(80)      :: fname

     
     real (kind=prec)   :: fdiffTol    
     real (kind=prec)   :: lambdaTol
     real (kind=prec)   :: k
     real (kind=prec)   :: c
     integer            :: nCGmax
     real (kind=prec)   :: alpha_1
     real (kind=prec)   :: startdm
     real (kind=prec)   :: gamma
  end type DCGiterControl_t
    type(DCGiterControl_t), private, save :: DCGiterControl
    
! iteration control for CG solver
  type  :: iterControl_t
     ! maximum number of iterations in one call to iterative solver
     integer					:: maxIt
     ! convergence criteria: return from solver if relative error < tol
      real(kind=prec) 			:: tol
     ! actual number of iterations before return
     integer					:: niter
     ! relative error for each iteration
      real(kind=prec) , pointer, dimension(:)	:: rerr
     ! logical variable indicating if algorithm "failed"
     logical					:: failed = .false.
  end type iterControl_t
  


public  :: DCGsolver

logical                                  :: keep_solution

Contains
!**********************************************************************
subroutine set_DCGiterControl(DCGiterControl)
   !  at present this is a simple routine to setup for
   !  iteration control; allowing parameters to be read from
   !   a file would be desireable.

   type(DCGiterControl_t), intent(inout)	:: DCGiterControl

     ! maximum number of iterations in one call to iterative solver
     DCGiterControl%maxIter = 3
     ! convergence criteria: return from solver if rms < rmsTol
     DCGiterControl%rmsTol  = 1.05
     ! initial value of lambda 
     DCGiterControl%lambda = 50.     
     ! model and data output file name
     DCGiterControl%fname = 'Modular'     
     
     
     ! inversion stalls when abs(rms - rmsPrev) < fdiffTol (2e-3 works well)
     DCGiterControl%fdiffTol = 2.0e-3
     ! exit if lambda < lambdaTol approx. 1e-4
     DCGiterControl%lambdaTol = 1.0e-8
     ! set lambda_i = lambda_{i-1}/k when the inversion stalls
     DCGiterControl%k = 10.
     ! the factor that ensures sufficient decrease in the line search >=1e-4
     DCGiterControl%c = 1.0e-4
     ! restart CG every nCGmax iterations to ensure conjugacy
     DCGiterControl%nCGmax = 8
     ! the starting step for the line search
     DCGiterControl%alpha_1 = 20.
     ! maximum initial delta mHat (overrides alpha_1)
     DCGiterControl%startdm = 20.
     ! optional relaxation parameter (Renormalised Steepest Descent algorithm)
     DCGiterControl%gamma = 0.99


end subroutine set_DCGiterControl
!**********************************************************************
subroutine setIterControl(CGiter)
   !  at present this is a simple routine to setup for
   !  iteration control; allowing parameters to be read from
   !   a file would be desireable.

   type(iterControl_t), intent(inout)	:: CGiter
   
   CGiter%maxit = 20
   CGiter%tol = 10E-4
   CGiter%niter = 0
   allocate(CGiter%rerr(0:CGiter%maxit))

   end subroutine setIterControl
!**********************************************************************
subroutine read_DCGiterControl(DCGiterControl,rFile,fileExists)

	type(DCGiterControl_t), intent(inout)	:: DCGiterControl
    character(*), intent(in)		        :: rFile
	logical, intent(out), optional          :: fileExists
    integer									:: ios
	logical                             	:: exists
	character(80)							:: string

    ! Initialize inverse solver configuration

    call set_DCGiterControl(DCGiterControl)

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

    read (ioInvCtrl,'(a36,a80)') string,DCGiterControl%fname
    if (output_level > 2) then
       write (*,*)
       write (*,'(a36,a80)') string,DCGiterControl%fname
    end if
    DCGiterControl%fname = adjustl(DCGiterControl%fname)
    read (ioInvCtrl,'(a36,g15.7)') string,DCGiterControl%lambda
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,DCGiterControl%lambda
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,DCGiterControl%k
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,DCGiterControl%k
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,DCGiterControl%startdm
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,DCGiterControl%startdm
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,DCGiterControl%fdiffTol
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,DCGiterControl%fdiffTol
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,DCGiterControl%rmsTol
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,DCGiterControl%rmsTol
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,DCGiterControl%lambdaTol
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,DCGiterControl%lambdaTol
    end if
    read (ioInvCtrl,'(a36,i4)') string,DCGiterControl%maxIter
    if (output_level > 2) then
       write (*,'(a36,i4)') string,DCGiterControl%maxIter
       write (*,*)
    end if
    close(ioInvCtrl)

end subroutine read_DCGiterControl
!**********************************************************************
subroutine DCGsolver(d,m0,m,lambda,fname)
  
  ! Subroutine to solve the inverse problem in data space using conjugate gradients (CG)  
   
   type(dataVectorMTX_t), intent(inout)        ::d
   !  m0 is prior model parameter
   type(modelParam_t), intent(in)		       ::m0
   !  m is solution parameter ... on input m contains starting guess
   type(modelParam_t), intent(inout)	       ::m
   !  lambda is regularization parameter
   real(kind=prec) , intent(inout)             ::lambda
   ! Inversion control file name
   character(*), intent(in), optional        :: fname
   
!  local variables
   type(dataVectorMTX_t)        :: dHat, b,dx,d_Pred,res,Nres,JmHat,Jm0,d_Pred_m0
   type(modelParam_t)			:: mHat,CmJTd,Cm_mHat
   real(kind=prec)		  		:: value,rms_old,F,mNorm,rms,rmsPrev
   integer						:: iter, ndata,DS_iter,DCG_iter,nfunc
   character(100)       		:: file_name_suffix,logFile,mFile 
   character(100)               :: mHatFile, gradFile, dataFile, resFile
   type(iterControl_t)			:: CGiter
   character(3)        			:: iterChar
   integer                      :: i,j,iDt,k,ios
   logical                      :: ok

   

   
   
   
if (present(fname)) then
      call read_DCGiterControl(DCGiterControl,fname,ok)
      if (ok) then
         lambda = DCGiterControl%lambda
      end if
else
      call set_DCGiterControl(DCGiterControl)
end if
! initialize the CG control parameters
		call setIterControl(CGiter)
        
        
! initialize the output to log file
   logFile = trim(DCGiterControl%fname)//'_DCG.log'
   open (unit=ioLog,file=logFile,status='unknown',position='append',iostat=ios)
   
   mHat 	=m
   Cm_mHat  =m
   
   m=  multBy_Cm(mHat) 
   call linComb(ONE,m,ONE,m0,m)
   
   JmHat	=d
   Jm0      =d
   dx		=d
   b		=d
   d_Pred	=d
   res		=d
   d_Pred_m0=d

   
call deall_solnVectorMTX(eAll)     
call zero_dataVectorMTX(JmHat)
call zero_dataVectorMTX(b)


   write(*,'(a41,es8.1)') 'The initial damping parameter lambda is ',lambda
   write(ioLog,'(a41,es8.1)') 'The initial damping parameter lambda is ',lambda
   
   
! Compute the predicted data for the current model m
        call Calc_FWD(lambda,d,m,mHat,d_Pred,res,eAll,F,mNorm,rms)
        
        call printf('START',lambda,f,mNorm,rms)
        call printf('START',lambda,f,mNorm,rms,logFile)

        ! initial function call
        nfunc = 1
        
        write(iterChar,'(i3.3)') 0
        
        ! output initial model and responses for later reference
        if (output_level > 1) then
          mFile = trim(DCGiterControl%fname)//'_DCG_'//iterChar//'.rho'
          call write_modelParam(m,trim(mFile))
        end if
        if (output_level > 2) then
          dataFile = trim(DCGiterControl%fname)//'_DCG_'//iterChar//'.dat'
          call write_dataVectorMTX(dHat,trim(dataFile))
        end if       
        

  
       d_Pred_m0=d_Pred


DCG_iter=1
do

	! Compute the right hand side vector (b) for the CG solver.
	! b= (d-dPred)+ J(m-m0)
		    
#ifdef MPI
            JmHat=d
            call zero_dataVectorMTX(JmHat)
	        call Master_job_Jmult(mHat,m,JmHat,eAll) 
#else
	        call Jmult(mHat,m,JmHat,eAll)
#endif
	
	        b=d
	        call linComb(ONE,res,ONE,JmHat,b)
	        call normalize_dataVectorMTX(b,1)  
	        call CG_DS_standard(b,dx,m,d,lambda,CGiter) 
            call normalize_with_dataVecMTX(dx,d,1)
	 
            nfunc = nfunc + CGiter%niter
#ifdef MPI           
            call Master_job_JmultT(m,dx,mHat,eAll)              
#else
	        call JmultT(m,dx,mHat,eAll)
#endif
            nfunc = nfunc + CGiter%niter

	      Cm_mHat=  multBy_Cm(mHat) 
	      mHat=Cm_mHat
	        	     
	     call linComb_modelParam(ONE,m0,ONE,mHat,m)
	! Compute the predicted data for the current model m
         rmsPrev=rms
         call Calc_FWD(lambda,d,m,mHat,d_Pred,res,eAll,F,mNorm,rms)
         nfunc = nfunc + 1
         
	  write(*,'(a25,i5)') 'Completed DCG iteration ',DCG_iter
	  write(ioLog,'(a25,i5)') 'Completed DCG iteration ',DCG_iter
      call printf('with',lambda,f,mNorm,rms)
      call printf('with',lambda,f,mNorm,rms,logFile)
        
         write(iterChar,'(i3.3)') DCG_iter
  
   	  if (output_level > 1) then
   	    mFile = trim(DCGiterControl%fname)//'_DCG_'//iterChar//'.rho'
        call write_modelParam(m,trim(mFile))
      end if
   	  if (output_level > 2) then
   	    mHatFile = trim(DCGiterControl%fname)//'_DCG_'//iterChar//'.prm'
        call write_modelParam(mHat,trim(mHatFile))
      end if
   	  if (output_level > 2) then
   	    dataFile = trim(DCGiterControl%fname)//'_DCG_'//iterChar//'.dat'
        call write_dataVectorMTX(d_Pred,trim(dataFile))
      end if
      ! compute residual for output: res = d-dHat; do not normalize by errors
   	  if (output_level > 2) then
   	    resFile = trim(DCGiterControl%fname)//'_DCG_'//iterChar//'.res'
        call write_dataVectorMTX(res,trim(resFile))
      end if   
      

      
      
      !  test for convergence ...
      if((rms.lt.DCGiterControl%rmsTol).or.(DCG_iter.ge.DCGiterControl%maxIter)) then
         exit
      end if

      DCG_iter=DCG_iter+1
end do
d=d_Pred
   write(*,'(a25,i5,a25,i5)') 'DCG   iterations:',DCG_iter,' function evaluations:',nfunc
   write(ioLog,'(a25,i5,a25,i5)') 'DCG   iterations:',DCG_iter,' function evaluations:',nfunc

    ! cleaning up
   call deall_dataVectorMTX(JmHat)
   call deall_dataVectorMTX(Jm0)
   call deall_dataVectorMTX(b)
   call deall_dataVectorMTX(d_Pred)
   call deall_dataVectorMTX(res)   
   call deall_dataVectorMTX(d_Pred_m0)     
   
   call deall_modelParam(mHat)
   call deall_modelParam(Cm_mHat)

end subroutine DCGsolver
!**********************************************************************
subroutine CG_DS_standard(b,x,m,d,lambda,CGiter)


  type (dataVectorMTX_t), intent(in)	 	::b
  type (dataVectorMTX_t), intent(inout) 	::x
  type(modelParam_t),  intent(in)           ::m
  type (dataVectorMTX_t), intent(in)	    ::d
  real(kind=prec),     intent(in)           ::lambda
  type(iterControl_t), intent(inout)	    :: CGiter
  
  !Local
    type (dataVectorMTX_t)              	:: r,p,Ap
    real(kind=prec)					 	::alpha,beta,r_norm_pre,r_norm,b_norm,error
    integer                          	::i,j,k,ii,iDt
    
     

r=b
p=r
Ap=d
b_norm=dotProd(b,b)
call zero_dataVectorMTX(x)
r_norm=dotProd(r,r)

                
ii = 1
CGiter%rerr(ii) = r_norm/b_norm
 write(ioLog,'(a18)') 'Relative CG-error:'
 write(ioLog,'(a9,i5,a10,es12.5,a10,es12.5)') 'CG-Iter= ',ii,', error = ', r_norm/b_norm, ' Lambda= ', lambda
loop: do while ((CGiter%rerr(ii).gt.CGiter%tol).and.(ii.lt.CGiter%maxIt))

             
! Compute matrix-vector product A*p and save the result in Ap  
       call MultA_DS(p,m,d,lambda,Ap)   

         do i=1,x%nTx 
          do iDt=1, x%d(i)%nDt 
	          r%d(i)%data(iDt)%errorBar= .false.
	          p%d(i)%data(iDt)%errorBar= .false.
	          x%d(i)%data(iDt)%errorBar= .false.
	          Ap%d(i)%data(iDt)%errorBar= .false.
          end do
         end do  
                       
! Compute alpha: alpha= (r^T r) / (p^T Ap)    
       alpha = r_norm/dotProd(p,Ap)
       
! Compute new x: x = x + alpha*p         
       Call scMultAdd_dataVectorMTX(alpha,p,x)  
                                 
! Compute new r: r = r - alpha*Ap   
       Call scMultAdd_dataVectorMTX(-alpha,Ap,r) 

        
                
        r_norm_pre=r_norm
        r_norm=dotProd(r,r)
! Compute beta: beta= r_norm /r_norm_previous           
        beta=r_norm/r_norm_pre
   
! Compute new p: p = r + beta*p    
          call linComb(ONE,r,beta,p,p)
          
       ii=ii+1
       CGiter%rerr(ii) = r_norm/b_norm 
        write(ioLog,'(a9,i5,a10,es12.5,a10,es12.5)') 'CG-Iter= ',ii,', error = ', r_norm/b_norm, ' Lambda= ', lambda
  end do loop

CGiter%niter = ii

! deallocate the help vectors
    call deall_dataVectorMTX(r)
    call deall_dataVectorMTX(p)
    call deall_dataVectorMTX(Ap)
    
end subroutine CG_DS_standard
!**********************************************************************
subroutine MultA_DS(p,m,d,lambda,Ap,JJT_P,CmJTp)
   type(dataVectorMTX_t), intent(in)          ::p
   type(modelParam_t), intent(in)             ::m
   type (dataVectorMTX_t), intent(in)         ::d
   real(kind=prec), intent(in)                ::lambda    
   type(dataVectorMTX_t), intent(out)         ::Ap
   type(dataVectorMTX_t),optional, intent(inout)         ::JJT_P
   type(modelParam_t) ,optional,intent(out)         ::CmJTp
   
!Local parameters
   type(modelParam_t)                      ::JTp,CmJTp_temp
   type(dataVectorMTX_t)                      ::lambdaP,p_temp
   integer                                 ::i,j,k,iDt

JTp		=m
if(present(CmJTp)) then
CmJTp	=m
end if
CmJTp_temp=m
p_temp	=p
lambdaP	=p

!  Mult Cd^(-1/2) p 

         call normalize_with_dataVecMTX(p_temp,d,1)              
! Compute   J^T  Cd^(-1/2) p                   
#ifdef MPI
            call linComb(R_ZERO,d,ONE,p_temp,p_temp) 
            call Master_job_JmultT(m,p_temp,JTp,eAll)
#else
            call JmultT(m,p_temp,JTp,eAll)
#endif
! Compute  Cm  J^T  Cd^(-1/2) p 
            CmJTp_temp= multBy_Cm(JTp) 
            if(present(CmJTp)) then
            CmJTp	=CmJTp_temp
            end if            
! Compute J Cm  J^T  Cd^(-1/2) p = Ap 
Ap=d    
#ifdef MPI
            call Master_job_Jmult(CmJTp_temp,m,Ap,eAll)
#else
            call Jmult(CmJTp_temp,m,Ap,eAll)
#endif

!Normalize: Cd^(-1/2)*Ap
            call normalize_with_dataVecMTX(Ap,d,1)
            !goto 999
            if (present(JJT_P)) then
            JJT_P=Ap
            end if            
            call scMult_dataVectorMTX(lambda,p,lambdaP)

!Add Cd^(-1/2)*Ap*Cd^(-1/2) to lambda*p
         do i=1,lambdaP%nTx
           do iDt=1,lambdaP%d(i)%nDt
             lambdaP%d(i)%data(iDt)%errorBar= .false.
            end do
         end do

           call linComb_dataVectorMTX(ONE,Ap,ONE,lambdaP,Ap)      

999 continue        
!Deallocate help vectors
    call deall_modelParam(JTp)
    !call deall_modelParam(CmJTp)
    call deall_dataVectorMTX(p_temp)
    call deall_dataVectorMTX(lambdaP)
      
end subroutine MultA_DS
!**********************************************************************
subroutine printf(comment,lambda,f,mNorm,rms,logfile)

   ! Compute the full penalty functional F
   ! Also output the predicted data and the EM solution
   ! that can be used for evaluating the gradient
   ! Assuming that the model norm is already scaled by Nmodel

    character(*), intent(in)               :: comment
    real(kind=prec), intent(in)  :: lambda, f, mNorm, rms
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
	write(io_unit,'(a3,es12.5)',advance='no') ' f=',f
	write(io_unit,'(a4,es12.5)',advance='no') ' m2=',mNorm
	write(io_unit,'(a5,f11.6)',advance='no') ' rms=',rms
	write(io_unit,'(a8,es12.5)') ' lambda=',lambda

	! flush(io_unit): this has the effect of flushing the buffer
	if (present(logfile)) then
		close(io_unit)
		open (unit=ioLog,file=logfile,status='old',position='append',iostat=ios)
	end if

   end subroutine printf
!**********************************************************************
subroutine Calc_FWD(lambda,d,m,mHat,d_Pred,res,eAll,F,mNorm,rms)

   ! Compute the full penalty functional F
   ! Also output the predicted data and the EM solution
   ! that can be used for evaluating the gradient
!Input
   real(kind=prec),    intent(in)              :: lambda
   type(dataVectorMTX_t), intent(in)           :: d
   type(modelParam_t), intent(in)              :: m,mHat
!Output   
   real(kind=prec),    intent(out)             :: F, mNorm
   type(dataVectorMTX_t), intent(inout)        :: d_Pred,res
   type(solnVectorMTX_t),  intent(inout)       :: eAll
    real(kind=prec), intent(inout)             :: rms

   !  local variables
   type(dataVectorMTX_t)    :: Nres
   real(kind=prec) :: SS
   integer :: Ndata, Nmodel


   ! initialize d_Pred
   d_Pred = d

   !  compute predicted data for current model parameter m
   !   also sets up forward solutions for all transmitters in eAll
   !   (which is created on the fly if it doesn't exist)
#ifdef MPI
      call Master_Job_fwdPred(m,d_Pred,eAll)
#else
      call fwdPred(m,d_Pred,eAll)
#endif



   ! initialize res
   res = d

   ! compute residual: res = d-d_Pred
   call linComb(ONE,d,MinusONE,d_Pred,res)

   ! normalize residuals, compute sum of squares
   Nres=res
   call normalize_dataVectorMTX(Nres,2)
   SS = dotProd(res,Nres)
   Ndata = countData(res)

   ! compute the model norm
   mNorm = dotProd(mHat,mHat)
   Nmodel = countModelParam(mHat)

   ! penalty functional = sum of squares + scaled model norm
   F = SS/Ndata + (lambda * mNorm/Nmodel)

   ! if required, compute the Root Mean Squared misfit
   	RMS = sqrt(SS/Ndata)
	
!Clean
call deall_dataVectorMTX(Nres)


end subroutine Calc_FWD
!**********************************************************************  
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
   subroutine compute_RMS(m,mhat,d,d_pred,lambda,RMS,res,f,mNorm,ss)

   !Input
    type(modelParam_t), intent(in)           :: m
    type(modelParam_t), intent(in)           :: mhat
   type(dataVectorMTX_t), intent(in)              :: d,d_pred
   real(kind=prec),    intent(in)           :: lambda
   !Output   
   real(kind=prec),    intent(out)          :: RMS,f,mNorm,ss
   type(dataVectorMTX_t),  intent(out)              :: res
   
   !Local
    type(dataVectorMTX_t)                           :: Nres
	!real(kind=prec)                                 :: SS
	Integer                                         :: Ndata,Nmodel
	
    res=d
   call linComb(ONE,d,MinusONE,d_Pred,res)
   Nres=res
   call normalize_dataVectorMTX(Nres,2)
   SS = dotProd(res,Nres)
   Ndata = countData(res)

   ! compute the model norm
   mNorm = dotProd(mHat,mHat)
   Nmodel = countModelParam(mHat)

   ! penalty functional = sum of squares + scaled model norm
   F = SS/Ndata + (lambda * mNorm/Nmodel)

   ! if required, compute the Root Mean Squared misfit
   	RMS = sqrt(SS/Ndata)	
	
	
	
   end subroutine compute_RMS
!**********************************************************************
  subroutine un_normalize_with_dataVecMTX(d_in_out,d,N)
  
     type(dataVectorMTX_t), intent(in)              :: d
     type(dataVectorMTX_t), intent(inout)           :: d_in_out
     integer, optional, intent(in)                  :: N
     integer                                        :: i,j,k,iDt
     
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
!**********************************************************************
!**********************************************************************
end module DCG
