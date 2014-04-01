module DCG

	use math_constants
 	use utilities
    use senscomp
    use main
#ifdef MPI
	Use MPI_main
	use MPI_sub
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
real(kind=prec) :: Desired_rms,k_step_solution
Logical         :: search_min_model=.false.
type(modelParam_t),pointer, dimension(:) :: s_hat
type(modelParam_t),pointer, dimension(:,:) :: JTw_matrix
logical                             :: keep_solution
type (dataVectorMTX_t)              	:: x_previous,CSEM_d,MT_d
Integer                                 :: Ndata_CSEM,Ndata_MT,start_kstep
   ! type(EMsolnMTX_t),save    :: eAll

Contains
!**********************************************************************
   subroutine set_DCGiterControl(DCGiterControl)
   !  at present this is a simple routine to setup for
   !  iteration control; allowing parameters to be read from
   !   a file would be desireable.

   type(DCGiterControl_t), intent(inout)	:: DCGiterControl

     ! maximum number of iterations in one call to iterative solver
     DCGiterControl%maxIter = 200
     ! convergence criteria: return from solver if rms < rmsTol
     DCGiterControl%rmsTol  = 1.05
     ! initial value of lambda 
     DCGiterControl%lambda = 100.     
     ! model and data output file name
     DCGiterControl%fname = 'Modular_DCG'     
     
     
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
   CGiter%tol = 10E-8
   CGiter%niter = 0
   allocate(CGiter%rerr(0:CGiter%maxit))

   end subroutine setIterControl
!**********************************************************************

  subroutine DCGsolver(d,m0,m,lambda)
  
  ! Subroutine to solve the inverse problem in data space using conjugate gradients (CG)  
   
   type(dataVectorMTX_t), intent(inout)		       ::d
   !  m0 is prior model parameter
   type(modelParam_t), intent(in)		       ::m0
   !  m is solution parameter ... on input m contains starting guess
   type(modelParam_t), intent(inout)	       ::m
   !  lambda is regularization parameter
   real(kind=prec) , intent(inout)							   ::lambda
   
!  local variables
   type(dataVectorMTX_t)			:: dHat, b,dx,d_Pred,res,Nres,JmHat,Jm0,d_Pred_m0
   type(modelParam_t)			:: mHat,CmJTd,Cm_mHat
   real(kind=prec)		  		:: value,rms_old,F,mNorm,rms
   integer						:: iter, ndata,DS_iter,CG_iter
   character(100)       		:: file_name_suffix
   type(iterControl_t)			:: CGiter
   character(3)        			:: iterChar
   integer                          	::i,j,iDt,k

   

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
   start_kstep=2
call zero_dataVectorMTX(JmHat)
call zero_dataVectorMTX(b)

! initialize the CG control parameters
		call setIterControl(CGiter)
   
open(130,file='DCG.log')
! Compute the predicted data for the current model m
        call Calc_FWD(lambda,d,m,d_Pred,res,eAll,F,mNorm,rms)
        file_name_suffix='Start_DCG'
         call write_output_files(1,d_Pred,m,file_name_suffix) 
        call printf('DCG_Start ',lambda,rms,mNorm,F)
        Desired_rms=rms/2.0
        
        write(130,'(a20)',advance='no') trim('DCG_Start ')//':'
 		write(130,'(a5,f11.6)',advance='no') ' rms=',rms
	    write(130,'(a4,es12.6)',advance='no') ' m2=',mNorm
	    write(130,'(a3,es12.6)',advance='no') ' F=',F    
		write(130,'(a8,f11.6)') ' lambda=',lambda     
  
       d_Pred_m0=d_Pred


  open(110,file='Lanczos_CG_DS.dat')
  open(150,file='x_sub.dat')
do DS_iter=1,5
	! Compute the right hand side vector (b) for the CG solver.
	! b= (d-dPred)+ J(m-m0)
	
	        if (DS_iter .gt. 1 )then	    
#ifdef MPI
            JmHat=d
            call zero_dataVectorMTX(JmHat)
	        call Master_job_Jmult(mHat,m,JmHat,eAll) 
#else
	        call Jmult(mHat,m,JmHat,eAll)
#endif
	
	        end if
	        b=d
	        call linComb(ONE,res,ONE,JmHat,b)
	        call normalize_dataVectorMTX(b,1)  
	        
	         rms_old=rms
	         call CG_DS(b,dx,m,m0,d,lambda,CGiter,d_Pred_m0,rms,mhat)
           !  call Lanczos_DS (b,m,m0,d,d_Pred,lambda,mhat,res,CGiter,DS_iter,rms,Jm0)
	     !  call Multi_Trans_DS (b,dx,m,m0,d,lambda,mhat,res,CGiter,DS_iter,rms)
	      
          ! call Lanczos_CG_DS(b,dx,m,m0,d,lambda,CGiter,DS_iter,rms)
	         
			 if (rms .le. 1.05 .or. abs(rms -rms_old) .lt. 0.01 ) then
	           !goto 999
			 end if
			 
	         !goto 10
             goto 5
	        call normalize_with_dataVecMTX(dx,d,1)

	 
#ifdef MPI           
	                call Master_job_JmultT(m,dx,mHat,eAll)              
#else
	                call JmultT(m,dx,mHat,eAll)
#endif

	      Cm_mHat=  multBy_Cm(mHat) 
	      mHat=Cm_mHat
	       
5    continue  	     
	     call linComb_modelParam(ONE,m0,ONE,mHat,m)
	! Compute the predicted data for the current model m
   
	   rms_old=rms
         call Calc_FWD(lambda,d,m,d_Pred,res,eAll,F,mNorm,rms)
           if (abs(rms -rms_old) .lt. 0.2) then
             !lambda=lambda/2.0
           end if
		   !lambda=lambda/2.0
    ! Write output model and data files
    file_name_suffix='DCG_MPI'
    call write_output_files(DS_iter,d_Pred,m,file_name_suffix)     
     ! Print output Information on the screen
    write(iterChar,'(i3.3)') DS_iter
    call printf('DCG_Iter '//iterChar,lambda,rms,mNorm,F)
    
        write(130,'(a20)',advance='no') 'DCG_Iter '//iterChar//':'
 		write(130,'(a5,f11.6)',advance='no') ' rms=',rms
	    write(130,'(a4,es12.6)',advance='no') ' m2=',mNorm
	    write(130,'(a3,es12.6)',advance='no') ' F=',F    
		write(130,'(a8,f11.6)',advance='no') ' lambda=',lambda  
        write(130,'(a16,i5)') ' # of CG iter.=', CGiter%niter


   10 continue        
 ! Clean temp vectors
end do
999 continue
d=d_Pred
close(130)
close(110)
close(150) 
end subroutine DCGsolver
!****************************************************************************************
subroutine Lanczos_CG_DS(b,x,m,m0,d,lambda,CGiter,DS_iter,start_RMS)


  type (dataVectorMTX_t), intent(in)	 	::b
  type (dataVectorMTX_t), intent(inout) 	::x
  type(modelParam_t),  intent(inout)       ::m
  type(modelParam_t),  intent(in)       ::m0  
  type (dataVectorMTX_t), intent(inout)	    ::d
  real(kind=prec),     intent(inout)       ::lambda,start_RMS
  type(iterControl_t), intent(inout)	:: CGiter
  integer,     intent(in)               :: DS_iter
  character(3)         					::iterChar
     character(100)       		:: file_name_suffix
  
  !Local
    type (dataVectorMTX_t)              	:: r,v,p,v_previous,Av,v_hat,d_Pred,res,w,x_n,x_opt
    type (dataVectorMTX_t)              	:: u,alpha_u,Au,beta_u,v_temp
    type(modelParam_t)			        :: ATu,beta_v,mHat,Cm_mHat,CmJTv,m_temp
    real(kind=prec)					 	::alpha,beta,delta,sigma,omega,mue,mue_previous,omega_previous
    real(kind=prec),pointer, dimension(:) :: beta_vec,delta_vec
    type (dataVectorMTX_t),pointer, dimension(:) :: v_vec,x_vec,p_vec
    real(kind=prec)						:: b_norm,r_norm,lambda1,lambda2,delta_hat,beta_zero,F,mNorm,rms,rms_old,sigma11,sigma22,sigma33,rms_1,target_rms,rms_opt
    integer                          	::i,j,iDt,k,ii,i_cg,i_lambda,jj
     type(solnVectorMTX_t)              :: eAll_temp


     target_rms=start_rms/2.0
     if (target_rms .lt. 1.0 ) then
       target_rms=1.05
     end if  
 mHat=m
 Cm_mHat=m
 CmJTv=m
 d_Pred=b
res=b
w=b
v_temp=b
x_n=b
x_opt=b
rms_opt=100000.00
         eAll_temp=eAll
         m_temp=m
    call zero_dataVectorMTX(d_Pred)  
    call zero_dataVectorMTX(res) 
    call zero_dataVectorMTX(x_opt) 
    call zero(mhat)
  call zero(Cm_mHat)
     call zero(w)
     call zero(v_temp)
     call zero(CmJTv)
allocate(beta_vec(0:20),delta_vec(0:20),v_vec(-1:20),x_vec(1:20),p_vec(1:20))

do i=-1,20
    v_vec(i)=b
    call zero_dataVectorMTX(v_vec(i))  
end do    
    
do i=1,20
    x_vec(i)=b
    p_vec(i)=b
    call zero_dataVectorMTX(x_vec(i))  
    call zero_dataVectorMTX(p_vec(i))  
end do    
    

call zero_dataVectorMTX(x)    
r=b
beta=sqrt(dotProd(r,r))
beta_vec(0)=beta
v=b
call scMult(1/beta,r,v)
v_vec(0)=v
p=r  




call zero_dataVectorMTX(v_vec(-1))
sigma=beta_vec(0)
omega_previous=R_ZERO
mue_previous=ONE

Av=b
v_hat=b
b_norm=dotProd(b,b)



do ii=0,10

    ! begin Lanczos step
      !compute delta = <Av,v>
      write(110,*)'Lanczos step ', ii
      !call MultA_DS(v_vec(ii),m,d,R_ZERO,Av,CmJTp) 
     
      
  !       do i=1,b%nTx 
  !        do iDt=1, b%d(i)%nDt 
	 !         Av%d(i)%data(iDt)%errorBar= .false.
  !        end do
  !       end do 
  !        
  !       if (ii .gt. 0) then
		!     do i=1,b%nTx 
	 !         do iDt=1, b%d(i)%nDt 
		!          v_vec(ii-1)%d(i)%data(iDt)%errorBar= .false.
	 !         end do
	 !        end do 		 
		!   call linComb (ONE,Av,-beta_vec(ii),v_vec(ii-1),w)
		! else
		!    w=Av
  !       end if     
  !       
  !
		! alpha=(dotProd(v_vec(ii),w))
		! delta_vec(ii)=alpha
		! do i=1,b%nTx 
  !        do iDt=1, b%d(i)%nDt 
	 !         w%d(i)%data(iDt)%errorBar= .false.
  !        end do
  !       end do 
  !       
  !       
		!call linComb (ONE,w,-alpha,v_vec(ii),v_vec(ii+1))
		!  
		!    
		!        call zero(v_temp)
	 !		   do jj=0,ii-1
		!		  sigma11= dotProd(v_vec(jj),v_vec(ii+1))
		!		  sigma22= dotProd(v_vec(jj),v_vec(jj))
		!		  sigma33=sigma11/sigma22
		!		  
		!		 		do i=1,d%nTx 
		!		          do iDt=1, d%d(i)%nDt 
		!			          v_vec(jj)%d(i)%data(iDt)%errorBar= .false.
		!		          end do
		!		         end do
		!		    call linComb (ONE,v_temp,sigma33,v_vec(jj),v_temp)
  !              end do
		!		 
		!		 		do i=1,d%nTx 
		!		          do iDt=1, d%d(i)%nDt 
		!			          v_vec(ii+1)%d(i)%data(iDt)%errorBar= .false.
		!		          end do
		!		         end do
		!	         call linComb (ONE,v_vec(ii+1),-ONE,v_temp,v_vec(ii+1))         
  !       
  !      beta_vec(ii+1)=sqrt(dotProd(v_vec(ii+1),v_vec(ii+1)))
  !      call scMult(1.0/beta_vec(ii+1),v_vec(ii+1),v_vec(ii+1))      
         

      delta_vec(ii)=(dotProd(Av,v_vec(ii)))
      !compute v_hat = Av - delta v - beta v_previous
           do i=1,Av%nTx
            do iDt=1,Av%d(i)%nDt
             do j=1,Av%d(i)%data(iDt)%nSite
               do k=1,Av%d(i)%data(iDt)%nComp
                      v_hat%d(i)%data(iDt)%value(k,j)=  Av%d(i)%data(iDt)%value(k,j)-(delta_vec(ii)*v_vec(ii)%d(i)%data(iDt)%value(k,j))-(beta_vec(ii)*v_vec(ii-1)%d(i)%data(iDt)%value(k,j))
              end do
            end do
           end do                                       
           end do
        !compute beta       
         beta_vec(ii+1)=sqrt(dotProd(v_hat,v_hat))
         !v_previous=v
         call scMult(ONE/beta_vec(ii+1),v_hat,v_vec(ii+1))
            lambda2=3.0
		    Lambda1=10**(lambda2)
            rms=1000.
    do i_lambda=1,20         
       !rest the initial vectors for each lambda  
         
         call zero_dataVectorMTX(x_vec(i_lambda))   
         r=b
         sigma=beta_vec(0)
         p=r
         p_vec(i_lambda)=r   

         omega_previous=R_ZERO
         mue_previous=ONE
        ! begin the CG- iterates
        do i_cg=0,ii
                delta_hat=delta_vec(i_cg)+lambda1
                 mue=ONE/(delta_hat-(omega_previous/ mue_previous))     
                 mue_previous=mue
                 omega=(beta_vec(i_cg+1)*mue)**2
                 omega_previous=omega
                 sigma=-beta_vec(i_cg+1)*mue*sigma
                 !Update x, r and p
		           do i=1,b%nTx
		            do iDt=1,b%d(i)%nDt
		             do j=1,b%d(i)%data(iDt)%nSite
		               do k=1,b%d(i)%data(iDt)%nComp
		                      x_vec(i_lambda)%d(i)%data(iDt)%value(k,j)=  x_vec(i_lambda)%d(i)%data(iDt)%value(k,j)+(mue*p_vec(i_lambda)%d(i)%data(iDt)%value(k,j))
		                      r%d(i)%data(iDt)%value(k,j)=  sigma*v_vec(i_cg+1)%d(i)%data(iDt)%value(k,j)
		                      p_vec(i_lambda)%d(i)%data(iDt)%value(k,j)=  r%d(i)%data(iDt)%value(k,j)+(omega*p_vec(i_lambda)%d(i)%data(iDt)%value(k,j))
		              end do
		            end do
		           end do                                       
		        end do 
		        r_norm=sqrt(dotProd(r,r))
                write(110,'(i5,ES12.5,2x,ES12.5,2x,ES12.5,2x,ES12.5,2x,ES12.5,2x,ES12.5,2x,ES12.5)')i_cg,lambda1,mue,omega,sigma,delta_hat,beta_vec(i_cg+1),beta_vec(0)    

        end do ! CG ietrates

         !goto 999
            call zero_dataVectorMTX(x_n)  
            x_n=x_vec(i_lambda)
	        call normalize_with_dataVecMTX(x_n,d,1)

	 
#ifdef MPI           
	                call Master_job_JmultT(m,x_n,mHat,eAll)              
#else
	                call JmultT(m,x_n,mHat,eAll)
#endif

	      Cm_mHat=  multBy_Cm(mHat) 
	      mHat=Cm_mHat
	       
	     call zero(m)
	     call linComb_modelParam(ONE,m0,ONE,mHat,m)
	! Compute the predicted data for the current model m
         rms_old=rms

         call Calc_FWD(lambda,d,m,d_Pred,res,eAll,F,mNorm,rms)  
          
    
         
         write(110,'(2i5,ES12.5,2x,f11.6,2x,ES12.5,2x,ES12.5,2x,ES12.5)')ii,i_lambda,Lambda1,rms,f,abs(sigma)/beta_vec(0),r_norm/b_norm
         
         write(i_lambda+200,'(i5,2x,f11.6,2x,ES12.5,2x,ES12.5)')ii,Lambda1,abs(sigma),abs(sigma)/beta_vec(0)
         if (rms .lt. rms_opt) then
             rms_opt=rms
             x_opt=x_vec(i_lambda)
          end if   
         !if (rms .gt. rms_old) then 
         !       m=m_temp
         !       eAll=eAll_temp
         !       goto 120
         ! end if       
         !if (rms .lt. target_rms) then
         !        m=m_temp
         !       eAll=eAll_temp
         !       goto 999
         !end if
         
         write(iterChar,'(i3.3)') ii
         file_name_suffix='DCG_MPI_'//iterChar
         call write_output_files(i_lambda,d_Pred,m,file_name_suffix)          
        
         !call linComb_modelParam(ONE,m0,R_ZERO,mHat,m)
         !call Calc_FWD(lambda,d,m,d_Pred,res,eAll,F,mNorm,rms_1)  
         
                m=m_temp
                eAll=eAll_temp
               lambda2=lambda2-0.25
               Lambda1=10**(lambda2)
    end do !lambda
120 continue     
end do  ! Lanczos steps
999 continue  
    m=m_temp
    eAll=eAll_temp
    Lambda=Lambda1
    x=x_opt
    write(10,'(1i5)') i_lambda
    
   

    

deallocate(beta_vec,delta_vec,v_vec,x_vec,p_vec)
    
    
end subroutine Lanczos_CG_DS 
!****************************************************************************************
subroutine Multi_Trans_DS(b,x,m,m0,d,lambda,mhat,res,CGiter,DS_iter,rms)


  type (dataVectorMTX_t), intent(in)	 	::b
  type (dataVectorMTX_t), intent(inout) 	::x
  type(modelParam_t),  intent(inout)     ::m
   type(modelParam_t),  intent(in)      :: m0
  type(modelParam_t),  intent(inout)       ::mhat
  type (dataVectorMTX_t), intent(inout)	    ::d
   type(dataVectorMTX_t),intent(inout)			:: res
  real(kind=prec),     intent(inout)       ::lambda,rms
  type(iterControl_t), intent(inout)	:: CGiter
  integer,intent(in)					::DS_iter
  
  
  !Local 
    type (dataVectorMTX_t),pointer, dimension(:) :: u_hat
    type(modelParam_t)                       :: JTu
    Integer                                  :: k_step,i_k_step,nTx,ii,jj,INFO,counter,counter1,ll,kk,i_lambda,i,idt
    real(kind=prec)                          :: beta,sum,mu_LU, lambda_LU,lambda_temp,F,mNorm,rms_old,sigma11
     real(kind=prec),pointer, dimension(:,:)   ::sts 
    character(3)         					 ::iterChar,iterChar1
    character(100)       		:: file_name_suffix
    real(kind=prec)	,pointer, dimension(:,:)  :: sts_temp,L_matrix, U_matrix
    real(kind=prec)	,pointer, dimension(:)    :: b_sub,y_sub,x_sub
    real(kind=prec)	,pointer, dimension(:)    ::  AP,sTs_b,uTr
    type(modelParam_t)			              :: q_temp,mhat_temp,m_temp
    type(dataVectorMTX_t)				          ::Jm,u_hat_temp,b_u,Jm_ub,residual
      type(dataVectorMTX_t)			:: d_Pred

  k_step=5
  nTx=b%nTx
   lambda_temp=lambda   
mhat=m
q_temp=m
mhat_temp=m
  Jm=d
  u_hat_temp=b
  b_u=b   
  Jm_ub=b 
  residual=b
  JTu=m0
  d_Pred=d
call zero(q_temp)
call zero(mhat)
call zero(mhat_temp)   
call zero(b_u) 
call zero(Jm_ub) 
 call zero(residual) 
 
 
  call zero(u_hat_temp)   
  call zero(Jm) 
  
  
  allocate(u_hat(k_step+1))
  allocate(JTw_matrix(k_step,nTx))
  allocate(s_hat(nTx),sts(nTx,nTx))
  allocate(sts_temp(nTx,nTx),L_matrix(nTx,nTx),U_matrix(nTx,nTx))
  allocate(b_sub(nTx),y_sub(nTx),x_sub(nTx),uTr(nTx))
  allocate(AP((nTx)*((nTx)+1)/2),sTs_b(nTx))
  
  


  
	  do ii=1,nTx
	  	s_hat(ii)=m
	  	call zero(s_hat(ii))
	  end do
	  do ii=1,k_step+1
	  	  u_hat(ii)=b
	  	  call zero(u_hat(ii))
	  end do
	  u_hat(1)=b
	  
	  
	 do jj=1,k_step 
	  do ii=1,nTx
	  	  JTw_matrix(jj,ii)=m
	  	  call zero(JTw_matrix(jj,ii))
	  end do
	end do  
	  	  	
	  
sts=R_zero
sts_temp=R_zero
L_matrix=R_zero
U_matrix=R_zero
y_sub=R_zero
x_sub=R_zero
b_sub=R_zero
sTs_b=R_zero


! 1-normalize each sub data vector [ u_hat(1), u_hat(2),...,u_hat(nTx)] with its norm: i.e.  u_hat(1)/||u_hat(1)||
      !call normalize_with_dataVecMTX(u_hat(1),d,1) 
      do ii=1,nTx
	   beta=sqrt(dotProd(u_hat(1)%d(ii),u_hat(1)%d(ii)))
	   call scMult(1/beta,u_hat(1)%d(ii),u_hat(1)%d(ii))
	  end do  
	  
 do i_k_step=1, k_step 


     
! 2- Normlize u_hat with the data error Cd^(-1/2):
!	 Mult u_hat  Cd^(-1/2) 



       
! The data vector for all transmitters is ready to pass it to JmultT.  

!3- Compute JT u_hat: the parallel subroutine 'Master_job_JmultT' returens back model parameters vectors 
!   corresponding to each transmitter; model(1,...,nTx)
!   However, the serial version returens only one model vector.                   
!   Compute   J^T  Cd^(-1/2) u_hat    
	  do ii=1,nTx
	  	call zero(s_hat(ii))
	  end do
	  call zero(JTu)
	    !  call linComb(R_ZERO,d,ONE,u_hat(i_k_step),u_hat(i_k_step))           
#ifdef MPI
            call Master_job_JmultT(m,u_hat(i_k_step),JTu,eAll,s_hat)
#else
            call JmultT(m,u_hat(i_k_step),JTu,eAll)
#endif



     sum=R_Zero
	  do ii=1,nTx
	  	  JTw_matrix(i_k_step,ii)=s_hat(ii)
		  sum=sum+sqrt(dotProd(JTw_matrix(1,ii),JTw_matrix(1,ii)))
	  end do
	   write(6,*)'############ Initial Lmabda ################'
	   write(6,*) sum/nTx
	  
	  
	  

 !4-   Compute  Cm J^T  Cd^(-1/2) u_hat and save them in JTw_matrix
    do ii=1,nTx
     ! JTw_matrix(i_k_step,ii)= multBy_Cm(JTw_matrix(i_k_step,ii)) 
   end do
   
 !5 - Make the symetric matrix sTs for the current i_k_step 
        do kk=1,nTx
          do jj=1,nTx
            sts(kk,jj)=dotProd(JTw_matrix(i_k_step,kk),JTw_matrix(i_k_step,jj))
          end do
        end do     
  
 !6 - Add Lambda to diag(sTs)   
             sts_temp=R_Zero
             sts_temp=sts
			 do ii=1,nTx
			    sts_temp(ii,ii)=sts(ii,ii)+1000
			 end do
			   
			do jj=1,nTx
			  do ii=1,jj
               AP(ii + (jj-1)*jj/2) = sts_temp(ii,jj) 
			  end do
			end do  
!7-  Preforme Cholesky decomposition on (sTs +Lambda I) matrix 
            call DPPTRF( 'U', nTx, AP, INFO )
          do jj=1,nTx
           do ii=1,jj
              L_matrix(ii,jj) = AP(ii + (jj-1)*jj/2) 
			  end do
		  end do
			
             write(6,*)'############ (sTs ) matrix ################'
			  do ii=1,nTx
			    write(6,'(19f15.2)')(sts(ii,jj),jj=1,nTx)
			 end do 
			 
			 write(6,*)'############ (sTs +Lambda I) matrix ################'
			  do ii=1,nTx
			    write(6,'(19f15.2)')(sts_temp(ii,jj),jj=1,nTx)
			 end do 

			 write(6,*)'############### L matrix #############'
			 do ii=1,nTx
			    write(6,'(19f15.2)')(L_matrix(ii,jj),jj=1,nTx)
			 end do 
! 8- make the right hand side to solve the projected problem:
!    b_sub=uhat(1,...,nTx)*b
 
			  do ii=1,nTx
			    b_sub(ii)=(dotProd(u_hat(i_k_step)%d(ii),b%d(ii)))
			 end do  
			 
			 write(6,*)'############### right hand side (b_sub vector) #############'
			 do ii=1,nTx
			    write(6,*) b_sub(ii)
			 end do 			 
!9- solve the problem, save solution in b_sub
  call DPPTRS( 'U', nTx, 1, AP, b_sub, nTx, INFO )	
             write(6,*)'############### Solution for the projected system #############'		 
 			 do ii=1,nTx
			    write(6,*) b_sub(ii)
            end do  

! 10- Model update for the projected system:
!  m_hat=JT b= JT u_hat *b_sub
! Model update 
 call zero(mhat)
 call zero(q_temp)
		 do ii=1, nTx
			  call scMult_modelParam (b_sub(ii),JTw_matrix(i_k_step,ii),q_temp)
			  call linComb_modelParam(ONE,mhat,ONE,q_temp,mhat)
		 end do

! 11- smooth the model update m_hat		  
	!mhat= multBy_Cm(mhat)
              
 
! 12-  compute the next u_hat
		  ! beta * u_hat   = J mhat - u_hat *(sTs) *b_sub 
		  !     k+1     k+1                k             k
		
		!  12a - compute J m  
#ifdef MPI
		   call Master_job_Jmult(mhat,m,Jm,eAll)
#else
		   call Jmult(mhat,m,Jm,eAll)
#endif 
! normalize Jm with Cd^1/2

         !call normalize_with_dataVecMTX(Jm,d,1) 
         
		 do ii=1,b%nTx 
          do jj=1, b%d(ii)%nDt 
	          Jm%d(ii)%data(jj)%errorBar= .false.
          end do
         end do 	
		!  12b-  (sTs) *b_sub 
		
		       do ii=1,nTx
		         sum=R_zero
		          do jj=1,nTx
		            sum=sum+(sts(ii,jj)*b_sub(jj))
		          end do
		           sTs_b(ii)=sum
		       end do     
		 ! 12c-  u_hat * sTs_b
		  do ii=1,nTx
		     call scMult(sTs_b(ii),u_hat(i_k_step)%d(ii),u_hat_temp%d(ii))
		 end do

	
          
   	    do ii=1,b%nTx 
          do jj=1, b%d(ii)%nDt 
	          u_hat_temp%d(ii)%data(jj)%errorBar= .false.
          end do
         end do 
                 	
		  call linComb (ONE,Jm,MinusONE,u_hat_temp,u_hat(i_k_step+1))
		!  Call normalize_with_dataVecMTX(u_hat(i_k_step+1),d,1) 

	  do ii=1,nTx
	   beta=sqrt(dotProd(u_hat(i_k_step+1)%d(ii),u_hat(i_k_step+1)%d(ii)))
	   call scMult(1/beta,u_hat(i_k_step+1)%d(ii),u_hat(i_k_step+1)%d(ii))
	  end do  
	  

 
          		  
   
 end do 
 ! Check orthogonalty
  write(6,*)'############### Check orthogonalty #############'
 do ii=1,nTx 
  do kk=1,k_step
    do jj=1,k_step
      beta=(dotProd(u_hat(kk)%d(ii),u_hat(jj)%d(ii)))
      write(6,*)ii,kk,jj,beta
    end do
  end do
 end do 


end subroutine Multi_Trans_DS  
!****************************************************************************************
subroutine update_sys(k_step,m,d,beta_kstep_plus_1,u,T_matrix,q,Au_vec)

integer,intent(in)                      			    :: k_step
real(kind=prec),intent(inout),dimension(:,:)    	    :: T_matrix
type (modelParam_t),intent(inout), dimension(:) 	     :: q
type (dataVectorMTX_t),intent(inout), dimension(:) 	     :: u
real(kind=prec),intent(inout)                            :: beta_kstep_plus_1
type(dataVectorMTX_t), intent(inout), dimension(:)                      :: Au_vec

  type(modelParam_t),  intent(in)       ::m
  type (dataVectorMTX_t), intent(in)	    ::d
    
!Local
    real(kind=prec),pointer,dimension(:,:)  :: T_matrix_Temp
    type (modelParam_t),pointer, dimension(:) 	    :: q_temp
    type(dataVectorMTX_t)                      :: Au,w
    type (modelParam_t)             	    :: q_kstep_plus_1
    real(kind=prec)                         :: alpha,beta_kstep,sigma11
    integer                                 :: ii,i,iDt,jj
     

q_kstep_plus_1= q(1)
Au=d
w=d



call zero(q_kstep_plus_1)
call zero(Au)
call zero(w)




 
 do ii=k_step,k_step

		 call MultA_JJT_DS(u(ii),m,d,Au,q_kstep_plus_1)
		 Au_vec(ii)=Au
		 
	  do jj=1,d%nTx 
	     s_hat(jj)=multBy_Cm(s_hat(jj))
	  	 JTw_matrix(ii,jj)=s_hat(jj)
	  end do			 
		 
		 
	     do i=1,d%nTx 
          do iDt=1, d%d(i)%nDt 
	          Au%d(i)%data(iDt)%errorBar= .false.
          end do
         end do 
		 
		 
		 call linComb (ONE,Au,-beta_kstep_plus_1,u(ii-1),w)
		 alpha=(dotProd(u(ii),w))
		
		 do i=1,d%nTx 
          do iDt=1, d%d(i)%nDt 
	          w%d(i)%data(iDt)%errorBar= .false.
          end do
         end do 
         
         
		call linComb (ONE,w,-alpha,u(ii),u(ii+1))
		
	 		    do jj=1,ii
	 		         sigma11= dotProd(u(ii+1),u(jj))
	 		         
				 		do i=1,d%nTx 
				          do iDt=1, d%d(i)%nDt 
					          u(ii+1)%d(i)%data(iDt)%errorBar= .false.
					          u(jj)%d(i)%data(iDt)%errorBar= .false.
				          end do
				         end do
				         
			         call linComb (ONE,u(ii+1),-sigma11,u(jj),u(ii+1))
	 		    end do
! Update T_Matrix and q      
	 		    
	  T_matrix(k_step,k_step)=alpha
      T_matrix(k_step-1,k_step)=beta_kstep_plus_1
      T_matrix(k_step,k_step-1)=beta_kstep_plus_1
      q(k_step)=q_kstep_plus_1

        beta_kstep_plus_1=sqrt(dotProd(u(ii+1),u(ii+1)))
        call scMult(1.0/beta_kstep_plus_1,u(ii+1),u(ii+1))
              
        write(6,*) k_step, 'Beta= ',beta_kstep,  beta_kstep_plus_1          
        write(6,*) k_step, 'Alpha= ',alpha       
      	    		


        

            





 end do
		
		

 

         

      

      

 end subroutine update_sys
!****************************************************************************************
subroutine scale_dHat(scaling_factor_CSEM_temp,scaling_factor_MT_temp,d,dHat)
		real(kind=prec) ,intent(in)                                  ::scaling_factor_CSEM_temp,scaling_factor_MT_temp 
		type (dataVectorMTX_t), intent(inout)	 	                 :: dHat
		type (dataVectorMTX_t), intent(in)	 	                 :: d
		
		Integer                     :: ii,jj,nTx,iDt,j,k
		
		
              nTx=dHat%nTx
	   
	    	  do jj=1,nTx
			    do iDt=1, dHat%d(jj)%nDt 
			       if (txDict(jj)%Tx_type=='MT') then
					 dHat%d(jj)%data(iDt)%value=(scaling_factor_MT_temp)*dHat%d(jj)%data(iDt)%value
					 dHat%d(jj)%data(iDt)%error=scaling_factor_MT_temp*d%d(jj)%data(iDt)%error
				   elseif (txDict(jj)%Tx_type=='CSEM') then
				       	 do j=1,d%d(jj)%data(iDt)%nSite
	                        do k=1,d%d(jj)%data(iDt)%nComp	
							     IF (abs(dHat%d(jj)%data(iDt)%value(k,j)-d%d(jj)%data(iDt)%value(k,j)) .gt. d%d(jj)%data(iDt)%error(k,j)) then
				                     dHat%d(jj)%data(iDt)%value(k,j)=scaling_factor_CSEM_temp*dHat%d(jj)%data(iDt)%value(k,j)
								 end if
								 
							end do
                         end do							
				   end if
			     end do
			  end do


end subroutine scale_dHat
!****************************************************************************************
subroutine scale_q_Matrix(k_step,scaling_factor_CSEM_temp,scaling_factor_MT_temp,b,m,q,beta_1)
		integer,intent(in)                      			     :: k_step
		real(kind=prec) ,intent(in)                              ::scaling_factor_CSEM_temp,scaling_factor_MT_temp 
		type (modelParam_t),intent(inout), dimension(:) 	     :: q
		type (dataVectorMTX_t), intent(in)	 	                 :: b
		type(modelParam_t),  intent(in)                          :: m
		real(kind=prec),optional ,intent(in)                              :: beta_1
		real(kind=prec) temp_norm,temp_norm_MT,temp_norm_CSEM,temp_norm_total,scaling_factor_MT,scaling_factor_CSEM,scaling_factor_q
        type(modelParam_t)			::  JTd_MT,JTd_CSEM,JTd_total,JTd_temp
		Integer                     :: ii,jj,nTx
		
  
  nTx=b%nTx
  
JTd_MT=m
JTd_CSEM=m
JTd_total=m	
JTd_temp=m


     call zero(JTd_MT)   
     call zero(JTd_CSEM)   
     call zero (JTd_total)
     
	 
  do ii=1, 1
 	  do jj=1,nTx
	              call scMult(beta_1,JTw_matrix(ii,jj),JTd_temp)
				  if (txDict(jj)%Tx_type=='MT') then
				       write(50,*)ii,jj,temp_norm,txDict(jj)%period,'MT'
					   call linComb_modelParam(ONE,JTd_MT,ONE,JTd_temp,JTd_MT) ! get JTd_MT
				  elseif (txDict(jj)%Tx_type=='CSEM') then
				       write(50,*)ii,jj,temp_norm,'CSEM'
					   call linComb_modelParam(ONE,JTd_CSEM,ONE,JTd_temp,JTd_CSEM) ! get JTd_CSEM
                  end if	
				  
				  call linComb_modelParam(ONE,JTd_total,ONE,JTd_temp,JTd_total) ! get JTd_total
	  end do
  end do
   write(110,*) 'scale q 2'
  
	  temp_norm_MT = sqrt(dotProd(JTd_MT,JTd_MT))
	  temp_norm_CSEM = sqrt(dotProd(JTd_CSEM,JTd_CSEM))
	  temp_norm_total=sqrt(dotProd(JTd_total,JTd_total))

	  
	  			       scaling_factor_MT =    scaling_factor_MT_temp !0.5*(temp_norm_total/temp_norm_MT) !scaling_factor_MT_temp    !1.0/temp_norm_MT !scaling_factor_MT_temp !(scaling_factor*temp_norm_total)/temp_norm_MT !1.0/ sqrt(dotProd(JTd_MT,JTd_MT)) !scaling_factor !(scaling_factor*temp_norm_total)/temp_norm_MT    
                       scaling_factor_CSEM =   scaling_factor_CSEM_temp !0.5*(temp_norm_total/temp_norm_CSEM) ! scaling_factor_CSEM_temp !10.0/temp_norm_CSEM !scaling_factor_CSEM_temp !1.0 ! 0.5*temp_norm_total/temp_norm_CSEM !scaling_factor !(2.0-scaling_factor) !((1.0-scaling_factor)*temp_norm_total)/temp_norm_CSEM  !(1.0-scaling_factor) !1.0/ sqrt(dotProd(JTd_CSEM,JTd_CSEM)) !(1.0-scaling_factor) !

  
  
  do ii=1, k_step
     call zero(JTd_MT)   
     call zero(JTd_CSEM)   
     call zero (JTd_total)	
 	  do jj=1,nTx
 				 temp_norm = sqrt(dotProd(JTw_matrix(ii,jj),JTw_matrix(ii,jj)))
				 !call scMult(1.0/temp_norm,JTw_matrix(ii,jj),JTw_matrix(ii,jj))
				 
				  if (txDict(jj)%Tx_type=='MT') then
				       write(50,*)ii,jj,temp_norm,txDict(jj)%period,'MT'
					   call linComb_modelParam(ONE,JTd_MT,ONE,JTw_matrix(ii,jj),JTd_MT) ! get JTd_MT
				  elseif (txDict(jj)%Tx_type=='CSEM') then
				       write(50,*)ii,jj,temp_norm,'CSEM'
					   call linComb_modelParam(ONE,JTd_CSEM,ONE,JTw_matrix(ii,jj),JTd_CSEM) ! get JTd_CSEM
                  end if	
				  
				  call linComb_modelParam(ONE,JTd_total,ONE,JTw_matrix(ii,jj),JTd_total) ! get JTd_total
	  end do

	  
	  
	  temp_norm_MT = sqrt(dotProd(JTd_MT,JTd_MT))
	  temp_norm_CSEM = sqrt(dotProd(JTd_CSEM,JTd_CSEM))
	  temp_norm_total=sqrt(dotProd(JTd_total,JTd_total))
	  write(55,*)ii,temp_norm_MT,temp_norm_total,scaling_factor_MT_temp
	 

				!if (ii == 1) then	 
				
					  call scMult(scaling_factor_MT,JTd_MT,JTd_MT)
					  call scMult(scaling_factor_CSEM,JTd_CSEM,JTd_CSEM)
	  
	  temp_norm_MT = sqrt(dotProd(JTd_MT,JTd_MT))
	  temp_norm_CSEM = sqrt(dotProd(JTd_CSEM,JTd_CSEM))
	  temp_norm_total=sqrt(dotProd(JTd_total,JTd_total))
	  write(60,*)ii,temp_norm_MT,temp_norm_total,scaling_factor_MT_temp
	  
  
					  call linComb_modelParam(ONE,JTd_MT,ONE,JTd_CSEM,q(ii))
		              !scaling_factor_q =1.0/ sqrt(dotProd(q(ii),q(ii)))
		              !call scMult(scaling_factor_q,q(ii),q(ii))
	  
	  
end do
 
  
  
end subroutine scale_q_Matrix
!****************************************************************************************
subroutine Lanczos_DS(b,m,m0,d,d_Pred_Orig,lambda,mhat,res,CGiter,DS_iter,rms,Jm0)


  type (dataVectorMTX_t), intent(in)	 	 :: b
  type(modelParam_t),  intent(inout)     :: m
  type(modelParam_t),  intent(in)        :: m0
  type(modelParam_t),  intent(inout)     :: mhat
  type(dataVectorMTX_t),  intent(inout)	 :: res 
  type (dataVectorMTX_t), intent(inout)	 :: d,d_Pred_Orig
  type (dataVectorMTX_t), intent(inout)     ::	Jm0
  type(iterControl_t), intent(inout)	 :: CGiter
  integer,intent(in)					 :: DS_iter
  real(kind=prec),intent(inout)			 :: rms,lambda

  
  !Local

    type (modelParam_t),pointer, dimension(:)  :: q,v
    type (dataVectorMTX_t),pointer, dimension(:)  :: u,u_recent,Au_vec
    real(kind=prec)	,pointer, dimension(:,:)   :: T_matrix
    real(kind=prec)					 	       :: Jm0_norm,beta_1,beta_kstep_plus_1,scaler_term
    integer                          	       ::Ndata,i,iDt,i_lambda,k_step,i_search,i_sub_search,ii,jj,i_sub_big_search,i_sub_small_search,fwd_calls,max_k_steps,min_k_steps,nTx,i_alpha_serach

   type(dataVectorMTX_t)			:: JmHat,d_Pred,b_start,d_start,u_temp,Au,dHat,d_Pred_approx,d_Pred_approx_temp
   character(100)       		:: file_name_suffix
   type(modelParam_t)			:: m_c,m_l,m_r,m_start,m_temp,JTd_MT,JTd_CSEM,JTd_total
   type(dataVectorMTX_t)			:: res_c,res_r,res_l,Nres
   real(kind=prec)              :: tol_rms,SS,start_lambda,step_lambda,scaling_factor1,old_rms_alpha,F,mNorm,lambda_c,lambda_r,lambda_l,mNorm_c,mNorm_r,mNorm_l,rms_c,rms_r,rms_l,step_size,sub_step_size,rms_temp,lambda11 
   real(kind=prec)              :: old_rms_MT,res_norm,opt_rms,scaling_factor_CSEM_temp,scaling_factor_MT_temp,scaling_factor_MT_opt,scaling_factor_CSEM_opt,RMS_MT,RMS_CSEM
   logical                      :: go_right,go_left,central,make_big_step,make_small_step,update_proj_sys
   character(3)         		:: iterChar,iterChar1
   character(100)       		:: modelFile,dataFile
   type(dataVectorMTX_t)			:: w,u_kstep
   type(solnVectorMTX_t)        :: eAll_temp
       real :: T1,T2, Seconds 
	       logical                  :: found_opt_MT_scale

	  

	
	
min_k_steps=start_kstep
max_k_steps=50

nTx=b%nTx
    k_step=min_k_steps
    
m_start=m
b_start=b
d_start=d
dHat=b
res_c=b
d_Pred_approx=b
d_Pred_approx_temp=b
m_temp=m
update_proj_sys=.false.
eAll_temp=eAll
u_temp=d
    call zero(u_temp)
	call zero (dHat)
	call zero (d_Pred_approx)	
	call zero (d_Pred_approx_temp)
    call zero (res_c)

allocate(T_matrix(max_k_steps,max_k_steps))
allocate(q(max_k_steps),v(k_step+1),u_recent(3),u(1:max_k_steps+1),Au_vec(max_k_steps+1))
allocate(JTw_matrix(max_k_steps,nTx))

T_matrix=R_zero
do ii=1,max_k_steps
  q(ii)= m
  call zero(q(ii))
end do
do ii=1,max_k_steps+1
  u(ii)= b
  Au_vec(ii)=b
  call zero(u(ii))
  call zero(Au_vec(ii))
end do

do ii=1,3
  u_recent(ii)=b
end do
rms_temp=1000.00

i_search=1
i_sub_search=1
i_sub_big_search=1
i_sub_small_search=1
fwd_calls=0
Lambda_c=10**(lambda)





  1 continue  

! scaling_factor=1.0
! write(105,*)'Before', sqrt(dotProd(b_start,b_start))
! do jj=1,nTx 
! if (txDict(jj)%Tx_type=='MT') then
!   call scMult(scaling_factor,b_start%d(jj),b_start%d(jj))
! elseif (txDict(jj)%Tx_type=='CSEM') then
!    scaling_factor1=1-scaling_factor
!    call scMult(scaling_factor1,b_start%d(jj),b_start%d(jj))
! end if  
!end do
!  write(105,*)'After', sqrt(dotProd(b_start,b_start))
!  
  
!call Arnold_bidiag_JJT(b,m,m0,d,k_step,DS_iter,lambda,T_matrix,q,beta_1)
 if (update_proj_sys) then 
     call update_sys (k_step,m_start,d_start,beta_kstep_plus_1,u,T_matrix,q,Au_vec)
 else
      call bidiag_JJT (b_start,m_start,d_start,k_step,DS_iter,lambda,T_matrix,q,beta_1,beta_kstep_plus_1,u,Au_vec)
     ! call Arnold_bidiag_JJT(b_start,m_start,m0,d_start,k_step,DS_iter,lambda,T_matrix,q,beta_1,beta_kstep_plus_1,u)
end if
   
   write(iterChar1,'(i3.3)') DS_iter
    if (DS_iter .eq. 1 .and. k_step .eq. min_k_steps) then
       open(10,file='L_serach.dat',STATUS = 'unknown')
    else
       open(10,file='L_serach.dat',status='unknown',position='append')
   end if
   
   
!if (k_step==min_k_steps) then
!             write(10,*) 'Start RMSs'
!             call compute_RMS(m,d,d_Pred_Orig,Lambda_c,RMS_c,res_c)
!			 write(10,'(4f10.3)') rms_c,RMS_MT,RMS_CSEM,Lambda_c
!end if
			 

   
   

   
   
   
   
   
   
 
!   write(10,*)' Start Search for Residual',k_step 
!   lambda11=4.0
!   Lambda_c=10**(lambda11)
!   old_rms_alpha=10000.0 
!	call cpu_time  ( T1 )
! do i_lambda=1, 21	 
!    call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c)
!	scaler_term=-ONE*beta_kstep_plus_1*k_step_solution
!	call scMult_dataVectorMTX(scaler_term,u(k_step+1),res_c)
!	res_norm=(dotProd(res_c,res_c))/(dotProd(b,b))
!	write(10,'(i5,f10.3,x,E9.3)') i_lambda,Lambda_c,res_norm
!	if  (res_norm .gt. 0.01) then
!	   lambda11=lambda11-0.5
!     lambda=lambda11
!	 goto 380
!	end if 
!	
!   lambda11=lambda11-0.25
!   Lambda_c=10**lambda11
! end do
! 380 continue


!
!
!
! 		    write(10,*)' Start Search for Weighting MT' 
!				 scaling_factor_CSEM_temp=1.0
!				 scaling_factor_MT_temp=0.0
!				  old_rms_alpha=10000.0
!				  
!             do i_alpha_serach=1,21
!				     !call scale_q_Matrix(k_step,scaling_factor_CSEM_temp,scaling_factor_MT_temp,b,m,q)
!				     call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c,Au_vec,dHat)
!                   !call Calc_FWD(Lambda_c,d,m_c,d_Pred,res_c,eAll,F,mNorm_c,rms_c)                 					 
!                	call un_normalize_with_dataVecMTX(dHat,d,1)
!					 	do i=1,d%nTx 
!				          do iDt=1, d%d(i)%nDt 
!					          dHat%d(i)%data(iDt)%errorBar= .false.
!				          end do
!				         end do	
!		  call scale_dHat(scaling_factor_CSEM_temp,scaling_factor_MT_temp,dHat)
!		  call linComb (ONE,d_Pred_Orig,ONE,dHat,d_Pred_approx)
!		  call compute_RMS(m_c,d,d_Pred_approx,Lambda_c,RMS_c,res_c)
!				  if (rms_c .gt. old_rms_alpha) then
!						scaling_factor_MT_temp=scaling_factor_MT_temp-0.1
!					    call scale_q_Matrix(k_step,scaling_factor_CSEM_temp,scaling_factor_MT_temp,b,m,q)
!						goto 530
!                end if					 
!					 old_rms_alpha = rms_c
!					 write(10,'(i5,3f10.3)') i_alpha_serach,scaling_factor_MT_temp,rms_c,Lambda_c
!				     scaling_factor_MT_temp=scaling_factor_MT_temp+0.1
!				 end do  
!30 continue	
!                  call scale_q_Matrix(k_step,scaling_factor_CSEM_temp,scaling_factor_MT_temp,b,m,q)
!
!
! 		    write(10,*)' Start Search for Weighting CSEM' 
!				 scaling_factor_CSEM_temp=0.0
!				 !scaling_factor_MT_temp=1.0
!				  old_rms_alpha=10000.0
!				  
!             do i_alpha_serach=1,21
!				     !call scale_q_Matrix(k_step,scaling_factor_CSEM_temp,scaling_factor_MT_temp,b,m,q)
!				     call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c,Au_vec,dHat)
!                   !call Calc_FWD(Lambda_c,d,m_c,d_Pred,res_c,eAll,F,mNorm_c,rms_c)                 					 
!                	call un_normalize_with_dataVecMTX(dHat,d,1)
!					 	do i=1,d%nTx 
!				          do iDt=1, d%d(i)%nDt 
!					          dHat%d(i)%data(iDt)%errorBar= .false.
!				          end do
!				         end do	
!		  call scale_dHat(scaling_factor_CSEM_temp,scaling_factor_MT_temp,dHat)
!		  call linComb (ONE,d_Pred_Orig,ONE,dHat,d_Pred_approx)
!		  call compute_RMS(m_c,d,d_Pred_approx,Lambda_c,RMS_c,res_c)
!				  if (rms_c .gt. old_rms_alpha) then
!						scaling_factor_CSEM_temp=scaling_factor_CSEM_temp-0.1
!					    call scale_q_Matrix(k_step,scaling_factor_CSEM_temp,scaling_factor_MT_temp,b,m,q)
!						goto 520
!                end if					 
!					 old_rms_alpha = rms_c
!					 write(10,'(i5,3f10.3)') i_alpha_serach,scaling_factor_CSEM_temp,rms_c,Lambda_c
!				     scaling_factor_CSEM_temp=scaling_factor_CSEM_temp+0.1
!				 end do  
!20 continue	
!
!
!
!
!
!
!
!	   call cpu_time  ( T2 )
!     Seconds=T2-T1 
!    write(10,'(a11,f10.3,a8)') 'End: total=', Seconds, ' Seconds'
!	  
!	       call scale_q_Matrix(k_step,scaling_factor_CSEM_temp,scaling_factor_MT_temp,b,m,q)
!
	  
	  
	  
	        !if (k_step .lt. max_k_steps) then
	  		!		 write(10,*) 'Increase the size of the projected system by adding additional K step'
			!	     !call Calc_FWD(lambda,d,m_start,d_Pred,res,eAll,F,mNorm,rms)
			!	     eAll=eAll_temp
			!	     !fwd_calls=fwd_calls +1
		    !        k_step=k_step+1
		    !        update_proj_sys=.true.
		    !        close(10)
		    !        !rms_temp=rms_c
		    !        m_temp=m_c
		    !        goto 1
	        ! end if
			  
	  
	  
!  write(10,*)' Start Search for Min model norm with 5% of the target RMS' 
!  lambda11=4.0
!  Lambda_c=10**(lambda11)
!  old_rms_alpha=10000.0 
!	call cpu_time  ( T1 )
!do i_lambda=1, 21	 
!   call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c,Au_vec,dHat)
!	call un_normalize_with_dataVecMTX(dHat,d,1)
!					 	do i=1,d%nTx 
!				          do iDt=1, d%d(i)%nDt 
!					          dHat%d(i)%data(iDt)%errorBar= .false.
!				          end do
!				         end do	
!	      call linComb (ONE,d_Pred_Orig,ONE,dHat,d_Pred_approx)
!		  call compute_RMS(m_c,d,d_Pred_approx,Lambda_c,RMS_c)
!		  tol_rms=1.2*opt_rms
!     if (rms_c .le. tol_rms) then
!	   lambda11=lambda11+0.25
!	   Lambda_c=10**lambda11
!	   Lambda=Lambda11
!	    rms_c=old_rms_alpha
!	    goto 395
!	 end if
!	 old_rms_alpha = rms_c		  
!	write(10,'(i5,4f10.3)') i_lambda,Lambda_c,rms_c,sqrt(dotProd(mhat,mhat)),1.05*rms_c
!  lambda11=lambda11-0.25
!  Lambda_c=10**lambda11
!end do 
!
!395 continue
!	   call cpu_time  ( T2 )
!      Seconds=T2-T1 
!     write(10,'(a11,f10.3,a8)') 'End: total=', Seconds, ' Seconds'	  
!	  


     		
! 
! 
! 

!
!   call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c)
!   call Calc_FWD(Lambda_c,d,m_c,d_Pred,res_c,eAll,F,mNorm_c,rms_c)
  
!Lambda_c=10**(lambda)
!

!				     call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c,Au_vec,dHat)               					 
!                     call un_normalize_with_dataVecMTX(dHat,d,1)
!					 	do i=1,d%nTx 
!				          do iDt=1, d%d(i)%nDt 
!					          dHat%d(i)%data(iDt)%errorBar= .false.
!				          end do
!				         end do	
!		  call linComb (ONE,d_Pred_Orig,ONE,dHat,d_Pred_approx)
!		  call compute_RMS(m_c,d,d_Pred_approx,Lambda_c,RMS_c,res_c,RMS_MT,RMS_CSEM)
!		  write(10,'(a5,5f10.3)') 'Before Scaling',rms_c,RMS_MT,RMS_CSEM
!		  
!		  
!






 


                 
 
 

 



	
 




 



    !call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c,Au_vec,dHat)
 	!call un_normalize_with_dataVecMTX(dHat,d,1)
 	!				 	do i=1,d%nTx 
 	!			          do iDt=1, d%d(i)%nDt 
 	!				          dHat%d(i)%data(iDt)%errorBar= .false.
 	!			          end do
 	!			         end do	
    !      call scale_dHat(scaling_factor_CSEM_temp,scaling_factor_MT_temp,dHat)			 
 	!      call linComb (ONE,d_Pred_Orig,ONE,dHat,d_Pred_approx)
 	!	  call compute_RMS(m_c,d,d_Pred_approx,Lambda_c,RMS_c,res_c)
 	!	  mNorm_c=sqrt(dotProd(m_c,m_c))
 	!	  JmHat=dHat


!		    write(10,*)' Start Search for Weighting MT Full'
!         		 
!				 old_rms_alpha=10000.0
!				 step_lambda = 10 ** (1 / 3.0)
!				 scaling_factor_CSEM_temp=1.0
!				 scaling_factor_CSEM_opt=1.0
!				 scaling_factor_MT_temp=0.0
!				 
!               do i_alpha_serach=1,21
!				     call scale_q_Matrix(k_step,scaling_factor_CSEM_temp,scaling_factor_MT_temp,b,m,q)
!				      call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c)
!                     call Calc_FWD(Lambda_c,d,m_c,d_Pred,res_c,eAll,F,mNorm_c,rms_c)  
!                  if (rms_c .gt. old_rms_alpha) then
!						scaling_factor_MT_temp=scaling_factor_MT_temp-0.1
!						scaling_factor_MT_opt=scaling_factor_MT_temp
!						goto 54
!                  end if					 
!					 old_rms_alpha = rms_c
!					 write(10,'(i5,3f10.3)') i_alpha_serach,scaling_factor_MT_temp,rms_c,Lambda_c
!				     scaling_factor_MT_temp=scaling_factor_MT_temp+0.1
!				 end do  
! 54 continue	
!
!
!                  
                            !call scale_q_Matrix(k_step,scaling_factor_CSEM_opt,scaling_factor_MT_opt,b,m,q)

!					 Jm0_norm = sqrt(dotProd(Jm0,Jm0))
!					 write(10,*)'Jm0_norm =', Jm0_norm
!            Lambda_c=10000.0
!			lambda11=4
!			old_rms_alpha=10000.0 
!	do i_lambda=1, 21				
!			!call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c,Au_vec,dHat)
!            !JmHat=d
!            !call zero_dataVectorMTX(JmHat)
!	        !call Master_job_Jmult(mHat,m,JmHat,eAll) 	
!			!					do i=1,d%nTx 
!			!					  do iDt=1, d%d(i)%nDt 
!			!						  JmHat%d(i)%data(iDt)%errorBar= .false.
!			!					   end do
!			!					 end do	    
!			!	    call linComb (ONE,d_Pred_Orig,ONE,JmHat,d_Pred_approx)								 
!	        !         write(iterChar,'(i3.3)') 1
!  			!	     file_name_suffix='Hybrid_FWD_JmHat_'//iterChar
!		    !         call write_output_files(DS_iter,d_Pred_approx,m_c,file_name_suffix) 
!			!		 call compute_RMS(m_c,d,d_Pred_approx,Lambda_c,RMS_r,res_c)
!            !    	 !write(10,*)'RMS FWD_JmHatt=', rms_c,Lambda_c		 
!
!					 
!		    call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c,Au_vec,dHat)
!			call un_normalize_with_dataVecMTX(dHat,d,1)
!								do i=1,d%nTx 
!								  do iDt=1, d%d(i)%nDt 
!									  dHat%d(i)%data(iDt)%errorBar= .false.
!									  d_Pred_Orig%d(i)%data(iDt)%errorBar= .false.
!								   end do
!								 end do	
!				    call linComb (ONE,d_Pred_Orig,ONE,dHat,d_Pred_approx)					  
!					  
!	                 write(iterChar,'(i3.3)') 1
!  				     file_name_suffix='Hybrid_FWD_approx_'//iterChar
!		             call write_output_files(DS_iter,d_Pred_approx,m_c,file_name_suffix) 
!                     call compute_RMS(m_c,d,d_Pred_approx,Lambda_c,RMS_c,res_c)						 
!			if (rms_c .gt. old_rms_alpha) then
!			   lambda11=lambda11+0.25
!			   Lambda_c=10**lambda11
!				     call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c)
!                     call Calc_FWD(Lambda_c,d,m_c,d_Pred,res_c,eAll,F,mNorm_c,rms_c) 
!			   			
!			   !Check if the Desired_rms is reached
!			     if (rms_c .le. Desired_rms) then
!			    !d_Pred=d_Pred_approx
!				Lambda=Lambda11
!				rms=RMS_c
!				opt_rms=RMS_c
!				mNorm_c=sqrt(dotProd(mhat,mhat))
!		           		  m=m_c
!						  lambda=dlog10(lambda_c)
!						  res=res_c
!						  mNorm=mNorm_c
!			      call linComb_modelParam(ONE,m,MinusONE,m0,mhat)				
!				   write(10,*) '1- Exit this Iteration with:'
!				   write(10,8551)DS_iter, rms,mNorm,10**lambda,k_step,(k_step*2)
!				   write(10,*)
!				    Desired_rms=rms/2
!					start_kstep=2			 
!				   goto 999
!				  end if 
!				 if (k_step .lt. max_k_steps) then
!				      write(10,*) 'Increase the size of the projected system by adding additional K step'
!				     eAll=eAll_temp
!				     !fwd_calls=fwd_calls +1
!		             k_step=k_step+1
!		             update_proj_sys=.true.
!		             close(10)
!		             !rms_temp=rms_c
!		             m_temp=m_c
!		             goto 1
!				 else
!					d_Pred=d_Pred_approx
!					Lambda=Lambda11
!					rms=RMS_c
!					opt_rms=RMS_c
!					mNorm_c=sqrt(dotProd(mhat,mhat))
!							  m=m_c
!							  lambda=dlog10(lambda_c)
!							  res=res_c
!							  mNorm=mNorm_c
!					  call linComb_modelParam(ONE,m,MinusONE,m0,mhat)				
!					   write(10,*) '1- Exit this Iteration with:'
!					   write(10,8551)DS_iter, rms,mNorm,10**lambda,k_step,(k_step*2)
!					   write(10,*)
!						Desired_rms=rms/2
!						start_kstep=2			 
!					   goto 999				 
!                 end if			 
!			 end if
!			old_rms_alpha = rms_c
!			
!                      
!					  
!				 !    call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c)
!                 !    call Calc_FWD(Lambda_c,d,m_c,d_Pred,res_c,eAll,F,mNorm_c,rms_c)  
!				!  
!	             !   write(iterChar,'(i3.3)') 1
!  				 !   file_name_suffix='Hybrid_FWD_Full_'//iterChar
!		         !   call write_output_files(DS_iter,d_Pred,m_c,file_name_suffix) 		
!                 !    call compute_RMS(m_c,d,d_Pred,Lambda_c,RMS_c,res_c)					 
!					  write(10,'(3f10.3,i5)') Lambda_c,rms_c,sqrt(dotProd(mhat,mhat)),k_step					  
!
!		 lambda11=lambda11-0.25
!		 Lambda_c=10**lambda11
!end do 
!                      !Lambda_c=562.0
!				      call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c)
!                      call Calc_FWD(Lambda_c,d,m_c,d_Pred,res_c,eAll,F,mNorm_c,rms_c)  
!		           		  m=m_c
!						  lambda=dlog10(lambda_c)
!						  rms=rms_c
!						  res=res_c
!						  mNorm=mNorm_c
!			      call linComb_modelParam(ONE,m,MinusONE,m0,mhat)
!				  fwd_calls=fwd_calls +1
!				   write(10,*) '1- Exit this Iteration with:'
!				   write(10,8551)DS_iter, rms,mNorm,10**lambda,k_step,(k_step*2)+fwd_calls
!				   write(10,*)
!				    Desired_rms=rms/2
!					start_kstep=5
!
!				  !close(10)
!				  goto 999
!					  
					  
					  
					  
					  
					  
!			scaling_factor_MT_temp=1.0
!		    scaling_factor_CSEM_temp=1.0
!			scaling_factor_MT_opt=1.0
!		    scaling_factor_CSEM_opt=1.0	
!
!			
!		    write(10,*)' Start Search for Weighting CSEM approx',k_step 
!				 scaling_factor_CSEM_temp=0.0
!				 scaling_factor_MT_temp=2.0
!				 scaling_factor_MT_opt=1.0
!				  old_rms_alpha=10000.0	
!                  old_rms_MT=10000.0	
!                  found_opt_MT_scale=.false.				  
!            do i_alpha_serach=1,21
!				     call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c,Au_vec,dHat)               					 
!                     call un_normalize_with_dataVecMTX(dHat,d,1)
!					 	do i=1,d%nTx 
!				          do iDt=1, d%d(i)%nDt 
!					          dHat%d(i)%data(iDt)%errorBar= .false.
!				          end do
!				         end do			  
!		  call scale_dHat(scaling_factor_CSEM_temp,scaling_factor_MT_temp,d,dHat)
!		  call linComb (ONE,d_Pred_Orig,ONE,dHat,d_Pred_approx)
!		  call compute_RMS(m_c,d,d_Pred_approx,Lambda_c,RMS_c,res_c,RMS_MT,RMS_CSEM)
!		       
!			   if (.not. found_opt_MT_scale) then
!		         if (rms_MT .gt. old_rms_MT) then
!                  found_opt_MT_scale= .true.
!                  scaling_factor_MT_temp=scaling_factor_MT_temp+0.1
!				  scaling_factor_MT_opt=scaling_factor_MT_temp
!                 end if
!			   end if
!			   
!				  if (RMS_c .gt. old_rms_alpha ) then
!				   !write(10,'(i5,5f10.3)') i_alpha_serach,scaling_factor_CSEM_temp,rms_c,RMS_MT,RMS_CSEM,Lambda_c
!						scaling_factor_CSEM_temp=scaling_factor_CSEM_temp-0.1
!                        scaling_factor_CSEM_opt=scaling_factor_CSEM_temp
!						goto 520
!                  end if					 
!					 old_rms_alpha = RMS_c
!					 old_rms_MT = RMS_MT
!					 write(10,'(i5,6f10.3)') i_alpha_serach,scaling_factor_MT_temp,scaling_factor_CSEM_temp,rms_c,RMS_MT,RMS_CSEM,Lambda_c
!				     scaling_factor_CSEM_temp=scaling_factor_CSEM_temp+0.1
!					 
!					 if (.not. found_opt_MT_scale) then
!					 scaling_factor_MT_temp=scaling_factor_MT_temp-0.1
!					 end if
!					 
!				 end do  
! 520 continue	 
		
                     if (scaling_factor_CSEM_opt .gt. 2.0) scaling_factor_CSEM_opt=2.0
                     if (scaling_factor_CSEM_temp .gt. 2.0) scaling_factor_CSEM_temp=2.0	


					 
			
			



 
 
 
 
 
 				!     call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c)
                !     call Calc_FWD(Lambda_c,d,m_c,d_Pred,res_c,eAll,F,mNorm_c,rms_c,RMS_MT,RMS_CSEM) 
				!	 write(10,*)'Before scaling',rms_c,RMS_MT,RMS_CSEM 
                ! 
                !  call scale_q_Matrix(k_step,scaling_factor_CSEM_opt,scaling_factor_MT_opt,b,m,q,beta_1)	
				!   
 				!     call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c)
                !     call Calc_FWD(Lambda_c,d,m_c,d_Pred,res_c,eAll,F,mNorm_c,rms_c,RMS_MT,RMS_CSEM) 				   
				!   write(10,*)'After scaling',rms_c,RMS_MT,RMS_CSEM 
				 
				 
				 
				 
			

!		 write(10,*)' Start Search for Lambda approxomated',k_step 
!		 lambda11=4.0
!		 Lambda_c=10**(lambda11)
!		 old_rms_alpha=10000.0 
!			call cpu_time  ( T1 )
!
!			
!		do i_lambda=1, 21	 
!           call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c,Au_vec,dHat)
!			call un_normalize_with_dataVecMTX(dHat,d,1)
!								do i=1,d%nTx 
!								  do iDt=1, d%d(i)%nDt 
!									  dHat%d(i)%data(iDt)%errorBar= .false.
!									  d_Pred_Orig%d(i)%data(iDt)%errorBar= .false.
!								   end do
!								 end do	
!				     call linComb (ONE,d_Pred_Orig,ONE,dHat,d_Pred_approx)					  
!                    call compute_RMS(m_c,d,d_Pred_approx,Lambda_c,RMS_c,res_c,RMS_MT,RMS_CSEM)						 
!			if (rms_c .gt. old_rms_alpha) then
!			    !write(10,'(i5,6f10.3)') i_lambda,scaling_factor_MT_temp,scaling_factor_CSEM_temp,rms_c,RMS_MT,RMS_CSEM,Lambda_c	
!			   lambda11=lambda11+0.25
!			   Lambda_c=10**lambda11
!			   Lambda=Lambda11
!				rms_c=old_rms_alpha
!				opt_rms=old_rms_alpha
!				mNorm_c=sqrt(dotProd(mhat,mhat))
!						        !call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c,Au_vec,dHat)
!			                    !call un_normalize_with_dataVecMTX(dHat,d,1)
!								!do i=1,d%nTx 
!								!  do iDt=1, d%d(i)%nDt 
!								!	  dHat%d(i)%data(iDt)%errorBar= .false.
!								!  end do
!								! end do	
!								!  call scale_dHat(scaling_factor_CSEM_temp,scaling_factor_MT_temp,d,dHat)
!				                !  call linComb (ONE,d_Pred_Orig,ONE,dHat,d_Pred_approx)
!				                !  call compute_RMS(m_c,d,d_Pred_approx,Lambda_c,RMS_c,res_c,RMS_MT,RMS_CSEM)
!				goto 390
!			 end if
!			 old_rms_alpha = rms_c
!           write(10,'(i5,6f10.3)') i_lambda,scaling_factor_MT_temp,scaling_factor_CSEM_temp,rms_c,RMS_MT,RMS_CSEM,Lambda_c			 
!		 lambda11=lambda11-0.25
!		 Lambda_c=10**lambda11
!		end do 
!
!       390 continue
		
		
		
		
!		    write(10,*)' Start Search for Weighting MT approx' 
!				 scaling_factor_MT_temp=0.0
!				 scaling_factor_CSEM_temp=1.0
!                scaling_factor_CSEM_opt=1.0
!				 old_rms_alpha=10000.0
!				  
!            do i_alpha_serach=1,21
!					call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c,Au_vec,dHat)
!					call un_normalize_with_dataVecMTX(dHat,d,1)
!								do i=1,d%nTx 
!								  do iDt=1, d%d(i)%nDt 
!									  dHat%d(i)%data(iDt)%errorBar= .false.
!									  d_Pred_Orig%d(i)%data(iDt)%errorBar= .false.
!								   end do
!								 end do	
!					 call scale_dHat(scaling_factor_CSEM_temp,scaling_factor_MT_temp,d,dHat)
!				     call linComb (ONE,d_Pred_Orig,ONE,dHat,d_Pred_approx)					  
!                    call compute_RMS(m_c,d,d_Pred_approx,Lambda_c,RMS_c,res_c,RMS_MT,RMS_CSEM)		
!				  if (RMS_MT .gt. old_rms_alpha) then
!						scaling_factor_MT_temp=scaling_factor_MT_temp-0.1
!					    scaling_factor_MT_opt=scaling_factor_MT_temp
!						goto 530
!               end if					 
!					 old_rms_alpha = RMS_MT
!					 write(10,'(i5,6f10.3)') i_lambda,scaling_factor_MT_temp,scaling_factor_CSEM_temp,rms_c,RMS_MT,RMS_CSEM,Lambda_c
!				     scaling_factor_MT_temp=scaling_factor_MT_temp+0.1
!				 end do  
!530 continue	
                    ! if (scaling_factor_MT_opt .gt. 2.0) scaling_factor_MT_opt=2.0
                    ! if (scaling_factor_MT_temp .gt. 2.0) scaling_factor_MT_temp=2.0					 
		            !
		
		

					 
		
		
	                 
	



		 write(10,*)' Start Search for Lambda approxomated',k_step 
		 lambda11=4.0
		 Lambda_c=10**(lambda11)
		 old_rms_alpha=10000.0 
			call cpu_time  ( T1 )
			scaling_factor_MT_temp=1.0
		    scaling_factor_CSEM_temp=1.0
			
		do i_lambda=1, 21	 
            call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c,Au_vec,dHat)
			call un_normalize_with_dataVecMTX(dHat,d,1)
								do i=1,d%nTx 
								  do iDt=1, d%d(i)%nDt 
									  dHat%d(i)%data(iDt)%errorBar= .false.
									  d_Pred_Orig%d(i)%data(iDt)%errorBar= .false.
								   end do
								 end do	
					 call scale_dHat(scaling_factor_CSEM_temp,scaling_factor_MT_temp,d,dHat)
				     call linComb (ONE,d_Pred_Orig,ONE,dHat,d_Pred_approx)					  
                     call compute_RMS(m_c,d,d_Pred_approx,Lambda_c,RMS_c,res_c,RMS_MT,RMS_CSEM)						 
			if (rms_c .gt. old_rms_alpha) then
            lambda11=lambda11+0.25
            Lambda_c=10**lambda11
            !lambda=lambda11
            call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c,Au_vec,dHat)
			call un_normalize_with_dataVecMTX(dHat,d,1)
								do i=1,d%nTx 
								  do iDt=1, d%d(i)%nDt 
									  dHat%d(i)%data(iDt)%errorBar= .false.
									  d_Pred_Orig%d(i)%data(iDt)%errorBar= .false.
								   end do
								 end do	
					 call scale_dHat(scaling_factor_CSEM_temp,scaling_factor_MT_temp,d,dHat)
				     call linComb (ONE,d_Pred_Orig,ONE,dHat,d_Pred_approx)					  
                     call compute_RMS(m_c,d,d_Pred_approx,Lambda_c,RMS_c,res_c,RMS_MT,RMS_CSEM)					
                  write(iterChar,'(i3.3)') k_step
  				  file_name_suffix='Hybrid_'//iterChar
		          call write_output_files(DS_iter,d_Pred_approx,m_c,file_name_suffix)               
                goto 380
              end if  
            !    !write(10,'(i5,6f10.3)') i_lambda,scaling_factor_MT_temp,scaling_factor_CSEM_temp,rms_c,RMS_MT,RMS_CSEM,Lambda_c	
			!   lambda11=lambda11+0.25
			!   Lambda_c=10**lambda11
			!   Lambda=Lambda11
			!               call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c,Au_vec,dHat)
			!	rms_c=old_rms_alpha
			!	opt_rms=old_rms_alpha
			!	mNorm_c=sqrt(dotProd(mhat,mhat))
			!	goto 380
			! end if
			 old_rms_alpha = rms_c
            write(10,'(i5,6f10.3)') i_lambda,scaling_factor_MT_temp,scaling_factor_CSEM_temp,rms_c,RMS_MT,RMS_CSEM,Lambda_c			 
		 lambda11=lambda11-0.25
		 Lambda_c=10**lambda11
		end do 

        380 continue
       call cpu_time  ( T2 )
       Seconds=T2-T1 
      write(10,'(a11,f10.3,a8)') 'End: total=', Seconds, ' Seconds'
		
		
!    write(10,*)' Start Search for Lambda Full FWD' 
!  
!  lambda11=4.0
!  Lambda_c=10**(lambda11)
!  old_rms_alpha=10000.0	 
!	call cpu_time  ( T1 )
!do i_lambda=1, 21	 
!  call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c)
!  !call Calc_FWD(Lambda_c,d,m_c,d_Pred,res_c,eAll,F,mNorm_c,rms_c)
!				   
!    !if (rms_c .gt. old_rms_alpha) then
!	!  lambda11=lambda11+0.25
!	!  lambda=lambda11
!	!   Lambda_c=10**lambda11
!	!    rms_c=old_rms_alpha
!	!	opt_rms=old_rms_alpha
!	!  goto 400
!	!end if
!	! old_rms_alpha = rms_c
!   write(10,'(i5,6f10.3)') i_lambda,scaling_factor_MT_temp,scaling_factor_CSEM_temp,rms_c,RMS_MT,RMS_CSEM,Lambda_c	
!  
!  lambda11=lambda11-0.25
!  Lambda_c=10**lambda11
!end do 
!400 continue	
       call cpu_time  ( T2 )
       Seconds=T2-T1 
      write(10,'(a11,f10.3,a8)') 'End: total=', Seconds, ' Seconds'		

      
 !		    write(10,*)' Start Search for Weighting CSEM Full',k_step 
 !        		 
	!			 old_rms_alpha=10000.0
	!			 step_lambda = 10 ** (1 / 3.0)
	!			 scaling_factor_CSEM_temp=0.0
	!			 scaling_factor_MT_temp=1.0
	!			 scaling_factor_MT_opt=1.0
	!			 
 !              do i_alpha_serach=1,21
	!			     call scale_q_Matrix(k_step,scaling_factor_CSEM_temp,scaling_factor_MT_temp,b,m,q,beta_1)
	!			     call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c)
 !                    call Calc_FWD(Lambda_c,d,m_c,d_Pred,res_c,eAll,F,mNorm_c,rms_c,RMS_MT,RMS_CSEM)  
 !                 if (rms_c .gt. old_rms_alpha) then
	!					scaling_factor_CSEM_temp=scaling_factor_CSEM_temp-0.1
	!					scaling_factor_CSEM_opt=scaling_factor_CSEM_temp
	!					goto 53
 !                 end if					 
	!				 old_rms_alpha = rms_c
	!				 write(10,'(i5,6f10.3)') i_alpha_serach,scaling_factor_MT_temp,scaling_factor_CSEM_temp,rms_c,RMS_MT,RMS_CSEM,Lambda_c
	!			     scaling_factor_CSEM_temp=scaling_factor_CSEM_temp+0.1
	!			 end do  
 !53 continue      
 !     
 !      call scale_q_Matrix(k_step,scaling_factor_CSEM_opt,scaling_factor_MT_opt,b,m,q,beta_1)	
      
	
!
!   call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c,Au_vec,dHat)
!	!call scale_dHat(scaling_factor_CSEM_temp,scaling_factor_MT_temp,d,dHat)	
!	call un_normalize_with_dataVecMTX(dHat,d,1)
!					 	do i=1,d%nTx 
!				          do iDt=1, d%d(i)%nDt 
!					          dHat%d(i)%data(iDt)%errorBar= .false.
!				          end do
!				         end do			 
!	      call linComb (ONE,d_Pred_Orig,ONE,dHat,d_Pred_approx)
!		  call compute_RMS(m_c,d,d_Pred_approx,Lambda_c,RMS_c,res_c)
!		  mNorm_c=sqrt(dotProd(m_c,m_c))
!		  JmHat=dHat
	

	 
          !goto 999
 
                !call scale_q_Matrix(k_step,scaling_factor,b,m,q)
 

!call bidiag_1 (b,m,d,k_step,T_matrix,q,v,beta_1) 
!q=v


d_Pred	=d
call zero(d_Pred)
call zero(mhat)




!Start Line search using the projected system



if (Desired_rms .lt. 1.05) then
  Desired_rms=1.05
end if
  !Desired_rms=1.00
	   	   go_right=.true.
	       central=.true.
	       go_left=.true.
	       make_big_step=.false.
	       make_small_step=.false.
	       
	       step_size=(0.25)
10 continue
  Lambda_c=10**(lambda)
  Lambda_r=10**(Lambda+(step_size))
  Lambda_l=10**(Lambda-(step_size))


!Central Lambda
if (central) then
  call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c)
  call Calc_FWD(Lambda_c,d,m_c,d_Pred,res_c,eAll,F,mNorm_c,rms_c)
                  ! write(iterChar,'(i3.3)') fwd_calls
  				  ! file_name_suffix='Hybrid_'//iterChar
		          ! call write_output_files(DS_iter,d_Pred,m_c,file_name_suffix) 
  fwd_calls=fwd_calls +1
end if  
!
!right Lambda 
if (go_right) then
  call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_r,beta_1,m_r)
  call Calc_FWD(Lambda_r,d,m_r,d_Pred,res_r,eAll,F,mNorm_r,rms_r)
                  ! write(iterChar,'(i3.3)') fwd_calls
  				  ! file_name_suffix='Hybrid_'//iterChar
		          ! call write_output_files(DS_iter,d_Pred,m_r,file_name_suffix) 
  fwd_calls=fwd_calls +1
end if  

!left Lambda 
if (go_left) then
  call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_l,beta_1,m_l)
  call Calc_FWD(Lambda_l,d,m_l,d_Pred,res_l,eAll,F,mNorm_l,rms_l)
                  ! write(iterChar,'(i3.3)') fwd_calls
  				  ! file_name_suffix='Hybrid_'//iterChar
		          ! call write_output_files(DS_iter,d_Pred,m_l,file_name_suffix)     
  fwd_calls=fwd_calls +1
end if  

!rms_r=10000
!mNorm_r=10000
!rms_l=10000
!mNorm_l=10000
       
if (i_search .eq. 1 .and. i_sub_small_search .eq. 1 .and. i_sub_big_search .eq. 1  .and. i_sub_search .eq. 1 .and. k_step .eq. min_k_steps )then
   write(10,'(a32,i5,a13)')'########## Outer loop Iteration= ',DS_iter,'###########'
   write(10,'(a32,f10.4)') 'Start RMS for this Iteration  =', rms
   write(10,'(a32,f10.4)') 'Target RMS for this Iteration =', Desired_rms
   !write(10,*) 'Search minimum model norm     =', search_min_model
   write(10,'(a32,f10.4)') 'Step size                     =', step_size
   write(10,8510)
end if 
if (i_sub_big_search .eq. 2 .or. i_sub_small_search .eq. 2) then
   write(10,8511)
end if
   
 
if (i_sub_big_search .gt. 1 .or. i_sub_small_search .gt. 1)then
  write(10,8521) Lambda_r,rms_r,mNorm_r
  write(10,8531) Lambda_c,rms_c,mNorm_c
  write(10,8541) Lambda_l,rms_l,mNorm_l
else
  write(10,8520) Lambda_r,rms_r,mNorm_r
  write(10,8530) Lambda_c,rms_c,mNorm_c
  write(10,8540) Lambda_l,rms_l,mNorm_l  
end if  

if (rms_c .gt. rms .and. rms_r .gt. rms .and. rms_l .gt. rms .and. rms_l .gt. rms_c .and. rms_r .gt. rms_c) then
		 if(i_sub_big_search .eq. 1 .and. i_sub_small_search .eq. 1) then
			 write(10,*)'The RMS in all directions is higher than the start one'	 
		 end if	 
	 lambda=dlog10(lambda_c) 
     if (make_big_step) then
         if(i_sub_big_search .eq. 1) then
            write(10,*)'Try to make the step size bigger'
         end if   
    
            go_right=.true.
		    central=.false.
		    go_left=.true.
		    step_size=step_size*1.5
		    i_sub_big_search =i_sub_big_search+1		    
		    if (i_sub_big_search .lt. 10 ) then
		         write(10,*)'              ______________________________ '
		          goto 10
		    else
		         make_big_step=.false.
	             make_small_step=.false.
	             step_size=0.75  
	             i_sub_big_search=1 
		    end if
     elseif (make_small_step) then
		     if(i_sub_small_search .eq. 1) then
		            write(10,*)'Try to make the step size smaller'
		     end if       
    
            go_right=.true.
		    central=.false.
		    go_left=.true.
		    step_size=step_size*0.5			    
		    i_sub_small_search =i_sub_small_search+1 
		    if (i_sub_small_search .lt. 10 ) then
		          write(10,*)'              ______________________________ '
		          goto 10
		    else
		         make_big_step=.true.
	             make_small_step=.false.
	             step_size=0.75 
	             i_sub_small_search=1 
		    end if		             
     end if
else
    if(i_search .eq. 1 ) then
       if (rms_r .lt. rms_c ) then
	    write(10,*)' The RMS in RIGHT direction is smaller than the start one'
	    write(10,8560) step_size
	  elseif (rms_l .lt. rms_c) then
	    write(10,*)' The RMS in LEFT direction is smaller than the start one'
	    write(10,8560) step_size
	 else
	    write(10,*)' The RMS in CENTRAL is the smallest one'		
     end if
   end if    	     
8560    FORMAT(' Keep going in that direction with the same step size=',f10.3)	
     i_sub_small_search=1
     i_sub_big_search=1
end if
 
! If the central RMS is lower than the Desired one, start to look for the min. model norm.  
if (rms_c .lt. Desired_rms)then
        if ( i_sub_search.eq. 1 )then
                write(10,*)'Reached the target RMS: Search for the min. model norm'
        end if
         
		        ! Pick the smallest RMS and exit
		         if (rms_l .lt. rms_c ) then
		          		  m=m_l
						  lambda=dlog10(lambda_l)
						  rms=rms_l
						  res=res_l
						  mNorm=mNorm_l
				elseif (rms_r .lt. rms_c ) then
		          		  m=m_r
						  lambda=dlog10(lambda_r)
						  rms=rms_r
						  res=res_r
						  mNorm=mNorm_r
		        else
		           		  m=m_c
						  lambda=dlog10(lambda_c)
						  rms=rms_c
						  res=res_c
						  mNorm=mNorm_c
				end if
			      call linComb_modelParam(ONE,m,MinusONE,m0,mhat)
				  call Calc_FWD(lambda,d,m,d_Pred,res,eAll,F,mNorm,rms)
				  fwd_calls=fwd_calls +1
				   write(10,*) '1- Exit this Iteration with:'
				   write(10,8551)DS_iter, rms,mNorm,10**lambda,k_step,(k_step*2)+fwd_calls
				   write(10,*)
				    Desired_rms=rms/2
					start_kstep=2
					
				    search_min_model=.false.
				   file_name_suffix='Hybrid_end'
		           call write_output_files(DS_iter,d_Pred,m,file_name_suffix)  
				  !close(10)
				  goto 999
		  		         				 
		 		  
               
    ! pick the min. model norm within an RMS tolerance (1.1 of the Desired RMS) and exist search
        if (mNorm_r .lt. mNorm_c  )then
          ! Check if the RMS within the tolerance of the DESIRED RMS
           if (rms_r .lt. (Desired_rms*1.01)) then
                 write(10,*)'The RIGHT solution is less than the target RMS and has min. Model norm'
				  m=m_r
				  lambda=dlog10(lambda_r)
		   else
		  ! Search in direction of the mNorm_r for the best RMS
	   	         m=m_r
				 lambda=dlog10(lambda_c)        
		         step_size=step_size/2.0  
		      if ( i_sub_search.eq. 1 )then   
		         write(10,*)'Go RIGHT in direction of the min. model norm'
		      end if   
                 go_right=.true.
		         central=.false.
		         go_left=.false.
		         write(10,8510)
                   write(10,8520) Lambda_r,rms_r,mNorm_r
		         i_sub_search=i_sub_search+1
		         if (i_sub_search .lt. 10 ) then
		          goto 10
		         end if
		         
		          
		   end if
				  
		  
		elseif (mNorm_l .lt. mNorm_c )then
		  ! Check if the RMS within the tolerance of the DESIRED RMS
           if (rms_l .lt. (Desired_rms*1.01)) then
              write(10,*)'The LEFT solution is less than the target RMS and has min. Model norm '
			  m=m_l
			  lambda=dlog10(lambda_l)
		   else
		    ! Search in direction of the mNorm_l for the best RMS 
		      m=m_l
			  lambda=dlog10(lambda_c)
		      step_size=step_size/2
		      if ( i_sub_search.eq. 1 )then   
		      write(10,*)'Go LEFT in direction of the min. model norm'
		      end if
              go_right=.false.
		      central=.false.
		      go_left=.true.
               write(10,8540) Lambda_l,rms_l,mNorm_l
		      i_sub_search=i_sub_search+1
		      if (i_sub_search .lt. 5 ) then
		        goto 10
		      end if
	   	   
		   end if
		     	  
		else
		  write(10,*) 'The Central solution is within the target RMS tolerance and has min. Model norm '
		   m=m_c
		  lambda=dlog10(lambda_c)
		end if
		    		  
		  call linComb_modelParam(ONE,m,MinusONE,m0,mhat)
		  call Calc_FWD(lambda,d,m,d_Pred,res,eAll,F,mNorm,rms)
		  fwd_calls=fwd_calls +1
		   write(10,*) '1- Exit this Iteration with:'
		   write(10,8551)DS_iter, rms,mNorm,10**lambda,k_step,k_step*2+fwd_calls
		   write(10,*)
		   Desired_rms=rms/2
		    search_min_model=.false.
		   file_name_suffix='Hybrid_end'
           call write_output_files(DS_iter,d_Pred,m,file_name_suffix)  
		  !close(10)
		  goto 999
 end if
    
	  if (rms_r .lt. rms_c .and. abs(rms_c-rms_r) .gt. 0.001 )then
	   	  lambda=dlog10(lambda_r)
	   	  lambda_c=dlog10(lambda_r)
	   	   go_right=.true.
	       central=.false.
	       go_left=.false. 
	       rms_l=rms_c
	       mNorm_l=mNorm_c
	       m_l=m_c
	       rms_c=rms_r
	       mNorm_c=mNorm_r
	   	   m=m_r
	   	   m_c=m_r
	   	    write(10,*)'------------------------------------------------------------'
	  elseif (rms_l .lt. rms_c .and. abs(rms_c-rms_l) .gt. 0.001 )then

	      lambda=dlog10(lambda_l)
	      lambda_c=dlog10(lambda_l)
	       go_right=.false.
	       central=.false.
	       go_left=.true. 
	       rms_r=rms_c
	       mNorm_r=mNorm_c
	       m_r=m_c
	       rms_c=rms_l
	       mNorm_c=mNorm_l
	       m=m_l
	       m_c=m_l


	        write(10,*)'------------------------------------------------------------'
	  else
		  m=m_c
		  lambda=dlog10(lambda_c)
		  if (k_step .lt. max_k_steps) then
		    if (rms_c .gt. rms_temp )then
		           		  m=m_c
						  lambda=dlog10(lambda_c)
						  rms=rms_c
						  res=res_c
						  mNorm=mNorm_c
			      !call linComb_modelParam(ONE,m,MinusONE,m0,mhat)
				  call linComb_modelParam(ONE,m,MinusONE,m0,mhat)
				  !call Calc_FWD(lambda,d,m,d_Pred,res,eAll,F,mNorm,rms)
				  !fwd_calls=fwd_calls +1
				   write(10,*) '2- Exit this Iteration with:'
				   write(10,8551)DS_iter, rms,mNorm,10**lambda,k_step,k_step*2+fwd_calls
				  write(10,*)
				   file_name_suffix='Hybrid_end'
		           call write_output_files(DS_iter,d_Pred,m,file_name_suffix)  
					  close(10)
					goto 999
		    else			        
				      write(10,*) 'Increase the size of the projected system by adding additional K step'
				     !call Calc_FWD(lambda,d,m_start,d_Pred,res,eAll,F,mNorm,rms)
				     eAll=eAll_temp
				     !fwd_calls=fwd_calls +1
		             k_step=k_step+1
		             update_proj_sys=.true.
		             close(10)
		             !rms_temp=rms_c
		             m_temp=m_c
		             goto 1
            end if
             
           else    
		   		          m=m_c
						  lambda=dlog10(lambda_c)
						  rms=rms_c
						  res=res_c
						  mNorm=mNorm_c
			  call linComb_modelParam(ONE,m,MinusONE,m0,mhat)
			  !call Calc_FWD(lambda,d,m,d_Pred,res,eAll,F,mNorm,rms)
			  !fwd_calls=fwd_calls +1
			   write(10,*) '2- Exit this Iteration with:'
			   write(10,8551)DS_iter, rms,mNorm,10**lambda,k_step,k_step*2+fwd_calls
			  write(10,*)
			   file_name_suffix='Hybrid_end'
	           call write_output_files(DS_iter,d_Pred,m,file_name_suffix)  
				  close(10)
				goto 999
		   end if	  
	  end if
  
	  

  i_search=i_search+1
  
  if (i_search .gt. 300 ) then
  		  m=m_c
		  lambda=dlog10(lambda_c)
  		  call linComb_modelParam(ONE,m,MinusONE,m0,mhat)
		  call Calc_FWD(lambda,d,m,d_Pred,res,eAll,F,mNorm,rms)
		  fwd_calls=fwd_calls +1
		   write(10,*) 'Reached the max. number of line search iterations'	  
		   write(10,*) 'Exit this Iteration with:'
		   write(10,8551)DS_iter, rms,mNorm,10**lambda,k_step,k_step*2+fwd_calls
		   write(10,*)		  
		   file_name_suffix='Hybrid_end'
           call write_output_files(DS_iter,d_Pred,m,file_name_suffix)  
	  close(10)
	  goto 999
  else
      goto 10
  end if	  
	        
           
 999 continue   
 
              write(10,*) 'End RMSs'
              call compute_RMS(m,d,d_Pred,Lambda,RMS_c,res_c)
			  write(10,'(4f10.3)') rms_c,RMS_MT,RMS_CSEM,Lambda_c
 

    close(10)
   call deall_solnVectorMTX(eAll_temp)
   
8510  FORMAT('         :    Lambda       RMS     mNorm ')
8520  FORMAT(' Right   : ',3f10.3)
8530  FORMAT(' Central : ',3f10.3)
8540  FORMAT(' Left    : ',3f10.3)
8511  FORMAT('              :    Lambda       RMS     mNorm ')  
8521  FORMAT('      Right   : ',3f10.3)
8531  FORMAT('      Central : ',3f10.3)
8541  FORMAT('      Left    : ',3f10.3)  

8550  FORMAT(' RMS = ',f10.3, ' mNorm = ',f10.3,' Lambda = ', f10.3 )
8551  FORMAT(' Iter No: ',i5,' RMS = ',f10.3, ' mNorm = ',f10.3,' Lambda = ', f10.3, ' K_Steps = ', i5, '# FWD_Calls', i5 )
		 	    
end subroutine Lanczos_DS 
!**************************************************************************************** 
Subroutine bidiag_1 (b,m,d,k_step,T_matrix,q,v,beta_1)
  !In
  type (dataVectorMTX_t), intent(in)	 	::b
  type(modelParam_t),  intent(in)       ::m
  type (dataVectorMTX_t), intent(in)	    ::d
  Integer            , intent(in)       :: k_step
  !out
  real(kind=prec),intent(inout),dimension(:,:)  :: T_matrix
  type (modelParam_t),intent(inout), dimension(:) :: q,v

  real(kind=prec),intent(out)                     :: beta_1
  
 !Local
   type (dataVectorMTX_t),pointer, dimension(:) :: u
  real(kind=prec)	,pointer, dimension(:)  :: beta
  type (dataVectorMTX_t)                   	:: Av
  real(kind=prec)                           :: alpha,sigma11
  Integer                                   :: i,iDt,ii,jj
    real(kind=prec),pointer,dimension(:,:)  :: T_matrix_T,T_matrix_Temp

  

allocate(beta(k_step+1)) 
allocate (T_matrix_T(k_step+1,k_step+1),T_matrix_Temp(k_step+1,k_step+1))
allocate(u(0:k_step+1))
T_matrix_T=R_zero
T_matrix_Temp=R_zero

Av=b

   call zero(Av)

   
T_matrix=R_zero

do ii=1,k_step+1
   u(ii)=b
   call zero(u(ii))
end do

  do ii=1,k_step+1
   q(ii)=m
   v(ii)=m
   call zero(q(ii))
   call zero(v(ii))
end do

u(1)=b
beta(1)=sqrt(dotProd(u(1),u(1)))
call scMult(1/beta(1),u(1),u(1))
beta_1=beta(1)

! Compute v =Cm JT Cd^(-1/2) u 

           call normalize_with_dataVecMTX(u(1),d,1)          
#ifdef MPI
            call Master_job_JmultT(m,u(1),q(1),eAll)
#else
            call JmultT(m,u(1),q(1),eAll)
#endif
            
            q(1)= multBy_Cm(q(1))
            v(1)=q(1)

			alpha=sqrt(dotProd(v(1),v(1)))
			call scMult(1/alpha,v(1),v(1))

do ii=1,k_step
 

#ifdef MPI
            call Master_job_Jmult(v(ii),m,Av,eAll)
#else
            call Jmult(v(ii),m,Av,eAll)
#endif
            call normalize_with_dataVecMTX(Av,d,1) 
 
		     do i=1,b%nTx 
	          do iDt=1, b%d(i)%nDt 
		          Av%d(i)%data(iDt)%errorBar= .false.
	          end do
	         end do  
	            
	 		 call linComb (ONE,Av,-alpha,u(ii),u(ii+1))
	 		 ! Gram-Schmidt orthogonalization

	 		    do jj=1,ii-1
	 		         sigma11= dotProd(u(ii+1),u(jj))/ dotProd(u(jj),u(jj))
	 		         
				 		do i=1,b%nTx 
				          do iDt=1, b%d(i)%nDt 
					          u(ii+1)%d(i)%data(iDt)%errorBar= .false.
					          u(jj)%d(i)%data(iDt)%errorBar= .false.
				          end do
				         end do
				         
	 		         call linComb (ONE,u(ii+1),-sigma11,u(jj),u(ii+1))
	 		    end do
 
	 		 
 		            
             beta(ii+1)=sqrt(dotProd(u(ii+1),u(ii+1)))
             call scMult(1/beta(ii+1),u(ii+1),u(ii+1))
         

            call normalize_with_dataVecMTX(u(ii+1),d,1)          
#ifdef MPI
            call Master_job_JmultT(m,u(ii+1),q(ii+1),eAll)
#else
            call JmultT(m,u(ii+1),q(ii+1),eAll)
#endif

            q(ii+1)= multBy_Cm(q(ii+1))
            call linComb (ONE,q(ii+1),-beta(ii+1),v(ii),v(ii+1))
            !v(ii+1)= multBy_Cm(v(ii+1))
             
	 		 ! Gram-Schmidt orthogonalization

	 		    do jj=1,ii-1
	 		         sigma11= dotProd(v(ii+1),v(jj))/ dotProd(v(jj),v(jj)) 
	 		         call linComb (ONE,v(ii+1),-sigma11,v(jj),v(ii+1))
	 		    end do
	 		    
        if (ii .eq. 1 )then
          T_matrix(ii,ii)=alpha
        elseif (ii .le. k_step )then
          T_matrix(ii,ii)=alpha
          T_matrix(ii,ii-1)=beta(ii+1)
       end if
       
       
            alpha=sqrt(dotProd(v(ii+1),v(ii+1)))
            call scMult(1/alpha,v(ii+1),v(ii+1)) 
                       

      write(6,*) ii, 'Beta= ', beta(ii+1)/beta(1),beta(ii)
      write(6,*) ii, 'Alpha= ', alpha
      write(6,*)ii, 'Orth 1,ii=', dotProd(u(1),u(ii))
       
 end do
 
             T_matrix_Temp=T_matrix
             T_matrix_T= transpose(T_matrix)
             T_matrix=matmul(T_matrix_T,T_matrix_Temp)
 
			  write(6,*)'############### bidiag1 T matrix #############'
			  do ii=1,k_step
			    write(6,'(3f10.2)')(T_matrix_temp(ii,jj),jj=1,k_step)
			 end do
			 write(6,*)'############### bidiag1 T matrix #############'
			  do ii=1,k_step
			    write(6,'(3f10.2)')(T_matrix(ii,jj),jj=1,k_step)
			 end do
			  
 beta_1=beta(1)
 
end Subroutine bidiag_1
!**************************************************************************************** 
Subroutine bidiag_JJT_test (b,m,d,k_step,DS_iter,lambda,T_matrix,q,beta_1,beta_kstep_plus_1,u)
  !In
  type (dataVectorMTX_t), intent(in)	 	::b
  type(modelParam_t),  intent(in)       ::m
  type (dataVectorMTX_t), intent(in)	    ::d
  Integer            , intent(inout)       :: k_step
  Integer            , intent(in)       ::DS_iter
  real(kind=prec),     intent(in)		   :: 	lambda	
  !out
  real(kind=prec),intent(inout),dimension(:,:)  :: T_matrix
  type (modelParam_t),intent(inout), dimension(:) :: q
  type (dataVectorMTX_t), intent(inout), dimension(:)  ::u
  real(kind=prec),intent(out)                     :: beta_1,beta_kstep_plus_1
  
 !Local
  real(kind=prec)	,pointer, dimension(:)  :: beta
  type (dataVectorMTX_t)                   	:: Au,u_temp,w
  real(kind=prec)                           :: alpha,sigma11,alpha_beta,lambda_sub,sum,sum1,norm1
  Integer                                   :: i,iDt,ii,jj,INFO,kk,nTx
  real(kind=prec)	,pointer, dimension(:)    ::  AP,b_sub,b_Tx,x_sub,alpha_vec,sigma11_vec
  real(kind=prec),pointer,dimension(:,:) ::T_matrix_temp

   type(modelParam_t)                         :: q_temp,s_hat_temp
  
  !Temp vectros
     type (dataVectorMTX_t)                   	:: r_vec,v_vec,Aq_vec,Qr
     type (dataVectorMTX_t), pointer, dimension(:)  :: q_vec
  
     nTx=b%nTx
     allocate(s_hat(nTx))

	  
	allocate(beta(k_step+1),alpha_vec(k_step+1),q_vec(k_step+1),sigma11_vec(k_step+1))   

do ii=1,k_step+1
  q_vec(ii)= b
  call zero(q_vec(ii))
end do

r_vec= b
v_vec=b
Aq_vec=b
Qr=b
   call zero(r_vec)
   call zero(v_vec)
   call zero(Aq_vec)   

         q_vec(1)=b   
  		 norm1=sqrt(dotProd(b,b))
         call scMult(1.0/norm1,b,q_vec(1)) 
         call MultA_JJT_DS(q_vec(1),m,d,r_vec,q(ii))
		 alpha_vec(1)=dotProd(q_vec(1),r_vec)
		 
		 do i=1,b%nTx 
          do iDt=1, b%d(i)%nDt 
	          r_vec%d(i)%data(iDt)%errorBar= .false.
          end do
         end do 
		 call linComb (ONE,r_vec,-alpha_vec(1),q_vec(1),r_vec)
		 
		  beta(1)=sqrt(dotProd(r_vec,r_vec))
		  write(75,*) 1,beta(1)
		  
		  do ii=2,k_step
		     v_vec=q_vec(ii-1)
			 call scMult(1.0/beta(ii-1),r_vec,q_vec(ii)) 
			 
			 call MultA_JJT_DS(q_vec(ii),m,d,Aq_vec,q(ii))
			 
			 do i=1,b%nTx 
			  do iDt=1, b%d(i)%nDt 
				  Aq_vec%d(i)%data(iDt)%errorBar= .false.
			  end do
			 end do 	
			 call linComb (ONE,Aq_vec,-beta(ii-1),v_vec,r_vec)
			 alpha_vec(ii)=dotProd(q_vec(ii),r_vec)
			 
			 do i=1,b%nTx 
			  do iDt=1, b%d(i)%nDt 
				  r_vec%d(i)%data(iDt)%errorBar= .false.
			  end do
			 end do 			 
			  call linComb (ONE,r_vec,-alpha_vec(ii),q_vec(ii),r_vec)


				
              
				!	  do jj=1,ii
				!		   sigma11_vec(jj)= dotProd(q_vec(jj),r_vec)
				!	  end do
				!		!
				!	  call zero(Qr)  	
				!	  do jj=1,ii    
				!				do i=1,d%nTx 
				!				  do iDt=1, d%d(i)%nDt 
				!					  q_vec(jj)%d(i)%data(iDt)%errorBar= .false.
				!				  end do
				!				 end do
				!				call scMult(sigma11_vec(jj),q_vec(jj),q_vec(jj)) 
				!				call linComb (ONE,Qr,ONE,q_vec(jj),Qr)
				!	  end do
				!	  
				!				do i=1,d%nTx 
				!				  do iDt=1, d%d(i)%nDt 
				!					  Qr%d(i)%data(iDt)%errorBar= .false.
				!				  end do
				!				 end do
				!				 call linComb (ONE,r_vec,-ONE,Qr,r_vec)
								  
			  beta(ii)=sqrt(dotProd(r_vec,r_vec))
			  write(75,*) ii,beta(ii)
		  end do	  
		 
 


         goto 10

Au=b
w=b
u_temp=b
q_temp=m
s_hat_temp=m

   call zero(Au)
   call zero(w)
   call zero(u_temp)
   call zero(q_temp)   
   call zero(s_hat_temp)   
do ii=1,21
  u(ii)= b
  call zero(u(ii))
end do

beta=R_zero
   
   u(1)=b
 !  Mult Cd^(-1/2) u(1) 
         !call normalize_with_dataVecMTX(u(1),d,1)  
		 beta(1)=sqrt(dotProd(u(1),u(1)))
         call scMult(1.0/beta(1),u(1),u(1))
         beta_1=beta(1)
 

 

 
 do ii=1,k_step
 
	  do jj=1,nTx
	  	s_hat(jj)=m
	    call zero(s_hat(jj))
	  end do
	  
		 call MultA_JJT_DS(u(ii),m,d,Au,q(ii))
		 
	  do jj=1,nTx
	     s_hat_temp=multBy_Cm(s_hat(jj))
	  	 JTw_matrix(ii,jj)=s_hat_temp
	  end do		 
	  
	  		 alpha=sqrt(dotProd(q(ii),q(ii)))
            !call scMult(1.0/alpha,q(ii),q(ii))
	  
	     do i=1,b%nTx 
          do iDt=1, b%d(i)%nDt 
	          Au%d(i)%data(iDt)%errorBar= .false.
          end do
         end do 
		 
		 if (ii .gt. 1) then
		     do i=1,b%nTx 
	          do iDt=1, b%d(i)%nDt 
		          u(ii-1)%d(i)%data(iDt)%errorBar= .false.
	          end do
	         end do 		 
		   call linComb (ONE,Au,-beta(ii),u(ii-1),w)
		 else
		    w=Au
		 end if
		  
		 
		 alpha=(dotProd(u(ii),w))
		
		 do i=1,b%nTx 
          do iDt=1, b%d(i)%nDt 
	          w%d(i)%data(iDt)%errorBar= .false.
          end do
         end do 
         
         
		call linComb (ONE,w,-alpha,u(ii),u(ii+1))
		
	 		    do jj=1,ii
	 		         sigma11= dotProd(u(ii+1),u(jj))
	 		         
				 		do i=1,d%nTx 
				          do iDt=1, d%d(i)%nDt 
					          !u(ii+1)%d(i)%data(iDt)%errorBar= .false.
					          u(jj)%d(i)%data(iDt)%errorBar= .false.
				          end do
				         end do
				         
			         call linComb (ONE,u(ii+1),-sigma11,u(jj),u(ii+1))
	 		    end do
	 	!call zero(u_temp)	    		
        beta(ii+1)=sqrt(dotProd(u(ii+1),u(ii+1)))
        call scMult(1.0/beta(ii+1),u(ii+1),u(ii+1))

        
        if (ii .eq. 1 )then
          T_matrix(ii,ii)=alpha
          T_matrix(ii,ii+1)=beta(ii+1)
        elseif (ii .lt. k_step )then
          T_matrix(ii,ii)=alpha
          T_matrix(ii,ii-1)=beta(ii)
          T_matrix(ii,ii+1)=beta(ii+1)
       else
          T_matrix(ii,ii)=alpha
          T_matrix(ii,ii-1)=beta(ii)
       end if
      
	    write(75,*)beta(ii),beta(ii+1),alpha,beta(ii+1)*alpha

									!if (ii .ge. 2 ) then
                                  !
									!		  if (associated(AP)) then
									!			deallocate(AP,b_sub,b_Tx)
									!		  end if
									!		  allocate(AP((ii)*((ii)+1)/2),b_sub(ii),b_Tx(ii))
									!		   AP=R_zero
									!		   do jj=1,ii
									!			  do kk=1,jj
									!			   AP(kk + (jj-1)*jj/2) = T_matrix(kk,jj) 
									!			  end do
									!			end do  
									!			 
									!			!Preforme Cholesky decomposition on (sTs +Lambda I) matrix 
									!			call DPPTRF( 'U', ii, AP, INFO )
									!			!Solve the problem, save solution in b_sub
									!			 b_sub=R_zero
									!			 b_sub(1)=beta_1
									!			 call DPPTRS( 'U', ii, 1, AP, b_sub, ii, INFO )
									!			 write(75,*) ii,b_sub(ii),b_sub(ii)*beta(ii+1)
									!			 
									!		   sum1=R_zero
									!		   do jj=1,ii
									!		   sum=R_zero
									!			  do kk=1,ii
									!			   sum= sum+(T_matrix(jj,kk)* b_sub(kk))
									!			  end do
									!				if (jj==1 ) then
									!				  b_Tx(jj)=beta_1-sum
									!				else
  									!				  b_Tx(jj)=R_zero-sum
									!				end if  
									!				sum1=sum1+(b_Tx(jj)*b_Tx(jj))
									!			end do 		
									!			 write(85,*) ii,sqrt(sum1)
									!			 
									!end if
		                          !
        
        
        
        
      write(6,*) ii, 'Alpha= ',alpha,'qTq= ',dotProd(q(ii),q(ii))
      write(6,*) ii, 'Orth= ', dotProd(u(1),u(ii))
 end do
 
beta_kstep_plus_1=beta(k_step+1)

  10 continue 
 do ii=1,k_step
 write(20,'(a6,5f10.5)') 'Orth= ', (dotProd(u(ii),u(jj)),jj=1,k_step)
 end do

 
 beta_1=beta(1)
 
			  
 
 
 
 
 

	  

end Subroutine bidiag_JJT_test
!**************************************************************************************** 
Subroutine bidiag_JJT (b,m,d,k_step,DS_iter,lambda,T_matrix,q,beta_1,beta_kstep_plus_1,u,Au_vec)
  !In
  type (dataVectorMTX_t), intent(in)	 	::b
  type(modelParam_t),  intent(in)       ::m
  type (dataVectorMTX_t), intent(in)	    ::d
  Integer            , intent(inout)       :: k_step
  Integer            , intent(in)       ::DS_iter
  real(kind=prec),     intent(in)		   :: 	lambda	
  !out
  real(kind=prec),intent(inout),dimension(:,:)  :: T_matrix
  type (modelParam_t),intent(inout), dimension(:) :: q
  type (dataVectorMTX_t), intent(inout), dimension(:)  ::u
  real(kind=prec),intent(out)                     :: beta_1,beta_kstep_plus_1
  type (dataVectorMTX_t) ,intent(inout), dimension(:) ::                 	 Au_vec
  
 !Local
  real(kind=prec)	,pointer, dimension(:)  :: beta,w1
  type (dataVectorMTX_t)                   	:: Au,u_temp,w
  real(kind=prec)                           :: alpha,sigma11,sigma22,sigma33,alpha_beta,lambda_sub,sum,sum1
  Integer                                   :: i,iDt,ii,jj,INFO,kk,nTx
  real(kind=prec)	,pointer, dimension(:)    ::  AP,b_sub,b_Tx,x_sub,alpha_vec
  real(kind=prec),pointer,dimension(:,:) ::T_matrix_temp,v1

   type(modelParam_t)                         :: q_temp,s_hat_temp
  
  
  
     nTx=b%nTx
     allocate(s_hat(nTx),alpha_vec(k_step))

	  
	  
  
allocate(T_matrix_temp(k_step,k_step))  
allocate(AP((k_step)*((k_step)+1)/2),b_sub(k_step),b_Tx(k_step))

 
allocate(beta(k_step+1)) 



Au=b
w=b
u_temp=b
q_temp=m
s_hat_temp=m

   call zero(Au)
   call zero(w)
   call zero(u_temp)
   call zero(q_temp)   
   call zero(s_hat_temp)   


beta=R_zero
   
   u(1)=b
 !  Mult Cd^(-1/2) u(1) 
         !call normalize_with_dataVecMTX(u(1),d,1)  
		 beta(1)=sqrt(dotProd(u(1),u(1)))
         call scMult(1.0/beta(1),u(1),u(1))
         beta_1=beta(1)
 

 

 
 do ii=1,k_step
 
	  do jj=1,nTx
	  	s_hat(jj)=m
	    call zero(s_hat(jj))
	  end do
	  
		 call MultA_JJT_DS(u(ii),m,d,Au,q(ii))
		 
		 Au_vec(ii)=Au
		 
	  do jj=1,nTx
	     s_hat_temp=multBy_Cm(s_hat(jj))
	  	 JTw_matrix(ii,jj)=s_hat_temp
	  end do		 
	  
	  		 !alpha=sqrt(dotProd(q(ii),q(ii)))
            !call scMult(1.0/alpha,q(ii),q(ii))
	  
	     do i=1,b%nTx 
          do iDt=1, b%d(i)%nDt 
	          Au%d(i)%data(iDt)%errorBar= .false.
          end do
         end do 
		 
		 if (ii .gt. 1) then
		     do i=1,b%nTx 
	          do iDt=1, b%d(i)%nDt 
		          u(ii-1)%d(i)%data(iDt)%errorBar= .false.
	          end do
	         end do 		 
		   call linComb (ONE,Au,-beta(ii),u(ii-1),w)
		 else
		    w=Au
		 end if
		  
		 
		 alpha=(dotProd(u(ii),w))
		
		 do i=1,b%nTx 
          do iDt=1, b%d(i)%nDt 
	          w%d(i)%data(iDt)%errorBar= .false.
          end do
         end do 
         
         
		call linComb (ONE,w,-alpha,u(ii),u(ii+1))
		  
		    
		        call zero(u_temp)
	 		   do jj=1,ii-1
				  sigma11= dotProd(u(jj),u(ii+1))
				  sigma22= dotProd(u(jj),u(jj))
				  sigma33=sigma11/sigma22
				  
				 		do i=1,d%nTx 
				          do iDt=1, d%d(i)%nDt 
					          u(jj)%d(i)%data(iDt)%errorBar= .false.
				          end do
				         end do
				    call linComb (ONE,u_temp,sigma33,u(jj),u_temp)
                end do
				 
				 		do i=1,d%nTx 
				          do iDt=1, d%d(i)%nDt 
					          u(ii+1)%d(i)%data(iDt)%errorBar= .false.
				          end do
				         end do
			         call linComb (ONE,u(ii+1),-ONE,u_temp,u(ii+1))
			    
    		
        beta(ii+1)=sqrt(dotProd(u(ii+1),u(ii+1)))
        call scMult(1.0/beta(ii+1),u(ii+1),u(ii+1))

        
        if (ii .eq. 1 )then
          T_matrix(ii,ii)=alpha
          T_matrix(ii,ii+1)=beta(ii+1)
        elseif (ii .lt. k_step )then
          T_matrix(ii,ii)=alpha
          T_matrix(ii,ii-1)=beta(ii)
          T_matrix(ii,ii+1)=beta(ii+1)
       else
          T_matrix(ii,ii)=alpha
          T_matrix(ii,ii-1)=beta(ii)
       end if
      
	    write(75,'(3f10.3)')beta(ii+1)/beta(1),beta(ii+1)/alpha,alpha/beta(ii+1)

!SVD on T        
     if (ii .gt. 1) then

         allocate(T_matrix_temp(ii,ii),w1(ii),v1(ii,ii))
         do jj=1,ii
            do kk=1, ii 
             T_matrix_temp(jj,kk)=T_matrix(jj,kk)
            end do 
         end do   
          call SVDCMP(T_matrix_temp,ii,ii,ii,ii,w1,v1)
          call SVDSORT(T_matrix_temp,w1,v1,ii,ii)
          do jj=1,ii
          write(100,*)jj,w1(jj)
          end do
         do jj=1,ii
          write(110,*)(v1(jj,kk),kk=1,ii)
         end do
         write(120,'(i5,5f12.5)')ii,beta(ii+1),v1(ii,ii),beta(ii+1)*abs(v1(ii,ii)),w1(ii),10**log10(w1(ii-1))
          
          !lambda=log10(w1(ii))
          !if (w1(ii) .lt. w1(1)/2.0)   goto 10
          !if (beta(ii+1)*abs(v1(ii,ii)) .lt. 1.0) goto 10
         
         deallocate(w1,v1,T_matrix_temp)
    end if    
! End SVD
 end do
 
beta_kstep_plus_1=beta(k_step+1)

  10 continue 
 do ii=1,k_step
  do jj=1,k_step
     write(20,*) 'Orth= ',ii,jj, (dotProd(u(ii),u(jj)))
  end do
 end do

 
 beta_1=beta(1)
 
			  
 
 
 
 
 

	  

end Subroutine bidiag_JJT
!**************************************************************************************** 
Subroutine Arnold_bidiag_JJT (b,m,m0,d,k_step,DS_iter,lambda,T_matrix,q,beta_1,beta_kstep_plus_1,u)
  !In
  type (dataVectorMTX_t), intent(in)	 	::b
  type(modelParam_t),  intent(in)       ::m,m0
  type (dataVectorMTX_t), intent(in)	    ::d
  Integer            , intent(in)       :: k_step,DS_iter
  real(kind=prec),intent(in)		    :: lambda
  !out
  real(kind=prec),intent(inout),dimension(:,:)  :: T_matrix
  type (modelParam_t),intent(inout), dimension(:) :: q
  type (dataVectorMTX_t),intent(inout), dimension(:) :: u

  real(kind=prec),intent(out)                     :: beta_1,beta_kstep_plus_1
  
 !Local
  real(kind=prec)	,pointer, dimension(:)  :: beta
  type (dataVectorMTX_t)                   	:: Au,u_temp,u_temp1,Jm,r
  real(kind=prec)                           :: lambda_sub,alpha,norm_Au_1,norm_Au_end,alpha_beta,lambda_test,sum
  Integer                                   :: i,iDt,ii,jj,INFO,restart,kk
  type (modelParam_t)                       :: mHat,q_temp
  real(kind=prec)	,pointer, dimension(:)    ::  AP,b_sub,x_sub,sub_residual
  real(kind=prec),pointer,dimension(:,:) ::T_matrix_temp
  
allocate(T_matrix_temp(k_step+1,k_step+1))  
allocate(AP((k_step)*((k_step)+1)/2),b_sub(k_step),x_sub(k_step),sub_residual(k_step))
allocate(beta(k_step+1)) 



if  (DS_iter .eq. 1 ) then
  lambda_sub=10.0
else
  if (lambda .lt. 0.1 )then
    lambda_sub=0.1
  else
    lambda_sub=lambda
  end if  
end if

mHat=m
q_temp=m
call zero(mHat)
call zero(q_temp)

Au=b
u_temp=b
u_temp1=b
Jm=b
r=b
   call zero(Au)
   call zero(u_temp1)
   call zero(Jm)
   
T_matrix=R_zero

do ii=0,k_step+1
   u(ii)=b
   call zero(u(ii))
end do

  do ii=1,k_step
   q(ii)=m
   call zero(q(ii))
end do


u(1)=b
beta(1)=sqrt(dotProd(u(1),u(1)))
call scMult(1/beta(1),u(1),u(1))
beta_1=beta(1)

 do ii=1,k_step
         
		 call MultA_JJT_DS(u(ii),m,d,Au,q(ii))
		 norm_Au_1=sqrt(dotProd(Au,Au))
		 
		 call zero(u_temp)
		 do jj=1,ii
			  T_matrix(jj,ii)=dotProd(u(jj),Au)
			  do i=1,b%nTx 
	           do iDt=1, b%d(i)%nDt 
		          Au%d(i)%data(iDt)%errorBar= .false.
		          u(jj)%d(i)%data(iDt)%errorBar= .false.
	           end do
	         end do 
         
		     call linComb (ONE,Au,-T_matrix(jj,ii),u(jj),Au)
		 end do

		 T_matrix(ii+1,ii)=sqrt(dotProd(Au,Au))
		 
         call scMult(1/T_matrix(ii+1,ii),Au,u(ii+1))
         norm_Au_end=sqrt(dotProd(Au,Au))

         alpha_beta=T_matrix(ii+1,ii)
         call scMult(alpha_beta,u(ii+1),u_temp)
         
         
      write(6,*) ii, 'Beta= ', alpha_beta/beta_1,sqrt(dotProd(u_temp,u_temp))
      
      if (ii .gt. 1 )then
             
             T_matrix_temp=R_Zero
             T_matrix_temp=T_matrix
			 do jj=1,ii
			    T_matrix_temp(jj,jj)=T_matrix(jj,jj)+  lambda_sub
			 end do
			 
		   AP=R_zero
		   do jj=1,ii
			  do kk=1,jj
               AP(kk + (jj-1)*jj/2) = T_matrix_temp(kk,jj) 
			  end do
			end do  
			 
!7-  Preforme Cholesky decomposition on (sTs +Lambda I) matrix 
            call DPPTRF( 'U', ii, AP, INFO )
!9- solve the problem, save solution in b_sub
             x_sub=R_zero
             x_sub(1)=beta_1
             call DPPTRS( 'U', ii, 1, AP, x_sub, ii, INFO )
             


             write(6,*) 'Sub_norm',T_matrix(ii+1,ii)*abs(x_sub(ii)), T_matrix(ii+1,ii)*abs(x_sub(ii))/beta_1

              
      end if  
      

 end do
 beta_kstep_plus_1=T_matrix(k_step+1,k_step)
 
			  do ii=1,k_step
			    write(6,'(a6,5f10.2)')'Matrix',(T_matrix(ii,jj),jj=1,k_step)
			 end do

			 do ii=1,k_step
			 write(6,'(a6,5f10.5)') 'Orth= ', (dotProd(u(ii),u(jj)),jj=1,k_step)
			 end do          
 
end Subroutine Arnold_bidiag_JJT
!****************************************************************************************  
subroutine get_model (k_step,T_matrix,q,m0,mhat,b,lambda,beta_1,m,Au_vec,dHat)

Integer, intent(in)                          :: k_step
real(kind=prec)	,intent(in), dimension(:,:)  :: T_matrix
type(modelParam_t),  intent(in)              :: m0
type(modelParam_t),  intent(inout)           :: mhat
type (modelParam_t),intent(in), dimension(:) :: q
!type (dataVectorMTX_t),intent(in), dimension(:) :: u
type (dataVectorMTX_t),intent(in)               :: b
real(kind=prec)	,intent(in)                  :: lambda,beta_1
type (modelParam_t),intent(out)              :: m
type (dataVectorMTX_t),intent(inout),pointer, dimension(:),optional               :: Au_vec
type (dataVectorMTX_t),intent(inout),optional               :: dHat
!Local 

real(kind=prec)	,pointer, dimension(:,:)  ::  T_matrix_temp,L_matrix, U_matrix
real(kind=prec)	,pointer, dimension(:)    ::  b_sub,y_sub,x_sub
type (modelParam_t)                       :: q_temp
type (dataVectorMTX_t)                       ::dHat_temp
Integer                                   :: ii,jj,kk,INFO,i,Idt
real(kind=prec)                           :: lambda_LU,mu_LU,sum,sum1,b_Tx
  character(100)       		:: file_name_suffix
    real(kind=prec)	,pointer, dimension(:)    ::  AP
  
  allocate(T_matrix_temp(k_step,k_step),L_matrix(k_step+1,k_step+1),U_matrix(k_step+1,k_step+1),b_sub(k_step),y_sub(k_step),x_sub(k_step))
  allocate(AP((k_step)*((k_step)+1)/2))
T_matrix_temp=R_zero
L_matrix=R_zero
U_matrix=R_zero
y_sub=R_zero
x_sub=R_zero
b_sub=R_zero


q_temp=m0
dHat_temp=b

             T_matrix_temp=R_Zero
             do ii=1,k_step
               do jj=1,k_step
                 T_matrix_temp(ii,jj)=T_matrix(ii,jj)
               end do
             end do    
			 do ii=1,k_step
			    T_matrix_temp(ii,ii)=T_matrix(ii,ii)+lambda
			 end do
			 
		   do jj=1,k_step
			  do ii=1,jj
               AP(ii + (jj-1)*jj/2) = T_matrix_temp(ii,jj) 
			  end do
			end do  
			 
!7-  Preforme Cholesky decomposition on (sTs +Lambda I) matrix 
            call DPPTRF( 'U', k_step, AP, INFO )
!9- solve the problem, save solution in b_sub
             b_sub(1)=beta_1
			 

              write(150,*) lambda
               write(150,*)'### L Matrix ###'
              do ii=1,k_step
                  write(150,*)(T_matrix_temp(ii,jj),jj=1,k_step) 
              end do  
              write(150,*)'### RHS ###'
               do ii=1,k_step   
                  write(150,*)b_sub(ii)
               end do  
               
             call DPPTRS( 'U', k_step, 1, AP, b_sub, k_step, INFO )
              x_sub=b_sub
              

               write(150,*)'### Solution ###'
               do ii=1,k_step  
                  write(150,*)x_sub(ii)
              end do
               write(150,*)'###########'
              
 
 											   sum1=R_zero
											   do jj=1,k_step
											   sum=R_zero
													do kk=1,k_step
													   sum= sum+(T_matrix_temp(jj,kk)* x_sub(kk))
													end do
													if (jj==1 ) then
													  b_Tx=beta_1-sum
													else
  													  b_Tx=R_zero-sum
													end if  
													sum1=sum1+(b_Tx*b_Tx)
												end do 		
												 !write(25,*) lambda,sqrt(sum1)/beta_1	
												! write(10,*) lambda,sqrt(sum1)	
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  
			  goto 10
            
            
            			   
			 mu_LU=T_matrix_temp(1,1)
			 ! Preforme LU decomposition on (T_matrix + lambda I)
			  do ii=1,k_step
			      do jj=1,k_step
					   if (ii .eq. jj ) then
					     L_matrix(ii,jj)=ONE
					     U_matrix(ii,jj)= mu_LU
					   end if
			
					   if (ii .gt. 1 .and. jj .eq. ii-1) then
					      lambda_LU=T_matrix_temp(ii,jj)/mu_LU
					      L_matrix(ii,jj)=lambda_LU
					      mu_LU=T_matrix_temp(ii,ii)-lambda_LU*T_matrix_temp(ii,jj)
					   end if
					   if (ii .lt. k_step .and. jj .eq. ii+1) then
					      U_matrix(ii,jj)=T_matrix_temp(ii,jj)
					   end if
			    end do
			  end do
			  
			  write(6,*)'############### T matrix #############'
			  do ii=1,k_step
			    write(6,'(3f10.2)')(T_matrix_temp(ii,jj),jj=1,k_step)
			 end do

			  write(6,*)'############### L matrix #############'
			  do ii=1,k_step
			    write(6,'(3f10.2)')(L_matrix(ii,jj),jj=1,k_step)
			 end do
			   write(6,*)'############### U matrix #############'
			   do ii=1,k_step
			    write(6,'(3f10.2)')(U_matrix(ii,jj),jj=1,k_step)
			 end do
! Make the right hand side for the subproblem U^T (d- F(m)+J (m-m0))= U^T b
write(6,*)'############### right hand side #############'
      do ii=1,k_step
        !b_sub(ii)=dotProd(u(ii),b)
         !write(75,*)dotProd(u(ii),b)
     end do
     b_sub(1)=beta_1
	  write(75,*)b_sub(1),beta_1
	  
! Forward substitution: Ly=b solve for y
write(6,*)'############### Y solution #############'
y_sub(1)=b_sub(1)/L_matrix(1,1)
write(6,*)y_sub(1)
 do ii=2,k_step
 sum=R_zero
    do jj=1,ii-1
      sum=sum+L_matrix(ii,jj)*y_sub(jj)
    end do
    y_sub(ii)=(1/L_matrix(1,1))*(b_sub(ii)-sum)
     write(6,*)y_sub(ii)
 end do

! Backward substitution Ux=y solve for x
write(6,*)'############### X solution #############'
x_sub(k_step)=y_sub(k_step)/U_matrix(k_step,k_step)
write(6,*)x_sub(k_step)
 do ii=k_step-1,1,-1
    sum=R_zero
    do jj=ii+1,k_step
      sum=sum+U_matrix(ii,jj)*x_sub(jj)
    end do
     x_sub(ii)=(1/U_matrix(ii,ii))*(y_sub(ii)-sum)
     write(6,*)x_sub(ii)
 end do
 

 
 10 continue 
 call zero(mhat)

 if(present(dHat) ) then
  call zero(dHat)
   call zero(dHat_temp)
 end if
 k_step_solution=x_sub(k_step)
 
! Model update 
		 do ii=1, k_step
			 call scMult_modelParam (x_sub(ii),q(ii),q_temp)
			 call linComb_modelParam(ONE,mhat,ONE,q_temp,mhat)
			 if(present(dHat) ) then
				 call scMult_dataVectorMTX(x_sub(ii),Au_vec(ii),dHat_temp)
							do i=1,b%nTx 
							  do iDt=1, b%d(i)%nDt 
								  dHat_temp%d(i)%data(iDt)%errorBar= .false.
							  end do
							 end do
				 call linComb_dataVectorMTX(ONE,dHat,ONE,dHat_temp,dHat)
			  end if	 
		 end do
! add to m0	 
         
	     call linComb_modelParam(ONE,m0,ONE,mHat,m)

end subroutine get_model
!****************************************************************************************  

   

subroutine CG_MS(b,x,m,d,lambda,CGiter)


  type (modelParam_t), intent(in)	 	::b
  type (modelParam_t), intent(inout) 	::x
  type(modelParam_t),  intent(in)       ::m
  type (dataVectorMTX_t), intent(in)	    ::d
  real(kind=prec),     intent(inout)       ::lambda
  type(iterControl_t), intent(inout)	:: CGiter
  character(3)         					::iterChar
  
  !Local
    type (modelParam_t)              	:: r,p,Ap
    real(kind=prec)					 	::alpha,beta,r_norm_pre,r_norm,b_norm,error,delta_new,delta_zero,delta_old
    integer                          	::cg_iter,i,j,k,ii,iDt
    
     
 
r=b
p=r
Ap=m
b_norm=dotProd(b,b)
call zero(x)
r_norm=dotProd(r,r)

 1 continue 
ii = 0
      CGiter%rerr(ii) = r_norm/b_norm

loop: do while ((CGiter%rerr(ii).gt.CGiter%tol).and.(ii.lt.CGiter%maxIt))

             
! Compute matrix-vector product A*p and save the result in Ap  
       call MultA_MS(p,m,d,lambda,Ap)   

                       
! Compute alpha: alpha= (r^T r) / (p^T Ap)    
       alpha = r_norm/dotProd(p,Ap)
       
! Compute new x: x = x + alpha*p           
       call linComb(ONE,x,alpha,p,x)                       
! Compute new r: r = r - alpha*Ap   
       call linComb(ONE,r,-alpha,Ap,r) 
        
                
        r_norm_pre=r_norm
        r_norm=dotProd(r,r)
! Compute beta: beta= r_norm /r_norm_previous           
        beta=r_norm/r_norm_pre
   
! Compute new p: p = r + beta*p    
          call linComb(ONE,r,beta,p,p)
          
       ii=ii+1
       CGiter%rerr(ii) = r_norm/b_norm 
       write(6,*) 'CG-error',ii, r_norm/b_norm
  end do loop
if (CGiter%rerr(ii).gt.CGiter%tol) then
lambda=lambda*2
goto 1
end if

CGiter%niter = ii

! deallocate the help vectors
    call deall_modelParam(r)
    call deall_modelParam(p)
    call deall_modelParam(Ap)
    
end subroutine CG_MS
!****************************************************************************************
 
subroutine CG_DS(b,x,m,m0,d,lambda,CGiter,d_Pred_Orig,RMS_Full,mhat)


  type (dataVectorMTX_t), intent(in)	 	::b
  type (dataVectorMTX_t), intent(inout) 	::x
  type(modelParam_t),  intent(inout)       ::m
  type(modelParam_t),  intent(in)       ::m0  
  type(modelParam_t),  intent(out)       ::mhat    
  type (dataVectorMTX_t), intent(in)	    ::d,d_Pred_Orig
  real(kind=prec),     intent(inout)       ::lambda
  real(kind=prec),     intent(in)       :: RMS_Full
  type(iterControl_t), intent(inout)	:: CGiter
  character(3)         					::iterChar
  type (dataVectorMTX_t)                ::JJT_P,JJT_P_N,sum_data_vec,sum_data_vec_x
  
  !Local
    type (dataVectorMTX_t)              	:: r,p,Ap,d_appox1,d_appox2,res,lambdaP,Ap_tri,Ap_tri_old
    real(kind=prec)					 	::target_rms,alpha,beta,r_norm_pre,r_norm,b_norm,error,delta_new,delta_zero,delta_old,RMS,alpha_lambda,lambda_tri,lambda11,old_rms,alpha_lambda_old
    integer                          	::cg_iter,i,j,k,ii,iDt,i_lambda,iii,nTx
    logical                             :: compute_alpha
    real(kind=prec)					 	:: alpha_vec(20),F,mNorm,rms_opt,lambda_opt
    type (dataVectorMTX_t)                ::JJT_P_vec(20),sum_JJT_P_vec,p_vec(20),r_vec(20)
    type(modelParam_t)                   :: CmJTp(20),q_temp,m_temp,mHat_opt

    
    nTx=b%nTx
     allocate(s_hat(nTx))
     target_rms=RMS_Full/2.0
     if (target_rms .lt. 1.0 ) then
       target_rms=1.00
     end if  
    
    d_appox1=d
    d_appox2=d
    lambdaP=d
    Ap_tri=d
    Ap_tri_old=d
    sum_data_vec=d
    sum_data_vec_x=d
    sum_JJT_P_vec=d
    mhat=m0
    q_temp=m0
    m_temp=m0
    mHat_opt=m0
    call zero(m_temp)
    call zero(mHat_opt)
    do i=1,20
        CmJTp(i)=m0
        call zero(CmJTp(i))
    end do    
    call zero_dataVectorMTX(d_appox2)
    call zero_dataVectorMTX(d_appox1) 
    call zero_dataVectorMTX(lambdaP)     
    call zero_dataVectorMTX(Ap_tri)  
    call zero_dataVectorMTX(Ap_tri_old) 
    Ap=d
    b_norm=dotProd(b,b)
	    r=b
        p=r
        p_vec(1)=p
        r_vec(1)=r
	   call zero_dataVectorMTX(x)
	   
    r_norm=dotProd(r,r)


         !call linComb(R_ZERO,d,ONE,r,r) 
         !call linComb(R_ZERO,d,ONE,p,p)
         !call linComb(R_ZERO,d,ONE,x,x) 
         !call linComb(R_ZERO,d,ONE,Ap,Ap)  

 1 continue                 
ii = 1
compute_alpha=.true.
CGiter%rerr(ii) = r_norm/b_norm
 write(10,*) 'CG-error',ii, r_norm/b_norm
 	  do i=1,nTx
	  	s_hat(i)=m
	    call zero(s_hat(j))
      end do
      rms_opt=1000.0
loop: do while ((CGiter%rerr(ii).gt.CGiter%tol).and.(ii.le.CGiter%maxIt))

             
! Compute matrix-vector product A*p and save the result in Ap  
       
       call MultA_DS(p_vec(ii),m,d,lambda,Ap,JJT_P,CmJTp(ii))
       !call MultA_JJT_DS(p,m,d,Ap,CmJTp(ii+1))
       JJT_P_vec(ii)=JJT_P

                          !call zero_dataVectorMTX(sum_JJT_P_vec) 
                          !do iii=1, ii+1
                          !   Call scMultAdd_dataVectorMTX(ONE,JJT_P_vec(iii),sum_JJT_P_vec)  
                          !end do
                          
! try the approximated FWD solution
!if (ii== 0) then
   		 lambda11=4
		 lambda_tri=10**(lambda11)
         old_rms=10000.0 
 
                   do i_lambda=1, 20
                              call zero_dataVectorMTX(sum_data_vec)
                              call zero_dataVectorMTX(sum_data_vec_x)
                                  mhat=m0
                                  call zero(mhat)
                        do iii=1, ii
                            call scMult_dataVectorMTX(lambda_tri,p_vec(iii),lambdaP)
                            !call scMult_dataVectorMTX(lambda_tri,p,lambdaP)
                            JJT_P_N=JJT_P_vec(iii)
                            call un_normalize_with_dataVecMTX(JJT_P_N,d,1)
                             do i=1,lambdaP%nTx
                               do iDt=1,lambdaP%d(i)%nDt
                                 lambdaP%d(i)%data(iDt)%errorBar= .false.
                                 JJT_P_N%d(i)%data(iDt)%errorBar= .false.
                                 sum_data_vec%d(i)%data(iDt)%errorBar= .false.
                                 JJT_P_vec(iii)%d(i)%data(iDt)%errorBar= .false.
                                end do
                             end do
                          
                           call linComb_dataVectorMTX(ONE,JJT_P_vec(iii),ONE,lambdaP,Ap_tri)   
                           ! Compute alpha_tri: alpha= (r^T r) / (p^T Ap) 
                           alpha_lambda =dotProd(r_vec(iii),r_vec(iii))/dotProd(p_vec(iii),Ap_tri)
                           Call scMultAdd_dataVectorMTX(alpha_lambda,p_vec(iii),sum_data_vec_x) ! Get X
                           Call scMultAdd_dataVectorMTX(alpha_lambda,JJT_P_N,sum_data_vec)     ! Get JJT del m
                           call linComb_modelParam(ONE,mhat,alpha_lambda,CmJTp(iii),mhat)      ! Get mHat
                         end do  

                          d_appox1=sum_data_vec  
                             do i=1,lambdaP%nTx
                               do iDt=1,lambdaP%d(i)%nDt
                                 d_appox1%d(i)%data(iDt)%errorBar= .false.
                                end do
                             end do
                          
                          call linComb_dataVectorMTX(ONE,d_Pred_Orig,ONE,d_appox1,d_appox2) 
                          call compute_RMS(m,d,d_appox2,Lambda,RMS,res)	
                          write(110,'(a10,i5,a12,f10.5,a9,f10.5)') 'CG iter:', ii, 'RMS_approx=', RMS, 'Lambda=',lambda_tri
                          !call linComb_modelParam(ONE,m0,ONE,mHat,m_temp)
                          !call Calc_FWD(lambda,d,m_temp,d_appox2,res,eAll,F,mNorm,rms)
                          !call compute_RMS(m,d,d_appox2,Lambda,RMS,res)	
                          !write(110,'(a10,i5,a6,f10.5,a9,f10.5)') 'CG iter:', ii, 'RMS_Full=', RMS, 'Lambda=',lambda_tri
                          !if (rms .lt. target_rms) then
                          !        lambda=lambda_tri
                          !        write(110,*)
                          !        goto 999
                          !end if
                          if (rms .lt. rms_opt) then
                              rms_opt=rms
                              mHat_opt=mHat
                              lambda_opt=lambda_tri
                              !lambda=lambda_tri
                              compute_alpha=.true.
                              !alpha=alpha_lambda
                              !Ap=Ap_tri
                              !x=sum_data_vec_x
                              
                          end if  
                          if (rms .lt. target_rms) then
                                  lambda=lambda_tri
                                  write(110,*)
                                  goto 999
                          end if
                         

                          
                          
                          old_rms= rms

                                
                           
                            lambda11=lambda11-0.25
		                    lambda_tri=10**lambda11
                   end do
                   380  continue 
                   write(110,*)
                    call zero_dataVectorMTX(d_appox1) 
!end if

                   
                   
                   
                         
                          
       
       
       
         do i=1,x%nTx 
          do iDt=1, x%d(i)%nDt 
	          r%d(i)%data(iDt)%errorBar= .false.
	          p%d(i)%data(iDt)%errorBar= .false.
	          x%d(i)%data(iDt)%errorBar= .false.
	          Ap%d(i)%data(iDt)%errorBar= .false.
          end do
         end do  
           write(110,*)'Lambda=',lambda                    
! Compute alpha: alpha= (r^T r) / (p^T Ap)
   if (compute_alpha) then
       alpha = r_norm/dotProd(p,Ap)
   end if    
       alpha_vec(ii)=alpha
       
       
        call zero(mhat)
        q_temp=m0
         call zero(q_temp)
       
 		 do iii=1, ii
			 call scMult_modelParam (alpha_vec(iii),CmJTp(iii),q_temp)
			 call linComb_modelParam(ONE,mhat,ONE,q_temp,mhat)
         end do    
         !call linComb_modelParam(ONE,m0,ONE,mHat,m)
     
     
! Compute new x: x = x + alpha*p         
       Call scMultAdd_dataVectorMTX(alpha,p,x)  
                        
! Compute new r: r = r - alpha*Ap   
       Call scMultAdd_dataVectorMTX(-alpha,Ap,r) 
       r_vec(ii+1)=r

        
                
        r_norm_pre=r_norm
        r_norm=dotProd(r,r)
! Compute beta: beta= r_norm /r_norm_previous           
        beta=r_norm/r_norm_pre
   
! Compute new p: p = r + beta*p    
          call linComb(ONE,r,beta,p,p)
          p_vec(ii+1)=p

                          
       ii=ii+1
       CGiter%rerr(ii) = r_norm/b_norm 
       write(55,'(a10,i5,2x,es12.6,2x,f10.5)') 'CG-error',ii, r_norm/b_norm,lambda
	   write(65,*) keep_solution
       write(10,*) 'CG-error',ii, r_norm/b_norm

       !write(6,*) 'Beta_CG',ii, sqrt(beta)/alpha
  end do loop

999 continue 
    write(110,*)'DONE with this Iter 1', ii
!keep_solution= .true.
!x_previous=x
    mHat=mHat_opt
    lambda=lambda_opt
   write(110,*)'DONE with this Iter 2', ii
CGiter%niter = ii

! deallocate the help vectors
    
     !call deall_dataVectorMTX(r)
       write(110,*)'DONE with this Iter 2_1', ii
    call deall_dataVectorMTX(p)
           write(110,*)'DONE with this Iter 2_2', ii
    call deall_dataVectorMTX(Ap)
           write(110,*)'DONE with this Iter 2_3', ii
    call deall_dataVectorMTX(d_appox1)    
           write(110,*)'DONE with this Iter 2_4', ii
    call deall_dataVectorMTX(d_appox2)   
           write(110,*)'DONE with this Iter 2_5', ii
    call deall_dataVectorMTX(JJT_P)  
       write(110,*)'DONE with this Iter 3', ii
       
end subroutine CG_DS
!###################################################################################
subroutine MultA_MS(p,m,d,lambda,Ap)

   type(modelParam_t), intent(in)          ::p
   type(modelParam_t), intent(out)         ::Ap
   type(modelParam_t), intent(in)          ::m
   type (dataVectorMTX_t), intent(in)	       ::d
   real(kind=prec), intent(in)             ::lambda
!Local parameters
   type(dataVectorMTX_t)                      ::Jp
   type(modelParam_t)                      ::lambdaP,JTCdJp
   integer                                 ::i,j,k,iDt

Jp		=d
JTCdJp	=m
lambdaP	=m


! Compute   J p 
#ifdef MPI
            call Master_job_Jmult(p,m,Jp,eAll)
#else
            call Jmult(p,m,Jp,eAll)
#endif

! Compute Cd  J p 
         do i=1,Jp%nTx
            do iDt=1,Jp%d(i)%nDt
             do j=1,Jp%d(i)%data(iDt)%nSite
               do k=1,Jp%d(i)%data(iDt)%nComp
                      Jp%d(i)%data(iDt)%value(k,j)=  (Jp%d(i)%data(iDt)%value(k,j)/d%d(i)%data(iDt)%error(k,j)**2)
              end do
            end do
           end do                                       
        end do 

! Compute JT Cd  J p                   
#ifdef MPI
            call Master_job_JmultT(m,Jp,JTCdJp,eAll)
#else
            call JmultT(m,Jp,JTCdJp,eAll)
#endif

            !lambdaP= multBy_Cm_Inv(p)
            call scMult(lambda,p,lambdaP)
            

           call linComb(ONE,JTCdJp,ONE,lambdaP,Ap)  
               
        
        
!Deallocate help vectors
    call deall_modelParam(JTCdJp)
    call deall_modelParam(lambdaP)
    call deall_dataVectorMTX(Jp)

    
    
                
end subroutine MultA_MS
!###################################################################################



!###################################################################################
subroutine MultA_DS(p,m,d,lambda,Ap,JJT_P,CmJTp)
   type(dataVectorMTX_t), intent(in)          ::p
   type(dataVectorMTX_t), intent(out)         ::Ap
   type(dataVectorMTX_t),optional, intent(inout)         ::JJT_P
   type(modelParam_t), intent(in)          ::m
   type (dataVectorMTX_t), intent(in)	       ::d
   real(kind=prec), intent(in)             ::lambda
   type(modelParam_t) ,intent(out)         ::CmJTp
!Local parameters
   type(modelParam_t)                      ::JTp
   type(dataVectorMTX_t)                      ::lambdaP,p_temp
   integer                                 ::i,j,k,iDt

JTp		=m
CmJTp	=m
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
            CmJTp= multBy_Cm(JTp)        
! Compute J Cm  J^T  Cd^(-1/2) p = Ap 
Ap=d    
#ifdef MPI
            call Master_job_Jmult(CmJTp,m,Ap,eAll)
#else
            call Jmult(CmJTp,m,Ap,eAll)
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
!###################################################################################
subroutine MultA_JJT_DS(p,m,d,Ap,CmJTp)
   type(dataVectorMTX_t), intent(in)          ::p
   type(dataVectorMTX_t), intent(out)         ::Ap
   type(modelParam_t), intent(in)          ::m
   type (dataVectorMTX_t), intent(in)	       ::d
   type(modelParam_t) ,intent(out)         ::CmJTp
!Local parameters
   type(modelParam_t)                      ::JTp,CmJTp1
   type(dataVectorMTX_t)                      ::p_temp
   integer                                 ::i,j,k,iDt

JTp		=m
CmJTp	=m
p_temp	=p


!  Mult Cd^(-1/2) p  
         call normalize_with_dataVecMTX(p_temp,d,1)
          
! Compute   J^T   p                   
#ifdef MPI
            call linComb(R_ZERO,d,ONE,p_temp,p_temp)
            call Master_job_JmultT(m,p_temp,JTp,eAll,s_hat)

#else
            call JmultT(m,p_temp,JTp,eAll)
#endif
! Compute  Cm  J^T   p 
            CmJTp= multBy_Cm(JTp)     
! Compute J Cm  J^T   p = Ap 
Ap=d    
#ifdef MPI
            call Master_job_Jmult(CmJTp,m,Ap,eAll)
#else
            call Jmult(CmJTp,m,Ap,eAll)
#endif

            
!Normalize: C^(-1/2)*Ap
			call normalize_with_dataVecMTX(Ap,d,1)

    
        
        
!Deallocate help vectors
     call deall_modelParam(JTp)
    call deall_dataVectorMTX(p_temp)

    
    
                
end subroutine MultA_JJT_DS

!###################################################################################


   subroutine printf(comment,lambda,rms,mNorm,F)

   ! print some comments, rms, f and lambda
  character(*), intent(in)               :: comment
  real(kind=prec), intent(inout)  :: lambda, rms
  real(kind=prec), intent(in), optional:: mNorm
  real(kind=prec), intent(in), optional:: F
  
		write(*,'(a20)',advance='no') trim(comment)//':'
		write(*,'(a5,f11.6)',advance='no') ' rms=',rms
		if (present(mNorm)) then
	    write(*,'(a4,es12.6)',advance='no') ' m2=',mNorm
	    end if
		if (present(F)) then
	    write(*,'(a3,es12.6)',advance='no') ' F=',F
	    end if	    
		write(*,'(a8,f11.6)') ' lambda=',lambda

   end subroutine printf
 !**********************************************************************
   subroutine Calc_FWD(lambda,d,m,d_Pred,res,eAll,F,mNorm,rms,RMS_MT,RMS_CSEM)

   ! Compute the full penalty functional F
   ! Also output the predicted data and the EM solution
   ! that can be used for evaluating the gradient
!Input
   real(kind=prec),    intent(in)           :: lambda
   type(dataVectorMTX_t), intent(in)           :: d
   type(modelParam_t), intent(in)           :: m
!Output   
   real(kind=prec),    intent(out)          :: F, mNorm
   type(dataVectorMTX_t), intent(inout)        :: d_Pred,res
   type(solnVectorMTX_t),  intent(inout)          :: eAll
    real(kind=prec), intent(inout)              :: rms
	real(kind=prec),optional, intent(out)     :: RMS_MT,RMS_CSEM

   !  local variables
   type(dataVectorMTX_t)    :: Nres,MT_d_temp,CSEM_d_temp,MT_d_N,CSEM_d_N
	real(kind=prec)                                 :: RMS_total
   real(kind=prec) :: SS
   integer :: Ndata


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
   mNorm = sqrt(dotProd(m,m))

   ! penalty functional = sum of squares + scaled model norm
   F = SS + (lambda * mNorm)

   ! if required, compute the Root Mean Squared misfit
   	RMS = sqrt(SS/Ndata)
	
	    if (present(RMS_MT)) then
		               call split_dataVectorMTX(res)
					   MT_d_temp= MT_d
					   CSEM_d_temp= CSEM_d
					   call deall_dataVectorMTX(MT_d)
					   call deall_dataVectorMTX(CSEM_d)
					  call split_dataVectorMTX(Nres)  
					  MT_d_N= MT_d
					  CSEM_d_N=CSEM_d
					   write(40,*) 'Before del'
						 call deall_dataVectorMTX(MT_d)
						 call deall_dataVectorMTX(CSEM_d) 
						 write(40,*) 'after del'
				!Compute RMS value for each data set
						 RMS_MT   = sqrt(dotProd(MT_d_temp,MT_d_N)/Ndata_MT)						 
						 RMS_CSEM = sqrt(dotProd(CSEM_d_temp,CSEM_d_N)/Ndata_CSEM)						 
						 RMS_total= sqrt(dotProd(res,Nres)/(Ndata_CSEM+Ndata_MT)) 
	     end if
		 

   call deall_dataVectorMTX(Nres)

                  write(40,*) 'ENd Calc_FWD'

   end subroutine Calc_FWD
   !**********************************************************************
   subroutine compute_RMS(m,d,d_pred,lambda,RMS,res,RMS_MT,RMS_CSEM)

   !Input
    type(modelParam_t), intent(in)           :: m
   type(dataVectorMTX_t), intent(in)              :: d,d_pred
   real(kind=prec),    intent(in)           :: lambda
   !Output   
   real(kind=prec),    intent(out)          :: RMS
    real(kind=prec),optional,    intent(out)         :: RMS_MT,RMS_CSEM
   type(dataVectorMTX_t),  intent(out)              :: res
   
   !Local
    type(dataVectorMTX_t)                           :: Nres,MT_d_temp,CSEM_d_temp,MT_d_N,CSEM_d_N
	real(kind=prec)                                 :: SS,F,mNorm,RMS_total
	Integer                                         :: Ndata
	
    res=d
   call linComb(ONE,d,MinusONE,d_Pred,res)
   Nres=res
   call normalize_dataVectorMTX(Nres,2)
   SS = dotProd(res,Nres)
   Ndata = countData(res)

   ! compute the model norm
   mNorm = sqrt(dotProd(m,m))

   ! penalty functional = sum of squares + scaled model norm
   F = SS + (lambda * mNorm)

   ! if required, compute the Root Mean Squared misfit
   	RMS = sqrt(SS/Ndata)	
	
	   write(50,*) 'Before split1'
	   call split_dataVectorMTX(res)
       
	   	   write(50,*) 'After split1'
       if (present(RMS_MT)) then
		    if (Ndata_MT .ne. 0) then
				MT_d_temp= MT_d
				call deall_dataVectorMTX(MT_d)
				call deall_dataVectorMTX(CSEM_d)
					   write(50,*) 'Before split2'
				call split_dataVectorMTX(Nres)  
					   	   write(50,*) 'After split2'
				MT_d_N= MT_d
				call deall_dataVectorMTX(MT_d)
				call deall_dataVectorMTX(CSEM_d)
				RMS_MT   = sqrt(dotProd(MT_d_temp,MT_d_N)/Ndata_MT)	
				call deall_dataVectorMTX(MT_d_temp)
				call deall_dataVectorMTX(MT_d_N)				
			end if
	   end if
	   
	   call split_dataVectorMTX(res)	   
       if (present(RMS_CSEM)) then
            if (Ndata_CSEM .ne. 0) then	   
			   CSEM_d_temp= CSEM_d
			   call deall_dataVectorMTX(MT_d)
               call deall_dataVectorMTX(CSEM_d)
               call split_dataVectorMTX(Nres)  
               CSEM_d_N=CSEM_d
               call deall_dataVectorMTX(CSEM_d)
               call deall_dataVectorMTX(MT_d)			   
			   RMS_CSEM = sqrt(dotProd(CSEM_d_temp,CSEM_d_N)/Ndata_CSEM)
	        end if             
	    end if
	
   end subroutine compute_RMS
!**********************************************************************
  subroutine un_normalize_with_dataVecMTX(d_in_out,d,N,scale_MT)
  
     type(dataVectorMTX_t), intent(in)              :: d
     type(dataVectorMTX_t), intent(inout)           :: d_in_out
     integer, optional, intent(in)               :: N
	 real(kind=prec),optional,intent(in)			     :: scale_MT
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
 subroutine write_output_files(Iter_number,data,model,file_name_suffix)
 
   type(dataVectorMTX_t), intent(in)              :: data
   type(modelParam_t), intent(in)              :: model
   integer,            intent(in)              :: Iter_number
   character(100), intent(in)   			   :: file_name_suffix
   
   
   character(100)       		:: modelFile,dataFile
   type(iterControl_t)			:: CGiter
   character(3)        			:: iterChar
   
   
  	   	  write(iterChar,'(i3.3)') Iter_number
	   	  modelFile = trim(file_name_suffix)//'_'//iterChar//'.cpr'
	      call write_modelParam(model,trim(modelFile))
    

	       dataFile =trim(file_name_suffix)//'_'//iterChar//'.imp'
           call write_dataVectorMTX(data,trim(dataFile))
  end subroutine write_output_files
  
  !******************************************************************************
!  Given a matrix A, with logical dimensions M by N and physical dimensions MP by NP, this
!  routine computes its singular value decomposition, A = U.W.V'. The matrix U replaces 
!  A on output. The diagonal matrix of singular values W is output as a vector W. The matrix 
!  V (not the Transpose V') is output as V. M must be greater or equal to N; if it is smaller,
!  then A should be filled up to square with zero rows.
!******************************************************************************
!  Verified: working on 6 Feb 1998
!
! NOTES: 
! 1. When used to solve Linear Equations with n equations and n unknowns,
!     mp=np=m=n.
! 2. When finding inverses, n = soln and b = IdentityMatrix(n)
!
! Modifications to Original Program:
! 1. Equalties/Inequalities are replaced by Differences and compared with EPS
! 2. The Matrix U in U.W.V' is actually stored in "a"
!
! DEBUG:
! 1.  "IMPLICIT NONE" and "REAL(DBP)"
! 2.  Parameters might need to be larger
!******************************************************************************
SUBROUTINE SVDCMP(a,m,n,mp,np,w,v)
   implicit none
   INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)
   
   INTEGER :: m,n,mp,np, i,j,l,k, its, nm, jj
   REAL(dbp) :: g, sscale, anorm, s, f,h,x,z,y,c

   REAL(dbp), parameter :: EPS = 3.0d-15

   INTEGER, PARAMETER :: nmax = 1000        !Maximum anticipated value of N
   REAL(dbp) :: A(mp,np), w(np), v(np,np), rv1(nmax)


!   PRINT *, 'Enter the tolerance or precision required'
!   READ *, eps
   !eps=10E-12

   if (m.lt.n) then
   	PRINT *, 'You must augment A with extra zero rows'
      call exit(10)
   ENDIF 

            !Householder Reduction to bidiagonal form
!(see Forsythe,Malcolm,Moler, "Computer Methods for Mathematical Computations"
   g=0.0d0
   sscale = 0.0d0
   anorm = 0.0d0
   do i = 1,n
      l = i + 1
      rv1(i) = sscale*g
      g = 0.0d0
      s = 0.0d0
      sscale = 0.0d0
      if (i.le.m) then
         do k = i,m
            sscale = sscale + dABS(a(k,i))
         end do       ! k loop
!         if (sscale.ne.0.0d0) then
			if ( dabs(sscale-0.0d0).gt.EPS ) then
            do k = i,m
               a(k,i) = a(k,i) / sscale
               s = s + a(k,i)*a(k,i)
            end do    ! k loop
            f = a(i,i)
            g = - SIGN(SQRT(s),f)
            h = f*g - s
            a(i,i) = f - g
            if (i.ne.n) then
               do j = l,n
                  s = 0.0d0
                  do k = i,m
                     s = s + a(k,i)*a(k,j)
                  end do      ! k loop
                  f = s / h
                  do k = i, m 
                     a(k,j) = a(k,j) + f*a(k,i)
                  end do   ! k loop
               end do      ! j loop
            end if
            do k = i, m 
               a(k,i) = sscale * a(k,i)
            end do         ! k loop
         end if
      end if

      w(i) = sscale * g
      g = 0.0d0
      s = 0.0d0
      sscale = 0.0d0
      if ((i.le.m).AND.(i.ne.n)) then
         do k = l, n
            sscale = sscale + dABS(a(i,k))
         end do         ! k loop
!         if (sscale.ne.0.0d0) then
			if ( dabs(sscale-0.0d0).gt.EPS ) then
            do k = l, n
               a(i,k) = a(i,k) /sscale
               s = s + a(i,k) * a(i,k)
            end do      ! k loop 
            f = a(i,l) 
            g = - SIGN(SQRT(s),f)
            h = f * g - s
            a(i,l) = f - g
            do k = l, n
               rv1(k) = a(i,k) / h
            end do      ! k loop
            if (i.ne.m) then
               do j = l, m 
                  s = 0.0d0
                  do k = l, n 
                     s = s + a(j,k)*a(i,k)
                  end do   ! k loop
                  do k = l, n 
                     a(j,k) = a(j,k) + s*rv1(k)
                  end do   ! k loop
               end do      ! j loop
            end if
				do k = l, n
               a(i,k) = sscale * a(i,k)
           	end do
         end if
      end if
      anorm = MAX(anorm, (dABS(w(i)) + dABS(rv1(i))))
   end do

! Accumulation of right-hand Transformations
   do i = n, 1, -1 
      if (i.lt.n) then
!         if (g.ne.0.0d0) then
			if ( dabs(g-0.0d0).gt.EPS ) then
            do j = l, n       ! Double division to avoid possible overflow
               v(j,i) = (a(i,j) / a(i,l)) / g
            end do      ! j loop
            do j = l, n
               s = 0.0d0
               do k = l, n
                  s = s + a(i,k)*v(k,j)
               end do   ! k loop
               do k = l, n
                  v(k,j) = v(k,j) + s * v(k,i)
               end do   ! k loop
           	end do      ! j loop
         end if
         do j = l, n 
            v(i,j) = 0.0d0
            v(j,i) = 0.0d0
         end do
      end if
      v(i,i) = 1.0d0
      g = rv1(i)
      l = i
   end do

! Accumulation of left-hand Transformations
   do i = n, 1, -1
      l = 1 + i
      g = w(i)
      if (i.lt.n) then
         do j = l, n
            a(i,j) = 0.0d0
         end do
      end if
!      if (g.ne.0.0d0) then
      if ( dabs(g-0.0d0).gt.EPS ) then
         g = 1.0d0 / g
         if (i.ne.n) then
            do j = l,n 
               s = 0.0d0
               do k = l, m
                  s = s + a(k,i)*a(k,j)
               end do   ! k loop
               f = (s/a(i,i)) * g
               do k = i, m 
                  a(k,j) = a(k,j) + f * a(k,i)
               end do   ! k loop
            end do      ! j loop
         end if
         do j = i, m 
            a(j,i) = a(j,i) * g
         end do         ! j loop
      else
         do j = i, m
            a(j,i) = 0.0d0
         end do         ! j loop
      end if
      a(i,i) = a(i,i) + 1.0d0
   end do               ! i loop

! Diagonalization of the bidigonal form
   do k = n, 1, -1                  !Loop over singular values
      do its = 1,30                 !Loop over allowed iterations
         do l = k, 1, -1            !Test for splitting
            nm = l - 1              ! Note that rv1(1) is always zero
!           if ( (dABS(rv1(l))+anorm) .eq. anorm ) GO TO 2
!          	if ( (dABS(w(nm))+anorm) .eq. anorm ) GO TO 1
            if ( dabs((dABS(rv1(l))+anorm) - anorm).lt.eps ) GO TO 2
          	if ( dabs((dABS(w(nm))+anorm) - anorm).lt.eps ) GO TO 1
         end do      !  l loop

1        c = 0.0d0                  ! Cancellation of rv1(l), if l>1 :
         s = 1.0d0
         do i = l, k
            f = s * rv1(i)
!            if ( (dABS(f)+anorm) .ne. anorm ) then
            if ( dabs( (dABS(f)+anorm) - anorm) .GT. eps ) then

               g = w(i)
               h = SQRT(f*f + g*g)
               w(i) = h
               h = 1.0d0 / h
               c = g * h
               s = -f * h
               do j = 1, m
                  y = a(j,nm)
                  z = a(j,i)
                  a(j,nm) = (y*c) + (z*s)
                  a(j,i) = -(y*s) + (z*c)
               end do   ! j loop
            end if
         end do         ! i loop
2        z = w(k) 
         if (l .eq. k) then         ! convergence
				if (z .lt. 0.0d0) then  ! Singular value is made non-negative
               w(k) = -z
               do j = 1,n
                  v(j,k) = -v(j,k)
               end do         ! j loop
	    		end if
            GO TO 3
         end if
         if (its.eq.30) then
         	PRINT*, 'No Convergence in 30 iterations'
            call exit(10)
         ENDIF
         x = w(l)          ! Shift from bottom 2-by-2 minor
         nm = k - 1
         y = w(nm)
         g = rv1(nm)
         h = rv1(k)
         f = ( (y-z)*(y+z) + (g-h)*(g+h) ) / ( 2.0d0*h*y)
         g = SQRT(f*f + 1.0d0)
         f = ( (x-z)*(x+z) + h*((y/(f+SIGN(g,f))) - h) ) / x

! Next   QR Transformation
         c = 1.0d0
         s = 1.0d0
         do j = l, nm
          	i = j + 1
            g = rv1(i)
            y = w(i)
            h = s*g
            g = c*g
            z = SQRT(f*f + h*h)
            rv1(j) = z
            c = f/z
            s = h/z
            f = (x*c) + (g*s)
            g = -(x*s) + (g*c)
            h = y*s
            y = y*c
            do jj = 1, n
               x = v(jj,j)
               z = v(jj,i)
               v(jj,j) = (x*c) + (z*s)
               v(jj,i) = -(x*s) + (z*c)
           	end do
            z = SQRT(f*f + h*h)
            w(j) = z
!            if (z.ne.0.0d0) then
            if (  dabs(z-0.0d0).gt.eps  ) then
               z = 1.0d0 / z
               c = f*z
               s = h*z
            end if
            f = (g*c) + (y*s)
            x = -(g*s) + (y*c)
            do jj = 1, m
               y = a(jj,j)
               z = a(jj,i)
               a(jj,j) = (y*c) + (z*s)
               a(jj,i) = -(y*s) + (z*c)
            end do
         end do         ! j loop
         rv1(l) = 0.0d0
         rv1(k) = f
         w(k) = x
      end do            ! its loop
3  continue
   end do               ! k loop

   return
END SUBROUTINE svdcmp


subroutine SVDSORT(U, W, V,m,n)
    !' Sort U, V, W  by Singular Values in decending order
    ! '   [meaning W(0) will be greatest singular value]
     INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)   
    real(dbp):: U(m,n),w(m),V(m,n)
    Integer :: i,j,k,m,n   
    real(dbp):: S      


    do i = 1 , N - 1
        K = i                    !Find next highest singular value index k
        S = W(K)
        do J = i + 1, N
            If (W(J) .gt. S) Then
                K = J
                S = W(K)
            End If
        end do
        If (K .ne. i) Then
            W(K) = W(i)          !Swap W(k), W(i)
            W(i) = S
            do J = 1, N       !Swap V(Row i), V(Row k)
                S = V(J, i)
                V(J, i) = V(J, K)
                V(J, K) = S
            end do
            do J = 1, M       !Swap U(Row i), U(Row k)
                S = U(J, i)
                U(J, i) = U(J, K)
                U(J, K) = S
            end do
        End If
        end do
end subroutine SVDSORT        
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine split_dataVectorMTX(d_in)
   type(dataVectorMTX_t), intent(inout)           :: d_in
   
   !LOCAL
   integer                  :: iTx,Tx_MT_counter,Tx_CSEM_counter,counter_MT,counter_CSEM
   character(100)           ::  dataFile
	
	

   
 !Split d into two data: CSEM and MT 
   Tx_MT_counter=0
   Tx_CSEM_counter=0   
   
   do iTx=1,d_in%nTx
	 if (txDict(iTx)%Tx_type=='MT') then
	   Tx_MT_counter=Tx_MT_counter+1
	 elseif (txDict(iTx)%Tx_type=='CSEM') then
	   Tx_CSEM_counter=Tx_CSEM_counter+1
     end if	   
   enddo

if (Tx_CSEM_counter .ne. 0 .and. Tx_MT_counter .ne. 0 ) then   
           call deall_dataVectorMTX(MT_d)
           call deall_dataVectorMTX(CSEM_d)		   
		   call create_dataVectorMTX(Tx_MT_counter,MT_d)
		   call create_dataVectorMTX(Tx_CSEM_counter,CSEM_d)
		   
		   counter_MT=0
		   counter_CSEM=0
		   
		 do iTx=1,d_in%nTx
			  if (txDict(iTx)%Tx_type=='MT') then
				 counter_MT=counter_MT+1
				  call create_dataVector(d_in%d(iTx)%ndt,MT_d%d(counter_MT))
				  call copy_dataVector(MT_d%d(counter_MT), d_in%d(iTx))		 
			elseif (txDict(iTx)%Tx_type=='CSEM') then
				  counter_CSEM=counter_CSEM+1
				  call create_dataVector(d_in%d(iTx)%ndt,CSEM_d%d(counter_CSEM))
				  call copy_dataVector(CSEM_d%d(counter_CSEM), d_in%d(iTx))
			end if
		end do	
				  MT_d%allocated = .true.
				  CSEM_d%allocated = .true.
				  
				 Ndata_MT = countData(MT_d)
				 Ndata_CSEM = countData(CSEM_d)

elseif (Tx_MT_counter .ne. 0 ) then   
		   !call create_dataVectorMTX(Tx_MT_counter,MT_d)
           MT_d=d_in
           MT_d%allocated = .true.
           Ndata_MT = countData(MT_d)		 
elseif (Tx_CSEM_counter .ne. 0 ) then   
		   !call create_dataVectorMTX(Tx_CSEM_counter,CSEM_d)
		   CSEM_d=d_in
		   CSEM_d%allocated = .true.
           Ndata_CSEM = countData(CSEM_d)
end if
		 
end subroutine split_dataVectorMTX
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**********************************************************************
subroutine CG_DS_modified(b,x,m,m0,d,lambda,CGiter,d_Pred_Orig,RMS_Full,mhat)


  type (dataVectorMTX_t), intent(in)	 	::b
  type (dataVectorMTX_t), intent(inout) 	::x
  type(modelParam_t),  intent(inout)       ::m
  type(modelParam_t),  intent(in)       ::m0  
  type(modelParam_t),  intent(out)       ::mhat    
  type (dataVectorMTX_t), intent(in)	    ::d,d_Pred_Orig
  real(kind=prec),     intent(inout)       ::lambda
  real(kind=prec),     intent(in)       :: RMS_Full
  type(iterControl_t), intent(inout)	:: CGiter
  character(3)         					::iterChar
  type (dataVectorMTX_t)                ::JJT_P,JJT_P_N,sum_data_vec,sum_data_vec_x
  
  !Local
    type (dataVectorMTX_t)              	:: r,p,Ap,d_appox1,d_appox2,res,lambdaP,Ap_tri,Ap_tri_old
    real(kind=prec)					 	::f_opt,target_rms,alpha,beta,r_norm_pre,r_norm,b_norm,error,delta_new,delta_zero,delta_old,RMS,alpha_lambda,lambda_tri,lambda11,old_rms,alpha_lambda_old
    integer                          	::cg_iter,i,j,k,ii,iDt,i_lambda,iii,nTx
    logical                             :: compute_alpha
    real(kind=prec)					 	:: alpha_vec(20),F,mNorm,rms_opt,lambda_opt
    type (dataVectorMTX_t)                ::JJT_P_vec(20),sum_JJT_P_vec,p_vec(21),r_vec(21)
    type(modelParam_t)                   :: CmJTp(20),q_temp,m_temp,mHat_opt
    real(kind=prec)              ::mNorm_prev,ss,ss_opt

    
    nTx=b%nTx
     target_rms=RMS_Full/2.0
     if (target_rms .lt. 1.0 ) then
       target_rms=1.00
     end if  
    
    d_appox1=d
    d_appox2=d
    lambdaP=d
    Ap_tri=d
    Ap_tri_old=d
    sum_data_vec=d
    sum_data_vec_x=d
    sum_JJT_P_vec=d
    mhat=m0
    q_temp=m0
    m_temp=m0
    mHat_opt=m0
    call zero(m_temp)
    call zero(mHat_opt)
    do i=1,20
        CmJTp(i)=m0
        call zero(CmJTp(i))
    end do    
    call zero_dataVectorMTX(d_appox2)
    call zero_dataVectorMTX(d_appox1) 
    call zero_dataVectorMTX(lambdaP)     
    call zero_dataVectorMTX(Ap_tri)  
    call zero_dataVectorMTX(Ap_tri_old) 
    Ap=d
    b_norm=dotProd(b,b)
	    r=b
        p=r
        p_vec(1)=p
        r_vec(1)=r
	   call zero_dataVectorMTX(x)
	   
    r_norm=dotProd(r,r)


         !call linComb(R_ZERO,d,ONE,r,r) 
         !call linComb(R_ZERO,d,ONE,p,p)
         !call linComb(R_ZERO,d,ONE,x,x) 
         !call linComb(R_ZERO,d,ONE,Ap,Ap)  

 1 continue                 
ii = 1
compute_alpha=.true.
CGiter%rerr(ii) = r_norm/b_norm
 write(10,*) 'CG-error',ii, r_norm/b_norm
 
      rms_opt=1000.0
      f_opt=1000000000.0
      ss_opt=10E10
loop: do while ((CGiter%rerr(ii).gt.CGiter%tol).and.(ii.le.1))

             
! Compute matrix-vector product A*p and save the result in Ap  
       
       call MultA_DS(p_vec(ii),m,d,lambda,Ap,JJT_P,CmJTp(ii))
       !call MultA_JJT_DS(p,m,d,Ap,CmJTp(ii+1))
       JJT_P_vec(ii)=JJT_P

                          !call zero_dataVectorMTX(sum_JJT_P_vec) 
                          !do iii=1, ii+1
                          !   Call scMultAdd_dataVectorMTX(ONE,JJT_P_vec(iii),sum_JJT_P_vec)  
                          !end do
                          
! try the approximated FWD solution
if (ii .eq. 1) then
   		 lambda11=4
		 lambda_tri=10**(lambda11)
         old_rms=10000.0 
 
                   do i_lambda=1, 20
                              call zero_dataVectorMTX(sum_data_vec)
                              call zero_dataVectorMTX(sum_data_vec_x)
                                  mhat=m0
                                  call zero(mhat)
                        do iii=1, ii
                            call scMult_dataVectorMTX(lambda_tri,p_vec(iii),lambdaP)
                            !call scMult_dataVectorMTX(lambda_tri,p,lambdaP)
                            JJT_P_N=JJT_P_vec(iii)
                            call un_normalize_with_dataVecMTX(JJT_P_N,d,1)
                             do i=1,lambdaP%nTx
                               do iDt=1,lambdaP%d(i)%nDt
                                 lambdaP%d(i)%data(iDt)%errorBar= .false.
                                 JJT_P_N%d(i)%data(iDt)%errorBar= .false.
                                 sum_data_vec%d(i)%data(iDt)%errorBar= .false.
                                 JJT_P_vec(iii)%d(i)%data(iDt)%errorBar= .false.
                                end do
                             end do
                          
                           call linComb_dataVectorMTX(ONE,JJT_P_vec(iii),ONE,lambdaP,Ap_tri)   
                           ! Compute alpha_tri: alpha= (r^T r) / (p^T Ap) 
                           alpha_lambda =dotProd(r_vec(iii),r_vec(iii))/dotProd(p_vec(iii),Ap_tri)
                           Call scMultAdd_dataVectorMTX(alpha_lambda,p_vec(iii),sum_data_vec_x) ! Get X
                           Call scMultAdd_dataVectorMTX(alpha_lambda,JJT_P_N,sum_data_vec)     ! Get JJT del m
                           call linComb_modelParam(ONE,mhat,alpha_lambda,CmJTp(iii),mhat)      ! Get mHat
                         end do  

                          d_appox1=sum_data_vec  
                             do i=1,lambdaP%nTx
                               do iDt=1,lambdaP%d(i)%nDt
                                 d_appox1%d(i)%data(iDt)%errorBar= .false.
                                end do
                             end do
                          
                          call linComb_dataVectorMTX(ONE,d_Pred_Orig,ONE,d_appox1,d_appox2) 
                          mNorm_prev=mNorm
                          call compute_RMS(m,mhat,d,d_appox2,lambda_tri,RMS,res,f,mNorm,ss)	
                          write(110,'(a10,i5,a10,es12.6,a13,f10.5,a4,es12.6,a8,es12.6,a5,es12.6)') 'CG iter:', ii, ' Lambda= ',lambda_tri, ' RMS_approx= ', RMS, ' F= ',f, ' mNorm= ',mNorm, ' SS= ',ss
                          write(110,*)(f)/ss
                          !call linComb_modelParam(ONE,m0,ONE,mHat,m_temp)
                          !call Calc_FWD(lambda,d,m_temp,d_appox2,res,eAll,F,mNorm,rms)
                          !call compute_RMS(m,d,d_appox2,Lambda,RMS,res)	
                          !write(110,'(a10,i5,a6,f10.5,a9,f10.5)') 'CG iter:', ii, 'RMS_Full=', RMS, 'Lambda=',lambda_tri
                          !if (rms .lt. target_rms) then
                          !        lambda=lambda_tri
                          !        write(110,*)
                          !        goto 999
                          !end if
                          
                          if (ss .lt. ss_opt) then
                          !write(110,*)abs(mNorm_prev-mNorm),mNorm_prev
                           !if (abs(mNorm_prev-mNorm) .lt.10.0 ) then
                              ss_opt=ss
                              mHat_opt=mHat
                              lambda_opt=lambda_tri
                              lambda=lambda_tri
                              compute_alpha=.true.
                              !!alpha=alpha_lambda
                              Ap=Ap_tri
                              !x=sum_data_vec_x
                              
                          end if  
                          if (rms .lt. target_rms) then
                                  lambda=lambda_tri
                                  write(110,*)
                                  goto 999
                          end if
                          !

                          
                          
                          old_rms= rms

                                
                           
                            lambda11=lambda11-0.25
		                    lambda_tri=10**lambda11
                   end do
                   380  continue 
                   write(110,*)
                    call zero_dataVectorMTX(d_appox1) 
                   return
end if

                   
                   
                   
                         
                          
       
       
       
         do i=1,x%nTx 
          do iDt=1, x%d(i)%nDt 
	          r%d(i)%data(iDt)%errorBar= .false.
	          p%d(i)%data(iDt)%errorBar= .false.
	          x%d(i)%data(iDt)%errorBar= .false.
	          Ap%d(i)%data(iDt)%errorBar= .false.
          end do
         end do  
           write(110,*)'Lambda=',lambda                    
! Compute alpha: alpha= (r^T r) / (p^T Ap)
   if (compute_alpha) then
       alpha = r_norm/dotProd(p,Ap)
   end if    
       alpha_vec(ii)=alpha
       
       
        call zero(mhat)
        q_temp=m0
         call zero(q_temp)
       
 		 do iii=1, ii
			 call scMult_modelParam (alpha_vec(iii),CmJTp(iii),q_temp)
			 call linComb_modelParam(ONE,mhat,ONE,q_temp,mhat)
         end do    
         !call linComb_modelParam(ONE,m0,ONE,mHat,m)
     
     
! Compute new x: x = x + alpha*p         
       Call scMultAdd_dataVectorMTX(alpha,p,x)  
                        
! Compute new r: r = r - alpha*Ap   
       Call scMultAdd_dataVectorMTX(-alpha,Ap,r) 
       r_vec(ii+1)=r

        
                
        r_norm_pre=r_norm
        r_norm=dotProd(r,r)
! Compute beta: beta= r_norm /r_norm_previous           
        beta=r_norm/r_norm_pre
   
! Compute new p: p = r + beta*p    
          call linComb(ONE,r,beta,p,p)
          p_vec(ii+1)=p

                          
       ii=ii+1
       CGiter%rerr(ii) = r_norm/b_norm 
       write(55,'(a10,i5,2x,es12.6,2x,f10.5)') 'CG-error',ii, r_norm/b_norm,lambda
	   write(65,*) keep_solution
       write(10,*) 'CG-error',ii, r_norm/b_norm

       !write(6,*) 'Beta_CG',ii, sqrt(beta)/alpha
  end do loop

999 continue 
    write(110,*)'DONE with this Iter 1', ii
!keep_solution= .true.
!x_previous=x
    mHat=mHat_opt
    lambda=lambda_opt
   write(110,*)'DONE with this Iter 2', ii
CGiter%niter = ii

! deallocate the help vectors
    
     !call deall_dataVectorMTX(r)
       write(110,*)'DONE with this Iter 2_1', ii
    call deall_dataVectorMTX(p)
           write(110,*)'DONE with this Iter 2_2', ii
    call deall_dataVectorMTX(Ap)
           write(110,*)'DONE with this Iter 2_3', ii
    call deall_dataVectorMTX(d_appox1)    
           write(110,*)'DONE with this Iter 2_4', ii
    call deall_dataVectorMTX(d_appox2)   
           write(110,*)'DONE with this Iter 2_5', ii
    call deall_dataVectorMTX(JJT_P)  
       write(110,*)'DONE with this Iter 3', ii
       
end subroutine CG_DS_modified


!**********************************************************************

end module DCG
