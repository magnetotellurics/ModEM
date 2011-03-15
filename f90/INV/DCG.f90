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
real(kind=prec) :: Desired_rms
Logical         :: search_min_model=.false.
   ! type(EMsolnMTX_t),save    :: eAll

Contains

!**********************************************************************
   subroutine setIterControl(CGiter)
   !  at present this is a simple routine to setup for
   !  iteration control; allowing parameters to be read from
   !   a file would be desireable.

   type(iterControl_t), intent(inout)	:: CGiter
   
   CGiter%maxit = 20
   CGiter%tol = 10E-3
   CGiter%niter = 0
   allocate(CGiter%rerr(CGiter%maxit))

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
   type(dataVectorMTX_t)			:: dHat, b,dx,d_Pred,res,Nres,JmHat
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
   dx		=d
   b		=d
   d_Pred	=d
   res		=d
call zero_dataVectorMTX(JmHat)
call zero_dataVectorMTX(b)

! initialize the CG control parameters
		call setIterControl(CGiter)
   
open(10,file='DCG.log')
! Compute the predicted data for the current model m
        call Calc_FWD(lambda,d,m,d_Pred,res,eAll,F,mNorm,rms)
        file_name_suffix='Start_DCG'
         call write_output_files(1,d_Pred,m,file_name_suffix) 
        call printf('DCG_Start ',lambda,rms,mNorm,F)
        Desired_rms=rms/2.0
        
        write(10,'(a20)',advance='no') trim('DCG_Start ')//':'
 		write(10,'(a5,f11.6)',advance='no') ' rms=',rms
	    write(10,'(a4,es12.6)',advance='no') ' m2=',mNorm
	    write(10,'(a3,es12.6)',advance='no') ' F=',F    
		write(10,'(a8,f11.6)') ' lambda=',lambda     
        
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
	        
	         
	       !  call CG_DS(b,dx,m,d,lambda,CGiter)
	       call Lanczos_DS (b,m,m0,d,lambda,mhat,res,CGiter,DS_iter,rms)
	      ! call Multi_Trans_DS (b,dx,m,m0,d,lambda,mhat,res,CGiter,DS_iter,rms)
	      !   call Lanczos_CG_DS(b,dx,m,d,lambda,CGiter,DS_iter)
	         
	       !  goto 999
	       
	         goto 10
	        call normalize_with_dataVecMTX(dx,d,1)

	 
#ifdef MPI           
	                call Master_job_JmultT(m,dx,mHat,eAll)              
#else
	                call JmultT(m,dx,mHat,eAll)
#endif

	      Cm_mHat=  multBy_Cm(mHat) 
	      mHat=Cm_mHat
	       
	     
	     call linComb_modelParam(ONE,m0,ONE,mHat,m)
	! Compute the predicted data for the current model m
	   rms_old=rms
         call Calc_FWD(lambda,d,m,d_Pred,res,eAll,F,mNorm,rms)
           if (rms .gt. rms_old) then
            !  lambda=lambda+1
           end if
    ! Write output model and data files
    file_name_suffix='DCG_MPI'
    call write_output_files(DS_iter,d_Pred,m,file_name_suffix)     
     ! Print output Information on the screen
    write(iterChar,'(i3.3)') DS_iter
    call printf('DCG_Iter '//iterChar,lambda,rms,mNorm,F)
    
        write(10,'(a20)',advance='no') 'DCG_Iter '//iterChar//':'
 		write(10,'(a5,f11.6)',advance='no') ' rms=',rms
	    write(10,'(a4,es12.6)',advance='no') ' m2=',mNorm
	    write(10,'(a3,es12.6)',advance='no') ' F=',F    
		write(10,'(a8,f11.6)') ' lambda=',lambda     


   10 continue        
 ! Clean temp vectors
end do
999 continue
d=d_Pred
close(10)
 
end subroutine DCGsolver
!****************************************************************************************
subroutine Lanczos_CG_DS(b,x,m,d,lambda,CGiter,DS_iter)


  type (dataVectorMTX_t), intent(in)	 	::b
  type (dataVectorMTX_t), intent(inout) 	::x
  type(modelParam_t),  intent(in)       ::m
  type (dataVectorMTX_t), intent(inout)	    ::d
  real(kind=prec),     intent(in)       ::lambda
  type(iterControl_t), intent(inout)	:: CGiter
  integer,     intent(in)               :: DS_iter
  character(3)         					::iterChar
  
  !Local
    type (dataVectorMTX_t)              	:: r,v,p,v_previous,Av,v_hat
    type (dataVectorMTX_t)              	:: u,alpha_u,Au,beta_u
    type(modelParam_t)			        :: ATu,beta_v
    real(kind=prec)					 	::alpha,beta,delta,sigma,omega,mue,mue_previous,omega_previous
    real(kind=prec)						:: b_norm,r_norm,lambda1
    integer                          	::i,j,iDt,k,ii,i_cg

 
    
      
call zero_dataVectorMTX(x)    
r=b
beta=sqrt(dotProd(r,r))
v=b
call scMult(1/beta,r,v)
p=r  
v_previous=b
call zero_dataVectorMTX(v_previous)
sigma=beta
omega_previous=R_ZERO
mue_previous=ONE

Av=b
v_hat=b
b_norm=dotProd(b,b)


open(10,file='CG.dat')
do ii=0,3
    ! begin Lanczos step
      !compute delta = <Av,v>
      call MultA_DS(v,m,d,R_Zero,Av) 
      delta=(dotProd(Av,v))
      !compute v_hat = Av - delta v - beta v_previous
           do i=1,Av%nTx
            do iDt=1,Av%d(i)%nDt
             do j=1,Av%d(i)%data(iDt)%nSite
               do k=1,Av%d(i)%data(iDt)%nComp
                      v_hat%d(i)%data(iDt)%value(k,j)=  Av%d(i)%data(iDt)%value(k,j)-(delta*v%d(i)%data(iDt)%value(k,j))-(beta*v_previous%d(i)%data(iDt)%value(k,j))
              end do
            end do
           end do                                       
        end do
        !compute beta       
         beta=sqrt(dotProd(v_hat,v_hat))
         v_previous=v
         call scMult(1/beta,v_hat,v)
              ! begin the CG- iterates
                lambda1=lambda
                            call zero_dataVectorMTX(x)    
                            r=b
                            sigma=sqrt(dotProd(r,r))
                            p=r  
                            omega_previous=R_ZERO
                            mue_previous=ONE
        do i_cg=1,5
                delta=delta+lambda1
                 mue=1/((delta-omega_previous)/ mue_previous)     
                 mue_previous=mue
                 omega=(beta*mue)**2
                 omega_previous=omega
                 sigma=-beta*mue*sigma
                 !Update x, r and p
		           do i=1,x%nTx
		            do iDt=1,x%d(i)%nDt
		             do j=1,x%d(i)%data(iDt)%nSite
		               do k=1,x%d(i)%data(iDt)%nComp
		                      x%d(i)%data(iDt)%value(k,j)=  x%d(i)%data(iDt)%value(k,j)+(mue*p%d(i)%data(iDt)%value(k,j))
		                      r%d(i)%data(iDt)%value(k,j)=  sigma*v%d(i)%data(iDt)%value(k,j)
		                      p%d(i)%data(iDt)%value(k,j)=  r%d(i)%data(iDt)%value(k,j)+(omega*p%d(i)%data(iDt)%value(k,j))
		              end do
		            end do
		           end do                                       
		        end do 
		        r_norm=sqrt(dotProd(r,r))
		        write(6,*) 'CG-error, r_norm',i_cg, r_norm/b_norm    
                write(6,*) 'CG-error, Beta',i_cg, beta
                 write(6,*) 'CG-error, Alpha',i_cg, alpha
                 write(6,*) 'CG-error, Sigma',i_cg, sigma
                write(10,*)'CG-error, r_norm: ',i_cg, r_norm   
                write(10,*) 'CG-error, beta: ',i_cg, beta
                write(10,*) 'CG-error, alpha: ',i_cg, alpha
                write(10,*) 'CG-error, Sigma',i_cg, sigma
                
                
               ! lambda1=lambda1/2
     end do


end do  
close(10)
    
    
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
    type(modelParam_t),pointer, dimension(:) :: s_hat
    type(modelParam_t),pointer, dimension(:,:) :: JTw_matrix
    type(modelParam_t)                       :: JTu
    Integer                                  :: k_step,i_k_step,nTx,ii,jj,INFO,counter,counter1,ll,kk,i_lambda
    real(kind=prec)                          :: beta,sum,mu_LU, lambda_LU,lambda_temp,F,mNorm,rms_old
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
  allocate(s_hat(nTx),sts(k_step*nTx,k_step*nTx))
  allocate(sts_temp(k_step*nTx,k_step*nTx),L_matrix(k_step*nTx,k_step*nTx),U_matrix(k_step*nTx,k_step*nTx))
  allocate(b_sub(k_step*nTx),y_sub(k_step*nTx),x_sub(k_step*nTx),uTr(k_step*nTx))
  allocate(AP((k_step*nTx)*((k_step*nTx)+1)/2),sTs_b(k_step*nTx))
  
  


  
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

       
	   beta=sqrt(dotProd(u_hat(1),u_hat(1)))
	   call scMult(1/beta,u_hat(1),u_hat(1))
       write(6,*) 'Beta', beta/sqrt(dotProd(b,b))	
 do i_k_step=1, k_step 

! 1-normalize each sub data vector [ u_hat(1), u_hat(2),...,u_hat(nTx)] with its norm: i.e.  u_hat(1)/||u_hat(1)||

      do ii=1,nTx
	   beta=sqrt(dotProd(u_hat(i_k_step)%d(ii),u_hat(i_k_step)%d(ii)))
	   call scMult(1/beta,u_hat(i_k_step)%d(ii),u_hat(i_k_step)%d(ii))
	  end do       
! 2- Normlize u_hat with the data error Cd^(-1/2):
!	 Mult u_hat  Cd^(-1/2) 

         call normalize_with_dataVecMTX(u_hat(i_k_step),d,1) 

       
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




	  do ii=1,nTx
	  	  JTw_matrix(i_k_step,ii)=s_hat(ii)
	  end do

 !4-   Compute  Cm J^T  Cd^(-1/2) u_hat 
    do ii=1,nTx
      JTw_matrix(i_k_step,ii)= multBy_Cm(JTw_matrix(i_k_step,ii)) 
   end do
   
 !5 - Make the symetric matrix sTs  
   counter1=0
   do ll=1,i_k_step  
     do ii=1,nTx
        counter1=counter1+1
        counter=0
        do kk=1,i_k_step
          do jj=1,nTx
            counter=counter+1
            sts(counter1,counter)=dotProd(JTw_matrix(ll,ii),JTw_matrix(kk,jj))
          end do
        end do  
     end do
   end do      
  
 !6 - Add Lambda to diag(sTs)   
             sts_temp=R_Zero
             sts_temp=sts
			 do ii=1,nTx*i_k_step
			    sts_temp(ii,ii)=sts(ii,ii)+100
			 end do
			   
			do jj=1,nTx*i_k_step
			  do ii=1,jj
               AP(ii + (jj-1)*jj/2) = sts_temp(ii,jj) 
			  end do
			end do  
!7-  Preforme Cholesky decomposition on (sTs +Lambda I) matrix 
            call DPPTRF( 'U', nTx*i_k_step, AP, INFO )
          do jj=1,nTx*i_k_step 
           do ii=1,jj
              L_matrix(ii,jj) = AP(ii + (jj-1)*jj/2) 
			  end do
		  end do
			
            write(6,*)'############ (sTs ) matrix ################'
			  do ii=1,nTx*i_k_step
			   if (i_k_step .eq. 1) then
			    write(6,'(4f10.2)')(sts(ii,jj),jj=1,nTx*i_k_step)
			     elseif (i_k_step .eq. 2) then 
			     write(6,'(8f10.2)')(sts(ii,jj),jj=1,nTx*i_k_step)
			     elseif (i_k_step .eq. 3) then 
			     write(6,'(12f10.2)')(sts(ii,jj),jj=1,nTx*i_k_step)
			     end if
			      
			 end do 
			 write(6,*)'############ (sTs +Lambda I) matrix ################'
			  do ii=1,nTx*i_k_step
			   if (i_k_step .eq. 1) then
			    write(6,'(4f10.2)')(sts_temp(ii,jj),jj=1,nTx*i_k_step)
			   elseif (i_k_step .eq. 2) then 
			    write(6,'(8f10.2)')(sts_temp(ii,jj),jj=1,nTx*i_k_step)
			   elseif (i_k_step .eq. 3) then 
			    write(6,'(12f10.2)')(sts_temp(ii,jj),jj=1,nTx*i_k_step)			    
			   end if 
			   
			 end do 

			 write(6,*)'############### L matrix #############'
			 do ii=1,nTx*i_k_step
			   if (i_k_step .eq. 1) then
			    write(6,'(4f10.2)')(L_matrix(ii,jj),jj=1,nTx*i_k_step)
			   elseif (i_k_step .eq. 2) then
			   write(6,'(8f10.2)')(L_matrix(ii,jj),jj=1,nTx*i_k_step)
			   elseif (i_k_step .eq. 3) then
			   write(6,'(12f10.2)')(L_matrix(ii,jj),jj=1,nTx*i_k_step)			   
			  end if 
			 end do 
! 8- make the right hand side to solve the projected problem:
!    b_sub=uhat(1,...,nTx)*b
 
            counter=0
            do jj=1,i_k_step
			  do ii=1,nTx
			   counter=counter+1
			    b_sub(counter)=(dotProd(u_hat(jj)%d(ii),b%d(ii)))
			 end do
		   end do	  
			 
			 write(6,*)'############### right hand side (b_sub vector) #############'
			 do ii=1,nTx*i_k_step
			    write(6,*) b_sub(ii)
			 end do 			 
!9- solve the problem, save solution in b_sub
  call DPPTRS( 'U', nTx*i_k_step, 1, AP, b_sub, nTx*i_k_step, INFO )	
             write(6,*)'############### Solution for the projected system #############'		 
 			 do ii=1,nTx*i_k_step
			    write(6,*) b_sub(ii)
            end do  

! 10- Model update for the projected system:
!  m_hat=JT b= JT u_hat *b_sub
! Model update 
 call zero(mhat)
 call zero(q_temp)
    counter=0
      do jj=1,i_k_step
		 do ii=1, nTx
		  counter=counter+1
			 call scMult_modelParam (b_sub(counter),JTw_matrix(jj,ii),q_temp)
			 call linComb_modelParam(ONE,mhat,ONE,q_temp,mhat)
			 call scMult(b_sub(counter),u_hat(jj)%d(ii),b_u%d(ii))
		 end do
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

         call normalize_with_dataVecMTX(Jm,d,1) 
         
		 do ii=1,b%nTx 
          do jj=1, b%d(ii)%nDt 
	          Jm%d(ii)%data(jj)%errorBar= .false.
          end do
         end do 	
		!  12b-  (sTs) *b_sub 
		
		       do ii=1,nTx*i_k_step
		         sum=R_zero
		          do jj=1,nTx*i_k_step
		            sum=sum+(sts_temp(ii,jj)*b_sub(jj))
		          end do
		           sTs_b(ii)=sum
		       end do     
		 ! 12c-  u_hat * sTs_b
		counter=0 
		 do jj=1,i_k_step
		  do ii=1,nTx
		    counter=counter+1
		     call scMult(sTs_b(counter),u_hat(jj)%d(ii),u_hat_temp%d(ii))
		 end do
		end do

	
          
   	    do ii=1,b%nTx 
          do jj=1, b%d(ii)%nDt 
	          u_hat_temp%d(ii)%data(jj)%errorBar= .false.
          end do
         end do 
                 	
		  call linComb (ONE,Jm,MinusONE,u_hat_temp,u_hat(i_k_step+1))
	      beta=sqrt(dotProd(u_hat(i_k_step+1),u_hat(i_k_step+1)))
          call scMult(1/beta,u_hat(i_k_step+1),u_hat(i_k_step+1))
                write(6,*) 'Beta', beta/sqrt(dotProd(b,b))		  
		  
		  
          !beta=sqrt(dotProd(u_hat(i_k_step+1),u_hat(i_k_step+1)))
          !call scMult(1/beta,u_hat(i_k_step+1),u_hat(i_k_step+1))
          !write(6,*) 'Beta', sqrt(beta/dotProd(b,b))
          
          	
    goto 1
		  
		  !beta=sqrt(dotProd(u_hat(i_k_step+1),u_hat(i_k_step+1)))
		  !call scMult(1/beta,u_hat(i_k_step+1),u_hat(i_k_step+1))
		  
 
          
          !  r = dd-(dHat+nu*u*b);
          ! u(k+1) = residual-u*u'*residual;
          
		counter=0 
		 do jj=1,i_k_step
		  do ii=1,nTx
		    counter=counter+1
	        beta= (dotProd(u_hat(jj)%d(ii),residual%d(ii)))
	        call scMult(beta,u_hat(jj)%d(ii),u_hat_temp%d(ii))
	      end do  
	  end do 
	  
 	    do ii=1,b%nTx 
          do jj=1, b%d(ii)%nDt 
	          residual%d(ii)%data(jj)%errorBar= .false.
	          u_hat_temp%d(ii)%data(jj)%errorBar= .false.
          end do
         end do 
         
          call linComb (ONE,residual,MinusONE,u_hat_temp,u_hat(i_k_step+1))
          
         ! u_hat(i_k_step+1)= u_hat(i_k_step)
          beta=sqrt(dotProd(u_hat(i_k_step+1),u_hat(i_k_step+1)))
          call scMult(1/beta,u_hat(i_k_step+1),u_hat(i_k_step+1))
          write(6,*) 'Beta', sqrt(beta/dotProd(b,b))    
 1 continue           
 end do 
 ! Check orthogonalty



 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 

 i_k_step= k_step
 write(iterChar,'(i3.3)') DS_iter
 call printf('Start1'//iterChar,lambda,rms,mNorm,F)
 ! Lambda search
lambda=1000
rms=10000

 do i_lambda=1, 20
 
 !6 - Add Lambda to diag(sTs)   
             sts_temp=R_Zero
             sts_temp=sts
			 do ii=1,nTx*i_k_step
			    sts_temp(ii,ii)=sts(ii,ii)+lambda
			 end do
			   
			do jj=1,nTx*i_k_step
			  do ii=1,jj
               AP(ii + (jj-1)*jj/2) = sts_temp(ii,jj) 
			  end do
			end do  
!7-  Preforme Cholesky decomposition on (sTs +Lambda I) matrix 
            call DPPTRF( 'U', nTx*i_k_step, AP, INFO )
          do jj=1,nTx*i_k_step 
           do ii=1,jj
              L_matrix(ii,jj) = AP(ii + (jj-1)*jj/2) 
			  end do
		  end do
			
            write(6,*)'############ (sTs ) matrix ################'
			  do ii=1,nTx*i_k_step
			   if (i_k_step .eq. 1) then
			    write(6,'(4f10.2)')(sts(ii,jj),jj=1,nTx*i_k_step)
			     elseif (i_k_step .eq. 2) then 
			     write(6,'(8f10.2)')(sts(ii,jj),jj=1,nTx*i_k_step)
			     elseif (i_k_step .eq. 3) then 
			     write(6,'(12f10.2)')(sts(ii,jj),jj=1,nTx*i_k_step)
			     end if
			      
			 end do 
			 write(6,*)'############ (sTs +Lambda I) matrix ################'
			  do ii=1,nTx*i_k_step
			   if (i_k_step .eq. 1) then
			    write(6,'(4f10.2)')(sts_temp(ii,jj),jj=1,nTx*i_k_step)
			   elseif (i_k_step .eq. 2) then 
			    write(6,'(8f10.2)')(sts_temp(ii,jj),jj=1,nTx*i_k_step)
			   elseif (i_k_step .eq. 3) then 
			    write(6,'(12f10.2)')(sts_temp(ii,jj),jj=1,nTx*i_k_step)			    
			   end if 
			   
			 end do 

			 write(6,*)'############### L matrix #############'
			 do ii=1,nTx*i_k_step
			   if (i_k_step .eq. 1) then
			    write(6,'(4f10.2)')(L_matrix(ii,jj),jj=1,nTx*i_k_step)
			   elseif (i_k_step .eq. 2) then
			   write(6,'(8f10.2)')(L_matrix(ii,jj),jj=1,nTx*i_k_step)
			   elseif (i_k_step .eq. 3) then
			   write(6,'(12f10.2)')(L_matrix(ii,jj),jj=1,nTx*i_k_step)			   
			  end if 
			 end do 
! 8- make the right hand side to solve the projected problem:
!    b_sub=uhat(1,...,nTx)*b
 
            counter=0
            do jj=1,i_k_step
			  do ii=1,nTx
			   counter=counter+1
			    b_sub(counter)=(dotProd(u_hat(jj)%d(ii),b%d(ii)))
			 end do
		   end do	  
			 
			 write(6,*)'############### right hand side (b_sub vector) #############'
			 do ii=1,nTx*i_k_step
			    write(6,*) b_sub(ii)
			 end do 			 
!9- solve the problem, save solution in b_sub
  call DPPTRS( 'U', nTx*i_k_step, 1, AP, b_sub, nTx*i_k_step, INFO )	
             write(6,*)'############### Solution for the projected system #############'		 
 			 do ii=1,nTx*i_k_step
			    write(6,*) b_sub(ii)
            end do  

! 10- Model update for the projected system:
!  m_hat=JT b= JT u_hat *b_sub
! Model update 
 call zero(mhat)
 call zero(q_temp)
    counter=0
      do jj=1,i_k_step
		 do ii=1, nTx
		  counter=counter+1
			 call scMult_modelParam (b_sub(counter),JTw_matrix(jj,ii),q_temp)
			 call linComb_modelParam(ONE,mhat,ONE,q_temp,mhat)
		 end do
	end do	 
! 11- smooth the model update m_hat		  
	!mhat= multBy_Cm(mhat)
 
         m_temp=m
 	     call linComb_modelParam(ONE,m0,ONE,mHat,m)

	          
	! Compute the predicted data for the current model m
	     rms_old=rms
         call Calc_FWD(lambda,d,m,d_Pred,res,eAll,F,mNorm,rms)
         
         if (rms .gt. rms_old )then
         call printf('Hi_DCG_Iter_sub'//iterChar,lambda,rms,mNorm,F)
         lambda=lambda*2
         m=m_temp
         goto 99
        end if
         
           
    ! Write output model and data files
    file_name_suffix='Multi_Trans_lambda'
    call write_output_files(i_lambda,d_Pred,m,file_name_suffix)  
    write(iterChar,'(i3.3)') i_lambda
    call printf('DCG_Iter_sub'//iterChar,lambda,rms,mNorm,F)
 
 lambda=lambda/2

 
 
 
 end do

  99 continue 
 call Calc_FWD(lambda,d,m,d_Pred,res,eAll,F,mNorm,rms)


 

 
 
end subroutine Multi_Trans_DS  
!****************************************************************************************
subroutine update_sys(k_step,m,d,beta_kstep_plus_1,u,T_matrix,q)

integer,intent(in)                      			    :: k_step
real(kind=prec),intent(inout),dimension(:,:)    	    :: T_matrix
type (modelParam_t),intent(inout), dimension(:) 	     :: q
type (dataVectorMTX_t),intent(inout), dimension(:) 	     :: u
real(kind=prec),intent(inout)                            :: beta_kstep_plus_1

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
subroutine Lanczos_DS(b,m,m0,d,lambda,mhat,res,CGiter,DS_iter,rms)


  type (dataVectorMTX_t), intent(in)	 	 :: b
  type(modelParam_t),  intent(inout)     :: m
  type(modelParam_t),  intent(in)        :: m0
  type(modelParam_t),  intent(inout)     :: mhat
  type(dataVectorMTX_t),  intent(inout)	 :: res 
  type (dataVectorMTX_t), intent(inout)	 :: d
  type(iterControl_t), intent(inout)	 :: CGiter
  integer,intent(in)					 :: DS_iter
  real(kind=prec),intent(inout)			 :: rms,lambda

  
  !Local

    type (modelParam_t),pointer, dimension(:)  :: q,v
    type (dataVectorMTX_t),pointer, dimension(:)  :: u,u_recent
    real(kind=prec)	,pointer, dimension(:,:)   :: T_matrix
    real(kind=prec)					 	       :: beta_1,beta_kstep_plus_1
    integer                          	       :: k_step,i_search,i_sub_search,ii,jj,i_sub_big_search,i_sub_small_search,fwd_calls,min_k_steps

   type(dataVectorMTX_t)			:: d_Pred,b_start,d_start
   character(100)       		:: file_name_suffix
   type(modelParam_t)			:: m_c,m_l,m_r,m_start,m_temp
   type(dataVectorMTX_t)			:: res_c,res_r,res_l
   real(kind=prec)              :: F,mNorm,lambda_c,lambda_r,lambda_l,mNorm_c,mNorm_r,mNorm_l,rms_c,rms_r,rms_l,step_size,sub_step_size,rms_temp  
   logical                      :: go_right,go_left,central,make_big_step,make_small_step,update_proj_sys
   character(3)         		:: iterChar,iterChar1,k_char
   character(100)       		:: modelFile,dataFile
   type(dataVectorMTX_t)			:: w,u_kstep
   type(solnVectorMTX_t)        :: eAll_temp


   
   
min_k_steps=2

    k_step=min_k_steps
    
m_start=m
b_start=b
d_start=d
m_temp=m
update_proj_sys=.false.
eAll_temp=eAll

    


allocate(T_matrix(20,20))
allocate(q(20),v(k_step+1),u_recent(3),u(1:21))

T_matrix=R_zero
do ii=1,20
  q(ii)= m
  call zero(q(ii))
end do
do ii=1,21
  u(ii)= b
  call zero(u(ii))
end do

do ii=1,3
  u_recent(ii)=b
end do
rms_temp=1000.00

i_search=1
i_sub_search=1
i_sub_big_search=1
i_sub_small_search=1
fwd_calls=1


  1 continue  


!call Arnold_bidiag_JJT(b,m,m0,d,k_step,DS_iter,lambda,T_matrix,q,beta_1)
 if (update_proj_sys) then 
     call update_sys (k_step,m_start,d_start,beta_kstep_plus_1,u,T_matrix,q)
 else
      call bidiag_JJT (b_start,m_start,d_start,k_step,DS_iter,lambda,T_matrix,q,beta_1,beta_kstep_plus_1,u)
     ! call Arnold_bidiag_JJT(b_start,m_start,m0,d_start,k_step,DS_iter,lambda,T_matrix,q,beta_1,beta_kstep_plus_1,u)
 end if
 
 do ii=1,10
   	   	write(iterChar,'(i3.3)') DS_iter
   	   	write(k_char,'(i3.3)') ii
	   	  modelFile = 'Q_matrix_'//iterChar// '_' //k_char//'.q'
	     ! call write_modelParam(q(ii),trim(modelFile))    
 end do
 
 do ii=1,10
   	   	write(iterChar,'(i3.3)') DS_iter
   	   	write(k_char,'(i3.3)') ii
	   	  modelFile = 'U_matrix_'//iterChar// '_' //k_char//'.u'
	     ! call write_dataVectorMTX(u(ii),trim(modelFile))
 end do

!call bidiag_1 (b,m,d,k_step,T_matrix,q,v,beta_1) 
!q=v
 write(6,*)DS_iter,' #### T Matrix ########'
			  do ii=1,k_step
			    write(6,*)(T_matrix(ii,jj),jj=1,k_step)
			 end do

d_Pred	=d
call zero(d_Pred)
call zero(mhat)

write(6,*),DS_iter,' Final Ksteps ', k_step


!Start Line search using the projected system



if (Desired_rms .lt. 1.05) then
  Desired_rms=1.05
end if

	   	   go_right=.true.
	       central=.true.
	       go_left=.true.
	       make_big_step=.false.
	       make_small_step=.true.
	       
	       step_size=(0.75/1)
10 continue
  Lambda_c=10**(lambda)
  Lambda_r=10**(Lambda+(step_size))
  Lambda_l=10**(Lambda-(step_size))


!Central Lambda
if (central) then
  !m_c=m_start
  call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_c,beta_1,m_c)
  call Calc_FWD(Lambda_c,d,m_c,d_Pred,res_c,eAll,F,mNorm_c,rms_c)
                   write(iterChar,'(i3.3)') fwd_calls
  				   file_name_suffix='Hybrid_'//iterChar
		           call write_output_files(DS_iter,d_Pred,m_c,file_name_suffix) 
  fwd_calls=fwd_calls +1
end if  
!right Lambda 
if (go_right) then
  !m_r=m_start
  call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_r,beta_1,m_r)
  call Calc_FWD(Lambda_r,d,m_r,d_Pred,res_r,eAll,F,mNorm_r,rms_r)
                   write(iterChar,'(i3.3)') fwd_calls
  				   file_name_suffix='Hybrid_'//iterChar
		           call write_output_files(DS_iter,d_Pred,m_r,file_name_suffix) 
  fwd_calls=fwd_calls +1
end if  
!left Lambda 
if (go_left) then
  !m_l=m_start
  call get_model(k_step,T_matrix,q,m0,mhat,b,Lambda_l,beta_1,m_l)
  call Calc_FWD(Lambda_l,d,m_l,d_Pred,res_l,eAll,F,mNorm_l,rms_l)
                   write(iterChar,'(i3.3)') fwd_calls
  				     file_name_suffix='Hybrid_'//iterChar
		           call write_output_files(DS_iter,d_Pred,m_l,file_name_suffix)     
  fwd_calls=fwd_calls +1
end if  


    write(iterChar1,'(i3.3)') DS_iter
    if (DS_iter .eq. 1 .and. k_step .eq. min_k_steps) then
       open(10,file='L_serach.dat',STATUS = 'unknown')
    else
       open(10,file='L_serach.dat',status='unknown',position='append')
   end if
       
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
  
if (rms_c .gt. rms .and. rms_r .gt. rms .and. rms_l .gt. rms) then
			 write(10,*)'The RMS in all directions is higher than the start one'	
			 write(10,*)'Using the start solution--> Recompute' 
			 write(10,*) 'Increase the size of the projected system by adding additional K step'
			 if (k_step .lt. 10) then
				     !call Calc_FWD(lambda,d,m_start,d_Pred,res,eAll,F,mNorm,rms)
				     eAll=eAll_temp
				     !fwd_calls=fwd_calls +1
		             k_step=k_step+1
		             update_proj_sys=.true.
		             close(10)
		             lambda=dlog10(lambda_c)
		             m_temp=m_c
		             goto 1
		     end if        
			 
elseif (rms_c .gt. rms .and. rms_r .gt. rms .and. rms_l .gt. rms .and. rms_l .gt. rms_c .and. rms_r .gt. rms_c) then
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
				elseif (rms_r .lt. rms_c ) then
		          		  m=m_r
						  lambda=dlog10(lambda_r)
		        else
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
				  close(10)
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
		  close(10)
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
		  if (k_step .lt. 10) then
		    if (rms_c .gt. rms_temp )then
		          m=m_temp
				  call linComb_modelParam(ONE,m,MinusONE,m0,mhat)
				  call Calc_FWD(lambda,d,m,d_Pred,res,eAll,F,mNorm,rms)
				  fwd_calls=fwd_calls +1
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
			  call linComb_modelParam(ONE,m,MinusONE,m0,mhat)
			  call Calc_FWD(lambda,d,m,d_Pred,res,eAll,F,mNorm,rms)
			  fwd_calls=fwd_calls +1
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
  
  if (i_search .gt. 30 ) then
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
Subroutine bidiag_JJT (b,m,d,k_step,DS_iter,lambda,T_matrix,q,beta_1,beta_kstep_plus_1,u)
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
  real(kind=prec)                           :: alpha,sigma11,alpha_beta,lambda_sub
  Integer                                   :: i,iDt,ii,jj,INFO,kk
  real(kind=prec)	,pointer, dimension(:)    ::  AP,b_sub,x_sub
  real(kind=prec),pointer,dimension(:,:) ::T_matrix_temp
  
allocate(T_matrix_temp(k_step,k_step))  
allocate(AP((k_step)*((k_step)+1)/2),b_sub(k_step),x_sub(k_step))
  
allocate(beta(k_step+1)) 



Au=b
w=b
u_temp=b
   call zero(Au)
   call zero(w)
   call zero(u_temp)
   

do ii=1,21
  u(ii)= b
  call zero(u(ii))
end do

beta=R_zero
u(1)=b
beta(1)=sqrt(dotProd(u(1),u(1)))
call scMult(1.0/beta(1),u(1),u_temp)
u(1)=u_temp
beta_1=beta(1)
 
 do ii=1,k_step

		 call MultA_JJT_DS(u(ii),m,d,Au,q(ii))
		 
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
            

        
  
        
        
        
        
      write(6,*) ii, 'Beta= ', beta(ii+1)/beta(1), beta(ii+1)
      write(6,*) ii, 'Orth= ', dotProd(u(1),u(ii))
 end do
 
beta_kstep_plus_1=beta(k_step+1)

  10 continue 
 do ii=1,k_step
 write(6,'(a6,5f10.5)') 'Orth= ', (dotProd(u(ii),u(jj)),jj=1,k_step)
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
subroutine get_model (k_step,T_matrix,q,m0,mhat,b,lambda,beta_1,m)

Integer, intent(in)                          :: k_step
real(kind=prec)	,intent(in), dimension(:,:)  :: T_matrix
type(modelParam_t),  intent(in)              :: m0
type(modelParam_t),  intent(inout)           :: mhat
type (modelParam_t),intent(in), dimension(:) :: q
!type (dataVectorMTX_t),intent(in), dimension(:) :: u
type (dataVectorMTX_t),intent(in)               :: b
real(kind=prec)	,intent(in)                  :: lambda,beta_1
type (modelParam_t),intent(inout)              :: m

!Local 

real(kind=prec)	,pointer, dimension(:,:)  ::  T_matrix_temp,L_matrix, U_matrix
real(kind=prec)	,pointer, dimension(:)    ::  b_sub,y_sub,x_sub
type (modelParam_t)                       :: q_temp
Integer                                   :: ii,jj,INFO
real(kind=prec)                           :: lambda_LU,mu_LU,sum
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
             call DPPTRS( 'U', k_step, 1, AP, b_sub, k_step, INFO )
              x_sub=b_sub
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
         !write(6,*)dotProd(u(ii),b)
     end do
     b_sub(1)=beta_1
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
! Model update 
		 do ii=1, k_step
			 call scMult_modelParam (x_sub(ii),q(ii),q_temp)
			 call linComb_modelParam(ONE,mhat,ONE,q_temp,mhat)
		 end do
 
		 !mhat= multBy_Cm(mhat)
! add to m0	 
         
	     call linComb_modelParam(ONE,m,ONE,mHat,m)

end subroutine get_model
!****************************************************************************************  

   

subroutine CG_MS(b,x,m,d,lambda,CGiter)


  type (modelParam_t), intent(in)	 	::b
  type (modelParam_t), intent(inout) 	::x
  type(modelParam_t),  intent(in)       ::m
  type (dataVectorMTX_t), intent(in)	    ::d
  real(kind=prec),     intent(in)       ::lambda
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

ii = 1
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

CGiter%niter = ii

! deallocate the help vectors
    call deall_modelParam(r)
    call deall_modelParam(p)
    call deall_modelParam(Ap)
    
end subroutine CG_MS
!****************************************************************************************
 
subroutine CG_DS(b,x,m,d,lambda,CGiter)


  type (dataVectorMTX_t), intent(in)	 	::b
  type (dataVectorMTX_t), intent(inout) 	::x
  type(modelParam_t),  intent(in)       ::m
  type (dataVectorMTX_t), intent(in)	    ::d
  real(kind=prec),     intent(in)       ::lambda
  type(iterControl_t), intent(inout)	:: CGiter
  character(3)         					::iterChar
  
  !Local
    type (dataVectorMTX_t)              	:: r,p,Ap
    real(kind=prec)					 	::alpha,beta,r_norm_pre,r_norm,b_norm,error,delta_new,delta_zero,delta_old
    integer                          	::cg_iter,i,j,k,ii,iDt
    
     

r=b
p=r
Ap=d
b_norm=dotProd(b,b)
call zero_dataVectorMTX(x)
r_norm=dotProd(r,r)


         !call linComb(R_ZERO,d,ONE,r,r) 
         !call linComb(R_ZERO,d,ONE,p,p)
         !call linComb(R_ZERO,d,ONE,x,x) 
         !call linComb(R_ZERO,d,ONE,Ap,Ap)  
                
ii = 1
CGiter%rerr(ii) = r_norm/b_norm
 write(10,*) 'CG-error',ii, r_norm/b_norm
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
       write(6,*) 'CG-error',ii, r_norm/b_norm
      write(10,*) 'CG-error',ii, r_norm/b_norm

       !write(6,*) 'Beta_CG',ii, sqrt(beta)/alpha
  end do loop

CGiter%niter = ii

! deallocate the help vectors
    call deall_dataVectorMTX(r)
    call deall_dataVectorMTX(p)
    call deall_dataVectorMTX(Ap)
    
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
subroutine MultA_DS(p,m,d,lambda,Ap)
   type(dataVectorMTX_t), intent(in)          ::p
   type(dataVectorMTX_t), intent(out)         ::Ap
   type(modelParam_t), intent(in)          ::m
   type (dataVectorMTX_t), intent(in)	       ::d
   real(kind=prec), intent(in)             ::lambda
!Local parameters
   type(modelParam_t)                      ::JTp,CmJTp
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

            call scMult_dataVectorMTX(lambda,p,lambdaP)
            
!Normalize: Cd^(-1/2)*Ap
            call normalize_with_dataVecMTX(Ap,d,1)

!Add Cd^(-1/2)*Ap*Cd^(-1/2) to lambda*p
         do i=1,lambdaP%nTx
           do iDt=1,lambdaP%d(i)%nDt
             lambdaP%d(i)%data(iDt)%errorBar= .false.
            end do
         end do
          
           call linComb_dataVectorMTX(ONE,Ap,ONE,lambdaP,Ap)      
        
        
!Deallocate help vectors
    call deall_modelParam(JTp)
    call deall_modelParam(CmJTp)
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
            call Master_job_JmultT(m,p_temp,JTp,eAll)

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
   subroutine Calc_FWD(lambda,d,m,d_Pred,res,eAll,F,mNorm,rms)

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

   !  local variables
   type(dataVectorMTX_t)    :: Nres
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
end module DCG
