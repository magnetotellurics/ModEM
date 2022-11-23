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
     type(DampingFactorANI) :: lambdaANI
     ! exit if lambda < lambdaTol approx. 1e-4
     real (kind=prec)   :: lambdaTol
     ! set lambda_i = lambda_{i-1}/k when the inversion stalls
     real (kind=prec)   :: k
     real (kind=prec)   :: p
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
  real (kind=prec), private, save :: lambda_avrg

Contains

!**********************************************************************
   subroutine set_NLCGiterControl(iterControl)
   !  at present this is a simple routine to setup for
   !  iteration control; allowing parameters to be read from
   !   a file would be desireable.

   type(NLCGiterControl_t), intent(inout)	:: iterControl

     ! maximum number of iterations in one call to iterative solver
     iterControl%maxIter = 600
     ! convergence criteria: return from solver if rms < rmsTol
     iterControl%rmsTol  = 1.05
     ! inversion stalls when abs(rms - rmsPrev) < fdiffTol (2e-3 works well)
     iterControl%fdiffTol = 2.0e-3
     ! initial value of lambda (will not override the NLCG input argument)
     iterControl%lambda = 10.
     iterControl%lambdaANI%xy = 0.
     iterControl%lambdaANI%yz = 0.
     iterControl%lambdaANI%xz = 0.
     ! exit if lambda < lambdaTol approx. 1e-4
     iterControl%lambdaTol = 1.0e-8
     ! set lambda_i = lambda_{i-1}/k when the inversion stalls
     iterControl%k = 10.
     iterControl%p = 1.
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
    iterControl%fname = adjustl(iterControl%fname)
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
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%lambdaANI%xy
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%lambdaANI%xy
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%lambdaANI%yz
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%lambdaANI%yz
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%lambdaANI%xz
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%lambdaANI%xz
    !end if
    !read (ioInvCtrl,'(a36,g15.7)') string,iterControl%p
    !if (output_level > 2) then
    !   write (*,'(a36,g15.7)') string,iterControl%p
       write(*,*)
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
   
   !!
   real(kind=prec) :: aniNorm

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

   !!!
   call cal_anisCsNorm(m,aniNorm,iterControl%lambdaANI)
   
   lambda_avrg = (iterControl%lambdaANI%xy+iterControl%lambdaANI%yz+iterControl%lambdaANI%xz)/3.d0
   
   ! penalty functional = sum of squares + scaled model norm
   F = SS/Ndata + (lambda/(1+lambda_avrg) * mNorm/Nmodel) + lambda/(1+lambda_avrg) * aniNorm/Nmodel

   ! scale mNorm for output
   mNorm = mNorm/Nmodel
   aniNorm = aniNorm/Nmodel
   write(1212,*) lambda/(1+lambda_avrg) * mNorm,lambda/(1+lambda_avrg) * aniNorm
   
   ! if required, compute the Root Mean Squared misfit
   if (present(RMS)) then
   	RMS = sqrt(SS/Ndata)
   end if

   call deall_dataVectorMTX(res)
   call deall_dataVectorMTX(Nres)
   call deall_modelParam(m)
   call deall_modelParam(JTd)

   end subroutine func

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
   type(dataVectorMTX_t), intent(in)              :: d
   type(modelParam_t), intent(in)           :: m0
   type(modelParam_t), intent(in)           :: mHat
   type(modelParam_t), intent(inout)          :: grad
   type(dataVectorMTX_t), intent(inout)              :: dHat
   type(solnVectorMTX_t), intent(inout)            :: eAll

   !  local variables
   real(kind=prec)       :: Ndata,Nmodel
   type(dataVectorMTX_t)    :: res
   type(modelParam_t) :: m,JTd,CmJTd

   ! integer :: j, Ny, NzEarth

   ! compute the smoothed model parameter vector
   call CmSqrtMult(mHat,m)

   ! overwriting the input with output
   call linComb(ONE,m,ONE,m0,m)

   ! initialize res
   res = d

   ! compute residual: res = (d-dHat)/Ndata
   call linComb(ONE,d,MinusONE,dHat,res)

   ! multiply by J^T
   call CdInvMult(res)

#ifdef MPI
        call Master_job_JmultT(m,res,JTd,eAll)
#else
        call JmultT(m,res,JTd,eAll)
#endif
   
   call AnalyzeAniModel(m0,JTd) ! analyze which anisotropic (including isotropic) model we want to obtain through inversion 

   !call write_modelParam_ws(JTd,'JTd_MT3DANI.txt') !debug
   !pause !debug   
   
   call CmSqrtMult(JTd,CmJTd)

   ! initialize grad
   grad = m

   ! compute the number of data and model parameters for scaling
   Ndata = countData(res)
   Nmodel = countModelParam(mHat)

   lambda_avrg = (iterControl%lambdaANI%xy+iterControl%lambdaANI%yz+iterControl%lambdaANI%xz)/3.d0
   
   ! multiply by 2 (to be consistent with the formula)
   ! and add the gradient of the model norm
   call linComb(MinusTWO/Ndata,CmJTd,TWO*lambda/(1+lambda_avrg)/Nmodel,mHat,grad)

   call cal_anisCsVec(m,JTd,iterControl%lambdaANI) ! Kong
   call CmSqrtMult(JTd,CmJTd) ! Kong
   call linComb(ONE,grad,TWO*lambda/(1+lambda_avrg)/Nmodel,CmJTd,grad)
   
   call deall_dataVectorMTX(res)
   call deall_modelParam(m)
   call deall_modelParam(JTd)
   call deall_modelParam(CmJTd)
   !call deall(eAll)
   end subroutine gradient

!**********************************************************************
   subroutine update_damping_parameter(lambda,mHat,F,grad,m0)

   real(kind=prec), intent(inout)  :: lambda
   type(modelParam_t), intent(in)              :: mHat,m0
   real(kind=prec), intent(inout)  :: F
   type(modelParam_t), intent(inout)             :: grad

   real(kind=prec) :: SS, mNorm, Nmodel, aniNorm
	 type(modelParam_t)          :: dSS, m, JTd

   ! compute the model norm
   mNorm = dotProd(mHat,mHat)
   Nmodel = countModelParam(mHat)
   
   ! compute the smoothed model parameter vector
   call CmSqrtMult(mHat,m)
   ! overwriting the input with output
   call linComb(ONE,m,ONE,m0,m)
   call cal_anisCsNorm(m,aniNorm,iterControl%lambdaANI)
   JTd = m
   call cal_anisCsVec(m,JTd,iterControl%lambdaANI) ! Kong
   call CmSqrtMult(JTd,m) ! Kong
   
   lambda_avrg = (iterControl%lambdaANI%xy+iterControl%lambdaANI%yz+iterControl%lambdaANI%xz)/3.d0
   
   ! (scaled) sum of squares = penalty functional - scaled model norm
   SS = F - (lambda/(1+lambda_avrg) * mNorm/Nmodel) - lambda/(1+lambda_avrg) * aniNorm/Nmodel

   ! initialize
   dSS = mHat

   ! subtract the model norm derivative from the gradient of the penalty functional
   call linComb(ONE,grad,MinusTWO*lambda/(1+lambda_avrg)/Nmodel,mHat,dSS)
   call linComb(ONE,dSS,MinusTWO*lambda/(1+lambda_avrg)/Nmodel,m,dSS)
   
   ! update the damping parameter lambda
   lambda = lambda/iterControl%k
   iterControl%lambdaANI%xy = iterControl%lambdaANI%xy/iterControl%p
   iterControl%lambdaANI%yz = iterControl%lambdaANI%yz/iterControl%p
   iterControl%lambdaANI%xz = iterControl%lambdaANI%xz/iterControl%p
 
   lambda_avrg = (iterControl%lambdaANI%xy+iterControl%lambdaANI%yz+iterControl%lambdaANI%xz)/3.d0
   
   ! penalty functional = (scaled) sum of squares + scaled model norm
   F = SS + (lambda/(1+lambda_avrg) * mNorm/Nmodel) + lambda/(1+lambda_avrg) * aniNorm/Nmodel/iterControl%p

	 ! add the model norm derivative to the gradient of the penalty functional
   call linComb(ONE,dSS,TWO*lambda/(1+lambda_avrg)/Nmodel,mHat,grad)
   call linComb(ONE,grad,TWO*lambda/(1+lambda_avrg)/Nmodel/iterControl%p,m,grad)

   call deall_modelParam(dSS)
   call deall_modelParam(m)
   call deall_modelParam(JTd)
   
   end subroutine update_damping_parameter

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

   if (present(fname)) then
      call read_NLCGiterControl(iterControl,fname,ok)
      if (ok) then
         lambda = iterControl%lambda
      end if
   else
      call set_NLCGiterControl(iterControl)
   end if

   ! initialize the output to log file
   logFile = trim(iterControl%fname)//'_NLCG.log'
   open (unit=ioLog,file=logFile,status='unknown',position='append',iostat=ios)
   open(1212,file=trim(iterControl%fname)//'_NLCG_Norms.dat')
   write(1212,*) ' model norm          ','    anisotropy norm'

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

   ! output (smoothed) initial model and responses for later reference
   call CmSqrtMult(mHat,m_minus_m0)
   call linComb(ONE,m_minus_m0,ONE,m0,m)
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
   if (output_level > 3) then
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
      		call update_damping_parameter(lambda,mHat,value,grad,m0)
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
   close(1212,iostat=ios)

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
   type(dataVectorMTX_t), intent(in)		       :: d
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
   type(dataVectorMTX_t), intent(in)		       :: d
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
        if(isnan(f)) then
          f = f_1
          starting_guess = .true.
          exit
        endif
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
