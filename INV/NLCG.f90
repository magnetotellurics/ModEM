module NLCG

! inherits SensComp, DataIO and all modules they use. Also Main_MPI and Sub_MPI

use invcore


implicit none

public  :: NLCGsolver

! iteration control for the NLCG solver is initialized once
! and saved in the module to be used by most subroutines

  type  :: NLCGiterControl_t
     ! maximum number of iterations in one call to iterative solver
     integer            :: maxIter
     ! convergence criteria: return from solver if rms < rmsTol
     real (kind=prec)   :: rmsTol
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
     ! the factor that ensures culvature condition in the line search
     real (kind=prec)   :: c2
     ! restart CG every nCGmax iterations to ensure conjugacy
     integer            :: nCGmax
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

    type(NLCGiterControl_t), intent(inout)  :: iterControl
    character(*), intent(in)                :: rFile
    logical, intent(out), optional          :: fileExists
    integer                                 :: ios
    logical                                 :: exists
    character(80)                           :: string

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
       write (*,*)
    end if

    close(ioInvCtrl)

   end subroutine read_NLCGiterControl

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
   type(dataVectorMTX_t), intent(inout)         :: d
   !  lambda is regularization parameter
   real(kind=prec), intent(inout)               :: lambda
   !  m0 is prior model parameter
   type(modelParam_t), intent(in)               :: m0
   !  m is solution parameter ... on input m contains starting guess
   type(modelParam_t), intent(inout)            :: m
   !  fname is a string that specifies the control file
   character(*), intent(in), optional           :: fname
   !  initial step size in the line search direction in model units
   real(kind=prec)                              :: startdm
   !  flavor is a string that specifies the algorithm to use
   character(80)                                :: flavor = 'Cubic'

   !  local variables
   type(dataVectorMTX_t)                        :: dHat, res
   type(modelParam_t)                           :: mHat, m_minus_m0
   type(modelParam_t)                           :: grad, g, h, gPrev
   !type(NLCGiterControl_t)			            :: iterControl
   real(kind=prec)                              :: value, valuePrev, rms
   real(kind=prec)                              :: rmsPrev, alpha, beta
   real(kind=prec)                              :: gnorm, mNorm, Nmodel
   real(kind=prec)                              :: grad_dot_h, g_dot_g
   real(kind=prec)                              :: g_dot_gPrev,g_dot_h
   real(kind=prec)                              :: gPrev_dot_gPrev 
   real(kind=prec)                              :: h_dot_g, h_dot_gPrev
   integer                                      :: iter, nCG, nLS, nfunc, ios
   logical                                      :: ok
   character(3)                                 :: iterChar
   character(100)                               :: mFile, mHatFile, gradFile
   character(100)                               :: dataFile, resFile, logFile
   type(solnVectorMTX_t)                        :: eAll

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
   write(*,'(a42,es12.5)') '    GRAD: initial norm of the gradient is',gnorm
   write(ioLog,'(a42,es12.5)') '     GRAD: initial norm of the gradient is',gnorm
   if (gnorm < TOL6) then
      call errStop('Problem with your gradient computations: first gradient is zero')
   else !if (alpha * gnorm > startdm) then
      alpha = startdm / gnorm
      write(*,'(a39,es12.5)') 'The initial value of alpha updated to ',alpha
      write(ioLog,'(a39,es12.5)') 'The initial value of alpha updated to ',alpha
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
	  case ('Wolfe')
              call lineSearchWolfe(lambda,d,m0,h,alpha,mHat,value,grad,rms,nLS,dHat,eAll)
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
            write(*,'(a34,es12.5)') 'The norm of the last gradient is ',gnorm
            write(ioLog,'(a34,es12.5)') 'The norm of the last gradient is ',gnorm
            !alpha = min(iterControl%alpha_1,startdm/gnorm)
      		alpha = min(ONE,startdm)/gnorm
      		write(*,'(a48,es12.5)') 'The value of line search step alpha updated to ',alpha
            write(ioLog,'(a48,es12.5)') 'The value of line search step alpha updated to ',alpha
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
      !if ((beta.lt.R_ZERO).or.(g_dot_g + beta*g_dot_h .le. R_ZERO)&
      !    .and.(nCG .ge. iterControl%nCGmax)) then  !PR+
      if ((g_dot_g + beta*g_dot_h .le. R_ZERO)&
          .and.(nCG .ge. iterControl%nCGmax)) then  !PR
          ! restart
          write(*,'(a45)') 'Restarting NLCG to restore orthogonality'
          write(ioLog,'(a45)') 'Restarting NLCG to restore orthogonality'
          nCG = 0
          beta = R_ZERO
      else
          nCG = nCG + 1
      end if
      call linComb(ONE,g,beta,h,h)
   end do

   ! multiply by C^{1/2} and add m_0
   call CmSqrtMult(mHat,m_minus_m0)
   call linComb(ONE,m_minus_m0,ONE,m0,m)
   d = dHat
   write(*,'(a25,i5,a25,i5)') 'NLCG iterations:',iter,' function evaluations:',nfunc
   write(ioLog,'(a25,i5,a25,i5)') 'NLCG iterations:',iter,' function evaluations:',nfunc
   close(ioLog,iostat=ios)

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
    ! The step size alpha should not be adjusted manually at all!
    ! Even when it is too small or too close to the previous try,
    ! adjusting it won't result in an improvement (it's better to exit
    ! the line search in that case, if anything)...
    ! if ((alpha_i - alpha < eps).or.(alpha < k*alpha_i)) then
    !     alpha = alpha_i/TWO ! reset alpha to ensure progress
    ! end if
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

  !**********************************************************************
  subroutine lineSearchWolfe(lambda,d,m0,h,alpha,mHat,f,grad, &
   & rms,niter,dHat,eAll)
 
   ! Contributed by Hao Dong
   
   ! Note: inexact line searches ultimately fit into two catalogs - 
   ! i.e. Armijo–Goldstein (back-tracking) and Wolfe conditions
   ! the latter also contains a few different criterias: (1) "Armijo 
   ! rule", and (2) "curvature condition", the latter requires an extra 
   ! gradient calculation at x_k + alpha_k*d_k
   ! (1) and modified (2) will form the "strong Wolfe" condition which is 
   ! the line search "convergence criteria" here

   ! Original Line search code here (by Anna) was based on the Numerical 
   ! Recipes and on the text by Michael Ferris, Chapter 3, p 59. 
   ! which is about the sufficient decrease (Armijo) condition 
   ! (ignoring the curvature condition).
   ! first interpolate using a quadratic approximation; if the
   ! solution does not satisfy the condition, then backtrack using
   ! cubic interpolation. 
   !
   ! The initial step size is set outside of this routine (in the LBFGS)
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
   ! the solution neet to satisfy the Armijo's rule only:
   !     f(alpha) < f(0) + c alpha f'(0).
   ! where c is a positive small value close to zero (0 < c <= 0.5)
   !
   ! This good part is this strategy only requires one gradient
   ! evaluation and is very efficient when computing gradients is
   ! expensive (as with EM). 
   !
   !
   ! but hey, penalty function evaluation is equally (if not more) expensive
   ! , and we need to calculate the gradient function for the next search 
   ! direction anyway - sooner or later, so here is my idea: 
   !
   ! when we have evaluated the penalty function for the initial guess (f_1)
   ! and quadratic interpolation (f)
   ! 1) if f < f_1, we just continue to calculate the gradient at the
   !    quadratic interpolation point, and use it to test the Wolfe condition
   !    if it satisfies, then the calculation is the same as Anna's scheme 
   !    2 penalty funtion evaluations and 1 gradient = 2 FWDs and 1 TRN
   ! 2) if f_1 < f, we do the same as 1), but only calculate grad at initial
   !    guess, still the calculation will be the same if everything works
   !    i.e. 2 FWDs and 1 TRN
   ! 3) if none of the 1) and 2) satisfied, the scheme falls back to the 
   !    More-Thunente scheme, the reason I didn't use a backtracking scheme
   !    like in Anna's Cubic subroutine, is the fact that a step length
   !    greater than alpha_1 may be needed to satisfy Wolfe's condition 
   !    (so "back"tracking may not be enough)
    
   ! in practise, this would almost never go to step 3 in the first 100
   ! iterations 
   !
   ! following Anna's idea, the major intention here is to use a cheap line
   ! search scheme (with only 3 forward-like-calculations) to quickly skip 
   ! to a small overall penalty function level and only to start bracketing
   ! when the quadratic interpolation doesn't work, in which case cubic 
   ! probably won't work either...
   ! 
   ! the real motivation, however, is try to implement Strong Wolfe condition 
   ! to prepare for L-BFGS method -
   ! as the standard Armijio backtracking method, which ignores the curvature 
   ! condition, cannot guarantee the stable converge of quasi-Newton-ish method

   real(kind=prec), intent(in)               :: lambda ! lagrange multiplier
   type(dataVectorMTX_t), intent(in)         :: d
   type(modelParam_t), intent(in)            :: m0 ! current model
   type(modelParam_t), intent(in)            :: h  ! search direction
   real(kind=prec), intent(inout)            :: alpha ! step size
   type(modelParam_t), intent(inout)         :: mHat 
   real(kind=prec), intent(inout)            :: f  ! penalty function
   type(modelParam_t), intent(inout)         :: grad ! function gradient
   real(kind=prec), intent(out)              :: rms
   integer, intent(out)                      :: niter
   type(dataVectorMTX_t), intent(out)        :: dHat
   type(solnVectorMTX_t), intent(inout)      :: eAll
   ! local variables
   real(kind=prec)                 :: alpha_1,alpha_i,alpha_j,mNorm
   real(kind=prec)                 :: alpha_l,alpha_r !left/right bound 
   logical                         :: starting_guess
   integer                         :: ibracket
   real(kind=prec)                 :: eps,k,c,c2,a,b,q1,q2,q3
   real(kind=prec)                 :: g_0,g_1,f_0,f_1,f_i,f_j,rms_1,mNorm_1
   type(modelParam_t)              :: mHat_0,mHat_1
   type(dataVectorMTX_t)           :: dHat_1
   type(solnVectorMTX_t)           :: eAll_1
   character(100)                  :: logFile

   ! parameters 
   c = iterControl%c ! ensures sufficient decrease
   c2 = iterControl%c2 ! ensures culvature condition
   !k = iterControl%alpha_k
   !eps = iterControl%alpha_tol
   logFile = trim(iterControl%fname)//'_NLCG.log'

   ! initialize the line search
   niter = 0
   mHat_0 = mHat
   f_0 = f
   starting_guess = .false. 

   ! rescale the search direction
   ! h_dot_h = dotProd(h,h)
   ! h = scMult_modelParam(ONE/sqrt(h_dot_h),h)

   ! g_0 is the directional derivative of our line search function 
   ! f'(0) = (df/dm).dot.h
   g_0 = dotProd(grad,h)

   ! setup the lower and upper boundary of alpha
   alpha_l = R_ZERO
   alpha_r = (f_0-(f_0*0.1))/(-g_0)/c2

   ! alpha_1 is the initial step size, which is set in LBFGS
   alpha_1 = alpha
   ! compute the trial mHat, f, dHat, eAll, rms
   mHat_1 = mHat_0
   ! mHat_1 = mHat_0 + dir*step
   call linComb(ONE,mHat_0,alpha_1,h,mHat_1)
   call func(lambda,d,m0,mHat_1,f_1,mNorm_1,dHat_1,eAll_1,rms_1)
   call printf('STARTLS',lambda,alpha_1,f_1,mNorm_1,rms_1)
   call printf('STARTLS',lambda,alpha_1,f_1,mNorm_1,rms_1,logFile)
   niter = niter + 1

   if (f_1 - f_0 >= LARGE_REAL) then
   ! oops, we are pushing too far away
       write(ioLog,'(a40)') 'Try a smaller starting value of alpha'
       write(*,'(a40)') 'Try a smaller starting value of alpha'
       write(ioLog,'(a10)') 'Exiting...'
       write(*,'(a10)') 'Exiting...'
       stop
   end if


   ! quadratic interpolation of a parabola
   a = TWO*(f_1 - f_0 - g_0*alpha_1)
   b = g_0*(alpha_1**2)
   alpha = - b/a
   ! if the curvature is -ve, there is no minimum; take the initial guess
   if (a < 0) then ! upside down parabola 
       starting_guess = .true.
       alpha = alpha_1
       dHat = dHat_1
       eAll = eAll_1
       mHat = mHat_1
       rms = rms_1
       f = f_1
       ! compute the gradient and exit
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
   call linComb(ONE,mHat_0,alpha,h,mHat)
   call func(lambda,d,m0,mHat,f,mNorm,dHat,eAll,rms)
   call printf('QUADLS',lambda,alpha,f,mNorm,rms)
   call printf('QUADLS',lambda,alpha,f,mNorm,rms,logFile)
   niter = niter + 1
   ! check whether the solution satisfies the sufficient decrease condition
   ! Strong Wolfe's condition needs the gradient 
   ! well, we are going to calculate it anyway - so why don't we do it now?
   if (f <= f_1) then 
       call gradient(lambda,d,m0,mHat,grad,dHat,eAll)
       g_1 = dotProd(grad, h)
       write(*,'(a29,es12.5)',advance='no') '    GRAD: computed, with g0=',g_0
       write(*,'(a4,es12.5)') ' g1=',g_1
       write(ioLog,'(a29,es12.5)',advance='no') '    GRAD: computed, with g0=',g_0
       write(ioLog,'(a4,es12.5)') ' g1=',g_1
       ! Note we test the Strong Wolfe's condition: Armijio's rule and
       ! curvature condition 
       if ((f <= f_0 + c * alpha * g_0).and.(abs(g_1) <= c2*abs(g_0))) then
           write(*,'(a53)') 'Strong Wolfe Condition satisfied, exiting line search'
           write(ioLog,'(a53)') 'Strong Wolfe Condition satisfied, exiting line search'
           call deall_dataVectorMTX(dHat_1)
           call deall_modelParam(mHat_0)
           call deall_modelParam(mHat_1)
           call deall_solnVectorMTX(eAll_1)
           return
       else 
           !ooops, we missed the Strong Wolfe's condtion (for one reason or 
           !the other 
           if ((alpha_r-alpha_l)*g_1<0) then
               ! update the left boundary for alpha
               alpha_l = alpha
           else
               ! update the right boundary for alpha
               alpha_r = alpha
           endif
       endif
   else 
       if (f_1 < f_0) then! is the initial making any progress?
           ! Test if the initial guess is good for Strong Wolfe condition 
           call gradient(lambda,d,m0,mHat_1,grad,dHat_1,eAll_1)
           g_1 = dotProd(grad, h) 
           write(*,'(a29,es12.5)',advance='no') '    GRAD: computed, with g0=',g_0
           write(*,'(a4,es12.5)') ' g1=',g_1
           write(ioLog,'(a29,es12.5)',advance='no') '    GRAD: computed, with g0=',g_0
           write(ioLog,'(a4,es12.5)') ' g1=',g_1
           if ((f_1 <= f_0 + c * alpha_1 * g_0).and.(abs(g_1) <= c2*abs(g_0))) then 
               write(*,'(a53)') 'Strong Wolfe Condition satisfied, exiting line search'
               write(ioLog,'(a53)') 'Strong Wolfe Condition satisfied, exiting line search'
               starting_guess = .true.
               alpha = alpha_1
               dHat = dHat_1
               eAll = eAll_1
               mHat = mHat_1
               rms = rms_1
               f= f_1
               call deall_dataVectorMTX(dHat_1)
               call deall_modelParam(mHat_0)
               call deall_modelParam(mHat_1)
               call deall_solnVectorMTX(eAll_1)
               return
           endif
       else
           !ooops, we missed the Strong Wolfe's condtion (for one reason or 
           !the other 
           if ((alpha_r-alpha_l)*g_1<0) then
               ! update the left boundary for alpha
               alpha_l = alpha_1
           else
               ! update the right boundary for alpha
               alpha_r = alpha_1
           endif
       endif
   endif
   ! someone has to spread the bad news
   write(*,'(a50)') '!======Strong Wolfe Condition NOT satisfied!======'
   write(ioLog,'(a50)') '!======Strong Wolfe Condition NOT satisfied!======'

   if (f > f_0) then
   ! this should not happen, but in practice it is possible to end up with
   ! a function increase at this point (e.g. in the current global code).
   ! Not that rare actually - even in normal MT.
   !
   ! Most likely, this is due to an inaccuracy in the gradient computations.
   ! In this case, we avoid an infinite loop by exiting line search.
   ! It is also possible that both f_1 and f are worse than the starting value!
   ! Then, take whichever is smaller. Ideally, want to decrease the tolerance
   ! for gradient computations if this happens.
   ! 
   ! it is unlikely, but also possible that we used a far-too-large alpha value
       if (f_1 < f) then ! pick the less-bad solution
           starting_guess = .true.
       endif
       write(*,'(a75)') 'Unable to fit a quadratic due to bad gradient estimate, exiting line search'
       write(ioLog,'(a75)') 'Unable to fit a quadratic due to bad gradient estimate, exiting line search'
   else ! /inhales
       ! sectioning and fit cubic (initialize)
       write(*,'(a50)') '!========Now sectioning with brute force========'
       write(ioLog,'(a50)') '!========Now sectioning with brute force========'
       write(6,*) 'alpha_l=',alpha_l,'alpha_r=',alpha_r
       write(ioLog,*) 'alpha_l=',alpha_l,'alpha_r=',alpha_r
       alpha_i = alpha_1 !START
       f_i = f_1
       alpha_j = alpha   !QUADRATIC
       f_j = f
       ibracket = 0;
       fit_cubic: do
       ! compute the minimizer
           q1 = f_i - f_0 - g_0 * alpha_i
           q2 = f_j - f_0 - g_0 * alpha_j
           q3 = alpha_i**2 * alpha_j**2 * (alpha_j - alpha_i)
           a = (alpha_i**2 * q2 - alpha_j**2 * q1)/q3
           b = (alpha_j**3 * q1 - alpha_i**3 * q2)/q3
           if ((b*b-3*a*g_0)<0) then ! failed to fit cubic
               write(*,'(a40)') 'SQRT of negative value in CUBIC INTERP!'
               write(ioLog,'(a40)') 'SQRT of negative value in CUBIC INTERP!'
               write(*,'(a35)') 'using default value to bracket...'
               write(ioLog,'(a35)') 'using default value to bracket...'
               alpha = sqrt(alpha_l*alpha_r)
           else
               if (b<=R_ZERO) then ! fit cubic 
                   alpha = (- b + sqrt(b*b - 3.0*a*g_0))/(3.0*a)
               else
                   alpha = -g_0/(b+sqrt(b*b - 3.0*a*g_0))
               endif
           endif
           ! if (alpha < 0.1*alpha_i) then !avoid too small a step 
           !     alpha = 0.1 * alpha_i
           !     write(*,'(a30)') 'searching step too small...'
           !     write(ioLog,'(a30)') 'searching step too small...'
           !     write(*,'(a21)') 'use default value...'
           !     write(ioLog,'(a21)') 'use default value...'
           ! endif
           ! if alpha is too close or too much smaller than its predecessor
           ! if ((alpha_j - alpha < eps).or.(alpha < k*alpha_j)) then
           !     alpha = alpha_j/TWO ! reset alpha to ensure progress
           ! end if
           ! compute the penalty functional
           call linComb(ONE,mHat_0,alpha,h,mHat)
           call func(lambda,d,m0,mHat,f,mNorm,dHat,eAll,rms)
           call printf('CUBICLS',lambda,alpha,f,mNorm,rms)
           call printf('CUBICLS',lambda,alpha,f,mNorm,rms,logFile)
           niter = niter + 1
           ibracket = ibracket + 1
           call gradient(lambda,d,m0,mHat,grad,dHat,eAll)
           g_1 = dotProd(grad, h) 
           write(*,'(a29,es12.5)',advance='no') '    GRAD: computed, with g0=',g_0
           write(*,'(a4,es12.5)') ' g1=',g_1
           write(ioLog,'(a29,es12.5)',advance='no') '    GRAD: computed, with g0=',g_0
           write(ioLog,'(a4,es12.5)') ' g1=',g_1
       ! check whether the solution satisfies the sufficient decrease condition
       ! if (f < f_0 + c * alpha * g_0) then
           if ((f <= f_0 + c * alpha * g_0).and.(abs(g_1) <= c2*abs(g_0))) then
               write(*,'(a53)') 'Strong Wolfe Condition satisfied, exiting line search'
               write(ioLog,'(a53)') 'Strong Wolfe Condition satisfied, exiting line search'
               call deall_dataVectorMTX(dHat_1)
               call deall_modelParam(mHat_0)
               call deall_modelParam(mHat_1)
               call deall_solnVectorMTX(eAll_1)
               return
           else
               if ((alpha_r-alpha_l)*g_1<0) then
                   ! update the left boundary for alpha
                   alpha_l = alpha
               else
                   ! update the right boundary for alpha
                   alpha_r = alpha
               endif
               write(6,*) 'alpha_l=',alpha_l,'alpha_r=',alpha_r
               write(ioLog,*) 'alpha_l=',alpha_l,'alpha_r=',alpha_r
           endif
       ! if not, iterate, using the two most recent values of f & alpha
           alpha_i = alpha_j
           f_i = f_j
           alpha_j = alpha
           f_j = f
           if (ibracket >= 10) then
               write(*,'(a69)') 'Warning: exiting bracketing since the it has iterated too many times!'
               write(ioLog,'(a69)') 'Warning: exiting bracketing since the it has iterated too many times!'
               if (f < f_1) then
                   starting_guess = .false.
               else 
                   starting_guess = .true.
                   call gradient(lambda,d,m0,mHat_1,grad,dHat_1,eAll_1)
               endif
               exit
           endif
       ! check that the function still decreases to avoid infinite loops in case of a bug
           if (abs(f_j - f_i) < TOL6) then
               write(*,'(a69)') 'Warning: exiting cubic search since the function no longer decreases!'
               write(ioLog,'(a69)') 'Warning: exiting cubic search since the function no longer decreases!'
               if (f < f_1) then
                   starting_guess = .false.
               else 
                   starting_guess = .true.
                   call gradient(lambda,d,m0,mHat_1,grad,dHat_1,eAll_1)
               endif
               exit
           endif
       end do fit_cubic
   end if
   if (starting_guess) then
       alpha = alpha_1
       dHat = dHat_1
       eAll = eAll_1
       mHat = mHat_1
       rms = rms_1
       f= f_1
   endif
   call deall_dataVectorMTX(dHat_1)
   call deall_modelParam(mHat_0)
   call deall_modelParam(mHat_1)
   call deall_solnVectorMTX(eAll_1)

  end subroutine lineSearchWolfe
  !**********************************************************************

end module NLCG
