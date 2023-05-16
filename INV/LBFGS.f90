module LBFGS

! inherits SensComp, DataIO and all modules they use. Also Main_MPI and Sub_MPI

use invcore

implicit none

public  :: LBFGSsolver

! iteration control for the LBFGS solver is initialized once
! and saved in the module to be used by most subroutines

  type  :: LBFGSiterControl_t
     ! maximum number of iterations in one call to iterative solver
     integer            :: maxIter
     ! convergence criteria: return from solver if rms < rmsTol
     real (kind=prec)   :: rmsTol
     ! the condition to identify when the inversion stalls
     real (kind=prec)   :: fdiffTol
     ! initial value of lambda (will not override the LBFGS input argument)
     real (kind=prec)   :: lambda
     ! exit if lambda < lambdaTol approx. 1e-4
     real (kind=prec)   :: lambdaTol
     ! set lambda_i = lambda_{i-1}/k when the inversion stalls
     real (kind=prec)   :: k
     ! the factor that ensures sufficient decrease in the line search
     real (kind=prec)   :: c
     ! the factor that ensures culvature condition in the line search
     real (kind=prec)   :: c2
     ! restart quasi Newton nQNmax iterations to ensure sufficient decend 
     integer            :: nQNmax ! just for book keeping 
     ! restart quasi Newton if orthogonality is lost (not necessarily needed)
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
  end type LBFGSiterControl_t

  type  :: modelParam_array_t ! container type
      type(modelParam_t), pointer             :: m
  end type

  type  :: LBFGSiterCache_t
     ! hard coded here to avoid too much memory use, or too few saves
     ! maxCache should be between 3 and 20
     integer                                    :: maxCache = 5
     integer                                    :: nCache
     type(modelParam_array_t), allocatable      :: deltaM(:),deltaG(:)
  end type LBFGSiterCache_t

  type(LBFGSiterControl_t), private, save :: iterControl
  type(LBFGSiterCache_t), private, save :: iterCache

Contains

!**********************************************************************
   subroutine set_LBFGSiterControl(iterControl)
   !  at present this is a simple routine to setup for
   !  iteration control; allowing parameters to be read from
   !   a file would be desireable.

   type(LBFGSiterControl_t), intent(inout)  :: iterControl

     ! maximum number of iterations in one call to iterative solver
     iterControl%maxIter = 300
     ! convergence criteria: return from solver if rms < rmsTol
     iterControl%rmsTol  = 1.05
     ! inversion stalls when abs(rms - rmsPrev) < fdiffTol (2e-3 works well)
     iterControl%fdiffTol = 2.0e-3
     ! initial value of lambda (will not override the LBFGS input argument)
     iterControl%lambda = 10.
     ! exit if lambda < lambdaTol approx. 1e-4
     iterControl%lambdaTol = 1.0e-4 ! makes no sense to use too small a value
     ! set lambda_i = lambda_{i-1}/k when the inversion stalls
     iterControl%k = 10.
     ! the factor that ensures sufficient decrease in the line search >=1e-4
     iterControl%c = 1e-4 
     ! the factor that ensures culvature condition in the line search c<c2<1
     iterControl%c2 = 0.9 ! use a value larger than 0.5
     ! restart QN every nQNmax iterations to ensure sufficient descend 
     iterControl%nQNmax = 8
     ! the starting step for the line search
     iterControl%alpha_1 = 20.
     ! maximum initial delta mHat (overrides alpha_1)
     iterControl%startdm = 20.
     ! optional relaxation parameter (Renormalised Steepest Descent algorithm)
     iterControl%gamma = 0.99
     ! model and data output file name
     iterControl%fname = 'YES_I_AM_LAZY'

   end subroutine set_LBFGSiterControl


   ! **************************************************************************
   ! * read_LBFGSiterControl reads the inverse solver configuration from file

   subroutine read_LBFGSiterControl(iterControl,rFile,fileExists)

    type(LBFGSiterControl_t), intent(inout)  :: iterControl
    character(*), intent(in)                :: rFile
    logical, intent(out), optional          :: fileExists
    integer                                 :: ios
    logical                                 :: exists
    character(80)                           :: string

    ! Initialize inverse solver configuration

    call set_LBFGSiterControl(iterControl)

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

   end subroutine read_LBFGSiterControl

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
   subroutine LBFGSsolver(d,lambda,m0,m,fname)

   ! computes inverse solution minimizing penalty functional
   !  for fixed value of regularization parameter, using
   !  limited memory Broyden-Fletcher-Goldfarb-Shanno method, 
   !  as described by Nocedal, 1980
   !  Various flavours of the algorithm and of the line search
   !  can be called from this routine
   !  NOTE: Quasi-Newtonian methods should only be used with line searches
   !  which follows the (Strong/Weak) Wolfe condition, i.e. the 
   !  curvature condition should be satisfed to ensure stable Quasi-Newton
   !  iterations. 
   !   
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
   character(80)                                :: flavor = 'Wolfe2'

   !  local variables
   type(dataVectorMTX_t)                        :: dHat, res
   type(modelParam_t)                           :: mHat, m_minus_m0
   type(modelParam_t)                           :: grad, g, h, gPrev
   type(modelParam_t)                           :: dM, dG
   type(LBFGSiterCache_t)                       :: saved

   real(kind=prec)                              :: value, valuePrev, rms
   real(kind=prec)                              :: rmsPrev, alpha, alphaPrev
   real(kind=prec)                              :: beta = 1.0
   real(kind=prec)                              :: gnorm, mNorm, Nmodel
   real(kind=prec)                              :: grad_dot_h !, g_dot_g
   integer                                      :: iter, flag, nLS, nfunc, ios
   integer                                      :: nQN, nQNmax, precType = 1
   logical                                      :: ok
   character(3)                                 :: iterChar
   character(100)                               :: mFile, mHatFile, gradFile
   character(100)                               :: dataFile, resFile, logFile
   type(solnVectorMTX_t)                        :: eAll

   if (present(fname)) then
      call read_LBFGSiterControl(iterControl,fname,ok)
      if (ok) then
         lambda = iterControl%lambda
      end if
   else
      call set_LBFGSiterControl(iterControl)
   end if

   ! initialize the output to log file
   logFile = trim(iterControl%fname)//'_LBFGS.log'
   open (unit=ioLog,file=logFile,status='unknown',position='append',iostat=ios)

   ! initialize the line search
   alpha = iterControl%alpha_1
   startdm = iterControl%startdm
   nQNmax = iterControl%nQNmax 

   write(*,'(a41,es8.1)') 'The initial damping parameter lambda is ',lambda
   write(*,'(a55,f12.6)') 'The initial line search step size (in model units) is ',startdm

   write(ioLog,'(a41,es8.1)') 'The initial damping parameter lambda is ',lambda
   write(ioLog,'(a55,f12.6)') 'The initial line search step size (in model units) is ',startdm


   ! starting from the prior hardcoded by setting mHat = 0 and m = m0
   ! m = m0
   ! mHat = m0
   ! call zero(mHat)

   ! starting model contains the rough deviations from the prior
   ! m_tilde 
   mHat = m

   !  compute the penalty functional and predicted data
   call func(lambda,d,m0,mHat,value,mNorm,dHat,eAll,rms)
   call printf('START',lambda,alpha,value,mNorm,rms)
   call printf('START',lambda,alpha,value,mNorm,rms,logFile)
   ! initial function call
   nfunc = 1

   write(iterChar,'(i3.3)') 0
   ! output (smoothed) initial model and responses for later reference
   ! 1. m - m_0 = C_m^(1/2) * \tilde{m} 
   ! 2. m = m_0 + (m - m_0)
   call CmSqrtMult(mHat,m_minus_m0)
   call linComb(ONE,m_minus_m0,ONE,m0,m)
   if (output_level > 1) then
     mFile = trim(iterControl%fname)//'_LBFGS_'//iterChar//'.rho'
     call write_modelParam(m,trim(mFile))
   end if
   if (output_level > 2) then
     dataFile = trim(iterControl%fname)//'_LBFGS_'//iterChar//'.dat'
     call write_dataVectorMTX(dHat,trim(dataFile))
   end if

   ! compute gradient of the full penalty functional
   call gradient(lambda,d,m0,mHat,grad,dHat,eAll)
   if (output_level > 3) then
     gradFile = trim(iterControl%fname)//'_LBFGS_'//iterChar//'.grt'
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

   ! initialize QN cache for deltaM and deltaG:
   call init_LBFGSiterCache(saved)
   ! g = - grad; h = g
   nQN = 0
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
      ! grad_dot_h = dotProd(grad,h)
      ! at the end of line search, set mHat to the new value
      ! mHat = mHat + alpha*h  and evaluate gradient at new mHat
      ! data and solnVector only needed for output
      write(*,'(a23)') 'Starting line search...'
      write(ioLog,'(a23)') 'Starting line search...'
      select case (flavor)
      case ('Cubic')
          call lineSearchCubic(lambda,d,m0,h,alpha,mHat,value,grad,rms,nLS,dHat,eAll)
          beta = 1.0
      case ('Quadratic')
          call lineSearchQuadratic(lambda,d,m0,h,alpha,mHat,value,grad,rms,nLS,dHat,eAll)
          beta = 1.0
      case ('Wolfe') ! wolfe with 2 FWD and 1 TRN per iteration
          call lineSearchWolfe(lambda,d,m0,h,alpha,mHat,value,grad,rms,nLS,dHat,eAll,flag)
          beta = 1.0
      case ('Wolfe2')! wolfe with 1 FWD and 1 TRN per iteration
          call lineSearchWolfe2(lambda,d,m0,h,alpha,mHat,value,grad,rms,nLS,dHat,eAll,flag)
          ! adjustment for beta, which controls the lambda update criterial
          beta = 0.3
      case default
          call errStop('Unknown line search requested in LBFGS')
      end select
      nfunc = nfunc + nLS
      ! save the previous (-) grad
      gPrev = g
      ! g = -grad here
      call linComb(MinusONE,grad,R_ZERO,grad,g)
      ! save the previous alpha 
      alphaPrev = alpha
      ! update the starting step for the next line search
      ! alpha = 2*(value - valuePrev)/grad_dot_h ! used in CG
      alpha = ONE 
      ! alpha = 1 should always be sufficient descend as the H_k^0 is 
      ! scaled with dm_dot_dg/dg_dot_dg, which ensures that the search
      ! direction is also scaled. 
      ! adjust the starting step to ensure superlinear convergence properties
      ! alpha = (ONE+0.01)*alpha ! used in CG
      ! alpha = max(ONE, alpha) ! used in CG
      write(*,'(a25,i5)') 'Completed LBFGS iteration ',iter
      write(ioLog,'(a25,i5)') 'Completed LBFGS iteration ',iter
      Nmodel = countModelParam(mHat)
      mNorm = dotProd(mHat,mHat)/Nmodel
      call printf('with',lambda,alpha,value,mNorm,rms)
      call printf('with',lambda,alpha,value,mNorm,rms,logFile)

      ! write out the intermediate model solution and responses
      call CmSqrtMult(mHat,m_minus_m0)
      call linComb(ONE,m_minus_m0,ONE,m0,m)
      write(iterChar,'(i3.3)') iter
      if (output_level > 1) then
        mFile = trim(iterControl%fname)//'_LBFGS_'//iterChar//'.rho'
        call write_modelParam(m,trim(mFile))
      end if
      if (output_level > 2) then
        mHatFile = trim(iterControl%fname)//'_LBFGS_'//iterChar//'.prm'
        call write_modelParam(mHat,trim(mHatFile))
      end if
      if (output_level > 2) then
        dataFile = trim(iterControl%fname)//'_LBFGS_'//iterChar//'.dat'
        call write_dataVectorMTX(dHat,trim(dataFile))
      end if
      ! compute residual for output: res = d-dHat; do not normalize by errors
      if (output_level > 2) then
        res = d
        call linComb(ONE,d,MinusONE,dHat,res)
        resFile = trim(iterControl%fname)//'_LBFGS_'//iterChar//'.res'
        call write_dataVectorMTX(res,trim(resFile))
      end if

      ! if alpha is too small, we are not making progress: update lambda
      ! the default criteria is on rms only
      ! if (abs(rmsPrev - rms) < iterControl%fdiffTol * beta) then
      ! I would recommend using this (object function) instead 
      if ((valuePrev-value)/value < iterControl%fdiffTol * beta) then
          ! update lambda, penalty functional and gradient
          call update_damping_parameter(lambda,mHat,value,grad)
          if (lambda < iterControl%lambdaTol) then
              write(*,'(a55)') 'Unable to get out of a local minimum. Exiting...'
              write(ioLog,'(a55)') 'Unable to get out of a local minimum. Exiting...'
              exit
          endif
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
          ! restart
          write(*,'(a55)') 'Restarting LBFGS with the damping parameter updated'
          call printf('to',lambda,alpha,value,mNorm,rms)
          write(ioLog,'(a55)') 'Restarting LBFGS with the damping parameter updated'
          call printf('to',lambda,alpha,value,mNorm,rms,logFile)
          h = g
          nQN = 0
          ! reset the saved dM and dG cache also
          saved%nCache = 0
          cycle  ! skip to descend direction
      endif
      ! Wolfe line search failed (for some reason) 
      if ((flag.eq.-1) .and. (nQN .ge. nQNmax) ) then ! reset the QN search direction
          ! use the gradient at current model
          ! need to reset the Hessian cache!
          gnorm = sqrt(dotProd(grad,grad))
          alpha = min(ONE,startdm)/gnorm
          ! g = - grad
          call linComb(MinusONE,grad,R_ZERO,grad,g)
          write(*,'(a42)') 'Restarting LBFGS with Hessian cache reset'
          call printf('to',lambda,alpha,value,mNorm,rms)
          write(ioLog,'(a42)') 'Restarting LBFGS with Hessian cache reset'
          call printf('to',lambda,alpha,value,mNorm,rms,logFile)
          h = g
          nQN = 0
          ! reset the saved dM and dG cache also
          saved%nCache = 0
          cycle  ! skip to descend direction 
      else
          nQN = nQN + 1
      endif
      ! dM = alphaPrev * h
      call linComb(alphaPrev, h, R_ZERO, h, dM)
      ! dG = grad - gradPrev
      ! but note that  g = -grad and gPrev = - gradPrev here
      call linComb(MinusONE, g, ONE, gPrev, dG)
      ! save dM and dG in our object...
      call update_LBFGSiterCache(saved,dM,dG) 
      call applyPrecond(grad,h,precType,ONE,ONE,saved)
      call linComb(MinusONE, h, R_ZERO, grad, h) 
      ! call update_Hessian(h,saved,grad)
      write(*,*) 'Hessian updated with results from previous ', saved%nCache, ' iteration(s)'
      write(ioLog,*) 'Hessian updated with results from previous ', saved%nCache, ' iteration(s)'
   end do
   ! multiply by C^{1/2} and add m_0
   call CmSqrtMult(mHat,m_minus_m0)
   call linComb(ONE,m_minus_m0,ONE,m0,m)
   d = dHat
   write(*,'(a25,i5,a25,i5)') 'LBFGS iterations:',iter,' function evaluations:',nfunc
   write(ioLog,'(a25,i5,a25,i5)') 'LBFGS iterations:',iter,' function evaluations:',nfunc
   close(ioLog,iostat=ios)

   ! cleaning up
   call deall_LBFGSiterCache(saved) 
   call deall_dataVectorMTX(dHat)
   call deall_dataVectorMTX(res)
   call deall_modelParam(mHat)
   call deall_modelParam(dM)
   call deall_modelParam(dG)
   call deall_modelParam(m_minus_m0)
   call deall_modelParam(grad)
   call deall_modelParam(g)
   call deall_modelParam(h)
   call deall_modelParam(gPrev)
   call deall_solnVectorMTX(eAll)

   end subroutine LBFGSsolver

  !**********************************************************************
  ! NOTE: one should not use this for LBFGS as LBFGS requires the 
  ! curvature condition (although you *can*, use at your own risk )
  subroutine lineSearchQuadratic(lambda,d,m0,h,alpha,mHat,f,grad, &
 &           rms,niter,dHat,eAll,gamma)

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
   ! To the best of my knowledge, it is not useful for LBFGS.

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
   logFile = trim(iterControl%fname)//'_LBFGS.log'

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

   ! alpha_1 is the initial step size, which is set in LBFGS
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
      STOP
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
  ! NOTE: one should not use this for LBFGS as LBFGS requires the 
  ! curvature condition (although you *can*, use at your own risk )
  subroutine lineSearchCubic(lambda,d,m0,h,alpha,mHat,f,grad, &
   & rms,niter,dHat,eAll,gamma)

   ! Line search that is based on the Numerical Recipes and on the
   ! text by Michael Ferris, Chapter 3, p 59. We only test the sufficient
   ! decrease (Armijo) condition (ignoring the curvature condition).
   ! We first interpolate using a quadratic approximation; if the
   ! solution does not satisfy the condition, we backtrack using
   ! cubic interpolation. This strategy only requires one gradient
   ! evaluation and is very efficient when computing gradients is
   ! expensive.
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
   ! Our solution has to satisfy the sufficient decrease condition
   !     f(alpha) < f(0) + c alpha f'(0).
   !
   ! The optional relaxation parameter gamma is needed for algorithms
   ! like the Renormalised Steepest Descent (RSD). See the dynamical
   ! systems in optimisation research (Pronzato et al [2000, 2001]).
   ! To the best of my knowledge, it is not useful for LBFGS.

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
   logFile = trim(iterControl%fname)//'_LBFGS.log'

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

   ! alpha_1 is the initial step size, which is set in LBFGS
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
		STOP
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
  ! line search with 2 FWD and 1 TRN with strong Wolfe Condition
  subroutine lineSearchWolfe(lambda,d,m0,h,alpha,mHat,f,grad, &
   & rms,niter,dHat,eAll,flag)
 
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
   ! but hey, penalty function evaluation is equally expensive, and we need 
   ! to calculate the gradient function for the next search direction anyway 
   ! - sooner or later, so here is my idea: 
   ! when we have evaluated the penalty function for the initial guess (f_1)
   ! and quadratic interpolation (f)
   ! 1) if f < f_1, we just continue to calculate the gradient at the
   !    quadratic interpolation point, and use it to test the Wolfe condition
   !    if it satisfies, then the calculation is the same as Anna's scheme 
   !    (2 penalty funtion evaluations and 1 gradient).
   ! 2) if f_1 < f, we do the same as 1), but only calculate grad at initial
   !    guess, still the calculation will be the same if everything works
   ! 3) if neither of the 1) and 2) is satisfied, the scheme falls back to the 
   !    cubic interpolation and sectioning with More-Thuente scheme
    
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
   ! condition, cannot guarantee the stable converge of Quasi-Newtonian method

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
   integer, intent(inout), optional          :: flag
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
   logFile = trim(iterControl%fname)//'_LBFGS.log'

   ! initialize the line search
   niter = 0
   mHat_0 = mHat
   f_0 = f
   starting_guess = .false. 

   if (present(flag)) then
      flag = 1 ! Wolfe condition not satisfied
   end if
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
       STOP
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
       if (present(flag)) then
           flag = -1 ! quadratic has no minimum - need restart
       end if
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
           if (present(flag)) then
               flag = 0 ! successful line search
           end if
           call deall_dataVectorMTX(dHat_1)
           call deall_modelParam(mHat_0)
           call deall_modelParam(mHat_1)
           call deall_solnVectorMTX(eAll_1)
           return
       else 
           !ooops, we missed the Strong Wolfe's condition (for one reason or 
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
               if (present(flag)) then
                   flag = 0 ! successful line search
               end if
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
           !ooops, we missed the Strong Wolfe's condition (for one reason or 
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

   if (f_1 > f_0) then
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
       if (f_1 < f) then ! pick the less-bad solution
           starting_guess = .true.
       endif
       write(*,'(a75)') 'Unable to fit a quadratic due to bad gradient estimate, exiting line search'
       write(ioLog,'(a75)') 'Unable to fit a quadratic due to bad gradient estimate, exiting line search'
       if (present(flag)) then
           flag = -1 ! bad gradient - need restart 
       end if
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
       ! interpolate with cubic 
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
                   write(*,'(a35)') 'b > 0, failed to fit CUBIC INTERP!'
                   write(ioLog,'(a35)') 'b > 0, failed to fit CUBIC INTERP!'
                   write(*,'(a35)') 'using default value to interpolate...'
                   write(ioLog,'(a35)') 'using default value to interpolate...'
                   alpha = -g_0/(b+sqrt(b*b - 3.0*a*g_0))
               endif
           endif
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
               if (present(flag)) then
                   flag = 0 ! successful line search
               end if
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
           if (ibracket >= 3) then
               write(*,'(a69)') 'Warning: exiting bracketing since the it has iterated too many times!'
               write(ioLog,'(a69)') 'Warning: exiting bracketing since the it has iterated too many times!'
               if (present(flag)) then
                   flag = 1 ! Wolfe condition not satisfied
               end if
               if (f < f_1) then
                   starting_guess = .false.
               else 
                   starting_guess = .true.
                   call gradient(lambda,d,m0,mHat,grad,dHat,eAll)
               endif
               exit
           endif
       ! check that the function still decreases to avoid infinite loops in case of a bug
           if (abs(f_j - f_i) < TOL6) then
               write(*,'(a69)') 'Warning: exiting cubic search since the function no longer decreases!'
               write(ioLog,'(a69)') 'Warning: exiting cubic search since the function no longer decreases!'
               if (present(flag)) then
                   flag = 1 ! Wolfe condition not satisfied
               end if
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

  ! line search with 1 FWD and 1 TRN with strong Wolfe Condition
  subroutine lineSearchWolfe2(lambda,d,m0,h,alpha,mHat,f,grad, &
   & rms,niter,dHat,eAll,flag,logid)
 
   ! Note: inexact line searches ultimately fit into two catalogs - 
   ! i.e. Armijo–Goldstein (back-tracking) and Wolfe conditions
   ! the latter also contains a few different criterias: (1) "Armijo 
   ! rule", and (2) "curvature condition", the latter requires an extra 
   ! gradient calculation at x_k + alpha_k*d_k
   ! (1) and modified (2) will form the "strong Wolfe" condition which is 
   ! the line search "convergence criteria" here

   ! see lineSearchWolfe subroutine for details of Anna's line search method
   ! (Armijo)
   ! 
   ! the really good part for LBFGS is, its search directions are sort of 
   ! "normalized". so there is a pretty good chance that the first guess (one)
   ! is already a good search step that satisfies Wolfe's rule.
   ! For most of the time, it is not even necessary to do a line "search".
   ! Also, if one examines the LBFGS line search (as in lineSearchWolfe), it
   ! is quite often that the f at quadratic interpolation does not improve 
   ! too much when compared with the first guess f_1. 
   ! 
   ! So here is my idea: 
   !
   ! when we have evaluated the penalty function for the initial guess (f_1), 
   ! we immediately calculate the gradient (g_1) at the initial guess.  
   ! 1) if f_1 and g_1 satisfy the Wolfe's rule, then we just proceed to the
   !    next iteration. The calculation cost is 1 penalty funtion evaluation
   !    and 1 gradient calculation. Everyone is happy(?)
   ! 
   ! 2) if f_1 and g_1 does not satify Wolfe's rule, we use a quadratic interp
   !    to find a minimum (f), and calculate the gradient (g) at the point. 
   !    if f and g satisfy the rule, then we proceed to the next iteration.
   !    The cost is 2 penalty function evaluations and 2 gradient calculations.
   !    this is actually worse than Anna's scheme (two f and one grad). 
   !    but don't worry yet, normally this should be a rare case (<5 %)
   !
   ! 3) if neither of the 1) nor 2) is satisfied, the scheme falls back to the 
   !    bracketing and sectioning with my variant of More-Thuente
   !    scheme
   !
   ! following Anna's idea, the major intention here is to use a cheap line
   ! search scheme (with merely 2 forward-like-calculations) to quickly skip 
   ! to a small overall penalty function level and only to start bracketing
   ! when the quadratic interpolation doesn't work, in which case cubic 
   ! probably won't work either...
   ! 
   ! in practice, wolfe2 should be slower (as the line search is not exact)
   ! in convergence (but still faster in time as the line search scheme is 
   ! 1/3 faster. 
   ! this may lead to quicker stop (or the update of the damping factor), as
   ! the descend may not meet the requirement of the "enough progress". 
   ! so I made a small modification in beta to allow the LBFGS to stall a
   ! little longer in each stage. 

   implicit none
   real(kind=prec), intent(in)               :: lambda ! lagrange multiplier
   type(dataVectorMTX_t), intent(in)         :: d  ! data vector
   type(modelParam_t), intent(in)            :: m0 ! current model
   type(modelParam_t), intent(in)            :: h  ! search direction
   real(kind=prec), intent(inout)            :: alpha ! step size
   type(modelParam_t), intent(inout)         :: mHat  ! next model
   real(kind=prec), intent(inout)            :: f  ! penalty function
   type(modelParam_t), intent(inout)         :: grad ! previous(next) gradient
   real(kind=prec), intent(out)              :: rms
   integer, intent(out)                      :: niter
   type(dataVectorMTX_t), intent(out)        :: dHat
   type(solnVectorMTX_t), intent(inout)      :: eAll
   integer, intent(inout), optional          :: flag
   character(100),intent(in), optional       :: logid

   ! local variables
   type(modelParam_t)              :: mHat_0, mHat_1, grad_1 ! initial grad
   type(dataVectorMTX_t)           :: dHat_1
   type(solnVectorMTX_t)           :: eAll_1
   ! to be deallocated, the above 5 variables need

   real(kind=prec)                 :: alpha_1,alpha_0, mNorm
   real(kind=prec)                 :: alpha_i,alpha_j ! brackets
   real(kind=prec)                 :: f_i,f_j ! function values at brackets
   real(kind=prec)                 :: g_i,g_j ! derivativess at brackets
   real(kind=prec)                 :: alphaPrev,alphaNext, fPrev,gPrev
   real(kind=prec)                 :: smin,smax ! the min/max step lengths
   real(kind=prec)                 :: left, right ! bounds for bracket
   logical                         :: starting_guess
   integer                         :: istrapped, nbracket 
   real(kind=prec)                 :: eps,k,c,c2 ! need to clear
   real(kind=prec)                 :: g_0,g_1,g,f_0,f_1,rms_1
   character(100)                  :: logFile

   ! parameters 
   c = iterControl%c ! ensures sufficient decrease
   c2 = iterControl%c2 ! ensures culvature condition
   !k = iterControl%alpha_k
   !eps = iterControl%alpha_tol
   if (.not.present(logid)) then
       logFile = trim(iterControl%fname)//'_LBFGS.log'
   else 
       logFile = logid
   endif

   ! initialize the line search
   niter = 0
   mHat_0 = mHat
   f_0 = f
   starting_guess = .false. 

   if (present(flag)) then
      flag = -1 ! Wolfe condition not satisfied
   end if
   ! g_0 is the directional derivative of our line search function 
   ! f'(0) = (df/dm).dot.h
   g_0 = dotProd(grad,h)
   if (g_0 >= R_ZERO) then
   ! quickly, blame the gradient calculation while you can (?) 
       write(ioLog,'(a45)') 'UNABLE TO PROCEED DUE TO A BAD GRADIENT(>0)'
       write(*,'(a45)') 'UNABLE TO PROCEED DUE TO A BAD GRADIENT(>0)'
       write(ioLog,'(a10)') 'Exiting...'
       write(*,'(a10)') 'Exiting...'
       STOP
   end if
   ! ====================================================================== !
   ! evaluate the functional at the first guess
   ! ====================================================================== !

   ! alpha_0 is the initial step size, which is set in LBFGS
   alpha_0 = alpha
   alpha_1 = alpha
   mHat_1 = mHat_0
   ! mHat_1 = mHat_0 + dir*step
   call linComb(ONE,mHat_0,alpha_1,h,mHat_1)
   ! compute the trial m, f, dHat, eAll, rms
   call func(lambda,d,m0,mHat_1,f_1,mNorm,dHat_1,eAll_1,rms_1)
   call printf('STARTLS',lambda,alpha_1,f_1,mNorm,rms_1)
   call printf('STARTLS',lambda,alpha_1,f_1,mNorm,rms_1,logFile)
   niter = niter + 1

   if (f_1 - f_0 >= LARGE_REAL) then
   ! oops, we are pushing too far away
       write(ioLog,'(a40)') 'Try a smaller starting value of alpha'
       write(*,'(a40)') 'Try a smaller starting value of alpha'
       write(ioLog,'(a10)') 'Exiting...'
       write(*,'(a10)') 'Exiting...'
       STOP
   end if
   ! calculate the gradient at the first guess
   call gradient(lambda,d,m0,mHat_1,grad,dHat_1,eAll_1)
   g_1 = dotProd(grad, h)
   grad_1 = grad
   write(*,'(a29,es12.5)',advance='no') &
       '    GRAD: computed, with g0=',g_0
   write(*,'(a4,es12.5)') ' g1=',g_1
   write(ioLog,'(a29,es12.5)',advance='no') &
       '    GRAD: computed, with g0=',g_0
   write(ioLog,'(a4,es12.5)') ' g1=',g_1
   ! ====================================================================== !
   ! test if the initial guess satisfies the Wolfe's condition
   ! ====================================================================== !
   if ((f_1.lt.f_0 + c * alpha_1 * g_0).and.(abs(g_1).lt.c2*abs(g_0))) then 
       istrapped = 0
       starting_guess = .true.
   else if (g_1 .gt. R_ZERO) then ! case 2
       ! we already trapped the minimum (change of derivative sign)
       istrapped = 2
   else if (f_1 .ge. f_0) then ! case 1 
       ! we have trapped the minimum (f_1 > f_0) 
       istrapped = 1
   else 
       ! sadly, we have not trapped the minimum
       ! write(ioLog, *) 'bracketing not successful (yet)'
       ! write(ioLog, *) 'f_1 = ', f_1
       ! write(ioLog, *) 'f_0 + c *alpha_1*g_0 = ', f_0 + c* alpha_1*g_0
       ! write(ioLog, *) 'g_0 = ', g_0
       ! write(ioLog, *) 'g_1 = ', g_1
       istrapped = -1
   endif

   if (istrapped.eq.-1) then  ! we have not yet trapped the minimum
   ! ====================================================================== !
   ! setup the bracketing parameters, for now, to prepare for the worst!
   ! ====================================================================== !
       alpha_i = R_ZERO
       f_i = f_0
       g_i = g_0
       alpha_j = alpha_1
       f_j = f_1
       g_j = g_1
       nbracket = 0
       ! setup the min and max step size
       smin = R_ZERO
       smax = (f_i-(f_i*0.98))/(-g_i*c)
   ! let's see what cards do we have in our hands...
       alphaPrev = R_ZERO ! original point
       fPrev = f_0
       gPrev = g_0
       alpha = alpha_1    ! starting guess
       f = f_1
       g = g_1
   ! ====================================================================== !
   ! bracketing session: try to find an interval that contains the minimizer
   ! ====================================================================== !
       bracket_session: do 
       ! now test if we have located the bracket
       ! the first half of Wolfe condition
           if ((f <= f_0 + c * alpha * g_0).and.(abs(g) <= c2*abs(g_0))) then 
               ! surprise! we have the proper step lengh = alpha now
               istrapped = 0
               exit ! no need to go on 
           else if (g .gt. R_ZERO) then ! case 2
               ! congratulations, we have find the bracket
               ! i.e. change of derivative sign
               alpha_i = alpha
               f_i = f
               g_i = g
               alpha_j = alphaPrev
               f_j = fPrev
               g_j = gPrev
               istrapped = 2
               exit ! finishing bracketing session
           else if (f .ge. fPrev) then ! case 1
               ! congratulations, we have find the bracket
               ! minimizer should be between alpha and alphaPrev
               ! i.e. f is increasing comparing with f_0
               ! this actually assumes that gPrev < 0
               alpha_i = alphaPrev
               f_i = fPrev
               g_i = gPrev
               alpha_j = alpha
               f_j = f
               g_j = g
               istrapped = 1 
               exit ! finishing bracketing session
           else if (nbracket.ge.2) then ! tried too many times...
               ! by this point the alpha should be quite large
               ! we are already quite far from the f_0 and g_0 
               ! the functional space cannot be considered quadratic 
               ! anymore - need to reset Hessian cache, if any
               alpha_i = R_ZERO
               f_i = f_0
               g_i = g_0
               alpha_j = alpha
               f_j = f
               g_j = g
               istrapped = -1
               exit ! finishing bracketing session
           endif
           ! if none of the above satisifies, update alpha
           if (2*alpha - alphaPrev < smax) then
               left = alpha + (alpha - alphaPrev)
               right = min(smax , alpha+3.0*(alpha-alphaPrev))
               if (nbracket.eq.0) then
                   if (g.ge.gPrev) then !g is not very helpful
                       call pickAlphaQuadratic(left,right,alphaPrev,fPrev,&
                           gPrev,alpha,f,alphaNext)
                   else 
                       call pickAlphaSecant(left,right,alphaPrev,fPrev,&
                           gPrev,alpha,g,alphaNext)
                   endif
               else
                   call pickAlphaCubic(left,right,alphaPrev,fPrev,&
                       gPrev,alpha,f,g,alphaNext)
               endif
           else 
               alphaNext = smax
           endif
           alphaPrev = alpha
           alpha = alphaNext
           ! now store the previous position
           fPrev = f
           gPrev = g
           ! evaluate the function and derivetive at the new alpha
           call linComb(ONE,mHat_0,alpha,h,mHat)
           call func(lambda,d,m0,mHat,f,mNorm,dHat,eAll,rms)
           call printf('CUBICLS',lambda,alpha,f,mNorm,rms)
           call printf('CUBICLS',lambda,alpha,f,mNorm,rms,logFile)
           niter = niter + 1
           nbracket = nbracket +1
           call gradient(lambda,d,m0,mHat,grad,dHat,eAll)
           g = dotProd(grad, h) 
           if (f.lt.f_1) then
               alpha_1 = alpha
               dHat_1 = dHat
               eAll_1 = eAll
               mHat_1 = mHat
               grad_1 = grad
               f_1 = f
               rms_1 = rms
           endif
           write(*,'(a29,es12.5)',advance='no') '    GRAD: computed, with g0=',g_0
           write(*,'(a4,es12.5)') ' g=',g
           write(ioLog,'(a29,es12.5)',advance='no') '    GRAD: computed, with g0=',g_0
           write(ioLog,'(a4,es12.5)') ' g=',g
       end do bracket_session
   else ! we have (somewhat) trapped the minimum
       alpha_i = R_ZERO
       f_i = f_0
       g_i = g_0
       alpha_j = alpha_1
       f_j = f_1
       g_j = g_1
   endif
               
   if (istrapped.eq.1) then
       write(*,'(a45)') '!======bracketing successful (case 1) ======'
       write(ioLog,'(a45)') '!=======bracketing successful (case 1) ======'
   else if (istrapped.eq.2) then
       write(*,'(a45)') '!=======bracketing successful (case 2) ======'
       write(ioLog,'(a45)') '!=======bracketing successful (case 2) ======'
   else if (istrapped.eq.-1) then
       write(*,'(a50)') '!=======bracketing failed (case -1) ======'
       write(ioLog,'(a50)') '!=======bracketing failed (case -1)======'
   else if (istrapped.eq.0) then
       ! say nothing here
       ! write(*,'(a50)') '!======= good alpha found (case 0) ======'
       ! write(ioLog,'(a50)') '!======= good alpha found (case 0) ======'
   else
       write(*,'(a50)') '!=========== it is a TRAP! =============='
       write(ioLog,'(a50)') '!=========== it is a TRAP! =============='
       stop
   endif
   nbracket = 0
   ! ====================================================================== !
   ! sectioning session: try to find a good searching step within brackets
   ! ====================================================================== !
   section_session: do
       if (istrapped .eq. 0) then ! we already found a good step
           exit ! no need to go on 
       endif
       ! firstly reduce the interval to avoid infinity loop
       left = alpha_i + min(0.1,c2)*(alpha_j - alpha_i)
       right = alpha_j - 0.1*(alpha_j - alpha_i)
       if ((istrapped .eq. 1).and.(nbracket.eq.0)) then ! try quadratic 
           call pickAlphaQuadratic(left,right,alpha_i,f_i,&
               g_i,alpha_j,f_j,alphaNext)
           alphaPrev = alpha
           alpha = alphaNext
           ! compute the penalty functional
           call linComb(ONE,mHat_0,alpha,h,mHat)
           call func(lambda,d,m0,mHat,f,mNorm,dHat,eAll,rms)
           call printf('QUADLS',lambda,alpha,f,mNorm,rms)
           call printf('QUADLS',lambda,alpha,f,mNorm,rms,logFile)
       else if ((istrapped .eq. 2).and.(nbracket.eq.0)) then ! try quadratic 
           call pickAlphaSecant(left,right,alpha_i,f_i,&
               g_i,alpha_j,g_j,alphaNext)
           alphaPrev = alpha
           alpha = alphaNext
           ! compute the penalty functional
           call linComb(ONE,mHat_0,alpha,h,mHat)
           call func(lambda,d,m0,mHat,f,mNorm,dHat,eAll,rms)
           call printf('QUADLS',lambda,alpha,f,mNorm,rms)
           call printf('QUADLS',lambda,alpha,f,mNorm,rms,logFile)
       else if ((istrapped .eq. -1).and.(nbracket.eq.0)) then ! jump 
           call pickAlphaSecant(left,1.0D+2,alpha_i,f_i,&
               g_i,alpha_j,g_j,alphaNext)
           alphaPrev = alpha
           alpha = alphaNext
           ! compute the penalty functional
           call linComb(ONE,mHat_0,alpha,h,mHat)
           call func(lambda,d,m0,mHat,f,mNorm,dHat,eAll,rms)
           call printf('QUADLS',lambda,alpha,f,mNorm,rms)
           call printf('QUADLS',lambda,alpha,f,mNorm,rms,logFile)
       else ! cubic 
           call pickAlphaCubic(left,right,alpha_i,f_i,&
               g_i,alpha_j,f_j,g_j,alphaNext)
           alphaPrev = alpha
           alpha = alphaNext
           ! compute the penalty functional
           call linComb(ONE,mHat_0,alpha,h,mHat)
           call func(lambda,d,m0,mHat,f,mNorm,dHat,eAll,rms)
           call printf('CUBICLS',lambda,alpha,f,mNorm,rms)
           call printf('CUBICLS',lambda,alpha,f,mNorm,rms,logFile)
       endif
       niter = niter + 1
       nbracket = nbracket + 1
       ! firstly store the previous values
       fPrev = f
       gPrev = g
       !calculatie gradient to test the Wolfe condition
       call gradient(lambda,d,m0,mHat,grad,dHat,eAll)
       g = dotProd(grad, h) 
       if (f.lt.f_1) then
           alpha_1 = alpha
           dHat_1 = dHat
           eAll_1 = eAll
           mHat_1 = mHat
           grad_1 = grad
           f_1 = f
           rms_1 = rms
       endif
       write(*,'(a29,es12.5)',advance='no') '    GRAD: computed, with g0=',g_0
       write(*,'(a4,es12.5)') ' g=',g
       write(ioLog,'(a29,es12.5)',advance='no') '    GRAD: computed, with g0=',g_0
       write(ioLog,'(a4,es12.5)') ' g=',g
! ======================================================================= !
! check if the cubic interpolation satisfies the condition
! ======================================================================= !
       if ((f <= f_0 + c * alpha * g_0).and.(abs(g) <= c2*abs(g_0))) then 
           istrapped = 0
           exit ! no need to go on 
       elseif ((f > f_0 + alpha*c*g_0).or.f > f_i) then!update the interval j
           alpha_j = alpha
           f_j = f
           g_j = g
       elseif ((g .gt. R_ZERO))then!update the interval j
           alpha_j = alpha
           f_j = f
           g_j = g
       else ! update the interval i
           alpha_i = alpha
           f_i = f
           g_i = g
           if ((alpha_j-alpha_i)*g >= 0) then
               alpha_j = alphaPrev
               f_j = fPrev
               g_j = gPrev
           endif
       endif
       if (abs((alpha_j-alpha_i)*g_i) .le. (f_0*c-(f_0-f_i))) then
           write(*,'(a65)') 'WARNING: no alpha that satisfies Wolfe condition can be found'
           write(ioLog,'(a65)') 'WARNING: no alpha that satisfies Wolfe condition can be found'
           istrapped = -1
           if (f_1 < f) then
               starting_guess = .true.
           endif
           exit ! no need to go on 
       endif
       if (abs(alpha_j - alpha_i) .le. 1e-3*alpha_0) then
           ! we didn't find an accetable point
           write(*,'(a69)') 'WARNING: exiting sectioning since the section interval is too small!'
           write(ioLog,'(a69)') 'WARNING: exiting sectioning since the section interval is too small!'
           istrapped = -1
           if (f_1 < f) then
               starting_guess = .true.
           endif
           exit
       endif
       if (abs(alphaPrev-alpha)/abs(alpha_i-alpha_j) <= 0.01) then
       ! no good minimizer possible for Wolfe condition
           write(*,'(a55)') 'WARNING: exiting as the step difference is too small!'
           write(ioLog,'(a55)') 'WARNING: exiting as the step difference is too small!'
           istrapped = -1
           if (f_1 < f) then
               starting_guess = .true.
           endif
           exit ! no need to go on 
       endif
       if (nbracket .ge. 3) then
       write(*,'(a43)') 'WARNING: maximum sectioning number reached!'
           write(ioLog,'(a43)') 'WARNING: maximum sectioning number reached!'
           istrapped = -1
           if (f_1 < f) then
               starting_guess = .true.
           endif
           exit
       endif
   end do section_session

   if (istrapped.eq.0) then
       if (present(flag)) then
           flag = 0 ! Wolfe condition satisfied, just go ahead
       endif
       write(*,'(a47)') 'Wolfe Condition satisfied, exiting line search'
       write(ioLog,'(a47)') 'Wolfe Condition satisfied, exiting line search'
   else 
       if (present(flag)) then
           flag = -1 ! Wolfe condition not satisfied, need to restart
       endif
       write(*,'(a40)') 'Wolfe Condition NOT satisfied, abort...'
       write(ioLog,'(a40)') 'Wolfe Condition NOT satisfied, abort...'
   endif

   if (starting_guess) then
       if (istrapped.ne.0) then
           write(6, *) 'recalling the best model so far...'
           write(ioLog, *) 'recalling the best model so far...'
       endif
       alpha = alpha_1
       dHat = dHat_1
       eAll = eAll_1
       mHat = mHat_1
       rms = rms_1
       f= f_1
       grad = grad_1
   endif
   call deall_dataVectorMTX(dHat_1)
   call deall_modelParam(mHat_0)
   call deall_modelParam(mHat_1)
   call deall_modelParam(grad_1)
   call deall_solnVectorMTX(eAll_1)
  end subroutine lineSearchWolfe2

  !**********************************************************************
  subroutine pickAlphaCubic(xi,xj,x1,f1,g1,x2,f2,g2,xc)
  ! this subroutine finds a minimizer x within the given interval 
  ! [xi, xj] with a cubic polynomial, which interpolates f and f' 
  ! at a1 and a2
  ! we assumes that f(x1) = f1, f(x2) = f2, f'(x1) = g1, f'(x2) = g2
  ! why, we assumes the polynomial to be:
  ! y = c3*x^3 + c2*x^2 + c1*x^1 + c0*x^0
    implicit none
    real(kind=prec), intent(in)        :: xi, xj
    real(kind=prec), intent(in)        :: x1, f1, g1
    real(kind=prec), intent(in)        :: x2, f2, g2
    real(kind=prec), intent(out)       :: xc
  ! local variables
    real(kind=prec)                    :: fs1,fs2,fleft,fright
    real(kind=prec)                    :: left,right,fc
    real(kind=prec)                    :: c3, c2, c1, c0 ! coefficients
    complex(kind=prec)                 :: s1, s2
  ! here we do something like 'un-scale', or 'preconditioning'
  ! imagine we do the linesearch in preconditioned (z) space
  ! i.e. x1 -> 0 x2 -> x2 - x1
  ! firstly find the coefficients of the cubic system
    c3 = (g1+g2)*(x2-x1)-2*(f2-f1)          ! 3rd order
    c2 = 3*(f2-f1) -(2*g1+g2)*(x2-x1)        ! 2nd order
    c1 = (x2-x1)*g1                          ! 1st order 
    c0 = f1  ! 0th order, but not really useful in finding the minimum
    
  ! 'preconditioning'
    left = (xi - x1)/(x2 - x1)
    right = (xj - x1)/(x2 - x1)
    if (left .gt. right) then 
        ! swap the two if needed
        xc = right
        right = left
        left = xc
    endif
  ! now find the solutions of f'(x) = 0
  ! should be degree of 2 and have 2 solutions(?)
  ! FIXME, this (dcmplx) is not right, because kind=prec may not be double 
    call roots2(dcmplx(3.0*c3,R_ZERO),dcmplx(2.0*c2,R_ZERO),&
        dcmplx(c1,R_ZERO),s1,s2)
  ! we calculate polynomial values at left, s1, s2 and right
  ! throw them out if they are complex
    if (abs(IMAG(s1)).gt.1e-10) then
        fs1 = 1e10
  ! throw them out if they are out of boundary 
    else if ((REAL(s1).le.left)) then 
        fs1 = 1e10
    else if ((REAL(s1).gt.right)) then 
        fs1 = 1e10
    else
        fs1 = c3*real(s1)**3.0+c2*real(s1)**2.0+c1*real(s1)+c0
    endif
    if (abs(IMAG(s2)).gt.1e-10) then
        fs2 = 1e10
  ! throw them out if they are out of boundary 
    else if ((REAL(s2).le.left)) then 
        fs2 = 1e10
    else if ((REAL(s2).gt.right)) then 
        fs2 = 1e10
    else
        fs2 = c3*real(s2)**3.0+c2*real(s2)**2.0+c1*real(s2)+c0
    endif
    fleft = c3*left**3.0+c2*left**2.0+c1*left+c0
    fright = c3*right**3.0+c2*right**2.0+c1*right+c0
    ! now compare them! 
    fc = 1e9
    xc = R_ZERO
    if (fleft .lt. fc) then
        fc = fleft
        xc = left
    endif
    if (fright .lt. fc) then
        fc = fright
        xc = right
    endif
    if (fs1 .lt. fc) then
        fc = fs1
        xc = real(s1)
    endif
    if (fs2 .lt. fc) then
        fc = fs2
        xc = real(s2)
    endif
    xc = x1 + xc * (x2 - x1)
    return
  end subroutine pickAlphaCubic

  subroutine roots2(a,b,c,s1,s2)
  ! solve degree 2 polynomial 
  ! y = ax^2 + bx + c
    implicit none
    complex(kind=prec),intent(in)     :: a,b,c
    complex(kind=prec),intent(out)    :: s1,s2
  ! local variables
    real(kind=prec)                   :: r,sx,sy,ux,uy,vx,vy,wx,wy
  
     ux = REAL(b)*REAL(b) - IMAG(b)*IMAG(b) - 4.0 * REAL(a)*REAL(c) &
         + 4.0 * IMAG(a)*IMAG(c)
     uy = 2.0*REAL(b)*IMAG(b) - 4.0 * REAL(a)*IMAG(c) &
         - 4.0*IMAG(a)*REAL(c)
     r = sqrt(ux*ux + uy*uy)
     vx = sqrt((r+ux)/2.0)
     vy = sqrt((r-ux)/2.0)
     if (uy<R_ZERO) then
         vy = -vy
     endif
     wx = (-REAL(b) - vx)/2.0
     wy = (-IMAG(b) - vy)/2.0
     ux = (-REAL(b) + vx)/2.0
     uy = (-IMAG(b) + vy)/2.0
     r = REAL(a)*REAL(a) + IMAG(a)*IMAG(a)
     sx = (REAL(a)*wx + IMAG(a)*wy)/r
     sy = (REAL(a)*wy - IMAG(a)*wx)/r
     s1 = CMPLX(sx,sy)
     sx = (REAL(a)*ux + IMAG(a)*uy)/r
     sy = (REAL(a)*uy - IMAG(a)*ux)/r
     s2 = CMPLX(sx,sy)
  end subroutine roots2

  !**********************************************************************
  subroutine pickAlphaQuadratic(xi,xj,x1,f1,g1,x2,f2,xc)
  ! this subroutine finds a minimizer x within the given interval 
  ! [xi, xj] with a quadartic polynomial, which interpolates f and f'
  ! at a1 and a2
  ! we assumes that f(x1) = f1, f(x2) = f2, f'(x1) = g1
  ! why, we assumes the polynomial to be:
  ! y = c2*x^2 + c1*x^1 + c0*x^0
    implicit none
    real(kind=prec), intent(in)        :: xi, xj
    real(kind=prec), intent(in)        :: x1, f1, g1
    real(kind=prec), intent(in)        :: x2, f2
    real(kind=prec), intent(out)       :: xc
  ! local variables
    real(kind=prec)                    :: fs,fleft,fright
    real(kind=prec)                    :: left,right,fc
    real(kind=prec)                    :: c2, c1, c0 ! coefficients
    real(kind=prec)                    :: s
  ! here we do something like 'un-scale', or 'preconditioning'
  ! imagine we do the linesearch in preconditioned (z) space
  ! i.e. x1 -> 0 x2 -> x2 - x1
  ! firstly find the coefficients of the quadratic system
    c2 = ((f1 - f2) + g1*( x2- x1)) ! 2nd order
    c2 = -c2 / ((x2 - x1)*(x2 - x1)) 
    c1 = g1  ! 1st order
    c0 = f1  ! 0th order, but not really useful in finding the minimum

  ! 'preconditioning'
    left = (xi - x1)/(x2-x1)
    right =(xj - x1)/(x2-x1)
    if (left .gt. right) then 
        ! swap the two if needed
        xc = right
        right = left
        left = xc
    endif
  ! now find the solutions of f'(x) = 0
    s = -c1/(TWO*c2)
  ! we calculate polynomial values at left, s and right
  ! throw it out if it is out of boundary [left right] 
    if (s.lt.left) then 
        fs = 1e10
    else if (s.gt.right) then 
        fs = 1e10
    else
        fs = c2*s*s+c1*s+c0
    endif
    fleft = c2*left*left+c1*left+c0
    fright = c2*right*right+c1*right+c0
    ! now compare them! 
    fc = 1e9
    xc = R_ZERO
    if (fleft .lt. fc) then
        fc = fleft
        xc = left
    endif
    if (fright .lt. fc) then
        fc = fright
        xc = right
    endif
    if (fs .lt. fc) then
        fc = fs
        xc = s
    endif
    xc = x1 + xc *(x2 - x1)
  end subroutine pickAlphaQuadratic

  !**********************************************************************
  subroutine pickAlphaSecant(xi,xj,x1,f1,g1,x2,g2,xc)
  ! this subroutine finds a minimizer x within the given interval 
  ! [xi, xj] with a quadartic polynomial, which interpolates 
  ! f1, g1 and g2 
  ! we assumes that f(x1) = f1, f'(x2) = g2, f'(x1) = g1
  ! why, we assumes the polynomial to be:
  ! y = c2*x^2 + c1*x^1 + c0*x^0
    implicit none
    real(kind=prec), intent(in)        :: xi, xj
    real(kind=prec), intent(in)        :: x1, f1, g1
    real(kind=prec), intent(in)        :: x2, g2
    real(kind=prec), intent(out)       :: xc
  ! local variables
    real(kind=prec)                    :: fs,fleft,fright
    real(kind=prec)                    :: left,right,fc
    real(kind=prec)                    :: c2, c1, c0 ! coefficients
    real(kind=prec)                    :: s
  ! here we do something like 'un-scale', or 'preconditioning'
  ! imagine we do the linesearch in preconditioned (z) space
  ! i.e. x1 -> 0 x2 -> x2 - x1
  ! firstly find the coefficients of the quadratic system
    c2 = 0.5*((g1 - g2)/( x1- x2)) ! 2nd order
    c1 = g1  ! 1st order
    c0 = f1  ! 0th order, but not really useful in finding the minimum
    
  ! 'preconditioning'
    left = (xi - x1)/(x2-x1)
    right =(xj - x1)/(x2-x1)
    if (left .gt. right) then 
        ! swap the two if needed
        xc = right
        right = left
        left = xc
    endif
  ! now find the solutions of f'(x) = 0
    s = -c1/(TWO*c2)
  ! we calculate polynomial values at left, s and right
  ! throw it out if it is out of boundary [left right] 
    if (s.lt.left) then 
        fs = 1e10
    else if (s.gt.right) then 
        fs = 1e10
    else
        fs = c2*s*s+c1*s+c0
    endif
    fleft = c2*left*left+c1*left+c0
    fright = c2*right*right+c1*right+c0
    ! now compare them! 
    fc = 1e9
    xc = R_ZERO
    if (fleft .lt. fc) then
        fc = fleft
        xc = left
    endif
    if (fright .lt. fc) then
        fc = fright
        xc = right
    endif
    if (fs .lt. fc) then
        fc = fs
        xc = s
    endif
    xc = x1 + xc *(x2 - x1)
  end subroutine pickAlphaSecant
  !**********************************************************************
  subroutine init_LBFGSiterCache(cache,maxsave) 
     ! initialize the stored queue of deltaM and deltaG
     implicit none
     type(LBFGSiterCache_t),intent(inout)        :: cache 
     integer, intent(in), optional               :: maxsave
     ! local variables 
     integer                                     :: i
     type(modelParam_t)                          :: m
     
     if (present(maxsave)) then
         cache%maxCache = maxsave
     endif
     if (.NOT.allocated(cache%deltaM)) then
         allocate(cache%deltaM(cache%maxCache))
     endif
     if (.NOT.allocated(cache%deltaG)) then
         allocate(cache%deltaG(cache%maxCache))
     endif
     do  i=1,cache%maxCache ! allocate them one by one
         allocate( cache%deltaM(i)%m, source = m )
         allocate( cache%deltaG(i)%m, source = m )
     end do
     ! no saves by far
     cache%nCache = 0
     return
  end subroutine init_LBFGSiterCache

  subroutine deall_LBFGSiterCache(cache) 
     ! deallocate the stored queue of deltaM and deltaG
     implicit none
     type(LBFGSiterCache_t),intent(inout)        :: cache 
     ! local variables 
     integer                                     :: i
     do  i=1,cache%maxCache
         call deall_modelParam(cache%deltaM(i)%m)
         call deall_modelParam(cache%deltaG(i)%m)
     end do
     return
  end subroutine deall_LBFGSiterCache
  
  !**********************************************************************
  subroutine update_LBFGSiterCache(cache,dM,dG,Bs) 
      ! update the stored queue of deltaM and deltaG
      ! note the index here 1->n is new -> old 
      implicit none
      type(LBFGSiterCache_t),intent(inout)        :: cache
      type(modelParam_t), intent(in)              :: dM, dG
      type(modelParam_t), intent(in),optional     :: Bs
      ! local variables 
      integer                                     :: i
      real(kind=prec)                             :: sBs, tau, phi
      cache%nCache = cache%nCache + 1
      if (cache%nCache.gt.cache%maxCache) then
          cache%nCache = cache%maxCache
      endif
      do i = 2, cache%nCache
          cache%deltaM(i)%m = cache%deltaM(i-1)%m
          cache%deltaG(i)%m = cache%deltaG(i-1)%m
      end do
      if (present(Bs)) then
          sBs = dotProd(dM,Bs)
          tau = dotProd(dM,dG)
          tau = abs(tau/sBs)
      else
          tau = ONE
      endif
      if (tau .lt. 0.2)then 
          phi = 0.8/(ONE-tau)
      elseif (tau .gt. ONE + 3.0) then
          phi = 3.0 / (tau - ONE)
      else
          phi = ONE
      endif
      cache%deltaG(1)%m = dG
      if (present(Bs)) then
          call linComb(phi, dG, (ONE-phi), Bs, cache%deltaG(1)%m)
      endif
      cache%deltaM(1)%m = dM
      return
  end subroutine update_LBFGSiterCache

  subroutine applyPrecond(r,rPrime,option,s0,s1,cache)
      ! apply precondition C to r, essentially this calculates
      !     rPrime = C^-1 * r
      ! option = 0 --> Diagonal Matrix
      ! option = 1 --> Hessian (with BFGS update)
      implicit none
      type(modelParam_t), intent(in)            :: r
      type(modelParam_t), intent(out)           :: rPrime
      integer,intent(in),optional               :: option
      real (kind=prec),intent(in),optional      :: s0,s1
      type(LBFGSiterCache_t),intent(in),optional:: cache
      ! note that 
      ! cache%dM_k = m_k+1 - m_k
      ! cache%dG_k = g_k+1 - g_k

      ! local variables 
      real (kind=prec),dimension(20)            :: a,p
      real (kind=prec)                          :: dg_dot_dm, dm_dot_q
      real (kind=prec)                          :: dg_dot_dg, b
      real (kind=prec)                          :: gamma1,lambda1
      type(modelParam_t)                        :: q,z
      integer                                   :: i, nstore, opt
      

      if (.not.present(option)) then 
          opt = 1
      else
          opt = option
      endif
      if (.not.present(s0)) then 
          gamma1 = ONE
      else
          gamma1 = s0
      endif
      if (.not.present(s1)) then 
          lambda1 = ONE
      else
          lambda1 = s1
      endif
      if (.not.present(cache)) then 
          nstore = 0
      else
          nstore = cache%nCache
      endif
      !if (nstore .eq. 0) then
      !    opt = 0
      ! endif
      select case (opt)
      case (0) ! diagonal preconditoning
          q = r
          if (nstore .ge. 1) then
             dg_dot_dm = dotProd(cache%deltaG(1)%m,cache%deltaM(1)%m)
             dg_dot_dg = dotProd(cache%deltaG(1)%m,cache%deltaG(1)%m)
             gamma1 = dg_dot_dm/dg_dot_dg 
          endif
          call linComb(gamma1,q,R_ZERO,q,z)
      case (1) ! preconditioning with BFGS H_k, with classic 'two loops' theme
          q = r
          do i = 1,nstore
            dg_dot_dm = dotProd(cache%deltaG(i)%m, cache%deltaM(i)%m)
            p(i) = ONE / dg_dot_dm
            dm_dot_q = dotProd(cache%deltaM(i)%m,q)
            a(i) = p(i) * dm_dot_q
            call linComb(ONE,q,-a(i),cache%deltaG(i)%m,q)
          end do
          ! now calculate the scale for the initial Hessian
          ! with most recent dM and dG
          if (nstore .ge. 1) then
             dg_dot_dm = dotProd(cache%deltaG(1)%m,cache%deltaM(1)%m)
             dg_dot_dg = dotProd(cache%deltaG(1)%m,cache%deltaG(1)%m)
             gamma1 = dg_dot_dm/dg_dot_dg 
          endif
          ! z =  Hessian^-1 * r
          call linComb(gamma1,q,R_ZERO,q,z)
          do i = nstore, 1, -1 !nstore:-1:1
             b = dotProd(cache%deltaG(i)%m,z)
             b = p(i) * b
             call linComb(ONE, z,  (a(i)-b), cache%deltaM(i)%m, z) 
          end do
      case default
          call errStop('Unknown preconditoner requested.')
      end select
      rPrime = z
      call deall_modelParam(q)
      call deall_modelParam(z)
      
  end subroutine applyPrecond

  !**********************************************************************
  subroutine update_Hessian(h,cache,grad)
      ! update the QN search direction (h) with stored list of d_M and d_g
      ! the (inverse) Hessian approximation is actually calculated 
      ! on-the-fly
      type(modelParam_t), intent(inout)         :: h 
      type(LBFGSiterCache_t),intent(in)         :: cache
      type(modelParam_t), intent(in)            :: grad
      ! local variables 
      real (kind=prec),dimension(20)            :: a,p
      real (kind=prec)                          :: dg_dot_dm, dm_dot_q
      real (kind=prec)                          :: dg_dot_dg, s, b
      type(modelParam_t)                        :: q,r
      integer                                   :: i, nstore

      q = grad
      nstore = cache%nCache
      do i = 1,nstore
          dg_dot_dm = dotProd(cache%deltaG(i)%m, cache%deltaM(i)%m)
          p(i) = 1 / dg_dot_dm
          dm_dot_q = dotProd(cache%deltaM(i)%m,q)
          a(i) = p(i) * dm_dot_q
          call linComb(ONE,q,-a(i),cache%deltaG(i)%m,q)
      end do
      ! now calculate the scale for the initial Hessian
      ! with most recent dM and dG
      dg_dot_dm = dotProd(cache%deltaG(1)%m,cache%deltaM(1)%m)
      dg_dot_dg = dotProd(cache%deltaG(1)%m,cache%deltaG(1)%m)
      s = dg_dot_dm/dg_dot_dg 
      call linComb(s,q,R_ZERO,q,r)
      ! r =  - Hessian * Grad
      do i = nstore, 1, -1 !nstore:-1:1
          b = dotProd(cache%deltaG(i)%m,r)
          b = p(i) * b
          call linComb(ONE, r,  (a(i)-b), cache%deltaM(i)%m, r) 
      end do
      ! h = -r
      call linComb(MinusONE, r, R_ZERO, h, h)
      call deall_modelParam(q)
      call deall_modelParam(r)
      return
  end subroutine update_Hessian

end module LBFGS
