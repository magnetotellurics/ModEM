module LBFGS
! this is an attempt to implement inversion search without Gary's tranformed  
! model space idea, but still use the transformed data space (because I 
! don't see why we should not!)
! the motivation is, unlike the 1D AR, the laplacian roughing operator seems 
! to be make the once-well-conditioned transformed penalty function space 
! quite rough and ill-conditioned, probably because the operator is "too rough"
! (it is applied for twice in the penalty function, as the laplacian is only 
! "half" of the Cm operator). 
! 
! so - here are three ideas: 
! 1. make a less-rough operator, say, laplace^(1/2), which seems not obvious
!    how I can implement it (you can ask Gary for details)
! 2. retreat to the original model space - and find a proper (simple)
!    preconditioner for the non-linear optimization (like what everyone else 
!    has been doing).
! 3. find a cool preconditioner for the transformed model space and laplacian 
!    operator - not sure if this is a good idea as this would make the once
!    simplified gradient and penalty function evaluation complicated again...
!
! any of the above requires a considerable modification of the inversion
! scheme! not sure exactly how much modification it would take - 
! but in theory it is not that substential 
! btw, it's April the 1st of 2019 today - maybe it hints something! 

! basicly what I want to do is to minimize the penalty functional 
! Phi = (d - f(m))^T C_d^-1 (d - f(m)) + lambda*(m-m_0)^T C_m^-1 (m-m_0)
! we still want to rescale data to eliminate the C_d^-1
! by setting d_hat = C_d^-0.5 d and 
! Phi = (d_hat-F_hat(m))^T(d_hat-F_hat(m)) + lambda*(m-m_0)C_m^-1(m-m_0)

! 
! inherits SensComp, DataIO and all modules they use. Also Main_MPI and Sub_MPI

use invcore

implicit none

public  :: LBFGSsolver 

! iteration control for general solver is initialized once
! and saved in the module to be used by most subroutines

  type  :: SolverIterControl_t
     ! maximum number of iterations in one call to iterative solver
     integer            :: maxIter
     ! convergence criteria: return from solver if rms < rmsTol
     real (kind=prec)   :: rmsTol
     ! the condition to identify when the inversion stalls
     real (kind=prec)   :: fdiffTol
     ! initial value of lambda (will not override the solver input argument)
     real (kind=prec)   :: lambda
     ! exit if lambda < lambdaTol approx. 1e-4
     real (kind=prec)   :: lambdaTol
     ! set lambda_i = lambda_{i-1}/k when the inversion stalls
     real (kind=prec)   :: k
     ! the factor that ensures sufficient decrease in the line search
     real (kind=prec)   :: c
     ! the factor that ensures culvature condition in the line search
     real (kind=prec)   :: c2
     ! restart quasi Newton nrestart iterations to ensure sufficient decend 
     integer            :: nrestart ! just for book keeping 
     ! restart CG if orthogonality is lost (not necessarily needed)
     real (kind=prec)   :: delta ! 0.5
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
     character(80)      :: fname
  end type SolverIterControl_t

  type  :: modelParam_array_t ! container type
      type(modelParam_t), pointer             :: m
  end type

  type  :: LBFGSiterCache_t
     ! hard coded here to avoid too much memory use, or too few saves
     ! maxCache should be between 3 and 20
     ! 3 <--- cheap  (in case of memory cost)  expensive ---> 20
     ! 3 <--- crude (in case of Hessian matrix) accurate ---> 20
     integer                                    :: maxCache = 5
     integer                                    :: nCache
     type(modelParam_array_t), allocatable      :: deltaM(:),deltaG(:)
  end type LBFGSiterCache_t

  type(SolverIterControl_t), private, save :: iterControl
  type(LBFGSiterCache_t), private, save :: iterCache

Contains

!**********************************************************************
   subroutine set_SolverIterControl(iterControl)
   !  at present this is a simple routine to setup for
   !  iteration control; allowing parameters to be read from
   !   a file would be desireable.

   type(SolverIterControl_t), intent(inout)  :: iterControl

     ! maximum number of iterations in one call to iterative solver
     iterControl%maxIter = 300
     ! convergence criteria: return from solver if rms < rmsTol
     iterControl%rmsTol  = 1.05
     ! inversion stalls when abs(rms - rmsPrev) < fdiffTol (2e-3 works well)
     iterControl%fdiffTol = 2.0e-3
     ! initial value of lambda (will not override the solver input argument)
     iterControl%lambda = 10.
     ! exit if lambda < lambdaTol approx. 1e-4
     iterControl%lambdaTol = 1.0e-8
     ! set lambda_i = lambda_{i-1}/k when the inversion stalls
     iterControl%k = 10.
     ! the factor that ensures sufficient decrease in the line search >=1e-4
     iterControl%c = 1.0e-4
     ! the factor that ensures culvature condition in the line search c<c2<1
     iterControl%c2 = 0.9 ! use a value larger than 0.5
     ! restart QN every nrestart iterations to ensure sufficient descend 
     iterControl%nrestart = 8
     ! the starting step for the line search
     iterControl%alpha_1 = 20.
     ! maximum initial delta mHat (overrides alpha_1)
     iterControl%startdm = 20.
     ! optional relaxation parameter (Renormalised Steepest Descent algorithm)
     iterControl%gamma = 0.99
     ! model and data output file name
     iterControl%fname = 'YES_I_AM_LAZY'

   end subroutine set_SolverIterControl


   ! **************************************************************************
   ! * read_SolverIterControl reads the inverse solver configuration from file

   subroutine read_SolverIterControl(iterControl,rFile,fileExists)

    type(SolverIterControl_t), intent(inout)  :: iterControl
    character(*), intent(in)                :: rFile
    logical, intent(out), optional          :: fileExists
    integer                                 :: ios
    logical                                 :: exists
    character(80)                           :: string

    ! Initialize inverse solver configuration

    call set_SolverIterControl(iterControl)

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

   end subroutine read_SolverIterControl

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
   subroutine update_damping_parameter2(lambda,m,m0,F,grad)

   real(kind=prec), intent(inout)              :: lambda
   type(modelParam_t), intent(in)              :: m
   type(modelParam_t), intent(in)              :: m0
   real(kind=prec), intent(inout)              :: F
   type(modelParam_t), intent(inout)           :: grad
   real(kind=prec)                             :: SS, mNorm, Nmodel
   type(modelParam_t)                          :: mHat
   type(modelParam_t)                          :: dSS

   ! initialize mHat
   mHat = m
   dSS = mHat
   ! now calculate the model normal (again)
   call linComb(ONE, m, MinusONE, m0, dSS)
   ! mHat = C_m(-1)(m-m0)
   call CmInvMult(dSS,mHat)

   ! compute the model norm
   mNorm = dotProd(mHat,dSS)
   Nmodel = countModelParam(mHat)

   ! (scaled) sum of squares = penalty functional - scaled model norm
   SS = F - (lambda * mNorm/Nmodel)
   ! write(6,*) 'F =', F
   ! write(6,*) 'lambda* mNorm/Nmodel = ', lambda *mNorm/Nmodel
   ! write(6,*) 'SS = ',SS

   ! subtract the model norm derivative from the gradient of the penalty functional
   ! (scaled) derivative of sum of squares = grad - scaled model derivative
   call linComb(ONE,grad,MinusTWO*lambda/Nmodel,mHat,dSS)

   ! update the damping parameter lambda
   lambda = lambda/iterControl%k

   ! penalty functional = (scaled) sum of squares + scaled model norm
   F = SS + (lambda * mNorm/Nmodel)
   ! write(6,*) 'Fnew =', F

   ! add the model norm derivative to the gradient of the penalty functional
   call linComb(ONE,dSS,TWO*lambda/Nmodel,mHat,grad)

   call deall_modelParam(dSS)
   call deall_modelParam(mHat)

   end subroutine update_damping_parameter2

!**********************************************************************

   subroutine LBFGSsolver(d,lambda,m0,m,fname)

   ! computes inverse solution minimizing penalty functional
   !  for fixed value of regularization parameter, using
   !  limited memory Broyden-Fletcher-Goldfarb-Shanno method, 
   !  as described by Nocedal, 1980
   !  Various flavours of the algorithm and of the line search
   !  can be called from this routine
   !  NOTE: Quasi-Newtonish methods should only be used with line searches
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
   character(80)                                :: flavor = 'Wolfe'

   !  local variables
   type(dataVectorMTX_t)                        :: dHat, res
   type(modelParam_t)                           :: mHat, m_minus_m0
   type(modelParam_t)                           :: grad, gradPrev, z, h
   type(modelParam_t)                           :: dM, dG
   type(LBFGSiterCache_t)                       :: saved
   ! to be deallocated, the above 12 variables need

   real(kind=prec)                              :: value, valuePrev, rms
   real(kind=prec)                              :: rmsPrev, alpha, alphaPrev
   ! real(kind=prec)                              :: beta
   real(kind=prec)                              :: gnorm, mNorm, Nmodel
   integer                                      :: iter, nQN, nLS, nfunc, ios
   integer                                      :: flag, precType = 1
   logical                                      :: ok
   character(3)                                 :: iterChar
   character(100)                               :: mFile, mHatFile, gradFile
   character(100)                               :: dataFile, resFile, logFile
   type(solnVectorMTX_t)                        :: eAll

   if (present(fname)) then
      call read_SolverIterControl(iterControl,fname,ok)
      if (ok) then
         lambda = iterControl%lambda
      end if
   else
      call set_SolverIterControl(iterControl)
   end if

   ! initialize the output to log file
   logFile = trim(iterControl%fname)//'_LBFGS.log'
   open (unit=ioLog,file=logFile,status='unknown',position='append',iostat=ios)

   ! initialize the line search
   alpha = iterControl%alpha_1
   startdm = iterControl%startdm

   write(*,'(a41,es8.1)') 'The initial damping parameter lambda is ',lambda
   write(*,'(a55,f12.6)') 'The initial line search step size (in model units) is ',startdm

   write(ioLog,'(a41,es8.1)') 'The initial damping parameter lambda is ',lambda
   write(ioLog,'(a55,f12.6)') 'The initial line search step size (in model units) is ',startdm

   ! starting model
   mHat = m ! in smoothed (untransformed) domain

   !  compute the penalty functional and predicted data
   call func2(lambda,d,m0,m,value,mNorm,dHat,eAll,rms)
   call printf('START',lambda,alpha,value,mNorm,rms)
   call printf('START',lambda,alpha,value,mNorm,rms,logFile)
   ! initial function call
   nfunc = 1
   write(iterChar,'(i3.3)') 0

   ! output initial model (0th) and responses for later reference
   if (output_level > 1) then
     mFile = trim(iterControl%fname)//'_LBFGS_'//iterChar//'.rho'
     call write_modelParam(m,trim(mFile))
   end if
   if (output_level > 2) then
     dataFile = trim(iterControl%fname)//'_LBFGS_'//iterChar//'.dat'
     call write_dataVectorMTX(dHat,trim(dataFile))
   end if

   ! compute gradient of the full penalty functional
   call gradient2(lambda,d,m0,m,grad,dHat,eAll)
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
   nQN = 0
   iter = 0
   ! z = C^-1*grad
   call applyPrecond(grad,z,PrecType,alpha,lambda,saved)
   ! h = -z
   call linComb(MinusONE,z,R_ZERO,z,h) 
   alpha = ONE
   do
      !  test for convergence ...
      if((rms.lt.iterControl%rmsTol).or.(iter.ge.iterControl%maxIter)) then
         exit
      end if
      iter = iter + 1
      ! save the values of the functional and the directional derivative
      rmsPrev = rms
      valuePrev = value
      ! save the previous grad
      gradPrev = grad
      ! at the end of line search, set m to the new value
      ! m = m + alpha*h  and evaluate gradient at new mHat
      ! data and solnVector only needed for output
      write(*,'(a23)') 'Starting line search...'
      write(ioLog,'(a23)') 'Starting line search...'
      select case (flavor)
      case ('Armijo')
          call lineSearchArmijo(lambda,d,m0,h,alpha,m,value,grad,rms,nLS,dHat,eAll,flag,logfile)
      case ('Wolfe')
          call lineSearchWolfe(lambda,d,m0,h,alpha,m,value,grad,rms,nLS,dHat,eAll,flag,logfile)
      case default
          call errStop('Unknown line search requested in LBFGS')
      end select
      nfunc = nfunc + nLS
      ! save the previous alpha 
      alphaPrev = alpha
      ! update the starting step for the next line search
      alpha = ONE 
      ! alpha = 1 should always be sufficient descend as the H_k^0 is 
      ! scaled with dm_dot_dg/dg_dot_dg, which ensures that the search
      ! direction is also scaled. 

      call linComb(ONE, m, MinusONE, m0, m_minus_m0) 
      call CmInvMult(m_minus_m0,mHat)
      Nmodel = countModelParam(mHat)
      mNorm = dotProd(m_minus_m0,mHat)/Nmodel
      write(*,'(a25,i5)') 'Completed LBFGS iteration ',iter
      write(ioLog,'(a25,i5)') 'Completed LBFGS iteration ',iter
      call printf('with',lambda,alpha,value,mNorm,rms)
      call printf('with',lambda,alpha,value,mNorm,rms,logFile)

      ! write out the intermediate model solution and responses
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
      if (output_level > 3) then
        res = d
        call linComb(ONE,d,MinusONE,dHat,res)
        resFile = trim(iterControl%fname)//'_LBFGS_'//iterChar//'.res'
        call write_dataVectorMTX(res,trim(resFile))
      end if

      if ((flag.eq.-1).and.(nQN.ne.0)) then ! need to reset the Hessian cache!
          ! use the gradient at current model
          gnorm = sqrt(dotProd(grad,grad))
          alpha = min(ONE,startdm/gnorm)
          saved%nCache = 0
          ! z = C^-1*grad
          call applyPrecond(grad,z,precType,alpha,lambda,saved)
          ! h = -z
          call linComb(MinusONE, z, R_ZERO, z, h)
          write(*,'(a42)') 'Restarting LBFGS with Hessian cache reset'
          write(ioLog,'(a42)') 'Restarting LBFGS with Hessian cache reset'
          alpha = ONE
          nQN = 0
          ! reset the saved dM and dG cache also
          cycle  ! skip to descend direction 
      elseif ((abs(rmsPrev - rms) < iterControl%fdiffTol).or.&
      ! if alpha is too small, we are not making progress: update lambda
          ((valuePrev-value)/value <= 0.0005)) then
          ! update lambda, penalty functional and gradient
          call update_damping_parameter2(lambda,m,m0,value,grad)
          if (lambda < iterControl%lambdaTol) then
              write(*,'(a55)') 'Unable to get out of a local minimum. Exiting...'
              write(ioLog,'(a55)') 'Unable to get out of a local minimum. Exiting...'
              exit
          endif
          ! update alpha
          gnorm = sqrt(dotProd(grad,grad))
          write(*,'(a34,es12.5)') 'The norm of the last gradient is ',gnorm
          write(ioLog,'(a34,es12.5)') 'The norm of the last gradient is ',gnorm
          alpha = min(ONE,startdm)/gnorm
          write(*,'(a48,es12.5)') 'The value of line search step alpha updated to ',alpha
          write(ioLog,'(a48,es12.5)') 'The value of line search step alpha updated to ',alpha
          ! restart
          write(*,'(a55)') 'Restarting LBFGS with the damping parameter updated'
          call printf('to',lambda,alpha,value,mNorm,rms)
          write(ioLog,'(a55)') 'Restarting LBFGS with the damping parameter updated'
          call printf('to',lambda,alpha,value,mNorm,rms,logFile)
          ! reset the saved dM and dG cache also
          saved%nCache = 0
          ! z = C^-1*grad
          call applyPrecond(grad,z,precType,alpha,lambda,saved)
          ! h = -z
          call linComb(MinusONE, z, R_ZERO, z, h)
          alpha = ONE
          nQN = 0
          ! reset the saved dM and dG cache also
          cycle  ! skip to descend direction
      ! Wolfe line search failed (for some reason) 
      endif
      ! else 
      nQN = nQN + 1
      ! dM = alphaPrev * h
      call linComb(alphaPrev, h, R_ZERO, h, dM)
      ! dG = grad - gradPrev
      call linComb(MinusONE, gradPrev, ONE, grad, dG)
      ! save dM and dG in our object...
      call update_LBFGSiterCache(saved,dM,dG) 
      ! z = C^-1*grad
      call applyPrecond(grad,z,precType,ONE,ONE,saved)
      ! h = -z
      call linComb(MinusONE, z, R_ZERO, z, h) 
      alpha = ONE
      write(*,*) 'Hessian updated with results from previous ', saved%nCache, ' iteration(s)'
      write(ioLog,*) 'Hessian updated with results from previous ', saved%nCache, ' iteration(s)'
   end do
   d = dHat
   write(*,'(a25,i5,a25,i5)') 'LBFGS iterations:',iter,' function evaluations:',nfunc
   write(ioLog,'(a25,i5,a25,i5)') 'LBFGS iterations:',iter,' function evaluations:',nfunc
   close(ioLog,iostat=ios)

   ! /active lightsaber
   call deall_LBFGSiterCache(saved) 
   call deall_dataVectorMTX(dHat)
   call deall_dataVectorMTX(res)
   call deall_modelParam(mHat)
   call deall_modelParam(dM)
   call deall_modelParam(dG)
   call deall_modelParam(m_minus_m0)
   call deall_modelParam(grad)
   call deall_modelParam(gradPrev)
   call deall_modelParam(z)
   call deall_modelParam(h)
   call deall_solnVectorMTX(eAll)

   end subroutine LBFGSsolver

  !**********************************************************************
  subroutine lineSearchWolfe(lambda,d,m0,h,alpha,m,f,grad, &
   & rms,niter,dHat,eAll,flag,logid)
 
   ! Note: inexact line searches ultimately fit into two catalogs - 
   ! i.e. Armijo–Goldstein (back-tracking) and Wolfe conditions
   ! the latter also contains a few different criterias: (1) "Armijo 
   ! rule", and (2) "curvature condition", the latter requires an extra 
   ! gradient calculation at x_k + alpha_k*d_k
   ! (1) and modified (2) will form the "strong Wolfe" condition which is 
   ! the line search "convergence criteria" here
   ! 
   ! the really good part for LBFGS is, its search directions are sort of 
   ! "normalized". so there is a pretty good chance that the first guess (one)
   ! is already a good search step that satisfies Wolfe's rule.
   ! For most of the time, it is not even necessary to do a line "search".
   ! Also, if one examines the LBFGS line search, it
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
   !    TODO: I should do the same for NLCG
   !
   ! following Anna's idea, the major intention here is to use a cheap line
   ! search scheme (with merely 2 forward-like-calculations) to quickly skip 
   ! to a small overall penalty function level and only to start bracketing
   ! when the quadratic interpolation doesn't work, in which case cubic 
   ! probably won't work either...
   ! 

   implicit none
   real(kind=prec), intent(in)               :: lambda ! lagrange multiplier
   type(dataVectorMTX_t), intent(in)         :: d  ! data vector
   type(modelParam_t), intent(in)            :: m0 ! prior model
   type(modelParam_t), intent(in)            :: h  ! search direction
   real(kind=prec), intent(inout)            :: alpha ! step size
   type(modelParam_t), intent(inout)         :: m  ! next model
   real(kind=prec), intent(inout)            :: f  ! penalty function
   type(modelParam_t), intent(inout)         :: grad ! previous(next) gradient
   real(kind=prec), intent(out)              :: rms
   integer, intent(out)                      :: niter
   type(dataVectorMTX_t), intent(out)        :: dHat
   type(solnVectorMTX_t), intent(inout)      :: eAll
   integer, intent(inout), optional          :: flag
   character(100),intent(in), optional       :: logid

   ! local variables
   type(modelParam_t)              :: m_0, m_1 !first guess 
   type(modelParam_t)              :: grad_1
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
   real(kind=prec)                 :: eps,k,c,c2 
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
   m_0 = m
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
   m_1 = m_0
   ! m_1 = m_0 + dir*step
   call linComb(ONE,m_0,alpha_1,h,m_1)
   ! compute the trial m, f, dHat, eAll, rms
   call func2(lambda,d,m0,m_1,f_1,mNorm,dHat_1,eAll_1,rms_1)
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
   call gradient2(lambda,d,m0,m_1,grad,dHat_1,eAll_1)
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
               ! we are already far from the f_0 and g_0 
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
           call linComb(ONE,m_0,alpha,h,m)
           call func2(lambda,d,m0,m,f,mNorm,dHat,eAll,rms)
           call printf('CUBICLS',lambda,alpha,f,mNorm,rms)
           call printf('CUBICLS',lambda,alpha,f,mNorm,rms,logFile)
           niter = niter + 1
           nbracket = nbracket +1
           call gradient2(lambda,d,m0,m,grad,dHat,eAll)
           g = dotProd(grad, h) 
           if (f.lt.f_1) then
               alpha_1 = alpha
               dHat_1 = dHat
               eAll_1 = eAll
               m_1 = m
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
       !if (f_j .gt. 1.5*f_i) then
       !    right = alpha_j - 0.618*(alpha_j - alpha_i)
       !else
       right = alpha_j - 0.1*(alpha_j - alpha_i)
       !endif
       if ((istrapped .eq. 1).and.(nbracket.eq.0)) then ! try quadratic 
           call pickAlphaQuadratic(left,right,alpha_i,f_i,&
               g_i,alpha_j,f_j,alphaNext)
           alphaPrev = alpha
           alpha = alphaNext
           ! compute the penalty functional
           call linComb(ONE,m_0,alpha,h,m)
           call func2(lambda,d,m0,m,f,mNorm,dHat,eAll,rms)
           call printf('QUADLS',lambda,alpha,f,mNorm,rms)
           call printf('QUADLS',lambda,alpha,f,mNorm,rms,logFile)
       else if ((istrapped .eq. 2).and.(nbracket.eq.0)) then ! try quadratic(secant)
           call pickAlphaSecant(left,right,alpha_i,f_i,&
               g_i,alpha_j,g_j,alphaNext)
           alphaPrev = alpha
           alpha = alphaNext
           ! compute the penalty functional
           call linComb(ONE,m_0,alpha,h,m)
           call func2(lambda,d,m0,m,f,mNorm,dHat,eAll,rms)
           call printf('SCNTLS',lambda,alpha,f,mNorm,rms)
           call printf('SCNTLS',lambda,alpha,f,mNorm,rms,logFile)
       else if ((istrapped .eq. -1).and.(nbracket.eq.0)) then ! oh, jump 
           call pickAlphaSecant(left,1.0D+2,alpha_i,f_i,&
               g_i,alpha_j,g_j,alphaNext)
           alphaPrev = alpha
           alpha = alphaNext
           ! compute the penalty functional
           call linComb(ONE,m_0,alpha,h,m)
           call func(lambda,d,m0,m,f,mNorm,dHat,eAll,rms)
           call printf('QUADLS',lambda,alpha,f,mNorm,rms)
           call printf('QUADLS',lambda,alpha,f,mNorm,rms,logFile)
       else ! cubic 
           call pickAlphaCubic(left,right,alpha_i,f_i,&
               g_i,alpha_j,f_j,g_j,alphaNext)
           alphaPrev = alpha
           alpha = alphaNext
           ! compute the penalty functional
           call linComb(ONE,m_0,alpha,h,m)
           call func2(lambda,d,m0,m,f,mNorm,dHat,eAll,rms)
           call printf('CUBICLS',lambda,alpha,f,mNorm,rms)
           call printf('CUBICLS',lambda,alpha,f,mNorm,rms,logFile)
       endif
       niter = niter + 1
       nbracket = nbracket + 1
       ! firstly store the previous values
       fPrev = f
       gPrev = g
       !calculatie gradient to test the Wolfe condition
       call gradient2(lambda,d,m0,m,grad,dHat,eAll)
       g = dotProd(grad, h) 
       if (f.lt.f_1) then
           alpha_1 = alpha
           dHat_1 = dHat
           eAll_1 = eAll
           m_1 = m
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
       if (nbracket .ge. 2) then
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
       m = m_1
       rms = rms_1
       f= f_1
       grad = grad_1
   endif
   call deall_dataVectorMTX(dHat_1)
   call deall_modelParam(m_0)
   call deall_modelParam(m_1)
   call deall_modelParam(grad_1)
   call deall_solnVectorMTX(eAll_1)
  end subroutine lineSearchWolfe

  !**********************************************************************
  subroutine lineSearchArmijo(lambda,d,m0,h,alpha,m,f,grad, &
   & rms,niter,dHat,eAll,flag,logid)
 
   ! Note: inexact line searches ultimately fit into two catalogs - 
   ! i.e. Armijo–Goldstein (back-tracking) and Wolfe conditions
   ! the latter also contains a few different criterias: (1) "Armijo 
   ! rule", and (2) "curvature condition", the latter requires an extra 
   ! gradient calculation at x_k + alpha_k*d_k
   ! (1) and modified (2) will form the "strong Wolfe" condition which is 
   ! the line search "convergence criteria" here

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
   !    cubic interpolation and sectioning with my variant of More-Thuente
   !    scheme
   !    TODO: I should make a more clear version of this, and stuff it 
   !          into NLCG also (not sure if Anna will like that)
   !
   ! following Anna's idea, the major intention here is to use a cheap line
   ! search scheme (with merely 2 forward-like-calculations) to quickly skip 
   ! to a small overall penalty function level and only to start bracketing
   ! when the quadratic interpolation doesn't work, in which case cubic 
   ! probably won't work either...
   ! 
   implicit none
   real(kind=prec), intent(in)               :: lambda ! lagrange multiplier
   type(dataVectorMTX_t), intent(in)         :: d
   type(modelParam_t), intent(in)            :: m0 ! current model
   type(modelParam_t), intent(in)            :: h  ! search direction
   real(kind=prec), intent(inout)            :: alpha ! step size
   type(modelParam_t), intent(inout)         :: m 
   real(kind=prec), intent(inout)            :: f  ! penalty function
   type(modelParam_t), intent(inout)         :: grad ! function gradient
   real(kind=prec), intent(out)              :: rms
   integer, intent(out)                      :: niter
   type(dataVectorMTX_t), intent(out)        :: dHat
   type(solnVectorMTX_t), intent(inout)      :: eAll
   integer, intent(inout), optional          :: flag
   character(100),intent(in), optional       :: logid

   ! local variables
   real(kind=prec)                 :: alpha_0, alpha_1,mNorm,alphaNext
   real(kind=prec)                 :: alpha_i,alpha_j ! brackets
   real(kind=prec)                 :: f_i,f_j ! function value at brackets
   real(kind=prec)                 :: g_i,g_j ! grad at brackets
   real(kind=prec)                 :: alphaPrev,fPrev,gPrev ! previous one
   real(kind=prec)                 :: smin,smax ! the min/max step lengths
   real(kind=prec)                 :: left, right ! bounds for bracket
   logical                         :: starting_guess
   integer                         :: istrapped,nbracket
   real(kind=prec)                 :: eps,k,c,c2,a,b,q1,q2,q3
   real(kind=prec)                 :: g_0,g_1,g,f_0,f_1,rms_1
   type(modelParam_t)              :: m_0,m_1,grad_1 ! initial grad
   type(dataVectorMTX_t)           :: dHat_1
   type(solnVectorMTX_t)           :: eAll_1
   character(100)                  :: logFile

   ! parameters 
   c = iterControl%c ! ensures sufficient decrease
   ! c2 = iterControl%c2 ! ensures culvature condition

   if (.not.present(logid)) then
       logFile = trim(iterControl%fname)//'_NLCG.log'
   else 
       logFile = logid
   endif

   ! initialize the line search
   niter = 0
   m_0 = m
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
   ! alpha_0 is the initial step size, which is set in LBFGS
   alpha_1 = alpha
   alpha_0 = alpha

   ! compute the trial m, f, dHat, eAll, rms
   m_1 = m_0
   ! m_1 = m_0 + dir*step
   call linComb(ONE,m_0,alpha_1,h,m_1)
   call func2(lambda,d,m0,m_1,f_1,mNorm,dHat_1,eAll_1,rms_1)
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

   call gradient2(lambda,d,m0,m_1,grad_1,dHat_1,eAll_1)
   g_1 = dotProd(grad_1, h)
   write(*,'(a29,es12.5)',advance='no') '    GRAD: computed, with g0=',g_0
   write(*,'(a4,es12.5)') ' g1=',g_1
   write(ioLog,'(a29,es12.5)',advance='no') '    GRAD: computed, with g0=',g_0
   write(ioLog,'(a4,es12.5)') ' g1=',g_1


   ! ====================================================================== !
   ! test if the initial guess satisfies the Armijo condition
   ! ====================================================================== !
   if ((f_1.lt.f_0 + c * alpha_1 * g_0)) then 
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
           if (f <= f_0 + c * alpha * g_0) then 
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
               ! we are already far from the f_0 and g_0 
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
           call linComb(ONE,m_0,alpha,h,m)
           call func2(lambda,d,m0,m,f,mNorm,dHat,eAll,rms)
           call printf('CUBICLS',lambda,alpha,f,mNorm,rms)
           call printf('CUBICLS',lambda,alpha,f,mNorm,rms,logFile)
           niter = niter + 1
           nbracket = nbracket +1
           call gradient2(lambda,d,m0,m,grad,dHat,eAll)
           g = dotProd(grad, h) 
           if (f.lt.f_1) then
               alpha_1 = alpha
               dHat_1 = dHat
               eAll_1 = eAll
               m_1 = m
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
       !if (f_j .gt. 1.5*f_i) then
       !    right = alpha_j - 0.618*(alpha_j - alpha_i)
       !else
       right = alpha_j - 0.1*(alpha_j - alpha_i)
       !endif
       if ((istrapped .eq. 1).and.(nbracket.eq.0)) then ! try quadratic 
           call pickAlphaQuadratic(left,right,alpha_i,f_i,&
               g_i,alpha_j,f_j,alphaNext)
           alphaPrev = alpha
           alpha = alphaNext
           ! compute the penalty functional
           call linComb(ONE,m_0,alpha,h,m)
           call func2(lambda,d,m0,m,f,mNorm,dHat,eAll,rms)
           call printf('QUADLS',lambda,alpha,f,mNorm,rms)
           call printf('QUADLS',lambda,alpha,f,mNorm,rms,logFile)
       else if ((istrapped .eq. 2).and.(nbracket.eq.0)) then ! try quadratic 
           call pickAlphaSecant(left,right,alpha_i,f_i,&
               g_i,alpha_j,g_j,alphaNext)
           alphaPrev = alpha
           alpha = alphaNext
           ! compute the penalty functional
           call linComb(ONE,m_0,alpha,h,m)
           call func2(lambda,d,m0,m,f,mNorm,dHat,eAll,rms)
           call printf('QUADLS',lambda,alpha,f,mNorm,rms)
           call printf('QUADLS',lambda,alpha,f,mNorm,rms,logFile)
       else if ((istrapped .eq. -1).and.(nbracket.eq.0)) then ! jump 
           call pickAlphaSecant(left,1.0D+2,alpha_i,f_i,&
               g_i,alpha_j,g_j,alphaNext)
           alphaPrev = alpha
           alpha = alphaNext
           ! compute the penalty functional
           call linComb(ONE,m_0,alpha,h,m)
           call func2(lambda,d,m0,m,f,mNorm,dHat,eAll,rms)
           call printf('QUADLS',lambda,alpha,f,mNorm,rms)
           call printf('QUADLS',lambda,alpha,f,mNorm,rms,logFile)
       else ! cubic 
           call pickAlphaCubic(left,right,alpha_i,f_i,&
               g_i,alpha_j,f_j,g_j,alphaNext)
           alphaPrev = alpha
           alpha = alphaNext
           ! compute the penalty functional
           call linComb(ONE,m_0,alpha,h,m)
           call func2(lambda,d,m0,m,f,mNorm,dHat,eAll,rms)
           call printf('CUBICLS',lambda,alpha,f,mNorm,rms)
           call printf('CUBICLS',lambda,alpha,f,mNorm,rms,logFile)
       endif
       niter = niter + 1
       nbracket = nbracket + 1
       ! firstly store the previous values
       fPrev = f
       gPrev = g
       !calculatie gradient to test the Wolfe condition
       call gradient2(lambda,d,m0,m,grad,dHat,eAll)
       g = dotProd(grad, h) 
       if (f.lt.f_1) then
           alpha_1 = alpha
           dHat_1 = dHat
           eAll_1 = eAll
           m_1 = m
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
       if (f <= f_0 + c * alpha * g_0) then 
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
       write(*,'(a47)') 'Armijo Condition satisfied, exiting line search'
       write(ioLog,'(a47)') 'Armijo Condition satisfied, exiting line search'
   else 
       if (present(flag)) then
           flag = -1 ! Wolfe condition not satisfied, need to restart
       endif
       write(*,'(a40)') 'Armijo Condition NOT satisfied, abort...'
       write(ioLog,'(a40)') 'Armijo Condition NOT satisfied, abort...'
   endif

   if (starting_guess) then
       if (istrapped.ne.0) then
           write(6, *) 'recalling the best model so far...'
           write(ioLog, *) 'recalling the best model so far...'
       endif
       alpha = alpha_1
       dHat = dHat_1
       eAll = eAll_1
       m = m_1
       rms = rms_1
       f= f_1
       grad = grad_1
   endif
   call deall_dataVectorMTX(dHat_1)
   call deall_modelParam(m_0)
   call deall_modelParam(m_1)
   call deall_modelParam(grad_1)
   call deall_solnVectorMTX(eAll_1)

  end subroutine lineSearchArmijo

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
    
    !write(6,*) 'xi = ', xi, 'xj = ', xj, 'x1 = ', x1, 'x2 = ', x2 
    !write(ioLog,*) 'xi = ', xi, 'xj = ', xj, 'x1 = ', x1, 'x2 = ', x2 
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
    ! write(6,*) 'sleft = ', left, 's1 = ', real(s1), 's2 = ', real(s2), 'sright = ', right
    ! write(ioLog,*) 'sleft = ', left, 's1 = ', real(s1), 's2 = ', real(s2), 'sright = ', right
    ! write(6,*) 'fleft = ', fleft, 'fs1 = ', fs1, 'fs2 = ', fs2, 'fright = ', fright
    ! write(ioLog,*) 'fleft = ', fleft, 'fs1 = ', fs1, 'fs2 = ', fs2, 'fright = ', fright
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
    

  ! write(ioLog,*) 'x1= ', x1, 'x2= ', x2, 'f1= ', f1, 'f2= ', f2, 'g1= ',g1 
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
    ! write(6,*) 'left = ', left, 's = ', s, 'right = ', right
    ! write(ioLog,*) 'left = ', left, 's = ', s, 'right = ', right
    ! write(6,*) 'fleft = ', fleft, 'fs = ', fs, 'fright = ', fright
    ! write(ioLog,*) 'fleft = ', fleft, 'fs = ', fs, 'fright = ', fright
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
    

  ! write(ioLog,*) 'x1= ', x1, 'x2= ', x2, 'f1= ', f1, 'f2= ', f2, 'g1= ',g1 
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
    ! write(6,*) 'left = ', left, 's = ', s, 'right = ', right
    ! write(ioLog,*) 'left = ', left, 's = ', s, 'right = ', right
    ! write(6,*) 'fleft = ', fleft, 'fs = ', fs, 'fright = ', fright
    ! write(ioLog,*) 'fleft = ', fleft, 'fs = ', fs, 'fright = ', fright
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

end module LBFGS
