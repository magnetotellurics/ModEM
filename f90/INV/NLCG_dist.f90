module NLCG_dist

  use ModEM_timers
  use invcore
  use DataSpace, only: countData
  use utilities

  implicit none

  public :: NLCGsolver_dist

  type :: NLCGDistControl_t
     integer            :: maxIter
     real(kind=prec)    :: rmsTol
     real(kind=prec)    :: fdiffTol
     real(kind=prec)    :: lambda
     real(kind=prec)    :: lambdaTol
     real(kind=prec)    :: k
     real(kind=prec)    :: c
     real(kind=prec)    :: c2
     integer            :: nCGmax
     real(kind=prec)    :: alpha_1
     real(kind=prec)    :: startdm
     real(kind=prec)    :: gamma
     character(80)      :: fname
     integer            :: nskip
     ! Distortion-specific control parameters
     real(kind=prec)    :: nu
     integer            :: nDistIter
     real(kind=prec)    :: alpha_C
     integer            :: nSigmaIter ! not used (yet)
     character(80)      :: distFname
     character(80)      :: init_distFname
  end type NLCGDistControl_t

  type(NLCGDistControl_t), private, save :: distControl

Contains

  subroutine set_NLCGDistControl(ctrl)
    type(NLCGDistControl_t), intent(inout) :: ctrl
    ctrl%maxIter = 300
    ctrl%rmsTol  = 1.05
    ctrl%fdiffTol = 2.0e-3
    ctrl%lambda = 1.
    ctrl%lambdaTol = 1.0e-4
    ctrl%k = 10.
    ctrl%c = 1.0e-4
    ctrl%c2 = 0.9
    ctrl%nCGmax = 8
    ctrl%alpha_1 = 20.
    ctrl%startdm = 20.
    ctrl%gamma = 0.99
    ctrl%fname = 'Modular'
    ctrl%nskip = 1
    ctrl%nu = 10.0
    ctrl%nDistIter = 1
    ctrl%alpha_C = 0.01
    ctrl%nSigmaIter = 5
    ctrl%distFname = 'distortion.dat'
    ctrl%init_distFname = ''
  end subroutine set_NLCGDistControl

  subroutine read_NLCGDistControl(ctrl, rFile, fileExists)
    type(NLCGDistControl_t), intent(inout) :: ctrl
    character(*), intent(in) :: rFile
    logical, intent(out), optional :: fileExists
    integer :: ios
    logical :: exists
    character(80) :: string

    call set_NLCGDistControl(ctrl)

    inquire(FILE=rFile, EXIST=exists)
    if (present(fileExists)) fileExists = exists
    if (.not. exists) return

    write(*,*) 'Reading inverse configuration from file ', trim(rFile)
    open(unit=ioInvCtrl, file=rFile, status='old', iostat=ios)
    if (ios /= 0) write(0,*) 'Error opening file: ', rFile

    read(ioInvCtrl, '(a36,a80)') string, ctrl%fname
    if (output_level > 2) write(*,'(a36,a80)') string, ctrl%fname
    ctrl%fname = adjustl(ctrl%fname)
    read(ioInvCtrl, '(a36,g15.7)') string, ctrl%lambda
    if (output_level > 2) write(*,'(a36,g15.7)') string, ctrl%lambda
    read(ioInvCtrl, '(a36,g15.7)') string, ctrl%k
    if (output_level > 2) write(*,'(a36,g15.7)') string, ctrl%k
    read(ioInvCtrl, '(a36,g15.7)') string, ctrl%startdm
    if (output_level > 2) write(*,'(a36,g15.7)') string, ctrl%startdm
    read(ioInvCtrl, '(a36,g15.7)') string, ctrl%fdiffTol
    if (output_level > 2) write(*,'(a36,g15.7)') string, ctrl%fdiffTol
    read(ioInvCtrl, '(a36,g15.7)') string, ctrl%rmsTol
    if (output_level > 2) write(*,'(a36,g15.7)') string, ctrl%rmsTol
    read(ioInvCtrl, '(a36,g15.7)') string, ctrl%lambdaTol
    if (output_level > 2) write(*,'(a36,g15.7)') string, ctrl%lambdaTol
    read(ioInvCtrl, '(a36,i4)') string, ctrl%maxIter
    if (output_level > 2) write(*,'(a36,i4)') string, ctrl%maxIter
    ! Distortion-specific parameters
    read(ioInvCtrl, '(a36,g15.7)', iostat=ios) string, ctrl%nu
    if (output_level > 2) write(*,'(a36,g15.7)') string, ctrl%nu
    read(ioInvCtrl, '(a36,i4)', iostat=ios) string, ctrl%nSigmaIter
    if (output_level > 2) write(*,'(a36,i4)') string, ctrl%nSigmaIter
    read(ioInvCtrl, '(a36,i4)', iostat=ios) string, ctrl%nDistIter
    if (output_level > 2) write(*,'(a36,i4)') string, ctrl%nDistIter
    read(ioInvCtrl, '(a36,g15.7)', iostat=ios) string, ctrl%alpha_C
    if (output_level > 2) write(*,'(a36,g15.7)') string, ctrl%alpha_C
    read(ioInvCtrl, '(a36,a80)', iostat=ios) string, ctrl%distFname
    if (output_level > 2) write(*,'(a36,a80)') string, ctrl%distFname
    ctrl%distFname = adjustl(ctrl%distFname)
    read(ioInvCtrl, '(a36,a80)', iostat=ios) string, ctrl%init_distFname
    if (output_level > 2) write(*,'(a36,a80)') string, ctrl%init_distFname
    ctrl%init_distFname = adjustl(ctrl%init_distFname)

    close(ioInvCtrl)
  end subroutine read_NLCGDistControl

  subroutine NLCGsolver_dist(d, lambda, nu, m0, m, distC, fname)
    type(dataVectorMTX_t), intent(inout) :: d
    real(kind=prec), intent(inout) :: lambda
    real(kind=prec), intent(in) :: nu
    type(modelParam_t), intent(in) :: m0
    type(modelParam_t), intent(inout) :: m
    type(distortionParam_t), intent(inout) :: distC
    character(*), intent(in), optional :: fname

    type(dataVectorMTX_t) :: dHat, res
    type(modelParam_t) :: mHat, m_minus_m0
    type(modelParam_t) :: gradM, g, h, gPrev
    type(distortionParam_t) :: gradC, gC, hC, gPrevC, distC_work
    type(solnVectorMTX_t) :: eAll
    real(kind=prec) :: value, valuePrev, rms, rmsPrev
    real(kind=prec) :: alpha, beta, gnorm, mNorm, distNorm, Ndata, Nmodel
    real(kind=prec) :: grad_dot_h, g_dot_g, g_dot_gPrev, gPrev_dot_gPrev
    real(kind=prec) :: g_dot_h
    integer :: iter, nCG, nLS, nfunc, ios, i, j
    logical :: ok
    character(3) :: iterChar
    character(100) :: mFile, mHatFile, gradFile, dataFile, resFile, logFile
    character(100) :: distFile

    if (present(fname)) then
       call read_NLCGDistControl(distControl, fname, ok)
       if (ok) lambda = distControl%lambda
    else
       call set_NLCGDistControl(distControl)
    end if
    distControl%nu = nu

    logFile = trim(distControl%fname)//'_NLCG_DIST.log'
    open(unit=ioLog, file=logFile, status='unknown', position='append', iostat=ios)

    alpha = distControl%alpha_1

    write(*,'(a41,es8.1)') 'The initial damping parameter lambda is ',lambda
    write(ioLog,'(a41,es8.1)') 'The initial damping parameter lambda is ',lambda
    write(*,'(a44,es8.1)') 'The distortion regularization parameter nu is ',distControl%nu
    write(ioLog,'(a44,es8.1)') 'The distortion regularization parameter nu is ',distControl%nu

    mHat = m

    ! Initialize distortion parameters
    if (distControl%init_distFname /= '' .and. distControl%init_distFname /= 'n') then
       call read_distortionParam(distControl%init_distFname, distC)
       write(*,*) 'Read initial distortion from ', trim(distControl%init_distFname)
    else
       call identity_distortionParam(distC)
    end if

    ! Create distortion gradient and work variables
    call create_distortionParam(distC%nSites, gradC)
    call create_distortionParam(distC%nSites, gC)
    call create_distortionParam(distC%nSites, hC)
    call create_distortionParam(distC%nSites, gPrevC)
    call create_distortionParam(distC%nSites, distC_work)

    ! First forward pass and function evaluation
    call func_dist(lambda, distControl%nu, d, m0, mHat, distC, value, mNorm, distNorm, dHat, eAll, rms)
    call printf('START', lambda, alpha, value, mNorm, rms)
    call printf('START', lambda, alpha, value, mNorm, rms, logFile)
    nfunc = 1

    write(iterChar, '(i3.3)') 0
    ! output (smoothed) initial model and responses for later reference
    call CmSqrtMult(mHat,m_minus_m0)
    call linComb(ONE,m_minus_m0,ONE,m0,m)
    if (output_level > 1) then
      mFile = trim(distControl%fname)//'_NLCG_DIST_'//iterChar//'.rho'
      call write_modelParam(m,trim(mFile))
    end if
    if (output_level > 2) then
      dataFile = trim(distControl%fname)//'_NLCG_DIST_'//iterChar//'.dat'
      call write_dataVectorMTX(dHat,trim(dataFile))
    end if

    ! First gradient (both model and distortion)
    call gradient_dist(lambda, distControl%nu, d, m0, mHat, distC, gradM, gradC, dHat, eAll)

    gnorm = sqrt(dotProd(gradM, gradM))
    write(*,'(a42,es12.5)') 'GRAD: initial norm of model gradient is', gnorm
    write(ioLog,'(a42,es12.5)') 'GRAD: initial norm of model gradient is', gnorm
    if (gnorm < TOL6) then
       call errStop('Problem with gradient computations: first gradient is zero')
    else
       alpha = distControl%startdm / gnorm
       write(*,'(a39,es12.5)') 'The initial value of alpha updated to ', alpha
       write(ioLog,'(a39,es12.5)') 'The initial value of alpha updated to ', alpha
    end if

    nCG = 0
    iter = 0
    g = gradM
    call linComb(MinusONE, gradM, R_ZERO, gradM, g)
    h = g

    call ModEM_timers_create("NLCG_DIST Iteration", .false.)

    ! Main alternating loop
    do
       call ModEM_timers_start("NLCG_DIST Iteration")

       if ((rms < distControl%rmsTol) .or. (iter >= distControl%maxIter)) exit

       iter = iter + 1

       rmsPrev = rms
       valuePrev = value
       grad_dot_h = dotProd(gradM, h)

       ! σ-update (outer NLCG sweep)
       write(*,'(a30)') 'Part 1: Model Sigma update '
       write(ioLog,'(a30)') 'Part 1: Model Sigma update '

       ! Line search for conductivity
       call lineSearchDistCubic(lambda, distControl%nu, d, m0, mHat, distC, &
            h, alpha, value, gradM, rms, nLS, dHat, eAll)
       nfunc = nfunc + nLS

       ! Update alpha for next sigma line search
       alpha = 2.0_prec * (value - valuePrev) / grad_dot_h
       alpha = (ONE + 0.01_prec) * alpha

       !write(*,'(a25,i5)') 'Completed sigma update, iter ', iter
       !write(ioLog,'(a25,i5)') 'Completed sigma update, iter ', iter

       ! C-update (gradient descent steps)
       write(*,'(a35)') '   Part 2: Distortion Matrix update'
       write(ioLog,'(a35)') '   Part 2: Distortion Matrix update'

       ! Re-compute the smoothed model m from updated mHat
       call CmSqrtMult(mHat, m_minus_m0)
       call linComb(ONE, m_minus_m0, ONE, m0, m)

       call copy_distortionParam(distC, distC_work)
       Ndata = countData(d)
       do i = 1, distControl%nDistIter
          ! Simple Gradient descent: C_j ← C_j - alpha_C * gradC_j
          ! should not be too hard as the distorsion matrices are 2x2xnSites
          call compute_distortion_gradient(d, eAll, m, distC_work, distControl%nu, Ndata, gradC)
          do j = 1, distC_work%nSites
             distC_work%C(1,1,j) = distC_work%C(1,1,j) - distControl%alpha_C * gradC%C(1,1,j)
             distC_work%C(1,2,j) = distC_work%C(1,2,j) - distControl%alpha_C * gradC%C(1,2,j)
             distC_work%C(2,1,j) = distC_work%C(2,1,j) - distControl%alpha_C * gradC%C(2,1,j)
             distC_work%C(2,2,j) = distC_work%C(2,2,j) - distControl%alpha_C * gradC%C(2,2,j)
          end do
       end do
       distC = distC_work

       ! Re-evaluate function value after C-update
       call func_dist(lambda, distControl%nu, d, m0, mHat, distC, value, mNorm, distNorm, dHat, eAll, rms)
       nfunc = nfunc + 1

      ! Recompute sigma gradient at the updated (mHat, distC) state.
      gPrev = g
      call gradient_dist(lambda, distControl%nu, d, m0, mHat, distC, gradM, gradC, dHat, eAll)
      call linComb(MinusONE, gradM, R_ZERO, gradM, g)

       call ModEM_timers_stop("NLCG_DIST Iteration", .false.)
       call ModEM_timers_print("NLCG_DIST Iteration", ioLog)


       write(*,'(a30,i5)') 'Completed NLCG_DIST iteration ',iter
       write(ioLog,'(a30,i5)') 'Completed NLCG_DIST iteration ',iter

       Nmodel = countModelParam(mHat)
       mNorm = dotProd(mHat, mHat) / Nmodel
       call printf('DIST', lambda, alpha, value, mNorm, rms)
       call printf('DIST', lambda, alpha, value, mNorm, rms, logFile)

       ! Write intermediate output
       if (mod(iter, distControl%nskip) == 0) then
          call CmSqrtMult(mHat, m_minus_m0)
          call linComb(ONE, m_minus_m0, ONE, m0, m)
          write(iterChar, '(i3.3)') iter
          if (output_level > 1) then
             mFile = trim(distControl%fname)//'_NLCG_DIST_'//iterChar//'.rho'
             call write_modelParam(m, trim(mFile))
          end if
          if (output_level > 2) then
             mHatFile = trim(distControl%fname)//'_NLCG_DIST_'//iterChar//'.prm'
             call write_modelParam(mHat, trim(mHatFile))
          end if
          if (output_level > 2) then
             dataFile = trim(distControl%fname)//'_NLCG_DIST_'//iterChar//'.dat'
             call write_dataVectorMTX(dHat, trim(dataFile))
          end if
          if (output_level > 2) then
             res = d
             call linComb(ONE, d, MinusONE, dHat, res)
             resFile = trim(distControl%fname)//'_NLCG_DIST_'//iterChar//'.res'
             call write_dataVectorMTX(res, trim(resFile))
             call deall_dataVectorMTX(res)
          end if
          ! Write distortion matrices
          distFile = trim(distControl%fname)//'_'//trim(distControl%distFname)
          call write_distortionParam(trim(distFile), distC)
       end if

       ! Check for stall: update lambda
       if (abs(rmsPrev - rms) < distControl%fdiffTol) then
          lambda = lambda / distControl%k
          if (lambda < distControl%lambdaTol) then
             write(*,'(a55)') 'Unable to get out of a local minimum. Exiting...'
             write(ioLog,'(a55)') 'Unable to get out of a local minimum. Exiting...'
             exit
          end if
          gnorm = sqrt(dotProd(gradM, gradM))
          alpha = min(ONE, distControl%startdm) / gnorm
          call linComb(MinusONE, gradM, R_ZERO, gradM, g)
          write(*,'(a55)') 'Restarting NLCG_DIST with the damping parameter updated'
          call printf('to', lambda, alpha, value, mNorm, rms)
          write(ioLog,'(a55)') 'Restarting NLCG_DIST with the damping parameter updated'
          call printf('to', lambda, alpha, value, mNorm, rms, logFile)
          h = g
          nCG = 0
          cycle
       end if

       ! CG beta (Polak-Ribiere) for model gradient
       g_dot_g = dotProd(g, g)
       g_dot_gPrev = dotProd(g, gPrev)
       gPrev_dot_gPrev = dotProd(gPrev, gPrev)
       g_dot_h = dotProd(g, h)

       beta = (g_dot_g - g_dot_gPrev) / gPrev_dot_gPrev

       if ((g_dot_g + beta * g_dot_h <= R_ZERO) .and. (nCG >= distControl%nCGmax)) then
          write(*,'(a45)') 'Restarting NLCG to restore orthogonality'
          write(ioLog,'(a45)') 'Restarting NLCG to restore orthogonality'
          nCG = 0
          beta = R_ZERO
       else
          nCG = nCG + 1
       end if
       call linComb(ONE, g, beta, h, h)
    end do

    ! Final output
    call CmSqrtMult(mHat, m_minus_m0)
    call linComb(ONE, m_minus_m0, ONE, m0, m)
    d = dHat
    write(*,'(a25,i5,a25,i5)') 'NLCG_DIST iterations:', iter, ' function evaluations:', nfunc
    write(ioLog,'(a25,i5,a25,i5)') 'NLCG_DIST iterations:', iter, ' function evaluations:', nfunc
    close(ioLog, iostat=ios)

    ! Write final distortion
    distFile = trim(distControl%fname)//'_'//trim(distControl%distFname)
    call write_distortionParam(trim(distFile), distC)

    ! Cleanup
    call deall_dataVectorMTX(dHat)
    call deall_dataVectorMTX(res)
    call deall_modelParam(mHat)
    call deall_modelParam(m_minus_m0)
    call deall_modelParam(gradM)
    call deall_modelParam(g)
    call deall_modelParam(h)
    call deall_modelParam(gPrev)
    call deall_solnVectorMTX(eAll)
    call deall_distortionParam(gradC)
    call deall_distortionParam(gC)
    call deall_distortionParam(hC)
    call deall_distortionParam(gPrevC)
    call deall_distortionParam(distC_work)
    call ModEM_timers_destory('NLCG_DIST Iteration')

  end subroutine NLCGsolver_dist

  !**********************************************************************
  subroutine lineSearchDistCubic(lambda, nu, d, m0, mHat, distC, h, alpha, f, grad, &
       rms, niter, dHat, eAll, gamma)

    real(kind=prec), intent(in) :: lambda, nu
    type(dataVectorMTX_t), intent(in) :: d
    type(modelParam_t), intent(in) :: m0
    type(modelParam_t), intent(inout) :: mHat
    type(distortionParam_t), intent(in) :: distC
    type(modelParam_t), intent(in) :: h
    real(kind=prec), intent(inout) :: alpha
    real(kind=prec), intent(inout) :: f
    type(modelParam_t), intent(inout) :: grad
    real(kind=prec), intent(out) :: rms
    integer, intent(out) :: niter
    type(dataVectorMTX_t), intent(out) :: dHat
    type(solnVectorMTX_t), intent(inout) :: eAll
    real(kind=prec), intent(in), optional :: gamma

    real(kind=prec) :: alpha_1, alpha_i, alpha_j, mNorm, distNorm
    logical :: starting_guess, relaxation
    real(kind=prec) :: eps, k, c, a, b, q1, q2, q3
    real(kind=prec) :: g_0, f_0, f_1, f_i, f_j, rms_1, mNorm_1, distNorm_1
    type(modelParam_t) :: mHat_0, mHat_1
    type(distortionParam_t) :: gradC_dummy
    type(dataVectorMTX_t) :: dHat_1
    type(solnVectorMTX_t) :: eAll_1
    character(100) :: logFile

    c = distControl%c
    logFile = trim(distControl%fname)//'_NLCG_DIST.log'

    niter = 0
    mHat_0 = mHat
    f_0 = f
    starting_guess = .false.

    g_0 = dotProd(grad, h)
    alpha_1 = alpha

    if (present(gamma)) then
       relaxation = .true.
    else
       relaxation = .false.
    end if

    mHat_1 = mHat_0
    call linComb(ONE, mHat_0, alpha_1, h, mHat_1)
    call func_dist(lambda, nu, d, m0, mHat_1, distC, f_1, mNorm_1, distNorm_1, dHat_1, eAll_1, rms_1)
    call printf('STARTLS', lambda, alpha, f_1, mNorm_1, rms_1)
    call printf('STARTLS', lambda, alpha, f_1, mNorm_1, rms_1, logFile)
    niter = niter + 1

    if (f_1 - f_0 >= LARGE_REAL) then
       call errStop("Try a smaller starting value of alpha ('Initial search step in model units' in InvCtrl file). Exiting..")
    end if

    a = (f_1 - f_0 - g_0*alpha_1) / (alpha_1**2)
    b = g_0
    if (a < 0) then
       starting_guess = .true.
       alpha = alpha_1
       dHat = dHat_1
       eAll = eAll_1
       mHat = mHat_1
       rms = rms_1
       f = f_1
       if (relaxation) then
          call linComb(ONE, mHat_0, gamma*alpha, h, mHat)
          call func_dist(lambda, nu, d, m0, mHat, distC, f, mNorm, distNorm, dHat, eAll, rms)
          call printf('RELAX', lambda, gamma*alpha, f, mNorm, rms)
          call printf('RELAX', lambda, gamma*alpha, f, mNorm, rms, logFile)
       end if
       call create_distortionParam(distC%nSites, gradC_dummy)
       call gradient_dist(lambda, nu, d, m0, mHat, distC, grad, gradC_dummy, dHat, eAll)
       call deall_distortionParam(gradC_dummy)
       write(*,'(a45)') 'Quadratic has no minimum, exiting line search'
       write(ioLog,'(a45)') 'Quadratic has no minimum, exiting line search'
       call deall_dataVectorMTX(dHat_1)
       call deall_modelParam(mHat_0)
       call deall_modelParam(mHat_1)
       call deall_solnVectorMTX(eAll_1)
       return
    end if

    alpha = -b / (TWO*a)
    call linComb(ONE, mHat_0, alpha, h, mHat)
    call func_dist(lambda, nu, d, m0, mHat, distC, f, mNorm, distNorm, dHat, eAll, rms)
    call printf('QUADLS', lambda, alpha, f, mNorm, rms)
    call printf('QUADLS', lambda, alpha, f, mNorm, rms, logFile)
    niter = niter + 1

    if (f < f_0 + c * alpha * g_0) then
       if (f_1 < f) then
          starting_guess = .true.
          alpha = alpha_1
          dHat = dHat_1
          eAll = eAll_1
          mHat = mHat_1
          rms = rms_1
          f = f_1
       end if
       if (relaxation) then
          call linComb(ONE, mHat_0, gamma*alpha, h, mHat)
          call func_dist(lambda, nu, d, m0, mHat, distC, f, mNorm, distNorm, dHat, eAll, rms)
          call printf('RELAX', lambda, gamma*alpha, f, mNorm, rms)
          call printf('RELAX', lambda, gamma*alpha, f, mNorm, rms, logFile)
       end if
       call create_distortionParam(distC%nSites, gradC_dummy)
       call gradient_dist(lambda, nu, d, m0, mHat, distC, grad, gradC_dummy, dHat, eAll)
       call deall_distortionParam(gradC_dummy)
       write(*,'(a60)') 'Sufficient decrease condition satisfied, exiting line search'
       write(ioLog,'(a60)') 'Sufficient decrease condition satisfied, exiting line search'
       call deall_dataVectorMTX(dHat_1)
       call deall_modelParam(mHat_0)
       call deall_modelParam(mHat_1)
       call deall_solnVectorMTX(eAll_1)
       return
    end if

    if (f > f_0) then
       write(*,'(a75)') 'Unable to fit a quadratic due to bad gradient estimate, exiting line search'
       write(ioLog,'(a75)') 'Unable to fit a quadratic due to bad gradient estimate, exiting line search'
    else
       alpha_i = alpha_1; f_i = f_1
       alpha_j = alpha;   f_j = f
       fit_cubic: do
          q1 = f_i - f_0 - g_0 * alpha_i
          q2 = f_j - f_0 - g_0 * alpha_j
          q3 = alpha_i**2 * alpha_j**2 * (alpha_j - alpha_i)
          a = (alpha_i**2 * q2 - alpha_j**2 * q1) / q3
          b = (alpha_j**3 * q1 - alpha_i**3 * q2) / q3
          alpha = (-b + sqrt(b*b - 3*a*g_0)) / (3*a)
          call linComb(ONE, mHat_0, alpha, h, mHat)
          call func_dist(lambda, nu, d, m0, mHat, distC, f, mNorm, distNorm, dHat, eAll, rms)
          call printf('CUBICLS', lambda, alpha, f, mNorm, rms)
          call printf('CUBICLS', lambda, alpha, f, mNorm, rms, logFile)
          niter = niter + 1
          if (f < f_0 + c * alpha * g_0) exit
          alpha_i = alpha_j; f_i = f_j
          alpha_j = alpha;   f_j = f
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
    if (starting_guess) then
       alpha = alpha_1
       dHat = dHat_1
       eAll = eAll_1
       mHat = mHat_1
       rms = rms_1
       f = f_1
    end if
    if (relaxation) then
       call linComb(ONE, mHat_0, gamma*alpha, h, mHat)
       call func_dist(lambda, nu, d, m0, mHat, distC, f, mNorm, distNorm, dHat, eAll, rms)
       call printf('RELAX', lambda, gamma*alpha, f, mNorm, rms)
       call printf('RELAX', lambda, gamma*alpha, f, mNorm, rms, logFile)
    end if
    call create_distortionParam(distC%nSites, gradC_dummy)
    call gradient_dist(lambda, nu, d, m0, mHat, distC, grad, gradC_dummy, dHat, eAll)
    call deall_distortionParam(gradC_dummy)
    write(*,'(a39)') 'Gradient computed, line search finished'
    write(ioLog,'(a39)') 'Gradient computed, line search finished'

    call deall_dataVectorMTX(dHat_1)
    call deall_modelParam(mHat_0)
    call deall_modelParam(mHat_1)
    call deall_solnVectorMTX(eAll_1)

  end subroutine lineSearchDistCubic

end module NLCG_dist
