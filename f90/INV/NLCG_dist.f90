module NLCG_dist
! nonlinear conjugate gradient solver for joint inversion of conductivity 
! and distortion parameters
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
     integer            :: nCreset    ! reset alpha_C every nCreset iterations
     real(kind=prec)    :: alphaMinScale ! lower (relative) bound for alpha 
     real(kind=prec)    :: alphaMaxScale ! upper (relative) bound for alpha 
     logical            :: useBBDistInit ! Barzilai-Borwein initialization
     real(kind=prec)    :: lambdaClock ! lock C updates when lambda < this
     integer            :: nClockMinIter  ! only lock after this many iters
     character(80)      :: distFname
     character(80)      :: init_distFname
  end type NLCGDistControl_t

  type(NLCGDistControl_t), private, save :: distControl

Contains

   subroutine update_damping_parameter_dist(lambda, mHat, distNorm, F, gradM)
      real(kind=prec), intent(inout) :: lambda
      type(modelParam_t), intent(in) :: mHat
      real(kind=prec), intent(in) :: distNorm
      real(kind=prec), intent(inout) :: F
      type(modelParam_t), intent(inout) :: gradM

      real(kind=prec) :: lambdaPrev, SSplusDist, mNorm, Nmodel
      type(modelParam_t) :: dSS

      mNorm = dotProd(mHat, mHat)
      Nmodel = countModelParam(mHat)

      lambdaPrev = lambda
      SSplusDist = F - (lambdaPrev * mNorm / Nmodel)

      dSS = mHat
      call linComb(ONE, gradM, MinusTWO * lambdaPrev / Nmodel, mHat, dSS)

      lambda = lambda / distControl%k
      F = SSplusDist + (lambda * mNorm / Nmodel)

      call linComb(ONE, dSS, TWO * lambda / Nmodel, mHat, gradM)
      call deall_modelParam(dSS)
   end subroutine update_damping_parameter_dist

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
    ! below are distortion-specific parameters
    ctrl%nu = 10.0 ! note this is a large value, so the distortion 
                   ! penalty is large (typical for data without strong    
                    ! distortion).
    ctrl%nDistIter = 1 ! number of distortion line searches
    ctrl%alpha_C   = 0.1 ! initial step size for distortion line search
    ctrl%nCreset   = 5   ! periodic reset alpha_C to initial value
    ctrl%alphaMinScale = 0.01_prec
    ctrl%alphaMaxScale = 10.0_prec
    ctrl%useBBDistInit = .true. ! defaulted to true
    ctrl%lambdaClock = 10.0     ! simple first-stage C-lock threshold
    ctrl%nClockMinIter  = 10    ! avoid locking too early
    ctrl%distFname = 'distortion.dat'
    ctrl%init_distFname = ''
  end subroutine set_NLCGDistControl

   subroutine read_NLCGDistControl(ctrl, rFile, hasCtrlFile)
    type(NLCGDistControl_t), intent(inout) :: ctrl
    character(*), intent(in) :: rFile
    logical, intent(out), optional :: hasCtrlFile
    integer :: ios
    logical :: exists
    character(80) :: string

    call set_NLCGDistControl(ctrl)

    inquire(FILE=rFile, EXIST=exists)
    ! see if we have a control file. 
    if (present(hasCtrlFile)) hasCtrlFile = exists
    ! If not, just use the defaults.
    if (.not. exists) then 
        return
    endif

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
    ! Distortion related parameters
    read(ioInvCtrl, '(a36,g15.7)', iostat=ios) string, ctrl%nu
    if (output_level > 2) write(*,'(a36,g15.7)') string, ctrl%nu
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

    ! Optional advanced (reads: not very useful) controls
    ! read(ioInvCtrl, '(a36,i4)', iostat=ios) string, ctrl%nCreset
    ! if (output_level > 2) write(*,'(a36,i4)') string, ctrl%nCreset
    ! read(ioInvCtrl, '(a36,g15.7)', iostat=ios) string, ctrl%lambdaClock
    ! if (output_level > 2 .and. ios == 0) write(*,'(a36,g15.7)') string, ctrl%lambdaClock
    ! read(ioInvCtrl, '(a36,i4)', iostat=ios) string, ctrl%nClockMinIter
    ! if (output_level > 2 .and. ios == 0) write(*,'(a36,i4)') string, ctrl%nClockMinIter
    ! read(ioInvCtrl, '(a36,g15.7)', iostat=ios) string, ctrl%alphaMinScale
    ! if (output_level > 2 .and. ios == 0) write(*,'(a36,g15.7)') string, ctrl%alphaMinScale
    ! read(ioInvCtrl, '(a36,g15.7)', iostat=ios) string, ctrl%alphaMaxScale
    ! if (output_level > 2 .and. ios == 0) write(*,'(a36,g15.7)') string, ctrl%alphaMaxScale
    ! read(ioInvCtrl, '(a36,l4)', iostat=ios) string, ctrl%useBBDistInit
    ! if (output_level > 2 .and. ios == 0) write(*,'(a36,l4)') string, ctrl%useBBDistInit

    close(ioInvCtrl)
  end subroutine read_NLCGDistControl

  ! ---------------------------------------------------------------------------
  ! lineSearchDistQuad 
  ! simple quadratic fit using phi(0), phi(alpha_1), and the directional 
  ! derivative gradC(0). The accepted step is the quadratic minimizer, 
  ! evaluated once to return the actual objective value.
  !
  ! The EM fields eAll and conductivity sigma is unchanged — only C is updated.
  ! Therefore each trial objective evaluation calls phi_C_only, which requires
  ! only cheap per-site matrix products (no EM solve, no adjoint solve).
  ! 
  ! The search direction is the negative C-gradient: h_C = -gradC (steepest
  ! descent), and the step is initialised to alpha_C_init.  At each trial:
  !   C_trial = C + alpha * h_C
  !   phi_trial = phi_C_only(C_trial, ...)
  ! The directional derivative of phi w.r.t. alpha along h_C is:
  !   g_0 = <gradC, h_C> = -<gradC, gradC>  (always negative)
  ! ---------------------------------------------------------------------------
  subroutine lineSearchDistQuad(d, eAll, sigma, mHat, distC_in, gradC, &
                                lambda, nu, Ndata, mNorm_over_Nmodel,   &
                                alpha_C_init, phi_0,                    &
                                distC_out, phi_out, alpha_out, nLS)

    use DistortGradient, only: phi_C_only

    type(dataVectorMTX_t),   intent(in)    :: d
    type(solnVectorMTX_t),   intent(in)    :: eAll
    type(modelParam_t),      intent(in)    :: sigma
    type(modelParam_t),      intent(in)    :: mHat
    type(distortionParam_t), intent(in)    :: distC_in
    type(distortionParam_t), intent(in)    :: gradC
    real(kind=prec),         intent(in)    :: lambda, nu
    real(kind=prec),         intent(in)    :: Ndata
    real(kind=prec),         intent(in)    :: mNorm_over_Nmodel
    real(kind=prec),         intent(in)    :: alpha_C_init
    real(kind=prec),         intent(in)    :: phi_0
    type(distortionParam_t), intent(inout) :: distC_out
    real(kind=prec),         intent(out)   :: phi_out
    real(kind=prec),         intent(out)   :: alpha_out
    integer,                 intent(out)   :: nLS

    integer :: iSite
    real(kind=prec) :: alpha_1, alpha_q, g_0, phi_1, phi_q, a_quad
    real(kind=prec) :: alpha_lo, alpha_hi
    type(distortionParam_t) :: distC_trial
    character(100) :: logFile

    logFile = trim(distControl%fname)//'_NLCG_DIST.log'
    g_0 = -dotProd_distortionParam(gradC, gradC)

    if (abs(g_0) < R_TINY) then
       call copy_distortionParam(distC_in, distC_out)
       phi_out = phi_0
       alpha_out = alpha_C_init
       nLS = 0
       return
    end if

    alpha_1 = max(alpha_C_init, R_TINY)
    call create_distortionParam(distC_in%nSites, distC_trial)

    write(*,'(a16,a7,es12.5)') '   DistLS: init ', ' phi0=', phi_0
    write(ioLog,'(a16,a7,es12.5)') '   DistLS: init ', ' phi0=', phi_0

    do iSite = 1, distC_in%nSites
       distC_trial%C(1,1,iSite) = distC_in%C(1,1,iSite) - alpha_1 * gradC%C(1,1,iSite)
       distC_trial%C(1,2,iSite) = distC_in%C(1,2,iSite) - alpha_1 * gradC%C(1,2,iSite)
       distC_trial%C(2,1,iSite) = distC_in%C(2,1,iSite) - alpha_1 * gradC%C(2,1,iSite)
       distC_trial%C(2,2,iSite) = distC_in%C(2,2,iSite) - alpha_1 * gradC%C(2,2,iSite)
    end do
    distC_trial%siteIndex = distC_in%siteIndex

    phi_1 = phi_C_only(d, eAll, sigma, distC_trial, nu, Ndata, lambda, mNorm_over_Nmodel)
    nLS = 1
    write(*,'(a16,i3,a8,es12.5,a7,es12.5)') '   DistLS: trial', 1, ' alpha=', alpha_1, ' phi=', phi_1
    write(ioLog,'(a16,i3,a8,es12.5,a7,es12.5)') '   DistLS: trial', 1, ' alpha=', alpha_1, ' phi=', phi_1

    a_quad = (phi_1 - phi_0 - g_0 * alpha_1) / max(alpha_1**2, R_TINY)
    if (a_quad <= R_TINY) then
       call copy_distortionParam(distC_trial, distC_out)
       phi_out = phi_1
       alpha_out = alpha_1
       write(*,'(a54)') '   DistLS: nonpositive curvature, using trial step'
       write(ioLog,'(a54)') '   DistLS: nonpositive curvature, using trial step'
       call deall_distortionParam(distC_trial)
       return
    end if

    alpha_q = -g_0 / (TWO * a_quad)
    alpha_lo = max(distControl%alphaMinScale * alpha_1, R_TINY)
    alpha_hi = max(distControl%alphaMaxScale * alpha_1, alpha_lo)
    alpha_q = min(alpha_hi, max(alpha_lo, alpha_q))

    do iSite = 1, distC_in%nSites
       distC_out%C(1,1,iSite) = distC_in%C(1,1,iSite) - alpha_q * gradC%C(1,1,iSite)
       distC_out%C(1,2,iSite) = distC_in%C(1,2,iSite) - alpha_q * gradC%C(1,2,iSite)
       distC_out%C(2,1,iSite) = distC_in%C(2,1,iSite) - alpha_q * gradC%C(2,1,iSite)
       distC_out%C(2,2,iSite) = distC_in%C(2,2,iSite) - alpha_q * gradC%C(2,2,iSite)
    end do
    distC_out%siteIndex = distC_in%siteIndex

    phi_q = phi_C_only(d, eAll, sigma, distC_out, nu, Ndata, lambda, mNorm_over_Nmodel)
    nLS = 2
    write(*,'(a16,i3,a8,es12.5,a7,es12.5)') '   DistLS: trial', 2, ' alpha=', alpha_q, ' phi=', phi_q
    write(ioLog,'(a16,i3,a8,es12.5,a7,es12.5)') '   DistLS: trial', 2, ' alpha=', alpha_q, ' phi=', phi_q

    if (phi_1 < phi_q) then
       call copy_distortionParam(distC_trial, distC_out)
       phi_out = phi_1
       alpha_out = alpha_1
       write(*,'(a35)') '   DistLS: using initial trial step'
       write(ioLog,'(a35)') '   DistLS: using initial trial step'
    else
       phi_out = phi_q
       alpha_out = alpha_q
       write(*,'(a36)') '   DistLS: using quadratic minimizer'
       write(ioLog,'(a36)') '   DistLS: using quadratic minimizer'
    end if

    call deall_distortionParam(distC_trial)

  end subroutine lineSearchDistQuad

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
    real(kind=prec) :: alpha_C          ! warm-started C line search step
    real(kind=prec) :: alpha_C_try      ! initial alpha passed to C line search (warm-start or BB)
    real(kind=prec) :: alpha_BB
    real(kind=prec) :: ss_bb, sy_bb, s_ij, y_ij
    real(kind=prec) :: mNorm_over_Nmodel ! lambda*||mHat||^2/Nmodel for phi_C_only
    real(kind=prec) :: grad_dot_h, g_dot_g, g_dot_gPrev, gPrev_dot_gPrev
    real(kind=prec) :: g_dot_h
    integer :: iter, nCG, nLS, nLSC, nfunc, ios, i, j
    logical :: hasCtrlFile, lockC, hasHistory
    character(3) :: iterChar
    character(100) :: mFile, mHatFile, gradFile, dataFile, resFile, logFile
    character(100) :: distFile
    real(kind=prec) :: lambda0
    real(kind=prec) :: lambdaRatio

    distControl%nu = nu
    if (present(fname)) then
       call read_NLCGDistControl(distControl, fname, hasCtrlFile)
       if (hasCtrlFile) lambda = distControl%lambda
    else
       call set_NLCGDistControl(distControl)
    end if

    lambda0 = max(lambda, R_TINY)
    lockC = .false.
    hasHistory = .false.

    logFile = trim(distControl%fname)//'_NLCG_DIST.log'
    open(unit=ioLog, file=logFile, status='unknown', position='append', iostat=ios)

    alpha   = distControl%alpha_1
    alpha_C = distControl%alpha_C   ! warm-started; may grow/shrink each outer iteration

    write(*,'(a41,es8.1)') 'The initial damping parameter lambda is ',lambda
    write(ioLog,'(a41,es8.1)') 'The initial damping parameter lambda is ',lambda
    write(*,'(a47,es8.1)') ' The distortion regularization parameter nu is ',distControl%nu
    write(ioLog,'(a47,es8.1)') ' The distortion regularization parameter nu is ',distControl%nu
    if (distControl%lambdaClock > R_ZERO) then
       write(*,'(a40,es12.5)') 'C updates will be locked when lambda < ', distControl%lambdaClock
       write(ioLog,'(a40,es12.5)') 'C updates will be locked when lambda < ', distControl%lambdaClock
    end if

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
    write(*,'(a44,es12.5)') '     GRAD: initial norm of model gradient is', gnorm
    write(ioLog,'(a44,es12.5)') '     GRAD: initial norm of model gradient is', gnorm
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

       ! sigma-update (one NLCG sweep)
       write(*,'(a30)') 'Part 1: Model Sigma update '
       write(ioLog,'(a30)') 'Part 1: Model Sigma update '

       ! Line search for conductivity
       call lineSearchCubic(lambda, distControl%nu, d, m0, mHat, distC, &
            h, alpha, value, gradM, rms, nLS, dHat, eAll)
       nfunc = nfunc + nLS

       ! Update alpha for next sigma line search
       alpha = 2.0_prec * (value - valuePrev) / grad_dot_h
       alpha = (ONE + 0.01_prec) * alpha

       !write(*,'(a25,i5)') 'Completed sigma update, iter ', iter
       !write(ioLog,'(a25,i5)') 'Completed sigma update, iter ', iter

      ! C-update: hybrid quadratic + Armijo line search (no EM/adjoint solves)
       write(*,'(a35)') '   Part 2: Distortion Matrix update'
       write(ioLog,'(a35)') '   Part 2: Distortion Matrix update'

       lambdaRatio = lambda / lambda0
       if (.not. lockC) then
          if (iter >= distControl%nClockMinIter) then
             if (distControl%lambdaClock > R_ZERO .and. lambda < distControl%lambdaClock) then 
                lockC = .true.
                write(*,'(a38,i4,a10,es12.5,a16,es12.5)') '   &
                    DistLS: locking C-updates at iter ', iter, ', &
                    lambda=', lambda, ', lambda/lambda0=', lambdaRatio
                write(ioLog,'(a38,i4,a10,es12.5,a16,es12.5)') ' &
                    DistLS: locking C-updates at iter  ', iter, ', &
                    lambda=', lambda, ', lambda/lambda0=', lambdaRatio
             end if
          end if
       end if
       if (.not. lockC) then
          ! Periodically reset the warm-start step to avoid persistent shrinkage.
          if (distControl%nCreset > 0) then
             if (mod(iter, distControl%nCreset) == 0) then
                alpha_C = distControl%alpha_C
                hasHistory = .false.
                write(*,'(a39,es12.5)') '  DistLS: reset alpha_C warm-start to ', alpha_C
                write(ioLog,'(a39,es12.5)') '  DistLS: reset alpha_C warm-start to ', alpha_C
             end if
          end if
       end if


       ! Re-compute the smoothed model m from updated mHat
       call CmSqrtMult(mHat, m_minus_m0)
       call linComb(ONE, m_minus_m0, ONE, m0, m)

       Ndata = countData(d)
       Nmodel = countModelParam(mHat)
       mNorm_over_Nmodel = dotProd(mHat, mHat) / Nmodel

       if (.not. lockC) then
          do i = 1, distControl%nDistIter
             ! save current C before update; used for history.
             call copy_distortionParam(distC, gPrevC)

             ! gradient at current distC (eAll locked)
             call compDistGrad(d, eAll, m, distC, distControl%nu, Ndata, gradC)

             alpha_C_try = alpha_C
             if (distControl%useBBDistInit .and. hasHistory) then
                ! Barzilai-Borwein step size estimation for C-update
                ! after the first iteration, use history to compute a BB 
                ! step size. Cheap so it almost feels guilty not to use it  
                ! Although not a guarantee of convergence.
                ss_bb = R_ZERO
                sy_bb = R_ZERO
                do j = 1, distC%nSites
                   ! Compute s_ij = distC - hC and y_ij = gradC - gC for
                   ! each element of the 2x2 distortion matrix
                   s_ij = distC%C(1,1,j) - hC%C(1,1,j)
                   y_ij = gradC%C(1,1,j) - gC%C(1,1,j)
                   ss_bb = ss_bb + s_ij*s_ij
                   sy_bb = sy_bb + s_ij*y_ij

                   s_ij = distC%C(1,2,j) - hC%C(1,2,j)
                   y_ij = gradC%C(1,2,j) - gC%C(1,2,j)
                   ss_bb = ss_bb + s_ij*s_ij
                   sy_bb = sy_bb + s_ij*y_ij

                   s_ij = distC%C(2,1,j) - hC%C(2,1,j)
                   y_ij = gradC%C(2,1,j) - gC%C(2,1,j)
                   ss_bb = ss_bb + s_ij*s_ij
                   sy_bb = sy_bb + s_ij*y_ij

                   s_ij = distC%C(2,2,j) - hC%C(2,2,j)
                   y_ij = gradC%C(2,2,j) - gC%C(2,2,j)
                   ss_bb = ss_bb + s_ij*s_ij
                   sy_bb = sy_bb + s_ij*y_ij
                end do

                if (sy_bb > R_TINY .and. ss_bb > R_TINY) then
                   alpha_BB = ss_bb / sy_bb
                   ! some safeguards to avoid too small or too large step sizes
                   alpha_BB = max(distControl%alphaMinScale * distControl%alpha_C, alpha_BB)
                   alpha_BB = min(distControl%alphaMaxScale * distControl%alpha_C, alpha_BB)
                   alpha_C_try = alpha_BB
                   write(*,'(a24,es12.5,a15,es12.5,a16,es12.5,a1)') &
                       ' DistLS: BB init alpha=', alpha_C_try, ' (ss=', &
                        ss_bb, ', sy=', sy_bb, ')'
                   write(ioLog,'(a24,es12.5,a15,es12.5,a16,es12.5,a1)') &
                       '   DistLS: BB init alpha=', alpha_C_try, ' (ss=', &
                        ss_bb, ', sy=', sy_bb, ')'
                end if
             end if

             ! quadratic line search along -gradC direction
             ! hard coded to use the quadratic line search (instead of Armijo backtracking)
             call lineSearchDistQuad(d, eAll, m, mHat, distC, gradC,        &
                                   lambda, distControl%nu, Ndata,           &
                                   mNorm_over_Nmodel,                       &
                                   alpha_C_try, value,                      &
                                   distC_work, value, alpha_C, nLSC)

             ! Update history: store pre-step C and current gradient.
             call copy_distortionParam(gPrevC, hC)
             call copy_distortionParam(gradC, gC)
             hasHistory = .true.

             call copy_distortionParam(distC_work, distC)
             nfunc = nfunc + nLSC

             ! warm-start: slightly enlarge alpha_C for next outer iteration
             ! (mirrors the NLCG alpha update heuristic — limited to initial value)
             alpha_C = min(distControl%alpha_C, alpha_C * 2.0_prec)
          end do
       else
          write(*,'(a40)') '   DistLS: C-update skipped (locked)    '
          write(ioLog,'(a40)') '   DistLS: C-update skipped (locked)    '
       end if

       ! After all C-update line search.
       ! Recompute rms and dHat using locked EM fields (no new forward solve).
       ! mHat has not changed, so EM solution remains valid.
       call reComputedHat_dist(d, eAll, m, distC, dHat, rms)
       
       ! Recompute distortion regularization penalty (very cheap)
       distNorm = psi_C(distC)

       ! Important: Recompute sigma gradient at the updated (mHat, distC) state.
       ! This is used for the next outer iteration.
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

       ! Check for stall: update lambda (same RMS-only criterion as NLCG)
       if (abs(rmsPrev - rms) < distControl%fdiffTol) then
          call update_damping_parameter_dist(lambda, mHat, distNorm, value, gradM)

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
  subroutine lineSearchCubic(lambda, nu, d, m0, mHat, distC, h, alpha, f, grad, &
       rms, niter, dHat, eAll, gamma)

    real(kind=prec), intent(in) :: lambda, nu
    type(dataVectorMTX_t), intent(in) :: d
    type(modelParam_t), intent(in) :: m0
    type(modelParam_t), intent(inout) :: mHat
    type(distortionParam_t), intent(in) :: distC
    type(modelParam_t), intent(in) :: h
    real(kind=prec), intent(inout) :: alpha
    real(kind=prec), intent(inout) :: f
    type(modelParam_t), intent(in) :: grad
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
       call errStop("Linear step too large - try a smaller starting alphaC. Exiting..")
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
    write(*,'(a46)') 'Line search finished; gradient deferred to outer loop'
    write(ioLog,'(a46)') 'Line search finished; gradient deferred to outer loop'

    call deall_dataVectorMTX(dHat_1)
    call deall_modelParam(mHat_0)
    call deall_modelParam(mHat_1)
    call deall_solnVectorMTX(eAll_1)

  end subroutine lineSearchCubic

end module NLCG_dist
