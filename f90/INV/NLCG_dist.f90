module NLCG_dist
! nonlinear conjugate gradient solver for joint inversion of conductivity 
! and distortion parameters
  use ModEM_timers
  use invcore
  use DistortGradient, only: phi_C_only
  use transmitters, only: txDict
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
      integer            :: nSigmaIter ! sigma-only NLCG iterations per block
      real(kind=prec)    :: alpha_C
      real(kind=prec)    :: alphaMinScale ! lower (relative) bound for alpha 
      real(kind=prec)    :: alphaMaxScale ! upper (relative) bound for alpha 
      logical            :: useBBDistInit ! Barzilai-Borwein initialization
      real(kind=prec)    :: lambdaClock ! lock C updates when lambda < this
      integer            :: nClockMinIter  ! only lock after this many iters
      integer            :: nCUpdatesMax ! max number of C line-search updates (0 = no limit)
      integer            :: lambdaCooldown ! min iterations between lambda updates (0 = disabled)
      integer            :: nCPeriods ! number of shortest periods used to update C (0 = all)
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
    ctrl%nu = 3.0  ! note 10 is a large value, so the distortion 
                   ! penalty is large (typical for data without strong    
                   ! distortion). use 1-3 for data with visible distortion
    ctrl%nDistIter = 1 ! number of distortion line searches per block
    ctrl%nSigmaIter = 5 ! sigma-only NLCG iterations per block (C fixed)
    ctrl%alpha_C   = 0.1 ! initial step size for distortion line search
    ctrl%alphaMinScale = 0.01_prec
    ctrl%alphaMaxScale = 10.0_prec
    ctrl%useBBDistInit = .true. ! defaulted to true
    ctrl%lambdaClock = 10.0     ! simple first-stage C-lock threshold
    ctrl%nClockMinIter  = 10    ! avoid locking too early
    ctrl%nCUpdatesMax = 0       ! 0 = no limit on the number of C updates
    ctrl%lambdaCooldown = 0     ! 0 = no cooldown between lambda updates
    ctrl%nCPeriods = 0          ! 0 = use all periods for the C update
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
    ! Distortion related parameters.  The minimal user-facing set is
    ! nu (distortion damping) followed by the optional output/restart
    ! file names.  The remaining lines are advanced controls with
    ! defaults; they may be omitted entirely.
    read(ioInvCtrl, '(a36,g15.7)', iostat=ios) string, ctrl%nu
    if (ios == 0) then
       if (output_level > 2) write(*,'(a36,g15.7)') string, ctrl%nu
    else
       ctrl%nu = 3.0_prec
    end if

    read(ioInvCtrl, '(a36,a80)', iostat=ios) string, ctrl%distFname
    if (ios == 0) then
       ctrl%distFname = adjustl(ctrl%distFname)
       if (output_level > 2) write(*,'(a36,a80)') string, ctrl%distFname
    else
       ctrl%distFname = 'distortion.dat'
    end if

    read(ioInvCtrl, '(a36,a80)', iostat=ios) string, ctrl%init_distFname
    if (ios == 0) then
       ctrl%init_distFname = adjustl(ctrl%init_distFname)
       if (len_trim(ctrl%init_distFname) == 0) ctrl%init_distFname = ''
       if (output_level > 2) write(*,'(a36,a80)') string, ctrl%init_distFname
    else
       ctrl%init_distFname = ''
    end if

    ! Advanced controls (optional).  Defaults are set in
    ! set_NLCGDistControl; these lines are only needed for tuning.
    read(ioInvCtrl, '(a36,i4)', iostat=ios) string, ctrl%nSigmaIter
    if (ios == 0) then
       if (output_level > 2) write(*,'(a36,i4)') string, ctrl%nSigmaIter
    else
       ctrl%nSigmaIter = 5
    end if

    read(ioInvCtrl, '(a36,i4)', iostat=ios) string, ctrl%nDistIter
    if (ios == 0) then
       if (output_level > 2) write(*,'(a36,i4)') string, ctrl%nDistIter
    else
       ctrl%nDistIter = 1
    end if

    read(ioInvCtrl, '(a36,g15.7)', iostat=ios) string, ctrl%alpha_C
    if (ios == 0) then
       if (output_level > 2) write(*,'(a36,g15.7)') string, ctrl%alpha_C
    else
       ctrl%alpha_C = 0.1_prec
    end if

    ! C-update lock: stop updating the distortion matrix once lambda has
    ! dropped below lambdaClock (only considered after nClockMinIter model
    ! iterations).  A lambdaClock of 0 disables locking.
    read(ioInvCtrl, '(a36,g15.7)', iostat=ios) string, ctrl%lambdaClock
    if (ios == 0) then
       if (output_level > 2) write(*,'(a36,g15.7)') string, ctrl%lambdaClock
    else
       ctrl%lambdaClock = 10.0_prec
    end if

    read(ioInvCtrl, '(a36,i4)', iostat=ios) string, ctrl%nClockMinIter
    if (ios == 0) then
       if (output_level > 2) write(*,'(a36,i4)') string, ctrl%nClockMinIter
    else
       ctrl%nClockMinIter = 10
    end if

    ! Optional hard cap on the number of C updates.  Once the cap is
    ! reached the distortion matrix is locked and the remaining iterations
    ! are sigma-only.  0 disables the cap.
    read(ioInvCtrl, '(a36,i4)', iostat=ios) string, ctrl%nCUpdatesMax
    if (ios == 0) then
       if (output_level > 2) write(*,'(a36,i4)') string, ctrl%nCUpdatesMax
    else
       ctrl%nCUpdatesMax = 0
    end if

    ! Optional cooldown between lambda reductions.  A value of 0 disables
    ! the cooldown; otherwise lambda is updated at most once per
    ! lambdaCooldown model iterations.
    read(ioInvCtrl, '(a36,i4)', iostat=ios) string, ctrl%lambdaCooldown
    if (ios == 0) then
       if (output_level > 2) write(*,'(a36,i4)') string, ctrl%lambdaCooldown
    else
       ctrl%lambdaCooldown = 0
    end if

    ! Optional restriction of the distortion update to the shortest periods.
    ! C is frequency independent, but the near-surface galvanic distortion is
    ! best constrained by the short-period responses.  A value of 0 uses all
    ! periods (current behaviour); N > 0 uses only the N shortest periods.
    read(ioInvCtrl, '(a36,i4)', iostat=ios) string, ctrl%nCPeriods
    if (ios == 0) then
       if (output_level > 2) write(*,'(a36,i4)') string, ctrl%nCPeriods
    else
       ctrl%nCPeriods = 0
    end if

    ! read(ioInvCtrl, '(a36,g15.7)', iostat=ios) string, ctrl%alphaMinScale
    ! if (output_level > 2 .and. ios == 0) write(*,'(a36,g15.7)') string, ctrl%alphaMinScale
    ! read(ioInvCtrl, '(a36,g15.7)', iostat=ios) string, ctrl%alphaMaxScale
    ! if (output_level > 2 .and. ios == 0) write(*,'(a36,g15.7)') string, ctrl%alphaMaxScale
    ! read(ioInvCtrl, '(a36,l4)', iostat=ios) string, ctrl%useBBDistInit
    ! if (output_level > 2 .and. ios == 0) write(*,'(a36,l4)') string, ctrl%useBBDistInit

    close(ioInvCtrl)
  end subroutine read_NLCGDistControl

  ! ---------------------------------------------------------------------------
  ! find_C_period_cutoff returns the nPer-th shortest transmitter period.  The
  ! C update then uses only data blocks with period <= periodMax.  If nPer <= 0
  ! or nPer covers all transmitters, periodMax is effectively infinite.
  ! ---------------------------------------------------------------------------
  subroutine find_C_period_cutoff(d, nPer, periodMax)
    type(dataVectorMTX_t), intent(in) :: d
    integer, intent(in) :: nPer
    real(kind=prec), intent(out) :: periodMax

    real(kind=prec), allocatable :: periods(:)
    real(kind=prec) :: tmp
    integer :: j, k, iMin, n

    periodMax = huge(1.0_prec)
    if (nPer <= 0) return

    n = d%nTx
    if (n <= 0) return
    allocate(periods(n))
    do j = 1, n
       periods(j) = txDict(d%d(j)%tx)%period
    end do

    ! Selection sort (nTx is small: one entry per period for MT).
    do j = 1, n - 1
       iMin = j
       do k = j + 1, n
          if (periods(k) < periods(iMin)) iMin = k
       end do
       if (iMin /= j) then
          tmp = periods(j)
          periods(j) = periods(iMin)
          periods(iMin) = tmp
       end if
    end do

    if (nPer < n) periodMax = periods(nPer)
    deallocate(periods)
  end subroutine find_C_period_cutoff

  ! ---------------------------------------------------------------------------
  ! countData_period counts the existing data values in all blocks with
  ! period <= periodMax.  This is the normalisation used by the C objective
  ! and C gradient when the distortion update is restricted to short periods.
  ! ---------------------------------------------------------------------------
  function countData_period(d, periodMax) result(ndata)
    type(dataVectorMTX_t), intent(in) :: d
    real(kind=prec), intent(in) :: periodMax
    integer :: ndata, j, i

    ndata = 0
    do j = 1, d%nTx
       if (txDict(d%d(j)%tx)%period > periodMax) cycle
       do i = 1, d%d(j)%nDt
          ndata = ndata + count(d%d(j)%data(i)%exist)
       end do
    end do
  end function countData_period

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
                                distC_out, phi_out, alpha_out, nLS,     &
                                periodMax)

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
    real(kind=prec),         intent(in), optional :: periodMax

    integer :: iSite, nBT
    real(kind=prec) :: alpha_1, alpha_q, g_0, phi_1, phi_q, a_quad
    real(kind=prec) :: alpha_lo, alpha_hi, alpha_try, phi_try
    real(kind=prec) :: c_armijo, alpha_min
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

    c_armijo = distControl%c
    alpha_1 = max(alpha_C_init, R_TINY)
    alpha_min = max(distControl%alphaMinScale * alpha_1, R_TINY)
    call create_distortionParam(distC_in%nSites, distC_trial)

    if (output_level > 3) then
       write(*,'(a16,a7,es12.5)') '   DistLS: init ', ' phi0=', phi_0
       write(ioLog,'(a16,a7,es12.5)') '   DistLS: init ', ' phi0=', phi_0
    end if

    nLS = 0

    ! Trial 1: initial step.
    call linComb_distortionParam(ONE, distC_in, -alpha_1, gradC, distC_trial)
    phi_1 = phi_C_only(d, eAll, sigma, distC_trial, nu, Ndata, lambda, mNorm_over_Nmodel, periodMax)
    nLS = nLS + 1
    if (output_level > 3) then
       write(*,'(a16,i3,a8,es12.5,a7,es12.5)') '   DistLS: trial', nLS, ' alpha=', alpha_1, ' phi=', phi_1
       write(ioLog,'(a16,i3,a8,es12.5,a7,es12.5)') '   DistLS: trial', nLS, ' alpha=', alpha_1, ' phi=', phi_1
    end if

    ! Quadratic fit and minimizer.  Keep the better of the initial and
    ! quadratic trials as the candidate, then enforce Armijo on that candidate.
    a_quad = (phi_1 - phi_0 - g_0 * alpha_1) / max(alpha_1**2, R_TINY)
    alpha_try = alpha_1
    phi_try = phi_1

    if (a_quad > R_TINY) then
       alpha_q = -g_0 / (TWO * a_quad)
       alpha_lo = max(distControl%alphaMinScale * alpha_1, R_TINY)
       alpha_hi = max(distControl%alphaMaxScale * alpha_1, alpha_lo)
       alpha_q = min(alpha_hi, max(alpha_lo, alpha_q))

       call linComb_distortionParam(ONE, distC_in, -alpha_q, gradC, distC_trial)
       phi_q = phi_C_only(d, eAll, sigma, distC_trial, nu, Ndata, lambda, mNorm_over_Nmodel, periodMax)
       nLS = nLS + 1
       if (output_level > 3) then
          write(*,'(a16,i3,a8,es12.5,a7,es12.5)') '   DistLS: trial', nLS, ' alpha=', alpha_q, ' phi=', phi_q
          write(ioLog,'(a16,i3,a8,es12.5,a7,es12.5)') '   DistLS: trial', nLS, ' alpha=', alpha_q, ' phi=', phi_q
       end if

       if (phi_q < phi_try) then
          alpha_try = alpha_q
          phi_try = phi_q
       end if
    end if

    ! Accept the best candidate only if it satisfies Armijo.
    if (phi_try <= phi_0 + c_armijo * alpha_try * g_0) then
       call linComb_distortionParam(ONE, distC_in, -alpha_try, gradC, distC_out)
       phi_out = phi_try
       alpha_out = alpha_try
       if (output_level > 3) then
          write(*,'(a44)') '   DistLS: best trial satisfies Armijo'
          write(ioLog,'(a44)') '   DistLS: best trial satisfies Armijo'
       end if
       goto 100
    end if

    ! Backtracking from the best candidate.
    do nBT = 1, 10
       alpha_try = alpha_try * 0.5_prec
       if (alpha_try < alpha_min) exit

       call linComb_distortionParam(ONE, distC_in, -alpha_try, gradC, distC_trial)
       phi_try = phi_C_only(d, eAll, sigma, distC_trial, nu, Ndata, lambda, mNorm_over_Nmodel, periodMax)
       nLS = nLS + 1
       if (output_level > 3) then
          write(*,'(a16,i3,a8,es12.5,a7,es12.5)') '   DistLS: trial', nLS, ' alpha=', alpha_try, ' phi=', phi_try
          write(ioLog,'(a16,i3,a8,es12.5,a7,es12.5)') '   DistLS: trial', nLS, ' alpha=', alpha_try, ' phi=', phi_try
       end if

       if (phi_try <= phi_0 + c_armijo * alpha_try * g_0) then
          call linComb_distortionParam(ONE, distC_in, -alpha_try, gradC, distC_out)
          phi_out = phi_try
          alpha_out = alpha_try
          if (output_level > 3) then
             write(*,'(a44)') '   DistLS: backtracked step satisfies Armijo'
             write(ioLog,'(a44)') '   DistLS: backtracked step satisfies Armijo'
          end if
          goto 100
       end if
    end do

    ! No acceptable step: keep the current C.
    call copy_distortionParam(distC_in, distC_out)
    phi_out = phi_0
    alpha_out = alpha_C_init
    write(*,'(a52)') '   DistLS: no Armijo step found, keeping current C'
    write(ioLog,'(a52)') '   DistLS: no Armijo step found, keeping current C'

100 continue
    if (output_level <= 3) then
       write(*,'(a26,es12.5,a10,es12.5)') '   DistLS: accepted alpha=', alpha_out, ' phi=', phi_out
       write(ioLog,'(a26,es12.5,a10,es12.5)') '   DistLS: accepted alpha=', alpha_out, ' phi=', phi_out
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
    type(distortionParam_t) :: distC_prev, distC_diff
    type(solnVectorMTX_t) :: eAll
    real(kind=prec) :: value, valuePrev, rms, rmsPrev
    real(kind=prec) :: alpha, beta, gnorm, mNorm, distNorm, Ndata_C, Nmodel
    real(kind=prec) :: alpha_C          ! warm-started C line search step
    real(kind=prec) :: alpha_C_try      ! initial alpha passed to C line search (warm-start or BB)
    real(kind=prec) :: alpha_BB
    real(kind=prec) :: ss_bb, sy_bb, s_ij, y_ij
    real(kind=prec) :: mNorm_over_Nmodel ! lambda*||mHat||^2/Nmodel for phi_C_only
    real(kind=prec) :: grad_dot_h, g_dot_g, g_dot_gPrev, gPrev_dot_gPrev
    real(kind=prec) :: g_dot_h
    real(kind=prec) :: value_before_C, value_after_C, rms_before_C, rms_after_C
    real(kind=prec) :: value_before_sigma, value_after_sigma, rms_before_sigma, rms_after_sigma
    real(kind=prec) :: dF_C, dF_sigma, share_C, dPhi_C
    real(kind=prec) :: phi_C0, phi_C_out, phi_C_before, phi_C_after
    real(kind=prec) :: periodMax
    real(kind=prec) :: gradC_norm, gradM_norm
    real(kind=prec) :: dC_norm, dC_max
    integer :: iter, nCG, nLS, nLSC, nfuncSigma, nfuncC, ios, i, j, inner, iSite
    integer :: nCUpdates, lastLambdaUpdateIter
    logical :: hasCtrlFile, lockC, hasHistory, stopInversion
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

    ! Optionally restrict the C update to the shortest periods.  periodMax is
    ! effectively infinite when nCPeriods = 0, preserving the previous
    ! all-period behaviour.
    call find_C_period_cutoff(d, distControl%nCPeriods, periodMax)
    Ndata_C = countData_period(d, periodMax)
    if (distControl%nCPeriods > 0 .and. Ndata_C <= 0) then
       call errStop('No data found in the selected C periods')
    end if

    lambda0 = max(lambda, R_TINY)
    lockC = .false.
    hasHistory = .false.
    stopInversion = .false.
    nCUpdates = 0
    lastLambdaUpdateIter = -distControl%lambdaCooldown

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
    write(*,'(a46,i4)') ' Distortion sigma iterations per block: ', distControl%nSigmaIter
    write(ioLog,'(a46,i4)') ' Distortion sigma iterations per block: ', distControl%nSigmaIter
    write(*,'(a46,i4)') ' Distortion steps per block: ', distControl%nDistIter
    write(ioLog,'(a46,i4)') ' Distortion steps per block: ', distControl%nDistIter
    write(*,'(a46,es12.5)') ' Initial distortion step size: ', distControl%alpha_C
    write(ioLog,'(a46,es12.5)') ' Initial distortion step size: ', distControl%alpha_C
    write(*,'(a46,i4)') ' Minimum iterations before C lock: ', distControl%nClockMinIter
    write(ioLog,'(a46,i4)') ' Minimum iterations before C lock: ', distControl%nClockMinIter
    write(*,'(a46,i4)') ' Maximum number of C updates: ', distControl%nCUpdatesMax
    write(ioLog,'(a46,i4)') ' Maximum number of C updates: ', distControl%nCUpdatesMax
    write(*,'(a46,i4)') ' Lambda update cooldown: ', distControl%lambdaCooldown
    write(ioLog,'(a46,i4)') ' Lambda update cooldown: ', distControl%lambdaCooldown
    if (distControl%nCPeriods > 0) then
       write(*,'(a46,i4)') ' Distortion periods used for C: ', distControl%nCPeriods
       write(ioLog,'(a46,i4)') ' Distortion periods used for C: ', distControl%nCPeriods
       write(*,'(a46,es12.5)') ' Distortion period cutoff (s): ', periodMax
       write(ioLog,'(a46,es12.5)') ' Distortion period cutoff (s): ', periodMax
    end if
    write(*,'(a46,a)') ' Distortion output file: ', trim(distControl%distFname)
    write(ioLog,'(a46,a)') ' Distortion output file: ', trim(distControl%distFname)

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
    call create_distortionParam(distC%nSites, distC_prev)
    call create_distortionParam(distC%nSites, distC_diff)

    ! First forward pass and function evaluation
    call func_dist(lambda, distControl%nu, d, m0, mHat, distC, value, mNorm, distNorm, dHat, eAll, rms)
    call printf('START', lambda, alpha, value, mNorm, rms)
    call printf('START', lambda, alpha, value, mNorm, rms, logFile)
    nfuncSigma = 1
    nfuncC = 0

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

    ! Initialise the CG loop counters; the first gradient (and its alpha) is
    ! computed at the top of Phase 2 below, using the post-C-update gradient
    ! that the first sigma line search actually employs.
    nCG = 0
    iter = 0

    call ModEM_timers_create("NLCG_DIST Iteration", .false.)

    ! Main block loop: C-update first, then nSigmaIter sigma-only iterations
    blockLoop: do
       if ((rms < distControl%rmsTol) .or. (iter >= distControl%maxIter)) exit blockLoop

       value_before_C = value
       rms_before_C = rms

       ! === Phase 1: C-update at current (mHat, eAll) state ===
       write(*,'(a40)') ' Block: Distortion Matrix update (Phase 1)'
       write(ioLog,'(a40)') ' Block: Distortion Matrix update (Phase 1)'

       ! Optional hard cap on the number of C updates.  Once the cap is
       ! reached, lock C and leave the remaining iterations to sigma.
       if (distControl%nCUpdatesMax > 0 .and. nCUpdates >= distControl%nCUpdatesMax) then
          if (.not. lockC) then
             lockC = .true.
             write(*,'(a46,i4,a15,i4,a1)') '   DistLS: locking C-updates after ', &
                  nCUpdates, ' updates (max ', distControl%nCUpdatesMax, ')'
             write(ioLog,'(a46,i4,a15,i4,a1)') '   DistLS: locking C-updates after ', &
                  nCUpdates, ' updates (max ', distControl%nCUpdatesMax, ')'
          end if
       end if

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
          ! Re-compute smoothed model m from current mHat
          call CmSqrtMult(mHat, m_minus_m0)
          call linComb(ONE, m_minus_m0, ONE, m0, m)

          mNorm_over_Nmodel = dotProd(mHat, mHat) / countModelParam(mHat)

          call copy_distortionParam(distC, distC_prev)

          ! Reference C objective before the update, used for the optional
          ! short-period diagnostic when nCPeriods > 0.
          phi_C_before = R_ZERO
          if (distControl%nCPeriods > 0) then
             phi_C_before = phi_C_only(d, eAll, m, distC, distControl%nu, &
                                       Ndata_C, lambda, mNorm_over_Nmodel, periodMax)
          end if

          do i = 1, distControl%nDistIter
             call copy_distortionParam(distC, gPrevC)
             call compDistGrad(d, eAll, m, distC, distControl%nu, Ndata_C, gradC, periodMax)
             gradC_norm = sqrt(dotProd_distortionParam(gradC, gradC))
             if (output_level > 3) then
                write(*,'(a30,es12.5)') '   DistDiag: ||gradC|| = ', gradC_norm
                write(ioLog,'(a30,es12.5)') '   DistDiag: ||gradC|| = ', gradC_norm
             end if

             alpha_C_try = alpha_C
             if (distControl%useBBDistInit .and. hasHistory) then
                ss_bb = R_ZERO
                sy_bb = R_ZERO
                do j = 1, distC%nSites
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
                   alpha_BB = max(distControl%alphaMinScale * distControl%alpha_C, alpha_BB)
                   alpha_BB = min(distControl%alphaMaxScale * distControl%alpha_C, alpha_BB)
                   alpha_C_try = alpha_BB
                   if (output_level > 3) then
                      write(*,'(a24,es12.5,a15,es12.5,a16,es12.5,a1)') &
                          ' DistLS: BB init alpha=', alpha_C_try, ' (ss=', ss_bb, ', sy=', sy_bb, ')'
                      write(ioLog,'(a24,es12.5,a15,es12.5,a16,es12.5,a1)') &
                          '   DistLS: BB init alpha=', alpha_C_try, ' (ss=', ss_bb, ', sy=', sy_bb, ')'
                   end if
                end if
             end if

             ! The C line search minimises the C-only objective (short-period
             ! data when nCPeriods > 0), not the full F.  phi_C_out is therefore
             ! kept separate and the full F is recomputed after the update.
             phi_C0 = phi_C_only(d, eAll, m, distC, distControl%nu, Ndata_C, &
                                 lambda, mNorm_over_Nmodel, periodMax)
             call lineSearchDistQuad(d, eAll, m, mHat, distC, gradC,        &
                                   lambda, distControl%nu, Ndata_C,         &
                                   mNorm_over_Nmodel,                       &
                                   alpha_C_try, phi_C0,                     &
                                   distC_work, phi_C_out, alpha_C, nLSC,    &
                                   periodMax)

             call copy_distortionParam(gPrevC, hC)
             call copy_distortionParam(gradC, gC)
             hasHistory = .true.
             call copy_distortionParam(distC_work, distC)
             nfuncC = nfuncC + nLSC
             alpha_C = min(distControl%alpha_C, alpha_C * 2.0_prec)
          end do
          nCUpdates = nCUpdates + distControl%nDistIter

          ! Recompute dHat and rms with new C (cheap, locked eAll).  The full
          ! objective is rebuilt from the full RMS, model penalty, and C
          ! regularisation; this stays correct even when the C line search
          ! used only the short-period data.
          call reComputedHat_dist(d, eAll, m, distC, dHat, rms)
          value = rms**2 + lambda * mNorm_over_Nmodel &
                + distControl%nu * psi_C(distC)

          if (distControl%nCPeriods > 0) then
             phi_C_after = phi_C_only(d, eAll, m, distC, distControl%nu, &
                                      Ndata_C, lambda, mNorm_over_Nmodel, periodMax)
             dPhi_C = phi_C_before - phi_C_after
          else
             dPhi_C = R_ZERO
          end if

          value_after_C = value
          rms_after_C = rms
          dF_C = value_before_C - value_after_C
          call linComb_distortionParam(ONE, distC, MinusONE, distC_prev, distC_diff)
          dC_norm = sqrt(dotProd_distortionParam(distC_diff, distC_diff))
          dC_max = R_ZERO
          do iSite = 1, distC%nSites
             do j = 1, 2
                dC_max = max(dC_max, abs(distC_diff%C(1,j,iSite)), abs(distC_diff%C(2,j,iSite)))
             end do
          end do
          if (output_level > 3) then
             write(*,'(a20,es12.5,a20,es12.5,a20,es12.5)') &
                  ' DistDiag: dF_C=', dF_C, ' rms_before=', rms_before_C, ' rms_after=', rms_after_C
             write(ioLog,'(a20,es12.5,a20,es12.5,a20,es12.5)') &
                  ' DistDiag: dF_C=', dF_C, ' rms_before=', rms_before_C, ' rms_after=', rms_after_C
             write(*,'(a20,es12.5,a20,es12.5)') ' DistDiag: ||dC||=', dC_norm, ' max|dC|=', dC_max
             write(ioLog,'(a20,es12.5,a20,es12.5)') ' DistDiag: ||dC||=', dC_norm, ' max|dC|=', dC_max
          end if
       else
          if (output_level > 3) then
             write(*,'(a40)') '   DistLS: C-update skipped (locked)    '
             write(ioLog,'(a40)') '   DistLS: C-update skipped (locked)    '
          end if
          value_after_C = value_before_C
          rms_after_C = rms_before_C
          dF_C = R_ZERO
          dPhi_C = R_ZERO
          phi_C_after = R_ZERO
          dC_norm = R_ZERO
          dC_max = R_ZERO
       end if

       value_before_sigma = value
       rms_before_sigma = rms

       ! === Phase 2: nSigmaIter sigma-only NLCG iterations (C fixed) ===
       if (.not. lockC) then
          ! C was updated: recompute the sigma gradient at the new C and
          ! restart CG from steepest descent.
          call gradient_dist(lambda, distControl%nu, d, m0, mHat, distC, gradM, gradC, dHat, eAll, needGradC=.false.)
          call linComb(MinusONE, gradM, R_ZERO, gradM, g)
          h = g
          nCG = 0
          gnorm = sqrt(dotProd(gradM, gradM))
          gradM_norm = gnorm
          if (output_level > 3) then
             write(*,'(a30,es12.5)') '   DistDiag: ||gradM|| = ', gradM_norm
             write(ioLog,'(a30,es12.5)') '   DistDiag: ||gradM|| = ', gradM_norm
          end if
          ! On the first outer cycle report the post-C-update gradient (the one
          ! the first sigma line search actually uses) and initialise alpha.
          if (iter == 0) then
             write(*,'(a44,es12.5)') '     GRAD: initial norm of model gradient is', gnorm
             write(ioLog,'(a44,es12.5)') '     GRAD: initial norm of model gradient is', gnorm
             if (gnorm < TOL6) then
                call errStop('Problem with gradient computations: first gradient is zero')
             else
                alpha = distControl%startdm / gnorm
                write(*,'(a39,es12.5)') 'The initial value of alpha updated to ', alpha
                write(ioLog,'(a39,es12.5)') 'The initial value of alpha updated to ', alpha
             end if
          else
             ! Warm-start the sigma step from the previous estimate instead of
             ! resetting to a full startdm/gnorm step after every C update.  A
             ! C update changes the model gradient and the objective landscape;
             ! forcing a full-length step here makes the first sigma trial
             ! overshoot badly.  Keep the previous step if it is shorter, and
             ! never exceed the configured initial step length.
             alpha = min(alpha, distControl%startdm / max(gnorm, R_TINY))
             alpha = max(alpha, R_TINY)
             if (output_level > 3) then
                write(*,'(a40,es12.5,a20,es12.5)') &
                     ' DistDiag: sigma warm-start alpha=', alpha, ' step_norm=', alpha*gnorm
                write(ioLog,'(a40,es12.5,a20,es12.5)') &
                     ' DistDiag: sigma warm-start alpha=', alpha, ' step_norm=', alpha*gnorm
             end if
          end if
          if (output_level > 3) then
             write(*,'(a39)') 'CG restarted after C-update (steepest descent)'
             write(ioLog,'(a39)') 'CG restarted after C-update (steepest descent)'
          end if
       else
          ! C is locked: the model state, gradient, and CG direction are
          ! unchanged from the end of the previous sigma block.  Keep them
          ! instead of resetting to steepest descent.
          gnorm = sqrt(dotProd(gradM, gradM))
          gradM_norm = gnorm
          if (output_level > 3) then
             write(*,'(a30,es12.5)') '   DistDiag: ||gradM|| = ', gradM_norm
             write(ioLog,'(a30,es12.5)') '   DistDiag: ||gradM|| = ', gradM_norm
             write(*,'(a52)') 'C locked: continuing CG without reset'
             write(ioLog,'(a52)') 'C locked: continuing CG without reset'
          end if
       end if
       write(*,'(a36,i3,a16,i3)') 'Phase 2: Inner sigma loop (up to ', distControl%nSigmaIter, ' iterations)'
       write(ioLog,'(a36,i3,a16,i3)') 'Phase 2: Inner sigma loop (up to ', distControl%nSigmaIter, ' iterations)'

       sigmaLoop: do inner = 1, distControl%nSigmaIter
          if ((rms < distControl%rmsTol) .or. (iter >= distControl%maxIter)) exit sigmaLoop

          iter = iter + 1
          call ModEM_timers_start("NLCG_DIST Iteration")

          rmsPrev = rms
          valuePrev = value
          grad_dot_h = dotProd(gradM, h)

          ! Sigma line search with current (fixed) distC
          write(*,'(a30,i5)') ' Sigma NLCG iteration ',iter
          write(ioLog,'(a30,i5)') ' Sigma NLCG iteration ',iter
          call lineSearchCubic(lambda, distControl%nu, d, m0, mHat, distC, &
               h, alpha, value, gradM, rms, nLS, dHat, eAll)
          nfuncSigma = nfuncSigma + nLS

           ! Update alpha for next sigma line search
           alpha = 2.0_prec * (value - valuePrev) / grad_dot_h
           alpha = (ONE + 0.01_prec) * alpha

           ! Stall check (before CG beta, so lambda update resets CG correctly)
           distNorm = psi_C(distC)
           if (abs(rmsPrev - rms) < distControl%fdiffTol) then
              if (distControl%lambdaCooldown > 0 .and. &
                  (iter - lastLambdaUpdateIter) < distControl%lambdaCooldown) then
                 if (output_level > 3) then
                    write(*,'(a48,i4,a16,i4)') '   Lambda update deferred: cooldown ', &
                         distControl%lambdaCooldown, ' iter since last=', iter - lastLambdaUpdateIter
                    write(ioLog,'(a48,i4,a16,i4)') '   Lambda update deferred: cooldown ', &
                         distControl%lambdaCooldown, ' iter since last=', iter - lastLambdaUpdateIter
                 end if
              else
                 call update_damping_parameter_dist(lambda, mHat, distNorm, value, gradM)
                 lastLambdaUpdateIter = iter
                 if (lambda < distControl%lambdaTol) then
                    write(*,'(a55)') 'Unable to get out of a local minimum. Exiting...'
                    write(ioLog,'(a55)') 'Unable to get out of a local minimum. Exiting...'
                    stopInversion = .true.
                    exit sigmaLoop
                 end if
                 gnorm = sqrt(dotProd(gradM, gradM))
                 alpha = min(ONE, distControl%startdm) / gnorm
                 call linComb(MinusONE, gradM, R_ZERO, gradM, g)
                 write(*,'(a55)') 'Restarting NLCG with damping parameter updated'
                 write(ioLog,'(a55)') 'Restarting NLCG with damping parameter updated'
                 Nmodel = countModelParam(mHat)
                 mNorm = dotProd(mHat, mHat) / Nmodel
                 call printf('to', lambda, alpha, value, mNorm, rms)
                 call printf('to', lambda, alpha, value, mNorm, rms, logFile)
                 call ModEM_timers_stop("NLCG_DIST Iteration", .false.)
                 call ModEM_timers_print("NLCG_DIST Iteration", ioLog)
                 h = g
                 nCG = 0
                 cycle sigmaLoop
              end if
           end if

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

           ! Sigma gradient at new mHat (C still fixed)
           gPrev = g
           call gradient_dist(lambda, distControl%nu, d, m0, mHat, distC, gradM, gradC, dHat, eAll, needGradC=.false.)
           call linComb(MinusONE, gradM, R_ZERO, gradM, g)

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

           ! Timer/printf at end of inner loop
           Nmodel = countModelParam(mHat)
           mNorm = dotProd(mHat, mHat) / Nmodel
           call printf('with', lambda, alpha, value, mNorm, rms)
           call printf('with', lambda, alpha, value, mNorm, rms, logFile)

           call ModEM_timers_stop("NLCG_DIST Iteration", .false.)
           call ModEM_timers_print("NLCG_DIST Iteration", ioLog)

        end do sigmaLoop

        value_after_sigma = value
        rms_after_sigma = rms
        dF_sigma = value_before_sigma - value_after_sigma
        if (abs(dF_C) + abs(dF_sigma) > R_TINY) then
           share_C = dF_C / (abs(dF_C) + abs(dF_sigma))
        else
           share_C = R_ZERO
        end if
        if (output_level > 3) then
           write(*,'(a20,es12.5,a20,es12.5,a20,es12.5)') &
                ' DistDiag: dF_sigma=', dF_sigma, ' rms_before=', rms_before_sigma, ' rms_after=', rms_after_sigma
           write(ioLog,'(a20,es12.5,a20,es12.5,a20,es12.5)') &
                ' DistDiag: dF_sigma=', dF_sigma, ' rms_before=', rms_before_sigma, ' rms_after=', rms_after_sigma
           if (distControl%nCPeriods > 0) then
              write(*,'(a20,es12.5,a20,es12.5)') &
                   ' DistDiag: dPhi_C=', dPhi_C, ' phi_C_after=', phi_C_after
              write(ioLog,'(a20,es12.5,a20,es12.5)') &
                   ' DistDiag: dPhi_C=', dPhi_C, ' phi_C_after=', phi_C_after
           end if
        end if
        if (distControl%nCPeriods > 0) then
           write(*,'(a18,es11.4,a12,es11.4,a12,es11.4,a12,es11.4,a10,es11.4)') &
                ' DistDiag: dF_C=', dF_C, ' dF_sigma=', dF_sigma, ' dPhi_C=', dPhi_C, &
                ' ||dC||=', dC_norm, ' rms=', rms
           write(ioLog,'(a18,es11.4,a12,es11.4,a12,es11.4,a12,es11.4,a10,es11.4)') &
                ' DistDiag: dF_C=', dF_C, ' dF_sigma=', dF_sigma, ' dPhi_C=', dPhi_C, &
                ' ||dC||=', dC_norm, ' rms=', rms
        else
           write(*,'(a18,es11.4,a12,es11.4,a12,es11.4,a12,es11.4,a10,es11.4)') &
                ' DistDiag: dF_C=', dF_C, ' dF_sigma=', dF_sigma, ' share_C=', share_C, &
                ' ||dC||=', dC_norm, ' rms=', rms
           write(ioLog,'(a18,es11.4,a12,es11.4,a12,es11.4,a12,es11.4,a10,es11.4)') &
                ' DistDiag: dF_C=', dF_C, ' dF_sigma=', dF_sigma, ' share_C=', share_C, &
                ' ||dC||=', dC_norm, ' rms=', rms
        end if

        if (stopInversion .or. (rms < distControl%rmsTol) .or. (iter >= distControl%maxIter)) exit blockLoop
    end do blockLoop

    ! Final output
    call CmSqrtMult(mHat, m_minus_m0)
    call linComb(ONE, m_minus_m0, ONE, m0, m)
    d = dHat
    write(*,'(a25,i5,a30,i5)') 'NLCG_DIST iterations:', iter, ' sigma function evaluations:', nfuncSigma
    write(ioLog,'(a25,i5,a30,i5)') 'NLCG_DIST iterations:', iter, ' sigma function evaluations:', nfuncSigma
    write(*,'(a30,i5,a8,i5,a1)') '       C function evaluations:', nfuncC, ' (total:', nfuncSigma + nfuncC, ')'
    write(ioLog,'(a30,i5,a8,i5,a1)') '       C function evaluations:', nfuncC, ' (total:', nfuncSigma + nfuncC, ')'
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
    call deall_distortionParam(distC_prev)
    call deall_distortionParam(distC_diff)
    call ModEM_timers_destroy('NLCG_DIST Iteration')

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
