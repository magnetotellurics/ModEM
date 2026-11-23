module NLCG_dist
! Nonlinear conjugate gradient solver for the joint inversion of the
! conductivity model and frequency-independent galvanic distortion matrices
! (Full_Impedance_Dist data, D = C Z).
!
! The distortion matrices are eliminated exactly: after every forward solve
! (including every line-search trial) C is replaced by the minimiser of the
! penalty functional for the current EM fields, which is a small per-site
! linear least-squares problem (see DistortGradient).  NLCG then minimises the
! reduced functional
!   Fr(mHat) = min_C [ SS/Ndata + lambda ||mHat||^2/Nmodel + nu psi(C) ],
! whose gradient is the conductivity gradient evaluated at the optimal C.
! The lambda strategy is the same as in the standard NLCG.
!
! User controls (in addition to the standard NLCG control lines):
!   nu                     dimensionless distortion damping factor
!   distortion output file
!   reference distortion   C0 in psi(C) = 1/2 sum_j w_j ||C_j - C0_j||^2
!                          (identity if not given)
!   distortion type        full (2x2 C, default) or diag (C = diag(gx,gy))

  use ModEM_timers
  use invcore
  use DistortGradient, only: DIST_FULL, DIST_DIAG, init_dist_penalty, &
                             deall_dist_penalty, psi_dist
  use utilities

  implicit none

  private
  public :: NLCGsolver_dist

  type :: NLCGDistControl_t
     ! standard NLCG controls
     integer            :: maxIter
     real(kind=prec)    :: rmsTol
     real(kind=prec)    :: fdiffTol
     real(kind=prec)    :: lambda
     real(kind=prec)    :: lambdaTol
     real(kind=prec)    :: k
     real(kind=prec)    :: c
     integer            :: nCGmax
     real(kind=prec)    :: alpha_1
     real(kind=prec)    :: startdm
     character(80)      :: fname
     integer            :: nskip
     ! distortion controls
     real(kind=prec)    :: nu
     character(80)      :: distFname
     character(80)      :: refDistFname
     integer            :: distType
  end type NLCGDistControl_t

  type(NLCGDistControl_t), save :: distControl

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
    ctrl%nCGmax = 8
    ctrl%alpha_1 = 20.
    ctrl%startdm = 20.
    ctrl%fname = 'Modular'
    ctrl%nskip = 1
    ctrl%nu = 0.1_prec
    ctrl%distFname = 'distortion.dat'
    ctrl%refDistFname = ''
    ctrl%distType = DIST_FULL
  end subroutine set_NLCGDistControl

  ! ---------------------------------------------------------------------------
  ! read_NLCGDistControl reads the standard NLCG control lines (positional,
  ! as in NLCG) followed by optional distortion lines identified by their
  ! label (the text before the first ':'), in any order:
  !   Distortion damping factor nu       : 0.1
  !   Distortion output filename         : distortion.dat
  !   Reference distortion filename      : (empty = identity)
  !   Distortion type (full|diag)        : full
  ! ---------------------------------------------------------------------------
  subroutine read_NLCGDistControl(ctrl, rFile, hasCtrlFile)
    type(NLCGDistControl_t), intent(inout) :: ctrl
    character(*), intent(in) :: rFile
    logical, intent(out) :: hasCtrlFile
    integer :: ios, iColon
    character(80) :: string
    character(256) :: line, label, val

    call set_NLCGDistControl(ctrl)

    inquire(FILE=rFile, EXIST=hasCtrlFile)
    if (.not. hasCtrlFile) return

    write(*,*) 'Reading inverse configuration from file ', trim(rFile)
    open(unit=ioInvCtrl, file=rFile, status='old', iostat=ios)
    if (ios /= 0) call errStop('Error opening file: '//trim(rFile))

    read(ioInvCtrl, '(a36,a80)') string, ctrl%fname
    ctrl%fname = adjustl(ctrl%fname)
    read(ioInvCtrl, '(a36,g15.7)') string, ctrl%lambda
    read(ioInvCtrl, '(a36,g15.7)') string, ctrl%k
    read(ioInvCtrl, '(a36,g15.7)') string, ctrl%startdm
    read(ioInvCtrl, '(a36,g15.7)') string, ctrl%fdiffTol
    read(ioInvCtrl, '(a36,g15.7)') string, ctrl%rmsTol
    read(ioInvCtrl, '(a36,g15.7)') string, ctrl%lambdaTol
    read(ioInvCtrl, '(a36,i4)') string, ctrl%maxIter

    do
       read(ioInvCtrl, '(a)', iostat=ios) line
       if (ios /= 0) exit
       if (len_trim(line) == 0) cycle
       iColon = index(line, ':')
       if (iColon == 0) then
          write(0,*) 'Warning: ignoring control line without label: ', &
    &        trim(line)
          cycle
       end if
       label = lowercase(adjustl(line(1:iColon-1)))
       val = adjustl(line(iColon+1:))

       ! try using a label-based parameter (should do that with the entire 
       ! modem control files)
       if (index(label, 'damping factor nu') > 0) then
          read(val, *, iostat=ios) ctrl%nu
          if (ios /= 0 .or. ctrl%nu < R_ZERO) &
               call errStop('Invalid distortion damping factor nu: '//&
    &          trim(val))
       elseif (index(label, 'distortion output') > 0) then
          if (len_trim(val) > 0) ctrl%distFname = val(1:80)
       elseif (index(label, 'reference distortion') > 0 .or. &
               index(label, 'initial distortion') > 0) then
          ctrl%refDistFname = val(1:80)
       elseif (index(label, 'distortion type') > 0) then
          select case (trim(lowercase(val)))
          case ('full', '')
             ! full matrix of C
             ctrl%distType = DIST_FULL
          case ('diag')
             ! one gain per electric field component (static shift/split)
             ctrl%distType = DIST_DIAG
          case default
             call errStop('Unknown distortion type (try full or diag): '//&
    &                trim(val))
          end select
       else
          write(0,*) 'Warning: unknown NLCG_DIST control: ', &
               trim(line(1:iColon-1))
       end if
    end do

    close(ioInvCtrl)
  end subroutine read_NLCGDistControl

  function lowercase(s) result(r)
    character(*), intent(in) :: s
    character(len(s)) :: r
    integer :: i, ic
    r = s
    do i = 1, len(s)
       ic = iachar(s(i:i))
       if (ic >= iachar('A') .and. ic <= iachar('Z')) r(i:i) = achar(ic + 32)
    end do
  end function lowercase

  ! ---------------------------------------------------------------------------
  subroutine NLCGsolver_dist(d, lambda, m0, m, distC, fname)
    type(dataVectorMTX_t), intent(inout) :: d
    real(kind=prec), intent(inout) :: lambda
    type(modelParam_t), intent(in) :: m0
    type(modelParam_t), intent(inout) :: m
    type(distortionParam_t), intent(inout) :: distC
    character(*), intent(in), optional :: fname

    type(dataVectorMTX_t) :: dHat, res
    type(modelParam_t) :: mHat, m_minus_m0, grad, g, h, gPrev
    type(solnVectorMTX_t) :: eAll
    type(distortionParam_t) :: distRef
    real(kind=prec) :: value, valuePrev, rms, rmsPrev, alpha, beta, gnorm
    real(kind=prec) :: mNorm, distNorm, Nmodel, nu, startdm, lambdaPrev
    real(kind=prec) :: grad_dot_h, g_dot_g, g_dot_gPrev, gPrev_dot_gPrev, &
   &        g_dot_h
    integer :: iter, nCG, nLS, nfunc, ios
    logical :: hasCtrlFile
    character(3) :: iterChar
    character(100) :: mFile, mHatFile, dataFile, resFile, logFile, distFile

    hasCtrlFile = .false.
    if (present(fname)) then
       call read_NLCGDistControl(distControl, fname, hasCtrlFile)
    else
       call set_NLCGDistControl(distControl)
    end if
    if (hasCtrlFile) lambda = distControl%lambda
    nu = distControl%nu
    startdm = distControl%startdm
    alpha = distControl%alpha_1

    logFile = trim(distControl%fname)//'_NLCG_DIST.log'
    open(unit=ioLog, file=logFile, status='unknown', position='append', &
   &        iostat=ios)
    distFile = trim(distControl%fname)//'_'//trim(distControl%distFname)

    ! reference distortion C0 (identity unless a file is given) and the
    ! data-based penalty weights
    call copy_distortionParam(distC, distRef)
    call identity_distortionParam(distRef)
    if (len_trim(distControl%refDistFname) > 0) then
       call read_distortionParam(trim(distControl%refDistFname), distRef)
    end if
    call init_dist_penalty(d, distRef, distControl%distType)
    call copy_distortionParam(distRef, distC)

    call write_header(ioLog)
    if (output_level > 0) call write_header(6)

    mHat = m
    call func_dist(lambda, nu, d, m0, mHat, distC, value, mNorm, distNorm, &
   &        dHat, eAll, rms)
    call printf('START', lambda, alpha, value, mNorm, rms)
    call printf('START', lambda, alpha, value, mNorm, rms, logFile)
    call print_dist_state('INIT', distC)
    nfunc = 1

    write(iterChar,'(i3.3)') 0
    call CmSqrtMult(mHat, m_minus_m0)
    call linComb(ONE, m_minus_m0, ONE, m0, m)
    if (output_level > 1) then
       mFile = trim(distControl%fname)//'_NLCG_DIST_'//iterChar//'.rho'
       call write_modelParam(m, trim(mFile))
    end if
    if (output_level > 2) then
       dataFile = trim(distControl%fname)//'_NLCG_DIST_'//iterChar//'.dat'
       call write_dataVectorMTX(dHat, trim(dataFile))
    end if
    call write_distortionParam(trim(distFile), distC)

    call gradient_dist(lambda, d, m0, mHat, distC, grad, dHat, eAll)
    gnorm = sqrt(dotProd(grad, grad))
    write(*,'(a42,es12.5)') '    GRAD: initial norm of the gradient is', gnorm
    write(ioLog,'(a42,es12.5)') '    GRAD: initial norm of the gradient is', &
   &        gnorm
    if (gnorm < TOL6) call errStop('Problem with gradient computations: first gradient is zero')
    alpha = startdm / gnorm

    nCG = 0
    iter = 0
    call linComb(MinusONE, grad, R_ZERO, grad, g)
    h = g
    call ModEM_timers_create("NLCG_DIST Iteration", .false.)

    do
       if ((rms < distControl%rmsTol) .or. (iter >= distControl%maxIter)) exit
       call ModEM_timers_start("NLCG_DIST Iteration")
       iter = iter + 1
       rmsPrev = rms
       valuePrev = value
       grad_dot_h = dotProd(grad, h)

       write(*,'(a23)') 'Starting line search...'
       write(ioLog,'(a23)') 'Starting line search...'
       call lineSearchCubic(lambda, nu, d, m0, mHat, distC, h, alpha, value, &
            grad, rms, nLS, dHat, eAll)
       nfunc = nfunc + nLS

       ! gradient of the reduced functional at the accepted (mHat, C)
       call gradient_dist(lambda, d, m0, mHat, distC, grad, dHat, eAll)
       gPrev = g
       call linComb(MinusONE, grad, R_ZERO, grad, g)

       if (value < valuePrev) then
          alpha = 2*(value - valuePrev)/grad_dot_h
          alpha = (ONE + 0.01_prec)*alpha
       end if
       write(*,'(a30,i5)') 'Completed NLCG_DIST iteration ', iter
       write(ioLog,'(a30,i5)') 'Completed NLCG_DIST iteration ', iter
       call ModEM_timers_stop("NLCG_DIST Iteration", .false.)
       call ModEM_timers_print("NLCG_DIST Iteration", ioLog)

       Nmodel = countModelParam(mHat)
       mNorm = dotProd(mHat, mHat)/Nmodel
       call printf('with', lambda, alpha, value, mNorm, rms)
       call printf('with', lambda, alpha, value, mNorm, rms, logFile)
       call print_dist_state('with', distC)

       if (mod(iter, distControl%nskip) == 0) then
          call CmSqrtMult(mHat, m_minus_m0)
          call linComb(ONE, m_minus_m0, ONE, m0, m)
          write(iterChar,'(i3.3)') iter
          if (output_level > 1) then
             mFile = trim(distControl%fname)//'_NLCG_DIST_'//iterChar//'.rho'
             call write_modelParam(m, trim(mFile))
          end if
          if (output_level > 2) then
             mHatFile = trim(distControl%fname)//'_NLCG_DIST_'//iterChar//'.prm'
             call write_modelParam(mHat, trim(mHatFile))
             dataFile = trim(distControl%fname)//'_NLCG_DIST_'//iterChar//'.dat'
             call write_dataVectorMTX(dHat, trim(dataFile))
             res = d
             call linComb(ONE, d, MinusONE, dHat, res)
             resFile = trim(distControl%fname)//'_NLCG_DIST_'//iterChar//'.res'
             call write_dataVectorMTX(res, trim(resFile))
          end if
          call write_distortionParam(trim(distFile), distC)
       end if

       if (abs(rmsPrev - rms) < distControl%fdiffTol) then
          ! The optimal C depends on nu and the EM fields only, not on lambda,
          ! so only the model-norm terms of F and its gradient change.
          lambdaPrev = lambda
          lambda = lambda / distControl%k
          if (lambda < distControl%lambdaTol) then
             write(*,'(a55)') 'Unable to get out of a local minimum. Exiting...'
             write(ioLog,'(a55)') 'Unable to get out of a local minimum. Exiting...'
             exit
          end if
          value = value + (lambda - lambdaPrev)*mNorm
          call linComb(ONE, grad, TWO*(lambda - lambdaPrev)/Nmodel, mHat, grad)
          gnorm = sqrt(dotProd(grad, grad))
          alpha = min(ONE, startdm)/gnorm
          call linComb(MinusONE, grad, R_ZERO, grad, g)
          write(*,'(a55)') 'Restarting NLCG with the damping parameter updated'
          write(ioLog,'(a55)') 'Restarting NLCG with the damping parameter updated'
          call printf('to', lambda, alpha, value, mNorm, rms)
          call printf('to', lambda, alpha, value, mNorm, rms, logFile)
          h = g
          nCG = 0
          cycle
       end if

       g_dot_g = dotProd(g, g)
       g_dot_gPrev = dotProd(g, gPrev)
       gPrev_dot_gPrev = dotProd(gPrev, gPrev)
       g_dot_h = dotProd(g, h)
       beta = (g_dot_g - g_dot_gPrev)/gPrev_dot_gPrev
       if ((g_dot_g + beta*g_dot_h <= R_ZERO) .and. (nCG >= distControl%nCGmax)) then
          write(*,'(a45)') 'Restarting NLCG to restore orthogonality'
          write(ioLog,'(a45)') 'Restarting NLCG to restore orthogonality'
          nCG = 0
          beta = R_ZERO
       else
          nCG = nCG + 1
       end if
       call linComb(ONE, g, beta, h, h)
    end do

    call CmSqrtMult(mHat, m_minus_m0)
    call linComb(ONE, m_minus_m0, ONE, m0, m)
    d = dHat
    write(*,'(a25,i5,a25,i5)') 'NLCG_DIST iterations:', iter, &
   &        ' function evaluations:', nfunc
    write(ioLog,'(a25,i5,a25,i5)') 'NLCG_DIST iterations:', iter, &
   &        ' function evaluations:', nfunc
    call write_distortionParam(trim(distFile), distC)
    close(ioLog, iostat=ios)

    call deall_dataVectorMTX(dHat)
    call deall_dataVectorMTX(res)
    call deall_modelParam(mHat)
    call deall_modelParam(m_minus_m0)
    call deall_modelParam(grad)
    call deall_modelParam(g)
    call deall_modelParam(h)
    call deall_modelParam(gPrev)
    call deall_solnVectorMTX(eAll)
    call deall_distortionParam(distRef)
    call deall_dist_penalty()
    call ModEM_timers_destroy('NLCG_DIST Iteration')

  contains

    subroutine write_header(iu)
      integer, intent(in) :: iu
      write(iu,'(a41,es10.3)') 'The initial damping parameter lambda is ', &
   &          lambda
      write(iu,'(a41,es10.3)') 'The distortion damping factor nu is     ', nu
      select case (distControl%distType)
      case (DIST_DIAG)
         write(iu,'(a41,a)') 'Distortion type:                        ', &
   &         ' diag (C = diag(gx,gy))'
      case default
         write(iu,'(a41,a)') 'Distortion type:                        ', &
   &         ' full tensor'
      end select
      if (len_trim(distControl%refDistFname) > 0) then
         write(iu,'(a42,a)') 'Reference distortion read from           ', &
              trim(distControl%refDistFname)
      else
         write(iu,'(a41,a)') 'Reference distortion:                   ', &
   &          ' C = I'
      end if
      write(iu,'(a42,a)') 'Distortion output file:                  ',    &
   &          trim(distFile)
    end subroutine write_header

  end subroutine NLCGsolver_dist

  !**********************************************************************
  ! one-line summary of the distortion state: penalty, mean and spread of
  ! the site log-gain (1/2 log|det C|) and the array-mean C
  subroutine print_dist_state(comment, distC)
    character(*), intent(in) :: comment
    type(distortionParam_t), intent(in) :: distC
    real(kind=prec) :: lg, s1, s2, Cbar(2,2), det
    integer :: j, n

    s1 = R_ZERO; s2 = R_ZERO; Cbar = R_ZERO; n = 0
    do j = 1, distC%nSites
       det = abs(distC%C(1,1,j)*distC%C(2,2,j) - distC%C(1,2,j)*distC%C(2,1,j))
       Cbar = Cbar + distC%C(:,:,j)
       if (det > R_TINY) then
          lg = 0.5_prec*log(det)
          s1 = s1 + lg; s2 = s2 + lg*lg; n = n + 1
       end if
    end do
    if (distC%nSites > 0) Cbar = Cbar / distC%nSites
    if (n > 0) then
       s1 = s1/n; s2 = sqrt(max(s2/n - s1*s1, R_ZERO))
    end if
    write(*,'(a10,a5,es10.3,a12,f8.4,a7,f7.4,a7,4f7.3)') & 
   &      trim(comment)//' DIST:', &
   &     'psi= ', psi_dist(distC), ' log(Gain)=', s1, ' sd=', s2, &
   &     'Cbar=', Cbar(1,1), Cbar(1,2), Cbar(2,1), Cbar(2,2)
    write(ioLog,'(a10,a5,es10.3,a12,f8.4,a7,f7.4,a7,4f7.3)')  &
   &      trim(comment)//' DIST:', &
   &     'psi=', psi_dist(distC), ' log(Gain)=', s1, ' sd=', s2, &
   &     'Cbar=', Cbar(1,1), Cbar(1,2), Cbar(2,1), Cbar(2,2)
    flush(ioLog)
  end subroutine print_dist_state

  !**********************************************************************
  ! Cubic line search on the reduced functional (same strategy as the
  ! standard NLCG lineSearchCubic).  Every trial re-eliminates C; the C of
  ! the accepted trial is returned in distC.  A step that increases the
  ! functional is never accepted: the starting point is restored and the
  ! step length reduced instead.
  subroutine lineSearchCubic(lambda, nu, d, m0, mHat, distC, h, alpha, f, grad, &
       rms, niter, dHat, eAll)
    real(kind=prec), intent(in) :: lambda, nu
    type(dataVectorMTX_t), intent(in) :: d
    type(modelParam_t), intent(in) :: m0
    type(modelParam_t), intent(inout) :: mHat
    type(distortionParam_t), intent(inout) :: distC
    type(modelParam_t), intent(in) :: h
    real(kind=prec), intent(inout) :: alpha
    real(kind=prec), intent(inout) :: f
    type(modelParam_t), intent(in) :: grad
    real(kind=prec), intent(out) :: rms
    integer, intent(out) :: niter
    type(dataVectorMTX_t), intent(inout) :: dHat
    type(solnVectorMTX_t), intent(inout) :: eAll

    real(kind=prec) :: alpha_1, alpha_i, alpha_j, mNorm, distNorm
    real(kind=prec) :: c, a, b, q1, q2, q3, g_0, f_0, f_1, f_i, f_j, rms_1
    real(kind=prec) :: mNorm_1, distNorm_1
    type(modelParam_t) :: mHat_0, mHat_1
    type(dataVectorMTX_t) :: dHat_1
    type(solnVectorMTX_t) :: eAll_1
    type(distortionParam_t) :: distC_0, distC_1
    logical :: useFirst, fail
    character(100) :: logFile

    c = distControl%c
    logFile = trim(distControl%fname)//'_NLCG_DIST.log'
    niter = 0
    mHat_0 = mHat
    f_0 = f
    call copy_distortionParam(distC, distC_0)
    call copy_distortionParam(distC, distC_1)
    g_0 = dotProd(grad, h)
    alpha_1 = alpha
    useFirst = .false.
    fail = .false.

    mHat_1 = mHat_0
    call linComb(ONE, mHat_0, alpha_1, h, mHat_1)
    call func_dist(lambda, nu, d, m0, mHat_1, distC_1, f_1, mNorm_1, distNorm_1, dHat_1, eAll_1, rms_1)
    call printf('STARTLS', lambda, alpha, f_1, mNorm_1, rms_1)
    call printf('STARTLS', lambda, alpha, f_1, mNorm_1, rms_1, logFile)
    niter = niter + 1

    a = (f_1 - f_0 - g_0*alpha_1)/(alpha_1**2)
    b = g_0
    if (a < 0) then
       ! the quadratic has no minimum: accept the first trial if it decreases f
       if (f_1 < f_0) then
          useFirst = .true.
       else
          fail = .true.
       end if
       write(*,'(a45)') 'Quadratic has no minimum, exiting line search'
       write(ioLog,'(a45)') 'Quadratic has no minimum, exiting line search'
    else
       alpha = -b/(TWO*a)
       call linComb(ONE, mHat_0, alpha, h, mHat)
       call copy_distortionParam(distC_0, distC)
       call func_dist(lambda, nu, d, m0, mHat, distC, f, mNorm, distNorm, dHat, eAll, rms)
       call printf('QUADLS', lambda, alpha, f, mNorm, rms)
       call printf('QUADLS', lambda, alpha, f, mNorm, rms, logFile)
       niter = niter + 1

       if (.not. (f < f_0 + c*alpha*g_0)) then
          if (f <= f_0) then
             alpha_i = alpha_1; f_i = f_1
             alpha_j = alpha;   f_j = f
             fit_cubic: do
                if (niter >= 5) exit fit_cubic
                q1 = f_i - f_0 - g_0*alpha_i
                q2 = f_j - f_0 - g_0*alpha_j
                q3 = alpha_i**2*alpha_j**2*(alpha_j - alpha_i)
                a = (alpha_i**2*q2 - alpha_j**2*q1)/q3
                b = (alpha_j**3*q1 - alpha_i**3*q2)/q3
                if (b*b - 3*a*g_0 < R_ZERO .or. abs(a) < R_TINY) exit fit_cubic
                alpha = (-b + sqrt(b*b - 3*a*g_0))/(3*a)
                call linComb(ONE, mHat_0, alpha, h, mHat)
                call copy_distortionParam(distC_0, distC)
                call func_dist(lambda, nu, d, m0, mHat, distC, f, mNorm, distNorm, dHat, eAll, rms)
                call printf('CUBICLS', lambda, alpha, f, mNorm, rms)
                call printf('CUBICLS', lambda, alpha, f, mNorm, rms, logFile)
                niter = niter + 1
                if (f < f_0 + c*alpha*g_0) exit fit_cubic
                alpha_i = alpha_j; f_i = f_j
                alpha_j = alpha;   f_j = f
                if (abs(f_j - f_i) < TOL8) exit fit_cubic
             end do fit_cubic
          end if
       end if
       if (f_1 < f) useFirst = .true.
       if ((.not. useFirst) .and. (f >= f_0)) fail = .true.
       if (useFirst .and. (f_1 >= f_0)) then
          useFirst = .false.
          fail = .true.
       end if
    end if

    if (useFirst) then
       alpha = alpha_1
       dHat = dHat_1
       eAll = eAll_1
       mHat = mHat_1
       rms = rms_1
       f = f_1
       call copy_distortionParam(distC_1, distC)
    else if (fail) then
       ! keep the starting point; the fields must be recomputed (only on
       ! rejection) so that the next gradient matches mHat_0
       write(*,'(a60)') 'Line search failed to decrease f; keeping previous model'
       write(ioLog,'(a60)') 'Line search failed to decrease f; keeping previous model'
       mHat = mHat_0
       call copy_distortionParam(distC_0, distC)
       call func_dist(lambda, nu, d, m0, mHat, distC, f, mNorm, distNorm, dHat, eAll, rms)
       niter = niter + 1
       alpha = 0.25_prec*min(alpha, alpha_1)
    end if

    call deall_dataVectorMTX(dHat_1)
    call deall_modelParam(mHat_0)
    call deall_modelParam(mHat_1)
    call deall_solnVectorMTX(eAll_1)
    call deall_distortionParam(distC_0)
    call deall_distortionParam(distC_1)
  end subroutine lineSearchCubic

end module NLCG_dist
