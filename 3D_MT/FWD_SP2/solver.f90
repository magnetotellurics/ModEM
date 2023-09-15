! *****************************************************************************
! module containing iterative equation solvers. Uses operators and
! pre-conditioners defined in SG3DFWC1 to solve equations for divergence
! correction, induction operator. 
! Source code is completely general; only the
! module interface is  specific to implementation of operators in SG3DFWC1
! 
! modified to translate cvectors into normal arrays to be used with sparse 
! matrix operation
module solver

   use math_constants   ! math/ physics constants
   use utilities, only: isnan
   use spoptools        ! for sparse-matrix operations
#ifdef CUDA
   use cudaFortMap      ! c GPU api bindings for fortran
#endif
   !use griddef	! staggered grid definitions
   !use sg_scalar
   !use sg_vector
   implicit none

  ! solverControl_t as a data type controls diagnostic tools for joint forward
  ! modeling-inversion scheme
  type  :: solverControl_t
     ! maximum number of iterations in one call to iterative solver
     integer                                               :: maxIt
     ! convergence criteria: return from solver if relative error < tol
     real (kind=prec)                                      :: tol
     ! actual number of iterations before return
     integer                                               :: niter
     ! relative error for each iteration
     real (kind=prec), pointer, dimension(:)               :: rerr
     ! logical variable indicating if algorithm "failed"
     logical                                               :: failed = .false.
  end type solverControl_t


Contains


! *****************************************************************************
! Solver contains subroutines for: 
! a) PCG- a quasi-generic pre-conditioned congugate gradient, and 
! b) QMR - Quasi-Minimal Residual method (pre-conditioned, no look-ahead)
! c) BICG - bicgstab Stablilized Biconjugate Gradient method
! d) TFQMR - Transpose-free QMR (don't use it with CC-DC, as TFQMR doesn't 
!            work well with frequent interuption) 
! *****************************************************************************

subroutine PCG(b,x,PCGiter)
  ! Purpose: a quasi-generic pre-conditioned conjugate gradient
  ! routine.
  ! solves Ax = b
  ! 
  ! modified to translate cvectors/scalars into normal arrays to be used with 
  ! sparse matrix operation. That said, the implementation is still the old way
  ! used in matrix free method:
  ! the real A is already defined in the DivCgrad function for divergence 
  ! correction.
    use modeloperator3d, only: A => DivCgrad, Minv => DivCgradILU

  implicit none
  complex (kind=prec), intent(in), dimension(:)    :: b
  complex (kind=prec), intent(inout),dimension(:)  :: x
  type (solverControl_t), intent(inout)            :: PCGiter

  ! local variables
  complex(kind=prec),allocatable,dimension(:)  :: r,s,p,q
  complex(kind=prec)    :: beta,alpha,delta,deltaOld
  real(kind=prec)       :: bnorm, rnorm
  integer               :: i
  integer               :: xsize


  xsize=size(x,1)
  ! Allocation of r, z, p, q
  ! residual
  allocate(r(xsize))
  allocate(s(xsize))
  allocate(p(xsize))
  allocate(q(xsize))

  ! r = b-A*x
  call A(x,r) ! in DivCorr, x should be zero (no initial guess)
  r = b - r   ! hence r should be indentical to b
  bnorm = sqrt(dot_product(b,b))
  if (isnan(bnorm)) then
  ! this usually means an inadequate model, in which case Maxwell's fails
      write(0,*) 'Error: b in PCG contains NaNs; exiting...'
      stop
  endif
  rnorm = sqrt(dot_product(r,r))
  i = 1
  PCGiter%rerr(i) =real(rnorm/bnorm)

  loop: do while ((PCGiter%rerr(i).gt.PCGiter%tol).and.(i.lt.PCGiter%maxIt))
     ! precondiction first
     call Minv(r,s)
     delta = dot_product(r,s)
     if(i.eq.1) then
        p = s
     else
        beta = delta/deltaOld
        p = s + beta*p
     end if
     call A(p,q)
     alpha = delta/dot_product(p,q)
     x = x + p * alpha
     r = r - q * alpha
     deltaOld = delta
     i = i + 1
     rnorm = sqrt(dot_product(r,r))
     if (isnan(rnorm)) then
         write(0,*) 'Error: residual in PCG contains NaNs; exiting...'
         stop
     endif
     PCGiter%rerr(i) = real(rnorm/bnorm)
  end do loop

  PCGiter%niter = i

  ! explicitly deallocate the temperary work arrays
  ! seems not necesary in fortran
  deallocate(r)
  deallocate(s)
  deallocate(p)
  deallocate(q)

end subroutine PCG ! PCG


! *****************************************************************************
subroutine QMR(b,x,QMRiter)
  ! Purpose ... a quasi-minimal residual method routine, set up for solving
  ! A x = b using routines in  mult_A. Actual code is generic, but interface
  ! is fairly specific

  ! interface...........
  ! redefining some of the interfaces for our convenience (locally)
  ! generic routines for vector operations for edge/ face nodes
  ! in a staggered grid
  ! 
  ! modified to be used with sparse matrix operators
  ! Again this is not that generic as the A is already defined in Mult_Aii... 
  use modeloperator3d, only: A => mult_Aii, M1solve => PC_Lsolve,        &
     &                         M2solve => PC_Usolve

  implicit none
  !  b is right hand side
  complex (kind=prec), intent(in), dimension(:)    :: b
  !  solution vector is x ... on input is provided with the initial
  !  guess, on output is the most recent iterate
  complex (kind=prec), intent(inout),dimension(:)  :: x
  type (solverControl_t), intent(inout) :: QMRiter

   ! local variables
  complex(kind=prec),allocatable,dimension(:)  :: XMIN,R,VT
  complex(kind=prec),allocatable,dimension(:)  :: Y,Z,WT,V,W,YT,ZT
  complex(kind=prec),allocatable,dimension(:)  :: P,Q,PT,D,S
  complex (kind=prec)          :: ETA,PDE,EPSIL,RDE,BETA,DELTA,RHO
  complex (kind=prec)          :: PSI,RHO1,GAMM,GAMM1,THET,THET1,TM2
  real (kind=prec)             :: bnorm,rnorm,rnormin
  complex (kind=prec)          :: rhoInv,psiInv
  integer                      :: iter, xsize
  logical                      :: adjoint, ilu_adjt

  xsize=size(x,1)
  ! Allocate work arrays
  allocate(XMIN(xsize))
  allocate(R(xsize))
  allocate(VT(xsize))
  allocate(Y(xsize))
  allocate(Z(xsize))
  allocate(WT(xsize))
  allocate(V(xsize))
  allocate(W(xsize))
  allocate(YT(xsize))
  allocate(ZT(xsize))
  allocate(P(xsize))
  allocate(Q(xsize))
  allocate(PT(xsize))
  allocate(D(xsize))
  allocate(S(xsize))

  ! NOTE: this iterative solver is QMR without look-ahead
  ! patterned after the scheme given on page 24 of Barrett et al.
  ! "Templates for the solution of linear systems of equations:
  ! Building blocks for iterative methods"
  ! Note that there are a couple of small differences, due to
  ! the fact that our system is complex

  ! R=b-Ax 
  adjoint = .false.
  call A(x,adjoint,R)
  ! R= b - Ax, for inital guess x, that has been inputted to the routine
  R = b - R

  ! Norm of rhs, residual
  bnorm = CDSQRT(dot_product(b, b))
  rnorm = CDSQRT(dot_product(R, R))
  rnormin = rnorm
  XMIN = x

  ! this usually means an inadequate model, in which case Maxwell's fails
  if (isnan(abs(bnorm))) then
      write(0,*) 'Error: b in QMR contains NaNs; exiting...'
      stop
  endif

  !  iter is iteration counter
  iter = 1
  QMRiter%rerr(iter) = real(rnorm/bnorm)
  ! write(6,*) 'initial residual:', QMRiter%rerr(iter) 

  ! L
  VT = R
  ilu_adjt = .false.
  call M1solve(VT,ilu_adjt,Y)
!  Y = VT
  RHO = CDSQRT(dot_product(Y,Y))
  ! U
  WT = R
  ilu_adjt = .true.
  call M2solve(WT,ilu_adjt,Z)
!  Z = WT
  PSI  = CDSQRT(dot_product(Z,Z))
  GAMM = C_ONE
  ETA  = C_MinusONE

  ! the do loop goes on while the relative error is greater than the tolerance
  ! and the iterations are less than maxIt
  loop: do while ((QMRiter%rerr(iter).gt.QMRiter%tol).and.&
       (iter.lt.QMRiter%maxIt))
      if ((RHO.eq.C_ZERO).or.(PSI.eq.C_ZERO)) then
        QMRiter%failed = .true.
        write(0,*) 'QMR FAILED TO CONVERGE : RHO'
        write(0,*) 'QMR FAILED TO CONVERGE : PSI'
        stop
      endif

      rhoInv = (1/RHO)*cmplx(1.0, 0.0, 8)
      psiInv = (1/PSI)*cmplx(1.0, 0.0, 8)
      V = VT * rhoInv
      W = WT * psiInv
      Y = Y * rhoinv
      Z = Z * psiinv

      DELTA = dot_product(Z,Y)
      if(DELTA.eq.C_ZERO) then
        QMRiter%failed = .true.
        write(0,*) 'QMR FAILS TO CONVERGE : DELTA'
        exit
      endif

      ilu_adjt = .false.
      call M2solve(Y,ilu_adjt,YT)
!      YT = Y
      ilu_adjt = .true.
      Call M1solve(Z,ilu_adjt,ZT)
!      ZT = Z

      if (iter.eq.1) then
        P = YT
        Q = ZT
      else
      ! these calculations are only done when iter greater than 1
        PDE = PSI*DELTA/EPSIL
        RDE = RHO*CONJG(DELTA/EPSIL)
        P = YT - PDE * P
        Q = ZT - RDE * Q
      endif

      adjoint = .false.
      Call A(P, adjoint, PT)
      EPSIL = dot_product(Q,PT)
      if (EPSIL.eq.C_ZERO) then
        QMRiter%failed = .true.
        write(0,*) 'QMR FAILED TO CONVERGE : EPSIL'
        exit
      endif

      BETA = EPSIL/DELTA
      if (BETA.eq.C_ZERO) then
        QMRiter%failed = .true.
        write(0,*) 'QMR FAILED TO CONVERGE : BETA'
        exit
      endif
      VT = PT - BETA * V

      RHO1 = RHO
      ilu_adjt = .false.
      Call M1solve(VT, ilu_adjt, Y)
!      Y = VT
      RHO = CDSQRT(dot_product(Y,Y))

      adjoint = .true.
      Call A(Q, adjoint, WT)
      WT = WT - conjg(BETA) * W

      ilu_adjt = .true.
      Call M2solve(WT,ilu_adjt,Z)
!      Z = WT
      PSI = CDSQRT(dot_product(Z,Z))

      if (iter.gt.1) then
        THET1 = THET
      endif
      THET = RHO/(GAMM*CDABS(BETA))
      GAMM1 = GAMM
      GAMM = C_ONE/CDSQRT(C_ONE + THET*THET)
      if (GAMM.eq.C_ZERO) then
        write(0,*) 'QMR FAILS TO CONVERGE : GAMM'
        exit
      endif

      ETA = -ETA*RHO1*GAMM*GAMM/(BETA*GAMM1*GAMM1)
      if (iter.eq.1) then
        D = ETA * P
        S = ETA * PT
      else
        TM2 = THET1*THET1*GAMM*GAMM
        D =  ETA*P + TM2*D
        S =  ETA*PT + TM2*S
      endif

      x = x + D
      R = R - S
      rnorm = CDSQRT(dot_product(R,R))
      if (rnorm .lt. rnormin) then! update the best solution so far
         rnormin = rnorm
         XMIN = X
      end if
      iter = iter + 1

      ! Keeping track of errors
      ! QMR book-keeping between divergence correction calls
      QMRiter%rerr(iter) = real(rnorm/bnorm)
      !write(6,*) 'iter # ', iter, 'residual:', QMRiter%rerr(iter)
  end do loop

  QMRiter%niter = iter
  ! x = XMIN ! use the last instead of the best
  QMRiter%rerr(iter) = real(rnormin/bnorm)

  ! deallocate all the work arrays
  ! seems not necesary here
  deallocate(XMIN)
  deallocate(R)
  deallocate(VT)
  deallocate(Y)
  deallocate(Z)
  deallocate(WT)
  deallocate(V)
  deallocate(W)
  deallocate(YT)
  deallocate(ZT)
  deallocate(P)
  deallocate(Q)
  deallocate(PT)
  deallocate(D)
  deallocate(S)

end subroutine QMR ! qmr

! *****************************************************************************
subroutine TFQMR(b,x,KSPiter,adjt)
  ! a transpose-free version of Quasi-Minimum Residue Algorithm,
  ! set up for solving
  ! A x = b using routines in  mult_Aii.
  ! solves for the interior (edge) field
  ! see: 
  ! Freund, Roland, A transpose-free quasi-minimal residual algorithm for
  ! non-Hermitian linear systems, SIAM J. Sci. Comp., 14 (1993), 470--482.
  !
  ! to me TFQMR is something like a middle ground between the Conjugate 
  ! Gradient and the Minimal Residual methods - the good part (like BiCG) 
  ! is this does not need the y = A^T*x operation like the original QMR
  ! therefore it is a good choice for the SP2 framework - as the CCGD 
  ! matrix is no longer symmetric 
  ! 
  ! the TFQMR also converges smoothly, which seems more stable comparing 
  ! with the BiCG (which is literally like a roller-coaster)
  ! the code is therefore easier to coupe with when the mixed precision
  ! method is considered. 
  ! 
  ! the downside of TFQMR, is that it doesn't work well when interupted 
  ! frequently (needs a long Krylov subspace) 
  ! so it won't work well with CC-DC - you have been warned
  ! 
  ! modified from my matlab version of TFQMR...
  ! so the naming might sound a little different from conventional ones
  ! also added the optional adjoint to solve adjoint system A^Tx = b 
  ! 
  ! NOTE: like BICG, TFQMR performs two sub line searches within a 
  !      iteration, but here we only store the relerr for the second sub
  !      just to be compatitive with QMR
  ! 
  ! interface...........
  ! redefining some of the interfaces for our convenience (locally)
  ! generic routines for vector operations for edge/ face nodes
  ! in a staggered grid
  !
  ! NOTE: this has not been extensively tested 
  ! if you have time reading this, test it!
  use modeloperator3d, only: A => mult_Aii, M1solve => PC_Lsolve,        &
     &                                      M2solve => PC_Usolve

  implicit none
  !  b is right hand side
  complex (kind=prec), intent(in), dimension(:)    :: b
  !  solution vector is x ... on input is provided with the initial
  !  guess, on output is the iterate with smallest residual. 
  !
  complex (kind=prec), intent(inout),dimension(:)  :: x
  type (solverControl_t), intent(inout)            :: KSPiter
  logical,intent(in),optional                      :: adjt

  ! local variables
  complex (kind=prec),allocatable,dimension(:)  :: R, R0, V, W, D
  complex (kind=prec),allocatable,dimension(:)  :: AX, AY, AD, Y, YP1
  complex (kind=prec),allocatable,dimension(:)  :: PY, PY1 ! preconditioned Y
  complex (kind=prec),allocatable,dimension(:)  :: xhalf,xmin
  real    (kind=prec)                           :: rnorm, rnorm1, bnorm,rnormin
  real    (kind=prec)                           :: xnorm, dnorm, btol
  real    (kind=prec)                           :: THETA, TAU, C
  complex (kind=prec)                           :: RHO, RHO1, ALPHA, BETA
  complex (kind=prec)                           :: ETA, SIGMA
  integer                                       :: iter, xsize, imin
  integer                                       :: maxiter
  logical                                       :: adjoint, ilu_adjt, converged

  ! stagnant check
  integer                                       :: maxstagsteps, stag
  integer                                       :: restarted, maxrestarts
  integer                                       :: last_restart, interval
  logical                                       :: restart
  real (kind=prec)                              :: eps = R_tiny
 
  if (present(adjt)) then
      adjoint = adjt
      ilu_adjt = adjt
  else
      adjoint = .false.
      ilu_adjt = .false.
  endif
  xsize = size(x,1)
  ! Norm of rhs
  bnorm = SQRT(dot_product(b, b))
  if (isnan(bnorm)) then
      write(0,*) 'Error: b in TFQMR contains NaNs; exiting...'
      stop
  else if ( bnorm .eq. 0.0) then ! zero rhs -> zero solution
      write(0,*) 'Warning: b in TFQMR has all zeros, returning zero solution'
      x = b 
      KSPiter%niter=1
      KSPiter%failed=.false.
      KSPiter%rerr=0.0
      return
  endif
  ! allocate the local variables
  allocate(R(xsize))
  ! now calculate the (original) residual
  call A(x,adjoint,R)
  ! R= b - Ax, for inital guess x, that has been input to the routine
  R = b - R
  ! Norm of residual
  rnorm = SQRT(dot_product(R, R))
  btol = KSPiter%tol * bnorm
  if ( rnorm .le. btol ) then ! the first guess is already good enough
     ! returning
      KSPiter%niter=1
      KSPiter%failed=.false.
      KSPiter%rerr(1)=real(rnorm/bnorm)
      deallocate(R)
      return 
  else
      ! allocate the rest here
      allocate(xhalf(xsize))
      allocate(xmin(xsize))
      allocate(AX(xsize))
      allocate(Y(xsize))
      allocate(YP1(xsize))
      allocate(AY(xsize))
      allocate(R0(xsize))
      allocate(D(xsize))
      allocate(AD(xsize))
      allocate(PY(xsize))
      allocate(PY1(xsize))
      allocate(W(xsize))
      allocate(V(xsize))
  end if 
!================= Now start configuring the iteration ===================!
  ! and store the best residual so far
  rnormin = rnorm
  xmin = x
  imin = 1
  KSPiter%rerr(1) = real(rnormin/bnorm)
  ! write(6,'(A, ES12.6)') ' initial relative residual = ', rnorm/bnorm
  converged = .false.
  maxiter = KSPiter%maxit 
  ! parameters for stagnant check
  stag = 0
  maxstagsteps = 5
  maxrestarts = 3
  last_restart = 0
  restarted = 0
  restart = .TRUE.
  ! interval = 120 !hard coded here
!============================== looooops! ================================!
  do iter= 1, maxIter
      ! if (mod((iter-last_restart), interval).eq.1) then
      !     restart = .true.
      ! end if
      if (restart) then
          ! write(6,'(A8,I4,A22,ES10.4)') 'iter: ', iter,'  &
          !         (re)started relres = ', rnorm/bnorm
          ! store the first residual
          R0 = R
          W = R0
          Y = W
          ! L
          call M1solve(Y,ilu_adjt,PY1)
          ! U
          call M2solve(PY1,ilu_adjt,PY)
          ! A*Y = AY
          call A(PY,adjoint,AY)
          ! V = AY
          V = AY
          D = C_ZERO
          AD = D
          THETA = 0
          ETA = C_ZERO
          TAU = rnorm
          RHO = dot_product(R, R)
          RHO1 = RHO
          restart = .FALSE.
      end if
      alpha = RHO / dot_product(R0, V)
      YP1 = Y - alpha * V
      ! =================== first half of the TFQMR iteration ===========!
      W = W - ALPHA * AY
      ! SIGMA 0 should be zero
      SIGMA = (THETA*THETA/ALPHA)*ETA
      D = PY + SIGMA * D
      AD = AY + SIGMA * AD
      THETA = SQRT(dot_product(W, W))/TAU
      C = 1.0D0 / SQRT(1.0D0 + THETA * THETA)
      TAU = TAU * THETA * C ! = norm(r)/norm(b)
      ETA = C*C * ALPHA
      xnorm = SQRT(dot_product(x,x))
      dnorm = SQRT(dot_product(D,D))
      if (abs(ETA)*dnorm .lt. eps*xnorm) then
          stag = stag + 1
      else 
          stag = 0
      end if
      ! update the xhalf - first half of the iteration
      xhalf = x + eta*d
      ! update the residual
      R = R - ETA*AD
      rnorm = SQRT(dot_product(R,R))
      KSPiter%rerr(iter)=real(rnorm/bnorm)
      ! stagnant check
      ! write(6,*) 'iter # ',iter,'first half relres= ', KSPiter%rerr(iter)
      if (rnorm.le.btol) then
          ! double check if the residual is really less than tol 
          call A(xhalf,adjoint,AX)
          R = b - AX
          rnorm = SQRT(dot_product(R,R))
          if (rnorm .le. btol) then
              x = xhalf
              KSPiter%rerr(iter)=real(rnorm/bnorm)
              KSPiter%failed = .false.
              KSPiter%niter = iter
              converged = .true.
              exit
          end if
      end if
      if (stag .ge. maxstagsteps) then
          stag = 0 ! bail out
          ! try restarting
          restarted = restarted + 1
          if (restarted .gt. maxrestarts) then
              ! stagnant - exiting
              converged = .false.
              KSPiter%failed = .true.
              KSPiter%niter = iter
              exit 
          else
              x = xhalf
              last_restart = iter
              restart = .TRUE.
              continue
          end if
      end if
      if (rnorm .lt. rnormin) then
          ! store the best solution so far
          rnormin = rnorm
          xmin = xhalf
          imin = iter
      end if
      ! L
      call M1solve(YP1,ilu_adjt,PY1)
      ! U
      call M2solve(PY1,ilu_adjt,PY)
      ! update AY
      call A(PY,adjoint,AY)
      Y = YP1
      ! ================== second half of the TFQMR iteration ===========!
      W = W - ALPHA * AY
      SIGMA = (THETA*THETA/ALPHA)*ETA
      D = PY + SIGMA * D
      AD = AY + SIGMA * AD
      THETA = SQRT(dot_product(W, W))/TAU
      C = 1.0D0 / SQRT(1.0D0 + THETA * THETA)
      TAU = TAU * THETA * C
      ETA = C*C * ALPHA
      ! stag check
      xnorm = SQRT(dot_product(xhalf,xhalf))
      dnorm = SQRT(dot_product(D,D))
      if (abs(ETA)*dnorm .lt. eps*xnorm) then
          stag = stag + 1
      else 
          stag = 0
      end if
      ! update the x - second half of the iteration
      x = xhalf + ETA*D
      ! update the residual
      R = R - ETA*AD
      rnorm = SQRT(dot_product(R,R))
      KSPiter%rerr(iter)=real(rnorm/bnorm)
      ! write(6,*) 'iter # ',iter,'second half relres= ', KSPiter%rerr(iter)
      if (rnorm.lt.btol) then
          ! double check if the residual is really less than tol 
          call A(x,adjoint,AX)
          R = b - AX
          rnorm = SQRT(dot_product(R,R))
          if (rnorm .le. btol) then
              KSPiter%rerr(iter)=real(rnorm/bnorm)
              KSPiter%failed = .false.
              KSPiter%niter = iter
              converged = .true.
              exit
          end if
      end if
      if (stag .ge. maxstagsteps) then
          stag = 0 ! bail out
          ! try restarting
          restarted = restarted + 1
          if (restarted .gt. maxrestarts) then
              ! stagnant - exiting
              converged = .false.
              KSPiter%failed = .true.
              KSPiter%niter = iter
              exit 
          else
              last_restart = iter
              restart = .TRUE.
              continue
          end if
      end if
      if (rnorm .lt. rnormin) then
          ! store the best solution so far
          rnormin = rnorm
          xmin = x
          imin = iter
      end if
      ! update the RHO
      RHO = dot_product(R0,W)
      BETA = RHO / RHO1
      ! store the previous RHO
      RHO1 = RHO
      ! update Y
      YP1 = W + BETA*Y
      ! L
      call M1solve(YP1,ilu_adjt,PY1)
      ! U
      call M2solve(PY1,ilu_adjt,PY)
      ! partial update of V
      V = BETA*(AY+BETA*V)
      ! update AY
      call A(PY,adjoint,AY)
      ! second part of update for V
      V = AY + V
      Y = YP1
  end do
 
  if (.not. converged) then 
      ! it should be noted that this is the way my matlab version works
      ! the TFQMR will return the 'best' (smallest residual) iteration
      ! KSPiter%niter=imin ! comment this line
      KSPiter%niter = maxiter
      ! KSPiter%rerr(KSPiter%maxit) = KSPiter%rerr(imin)  ! and this line
      ! to use the last iteration result instead of the 'best'
  end if
  deallocate(xhalf)
  deallocate(xmin)
  deallocate(AX)
  deallocate(R)
  deallocate(R0)
  deallocate(Y)
  deallocate(YP1)
  deallocate(PY)
  deallocate(PY1)
  deallocate(AY)
  deallocate(V)
  deallocate(W)
  deallocate(D)
  deallocate(AD)
end subroutine TFQMR ! tfqmr 

! *****************************************************************************
subroutine BICG(b,x,KSPiter,adjt)
  ! Stablized version of BiConjugate Gradient, set up for solving
  ! A x = b using routines in  mult_Aii.
  ! solves for the interior (edge) field
  !
  ! modified from my matlab version of BICGstab...
  ! so the naming might sound a little different from conventional ones
  ! also added the optional adjoint to solve adjoint system A^Tx = b
  ! 
  ! NOTE: BICG actually performs two sub (or half) line searches within a 
  !       iteration, but here we only store the relerr for the second sub
  !       just to be compatitive with QMR
  !
  ! 
  ! interface...........
  ! redefining some of the interfaces for our convenience (locally)
  ! generic routines for vector operations for edge/ face nodes
  ! in a staggered grid
  !
  ! NOTE: this has not been extensively tested - I believe it feels a 
  ! little unstable (dispite the name)...
  ! if you have time reading this, test it!
  use modeloperator3d, only: A => mult_Aii, M1solve => PC_Lsolve,        &
     &                                      M2solve => PC_Usolve

  implicit none
  !  b is right hand side
  complex (kind=prec), intent(in), dimension(:)    :: b
  !  solution vector is x ... on input is provided with the initial
  !  guess, on output is the iterate with smallest residual. 
  !
  complex (kind=prec), intent(inout),dimension(:)  :: x
  type (solverControl_t), intent(inout)            :: KSPiter
  logical,intent(in),optional                      :: adjt

  ! local variables
  complex (kind=prec),allocatable,dimension(:)  :: R, RT, V,T
  complex (kind=prec),allocatable,dimension(:)  :: P,PT,PH,S,ST,SH,AX
  complex (kind=prec),allocatable,dimension(:)  :: xhalf,xmin
  real    (kind=prec)                           :: rnorm, bnorm, rnormin, btol
  complex (kind=prec)                           :: RHO, ALPHA, BETA, OMEGA
  complex (kind=prec)                           :: RTV,TT,RHO1
  integer                                       :: iter, xsize, imin
  integer                                       :: maxiter
  logical                                       :: adjoint, ilu_adjt, converged
 
  if (present(adjt)) then
      adjoint = adjt
      ilu_adjt = adjt
  else
      adjoint = .false.
      ilu_adjt = .false.
  endif
  xsize = size(x,1)
  ! Norm of rhs
  bnorm = SQRT(dot_product(b, b))
  if (isnan(bnorm)) then
  ! this usually means an inadequate model, in which case Maxwell's fails
      write(0,*) 'Error: b in BICG contains NaNs; exiting...'
      stop
  else if ( bnorm .eq. 0.0) then ! zero rhs -> zero solution
      write(0,*) 'Warning: b in BICG has all zeros, returning zero solution'
      x = b 
      KSPiter%niter=1
      KSPiter%failed=.false.
      KSPiter%rerr=0.0
      return
  endif
  ! allocate the local variables
  allocate(R(xsize))
  ! now calculate the (original) residual
  call A(x,adjoint,R)
  ! R= b - Ax, for inital guess x, that has been inputted to the routine
  R = b - R
  ! Norm of residual
  rnorm = CDSQRT(dot_product(R, R))
  btol = KSPiter%tol * bnorm
  if ( rnorm .le. btol ) then ! the first guess is already good enough
     ! returning
      KSPiter%niter=1
      KSPiter%failed=.false.
      KSPiter%rerr(1)=real(rnorm/bnorm)
      deallocate(R)
      return 
  else
      ! allocate the rest here
      ! allocate(xhalf(xsize))
      allocate(xmin(xsize))
      allocate(AX(xsize))
      allocate(RT(xsize))
      allocate(P(xsize))
      allocate(PT(xsize))
      allocate(PH(xsize))
      allocate(S(xsize))
      allocate(ST(xsize))
      allocate(SH(xsize))
      allocate(V(xsize))
      allocate(T(xsize))
  end if 
!================= Now start configuring the iteration ===================!
  ! the adjoint (shadow) residual
  rnormin = rnorm
  KSPiter%rerr(1) = real(rnormin/bnorm)
  ! write(6,*) 'initial residual',  KSPiter%rerr(1)
  converged = .false.
  maxiter = KSPiter%maxit 
  imin = 0
  RHO = C_ONE
  OMEGA = C_ONE
  RT = R
  xmin = x
  imin = 1
!============================== looooops! ================================!
  do iter= 1, maxiter
      RHO1 = RHO
      RHO = dot_product(RT,R)
      if (RHO .eq. 0.0) then
          KSPiter%failed = .true.
          exit
      end if 
      if (iter .eq. 1) then
          P = R
      else 
          BETA = (RHO/RHO1)*(ALPHA/OMEGA) 
          if (BETA .eq. 0.0) then
              KSPiter%failed = .true.
              exit
          end if
          P= R + BETA * (P - OMEGA * V);
      end if 
      ! first half of the iteration
      ! L
      call M1solve(P,ilu_adjt,PT)
      ! U
      call M2solve(PT,ilu_adjt,PH)
!      PH = P
      call A(PH,adjoint,V)
      RTV = dot_product(RT,V)
      if (RTV.eq.0.0) then
          KSPiter%failed = .true.
          exit
      end if
      ALPHA = RHO / RTV
      if (ALPHA.eq.0.0) then
          KSPiter%failed = .true.
          exit
      end if
      x = x + ALPHA * PH ! the first half 
      S = R - ALPHA*V  !residual for the 0.5 x
      ! second half of the iteration
      ! L
      call M1solve(S,ilu_adjt,ST)
      ! U
      call M2solve(ST,ilu_adjt,SH)
!      SH = S
      call A(SH,adjoint,T)
      TT = dot_product(T,T)
      if (TT.eq.0.0) then
          KSPiter%failed = .true.
          exit
      end if
      OMEGA = dot_product(T,S)/TT
      if (OMEGA.eq.0.0) then
          KSPiter%failed = .true.
          exit
      end if
      x = x + OMEGA * SH  ! the second half 
      R = S - OMEGA * T  !residual for the 1.0 x
      rnorm = CDSQRT(dot_product(R,R))
      KSPiter%rerr(iter) = real(rnorm / bnorm)
      ! write(6,*) 'iter # ',iter,' x residual: ', KSPiter%rerr(2*iter)
      if (rnorm.lt.btol) then
          KSPiter%failed = .false.
          KSPiter%niter = iter
          converged = .true.
          exit
      end if
      if (rnorm .lt. rnormin) then
          rnormin = rnorm
          xmin = x
          imin = iter
      end if
  end do
 
  if (.not. converged) then 
      ! it should be noted that this is the way my matlab version works
      ! the bicg will return the 'best' (smallest residual) iteration
      ! x = xmin;  ! comment this line 
      KSPiter%niter=maxiter
      ! KSPiter%rerr(KSPiter%maxit) = KSPiter%rerr(imin)  ! and this line
      ! to use the last iteration result instead of the 'best'
  end if

  ! deallocate(xhalf)
  deallocate(xmin)
  deallocate(AX)
  deallocate(R)
  deallocate(RT)
  deallocate(P)
  deallocate(PT)
  deallocate(PH)
  deallocate(S)
  deallocate(ST)
  deallocate(SH)
  deallocate(V)
  deallocate(T)
end subroutine BICG ! BICG
  



end module solver ! spsolver
