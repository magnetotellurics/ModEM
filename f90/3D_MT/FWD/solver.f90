! *****************************************************************************
! module containing iterative equation solvers. Uses operators and
! pre-conditioners defined in SG3DFWC1 to solve equations for divergence
! correction, induction operator. Source code is completely general; only the
! module interface is  specific to implementation of operators in SG3DFWC1
module solver

   use math_constants	! math/ physics constants
   use utilities, only: isnan
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
     real (kind=prec)                              :: tol
     ! actual number of iterations before return
     integer                                               :: niter
     ! relative error for each iteration
     real (kind=prec), pointer, dimension(:)   :: rerr
     ! logical variable indicating if algorithm "failed"
     logical                                               :: failed = .false.
  end type solverControl_t

Contains


! *****************************************************************************
! Solver contains subroutines for: a) PCG- a quasi-generic pre-conditioned
! congugate gradient, and b) QMR - Quasi-Minimal Residual method
! (pre-conditioned, no look-ahead)
! *****************************************************************************

subroutine PCG(b,x, PCGiter)
  ! Purpose: a quasi-generic pre-conditioned conjugate gradient
  ! routine, set up for solving DivCgrad phiSOL = phiEMrhs using
  ! routines in divCor.  Actual code is generic, but the interface
  ! is not

   ! inteface..............
   ! redefining some of the interfaces for our convenience (locally)
   ! generic routines for scalar operations on corner/ center nodes
   ! on a staggered grid
    use sg_scalar

  ! routines for divergence correction
  use modeloperator3D, only:  A => DivCgrad, Minv => DivCgradILU

  implicit none
  type (cscalar), intent(in)	        :: b
  type (cscalar), intent(inout)	        :: x
  type (solverControl_t), intent(inout) 	:: PCGiter

  ! local variables
  type (cscalar)	:: r,s,p,q
  complex(kind=prec)	:: beta,alpha,delta,deltaOld
  complex(kind=prec)	:: bnorm, rnorm
  integer		:: i

  if (.not.b%allocated) then
      write(0,*) 'b in PCG not allocated yet'
      stop
  end if

  if (.not.x%allocated) then
      write(0,*) 'x in PCG not allocated yet'
      stop
  end if

  ! Allocation of r, z, p, q
  Call create(x%grid, r, x%gridType)
  Call create(x%grid, s, x%gridType)
  Call create(x%grid, p, x%gridType)
  Call create(x%grid, q, x%gridType)

  Call A(x,r)
  Call linComb(C_ONE,b,C_MinusOne,r,r)
  bnorm = dotProd(b,b)
  rnorm = dotProd(r,r)
  i = 1
  PCGiter%rerr(i) = real(rnorm/bnorm)

  loop: do while ((PCGiter%rerr(i).gt.PCGiter%tol).and.(i.lt.PCGiter%maxIt))
     Call Minv(r,s)
     delta = dotProd(r,s)
     if(i.eq.1) then
        p = s
     else
        beta = delta/deltaOld
        Call linComb(C_ONE,s,beta,p,p)
     end if

     Call A(p,q)
     alpha = delta/dotProd(p,q)
     Call scMultAdd(alpha,p,x)
     Call scMultAdd(-alpha,q,r)
     deltaOld = delta
     i = i + 1
     rnorm = dotProd(r,r)
     PCGiter%rerr(i) = real(rnorm/bnorm)

  end do loop

  PCGiter%niter = i

  ! deallocate all the work arrays
  Call deall(r)
  Call deall(s)
  Call deall(p)
  Call deall(q)

end subroutine PCG ! PCG


! *****************************************************************************
subroutine QMR(b,x, QMRiter)
  ! Purpose ... a quasi-minimal residual method routine, set up for solving
  ! A x = b using routines in  multA. Actual code is generic, but interface
  ! is fairly specific

  ! interface...........
  ! redefining some of the interfaces for our convenience (locally)
  ! generic routines for vector operations for edge/ face nodes
  ! in a staggered grid
  use sg_vector

  ! routines for solving Maxwell's equation
  use modeloperator3D, only: A => multA_N, M1solve, M2solve

  implicit none
  !  b is right hand side
  type (cvector), intent(in)      	:: b
  !  solution vector is x ... on input is provided with the initial
  !   guess, on output is the most recent iterate
  type (cvector), intent(inout)   	:: x
  type (solverControl_t), intent(inout)	:: QMRiter

   ! local variables
  type (cvector)      	    :: AX,R,VT
  type (cvector)	    :: Y,Z,WT,V,W,YT,ZT,P,Q,PT,D,S
  logical                   :: adjoint, ilu_adjt
  complex (kind=prec)          :: ETA,PDE,EPSIL,RDE,BETA,DELTA,RHO
  complex (kind=prec)          :: PSI,RHO1,GAMM,GAMM1,THET,THET1,TM2
  complex (kind=prec)          :: bnorm,rnorm
  complex (kind=prec)          :: rhoInv,psiInv
  integer                   :: iter

  if (.not.b%allocated) then
      write(0,*) 'Error: b in QMR not allocated yet'
      stop
  end if

  if (.not.x%allocated) then
      write(0,*) 'Error: x in QMR not allocated yet'
      stop
  end if

  ! Allocate work arrays
  Call create(x%grid, AX, x%gridType)
  Call create(x%grid, R, x%gridType)
  Call create(x%grid, VT, x%gridType)
  Call create(x%grid, Y,x%gridType)
  Call create(x%grid, Z,x%gridType)
  Call create(x%grid, WT,x%gridType)
  Call create(x%grid, V,x%gridType)
  Call create(x%grid, W,x%gridType)
  Call create(x%grid, YT,x%gridType)
  Call create(x%grid, ZT,x%gridType)
  Call create(x%grid, P,x%gridType)
  Call create(x%grid, Q,x%gridType)
  Call create(x%grid, PT,x%gridType)
  Call create(x%grid, D,x%gridType)
  Call create(x%grid, S,x%gridType)

  ! NOTE: this iterative solver is QMR without look-ahead
  ! patterned after the scheme given on page 24 of Barrett et al.
  ! "Templates for the solution of linear systems of equations:
  ! Building blocks for iterative methods"
  ! Note that there are a couple of small differences, due to
  ! the fact that our system is complex (agrees with
  ! matlab6 version of qmr)

  adjoint = .false.
  ! R is Ax
  Call A(x, adjoint, R)
  ! b - Ax, for inital guess x, that has been inputted to the routine
  Call linComb(C_ONE,b,C_MinusOne,R,R)

  ! Norm of rhs, residual
  bnorm = CDSQRT(dotProd(b, b))
  rnorm = CDSQRT(dotProd(R, R))
!  print *, "bnorm,rnorm=",bnorm,rnorm
  ! this usually means an inadequate model, in which case Maxwell's fails
  if (isnan(abs(bnorm))) then
      write(0,*) 'Error: b in QMR contains NaNs; exiting...'
      stop
  endif

  !  iter is iteration counter
  iter = 1
  QMRiter%rerr(iter) = real(rnorm/bnorm)

  VT = R
  ilu_adjt = .false.
  Call M1solve(VT,ilu_adjt,Y)
  RHO = CDSQRT(dotProd(Y,Y))

  WT = R
  ilu_adjt = .true.
  Call M2solve(WT,ilu_adjt,Z)
  PSI  = CDSQRT(dotProd(Z,Z))
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
        exit
      endif

      rhoInv = (1/RHO)*cmplx(1.0, 0.0, 8)
      psiInv = (1/PSI)*cmplx(1.0, 0.0, 8)
      Call scMult(rhoInv, VT, V)
      Call scMult(psiInv, WT, W)
      Call scMult(rhoInv, Y, Y)
      Call scMult(psiInv, Z, Z)

      DELTA = dotProd(Z,Y)
      if(DELTA.eq.C_ZERO) then
        QMRiter%failed = .true.
        write(0,*) 'QMR FAILS TO CONVERGE : DELTA'
        exit
      endif

      ilu_adjt = .false.
      Call M2solve(Y,ilu_adjt,YT)
      ilu_adjt = .true.
      Call M1solve(Z,ilu_adjt,ZT)

      if (iter.eq.1) then
        P = YT
        Q = ZT
      else
      	! these calculations are only done when iter greater than 1
        PDE = -PSI*DELTA/EPSIL
        RDE = -RHO*CONJG(DELTA/EPSIL)
        Call linComb(C_ONE,YT,PDE,P,P)
        Call linComb(C_ONE,ZT,RDE,Q,Q)
      endif

      adjoint = .false.
      Call A(P, adjoint, PT)
      EPSIL = dotProd(Q,PT)
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
      Call linComb(C_ONE,PT,-BETA,V,VT)

      RHO1 = RHO
      ilu_adjt = .false.
      Call M1solve(VT, ilu_adjt, Y)
      RHO = CDSQRT(dotProd(Y,Y))

      adjoint = .true.
      Call A(Q, adjoint, WT)
      Call scMultAdd(-conjg(BETA),W,WT)

      ilu_adjt = .true.
      Call M2solve(WT,ilu_adjt,Z)
      PSI = CDSQRT(dotProd(Z,Z))

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
        Call scMult(ETA, P, D)
        Call scMult(ETA, PT, S)
      else
        TM2 = THET1*THET1*GAMM*GAMM
        Call linComb(ETA,P,TM2,D,D)
        Call linComb(ETA,PT,TM2,S,S)
      endif

      Call scMultAdd(C_ONE,D,x)
      Call scMultAdd(C_MinusONE,S,R)
      ! A new AX
      rnorm = CDSQRT(dotProd(R,R))
!      print *,"iter=",iter,' rnorm=',rnorm
      iter = iter + 1

      ! Keeping track of errors
      ! QMR book-keeping between divergence correction calls
      QMRiter%rerr(iter) = real(rnorm/bnorm)

  end do loop

  QMRiter%niter = iter

  ! deallocate all the work arrays
  Call deall(AX)
  Call deall(R)
  Call deall(VT)
  Call deall(Y)
  Call deall(Z)
  Call deall(WT)
  Call deall(V)
  Call deall(W)
  Call deall(YT)
  Call deall(ZT)
  Call deall(P)
  Call deall(Q)
  Call deall(PT)
  Call deall(D)
  Call deall(S)

end subroutine qmr ! qmr


end module solver ! solver
