!**********************************************************************
! driver modules for solving the forward EM problem, including setup and
! solver

module EMsolve3D
  use sg_boundary			! work between different data types
  					! (between boundary conditions and
					! complex vectors)
  use sg_diff_oper			 ! generic differential operators
  use sg_sparse_vector, only: add_scvector
  use modelOperator3d                   ! quasi-static Maxwell operator module
  use solnspace

  implicit none
  public	:: deallSolverDiag, deallEMsolveControl
  public        :: createSolverDiag, getEMsolveDiag, setEMsolveControl
  private	:: SdivCorr

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
  
  type :: emsolve_control
    ! Values of solver control parameters, e.g., read in from file
    !   plus other information on how the solver is to be initialized, called, etc.
    !  idea is that this is the public access version of this info, which is
    !   copied into private version for actual solver control
    integer                   ::      IterPerDivCor, MaxDivCor, MaxIterDivCor
    real(kind = 8)            ::       tolEMfwd, tolEMadj, tolDivCor
    logical                   ::      E0fromFile
    logical                   ::      UseDefaults
    logical                   ::      read_E0_from_File=.false.
    character (len=80)        ::      E0fileName
    integer                   ::      ioE0
  end type emsolve_control

  type :: emsolve_diag
    ! Solver diagnostic arrays, computed during run of forward solver.
    !  idea is that this is the public access version of this info, which is
    !   copied from the private version in module em_solve where this info is
    !   initially stored
    logical           ::      diagOut
    character (len=80)        :: fn_diagn
    integer                   :: ioDiag
    integer           ::              nIterTotal, nDivCor
    real(kind = 8), pointer, dimension(:)      ::      EMrelErr
    real(kind = 8), pointer, dimension(:,:)    ::      divJ
    real(kind = 8), pointer, dimension(:,:)    ::      DivCorRelErr
  end type emsolve_diag


  ! Default solver control parameters
  ! number of QMR iterations for each call to divergence correction:
  integer, parameter    ::              IterPerDivCorDef = 500
  ! maximum number of divergence correction calls allowed
  integer, parameter    ::              MaxDivCorDef = 200
  ! maximum number of PCG iterations for divergence correction
  integer, parameter    ::              MaxIterDivCorDef = 10
  ! misfit tolerance for convergence of EMsolve algorithm
  real(kind=prec), parameter       ::      tolEMDef = 1E-7
  ! misfit tolerance for convergence of divergence correction solver
  real(kind=prec), parameter       ::      tolDivCorDef = 1E-5

  save

  type(timer_t), private :: timer

  ! Actual values of control parameters must be set before first use,
  !     by call to setEMsolveControl
  !  of em_solve; are saved between calls, private to this module
  integer,  private        ::      IterPerDivCor, MaxDivCor, MaxIterDivCor
  integer,  private        ::      MaxIterTotal ! = MaxDivCor*IterPerDivCor
  real(kind=prec), private   ::      tolEMfwd, tolEMadj, tolDivCor

  ! EMsolve diagnostics: these are computed during execution of em_solve
  !   can be retrieved by call to getEmsolveDiag
  integer, private        ::      nIterTotal, nDivCor
  logical, private 		::	failed
  ! nIterTotal keeps tally on number of iterations so far
  ! nDivCor keeps tally on number of divergence correction so far
  real(kind=prec), pointer, dimension(:), private	::	EMrelErr
  real(kind=prec), pointer, dimension(:,:), private	::	divJ
  real(kind=prec), pointer, dimension(:,:), private	::	DivCorRelErr

Contains

  !**********************************************************************
  subroutine Mackie_solve(bRHS,omega,eSol,imode)
    ! redefine some of the interfaces (locally) for our convenience
    use sg_vector !, only: copy => copy_cvector, &
         !scMult => scMultReal_cvector
    ! generic routines for vector operations on the edge/face nodes
    ! in a staggered grid
    ! in copy, remember the order is copy(new, old) i.e new = old

    implicit none
    !  INPUTS:
	integer imode
    type (RHS_t), intent(in)		:: bRHS
    real(kind=prec), intent(in)	:: omega
    !  OUTPUTS:
    !     eSol must be allocated before calling this routine
    type (cvector), intent(inout)	:: eSol

    ! LOCAL VARIABLES
    logical				:: converged,trans,ltemp
    integer				:: status, iter
    complex(kind=prec)         	:: iOmegaMuInv
    type (cvector)			:: b,temp
    type (cscalar)			:: phi0
    type (cboundary)             	:: tempBC
    type (solverControl_t)			:: PBICGiter


    !  Zero solver diagnostic variables
    nIterTotal = 0
    nDivCor = 0
    EMrelErr = R_ZERO
    divJ = R_ZERO
    DivCorRelErr = R_ZERO
    failed = .false.

    trans = (bRHS%adj .eq. TRN)

    if (.not.eSol%allocated) then
       write(0,*) 'eSol in EMsolve not allocated yet'
       stop
    endif

    ! allocate/initialize local data structures
    Call create_cvector(bRHS%grid, b, eSol%gridType)
    Call create_cvector(bRHS%grid, temp, eSol%gridType)
    Call create_cboundary(bRHS%grid, tempBC)


    if(bRHS%nonzero_Source) then
       call create_cscalar(bRHS%grid,phi0,CORNER)
    endif

    
    if(trans) then
       !  In this case boundary conditions do not enter into forcing
       !    for reduced (interior node) linear system; solution on
       !    boundary is determined after solving for interior nodes
       if(bRHS%nonZero_Source) then    ! comb%b(k)%nonzero_source = .true.
          if(bRHS%sparse_Source) then   !comb%b(k)%sparse_Source = .false.
             ! Note: C_ONE = (1,0) (double complex)
	         call add_scvector(C_ONE,bRHS%sSparse,b)
          else
              b = bRHS%s  ! run this code
          endif
       else
          call zero(Esol)
          write(0,*) 'Warning: no sources for adjoint problem'
          write(0,*) 'Solution is identically zero'
          ! just copy input BC into boundary nodes of solution and return
          if(bRHS%nonzero_BC) then
             Call setBC(bRHS%bc, eSol)
          else
             Call setBC(tempBC, eSol)
          endif
          return
       endif

       ! use this processing
       call diagDiv(b,V_E,temp)
       call Div(temp,phi0) !debug  
            
    else
       ! In the usual forward model case BC do enter into forcing
       !   First compute contribution of BC term to RHS of reduced interior
       !    node system of equations : - A_IB*b
       if (bRHS%nonzero_BC) then
          !   copy from rHS structure into zeroed complex edge vector temp
          Call setBC(bRHS%bc, temp)

          !   Then multiply by curl_curl operator (use MultA_N ...
          !     Note that MultA_N already multiplies by volume weights
	  !     required to symetrize problem, so the result is V*A_IB*b)
          ltemp = .false.
          Call calcb(temp, ltemp, b)   ! debug ”“∂ÀœÓ
 
       endif

    endif


    if(bRHS%nonzero_Source) then
       iOmegaMuInv = ISIGN/cmplx(0.0,omega*MU_0,prec)
       call scMult(iOmegaMuInv,phi0,phi0)             
    endif
       
    ! Need to make sure first guess is zero on boundaries
    ! tempBC has all zeros on the boundaries
!    Call setBC(tempBC, eSol)


    ! resetting
    nIterTotal = 0
    nDivCor = 0
    call reset_time(timer)

    ! Initialize iteration control/diagnostic structure for QMR, PCG
    if (trans) then
       PBICGiter%tol = tolEMadj
    else
      if (bRHS%nonzero_BC) then
        PBICGiter%tol = tolEMfwd
      else
        PBICGiter%tol = tolEMadj
       end if
    end if


    PBICGiter%niter = 0
    PBICGiter%maxIt = IterPerDivCor
    allocate(PBICGiter%rerr(IterPerDivCor), STAT=status)
    PBICGiter%rerr = 0.0

    converged = .false.
    failed = .false.

    !if(bRHS%nonzero_Source) then
    !   nDivCor = 1
    !   Call SdivCorr(eSol,phi0)
    !endif
    
    loop: do while ((.not.converged).and.(.not.failed))  

	  Call PBICGSTAB(b,eSol,imode,PBICGiter,failed,converged,bRHS%nonzero_Source,phi0)

      failed = failed .or. PBICGiter%failed

      !  update diagnostics output from PBICGSTAB
      nIterTotal = nIterTotal + 1
      EMrelErr(nIterTotal) = PBICGiter%rerr(PBICGiter%niter)

      if(.false.) then
        nDivCor = nDivCor+1
        if( nDivCor < MaxDivCor) then
          ! do divergence correction
          if(bRHS%nonzero_Source) then
             Call SdivCorr(eSol,phi0)
          else
             Call SdivCorr(eSol)
          endif
        else
          ! max number of divergence corrections exceeded; convergence failed
          failed = .true.
        endif
	    endif

    end do loop


    if (output_level > 1) then
       write (*,'(a12,a20,i8,g15.7)') node_info, 'finished solving:', nIterTotal, EMrelErr(nIterTotal)
	   write (*,'(a12,a22,f12.6)')    node_info, ' time taken (mins) ', elapsed_time(timer)/60.0
    end if

!    if(bRHS%nonzero_Source) then
!       !Call setBC(tempBC, eSol)
!    else
!       Call setBC(bRHS%bc, eSol)
!    endif
    

    ! deallocate local temporary arrays
    Call deall(phi0)
    Call deall(b)
    Call deall(temp)
    Call deall(tempBC)
    deallocate(PBICGiter%rerr, STAT=status)

  end subroutine Mackie_solve

!****************************
subroutine PBICGSTAB(b,x,imode,PBICGiter,failed,converged,nzSrc,phi0)

  implicit none
  integer imode
  type (cvector), intent(in)      	:: b
  type (cvector), intent(inout)   	:: x
  type (cscalar), intent(in)		:: phi0
  type (solverControl_t)			:: PBICGiter
  ! local variables
  type (cvector) :: R,RTLD,P,V,T,PHAT,SHAT 
  complex (kind=prec) :: rho,temp,alpha,beta,womega,tt,st,c3
  complex (kind=prec)          :: bnorm,rnorm0,rnorm
  integer                   :: iter
  logical                   :: adjoint,failed,converged,nzSrc
  integer ipotloop,ipotmax,icount
  real (kind=prec) :: erre,erre1,erre2,errchk,errlast,errend,dcabs1

  ! Allocate work arrays
  Call create(x%grid, R, x%gridType)
  Call create(x%grid, RTLD, x%gridType)
  Call create(x%grid, P, x%gridType)
  Call create(x%grid, V,x%gridType)
  Call create(x%grid, T,x%gridType)
  Call create(x%grid, PHAT,x%gridType)
  Call create(x%grid, SHAT,x%gridType)

!    call writeEsol(x,107) !debug
!    call writeEsol(b,108) !debug
!	pause !debug

    ipotloop=0
    ipotmax=2
50  ipotmax=ipotmax+2
    if(ipotmax.gt.2) ipotmax=2

    call copy_cvector(V,b) !V=b

    call MInv_multB(V,R) !R=M^-1*b

    bnorm = CDSQRT(dotProd(R, R))
	   
    adjoint = .false.  
    call multA_N(x, adjoint, P) !P=A*x0

    call linComb(C_ONE,V,C_MinusOne,P,V) !V=b-P

    call MInv_multB(V,R) !R=M^-1*(b-P)

    call copy_cvector(RTLD,R) !RTLD=R

    if(ipotloop.eq.0) then
      rnorm0=CDSQRT(dotProd(R, R))
    endif

    !if(imode.eq.1)then
    !  print*,  &
    !     'iteration count and convergence for hy polarization'
    !else
    !  print*,  &
    !     'iteration count and convergence for hx polarization'
    !endif
    !print*,' '

	rho=C_ONE
	alpha=C_ONE
	womega=C_ONE

    call zero(V) 
	call zero(P) 

	do icount=1,PBICGiter%maxIt
      
	  temp=zdotu_cvector_f(RTLD, R)

      if (1e9*dcabs1(temp) .le. eps()) then
        call errStop('convergence fails')
      end if

      beta=alpha*(temp/rho)/womega
      rho=temp

      call linComb(-womega,V,C_One,P,P) !P=P-w*V

	  call scMult(beta,P,P)

	  call linComb(C_One,R,C_One,P,P) !P=R+P

      adjoint = .false.  
      call multA_N(P, adjoint, PHAT) !PHAT=A*P

      call MInv_multB(PHAT,V) !V=M^-1*PHAT

	  temp=zdotu_cvector_f(RTLD, V)

	  alpha=rho/temp

	  call linComb(C_One,R,-alpha,V,R)

      adjoint = .false.  
      call multA_N(R, adjoint, SHAT) !SHAT=A*R

      call MInv_multB(SHAT,T) !T=M^-1*SHAT

      tt=zdotu_cvector_f(T, T)

	  st=zdotu_cvector_f(R, T)
      
      womega=st/tt

	  call scMultAdd(womega,R,x)

      call linComb(alpha,P,C_ONE,x,x)

	  call linComb(-womega,T,C_ONE,R,R)

      rnorm=CDSQRT(dotProd(R, R))

      erre=rnorm/rnorm0
      erre1=rnorm/bnorm
      PBICGiter%niter=PBICGiter%niter+1
	  PBICGiter%rerr(PBICGiter%niter)=erre1
      !write(6,'(i5,2(a,1pg11.3))')  &
      !      PBICGiter%niter,'  r/r0 =',erre,'  r/b =',erre1
      if(PBICGiter%niter.eq.PBICGiter%maxIt) then
        if(erre1.le.PBICGiter%tol) then
          converged=.true.
        else
          converged=.false.
		  PBICGiter%failed=.true.
		endif
        !write(6,*)' '
        goto 999
      end if
      
      if((icount.ge.ipotmax)) then

	    if(.true.) then
		  nDivCor = nDivCor+1	  
		  if( nDivCor < MaxDivCor) then
		    if(nzSrc) then
		      call SdivCorr(x,phi0)
		      !print*,'Pause in SensDivCorr' !debug
              !pause !debug
		    else
		      call SdivCorr(x)
		    endif
		  else
		    failed = .true.
		  endif 	  
		  
          !print*
          ipotloop=ipotloop+1
          go to 50
        endif

      end if

      if(PBICGiter%niter.lt.PBICGiter%maxIt) then
        if(erre1.le.PBICGiter%tol) then
          converged=.true.
          !write(6,*)' '
          goto 999
        endif
      end if
      errlast=erre1

	end do

999 continue

  Call deall(R)
  Call deall(RTLD)
  Call deall(P)
  Call deall(V)
  Call deall(T)
  Call deall(PHAT)
  Call deall(SHAT)
  
end subroutine PBICGSTAB


  !**************************************
  subroutine SdivCorr(inE, phi0)
    ! Purpose: driver routine to compute divergence correction for input electric
    ! field vector inE output corrected ! electric field in outE
    !  Optional argument phi0 is scaled divergence of source term
    !   to be subtracted from current divergence

    implicit none
    type (cvector), intent(inout)	:: inE
    type (cscalar), intent(in), optional	:: phi0

    !  local variables
	  integer job
    type (solverControl_t)			:: PCGiter
    type (rscalar)		        :: phiSol, RephiRHS, ImphiRHS
    complex (kind=prec)        	:: c2
    integer				:: status
    character (len=80)              	:: Desc = ''
    logical				:: SourceTerm

    SourceTerm = present(phi0)

    ! initialize PCGiter (maximum iterations allowed per set of diveregence
    ! correction, error tolerence, and relative error book keeping)
    PCGiter%maxIt = MaxIterDivCor
    PCGiter%tol = tolDivCor
    allocate(PCGiter%rerr(PCGiter%maxIt), STAT = status)
    PCGiter%rerr = 0.0

    Desc = CORNER
    ! alocating phiSol, phiRHS
    Call create_rscalar(inE%grid, phiSol, Desc)
    Call create_rscalar(inE%grid, RephiRHS, Desc)
	  Call create_rscalar(inE%grid, ImphiRHS, Desc)

    ! compute divergence of magnetic field
    Call DivC(inE, RephiRHS, ImphiRHS)
   

    !  If source term is present, subtract from divergence of H
    if(SourceTerm) then
       call subtractr_rcscalar(RephiRHS,phi0,RephiRHS)  ! problems
       call subtracti_rcscalar(ImphiRHS,phi0,ImphiRHS)
    endif


! ======================Real Part Correction==============================
    ! compute the size of current Divergence before (using dot product)
    divJ(1,nDivCor) = dotProd(RephiRHS,RephiRHS)+dotProd(ImphiRHS,ImphiRHS)

    ! point-wise multiplication with volume weights centered on corner nodes
    Call dMult_rscalar(V_N,RephiRHS,RephiRHS)


    ! PCG is a generic pre-conditioned CG algorithm
    Call PCG(RephiRHS,phiSol,PCGiter)
    DivCorRelErr(:,nDivCor) = PCGiter%rerr


    !if (output_level > 2) then
    !   write (*,'(a12,a41,i5,g15.7)') node_info, &
    !      'finished Real Part divergence correction:', &
	!	   PCGiter%niter, PCGiter%rerr(PCGiter%niter)
    !end if

    ! compute gradient of phiSol (Divergence correction for inE)
	  job=0
    Call Grad(phiSol,inE,job)

! ======================Imaginary Part Correction==============================
    ! point-wise multiplication with volume weights centered on corner nodes
    Call dMult_rscalar(V_N,ImphiRHS,ImphiRHS)

    ! PCG is a generic pre-conditioned CG algorithm
	  call zero(phiSol)
    Call PCG(ImphiRHS,phiSol,PCGiter)
    DivCorRelErr(:,nDivCor) = PCGiter%rerr

    !if (output_level > 2) then
    !   write (*,'(a12,a41,i5,g15.7)') node_info, &
    !      'finished Imag Part divergence correction:', &
	!	   PCGiter%niter, PCGiter%rerr(PCGiter%niter)
    !end if

    ! compute gradient of phiSol (Divergence correction for inE)
	  job=1
    Call Grad(phiSol,inE,job)

    ! divergence of the corrected output electrical field
    Call DivC(inE, RephiRHS, ImphiRHS)


    !  If source term is present, subtract from divergence of currents
    if(SourceTerm) then
       call subtractr_rcscalar(RephiRHS,phi0,RephiRHS)
       call subtracti_rcscalar(ImphiRHS,phi0,ImphiRHS)
    endif

    ! as in WS code, compute the size of current Divergence after
    ! (using the dot product)
    divJ(2,nDivCor) = dotProd(RephiRHS,RephiRHS)+dotProd(ImphiRHS,ImphiRHS)

    ! output level defined in basic file_units module
    !if (output_level > 3) then
    !   write(*,'(a12,a48,g15.7)') node_info, 'divergence of magnetic filed before correction: ', divJ(1, nDivCor)
    !   write(*,'(a12,a48,g15.7)') node_info, 'divergence of magnetic filed  after correction: ', divJ(2, nDivCor)
    !end if

    ! deallocate the temporary work arrays
    Call deall(phiSol)
    Call deall(RephiRHS)
	  Call deall(ImphiRHS)
    deallocate(PCGiter%rerr, STAT = status)

  end subroutine SdivCorr



!***************************
subroutine PCG(b,x, PCGiter)
  use sg_scalar


  implicit none
  type (rscalar), intent(in)	        :: b
  type (rscalar), intent(inout)	        :: x
  type (solverControl_t), intent(inout) 	:: PCGiter
  ! local variables
  integer i
  real (kind=prec)          :: bnorm,rnorm,b1,c1,d1,e1,l1
  type (rscalar)  :: R,Z,Q,P,S,D
  logical :: debug=.false.

  if (.not.b%allocated) then
      write(0,*) 'b in PCG not allocated yet'
      stop
  end if

  if (.not.x%allocated) then
      write(0,*) 'x in PCG not allocated yet'
      stop
  end if

  b1=R_ZERO
  c1=R_ZERO
  d1=R_ZERO
  e1=R_ZERO
  l1=R_ZERO

  call create_rscalar(x%grid, R, x%gridType)
  call create_rscalar(x%grid, Z, x%gridType)
  call create_rscalar(x%grid, Q, x%gridType)
  call create_rscalar(x%grid, P, x%gridType)
  call create_rscalar(x%grid, S, x%gridType)
  call create_rscalar(x%grid, D, x%gridType)

  call DivA(x,R) !R=A*x0

  call linComb(ONE,b,MinusONE,R,R) !R=b-A*x0

  bnorm = sqrt(dotProd(b,b))
  rnorm = sqrt(dotProd(r,r))
  i = 1
  PCGiter%rerr(i) = rnorm/bnorm

  loop: do while ((PCGiter%rerr(i).gt.PCGiter%tol).and.(i.lt.PCGiter%maxIt))
! do i=1,10
   
    call DivMInvB(R,Z) !M*z=R

	call DivA(Z,Q) !Q=A*Z

    c1=zdotu_rscalar_f(Q,S)


    if((dabs(e1).eq.R_ZERO).or.(i.eq.1)) then
      b1=R_ZERO
    else
      b1=-c1/e1
    end if

    call linComb(MinusONE,Z,b1,D,D)

	call linComb(ONE,Q,b1,P,P)

	call DivMInvB(P,S)

	e1=zdotu_rscalar_f(P,S)
    d1=zdotu_rscalar_f(Z,P)

    if(dabs(e1).eq.R_ZERO) then
      l1=R_ZERO
    else
      l1=-d1/e1
    end if

	call linComb(ONE,x,l1,D,x)
  	
	call linComb(ONE,R,l1,P,R)

     i = i + 1
     rnorm = sqrt(dotProd(r,r))
     PCGiter%rerr(i) = rnorm/bnorm

  end do loop

  PCGiter%niter = i

  Call deall(R)
  Call deall(Z)
  Call deall(Q)
  Call deall(P)
  Call deall(S)
  Call deall(D)
end subroutine PCG


  !**********************************************************************
  ! setEMsolveControl sets actual solver control parameters, using info
  !  in structure solverControl, and allocates diagnostic arrays
  subroutine setEMsolveControl(solverControl,tolEM)

     type (emsolve_control), intent(in)	::	solverControl
     real (8), intent(in), optional     ::  tolEM

     if(solverControl%UseDefaults) then
        IterPerDivCor = IterPerDivCorDef
        MaxDivCor = MaxDivCorDef
        MaxIterTotal = MaxDivCor*IterPerDivCor
        MaxIterDivCor = MaxIterDivCorDef
        tolEMfwd = tolEMDef
        tolEMadj = tolEMDef
        tolDivCor = tolDivCorDef
     else
        IterPerDivCor = solverControl%IterPerDivCor
        MaxDivCor = solverControl%MaxDivCor
        MaxIterTotal = MaxDivCor*IterPerDivCor
        MaxIterDivCor = solverControl%MaxIterDivCor
        tolEMfwd = solverControl%tolEMfwd
        tolEMadj = solverControl%tolEMadj
        tolDivCor = solverControl%tolDivCor
     endif

     if (present(tolEM)) then
        tolEMfwd = tolEM
        tolEMadj = tolEM
     endif

     !  first check to see if diagnostic arrays are allocated
     !     ... if so deallocate first
     if(associated(EMrelErr)) then
        deallocate(EMrelErr)
     endif
     if(associated(divJ)) then
        deallocate(divJ)
     endif
     if(associated(DivCorRelErr)) then
        deallocate(DivCorRelErr)
     endif
     !   then allocate all arrays
     allocate(EMrelErr(MaxIterTotal))
     allocate(divJ(2,MaxDivCor))
     allocate(DivCorRelErr(MaxIterDivCor,MaxDivCor))

  end subroutine setEMsolveControl

   ! ***************************************************************************
   ! * readEMsolveControl reads the EM solver configuration from file
   subroutine readEMsolveControl(solverControl,rFile,fileExists,tolEM)

	type(emsolve_control), intent(inout)	:: solverControl
    character(*), intent(in)		        :: rFile
	logical, intent(out), optional          :: fileExists
	real(8), intent(in), optional           :: tolEM
    integer									:: ios
	logical                             	:: exists
	character(80)							:: string
	integer									:: istat

    ! Initialize inverse solver configuration

    inquire(FILE=rFile,EXIST=exists)
    if (present(fileExists)) then
       fileExists = exists
    end if

    if (.not. exists) then
       solverControl%UseDefaults = .true.
       if (present(tolEM)) then
          call setEMsolveControl(solverControl,tolEM)
       else
          call setEMsolveControl(solverControl)
       end if
       return
    else
       solverControl%UseDefaults = .false.
       write(*,*) node_info,'Reading EM solver configuration from file ',trim(rFile)
    end if

    !open (unit=ioFwdCtrl,file=rFile,status='old',iostat=ios)
     open (unit=ioFwdCtrl,file=rFile,form='formatted',status='old',iostat=ios)
    if(ios/=0) then
       write(0,*) node_info,'Error opening file: ', rFile
    end if

    ! This is the list of options specified in the startup file

    read (ioFwdCtrl,'(a48,i5)') string,solverControl%IterPerDivCor
    if (output_level > 2) then
       write (*,*)
       write (*,'(a12,a48,i5)') node_info,string,solverControl%IterPerDivCor
    end if
    read (ioFwdCtrl,'(a48,i5)') string,solverControl%MaxDivCor
    if (output_level > 2) then
       write (*,'(a12,a48,i5)') node_info,string,solverControl%MaxDivCor
    end if
    read (ioFwdCtrl,'(a48,i5)') string,solverControl%MaxIterDivCor
    if (output_level > 2) then
       write (*,'(a12,a48,i5)') node_info,string,solverControl%MaxIterDivCor
    end if
    read (ioFwdCtrl,'(a48,g15.7)') string,solverControl%tolEMfwd
    if (output_level > 2) then
       write (*,'(a12,a48,g15.7)') node_info,string,solverControl%tolEMfwd
    end if
    read (ioFwdCtrl,'(a48,g15.7)') string,solverControl%tolEMadj
    if (output_level > 2) then
       write (*,'(a12,a48,g15.7)') node_info,string,solverControl%tolEMadj
    end if
    read (ioFwdCtrl,'(a48,g15.7)') string,solverControl%tolDivCor
    if (output_level > 2) then
       write (*,'(a12,a48,g15.7)') node_info,string,solverControl%tolDivCor
    end if



! Check if there is addtional line.
! if yes, it is corresponde to the larger E field solution.

      read(ioFwdCtrl,'(a48)',advance='no',iostat=istat) string
      read(ioFwdCtrl,'(a80)',iostat=istat) solverControl%E0fileName      
   if (istat .eq. 0 ) then
     if (index(string,'#')>0) then
      ! This is a comment line
       solverControl%read_E0_from_File=.false.
     else
	     if (output_level > 2) then
	       write (*,'(a12,a48,a80)') node_info,string,solverControl%E0fileName
	     end if
       solverControl%read_E0_from_File=.true.
       solverControl%ioE0=ioE
     end if

 else
     solverControl%read_E0_from_File=.false.
  end if

    close(ioFwdCtrl)

    call setEMsolveControl(solverControl)

   end subroutine readEMsolveControl

  !**********************************************************************
  !   deallEMsolveControl deallocate
  subroutine  deallEMsolveControl()

     integer istat

     deallocate(EMrelErr, STAT=istat)
     deallocate(divJ, STAT=istat)
     deallocate(DivCorRelErr, STAT=istat)

  end subroutine deallEMsolveControl

  !**********************************************************************
  ! getEMsolveDiag retrieves solver diagnositics
  subroutine getEMsolveDiag(solverDiagnostics)

     type (emsolve_diag), intent(inout)    ::      solverDiagnostics

     solverDiagnostics%nIterTotal = nIterTotal
     solverDiagnostics%nDivCor = nDivCor
     solverDiagnostics%EMrelErr = EMrelErr
     solverDiagnostics%divJ = divJ
     solverDiagnostics%DivCorRelErr = DivCorRelErr

  end subroutine getEMsolveDiag

  !***************************************************************************
  ! * createSolverDiag initializes emsolve_diag structure
  subroutine createSolverDiag(solverParams,solverDiag)

    implicit none
    type (emsolve_control), intent(in)  :: solverParams
    type (emsolve_diag), intent(inout)  :: solverDiag
    integer                             :: maxIterTotal

    maxIterTotal = solverParams%MaxDivCor*solverParams%IterPerDivCor
    allocate(solverDiag%EMrelErr(maxIterTotal))
    allocate(solverDiag%divJ(2,solverParams%MaxDivCor))
    allocate(solverDiag%DivCorRelErr(solverParams%MaxIterDivCor,  &
        solverParams%MaxDivCor))

  end subroutine createSolverDiag

  !***************************************************************************
  ! * deallSolverDiag deallocates emsolve_diag structure
  subroutine deallSolverDiag(solverDiag)
    type (emsolve_diag), intent(inout)  :: solverDiag

    deallocate(solverDiag%EMrelErr)
    deallocate(solverDiag%divJ)
    deallocate(solverDiag%DivCorRelErr)

  end subroutine deallSolverDiag


!=======================================================================

      real(kind=prec) function eps()

        real(kind=prec) epst,epsp1
      
        epst=1.d0
5       eps=epst
        epst=eps/2.d0
        epsp1=epst+1.d0
        if(epsp1 .gt. 1.d0) goto 5

        return
      end function eps

!=======================================================================

end module EMsolve3D
