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
  use solver				! generic solvers
  use solnspace

  implicit none
  public	:: FWDSolve3D ,ModelOperatorSetup, ModelOperatorCleanup
  public	:: deallSolverDiag, deallEMsolveControl
  public        :: createSolverDiag, getEMsolveDiag, setEMsolveControl
  private	:: SdivCorr

  type :: emsolve_control
    ! Values of solver control parameters, e.g., read in from file
    !   plus other information on how the solver is to be initialized, called, etc.
    !  idea is that this is the public access version of this info, which is
    !   copied into private version for actual solver control
    integer                   ::      IterPerDivCor, MaxDivCor, MaxIterDivCor
    real(kind = 8)            ::      tolEMfwd, tolEMadj, tolDivCor
    logical                   ::      E0fromFile
    logical                   ::      UseDefaults
    logical                   ::      read_E0_from_File=.false.
    character (len=80)        ::      E0fileName
    integer                   ::      ioE0
    character (len=80)        ::      AirLayersMethod
    integer                   ::      AirLayersNz
    real(kind = 8)            ::      AirLayersMaxHeight, AirLayersAlpha, AirLayersMinTopDz
    real(kind = 8), pointer, dimension(:)   :: AirLayersDz
    logical                   ::      AirLayersPresent=.false.
    character (len=10)        ::      solver_name="QMR"

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
  integer, parameter    ::              IterPerDivCorDef = 40
  ! maximum number of divergence correction calls allowed
  integer, parameter    ::              MaxDivCorDef = 20
  ! maximum number of PCG iterations for divergence correction
  integer, parameter    ::              MaxIterDivCorDef = 100
  ! misfit tolerance for convergence of EMsolve algorithm
  real(kind=prec), parameter       ::   tolEMDef = 1E-7
  ! misfit tolerance for convergence of divergence correction solver
  real(kind=prec), parameter       ::   tolDivCorDef = 1E-5
  !Solver name, by default we use QMR
  character (len=10)               ::   solver_name="QMR"

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
! Solves the forward EM problem;
!
! If bRHS%adj = 'TRN' solves transposed problem  A^T x = b
!
! Note [AK 2018-05-10]: any physical source has already been pre-multiplied
! by [- ISIGN i\omega\mu_0] to yield [- ISIGN i\omega\mu_0 j] on input to this
! routine. Note that this also holds for the secondary field formulation,
! where j = dsigma * e, as well as for the tidal forcing, where j = sigma (v x B).
! However, we still want to pre-compute the source RHS outside of this routine, for
! generality: specifically, Jmult supplies an interior source on the RHS that is
! not physical and is not pre-multiplied by that factor (except in Pmult). So it's
! cleaner to pass on the complete interior forcing in bRHS.
! For divergence correction, we divide by [+ ISIGN i\omega\mu_0] to get [- Div(j)].
! The plus sign is needed because we're taking the divergence of
!  curl(curl(E)) + ISIGN i\omega\mu_0 sigma E = f - curl(curl(b))
! Terms 1 and 4 cancel, leaving Div(sigma E) - Div(f)/(+ ISIGN i\omega\mu_0) = 0.
! For a physical source j, this is equivalent to Div(sigma E) + Div(j) = 0; but
! the divergence correction may be applied also for non-physical sources, such as
! in Jmult ('FWD') and JmultT ('TRN').

  subroutine FWDsolve3D(bRHS,omega,eSol)

    ! redefine some of the interfaces (locally) for our convenience
    use sg_vector !, only: copy => copy_cvector, &
         !scMult => scMultReal_cvector
    ! generic routines for vector operations on the edge/face nodes
    ! in a staggered grid
    ! in copy, remember the order is copy(new, old) i.e new = old

    implicit none
    !  INPUTS:
    type (RHS_t), intent(in)		:: bRHS
    real(kind=prec), intent(in)	:: omega
    !  OUTPUTS:
    !     eSol must be allocated before calling this routine
    type (cvector), intent(inout)	:: eSol

    ! LOCAL VARIABLES
    logical				:: converged,trans,ltemp
    integer				:: status, iter
    complex(kind=prec)			:: i_omega_mu
    complex(kind=prec)         	:: iOmegaMuInv
    type (cvector)			:: b,temp
    type (cscalar)			:: phi0
    type (cboundary)             	:: tempBC
    type (solverControl_t)			:: QMRiter


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
    ! this is just a work array, at a given instance only single frequency and
    ! mode is being used
    Call create_cboundary(bRHS%grid, tempBC)

    if(bRHS%nonzero_Source) then
       call create_cscalar(bRHS%grid,phi0,CORNER)
    endif

    i_omega_mu = cmplx(0.,1.0d0*ISIGN*MU_0*omega,kind=prec)


    ! Using boundary condition and sources from rHS data structure
    ! construct vector b (defined only on interior nodes) for rHS of
    ! reduced (interior nodes only) linear system of equations

    if(trans) then
       !  In this case boundary conditions do not enter into forcing
       !    for reduced (interior node) linear system; solution on
       !    boundary is determined after solving for interior nodes
       if(bRHS%nonZero_Source) then
          if(bRHS%sparse_Source) then
             ! Note: C_ONE = (1,0) (double complex)
	     call add_scvector(C_ONE,bRHS%sSparse,b)
          else
              b = bRHS%s
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

       !  divide by volume weights before computing divergence of sources
       call diagDiv(b,V_E,temp)
       call Div(temp,phi0)
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
          Call MultA_N(temp, ltemp, b)

          !  change sign of result; obtained [- V_E A_IB b]
          Call scMult(MinusOne,b,b)
       endif

       ! Add internal sources if appropriate: Note that these must be multiplied
       !  explicitly by volume weights; in the interior of the domain,
       ! [V_E^{-1} C_II^T V_F] C_II e + i\omega\mu_0 \sigma e = - i\omega\mu_0 j - [V_E^{-1} C_II^T V_F] C_IB b
       ! Multiply by V_E throughout to obtain a symmetric system:
       !    V_E A_II e = - i\omega\mu_0 V_E j - V_E A_IB b
       ! where A_II = V_E^{-1} C_II^T V_F C_II + i\omega\mu_0 \sigma, and A_IB = V_E^{-1} C_II^T V_F C_IB.
       if (bRHS%nonzero_Source) then
          if (bRHS%sparse_Source) then
             ! temp  = bRHS%sSparse
             call zero(temp)
             call add_scvector(C_ONE,bRHS%sSparse,temp)
             call Div(temp,phi0)
          else
	     temp = bRHS%s
          endif

	  ! At this point, temp = - ISIGN * i\omega\mu_0 j
	  ! Now Div(f) - will later divide by i_omega_mu to get the general divergence correction
	  call Div(temp,phi0)

          call diagMult(V_E,temp,temp)
	  ! Now temp stores [-i\omega\mu_0 V_E j], and b stores [-V_E A_IB b]

          if(bRHS%nonzero_BC) then
             Call add(temp,b,b)
          else
             Call copy_cvector(b,temp)
          endif
       endif
    endif

    if(bRHS%nonzero_Source) then
       call scMult(C_ONE/i_omega_mu,phi0,phi0)
    endif

    ! Need to make sure first guess is zero on boundaries
    ! tempBC has all zeros on the boundaries
    Call setBC(tempBC, eSol)

    ! Outer part of QMR loop ... alternates between Calls to QMR
    ! and Calls to divcor  ... this will be part of EMsolve
    !
    ! eSol = current best solution
    ! b = rHS
    !
    ! at present we don't really have the option to skip
    ! the divergence correction.  Not sure how/if this should
    ! be done.

    ! resetting
    nIterTotal = 0
    nDivCor = 0
    call reset_time(timer)

    ! Initialize iteration control/diagnostic structure for QMR, PCG
    if (trans) then
       QMRiter%tol = tolEMadj
    else
      if (bRHS%nonzero_BC) then
        QMRiter%tol = tolEMfwd
      else
        QMRiter%tol = tolEMadj
       end if
    end if


    QMRiter%niter = 0
    QMRiter%maxIt = IterPerDivCor
    allocate(QMRiter%rerr(IterPerDivCor), STAT=status)
    QMRiter%rerr = 0.0

    converged = .false.
    failed = .false.
    !  idea to test: for non-zero source START with divergence
    !    correction
    if(bRHS%nonzero_Source) then
       Call copy_cvector(temp, eSol)
       nDivCor = 1
       Call SdivCorr(temp,eSol,phi0)
    endif
    loop: do while ((.not.converged).and.(.not.failed))

    if (trim(solver_name) .eq. "QMR") then
       Call QMR(b, eSol,QMRiter)
    elseif (trim(solver_name) .eq. "BICG") then
       Call BICG(b, eSol,QMRiter)
    else
       Write(*,*)"Unknown Solver Method"
    end if


       ! algorithm is converged when the relative error is less than tolerance
       ! (in which case QMRiter%niter will be less than QMRiter%maxIt)
       converged = QMRiter%niter .lt. QMRiter%maxIt

       ! there are two ways of failing: 1) QMR did not work or
       !        2) total number of divergence corrections exceeded
       failed = failed .or. QMRiter%failed

       !  update diagnostics output from QMR
       do iter = 1,QMRiter%niter
           EMrelErr(nIterTotal+iter) = QMRiter%rerr(iter)
       enddo
       nIterTotal = nIterTotal + QMRiter%niter

       nDivCor = nDivCor+1
       if( nDivCor < MaxDivCor) then
          ! do divergence correction
          Call copy_cvector(temp, eSol)
          if(bRHS%nonzero_Source) then
             Call SdivCorr(temp,eSol,phi0)
          else
             Call SdivCorr(temp,eSol)
          endif
       else
          ! max number of divergence corrections exceeded; convergence failed
          failed = .true.
       endif

    end do loop

    if (output_level > 1) then
       write (*,'(a12,a20,i8,g15.7)') node_info, 'finished solving:', nIterTotal, EMrelErr(nIterTotal)
	   write (*,'(a12,a22,f12.6)')    node_info, ' time taken (mins) ', elapsed_time(timer)/60.0
    end if

    !  After solving symetrized system, need to do different things for
    !   transposed, standard cases
    if(trans) then
    !   Multiply solution on interior nodes by volume weights
       ! eSol = V_E*eSol
        call diagMult(V_E,eSol,eSol)
       ! then compute solution on boundary nodes: first  A_IB^T eSol
       call AdjtBC(eSol, tempBC)
       ! then b - A_IB^T eSol, where b is input boundary values (if any)
       !  C_MinusONE = (-1,0) and C_ONE = (1,0) (double complex) are defined in
       !    SG_Basics/math_constants.f90
       !   tempBC = rhs%bc - tempBC
       if(bRHS%nonzero_BC) then
           call linComb_cboundary(C_MinusONE,tempBC,C_ONE,bRHS%bc,tempBC)
       else
           call scMult_cboundary(C_MinusONE, tempBC, tempBC)
       endif
       !  and copy result into boundary nodes of eSol
       Call setBC(tempBC, eSol)
    else
        ! just copy input BC into boundary nodes of solution
       if(bRHS%nonzero_BC) then
          Call setBC(bRHS%bc, eSol)
       else
          Call setBC(tempBC, eSol)
       endif
    endif



    ! deallocate local temporary arrays
    Call deall(phi0)
    Call deall(b)
    Call deall(temp)
    Call deall(tempBC)
    deallocate(QMRiter%rerr, STAT=status)

  end subroutine FWDsolve3D

!**********************************************************************
! solver_divcorr contains the subroutine that solves the divergence correction
! using pre-conditioned conjugate gradient.
! Taking the divergence of curl(curl(E)) + i\omega\mu_0 sigma E = f, where f
! is a general RHS, and since div(curl(xxx)) is identically zero, we get
! i\omega\mu_0 div(sigma E) = div(f) => div(sigma E) - div(f)/(i\omega\mu_0) = 0.
! This is what we want for a general adjoint solution. However, note that for
! a physical current source j, we would have f = - i\omega\mu_0 j,
! and would there have div(sigma E) + div(J) = 0, where J is the physical source
! consistent with curl(H) = sigma E + J.
subroutine SdivCorr(inE,outE,phi0)
  ! Purpose: driver routine to compute divergence correction for input electric
  ! field vector inE output corrected ! electric field in outE
  !  Optional argument phi0 is scaled divergence of source term
  !   to be subtracted from current divergence

  use sg_scalar ! mult => diagMult_scalar, linComb => linComb_cscalar
  ! rename routines for linear algebra operations; change to apply
  !  PCG to different problem

  implicit none
  type (cvector), intent(in)	:: inE
  type (cvector), intent(inout)	:: outE
  type (cscalar), intent(in), optional	:: phi0

  !  local variables
  type (solverControl_t)			:: PCGiter
  type (cscalar)		        :: phiSol, phiRHS
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
  Call create_cscalar(inE%grid, phiSol, Desc)
  Call create_cscalar(inE%grid, phiRHS, Desc)

  ! compute divergence of currents for input electric field
  Call DivC(inE, phiRHS)

  !  If source term is present, subtract from divergence of currents
  if(SourceTerm) then
     call subtract(phiRHS,phi0,phiRHS)
  endif

  ! compute the size of current Divergence before (using dot product)
  divJ(1,nDivCor) = sqrt(dotProd(phiRHS,phiRHS))

  ! point-wise multiplication with volume weights centered on corner nodes
  Call diagMult(V_N,phiRHS,phiRHS)

  ! PCG is a generic pre-conditioned CG algorithm
  Call PCG(phiRHS,phiSol,PCGiter)
  DivCorRelErr(:,nDivCor) = PCGiter%rerr

  if (output_level > 2) then
     write (*,'(a12,a32,i5,g15.7)') node_info, &
        'finished divergence correction:', PCGiter%niter, PCGiter%rerr(PCGiter%niter)
  end if

  ! compute gradient of phiSol (Divergence correction for inE)
  Call Grad(phiSol,outE)

  ! subtract Divergence correction from inE
  Call subtract(inE, outE, outE)

  ! divergence of the corrected output electrical field
  Call DivC(outE,phiRHS)

  !  If source term is present, subtract from divergence of currents
  if(SourceTerm) then
     call subtract(phiRHS,phi0,phiRHS)
  endif

  ! as in WS code, compute the size of current Divergence after
  ! (using the dot product)
  divJ(2,nDivCor) = sqrt(dotProd(phiRHS,phiRHS))

  ! output level defined in basic file_units module
  if (output_level > 2) then
     write(*,'(a12,a47,g15.7)') node_info, 'divergence of currents before correction: ', divJ(1, nDivCor)
     write(*,'(a12,a47,g15.7)') node_info, 'divergence of currents  after correction: ', divJ(2, nDivCor)
  end if

  ! deallocate the temporary work arrays
  Call deall(phiSol)
  Call deall(phiRHS)
  deallocate(PCGiter%rerr, STAT = status)

end subroutine SdivCorr ! SdivCorr

  !**********************************************************************
  ! Bundles the Inits that are used for an EM problem. These Inits can be
  ! used separately as well.
  subroutine ModelOperatorSetUp()

    ! Initialize EM equation operator arrays
    Call AdiagInit()
    Call DiluInit()
    Call DivCorrInit()

    ! Set up model operators
    ! Set up operator arrays that only need grid geometry information
    ! discretization of del X del X E
    Call CurlcurleSetUp()

  end subroutine ModelOperatorSetUp

  !**********************************************************************
  ! Deallocate the model operators after an EM problem is finished
  subroutine ModelOperatorCleanUp()

    ! Deallocate EM equation operator arrays
    call deall_Adiag()
	call DeallocateDilu()
	call Deallocate_DivCorr()

    ! Deallocate model operators arrays
 	call CurlcurleCleanUp()

  end subroutine ModelOperatorCleanUp

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
        solver_name="QMR"
     else
        IterPerDivCor = solverControl%IterPerDivCor
        MaxDivCor = solverControl%MaxDivCor
        MaxIterTotal = MaxDivCor*IterPerDivCor
        MaxIterDivCor = solverControl%MaxIterDivCor
        tolEMfwd = solverControl%tolEMfwd
        tolEMadj = solverControl%tolEMadj
        tolDivCor = solverControl%tolDivCor
        solver_name=solverControl%solver_name
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
    character (len=1000) :: line_text
 	character(len=100),dimension(10) :: args
 	character(len=100)              :: search_string
	integer                       :: nargs
    real(kind = 8)            :: realval
    integer					  :: intval
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

solverControl%IterPerDivCor = IterPerDivCorDef
solverControl%MaxDivCor     = IterPerDivCorDef
solverControl%MaxIterDivCor = MaxIterDivCorDef
solverControl%tolEMfwd      = tolEMDef
solverControl%tolEMadj      = tolEMDef
solverControl%tolDivCor     = tolDivCorDef
solverControl%solver_name   = "QMR"
do
   read (ioFwdCtrl,"(a)",iostat=ierr) line_text ! Read line into character variable
   line_text=trim(line_text)
   if (ierr /= 0) exit                          ! End of the file

    if(index(line_text, "Number of QMR iters per divergence correction") .ne. 0)then
        call parse(line_text,":",args,nargs)
         READ(args(2),*)intval
         solverControl%IterPerDivCor=intval
         if (output_level > 2) then
           write (*,'(a12,a,i5)') node_info,"Number of QMR iters per divergence correction: ",solverControl%IterPerDivCor
          end if
   elseif(index(line_text, "Maximum number of divergence correction calls") .ne. 0)then
        call parse(line_text,":",args,nargs)
         READ(args(2),*)intval
         solverControl%MaxDivCor=intval

        if (output_level > 2) then
           write (*,'(a12,a,i5)') node_info,"Maximum number of divergence correction calls: ",solverControl%MaxDivCor
        end if
    elseif(index(line_text, "Maximum number of divergence correction iters") .ne. 0)then
        call parse(line_text,":",args,nargs)
         READ(args(2),*)intval
         solverControl%MaxIterDivCor=intval

        if (output_level > 2) then
           write (*,'(a12,a,i5)') node_info,"Maximum number of divergence correction iters: ",solverControl%MaxIterDivCor
        end if
    elseif(index(line_text, "Misfit tolerance for EM forward solver") .ne. 0)then
        call parse(line_text,":",args,nargs)
        READ(args(2),*)realval
         solverControl%tolEMfwd=realval

        if (output_level > 2) then
           write (*,'(a12,a,g15.7)') node_info,"Misfit tolerance for EM forward solver: ",solverControl%tolEMfwd
        end if
    elseif(index(line_text, "Misfit tolerance for EM adjoint solver") .ne. 0)then
        call parse(line_text,":",args,nargs)
        READ(args(2),*)realval
         solverControl%tolEMadj=realval

        if (output_level > 2) then
           write (*,'(a12,a,g15.7)') node_info,"Misfit tolerance for EM adjoint solver: ",solverControl%tolEMadj
        end if
    elseif(index(line_text, "Misfit tolerance for divergence correction") .ne. 0)then
        call parse(line_text,":",args,nargs)
        READ(args(2),*)realval
         solverControl%tolDivCor=realval

        if (output_level > 2) then
           write (*,'(a12,a,g15.7)') node_info,"Misfit tolerance for divergence correction: ",solverControl%tolDivCor
        end if
    elseif(index(line_text, "Optional EM solution file name for nested BC") .ne. 0)then
        call parse(line_text,":",args,nargs)
         solverControl%E0fileName=args(2)
         inquire(FILE=solverControl%E0fileName,EXIST=exists)
            if (.not.exists) then
                        write(*,*) node_info,'Nested E-field solution file not found and will not be used. '
                        solverControl%read_E0_from_File=.false.
            else
                        if (output_level > 2) then
                            write (*,'(a12,a,a80)') node_info,"Optional EM solution file name for nested BC: ",adjustl(solverControl%E0fileName)
                        end if
                        solverControl%read_E0_from_File=.true.
                        solverControl%ioE0=ioE
            end if
    elseif(index(line_text, "Air layers mirror|fixed height|read from file") .ne. 0)then
        solverControl%AirLayersPresent = .true.
        call parse(line_text,":",args,nargs)
        solverControl%AirLayersMethod=trim(args(2))

        if (index(solverControl%AirLayersMethod,'mirror')>0) then
            read(ioFwdCtrl,'(a48)',advance='no',iostat=istat) string
            read(ioFwdCtrl,*,iostat=istat) solverControl%AirLayersNz,solverControl%AirLayersAlpha,solverControl%AirLayersMinTopDz
        else if (index(solverControl%AirLayersMethod,'fixed height')>0) then
            read(ioFwdCtrl,'(a48)',advance='no',iostat=istat) string
            read(ioFwdCtrl,*,iostat=istat) solverControl%AirLayersNz,solverControl%AirLayersMaxHeight
        else if (index(solverControl%AirLayersMethod,'read from file')>0) then
            read(ioFwdCtrl,'(a48)',advance='no',iostat=istat) string
            read(ioFwdCtrl,'(i3)',advance='no',iostat=istat) solverControl%AirLayersNz
            allocate(solverControl%AirLayersDz(solverControl%AirLayersNz), STAT=istat)
            if (solverControl%AirLayersNz > 0) then
                read(ioFwdCtrl,*,iostat=istat) solverControl%AirLayersDz
            end if
        else
            solverControl%AirLayersPresent = .false.
            call warning('Unknown air layers method option in readEMsolveControl')
        end if

        if (solverControl%AirLayersNz <= 0) then
            write(*,*) node_info,'Problem reading the air layers. Resort to defaults '
            solverControl%AirLayersPresent = .false.
        end if
    elseif(index(line_text, "Solver QMR|BICG") .ne. 0)then
        call parse(line_text,":",args,nargs)
         solverControl%solver_name=args(2)

        if (output_level > 2) then
           write (*,'(a12,a,a)') node_info,"Solver QMR|BICG: ",solverControl%solver_name
        end if
    else
        if (output_level > 2) then
         write (*,'(a12,a,a)') node_info,"Unknown line in file: ",rFile
        end if
    end if

end do


close(ioFwdCtrl)

    call setEMsolveControl(solverControl)


   end subroutine readEMsolveControl


  !**********************************************************************
  !   deallEMsolveControl deallocate
  subroutine  deallEMsolveControl(solverControl)
     type(emsolve_control), intent(inout),optional    :: solverControl

     integer istat

     if (present(solverControl)) then
        deallocate(solverControl%AirLayersDz, STAT=istat)
     end if

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

end module EMsolve3D
