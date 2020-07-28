!**********************************************************************
! driver modules for solving the forward EM problem, including setup and
! solver
! this module is modified to use sparse matrix

module EMsolve3D
  use sg_boundary! work between different data types
 ! (between boundary conditions and
! complex vectors)
  use sg_sparse_vector, only: add_scvector
  use modelOperator3D  ! Maxwell operator module for sp
  use vectranslate     ! translate back and forth between Cvec and vec
  use solver           ! generic solvers rewrite for sp
  use solnspace

  implicit none
  public        :: FWDSolve3D
  public        :: deallSolverDiag, deallEMsolveControl
  public        :: createSolverDiag, getEMsolveDiag, setEMsolveControl
  private       :: SdivCorr

  type :: emsolve_control
    ! Values of solver control parameters, e.g., read in from file
    ! plus other information on how the solver is to be initialized, called, etc.
    ! idea is that this is the public access version of this info, which is
    ! copied into private version for actual solver control
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
  end type emsolve_control

  type :: emsolve_diag
    ! Solver diagnostic arrays, computed during run of forward solver.
    !  idea is that this is the public access version of this info, which is
    !   copied from the private version in module em_solve where this info is
    !   initially stored
    logical                   :: diagOut
    character (len=80)        :: fn_diagn
    integer                   :: ioDiag
    integer                   :: nIterTotal, nDivCor
    real(kind = 8), pointer, dimension(:)      ::      EMrelErr
    real(kind = 8), pointer, dimension(:,:)    ::      divJ
    real(kind = 8), pointer, dimension(:,:)    ::      DivCorRelErr
  end type emsolve_diag

  ! Default solver control parameters
  ! number of iterations for each call to divergence correction:
  integer, parameter    ::              IterPerDivCorDef = 40
  ! maximum number of divergence correction calls allowed
  integer, parameter    ::              MaxDivCorDef = 20
  ! maximum number of PCG iterations for divergence correction
  integer, parameter    ::              MaxIterDivCorDef = 100
  ! misfit tolerance for convergence of EMsolve algorithm
  real(kind=prec), parameter       ::      tolEMDef = 1E-8
  ! misfit tolerance for convergence of divergence correction solver
  real(kind=prec), parameter       ::      tolDivCorDef = 1E-7

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
  integer, private   :: nIterTotal, nDivCor
  logical, private   :: failed
  ! nIterTotal keeps tally on number of iterations so far
  ! nDivCor keeps tally on number of divergence correction so far
  real(kind=prec), pointer, dimension(:), private	::	EMrelErr
  real(kind=prec), pointer, dimension(:,:), private	::	divJ
  real(kind=prec), pointer, dimension(:,:), private	::	DivCorRelErr

Contains

!**********************************************************************
! main subroutine to Solve the forward EM problem;
! modified to use the sparse matrix data structure defined in
! sp modelOperator3D module
! If bRHS%adj = 'TRN' solves transposed problem  A^T x = b
  subroutine FWDsolve3D(bRHS,omega,eSol,comm_local)

    ! redefine some of the interfaces (locally) for our convenience
    use sg_vector !, only: copy => copy_cvector, &
    use vectranslate     ! translate back and forth between Cvec and vec
    use solver
    use spoptools
    ! generic routines for vector operations on the edge/face nodes
    ! in a staggered grid
    ! in cvec copy, remember the order is copy(new, old) i.e new = old

    implicit none
    !  INPUTS:
    type (RHS_t), intent(in)      :: bRHS
    real(kind=prec), intent(in)   :: omega
    !dummy parameter for compatibiliy
    integer, intent(in),optional    :: comm_local 
    !  OUTPUTS:
    !  eSol must be allocated before calling this routine
    type (cvector), intent(inout) :: eSol

    ! LOCAL VARIABLES
    logical                     :: converged,trans
    integer                     :: iter, fid
    integer                     :: Ne,Nei,Nni,Nn,i
    complex(kind=prec)          :: iOmegaMuInv
    ! e(lectric field) s(ource) b(rhs) phi0(div(s))
    complex(kind=prec), pointer, dimension (:) :: e,s,b
    complex(kind=prec), allocatable, dimension (:) :: ei,phi0,phii
    complex(kind=prec), allocatable, dimension (:) :: temp,etemp
    character(80)                                  :: cfile
    !  band-aid cvector ...
    type (cvector)              :: tvec
    !  diagnostic structure for Krylov Subspace Solvers(KSS)
    type (solverControl_t)      :: KSSiter
    !  initialize solver diagnostic variables
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
    else
    !   determine the edge numbers of the mesh
    !   need to write a interface for these
    !   since these will be private after debugging
        Nei = size(EDGEi,1)
        Ne = size(EDGEb,1)+Nei
    end if
    ! allocate/initialize local data structures
    ! cboundary is a quite complex type...
    ! *essentially it should be e(EDGEb)
    ! for now we don't have an interface to deal with cboundary
    ! so just use cvectors to deliver the value...
    call create_cvector(bRHS%grid, tvec, eSol%gridtype)
    allocate(e(Ne))
    allocate(ei(Nei))
    allocate(s(Ne))
    allocate(b(Nei))
    allocate(temp(Ne))
    ! at this point e should be all zeros if there's no initial guess
    call getCVector(eSol,e)
    if(bRHS%nonZero_Source) then ! source (TRN)
        !   this is for *all* cells
        Nni = size(NODEi,1)
        Nn  = size(NODEb,1) + Nni
        allocate(phi0(Nn))
    endif
    ! Using boundary condition and sources from rHS data structure
    ! construct vector b (defined only on interior nodes) for rHS of
    ! reduced (interior nodes only) linear system of equations

    if(trans) then ! TRN, trans=.true.
       !  In this case boundary conditions do not enter into forcing
       !    for reduced (interior node) linear system; solution on
       !    boundary is determined after solving for interior nodes
       if (bRHS%nonZero_Source) then
          if (bRHS%sparse_Source) then
             ! sparse source
             ! not sure how to do it efficiently with normal array
             ! for now it is just a walkaround, probably not going to
             ! be used by most
              call add_scvector(C_ONE,bRHS%sSparse,tvec)
              call getVector(tvec,s)
          else
             ! normal source
              call getVector(bRHS%s,s)
          endif
       else
          ! doesn't need to tamper with this part for sparse matrix
          ! let it go with cvectors...
          call zero(eSol)
          write(0,*) 'Warning: no sources for adjoint problem'
          write(0,*) 'Solution is identically zero'
          if(bRHS%nonzero_BC) then
          ! just copy input BC into boundary nodes of solution and return
             Call setBC(bRHS%bc, eSol)
          endif ! otherwise the eSol should be all zeros
          return
       endif
       !  note that Div is formed from inner edges to all nodes
       b = s(EDGEi) ! taking only the interior edges
       call Div(b,phi0)
    else ! trans = .false.
       ! In the usual forward model case BC do enter into forcing
       !   First compute contribution of BC term to RHS of reduced interior
       !    node system of equations : - A_IB*b
       if (bRHS%nonzero_BC) then
          !   copy from rHS structure into zeroed complex edge vector
          !   note that bRHS%bc is a cboundary type
          Call setBC(bRHS%bc, tvec) ! setBC -> copy_bcvector
          !   get info form BC
          call getVector(tvec,s)
          !   but only the boundary parts
          e(EDGEb) = s(EDGEb) ! note the edgeb is not changed since
          !   Then multiply by curl_curl operator (use Mult_Aib ...
          !   Note that Mult_Aib already multiplies by volume weights
          !   required to symetrize problem, so the result is V*A_IB*b)
          !   essentially b = A(i,b)*e(b)
          Call Mult_Aib(e(EDGEb), trans, b)
       endif
       ! Add internal sources if appropriate: Note that these must be multiplied
       !  explictly by volume weights
       if (bRHS%nonzero_Source) then
          if (bRHS%sparse_Source) then
             ! sparse source
             call zero(tvec)
             call add_scvector(C_ONE,bRHS%sSparse,tvec)
             call getVector(tvec,s)
          else
             ! normal source
             call getVector(bRHS%s, s)
          endif
          call Div(s,phi0)
          temp = Vedge*s
          if(bRHS%nonzero_BC) then
             b = temp(EDGEi) - b
          else
             b = temp(EDGEi)
          endif
       else
       !   no source
           b = -b
       endif
    endif
    if(bRHS%nonzero_Source) then
       ! leave it as it was
       iOmegaMuInv = ISIGN/cmplx(0.0,omega*MU_0,prec) ! original
       phi0 = phi0 * iOmegaMuInv
    endif
    ! Outer part of KSS loop ... alternates between Calls to KSS solver
    ! and Calls to divcor  ... this will be part of EMsolve
    !
    ! e = current best solution (only on interior edges)
    ! b = rHS
    !
    ! resetting
    nIterTotal = 0
    nDivCor = 0
    call reset_time(timer)
    ! Initialize iteration control/diagnostic structure for KSS
    if (trans) then
       KSSiter%tol = tolEMadj
    else
      if (bRHS%nonzero_BC) then
        KSSiter%tol = tolEMfwd
      else
        KSSiter%tol = tolEMadj
      end if
    end if

    KSSiter%niter = 0
    KSSiter%maxIt = IterPerDivCor
    allocate(KSSiter%rerr(IterPerDivCor))
    KSSiter%rerr = 0.0
    converged = .false.
    failed = .false.
    !   just take the interior elements
    ei = e(EDGEi)
    !  idea to test: for non-zero source START with divergence
    !   correction
    if(bRHS%nonzero_Source) then
       nDivCor = 1
       Call SdivCorr(ei,phi0)
    endif
    loop: do while ((.not.converged).and.(.not.failed))
       ! Call BICG(b, ei, KSSiter)
       Call QMR(b, ei, KSSiter)
       ! algorithm is converged when the relative error is less than tolerance
       ! (in which case KSSiter%niter will be less than KSSiter%maxIt)
       converged = KSSiter%niter .lt. KSSiter%maxIt
       ! there are two ways of failing:
       !    1) the specific KSS did not work or
       !    2) total number of divergence corrections exceeded
       failed = failed .or. KSSiter%failed
       ! update diagnostics output from KSS
       do iter = 1,KSSiter%niter
           EMrelErr(nIterTotal+iter) = KSSiter%rerr(iter)
       end do
       if (KSSiter%niter.eq.0) then ! in case a initial guess is good enough
           KSSiter%niter = 1
           EMrelErr(KSSiter%niter) = KSSiter%rerr(1)
       endif
       nIterTotal = nIterTotal + KSSiter%niter
       nDivCor = nDivCor+1
       if( nDivCor < MaxDivCor) then
          ! do divergence correction
          if(bRHS%nonzero_Source) then
             Call SdivCorr(ei,phi0)
          else
             Call SdivCorr(ei)
          endif
       else
          ! max number of divergence corrections exceeded; convergence failed
          failed = .true.
       endif
       if (output_level > 3) then
           write (6,*) 'iter: ', nIterTotal, ' residual: ',                 &
   &    EMrelErr(nIterTotal)
       end if
    end do loop
    if (output_level > 2) then
       write (*,'(a12,a20,i8,g15.7)') node_info, 'finished solving:',         &
   &    nIterTotal, EMrelErr(nIterTotal)
       write (*,'(a12,a22,f12.6)')    node_info, ' time taken (mins) ',       &
   &    elapsed_time(timer)/60.0
    end if
    e(EDGEi) = ei
    !  After solving symetrized system, need to do different things for
    !   transposed, standard cases
    if(trans) then ! trans = .true.
       !   compute solution on boundary nodes: first  A_IB^T eSol
!       write(fid) Nei
!       write(fid) ei
!       close(fid)
       call Mult_Aib(ei ,trans, s)
       !   Multiply solution on interior nodes by volume weights
       !   but after filling the solution on boundary
       e = Vedge*e
       ! then b - A_IB^T eSol, where b is input boundary values (if any)
       if(bRHS%nonzero_BC) then
           e(EDGEb) = e(EDGEb) - s(EDGEb)
       else
           e(EDGEb) = -s(EDGEb)
       endif
    else ! trans = .false.
       ! just copy input BC into boundary nodes of solution
       if(.not.bRHS%nonzero_BC) then
           e(EDGEb) = 0
       endif
    endif
    call setVector(e,eSol)
    ! deallocate local temporary arrays
    deallocate(b)
    deallocate(s)
    deallocate(e)
    deallocate(ei)
    deallocate(temp)
    if(bRHS%nonzero_Source) then
        deallocate(phi0)
    end if
    Call deall(tvec)
    deallocate(KSSiter%rerr)

  end subroutine FWDsolve3D

!**********************************************************************
! solver_divcorrr contains the subroutine that would solve the divergence
! correction. Solves the divergene correction using pre-conditioned
! conjuagte gradient
subroutine SdivCorr(inE,phi0)
  ! Purpose: driver routine to compute divergence correction for input/output
  ! electric field vector inE
  ! Optional argument phi0 is scaled divergence of source term
  ! to be subtracted from current divergence

  use modelOperator3D
  ! rename routines for linear algebra operations; change to apply
  ! PCG to different problem

  implicit none

  complex(kind=prec), intent(inout),        dimension (:) :: inE
  complex(kind=prec), intent(in), optional, dimension (:) :: phi0

  !  local variables
  type (solverControl_t)                         :: PCGiter
  complex(kind=prec), allocatable, dimension (:) :: phiSol, phiRHS, phiAll
  complex(kind=prec), allocatable, dimension (:) :: cE
  integer                                        :: Nni,Nn,Nei,fid
  logical                                        :: SourceTerm

  SourceTerm = present(phi0)

  ! initialize PCGiter (maximum iterations allowed per set of diveregence
  ! correction, error tolerence, and relative error book keeping)
  PCGiter%maxIt = MaxIterDivCor
  PCGiter%tol = tolDivCor
  allocate(PCGiter%rerr(PCGiter%maxIt))
  PCGiter%rerr = 0.0

  ! alocating phiSol, phiRHS, phiAll
  Nni = size(NODEi)
  Nn = size(NODEb)+Nni
  Nei = size(inE)
  allocate(phiSol(Nni))
  allocate(phiRHS(Nni))
  allocate(phiAll(Nn))
  allocate(cE(Nei))
  phiSol=C_ZERO ! some initializations
  phiAll=C_ZERO

  ! compute divergence of currents for input electric field
  ! phi = V*Div*sigma*e
  call DivC(inE,phiRHS)
  !  If source term is present, subtract from divergence of currents
  if(SourceTerm) then
     phiRHS = phiRHS - phi0(NODEi)
  endif
  ! compute current Divergence before correction (using dot product)
  ! note this is the phi with volumn weights - not really just divJ
  divJ(1,nDivCor) = sqrt(dot_product(phiRHS,phiRHS))
  if (divJ(1,nDivCor).eq.0.0)  then
      ! initial current divergence is already zero
      ! no point to do the correction
      divJ(2,nDivCor) = divJ(1,nDivCor)
      if (output_level > 2) then
          write (*,'(a12,a25,g15.7,a25)')  node_info,                   &
    &           'current divergence is: ', divJ(2,nDivCor),             &
    &           'no need to do correction'
      endif
      return
  endif
  ! PCG is a (sort of) generic pre-conditioned CG algorithm
  Call PCG(phiRHS,phiSol,PCGiter)
  DivCorRelErr(:,nDivCor) = PCGiter%rerr
  if (output_level > 2) then
     write (*,'(a12,a32,i5,g15.7)') node_info,                            &
        'finished divergence correction:', PCGiter%niter,                 &
         PCGiter%rerr(PCGiter%niter)
  end if
  ! compute gradient of phiSol (Divergence correction for inE)
  phiAll(NODEi) = phiSol
  call Grad(phiAll,cE)

  ! subtract Divergence correction from inE
  inE = inE - cE

  ! divergence of the corrected output electrical field
  Call DivC(inE,phiRHS)

  !  If source term is present, subtract from divergence of currents
  if(SourceTerm) then
     phiRHS = phiRHS - phi0(NODEi)
  endif

  ! as in WS code, compute the current Divergence after correction
  ! (using the dot product)
  divJ(2,nDivCor) = sqrt(dot_product(phiRHS,phiRHS))

  ! output level defined in basic file_units module
  if (output_level > 3) then
     write(*,'(a12,a47,g15.7)') node_info, 'divergence of currents before correction: ', divJ(1, nDivCor)
     write(*,'(a12,a47,g15.7)') node_info, 'divergence of currents  after correction: ', divJ(2, nDivCor)
  end if

  ! deallocate the temporary work arrays
  deallocate(phiRHS)
  deallocate(phiSol)
  deallocate(phiAll)
  deallocate(cE)
  deallocate(PCGiter%rerr)

end subroutine SdivCorr ! SdivCorr

  !**********************************************************************
  ! Bundles the Inits that are used for an EM problem. These Inits can be
  ! used separately as well.
  subroutine ModelOperatorSetUp()

    ! COPIED OVER FROM MATRIX FREE VERSION FOR CONSISTENCY. CURRENTLY NOT USED
    ! BUT NECESSARY FOR REUSE OF THE SAME FORWARD SOLVER MODULE

  end subroutine ModelOperatorSetUp

  !**********************************************************************
  ! Deallocate the model operators after an EM problem is finished
  subroutine ModelOperatorCleanUp()

    ! COPIED OVER FROM MATRIX FREE VERSION FOR CONSISTENCY. CURRENTLY NOT USED
    ! BUT NECESSARY FOR REUSE OF THE SAME FORWARD SOLVER MODULE


  end subroutine ModelOperatorCleanUp

  !**********************************************************************
  ! setEMsolveControl sets actual solver control parameters, using info
  !  in structure solverControl, and allocates diagnostic arrays
  subroutine setEMsolveControl(solverControl,tolEM)

     type (emsolve_control), intent(in) ::  solverControl
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

    type(emsolve_control), intent(inout)    :: solverControl
    character(*), intent(in)                :: rFile
    logical, intent(out), optional          :: fileExists
    real(8), intent(in), optional           :: tolEM
    integer                                 :: ios
	logical                             :: exists
    character(80)                           :: string
    integer                                 :: istat

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



    ! Check if there is an additional line.
    ! if yes, it corresponds to the larger E field solution.

    read(ioFwdCtrl,'(a48)',advance='no',iostat=istat) string
    read(ioFwdCtrl,'(a80)',iostat=istat) solverControl%E0fileName
    if (istat .eq. 0 ) then
        inquire(FILE=solverControl%E0fileName,EXIST=exists)
        if (index(string,'#')>0) then
            ! This is a comment line
            solverControl%read_E0_from_File=.false.
        else if (.not.exists) then
            write(*,*) node_info,'Nested E-field solution file not found and will not be used. '
            solverControl%read_E0_from_File=.false.
        else
            if (output_level > 2) then
                write (*,'(a12,a48,a80)') node_info,string,adjustl(solverControl%E0fileName)
            end if
            solverControl%read_E0_from_File=.true.
            solverControl%ioE0=ioE
        end if

    else
        solverControl%read_E0_from_File=.false.
    end if

    ! Now keep on reading for the air layers info
    ! If the number of air layers conflicts with that from the model file, we
    ! update the grid to use the controls
    ! Method options are: mirror; fixed height; read from file
    ! For backwards compatibility, defaults to what was previously hardcoded

    read(ioFwdCtrl,'(a48)',advance='no',iostat=istat) string
    read(ioFwdCtrl,'(a80)',iostat=istat) solverControl%AirLayersMethod
    !if (output_level > 2) then
    !    write (*,'(a12,a48,a80)') node_info,string,adjustl(solverControl%AirLayersMethod)
    !end if
    if (istat .eq. 0 ) then
        solverControl%AirLayersPresent = .true.
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
    else
        solverControl%AirLayersPresent = .false.
    end if

    if (solverControl%AirLayersNz <= 0) then
        write(*,*) node_info,'Problem reading the air layers. Resort to defaults '
        solverControl%AirLayersPresent = .false.
    end if

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
