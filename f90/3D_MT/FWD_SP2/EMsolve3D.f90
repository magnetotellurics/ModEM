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
  use utilities

  implicit none
  public        :: FWDSolve3D
  public        :: deallSolverDiag, deallEMsolveControl
  public        :: createSolverDiag, getEMsolveDiag, setEMsolveControl
  private       :: SdivCorr

#ifdef FG
  interface FWDSolve3D
     MODULE PROCEDURE FWDSolve3D
     MODULE PROCEDURE FWDSolve3Dfg
  end interface
#endif

  type :: emsolve_control
    ! Values of solver control parameters, e.g., read in from file
    ! plus other information on how the solver is to be initialized, 
    ! called, etc.
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
    real(kind = 8)            ::      AirLayersMaxHeight, AirLayersAlpha
    real(kind = 8)            ::      AirLayersMinTopDz
    real(kind = 8), pointer, dimension(:)   :: AirLayersDz
    logical                   ::      AirLayersPresent=.false.
    character (len=10)        ::      solver_name="BICG"
    character (len=50) , public      ::   get_1D_from="Geometric_mean"
  end type emsolve_control

  type :: emsolve_diag
    ! Solver diagnostic arrays, computed during run of forward solver.
    ! idea is that this is the public access version of this info, which is
    ! copied from the private version in module em_solve where this info is
    ! initially stored
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
  integer, parameter    ::              IterPerDivCorDef = 150
  ! maximum number of divergence correction calls allowed
  integer, parameter    ::              MaxDivCorDef = 8
  ! maximum number of PCG iterations for divergence correction
  integer, parameter    ::              MaxIterDivCorDef = 100
  ! misfit tolerance for convergence of EMsolve algorithm
  real(kind=prec), parameter       ::      tolEMDef = 1E-10
  ! misfit tolerance for convergence of divergence correction solver
  real(kind=prec), parameter       ::      tolDivCorDef = 1E-7
  !Solver name, by default we use BICG
  character (len=10)  		   ::   solver_name="BICG"
  character (len=50) , public      ::   get_1D_from="Geometric_mean"


  save

  type(timer_t), private :: timer

  ! Actual values of control parameters must be set before first use,
  !     by call to setEMsolveControl
  !  of em_solve; are saved between calls, private to this module
  integer,  private        ::      IterPerDivCor, MaxDivCor, MaxIterDivCor
  integer,  private        ::      MaxIterTotal ! = MaxDivCor*IterPerDivCor
  real(kind=prec), private ::      tolEMfwd, tolEMadj, tolDivCor

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
!
! below is Anna's comment copied from the MF equivalent subroutine
!
! Note [AK 2018-05-10]:
! Any physical source has already been pre-multiplied by
!       [- ISIGN i\omega\mu_0]
! to yield
!       [- ISIGN i\omega\mu_0 j]
! on input to this routine. Note that this also holds for the secondary field
! formulation, where
!       j = dsigma * e,
! as well as for the tidal forcing, where
!       j = sigma (v x B).
! However, we still want to pre-compute the source RHS outside of this
! routine, for "generality".
! Specifically, Jmult supplies an interior source on the RHS that is
! not physical and is not pre-multiplied by that factor (except in Pmult).
! So it's cleaner to pass on the complete interior forcing in bRHS.
! For divergence correction, we divide by [+ ISIGN i\omega\mu_0] to get
!       i[- Div(j)].
! The plus sign is needed because we're taking the divergence of
!       curl(curl(E)) + ISIGN i\omega\mu_0 sigma E = f - curl(curl(b))
! Terms 1 and 4 vanishes, leaving:
!       Div(sigma E) - Div(f)/(+ ISIGN i\omega\mu_0) =  0.
! For a physical source j, this is equivalent to Div(sigma E) + Div(j) = 0;
! but the divergence correction may be applied also for non-physical sources,
! such  as in Jmult ('FWD') and JmultT ('TRN').
  subroutine FWDSolve3D(bRHS,omega,eSol,device_id)

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
    integer,optional,intent(in)   :: device_id ! use GPU device
    !  OUTPUTS:
    !  eSol must be allocated before calling this routine
    type (cvector), intent(inout) :: eSol

    ! LOCAL VARIABLES
    logical                     :: converged,trans
    integer                     :: iter, fid, ierr
    integer                     :: Ne,Nei,Nni,Nn,i
    complex(kind=prec)          :: iOmegaMuInv
    ! e(lectric field) s(ource) b(rhs) phi0(div(s))
    complex(kind=prec), pointer, dimension (:)     :: e,s,b
    complex(kind=prec), allocatable, dimension (:) :: ei,si,phi0
    complex(kind=prec), allocatable, dimension (:) :: temp, stemp
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
    iOmegaMuInv = C_ONE/cmplx(0.0,1.0d0*ISIGN*omega*MU_0,kind=prec) 
    if (.not.eSol%allocated) then
       call errStop('eSol in EMsolve not allocated yet')
    else
    !   determine the edge numbers of the mesh
    !   need to write a interface for these
    !   since these will be private after debugging
        Nei = size(EDGEi,1)
        Ne = size(EDGEb,1)+Nei
        if (output_level > 3) then
            write(*,'(a36,i8,a4,i8)') 'FWDSolve3D model grid #edges: Nei=', &
                Nei,' Ne=',Ne
        end if
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
    allocate(si(Nei))
    allocate(b(Nei))
    allocate(temp(Ne))
    allocate(stemp(Nei))
    ! at this point e should be all zeros if there's no initial guess
    call getCVector(eSol,e)
    if(bRHS%nonZero_Source) then ! source (TRN)
        !   this is for *all* nodes
        Nni = size(NODEi,1)
        Nn  = size(NODEb,1) + Nni
        if (output_level > 3) then
            write(*,'(a36,i8,a4,i8)') 'FWDSolve3D source grid #nodes: Nni=',&
   &                Nni,' Nn=',Nn
        end if
    ! uncomment the following line to try divergence correction in CCGD
    !    allocate(phi0(Nn)) ! make sure you *WANT* to do this, first!
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
              call getVector(tvec,s) !s is of size nEdge (all edges)
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
       ! NOTE that here we DO NOT divide the source by volume weights before
       ! the Div as the divcorr operation in SP is using VDiv instead of Div
       si = s(EDGEi) ! taking only the interior edges
       ! note that Div is formed from inner edges to all nodes
       ! uncomment the following line, to do divergence correction
       ! call Div(si,phi0)
       ! divide by iOmegaMu and Volume weight to get the source term j
       si = si * iOmegaMuInv / Vedge(EDGEi)
       ! calculate the modification term V_E GD_II j
       call RMATxCVEC(GDii,si,stemp)
       ! now i\omega\mu_0 V_E j + V_E GD_II j
       b = s(EDGEi) + stemp * Vedge(EDGEi)
    else ! trans = .false.
       ! In the usual forward model case BC does enter into forcing
       ! First compute contribution of BC term to RHS of reduced interior
       ! node system of equations : - A_IB*b
       if (bRHS%nonzero_BC) then
          !   copy from rHS structure into zeroed complex edge vector
          !   note that bRHS%bc is a cboundary type
          Call setBC(bRHS%bc, tvec) ! setBC -> copy_bcvector
          !   get info form BC
          call getVector(tvec,s)
          !   but only the boundary parts
          e(EDGEb) = s(EDGEb)
          !   Then multiply by curl_curl operator (use Mult_Aib ...
          !   Note that Mult_Aib is already multiplied by volume weights
          !   required to symmetrize the problem, so the result is V*A_IB*b)
          !   essentially b = A(i,b)*e(b)
          Call Mult_Aib(e(EDGEb), trans, b)
       endif
       ! Add internal sources if appropriate:
       ! Note that these must be multiplied explictly by volume weights
       !     [V_E^{-1} C_II^T V_F] C_II e + i\omega\mu_0\sigma e
       !   = - i\omega\mu_0  j - V_E^{-1} G_IA \Lambda D_AI j 
       !     - [V_E^{-1} C_II^T V_F] C_IB b
       ! here we multiply by V_E throughout to obtain a symmetric system:
       !      V_E A_II e = - i\omega\mu_0 V_E j - V_E GD_II j - V_E A_IB b
       ! where
       !      A_II = V_E^{-1} C_II^T V_F C_II + i\omega\mu_0 \sigma,
       ! while
       !      A_IB = V_E^{-1} C_II^T V_F C_IB, 
       ! and 
       !      GD_II = V_E{-1} G_IA \Lambda D_AI.  
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
          ! At this point, s = - ISIGN * i\omega\mu_0 j
          ! Now Div(s) - will later divide by i_omega_mu to get the general
          ! divergence correction (j)
          temp = s*Vedge
          ! uncomment the following line to do divergence correction 
          ! call Div(temp(EDGEi), phi0)
          ! now temp = - ISIGN * i\omega\mu_0 V_E j
          ! divide by iOmegaMu to get the source term (j)
          si = s(EDGEi) * iOmegaMuInv 
          ! calculate the modification term GD_II j
          call RMATxCVEC(GDii,si,stemp)
          ! i\omega\mu_0 V_E j + V_E GD_II j
          stemp = temp(EDGEi) + stemp * Vedge(EDGEi)
          ! now add the V_E A_IB b term
          if(bRHS%nonzero_BC) then
             b = stemp - b
          else
             b = stemp
          endif
      else ! there is no source
          b = -b
      endif 
    endif
    ! uncomment the following 3 lines to do divergence correction 
    ! if (bRHS%nonzero_Source) then
    !     phi0 = phi0 * iOmegaMuInv ! 1/i_omega_mu
    ! endif
    ! Outer part of KSS loop ... alternates between Calls to KSS solver
    ! and Calls to divcor  ... this will be part of EMsolve
    !
    ! e = current best solution (only on interior edges)
    ! b = rHS
    ! resetting
    nIterTotal = 0
    nDivCor = 0
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
#if defined(CUDA) || defined(HIP)
    ! FIXME this is now hard coded here
    ! need a more elegant way to deal with it
    KSSiter%maxIt = maxIterTotal
    MaxDivCor = 1
#else
    KSSiter%maxIt = IterPerDivCor
    call reset_time(timer)
#endif
    allocate(KSSiter%rerr(KSSiter%maxIt))
    KSSiter%rerr = 0.0
    converged = .false.
    failed = .false.
    !   just take the interior elements
    !   Note: e here can be used for some initial guess
    ei = e(EDGEi)
    !  idea to test: for non-zero source START with divergence
    !   correction
    ! if(bRHS%nonzero_Source) then
    !   nDivCor = 1
    !   Call SdivCorr(ei,phi0)
    ! endif
    loop: do while ((.not.converged).and.(.not.failed))
#if defined(CUDA) || defined(HIP)
       if (device_id.ge.0) then
           ! before start, need to tell if the device is available
           ierr = cf_hookDev(device_id)
           call reset_time(timer)
           if (trim(solver_name) .eq. 'PCG') then
             write(*,*) '[WARNING] GPU PCG is not yet implemented'
             write(*,*) '[WARNING] Fall back to CPU version of PCG'
             write(*,*) 'I am using PCG with initial relative error ', &
    &                KSSiter%rerr(1)
             Call PCG(b, ei, KSSiter)
           elseif (trim(solver_name) .eq. 'QMR') then
             write(*,*) '[WARNING] GPU QMR is not yet implemented '
             write(*,*) '[WARNING] Fall back to CPU version of QMR'
             write(*,*) 'I am using QMR with initial relative error ', &
    &                KSSiter%rerr(1)
             Call QMR(b, ei, KSSiter)
           elseif (trim(solver_name) .eq. 'TFQMR') then
             write(*,*) '[WARNING] GPU TFQMR is not yet implemented ', &
    &                KSSiter%rerr(1)
             write(*,*) '[WARNING] Fall back to CPU version of TFQMR'
             write(*,*) 'I am using TFQMR with initial relative error ', &
    &                KSSiter%rerr(1)
             Call TFQMR(b, ei, KSSiter)
           elseif (trim(solver_name) .eq. 'BICG') then
             write(*,*) 'I am using BICG with initial relative error ', &
    &                KSSiter%rerr(1)
             Call BiCG(b, ei, KSSiter, device_id)
           else
             write(1,*) 'ERROR: Unknown Forward Solver Method: ', &
    &                trim(solver_name)
             write(6,*) '[WARNING] Fall back to CPU version of TFQMR'
             write(6,*) 'I am using TFQMR with initial relative error ', &
    &                KSSiter%rerr(1)
             Call TFQMR(b, ei, KSSiter)
           end if
       else
           write(6,*) '[WARNING] could not find a valid GPU...'
           if (trim(solver_name) .eq. 'PCG') then
             write(6,*) '[WARNING] Fall back to CPU version of PCG'
             write(6,*) 'I am using PCG with initial relative error ', &
    &                KSSiter%rerr(1)
             Call PCG(b, ei, KSSiter)
           elseif (trim(solver_name) .eq. 'QMR') then
             write(6,*) '[WARNING] Fall back to CPU version of QMR'
             write(6,*) 'I am using QMR with initial relative error ', &
    &                KSSiter%rerr(1)
             Call QMR(b, ei, KSSiter)
           elseif (trim(solver_name) .eq. 'TFQMR') then
             write(6,*) '[WARNING] Fall back to CPU version of TFQMR'
             write(6,*) 'I am using TFQMR with initial relative error ', &
    &                KSSiter%rerr(1)
             Call TFQMR(b, ei, KSSiter)
           elseif (trim(solver_name) .eq. 'BICG') then
             write(6,*) '[WARNING] Fall back to CPU version of BICG'
             write(6,*) 'I am using BICG with initial relative error ', &
    &                KSSiter%rerr(1)
             Call BiCG(b, ei, KSSiter)
           else
             write(1,*) 'ERROR: Unknown Forward Solver Method: ', &
    &                trim(solver_name)
             write(6,*) '[WARNING] Fall back to CPU version of TFQMR'
             write(6,*) 'I am using TFQMR with initial relative error ', &
    &                KSSiter%rerr(1)
             Call TFQMR(b, ei, KSSiter)
           end if
       end if
#else
      if (trim(solver_name) .eq. 'PCG') then
        write(6,*) 'I am using PCG with initial relative error ',   &
    &           KSSiter%rerr(1)
        Call PCG(b, ei, KSSiter)
      elseif (trim(solver_name) .eq. 'QMR') then
        write(6,*) 'I am using QMR with initial relative error ',   &
    &           KSSiter%rerr(1)
        Call QMR(b, ei, KSSiter)
      elseif (trim(solver_name) .eq. 'TFQMR') then
        write(6,*) 'I am using TFQMR with initial relative error ', &
    &           KSSiter%rerr(1)
        Call TFQMR(b, ei, KSSiter)
      elseif (trim(solver_name) .eq. 'BICG') then
        write(6,*) 'I am using BICG with initial relative error ', &
    &           KSSiter%rerr(1)
        Call BiCG(b, ei, KSSiter)
      else
        write(1,*) 'ERROR: Unknown Forward Solver Method: ', trim(solver_name)
        write(6,*) '[WARNING] Fall back to TFQMR'
        write(6,*) 'I am using TFQMR with initial relative error ', &
    &                KSSiter%rerr(1)
        Call TFQMR(b, ei, KSSiter)
      end if
#endif
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
          ! uncomment the following lines to try divergence correction after
          ! solving the system matrix...
          ! if(bRHS%nonzero_Source) then
          !    Call SdivCorr(ei,phi0)     
          ! else
          !    Call SdivCorr(ei)
          ! endif
       else
          ! max number of divergence corrections exceeded; convergence failed
          failed = .true.
       endif
       if (output_level > 3) then
           write (6,*) 'iter: ', nIterTotal, ' residual: ',             &
   &    EMrelErr(nIterTotal)
       end if
    end do loop
    if (output_level > 2) then
       write (*,'(a12,a20,i8,g15.7)') node_info, 'finished solving:',     &
   &            nIterTotal, EMrelErr(nIterTotal)
       write (*,'(a12,a22,f12.6)')    node_info, 'solving time (sec): ',  &
   &            elapsed_time(timer)
    end if
    e(EDGEi) = ei
    !  After solving symetrized system, need to do different things for
    !   transposed, standard cases
    if(trans) then ! trans = .true.
       !   compute solution on boundary nodes: first  A_IB^T eSol
       call Mult_Aib(ei ,trans, s)
       !   Multiply solution on interior nodes by volume weights
       !   but after filling the solution on boundary
       temp = Vedge*e
       e = temp
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
    deallocate(e)
    deallocate(ei)
    deallocate(s)
    deallocate(si)
    deallocate(b)
    deallocate(temp)
    deallocate(stemp)
    ! uncomment the following lines for divergence correction
    ! if(bRHS%nonzero_Source) then
    !     deallocate(phi0)
    ! end if
    Call deall(tvec)
    deallocate(KSSiter%rerr)

    ! Release GPU lock and cleanup shared memory on exit
#if defined(CUDA) || defined(HIP)
    if (device_id >= 0) then
        call cf_releaseDev(device_id)
    end if
    call cf_cleanupLock()
#endif

  end subroutine FWDsolve3D

#if defined(MPI) && defined(FG)
  ! fine-grained parallel version
  subroutine FWDsolve3Dfg(bRHS,omega,eSol,device_id,comm_local)
!----------------------------------------------------------------------
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
     type (RHS_t), intent(in)           :: bRHS
     real(kind=prec), intent(in)        :: omega
     integer, intent(in)                :: device_id 
     integer, intent(in)                :: comm_local 
     !  OUTPUTS:
     !  eSol must be allocated before calling this routine
     type (cvector), intent(inout)      :: eSol
!----------------------------------------------------------------------
!               variables related to parallel computation
!----------------------------------------------------------------------
     integer, dimension(:), allocatable :: idx
     integer, dimension(:), pointer     :: isizes, isubs, displs
     integer, dimension(3)              :: iedges
     integer                            :: nrow, ncol, nnz ! spMAT parameter
     integer                            :: block_size
     integer                            :: rank_local,size_local
     integer                            :: ierr
     integer                            :: Nsub,Ntotal,istart,iend,csize
     integer, allocatable, dimension(:) :: ilocal,jlocal
     real(kind=prec)   , dimension (:), allocatable :: vlocal
     complex(kind=prec), dimension (:), allocatable :: clocal
     complex(kind=prec), dimension (:), allocatable :: xlocal
     complex(kind=prec), dimension (:), allocatable :: blocal
     real(kind=prec)                    :: ptol=1e-2
     real                               :: normu
!----------------------------------------------------------------------
     ! LOCAL VARIABLES
     logical                                        :: converged,trans
     integer                                        :: iter, fid
     integer                                        :: Ne,Nei,Neb
     integer                                        :: Nn,Nni,Nnb,i
     complex(kind=prec)                             :: iOmegaMuInv
     ! e(lectric field) s(ource) b(rhs) phi0(div(s))
     complex(kind=prec), pointer, dimension (:)     :: e,s
     complex(kind=prec), allocatable, dimension (:) :: ei,b,si,phi0,phii
     complex(kind=prec), allocatable, dimension (:) :: temp,stemp
     
     character(80)                                  :: cfile
     ! band-aid cvector ...
     type (cvector)                                 :: tmpvec
     ! diagnostic structure for Krylov Subspace Solvers(KSS)
     type (solverControl_t)                         :: KSSiter
     ! initialize solver diagnostic variables
     ! now initialize solver diagnostic variables
     nIterTotal = 0
     nDivCor = 0
     EMrelErr = R_ZERO
     divJ = R_ZERO
     DivCorRelErr = R_ZERO
     failed = .false.
     iOmegaMuInv = C_ONE/cmplx(0.0,1.0d0*ISIGN*omega*MU_0,kind=prec)
     ! now see how many workers do we have
     call MPI_COMM_RANK(comm_local,rank_local,ierr)
     call MPI_COMM_SIZE(comm_local,size_local,ierr)
     !------------------------------------------------------------------------
     !     check the input variables
     !------------------------------------------------------------------------
     if (rank_local.eq.0) then ! leader
         trans = (bRHS%adj .eq. TRN)
         if (.not.eSol%allocated) then
             call errStop('eSol in EMsolve not allocated yet')
         end if 
         ! determine the edge numbers of the mesh
         ! need to write a interface for these
         ! since these will be private after debugging
         Nei = size(EDGEi,1)
         Neb = size(EDGEb,1)
         Ne = Nei+Neb
         if(bRHS%nonZero_Source) then ! source (TRN)
             !   this is for *all* nodes
             Nni = size(NODEi,1)
             Nnb = size(NODEb,1)
             Nn = Nni+Nnb
             if (output_level > 3) then
                 write(*,'(a36,i8,a4,i8)')                         &
    &              'FWDSolve3D source grid #nodes:           Nni=',&
    &              Nni,' Nn=',Nn
             end if
         ! uncomment the following line to try divergence correction in CCGD
         ! allocate(phi0(Nn)) ! make sure you *WANT* to do this, first!
         endif
         iedges(1) = eSol%Nx*(eSol%Ny-1)*(eSol%Nz-1)
         iedges(2) = eSol%Ny*(eSol%Nx-1)*(eSol%Nz-1)
         iedges(3) = eSol%Nz*(eSol%Nx-1)*(eSol%Ny-1)
     end if
     ! common part for leader and workers
     ! broadcast those parameters to all workers
     call MPI_BCAST(trans,1, MPI_LOGICAL,0, comm_local,ierr)
     call MPI_BCAST(Nei,1, MPI_INTEGER,0, comm_local,ierr)
     call MPI_BCAST(Neb,1, MPI_INTEGER,0, comm_local,ierr)
     call MPI_BCAST(Ne ,1, MPI_INTEGER,0, comm_local,ierr)
     call MPI_BCAST(Nni,1, MPI_INTEGER,0, comm_local,ierr)
     call MPI_BCAST(Nnb,1, MPI_INTEGER,0, comm_local,ierr)
     call MPI_BCAST(Nn ,1, MPI_INTEGER,0, comm_local,ierr)
     call MPI_BCAST(iedges ,3, MPI_INTEGER,0, comm_local,ierr)
     ! calculate the sizes of each local row matrix
     ! as well as the size of block preconditioners
     ! *note* every process knows about this
     call calc_dist(size_local,iedges,isizes,isubs)
     ! now setup the displacement to send to each process 
     allocate(displs(size_local))
     displs = 0
     do i=2,size_local
         displs(i) = sum(isizes(1:i-1))
     end do
     ! output information for debug
     Ntotal = size(isubs)
     csize = 0
     Nsub = 0
     istart = 1
     do i=1,Ntotal
         csize = csize + isubs(i)
         if (csize.le.sum(isizes(1:rank_local))) then ! go on
             istart = istart + 1
         elseif (csize.le.sum(isizes(1:rank_local+1))) then 
             ! count sub-blocks in the local block
             Nsub = Nsub + 1
         else 
             exit
         end if
     end do
     iend = istart+Nsub-1
     if (output_level > 3) then
         ! for debug
         write(6,*) 'number of sub blocks =', Nsub, rank_local
         write(6,*) 'from #', istart, 'to #', iend, rank_local
     end if
     if (rank_local .eq. 0) then !leader
         ! leader does all the works to allocate/initialize 
         ! local data structures
         if (output_level > 3) then
             ! for debug
             write(6,*) 'system matrix iedges =', iedges
             write(6,*) 'system matrix isubs =', isubs
             write(6,*) 'system matrix isizes =', isizes
         endif 
         ! cboundary is a quite complex type...
         ! *essentially it should be e(EDGEb)
         ! for now we don't have an interface to deal with cboundary
         ! so just use cvectors to deliver the value...
         call create_cvector(bRHS%grid, tmpvec, eSol%gridtype)
         allocate(e(Ne))
         allocate(ei(Nei))
         allocate(s(Ne))
         allocate(si(Nei))
         allocate(b(Nei))
         allocate(temp(Ne))
         allocate(stemp(Nei))
         ! at this point e should be all zeros if there's no initial guess
         call getCVector(eSol,e)
         ! Using boundary condition and sources from rHS data structure
         ! construct vector b (defined only on interior nodes) for rHS of
         ! reduced (interior nodes only) linear system of equations
         if (trans) then ! TRN, trans=.true.
             !  In this case boundary conditions do not enter into forcing
             !  for reduced (interior node) linear system; solution on
             !  boundary is determined after solving for interior nodes
             if (bRHS%nonZero_Source) then
                 if (bRHS%sparse_Source) then
                     ! sparse source
                     ! not sure how to do it efficiently with normal array
                     ! for now it is just a walkaround, probably not going to
                     ! be used by most
                     call add_scvector(C_ONE,bRHS%sSparse,tmpvec)
                     call getVector(tmpvec,s)
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
                 if (bRHS%nonzero_BC) then
                     ! just copy input BC into boundary nodes of solution 
                     ! and return
                     Call setBC(bRHS%bc, eSol)
                 endif ! otherwise the eSol should be all zeros
                 return
             endif
             ! NOTE that here we DO NOT divide the source by volume weights
             ! before the Div as the divcorr operation in SP is using VDiv
             ! instead of Div
             si = s(EDGEi) ! taking only the interior edges
             ! note that Div is formed from inner edges to all nodes
             ! uncomment the following line, to do divergence correction
             ! call Div(si,phi0)
             ! divide by iOmegaMu and Volume weight to get the source term j
             si = si * iOmegaMuInv / Vedge(EDGEi)
             ! calculate the modification term V_E GD_II j
             call RMATxCVEC(GDii,si,stemp)
             ! now i\omega\mu_0 V_E j + V_E GD_II j
             b = s(EDGEi) + stemp * Vedge(EDGEi)
         else ! trans = .false.
             ! In the usual forward model case BC does enter into forcing
             ! First compute contribution of BC term to RHS of reduced 
             ! interior node system of equations : - A_IB*b
             if (bRHS%nonzero_BC) then
                 !   copy from rHS structure into zeroed complex edge vector
                 !   note that bRHS%bc is a cboundary type
                 Call setBC(bRHS%bc, tmpvec) ! setBC -> copy_bcvector
                 !   get info form BC
                 call getVector(tmpvec,s)
                 !   but only the boundary parts
                 e(EDGEb) = s(EDGEb)
                 !   Then multiply by curl_curl operator (use Mult_Aib ...
                 !   Note that Mult_Aib already multiplies by volume weights
                 !   required to symetrize problem, so the result is V*A_IB*b)
                 !   essentially b = A(i,b)*e(b)
                 Call Mult_Aib(e(EDGEb), trans, b)
             endif
             ! Add internal sources if appropriate: 
             ! Note that these must be multiplied explictly by volume weights
             !     [V_E^{-1} C_II^T V_F] C_II e + i\omega\mu_0\sigma e
             !   = - i\omega\mu_0  j - V_E^{-1} G_IA \Lambda D_AI j
             !     - [V_E^{-1} C_II^T V_F] C_IB b
             ! here we multiply by V_E throughout to obtain a symmetric system:
             !      V_E A_II e = - i\omega\mu_0 V_E j - V_E GD_II j - V_E A_IB b
             ! where
             !      A_II = V_E^{-1} C_II^T V_F C_II + i\omega\mu_0 \sigma,
             ! while
             !      A_IB = V_E^{-1} C_II^T V_F C_IB,
             ! and
             !      GD_II = V_E{-1} G_IA \Lambda D_AI.
             if (bRHS%nonzero_Source) then
                 if (bRHS%sparse_Source) then
                     ! sparse source
                     call zero(tmpvec)
                     call add_scvector(C_ONE,bRHS%sSparse,tmpvec)
                     call getVector(tmpvec,s)
                 else
                     ! normal source
                     call getVector(bRHS%s, s)
                 endif
                 ! At this point, s = - ISIGN * i\omega\mu_0 j
                 ! Now Div(s) - will later divide by i_omega_mu to get the 
                 ! general divergence correction (j)
                 temp = s*Vedge
                 ! now temp = - ISIGN * i\omega\mu_0 V_E j
                 ! divide by iOmegaMu to get the source (j)
                 si = s(EDGEi) * iOmegaMuInv
                 ! calculate the modification term GD_II j
                 call RMATxCVEC(GDii,si,stemp)
                 ! i\omega\mu_0 V_E j + V_E GD_II j
                 stemp = temp(EDGEi) + stemp * Vedge(EDGEi)
                 ! now add the V_E A_IB b term
                 if(bRHS%nonzero_BC) then
                     b = stemp - b
                 else
                     b = stemp
                 endif
             else! there is no source
                 b = -b
             endif
         endif ! trans 
     endif ! if rank == 0
     !-------------------------------------------------------------------------
     !    now start to decompose the system matrix
     !-------------------------------------------------------------------------
     if (rank_local.eq.0) then ! leader
         if (output_level > 2) then
             write (*,'(a12,a25,a15,i8)') node_info,'sending sub-systems to', &
   &                 ' workers...', rank_local
         endif
         ! firstly leader split the system matrix AAii (real)
         ! now calculate the rows that should be stored locally...
         ! AAii (Nei x Nei)
         do i = 2,size_local !now send those to your fellow workers
             call splitMAT(AAii,i-1,size_local,Alocal,isizes)
             nrow = size(Alocal%row)
             ncol = size(Alocal%col)
             nnz = ncol
             ! send the info to workers
             call MPI_SEND(nrow,1,MPI_INTEGER,i-1, 1,comm_local,ierr)
             call MPI_SEND(ncol,1,MPI_INTEGER,i-1, 1,comm_local,ierr)
             call MPI_SEND(nnz,1,MPI_INTEGER,i-1, 1,comm_local,ierr)
             ! now send the local matrix (in CSR format) that will be 
             ! stored in that worker process
             call MPI_SEND(Alocal%row,nrow, MPI_INTEGER,i-1,2, comm_local,ierr)
             call MPI_SEND(Alocal%col,ncol, MPI_INTEGER,i-1,2, comm_local,ierr)
             call MPI_SEND(Alocal%val,nnz, MPI_DOUBLE,i-1,2,         &
    &             comm_local, ierr)
             call deall_spMatCSR(Alocal) ! and release the temp sp matrix
         end do
         ! now deal with the local rows in leader
         call splitMAT(AAii,rank_local,size_local,Alocal,isizes)
         ! this is required as the original splitMAT is designed for PETSc
         ! (c index that starts from 0)
         Alocal%col = Alocal%col+1
         Alocal%row = Alocal%row+1
     else !workers
         if (output_level > 2) then
             write (*,'(a12,a25,a15,i8)') node_info,'receiving sub-systems', &
    &              ' from leader...', rank_local
         endif
         ! get the info from leader
         call MPI_RECV(nrow,1, MPI_INTEGER,0,1, comm_local,            &
    &        MPI_STATUS_IGNORE,ierr)
         call MPI_RECV(ncol,1, MPI_INTEGER,0,1, comm_local,            &
    &        MPI_STATUS_IGNORE,ierr)
         call MPI_RECV(nnz,1, MPI_INTEGER,0,1, comm_local,            &
    &        MPI_STATUS_IGNORE,ierr)
         ! now get the local matrix (in CSR format) 
         allocate(ilocal(nrow))
         allocate(jlocal(ncol))
         allocate(vlocal(nnz))
         call MPI_RECV(ilocal,nrow, MPI_INTEGER,0,2, comm_local,       &
    &        MPI_STATUS_IGNORE,ierr)
         call MPI_RECV(jlocal,ncol, MPI_INTEGER,0,2, comm_local,       & 
    &        MPI_STATUS_IGNORE,ierr)
         call MPI_RECV(vlocal,nnz, MPI_DOUBLE,0,2, comm_local,&
    &        MPI_STATUS_IGNORE,ierr)
         ! now assemble the matrix manually
         call create_spMatCSR_Real(isizes(rank_local+1),Nei,nnz,Alocal)
         ! this is required as the original splitMAT is designed for PETSc
         ! (c index that starts from 0)
         Alocal%row=ilocal+1
         Alocal%col=jlocal+1
         Alocal%val=vlocal
         Alocal%allocated=.true.
         deallocate(ilocal)
         deallocate(jlocal)
         deallocate(vlocal)
     end if
     !-------------------------------------------------------------------------
     !    now start to decompose the L/U matrix
     !-------------------------------------------------------------------------
     if (rank_local.eq.0) then ! leader
         ! firstly leader split the system matrix L/U (complex)
         ! now calculate the blocks that should be stored locally...
         ! L/U (Nei x Nei)
         do i = 2,size_local !now send those to your fellow workers
             ! L matrix
             call splitBlkMAT(L,i-1,size_local,Llocal,isizes)
             nrow = size(Llocal%row)
             ncol = size(Llocal%col)
             nnz = ncol
             ! send the info to workers
             call MPI_SEND(nrow,1,MPI_INTEGER,i-1, 1,comm_local,ierr)
             call MPI_SEND(ncol,1,MPI_INTEGER,i-1, 1,comm_local,ierr)
             call MPI_SEND(nnz,1,MPI_INTEGER,i-1, 1,comm_local,ierr)
             ! now send the local matrix (in CSR format) that will be 
             ! stored in that worker process
             call MPI_SEND(Llocal%row,nrow, MPI_INTEGER,i-1,2, comm_local,ierr)
             call MPI_SEND(Llocal%col,ncol, MPI_INTEGER,i-1,2, comm_local,ierr)
             call MPI_SEND(Llocal%val,nnz, MPI_DOUBLE_COMPLEX,i-1,2,         &
    &             comm_local, ierr)
             call deall_spMatCSR(Llocal) ! and release the temp sp matrix
             ! U matrix
             call splitBlkMAT(U,i-1,size_local,Ulocal,isizes)
             nrow = size(Ulocal%row)
             ncol = size(Ulocal%col)
             nnz = ncol
             ! send the info to workers
             call MPI_SEND(nrow,1,MPI_INTEGER,i-1, 1,comm_local,ierr)
             call MPI_SEND(ncol,1,MPI_INTEGER,i-1, 1,comm_local,ierr)
             call MPI_SEND(nnz,1,MPI_INTEGER,i-1, 1,comm_local,ierr)
             ! now send the local matrix (in CSR format) that will be 
             ! stored in that worker process
             call MPI_SEND(Ulocal%row,nrow, MPI_INTEGER,i-1,2, comm_local,ierr)
             call MPI_SEND(Ulocal%col,ncol, MPI_INTEGER,i-1,2, comm_local,ierr)
             call MPI_SEND(Ulocal%val,nnz, MPI_DOUBLE_COMPLEX,i-1,2,         &
    &             comm_local, ierr)
             call deall_spMatCSR(Ulocal) ! and release the temp sp matrix
         end do
         ! now deal with the local rows in leader
         call splitBlkMAT(L,rank_local,size_local,Llocal,isizes)
         call splitBlkMAT(U,rank_local,size_local,Ulocal,isizes)
     else !workers
         ! L 
         ! get the info from leader
         call MPI_RECV(nrow,1, MPI_INTEGER,0,1, comm_local,            &
    &        MPI_STATUS_IGNORE,ierr)
         call MPI_RECV(ncol,1, MPI_INTEGER,0,1, comm_local,            &
    &        MPI_STATUS_IGNORE,ierr)
         call MPI_RECV(nnz,1, MPI_INTEGER,0,1, comm_local,            &
    &        MPI_STATUS_IGNORE,ierr)
         ! now get the local matrix (in CSR format) 
         allocate(ilocal(nrow))
         allocate(jlocal(ncol))
         allocate(clocal(nnz))
         call MPI_RECV(ilocal,nrow, MPI_INTEGER,0,2, comm_local,       &
    &        MPI_STATUS_IGNORE,ierr)
         call MPI_RECV(jlocal,ncol, MPI_INTEGER,0,2, comm_local,       & 
    &        MPI_STATUS_IGNORE,ierr)
         call MPI_RECV(clocal,nnz, MPI_DOUBLE_COMPLEX,0,2, comm_local,&
    &        MPI_STATUS_IGNORE,ierr)
         ! now assemble the matrix manually
         call create_spMatCSR_Cmplx(isizes(rank_local+1),isizes(rank_local+1),&
    &        nnz,Llocal)
         Llocal%row=ilocal
         Llocal%col=jlocal
         Llocal%val=clocal
         Llocal%allocated=.true.
         ! setup the triangular flags
         Llocal%lower=.true.
         deallocate(ilocal)
         deallocate(jlocal)
         deallocate(clocal)
         ! U 
         ! get the info from leader
         call MPI_RECV(nrow,1, MPI_INTEGER,0,1, comm_local,            &
    &        MPI_STATUS_IGNORE,ierr)
         call MPI_RECV(ncol,1, MPI_INTEGER,0,1, comm_local,            &
    &        MPI_STATUS_IGNORE,ierr)
         call MPI_RECV(nnz,1, MPI_INTEGER,0,1, comm_local,            &
    &        MPI_STATUS_IGNORE,ierr)
         ! now get the local matrix (in CSR format) 
         allocate(ilocal(nrow))
         allocate(jlocal(ncol))
         allocate(clocal(nnz))
         call MPI_RECV(ilocal,nrow, MPI_INTEGER,0,2, comm_local,       &
    &        MPI_STATUS_IGNORE,ierr)
         call MPI_RECV(jlocal,ncol, MPI_INTEGER,0,2, comm_local,       & 
    &        MPI_STATUS_IGNORE,ierr)
         call MPI_RECV(clocal,nnz, MPI_DOUBLE_COMPLEX,0,2, comm_local,&
    &        MPI_STATUS_IGNORE,ierr)
         ! now assemble the matrix manually
         call create_spMatCSR_Cmplx(isizes(rank_local+1),isizes(rank_local+1),&
    &        ncol,Ulocal)
         Ulocal%row=ilocal
         Ulocal%col=jlocal
         Ulocal%val=clocal
         Ulocal%allocated=.true.
         ! setup the triangular flags
         Ulocal%upper=.true.
         deallocate(ilocal)
         deallocate(jlocal)
         deallocate(clocal)
     end if
     ! for all processes
     ! do the transopose
     call CMATtrans(Llocal,LHlocal)
     call CMATtrans(Ulocal,UHlocal)
     !-------------------------------------------------------------------------
     !    also split the RHS vectors
     !-------------------------------------------------------------------------
     ! send (part of) b to all processes
     allocate(blocal(isizes(rank_local+1)))
     call MPI_SCATTERV(b, isizes, displs, MPI_DOUBLE_COMPLEX, blocal,&
    &    isizes(rank_local+1), MPI_DOUBLE_COMPLEX, 0, comm_local, ierr)
     ! send (part of) iOmegaMuSigma to all processes
     allocate(VomegaMuSigLoc(isizes(rank_local+1)))
     call MPI_SCATTERV(VomegaMuSig, isizes, displs, MPI_DOUBLE, &
    &    VomegaMuSigLoc ,isizes(rank_local+1), MPI_DOUBLE, 0, &
    &    comm_local, ierr)
     ! send (part of) x to all processes
     if (rank_local.eq.0) then ! leader 
         ! just take the interior elements
         ei = e(EDGEi)
         ! idea to test: for non-zero source START with divergence
         !               correction
         ! if (bRHS%nonzero_Source) then
         !     nDivCor = 1
         !     Call SdivCorr(ei,phi0)
         ! endif
         ! also send the local x vector
     endif
     ! also send (part of) ei, or x to all processes
     allocate(xlocal(isizes(rank_local+1)))
     call MPI_SCATTERV(ei, isizes, displs, MPI_DOUBLE_COMPLEX, xlocal,&
    &    isizes(rank_local+1), MPI_DOUBLE_COMPLEX, 0, comm_local, ierr)
     !------------------------------------------------------------------------!
     ! Outer part of KSS loop ... alternates between Calls to KSS solver
     ! and Calls to divcor  ... this will be part of EMsolve
     ! these part should be thread safe
     !
     ! e = current best solution (only on interior edges)
     ! b = rHS
     !
     ! at present we don't really have the option to skip
     ! the divergence correction.  Not sure how/if this should
     ! be done.
     ! resetting
     nIterTotal = 0
     nDivCor = 0
     call reset_time(timer)
     ! Initialize iteration control/diagnostic structure for KSS
     KSSiter%niter = 0
#if defined(CUDA) || defined(HIP)
     ! FIXME this is now hard coded here
     ! need a more elegant way to deal with it
     KSSiter%maxIt = maxIterTotal
     MaxDivCor = 1
#else
     KSSiter%maxIt = IterPerDivCor
     call reset_time(timer)
#endif
     ! distribute the structure to all processes
     if (rank_local.eq.0) then
         if (trans) then
             KSSiter%tol = tolEMadj
         else
             if (bRHS%nonzero_BC) then
                 KSSiter%tol = tolEMfwd
             else
                 KSSiter%tol = tolEMadj
             end if
         end if
         call MPI_BCAST(KSSiter%tol,1, MPI_DOUBLE,0, comm_local,ierr)
         call MPI_BCAST(KSSiter%maxIt,1, MPI_INTEGER,0, comm_local,ierr)
         call MPI_BCAST(MaxDivCor,1, MPI_INTEGER,0, comm_local,ierr)
         call MPI_BCAST(IterPerDivCor,1, MPI_INTEGER,0, comm_local,ierr)
     else
         call MPI_BCAST(KSSiter%tol,1, MPI_DOUBLE,0, comm_local,ierr)
         call MPI_BCAST(KSSiter%maxIt,1, MPI_INTEGER,0, comm_local,ierr)
         call MPI_BCAST(MaxDivCor,1, MPI_INTEGER,0, comm_local,ierr)
         call MPI_BCAST(IterPerDivCor,1, MPI_INTEGER,0, comm_local,ierr)
     end if
     allocate(KSSiter%rerr(KSSiter%maxIt))
     KSSiter%rerr = 0.0
     converged = .false.
     failed = .false.
     if (output_level > 3) then
         write (*,'(a12,a30)') node_info,'finished transmitting data...'
     endif
!------------------------------------------------------------------------------
!     iteration starts
!------------------------------------------------------------------------------
     loop: do while ((.not.converged).and.(.not.failed))
         ! fine-grained parallel version 
#if defined(CUDA) || defined(HIP)
         if ((device_id.ge.0) .and. (size_local .gt. 1)) then
             ! before start, need to tell if the device is available
             ierr = cf_hookDev(device_id)
             if (trim(solver_name) .ne. 'BICG') then
                 write(6,*) '[WARNING] the GPU ', solver_name, &
    &                    'algorithm is not yet implemented'
                 write(6,*) '[WARNING] Fall back to BICG GPU solver'
                 write(6,*) 'I am using BICG with initial relative error ', &
    &                    KSSiter%rerr(1)
             endif
             call reset_time(timer)
             call BiCG(blocal,xlocal,KSSiter,comm_local,device_id)
         else
             write(6,*) '[WARNING] could not find enough/valid CUDA or HIP GPUs '
             write(6,*) '[WARNING] fall back to CPU FG '
             call reset_time(timer)
             call BiCG(blocal,xlocal,KSSiter,comm_local)
         endif
#else
         ! cpu only, currently only BICG is implemented
         if (trim(solver_name) .ne. 'BICG') then
             write(*,*) '[WARNING] the Fine Grain ', solver_name,           &
    &                'algorithm is not yet implemented'
             write(*,*) '[WARNING] Fall back to BICG FG solver'
             write(*,*) 'I am using BICG with initial relative error ',     &
    &                KSSiter%rerr(1)
         endif
         call reset_time(timer)
         call BiCG(blocal,xlocal,KSSiter,comm_local)
#endif
         ! algorithm is converged when the relative error is less than 
         ! tolerance
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
         if (nDivCor < MaxDivCor) then
             if (rank_local.eq.0) then ! leader
                 if (output_level > 3) then
                     ! note that the relative residual is "preconditioned"
                     ! (i.e. after applying the ILU preconditioner)
                     ! which is different from what you are experiencing 
                     ! in the SP/SP2/MF version
                     write(6,*) 'iter: ',nIterTotal,' residual: ',       &
    &                     EMrelErr(nIterTotal)
                 endif
             endif
         else
             ! max number of divergence corrections exceeded; 
             ! convergence failed
             failed = .true.
         endif
     end do loop
     if (rank_local.eq.0) then ! leader
         if (output_level > 2) then
             write (*,'(a12,a20,i8,g15.7)') node_info,'finished solving:',& 
    &             nIterTotal, EMrelErr(nIterTotal)
             write (*,'(a12,a22,f12.6)') node_info,'solving time (sec): ',&
    &             elapsed_time(timer)
         end if
     end if
     ! now retrieve the full ei from the xlocal in different processes
     call MPI_GATHERV(xlocal, isizes(rank_local+1), MPI_DOUBLE_COMPLEX, ei, &
 &           isizes, displs, MPI_DOUBLE_COMPLEX, 0, comm_local, ierr)
     ! common part for leader and workers
     call deall_spMatCSR(Alocal) ! release the temp sp matrix
     call deall_spMatCSR(Llocal) ! release the temp sp matrix
     call deall_spMatCSR(Ulocal) ! release the temp sp matrix
     call deall_spMatCSR(LHlocal) ! release the temp sp matrix
     call deall_spMatCSR(UHlocal) ! release the temp sp matrix
     if (rank_local.eq.0) then ! leader prepare the results...
         e(EDGEi) = ei
         !  After solving symetrized system, need to do different things for
         !   transposed, standard cases
         if(trans) then ! trans = .true.
            !   compute solution on boundary nodes: first  A_IB^T eSol
            call Mult_Aib(ei ,trans, s)
            !   Multiply solution on interior nodes by volume weights
            !   but after filling the solution on boundary
            temp = Vedge*e
            e = temp
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
     else !while others sit back and watch 
         if (output_level > 3) then
             write(*,'(a36,i8,a4,i8)') 'Waiting for Leader to finish: rank#', &
   &             rank_local
         end if
     end if

    ! Release GPU lock so other processes can hook on
#if defined(CUDA) || defined(HIP)
    if (device_id >= 0 .and. size_local > 1) then
        call cf_releaseDev(device_id)
    end if
#endif

     call MPI_BARRIER(comm_local,ierr)
    ! Cleanup shared-memory lock name (only rank 0)
#if defined(CUDA) || defined(HIP)
    if (rank_local .eq. 0) then
        call cf_cleanupLock()
    end if
#endif
     if (rank_local .eq. 0) then ! Leader 
         ! deallocate local temporary arrays
         deallocate(b)
         deallocate(s)
         deallocate(si)
         deallocate(e)
         deallocate(ei)
         deallocate(temp)
         deallocate(stemp)
         call deall(tmpvec)
     endif
     ! deallocate(phi0)
     deallocate(isizes)
     deallocate(isubs)
     deallocate(displs)
     deallocate(xlocal)
     deallocate(blocal)
     deallocate(VomegaMuSigLoc) 
     deallocate(KSSiter%rerr)

  end subroutine FWDsolve3Dfg
#endif

!**********************************************************************
! solver_divcorrr contains the subroutine that would solve the divergence
! correction. Solves the divergene correction using pre-conditioned
! conjugate gradient
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
  ! all nodes -> inner edges
  call Grad(phiAll,cE)

  ! subtract Divergence correction from inE
  inE = inE - cE

  ! divergence of the corrected output electrical field
  ! inner Edges -> inner Nodes
  Call DivC(inE,phiRHS)

  !  If source term is present, subtract from divergence of currents
  if(SourceTerm) then
     phiRHS = phiRHS - phi0(NODEi)
  endif

  ! as in WS code, compute the current Divergence after correction
  ! (using the dot product)
  divJ(2,nDivCor) = sqrt(dot_product(phiRHS,phiRHS))

  ! output level defined in basic file_units module
  ! write(6,*) divJ(1,nDivCor), divJ(2,nDivCor)
  if (output_level > 2) then
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
        solver_name="BICG"
        get_1D_from="Geometric_mean"
     else
        IterPerDivCor = solverControl%IterPerDivCor
        MaxDivCor = solverControl%MaxDivCor
        MaxIterTotal = MaxDivCor*IterPerDivCor
        MaxIterDivCor = solverControl%MaxIterDivCor
        tolEMfwd = solverControl%tolEMfwd
        tolEMadj = solverControl%tolEMadj
        tolDivCor = solverControl%tolDivCor
        solver_name=solverControl%solver_name
        get_1D_from=solverControl%get_1D_from
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


    read(ioFwdCtrl,'(a48)',advance='no',iostat=istat) string
    read(ioFwdCtrl,'(a10)',iostat=istat) solverControl%solver_name
    if (istat .ne. 0) then
       solverControl%solver_name = 'QMR' ! default
    elseif (output_level > 2) then
       write (*,'(a12,a48,a)') node_info,string,adjustl(solverControl%solver_name)
    end if

    ! For any secondary field calculation approach...
    read(ioFwdCtrl,'(a48)',advance='no',iostat=istat) string
    read(ioFwdCtrl,'(a50)',iostat=istat) solverControl%get_1D_from
    if (istat .ne. 0) then
       solverControl%get_1D_from = 'Geometric_mean' ! default
    elseif (output_level > 2) then
       write (*,'(a12,a48,a)') node_info,string,adjustl(solverControl%get_1D_from)
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
