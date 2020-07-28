!   this module implements model operators as sparse matrix multiplication
!   This will replace modelOperator3D in the original ModEM implementation
!   Through extension and modifcation of this module (need to work out how
!    best to do this, alternative forward problem implementations
!    (spherical coordinates, multiresolution, etc.) can be developed
!    NOTE:  many variants might require minimum modifications (e.g., use
!      different modules, with the same method names)
!      Others may require a complete rewriting of this module

module modelOperator3D
   use vecTranslate
   use spOpTools
   use spOpTopology_SG
   use MetricElements_SG
   use ModelSpace
   ! temporarily use the following modules here
   use sg_boundary ! seems ironic in a sparse matrix module
   use nestedem
   use boundary_ws

   implicit none
   save
 !  private
   !   make all public for initial testing
   public

   !!!!!!!>>>>>>>>> conductivities on cells and edges (private)
   !   some of these might need to be public
   type(rscalar)           :: sigma_C ! cells - needed for boundary conditions
   type(grid_t), target    :: mGrid   ! THE model grid
   real(kind=prec)         :: omega   ! THE (active) frequency
   integer, allocatable, dimension(:) :: EDGEi, EDGEb
   integer, allocatable, dimension(:) :: NODEi, NODEb

   !   curl-curl operator
   type(spMatCSR_Real)     :: CCii,CCib  !   sparse matrix representation
                                         !    of curl-curl operator
                                         !  interior/interior and int./bdry
   real(kind=prec),allocatable  :: VomegaMuSig(:) !    diagonal component of
                                                  !  symetrized operator
   type(spMatCSR_Cmplx)         :: L, U    !   upper and lower triangular
   type(spMatCSR_Cmplx)         :: LH,UH   !  matrices for preconditioner
                                           !   based on interior edges only
                                           !   LH, UH are Hermitian Conjugate
                                           !   transpose of L and  U
   !   divergence correction operators
   type(spMatCSR_Real)         ::  VDiv  ! div : edges->nodes (interior only)
   type(spMatCSR_Real)         ::  VDsG  !   operator for div correction
   type(spMatCSR_Real)         ::  VDs   !     divergence of current operator
   type(spMatCSR_Real)         ::  VDsG_L  !  preconditioner for div cor
   type(spMatCSR_Real)         ::  VDsG_U  !  preconditioner for div cor
   complex(kind=prec),allocatable,dimension(:)   :: tempPhi

   public       :: ModelDataInit,ModelDataCleanup,UpdateFreq,   &
                   UpdateCond, UpdateFreqCond, DivCgrad, Grad, Div, &
                   DivC , Mult_Aii, Mult_Aib, PC_Lsolve, PC_Usolve, DivCgradILU

Contains

   subroutine ModelDataInit(inGrid)
   !**********************************************************************
   !   this should do all of the setup that can be done once and for all
   !    once grid is defined.  The setup done here is more extensive than
   !     in the original ModEM matrix free implementation.  This includes
   !     setting all metric elements and basic operator topologies
   !**********************************************************************

      implicit none
      !  INPUTS:
      type (grid_t), intent(in)             :: inGrid
      integer                               :: nz,nInterior

      !   copy inGrid to mGrid :  do we really need a copy?
      call copy_grid(mGrid,inGrid)

      !   set sparse matrices for curl (T) and grad (G)
      !    operator topolgies; these sparse matrices are stored
      !    in module spOpTopology
      call setCurlTopology(mGrid)
      call setGradTopology(mGrid)
      call boundaryIndex(EDGE,mGrid,EDGEb,EDGEi)
      nInterior = size(EDGEi)
      !   find indicies (in vector of all) of boundary and interior edges
      !   allocate for diagonal part of curl-curl operator
      !     (maybe this should just be for interior edges)
      !      here for all edges
      allocate(VomegaMuSig(nInterior))

      !    set metric elements
      !     these are stored as real arrays (all elements, including
      !         on boundaries) in module MetricElements
      call setFaceArea(mGrid)
      call setDualFaceArea(mGrid)
      call setEdgeLength(mGrid)
      call setDualEdgeLength(mGrid)
      call setVnode(mGrid)
      call setVedge(mGrid)

      ! Use the grid (which, potentially, maybe have been updated!) to set up
      !   all the grid length, surface and volume elements stored in GridCalc.
      ! Want to initialize them here in case the grid gets updated along the way.
      ! The reason for storing them in GridCalc is that they are also used
      !   by ModelMap, EMfieldInterp, nestedEM
      ! THIS IS CLEARLY A CONCEPTUAL DUPLICATE OF MetricElements but including
      ! it here for now for compatibility - needed for EMfieldInterp
      Call EdgeLength(mGrid, l_E)
      Call EdgeArea(mGrid, S_E)
      Call EdgeVolume(mGrid, V_E, l_E, S_E)
      Call FaceLength(mGrid, l_F)
      Call FaceArea(mGrid, S_F)
      Call FaceVolume(mGrid, V_F, l_F, S_F)
      Call CellVolume(mGrid, V_C)
      Call NodeVolume(mGrid, V_N) ! used for divergence correction

      ! set a default omega
      omega = 0.0
      !  specific model operators
      call CurlCurlSetup()
      call DivCorInit()

   end subroutine ModelDataInit
   !**********************************************************************
   subroutine ModelDataCleanUp()

      ! Deallocated the grid
      call deall_grid(mGrid)
      !    and the curl and grad topology matrices
      call deall_spMatCSR(T)
      call deall_spMatCSR(G)
      call CurlCurlCleanup()
      call deall_PC()
      call deall_DivCor()

      !    interior and edge indicies
      deallocate(EDGEi)
      deallocate(EDGEb)
      deallocate(NODEi)
      deallocate(NODEb)
      !    and the metric elements
      call deall_MetricElements()

      ! and the grid elements stored in GridCalc
      call deall_rvector(l_E)
      call deall_rvector(S_E)
      call deall_rvector(V_E)
      call deall_rvector(l_F)
      call deall_rvector(S_F)
      call deall_rvector(V_F)
      call deall_rscalar(V_C)
      call deall_rscalar(V_N)

      ! and the edge conductivities
      if(allocated(VomegaMuSig)) then
         deallocate(VomegaMuSig)
      endif

      ! and the cell conductivities
      ! note that sigma_C is only needed to set up boundary conditions
      call deall_rscalar(sigma_C)

   end subroutine ModelDataCleanUp
   ! **************************************************************************
   ! * UpdateFreq updates the frequency that is currently being use
   subroutine UpdateFreq(inOmega)

      real (kind=prec), intent (in)             :: inOmega

      call updateOmegaMuSig(inOmega)
      call PC_setup()
      ! this is also needed as new omega is also used for div correction
      call DivCorSetup()

   end subroutine UpdateFreq  ! UpdateFreq
   ! ***************************************************************************
   ! * UpdateCond updates the conductivity values on the edges
   subroutine UpdateCond(CondParam)

       implicit none
       type(modelParam_t), intent(in)      :: CondParam      ! input conductivity
       type(rvector)                       :: sigTemp
       real(kind=prec), pointer, dimension(:)    :: sigVec

       ! COPIED OVER FROM MATRIX FREE VERSION FOR CONSISTENCY. CURRENTLY NOT USED
       ! BUT NECESSARY FOR REUSE OF THE SAME FORWARD SOLVER MODULE
       if(associated(sigVec)) then
           nullify(sigVec)
       endif
       call create(mGrid,sigTemp,EDGE)
       call ModelParamToEdge(CondParam,sigTemp)
       call getVector(sigTemp,sigVec)
       call deall(sigTemp)
       omega = ONE ! setup an (arbitary) working omega
       VomegaMuSig = MU_0*omega*sigVec(EDGEi)*Vedge(EDGEi)
       deallocate(sigVec)
      ! TEMPORARY; REQUIRED FOR BOUNDARY CONDITIONS
      !  set static array for cell conductivities
      !  this stores conductivity values in a module structure
      !  that is readily accesible to boundary condition routines
      !  rvector sigma_C is created if it is not yet allocated
       call ModelParamToCell(CondParam, sigma_C)

   end subroutine UpdateCond  ! UpdateCond
   ! ***************************************************************************
   ! * UpdateFreqCond updates fequency and conductivity  on edges
   subroutine UpdateFreqCond(inOmega,CondParam)

      implicit none
      type(modelParam_t), intent(in)      :: CondParam      ! input conductivity
      real (kind=prec), intent (in)       :: inOmega

      call updateOmegaMuSig(inOmega,CondParam)
      call PC_setup()
      call DivCorSetup()

      ! TEMPORARY; REQUIRED FOR BOUNDARY CONDITIONS
      !  set static array for cell conductivities
      !  this stores conductivity values in a module structure
      !  that is readily accesible to boundary condition routines
      !  rvector sigma_C is created if it is not yet allocated
      call ModelParamToCell(CondParam, sigma_C)

   end subroutine UpdateFreqCond  !
   ! **************************************************************************
   subroutine  updateOmegaMuSig(inOmega,m)
   !   updates VomegaMuSig; private to this module

      real (kind=prec), intent (in)       :: inOmega
      type(modelParam_t), intent(in),optional :: m
      type(rvector)      :: sigTemp
      real(kind=prec), pointer, dimension(:)    :: sigVec

      if(present(m)) then
         call create(mGrid,sigTemp,EDGE)
         call ModelParamToEdge(m,sigTemp)
         if(associated(sigVec)) then
             nullify(sigVec)
         endif
         call getVector(sigTemp,sigVec)
         call deall(sigTemp)
         VomegaMuSig = MU_0*inOmega*sigVec(EDGEi)*Vedge(EDGEi)
         deallocate(sigVec)
         omega = inOmega
      else
         if(omega.gt.0) then
            VomegaMuSig = VomegaMuSig/omega
         endif
         VomegaMuSig = VomegaMuSig*inOmega
         omega = inOmega
      endif
   end subroutine updateOmegaMuSig

!**********************************************************************
! Sets boundary conditions. Currently a wrapper for BC_x0_WS.
! Uses input 3D conductivity in cells sigma_C, that has to be initialized
! by updateCond before calling this routine. Also uses mGrid set by
! ModelDataInit. Uses omega, which is set by updateFreq.
! We always run this after setting the private variable omega, anyway.
! BOUNDARY CONDITIONS SHOULD BE DONE DIFFERENTLY FOR SPARSE MATRIX VERSION.
! FOR NOW, COPYING THE LATEST FIX FROM MATRIX FREE VERSION. [AK]
  Subroutine ComputeBC(iTx,imode,E0,BC)

    !  Input mode, period
    integer, intent(in)     :: imode
    integer, intent(in)     :: iTx

    ! local variable
    real(kind=prec) :: period

    ! Output electric field first guess (for iterative solver)
    type(cvector), intent(inout)    :: E0
    ! Output boundary conditions
    type(cboundary), intent(inout)  :: BC

    period = (2*PI)/omega ! period is seconds

    ! Compute E0 using Weerachai 2D approach; can get BC from that
    call BC_x0_WS(imode,period,mGrid,sigma_C,E0,BC)
    
    !call getBC(E0,BC)
   
    ! Cell conductivity array is no longer needed
    ! NOT TRUE: needed for imode=2
    ! call deall_rscalar(sigma_C)

  end subroutine ComputeBC

!**********************************************************************
! Sets boundary conditions. Currently a wrapper for BC_x0_WS.
! Uses input 3D conductivity in cells sigma_C, that has to be initialized
! by updateCond before calling this routine. Also uses mGrid set by
! ModelDataInit. Uses omega, which is set by updateFreq.
! We always run this after setting the private variable omega, anyway.
!  Subroutine SetBound(imode,E0,BC,iTx)
!
!    !  Input mode, period
!    integer, intent(in)		:: imode
!    integer, intent(in)		:: iTx
!
!    ! local variable
!    real(kind=prec)	:: period
!
!    ! Output electric field first guess (for iterative solver)
!    type(cvector), intent(inout)	:: E0
!    ! Output boundary conditions
!    type(cboundary), intent(inout)	:: BC
!
!    period = (2*PI)/omega ! period is seconds
!
!    if (BC_AND_E0_FROM_FILE) then
!       ! we are going to make a huge assumption here: nPol == 1 always for this case
!       !  and of course transmitters are in same order always
!       E0 = E0_from_file(iTx)
!       call getBC(E0,BC)
!       !   do we now need to set boundary edges of E0 == 0?
!    else
!       if (BC%read_E_from_file) then
!
!          call BC_x0_WS(imode,period,mGrid,sigma_C,E0,BC)
!
!          ! The BC are already computed from a larger grid for all transmitters and modes and stored in BC_from_file.
!          ! Overwrite BC with BC_from_file.
!          ! Note: Right now we are using the same period layout for both grid.
!          ! This why, it is enough to know the period and mode index to pick up the BC from BC_from_file vector.
!          BC = BC_from_file((iTx*2)-(2-imode))
!       else
!         ! Compute the BC using Weerachai 2D approach
!          call BC_x0_WS(imode,period,mGrid,sigma_C,E0,BC)
!       end if
!   end if
!
!
!    ! Cell conductivity array is no longer needed
!    ! NOT TRUE: needed for imode=2
!    ! call deall_rscalar(sigma_C)
!
!  end subroutine SetBound

!*****************************************************************************
   subroutine CurlCurlSetUp()
   !   using existing curl operator, create sparse matrix CC
   !   Note: this is the symmetric form, multiplied by edge volume
   !   elements
      type(spMatCSR_Real)       ::   Temp
      type(spMatCSR_Real)       ::   CC
      type(spMatCSR_Real)       ::   Ttrans
      integer          :: m,n,nz
      real(kind=prec),allocatable,dimension(:)  :: Dtemp
      integer   :: fid
      character*80  ::  cfile

      m = T%nRow
      n = T%nCol
      nz = T%row(T%nRow+1)-1
      allocate(Dtemp(m))
      call create_spMatCSR(m,n,nz,Temp)
      call create_spMatCSR(n,m,nz,Ttrans)
      call create_spMatCSR(m,n,nz,CC)
      call RMATxDIAG(T,EdgeL,Temp)
      Dtemp = DualEdgeL/FaceA
      call DIAGxRMAT(Dtemp,Temp,CC)
      call RMATtrans(T,Ttrans)
      call RMATxRMAT(Ttrans,CC,Temp)
      call DIAGxRMAT(EdgeL,Temp,CC)
      call subMatrix_Real(CC,EDGEi,EDGEi,CCii)
      call subMatrix_Real(CC,EDGEi,EDGEb,CCib)
      call deall_spMatCSR(Temp)
      call deall_spMatCSR(Ttrans)
      call deall_spMatCSR(CC)
      deallocate(Dtemp)
      return
   end subroutine
!*****************************************************************************
   subroutine CurlCurlCleanup()
   !   probably don't need a separate routine for this
      call deall_spMatCSR(CCii)
      call deall_spMatCSR(CCib)
   end subroutine
!*****************************************************************************
   subroutine Mult_Aii(x,adjt,y)
   !   implement the sparse matrix multiply for curl-curl operator
   !   for interior elements
   !   assume output y is already allocated
   complex(kind=prec),intent(in), dimension(:)    :: x
   complex(kind=prec),intent(inout), dimension(:)    :: y
   logical adjt

      call RMATxCVEC(CCii,x,y)
      if(adjt) then
         y = y-cmplx(0,1,8)*ISIGN*VomegaMuSig*x
      else
         y = y+cmplx(0,1,8)*ISIGN*VomegaMuSig*x
      endif

   end subroutine
!*****************************************************************************
   subroutine Mult_Aib(x,adjt,y)
   !   implement the sparse matrix multiply for curl-curl operator
   !   for interior/boundary elements
   !   assume output y is already allocated
   complex(kind=prec),intent(in), dimension(:)       :: x
   logical,intent(in)                                :: adjt
   complex(kind=prec),intent(inout), dimension(:)    :: y
   complex(kind=prec),allocatable, dimension(:)      :: yb
   type(spMatCSR_Real)                               :: CCibt


      if(adjt) then
          call RMATtrans(CCib,CCibt)
          allocate(yb(size(EDGEb)))
          call RMATxCVEC(CCibt,x,yb)
          y(EDGEb)=yb;
          call deall_spMATcsr(CCibt)
          deallocate(yb)
      else
          call RMATxCVEC(CCib,x,y)
      endif
   end subroutine
!*****************************************************************************
   subroutine PC_setupGS()
   !   implement the sparse matrix multiply for curl-curl operator
   !    call after setting CC and VomegaMuSigma for current problem
   !   Initially try the GS preconditioner implemented in ModeMM
   !     by Maxim: L = tril(A), U = diag(A)\triu(A)
   !    will need to do something different to make parallel.

      type(spMatCSR_Real)       ::   Temp
      type(spMatCSR_Cmplx)       ::   TempC
      complex(kind=prec),allocatable,dimension(:)    ::   D
      integer                  :: i
      call upperTri_Real(CCii,Temp)
      call CSR_R2Cdiag(Temp,VomegaMuSig,TempC)
      call diag_Cmplx(TempC,D)
      do i = 1,size(D)
         D(i) = 1.0_dp/D(i)
      enddo
      call DIAGxCMAT(D,TempC,U)
      ! well the flag is only set to Temp... manually set it here
      U%upper=.TRUE.
      deallocate(D)
      call lowerTri_Real(CCii,Temp)
      call CSR_R2Cdiag(Temp,VomegaMuSig,L)
      ! well the flag is only set to Temp... manually set it here
      L%lower=.TRUE.
      call deall_spMatCSR(TempC)
      call deall_spMatCSR(Temp)
      call CMATtrans(L,LH)
      call CMATtrans(U,UH)
      return
   end subroutine
!*****************************************************************************
   subroutine PC_setupDILU()
   !   implement the sparse matrix multiply for curl-curl operator
   !    call after setting CC and VomegaMuSigma for current problem
   !   Initially try the GS preconditioner implemented in ModeMM
   !     by Maxim: L = tril(A), U = diag(A)\triu(A)
   !    will need to do something different to make parallel.

      type(spMatCSR_Cmplx)       ::   TempC
      complex(kind=prec),allocatable,dimension(:)    ::   D
      integer                  :: i

      call CSR_R2Cdiag(CCii,VomegaMuSig,TempC)
      call Dilu_Cmplx(TempC,L,U)
      call CMATtrans(L,LH)
      call CMATtrans(U,UH)
      call deall_spMatCSR(TempC)
      return
   end subroutine
!*****************************************************************************
   subroutine PC_setup()
   !   block DILU preconditioner for CC operator, should be
   !    comparable to what is implemented for matrix-free verson

       implicit none
       integer, allocatable, dimension(:)               :: ix,iy,iz
       real(kind=prec), allocatable, dimension(:)       :: d
       integer                                          :: nx,ny,na,nz
       integer                                          :: nEdge,nEdgeT,n,j
       type(spMatCSR_real)                ::  CCxx
       type(spMatCSR_Cmplx)               ::  Axx
       type(spMatCSR_Cmplx),pointer       ::  Lblk(:),Ublk(:)

    !   find indicies of x, y, z elements

      !   this generates indicies (in list of interior edges)
      !     for x, y, z edges
      nEdgeT = 0
      call setLimits(XEDGE,mGrid,nx,ny,nz)
      nEdge = nx*(ny-2)*(nz-2)
      allocate(ix(nEdge))
      ix = (/ (j,j=nEdgeT+1,nEdgeT+nEdge) /)

      nEdgeT = nEdgeT+nEdge
      call setLimits(YEDGE,mGrid,nx,ny,nz)
      nEdge = (nx-2)*ny*(nz-2)
      allocate(iy(nEdge))
      iy = (/ (j,j=nEdgeT+1,nEdgeT+nEdge) /)

      nEdgeT = nEdgeT+nEdge
      call setLimits(ZEDGE,mGrid,nx,ny,nz)
      nEdge = (nx-2)*(ny-2)*nz
      allocate(iz(nEdge))
      iz = (/ (j,j=nEdgeT+1,nEdgeT+nEdge) /)

    !    construct submatrices for x, y, z components
      allocate(Lblk(3))
      allocate(Ublk(3))
      call SubMatrix_Real(CCii,ix,ix,CCxx)
      n = size(ix)
      allocate(d(n))
      d = VomegaMuSig(ix)
      call CSR_R2Cdiag(CCxx,d,Axx)
      call dilu_Cmplx(Axx,Lblk(1),Ublk(1))
      deallocate(d)

      call SubMatrix_Real(CCii,iy,iy,CCxx)
      n = size(iy)
      allocate(d(n))
      d = VomegaMuSig(iy)
      call CSR_R2Cdiag(CCxx,d,Axx)
      call dilu_Cmplx(Axx,Lblk(2),Ublk(2))
      deallocate(d)

      call SubMatrix_Real(CCii,iz,iz,CCxx)
      n = size(iz)
      allocate(d(n))
      d = VomegaMuSig(iz)
      call CSR_R2Cdiag(CCxx,d,Axx)
      call dilu_Cmplx(Axx,Lblk(3),Ublk(3))
      deallocate(d)

   !  could merge into a single LT and UT matrix, or solve systems
   !      individually
      call BlkDiag_Cmplx(Lblk,L)
      call BlkDiag_Cmplx(Ublk,U)
      call CMATtrans(L,LH)
      call CMATtrans(U,UH)
      deallocate(ix)
      deallocate(iy)
      deallocate(iz)
      call deall_spMatCSR(CCxx)
      call deall_spMatCSR(Axx)
      do j=1,3
          call deall_spMatCSR(Lblk(j))
          call deall_spMatCSR(Ublk(j))
      enddo
      deallocate(Lblk)
      deallocate(Ublk)
      return
   end subroutine
!*****************************************************************************
   subroutine deall_PC()
   !   implement the sparse matrix multiply for curl-curl operator
      call deall_spMatCSR(L)
      call deall_spMatCSR(U)
      call deall_spMatCSR(LH)
      call deall_spMatCSR(UH)
      return
   end subroutine
!*****************************************************************************
   subroutine PC_Lsolve(x,adjt,y)
   !   implement the sparse matrix solve for curl-curl operator
      complex(kind=prec), intent(in), dimension(:) :: x
      logical, intent(in)             :: adjt
      complex(kind=prec), intent(inout), dimension(:) :: y
      if(adjt) then
         call UTsolve_Cmplx(LH,x,y)
      else
         call LTsolve_Cmplx(L,x,y)
      endif
   end subroutine
!*****************************************************************************
   subroutine PC_Usolve(x,adjt,y)
   !   implement the sparse matrix solve for curl-curl operator
      complex(kind=prec), intent(in), dimension(:) :: x
      logical, intent(in)             :: adjt
      complex(kind=prec), intent(inout), dimension(:) :: y
      if(adjt) then
         call LTsolve_Cmplx(UH,x,y)
      else
         call UTsolve_Cmplx(U,x,y)
      endif
   end subroutine
!*****************************************************************************
!    Divergence correction routines
!*****************************************************************************
   subroutine DivCorInit()
   !  These parts of the setup can be done once and for all
   !    without knowing conductivity
      implicit none
      type(spMatCSR_Real)                         ::   Temp,Temp2
      real(kind=prec), allocatable, dimension(:)  ::   d
      integer, allocatable, dimension(:)          ::  allNodes
      integer                                     :: i,m

      !   set indicies for interior and boundary nodes
      call boundaryIndex(CORNER,mGrid,NODEb,NODEi)

      !   (1) first construct VDiv operator
      !   transpose of topology
      call RMATtrans(G,Temp)
      !  pre-multiply by dual face area
      call RMATxDIAG(Temp,DualFaceA,Temp2)
      !   select out interior nodes, edges
      call subMatrix_Real(Temp2,NODEi,EDGEi,VDiv)
      call deall_spMatCSR(Temp)
      call deall_spMatCSR(Temp2)

      !   (2) next turn G into actual gradient (not just topology,
      !       all nodes-> interior edges (not clear this is what we want!)
      allocate(d(G%nRow))
      do i=1,G%nRow
         d(i) = 1./EdgeL(i)
      enddo
      allocate(allNodes(G%nCol))
      do i=1,G%nCol
         allNodes(i) = i
      enddo
      call DIAGxRMAT(d,G,Temp)
      call subMatrix_Real(Temp,EDGEi,allNodes,G)
      call deall_spMatCSR(Temp)
      deallocate(allNodes)
      deallocate(d)
      !   this is a temporary array needed with this implementation of
      !    of the preconditioner --- allocate once, and reuse, clean up
      !     at end
      !
      !   seems wrong -
      !   this should be of size NODEi to be used with VDsG (size(NODEi)
      !   by size(NODEi))
      m = size(NODEi)
      allocate(tempPhi(m))
      return
   end subroutine
!*****************************************************************************
   subroutine DivCorSetup()
   !  to complete setup conductivity is required
   !   DivCorInit has to be called before this routine

      type(spMatCSR_Real)       ::   Temp
      real(kind=prec), allocatable, dimension(:)   :: d
      integer, allocatable, dimension(:)   ::  allNodes
      integer    :: n,i

      ! Construct VDs .. multiply VDiv by Conductivity on edges; can
      !    use VomegaMuSig
      n  = VDiv%nCol
      allocate(d(n))
      d = VomegaMuSig/(Mu_0*omega)
      d = d/Vedge(EDGEi)
      call RMATxDIAG(VDiv,d,VDs)
      ! Construct VDsG   ...  symmetric operator for divergence correction
      !      solver
      allocate(allNodes(G%nRow))
      do i=1,G%nRow
         allNodes(i) = i
      enddo
      call subMatrix_Real(G,allNodes,NODEi,Temp)
      call RMATxRMAT(VDs,Temp,VDsG)
      ! Setup preconditioner
      call Dilu_Real(VDsG,VDsG_L,VDsG_U)
      !call CholInc_real(VDsG,VDsG_L)
      !call RMATtrans(VDsG_L,VDsG_U)
      call deall_spMatCSR(Temp)
      deallocate(d)
      deallocate(allNodes)
      return
   end subroutine
!*****************************************************************************
   subroutine deall_DivCor()
   !   implement the sparse matrix multiply for curl-curl operator
      call deall_spMatCSR(VDiv)
      call deall_spMatCSR(VDs)
      call deall_spMatCSR(VDsG)
      call deall_spMatCSR(VDsG_L)
      call deall_spMatCSR(VDsG_U)
      deallocate(tempPhi)
   end subroutine
!*****************************************************************************
   subroutine DivCgrad(inPhi,outPhi)
   !   implement the sparse matrix multiply for divergence correction op
      complex(kind=prec),intent(in)    :: inPhi(:)
      complex(kind=prec),intent(inout)    :: outPhi(:)

      call RMATxCVEC(VDsG, inPhi, outPhi)
      return

   end subroutine
!*****************************************************************************
   subroutine DivCgradILU(inPhi,outPhi)
   !   implement preconditioner for divergence correction op
      implicit none
      complex(kind=prec),intent(in)       :: inPhi(:)
      complex(kind=prec),intent(inout)    :: outPhi(:)

      call LTsolve_Real(VDsG_L,inPhi,tempPhi)
      call UTsolve_Real(VDsG_U,tempPhi,outPhi)
      return
   end subroutine
!*****************************************************************************
   subroutine DivC(inE,outPhi)
   !   implement multiplication by DivC
      implicit none
      complex(kind=prec),intent(in)    :: inE(:)
      complex(kind=prec),intent(inout)    :: outPhi(:)
        call RMATxCVEC(VDs,inE,outPhi)
      return
   end subroutine
!*****************************************************************************
   subroutine Grad(inPhi,outE)
   !   implement multiplication by E=grad(phi) from all nodes -> inner edges
      implicit none
      complex(kind=prec),intent(in)    :: inPhi(:)
      complex(kind=prec),intent(inout)    :: outE(:)
        call RMATxCVEC(G,inPhi,outE)
      return
   end subroutine

!*****************************************************************************
   subroutine Div(inE,outPhi)
   !   implement multiplication by phi=div(E) from inner edges -> all nodes
      implicit none
      complex(kind=prec),intent(in)       :: inE(:)
      complex(kind=prec),intent(inout)    :: outPhi(:)
      type(spMatCSR_Real)                 :: D
      call RMATtrans(G,D)
      call RMATxCVEC(D,inE,outPhi)
      call deall_spMatCSR(D)
      return
   end subroutine
end module modelOperator3D
