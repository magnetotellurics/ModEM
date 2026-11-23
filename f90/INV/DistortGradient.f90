module DistortGradient
! Distortion matrix estimation for the joint conductivity/distortion inversion
! (NLCG_DIST).  
! For a fixed conductivity model, the distorted impedances C_j Z_jn are 
! linear in the distortion matrices C_j, so the optimal C for the current fields 
! can be obtained exactly from trivial per-site least-squares problems (with no 
! additional EM solves).
!
! The distortion penalty is
!   psi(C) = 1/2 sum_j w_j ||C_j - C0_j||_F^2
! with C0 the reference distortion (identity by default) and w_j the mean
! diagonal of the data-misfit Hessian with respect to C_j, computed once from
! the observed impedances and their errors.  With these weights the damping
! factor nu is dimensionless: the estimated C_j moves roughly a fraction
! 1/(1+nu) of the way from C0_j towards its data-only fit.
!
! Two distortion types are supported:
!   full : general real 2x2 C_j (4 parameters per site)
!   diag : C_j = diag(gx_j, gy_j), one gain per electric field component
!          (2 parameters per site; the static-shift/splitting part of C)

  use DistortionParam
  use DataSpace
  use dataTypes
  use dataFunc, only: dataResp
  use SolnSpace, only: solnVector_t, solnVectorMTX_t
  use ModelSpace, only: modelParam_t
  use math_constants
  use utilities

  implicit none

  private

  public :: DIST_FULL, DIST_DIAG
  public :: init_dist_penalty, deall_dist_penalty, psi_dist, solveDistExact

  integer, parameter :: DIST_FULL = 1
  integer, parameter :: DIST_DIAG = 2

  ! penalty weights w_j, reference distortion C0 and distortion type
  real(kind=prec), allocatable, save :: distW(:)
  real(kind=prec), allocatable, save :: distC0(:,:,:)
  integer, save                      :: distType = DIST_FULL

Contains

  ! ---------------------------------------------------------------------------
  ! component indices (Re part) of Z(row,col) in the 8-component impedance
  ! vector; the Im part is the next index
  ! ---------------------------------------------------------------------------
  pure integer function iZre(iRow, iCol)
    integer, intent(in) :: iRow, iCol
    iZre = 4*(iRow-1) + 2*(iCol-1) + 1
  end function iZre

  function find_dist_site(dist, rxIdx) result(idx)
    type(distortionParam_t), intent(in) :: dist
    integer, intent(in) :: rxIdx
    integer :: idx, k
    idx = 0
    if (.not. allocated(dist%siteIndex)) return
    do k = 1, dist%nSites
       if (dist%siteIndex(k) == rxIdx) then
          idx = k
          return
       end if
    end do
  end function find_dist_site

  ! ---------------------------------------------------------------------------
  ! init_dist_penalty sets the distortion type, the reference distortion C0
  ! (copied from Cref) and the penalty weights w_j.  w_j is the mean diagonal
  ! of the data-misfit Hessian with respect to the entries of C_j, evaluated
  ! with the observed impedances (2/Ndata * sum |Z|^2/err^2).  It only depends
  ! on the observed data, so it is fixed during the inversion.
  ! ---------------------------------------------------------------------------
  subroutine init_dist_penalty(d_obs, Cref, dType)
    type(dataVectorMTX_t), intent(in) :: d_obs
    type(distortionParam_t), intent(in) :: Cref
    integer, intent(in) :: dType

    integer :: j, i, iSite, iDistSite, iRow, iCol, k, nAvail, ic, kc
    real(kind=prec) :: wRe, wIm, wbar, Ndata

    call deall_dist_penalty()
    distType = dType
    allocate(distW(Cref%nSites), distC0(2,2,Cref%nSites))
    distC0 = Cref%C
    ! a diagonal model can only use the diagonal of the reference
    if (distType == DIST_DIAG) then
       distC0(1,2,:) = R_ZERO
       distC0(2,1,:) = R_ZERO
    end if
    distW = R_ZERO
    Ndata = real(countData(d_obs), kind=prec)

    do j = 1, d_obs%nTx
       do i = 1, d_obs%d(j)%nDt
          if (d_obs%d(j)%data(i)%dataType /= Full_Impedance_Dist) cycle
          associate (blk => d_obs%d(j)%data(i))
          do iSite = 1, blk%nSite
             iDistSite = find_dist_site(Cref, blk%rx(iSite))
             if (iDistSite < 1) cycle
             ! d^2/dC(iRow,k)^2 of the misfit of data element (iRow,iCol)
             ! is 2/N |Z(k,iCol)|^2 / err(iRow,iCol)^2
             do iRow = 1, 2
                do iCol = 1, 2
                   ic = iZre(iRow,iCol)
                   if (.not. blk%exist(ic,iSite)) cycle
                   wRe = ONE; wIm = ONE
                   if (blk%errorBar) then
                      wRe = ONE / blk%error(ic,iSite)**2
                      wIm = ONE / blk%error(ic+1,iSite)**2
                   end if
                   do k = 1, 2
                      kc = iZre(k,iCol)
                      if (.not. blk%exist(kc,iSite)) cycle
                      distW(iDistSite) = distW(iDistSite) + (TWO/Ndata) * &
                           (wRe * blk%value(kc,iSite)**2 + wIm * blk%value(kc+1,iSite)**2)
                   end do
                end do
             end do
          end do
          end associate
       end do
    end do

    ! mean curvature over the free parameters of C_j: 4 for full, 2 for
    ! diag (the two diagonal entries).  This keeps nu dimensionless and
    ! comparable between distortion types.
    distW = distW / 4.0_prec
    ! sites without usable data get the mean weight
    nAvail = count(distW > R_TINY)
    if (nAvail > 0) then
       wbar = sum(distW, mask=distW > R_TINY) / nAvail
    else
       wbar = ONE
    end if
    where (distW <= R_TINY) distW = wbar

  end subroutine init_dist_penalty

  subroutine deall_dist_penalty()
    if (allocated(distW)) deallocate(distW)
    if (allocated(distC0)) deallocate(distC0)
    distType = DIST_FULL
  end subroutine deall_dist_penalty

  ! ---------------------------------------------------------------------------
  ! psi_dist = 1/2 sum_j w_j ||C_j - C0_j||_F^2
  ! ---------------------------------------------------------------------------
  function psi_dist(distC) result(psi)
    type(distortionParam_t), intent(in) :: distC
    real(kind=prec) :: psi
    integer :: iSite

    psi = R_ZERO
    if (.not. allocated(distW)) call errStop('psi_dist: distortion penalty not initialised')
    do iSite = 1, distC%nSites
       psi = psi + distW(iSite) * sum((distC%C(:,:,iSite) - distC0(:,:,iSite))**2)
    end do
    psi = 0.5_prec * psi
  end function psi_dist

  ! ---------------------------------------------------------------------------
  ! solveDistExact replaces distC by the minimiser of
  !   SS(C)/Ndata + nu * psi(C)
  ! for the locked EM solution eAll.  Only Full_Impedance_Dist data depend on
  ! C.  Because the result is the exact minimiser, the conductivity gradient
  ! evaluated at this C is the gradient of the reduced objective
  ! min_C F(sigma, C).
  !
  ! full: the problem separates into the two rows of C_j (row i of C only
  !       affects row i of the distorted impedance): a 2x2 system per row.
  ! diag: C_j = diag(gx,gy): one scalar problem per row (gx for row 1,
  !       gy for row 2).
  ! ---------------------------------------------------------------------------
  subroutine solveDistExact(d_obs, eAll, sigma, distC, nu)
    type(dataVectorMTX_t), intent(in) :: d_obs
    type(solnVectorMTX_t), intent(in) :: eAll
    type(modelParam_t), intent(in) :: sigma
    type(distortionParam_t), intent(inout) :: distC
    real(kind=prec), intent(in) :: nu

    real(kind=prec), allocatable :: H(:,:,:,:), b(:,:,:), Hg(:,:), bg(:,:)
    real(kind=prec) :: Resp(8), z(2), dval, w, A(2,2), rhs(2), det, Ndata
    integer :: j, i, iSite, iDistSite, iRow, iCol, k, part, nS, ic

    nS = distC%nSites
    if (nS <= 0) return
    if (.not. allocated(distW)) call errStop('solveDistExact: distortion penalty not initialised')
    Ndata = real(countData(d_obs), kind=prec)

    allocate(H(2,2,2,nS), b(2,2,nS), Hg(2,nS), bg(2,nS))
    H = R_ZERO; b = R_ZERO; Hg = R_ZERO; bg = R_ZERO

    ! normal equations: for data element (iRow,iCol) and part (Re/Im)
    !   d = C(iRow,1) Z(1,iCol) + C(iRow,2) Z(2,iCol) = c . z
    do j = 1, d_obs%nTx
       do i = 1, d_obs%d(j)%nDt
          if (d_obs%d(j)%data(i)%dataType /= Full_Impedance_Dist) cycle
          associate (blk => d_obs%d(j)%data(i))
          do iSite = 1, blk%nSite
             iDistSite = find_dist_site(distC, blk%rx(iSite))
             if (iDistSite < 1) cycle
             call dataResp(eAll%solns(j), sigma, Full_Impedance, &
                  blk%rx(iSite), Resp, blk%orient(iSite))
             do iRow = 1, 2
                do iCol = 1, 2
                   do part = 0, 1
                      ic = iZre(iRow,iCol) + part
                      if (.not. blk%exist(ic,iSite)) cycle
                      z(1) = Resp(iZre(1,iCol) + part)
                      z(2) = Resp(iZre(2,iCol) + part)
                      dval = blk%value(ic,iSite)
                      w = ONE
                      if (blk%errorBar) then
                         if (blk%error(ic,iSite)**2 <= R_TINY) &
                              call errStop('data error too small in solveDistExact')
                         w = ONE / blk%error(ic,iSite)**2
                      end if
                      w = (TWO / Ndata) * w
                      do k = 1, 2
                         H(:,k,iRow,iDistSite) = H(:,k,iRow,iDistSite) + w * z * z(k)
                      end do
                      b(:,iRow,iDistSite) = b(:,iRow,iDistSite) + w * dval * z
                      ! diag: predicted value is g * Z(iRow,iCol),
                      ! accumulated per row
                      Hg(iRow,iDistSite) = Hg(iRow,iDistSite) + w * z(iRow)**2
                      bg(iRow,iDistSite) = bg(iRow,iDistSite) + w * dval * z(iRow)
                   end do
                end do
             end do
          end do
          end associate
       end do
    end do

    do iSite = 1, nS
       select case (distType)
       case (DIST_DIAG)
          ! psi_j = 1/2 w_j ||diag(gx,gy) - C0||^2; one scalar solve per row
          distC%C(:,:,iSite) = R_ZERO
          distC%C(1,1,iSite) = (bg(1,iSite) + nu*distW(iSite)*distC0(1,1,iSite)) / &
                               (Hg(1,iSite) + nu*distW(iSite))
          distC%C(2,2,iSite) = (bg(2,iSite) + nu*distW(iSite)*distC0(2,2,iSite)) / &
                               (Hg(2,iSite) + nu*distW(iSite))
       case default
          do iRow = 1, 2
             A = H(:,:,iRow,iSite)
             A(1,1) = A(1,1) + nu * distW(iSite)
             A(2,2) = A(2,2) + nu * distW(iSite)
             rhs = b(:,iRow,iSite) + nu * distW(iSite) * distC0(iRow,:,iSite)
             det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
             if (abs(det) <= R_TINY) call errStop('singular distortion system in solveDistExact')
             distC%C(iRow,1,iSite) = ( A(2,2)*rhs(1) - A(1,2)*rhs(2)) / det
             distC%C(iRow,2,iSite) = (-A(2,1)*rhs(1) + A(1,1)*rhs(2)) / det
          end do
       end select
    end do

    deallocate(H, b, Hg, bg)
  end subroutine solveDistExact

end module DistortGradient
