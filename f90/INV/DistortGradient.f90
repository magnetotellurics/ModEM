module DistortGradient

  use DistortionParam
  use DataSpace
  use dataTypes
  use dataFunc, only: dataResp
  use SolnSpace, only: solnVector_t, solnVectorMTX_t
  use ModelSpace, only: modelParam_t
  use fields_orientation, only: orient_t
  use math_constants
  use utilities

  implicit none

  public :: compDistGrad, phi_C_only

Contains

  ! compute the gradient of the data-misfit functional for each 
  ! distortion parameters C_j
  subroutine compDistGrad(d_obs, eAll, sigma, distC, nu, Ndata, gradC)
    type(dataVectorMTX_t), intent(in) :: d_obs
    type(solnVectorMTX_t), intent(in) :: eAll
    type(modelParam_t), intent(in) :: sigma
    type(distortionParam_t), intent(in) :: distC
    real(kind=prec), intent(in) :: nu, Ndata
    type(distortionParam_t), intent(inout) :: gradC

    integer :: j, i, iSite, iTx, iDt, nSites, nComp, nFunc, iFunc, iDistSite, iRx
    integer :: iCompRe, iCompIm, iRow, iCol
    real(kind=prec) :: rRe, rIm, err2Re, err2Im
    complex(kind=prec) :: Z_jn(2,2), A_bar(2,2), A_jn(2,2), A_w(2,2), W(2,2)
    real(kind=prec) :: Resp_undist(8)
    logical :: existsRe, existsIm

    call zero_distortionParam(gradC)

    ! loop over transmitters, data types, and sites to accumulate the gradient
    do j = 1, d_obs%nTx
       iTx = d_obs%d(j)%tx
       do i = 1, d_obs%d(j)%nDt
          iDt = d_obs%d(j)%data(i)%dataType
          if (iDt /= Full_Impedance_Dist) cycle
          nSites = d_obs%d(j)%data(i)%nSite
          nComp = d_obs%d(j)%data(i)%nComp
          nFunc = nComp / 2
          do iSite = 1, nSites
             iRx = d_obs%d(j)%data(i)%rx(iSite)
             iDistSite = find_dist_site(distC, iRx)
             if (iDistSite < 1) cycle

             ! Get undistorted Z_jn from the forward EM solution
             call dataResp(eAll%solns(j), sigma, Full_Impedance, &
                  iRx, Resp_undist, &
                  d_obs%d(j)%data(i)%orient(iSite))
             Z_jn(1,1) = cmplx(Resp_undist(1), Resp_undist(2), prec)
             Z_jn(1,2) = cmplx(Resp_undist(3), Resp_undist(4), prec)
             Z_jn(2,1) = cmplx(Resp_undist(5), Resp_undist(6), prec)
             Z_jn(2,2) = cmplx(Resp_undist(7), Resp_undist(8), prec)

             ! Raw residual seed: A_jn = C_j * Z_jn
             A_jn(1,1) = distC%C(1,1,iDistSite)*Z_jn(1,1) + distC%C(1,2,iDistSite)*Z_jn(2,1)
             A_jn(1,2) = distC%C(1,1,iDistSite)*Z_jn(1,2) + distC%C(1,2,iDistSite)*Z_jn(2,2)
             A_jn(2,1) = distC%C(2,1,iDistSite)*Z_jn(1,1) + distC%C(2,2,iDistSite)*Z_jn(2,1)
             A_jn(2,2) = distC%C(2,1,iDistSite)*Z_jn(1,2) + distC%C(2,2,iDistSite)*Z_jn(2,2)

             ! Build weighted residual matrix A_w with component-wise errors:
             !   A_w(row,col) = (Re residual / errRe^2) + i (Im residual / errIm^2)
             ! Missing components (exist=.false.) contribute zero.
             A_w = C_ZERO
             do iFunc = 1, nFunc
                iRow = (iFunc + 1) / 2
                iCol = mod(iFunc + 1, 2) + 1
                iCompRe = 2*iFunc - 1
                iCompIm = 2*iFunc

                existsRe = d_obs%d(j)%data(i)%exist(iCompRe, iSite)
                existsIm = d_obs%d(j)%data(i)%exist(iCompIm, iSite)

                rRe = R_ZERO
                rIm = R_ZERO

                if (existsRe) then 
                    rRe = real(A_jn(iRow,iCol)) - d_obs%d(j)%data(i)%value(iCompRe, iSite)
                endif
                if (existsIm) then
                    rIm = aimag(A_jn(iRow,iCol)) - d_obs%d(j)%data(i)%value(iCompIm, iSite)
                endif

                if (d_obs%d(j)%data(i)%errorBar) then
                   if (existsRe) then
                      err2Re = d_obs%d(j)%data(i)%error(iCompRe, iSite)**2
                      if (err2Re <= R_TINY) then 
                          call errStop('data error too small in compDistGrad (Re)')
                      endif
                      rRe = rRe / err2Re
                   end if
                   if (existsIm) then
                      err2Im = d_obs%d(j)%data(i)%error(iCompIm, iSite)**2
                      if (err2Im <= R_TINY) then 
                          call errStop('data error too small in compDistGrad (Im)')
                      endif
                      rIm = rIm / err2Im
                   end if
                end if

                A_w(iRow, iCol) = cmplx(rRe, rIm, prec)
             end do

             ! Ā = conj(weighted residual)
             A_bar = conjg(A_w)

             ! ∂φ_d/∂C_j = Σ_n Ā_jn · Z_jn^T
             W(1,1) = A_bar(1,1)*Z_jn(1,1) + A_bar(1,2)*Z_jn(1,2)
             W(1,2) = A_bar(1,1)*Z_jn(2,1) + A_bar(1,2)*Z_jn(2,2)
             W(2,1) = A_bar(2,1)*Z_jn(1,1) + A_bar(2,2)*Z_jn(1,2)
             W(2,2) = A_bar(2,1)*Z_jn(2,1) + A_bar(2,2)*Z_jn(2,2)

             gradC%C(1,1,iDistSite) = gradC%C(1,1,iDistSite) + real(W(1,1))
             gradC%C(1,2,iDistSite) = gradC%C(1,2,iDistSite) + real(W(1,2))
             gradC%C(2,1,iDistSite) = gradC%C(2,1,iDistSite) + real(W(2,1))
             gradC%C(2,2,iDistSite) = gradC%C(2,2,iDistSite) + real(W(2,2))
           end do
        end do
     end do

     ! Scale by (2/Ndata) for consistency with the total penalty functional
     do iSite = 1, distC%nSites
        gradC%C(:,:,iSite) = (TWO / Ndata) * gradC%C(:,:,iSite)
     end do

    ! Add regularization gradient: ν · (C_j - I)
    do iSite = 1, distC%nSites
       gradC%C(1,1,iSite) = gradC%C(1,1,iSite) + nu * (distC%C(1,1,iSite) - ONE)
       gradC%C(1,2,iSite) = gradC%C(1,2,iSite) + nu * distC%C(1,2,iSite)
       gradC%C(2,1,iSite) = gradC%C(2,1,iSite) + nu * distC%C(2,1,iSite)
       gradC%C(2,2,iSite) = gradC%C(2,2,iSite) + nu * (distC%C(2,2,iSite) - ONE)
    end do

  end subroutine compDistGrad

  function phi_C_only(d_obs, eAll, sigma, distC, nu, Ndata, lambda, &
                      mNorm_over_Nmodel) result(phi)
  ! phi_C_only evaluates the total objective regarding C only, reusing locked 
  ! eAll solutions from the sigma update step.
  !
  ! As eAll and sigma are already calculated and fixed (only C changed), no full 
  ! forward calculation solve is required in this function
  ! The data-misfit term is computed by re-evaluating the distorted predicted
  ! data as D_hat_jn = C_j * Z_jn (where Z_jn comes from calling dataResp on
  ! eAll) and forming the weighted sum-of-squares.
  ! The distortion regularisation term psi(C) = 0.5 * sum_j ||C_j - I||^2_F
  ! is also included.

    type(dataVectorMTX_t), intent(in) :: d_obs ! observed data
    type(solnVectorMTX_t), intent(in) :: eAll ! locked EM solutions
    type(modelParam_t),    intent(in) :: sigma ! locked conductivity model
    type(distortionParam_t), intent(in) :: distC ! distortion matrices
    real(kind=prec), intent(in) :: nu, Ndata, lambda ! regularisation weights
    real(kind=prec), intent(in) :: mNorm_over_Nmodel ! pre-computed 
                                  ! lambda*||mHat||^2/Nmodel term (scalar)
    real(kind=prec) :: phi ! F_C = SS/Ndata + lambda*mNorm + nu*psi(C)

    integer :: j, i, iSite, iDt, nSites, nComp, nFunc, iFunc, iDistSite, iRx
    integer :: iCompRe, iCompIm, iRow, iCol
    real(kind=prec) :: SS, psiC, rRe, rIm, err2Re, err2Im
    real(kind=prec) :: Resp_undist(8)
    complex(kind=prec) :: Z_jn(2,2), A_jn(2,2)
    logical :: existsRe, existsIm

    SS   = R_ZERO
    psiC = R_ZERO

    ! data-misfit contribution
    do j = 1, d_obs%nTx
       do i = 1, d_obs%d(j)%nDt
          iDt = d_obs%d(j)%data(i)%dataType
          if (iDt /= Full_Impedance_Dist) cycle
          nSites = d_obs%d(j)%data(i)%nSite
          nComp  = d_obs%d(j)%data(i)%nComp
          nFunc  = nComp / 2
          do iSite = 1, nSites
             iRx = d_obs%d(j)%data(i)%rx(iSite)
             iDistSite = find_dist_site(distC, iRx)
             if (iDistSite < 1) cycle

             ! undistorted Z from locked eAll
             call dataResp(eAll%solns(j), sigma, Full_Impedance, &
                  iRx, Resp_undist, d_obs%d(j)%data(i)%orient(iSite))
             Z_jn(1,1) = cmplx(Resp_undist(1), Resp_undist(2), prec)
             Z_jn(1,2) = cmplx(Resp_undist(3), Resp_undist(4), prec)
             Z_jn(2,1) = cmplx(Resp_undist(5), Resp_undist(6), prec)
             Z_jn(2,2) = cmplx(Resp_undist(7), Resp_undist(8), prec)

             ! Predicted distorted response A_jn = C_j * Z_jn
             A_jn(1,1) = distC%C(1,1,iDistSite)*Z_jn(1,1) + distC%C(1,2,iDistSite)*Z_jn(2,1)
             A_jn(1,2) = distC%C(1,1,iDistSite)*Z_jn(1,2) + distC%C(1,2,iDistSite)*Z_jn(2,2)
             A_jn(2,1) = distC%C(2,1,iDistSite)*Z_jn(1,1) + distC%C(2,2,iDistSite)*Z_jn(2,1)
             A_jn(2,2) = distC%C(2,1,iDistSite)*Z_jn(1,2) + distC%C(2,2,iDistSite)*Z_jn(2,2)

             ! Component-wise residual and weighting, consistent with CdInvMult:
             !   SS = sum_k (res_k^2 / err_k^2), counting only existing components.
             do iFunc = 1, nFunc
                iRow = (iFunc + 1) / 2
                iCol = mod(iFunc + 1, 2) + 1
                iCompRe = 2*iFunc - 1
                iCompIm = 2*iFunc

                existsRe = d_obs%d(j)%data(i)%exist(iCompRe, iSite)
                existsIm = d_obs%d(j)%data(i)%exist(iCompIm, iSite)

                if (existsRe) then
                   rRe = real(A_jn(iRow,iCol)) - d_obs%d(j)%data(i)%value(iCompRe, iSite)
                   if (d_obs%d(j)%data(i)%errorBar) then
                      err2Re = d_obs%d(j)%data(i)%error(iCompRe, iSite)**2
                      if (err2Re <= R_TINY) then
                          call errStop('data error too small in phi_C_only (Re)')
                      endif
                      SS = SS + (rRe*rRe) / err2Re
                   else
                      SS = SS + rRe*rRe
                   end if
                end if

                if (existsIm) then
                   rIm = aimag(A_jn(iRow,iCol)) - d_obs%d(j)%data(i)%value(iCompIm, iSite)
                   if (d_obs%d(j)%data(i)%errorBar) then
                      err2Im = d_obs%d(j)%data(i)%error(iCompIm, iSite)**2
                      if (err2Im <= R_TINY) then
                          call errStop('data error too small in phi_C_only (Im)')
                      endif
                      SS = SS + (rIm*rIm) / err2Im
                   else
                      SS = SS + rIm*rIm
                   end if
                end if
             end do
          end do
       end do
    end do

    ! distortion regularisation
    do iSite = 1, distC%nSites
       psiC = psiC &
            + (distC%C(1,1,iSite) - ONE)**2 &
            + distC%C(1,2,iSite)**2 &
            + distC%C(2,1,iSite)**2 &
            + (distC%C(2,2,iSite) - ONE)**2
    end do
    psiC = psiC * 0.5_prec

    ! normalized total objective function 
    phi = SS / Ndata + lambda * mNorm_over_Nmodel + nu * psiC

  end function phi_C_only 

  function find_dist_site(dist, rxIdx) result(idx)
  ! a simple untility to find the index of a distortion site given 
  ! the receiver index
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

end module DistortGradient
