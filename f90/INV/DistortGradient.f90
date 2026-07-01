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

  public :: compute_distortion_gradient

Contains

  subroutine compute_distortion_gradient(d_obs, eAll, sigma, distC, nu, Ndata, gradC)
    type(dataVectorMTX_t), intent(in) :: d_obs
    type(solnVectorMTX_t), intent(in) :: eAll
    type(modelParam_t), intent(in) :: sigma
    type(distortionParam_t), intent(in) :: distC
    real(kind=prec), intent(in) :: nu, Ndata
    type(distortionParam_t), intent(inout) :: gradC

    integer :: j, i, iSite, iTx, iDt, nSites, nComp, nFunc, iFunc, iComp, iDistSite, iRx
    real(kind=prec) :: var
    complex(kind=prec) :: Z_jn(2,2), A_bar(2,2), A_jn(2,2), W(2,2)
    real(kind=prec) :: Resp_undist(8)
    logical :: exists

    call zero_distortionParam(gradC)

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

             ! Compute A_jn = C_j * Z_jn - D_obs
             A_jn(1,1) = distC%C(1,1,iDistSite)*Z_jn(1,1) + distC%C(1,2,iDistSite)*Z_jn(2,1)
             A_jn(1,2) = distC%C(1,1,iDistSite)*Z_jn(1,2) + distC%C(1,2,iDistSite)*Z_jn(2,2)
             A_jn(2,1) = distC%C(2,1,iDistSite)*Z_jn(1,1) + distC%C(2,2,iDistSite)*Z_jn(2,1)
             A_jn(2,2) = distC%C(2,1,iDistSite)*Z_jn(1,2) + distC%C(2,2,iDistSite)*Z_jn(2,2)

             iComp = 1
             do iFunc = 1, nFunc
                exists = d_obs%d(j)%data(i)%exist(iComp, iSite)
                if (exists) then
                   A_jn((iFunc+1)/2, mod(iFunc+1,2)+1) = &
                        A_jn((iFunc+1)/2, mod(iFunc+1,2)+1) - &
                        cmplx(d_obs%d(j)%data(i)%value(iComp,iSite), &
                        d_obs%d(j)%data(i)%value(iComp+1,iSite), prec)
                end if
                iComp = iComp + 2
              end do

             ! Normalize by variance: N_jn = A_jn / σ²_jn
             if (d_obs%d(j)%data(i)%errorBar) then
                var = d_obs%d(j)%data(i)%error(1,iSite)**2
                if (var > R_TINY) then
                   A_jn = A_jn / var
                end if
             end if

             ! Ā = conj(A_normalized)
             A_bar = conjg(A_jn)

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

  end subroutine compute_distortion_gradient

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

end module DistortGradient
