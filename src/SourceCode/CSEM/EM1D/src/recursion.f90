!---------------------------------------------------------
! integrand 1 for Ex for long wire source, receiver above source
!---------------------------------------------------------
function i1abvExwire(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  !put everything together (except factor -1/4pi, which is multiplied at the very end)
  res = - kappa * iabvA1TE(kappa)

endfunction i1abvExwire


!---------------------------------------------------------
! integrand 1 for Ex for long wire source, receiver below source
!---------------------------------------------------------
function i1blwExwire(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = - kappa * iblwA1TE(kappa)

endfunction i1blwExwire


!---------------------------------------------------------
! integrand 2 for Ex for long wire source, receiver above source
!---------------------------------------------------------
function i2abvExwire(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  !internal variables
  complex(kind=real64)            :: resTE, resTM  !results from TE and TM parts

  !TE part
  resTE = iabvA1TE(kappa)
  !TM part
  resTM = iabvA1TM(kappa)

  res = resTE - resTM

endfunction i2abvExwire


!---------------------------------------------------------
! integrand 2 for Ex for long wire source, receiver below source
!---------------------------------------------------------
function i2blwExwire(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  !internal variables
  complex(kind=real64)            :: resTE, resTM  !results from TE and TM parts

  !TE part
  resTE = iblwA1TE(kappa)
  !TM part
  resTM = iblwA1TM(kappa)

  res = resTE - resTM

endfunction i2blwExwire


!---------------------------------------------------------
! integrand 1 for Hx for long wire source, rec. above source
!---------------------------------------------------------
function idabvHxwire(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  !internal variables
  complex(kind=real64)            :: resTE, resTM  !results from TE and TM parts

  !TE part
  resTE = iabvD1TE(kappa)
  !TM part
  resTM = iabvD1TM(kappa)

  res = resTE - resTM

endfunction idabvHxwire


!---------------------------------------------------------
! integrand 1 for Hx for long wire source, rec. below source
!---------------------------------------------------------
function idblwHxwire(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  !internal variables
  complex(kind=real64)            :: resTE, resTM  !results from TE and TM parts

  !TE part
  resTE = iblwD1TE(kappa)
  !TM part
  resTM = iblwD1TM(kappa)

  res = resTE - resTM

endfunction idblwHxwire


!---------------------------------------------------------
! integrand A0 TE for receivers above source
!---------------------------------------------------------
function iabvA0TE(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa * iabvA1TE(kappa)

endfunction iabvA0TE


!---------------------------------------------------------
! integrand A0 TM for receivers above source
!---------------------------------------------------------
function iabvA0TM(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa * iabvA1TM(kappa)

endfunction iabvA0TM


!---------------------------------------------------------
! integrand A1 TE for receivers above source
!---------------------------------------------------------
function iabvA1TE(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  !internal variables
  complex(kind=real64)            :: refETE    !temp result
  complex(kind=real64)            :: rt11EH    !temp result
  complex(kind=real64)            :: ra11      !temp result
  complex(kind=real64)            :: rdnsrc11  !temp result
  complex(kind=real64)            :: tu11      !contribution of homogeneous field

  !  total reflection responses for "up-down symmetric" media
  call set_fresnelcoef_TE(kappa,pvert)

  call get_rt_abv(rt11EH,ra11,rdnsrc11,kappa,tu11,set_trans_abvTE)
  refETE = (1. + ra11) * (1. + rdnsrc11) * rt11EH

  !put everything together (except factor -1/4pi, which is multiplied at the very end)
  res = (dmu0 / sqrt(pvert(ilayrec)*pvert(ilaysrc))) * refETE

endfunction iabvA1TE


!---------------------------------------------------------
! integrand A1 TM for receivers above source
!---------------------------------------------------------
function iabvA1TM(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  !internal variables
  complex(kind=real64)            :: refETM    !temp result
  complex(kind=real64)            :: rt22EH    !temp result
  complex(kind=real64)            :: ra22      !temp result
  complex(kind=real64)            :: rdnsrc22  !temp result
  complex(kind=real64)            :: tu22      !contribution of homogeneous field

  !  total reflection responses for "up-down symmetric" media
  call set_fresnelcoef_TM(kappa,pvert)
  call get_rt_abv(rt22EH,ra22,rdnsrc22,kappa,tu22,set_trans_abvTM)
  refETM = (1. + ra22) * (1.+rdnsrc22) * rt22EH

  !put everything together (except factor -1/4pi, which is multiplied at the very end)
  res = (sqrt(pvert(ilayrec)*pvert(ilaysrc) / (epsh(ilayrec)*epsh(ilaysrc)))) * refETM

endfunction iabvA1TM


!---------------------------------------------------------
! integrand A0 TE for receivers below source
!---------------------------------------------------------
function iblwA0TE(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa * iblwA1TE(kappa)

endfunction iblwA0TE


!---------------------------------------------------------
! integrand A0 TM for receivers below source
!---------------------------------------------------------
function iblwA0TM(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa * iblwA1TM(kappa)

endfunction iblwA0TM


!---------------------------------------------------------
! integrand A1 TE for receivers below source
!---------------------------------------------------------
function iblwA1TE(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  !internal variables
  complex(kind=real64)            :: refETE    !temp result
  complex(kind=real64)            :: rt11EH    !temp result
  complex(kind=real64)            :: rb11      !temp result
  complex(kind=real64)            :: rupsrc11  !temp result
  complex(kind=real64)            :: td11      !contribution of homogeneous field

  !  total reflection responses for "up-down symmetric" media
  call set_fresnelcoef_TE(kappa,pvert)
  call get_rt_blw(rt11EH,rb11,rupsrc11,kappa,td11,set_trans_blwTE)
  refETE = (1. + rb11) * (1.+rupsrc11) * rt11EH

  !put everything together (except factor -1/4pi, which is multiplied at the very end)
  res = (dmu0 / sqrt(pvert(ilayrec)*pvert(ilaysrc))) * refETE

endfunction iblwA1TE


!---------------------------------------------------------
! integrand A1 TM for receivers below source
!---------------------------------------------------------
function iblwA1TM(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  !internal variables
  complex(kind=real64)            :: refETM    !temp result
  complex(kind=real64)            :: rt22EH    !temp result
  complex(kind=real64)            :: rb22      !temp result
  complex(kind=real64)            :: rupsrc22  !temp result
  complex(kind=real64)            :: td22      !contribution of homogeneous field

  !  total reflection responses for "up-down symmetric" media
  call set_fresnelcoef_TM(kappa,pvert)
  call get_rt_blw(rt22EH,rb22,rupsrc22,kappa,td22,set_trans_blwTM)
  refETM = (1. + rb22) * (1.+rupsrc22) * rt22EH

  !put everything together (except factor -1/4pi, which is multiplied at the very end)
  res = (sqrt(pvert(ilayrec)*pvert(ilaysrc) / (epsh(ilayrec)*epsh(ilaysrc)))) * refETM

endfunction iblwA1TM


!---------------------------------------------------------
! integrand Dz1 TM for receivers above source
!---------------------------------------------------------
function iabvDz1TM(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa**2 * iabvD1TM(kappa)

endfunction iabvDz1TM


!---------------------------------------------------------
! integrand Dz1 TM for receivers below source
!---------------------------------------------------------
function iblwDz1TM(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa**2 * iblwD1TM(kappa)

endfunction iblwDz1TM


!---------------------------------------------------------
! integrand D0 TE for receivers above source
!---------------------------------------------------------
function iabvD0TE(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa * iabvD1TE(kappa)

endfunction iabvD0TE


!---------------------------------------------------------
! integrand D0 TM for receivers above source
!---------------------------------------------------------
function iabvD0TM(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa * iabvD1TM(kappa)

endfunction iabvD0TM


!---------------------------------------------------------
! integrand D1 TE for receivers above source
!---------------------------------------------------------
function iabvD1TE(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  !internal variables
  complex(kind=real64)            :: refHTE    !temp result
  complex(kind=real64)            :: rt11EH    !temp result
  complex(kind=real64)            :: ra11      !temp result
  complex(kind=real64)            :: rdnsrc11  !temp result
  complex(kind=real64)            :: tu11      !contribution of homogeneous field

  !  total reflection responses for "up-down symmetric" media
  call set_fresnelcoef_TE(kappa,pvert)
  call get_rt_abv(rt11EH,ra11,rdnsrc11,kappa,tu11,set_trans_abvTE)
  refHTE = (1. - ra11) * (1.+rdnsrc11) * rt11EH

  !put everything together (except factor -1/4pi, which is multiplied at the very end)
  res = - (sqrt(pvert(ilayrec) /pvert(ilaysrc))) * refHTE

endfunction iabvD1TE


!---------------------------------------------------------
! integrand D1 TM for receivers above source
!---------------------------------------------------------
function iabvD1TM(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  !internal variables
  complex(kind=real64)            :: refHTM    !temp result
  complex(kind=real64)            :: rt22EH    !temp result
  complex(kind=real64)            :: ra22      !temp result
  complex(kind=real64)            :: rdnsrc22  !temp result
  complex(kind=real64)            :: tu22      !contribution of homogeneous field

  !  total reflection responses for "up-down symmetric" media
  call set_fresnelcoef_TM(kappa,pvert)
  call get_rt_abv(rt22EH,ra22,rdnsrc22,kappa,tu22,set_trans_abvTM)
  refHTM = (1. - ra22) * (1.+rdnsrc22) * rt22EH

  !put everything together (except factor -1/4pi, which is multiplied at the very end)
  res = - (sqrt(epsh(ilayrec)*pvert(ilaysrc) / (epsh(ilaysrc)*pvert(ilayrec)))) * refHTM

endfunction iabvD1TM


!---------------------------------------------------------
! integrand D0 TE for receivers below source
!---------------------------------------------------------
function iblwD0TE(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa *  iblwD1TE(kappa)

endfunction iblwD0TE


!---------------------------------------------------------
! integrand D0 TM for receivers below source
!---------------------------------------------------------
function iblwD0TM(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa * iblwD1TM(kappa)

endfunction iblwD0TM


!---------------------------------------------------------
! integrand D1 TE for receivers below source
!---------------------------------------------------------
function iblwD1TE(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  !internal variables
  complex(kind=real64)            :: refHTE    !temp result
  complex(kind=real64)            :: rt11EH    !temp result
  complex(kind=real64)            :: rb11      !temp result
  complex(kind=real64)            :: rupsrc11  !temp result
  complex(kind=real64)            :: td11      !contribution of homogeneous field

  !  total reflection responses for "up-down symmetric" media
  call set_fresnelcoef_TE(kappa,pvert)
  call get_rt_blw(rt11EH,rb11,rupsrc11,kappa,td11,set_trans_blwTE)
  refHTE = (rb11 - 1.) * (1.+rupsrc11) * rt11EH

  !put everything together
  res = - (sqrt(pvert(ilayrec) / pvert(ilaysrc))) * refHTE

endfunction iblwD1TE


!---------------------------------------------------------
! integrand D1 TM for receivers below source
!---------------------------------------------------------
function iblwD1TM(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  !internal variables
  complex(kind=real64)            :: refHTM    !temp result
  complex(kind=real64)            :: rt22EH    !temp result
  complex(kind=real64)            :: rb22      !temp result
  complex(kind=real64)            :: rupsrc22  !temp result
  complex(kind=real64)            :: td22      !contribution of homogeneous field

  !  total reflection responses for "up-down symmetric" media
  call set_fresnelcoef_TM(kappa,pvert)
  call get_rt_blw(rt22EH,rb22,rupsrc22,kappa,td22,set_trans_blwTM)
  refHTM = (rb22 - 1.) * (1.+rupsrc22) * rt22EH

  !put everything together
  res = - (sqrt(epsh(ilayrec)*pvert(ilaysrc) / (pvert(ilayrec)*epsh(ilaysrc)))) * refHTM

endfunction iblwD1TM


!---------------------------------------------------------
! integrand Az1 TE for receivers above source
!---------------------------------------------------------
function iabvAz1TE(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa**2 * iabvA1TE(kappa)

endfunction iabvAz1TE


!---------------------------------------------------------
! integrand Az1 TE for receivers below source
!---------------------------------------------------------
function iblwAz1TE(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa**2 * iblwA1TE(kappa)

endfunction iblwAz1TE


!---------------------------------------------------------
! integrand B0 TE for receivers above source
!---------------------------------------------------------
function iabvB0TE(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa * iabvB1TE(kappa)

endfunction iabvB0TE


!---------------------------------------------------------
! integrand B0 TM for receivers above source
!---------------------------------------------------------
function iabvB0TM(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa * iabvB1TM(kappa)

endfunction iabvB0TM


!---------------------------------------------------------
! integrand B1 TE for receivers above source
!---------------------------------------------------------
function iabvB1TE(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  !internal variables
  complex(kind=real64)            :: refETE    !temp result
  complex(kind=real64)            :: rt11EH    !temp result
  complex(kind=real64)            :: ra11      !temp result
  complex(kind=real64)            :: rdnsrc11  !temp result
  complex(kind=real64)            :: tu11      !contribution of homogeneous field

  !  total reflection responses for "up-down symmetric" media
  call set_fresnelcoef_TE(kappa,pvert)
  call get_rt_abv(rt11EH,ra11,rdnsrc11,kappa,tu11,set_trans_abvTE)
  refETE = (1. + ra11) * (rdnsrc11 - 1.) * rt11EH

  !put everything together (except factor -1/4pi, which is multiplied at the very end)
  res = sqrt(pvert(ilaysrc) / pvert(ilayrec)) * refETE

endfunction iabvB1TE


!---------------------------------------------------------
! integrand B1 TM for receivers above source
!---------------------------------------------------------
function iabvB1TM(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  !internal variables
  complex(kind=real64)            :: refETM    !temp result
  complex(kind=real64)            :: rt22EH    !temp result
  complex(kind=real64)            :: ra22      !temp result
  complex(kind=real64)            :: rdnsrc22  !temp result
  complex(kind=real64)            :: tu22      !contribution of homogeneous field

  !  total reflection responses for "up-down symmetric" media
  call set_fresnelcoef_TM(kappa,pvert)
  call get_rt_abv(rt22EH,ra22,rdnsrc22,kappa,tu22,set_trans_abvTM)
  refETM = (1. + ra22) * (rdnsrc22 - 1.) * rt22EH

  !put everything together (except factor -1/4pi, which is multiplied at the very end)
  res = sqrt(pvert(ilayrec)*epsh(ilaysrc) / (epsh(ilayrec)*pvert(ilaysrc))) * refETM

endfunction iabvB1TM


!---------------------------------------------------------
! integrand B0 TE for receivers below source
!---------------------------------------------------------
function iblwB0TE(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa * iblwB1TE(kappa)

endfunction iblwB0TE


!---------------------------------------------------------
! integrand B0 TM for receivers below source
!---------------------------------------------------------
function iblwB0TM(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa * iblwB1TM(kappa)

endfunction iblwB0TM


!---------------------------------------------------------
! integrand B1 TE for receivers below source
!---------------------------------------------------------
function iblwB1TE(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  !internal variables
  complex(kind=real64)            :: refETE    !temp result
  complex(kind=real64)            :: rt11EH    !temp result
  complex(kind=real64)            :: rb11      !temp result
  complex(kind=real64)            :: rupsrc11  !temp result
  complex(kind=real64)            :: td11      !contribution of homogeneous field

  !  total reflection responses for "up-down symmetric" media
  call set_fresnelcoef_TE(kappa,pvert)
  call get_rt_blw(rt11EH,rb11,rupsrc11,kappa,td11,set_trans_blwTE)
  refETE = (1. + rb11) * (1. - rupsrc11) * rt11EH

  !put everything together (except factor -1/4pi, which is multiplied at the very end)
  res = sqrt(pvert(ilaysrc) / pvert(ilayrec)) * refETE

endfunction iblwB1TE


!---------------------------------------------------------
! integrand B1 TM for receivers below source
!---------------------------------------------------------
function iblwB1TM(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  !internal variables
  complex(kind=real64)            :: refETM    !temp result
  complex(kind=real64)            :: rt22EH    !temp result
  complex(kind=real64)            :: rb22      !temp result
  complex(kind=real64)            :: rupsrc22  !temp result
  complex(kind=real64)            :: td22      !contribution of homogeneous field

  !  total reflection responses for "up-down symmetric" media
  call set_fresnelcoef_TM(kappa,pvert)
  call get_rt_blw(rt22EH,rb22,rupsrc22,kappa,td22,set_trans_blwTM)
  refETM = (1. + rb22) * (1. - rupsrc22) * rt22EH

  !put everything together (except factor -1/4pi, which is multiplied at the very end)
  res = sqrt(pvert(ilayrec)*epsh(ilaysrc) / (epsh(ilayrec)*pvert(ilaysrc))) * refETM

endfunction iblwB1TM


!---------------------------------------------------------
! integrand C0 TE for receivers above source
!---------------------------------------------------------
function iabvC0TE(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa * iabvC1TE(kappa)

endfunction iabvC0TE


!---------------------------------------------------------
! integrand C0 TM for receivers above source
!---------------------------------------------------------
function iabvC0TM(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa * iabvC1TM(kappa)

endfunction iabvC0TM


!---------------------------------------------------------
! integrand C1 TE for receivers above source
!---------------------------------------------------------
function iabvC1TE(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  !internal variables
  complex(kind=real64)            :: refHTE    !temp result
  complex(kind=real64)            :: rt11EH    !temp result
  complex(kind=real64)            :: ra11      !temp result
  complex(kind=real64)            :: rdnsrc11  !temp result
  complex(kind=real64)            :: tu11      !contribution of homogeneous field

  !  total reflection responses for "up-down symmetric" media
  call set_fresnelcoef_TE(kappa,pvert)
  call get_rt_abv(rt11EH,ra11,rdnsrc11,kappa,tu11,set_trans_abvTE)
  refHTE = (1. - ra11) * (rdnsrc11 - 1.) * rt11EH

  !put everything together (except factor -1/4pi, which is multiplied at the very end)
  res = - (sqrt(pvert(ilayrec) * pvert(ilaysrc)) / dmu0) * refHTE

endfunction iabvC1TE


!---------------------------------------------------------
! integrand C1 TM for receivers above source
!---------------------------------------------------------
function iabvC1TM(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  !internal variables
  complex(kind=real64)            :: refHTM    !temp result
  complex(kind=real64)            :: rt22EH    !temp result
  complex(kind=real64)            :: ra22      !temp result
  complex(kind=real64)            :: rdnsrc22  !temp result
  complex(kind=real64)            :: tu22      !contribution of homogeneous field

  !  total reflection responses for "up-down symmetric" media
  call set_fresnelcoef_TM(kappa,pvert)
  call get_rt_abv(rt22EH,ra22,rdnsrc22,kappa,tu22,set_trans_abvTM)
  refHTM = (1. - ra22) * (rdnsrc22 - 1.) * rt22EH

  !put everything together (except factor -1/4pi, which is multiplied at the very end)
  res = - sqrt(epsh(ilayrec)*epsh(ilaysrc) / (pvert(ilayrec)*pvert(ilaysrc))) * refHTM

endfunction iabvC1TM


!---------------------------------------------------------
! integrand C0 TE for receivers below source
!---------------------------------------------------------
function iblwC0TE(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa *  iblwC1TE(kappa)

endfunction iblwC0TE


!---------------------------------------------------------
! integrand C0 TM for receivers below source
!---------------------------------------------------------
function iblwC0TM(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa * iblwC1TM(kappa)

endfunction iblwC0TM


!---------------------------------------------------------
! integrand C1 TE for receivers below source
!---------------------------------------------------------
function iblwC1TE(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  !internal variables
  complex(kind=real64)            :: refHTE    !temp result
  complex(kind=real64)            :: rt11EH    !temp result
  complex(kind=real64)            :: rb11      !temp result
  complex(kind=real64)            :: rupsrc11  !temp result
  complex(kind=real64)            :: td11      !contribution of homogeneous field

  !  total reflection responses for "up-down symmetric" media
  call set_fresnelcoef_TE(kappa,pvert)
  call get_rt_blw(rt11EH,rb11,rupsrc11,kappa,td11,set_trans_blwTE)
  refHTE = (rb11 - 1.) * (1. - rupsrc11) * rt11EH

  !put everything together
  res = - (sqrt(pvert(ilayrec) * pvert(ilaysrc)) / dmu0) * refHTE

endfunction iblwC1TE


!---------------------------------------------------------
! integrand C1 TM for receivers below source
!---------------------------------------------------------
function iblwC1TM(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  !internal variables
  complex(kind=real64)            :: refHTM    !temp result
  complex(kind=real64)            :: rt22EH    !temp result
  complex(kind=real64)            :: rb22      !temp result
  complex(kind=real64)            :: rupsrc22  !temp result
  complex(kind=real64)            :: td22      !contribution of homogeneous field

  !  total reflection responses for "up-down symmetric" media
  call set_fresnelcoef_TM(kappa,pvert)
  call get_rt_blw(rt22EH,rb22,rupsrc22,kappa,td22,set_trans_blwTM)
  refHTM = (rb22 - 1.) * (1. - rupsrc22) * rt22EH

  !put everything together
  res = - sqrt(epsh(ilayrec)*epsh(ilaysrc) / (pvert(ilayrec)*pvert(ilaysrc))) * refHTM

endfunction iblwC1TM


!---------------------------------------------------------
! integrand Cz1 TM for receivers above source - Ez for horiz. magn. dipole
!---------------------------------------------------------
function iabvCz1TM(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa**2 * iabvC1TM(kappa)

endfunction iabvCz1TM


!---------------------------------------------------------
! integrand Cz1 TM for receivers below source - Ez for horiz. magn. dipole
!---------------------------------------------------------
function iblwCz1TM(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa**2 * iblwC1TM(kappa)

endfunction iblwCz1TM


!---------------------------------------------------------
! integrand Bz1 TE for receivers above source - Hz for horiz. magn. dipole source
!---------------------------------------------------------
function iabvBz1TE(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa**2 * iabvB1TE(kappa)

endfunction iabvBz1TE


!---------------------------------------------------------
! integrand Bz1 TE for receivers below source - Hz for horiz. magn. dipole source
!---------------------------------------------------------
function iblwBz1TE(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa**2 * iblwB1TE(kappa)

endfunction iblwBz1TE


!---------------------------------------------------------
! integrand B1 TM VED for receivers above source - vertical electric dipole source
!---------------------------------------------------------
function iabvB1TMved(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa**2 * iabvB1TM(kappa)

endfunction iabvB1TMved


!---------------------------------------------------------
! integrand B1 TM VED for receivers below source - vertical electric dipole source
!---------------------------------------------------------
function iblwB1TMved(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa**2 * iblwB1TM(kappa)

endfunction iblwB1TMved


!---------------------------------------------------------
! integrand C1 TM VED for receivers above source - vertical electric dipole source
!---------------------------------------------------------
function iabvC1TMved(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa**2 * iabvC1TM(kappa)

endfunction iabvC1TMved


!---------------------------------------------------------
! integrand C1 TM VED for receivers below source - vertical electric dipole source
!---------------------------------------------------------
function iblwC1TMved(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa**2 * iblwC1TM(kappa)

endfunction iblwC1TMved


!---------------------------------------------------------
! integrand C0 TM VED for receivers above source - vertical electric dipole source
!---------------------------------------------------------
function iabvC0TMved(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa**3 * iabvC1TM(kappa)

endfunction iabvC0TMved


!---------------------------------------------------------
! integrand C0 TM VED for receivers below source - vertical electric dipole source
!---------------------------------------------------------
function iblwC0TMved(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa**3 * iblwC1TM(kappa)

endfunction iblwC0TMved


!---------------------------------------------------------
! integrand A1 TE VMD for receivers above source - vertical magnetic dipole source
!---------------------------------------------------------
function iabvA1TEvmd(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa**2 * iabvA1TE(kappa)

endfunction iabvA1TEvmd


!---------------------------------------------------------
! integrand A1 TE VMD for receivers below source - vertical magnetic dipole source
!---------------------------------------------------------
function iblwA1TEvmd(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa**2 * iblwA1TE(kappa)

endfunction iblwA1TEvmd


!---------------------------------------------------------
! integrand D1 TE VMD for receivers above source - vertical magnetic dipole source
!---------------------------------------------------------
function iabvD1TEvmd(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa**2 * iabvD1TE(kappa)

endfunction iabvD1TEvmd


!---------------------------------------------------------
! integrand D1 TE VMD for receivers below source - vertical magnetic dipole source
!---------------------------------------------------------
function iblwD1TEvmd(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa**2 * iblwD1TE(kappa)

endfunction iblwD1TEvmd


!---------------------------------------------------------
! integrand A0 TE VMD for receivers above source - vertical magnetic dipole source
!---------------------------------------------------------
function iabvA0TEvmd(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa**3 * iabvA1TE(kappa)

endfunction iabvA0TEvmd


!---------------------------------------------------------
! integrand A0 TE VMD for receivers below source - vertical magnetic dipole source
!---------------------------------------------------------
function iblwA0TEvmd(kappa) result(res)

  implicit none

  !external variable
  real(kind=real64)  :: kappa
  complex(kind=real64)            :: res

  res = kappa**3 * iblwA1TE(kappa)

endfunction iblwA0TEvmd


!---------------------------------------------------------
! "RT-reflection amplitude" TE or TM mode for receivers above source
!  the same recursions are used for TE and TM
!  RS 2010
!---------------------------------------------------------
subroutine get_rt_abv(rtEH,ra,rdnsrc,kappa,tu,transfunc)

  implicit none

  !external variables
  complex(kind=real64)  :: rtEH
  complex(kind=real64)  :: ra
  complex(kind=real64)  :: rdnsrc
  real(kind=real64)     :: kappa
  complex(kind=real64)  :: tu
  interface
    subroutine transfunc() !set_trans_abvTE or set_trans_abvTM
    endsubroutine transfunc
  end interface  

  !internal variables
  complex(kind=real64)           :: rupsrc
  complex(kind=real64)           :: ru, ru_up


  !transmission factors through all layers for receivers above source
  call transfunc()


  !"downward" reflection response from downward radiation of source
  !start in lower halfspace: no reflection response from there
  rdnsrc = refl_response_dn(ilaysrc,nlay-1,trans_below_src,ilaysrc,zeroc)

  !"upward" reflection response from upward radiation of source
  !start in upper halfspace: no reflection response from there
  rupsrc = refl_response_up(2,ilaysrc,trans_above_src,2,zeroc)


  !---------------------------------------------------------
  !specific for receivers above source
  !---------------------------------------------------------

  !"up" reflection response from above receivers
  ra = refl_response_up(2,ilayrec,trans_above_rec,2,zeroc)

  !"down" reflection response from between source and receivers
  !"reflection response" from the "layer" between source and just below next interface is zero
  !--> start computing response from just below next interface
  ru = refl_response_dn(ilayrec,ilaysrc-1,trans_below_rec,ilayrec,zeroc)


  !"up" transmission response from between source and receivers - iterate from receiver to source layer
  ru_up = zeroc
  tu = trans_response_up(ilayrec,ilaysrc,trans_below_rec,ilayrec,trans_below_rec(ilayrec),ru_up)


  !-----------------------------------------------------------------
  ! (nearly) total reflection response for "up-down symmetric" media
  !-----------------------------------------------------------------

  !valid for horizontal electric dipole or vertical magnetic dipole sources!!!
  rtEH = (1._real64/(1._real64 - ru*ra)) * tu * (1._real64 / (1._real64 - rdnsrc*rupsrc))


endsubroutine get_rt_abv


!---------------------------------------------------------
! "RT-reflection amplitude" TE or TM mode for receivers below source
!  the same recursions are used for TE and TM
!  RS 2010
!---------------------------------------------------------
subroutine get_rt_blw(rtEH,rb,rupsrc,kappa,td,transfunc)

  implicit none

  !external variables
  complex(kind=real64)  :: rtEH
  complex(kind=real64)  :: rb
  complex(kind=real64)  :: rupsrc
  real(kind=real64)     :: kappa
  complex(kind=real64)  :: td
  interface
    subroutine transfunc() !set_trans_blwTE or set_trans_blwTM
    endsubroutine transfunc
  end interface  

  !internal variables
  complex(kind=real64)           :: rdnsrc
  complex(kind=real64)           :: rd,rd_dn


  !transmission factors through all layers for receivers below source
  call transfunc()


  !"downward" reflection response from downward radiation of source
  !start in lower halfspace: no reflection response from there
  !also have to say which transmission coeff. to use
  rdnsrc = refl_response_dn(ilaysrc,nlay-1,trans_below_src,ilaysrc,zeroc)


  !"upward" reflection response from upward radiation of source
  !start in upper halfspace: no reflection response from there
  rupsrc = refl_response_up(2,ilaysrc,trans_above_src,2,zeroc)


  !---------------------------------------------------------
  !specific for receivers below source
  !---------------------------------------------------------

  !"up" reflection response from between source and receivers
  !"reflection response" from the "layer" between source and just above next interface is zero
  !--> start computing response from just above next interface
  rd = refl_response_up(ilaysrc+1,ilayrec,trans_above_rec,ilaysrc,zeroc)

  !"down" reflection response from below receivers
  rb = refl_response_dn(ilayrec,nlay-1,trans_below_rec,ilayrec,zeroc)


  !"down" transmission response from between source and receivers
  rd_dn = zeroc
  td = trans_response_dn(ilaysrc,ilayrec,trans_above_rec,ilaysrc,trans_above_rec(ilayrec),rd_dn)


  !valid for horizontal electric dipole or vertical magnetic dipole sources!!!
  rtEH = 1._real64/(1._real64 - rd*rb) * td * (1._real64 / (1._real64 - rupsrc*rdnsrc))


endsubroutine get_rt_blw


!---------------------------------------------------------
! EM1D function realval(intfunc,kappa)
!
! get real part of an integrand
!
! Rita Streich 2009
!---------------------------------------------------------
real(kind=real64) function realval(intfunc,kappa)

  implicit none

  !external variables
  complex(kind=real64),external  :: intfunc  !function to be integrated
  real(kind=real64),intent(in)   :: kappa    !point at which to evaluate function

  !internal variables
  complex(kind=real64)           :: funval    !complex function value

  !qintfunc defined globally in module!
  funval = intfunc(kappa)

  realval = real(funval,kind=real64)

endfunction realval


!---------------------------------------------------------
! EM1D function imagval(intfunc,kappa)
!
! get imaginary part of an integrand
!
! Rita Streich 2009
!---------------------------------------------------------
real(kind=real64) function imagval(intfunc,kappa)

  implicit none

  !external variables
  complex(kind=real64),external  :: intfunc  !function to be integrated
  real(kind=real64),intent(in)   :: kappa    !point at which to evaluate function

  !internal variables
  complex(kind=real64)           :: funval    !complex function value

  !qintfunc defined globally in module!
  funval = intfunc(kappa)

  imagval = aimag(funval)

endfunction imagval



!---------------------------------------------------------
! EM1D subroutine set_fresnelcoef_TE:
! set vertical wavenumbers and TE interface reflection and transmission coefficients
! RS 2010
!---------------------------------------------------------
subroutine set_fresnelcoef_TE(kappa,pvert)

  implicit none

  !external variable
  real(kind=real64),intent(in)   :: kappa   !horiz. wavenumber
  complex(kind=real64),dimension(:),pointer :: pvert

  !internal variables
  real(kind=real64)              :: kappasq !horiz. wavenumber squared
  real(kind=real64)              :: phsq    !horizontal slowness squared
  integer(kind=int32)            :: ilay    !layer counter


  select case (ikap)
  case (1:filtlen)

    pvert => pvertall1(:,ikap)
    rup => rupallTE(:,ikap)
    rdn => rdnallTE(:,ikap)
    tup => tupallTE(:,ikap)
    tdn => tdnallTE(:,ikap)

  case default
    pvert => pvertall1(:,0)
    rup => rupallTE(:,0)
    rdn => rdnallTE(:,0)
    tup => tupallTE(:,0)
    tdn => tdnallTE(:,0)

    kappasq = kappa**2
    phsq = kappasq / omegasq

    !vertical slowness for each layer
    pvert(:) = sqrt(epsmuh(:) - phsq)

    !reflection and transmission coefficients through interfaces
    do ilay = 1,nlay-1
      rup(ilay) = (pvert(ilay+1) - pvert(ilay)) / (pvert(ilay+1) + pvert(ilay))
      if(pvert(ilay).eq.pvert(ilay+1)) then
        tup(ilay) = 1._real64
      else
        tup(ilay) = 2._real64*sqrt(pvert(ilay)*pvert(ilay+1)) / (pvert(ilay+1) + pvert(ilay))
      endif
    enddo
    rdn = -rup
    tdn = tup

  end select

endsubroutine set_fresnelcoef_TE



!---------------------------------------------------------
! EM1D subroutine set_fresnelcoef_TM:
! set vertical wavenumbers and TM interface reflection and transmission coefficients
! RS 2010
!---------------------------------------------------------
subroutine set_fresnelcoef_TM(kappa,pvert)

  implicit none

  !external variable
  real(kind=real64),intent(in)   :: kappa   !horiz. wavenumber
  complex(kind=real64),dimension(:),pointer :: pvert

  !internal variables
  real(kind=real64)              :: kappasq !horiz. wavenumber squared
  real(kind=real64)              :: phsq    !horizontal slowness squared
  integer(kind=int32)            :: ilay    !layer counter

  select case (ikap)
  case (1:filtlen)

    pvert => pvertall2(:,ikap)
    rup => rupallTM(:,ikap)
    rdn => rdnallTM(:,ikap)
    tup => tupallTM(:,ikap)
    tdn => tdnallTM(:,ikap)

  case default

    pvert => pvertall2(:,0)
    rup => rupallTM(:,0)
    rdn => rdnallTM(:,0)
    tup => tupallTM(:,0)
    tdn => tdnallTM(:,0)

    kappasq = kappa**2
    phsq = kappasq / omegasq

    !vertical slowness for each layer
    pvert(:) = sqrt(epsmuh(:) - epsmuratio(:)*phsq)
    !fix against acausal propagation (OK?)
    pvert = cmplx(real(pvert,kind=real64),abs(aimag(pvert)),kind=real64)

    !reflection and transmission coefficients through interfaces
    do ilay = 1,nlay-1
      rup(ilay) = (epsh(ilay+1)*pvert(ilay) - epsh(ilay)*pvert(ilay+1)) / &
                  (epsh(ilay+1)*pvert(ilay) + epsh(ilay)*pvert(ilay+1))
       	!WARNING: minus is not in formula 121 of Løseth and Ursin, but setting a minus here
	!  makes field values reasonable for receivers in different layer than source
	!  (checked on 2-layer model: integrals A1TE and A1TM have to nearly cancel out)
      if(pvert(ilay).eq.pvert(ilay+1)) then
        tup(ilay) = 1._real64
      else
        tup(ilay) = -2._real64*sqrt(epsh(ilay)*epsh(ilay+1)*pvert(ilay)*pvert(ilay+1)) / &
                      (epsh(ilay+1)*pvert(ilay) + epsh(ilay)*pvert(ilay+1))
      endif
    enddo
    rdn = -rup
    tdn = tup

  end select

endsubroutine set_fresnelcoef_TM


!---------------------------------------------------------
! EM1D subroutine set_trans_abvTE
! set transmission coefficients through homogeneous regions inside layer stack
! for receivers above source
! RS 2010-2011
!RS 06/2011: pvert1 corresponds to TE, pvert2 to TM
!---------------------------------------------------------
subroutine set_trans_abvTE()

  implicit none

  !transmission through homogeneous layers above source
  trans_above_src(2:ilaysrc) = exp(jomega*pvert(2:ilaysrc) * dz_above_src )

  !transmission through homogeneous layers below source
  trans_below_src(ilaysrc:nlay-1) = exp(jomega*pvert(ilaysrc:nlay-1) * dz_below_src)

  !transmission through homogeneous layers above receivers
  trans_above_rec(2:ilayrec) = exp(jomega*pvert(2:ilayrec) * dz_above_rec)

  !transmission through homogeneous regions between source and receiver
  trans_below_rec(ilayrec:ilaysrc) = exp(jomega*pvert(ilayrec:ilaysrc) * dz_below_rec)

endsubroutine set_trans_abvTE


!---------------------------------------------------------
! EM1D subroutine set_trans_abvTE
! set transmission coefficients through homogeneous regions inside layer stack
! for receivers above source
! RS 2010-2011
!RS 06/2011: pvert1 corresponds to TE, pvert2, to TM
!---------------------------------------------------------
subroutine set_trans_abvTM()

  implicit none
  
  !transmission through homogeneous layers above source
  trans_above_src(2:ilaysrc) = exp(jomega*pvert(2:ilaysrc) * dz_above_src )

  !transmission through homogeneous layers below source
  trans_below_src(ilaysrc:nlay-1) = exp(jomega*pvert(ilaysrc:nlay-1) * dz_below_src)

  !transmission through homogeneous layers above receivers
  trans_above_rec(2:ilayrec) = exp(jomega*pvert(2:ilayrec) * dz_above_rec)

  !transmission through homogeneous regions between source and receiver
  trans_below_rec(ilayrec:ilaysrc) = exp(jomega*pvert(ilayrec:ilaysrc) * dz_below_rec)

endsubroutine set_trans_abvTM


!---------------------------------------------------------
! EM1D subroutine set_trans_blw
! set transmission coefficients through homogeneous regions inside layer stack
! for receivers below source
! RS 2010-2011
!RS 06/2011: pvert1 corresponds to TE, pvert2, to TM
!---------------------------------------------------------
subroutine set_trans_blwTE()

  implicit none

  !transmission through homogeneous layers above source
  trans_above_src(2:ilaysrc) = exp(jomega*pvert(2:ilaysrc) * dz_above_src )

  !transmission through homogeneous layers below source
  trans_below_src(ilaysrc:nlay-1) = exp(jomega*pvert(ilaysrc:nlay-1) * dz_below_src)

  !transmission through homogeneous regions between source and receiver
  trans_above_rec(ilaysrc:ilayrec) = exp(jomega*pvert(ilaysrc:ilayrec) * dz_above_rec)

  !transmission through homogeneous layers below receivers
  trans_below_rec(ilayrec:nlay-1) = exp(jomega*pvert(ilayrec:nlay-1) * dz_below_rec)

endsubroutine set_trans_blwTE

!---------------------------------------------------------
! EM1D subroutine set_trans_blw
! set transmission coefficients through homogeneous regions inside layer stack
! for receivers below source
! RS 2010-2011
!RS 06/2011: pvert1 corresponds to TE, pvert2, to TM
!---------------------------------------------------------
subroutine set_trans_blwTM()

  implicit none

  !transmission through homogeneous layers above source
  trans_above_src(2:ilaysrc) = exp(jomega*pvert(2:ilaysrc) * dz_above_src )

  !transmission through homogeneous layers below source
  trans_below_src(ilaysrc:nlay-1) = exp(jomega*pvert(ilaysrc:nlay-1) * dz_below_src)

  !transmission through homogeneous regions between source and receiver
  trans_above_rec(ilaysrc:ilayrec) = exp(jomega*pvert(ilaysrc:ilayrec) * dz_above_rec)

  !transmission through homogeneous layers below receivers
  trans_below_rec(ilayrec:nlay-1) = exp(jomega*pvert(ilayrec:nlay-1) * dz_below_rec)

endsubroutine set_trans_blwTM

!---------------------------------------------------------
! EM1D function refl_response_dn
! "downward" reflection response, eq. 70a in Loeseth & Ursin 2007
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function refl_response_dn(ilaytop,ilaybot,trans_blw,startidx,refstart)

  implicit none

  !external variables
  integer(kind=int32)   :: ilaytop,ilaybot  !layer indices: start recursion from ilaybot, go to ilaytop
  integer(kind=int32),intent(in)      :: startidx  !start index of transm. coeff. vector
  complex(kind=real64),dimension(startidx:),intent(in) :: trans_blw !transm. coeff. either with respect to source or receiver
  complex(kind=real64),intent(in) :: refstart  !start value for refl. response, input this so function can be called step by step

  !internal variable
  integer(kind=int32)   :: ilay    !layer counter

  refl_response_dn = refstart

  !from lowest non-halfspace layer to the layer containing the source
  do ilay=ilaybot,ilaytop,-1
    refl_response_dn = trans_blw(ilay) * &
               (rdn(ilay) + tup(ilay)*refl_response_dn*(1._real64/(1._real64 - rup(ilay)*refl_response_dn)) * tdn(ilay)) &
               * trans_blw(ilay)
  enddo

endfunction refl_response_dn


!---------------------------------------------------------
! EM1D function refl_response_up
! "upward" reflection response, eq. 71a in Loeseth & Ursin 2007
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function refl_response_up(ilaytop,ilaybot,trans_abv,startidx,refstart)

  implicit none

  !external variables
  integer(kind=int32)   :: ilaytop,ilaybot  !layer indices: start recursion from ilaytop, go to ilaybot
  integer(kind=int32),intent(in)      :: startidx  !start index of transm. coeff. vector
  complex(kind=real64),dimension(startidx:),intent(in) :: trans_abv !transm. coeff. either with respect to source or receiver
  complex(kind=real64),intent(in) :: refstart  !start value for refl. response, input this so function can be called step by step

  !internal variable
  integer(kind=int32)   :: ilay    !layer counter


  refl_response_up = refstart

  do ilay=ilaytop,ilaybot
    refl_response_up = trans_abv(ilay) * &
       (rup(ilay-1) + tdn(ilay-1)*refl_response_up*(1._real64/(1._real64 - rdn(ilay-1)*refl_response_up)) * tup(ilay-1)) &
       * trans_abv(ilay)
  enddo

endfunction refl_response_up


!---------------------------------------------------------
! EM1D function trans_response_up
! "upward" transmission response from between source and receivers, for receivers above source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function trans_response_up(ilaytop,ilaybot,trans_blw,startidx,trstart,ru_up) result(tu)

  implicit none

  !external variables
  integer(kind=int32)   :: ilaytop,ilaybot  !layer indices: start recursion from ilaytop, go to ilaybot
  integer(kind=int32),intent(in)      :: startidx  !start index of transm. coeff. vector
  complex(kind=real64),dimension(startidx:),intent(in) :: trans_blw !transm. coeff. with respect to receiver layer
  complex(kind=real64),intent(in) :: trstart  !start value for transm. response, input so function can be called step by step
  complex(kind=real64),intent(inout)  :: ru_up     !refl. response, input and return!


  !internal variable
  integer(kind=int32)   :: ilay    !layer counter


  !tu = trans_blw(ilaytop)
  !ru_up = 0._real64
  tu = trstart

  !"up" --> recursion from top to bottom
  do ilay=ilaytop+1,ilaybot
    tu = tu * (1._real64/(1._real64 - rdn(ilay-1)*ru_up)) * tup(ilay-1) * trans_blw(ilay)

    ru_up = trans_blw(ilay) * &
              (rup(ilay-1) + tdn(ilay-1)*ru_up*(1._real64/(1._real64 - rdn(ilay-1)*ru_up)) * tup(ilay-1) ) &
              * trans_blw(ilay)
  enddo

endfunction trans_response_up


!---------------------------------------------------------
! EM1D function trans_response_dn
! "downward" transmission response from between source and receivers, for receivers below source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function trans_response_dn(ilaytop,ilaybot,trans_abv,startidx,trstart,rd_dn) result(td)

  implicit none

  !external variables
  integer(kind=int32)   :: ilaytop,ilaybot  !layer indices: start recursion from ilaytop, go to ilaybot
  integer(kind=int32),intent(in)      :: startidx  !start index of transm. coeff. vector
  complex(kind=real64),dimension(startidx:),intent(in) :: trans_abv !transm. coeff. with respect to receiver layer
  complex(kind=real64),intent(in)     :: trstart   !start value for rtansm. response, input to enable step-by-step calling
  complex(kind=real64),intent(inout)  :: rd_dn   !refl. response, input and return!

  !internal variable
  integer(kind=int32)   :: ilay    !layer counter


  !td = trans_abv(ilaybot)
  !rd_dn = 0._real64
  td = trstart

  !"down" --> recursion from bottom to top
  do ilay=ilaybot-1,ilaytop,-1
    td = td * (1._real64/(1._real64 - rup(ilay)*rd_dn)) * tdn(ilay) * trans_abv(ilay)

    rd_dn = trans_abv(ilay) * &
              (rdn(ilay) + tup(ilay)*rd_dn*(1._real64/(1._real64 - rup(ilay)*rd_dn)) * tdn(ilay) ) &
              * trans_abv(ilay)
  enddo

endfunction trans_response_dn

