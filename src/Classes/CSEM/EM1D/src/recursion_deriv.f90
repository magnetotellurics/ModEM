!---------------------------------------------------------
! derivative of integrand A0 TM with respect to layer ilaym for receivers above source
! Bessel function not included here, but handled by zhankl
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iabvA0TMderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iabvA1TMderiv(kappa)
endfunction iabvA0TMderiv

complex(kind=real64) function iabvA0TMderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iabvA1TMderivh(kappa)
endfunction iabvA0TMderivh

complex(kind=real64) function iabvA0TMderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iabvA1TMderivv(kappa)
endfunction iabvA0TMderivv

!---------------------------------------------------------
! derivative of integrand A0 TM with respect to layer ilaym for receivers below source
! Bessel function not included here, but handled by zhankl
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iblwA0TMderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iblwA1TMderiv(kappa)
endfunction iblwA0TMderiv

complex(kind=real64) function iblwA0TMderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iblwA1TMderivh(kappa)
endfunction iblwA0TMderivh

complex(kind=real64) function iblwA0TMderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iblwA1TMderivv(kappa)
endfunction iblwA0TMderivv

!---------------------------------------------------------
! derivative with respect to layer ilaym of integrand A1 TM for receivers above source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iabvA1TMderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRaTM ! dRa / depsilon_m
  complex(kind=real64)            :: refRaTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRaTM = get_drefR_abv(1._real64,1._real64,set_trans_abvTM, fact_exp_tmiso,fact_isotm,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = sqrt(pvert(ilayrec)*pvert(ilaysrc) / (epsh(ilayrec)*epsh(ilaysrc))) * drefRaTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRaTM = iabvA1TM(kappa) !factor sqrt(pvert pvert epsr epsr ...) is contained in integral value
    termadd = 0.5_real64 * (-1./epsh(ilaym) + dmu0/(2.*pvert(ilaym)**2) ) * refRaTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res + termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res + termadd
  endif

endfunction iabvA1TMderiv


!---------------------------------------------------------
! derivative with respect to layer ilaym of integrand A1 TM for receivers above source
! VTI for epsh
! RS 2011
!---------------------------------------------------------
complex(kind=real64) function iabvA1TMderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRaTM ! dRa / depsilon_m
  complex(kind=real64)            :: refRaTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m
  real(kind=real64)               :: phsq     !horiz. slowness

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRaTM = get_drefR_abv(1._real64,1._real64,set_trans_abvTM, fact_exp_tmvtih,fact_vtitmh,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = sqrt(pvert(ilayrec)*pvert(ilaysrc) / (epsh(ilayrec)*epsh(ilaysrc))) * drefRaTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRaTM = iabvA1TM(kappa) !factor sqrt(pvert pvert epsr epsr ...) is contained in integral value
    phsq = (kappa / omega)**2
    termadd = 0.5_real64 * (-1./epsh(ilaym) + (dmu0 - phsq/epsv(ilaym))/(2.*pvert(ilaym)**2) ) * refRaTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res + termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res + termadd
  endif

endfunction iabvA1TMderivh


!---------------------------------------------------------
! derivative with respect to layer ilaym of integrand A1 TM for receivers above source
! VTI for epsv
! RS 2011
!---------------------------------------------------------
complex(kind=real64) function iabvA1TMderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRaTM ! dRa / depsilon_m
  complex(kind=real64)            :: refRaTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m
  real(kind=real64)               :: phsq     !horiz. slowness

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRaTM = get_drefR_abv(1._real64,1._real64,set_trans_abvTM, fact_exp_tmvtiv,fact_vtitmv,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = sqrt(pvert(ilayrec)*pvert(ilaysrc) / (epsh(ilayrec)*epsh(ilaysrc))) * drefRaTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRaTM = iabvA1TM(kappa) !factor sqrt(pvert pvert epsr epsr ...) is contained in integral value
    phsq = (kappa / omega)**2
    termadd = (epsh(ilaym)*phsq / (4.* (epsv(ilaym)*pvert(ilaym))**2)) * refRaTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res + termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res + termadd
  endif

endfunction iabvA1TMderivv


!---------------------------------------------------------
! derivative of integrand A1 TM with respect to layer ilaym for receivers below source
! Bessel function not included here, but handled by zhankl
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iblwA1TMderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRaTM ! dRa / depsilon_m
  complex(kind=real64)            :: refRaTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRaTM = get_drefR_blw(1._real64,1._real64,set_trans_blwTM, fact_exp_tmiso,fact_isotm,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = sqrt(pvert(ilayrec)*pvert(ilaysrc) / (epsh(ilayrec)*epsh(ilaysrc))) * drefRaTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRaTM = iblwA1TM(kappa) !factor sqrt(pvert pvert epsr epsr ...) is contained in integral value
    termadd = 0.5_real64 * (-1./epsh(ilaym) + dmu0/(2.*pvert(ilaym)**2) ) * refRaTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res + termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res + termadd
  endif

endfunction iblwA1TMderiv


!---------------------------------------------------------
! derivative of integrand A1 TM with respect to layer ilaym for receivers below source
! VTI for epsh
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iblwA1TMderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRaTM ! dRa / depsilon_m
  complex(kind=real64)            :: refRaTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m
  real(kind=real64)               :: phsq     !horiz. slowness

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRaTM = get_drefR_blw(1._real64,1._real64,set_trans_blwTM, fact_exp_tmvtih,fact_vtitmh,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = sqrt(pvert(ilayrec)*pvert(ilaysrc) / (epsh(ilayrec)*epsh(ilaysrc))) * drefRaTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRaTM = iblwA1TM(kappa) !factor sqrt(pvert pvert epsr epsr ...) is contained in integral value
    phsq = (kappa / omega)**2
    termadd = 0.5_real64 * (-1./epsh(ilaym) + (dmu0 - phsq/epsv(ilaym))/(2.*pvert(ilaym)**2) ) * refRaTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res + termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res + termadd
  endif

endfunction iblwA1TMderivh


!---------------------------------------------------------
! derivative of integrand A1 TM with respect to layer ilaym for receivers below source
! VTI for epsv
! RS 2011
!---------------------------------------------------------
complex(kind=real64) function iblwA1TMderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRaTM ! dRa / depsilon_m
  complex(kind=real64)            :: refRaTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m
  real(Kind=real64)               :: phsq     !horiz. slowness

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRaTM = get_drefR_blw(1._real64,1._real64,set_trans_blwTM, fact_exp_tmvtiv,fact_vtitmv,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = sqrt(pvert(ilayrec)*pvert(ilaysrc) / (epsh(ilayrec)*epsh(ilaysrc))) * drefRaTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRaTM = iblwA1TM(kappa) !factor sqrt(pvert pvert epsr epsr ...) is contained in integral value
    phsq = (kappa / omega)**2
    termadd = (epsh(ilaym)*phsq / (4.* (epsv(ilaym)*pvert(ilaym))**2)) * refRaTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res + termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res + termadd
  endif

endfunction iblwA1TMderivv


!---------------------------------------------------------
! derivative of integrand A0 TE with respect to layer ilaym for receivers above source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iabvA0TEderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iabvA1TEderiv(kappa)
endfunction iabvA0TEderiv

!---------------------------------------------------------
! derivative of integrand A0 TE with respect to layer ilaym for receivers below source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iblwA0TEderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iblwA1TEderiv(kappa)
endfunction iblwA0TEderiv


!---------------------------------------------------------
! derivative of integrand A1 TE with respect to layer ilaym for receivers above source
! also valid for VTI
! RS 2010-2011
!---------------------------------------------------------
complex(kind=real64) function iabvA1TEderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRaTE ! dRa / depsilon_m
  complex(kind=real64)            :: refRaTE  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TE(kappa,pvert)

  !  derivative of total reflection response
  drefRaTE = get_drefR_abv(1._real64,1._real64,set_trans_abvTE, fact_exp_teiso,fact_isote,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = (dmu0 / sqrt(pvert(ilayrec)*pvert(ilaysrc))) * drefRaTE

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRaTE = iabvA1TE(kappa) !factor sqrt(pvert pvert epsr epsr ...) is contained in integral value
    termadd = (dmu0 / (4.*pvert(ilaym)**2) ) * refRaTE
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res - termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res - termadd
  endif

endfunction iabvA1TEderiv


!---------------------------------------------------------
! derivative of integrand A1 TE with respect to layer ilaym for receivers below source
! also valid for VTI
! RS 2010-2011
!---------------------------------------------------------
complex(kind=real64) function iblwA1TEderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRaTE ! dRa / depsilon_m
  complex(kind=real64)            :: refRaTE  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TE(kappa,pvert)

  !  derivative of total reflection response
  drefRaTE = get_drefR_blw(1._real64,1._real64,set_trans_blwTE, fact_exp_teiso,fact_isote,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = (dmu0 / sqrt(pvert(ilayrec)*pvert(ilaysrc))) * drefRaTE

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRaTE = iblwA1TE(kappa) !factor sqrt(pvert pvert epsr epsr ...) is contained in integral value
    termadd = (dmu0 / (4.*pvert(ilaym)**2) ) * refRaTE
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res - termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res - termadd
  endif

endfunction iblwA1TEderiv


!---------------------------------------------------------
! derivative of integrand Az1 TE with respect to layer ilaym for receivers above source - for Hz
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iabvAz1TEderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iabvA1TEderiv(kappa)
endfunction iabvAz1TEderiv

!---------------------------------------------------------
! derivative of integrand Az1 TE with respect to layer ilaym for receivers below source - for Hz
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iblwAz1TEderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iblwA1TEderiv(kappa)
endfunction iblwAz1TEderiv


!---------------------------------------------------------
! derivative of integrand D0 TE with respect to layer ilaym for receivers above source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iabvD0TEderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iabvD1TEderiv(kappa)
endfunction iabvD0TEderiv

!---------------------------------------------------------
! derivative of integrand D0 TE with respect to layer ilaym for receivers below source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iblwD0TEderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iblwD1TEderiv(kappa)
endfunction iblwD0TEderiv


!---------------------------------------------------------
! derivative of integrand D1 TE with respect to layer ilaym for receivers above source
! also valid for VTI
! RS 2010-2011
!---------------------------------------------------------
complex(kind=real64) function iabvD1TEderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRdTE ! dRa / depsilon_m
  complex(kind=real64)            :: refRdTE  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TE(kappa,pvert)

  !  derivative of total reflection response
  drefRdTE = get_drefR_abv(-1._real64,1._real64,set_trans_abvTE, fact_exp_teiso,fact_isote,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = - sqrt(pvert(ilayrec)/pvert(ilaysrc)) * drefRdTE

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRdTE = iabvD1TE(kappa) !factor sqrt(pvert pvert epsr epsr ...) is contained in integral value
    termadd = (dmu0 / (4.*pvert(ilaym)**2)) * refRdTE !"minus" is contained in refRdTE
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res + termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res - termadd !"minus" is contained in refRdTE
  endif

endfunction iabvD1TEderiv


!---------------------------------------------------------
! derivative of integrand D1 TE with respect to layer ilaym for receivers below source
! also valid for VTI
! RS 2010-2011
!---------------------------------------------------------
complex(kind=real64) function iblwD1TEderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRdTE ! dRa / depsilon_m
  complex(kind=real64)            :: refRdTE  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TE(kappa,pvert)

  !  derivative of total reflection response
  drefRdTE = get_drefR_blw(-1._real64,1._real64,set_trans_blwTE, fact_exp_teiso,fact_isote,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = - sqrt(pvert(ilayrec)/pvert(ilaysrc)) * drefRdTE

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRdTE = iblwD1TE(kappa)
    termadd = (dmu0 / (4.*pvert(ilaym)**2)) * refRdTE !"minus" is contained in refRdTE
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res + termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res - termadd !"minus" is contained in refRdTE
  endif

endfunction iblwD1TEderiv


!---------------------------------------------------------
! derivative of integrand D0 TM with respect to layer ilaym for receivers above source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iabvD0TMderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iabvD1TMderiv(kappa)
endfunction iabvD0TMderiv

complex(kind=real64) function iabvD0TMderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iabvD1TMderivh(kappa)
endfunction iabvD0TMderivh

complex(kind=real64) function iabvD0TMderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iabvD1TMderivv(kappa)
endfunction iabvD0TMderivv

!---------------------------------------------------------
! derivative of integrand D0 TM with respect to layer ilaym for receivers below source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iblwD0TMderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iblwD1TMderiv(kappa)
endfunction iblwD0TMderiv

complex(kind=real64) function iblwD0TMderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iblwD1TMderivh(kappa)
endfunction iblwD0TMderivh

complex(kind=real64) function iblwD0TMderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iblwD1TMderivv(kappa)
endfunction iblwD0TMderivv

!---------------------------------------------------------
! derivative of integrand D1 TM with respect to layer ilaym for receivers above source
! RS 2010-2011
!---------------------------------------------------------
complex(kind=real64) function iabvD1TMderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRdTM ! dRa / depsilon_m
  complex(kind=real64)            :: refRdTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRdTM = get_drefR_abv(-1._real64,1._real64,set_trans_abvTM, fact_exp_tmiso,fact_isotm,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = - sqrt(epsh(ilayrec)*pvert(ilaysrc) / (pvert(ilayrec)*epsh(ilaysrc))) * drefRdTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRdTM = iabvD1TM(kappa) !contains a minus
    termadd = (1./(2.*epsh(ilaym)) - dmu0/(4.*pvert(ilaym)**2)) * refRdTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res + termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res - termadd
  endif

endfunction iabvD1TMderiv


!---------------------------------------------------------
! derivative of integrand D1 TM with respect to layer ilaym for receivers above source
! VTI for epsh
! RS 2010-2011
!---------------------------------------------------------
complex(kind=real64) function iabvD1TMderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRdTM ! dRa / depsilon_m
  complex(kind=real64)            :: refRdTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m
  real(kind=real64)               :: phsq     !horiz. slowness

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRdTM = get_drefR_abv(-1._real64,1._real64,set_trans_abvTM, fact_exp_tmvtih,fact_vtitmh,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = - sqrt(epsh(ilayrec)*pvert(ilaysrc) / (pvert(ilayrec)*epsh(ilaysrc))) * drefRdTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRdTM = iabvD1TM(kappa) !contains a minus
    phsq = (kappa / omega)**2
    termadd = (1./(2.*epsh(ilaym)) - (dmu0 - phsq/epsv(ilaym))/(4.*pvert(ilaym)**2)) * refRdTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res + termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res - termadd
  endif

endfunction iabvD1TMderivh


!---------------------------------------------------------
! derivative of integrand D1 TM with respect to layer ilaym for receivers above source
! VTI for epsv
! RS 2010-2011
!---------------------------------------------------------
complex(kind=real64) function iabvD1TMderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRdTM ! dRa / depsilon_m
  complex(kind=real64)            :: refRdTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m
  real(kind=real64)               :: phsq     !horiz. slowness

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRdTM = get_drefR_abv(-1._real64,1._real64,set_trans_abvTM, fact_exp_tmvtiv,fact_vtitmv,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = - sqrt(epsh(ilayrec)*pvert(ilaysrc) / (pvert(ilayrec)*epsh(ilaysrc))) * drefRdTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRdTM = iabvD1TM(kappa) !contains a minus
    phsq = (kappa / omega)**2
    !include minus here, then +/- termadd is the same as for epsh
    termadd = (- epsh(ilaym)*phsq / (4.* (epsv(ilaym)*pvert(ilaym))**2)) * refRdTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res + termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res - termadd
  endif

endfunction iabvD1TMderivv


!---------------------------------------------------------
! derivative of integrand D1 TM with respect to layer ilaym for receivers below source
! RS 2010-2011
!---------------------------------------------------------
complex(kind=real64) function iblwD1TMderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRdTM ! dRa / depsilon_m
  complex(kind=real64)            :: refRdTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRdTM = get_drefR_blw(-1._real64,1._real64,set_trans_blwTM, fact_exp_tmiso,fact_isotm,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = - sqrt(epsh(ilayrec)*pvert(ilaysrc) / (pvert(ilayrec)*epsh(ilaysrc))) * drefRdTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRdTM = iblwD1TM(kappa) !contains a minus
    termadd = (1./(2.*epsh(ilaym)) - dmu0/(4.*pvert(ilaym)**2)) * refRdTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res + termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res - termadd
  endif

endfunction iblwD1TMderiv


!---------------------------------------------------------
! derivative of integrand D1 TM with respect to layer ilaym for receivers below source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iblwD1TMderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRdTM ! dRa / depsilon_m
  complex(kind=real64)            :: refRdTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m
  real(kind=real64)               :: phsq     !horiz. slowness

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRdTM = get_drefR_blw(-1._real64,1._real64,set_trans_blwTM, fact_exp_tmvtih,fact_vtitmh,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = - sqrt(epsh(ilayrec)*pvert(ilaysrc) / (pvert(ilayrec)*epsh(ilaysrc))) * drefRdTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRdTM = iblwD1TM(kappa) !contains a minus
    phsq = (kappa / omega)**2
    termadd = (1./(2.*epsh(ilaym)) - (dmu0 - phsq/epsv(ilaym))/(4.*pvert(ilaym)**2)) * refRdTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res + termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res - termadd
  endif

endfunction iblwD1TMderivh


!---------------------------------------------------------
! derivative of integrand D1 TM with respect to layer ilaym for receivers below source
! RS 2010-2011
!---------------------------------------------------------
complex(kind=real64) function iblwD1TMderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRdTM ! dRa / depsilon_m
  complex(kind=real64)            :: refRdTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m
  real(kind=real64)               :: phsq     !horiz. slowness

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRdTM = get_drefR_blw(-1._real64,1._real64,set_trans_blwTM, fact_exp_tmvtiv,fact_vtitmv,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = - sqrt(epsh(ilayrec)*pvert(ilaysrc) / (pvert(ilayrec)*epsh(ilaysrc))) * drefRdTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRdTM = iblwD1TM(kappa) !contains a minus
    phsq = (kappa / omega)**2
    !include minus here, then +/- termadd is the same as for epsh
    termadd = (- epsh(ilaym)*phsq / (4.* (epsv(ilaym)*pvert(ilaym))**2)) * refRdTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res + termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res - termadd
  endif

endfunction iblwD1TMderivv


!---------------------------------------------------------
! derivative of integrand Dz1 TM with respect to layer ilaym for receivers above source - for Ez
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iabvDz1TMderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iabvD1TMderiv(kappa)
endfunction iabvDz1TMderiv

complex(kind=real64) function iabvDz1TMderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iabvD1TMderivh(kappa)
endfunction iabvDz1TMderivh

complex(kind=real64) function iabvDz1TMderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iabvD1TMderivv(kappa)
endfunction iabvDz1TMderivv

!---------------------------------------------------------
! derivative of integrand Dz1 TM with respect to layer ilaym for receivers below source - for Ez
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iblwDz1TMderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iblwD1TMderiv(kappa)
endfunction iblwDz1TMderiv

complex(kind=real64) function iblwDz1TMderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iblwD1TMderivh(kappa)
endfunction iblwDz1TMderivh

complex(kind=real64) function iblwDz1TMderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iblwD1TMderivv(kappa)
endfunction iblwDz1TMderivv

!---------------------------------------------------------
! derivative of integrand B0 TE with respect to layer ilaym for receivers above source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iabvB0TEderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iabvB1TEderiv(kappa)
endfunction iabvB0TEderiv

!---------------------------------------------------------
! derivative of integrand B0 TE with respect to layer ilaym for receivers below source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iblwB0TEderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iblwB1TEderiv(kappa)
endfunction iblwB0TEderiv


!---------------------------------------------------------
! derivative of integrand B1 TE with respect to layer ilaym for receivers above source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iabvB1TEderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRbTE ! dRa / depsilon_m
  complex(kind=real64)            :: refRbTE  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TE(kappa,pvert)

  !  derivative of total reflection response
  drefRbTE = get_drefR_abv(1._real64,-1._real64,set_trans_abvTE, fact_exp_teiso,fact_isote,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = sqrt(pvert(ilaysrc)/pvert(ilayrec)) * drefRbTE

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRbTE = iabvB1TE(kappa)
    termadd = (dmu0 / (4.*pvert(ilaym)**2)) * refRbTE
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res - termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res + termadd
  endif

endfunction iabvB1TEderiv


!---------------------------------------------------------
! derivative of integrand B1 TE with respect to layer ilaym for receivers below source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iblwB1TEderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRbTE ! dRa / depsilon_m
  complex(kind=real64)            :: refRbTE  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TE(kappa,pvert)

  !  derivative of total reflection response
  drefRbTE = get_drefR_blw(1._real64,-1._real64,set_trans_blwTE, fact_exp_teiso,fact_isote,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = sqrt(pvert(ilaysrc)/pvert(ilayrec)) * drefRbTE

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRbTE = iblwB1TE(kappa)
    termadd = (dmu0 / (4.*pvert(ilaym)**2)) * refRbTE
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res - termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res + termadd
  endif

endfunction iblwB1TEderiv


!---------------------------------------------------------
! derivative of integrand Bz1 TE with respect to layer ilaym for receivers above source - for Hz
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iabvBz1TEderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iabvB1TEderiv(kappa)
endfunction iabvBz1TEderiv

!---------------------------------------------------------
! derivative of integrand Bz1 TE with respect to layer ilaym for receivers below source - for Hz
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iblwBz1TEderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iblwB1TEderiv(kappa)
endfunction iblwBz1TEderiv


!---------------------------------------------------------
! derivative of integrand B0 TM with respect to layer ilaym for receivers above source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iabvB0TMderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iabvB1TMderiv(kappa)
endfunction iabvB0TMderiv

complex(kind=real64) function iabvB0TMderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iabvB1TMderivh(kappa)
endfunction iabvB0TMderivh

complex(kind=real64) function iabvB0TMderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iabvB1TMderivv(kappa)
endfunction iabvB0TMderivv

!---------------------------------------------------------
! derivative of integrand B0 TM with respect to layer ilaym for receivers below source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iblwB0TMderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iblwB1TMderiv(kappa)
endfunction iblwB0TMderiv

complex(kind=real64) function iblwB0TMderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iblwB1TMderivh(kappa)
endfunction iblwB0TMderivh

complex(kind=real64) function iblwB0TMderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iblwB1TMderivv(kappa)
endfunction iblwB0TMderivv


!---------------------------------------------------------
! derivative of integrand B1 TM with respect to layer ilaym for receivers above source
! RS 2010-2011
!---------------------------------------------------------
complex(kind=real64) function iabvB1TMderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRbTM ! dRa / depsilon_m
  complex(kind=real64)            :: refRbTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRbTM = get_drefR_abv(1._real64,-1._real64,set_trans_abvTM, fact_exp_tmiso,fact_isotm,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = sqrt(pvert(ilayrec)*epsh(ilaysrc) / (epsh(ilayrec)*pvert(ilaysrc))) * drefRbTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRbTM = iabvB1TM(kappa)
    termadd = (1./(2.*epsh(ilaym)) - dmu0/(4.*pvert(ilaym)**2)) * refRbTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res - termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res + termadd
  endif

endfunction iabvB1TMderiv


!---------------------------------------------------------
! derivative of integrand B1 TM with respect to layer ilaym for receivers above source
! RS 2010-2011
!---------------------------------------------------------
complex(kind=real64) function iabvB1TMderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRbTM ! dRa / depsilon_m
  complex(kind=real64)            :: refRbTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m
  real(kind=real64)   :: phsq  !horiz. slowness

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRbTM = get_drefR_abv(1._real64,-1._real64,set_trans_abvTM, fact_exp_tmvtih,fact_vtitmh,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = sqrt(pvert(ilayrec)*epsh(ilaysrc) / (epsh(ilayrec)*pvert(ilaysrc))) * drefRbTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRbTM = iabvB1TM(kappa)
    phsq = (kappa / omega)**2
    termadd = (1./(2.*epsh(ilaym)) - (dmu0 - phsq/epsv(ilaym))/(4.*pvert(ilaym)**2)) * refRbTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res - termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res + termadd
  endif

endfunction iabvB1TMderivh


!---------------------------------------------------------
! derivative of integrand B1 TM with respect to layer ilaym for receivers above source
! VTI for epsv
! RS 2010-2011
!---------------------------------------------------------
complex(kind=real64) function iabvB1TMderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRbTM ! dRa / depsilon_m
  complex(kind=real64)            :: refRbTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m
  real(Kind=real64)   :: phsq   !horiz. slowness

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRbTM = get_drefR_abv(1._real64,-1._real64,set_trans_abvTM, fact_exp_tmvtiv,fact_vtitmv,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = sqrt(pvert(ilayrec)*epsh(ilaysrc) / (epsh(ilayrec)*pvert(ilaysrc))) * drefRbTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRbTM = iabvB1TM(kappa)
    phsq = (kappa / omega)**2
    termadd = (- epsh(ilaym)*phsq / (4.* (epsv(ilaym) * pvert(ilaym))**2)) * refRbTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res - termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res + termadd
  endif

endfunction iabvB1TMderivv


!---------------------------------------------------------
! derivative of integrand B1 TM with respect to layer ilaym for receivers below source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iblwB1TMderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRbTM ! dRa / depsilon_m
  complex(kind=real64)            :: refRbTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRbTM = get_drefR_blw(1._real64,-1._real64,set_trans_blwTM, fact_exp_tmiso,fact_isotm,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = sqrt(pvert(ilayrec)*epsh(ilaysrc) / (epsh(ilayrec)*pvert(ilaysrc))) * drefRbTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRbTM = iblwB1TM(kappa)
    termadd = (1./(2.*epsh(ilaym)) - dmu0/(4.*pvert(ilaym)**2)) * refRbTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res - termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res + termadd
  endif

endfunction iblwB1TMderiv


!---------------------------------------------------------
! derivative of integrand B1 TM with respect to layer ilaym for receivers below source
! VTI for epsh
! RS 2011
!---------------------------------------------------------
complex(kind=real64) function iblwB1TMderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRbTM ! dRa / depsilon_m
  complex(kind=real64)            :: refRbTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m
  real(kind=real64)  :: phsq  !horiz. slowness

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRbTM = get_drefR_blw(1._real64,-1._real64,set_trans_blwTM, fact_exp_tmvtih,fact_vtitmh,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = sqrt(pvert(ilayrec)*epsh(ilaysrc) / (epsh(ilayrec)*pvert(ilaysrc))) * drefRbTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRbTM = iblwB1TM(kappa)
    phsq = (kappa / omega)**2
    termadd = (1./(2.*epsh(ilaym)) - (dmu0 - phsq/epsv(ilaym))/(4.*pvert(ilaym)**2)) * refRbTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res - termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res + termadd
  endif

endfunction iblwB1TMderivh


!---------------------------------------------------------
! derivative of integrand B1 TM with respect to layer ilaym for receivers below source
! VTI for epsv
! RS 2011
!---------------------------------------------------------
complex(kind=real64) function iblwB1TMderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRbTM ! dRa / depsilon_m
  complex(kind=real64)            :: refRbTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m
  real(kind=real64)  :: phsq  !horiz. slowness

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRbTM = get_drefR_blw(1._real64,-1._real64,set_trans_blwTM, fact_exp_tmvtiv,fact_vtitmv,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = sqrt(pvert(ilayrec)*epsh(ilaysrc) / (epsh(ilayrec)*pvert(ilaysrc))) * drefRbTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRbTM = iblwB1TM(kappa)
    phsq = (kappa / omega)**2
    termadd = (- epsh(ilaym)*phsq / (4.* (epsv(ilaym) * pvert(ilaym))**2)) * refRbTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res - termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res + termadd
  endif

endfunction iblwB1TMderivv


!---------------------------------------------------------
! derivative of integrand C0 TE with respect to layer ilaym for receivers above source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iabvC0TEderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iabvC1TEderiv(kappa)
endfunction iabvC0TEderiv

!---------------------------------------------------------
! derivative of integrand C0 TE with respect to layer ilaym for receivers below source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iblwC0TEderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iblwC1TEderiv(kappa)
endfunction iblwC0TEderiv


!---------------------------------------------------------
! derivative of integrand C1 TE with respect to layer ilaym for receivers above source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iabvC1TEderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRcTE ! dRa / depsilon_m
  complex(kind=real64)            :: refRcTE  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TE(kappa,pvert)

  !  derivative of total reflection response
  drefRcTE = get_drefR_abv(-1._real64,-1._real64,set_trans_abvTE, fact_exp_teiso,fact_isote,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = - (sqrt(pvert(ilayrec)*pvert(ilaysrc)) / dmu0) * drefRcTE

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRcTE = iabvC1TE(kappa) !contains a minus
    termadd = (dmu0 / (4.*pvert(ilaym)**2)) * refRcTE
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res + termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res + termadd
  endif

endfunction iabvC1TEderiv


!---------------------------------------------------------
! derivative of integrand C1 TE with respect to layer ilaym for receivers below source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iblwC1TEderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRcTE ! dRa / depsilon_m
  complex(kind=real64)            :: refRcTE  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TE(kappa,pvert)

  !  derivative of total reflection response
  drefRcTE = get_drefR_blw(-1._real64,-1._real64,set_trans_blwTE, fact_exp_teiso,fact_isote,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = - (sqrt(pvert(ilayrec)*pvert(ilaysrc)) / dmu0) * drefRcTE

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRcTE = iblwC1TE(kappa) !contains  aminus
    termadd = (dmu0 / (4.*pvert(ilaym)**2)) * refRcTE
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res + termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res + termadd
  endif

endfunction iblwC1TEderiv


!---------------------------------------------------------
! derivative of integrand C0 TM with respect to layer ilaym for receivers above source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iabvC0TMderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iabvC1TMderiv(kappa)
endfunction iabvC0TMderiv

complex(kind=real64) function iabvC0TMderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iabvC1TMderivh(kappa)
endfunction iabvC0TMderivh

complex(kind=real64) function iabvC0TMderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iabvC1TMderivv(kappa)
endfunction iabvC0TMderivv

!---------------------------------------------------------
! derivative of integrand C0 TM with respect to layer ilaym for receivers below source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iblwC0TMderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iblwC1TMderiv(kappa)
endfunction iblwC0TMderiv

complex(kind=real64) function iblwC0TMderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iblwC1TMderivh(kappa)
endfunction iblwC0TMderivh

complex(kind=real64) function iblwC0TMderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa * iblwC1TMderivv(kappa)
endfunction iblwC0TMderivv


!---------------------------------------------------------
! derivative of integrand C1 TM with respect to layer ilaym for receivers above source
! RS 2010-2011
!---------------------------------------------------------
complex(kind=real64) function iabvC1TMderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRcTM ! dRc / depsilon_m
  complex(kind=real64)            :: refRcTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRcTM = get_drefR_abv(-1._real64,-1._real64,set_trans_abvTM, fact_exp_tmiso,fact_isotm,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = - sqrt(epsh(ilayrec)*epsh(ilaysrc) / (pvert(ilayrec)*pvert(ilaysrc))) * drefRcTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRcTM = iabvC1TM(kappa) !contains a minus
    termadd = (1./(2.*epsh(ilaym)) - dmu0/(4.*pvert(ilaym)**2)) * refRcTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res + termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res + termadd
  endif

endfunction iabvC1TMderiv


!---------------------------------------------------------
! derivative of integrand C1 TM with respect to layer ilaym for receivers above source
! VTI for epsh
! RS 2010-2011
!---------------------------------------------------------
complex(kind=real64) function iabvC1TMderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRcTM ! dRc / depsilon_m
  complex(kind=real64)            :: refRcTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m
  real(kind=real64)   :: phsq   !horiz. slowness

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRcTM = get_drefR_abv(-1._real64,-1._real64,set_trans_abvTM, fact_exp_tmvtih,fact_vtitmh,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = - sqrt(epsh(ilayrec)*epsh(ilaysrc) / (pvert(ilayrec)*pvert(ilaysrc))) * drefRcTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRcTM = iabvC1TM(kappa) !contains a minus
    phsq = (kappa / omega)**2
    termadd = (1./(2.*epsh(ilaym)) - (dmu0 - phsq/epsv(ilaym))/(4.*pvert(ilaym)**2)) * refRcTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res + termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res + termadd
  endif

endfunction iabvC1TMderivh


!---------------------------------------------------------
! derivative of integrand C1 TM with respect to layer ilaym for receivers above source
! VTI for epsv
! RS 2010-2011
!---------------------------------------------------------
complex(kind=real64) function iabvC1TMderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRcTM ! dRc / depsilon_m
  complex(kind=real64)            :: refRcTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m
  real(kind=real64)   :: phsq   !horiz. slowness

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRcTM = get_drefR_abv(-1._real64,-1._real64,set_trans_abvTM, fact_exp_tmvtiv,fact_vtitmv,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = - sqrt(epsh(ilayrec)*epsh(ilaysrc) / (pvert(ilayrec)*pvert(ilaysrc))) * drefRcTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRcTM = iabvC1TM(kappa) !contains a minus
    phsq = (kappa / omega)**2
    termadd = (- epsh(ilaym)*phsq / (4.* (epsv(ilaym) * pvert(ilaym))**2)) * refRcTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res + termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res + termadd
  endif

endfunction iabvC1TMderivv


!---------------------------------------------------------
! derivative of integrand C1 TM with respect to layer ilaym for receivers below source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iblwC1TMderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRcTM ! dRc / depsilon_m
  complex(kind=real64)            :: refRcTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m
  !!!TEST
!!$  complex(kind=real64)  :: dT_depsm_num
!!$  real(kind=real64)     :: errnum
!!$  real(kind=real64)     :: resdif


  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRcTM = get_drefR_blw(-1._real64,-1._real64,set_trans_blwTM, fact_exp_tmiso,fact_isotm,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = - sqrt(epsh(ilayrec)*epsh(ilaysrc) / (pvert(ilayrec)*pvert(ilaysrc))) * drefRcTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRcTM = iblwC1TM(kappa) !contains a minus
    termadd = (1./(2.*epsh(ilaym)) - dmu0/(4.*pvert(ilaym)**2)) * refRcTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res + termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res + termadd
  endif

  !test: numerical derivative d_someting / depsilon_m
!!$  kappatmp = kappa
!!$  dT_depsm_num = dfridr_cplx(R_num,epsh(ilaym),epsh(ilaym)*0.5,errnum)
!!$  resdif = abs((res - dT_depsm_num) / dT_depsm_num)

endfunction iblwC1TMderiv


!---------------------------------------------------------
! derivative of integrand C1 TM with respect to layer ilaym for receivers below source
! VTI for epsh
! RS 2011
!---------------------------------------------------------
complex(kind=real64) function iblwC1TMderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRcTM ! dRc / depsilon_m
  complex(kind=real64)            :: refRcTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m
  real(kind=real64)   :: phsq   !horiz. slowness


  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRcTM = get_drefR_blw(-1._real64,-1._real64,set_trans_blwTM, fact_exp_tmvtih,fact_vtitmh,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = - sqrt(epsh(ilayrec)*epsh(ilaysrc) / (pvert(ilayrec)*pvert(ilaysrc))) * drefRcTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRcTM = iblwC1TM(kappa) !contains a minus
    phsq = (kappa / omega)**2
    termadd = (1./(2.*epsh(ilaym)) - (dmu0 - phsq/epsv(ilaym))/(4.*pvert(ilaym)**2)) * refRcTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res + termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res + termadd
  endif

endfunction iblwC1TMderivh


!---------------------------------------------------------
! derivative of integrand C1 TM with respect to layer ilaym for receivers below source
! VTI for epsv
! RS 2011
!---------------------------------------------------------
complex(kind=real64) function iblwC1TMderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: drefRcTM ! dRc / depsilon_m
  complex(kind=real64)            :: refRcTM  !total reflection response
  complex(kind=real64)            :: termadd  !additional term if source and/or receiver is in layer m
  real(kind=real64)   :: phsq   !horiz. slowness

  !set vertical wavenumbers and interface reflection coeff.
  call set_fresnelcoef_TM(kappa,pvert)

  !  derivative of total reflection response
  drefRcTM = get_drefR_blw(-1._real64,-1._real64,set_trans_blwTM, fact_exp_tmvtiv,fact_vtitmv,kappa)

  !entire integrand (except for Bessel function, which is handled by zhankl)
  res = - sqrt(epsh(ilayrec)*epsh(ilaysrc) / (pvert(ilayrec)*pvert(ilaysrc))) * drefRcTM

  !special contributions if layer m is the source or(and) receiver layer
  if ((ilaym.eq.ilayrec) .or. (ilaym.eq.ilaysrc)) then
    refRcTM = iblwC1TM(kappa) !contains a minus
    phsq = (kappa / omega)**2
    termadd = (- epsh(ilaym)*phsq / (4.* (epsv(ilaym) * pvert(ilaym))**2)) * refRcTM
    !layer m is receiver layer
    if (ilaym .eq. ilayrec) res = res + termadd
    !layer m is source layer
    if (ilaym .eq. ilaysrc) res = res + termadd
  endif

endfunction iblwC1TMderivv


!---------------------------------------------------------
! derivative of integrand Cz1 TM with respect to layer ilaym for receivers above source - for Ez
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iabvCz1TMderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iabvC1TMderiv(kappa)
endfunction iabvCz1TMderiv

complex(kind=real64) function iabvCz1TMderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iabvC1TMderivh(kappa)
endfunction iabvCz1TMderivh

complex(kind=real64) function iabvCz1TMderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iabvC1TMderivv(kappa)
endfunction iabvCz1TMderivv

!---------------------------------------------------------
! derivative of integrand Cz1 TM with respect to layer ilaym for receivers below source - for Ez
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iblwCz1TMderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iblwC1TMderiv(kappa)
endfunction iblwCz1TMderiv

complex(kind=real64) function iblwCz1TMderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iblwC1TMderivh(kappa)
endfunction iblwCz1TMderivh

complex(kind=real64) function iblwCz1TMderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iblwC1TMderivv(kappa)
endfunction iblwCz1TMderivv

!---------------------------------------------------------
! derivative of integrand B1 TM for VED with respect to layer ilaym for receivers above source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iabvB1TMvedderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iabvB1TMderiv(kappa)
endfunction iabvB1TMvedderiv

complex(kind=real64) function iabvB1TMvedderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iabvB1TMderivh(kappa)
endfunction iabvB1TMvedderivh

complex(kind=real64) function iabvB1TMvedderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iabvB1TMderivv(kappa)
endfunction iabvB1TMvedderivv

!---------------------------------------------------------
! derivative of integrand B1 TM for VED with respect to layer ilaym for receivers below source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iblwB1TMvedderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iblwB1TMderiv(kappa)
endfunction iblwB1TMvedderiv

complex(kind=real64) function iblwB1TMvedderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iblwB1TMderivh(kappa)
endfunction iblwB1TMvedderivh

complex(kind=real64) function iblwB1TMvedderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iblwB1TMderivv(kappa)
endfunction iblwB1TMvedderivv

!---------------------------------------------------------
! derivative of integrand C0 TM for VED with respect to layer ilaym for receivers above source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iabvC0TMvedderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**3 * iabvC1TMderiv(kappa)
endfunction iabvC0TMvedderiv

complex(kind=real64) function iabvC0TMvedderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**3 * iabvC1TMderivh(kappa)
endfunction iabvC0TMvedderivh

complex(kind=real64) function iabvC0TMvedderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**3 * iabvC1TMderivv(kappa)
endfunction iabvC0TMvedderivv

!---------------------------------------------------------
! derivative of integrand C0 TM for VED with respect to layer ilaym for receivers below source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iblwC0TMvedderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**3 * iblwC1TMderiv(kappa)
endfunction iblwC0TMvedderiv

complex(kind=real64) function iblwC0TMvedderivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**3 * iblwC1TMderivh(kappa)
endfunction iblwC0TMvedderivh

complex(kind=real64) function iblwC0TMvedderivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**3 * iblwC1TMderivv(kappa)
endfunction iblwC0TMvedderivv


!---------------------------------------------------------
! derivative of integrand D1 TE for VMD with respect to layer ilaym for receivers above source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iabvD1TEvmdderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iabvD1TEderiv(kappa)
endfunction iabvD1TEvmdderiv

!---------------------------------------------------------
! derivative of integrand D1 TE for VMD with respect to layer ilaym for receivers below source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iblwD1TEvmdderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**2 * iblwD1TEderiv(kappa)
endfunction iblwD1TEvmdderiv

!---------------------------------------------------------
! derivative of integrand A0 TE for VMD with respect to layer ilaym for receivers above source - for Hz
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iabvA0TEvmdderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**3 * iabvA1TEderiv(kappa)
endfunction iabvA0TEvmdderiv

!---------------------------------------------------------
! derivative of integrand A0 TE for VMD with respect to layer ilaym for receivers below source - for Hz
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function iblwA0TEvmdderiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = kappa**3 * iblwA1TEderiv(kappa)
endfunction iblwA0TEvmdderiv


!---------------------------------------------------------
! derivative of integrand (A0 TE) of first integral for Ex for long wire sources
!   with respect to layer ilaym for receivers above source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function i1abvExwirederiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = - kappa * iabvA1TEderiv(kappa)
endfunction i1abvExwirederiv

!---------------------------------------------------------
! derivative of integrand (A0 TE) of first integral for Ex for long wire sources
!   with respect to layer ilaym for receivers below source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function i1blwExwirederiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  res = - kappa * iblwA1TEderiv(kappa)
endfunction i1blwExwirederiv


!---------------------------------------------------------
! derivative of integrand for Ex for wire end points with respect to layer ilaym, receivers above source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function i2abvExwirederiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: resTE, resTM  !results from TE and TM parts

  !TE part
  resTE = iabvA1TEderiv(kappa)
  !TM part
  resTM = iabvA1TMderiv(kappa)

  res = resTE - resTM
endfunction i2abvExwirederiv

complex(kind=real64) function i2abvExwirederivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: resTE, resTM  !results from TE and TM parts

  !TE part
  resTE = iabvA1TEderiv(kappa)
  !TM part
  resTM = iabvA1TMderivh(kappa)

  res = resTE - resTM
endfunction i2abvExwirederivh

complex(kind=real64) function i2abvExwirederivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !TE part is zero
  !TM part
  res = - iabvA1TMderivv(kappa)
endfunction i2abvExwirederivv


!---------------------------------------------------------
! derivative of integrand for Ex for wire end points with respect to layer ilaym, receivers below source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function i2blwExwirederiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: resTE, resTM  !results from TE and TM parts

  !TE part
  resTE = iblwA1TEderiv(kappa)
  !TM part
  resTM = iblwA1TMderiv(kappa)

  res = resTE - resTM
endfunction i2blwExwirederiv

complex(kind=real64) function i2blwExwirederivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: resTE, resTM  !results from TE and TM parts

  !TE part
  resTE = iblwA1TEderiv(kappa)
  !TM part
  resTM = iblwA1TMderivh(kappa)

  res = resTE - resTM
endfunction i2blwExwirederivh

complex(kind=real64) function i2blwExwirederivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !TE part = 0
  !TM part
  res = - iblwA1TMderivv(kappa)
endfunction i2blwExwirederivv


!---------------------------------------------------------
! derivative of integrand for Hx and Hy for wire end points with respect to layer ilaym, receivers above source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function idabvHxwirederiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: resTE, resTM  !results from TE and TM parts
  !!!TEST
  !complex(kind=real64)  :: dT_depsm_num
  !real(kind=real64)     :: errnum
  !real(kind=real64)     :: resdif

  !TE part
  resTE = iabvD1TEderiv(kappa)
  !TM part
  resTM = iabvD1TMderiv(kappa)

  res = resTE - resTM

  !test: numerical derivative d_someting / depsilon_m
  !kappatmp = kappa
  !dT_depsm_num = dfridr_cplx(R_num,epsh(ilaym),epsh(ilaym)*0.5,errnum)
  !resdif = abs((res - dT_depsm_num) / dT_depsm_num)
endfunction idabvHxwirederiv

complex(kind=real64) function idabvHxwirederivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: resTE, resTM  !results from TE and TM parts

  !TE part
  resTE = iabvD1TEderiv(kappa)
  !TM part
  resTM = iabvD1TMderivh(kappa)

  res = resTE - resTM
endfunction idabvHxwirederivh

complex(kind=real64) function idabvHxwirederivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !TE part = 0
  !TM part
  res = - iabvD1TMderivv(kappa)
endfunction idabvHxwirederivv

!---------------------------------------------------------
! derivative of integrand for Hx and Hy for wire end points with respect to layer ilaym, receivers below source
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function idblwHxwirederiv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: resTE, resTM  !results from TE and TM parts

  !TE part
  resTE = iblwD1TEderiv(kappa)
  !TM part
  resTM = iblwD1TMderiv(kappa)

  res = resTE - resTM
endfunction idblwHxwirederiv

complex(kind=real64) function idblwHxwirederivh(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !internal variables
  complex(kind=real64)            :: resTE, resTM  !results from TE and TM parts

  !TE part
  resTE = iblwD1TEderiv(kappa)
  !TM part
  resTM = iblwD1TMderivh(kappa)

  res = resTE - resTM
endfunction idblwHxwirederivh

complex(kind=real64) function idblwHxwirederivv(kappa) result(res)
  implicit none
  !external variable
  real(kind=real64)               :: kappa    !horizontal wavenumber

  !TE part = 0
  !TM part
  res = - iblwD1TMderivv(kappa)
endfunction idblwHxwirederivv



!!!TEST
complex(kind=real64) function R_num(modpar)

  implicit none

  complex(kind=real64)  :: modpar
  real(kind=real64)     :: kappa

  complex(kind=real64)  :: eps_safe,epsmu_safe,omsq_epsmu_safe,rnum1


  !remember original model
  eps_safe = epsh(ilaym)
  epsmu_safe = epsmuv(ilaym)
  omsq_epsmu_safe = omsq_epsmuv(ilaym)

  !replace epsilon for layer ilaym
  epsh(ilaym) = modpar
  epsmuv(ilaym) = epsh(ilaym)*dmu0
  omsq_epsmuv(ilaym) = omegasq*epsmuv(ilaym)
  branchpt(1:nlay) = real(-jomega*epsmuv)
  branchpt(nlay+1:2*nlay) = real(-jomega*epsmuh)
  call sort_dbl(2*nlay,branchpt)
  !retrieve current wavenumber
  kappa = kappatmp

  !recompute refl. coeff. - only need this when checking parts of integrals
  !call set_fresnelcoef_TE(kappa,pvert)
  !transmission factors through all layers for receivers below or above source
  !call set_trans_blwTE()

!!!receivers below source
  !rnum1 = zeroc
  !R_num = trans_response_dn(ilaysrc,ilayrec,trans_above_rec,ilaysrc,trans_above_rec(ilayrec),rnum1)
  !R_num = refl_response_dn(ilaysrc,ilayrec-1,trans_above_rec,ilaysrc,zeroc) !used inside dtrans...
  !R_num = refl_response_up(ilaysrc+1,ilayrec,trans_above_rec,ilaysrc,zeroc)
  !R_num = refl_response_up(2,ilaysrc,trans_above_src,2,zeroc)
  !R_num = refl_response_dn(ilaysrc,nlay-1,trans_below_src,ilaysrc,zeroc)
  !R_num = refl_response_dn(ilayrec,nlay-1,trans_below_rec,ilayrec,zeroc)

!!!receivers above source
  !R_num = refl_response_dn(ilayrec,nlay-1,trans_below_rec,ilayrec,zeroc) !Ru
  !R_num = refl_response_up(2,ilayrec,trans_above_rec,2,zeroc)  !Ra
  !rnum1 = zeroc
  !R_num = trans_response_up(ilayrec,ilaysrc,trans_below_rec,ilayrec,trans_below_rec(ilayrec),rnum1)
  !R_num = refl_response_up(ilayrec+1,ilaysrc,trans_below_rec,ilayrec,zeroc) !used inside dtrans...


  R_num = idabvHxwire(kappa)


  !back to original
  epsh(ilaym) = eps_safe
  epsmuv(ilaym) = epsmu_safe
  omsq_epsmuv(ilaym) = omsq_epsmu_safe
  branchpt(1:nlay) = real(-jomega*epsmuv)
  branchpt(nlay+1:2*nlay) = real(-jomega*epsmuh)
  call sort_dbl(2*nlay,branchpt)

  !recompute refl. coeff.
  !call set_fresnelcoef_TE(kappa,pvert)
  !transmission factors through all layers for receivers below or above source
  !call set_trans_blwTE()

endfunction R_num


!---------------------------------------------------------
! derivative drefR(a,d,b,c) / depsilon_m for receivers above source, TE or TM
! with respect to layer ilaym
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function get_drefR_abv(sgn1,sgn2,transfunc, func_dexp, func_coef,kappa) result(drefR)

  implicit none

  !external variables
  real(kind=real64),intent(in)   :: sgn1,sgn2 !+1 or -1
  complex(kind=real64),external  :: func_dexp, func_coef !functions for little terms that differ between TE TM ISO VTI epsh epsv
  real(kind=real64)              :: kappa     !horiz. wavenumber
  interface
    subroutine transfunc() !set_trans_abvTE or set_trans_abvTM
    endsubroutine transfunc
  end interface  

  !internal variables
  complex(kind=real64)  :: rupsrc,rdnsrc,ra,ru,tu  !contributions to total response as for forward computation
  !there are 5 contributions (using chain rule for derivatives): drefRa/dTd * dTd/depsilon_m etc
  complex(kind=real64)  :: drefR_dTu, drefR_dRa, drefR_dRu, drefR_dRupsrc, drefR_dRdnsrc
  complex(kind=real64)  :: dTu_depsm, dRa_depsm, dRu_depsm, dRupsrc_depsm, dRdnsrc_depsm
  complex(kind=real64)  :: rsrcterm,ruaterm,ru_up

  !transmission factors through all layers for receivers above source
  call transfunc()

  !get all the components as for reflection response
  rdnsrc = refl_response_dn(ilaysrc,nlay-1,trans_below_src,ilaysrc,zeroc)
  rupsrc = refl_response_up(2,ilaysrc,trans_above_src,2,zeroc)
  ra = refl_response_up(2,ilayrec,trans_above_rec,2,zeroc)
  ru = refl_response_dn(ilayrec,ilaysrc-1,trans_below_rec,ilayrec,zeroc)
  ru_up = zeroc
  tu = trans_response_up(ilayrec,ilaysrc,trans_below_rec,ilayrec,trans_below_rec(ilayrec),ru_up)

  !get factors drefR/dTu etc.
  rsrcterm = 1.-rdnsrc*rupsrc
  ruaterm = 1.-ru*ra
  drefR_dTu = (1.+sgn1*ra) * (rdnsrc+sgn2) / (ruaterm * rsrcterm)
  drefR_dRa = tu*(rdnsrc+sgn2)*(ru+sgn1) / ( rsrcterm * ruaterm**2)
  drefR_dRu = tu*(1.+sgn1*ra)*ra*(rdnsrc+sgn2) / (ruaterm**2 * rsrcterm)
  drefR_dRdnsrc = tu*(1.+sgn1*ra)*(1.+sgn2*rupsrc) / (ruaterm * rsrcterm**2)
  drefR_dRupsrc = tu*(1.+sgn1*ra)*(rdnsrc+sgn2)*rdnsrc / (ruaterm * rsrcterm**2)

  !derivative dTu / depsilon_m - from receiver to source
  dTu_depsm = dtrans_response_up(ilayrec,ilaysrc,trans_below_rec,dz_below_rec,ilayrec,ilaysrc,kappa, func_dexp,func_coef)
  !derivative dRa / depsilon_m - from top to receiver
  dRa_depsm = drefl_response_up(1,ilayrec,trans_above_rec,dz_above_rec,2,ilayrec,kappa, func_dexp,func_coef)
  !derivative dRu / depsilon_m - recursion from source to receiver
  dRu_depsm = drefl_response_dn(ilayrec,ilaysrc,trans_below_rec,dz_below_rec,ilayrec,ilaysrc,kappa, func_dexp,func_coef)
  !derivative dRupsrc / depsilon_m - from top to source
  dRupsrc_depsm = drefl_response_up(1,ilaysrc,trans_above_src,dz_above_src,2,ilaysrc,kappa, func_dexp,func_coef)
  !derivative dRdnsrc / depsilon_m - from bottom to source
  dRdnsrc_depsm = drefl_response_dn(ilaysrc,nlay,trans_below_src,dz_below_src,ilaysrc,nlay-1,kappa, func_dexp,func_coef)

  !final result
  drefR = drefR_dTu*dTu_depsm + drefR_dRa*dRa_depsm + drefR_dRu*dRu_depsm + &
             drefR_dRupsrc*dRupsrc_depsm + drefR_dRdnsrc*dRdnsrc_depsm

endfunction get_drefR_abv


!---------------------------------------------------------
! derivative drefR(a,d,b,c) / depsilon_m for receivers below source, TE or TM
! with respect to layer ilaym
! RS 2010
!---------------------------------------------------------
complex(kind=real64) function get_drefR_blw(sgn1,sgn2,transfunc, func_dexp,func_coef,kappa) result(drefR)

  implicit none

  !external variables
  real(kind=real64),intent(in)   :: sgn1,sgn2 !+1 or -1
  !functions for computing little terms that differ between TE TM iso vti epsh epsv
  complex(kind=real64),external  :: func_dexp, func_coef
  real(kind=real64)              :: kappa     !horiz. wavenumber
  interface
    subroutine transfunc() !set_trans_blwTE or set_trans_blwTM
    endsubroutine transfunc
  end interface  

  !internal variables
  complex(kind=real64)  :: rupsrc,rdnsrc,rb,rd,td  !contributions to total response as for forward computation
  !there are 5 contributions (using chain rule for derivatives): drefR/dTd * dTd/depsilon_m etc
  complex(kind=real64)  :: drefR_dTd, drefR_dRb, drefR_dRd, drefR_dRupsrc, drefR_dRdnsrc
  complex(kind=real64)  :: dTd_depsm, dRb_depsm, dRd_depsm, dRupsrc_depsm, dRdnsrc_depsm
  complex(kind=real64)  :: rsrcterm,rdbterm,rd_dn


  !transmission factors through all layers for receivers below source
  call transfunc()

  !get all the components as for reflection response
  rdnsrc = refl_response_dn(ilaysrc,nlay-1,trans_below_src,ilaysrc,zeroc)
  rupsrc = refl_response_up(2,ilaysrc,trans_above_src,2,zeroc)
  rd = refl_response_up(ilaysrc+1,ilayrec,trans_above_rec,ilaysrc,zeroc)
  rb = refl_response_dn(ilayrec,nlay-1,trans_below_rec,ilayrec,zeroc)
  rd_dn = zeroc
  td = trans_response_dn(ilaysrc,ilayrec,trans_above_rec,ilaysrc,trans_above_rec(ilayrec),rd_dn)

  !get factors drefRa/dTd - these are simple ...
  rsrcterm = 1.-rupsrc*rdnsrc
  rdbterm = 1.-rd*rb
  drefR_dTd = (rb+sgn1) * (1.+sgn2*rupsrc) / (rdbterm * rsrcterm)
  drefR_dRb = td*(1.+sgn2*rupsrc)*(1.+sgn1*rd) / (rsrcterm * rdbterm**2 )
  drefR_dRd = td*(rb+sgn1)*rb*(1.+sgn2*rupsrc) / (rdbterm**2 * rsrcterm)
  drefR_dRupsrc = td*(rb+sgn1)*(rdnsrc+sgn2) / ( rdbterm * rsrcterm**2 )
  drefR_dRdnsrc = td*(rb+sgn1)*(1.+sgn2*rupsrc)*rupsrc / (rdbterm * rsrcterm**2)

  !derivative dTd / depsilon_m - from receiver to source
  dTd_depsm = dtrans_response_dn(ilaysrc,ilayrec,trans_above_rec,dz_above_rec,ilaysrc,ilayrec,kappa, func_dexp,func_coef)
  !derivative dRb / depsilon_m - from bottom to receiver
  dRb_depsm = drefl_response_dn(ilayrec,nlay,trans_below_rec,dz_below_rec,ilayrec,nlay-1,kappa, func_dexp,func_coef)
  !derivative dRd / depsilon_m - recursion from source to receiver
  dRd_depsm = drefl_response_up(ilaysrc,ilayrec,trans_above_rec,dz_above_rec,ilaysrc,ilayrec,kappa, func_dexp,func_coef)
  !derivative dRupsrc / depsilon_m - from top to source
  dRupsrc_depsm = drefl_response_up(1,ilaysrc,trans_above_src,dz_above_src,2,ilaysrc,kappa, func_dexp,func_coef)
  !derivative dRdnsrc / depsilon_m - from bottom to source
  dRdnsrc_depsm = drefl_response_dn(ilaysrc,nlay,trans_below_src,dz_below_src,ilaysrc,nlay-1,kappa, func_dexp,func_coef)

  !final result
  drefR = drefR_dTd*dTd_depsm + drefR_dRb*dRb_depsm + drefR_dRd*dRd_depsm + &
             drefR_dRupsrc*dRupsrc_depsm + drefR_dRdnsrc*dRdnsrc_depsm

endfunction get_drefR_blw


!---------------------------------------------------------
! derivative of "downgoing" transmission response between source and receiver
! TE or TM, generic for iso, vti epsh, epsv
! deriv. with respect to layer ilaym, i.e. dTd_TM / depsilon_m
! RS 2011
!---------------------------------------------------------
complex(kind=real64) function dtrans_response_dn(ilaytop,ilaybot,trans_abv,dz_abv,startidx,endidx,kappa, &
     func_dexp,func_coef) result(dtransT)

  implicit none

  !external variables
  integer(kind=int32)   :: ilaytop,ilaybot  !layer indices: start recursion from ilaybot, go to ilaytop
  integer(kind=int32),intent(in)      :: startidx,endidx  !start and end index of transm. coeff. vector
  complex(kind=real64),dimension(startidx:endidx),intent(in) :: trans_abv !transm. coeff. either with respect to src or rec
  real(kind=real64),dimension(startidx:endidx),intent(in)    :: dz_abv    !layer thicknesses
  real(kind=real64),intent(in)    :: kappa    !horiz. wavenumber
  complex(kind=real64),external   :: func_dexp,func_coef  !functions, different for TE TM iso vti epsh epsv

  !internal variables
  integer(kind=int32)   :: ilay    !layer counter
  complex(kind=real64)  :: refT,drefT   !reflection response and its derivative TE or TM
  complex(kind=real64)  :: transT !transmission response TE or TM
  integer(kind=int32)   :: ilaym1  !layer index for layer just above ilaym
  complex(kind=real64)  :: term1,term2,term3,denom,fact_coef,fact_dexp


  !initialize - need this?
  drefT = 0._real64
  refT = 0._real64
  dtransT = 0._real64
  transT = 0._real64
  ilaym1 = ilaym + 1

  !term stays zero if layer m is outside range in which this recursion applies
  if ((ilaym.lt.ilaytop) .or. (ilaym.gt.ilaybot)) return


  !derivative for layer ilaym
  !--> special (simplified) case if this is at the bottom of the layer stack considered (i.e. ilaym = receiver layer)
  fact_dexp = func_dexp(dz_abv(ilaym),kappa) / 2.
  if (ilaym .eq. ilaybot) then
    transT = trans_abv(ilaym)
    dtransT = trans_abv(ilaym) * fact_dexp
    drefT = 0._real64
    refT = 0._real64
  else
    !transmission response to layer just below ilaym, refTM is also returned
    refT = 0._real64
    transT = trans_response_dn(ilaym1,ilaybot,trans_abv,ilaysrc,trans_abv(ilaybot),refT)
    !contribution from deriv. of refl. and transm. coeff.
    fact_coef = func_coef(kappa)
    denom = (1.-rup(ilaym)*refT)
    term1 = - trans_abv(ilaym) * transT * tup(ilaym) * fact_coef / denom
    term2 = refT * tup(ilaym)**2 / denom
    dtransT = term1 * (rdn(ilaym) + term2)

    !reflection response from layers up to the one below ilaym
    !refTM = refl_response_dn(ilaym1,ilaybot-1,trans_abv,startidx,zeroc)
    !derivative of reflection response:
    !terms containing refTM
    term1 = trans_abv(ilaym)**2 * tup(ilaym)**2 * fact_coef
    term2 = 2. * refT * rdn(ilaym) / denom
    term3 = (refT * tup(ilaym) / denom)**2
    drefT = term1 * (1. - term2 - term3)

    !term from derivative of exponential factor (contains pvert_m) --> contains refTM for layer ilaym --> next iteration
    !reflection and transmission response for layers up to layer m
    transT = trans_response_dn(ilaym,ilaym1,trans_abv,ilaysrc,transT,refT)
    !contribution from derivative of transmission factor
    dtransT = dtransT + transT * fact_dexp
    !last term for deriv. of refl. response, needs refl. response to layer m
    drefT = drefT + refT * 2. * fact_dexp
  endif

  !derivative for the layer above ilaym
  !special case again since the interface refl./transm. coeff. still depend on eps(ilaym)
  if (ilaym.gt.ilaytop) then !nothing more to do if ilaym is top layer of the stack considered
    !FIRST: update dtransTM using results from previous layer
    !we already have transm. and refl. responses for layers up to ilaym, so no update needed here
    denom = 1.-rup(ilaym-1)*refT
    dtransT = (trans_abv(ilaym-1)*tdn(ilaym-1)/denom) * &
                (dtransT + (transT*rup(ilaym-1)/denom)*drefT )

    !now add special contribution from derivatives of refl./transm. coeff.
    fact_coef = func_coef(kappa)
    term1 = trans_abv(ilaym-1) * transT * tup(ilaym-1) * fact_coef / denom
    term2 = refT * tup(ilaym-1)**2 / denom
    dtransT = dtransT + term1 * (rdn(ilaym-1) + term2)

    !derivative of reflection response for this layer
    !FIRST add the term containing deriv. of refl. up to layer m
    drefT = drefT * (trans_abv(ilaym-1) * tup(ilaym-1) / denom)**2

    !now add special contribution arising from interface refl./transm. coeff. still containing epsilon_m
    term1 = - trans_abv(ilaym-1)**2 * tup(ilaym-1)**2 * fact_coef
    term2 = 2. * refT * rdn(ilaym-1) / denom
    term3 = (refT * tup(ilaym-1) / denom)**2
    drefT = drefT + term1 * (1. - term2 - term3)
  endif


  !iteration for remaining layers, from 2 layers above ilaym to top of stack considered
  do ilay=ilaym-2,ilaytop,-1
    !update transm. response, up to the layer below ilay, refTM is also updated
    transT = trans_response_dn(ilay+1,ilay+2,trans_abv,ilaysrc,transT,refT)

    !update derivative of transm. response
    denom = 1.-rup(ilay)*refT
    dtransT = (trans_abv(ilay)*tdn(ilay)/denom) * &
                (dtransT + (transT*rup(ilay)/denom)*drefT )

    !update derivative of refl. response
    drefT = drefT * (trans_abv(ilay) * tup(ilay) / denom)**2
  enddo

endfunction dtrans_response_dn


!---------------------------------------------------------
! derivative of "upgoing" TM transmission response between source and receiver
! TE or TM, generic for iso, vti epsh, epsv
! deriv. with respect to layer ilaym, i.e. dTd_TM / depsilon_m
! RS 2011
!---------------------------------------------------------
complex(kind=real64) function dtrans_response_up(ilaytop,ilaybot,trans_blw,dz_blw,startidx,endidx,kappa, &
       func_dexp,func_coef) result(dtransT)

  implicit none

  !external variables
  integer(kind=int32)   :: ilaytop,ilaybot  !layer indices: start recursion from ilaybot, go to ilaytop
  integer(kind=int32),intent(in)      :: startidx,endidx  !start and end index of transm. coeff. vector
  complex(kind=real64),dimension(startidx:endidx),intent(in) :: trans_blw !transm. coeff. either with respect to src or rec
  real(kind=real64),dimension(startidx:endidx),intent(in)    :: dz_blw    !layer thicknesses
  real(kind=real64),intent(in)    :: kappa    !horiz. wavenumber
  complex(kind=real64),external   :: func_dexp,func_coef  !functions, different for TE TM iso vti epsh epsv

  !internal variables
  integer(kind=int32)   :: ilay    !layer counter
  complex(kind=real64)  :: refT,drefT   !reflection response and its derivative
  complex(kind=real64)  :: transT !transmission response
  integer(kind=int32)   :: ilaym1  !layer index for layer just above ilaym
  complex(kind=real64)  :: term1,term2,term3,denom,fact_coef,fact_dexp


  !initialize - need this?
  drefT = 0._real64
  refT = 0._real64
  dtransT = 0._real64
  transT = 0._real64
  ilaym1 = ilaym - 1

  !term stays zero if layer m is outside range in which this recursion applies
  if ((ilaym.lt.ilaytop) .or. (ilaym.gt.ilaybot)) return


  !derivative for layer ilaym
  !--> special (simplified) case if this is at the top of the layer stack considered (i.e. ilaym = receiver layer)
  fact_dexp = func_dexp(dz_blw(ilaym),kappa) / 2.
  if (ilaym .eq. ilaytop) then
    transT = trans_blw(ilaym)
    dtransT = trans_blw(ilaym) * fact_dexp
    drefT = 0._real64
    refT = 0._real64
  else
    !transmission response to layer just above ilaym, refTM is also returned
    refT = 0._real64
    transT = trans_response_up(ilaytop,ilaym1,trans_blw,startidx,trans_below_rec(ilayrec),refT)

    !contribution from deriv. of refl. and transm. coeff.
    fact_coef = func_coef(kappa)
    denom = 1.-rdn(ilaym1)*refT
    term1 = - trans_blw(ilaym) * transT * tup(ilaym1) * fact_coef / denom
    term2 = tup(ilaym1)**2 * refT / denom
    dtransT = term1 * (rup(ilaym1) + term2)

    !derivative of reflection response:
    !terms containing refTM
    term1 = trans_blw(ilaym)**2 * tup(ilaym1)**2 *  fact_coef
    term2 = 2. * refT * rup(ilaym1) / denom
    term3 = (tdn(ilaym1) * refT / denom)**2
    drefT = term1 * (1. - term2 - term3)

    !term from derivative of exponential factor (contains pvert_m) --> contains refTM for layer ilaym --> next iteration
    !reflection and transmission response for layers down to layer m
    transT = trans_response_up(ilaym1,ilaym,trans_blw,startidx,transT,refT)
    !contribution from derivative of transmission factor
    dtransT = dtransT + transT * fact_dexp
    !last term for deriv. of refl. response, needs refl. response to layer m
    drefT = drefT + refT * 2. * fact_dexp
  endif

  !derivative for the layer below ilaym
  !special case again since the interface refl./transm. coeff. still depend on eps(ilaym)
  if (ilaym.lt.ilaybot) then !nothing more to do if ilaym is bottom layer of the stack considered
    !FIRST: update dtransTM using results from previous layer
    !we already have transm. and refl. responses for layers down to ilaym, so no update needed here
    denom = 1.-rdn(ilaym)*refT
    dtransT = (trans_blw(ilaym+1)*tup(ilaym)/denom) * &
                (dtransT + (transT*rdn(ilaym)/denom)*drefT )

    !now add special contribution from derivatives of refl./transm. coeff.
    fact_coef = func_coef(kappa)
    
    term1 = trans_blw(ilaym+1) * transT * tup(ilaym) * fact_coef / denom
    term2 = refT * tup(ilaym)**2 / denom
    dtransT = dtransT + term1 * (rup(ilaym) + term2)

    !derivative of reflection response for this layer
    !FIRST add the term containing deriv. of refl. up to layer m
    drefT = drefT * (trans_blw(ilaym+1) * tdn(ilaym) / denom)**2

    !now add special contribution arising from interface refl./transm. coeff. still containing epsilon_m
    term1 = - (trans_blw(ilaym+1) * tup(ilaym))**2 * fact_coef
    term2 = 2.* refT *rup(ilaym) / denom
    term3 = (tdn(ilaym) * refT / denom)**2
    drefT = drefT + term1 * (1. - term2 - term3)
  endif


  !iteration for remaining layers, from 2 layers above ilaym to top of stack considered
  do ilay=ilaym+2,ilaybot
    !update transm. response, down to the layer below ilay, refTM is also updated
    transT = trans_response_up(ilay-2,ilay-1,trans_blw,startidx,transT,refT)

    !update derivative of transm. response
    denom = 1.-rdn(ilay-1)*refT
    dtransT = (trans_blw(ilay)*tup(ilay-1)/denom) * &
                (dtransT + (transT*rdn(ilay-1)/denom)*drefT )

    !update derivative of refl. response
    drefT = drefT * (trans_blw(ilay) * tdn(ilay-1) / denom)**2
  enddo

endfunction dtrans_response_up


!---------------------------------------------------------
! derivative of "downgoing" TE or TM reflection response for receivers below source
! with respect to layer ilaym, e.g. dRbTM / depsilon_m
! TE or TM, generic for iso, vti epsh, epsv
! RS 2011
!---------------------------------------------------------
complex(kind=real64) function drefl_response_dn(ilaytop,ilaybot,trans_blw,dz_blw,startidx,endidx,kappa, &
      func_dexp,func_coef) result(drefT)

  implicit none

  !external variables
  integer(kind=int32)   :: ilaytop,ilaybot  !layer indices: start recursion from ilaybot, go to ilaytop
  integer(kind=int32),intent(in)      :: startidx,endidx  !start and end index of transm. coeff. vector
  complex(kind=real64),dimension(startidx:endidx),intent(in) :: trans_blw !transm. coeff. either with respect to src or rec
  real(kind=real64),dimension(startidx:endidx),intent(in)    :: dz_blw    !layer thicknesses
  real(kind=real64),intent(in)    :: kappa    !horiz. wavenumber
  complex(kind=real64),external   :: func_dexp,func_coef  !functions, different for TE TM iso vti epsh epsv

  !internal variables
  integer(kind=int32)   :: ilay    !layer counter
  complex(kind=real64)  :: refT   !reflection response from below
  integer(kind=int32)   :: ilaym1  !layer index for layer just below ilaym
  complex(kind=real64)  :: term1,term2,term3,denom,fact_isotm,fact_dexp


  !initialize
  drefT = 0._real64
  ilaym1 = ilaym + 1

  !term stays zero if layer m is outside range in which this recursion applies
  if ((ilaym.lt.ilaytop) .or. (ilaym.gt.ilaybot)) return

  !derivative of refl. response for layer ilaym
  !- depends on recursive refl. response for next lower layer
  !--> special case: if layer m is the lowest of the stack considered, refl. response for next layer does not exist
  !--> then derivative is zero
  if (ilaym .eq. ilaybot) then
    drefT = 0._real64
    refT = 0._real64
  else
    !reflection response from layers up to the one below ilaym
    refT = refl_response_dn(ilaym1,ilaybot-1,trans_blw,startidx,zeroc)

    !derivative of reflection response:
    !terms containing refTM
    fact_isotm = func_coef(kappa)
    denom = 1.-rup(ilaym)*refT
    term1 = trans_blw(ilaym)**2 * tup(ilaym)**2 * fact_isotm
    term2 = 2. * refT * rdn(ilaym) / denom
    term3 = (refT * tup(ilaym) / denom)**2
    drefT = term1 * (1. - term2 - term3)

    !term from derivative of exponential factor (contains pvert_m) --> contains refTM for layer ilaym --> next iteration
    refT = refl_response_dn(ilaym,ilaym,trans_blw,startidx,refT)
    
    fact_dexp = func_dexp(dz_blw(ilaym),kappa)
    drefT = drefT + refT * fact_dexp
  endif

  !derivative of refl response up to the layer above ilaym
  !special case again since the interface refl./transm. coeff. still depend on ilaym
  if (ilaym.gt.ilaytop) then !nothing more to do if ilaym is top layer of the stack considered
    !FIRST add the term containing deriv. of refl. up to layer m
    denom = 1.-rup(ilaym-1)*refT
    drefT = drefT * (trans_blw(ilaym-1) * tup(ilaym-1) / denom)**2

    !now add special contribution arising from interface refl./transm. coeff. still containing epsilon_m
    fact_isotm = func_coef(kappa)
    term1 = - trans_blw(ilaym-1)**2 * tup(ilaym-1)**2 * fact_isotm
    term2 = 2. * refT * rdn(ilaym-1) / denom
    term3 = (refT * tup(ilaym-1) / denom)**2
    drefT = drefT + term1 * (1. - term2 - term3)
  endif

  !iteration for remaining layers, from 2 layers above ilaym to top of stack considered
  do ilay=ilaym-2,ilaytop,-1
    !next iteration for refl. response --> recursion up to the layer below ilay
    refT = refl_response_dn(ilay+1,ilay+1,trans_blw,startidx,refT)

    !update derivative of refl. response
    denom = 1.-rup(ilay)*refT
    drefT = drefT * (trans_blw(ilay) * tup(ilay) / denom)**2
  enddo

endfunction drefl_response_dn


!---------------------------------------------------------
! derivative of "upgoing" reflection response
! with respect to layer ilaym, e.g. dRbTM / depsilon_m
! TE or TM, generic for iso, vti epsh, epsv
! RS 2011
!---------------------------------------------------------
complex(kind=real64) function drefl_response_up(ilaytop,ilaybot,trans_abv,dz_abv,startidx,endidx,kappa, &
     func_dexp,func_coef) result(drefT)

  implicit none

  !external variables
  integer(kind=int32)   :: ilaytop,ilaybot  !layer indices: start recursion from ilaybot, go to ilaytop
  integer(kind=int32),intent(in)      :: startidx,endidx !start and end index of transm. coeff. vector
  complex(kind=real64),dimension(startidx:endidx),intent(in) :: trans_abv !transm. coeff. either with respect to src or rec
  real(kind=real64),dimension(startidx:endidx),intent(in)    :: dz_abv    !layer thicknesses
  real(kind=real64),intent(in)    :: kappa    !horiz. wavenumber
  complex(kind=real64),external   :: func_dexp,func_coef  !functions, different for TE TM iso vti epsh epsv

  !internal variables
  integer(kind=int32)   :: ilay    !layer counter
  complex(kind=real64)  :: refT   !reflection response from below
  integer(kind=int32)   :: ilaym1  !layer index for layer just above ilaym
  complex(kind=real64)  :: denom   !denominator
  complex(kind=real64)  :: fact_coef,fact_dexp !terms pecific to TM isotropic case
  complex(kind=real64)  :: term1,term2,term3      !parts of derivatives

  !initialize
  drefT = 0._real64
  ilaym1 = ilaym - 1

  !term stays zero if layer m is outside range in which this recursion applies
  if ((ilaym.lt.ilaytop) .or. (ilaym.gt.ilaybot)) return

  !derivative of refl. response for layer ilaym
  !- depends on recursive refl. response for next higher layer
  !--> special case: if layer m is the top of the stack considered, refl. response for layer above does not exist
  !--> then derivative is zero
  if (ilaym .eq. ilaytop) then
    drefT = 0._real64
    refT = 0._real64
  else
    !reflection response from layers down to the one above ilaym
    refT = refl_response_up(ilaytop+1,ilaym1,trans_abv,startidx,zeroc)
    !terms containing refTM
    denom = 1. - rdn(ilaym1) * refT
    fact_coef = func_coef(kappa)
    term1 = (trans_abv(ilaym) * tup(ilaym1))**2 *  fact_coef
    term2 = 2. * refT * rup(ilaym1) / denom
    term3 = (tdn(ilaym1) * refT / denom)**2
    drefT = term1 * (1. - term2 - term3)

    !term from derivative of exponential factor (contains pvert_m) --> contains refTM for layer ilaym --> next iteration
    refT = refl_response_up(ilaym,ilaym,trans_abv,startidx,refT)
    
    fact_dexp = func_dexp(dz_abv(ilaym),kappa)
    drefT = drefT + fact_dexp * refT
  endif

  !derivative of refl response down to the layer below ilaym
  !special case again since the interface refl./transm. coeff. still depend on ilaym
  if (ilaym.lt.ilaybot) then !nothing more to do if ilaym is bottom layer of the stack considered
    !FIRST add the term containing deriv. of refl. up to layer m
    denom = 1. - rdn(ilaym) * refT
    drefT = drefT * (trans_abv(ilaym+1) * tdn(ilaym) / denom)**2

    !now add special contribution arising from interface refl./transm. coeff. still containing epsilon_m
    fact_coef = func_coef(kappa)
    term1 = - (trans_abv(ilaym+1) * tup(ilaym))**2 * fact_coef
    term2 = 2.* refT *rup(ilaym) / denom
    term3 = (tdn(ilaym) * refT / denom)**2
    drefT = drefT + term1 * (1. - term2 - term3)
  endif

  !iteration for remaining layers, from 2 layers below ilaym to bottom of stack considered
  do ilay=ilaym+2,ilaybot
    !next iteration for refl. response --> recursion down to the layer above ilay
    refT = refl_response_up(ilay-1,ilay-1,trans_abv,startidx,refT)

    !update derivative of refl. response
    denom = 1.-rdn(ilay-1)*refT
    drefT = drefT * (trans_abv(ilay) * tdn(ilay-1) / denom)**2
  enddo

endfunction drefl_response_up


!-----------------------------------------
!factors from derivative of exponential terms
!-----------------------------------------
!-----------------------------------------
!TM, isotropic
!TM and TE are equal in isotropic case, but TM uses pvert2 and TE uses pvert1
!for factor 2 in exponent --> divide by 2 when using for transmission response deriv.
!RS 2011
!----------------------------------------------
complex(kind=real64) function fact_exp_tmiso(dz,kappa)
  implicit none
  real(kind=real64)  :: dz  !layer thickness
  real(kind=real64)  :: kappa  !horiz. wavenumber, not needed here but input for consistency with anisotropic functions
  fact_exp_tmiso = j_om_mu * dz / pvert(ilaym)
endfunction fact_exp_tmiso

!----------------------------------------------
!TE, isotropic, also valid for VTI!
!for factor 2 in exponent --> divide by 2 when using for transmission response deriv.
!RS 2011
!----------------------------------------------
complex(kind=real64) function fact_exp_teiso(dz,kappa)
  implicit none
  real(kind=real64)  :: dz  !layer thickness
  real(kind=real64)  :: kappa  !horiz. wavenumber, not needed here but input for consistency with anisotropic functions
  fact_exp_teiso = j_om_mu * dz / pvert(ilaym)
endfunction fact_exp_teiso

!-----------------------------------------
!TM, VTI-anisotropic, for deriv. with respect to epsh
!for factor 2 in exponent --> divide by 2 when using for transmission response deriv.
!RS 2011
!----------------------------------------------
complex(kind=real64) function fact_exp_tmvtih(dz,kappa)
  implicit none
  !external variables
  real(kind=real64)  :: dz     !layer thickness
  real(kind=real64)  :: kappa  !horiz. wavenumber
  !internal variable
  real(kind=real64)  :: phsq   !horiz. slowness
  
  phsq = (kappa / omega)**2
  fact_exp_tmvtih = jomega * dz * (dmu0 - phsq/epsv(ilaym)) / pvert(ilaym)
endfunction fact_exp_tmvtih

!-----------------------------------------
!TM, VTI-anisotropic, for deriv. with respect to epsv
!for factor 2 in exponent --> divide by 2 when using for transmission response deriv.
!RS 2011
!----------------------------------------------
complex(kind=real64) function fact_exp_tmvtiv(dz,kappa)
  implicit none
  !external variables
  real(kind=real64)  :: dz     !layer thickness
  real(kind=real64)  :: kappa  !horiz. wavenumber
  !internal variable
  real(kind=real64)  :: phsq   !horiz. slowness
  
  phsq = (kappa / omega)**2
  fact_exp_tmvtiv = jomega * epsh(ilaym) * phsq * dz / (epsv(ilaym)**2 * pvert(ilaym))
endfunction fact_exp_tmvtiv


!-----------------------------------------
!"coefficient factors"
!-----------------------------------------
!-----------------------------------------
!TM, isotropic
!RS 2011
!----------------------------------------------
complex(kind=real64) function fact_isotm(kappa)
  implicit none
  real(kind=real64),intent(in)  :: kappa  !horiz. wavenumber
  fact_isotm = (2.*pvert(ilaym)**2 - dmu0*epsh(ilaym)) / (4. * pvert(ilaym)**2 * epsh(ilaym))
endfunction fact_isotm

!-----------------------------------------
!TE, isotropic, also valid for VTI
!RS 2011
!----------------------------------------------
complex(kind=real64) function fact_isote(kappa)
  implicit none
  real(kind=real64),intent(in)  :: kappa  !horiz. wavenumber
  fact_isote = dmu0 / (4. * pvert(ilaym)**2)
endfunction fact_isote

!-----------------------------------------
!TM, VTI-anisotropic for epsh
!RS 2011
!----------------------------------------------
complex(kind=real64) function fact_vtitmh(kappa)
  implicit none
  real(kind=real64),intent(in)  :: kappa  !horiz. wavenumber
  fact_vtitmh = 1./(4.*epsh(ilaym))
endfunction fact_vtitmh

!-----------------------------------------
!TM, VTI-anisotropic for epsv
!RS 2011
!----------------------------------------------
complex(kind=real64) function fact_vtitmv(kappa)
  implicit none
  real(kind=real64),intent(in)  :: kappa  !horiz. wavenumber
  real(Kind=real64)    :: phsq   !horiz. slowness
  phsq = (kappa / omega)**2
  fact_vtitmv = - phsq / (4.*epsv(ilaym) * (epsv(ilaym)*dmu0 - phsq))
endfunction fact_vtitmv

