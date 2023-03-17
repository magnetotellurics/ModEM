!****************************************************************
!  1D EM subroutine precomp_intvals_deriv_hed
!    compute integral values for derivatives for all integrals for HED source
!
!  Rita Streich 2010-2011
!****************************************************************
subroutine precomp_intvals_deriv_hed(refl_var,sz,zr,ilay, &
    funcA0TE,funcA0TM,funcA1TE,funcA1TM,funcDz1TM,funcD0TE,funcD0TM,funcD1TE,funcD1TM,funcAz1TE,funcDz1TMfwd, &
    funcA0TMv,funcA1TMv,funcDz1TMv,funcD0TMv,funcD1TMv)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and receiver depth
  integer(kind=int32),intent(in) :: ilay     !layer index, compute derivatives with respect to this layer
  !derivatives of integrals with respect to epsh or isotropic
  complex(kind=real64),external  :: funcA0TE,funcA0TM,funcA1TE,funcA1TM,funcDz1TM,funcD0TE,funcD0TM,funcD1TE,funcD1TM,funcAz1TE
  complex(kind=real64),external,optional   :: funcDz1TMfwd  !function only needed for derivatives in receiver layer
  !derivatives of TM integrals with respect to epsv in VTI case
  complex(kind=real64),external,optional  :: funcA0TMv,funcA1TMv,funcDz1TMv,funcD0TMv,funcD1TMv

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter
  integer(kind=int32),dimension(nrelmax)  :: iint  !indicates locations for integration results in output structure

  !ilaym is defined in refl_mod, visible nearly everywhere
  ilaym = ilay

  !receivers above or below source, not exactly at the same depth
  !TEST: the "first" of each related set of integrals is always well-behaved, 
  !  so the others with higher powers of kappa should be ok too???
  !--> would not need adaptive integration then any more??? - but result is not reliable
  if (sz .ne. zr) then

    !related integrals: IA1TEderiv, IA0TEderiv, IAz1TEderiv, same integrand just different kappa factors and Bessel function orders
    !basis is IA1TEderiv, contains no kappa, Bessel order 1
    !IA0TEderiv contains kappa^1, Bessel order 0
    !IAz1TEderiv contains kappa^2, Bessel order 1
    ibesord(1:3:2) = 1
    ibesord(2) = 0
    ijrel = 1
    ijrel(1,3) = 2 !factor kappa^2 for third integral
    iint(1) = 3    !integral IabvA1TEderiv
    iint(2) = 1    !integral IabvA0TEderiv
    iint(3) = 10   !integral IabvAz1TEderiv (for Hz)

    call precomp_intval_fht_rel(refl_var,3,ibesord,ijrel,funcA1TE,iint)

    !related integrals IA1TMderiv and IA0TMderiv
    !Bessel order as before (using first 2 elements only)
    !exponents for first 2 related transforms also as before
    iint(1) = 4    !integral IabvA1TMderiv
    iint(2) = 2    !integral IabvA0TMderiv
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcA1TM,iint)

    !related integrals ID1TEderiv and ID0TEderiv
    iint(1) = 8    !integral IabvD1TEderiv
    iint(2) = 6    !integral IabvD0TEderiv
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcD1TE,iint)

    !related integrals D_TMderiv
    iint(1) = 9    !integral iabvD1TMderiv
    iint(2) = 7    !integral iabvD0TMderiv
    iint(3) = 5    !integral iabvDz1TMderiv (for Ez)
    call precomp_intval_fht_rel(refl_var,3,ibesord,ijrel,funcD1TM,iint)
 
    !special in receiver layer: there is an eps in front of the integral, so we also need 
    ! derivative of term in front of integral times integral without derivative
    if (ilaym .eq. ilayrec) then
      ibesord = 1
      call precomp_intval_fht(refl_var,ibesord,funcDz1TMfwd,11)
    endif

    !VTI-anisotropy: derivatives of TM integrals for epsv
    if (aniso.ne.iso) then
      ibesord(2) = 0
      !related integrals IA1TMderivv and IA0TMderivv
      iint(1) = 12    !integral Iabv/blwA1TMderivv
      iint(2) = 13    !integral Iabv/blwA0TMderivv
      call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcA1TMv,iint)
 
      !related integrals D_TMderiv
      iint(1) = 14    !integral iabv/blwD1TMderivv
      iint(2) = 15    !integral iabv/blwD0TMderivv
      iint(3) = 16    !integral iabv/blwDz1TMderivv (for Ez)
      call precomp_intval_fht_rel(refl_var,3,ibesord,ijrel,funcD1TMv,iint)
    endif

  !receiver exactly at source depth
  else

    write(*,'(a)') 'WARNING: source depth = receiver depth, entering adaptive integration, this can be SLOW!'

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,funcA0TM,2,sz,zr)
    besorder = 1._real64
    call precomp_intval_adaptive(refl_var,besorder,funcDz1TM,5,sz,zr)

    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,funcD0TE,6,sz,zr)
    call precomp_intval_adaptive(refl_var,besorder,funcD0TM,7,sz,zr)

    if (ilaym .eq. ilayrec) then
      besorder = 1._real64
      call precomp_intval_adaptive(refl_var,besorder,funcDz1TMfwd,11,sz,zr)
    endif

    !anisotropic, badly behaved
    if (aniso.ne.iso) then
      besorder = 0._real64
      call precomp_intval_adaptive(refl_var,besorder,funcA0TMv,13,sz,zr)
      besorder = 1._real64
      call precomp_intval_adaptive(refl_var,besorder,funcDz1TMv,16,sz,zr)
      besorder = 0._real64
      call precomp_intval_adaptive(refl_var,besorder,funcD0TMv,15,sz,zr)
    endif

    !do spline interpolation here and not inside precomp_intval_adaptive
    ! so that there is no need to go back and forth between radii and their logarithms

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,2))
    call spline(refl_var%radlog,refl_var%intvalim(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,2))
    call spline(refl_var%radlog,refl_var%intvalre(:,5),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,5))
    call spline(refl_var%radlog,refl_var%intvalim(:,5),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,5))
    call spline(refl_var%radlog,refl_var%intvalre(:,6),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,6))
    call spline(refl_var%radlog,refl_var%intvalim(:,6),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,6))
    call spline(refl_var%radlog,refl_var%intvalre(:,7),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,7))
    call spline(refl_var%radlog,refl_var%intvalim(:,7),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,7))

    if (ilaym .eq. ilayrec) then
      call spline(refl_var%radlog,refl_var%intvalre(:,11),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,11))
      call spline(refl_var%radlog,refl_var%intvalim(:,11),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,11))
    endif
    
    !anisotropic, badly behaved
    if (aniso.ne.iso) then
      call spline(refl_var%radlog,refl_var%intvalre(:,13),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,13))
      call spline(refl_var%radlog,refl_var%intvalim(:,13),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,13))
      call spline(refl_var%radlog,refl_var%intvalre(:,16),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,16))
      call spline(refl_var%radlog,refl_var%intvalim(:,16),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,16))
      call spline(refl_var%radlog,refl_var%intvalre(:,15),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,15))
      call spline(refl_var%radlog,refl_var%intvalim(:,15),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,15))
    endif

    !compute well-behaved integrals by fast Hankel transform
    !related integrals: IA1TEderiv, IA0TEderiv, IAz1TEderiv
    ibesord(1:3:2) = 1
    ibesord(2) = 0
    ijrel = 1
    ijrel(1,3) = 2 !factor kappa^2 for third integral
    iint(1) = 3    !integral IabvA1TEderiv
    iint(2) = 1    !integral IabvA0TEderiv
    iint(3) = 10   !integral IabvAz1TEderiv (for Hz)
    call precomp_intval_fht_rel(refl_var,3,ibesord,ijrel,funcA1TE,iint)

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,funcA1TM,4)

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,funcD1TE,8)
    call precomp_intval_fht(refl_var,ibesord,funcD1TM,9)
    
    if (aniso.ne.iso) then
      ibesord = 1
      call precomp_intval_fht(refl_var,ibesord,funcA1TMv,12)
      call precomp_intval_fht(refl_var,ibesord,funcD1TMv,14)
    endif

  endif

endsubroutine precomp_intvals_deriv_hed


!****************************************************************
!  1D EM subroutine precomp_intvals_deriv_hed_Exy
!    compute integral values for derivatives for all integrals for HED source, Ex/Ey only
!
!  Rita Streich 2010-2011
!****************************************************************
subroutine precomp_intvals_deriv_hed_Exy(refl_var,sz,zr,ilay, &
    funcA0TE,funcA0TM,funcA1TE,funcA1TM, &
    funcA0TMv,funcA1TMv)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and receiver depth
  integer(kind=int32),intent(in) :: ilay     !layer index, compute derivatives with respect to this layer
  !derivatives of integrals with respect to epsh or isotropic
  complex(kind=real64),external  :: funcA0TE,funcA0TM,funcA1TE,funcA1TM
  !derivatives of TM integrals with respect to epsv in VTI case
  complex(kind=real64),external,optional  :: funcA0TMv,funcA1TMv

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter
  integer(kind=int32),dimension(nrelmax)  :: iint  !indicates locations for integration results in output structure

  !ilaym is defined in refl_mod, visible nearly everywhere
  ilaym = ilay

  !receivers above or below source, not exactly at the same depth
  !TEST: the "first" of each related set of integrals is always well-behaved, 
  !  so the others with higher powers of kappa should be ok too???
  !--> would not need adaptive integration then any more??? - but result is not reliable
  if (sz .ne. zr) then

    !related integrals: IA1TEderiv, IA0TEderiv, IAz1TEderiv, same integrand just different kappa factors and Bessel function orders
    !basis is IA1TEderiv, contains no kappa, Bessel order 1
    !IA0TEderiv contains kappa^1, Bessel order 0
    !IAz1TEderiv contains kappa^2, Bessel order 1
    ibesord = 1
    ibesord(2) = 0
    ijrel = 1
    iint(1) = 3    !integral IabvA1TEderiv
    iint(2) = 1    !integral IabvA0TEderiv

    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcA1TE,iint)

    !related integrals IA1TMderiv and IA0TMderiv
    !Bessel order as before (using first 2 elements only)
    !exponents for first 2 related transforms also as before
    iint(1) = 4    !integral IabvA1TMderiv
    iint(2) = 2    !integral IabvA0TMderiv
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcA1TM,iint)

    !VTI-anisotropy: derivatives of TM integrals for epsv
    if (aniso.ne.iso) then
      ibesord(2) = 0
      !related integrals IA1TMderivv and IA0TMderivv
      iint(1) = 12    !integral Iabv/blwA1TMderivv
      iint(2) = 13    !integral Iabv/blwA0TMderivv
      call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcA1TMv,iint)
 
    endif

  !receiver exactly at source depth
  else

    if (refl_var%infolevel .gt. output_final) &
      write(*,'(a)') 'WARNING: source depth = receiver depth, entering adaptive integration, this can be SLOW!'

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,funcA0TM,2,sz,zr)

    !anisotropic, badly behaved
    if (aniso.ne.iso) then
      besorder = 0._real64
      call precomp_intval_adaptive(refl_var,besorder,funcA0TMv,13,sz,zr)
    endif

    !do spline interpolation here and not inside precomp_intval_adaptive
    ! so that there is no need to go back and forth between radii and their logarithms

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,2))
    call spline(refl_var%radlog,refl_var%intvalim(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,2))

    !anisotropic, badly behaved
    if (aniso.ne.iso) then
      call spline(refl_var%radlog,refl_var%intvalre(:,13),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,13))
      call spline(refl_var%radlog,refl_var%intvalim(:,13),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,13))
    endif

    !compute well-behaved integrals by fast Hankel transform
    !related integrals: IA1TEderiv, IA0TEderiv, IAz1TEderiv
    ibesord = 1
    ibesord(2) = 0
    ijrel = 1
    iint(1) = 3    !integral IabvA1TEderiv
    iint(2) = 1    !integral IabvA0TEderiv
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcA1TE,iint)

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,funcA1TM,4)

    if (aniso.ne.iso) then
      ibesord = 1
      call precomp_intval_fht(refl_var,ibesord,funcA1TMv,12)
    endif

  endif

endsubroutine precomp_intvals_deriv_hed_Exy


!****************************************************************
!  1D EM subroutine precomp_intvals_deriv_hed_Ez
!    compute integral values for derivatives for all integrals for HED source, Ez only
!
!  Rita Streich 2010-2011
!****************************************************************
subroutine precomp_intvals_deriv_hed_Ez(refl_var,sz,zr,ilay, &
    funcDz1TM,funcDz1TMfwd, &
    funcDz1TMv)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and receiver depth
  integer(kind=int32),intent(in) :: ilay     !layer index, compute derivatives with respect to this layer
  !derivatives of integrals with respect to epsh or isotropic
  complex(kind=real64),external  :: funcDz1TM
  complex(kind=real64),external,optional   :: funcDz1TMfwd  !function only needed for derivatives in receiver layer
  !derivatives of TM integrals with respect to epsv in VTI case
  complex(kind=real64),external,optional  :: funcDz1TMv

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter

  !ilaym is defined in refl_mod, visible nearly everywhere
  ilaym = ilay

  !receivers above or below source, not exactly at the same depth
  !TEST: the "first" of each related set of integrals is always well-behaved, 
  !  so the others with higher powers of kappa should be ok too???
  !--> would not need adaptive integration then any more??? - but result is not reliable
  if (sz .ne. zr) then

    !related integrals: IA1TEderiv, IA0TEderiv, IAz1TEderiv, same integrand just different kappa factors and Bessel function orders
    !basis is IA1TEderiv, contains no kappa, Bessel order 1
    !IA0TEderiv contains kappa^1, Bessel order 0
    !IAz1TEderiv contains kappa^2, Bessel order 1
    ibesord = 1

    !related integrals D_TMderiv
    call precomp_intval_fht(refl_var,ibesord,funcDz1TM,5)

    !special in receiver layer: there is an eps in front of the integral, so we also need 
    ! derivative of term in front of integral times integral without derivative
    if (ilaym .eq. ilayrec) then
      call precomp_intval_fht(refl_var,ibesord,funcDz1TMfwd,11)
    endif

    !VTI-anisotropy: derivatives of TM integrals for epsv
    if (aniso.ne.iso) then
      call precomp_intval_fht(refl_var,ibesord,funcDz1TMv,16)
    endif

  !receiver exactly at source depth
  else

    write(*,'(a)') 'WARNING: source depth = receiver depth, entering adaptive integration, this can be SLOW!'

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    besorder = 1._real64

    !compute badly behaved integrals by adaptive integration
    call precomp_intval_adaptive(refl_var,besorder,funcDz1TM,5,sz,zr)

    if (ilaym .eq. ilayrec) then
      call precomp_intval_adaptive(refl_var,besorder,funcDz1TMfwd,11,sz,zr)
    endif

    !anisotropic, badly behaved
    if (aniso.ne.iso) then
      call precomp_intval_adaptive(refl_var,besorder,funcDz1TMv,16,sz,zr)
    endif

    !do spline interpolation here and not inside precomp_intval_adaptive
    ! so that there is no need to go back and forth between radii and their logarithms

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,5),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,5))
    call spline(refl_var%radlog,refl_var%intvalim(:,5),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,5))

    if (ilaym .eq. ilayrec) then
      call spline(refl_var%radlog,refl_var%intvalre(:,11),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,11))
      call spline(refl_var%radlog,refl_var%intvalim(:,11),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,11))
    endif
    
    !anisotropic, badly behaved
    if (aniso.ne.iso) then
      call spline(refl_var%radlog,refl_var%intvalre(:,16),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,16))
      call spline(refl_var%radlog,refl_var%intvalim(:,16),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,16))
    endif

  endif

endsubroutine precomp_intvals_deriv_hed_Ez


!****************************************************************
!  1D EM subroutine precomp_intvals_deriv_hed_Hxy
!    compute integral values for derivatives for all integrals for HED source, Hx/Hy only
!
!  Rita Streich 2010-2011
!****************************************************************
subroutine precomp_intvals_deriv_hed_Hxy(refl_var,sz,zr,ilay, &
    funcD0TE,funcD0TM,funcD1TE,funcD1TM, funcD0TMv,funcD1TMv)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and receiver depth
  integer(kind=int32),intent(in) :: ilay     !layer index, compute derivatives with respect to this layer
  !derivatives of integrals with respect to epsh or isotropic
  complex(kind=real64),external  :: funcD0TE,funcD0TM,funcD1TE,funcD1TM
  !derivatives of TM integrals with respect to epsv in VTI case
  complex(kind=real64),external,optional  :: funcD0TMv,funcD1TMv

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter
  integer(kind=int32),dimension(nrelmax)  :: iint  !indicates locations for integration results in output structure

  !ilaym is defined in refl_mod, visible nearly everywhere
  ilaym = ilay

  !receivers above or below source, not exactly at the same depth
  !TEST: the "first" of each related set of integrals is always well-behaved, 
  !  so the others with higher powers of kappa should be ok too???
  !--> would not need adaptive integration then any more??? - but result is not reliable
  if (sz .ne. zr) then

    !related integrals: IA1TEderiv, IA0TEderiv, IAz1TEderiv, same integrand just different kappa factors and Bessel function orders
    !basis is IA1TEderiv, contains no kappa, Bessel order 1
    !IA0TEderiv contains kappa^1, Bessel order 0
    !IAz1TEderiv contains kappa^2, Bessel order 1
    ibesord = 1
    ibesord(2) = 0
    ijrel = 1

    !related integrals ID1TEderiv and ID0TEderiv
    iint(1) = 8    !integral IabvD1TEderiv
    iint(2) = 6    !integral IabvD0TEderiv
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcD1TE,iint)

    !related integrals D_TMderiv
    iint(1) = 9    !integral iabvD1TMderiv
    iint(2) = 7    !integral iabvD0TMderiv
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcD1TM,iint)
 
    !VTI-anisotropy: derivatives of TM integrals for epsv
    if (aniso.ne.iso) then
      ibesord(2) = 0
 
      !related integrals D_TMderiv
      iint(1) = 14    !integral iabv/blwD1TMderivv
      iint(2) = 15    !integral iabv/blwD0TMderivv
      call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcD1TMv,iint)
    endif

  !receiver exactly at source depth
  else

    write(*,'(a)') 'WARNING: source depth = receiver depth, entering adaptive integration, this can be SLOW!'

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,funcD0TE,6,sz,zr)
    call precomp_intval_adaptive(refl_var,besorder,funcD0TM,7,sz,zr)

    !anisotropic, badly behaved
    if (aniso.ne.iso) then
      besorder = 0._real64
      call precomp_intval_adaptive(refl_var,besorder,funcD0TMv,15,sz,zr)
    endif

    !do spline interpolation here and not inside precomp_intval_adaptive
    ! so that there is no need to go back and forth between radii and their logarithms

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,6),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,6))
    call spline(refl_var%radlog,refl_var%intvalim(:,6),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,6))
    call spline(refl_var%radlog,refl_var%intvalre(:,7),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,7))
    call spline(refl_var%radlog,refl_var%intvalim(:,7),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,7))

    !anisotropic, badly behaved
    if (aniso.ne.iso) then
      call spline(refl_var%radlog,refl_var%intvalre(:,15),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,15))
      call spline(refl_var%radlog,refl_var%intvalim(:,15),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,15))
    endif

    !compute well-behaved integrals by fast Hankel transform
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,funcD1TE,8)
    call precomp_intval_fht(refl_var,ibesord,funcD1TM,9)
    
    if (aniso.ne.iso) then
      ibesord = 1
      call precomp_intval_fht(refl_var,ibesord,funcD1TMv,14)
    endif

  endif

endsubroutine precomp_intvals_deriv_hed_Hxy


!****************************************************************
!  1D EM subroutine precomp_intvals_deriv_hed_Hz
!    compute integral values for derivatives for all integrals for HED source, Hz only
!
!  Rita Streich 2010-2011
!****************************************************************
subroutine precomp_intvals_deriv_hed_Hz(refl_var,sz,zr,ilay, funcAz1TE)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and receiver depth
  integer(kind=int32),intent(in) :: ilay     !layer index, compute derivatives with respect to this layer
  !derivatives of integrals with respect to epsh or isotropic
  complex(kind=real64),external  :: funcAz1TE

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")

  !ilaym is defined in refl_mod, visible nearly everywhere
  ilaym = ilay

  !the integral _should be_ well-behaved, no need to do adaptive integration
  ibesord = 1
  call precomp_intval_fht(refl_var,ibesord,funcAz1TE,10)

endsubroutine precomp_intvals_deriv_hed_Hz


!****************************************************************
!  1D EM subroutine precomp_intvals_deriv_ved
!    compute integral values for derivatives for all integrals for VED source
!
!  Rita Streich 2010-2011
!****************************************************************
subroutine precomp_intvals_deriv_ved(refl_var,sz,zr,ilay, &
    funcB1TMved,funcC0TMved,funcC1TMved, funcB1TMfwd,funcC0TMfwd,funcC1TMfwd, funcB1TMvedv,funcC0TMvedv,funcC1TMvedv)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and receiver depth
  integer(kind=int32),intent(in) :: ilay     !layer index, compute derivatives with respect to this layer
  complex(kind=real64),external   :: funcB1TMved,funcC0TMved,funcC1TMved
  complex(kind=real64),external   :: funcB1TMfwd,funcC0TMfwd,funcC1TMfwd
  complex(kind=real64),external,optional   :: funcB1TMvedv,funcC0TMvedv,funcC1TMvedv

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter
  integer(kind=int32),dimension(nrelmax)  :: iint  !indicates locations for integration results in output structure


  !ilaym is defined in medium_mod, visible nearly everywhere
  ilaym = ilay

  !receivers not exactly at source depth
  if (sz .ne. zr) then

    !compute well-behaved integrals by fast Hankel transform
    !related integrals: iabvC1TMvedderiv and iabvC0TMvedderiv
    ibesord(1) = 1
    ibesord(2) = 0
    ijrel = 1
    iint(1) = 3    !integral iabvCz1TMderiv (for Hx and Hy)
    iint(2) = 2    !integral iabvC0TMvedderiv (for Ez)
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcC1TMved,iint)

    !integrals for Ex and Ey
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,funcB1TMved,1)

    !need special terms for ALL components in source layer
    if (ilaym .eq. ilaysrc) then

      ibesord = 1
      call precomp_intval_fht(refl_var,ibesord,funcB1TMfwd,4) !Ex and Ey

        !related integrals: iabvC1TMved and iabvC0TMved
        ibesord(2) = 0
        iint(1) = 6    !integral iabvC1TMved
        iint(2) = 5    !integral iabvC0TMved
        call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcC1TMfwd,iint)

    !special term for Ez needed in source AND receiver layer
    elseif (ilaym .eq. ilayrec) then
      ibesord = 0
      call precomp_intval_fht(refl_var,ibesord,funcC0TMfwd,5) !Ez
    endif
    
    !VTI: derivatives of integrals for epsv
    if (aniso.ne.iso) then
      !related integrals: iabvC1TMvedderiv and iabvC0TMvedderiv
      ibesord(1) = 1
      ibesord(2) = 0
      ijrel = 1
      iint(1) = 9    !integral iabvCz1TMderivv (for Hx and Hy)
      iint(2) = 8    !integral iabvC0TMvedderivv (for Ez)
      call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcC1TMvedv,iint)

      !integral for Ex and Ey
      ibesord = 1
      call precomp_intval_fht(refl_var,ibesord,funcB1TMvedv,7)
    endif

  !receiver exactly at source depth
  else

    !precompute radii at logarithmic spacing, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 1._real64
    call precomp_intval_adaptive(refl_var,besorder,funcB1TMved,1,sz,zr)
    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,funcC0TMved,2,sz,zr)

    !special term for Ez needed in source AND receiver layer
    if ((ilaym .eq. ilayrec) .or. (ilaym .eq. ilaysrc)) then
      besorder = 0._real64
      call precomp_intval_adaptive(refl_var,besorder,funcC0TMfwd,5,sz,zr)
    endif
    !need special terms for all components in source layer
    if (ilaym .eq. ilaysrc) then
      besorder = 1._real64
      call precomp_intval_adaptive(refl_var,besorder,funcB1TMfwd,4,sz,zr)
    endif

    !VTI: derivatives of integrals for epsv
    if (aniso.ne.iso) then
      besorder = 1._real64
      call precomp_intval_adaptive(refl_var,besorder,funcB1TMvedv,7,sz,zr)
      besorder = 0._real64
      call precomp_intval_adaptive(refl_var,besorder,funcC0TMvedv,8,sz,zr)
    endif

    !do spline interpolation here and not inside precomp_intval_adaptive
    ! so that there is no need to go back and forth between radii and their logarithms

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,1),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,1))
    call spline(refl_var%radlog,refl_var%intvalim(:,1),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,1))
    call spline(refl_var%radlog,refl_var%intvalre(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,2))
    call spline(refl_var%radlog,refl_var%intvalim(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,2))

    if (aniso.ne.iso) then
      call spline(refl_var%radlog,refl_var%intvalre(:,7),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,7))
      call spline(refl_var%radlog,refl_var%intvalim(:,7),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,7))
      call spline(refl_var%radlog,refl_var%intvalre(:,8),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,8))
      call spline(refl_var%radlog,refl_var%intvalim(:,8),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,8))
    endif

    !special term for Ez needed in source AND receiver layer
    if ((ilaym .eq. ilayrec) .or. (ilaym .eq. ilaysrc)) then
      call spline(refl_var%radlog,refl_var%intvalre(:,5),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,5))
      call spline(refl_var%radlog,refl_var%intvalim(:,5),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,5))
    endif
    !need special terms for all components in source layer
    if (ilaym .eq. ilaysrc) then
      call spline(refl_var%radlog,refl_var%intvalre(:,4),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,4))
      call spline(refl_var%radlog,refl_var%intvalim(:,4),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,4))

      ibesord = 1
      call precomp_intval_fht(refl_var,ibesord,funcC1TMfwd,6)
    endif

    !compute well-behaved integrals by fast Hankel transform
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,funcC1TMved,3)
    if (aniso.ne.iso) then
      call precomp_intval_fht(refl_var,ibesord,funcC1TMvedv,9)
    endif
    
  endif !receiver location above or below source

endsubroutine precomp_intvals_deriv_ved


!****************************************************************
!  1D EM subroutine precomp_intvals_deriv_ved_Exy
!    compute integral values for derivatives for all integrals for VED source, Ex/Ey only
!
!  Rita Streich 2010-2011
!****************************************************************
subroutine precomp_intvals_deriv_ved_Exy(refl_var,sz,zr,ilay, funcB1TMved, funcB1TMfwd, funcB1TMvedv)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and receiver depth
  integer(kind=int32),intent(in) :: ilay     !layer index, compute derivatives with respect to this layer
  complex(kind=real64),external   :: funcB1TMved
  complex(kind=real64),external   :: funcB1TMfwd
  complex(kind=real64),external,optional   :: funcB1TMvedv

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter


  !ilaym is defined in medium_mod, visible nearly everywhere
  ilaym = ilay

  !receivers not exactly at source depth
  if (sz .ne. zr) then

    !integrals for Ex and Ey
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,funcB1TMved,1)

    !need special terms for ALL components in source layer
    if (ilaym .eq. ilaysrc) then
      ibesord = 1
      call precomp_intval_fht(refl_var,ibesord,funcB1TMfwd,4) !Ex and Ey
    endif
    
    !VTI: derivatives of integrals for epsv
    if (aniso.ne.iso) then
      !integral for Ex and Ey
      ibesord = 1
      call precomp_intval_fht(refl_var,ibesord,funcB1TMvedv,7)
    endif

  !receiver exactly at source depth
  else

    !precompute radii at logarithmic spacing, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 1._real64
    call precomp_intval_adaptive(refl_var,besorder,funcB1TMved,1,sz,zr)

    !need special terms for all components in source layer
    if (ilaym .eq. ilaysrc) then
      besorder = 1._real64
      call precomp_intval_adaptive(refl_var,besorder,funcB1TMfwd,4,sz,zr)
    endif

    !VTI: derivatives of integrals for epsv
    if (aniso.ne.iso) then
      besorder = 1._real64
      call precomp_intval_adaptive(refl_var,besorder,funcB1TMvedv,7,sz,zr)
    endif

    !do spline interpolation here and not inside precomp_intval_adaptive
    ! so that there is no need to go back and forth between radii and their logarithms

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,1),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,1))
    call spline(refl_var%radlog,refl_var%intvalim(:,1),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,1))

    if (aniso.ne.iso) then
      call spline(refl_var%radlog,refl_var%intvalre(:,7),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,7))
      call spline(refl_var%radlog,refl_var%intvalim(:,7),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,7))
    endif

    !need special terms for all components in source layer
    if (ilaym .eq. ilaysrc) then
      call spline(refl_var%radlog,refl_var%intvalre(:,4),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,4))
      call spline(refl_var%radlog,refl_var%intvalim(:,4),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,4))
    endif

  endif !receiver location above or below source

endsubroutine precomp_intvals_deriv_ved_Exy


!****************************************************************
!  1D EM subroutine precomp_intvals_deriv_ved_Ez
!    compute integral values for derivatives for all integrals for VED source, Ez only
!
!  Rita Streich 2010-2011
!****************************************************************
subroutine precomp_intvals_deriv_ved_Ez(refl_var,sz,zr,ilay, funcC0TMved, funcC0TMfwd,funcC0TMvedv)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and receiver depth
  integer(kind=int32),intent(in) :: ilay     !layer index, compute derivatives with respect to this layer
  complex(kind=real64),external   :: funcC0TMved
  complex(kind=real64),external   :: funcC0TMfwd
  complex(kind=real64),external,optional   :: funcC0TMvedv

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter


  !ilaym is defined in medium_mod, visible nearly everywhere
  ilaym = ilay

  ibesord = 0

  !receivers not exactly at source depth
  if (sz .ne. zr) then
    call precomp_intval_fht(refl_var,ibesord,funcC0TMved,2) !integral iabvC0TMvedderiv (for Ez)

    if ((ilaym .eq. ilaysrc) .or. (ilaym .eq. ilayrec)) then
      !special term for Ez needed in source AND receiver layer
      call precomp_intval_fht(refl_var,ibesord,funcC0TMfwd,5) !integral iabvC0TMved (for Ez)
    endif
    
    !VTI: derivatives of integrals for epsv
    if (aniso.ne.iso) then
      call precomp_intval_fht(refl_var,ibesord,funcC0TMvedv,8)  !integral iabvC0TMvedderivv (for Ez)
    endif

  !receiver exactly at source depth
  else

    !precompute radii at logarithmic spacing, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    besorder = 0._real64

    !compute badly behaved integrals by adaptive integration
    call precomp_intval_adaptive(refl_var,besorder,funcC0TMved,2,sz,zr)

    !special term for Ez needed in source AND receiver layer
    if ((ilaym .eq. ilayrec) .or. (ilaym .eq. ilaysrc)) then
      call precomp_intval_adaptive(refl_var,besorder,funcC0TMfwd,5,sz,zr)
    endif

    !VTI: derivatives of integrals for epsv
    if (aniso.ne.iso) then
      call precomp_intval_adaptive(refl_var,besorder,funcC0TMvedv,8,sz,zr)
    endif

    !do spline interpolation here and not inside precomp_intval_adaptive
    ! so that there is no need to go back and forth between radii and their logarithms

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,2))
    call spline(refl_var%radlog,refl_var%intvalim(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,2))

    if (aniso.ne.iso) then
      call spline(refl_var%radlog,refl_var%intvalre(:,8),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,8))
      call spline(refl_var%radlog,refl_var%intvalim(:,8),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,8))
    endif

    !special term for Ez needed in source AND receiver layer
    if ((ilaym .eq. ilayrec) .or. (ilaym .eq. ilaysrc)) then
      call spline(refl_var%radlog,refl_var%intvalre(:,5),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,5))
      call spline(refl_var%radlog,refl_var%intvalim(:,5),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,5))
    endif
    
  endif !receiver location above or below source

endsubroutine precomp_intvals_deriv_ved_Ez


!****************************************************************
!  1D EM subroutine precomp_intvals_deriv_ved_Hxy
!    compute integral values for derivatives for all integrals for VED source
!
!  Rita Streich 2010-2011
!****************************************************************
subroutine precomp_intvals_deriv_ved_Hxy(refl_var,sz,zr,ilay, funcC1TMved, funcC1TMfwd, funcC1TMvedv)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and receiver depth
  integer(kind=int32),intent(in) :: ilay     !layer index, compute derivatives with respect to this layer
  complex(kind=real64),external   :: funcC1TMved
  complex(kind=real64),external   :: funcC1TMfwd
  complex(kind=real64),external,optional   :: funcC1TMvedv

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")


  !ilaym is defined in medium_mod, visible nearly everywhere
  ilaym = ilay

  ibesord = 1

  !all integrals well-behaved, no need to distinguish special case for source depth = receiver depth

    call precomp_intval_fht(refl_var,ibesord,funcC1TMved,3) !integral iabvCz1TMderiv (for Hx and Hy)

    !need special terms for ALL components in source layer
    if (ilaym .eq. ilaysrc) then
      call precomp_intval_fht(refl_var,ibesord,funcC1TMfwd,6) !integral iabvC1TMved (for Hx and Hy)
    endif
    
    !VTI: derivatives of integrals for epsv
    if (aniso.ne.iso) then
      call precomp_intval_fht(refl_var,ibesord,funcC1TMvedv,9)  !integral iabvCz1TMderivv (for Hx and Hy)
    endif

endsubroutine precomp_intvals_deriv_ved_Hxy


!****************************************************************
!  1D EM subroutine precomp_intvals_deriv_hmd
!    compute integral values for derivatives for all integrals for HMD source
!
!  Rita Streich 2010
!****************************************************************
subroutine precomp_intvals_deriv_hmd(refl_var,sz,zr,ilay, &
    funcB0TE,funcB0TM,funcB1TE,funcB1TM,funcCz1TM,funcC0TE,funcC0TM,funcC1TE,funcC1TM,funcBz1TE,funcCz1TMfwd, &
    funcB0TMv,funcB1TMv,funcCz1TMv,funcC0TMv,funcC1TMv)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and receiver depth
  integer(kind=int32),intent(in) :: ilay     !layer index, compute derivatives with respect to this layer
  complex(kind=real64),external   :: funcB0TE,funcB0TM,funcB1TE,funcB1TM,funcCz1TM,funcC0TE,funcC0TM,funcC1TE,funcC1TM,funcBz1TE
  complex(kind=real64),external   :: funcCz1TMfwd  !function only needed for derivatives in receiver layer
  complex(kind=real64),external,optional   :: funcB0TMv,funcB1TMv,funcCz1TMv,funcC0TMv,funcC1TMv !functions for epsv derivatives

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter
  integer(kind=int32),dimension(nrelmax)  :: iint  !indicates locations for integration results in output structure


  !ilaym is defined in medium_mod, visible nearly everywhere
  ilaym = ilay

  !receivers above or below source, not exactly at the same depth
  if (sz .ne. zr) then

    !related integrals: IB1TEderiv, IB0TEderiv, IBz1TEderiv, same integrand just different kappa factors and Bessel function orders
    !basis is IB1TEderiv, contains no kappa, Bessel order 1
    !IB0TEderiv contains kappa^1, Bessel order 0
    !IBz1TEderiv contains kappa^2, Bessel order 1
    ibesord(1:3:2) = 1
    ibesord(2) = 0
    ijrel = 1
    ijrel(1,3) = 2 !factor kappa^2 for third integral
    iint(1) = 3    !integral IabvB1TEderiv
    iint(2) = 1    !integral IabvB0TEderiv
    iint(3) = 10   !integral IabvBz1TEderiv (for Hz)

    call precomp_intval_fht_rel(refl_var,3,ibesord,ijrel,funcB1TE,iint)

    !related integrals IB1TMderiv and IB0TMderiv
    !Bessel order as before (using first 2 elements only)
    !exponents for first 2 related transforms also as before
    iint(1) = 4    !integral IabvB1TMderiv
    iint(2) = 2    !integral IabvB0TMderiv
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcB1TM,iint)

    !related integrals IC1TEderiv and IC0TEderiv
    iint(1) = 8    !integral IabvC1TEderiv
    iint(2) = 6    !integral IabvC0TEderiv
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcC1TE,iint)

    !related integrals C_TM
    iint(1) = 9    !integral iabvC1TMderiv
    iint(2) = 7    !integral iabvC0TMderiv
    iint(3) = 5    !integral iabvCz1TMderiv (for Ez)
    call precomp_intval_fht_rel(refl_var,3,ibesord,ijrel,funcC1TM,iint)

    !special in receiver layer: there is an eps in front of the integral, so we also need 
    ! derivative of term in front of integral times integral without derivative
    if (ilaym .eq. ilayrec) then
      ibesord = 1
      call precomp_intval_fht(refl_var,ibesord,funcCz1TMfwd,11)
    endif
    
    !VTI: derivatives for epsv
    if (aniso.ne.iso) then
      ibesord(2) = 0
      !related integrals IB1TMderivv and IB0TMderivv
      iint(1) = 12    !integral Iabv/blwB1TMderivv
      iint(2) = 13    !integral Iabv/blwB0TMderivv
      call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcB1TMv,iint)

      !related integrals C_TMv
      iint(1) = 14    !integral iabv/blwC1TMderivv
      iint(2) = 15    !integral iabv/blwC0TMderivv
      iint(3) = 16    !integral iabv/blwCz1TMderivv (for Ez)
      call precomp_intval_fht_rel(refl_var,3,ibesord,ijrel,funcC1TMv,iint)
    endif

  !receiver exactly at source depth
  else

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,funcB0TE,1,sz,zr)
    call precomp_intval_adaptive(refl_var,besorder,funcB0TM,2,sz,zr)

    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,funcC0TE,6,sz,zr)

    besorder = 1._real64
    call precomp_intval_adaptive(refl_var,besorder,funcBz1TE,10,sz,zr)

    !VTI: derivatives for epsv
    if (aniso.ne.iso) then
      besorder = 0._real64
      call precomp_intval_adaptive(refl_var,besorder,funcB0TMv,13,sz,zr)
    endif
    
    !do spline interpolation here and not inside precomp_intval_adaptive
    ! so that there is no need to go back and forth between radii and their logarithms

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,1),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,1))
    call spline(refl_var%radlog,refl_var%intvalim(:,1),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,1))
    call spline(refl_var%radlog,refl_var%intvalre(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,2))
    call spline(refl_var%radlog,refl_var%intvalim(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,2))
    call spline(refl_var%radlog,refl_var%intvalre(:,6),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,6))
    call spline(refl_var%radlog,refl_var%intvalim(:,6),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,6))
    call spline(refl_var%radlog,refl_var%intvalre(:,10),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,10))
    call spline(refl_var%radlog,refl_var%intvalim(:,10),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,10))

    if (aniso.ne.iso) then
      call spline(refl_var%radlog,refl_var%intvalre(:,13),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,13))
      call spline(refl_var%radlog,refl_var%intvalim(:,13),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,13))
    endif

    !compute well-behaved integrals by fast Hankel transform
    !related integrals C_TMderiv
    ibesord(1:3:2) = 1
    ibesord(2) = 0
    ijrel = 1
    ijrel(1,3) = 2 !factor kappa^2 for third integral
    iint(1) = 9    !integral iabvC1TMderiv
    iint(2) = 7    !integral iabvC0TMderiv
    iint(3) = 5    !integral iabvCz1TMderiv (for Ez)
    call precomp_intval_fht_rel(refl_var,3,ibesord,ijrel,funcC1TM,iint)

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,funcB1TE,3)
    call precomp_intval_fht(refl_var,ibesord,funcB1TM,4)

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,funcC1TE,8)
    
    if (ilaym .eq. ilayrec) then
      ibesord = 1
      call precomp_intval_fht(refl_var,ibesord,funcCz1TMfwd,11)
    endif

    if (aniso.ne.iso) then
    !related integrals C_TMderivv
      ibesord(1:3:2) = 1
      ibesord(2) = 0
      ijrel = 1
      ijrel(1,3) = 2 !factor kappa^2 for third integral
      iint(1) = 14    !integral iabvC1TMderivv
      iint(2) = 15    !integral iabvC0TMderivv
      iint(3) = 16    !integral iabvCz1TMderivv (for Ez)
      call precomp_intval_fht_rel(refl_var,3,ibesord,ijrel,funcC1TMv,iint)
      
      ibesord = 1
      call precomp_intval_fht(refl_var,ibesord,funcB1TMv,12)
    endif

  endif

endsubroutine precomp_intvals_deriv_hmd


!****************************************************************
!  1D EM subroutine precomp_intvals_deriv_hmd_Exy
!    compute integral values for derivatives for all integrals for HMD source, Ex/Ey only
!
!  Rita Streich 2010
!****************************************************************
subroutine precomp_intvals_deriv_hmd_Exy(refl_var,sz,zr,ilay, funcB0TE,funcB0TM,funcB1TE,funcB1TM, funcB0TMv,funcB1TMv)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and receiver depth
  integer(kind=int32),intent(in) :: ilay     !layer index, compute derivatives with respect to this layer
  complex(kind=real64),external   :: funcB0TE,funcB0TM,funcB1TE,funcB1TM
  complex(kind=real64),external,optional   :: funcB0TMv,funcB1TMv !functions for epsv derivatives

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter
  integer(kind=int32),dimension(nrelmax)  :: iint  !indicates locations for integration results in output structure


  !ilaym is defined in medium_mod, visible nearly everywhere
  ilaym = ilay

  !receivers above or below source, not exactly at the same depth
  if (sz .ne. zr) then

    !related integrals: IB1TEderiv, IB0TEderiv, IBz1TEderiv, same integrand just different kappa factors and Bessel function orders
    !basis is IB1TEderiv, contains no kappa, Bessel order 1
    !IB0TEderiv contains kappa^1, Bessel order 0
    !IBz1TEderiv contains kappa^2, Bessel order 1
    ibesord = 1
    ibesord(2) = 0
    ijrel = 1
    iint(1) = 3    !integral IabvB1TEderiv
    iint(2) = 1    !integral IabvB0TEderiv
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcB1TE,iint)

    !related integrals IB1TMderiv and IB0TMderiv
    !Bessel order as before (using first 2 elements only)
    !exponents for first 2 related transforms also as before
    iint(1) = 4    !integral IabvB1TMderiv
    iint(2) = 2    !integral IabvB0TMderiv
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcB1TM,iint)

    !VTI: derivatives for epsv
    if (aniso.ne.iso) then
      ibesord(2) = 0
      !related integrals IB1TMderivv and IB0TMderivv
      iint(1) = 12    !integral Iabv/blwB1TMderivv
      iint(2) = 13    !integral Iabv/blwB0TMderivv
      call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcB1TMv,iint)
    endif

  !receiver exactly at source depth
  else

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,funcB0TE,1,sz,zr)
    call precomp_intval_adaptive(refl_var,besorder,funcB0TM,2,sz,zr)

    !VTI: derivatives for epsv
    if (aniso.ne.iso) then
      besorder = 0._real64
      call precomp_intval_adaptive(refl_var,besorder,funcB0TMv,13,sz,zr)
    endif
    
    !do spline interpolation here and not inside precomp_intval_adaptive
    ! so that there is no need to go back and forth between radii and their logarithms

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,1),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,1))
    call spline(refl_var%radlog,refl_var%intvalim(:,1),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,1))
    call spline(refl_var%radlog,refl_var%intvalre(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,2))
    call spline(refl_var%radlog,refl_var%intvalim(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,2))

    if (aniso.ne.iso) then
      call spline(refl_var%radlog,refl_var%intvalre(:,13),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,13))
      call spline(refl_var%radlog,refl_var%intvalim(:,13),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,13))
    endif

    !compute well-behaved integrals by fast Hankel transform
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,funcB1TE,3)
    call precomp_intval_fht(refl_var,ibesord,funcB1TM,4)

    if (aniso.ne.iso) then
      ibesord = 1
      call precomp_intval_fht(refl_var,ibesord,funcB1TMv,12)
    endif

  endif

endsubroutine precomp_intvals_deriv_hmd_Exy


!****************************************************************
!  1D EM subroutine precomp_intvals_deriv_hmd_Ez
!    compute integral values for derivatives for all integrals for HMD source, Ez only
!
!  Rita Streich 2010
!****************************************************************
subroutine precomp_intvals_deriv_hmd_Ez(refl_var,sz,zr,ilay, funcCz1TM,funcCz1TMfwd, funcCz1TMv)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and receiver depth
  integer(kind=int32),intent(in) :: ilay     !layer index, compute derivatives with respect to this layer
  complex(kind=real64),external   :: funcCz1TM
  complex(kind=real64),external   :: funcCz1TMfwd  !function only needed for derivatives in receiver layer
  complex(kind=real64),external,optional   :: funcCz1TMv !functions for epsv derivatives

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")

  !ilaym is defined in medium_mod, visible nearly everywhere
  ilaym = ilay

  ibesord = 1

  !integrals are (probabbly) well-behaved, no need for adaptive integration --> no need to distinguish case sz == zr
  !receivers above or below source, not exactly at the same depth

    call precomp_intval_fht(refl_var,ibesord,funcCz1TM,5)  !integral iabvCz1TMderiv (for Ez)

    !special in receiver layer: there is an eps in front of the integral, so we also need 
    ! derivative of term in front of integral times integral without derivative
    if (ilaym .eq. ilayrec) then
      call precomp_intval_fht(refl_var,ibesord,funcCz1TMfwd,11)
    endif
    
    !VTI: derivatives for epsv
    if (aniso.ne.iso) then
      call precomp_intval_fht(refl_var,ibesord,funcCz1TMv,16)
    endif

endsubroutine precomp_intvals_deriv_hmd_Ez


!****************************************************************
!  1D EM subroutine precomp_intvals_deriv_hmd_Hxy
!    compute integral values for derivatives for all integrals for HMD source, Hx/Hy only
!
!  Rita Streich 2010
!****************************************************************
subroutine precomp_intvals_deriv_hmd_Hxy(refl_var,sz,zr,ilay, funcC0TE,funcC0TM,funcC1TE,funcC1TM, funcC0TMv,funcC1TMv)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and receiver depth
  integer(kind=int32),intent(in) :: ilay     !layer index, compute derivatives with respect to this layer
  complex(kind=real64),external   :: funcC0TE,funcC0TM,funcC1TE,funcC1TM
  complex(kind=real64),external,optional   :: funcC0TMv,funcC1TMv !functions for epsv derivatives

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter
  integer(kind=int32),dimension(nrelmax)  :: iint  !indicates locations for integration results in output structure


  !ilaym is defined in medium_mod, visible nearly everywhere
  ilaym = ilay

  !receivers above or below source, not exactly at the same depth
  if (sz .ne. zr) then

    ibesord = 1
    ibesord(2) = 0
    ijrel = 1

    !related integrals IC1TEderiv and IC0TEderiv
    iint(1) = 8    !integral IabvC1TEderiv
    iint(2) = 6    !integral IabvC0TEderiv
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcC1TE,iint)

    !related integrals C_TM
    iint(1) = 9    !integral iabvC1TMderiv
    iint(2) = 7    !integral iabvC0TMderiv
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcC1TM,iint)

    !VTI: derivatives for epsv
    if (aniso.ne.iso) then
      !related integrals C_TMv
      iint(1) = 14    !integral iabv/blwC1TMderivv
      iint(2) = 15    !integral iabv/blwC0TMderivv
      call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcC1TMv,iint)
    endif

  !receiver exactly at source depth
  else

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,funcC0TE,6,sz,zr)

    !do spline interpolation here and not inside precomp_intval_adaptive
    ! so that there is no need to go back and forth between radii and their logarithms

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,6),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,6))
    call spline(refl_var%radlog,refl_var%intvalim(:,6),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,6))

    !compute well-behaved integrals by fast Hankel transform
    !related integrals C_TMderiv
    ibesord = 1
    ibesord(2) = 0
    ijrel = 1
    iint(1) = 9    !integral iabvC1TMderiv
    iint(2) = 7    !integral iabvC0TMderiv
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcC1TM,iint)

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,funcC1TE,8)
    
    if (aniso.ne.iso) then
    !related integrals C_TMderivv
      ibesord = 1
      ibesord(2) = 0
      iint(1) = 14    !integral iabvC1TMderivv
      iint(2) = 15    !integral iabvC0TMderivv
      call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcC1TMv,iint)
    endif

  endif

endsubroutine precomp_intvals_deriv_hmd_Hxy


!****************************************************************
!  1D EM subroutine precomp_intvals_deriv_hmd_Hz
!    compute integral values for derivatives for all integrals for HMD source, Hz only
!
!  Rita Streich 2010
!****************************************************************
subroutine precomp_intvals_deriv_hmd_Hz(refl_var,sz,zr,ilay, funcBz1TE)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and receiver depth
  integer(kind=int32),intent(in) :: ilay     !layer index, compute derivatives with respect to this layer
  complex(kind=real64),external   :: funcBz1TE

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter


  !ilaym is defined in medium_mod, visible nearly everywhere
  ilaym = ilay

  !receivers above or below source, not exactly at the same depth
  if (sz .ne. zr) then

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,funcBz1TE,10) !integral IabvBz1TEderiv (for Hz)

  !receiver exactly at source depth
  else

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 1._real64
    call precomp_intval_adaptive(refl_var,besorder,funcBz1TE,10,sz,zr)

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,10),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,10))
    call spline(refl_var%radlog,refl_var%intvalim(:,10),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,10))

  endif

endsubroutine precomp_intvals_deriv_hmd_Hz


!****************************************************************
!  1D EM subroutine precomp_intvals_deriv_vmd
!    compute integral values for derivatives for all integrals for VMD source
!
!  Rita Streich 2010
!****************************************************************
subroutine precomp_intvals_deriv_vmd(refl_var,sz,zr,ilay, funcA1TEvmd,funcD1TEvmd,funcA0TEvmd)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and receiver depth
  integer(kind=int32),intent(in) :: ilay     !layer index, compute derivatives with respect to this layer
  complex(kind=real64),external  :: funcA1TEvmd,funcD1TEvmd,funcA0TEvmd !integral derivative functions

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter
  integer(kind=int32),dimension(nrelmax)  :: iint  !indicates locations for integration results in output structure


  !ilaym is defined in medium_mod, visible nearly everywhere
  ilaym = ilay

  !receivers above or below source, not at exactly the same depth
  if (sz .ne. zr) then

    !related integrals: iabvAz1TEderiv and iabvA0TEvmdderiv
    ibesord(1) = 1
    ibesord(2) = 0
    ijrel = 1
    iint(1) = 1    !integral iabv/blwAz1TEderiv (for Ex and Ey)
    iint(2) = 3    !integral iabv/blwA0TEvmdderiv (for Hz)
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,funcA1TEvmd,iint)

    !no integral for Ez - Ez for VMD source is zero

    !integral for Hx and Hy
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,funcD1TEvmd,2)

  !receiver exactly at source depth
  else

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 1._real64
    call precomp_intval_adaptive(refl_var,besorder,funcD1TEvmd,2,sz,zr)

    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,funcA0TEvmd,3,sz,zr)

    !do spline interpolation here and not inside precomp_intval_adaptive
    ! so that there is no need to go back and forth between radii and their logarithms

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,2))
    call spline(refl_var%radlog,refl_var%intvalim(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,2))
    call spline(refl_var%radlog,refl_var%intvalre(:,3),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,3))
    call spline(refl_var%radlog,refl_var%intvalim(:,3),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,3))


    !compute well-behaved integrals by fast Hankel transform
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,funcA1TEvmd,1)
  endif

endsubroutine precomp_intvals_deriv_vmd


!****************************************************************
!  1D EM subroutine precomp_intvals_deriv_vmd_Exy
!    compute integral values for derivatives for all integrals for VMD source, Ex/Ey only
!
!  Rita Streich 2010
!****************************************************************
subroutine precomp_intvals_deriv_vmd_Exy(refl_var,sz,zr,ilay, funcA1TEvmd)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and receiver depth
  integer(kind=int32),intent(in) :: ilay     !layer index, compute derivatives with respect to this layer
  complex(kind=real64),external  :: funcA1TEvmd !integral derivative functions

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")


  !ilaym is defined in medium_mod, visible nearly everywhere
  ilaym = ilay

  !integral is well-behaved, no need for special case sz == zr and adaptive integration

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,funcA1TEvmd,1)  !integral iabv/blwAz1TEderiv (for Ex and Ey)

endsubroutine precomp_intvals_deriv_vmd_Exy


!****************************************************************
!  1D EM subroutine precomp_intvals_deriv_vmd_Hxy
!    compute integral values for derivatives for all integrals for VMD source, Hx/Hy only
!
!  Rita Streich 2010
!****************************************************************
subroutine precomp_intvals_deriv_vmd_Hxy(refl_var,sz,zr,ilay, funcD1TEvmd)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and receiver depth
  integer(kind=int32),intent(in) :: ilay     !layer index, compute derivatives with respect to this layer
  complex(kind=real64),external  :: funcD1TEvmd !integral derivative functions

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter


  !ilaym is defined in medium_mod, visible nearly everywhere
  ilaym = ilay

  !receivers above or below source, not at exactly the same depth
  if (sz .ne. zr) then

    !integral for Hx and Hy
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,funcD1TEvmd,2)

  !receiver exactly at source depth
  else

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 1._real64
    call precomp_intval_adaptive(refl_var,besorder,funcD1TEvmd,2,sz,zr)

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,2))
    call spline(refl_var%radlog,refl_var%intvalim(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,2))
  endif

endsubroutine precomp_intvals_deriv_vmd_Hxy


!****************************************************************
!  1D EM subroutine precomp_intvals_deriv_vmd_Hz
!    compute integral values for derivatives for all integrals for VMD source, Hz only
!
!  Rita Streich 2010-2011
!****************************************************************
subroutine precomp_intvals_deriv_vmd_Hz(refl_var,sz,zr,ilay, funcA0TEvmd)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and receiver depth
  integer(kind=int32),intent(in) :: ilay     !layer index, compute derivatives with respect to this layer
  complex(kind=real64),external  :: funcA0TEvmd !integral derivative functions

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter


  !ilaym is defined in medium_mod, visible nearly everywhere
  ilaym = ilay

  !receivers above or below source, not at exactly the same depth
  if (sz .ne. zr) then

    ibesord = 0
    call precomp_intval_fht(refl_var,ibesord,funcA0TEvmd,3)  !integral iabv/blwA0TEvmdderiv (for Hz)

  !receiver exactly at source depth
  else

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,funcA0TEvmd,3,sz,zr)

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,3),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,3))
    call spline(refl_var%radlog,refl_var%intvalim(:,3),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,3))
  endif

endsubroutine precomp_intvals_deriv_vmd_Hz


!****************************************************************
!  1D EM subroutine precomp_intvals_deriv_wire
!    compute integral values for derivatives for all integrals along wires
!
!  Rita Streich 2010
!****************************************************************
subroutine precomp_intvals_deriv_wire(refl_var,sz,zr,ilay, &
  func1Exwire,funcD0TE,funcAz1TE,func2Exwire,funcD0TM,funcdHxwire,funcD0TMfwd, &
  func2Exwirev,funcD0TMv,funcdHxwirev)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and receiver depth
  integer(kind=int32),intent(in)  :: ilay     !layer index, compute derivatives with respect to this layer
  complex(kind=real64),external   :: func1Exwire,funcD0TE,funcAz1TE,func2Exwire,funcD0TM,funcdHxwire
  complex(kind=real64),external   :: funcD0TMfwd  !function only needed for derivatives in receiver layer
  complex(kind=real64),external,optional   :: func2Exwirev,funcD0TMv,funcdHxwirev !integral derivatives for epsv

  !internal variables
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter


  !ilaym is defined in medium_mod, visible nearly everywhere
  ilaym = ilay

  !receivers above or below source
  if (sz .ne. zr) then

    !integrals along wires
    !integral for Ex (and Ey), TE only
    ibesord = 0
    call precomp_intval_fht(refl_var,ibesord,func1Exwire,1)

    !integral for Hy (and Hx)
    ibesord = 0
    call precomp_intval_fht(refl_var,ibesord,funcD0TE,2)

    !integral for Hz
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,funcAz1TE,3)

    !integrals for end points
    !integral for Ex and Ey
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,func2Exwire,4)

    !integral for Ez
    ibesord = 0
    call precomp_intval_fht(refl_var,ibesord,funcD0TM,5)

    !integral for Hx and Hy
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,funcdHxwire,6)

    !special in receiver layer for Ez: there is an eps in front of the integral, so we also need 
    ! derivative of term in front of integral times integral without derivative
    if (ilaym .eq. ilayrec) then
      ibesord = 0
      call precomp_intval_fht(refl_var,ibesord,funcD0TMfwd,7)
    endif
    
    !VTI: epsv derivatives of integrals containing TM parts: end point integrals only!
    if (aniso.ne.iso) then
      ibesord = 1
      call precomp_intval_fht(refl_var,ibesord,func2Exwirev,8)
      ibesord = 0
      call precomp_intval_fht(refl_var,ibesord,funcD0TMv,9)
      ibesord = 1
      call precomp_intval_fht(refl_var,ibesord,funcdHxwirev,10)
    endif

  !receiver exactly at source depth
  else

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    !integrals along wire
    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,func1Exwire,1,sz,zr)
    call precomp_intval_adaptive(refl_var,besorder,funcD0TE,2,sz,zr)
    !integrals for end points
    besorder = 1._real64
    call precomp_intval_adaptive(refl_var,besorder,func2Exwire,4,sz,zr)
    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,funcD0TM,5,sz,zr)

    !special term for derivative of Ez in receiver layer
    if (ilaym .eq. ilayrec) then
      besorder = 0._real64
      call precomp_intval_adaptive(refl_var,besorder,funcD0TMfwd,7,sz,zr)
    endif

    if (aniso.ne.iso) then
      besorder = 1._real64
      call precomp_intval_adaptive(refl_var,besorder,func2Exwirev,8,sz,zr)
      besorder = 0._real64
      call precomp_intval_adaptive(refl_var,besorder,funcD0TMv,9,sz,zr)
    endif

    !get spline derivatives

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,1),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,1))
    call spline(refl_var%radlog,refl_var%intvalim(:,1),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,1))
    call spline(refl_var%radlog,refl_var%intvalre(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,2))
    call spline(refl_var%radlog,refl_var%intvalim(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,2))
    call spline(refl_var%radlog,refl_var%intvalre(:,4),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,4))
    call spline(refl_var%radlog,refl_var%intvalim(:,4),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,4))
    call spline(refl_var%radlog,refl_var%intvalre(:,5),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,5))
    call spline(refl_var%radlog,refl_var%intvalim(:,5),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,5))

    !special term for derivative of Ez in receiver layer
    if (ilaym .eq. ilayrec) then
      call spline(refl_var%radlog,refl_var%intvalre(:,7),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,7))
      call spline(refl_var%radlog,refl_var%intvalim(:,7),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,7))
    endif
    
    if (aniso.ne.iso) then
      call spline(refl_var%radlog,refl_var%intvalre(:,8),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,8))
      call spline(refl_var%radlog,refl_var%intvalim(:,8),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,8))
      call spline(refl_var%radlog,refl_var%intvalre(:,9),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,9))
      call spline(refl_var%radlog,refl_var%intvalim(:,9),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,9))

      ibesord = 1
      call precomp_intval_fht(refl_var,ibesord,funcdHxwirev,10)
    endif

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,funcAz1TE,3)
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,funcdHxwire,6)

  endif

endsubroutine precomp_intvals_deriv_wire


!****************************************************************
!  1D EM subroutine precomp_intvals_deriv_wire_Exy
!    compute integral values for derivatives for all integrals along wires, Ex/Ey only
!
!  Rita Streich 2011
!****************************************************************
subroutine precomp_intvals_deriv_wire_Exy(refl_var,sz,zr,ilay, func1Exwire,func2Exwire, func2Exwirev)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and receiver depth
  integer(kind=int32),intent(in)  :: ilay     !layer index, compute derivatives with respect to this layer
  complex(kind=real64),external   :: func1Exwire,func2Exwire
  complex(kind=real64),external,optional   :: func2Exwirev !integral derivatives for epsv

  !internal variables
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter


  !ilaym is defined in medium_mod, visible nearly everywhere
  ilaym = ilay

  !receivers above or below source
  if (sz .ne. zr) then

    !integrals along wires
    !integral for Ex (and Ey), TE only
    ibesord = 0
    call precomp_intval_fht(refl_var,ibesord,func1Exwire,1)

    !integrals for end points
    !integral for Ex and Ey
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,func2Exwire,4)

    !VTI: epsv derivatives of integrals containing TM parts: end point integrals only!
    if (aniso.ne.iso) then
      ibesord = 1
      call precomp_intval_fht(refl_var,ibesord,func2Exwirev,8)
    endif

  !receiver exactly at source depth
  else

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    !integrals along wire
    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,func1Exwire,1,sz,zr)
    !integrals for end points
    besorder = 1._real64
    call precomp_intval_adaptive(refl_var,besorder,func2Exwire,4,sz,zr)

    if (aniso.ne.iso) then
      besorder = 1._real64
      call precomp_intval_adaptive(refl_var,besorder,func2Exwirev,8,sz,zr)
    endif

    !get spline derivatives

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,1),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,1))
    call spline(refl_var%radlog,refl_var%intvalim(:,1),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,1))
    call spline(refl_var%radlog,refl_var%intvalre(:,4),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,4))
    call spline(refl_var%radlog,refl_var%intvalim(:,4),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,4))

    if (aniso.ne.iso) then
      call spline(refl_var%radlog,refl_var%intvalre(:,8),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,8))
      call spline(refl_var%radlog,refl_var%intvalim(:,8),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,8))
    endif

  endif

endsubroutine precomp_intvals_deriv_wire_Exy


!****************************************************************
!  1D EM subroutine precomp_intvals_deriv_wire_Ez
!    compute integral values for derivatives for all integrals along wires, Ez only
!
!  Rita Streich 2011
!****************************************************************
subroutine precomp_intvals_deriv_wire_Ez(refl_var,sz,zr,ilay, funcD0TM,funcD0TMfwd, funcD0TMv)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and receiver depth
  integer(kind=int32),intent(in)  :: ilay     !layer index, compute derivatives with respect to this layer
  complex(kind=real64),external   :: funcD0TM
  complex(kind=real64),external   :: funcD0TMfwd  !function only needed for derivatives in receiver layer
  complex(kind=real64),external,optional   :: funcD0TMv !integral derivatives for epsv

  !internal variables
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter


  !ilaym is defined in medium_mod, visible nearly everywhere
  ilaym = ilay

  !receivers above or below source
  if (sz .ne. zr) then

    !integrals for end points
    !integral for Ez
    ibesord = 0
    call precomp_intval_fht(refl_var,ibesord,funcD0TM,5)

    !special in receiver layer for Ez: there is an eps in front of the integral, so we also need 
    ! derivative of term in front of integral times integral without derivative
    if (ilaym .eq. ilayrec) then
      ibesord = 0
      call precomp_intval_fht(refl_var,ibesord,funcD0TMfwd,7)
    endif
    
    !VTI: epsv derivatives of integrals containing TM parts: end point integrals only!
    if (aniso.ne.iso) then
      ibesord = 0
      call precomp_intval_fht(refl_var,ibesord,funcD0TMv,9)
    endif

  !receiver exactly at source depth
  else

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    !integrals along wire
    !integrals for end points
    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,funcD0TM,5,sz,zr)

    !special term for derivative of Ez in receiver layer
    if (ilaym .eq. ilayrec) then
      besorder = 0._real64
      call precomp_intval_adaptive(refl_var,besorder,funcD0TMfwd,7,sz,zr)
    endif

    if (aniso.ne.iso) then
      besorder = 0._real64
      call precomp_intval_adaptive(refl_var,besorder,funcD0TMv,9,sz,zr)
    endif

    !get spline derivatives

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,5),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,5))
    call spline(refl_var%radlog,refl_var%intvalim(:,5),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,5))

    !special term for derivative of Ez in receiver layer
    if (ilaym .eq. ilayrec) then
      call spline(refl_var%radlog,refl_var%intvalre(:,7),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,7))
      call spline(refl_var%radlog,refl_var%intvalim(:,7),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,7))
    endif
    
    if (aniso.ne.iso) then
      call spline(refl_var%radlog,refl_var%intvalre(:,9),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,9))
      call spline(refl_var%radlog,refl_var%intvalim(:,9),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,9))
    endif

  endif

endsubroutine precomp_intvals_deriv_wire_Ez


!****************************************************************
!  1D EM subroutine precomp_intvals_deriv_wire_Hxy
!    compute integral values for derivatives for all integrals along wires, Hx/Hy only
!
!  Rita Streich 2011
!****************************************************************
subroutine precomp_intvals_deriv_wire_Hxy(refl_var,sz,zr,ilay, funcD0TE,funcdHxwire, funcdHxwirev)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and receiver depth
  integer(kind=int32),intent(in)  :: ilay     !layer index, compute derivatives with respect to this layer
  complex(kind=real64),external   :: funcD0TE,funcdHxwire
  complex(kind=real64),external,optional   :: funcdHxwirev !integral derivatives for epsv

  !internal variables
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter


  !ilaym is defined in medium_mod, visible nearly everywhere
  ilaym = ilay

  !receivers above or below source
  if (sz .ne. zr) then

    !integrals along wires
    !integral for Hy (and Hx)
    ibesord = 0
    call precomp_intval_fht(refl_var,ibesord,funcD0TE,2)

    !integrals for end points
    !integral for Hx and Hy
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,funcdHxwire,6)

    !VTI: epsv derivatives of integrals containing TM parts: end point integrals only!
    if (aniso.ne.iso) then
      ibesord = 1
      call precomp_intval_fht(refl_var,ibesord,funcdHxwirev,10)
    endif

  !receiver exactly at source depth
  else

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    !integrals along wire
    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,funcD0TE,2,sz,zr)

    !get spline derivatives

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,2))
    call spline(refl_var%radlog,refl_var%intvalim(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,2))
    
    if (aniso.ne.iso) then
      ibesord = 1
      call precomp_intval_fht(refl_var,ibesord,funcdHxwirev,10)
    endif

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,funcdHxwire,6)

  endif

endsubroutine precomp_intvals_deriv_wire_Hxy


!****************************************************************
!  1D EM subroutine precomp_intvals_deriv_wire_Hz
!    compute integral values for derivatives for all integrals along wires, Hz only
!
!  Rita Streich 2010
!****************************************************************
subroutine precomp_intvals_deriv_wire_Hz(refl_var,sz,zr,ilay, funcAz1TE)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and receiver depth
  integer(kind=int32),intent(in)  :: ilay     !layer index, compute derivatives with respect to this layer
  complex(kind=real64),external   :: funcAz1TE

  !internal variables
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")


  !ilaym is defined in medium_mod, visible nearly everywhere
  ilaym = ilay

  !no need for special case for receiver exactly at source depth
  !integral along wire for Hz
  ibesord = 1
  call precomp_intval_fht(refl_var,ibesord,funcAz1TE,3)

endsubroutine precomp_intvals_deriv_wire_Hz


