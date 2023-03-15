!****************************************************************
!  1D EM subroutine precomp_intval
!    compute integral values for one integral for logarithmically spaced radii
!    by fast Hankel transform
!
!  source and iReceiver depths are not needed here since they have already been absorbed
!    into dz vectors
!
!  Rita Streich 2009
!****************************************************************
subroutine precomp_intval_fht(refl_var,ibesord,intfunc,iint)

  implicit none

  !external variables
  type(refl_struct)               :: refl_var  !all stuff for 1D computations
  integer(kind=int32),dimension(1),intent(in)  :: ibesord   !bessel function order
  complex(kind=real64),external   :: intfunc   !integral function (for receivers above or below source)
  integer(kind=int32),intent(in)  :: iint      !where to write output

  !internal variables
  integer(kind=int32)             :: irad      !radius counter
  integer(kind=int32)             :: ierr      !error index
  integer(kind=int32)             :: NOFUN1    !number of function evaluations


  !use temp array for integral values because it has to be a 2D array
  call zhankl(refl_var%rmax,refl_var%nrad,NREL,TOL,NTOL,ibesord,intfunc,IJREL,ZWORK,refl_var%intvaltmp, &
              refl_var%radlog,NOFUN1,ierr)

  !take logarithm of radii - this makes spline interpolation more accurate
  do irad=1,refl_var%nrad
    refl_var%radlog(irad) = log10(refl_var%radlog(irad))
  enddo

  refl_var%intvalre(:,iint) = real(refl_var%intvaltmp(:,1))
  refl_var%intvalim(:,iint) = aimag(refl_var%intvaltmp(:,1))

  !spline derivatives
  call spline(refl_var%radlog,refl_var%intvalre(:,iint),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,iint))
  call spline(refl_var%radlog,refl_var%intvalim(:,iint),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,iint))

endsubroutine precomp_intval_fht


!****************************************************************
!  1D EM subroutine precomp_intval_fht_rel
!    compute integral values for one integral for logarithmically spaced radii
!    by fast Hankel transform - expoint related Hankel integrals!
!
!  Rita Streich 2010
!****************************************************************
subroutine precomp_intval_fht_rel(refl_var,nrel,ibesord,ijrel,intfunc,iint)

  implicit none

  !external variables
  type(refl_struct)               :: refl_var  !all stuff for 1D computations
  integer(kind=int32),intent(in)  :: nrel      !nr of related Hankel integrals
  integer(kind=int32),dimension(nrel),intent(in)  :: ibesord   !bessel function order
  integer(kind=int32),dimension(2,nrel),intent(in)  :: ijrel   !parameters for related Hankel integrals
  complex(kind=real64),external   :: intfunc   !integral function (for receivers above or below source)
  integer(kind=int32),dimension(nrel),intent(in)  :: iint      !where to write output

  !internal variables
  integer(kind=int32)             :: irad      !radius counter
  integer(kind=int32)             :: ierr      !error index
  integer(kind=int32)             :: NOFUN1    !number of function evaluations
  integer(kind=int32)             :: irel      !counter for related integrals
  integer(kind=int32)             :: iout      !output data index


  !use temp array for integral values because it has to be a 2D array
  call zhankl(refl_var%rmax,refl_var%nrad,NREL,TOL,NTOL,ibesord,intfunc,IJREL,ZWORK,refl_var%intvaltmp, &
              refl_var%radlog,NOFUN1,ierr)

  !take logarithm of radii - this makes spline interpolation more accurate
  do irad=1,refl_var%nrad
    refl_var%radlog(irad) = log10(refl_var%radlog(irad))
  enddo

  do irel=1,nrel
    iout = iint(irel)

    refl_var%intvalre(:,iout) = real(refl_var%intvaltmp(:,irel))
    refl_var%intvalim(:,iout) = aimag(refl_var%intvaltmp(:,irel))

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,iout),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,iout))
    call spline(refl_var%radlog,refl_var%intvalim(:,iout),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,iout))
  enddo

endsubroutine precomp_intval_fht_rel


!****************************************************************
!  1D EM subroutine precomp_intval_adaptive
!    compute integral values for one integral for logarithmically spaced radii
!    by adaptive integration
!
!  Rita Streich 2009
!****************************************************************
subroutine precomp_intval_adaptive(refl_var,besorder,intfunc,iint,sz,zr)

  implicit none

  !external variables
  type(refl_struct)               :: refl_var  !all stuff for 1D computations
  real(kind=real64),intent(in)    :: besorder  !bessel function order
  complex(kind=real64),external   :: intfunc   !integral kernel function
  integer(kind=int32),intent(in)  :: iint      !where to write output
  real(kind=real64),intent(in)    :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32)     :: irad     !counter for radii
  integer(kind=int32)     :: npieces  !number of sub-intervals for adaptive integration
  real(kind=real64)       :: r        !temp radius
  integer(kind=int32)     :: newint   !switch for adpative integration, no previous kernels are re-used
  integer(kind=int32)     :: ierr     !error index


  !compute integral value for each radius separately
  do irad = 1,refl_var%nrad

    r = refl_var%radlog(irad)
    npieces = getnpieces(r,sz,zr)

    ikap = -1 !indicate that we have to recompute all interface refl./transm. coeff.
    newint = 1
    CALL BESAUT(refl_var%intvalre(irad,iint),refl_var%intvalim(irad,iint),besorder,GAUSLO,GAUSHI,r,intfunc, &
                relerr,abserr,npieces,newint,0,branchpt,ierr)

  enddo

endsubroutine precomp_intval_adaptive



!****************************************************************
!  1D EM subroutine precomp_intvals_hed
!    compute integral values for all integrals for HED source
!
!  Rita Streich 2009
!****************************************************************
subroutine precomp_intvals_hed(refl_var,sz,zr)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter
  integer(kind=int32),dimension(nrelmax)  :: iint  !indicates locations for integration results in output structure

  !receivers above source
  if (sz .gt. zr) then


    !related integrals: IA1TE, IA0TE, IAz1TE, same integrand just different kappa factors and Bessel function orders
    !basis is IA1TE, contains no kappa, Bessel order 1
    !IA0TE contains kappa^1, Bessel order 0
    !IAz1TE contains kappa^2, Bessel order 1
    ibesord(1:3:2) = 1
    ibesord(2) = 0
    ijrel = 1
    ijrel(1,3) = 2 !factor kappa^2 for third integral
    iint(1) = 3    !integral IabvA1TE
    iint(2) = 1    !integral IabvA0TE
    iint(3) = 10   !integral IabvAz1TE (for Hz)

    call precomp_intval_fht_rel(refl_var,3,ibesord,ijrel,iabvA1TE,iint)

    !related integrals IA1TM and IA0TM
    !Bessel order as before (using first 2 elements only)
    !exponents for first 2 related transforms also as before
    iint(1) = 4    !integral IabvA1TM
    iint(2) = 2    !integral IabvA0TM
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iabvA1TM,iint)

    !related integrals ID1TE and ID0TE
    iint(1) = 8    !integral IabvD1TE
    iint(2) = 6    !integral IabvD0TE
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iabvD1TE,iint)

    !related integrals D_TM
    iint(1) = 9    !integral iabvD1TM
    iint(2) = 7    !integral iabvD0TM
    iint(3) = 5    !integral iabvDz1TM (for Ez)
    call precomp_intval_fht_rel(refl_var,3,ibesord,ijrel,iabvD1TM,iint)
 
  !iReceiver below source
  elseif (sz .lt. zr) then

    ibesord(1:3:2) = 1
    ibesord(2) = 0
    ijrel = 1
    ijrel(1,3) = 2 !factor kappa^2 for third integral
    iint(1) = 3    !integral IblwA1TE
    iint(2) = 1    !integral IblwA0TE
    iint(3) = 10   !integral IblwAz1TE (for Hz)

    call precomp_intval_fht_rel(refl_var,3,ibesord,ijrel,iblwA1TE,iint)

    iint(1) = 4    !integral IblwA1TM
    iint(2) = 2    !integral IblwA0TM
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iblwA1TM,iint)

    iint(1) = 8    !integral IblwD1TE
    iint(2) = 6    !integral IblwD0TE
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iblwD1TE,iint)

    iint(1) = 9    !integral iblwD1TM
    iint(2) = 7    !integral iblwD0TM
    iint(3) = 5    !integral iblwDz1TM (for Ez)
    call precomp_intval_fht_rel(refl_var,3,ibesord,ijrel,iblwD1TM,iint)

  !iReceiver exactly at source depth
  else
    if (refl_var%infolevel.ge.output_more) then
      write(*,'(a)') 'WARNING: source depth = iReceiver depth, entering adaptive integration, this can be SLOW!'
    endif

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,iabvA0TM,2,sz,zr)
    besorder = 1._real64
    call precomp_intval_adaptive(refl_var,besorder,iabvDz1TM,5,sz,zr)

    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,iabvD0TE,6,sz,zr)
    call precomp_intval_adaptive(refl_var,besorder,iabvD0TM,7,sz,zr)

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


    !compute well-behaved integrals by fast Hankel transform
    !related integrals: IA1TE, IA0TE, IAz1TE
    ibesord(1:3:2) = 1
    ibesord(2) = 0
    ijrel = 1
    ijrel(1,3) = 2 !factor kappa^2 for third integral
    iint(1) = 3    !integral IabvA1TE
    iint(2) = 1    !integral IabvA0TE
    iint(3) = 10   !integral IabvAz1TE (for Hz)
    call precomp_intval_fht_rel(refl_var,3,ibesord,ijrel,iabvA1TE,iint)

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iabvA1TM,4)

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iabvD1TE,8)
    call precomp_intval_fht(refl_var,ibesord,iabvD1TM,9)

  endif

endsubroutine precomp_intvals_hed


!****************************************************************
!  FD EM subroutine precomp_intvals_Exy_hed
!    compute integral values for integrals needed for Ex and Ey, HED source
!
!  Rita Streich 2009
!****************************************************************
subroutine precomp_intvals_Exy_hed(refl_var,sz,zr)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter
  integer(kind=int32),dimension(nrelmax)  :: iint  !indicates locations for integration results in output structure


  !receivers above source
  if (sz .gt. zr) then

    !related integrals: IA1TE, IA0TE, same integrand just different kappa factors and Bessel function orders
    !basis is IA1TE, contains no kappa, Bessel order 1
    !IA0TE contains kappa^1, Bessel order 0
    ibesord = 1
    ibesord(2) = 0
    ijrel = 1
    iint(1) = 3    !integral IabvA1TE - 3rd in refl_var structure
    iint(2) = 1    !integral IabvA0TE

    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iabvA1TE,iint)

    !related integrals IA1TM and IA0TM
    !Bessel order as before
    !exponents for first 2 related transforms also as before
    iint(1) = 4    !integral IabvA1TM
    iint(2) = 2    !integral IabvA0TM
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iabvA1TM,iint)

  !iReceiver below source
  elseif (sz .lt. zr) then

    ibesord = 1
    ibesord(2) = 0
    ijrel = 1
    iint(1) = 3    !integral IblwA1TE
    iint(2) = 1    !integral IblwA0TE

    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iblwA1TE,iint)

    iint(1) = 4    !integral IblwA1TM
    iint(2) = 2    !integral IblwA0TM
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iblwA1TM,iint)

  !iReceiver exactly at source depth
  else

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,iabvA0TM,2,sz,zr)

    !do spline interpolation here and not inside precomp_intval_adaptive
    ! so that there is no need to go back and forth between radii and their logarithms

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,2))
    call spline(refl_var%radlog,refl_var%intvalim(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,2))


    !compute well-behaved integrals by fast Hankel transform
    !related integrals: IA1TE, IA0TE
    ibesord(1) = 1
    ibesord(2) = 0
    ijrel = 1
    iint(1) = 3    !integral IabvA1TE
    iint(2) = 1    !integral IabvA0TE
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iabvA1TE,iint)

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iabvA1TM,4)

  endif

endsubroutine precomp_intvals_Exy_hed


!****************************************************************
!  FD EM subroutine precomp_intvals_Ez_hed
!    compute integral values for the integral needed for Ez, HED source
!
!  Rita Streich 2009
!****************************************************************
subroutine precomp_intvals_Ez_hed(refl_var,sz,zr)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter


  !receivers above source
  if (sz .gt. zr) then

    !integral for Ez
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iabvDz1TM,5)

  !iReceiver below source
  elseif (sz .lt. zr) then

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iblwDz1TM,5)

  !iReceiver exactly at source depth
  else

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 1._real64
    call precomp_intval_adaptive(refl_var,besorder,iabvDz1TM,5,sz,zr)

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,5),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,5))
    call spline(refl_var%radlog,refl_var%intvalim(:,5),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,5))

  endif

endsubroutine precomp_intvals_Ez_hed


!****************************************************************
!  FD EM subroutine precomp_intvals_Hxy_hed
!    compute integral values for all integrals for HED source, Hx and Hy
!
!  Rita Streich 2011
!****************************************************************
subroutine precomp_intvals_Hxy_hed(refl_var,sz,zr)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter
  integer(kind=int32),dimension(nrelmax)  :: iint  !indicates locations for integration results in output structure


  !receivers above source
  if (sz .gt. zr) then

    !related integrals: IA1TE, IA0TE, IAz1TE, same integrand just different kappa factors and Bessel function orders
    !basis is IA1TE, contains no kappa, Bessel order 1
    !IA0TE contains kappa^1, Bessel order 0
    !IAz1TE contains kappa^2, Bessel order 1
    ibesord = 1
    ibesord(2) = 0
    ijrel = 1

    !related integrals ID1TE and ID0TE
    iint(1) = 8    !integral IabvD1TE
    iint(2) = 6    !integral IabvD0TE
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iabvD1TE,iint)

    !related integrals D_TM
    iint(1) = 9    !integral iabvD1TM
    iint(2) = 7    !integral iabvD0TM
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iabvD1TM,iint)

  !iReceiver below source
  elseif (sz .lt. zr) then

    ibesord = 1
    ibesord(2) = 0
    ijrel = 1

    iint(1) = 8    !integral IblwD1TE
    iint(2) = 6    !integral IblwD0TE
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iblwD1TE,iint)

    iint(1) = 9    !integral iblwD1TM
    iint(2) = 7    !integral iblwD0TM
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iblwD1TM,iint)

  !iReceiver exactly at source depth
  else

    if (refl_var%infolevel .ge. output_more) then
      write(*,'(a)') 'WARNING: source depth = iReceiver depth, entering adaptive integration, this can be SLOW!'
    endif

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,iabvD0TE,6,sz,zr)
    call precomp_intval_adaptive(refl_var,besorder,iabvD0TM,7,sz,zr)

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


    !compute well-behaved integrals by fast Hankel transform
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iabvD1TE,8)
    call precomp_intval_fht(refl_var,ibesord,iabvD1TM,9)

  endif

endsubroutine precomp_intvals_Hxy_hed


!****************************************************************
!  FD EM subroutine precomp_intvals_Hz_hed
!    compute integral values for all integrals for HED source, Hz only
!
!  Rita Streich 2011
!****************************************************************
subroutine precomp_intvals_Hz_hed(refl_var,sz,zr)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")


  !receivers above source
  !include case that iReceiver is exactly at source depth here: integral is "well-behaved", adaptive integration not needed
  if (sz .ge. zr) then
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iabvAz1TE,10)

  !iReceiver below source
  else
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iblwAz1TE,10)
  endif

endsubroutine precomp_intvals_Hz_hed


!****************************************************************
!  1D EM subroutine precomp_intvals_ved
!    compute integral values for all integrals for VED source
!
!  Rita Streich 2009
!****************************************************************
subroutine precomp_intvals_ved(refl_var,sz,zr)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter
  integer(kind=int32),dimension(nrelmax)  :: iint  !indicates locations for integration results in output structure


  !receivers above source
  if (sz .gt. zr) then

    !compute well-behaved integrals by fast Hankel transform
    !related integrals: iabvC1TMved and iabvC0TMved
    ibesord(1) = 1
    ibesord(2) = 0
    ijrel = 1
    iint(1) = 3    !integral iabvC1TMved
    iint(2) = 2    !integral iabvC0TMved
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iabvC1TMved,iint)


    !integrals for Ex and Ey
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iabvB1TMved,1)

  !iReceiver below source
  elseif (sz .lt. zr) then

    !related integrals: iblwC1TMved and iblwC0TMved
    ibesord(1) = 1
    ibesord(2) = 0
    ijrel = 1
    iint(1) = 3    !integral iblwC1TMved
    iint(2) = 2    !integral iblwC0TMved
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iblwC1TMved,iint)

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iblwB1TMved,1)

  !iReceiver exactly at source depth
  else

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 1._real64
    call precomp_intval_adaptive(refl_var,besorder,iabvB1TMved,1,sz,zr)

    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,iabvC0TMved,2,sz,zr)


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


    !compute well-behaved integrals by fast Hankel transform
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iabvC1TMved,3)

  endif

endsubroutine precomp_intvals_ved


!****************************************************************
!  FD EM subroutine precomp_intvals_Exy_ved
!    compute integral values for the integral needed for Ex and Ey, VED source
!
!  Rita Streich 2009
!****************************************************************
subroutine precomp_intvals_Exy_ved(refl_var,sz,zr)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter


  !receivers above source
  if (sz .gt. zr) then

    !integral for Ex and Ey
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iabvB1TMved,1)

  !iReceiver below source
  elseif (sz .lt. zr) then

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iblwB1TMved,1)

  !iReceiver exactly at source depth
  else

    !logarithmic radii, same as returned from Fast Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 1._real64
    call precomp_intval_adaptive(refl_var,besorder,iabvB1TMved,1,sz,zr)


    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,1),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,1))
    call spline(refl_var%radlog,refl_var%intvalim(:,1),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,1))

  endif

endsubroutine precomp_intvals_Exy_ved


!****************************************************************
!  FD EM subroutine precomp_intvals_Ez_ved
!    compute integral values for the integral needed for Ez, VED source
!
!  Rita Streich 2009
!****************************************************************
subroutine precomp_intvals_Ez_ved(refl_var,sz,zr)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter


  !receivers above source
  if (sz .gt. zr) then

    !integral for Ez
    ibesord = 0
    call precomp_intval_fht(refl_var,ibesord,iabvC0TMved,2)

  !iReceiver below source
  elseif (sz .lt. zr) then

    ibesord = 0
    call precomp_intval_fht(refl_var,ibesord,iblwC0TMved,2)

  !iReceiver exactly at source depth
  else

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,iabvC0TMved,2,sz,zr)


    !do spline interpolation here and not inside precomp_intval_adaptive
    ! so that there is no need to go back and forth between radii and their logarithms

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,2))
    call spline(refl_var%radlog,refl_var%intvalim(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,2))

  endif

endsubroutine precomp_intvals_Ez_ved


!****************************************************************
!  FD EM subroutine precomp_intvals_Hxy_ved
!    compute integral values for all integrals for VED source, Hx and Hy
!
!  Rita Streich 2009-2011
!****************************************************************
subroutine precomp_intvals_Hxy_ved(refl_var,sz,zr)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")


  !receivers above source
  if (sz .ge. zr) then

    !integral for Hx and Hy
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iabvC1TMved,3)

  !iReceiver below source
  else
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iblwC1TMved,3)
  endif

endsubroutine precomp_intvals_Hxy_ved


!****************************************************************
!  FD EM subroutine precomp_intvals_hmd
!    compute integral values for all integrals for HMD source
!
!  Rita Streich 2009
!****************************************************************
subroutine precomp_intvals_hmd(refl_var,sz,zr)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter
  integer(kind=int32),dimension(nrelmax)  :: iint  !indicates locations for integration results in output structure


  !receivers above source
  if (sz .gt. zr) then

    !related integrals: IB1TE, IB0TE, IBz1TE, same integrand just different kappa factors and Bessel function orders
    !basis is IB1TE, contains no kappa, Bessel order 1
    !IB0TE contains kappa^1, Bessel order 0
    !IBz1TE contains kappa^2, Bessel order 1
    ibesord(1:3:2) = 1
    ibesord(2) = 0
    ijrel = 1
    ijrel(1,3) = 2 !factor kappa^2 for third integral
    iint(1) = 3    !integral IabvB1TE
    iint(2) = 1    !integral IabvB0TE
    iint(3) = 10   !integral IabvBz1TE (for Hz)

    call precomp_intval_fht_rel(refl_var,3,ibesord,ijrel,iabvB1TE,iint)

    !related integrals IB1TM and IB0TM
    !Bessel order as before (using first 2 elements only)
    !exponents for first 2 related transforms also as before
    iint(1) = 4    !integral IabvB1TM
    iint(2) = 2    !integral IabvB0TM
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iabvB1TM,iint)

    !related integrals IC1TE and IC0TE
    iint(1) = 8    !integral IabvC1TE
    iint(2) = 6    !integral IabvC0TE
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iabvC1TE,iint)

    !related integrals C_TM
    iint(1) = 9    !integral iabvC1TM
    iint(2) = 7    !integral iabvC0TM
    iint(3) = 5    !integral iabvCz1TM (for Ez)
    call precomp_intval_fht_rel(refl_var,3,ibesord,ijrel,iabvC1TM,iint)

  !iReceiver below source
  elseif (sz .lt. zr) then

    !related integrals: IB1TE, IB0TE, IBz1TE, same integrand just different kappa factors and Bessel function orders
    !basis is IB1TE, contains no kappa, Bessel order 1
    !IB0TE contains kappa^1, Bessel order 0
    !IBz1TE contains kappa^2, Bessel order 1
    ibesord(1:3:2) = 1
    ibesord(2) = 0
    ijrel = 1
    ijrel(1,3) = 2 !factor kappa^2 for third integral
    iint(1) = 3    !integral IblwB1TE
    iint(2) = 1    !integral IblwB0TE
    iint(3) = 10   !integral IblwBz1TE (for Hz)

    call precomp_intval_fht_rel(refl_var,3,ibesord,ijrel,iblwB1TE,iint)

    !related integrals IB1TM and IB0TM
    !Bessel order as before (using first 2 elements only)
    !exponents for first 2 related transforms also as before
    iint(1) = 4    !integral IblwB1TM
    iint(2) = 2    !integral IblwB0TM
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iblwB1TM,iint)

    !related integrals IC1TE and IC0TE
    iint(1) = 8    !integral IblwC1TE
    iint(2) = 6    !integral IblwC0TE
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iblwC1TE,iint)

    !related integrals C_TM
    iint(1) = 9    !integral iblwC1TM
    iint(2) = 7    !integral iblwC0TM
    iint(3) = 5    !integral iblwCz1TM (for Ez)
    call precomp_intval_fht_rel(refl_var,3,ibesord,ijrel,iblwC1TM,iint)

  !iReceiver exactly at source depth
  else

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,iabvB0TE,1,sz,zr)
    call precomp_intval_adaptive(refl_var,besorder,iabvB0TM,2,sz,zr)

    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,iabvC0TE,6,sz,zr)

    besorder = 1._real64
    call precomp_intval_adaptive(refl_var,besorder,iabvBz1TE,10,sz,zr)


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


    !compute well-behaved integrals by fast Hankel transform
    !related integrals C_TM
    ibesord(1:3:2) = 1
    ibesord(2) = 0
    ijrel = 1
    ijrel(1,3) = 2 !factor kappa^2 for third integral
    iint(1) = 9    !integral iabvC1TM
    iint(2) = 7    !integral iabvC0TM
    iint(3) = 5    !integral iabvCz1TM (for Ez)
    call precomp_intval_fht_rel(refl_var,3,ibesord,ijrel,iabvC1TM,iint)

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iabvB1TE,3)
    call precomp_intval_fht(refl_var,ibesord,iabvB1TM,4)

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iabvC1TE,8)

  endif

endsubroutine precomp_intvals_hmd


!****************************************************************
!  FD EM subroutine precomp_intvals_Exy_hmd
!    compute integral values for all integrals needed for Ex and Ey, HMD source
!
!  Rita Streich 2009
!****************************************************************
subroutine precomp_intvals_Exy_hmd(refl_var,sz,zr)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter
  integer(kind=int32),dimension(nrelmax)  :: iint  !indicates locations for integration results in output structure


  !receivers above source
  if (sz .gt. zr) then

    !related integrals: IB1TE, IB0TE, same integrand just different kappa factors and Bessel function orders
    !basis is IB1TE, contains no kappa, Bessel order 1
    !IB0TE contains kappa^1, Bessel order 0
    !IBz1TE contains kappa^2, Bessel order 1
    ibesord(1) = 1
    ibesord(2) = 0
    ijrel = 1
    iint(1) = 3    !integral IabvB1TE
    iint(2) = 1    !integral IabvB0TE

    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iabvB1TE,iint)

    !related integrals IB1TM and IB0TM
    !Bessel order as before (using first 2 elements only)
    !exponents for first 2 related transforms also as before
    iint(1) = 4    !integral IabvB1TM
    iint(2) = 2    !integral IabvB0TM
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iabvB1TM,iint)

  !iReceiver below source
  elseif (sz .lt. zr) then

    !related integrals: IB1TE, IB0TE, same integrand just different kappa factors and Bessel function orders
    !basis is IB1TE, contains no kappa, Bessel order 1
    !IB0TE contains kappa^1, Bessel order 0
    ibesord(1) = 1
    ibesord(2) = 0
    ijrel = 1
    iint(1) = 3    !integral IblwB1TE
    iint(2) = 1    !integral IblwB0TE

    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iblwB1TE,iint)

    !related integrals IB1TM and IB0TM
    !Bessel order as before (using first 2 elements only)
    !exponents for first 2 related transforms also as before
    iint(1) = 4    !integral IblwB1TM
    iint(2) = 2    !integral IblwB0TM
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iblwB1TM,iint)

  !iReceiver exactly at source depth
  else

    !precompute logarithmic radii, as returned from FHT
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,iabvB0TE,1,sz,zr)
    call precomp_intval_adaptive(refl_var,besorder,iabvB0TM,2,sz,zr)


    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,1),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,1))
    call spline(refl_var%radlog,refl_var%intvalim(:,1),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,1))
    call spline(refl_var%radlog,refl_var%intvalre(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,2))
    call spline(refl_var%radlog,refl_var%intvalim(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,2))


    !compute well-behaved integrals by fast Hankel transform
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iabvB1TE,3)
    call precomp_intval_fht(refl_var,ibesord,iabvB1TM,4)

  endif

endsubroutine precomp_intvals_Exy_hmd


!****************************************************************
!  FD EM subroutine precomp_intvals_Ez_hmd
!    compute integral values for the integral needed for Ez, HMD source
!
!  Rita Streich 2009
!****************************************************************
subroutine precomp_intvals_Ez_hmd(refl_var,sz,zr)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")


  !receivers above source
  if (sz .ge. zr) then

    !integral for Ez
    !the integral is "easy" --> don't need special case for iReceiver exactly at source depth
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iabvCz1TM,5)

  !iReceiver below source
  else

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iblwCz1TM,5)

  endif

endsubroutine precomp_intvals_Ez_hmd


!****************************************************************
!  FD EM subroutine precomp_intvals_Hxy_hmd
!    compute integral values for all integrals for HMD source, Hx/Hy
!
!  Rita Streich 2009-2011
!****************************************************************
subroutine precomp_intvals_Hxy_hmd(refl_var,sz,zr)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter
  integer(kind=int32),dimension(nrelmax)  :: iint  !indicates locations for integration results in output structure


  !receivers above source
  if (sz .gt. zr) then

    !related integrals: IB1TE, IB0TE, IBz1TE, same integrand just different kappa factors and Bessel function orders
    !basis is IB1TE, contains no kappa, Bessel order 1
    !IB0TE contains kappa^1, Bessel order 0
    !IBz1TE contains kappa^2, Bessel order 1
    ibesord = 1
    ibesord(2) = 0
    ijrel = 1

    !related integrals IC1TE and IC0TE
    iint(1) = 8    !integral IabvC1TE
    iint(2) = 6    !integral IabvC0TE
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iabvC1TE,iint)

    !related integrals C_TM
    iint(1) = 9    !integral iabvC1TM
    iint(2) = 7    !integral iabvC0TM
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iabvC1TM,iint)

  !iReceiver below source
  elseif (sz .lt. zr) then

    !related integrals: IB1TE, IB0TE, IBz1TE, same integrand just different kappa factors and Bessel function orders
    !basis is IB1TE, contains no kappa, Bessel order 1
    !IB0TE contains kappa^1, Bessel order 0
    !IBz1TE contains kappa^2, Bessel order 1
    ibesord = 1
    ibesord(2) = 0
    ijrel = 1

    !related integrals IC1TE and IC0TE
    iint(1) = 8    !integral IblwC1TE
    iint(2) = 6    !integral IblwC0TE
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iblwC1TE,iint)

    !related integrals C_TM
    iint(1) = 9    !integral iblwC1TM
    iint(2) = 7    !integral iblwC0TM
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iblwC1TM,iint)

  !iReceiver exactly at source depth
  else

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,iabvC0TE,6,sz,zr)

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
    !related integrals C_TM
    ibesord = 1
    ibesord(2) = 0
    ijrel = 1
    iint(1) = 9    !integral iabvC1TM
    iint(2) = 7    !integral iabvC0TM
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iabvC1TM,iint)

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iabvC1TE,8)

  endif

endsubroutine precomp_intvals_Hxy_hmd


!****************************************************************
!  FD EM subroutine precomp_intvals_Hz_hmd
!    compute integral values for all integrals for HMD source, Hz only
!
!  Rita Streich 2009-2011
!****************************************************************
subroutine precomp_intvals_Hz_hmd(refl_var,sz,zr)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter


  !receivers above source
  if (sz .gt. zr) then

    !IBz1TE contains kappa^2, Bessel order 1
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iabvBz1TE,10)

  !iReceiver below source
  elseif (sz .lt. zr) then

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iblwBz1TE,10)

  !iReceiver exactly at source depth
  else

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 1._real64
    call precomp_intval_adaptive(refl_var,besorder,iabvBz1TE,10,sz,zr)

    !do spline interpolation here and not inside precomp_intval_adaptive
    ! so that there is no need to go back and forth between radii and their logarithms

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,10),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,10))
    call spline(refl_var%radlog,refl_var%intvalim(:,10),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,10))

  endif

endsubroutine precomp_intvals_Hz_hmd


!****************************************************************
!  1D EM subroutine precomp_intvals_vmd
!    compute integral values for all integrals for VMD source
!
!  Rita Streich 2009
!****************************************************************
subroutine precomp_intvals_vmd(refl_var,sz,zr)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter
  integer(kind=int32),dimension(nrelmax)  :: iint  !indicates locations for integration results in output structure


  !receivers above source
  if (sz .gt. zr) then

    !related integrals: iabvA1TEvmd and iabvA0TEvmd
    ibesord(1) = 1
    ibesord(2) = 0
    ijrel = 1
    iint(1) = 1    !integral iabvA1TEvmd (for Ex and Ey)
    iint(2) = 3    !integral iabvA0TEvmd (for Hz)
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iabvA1TEvmd,iint)

    !no integral for Ez - Ez for VMD source is zero

    !integral for Hx and Hy
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iabvD1TEvmd,2)

  !iReceiver below source
  elseif (sz .lt. zr) then

    !related integrals: iblwA1TEvmd and iblwA0TEvmd
    ibesord(1) = 1
    ibesord(2) = 0
    ijrel = 1
    iint(1) = 1    !integral iblwA1TEvmd
    iint(2) = 3    !integral iblwA0TEvmd
    call precomp_intval_fht_rel(refl_var,2,ibesord,ijrel,iblwA1TEvmd,iint)

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iblwD1TEvmd,2)

  !iReceiver exactly at source depth
  else

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 1._real64
    call precomp_intval_adaptive(refl_var,besorder,iabvD1TEvmd,2,sz,zr)

    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,iabvA0TEvmd,3,sz,zr)

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
    call precomp_intval_fht(refl_var,ibesord,iabvA1TEvmd,1)

  endif

endsubroutine precomp_intvals_vmd


!****************************************************************
!  FD EM subroutine precomp_intvals_Exy_vmd
!    compute integral values for the integral needed for Ex and Ey, VMD source
!
!  Rita Streich 2009
!****************************************************************
subroutine precomp_intvals_Exy_vmd(refl_var,sz,zr)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")


  !receivers above source
  if (sz .ge. zr) then

    !integral for Ex and Ey
    !the integral is "easy" --> don't need special case for iReceiver exactly at source depth
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iabvA1TEvmd,1)

  !receivers below source
  else
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iblwA1TEvmd,1)
  endif

endsubroutine precomp_intvals_Exy_vmd


!****************************************************************
!  FD EM subroutine precomp_intvals_Hxy_vmd
!    compute integral values for all integrals for VMD source, Hx/Hy
!
!  Rita Streich 2009-2011
!****************************************************************
subroutine precomp_intvals_Hxy_vmd(refl_var,sz,zr)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter


  !receivers above source
  if (sz .gt. zr) then

    !integral for Hx and Hy
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iabvD1TEvmd,2)

  !iReceiver below source
  elseif (sz .lt. zr) then

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iblwD1TEvmd,2)

  !iReceiver exactly at source depth
  else

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 1._real64
    call precomp_intval_adaptive(refl_var,besorder,iabvD1TEvmd,2,sz,zr)

    !do spline interpolation here and not inside precomp_intval_adaptive
    ! so that there is no need to go back and forth between radii and their logarithms

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,2))
    call spline(refl_var%radlog,refl_var%intvalim(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,2))

  endif

endsubroutine precomp_intvals_Hxy_vmd


!****************************************************************
!  FD EM subroutine precomp_intvals_Hz_vmd
!    compute integral values for all integrals for VMD source, Hz only
!
!  Rita Streich 2009-2011
!****************************************************************
subroutine precomp_intvals_Hz_vmd(refl_var,sz,zr)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32),dimension(nrelmax)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter


  !receivers above source
  if (sz .gt. zr) then

    !integral for Hz
    ibesord = 0
    call precomp_intval_fht(refl_var,ibesord,iabvA0TEvmd,3)

  !iReceiver below source
  elseif (sz .lt. zr) then

    ibesord = 0
    call precomp_intval_fht(refl_var,ibesord,iblwA0TEvmd,3)

  !iReceiver exactly at source depth
  else

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,iabvA0TEvmd,3,sz,zr)

    !do spline interpolation here and not inside precomp_intval_adaptive
    ! so that there is no need to go back and forth between radii and their logarithms

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,3),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,3))
    call spline(refl_var%radlog,refl_var%intvalim(:,3),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,3))

  endif

endsubroutine precomp_intvals_Hz_vmd


!****************************************************************
!  FD EM subroutine precomp_intvals_wire
!    compute integral values for all integrals along wires
!
!  Rita Streich 2009
!****************************************************************
subroutine precomp_intvals_wire(refl_var,sz,zr)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter


  !receivers above source
  if (sz .gt. zr) then

    !integrals along wires
    !integral for Ex (and Ey)
    ibesord = 0
    call precomp_intval_fht(refl_var,ibesord,i1abvExwire,1)

    !integral for Hy (and Hx)
    ibesord = 0
    call precomp_intval_fht(refl_var,ibesord,iabvD0TE,2)

    !integral for Hz
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iabvAz1TE,3)

    !integrals for end points
    !integral for Ex and Ey
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,i2abvExwire,4)

    !integral for Ez
    ibesord = 0
    call precomp_intval_fht(refl_var,ibesord,iabvD0TM,5)

    !integral for Hx and Hy
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,idabvHxwire,6)

  !iReceiver below source
  elseif (sz .lt. zr) then

    ibesord = 0
    call precomp_intval_fht(refl_var,ibesord,i1blwExwire,1)
    call precomp_intval_fht(refl_var,ibesord,iblwD0TE,2)
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iblwAz1TE,3)

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,i2blwExwire,4)
    ibesord = 0
    call precomp_intval_fht(refl_var,ibesord,iblwD0TM,5)
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,idblwHxwire,6)

  !iReceiver exactly at source depth
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
    call precomp_intval_adaptive(refl_var,besorder,i1abvExwire,1,sz,zr)
    call precomp_intval_adaptive(refl_var,besorder,iabvD0TE,2,sz,zr)
    !integrals for end points
    besorder = 1._real64
    call precomp_intval_adaptive(refl_var,besorder,i2abvExwire,4,sz,zr)
    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,iabvD0TM,5,sz,zr)


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


    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,iabvAz1TE,3)
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,idabvHxwire,6)

  endif

endsubroutine precomp_intvals_wire


!****************************************************************
!  FD EM subroutine precomp_intvals_Exy_wire
!    compute integral values for Ex or Ey component for integrals along wires
!    and at wire end points
!
!  Rita Streich 2009
!****************************************************************
subroutine precomp_intvals_Exy_wire(refl_var,sz,zr)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter


  !receivers above source
  if (sz .gt. zr) then

    !integral for Ex and Ey along wires
    ibesord = 0
    call precomp_intval_fht(refl_var,ibesord,i1abvExwire,1)

    !integral for Ex and Ey, end point contributions
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,i2abvExwire,4)

  !iReceiver below source
  elseif (sz .lt. zr) then

    ibesord = 0
    call precomp_intval_fht(refl_var,ibesord,i1blwExwire,1)
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,i2blwExwire,4)

  !iReceiver exactly at source depth
  else

    !precompute logarithmic radii as returned from fast Hankel transform
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    !integrals along wire
    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,i1abvExwire,1,sz,zr)
    besorder = 1._real64
    call precomp_intval_adaptive(refl_var,besorder,i2abvExwire,4,sz,zr)

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

  endif

endsubroutine precomp_intvals_Exy_wire


!****************************************************************
! FD EM subroutine precomp_intvals_Ez_wire
!   compute integral values for the integral needed for Ez (contribution from wire end points)
!
! Rita Streich 2009
!****************************************************************
subroutine precomp_intvals_Ez_wire(refl_var,sz,zr)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter


  !receivers above source
  if (sz .gt. zr) then

    !integral for Ez
    ibesord = 0
    call precomp_intval_fht(refl_var,ibesord,iabvD0TM,5)

  !iReceiver below source
  elseif (sz .lt. zr) then

    ibesord = 0
    call precomp_intval_fht(refl_var,ibesord,iblwD0TM,5)

  !iReceiver exactly at source depth
  else

    !precompute logarithmic radii, start from largest, so that it's the same as for Hankel transforms
    fact = exp(logspace)
    refl_var%radlog(refl_var%nrad) = refl_var%rmax
    do irad = refl_var%nrad-1,1,-1
        refl_var%radlog(irad) = refl_var%radlog(irad+1) / fact
    enddo

    !compute badly behaved integrals by adaptive integration
    besorder = 0._real64
    call precomp_intval_adaptive(refl_var,besorder,iabvD0TM,5,sz,zr)


    !get spline derivatives

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,5),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,5))
    call spline(refl_var%radlog,refl_var%intvalim(:,5),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,5))

  endif

endsubroutine precomp_intvals_Ez_wire


!****************************************************************
!  FD EM subroutine precomp_intvals_Hxy_wire
!    compute integral values for all integrals along wires, Hx/Hy
!
!  Rita Streich 2009-2011
!****************************************************************
subroutine precomp_intvals_Hxy_wire(refl_var,sz,zr)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")
  real(kind=real64)             :: besorder  !Bessel function order (real)
  real(kind=real64)             :: fact      !multiplication factor for radii
  integer(kind=int32)           :: irad      !radius counter


  !receivers above source
  if (sz .gt. zr) then

    !integral for Hy and Hx along wire
    ibesord = 0
    call precomp_intval_fht(refl_var,ibesord,iabvD0TE,2)

    !integral for Hx and Hy, end points
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,idabvHxwire,6)

  !iReceiver below source
  elseif (sz .lt. zr) then

    ibesord = 0
    call precomp_intval_fht(refl_var,ibesord,iblwD0TE,2)
    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,idblwHxwire,6)

  !iReceiver exactly at source depth
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
    call precomp_intval_adaptive(refl_var,besorder,iabvD0TE,2,sz,zr)

    !get spline derivatives

    !take logarithm of radii - this makes spline interpolation more accurate
    do irad=1,refl_var%nrad
      refl_var%radlog(irad) = log10(refl_var%radlog(irad))
    enddo

    !spline derivatives
    call spline(refl_var%radlog,refl_var%intvalre(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivre(:,2))
    call spline(refl_var%radlog,refl_var%intvalim(:,2),refl_var%nrad,spl_endval,spl_endval,refl_var%spl_derivim(:,2))

    ibesord = 1
    call precomp_intval_fht(refl_var,ibesord,idabvHxwire,6)
  endif

endsubroutine precomp_intvals_Hxy_wire


!****************************************************************
!  FD EM subroutine precomp_intvals_Hz_wire
!    compute integral values for all integrals along wires, Hz only
!
!  Rita Streich 2009-2011
!****************************************************************
subroutine precomp_intvals_Hz_wire(refl_var,sz,zr)

  implicit none

  !external variables
  type(refl_struct)             :: refl_var  !stuff to remember throughout 1D computations
  real(kind=real64),intent(in)  :: sz,zr     !source and iReceiver depth

  !internal variables
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")


  ibesord = 1

  !receivers above source
  if (sz .ge. zr) then
    !integrals along wires
    !integral for Hz - no need for adaptive integration for this integral
    call precomp_intval_fht(refl_var,ibesord,iabvAz1TE,3)
  !iReceiver below source
  else
    call precomp_intval_fht(refl_var,ibesord,iblwAz1TE,3)
  endif

endsubroutine precomp_intvals_Hz_wire


