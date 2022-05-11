!!$!****************************************************************
!!$!  1D EM function interp_1val - inefficient and OBSOLETE!
!!$!    interpolate one integral value for one point
!!$!    assumes that values at logarithmically spaced radii and spline derivatives are available
!!$!
!!$!  Rita Streich 2009
!!$!**************************************************************
!!$function interp_1val(refl_var,intfunc,iint,r,sz,zr,ibesord,wellbehaved,rspl) result(val)
!!$
!!$  implicit none
!!$
!!$  !EXTERNAL VARIABLES
!!$  complex(kind=real64)   :: val         !interpolated integral value
!!$  type(refl_struct)      :: refl_var    !structure for 1D computations
!!$  complex(kind=real64),external   :: intfunc          !function to use for small (or zero) radius
!!$  integer(kind=int32),intent(in)  :: iint             !which part of "intvalues" to use
!!$  real(kind=real64),intent(in)    :: r                !radius
!!$  real(kind=real64),intent(in)    :: sz,zr            !source and receiver depth
!!$  integer(kind=int32),dimension(NREL),intent(in)  :: ibesord  !bessel function order
!!$  logical,intent(in)     :: wellbehaved  !indicates if we can use fast Hankel transform for sz = zr
!!$  real(kind=real64),intent(in)    :: rspl !minimum radius at which spline interpolation is allowed
!!$
!!$  !internal variables
!!$  integer(kind=int32)    :: ierr  !error index
!!$  integer(kind=int32)    :: npieces  !number of subintervals in adaptive integratoin (for zero radius)
!!$  real(kind=real64)      :: besorder !Bessel function order
!!$  integer(kind=int32)    :: newint   !indicates that integral kernels in adaptive integration need to be computed from scratch
!!$  real(kind=real64)      :: intvalre,intvalim    !real and imag part of interpolated integral value
!!$  complex(kind=real64),dimension(1,1)  :: valtmp   !temp integral value
!!$  real(kind=real64),dimension(1)     :: rout     !temp radius returned by FHT
!!$  integer(kind=int32)    :: NOFUN1   !Nr of function evaluations in FHT
!!$  real(kind=real64)      :: errtot   !estimates absolute rrrro in numerical integration
!!$  real(kind=real64)      :: rlog     !log10 of radius
!!$
!!$
!!$  !radius smaller than threshold for spline interpolation
!!$  if (r .lt. rspl) then
!!$
!!$    !zero radius, but source and receiver not at the same depth
!!$    if (r .eq. 0.d0) then
!!$
!!$      !Bessel function J1(0) = 0, so no need to enter adaptive integration if besorder =1
!!$      if (ibesord(1) .ne. 1) then
!!$
!!$        !the case that receiver is exactly at source point is NOT included - work this out later!!!
!!$
!!$        !get real part of integral
!!$        !realval will call intfunc
!!$        call qagic (realval, intfunc, 0._real64, 1, abserr, relerr, intvalre, errtot, NOFUN1, ierr )
!!$        !get imag part of integral
!!$        call qagic (imagval, intfunc, 0._real64, 1, abserr, relerr, intvalim, errtot, NOFUN1, ierr )
!!$
!!$        val = cmplx(intvalre,intvalim)
!!$
!!$      else
!!$        val = 0._real64
!!$      endif
!!$
!!$    !non-zero, small radius
!!$    else
!!$
!!$      !receiver not at source depth or "easy" integral
!!$      if (wellbehaved) then
!!$
!!$        !use temp array for integral values because it has to be a 2D array
!!$        call zhankl(r,1,NREL,TOL,NTOL,ibesord,intfunc,IJREL,ZWORK,valtmp,rout,NOFUN1,ierr)
!!$        val = valtmp(1,1)
!!$
!!$      !receiver exactly at source depth and badly behaved integral
!!$      else
!!$
!!$        npieces = getnpieces(r,sz,zr)
!!$        besorder = real(ibesord(1),kind=real64)
!!$        newint = 1
!!$        ikap = -1 !indicate that we have to recompute all interface refl./transm. coeff.
!!$        CALL BESAUT(intvalre,intvalim,besorder,GAUSLO,GAUSHI,r,intfunc,relerr,abserr,npieces,newint,0,branchpt,ierr)
!!$        val = cmplx(intvalre,intvalim)
!!$
!!$      endif !receiver depth relative to source depth
!!$
!!$    endif
!!$
!!$  !radius sufficiently large so we can apply spline interpolation
!!$  else
!!$    rlog = log10(r)
!!$    call splint(refl_var%radlog,refl_var%intvalre(:,iint),refl_var%spl_derivre(:,iint),refl_var%nrad,rlog,intvalre)
!!$    call splint(refl_var%radlog,refl_var%intvalim(:,iint),refl_var%spl_derivim(:,iint),refl_var%nrad,rlog,intvalim)
!!$    val = cmplx(intvalre,intvalim)
!!$
!!$  endif
!!$
!!$endfunction interp_1val


!****************************************************************
!  1D EM function splinterp_1val
!    spline interpolate one integral value for one point
!    assumes that values at logarithmically spaced radii and spline derivatives are available
!
!  Rita Streich 2009
!**************************************************************
function splinterp_1val(refl_var,iint,r) result(val)

  implicit none

  !EXTERNAL VARIABLES
  complex(kind=real64)   :: val         !interpolated integral value
  type(refl_struct)      :: refl_var    !structure for 1D computations
  integer(kind=int32),intent(in)  :: iint             !which part of "intvalues" to use
  real(kind=real64),intent(in)    :: r                !radius

  !internal variables
  real(kind=real64)      :: rlog     !log10 of radius
  real(kind=real64)      :: intvalre,intvalim    !real and imag part of interpolated integral value


  rlog = log10(r)
  call splint(refl_var%radlog,refl_var%intvalre(:,iint),refl_var%spl_derivre(:,iint),refl_var%nrad,rlog,intvalre)
  call splint(refl_var%radlog,refl_var%intvalim(:,iint),refl_var%spl_derivim(:,iint),refl_var%nrad,rlog,intvalim)
  val = cmplx(intvalre,intvalim)

endfunction splinterp_1val


!****************************************************************
!  1D EM function compute_1val
!    explicitly compute one integral value for one point near source
!
!  Rita Streich 2009
!**************************************************************
function compute_1val(intfunc,r,sz,zr,ibesord,wellbehaved) result(val)

  implicit none

  !EXTERNAL VARIABLES
  complex(kind=real64)   :: val         !interpolated integral value
  complex(kind=real64),external   :: intfunc          !function to use for small (or zero) radius
  real(kind=real64),intent(in)    :: r                !radius
  real(kind=real64),intent(in)    :: sz,zr            !source and receiver depth
  integer(kind=int32),dimension(NREL),intent(in)  :: ibesord  !bessel function order
  logical,intent(in)     :: wellbehaved  !indicates if we can use fast Hankel transform for sz = zr

  !internal variables
  integer(kind=int32)    :: ierr  !error index
  integer(kind=int32)    :: npieces  !number of subintervals in adaptive integratoin (for zero radius)
  real(kind=real64)      :: besorder !Bessel function order
  integer(kind=int32)    :: newint   !indicates that integral kernels in adaptive integration need to be computed from scratch
  real(kind=real64)      :: intvalre,intvalim    !real and imag part of interpolated integral value
  complex(kind=real64),dimension(1,1)  :: valtmp   !temp integral value
  real(kind=real64),dimension(1)     :: rout     !temp radius returned by FHT
  integer(kind=int32)    :: NOFUN1   !Nr of function evaluations in FHT


  !receiver not at source depth or "easy" integral
  if (wellbehaved) then

    !use temp array for integral values because it has to be a 2D array
    call zhankl(r,1,NREL,TOL,NTOL,ibesord,intfunc,IJREL,ZWORK,valtmp,rout,NOFUN1,ierr)
    val = valtmp(1,1)

  !receiver exactly at source depth and badly behaved integral
  else

    npieces = getnpieces(r,sz,zr)
    besorder = real(ibesord(1),kind=real64)
    newint = 1
    ikap = -1 !indicate that we have to recompute all interface refl./transm. coeff.
    CALL BESAUT(intvalre,intvalim,besorder,GAUSLO,GAUSHI,r,intfunc,relerr,abserr,npieces,newint,0,branchpt,ierr)
    val = cmplx(intvalre,intvalim)

  endif !receiver depth relative to source depth

endfunction compute_1val


!****************************************************************
!  1D EM function compute_1valr0
!    explicitly compute one integral value for one point at horizontal position of the source point
!
!  Rita Streich 2009
!**************************************************************
function compute_1valr0(intfunc) result(val)

  implicit none

  !EXTERNAL VARIABLES
  complex(kind=real64)   :: val         !interpolated integral value
  complex(kind=real64),external   :: intfunc          !function to use for small (or zero) radius
!!$  integer(kind=int32),dimension(NREL),intent(in)  :: ibesord  !bessel function order

  !internal variables
  integer(kind=int32)    :: ierr  !error index
  real(kind=real64)      :: intvalre,intvalim    !real and imag part of interpolated integral value
  integer(kind=int32)    :: NOFUN1   !Nr of function evaluations in FHT
  real(kind=real64)      :: errtot   !estimates absolute rrrro in numerical integration


  !Bessel function J1(0) = 0, so no need to enter adaptive integration if besorder =1
!!$  if (ibesord(1) .ne. 1) then

    !the case that receiver is exactly at source point is NOT included - work this out later!!!

    ikap = -1 !indicate that we have to recompute all interface refl./transm. coeff.
    !get real part of integral
    !realval will call intfunc
    call qagic (realval, intfunc, 0._real64, 1, abserr, relerr, intvalre, errtot, NOFUN1, ierr )
    !get imag part of integral
    call qagic (imagval, intfunc, 0._real64, 1, abserr, relerr, intvalim, errtot, NOFUN1, ierr )

    val = cmplx(intvalre,intvalim,kind=real64)

!!$  else
!!$    val = 0._real64
!!$  endif

endfunction compute_1valr0


