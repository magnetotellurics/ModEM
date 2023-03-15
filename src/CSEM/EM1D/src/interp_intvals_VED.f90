!------------------------------------------------------------
!  1D EM subroutine interp_intvals_ved_allcomp
!
!  get interpolated field values at iReceiver locations, VED source, 
!    same coordinates for all field components
!
!  Rita Streich 2009
!------------------------------------------------------------
subroutine interp_intvals_ved_allcomp(refl_var,src,ifreq,sz,zr,bgdat,Ex,Ey,Ez,Hx,Hy,Hz,omeps_srcv,omeps_recv, &
  funcB1TMved,funcC0TMved,funcC1TMved,ilay, funcB1TMfwd,funcC0TMfwd,funcC1TMfwd, &
  funcB1TMvedv,funcC0TMvedv,funcC1TMvedv,Exv,Eyv,Ezv,Hxv,Hyv,Hzv)

  implicit none

  !external variables
  type(refl_struct)               :: refl_var   !everything needed throughout 1D computations
  type(sorec),intent(in)          :: src        !source definition (we need the currents here)
  integer(kind=int32),intent(in)  :: ifreq      !frequency index
  real(kind=real64),intent(in)    :: sz,zr      !source and iReceiver depth
  type(backgrounddata) :: bgdat      !coordinate vectors and output EM fields
  complex(kind=real64),dimension(:),target     :: Ex,Ey,Ez  !electric field: nr of receivers x 3 components
  complex(kind=real64),dimension(:),target     :: Hx,Hy,Hz  !magnetic field: nr of receivers x 3 components
  complex(kind=real64),intent(in) :: omeps_srcv  !omega * epsilon in source layer
  complex(kind=real64),intent(in) :: omeps_recv  !omega * epsilon in iReceiver layer
  complex(kind=real64),external   :: funcB1TMved,funcC0TMved,funcC1TMved
  integer(kind=int32),intent(in)  :: ilay       !layer index for derivatives - leave at zero for forward modeling
  complex(kind=real64),external,optional   :: funcB1TMfwd,funcC0TMfwd,funcC1TMfwd
  complex(kind=real64),external,optional   :: funcB1TMvedv,funcC0TMvedv,funcC1TMvedv
  complex(kind=real64),dimension(:),target,optional   :: Exv,Eyv,Ezv  !electric field: nr of receivers x 3 components for epsv
  complex(kind=real64),dimension(:),target,optional   :: Hxv,Hyv,Hzv  !magnetic field: nr of receivers x 3 components for epsv

  !internal variables
  integer(kind=int32)   :: isrc    !source element counter
  real(kind=real64)     :: x,y,r   !temp source-iReceiver distances
  integer(kind=int32)   :: irec    !iReceiver counter
  real(kind=real64)     :: beta    !temp angle

  complex(kind=real64)          :: IB1TMved,IC0TMved,IC1TMved !interpolated integral values

  logical,dimension(nintVEDdvti)    :: wellbehaved  !indicates if Hankel integration can be used
  logical                       :: sz_eq_zr     !indicates if source and reveicer are at the same depth
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")

  complex(kind=real64)  :: Er,Hbeta             !temp field values in cylindrical coordinates
  real(kind=real64)     :: cosbeta,sinbeta      !cos(beta) and sin(beta), precompute for efficiency
  complex(kind=real64)  :: cur      !temp source current
  complex(kind=real64)  :: JEz      !source current times constants
  integer(kind=int32)   :: idx      !source element index
  integer(kind=int32)   :: recidx   !iReceiver index
  complex(kind=real64)  :: fact_ErHb,fact_Ezsrc,fact_Ezrec,fact_Ez !factors in front of integrals
  real(kind=real64)     :: isrec    !indicates if source and iReceiver are in the same layer
  !flag for computing epsv derivatives, not needed for forward computation, so "aniso" value cannot be used here
  logical               :: with_dvert
  complex(kind=real64),dimension(:),pointer :: Exrec,Eyrec,Ezrec,Hxrec,Hyrec,Hzrec !point to E and H for isotropic and Ev, Hv for VTI case

  !indicators for fast Hankel transform or adaptive integration
  wellbehaved = .true.
  sz_eq_zr = .false.
  if (sz.eq.zr) then
    wellbehaved(1:2) = .false.
    sz_eq_zr = .true.
  endif

  if (present(funcC1TMvedv)) then
    with_dvert = .true.
    if (sz_eq_zr) then
      wellbehaved(7:8) = .false.
    endif
  else
    with_dvert = .false.
  endif


  do isrc = refl_var%isrcstart,refl_var%isrcend

    !get source current times constant factors
    idx = refl_var%isrcperz(isrc)
    cur = conjg(src%cur(idx,ifreq))
    JEz = cur * src%ljz(idx) / dfourpi

    do irec=refl_var%irecstart,refl_var%irecend
      recidx = refl_var%irecperzExy(irec)

      y = bgdat%Exypos(recidx,2) - refl_var%ys(isrc)
      x = bgdat%Exypos(recidx,1) - refl_var%xs(isrc)

      r = sqrt(x**2 + y**2)

      beta = atan2(y,x)
      cosbeta = cos(beta)
      sinbeta = sin(beta)

      r_is_zero: if (r.eq.0._real64) then

        !quick & dirty: skip the point if iReceiver is right at source point
        if (sz_eq_zr) then
          if (refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle iReceiver right at source point yet!'

          !still add special contribution right at source point - here for forward modeling only!!!
          if (ilay .eq. 0)  Ez(recidx) = Ez(recidx) + JEz / (dci * omeps_srcv)

          !skip numerical integral evaluations - they would most likely fail anyway, need to look at this...
          cycle
        endif

        !we only have an Ez component for r=0
        IB1TMved = 0._real64
        IC0TMved = compute_1valr0(funcC0TMved)
        IC1TMved = 0._real64

      else  !r is not zero
        smallr: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .true.
          call prepare_refcoef(refl_var,r,ved,aniso)

          ibesord = 1
          IB1TMved = compute_1val(funcB1TMved,r,sz,zr,ibesord,wellbehaved(1))
          ibesord = 0
          IC0TMved = compute_1val(funcC0TMved,r,sz,zr,ibesord,wellbehaved(2))
          ibesord = 1
          IC1TMved = compute_1val(funcC1TMved,r,sz,zr,ibesord,wellbehaved(3))

        else !r larger than threshold radius for spline interpolation

          IB1TMved = splinterp_1val(refl_var,1,r)
          IC0TMved = splinterp_1val(refl_var,2,r)
          IC1TMved = splinterp_1val(refl_var,3,r)

        endif smallr

      endif r_is_zero
 
      !Er relative to horizontal source orientation (factor 4pi and dipole length absorbed into JEh)
      Er = JEz * dci * IB1TMved / omeps_srcv

      !Ex, Ey rotated back to original x-y coordinate system
      !don't take complex conjugate here, but at the very end!
      Ex(recidx) = Ex(recidx) + cosbeta*Er
      Ey(recidx) = Ey(recidx) + sinbeta*Er

      !Ez in original x-y coordinate system
      Ez(recidx) = Ez(recidx) - JEz * IC0TMved / (omeps_srcv * omeps_recv)

      Hbeta = JEz * dci * IC1TMved / omeps_srcv

      !Hx, Hy in original x-y coordinate system
      Hx(recidx) = Hx(recidx) - sinbeta*Hbeta
      Hy(recidx) = Hy(recidx) + cosbeta*Hbeta

      !no Hz for VED source
      
      dvert: if (with_dvert) then
        r_is_zerov: if (r.eq.0._real64) then

          !cycling already done for epsh
          !special contribution right at source point is ignored here and added later!!!
          !if (sz_eq_zr) cycle

          !we only have an Ez component for r=0
          IB1TMved = 0._real64
          IC0TMved = compute_1valr0(funcC0TMvedv)
          IC1TMved = 0._real64

        else  !r is not zero
          smallrv: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,ved,aniso)

            ibesord = 1
            IB1TMved = compute_1val(funcB1TMvedv,r,sz,zr,ibesord,wellbehaved(7))
            ibesord = 0
            IC0TMved = compute_1val(funcC0TMvedv,r,sz,zr,ibesord,wellbehaved(8))
            ibesord = 1
            IC1TMved = compute_1val(funcC1TMvedv,r,sz,zr,ibesord,wellbehaved(9))

          else !r larger than threshold radius for spline interpolation

            IB1TMved = splinterp_1val(refl_var,7,r)
            IC0TMved = splinterp_1val(refl_var,8,r)
            IC1TMved = splinterp_1val(refl_var,9,r)
          endif smallrv

        endif r_is_zerov

        !Er relative to horizontal source orientation (factor 4pi and dipole length absorbed into JEh)
        Er = JEz * dci * IB1TMved / omeps_srcv

        !Ex, Ey rotated back to original x-y coordinate system
        !don't take complex conjugate here, but at the very end!
        Exv(recidx) = Exv(recidx) + cosbeta*Er
        Eyv(recidx) = Eyv(recidx) + sinbeta*Er

        !Ez in original x-y coordinate system
        Ezv(recidx) = Ezv(recidx) - JEz * IC0TMved / (omeps_srcv * omeps_recv)

        Hbeta = JEz * dci * IC1TMved / omeps_srcv

        !Hx, Hy in original x-y coordinate system
        Hxv(recidx) = Hxv(recidx) - sinbeta*Hbeta
        Hyv(recidx) = Hyv(recidx) + cosbeta*Hbeta
      endif dvert

    enddo !irec
  enddo  !source elements


  !special terms for all components in source layer:
  !add to overall integral derivatives for isotropic case, to epsv derivatives for VTI
  deriv_ilaysrc: if (ilay .eq. ilaysrc) then

    if (with_dvert) then
      Exrec => Exv
      Eyrec => Eyv
      Ezrec => Ezv
      Hxrec => Hxv
      Hyrec => Hyv
      Hzrec => Hzv
    else
      Exrec => Ex
      Eyrec => Ey
      Ezrec => Ez
      Hxrec => Hx
      Hyrec => Hy
      Hzrec => Hz
    endif

    !check if receivers are in the same layer - in that case, add all special terms in one go
    if (ilay .eq. ilayrec) then
      isrec = 1._real64
    else
      isrec = 0._real64
    endif

    !indicators for fast Hankel transform or adaptive integration
    wellbehaved(4:6) = wellbehaved(1:3)

    do isrc = refl_var%isrcstart,refl_var%isrcend

      !get source current times constant factors
      idx = refl_var%isrcperz(isrc)
      cur = conjg(src%cur(idx,ifreq))
      JEz = cur * src%ljz(idx) / dfourpi
      fact_ErHb = - JEz * dci / (omeps_srcv*epsv(ilaysrc))
      fact_Ezsrc = JEz / (omeps_srcv * omeps_recv * epsv(ilaysrc))
      fact_Ezrec = isrec * JEz / (omeps_srcv * omeps_recv * epsv(ilayrec))
      fact_Ez = fact_Ezsrc + fact_Ezrec

      do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzExy(irec)

        y = bgdat%Exypos(recidx,2) - refl_var%ys(isrc)
        x = bgdat%Exypos(recidx,1) - refl_var%xs(isrc)

        r = sqrt(x**2 + y**2)

        beta = atan2(y,x)
        cosbeta = cos(beta)
        sinbeta = sin(beta)

        r_is_zero2: if (r.eq.0._real64) then

          !quick & dirty: skip the point if iReceiver is right at source point
          if (sz_eq_zr) then
            !still add special contribution right at source point - here for derivatives only
            Ezrec(recidx) = Ezrec(recidx) - JEz / (dci * omeps_srcv * epsv(ilaysrc))
            !skip numerical integral evaluations - they would most likely fail anyway, need to look at this...
            cycle
          endif

          !we only have an Ez component for r=0
          IB1TMved = 0._real64
          IC0TMved = compute_1valr0(funcC0TMfwd)
          IC1TMved = 0._real64

        else  !r is not zero
          smallr2: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            call prepare_refcoef(refl_var,r,ved,aniso)

            ibesord = 1
            IB1TMved = compute_1val(funcB1TMfwd,r,sz,zr,ibesord,wellbehaved(4))
            ibesord = 0
            IC0TMved = compute_1val(funcC0TMfwd,r,sz,zr,ibesord,wellbehaved(5))
            ibesord = 1
            IC1TMved = compute_1val(funcC1TMfwd,r,sz,zr,ibesord,wellbehaved(6))

          else !r larger than threshold radius for spline interpolation

            IB1TMved = splinterp_1val(refl_var,4,r)
            IC0TMved = splinterp_1val(refl_var,5,r)
            IC1TMved = splinterp_1val(refl_var,6,r)

          endif smallr2
        endif r_is_zero2


        Er = fact_ErHb * IB1TMved
        !Ex, Ey rotated back to original x-y coordinate system
        !don't take complex conjugate here, but at the very end!
        Exrec(recidx) = Exrec(recidx) + cosbeta*Er
        Eyrec(recidx) = Eyrec(recidx) + sinbeta*Er

        !Ez in original x-y coordinate system
        !if receciver layer = source layer, fact_Ez includes contributions for both source and iReceiver layer
        Ezrec(recidx) = Ezrec(recidx) + fact_Ez * IC0TMved

        Hbeta = fact_ErHb * IC1TMved
        !Hx, Hy in original x-y coordinate system
        Hxrec(recidx) = Hxrec(recidx) - sinbeta*Hbeta
        Hyrec(recidx) = Hyrec(recidx) + cosbeta*Hbeta
      enddo !irec
    enddo  !source elements

  !derivatives for iReceiver layer, receivers NOT in source layer
  elseif (ilay .eq. ilayrec) then

    if (with_dvert) then
      Ezrec => Ezv
    else
      Ezrec => Ez
    endif

    !indicators for fast Hankel transform or adaptive integration
    wellbehaved(5) = wellbehaved(2)

    do isrc = refl_var%isrcstart,refl_var%isrcend

      !get source current times constant factors
      idx = refl_var%isrcperz(isrc)
      cur = conjg(src%cur(idx,ifreq))
      JEz = cur * src%ljz(idx) / dfourpi
      fact_Ez = JEz / (omeps_srcv * omeps_recv * epsv(ilayrec))

      do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzExy(irec)

        y = bgdat%Exypos(recidx,2) - refl_var%ys(isrc)
        x = bgdat%Exypos(recidx,1) - refl_var%xs(isrc)
        r = sqrt(x**2 + y**2)

        r_is_zero3: if (r.eq.0._real64) then
          !quick & dirty: skip the point if iReceiver is right at source point
          if (sz_eq_zr) cycle
          IC0TMved = compute_1valr0(funcC0TMfwd)

        else  !r is not zero
          smallr3: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation
            !reflection coeff. for this radius
            call prepare_refcoef(refl_var,r,ved,aniso)

            ibesord = 0
            IC0TMved = compute_1val(funcC0TMfwd,r,sz,zr,ibesord,wellbehaved(5))
          else !r larger than threshold radius for spline interpolation

            IC0TMved = splinterp_1val(refl_var,5,r)
          endif smallr3
        endif r_is_zero3

        Ezrec(recidx) = Ezrec(recidx) + fact_Ez * IC0TMved
      enddo !irec
    enddo  !source elements

  endif deriv_ilaysrc

endsubroutine interp_intvals_ved_allcomp


!------------------------------------------------------------
!  1D EM subroutine interp_intvals_ved_Exy
!
!  get interpolated field values at iReceiver locations, VED source, 
!    Ex and / or Ey only
!
!  Rita Streich 2009
!------------------------------------------------------------
subroutine interp_intvals_ved_Exy(refl_var,src,ifreq,sz,zr,bgdat,Ex,Ey,omeps_srcv, &
  funcB1TMved,ilay, funcB1TMfwd, funcB1TMvedv,Exv,Eyv)

  implicit none

  !external variables
  type(refl_struct)               :: refl_var   !everything needed throughout 1D computations
  type(sorec),intent(in)          :: src        !source definition (we need the currents here)
  integer(kind=int32),intent(in)  :: ifreq      !frequency index
  real(kind=real64),intent(in)    :: sz,zr      !source and iReceiver depth
  type(backgrounddata) :: bgdat      !coordinate vectors and output EM fields
  complex(kind=real64),dimension(:),target     :: Ex,Ey  !electric field: nr of receivers x 3 components
  complex(kind=real64),intent(in) :: omeps_srcv  !omega * epsilon in source layer
  complex(kind=real64),external   :: funcB1TMved
  integer(kind=int32),intent(in)  :: ilay       !layer index for derivatives - leave at zero for forward modeling
  complex(kind=real64),external,optional   :: funcB1TMfwd
  complex(kind=real64),external,optional   :: funcB1TMvedv
  complex(kind=real64),dimension(:),target,optional   :: Exv,Eyv  !electric field: nr of receivers x 3 components for epsv

  !internal variables
  integer(kind=int32)   :: isrc    !source element counter
  real(kind=real64)     :: x,y,r   !temp source-iReceiver distances
  integer(kind=int32)   :: irec    !iReceiver counter
  real(kind=real64)     :: beta    !temp angle

  complex(kind=real64)          :: IB1TMved !interpolated integral values

  logical,dimension(nintVEDdvti)    :: wellbehaved  !indicates if Hankel integration can be used
  logical                       :: sz_eq_zr     !indicates if source and reveicer are at the same depth
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")

  complex(kind=real64)  :: Er                   !temp field values in cylindrical coordinates
  real(kind=real64)     :: cosbeta,sinbeta      !cos(beta) and sin(beta), precompute for efficiency
  complex(kind=real64)  :: cur      !temp source current
  complex(kind=real64)  :: JEz      !source current times constants
  integer(kind=int32)   :: idx      !source element index
  integer(kind=int32)   :: recidx   !iReceiver index
  complex(kind=real64)  :: fact_ErHb  !factors in front of integrals
  real(kind=real64)     :: isrec    !indicates if source and iReceiver are in the same layer
  !flag for computing epsv derivatives, not needed for forward computation, so "aniso" value cannot be used here
  logical               :: with_dvert
  complex(kind=real64),dimension(:),pointer :: Exrec,Eyrec !point to Ex,Ey for isotropic and Exv,Eyv for VTI case

  !indicators for fast Hankel transform or adaptive integration
  wellbehaved = .true.
  sz_eq_zr = .false.
  if (sz.eq.zr) then
    wellbehaved(1:2) = .false.
    sz_eq_zr = .true.
  endif

  if (present(funcB1TMvedv)) then
    with_dvert = .true.
    if (sz_eq_zr) then
      wellbehaved(7:8) = .false.
    endif
  else
    with_dvert = .false.
  endif


  do isrc = refl_var%isrcstart,refl_var%isrcend

    !get source current times constant factors
    idx = refl_var%isrcperz(isrc)
    cur = conjg(src%cur(idx,ifreq))
    JEz = cur * src%ljz(idx) / dfourpi


    !same positions for Ex and Ey
    exy_equalpos: if (bgdat%nExy.gt.0) then

    do irec=refl_var%irecstart,refl_var%irecend
      recidx = refl_var%irecperzExy(irec)

      y = bgdat%Exypos(recidx,2) - refl_var%ys(isrc)
      x = bgdat%Exypos(recidx,1) - refl_var%xs(isrc)

      r = sqrt(x**2 + y**2)

      beta = atan2(y,x)
      cosbeta = cos(beta)
      sinbeta = sin(beta)

      r_is_zero: if (r.eq.0._real64) then

        !quick & dirty: skip the point if iReceiver is right at source point
        if (sz_eq_zr) then
          if (refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle iReceiver right at source point yet!'
          !skip numerical integral evaluations - they would most likely fail anyway, need to look at this...
          cycle
        endif
        IB1TMved = 0._real64
      else  !r is not zero
        smallr: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .true.
          call prepare_refcoef(refl_var,r,ved,aniso)
          ibesord = 1
          IB1TMved = compute_1val(funcB1TMved,r,sz,zr,ibesord,wellbehaved(1))

        else !r larger than threshold radius for spline interpolation

          IB1TMved = splinterp_1val(refl_var,1,r)
        endif smallr
      endif r_is_zero
 
      !Er relative to horizontal source orientation (factor 4pi and dipole length absorbed into JEh)
      Er = JEz * dci * IB1TMved / omeps_srcv

      !Ex, Ey rotated back to original x-y coordinate system
      !don't take complex conjugate here, but at the very end!
      Ex(recidx) = Ex(recidx) + cosbeta*Er
      Ey(recidx) = Ey(recidx) + sinbeta*Er

      dvert: if (with_dvert) then
        r_is_zerov: if (r.eq.0._real64) then

          !cycling already done for epsh
          !special contribution right at source point is ignored here and added later!!!
          !if (sz_eq_zr) cycle
          IB1TMved = 0._real64
        else  !r is not zero
          smallrv: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,ved,aniso)
            ibesord = 1
            IB1TMved = compute_1val(funcB1TMvedv,r,sz,zr,ibesord,wellbehaved(7))

          else !r larger than threshold radius for spline interpolation
            IB1TMved = splinterp_1val(refl_var,7,r)
          endif smallrv
        endif r_is_zerov

        !Er relative to horizontal source orientation (factor 4pi and dipole length absorbed into JEh)
        Er = JEz * dci * IB1TMved / omeps_srcv

        !Ex, Ey rotated back to original x-y coordinate system
        !don't take complex conjugate here, but at the very end!
        Exv(recidx) = Exv(recidx) + cosbeta*Er
        Eyv(recidx) = Eyv(recidx) + sinbeta*Er
      endif dvert
    enddo !irec

    else

      have_ex: if (bgdat%nEx .gt. 0) then

    do irec=refl_var%irecstart,refl_var%irecend
      recidx = refl_var%irecperzExy(irec)

      y = bgdat%Expos(recidx,2) - refl_var%ys(isrc)
      x = bgdat%Expos(recidx,1) - refl_var%xs(isrc)

      r = sqrt(x**2 + y**2)

      beta = atan2(y,x)
      cosbeta = cos(beta)
      sinbeta = sin(beta)

      r_is_zeroex: if (r.eq.0._real64) then

        !quick & dirty: skip the point if iReceiver is right at source point
        if (sz_eq_zr) then
          if (refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle iReceiver right at source point yet!'
          !skip numerical integral evaluations - they would most likely fail anyway, need to look at this...
          cycle
        endif
        IB1TMved = 0._real64
      else  !r is not zero
        smallrex: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .true.
          call prepare_refcoef(refl_var,r,ved,aniso)
          ibesord = 1
          IB1TMved = compute_1val(funcB1TMved,r,sz,zr,ibesord,wellbehaved(1))

        else !r larger than threshold radius for spline interpolation

          IB1TMved = splinterp_1val(refl_var,1,r)
        endif smallrex
      endif r_is_zeroex
 
      !Er relative to horizontal source orientation (factor 4pi and dipole length absorbed into JEh)
      Er = JEz * dci * IB1TMved / omeps_srcv

      !Ex, Ey rotated back to original x-y coordinate system
      !don't take complex conjugate here, but at the very end!
      Ex(recidx) = Ex(recidx) + cosbeta*Er

      dvertex: if (with_dvert) then
        r_is_zerovex: if (r.eq.0._real64) then

          !cycling already done for epsh
          !special contribution right at source point is ignored here and added later!!!
          !if (sz_eq_zr) cycle
          IB1TMved = 0._real64
        else  !r is not zero
          smallrvex: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,ved,aniso)
            ibesord = 1
            IB1TMved = compute_1val(funcB1TMvedv,r,sz,zr,ibesord,wellbehaved(7))

          else !r larger than threshold radius for spline interpolation
            IB1TMved = splinterp_1val(refl_var,7,r)
          endif smallrvex
        endif r_is_zerovex

        !Er relative to horizontal source orientation (factor 4pi and dipole length absorbed into JEh)
        Er = JEz * dci * IB1TMved / omeps_srcv

        !Ex, Ey rotated back to original x-y coordinate system
        !don't take complex conjugate here, but at the very end!
        Exv(recidx) = Exv(recidx) + cosbeta*Er
      endif dvertex
    enddo !irec

      endif have_ex

      have_ey: if (bgdat%nEy .gt. 0) then

    do irec=refl_var%irecstart,refl_var%irecend
      recidx = refl_var%irecperzExy(irec)

      y = bgdat%Eypos(recidx,2) - refl_var%ys(isrc)
      x = bgdat%Eypos(recidx,1) - refl_var%xs(isrc)

      r = sqrt(x**2 + y**2)

      beta = atan2(y,x)
      cosbeta = cos(beta)
      sinbeta = sin(beta)

      r_is_zeroey: if (r.eq.0._real64) then

        !quick & dirty: skip the point if iReceiver is right at source point
        if (sz_eq_zr) then
          if (refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle iReceiver right at source point yet!'
          !skip numerical integral evaluations - they would most likely fail anyway, need to look at this...
          cycle
        endif
        IB1TMved = 0._real64
      else  !r is not zero
        smallrey: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .true.
          call prepare_refcoef(refl_var,r,ved,aniso)
          ibesord = 1
          IB1TMved = compute_1val(funcB1TMved,r,sz,zr,ibesord,wellbehaved(1))

        else !r larger than threshold radius for spline interpolation

          IB1TMved = splinterp_1val(refl_var,1,r)
        endif smallrey
      endif r_is_zeroey
 
      !Er relative to horizontal source orientation (factor 4pi and dipole length absorbed into JEh)
      Er = JEz * dci * IB1TMved / omeps_srcv

      !Ex, Ey rotated back to original x-y coordinate system
      !don't take complex conjugate here, but at the very end!
      Ey(recidx) = Ey(recidx) + sinbeta*Er

      dvertey: if (with_dvert) then
        r_is_zerovey: if (r.eq.0._real64) then

          !cycling already done for epsh
          !special contribution right at source point is ignored here and added later!!!
          !if (sz_eq_zr) cycle
          IB1TMved = 0._real64
        else  !r is not zero
          smallrvey: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,ved,aniso)
            ibesord = 1
            IB1TMved = compute_1val(funcB1TMvedv,r,sz,zr,ibesord,wellbehaved(7))

          else !r larger than threshold radius for spline interpolation
            IB1TMved = splinterp_1val(refl_var,7,r)
          endif smallrvey
        endif r_is_zerovey

        !Er relative to horizontal source orientation (factor 4pi and dipole length absorbed into JEh)
        Er = JEz * dci * IB1TMved / omeps_srcv

        !Ex, Ey rotated back to original x-y coordinate system
        !don't take complex conjugate here, but at the very end!
        Eyv(recidx) = Eyv(recidx) + sinbeta*Er
      endif dvertey
    enddo !irec

      endif have_ey
    endif exy_equalpos

  enddo  !source elements


  !special terms for all components in source layer:
  !add to overall integral derivatives for isotropic case, to epsv derivatives for VTI
  deriv_ilaysrc: if (ilay .eq. ilaysrc) then

    if (with_dvert) then
      Exrec => Exv
      Eyrec => Eyv
    else
      Exrec => Ex
      Eyrec => Ey
    endif

    !check if receivers are in the same layer - in that case, add all special terms in one go
    if (ilay .eq. ilayrec) then
      isrec = 1._real64
    else
      isrec = 0._real64
    endif

    !indicators for fast Hankel transform or adaptive integration
    wellbehaved(4:6) = wellbehaved(1:3)

    do isrc = refl_var%isrcstart,refl_var%isrcend

      !get source current times constant factors
      idx = refl_var%isrcperz(isrc)
      cur = conjg(src%cur(idx,ifreq))
      JEz = cur * src%ljz(idx) / dfourpi
      fact_ErHb = - JEz * dci / (omeps_srcv*epsv(ilaysrc))


    !same positions for Ex and Ey
    exy_equalpossrc: if (bgdat%nExy.gt.0) then

      do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzExy(irec)

        y = bgdat%Exypos(recidx,2) - refl_var%ys(isrc)
        x = bgdat%Exypos(recidx,1) - refl_var%xs(isrc)

        r = sqrt(x**2 + y**2)

        beta = atan2(y,x)
        cosbeta = cos(beta)
        sinbeta = sin(beta)

        r_is_zero2src: if (r.eq.0._real64) then

          !quick & dirty: skip the point if iReceiver is right at source point
          if (sz_eq_zr) cycle

          !we only have an Ez component for r=0
          IB1TMved = 0._real64
        else  !r is not zero
          smallr2src: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            call prepare_refcoef(refl_var,r,ved,aniso)

            ibesord = 1
            IB1TMved = compute_1val(funcB1TMfwd,r,sz,zr,ibesord,wellbehaved(4))
          else !r larger than threshold radius for spline interpolation
            IB1TMved = splinterp_1val(refl_var,4,r)
          endif smallr2src
        endif r_is_zero2src

        Er = fact_ErHb * IB1TMved
        !Ex, Ey rotated back to original x-y coordinate system
        !don't take complex conjugate here, but at the very end!
        Exrec(recidx) = Exrec(recidx) + cosbeta*Er
        Eyrec(recidx) = Eyrec(recidx) + sinbeta*Er
      enddo !irec

    else

      have_exsrc: if (bgdat%nEx .gt. 0) then

      do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzExy(irec)

        y = bgdat%Expos(recidx,2) - refl_var%ys(isrc)
        x = bgdat%Expos(recidx,1) - refl_var%xs(isrc)

        r = sqrt(x**2 + y**2)

        beta = atan2(y,x)
        cosbeta = cos(beta)
        sinbeta = sin(beta)

        r_is_zero2srcex: if (r.eq.0._real64) then

          !quick & dirty: skip the point if iReceiver is right at source point
          if (sz_eq_zr) cycle

          !we only have an Ez component for r=0
          IB1TMved = 0._real64
        else  !r is not zero
          smallr2srcex: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            call prepare_refcoef(refl_var,r,ved,aniso)

            ibesord = 1
            IB1TMved = compute_1val(funcB1TMfwd,r,sz,zr,ibesord,wellbehaved(4))
          else !r larger than threshold radius for spline interpolation
            IB1TMved = splinterp_1val(refl_var,4,r)
          endif smallr2srcex
        endif r_is_zero2srcex

        Er = fact_ErHb * IB1TMved
        Exrec(recidx) = Exrec(recidx) + cosbeta*Er
      enddo !irec

      endif have_exsrc

      have_eysrc: if (bgdat%nEy .gt. 0) then

      do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzExy(irec)

        y = bgdat%Eypos(recidx,2) - refl_var%ys(isrc)
        x = bgdat%Eypos(recidx,1) - refl_var%xs(isrc)

        r = sqrt(x**2 + y**2)

        beta = atan2(y,x)
        cosbeta = cos(beta)
        sinbeta = sin(beta)

        r_is_zero2srcey: if (r.eq.0._real64) then

          !quick & dirty: skip the point if iReceiver is right at source point
          if (sz_eq_zr) cycle

          !we only have an Ez component for r=0
          IB1TMved = 0._real64
        else  !r is not zero
          smallr2srcey: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            call prepare_refcoef(refl_var,r,ved,aniso)

            ibesord = 1
            IB1TMved = compute_1val(funcB1TMfwd,r,sz,zr,ibesord,wellbehaved(4))
          else !r larger than threshold radius for spline interpolation
            IB1TMved = splinterp_1val(refl_var,4,r)
          endif smallr2srcey
        endif r_is_zero2srcey

        Er = fact_ErHb * IB1TMved
        !Ex, Ey rotated back to original x-y coordinate system
        !don't take complex conjugate here, but at the very end!
        Eyrec(recidx) = Eyrec(recidx) + sinbeta*Er
      enddo !irec

      endif have_eysrc
    endif exy_equalpossrc

    enddo  !source elements
  endif deriv_ilaysrc

endsubroutine interp_intvals_ved_Exy


!------------------------------------------------------------
!  1D EM subroutine interp_intvals_ved_Ez
!
!  get interpolated field values at iReceiver locations, VED source, 
!    Ez only
!
!  Rita Streich 2009
!------------------------------------------------------------
subroutine interp_intvals_ved_Ez(refl_var,src,ifreq,sz,zr,bgdat,Ez,omeps_srcv,omeps_recv, &
  funcC0TMved,ilay, funcC0TMfwd, funcC0TMvedv,Ezv)

  implicit none

  !external variables
  type(refl_struct)               :: refl_var   !everything needed throughout 1D computations
  type(sorec),intent(in)          :: src        !source definition (we need the currents here)
  integer(kind=int32),intent(in)  :: ifreq      !frequency index
  real(kind=real64),intent(in)    :: sz,zr      !source and iReceiver depth
  type(backgrounddata) :: bgdat      !coordinate vectors and output EM fields
  complex(kind=real64),dimension(:),target     :: Ez  !electric field: nr of receivers x 3 components
  complex(kind=real64),intent(in) :: omeps_srcv  !omega * epsilon in source layer
  complex(kind=real64),intent(in) :: omeps_recv  !omega * epsilon in iReceiver layer
  complex(kind=real64),external   :: funcC0TMved
  integer(kind=int32),intent(in)  :: ilay       !layer index for derivatives - leave at zero for forward modeling
  complex(kind=real64),external,optional   :: funcC0TMfwd
  complex(kind=real64),external,optional   :: funcC0TMvedv
  complex(kind=real64),dimension(:),target,optional   :: Ezv  !electric field: nr of receivers x 3 components for epsv

  !internal variables
  integer(kind=int32)   :: isrc    !source element counter
  real(kind=real64)     :: x,y,r   !temp source-iReceiver distances
  integer(kind=int32)   :: irec    !iReceiver counter

  complex(kind=real64)          :: IC0TMved !interpolated integral values

  logical,dimension(nintVEDdvti)    :: wellbehaved  !indicates if Hankel integration can be used
  logical                       :: sz_eq_zr     !indicates if source and reveicer are at the same depth
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")

  complex(kind=real64)  :: cur      !temp source current
  complex(kind=real64)  :: JEz      !source current times constants
  integer(kind=int32)   :: idx      !source element index
  integer(kind=int32)   :: recidx   !iReceiver index
  complex(kind=real64)  :: fact_Ezsrc,fact_Ezrec,fact_Ez !factors in front of integrals
  real(kind=real64)     :: isrec    !indicates if source and iReceiver are in the same layer
  !flag for computing epsv derivatives, not needed for forward computation, so "aniso" value cannot be used here
  logical               :: with_dvert
  complex(kind=real64),dimension(:),pointer :: Ezrec !point to E and H for isotropic and Ev, Hv for VTI case

  !indicators for fast Hankel transform or adaptive integration
  wellbehaved = .true.
  sz_eq_zr = .false.
  if (sz.eq.zr) then
    wellbehaved(1:2) = .false.
    sz_eq_zr = .true.
  endif

  if (present(funcC0TMvedv)) then
    with_dvert = .true.
    if (sz_eq_zr) then
      wellbehaved(7:8) = .false.
    endif
  else
    with_dvert = .false.
  endif


  do isrc = refl_var%isrcstart,refl_var%isrcend

    !get source current times constant factors
    idx = refl_var%isrcperz(isrc)
    cur = conjg(src%cur(idx,ifreq))
    JEz = cur * src%ljz(idx) / dfourpi

    do irec=refl_var%irecstart,refl_var%irecend
      recidx = refl_var%irecperzEz(irec)

      y = bgdat%Ezpos(recidx,2) - refl_var%ys(isrc)
      x = bgdat%Ezpos(recidx,1) - refl_var%xs(isrc)
      r = sqrt(x**2 + y**2)

      r_is_zero: if (r.eq.0._real64) then

        !quick & dirty: skip the point if iReceiver is right at source point
        if (sz_eq_zr) then
          if (refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle iReceiver right at source point yet!'
          !still add special contribution right at source point - here for forward modeling only!!!
          if (ilay .eq. 0)  Ez(recidx) = Ez(recidx) + JEz / (dci * omeps_srcv)
          !skip numerical integral evaluations - they would most likely fail anyway, need to look at this...
          cycle
        endif
        IC0TMved = compute_1valr0(funcC0TMved)
      else  !r is not zero
        smallr: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .true.
          call prepare_refcoef(refl_var,r,ved,aniso)

          ibesord = 0
          IC0TMved = compute_1val(funcC0TMved,r,sz,zr,ibesord,wellbehaved(2))
        else !r larger than threshold radius for spline interpolation
          IC0TMved = splinterp_1val(refl_var,2,r)
        endif smallr
      endif r_is_zero
 
      !Ez in original x-y coordinate system
      Ez(recidx) = Ez(recidx) - JEz * IC0TMved / (omeps_srcv * omeps_recv)

      dvert: if (with_dvert) then
        r_is_zerov: if (r.eq.0._real64) then

          !cycling already done for epsh
          !special contribution right at source point is ignored here and added later!!!
          !if (sz_eq_zr) cycle

          IC0TMved = compute_1valr0(funcC0TMvedv)
        else  !r is not zero
          smallrv: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,ved,aniso)
            ibesord = 0
            IC0TMved = compute_1val(funcC0TMvedv,r,sz,zr,ibesord,wellbehaved(8))
          else !r larger than threshold radius for spline interpolation
            IC0TMved = splinterp_1val(refl_var,8,r)
          endif smallrv
        endif r_is_zerov

        !Ez in original x-y coordinate system
        Ezv(recidx) = Ezv(recidx) - JEz * IC0TMved / (omeps_srcv * omeps_recv)
      endif dvert
    enddo !irec
  enddo  !source elements


  !special terms for all components in source layer:
  !add to overall integral derivatives for isotropic case, to epsv derivatives for VTI
  deriv_ilaysrc: if (ilay .eq. ilaysrc) then

    if (with_dvert) then
      Ezrec => Ezv
    else
      Ezrec => Ez
    endif

    !check if receivers are in the same layer - in that case, add all special terms in one go
    if (ilay .eq. ilayrec) then
      isrec = 1._real64
    else
      isrec = 0._real64
    endif

    !indicators for fast Hankel transform or adaptive integration
    wellbehaved(4:6) = wellbehaved(1:3)

    do isrc = refl_var%isrcstart,refl_var%isrcend

      !get source current times constant factors
      idx = refl_var%isrcperz(isrc)
      cur = conjg(src%cur(idx,ifreq))
      JEz = cur * src%ljz(idx) / dfourpi
      fact_Ezsrc = JEz / (omeps_srcv * omeps_recv * epsv(ilaysrc))
      fact_Ezrec = isrec * JEz / (omeps_srcv * omeps_recv * epsv(ilayrec))
      fact_Ez = fact_Ezsrc + fact_Ezrec

      do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzEz(irec)

        y = bgdat%Ezpos(recidx,2) - refl_var%ys(isrc)
        x = bgdat%Ezpos(recidx,1) - refl_var%xs(isrc)
        r = sqrt(x**2 + y**2)

        r_is_zero2: if (r.eq.0._real64) then

          !quick & dirty: skip the point if iReceiver is right at source point
          if (sz_eq_zr) then
            !still add special contribution right at source point - here for derivatives only
            Ezrec(recidx) = Ezrec(recidx) - JEz / (dci * omeps_srcv * epsv(ilaysrc))
            cycle
          endif
          IC0TMved = compute_1valr0(funcC0TMfwd)
        else  !r is not zero
          smallr2: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            call prepare_refcoef(refl_var,r,ved,aniso)

            ibesord = 0
            IC0TMved = compute_1val(funcC0TMfwd,r,sz,zr,ibesord,wellbehaved(5))
          else !r larger than threshold radius for spline interpolation
            IC0TMved = splinterp_1val(refl_var,5,r)
          endif smallr2
        endif r_is_zero2

        Ezrec(recidx) = Ezrec(recidx) + fact_Ez * IC0TMved
      enddo !irec
    enddo  !source elements

  !derivatives for iReceiver layer, receivers NOT in source layer
  elseif (ilay .eq. ilayrec) then

    if (with_dvert) then
      Ezrec => Ezv
    else
      Ezrec => Ez
    endif

    !indicators for fast Hankel transform or adaptive integration
    wellbehaved(5) = wellbehaved(2)

    do isrc = refl_var%isrcstart,refl_var%isrcend

      !get source current times constant factors
      idx = refl_var%isrcperz(isrc)
      cur = conjg(src%cur(idx,ifreq))
      JEz = cur * src%ljz(idx) / dfourpi
      fact_Ez = JEz / (omeps_srcv * omeps_recv * epsv(ilayrec))

      do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzEz(irec)

        y = bgdat%Ezpos(recidx,2) - refl_var%ys(isrc)
        x = bgdat%Ezpos(recidx,1) - refl_var%xs(isrc)
        r = sqrt(x**2 + y**2)

        r_is_zero3: if (r.eq.0._real64) then
          !quick & dirty: skip the point if iReceiver is right at source point
          if (sz_eq_zr) cycle
          IC0TMved = compute_1valr0(funcC0TMfwd)

        else  !r is not zero
          smallr3: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation
            !reflection coeff. for this radius
            call prepare_refcoef(refl_var,r,ved,aniso)

            ibesord = 0
            IC0TMved = compute_1val(funcC0TMfwd,r,sz,zr,ibesord,wellbehaved(5))
          else !r larger than threshold radius for spline interpolation

            IC0TMved = splinterp_1val(refl_var,5,r)
          endif smallr3
        endif r_is_zero3

        Ezrec(recidx) = Ezrec(recidx) + fact_Ez * IC0TMved
      enddo !irec
    enddo  !source elements

  endif deriv_ilaysrc

endsubroutine interp_intvals_ved_Ez


!------------------------------------------------------------
!  1D EM subroutine interp_intvals_ved_Hxy
!
!  get interpolated field values at iReceiver locations, VED source, 
!    Hx and/or Hy only
!
!  Rita Streich 2009
!------------------------------------------------------------
subroutine interp_intvals_ved_Hxy(refl_var,src,ifreq,sz,zr,bgdat,Hx,Hy,omeps_srcv, &
  funcC1TMved,ilay, funcC1TMfwd, funcC1TMvedv,Hxv,Hyv)

  implicit none

  !external variables
  type(refl_struct)               :: refl_var   !everything needed throughout 1D computations
  type(sorec),intent(in)          :: src        !source definition (we need the currents here)
  integer(kind=int32),intent(in)  :: ifreq      !frequency index
  real(kind=real64),intent(in)    :: sz,zr      !source and iReceiver depth
  type(backgrounddata) :: bgdat      !coordinate vectors and output EM fields
  complex(kind=real64),dimension(:),target     :: Hx,Hy  !magnetic field: nr of receivers x 3 components
  complex(kind=real64),intent(in) :: omeps_srcv  !omega * epsilon in source layer
  complex(kind=real64),external   :: funcC1TMved
  integer(kind=int32),intent(in)  :: ilay       !layer index for derivatives - leave at zero for forward modeling
  complex(kind=real64),external,optional   :: funcC1TMfwd
  complex(kind=real64),external,optional   :: funcC1TMvedv
  complex(kind=real64),dimension(:),target,optional   :: Hxv,Hyv  !magnetic field: nr of receivers x 3 components for epsv

  !internal variables
  integer(kind=int32)   :: isrc    !source element counter
  real(kind=real64)     :: x,y,r   !temp source-iReceiver distances
  integer(kind=int32)   :: irec    !iReceiver counter
  real(kind=real64)     :: beta    !temp angle

  complex(kind=real64)          :: IC1TMved !interpolated integral values

  logical,dimension(nintVEDdvti)    :: wellbehaved  !indicates if Hankel integration can be used
  logical                       :: sz_eq_zr     !indicates if source and reveicer are at the same depth
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")

  complex(kind=real64)  :: Hbeta                !temp field values in cylindrical coordinates
  real(kind=real64)     :: cosbeta,sinbeta      !cos(beta) and sin(beta), precompute for efficiency
  complex(kind=real64)  :: cur      !temp source current
  complex(kind=real64)  :: JEz      !source current times constants
  integer(kind=int32)   :: idx      !source element index
  integer(kind=int32)   :: recidx   !iReceiver index
  complex(kind=real64)  :: fact_ErHb !factors in front of integrals
  real(kind=real64)     :: isrec    !indicates if source and iReceiver are in the same layer
  !flag for computing epsv derivatives, not needed for forward computation, so "aniso" value cannot be used here
  logical               :: with_dvert
  complex(kind=real64),dimension(:),pointer :: Hxrec,Hyrec !point to E and H for isotropic and Ev, Hv for VTI case

  !indicators for fast Hankel transform or adaptive integration
  wellbehaved = .true.
  sz_eq_zr = .false.
  if (sz.eq.zr) then
    wellbehaved(1:2) = .false.
    sz_eq_zr = .true.
  endif

  if (present(funcC1TMvedv)) then
    with_dvert = .true.
    if (sz_eq_zr) then
      wellbehaved(7:8) = .false.
    endif
  else
    with_dvert = .false.
  endif


  do isrc = refl_var%isrcstart,refl_var%isrcend

    !get source current times constant factors
    idx = refl_var%isrcperz(isrc)
    cur = conjg(src%cur(idx,ifreq))
    JEz = cur * src%ljz(idx) / dfourpi


    !same positions for Hx and Hy
    hxy_equalpos: if (bgdat%nHxy.gt.0) then

    do irec=refl_var%irecstart,refl_var%irecend
      recidx = refl_var%irecperzHxy(irec)

      y = bgdat%Hxypos(recidx,2) - refl_var%ys(isrc)
      x = bgdat%Hxypos(recidx,1) - refl_var%xs(isrc)
      r = sqrt(x**2 + y**2)

      beta = atan2(y,x)
      cosbeta = cos(beta)
      sinbeta = sin(beta)

      r_is_zero: if (r.eq.0._real64) then

        !quick & dirty: skip the point if iReceiver is right at source point
        if (sz_eq_zr) then
          if (refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle iReceiver right at source point yet!'
          cycle
        endif
        IC1TMved = 0._real64
      else  !r is not zero
        smallr: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .true.
          call prepare_refcoef(refl_var,r,ved,aniso)

          ibesord = 1
          IC1TMved = compute_1val(funcC1TMved,r,sz,zr,ibesord,wellbehaved(3))
        else !r larger than threshold radius for spline interpolation
          IC1TMved = splinterp_1val(refl_var,3,r)
        endif smallr
      endif r_is_zero
 
      Hbeta = JEz * dci * IC1TMved / omeps_srcv

      !Hx, Hy in original x-y coordinate system
      Hx(recidx) = Hx(recidx) - sinbeta*Hbeta
      Hy(recidx) = Hy(recidx) + cosbeta*Hbeta

      dvert: if (with_dvert) then
        r_is_zerov: if (r.eq.0._real64) then
          !cycling already done for epsh
          !if (sz_eq_zr) cycle
          IC1TMved = 0._real64
        else  !r is not zero
          smallrv: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,ved,aniso)

            ibesord = 1
            IC1TMved = compute_1val(funcC1TMvedv,r,sz,zr,ibesord,wellbehaved(9))
          else !r larger than threshold radius for spline interpolation
            IC1TMved = splinterp_1val(refl_var,9,r)
          endif smallrv
        endif r_is_zerov

        Hbeta = JEz * dci * IC1TMved / omeps_srcv

        !Hx, Hy in original x-y coordinate system
        Hxv(recidx) = Hxv(recidx) - sinbeta*Hbeta
        Hyv(recidx) = Hyv(recidx) + cosbeta*Hbeta
      endif dvert
    enddo !irec

    else

      have_hx: if (bgdat%nHx .gt. 0) then

    do irec=refl_var%irecstart,refl_var%irecend
      recidx = refl_var%irecperzHxy(irec)

      y = bgdat%Hxpos(recidx,2) - refl_var%ys(isrc)
      x = bgdat%Hxpos(recidx,1) - refl_var%xs(isrc)
      r = sqrt(x**2 + y**2)

      beta = atan2(y,x)
      cosbeta = cos(beta)
      sinbeta = sin(beta)

      r_is_zerohx: if (r.eq.0._real64) then

        !quick & dirty: skip the point if iReceiver is right at source point
        if (sz_eq_zr) then
          if (refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle iReceiver right at source point yet!'
          cycle
        endif
        IC1TMved = 0._real64
      else  !r is not zero
        smallrhx: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .true.
          call prepare_refcoef(refl_var,r,ved,aniso)

          ibesord = 1
          IC1TMved = compute_1val(funcC1TMved,r,sz,zr,ibesord,wellbehaved(3))
        else !r larger than threshold radius for spline interpolation
          IC1TMved = splinterp_1val(refl_var,3,r)
        endif smallrhx
      endif r_is_zerohx
 
      Hbeta = JEz * dci * IC1TMved / omeps_srcv

      !Hx, Hy in original x-y coordinate system
      Hx(recidx) = Hx(recidx) - sinbeta*Hbeta

      dverthx: if (with_dvert) then
        r_is_zerovhx: if (r.eq.0._real64) then
          !cycling already done for epsh
          !if (sz_eq_zr) cycle
          IC1TMved = 0._real64
        else  !r is not zero
          smallrvhx: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,ved,aniso)

            ibesord = 1
            IC1TMved = compute_1val(funcC1TMvedv,r,sz,zr,ibesord,wellbehaved(9))
          else !r larger than threshold radius for spline interpolation
            IC1TMved = splinterp_1val(refl_var,9,r)
          endif smallrvhx
        endif r_is_zerovhx

        Hbeta = JEz * dci * IC1TMved / omeps_srcv

        !Hx, Hy in original x-y coordinate system
        Hxv(recidx) = Hxv(recidx) - sinbeta*Hbeta
      endif dverthx
    enddo !irec

      endif have_hx

      have_hy: if (bgdat%nHy .gt. 0) then

    do irec=refl_var%irecstart,refl_var%irecend
      recidx = refl_var%irecperzHxy(irec)

      y = bgdat%Hypos(recidx,2) - refl_var%ys(isrc)
      x = bgdat%Hypos(recidx,1) - refl_var%xs(isrc)
      r = sqrt(x**2 + y**2)

      beta = atan2(y,x)
      cosbeta = cos(beta)
      sinbeta = sin(beta)

      r_is_zerohy: if (r.eq.0._real64) then

        !quick & dirty: skip the point if iReceiver is right at source point
        if (sz_eq_zr) then
          if (refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle iReceiver right at source point yet!'
          cycle
        endif
        IC1TMved = 0._real64
      else  !r is not zero
        smallrhy: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .true.
          call prepare_refcoef(refl_var,r,ved,aniso)

          ibesord = 1
          IC1TMved = compute_1val(funcC1TMved,r,sz,zr,ibesord,wellbehaved(3))
        else !r larger than threshold radius for spline interpolation
          IC1TMved = splinterp_1val(refl_var,3,r)
        endif smallrhy
      endif r_is_zerohy
 
      Hbeta = JEz * dci * IC1TMved / omeps_srcv

      !Hx, Hy in original x-y coordinate system
      Hy(recidx) = Hy(recidx) + cosbeta*Hbeta

      dverthy: if (with_dvert) then
        r_is_zerovhy: if (r.eq.0._real64) then
          !cycling already done for epsh
          !if (sz_eq_zr) cycle
          IC1TMved = 0._real64
        else  !r is not zero
          smallrvhy: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,ved,aniso)

            ibesord = 1
            IC1TMved = compute_1val(funcC1TMvedv,r,sz,zr,ibesord,wellbehaved(9))
          else !r larger than threshold radius for spline interpolation
            IC1TMved = splinterp_1val(refl_var,9,r)
          endif smallrvhy
        endif r_is_zerovhy

        Hbeta = JEz * dci * IC1TMved / omeps_srcv

        !Hx, Hy in original x-y coordinate system
        Hyv(recidx) = Hyv(recidx) + cosbeta*Hbeta
      endif dverthy
    enddo !irec

      endif have_hy
    endif hxy_equalpos


  enddo  !source elements


  !special terms for all components in source layer:
  !add to overall integral derivatives for isotropic case, to epsv derivatives for VTI
  deriv_ilaysrc: if (ilay .eq. ilaysrc) then

    if (with_dvert) then
      Hxrec => Hxv
      Hyrec => Hyv
    else
      Hxrec => Hx
      Hyrec => Hy
    endif

    !check if receivers are in the same layer - in that case, add all special terms in one go
    if (ilay .eq. ilayrec) then
      isrec = 1._real64
    else
      isrec = 0._real64
    endif

    !indicators for fast Hankel transform or adaptive integration
    wellbehaved(4:6) = wellbehaved(1:3)

    do isrc = refl_var%isrcstart,refl_var%isrcend

      !get source current times constant factors
      idx = refl_var%isrcperz(isrc)
      cur = conjg(src%cur(idx,ifreq))
      JEz = cur * src%ljz(idx) / dfourpi
      fact_ErHb = - JEz * dci / (omeps_srcv*epsv(ilaysrc))


    !same positions for Hx and Hy
    hxy_equalpossrc: if (bgdat%nHxy.gt.0) then

      do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzHxy(irec)

        y = bgdat%Hxypos(recidx,2) - refl_var%ys(isrc)
        x = bgdat%Hxypos(recidx,1) - refl_var%xs(isrc)
        r = sqrt(x**2 + y**2)

        beta = atan2(y,x)
        cosbeta = cos(beta)
        sinbeta = sin(beta)

        r_is_zero2: if (r.eq.0._real64) then

          !quick & dirty: skip the point if iReceiver is right at source point
          if (sz_eq_zr) cycle

          IC1TMved = 0._real64
        else  !r is not zero
          smallr2: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            call prepare_refcoef(refl_var,r,ved,aniso)

            ibesord = 1
            IC1TMved = compute_1val(funcC1TMfwd,r,sz,zr,ibesord,wellbehaved(6))
          else !r larger than threshold radius for spline interpolation
            IC1TMved = splinterp_1val(refl_var,6,r)
          endif smallr2
        endif r_is_zero2

        Hbeta = fact_ErHb * IC1TMved
        !Hx, Hy in original x-y coordinate system
        Hxrec(recidx) = Hxrec(recidx) - sinbeta*Hbeta
        Hyrec(recidx) = Hyrec(recidx) + cosbeta*Hbeta
      enddo !irec

    else

      have_hxsrc: if (bgdat%nHx .gt. 0) then

      do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzHxy(irec)

        y = bgdat%Hxpos(recidx,2) - refl_var%ys(isrc)
        x = bgdat%Hxpos(recidx,1) - refl_var%xs(isrc)
        r = sqrt(x**2 + y**2)

        beta = atan2(y,x)
        sinbeta = sin(beta)

        r_is_zero2hx: if (r.eq.0._real64) then

          !quick & dirty: skip the point if iReceiver is right at source point
          if (sz_eq_zr) cycle

          IC1TMved = 0._real64
        else  !r is not zero
          smallr2hx: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            call prepare_refcoef(refl_var,r,ved,aniso)

            ibesord = 1
            IC1TMved = compute_1val(funcC1TMfwd,r,sz,zr,ibesord,wellbehaved(6))
          else !r larger than threshold radius for spline interpolation
            IC1TMved = splinterp_1val(refl_var,6,r)
          endif smallr2hx
        endif r_is_zero2hx

        Hbeta = fact_ErHb * IC1TMved
        !Hx, Hy in original x-y coordinate system
        Hxrec(recidx) = Hxrec(recidx) - sinbeta*Hbeta
      enddo !irec

      endif have_hxsrc

      have_hysrc: if (bgdat%nHy .gt. 0) then

      do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzHxy(irec)

        y = bgdat%Hypos(recidx,2) - refl_var%ys(isrc)
        x = bgdat%Hypos(recidx,1) - refl_var%xs(isrc)
        r = sqrt(x**2 + y**2)

        beta = atan2(y,x)
        cosbeta = cos(beta)

        r_is_zero2hy: if (r.eq.0._real64) then

          !quick & dirty: skip the point if iReceiver is right at source point
          if (sz_eq_zr) cycle

          IC1TMved = 0._real64
        else  !r is not zero
          smallr2hy: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            call prepare_refcoef(refl_var,r,ved,aniso)

            ibesord = 1
            IC1TMved = compute_1val(funcC1TMfwd,r,sz,zr,ibesord,wellbehaved(6))
          else !r larger than threshold radius for spline interpolation
            IC1TMved = splinterp_1val(refl_var,6,r)
          endif smallr2hy
        endif r_is_zero2hy

        Hbeta = fact_ErHb * IC1TMved
        !Hx, Hy in original x-y coordinate system
        Hyrec(recidx) = Hyrec(recidx) + cosbeta*Hbeta
      enddo !irec

      endif have_hysrc
    endif hxy_equalpossrc

    enddo  !source elements
  endif deriv_ilaysrc

endsubroutine interp_intvals_ved_Hxy









