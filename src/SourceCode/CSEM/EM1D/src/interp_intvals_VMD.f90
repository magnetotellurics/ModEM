!------------------------------------------------------------
!  1D EM subroutine interp_intvals_vmd_allcomp
!
!  get interpolated field values at receiver locations, VMD source
!  NO change required for VTI since VMD fields do not depend on epsv
!  Rita Streich 2009-2011
!------------------------------------------------------------
subroutine interp_intvals_vmd_allcomp(refl_var,src,ifreq,sz,zr,bgdat,Ex,Ey,Ez,Hx,Hy,Hz,j_om_mu,ommusq, &
  funcA1TEvmd,funcD1TEvmd,funcA0TEvmd,ilay)

  implicit none

  !external variables
  type(refl_struct)               :: refl_var   !everything needed throughout 1D computations
  type(sorec),intent(in)          :: src        !source definition (we need the currents here)
  integer(kind=int32),intent(in)  :: ifreq      !frequency index
  real(kind=real64),intent(in)    :: sz,zr      !source and receiver depth
  type(backgrounddata) :: bgdat      !coordinate vectors and output EM fields
  complex(kind=real64),dimension(:)   :: Ex,Ey,Ez,Hx,Hy,Hz  !electric and magnetic field vectors
  complex(kind=real64),intent(in) :: j_om_mu    !j * omega * mu0
  real(kind=real64),intent(in)    :: ommusq     ! (omega * mu0) **2
  complex(kind=real64),external   :: funcA1TEvmd,funcD1TEvmd,funcA0TEvmd
  integer(kind=int32),intent(in)  :: ilay       !layer index for derivatives - leave at zero for forward modeling

  !internal variables
  integer(kind=int32)   :: isrc    !source element counter
  real(kind=real64)     :: x,y,r   !temp source-receiver distances
  integer(kind=int32)   :: irec    !receiver counter
  real(kind=real64)     :: beta    !temp angle

  complex(kind=real64)          :: IA1TEvmd,ID1TEvmd,IA0TEvmd !interpolated integral values

  logical,dimension(nintVMDdvti)    :: wellbehaved  !indicates if Hankel integration can be used
  logical                       :: sz_eq_zr     !indicates if source and reveicer are at the same depth
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")

  complex(kind=real64)  :: Hr,Ebeta             !temp field values in cylindrical coordinates
  real(kind=real64)     :: cosbeta,sinbeta      !cos(beta) and sin(beta), precompute for efficiency
  complex(kind=real64)  :: cur      !temp source current
  complex(kind=real64)  :: JMz      !source current times constants
  integer(kind=int32)   :: idx      !source element index
  integer(kind=int32)   :: recidx   !receiver index

  !indicators for fast Hankel transform or adaptive integration
  wellbehaved = .TRUE.
  sz_eq_zr = .false.
  if(sz.eq.zr) then
    wellbehaved(1:2) = .false.
    sz_eq_zr = .TRUE.
  endif


  do isrc = refl_var%isrcstart,refl_var%isrcend

    !get source current times constant factors
    idx = refl_var%isrcperz(isrc)
    cur = conjg(src%cur(idx,ifreq))
    JMz = - cur * src%akz(idx) * j_om_mu / dfourpi

    do irec=refl_var%irecstart,refl_var%irecend
      recidx = refl_var%irecperzExy(irec)

      y = bgdat%Exypos(recidx,2) - refl_var%ys(isrc)
      x = bgdat%Exypos(recidx,1) - refl_var%xs(isrc)
      r = sqrt(x**2 + y**2)

      beta = atan2(y,x)
      cosbeta = cos(beta)
      sinbeta = sin(beta)

      !check case r=0 only once...
      r_is_zero: if(r .eq. 0._real64) then

        !quick & dirty: skip the point if receiver is right at source point
        if(sz_eq_zr) then
          if(refl_var%infolevel.ge.output_more) then
            write(*,'(a)') 'WARNING: cannot handle receiver right at source point yet!'
          endif
          !still add special contribution right at source point
          if(ilay.eq.0) Hz(recidx) = Hz(recidx) - JMz / j_om_mu
          !skip numerical integral evaluations - they would most likely fail anyway, need to look at this...
          cycle
        endif

        IA1TEvmd = 0._real64
        ID1TEvmd = 0._real64
        IA0TEvmd = compute_1valr0(funcA0TEvmd)
      else !r is not zero
        smallr: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .TRUE.
          call prepare_refcoef(refl_var,r,vmd,aniso)

          ibesord = 1
          IA1TEvmd = compute_1val(funcA1TEvmd,r,sz,zr,ibesord,wellbehaved(1))
          ibesord = 1
          ID1TEvmd = compute_1val(funcD1TEvmd,r,sz,zr,ibesord,wellbehaved(2))
          ibesord = 0
          IA0TEvmd = compute_1val(funcA0TEvmd,r,sz,zr,ibesord,wellbehaved(3))

        else !r larger than threshold radius for spline interpolation

          IA1TEvmd = splinterp_1val(refl_var,1,r)
          ID1TEvmd = splinterp_1val(refl_var,2,r)
          IA0TEvmd = splinterp_1val(refl_var,3,r)
        endif smallr
      endif r_is_zero

      Hr = JMz / j_om_mu * ID1TEvmd

      !Hx, Hy in original x-y coordinate system
      Hx(recidx) = Hx(recidx) + cosbeta*Hr
      Hy(recidx) = Hy(recidx) + sinbeta*Hr

      Ebeta = - (JMz / j_om_mu) * IA1TEvmd

      !Ex, Ey in original x-y coordinate system
      Ex(recidx) = Ex(recidx) - sinbeta*Ebeta
      Ey(recidx) = Ey(recidx) + cosbeta*Ebeta

      !no Ez for VMD source

      !Hz in original x-y coordinate system
      Hz(recidx) = Hz(recidx) + (JMz/ommusq) * IA0TEvmd
    enddo !irec
  enddo  !source elements

endsubroutine interp_intvals_vmd_allcomp


!------------------------------------------------------------
!  1D EM subroutine interp_intvals_vmd_Exy, Ex and/or Ey only
!
!  get interpolated field values at receiver locations, VMD source
!  NO change required for VTI since VMD fields do not depend on epsv
!  Rita Streich 2009-2011
!------------------------------------------------------------
subroutine interp_intvals_vmd_Exy(refl_var,src,ifreq,sz,zr,bgdat,Ex,Ey,j_om_mu, funcA1TEvmd,ilay)

  implicit none

  !external variables
  type(refl_struct)               :: refl_var   !everything needed throughout 1D computations
  type(sorec),intent(in)          :: src        !source definition (we need the currents here)
  integer(kind=int32),intent(in)  :: ifreq      !frequency index
  real(kind=real64),intent(in)    :: sz,zr      !source and receiver depth
  type(backgrounddata) :: bgdat      !coordinate vectors and output EM fields
  complex(kind=real64),dimension(:)   :: Ex,Ey  !electric and magnetic field vectors
  complex(kind=real64),intent(in) :: j_om_mu    !j * omega * mu0
  complex(kind=real64),external   :: funcA1TEvmd
  integer(kind=int32),intent(in)  :: ilay       !layer index for derivatives - leave at zero for forward modeling

  !internal variables
  integer(kind=int32)   :: isrc    !source element counter
  real(kind=real64)     :: x,y,r   !temp source-receiver distances
  integer(kind=int32)   :: irec    !receiver counter
  real(kind=real64)     :: beta    !temp angle

  complex(kind=real64)          :: IA1TEvmd !interpolated integral values

  logical,dimension(nintVMDdvti)    :: wellbehaved  !indicates if Hankel integration can be used
  logical                       :: sz_eq_zr     !indicates if source and reveicer are at the same depth
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")

  complex(kind=real64)  :: Ebeta                !temp field values in cylindrical coordinates
  real(kind=real64)     :: cosbeta,sinbeta      !cos(beta) and sin(beta), precompute for efficiency
  complex(kind=real64)  :: cur      !temp source current
  complex(kind=real64)  :: JMz      !source current times constants
  integer(kind=int32)   :: idx      !source element index
  integer(kind=int32)   :: recidx   !receiver index

  !indicators for fast Hankel transform or adaptive integration
  wellbehaved = .TRUE.
  sz_eq_zr = .false.
  if(sz.eq.zr) then
    wellbehaved(1:2) = .false.
    sz_eq_zr = .TRUE.
  endif


  do isrc = refl_var%isrcstart,refl_var%isrcend

    !get source current times constant factors
    idx = refl_var%isrcperz(isrc)
    cur = conjg(src%cur(idx,ifreq))
    JMz = - cur * src%akz(idx) * j_om_mu / dfourpi


    !same positions for Ex and Ey
    exy_equalpos: if(bgdat%nExy.gt.0) then

    do irec=refl_var%irecstart,refl_var%irecend
      recidx = refl_var%irecperzExy(irec)

      y = bgdat%Exypos(recidx,2) - refl_var%ys(isrc)
      x = bgdat%Exypos(recidx,1) - refl_var%xs(isrc)
      r = sqrt(x**2 + y**2)

      beta = atan2(y,x)
      cosbeta = cos(beta)
      sinbeta = sin(beta)

      !check case r=0 only once...
      r_is_zero: if(r .eq. 0._real64) then

        !quick & dirty: skip the point if receiver is right at source point
        if(sz_eq_zr) then
          if(refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle receiver right at source point yet!'
          cycle
        endif

        IA1TEvmd = 0._real64
      else !r is not zero
        smallr: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .TRUE.
          call prepare_refcoef(refl_var,r,vmd,aniso)

          ibesord = 1
          IA1TEvmd = compute_1val(funcA1TEvmd,r,sz,zr,ibesord,wellbehaved(1))
        else !r larger than threshold radius for spline interpolation
          IA1TEvmd = splinterp_1val(refl_var,1,r)
        endif smallr
      endif r_is_zero

      Ebeta = - (JMz / j_om_mu) * IA1TEvmd

      !Ex, Ey in original x-y coordinate system
      Ex(recidx) = Ex(recidx) - sinbeta*Ebeta
      Ey(recidx) = Ey(recidx) + cosbeta*Ebeta
    enddo !irec

    else

      have_ex: if(bgdat%nEx .gt. 0) then

    do irec=refl_var%irecstart,refl_var%irecend
      recidx = refl_var%irecperzExy(irec)

      y = bgdat%Expos(recidx,2) - refl_var%ys(isrc)
      x = bgdat%Expos(recidx,1) - refl_var%xs(isrc)
      r = sqrt(x**2 + y**2)

      beta = atan2(y,x)
      sinbeta = sin(beta)

      !check case r=0 only once...
      r_is_zeroex: if(r .eq. 0._real64) then

        !quick & dirty: skip the point if receiver is right at source point
        if(sz_eq_zr) then
          if(refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle receiver right at source point yet!'
          cycle
        endif

        IA1TEvmd = 0._real64
      else !r is not zero
        smallrex: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .TRUE.
          call prepare_refcoef(refl_var,r,vmd,aniso)

          ibesord = 1
          IA1TEvmd = compute_1val(funcA1TEvmd,r,sz,zr,ibesord,wellbehaved(1))
        else !r larger than threshold radius for spline interpolation
          IA1TEvmd = splinterp_1val(refl_var,1,r)
        endif smallrex
      endif r_is_zeroex

      Ebeta = - (JMz / j_om_mu) * IA1TEvmd

      !Ex, Ey in original x-y coordinate system
      Ex(recidx) = Ex(recidx) - sinbeta*Ebeta
    enddo !irec

      endif have_ex

      have_ey: if(bgdat%nEy .gt. 0) then

    do irec=refl_var%irecstart,refl_var%irecend
      recidx = refl_var%irecperzExy(irec)

      y = bgdat%Eypos(recidx,2) - refl_var%ys(isrc)
      x = bgdat%Eypos(recidx,1) - refl_var%xs(isrc)
      r = sqrt(x**2 + y**2)

      beta = atan2(y,x)
      cosbeta = cos(beta)

      !check case r=0 only once...
      r_is_zeroey: if(r .eq. 0._real64) then

        !quick & dirty: skip the point if receiver is right at source point
        if(sz_eq_zr) then
          if(refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle receiver right at source point yet!'
          cycle
        endif

        IA1TEvmd = 0._real64
      else !r is not zero
        smallrey: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .TRUE.
          call prepare_refcoef(refl_var,r,vmd,aniso)

          ibesord = 1
          IA1TEvmd = compute_1val(funcA1TEvmd,r,sz,zr,ibesord,wellbehaved(1))
        else !r larger than threshold radius for spline interpolation
          IA1TEvmd = splinterp_1val(refl_var,1,r)
        endif smallrey
      endif r_is_zeroey

      Ebeta = - (JMz / j_om_mu) * IA1TEvmd

      !Ex, Ey in original x-y coordinate system
      Ey(recidx) = Ey(recidx) + cosbeta*Ebeta
    enddo !irec

      endif have_ey
    endif exy_equalpos

  enddo  !source elements

endsubroutine interp_intvals_vmd_Exy


!------------------------------------------------------------
!  1D EM subroutine interp_intvals_vmd_Hxy
!
!  get interpolated field values at receiver locations, VMD source
!  NO change required for VTI since VMD fields do not depend on epsv
!  Rita Streich 2009-2011
!------------------------------------------------------------
subroutine interp_intvals_vmd_Hxy(refl_var,src,ifreq,sz,zr,bgdat,Hx,Hy,j_om_mu, funcD1TEvmd,ilay)

  implicit none

  !external variables
  type(refl_struct)               :: refl_var   !everything needed throughout 1D computations
  type(sorec),intent(in)          :: src        !source definition (we need the currents here)
  integer(kind=int32),intent(in)  :: ifreq      !frequency index
  real(kind=real64),intent(in)    :: sz,zr      !source and receiver depth
  type(backgrounddata) :: bgdat      !coordinate vectors and output EM fields
  complex(kind=real64),dimension(:)   :: Hx,Hy  !electric and magnetic field vectors
  complex(kind=real64),intent(in) :: j_om_mu    !j * omega * mu0
  complex(kind=real64),external   :: funcD1TEvmd
  integer(kind=int32),intent(in)  :: ilay       !layer index for derivatives - leave at zero for forward modeling

  !internal variables
  integer(kind=int32)   :: isrc    !source element counter
  real(kind=real64)     :: x,y,r   !temp source-receiver distances
  integer(kind=int32)   :: irec    !receiver counter
  real(kind=real64)     :: beta    !temp angle

  complex(kind=real64)          :: ID1TEvmd !interpolated integral values

  logical,dimension(nintVMDdvti)    :: wellbehaved  !indicates if Hankel integration can be used
  logical                       :: sz_eq_zr     !indicates if source and reveicer are at the same depth
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")

  complex(kind=real64)  :: Hr                   !temp field values in cylindrical coordinates
  real(kind=real64)     :: cosbeta,sinbeta      !cos(beta) and sin(beta), precompute for efficiency
  complex(kind=real64)  :: cur      !temp source current
  complex(kind=real64)  :: JMz      !source current times constants
  integer(kind=int32)   :: idx      !source element index
  integer(kind=int32)   :: recidx   !receiver index

  !indicators for fast Hankel transform or adaptive integration
  wellbehaved = .TRUE.
  sz_eq_zr = .false.
  if(sz.eq.zr) then
    wellbehaved(1:2) = .false.
    sz_eq_zr = .TRUE.
  endif


  do isrc = refl_var%isrcstart,refl_var%isrcend

    !get source current times constant factors
    idx = refl_var%isrcperz(isrc)
    cur = conjg(src%cur(idx,ifreq))
    JMz = - cur * src%akz(idx) * j_om_mu / dfourpi


    !same positions for Hx and Hy
    hxy_equalpos: if(bgdat%nExy.gt.0) then

    do irec=refl_var%irecstart,refl_var%irecend
      recidx = refl_var%irecperzHxy(irec)

      y = bgdat%Hxypos(recidx,2) - refl_var%ys(isrc)
      x = bgdat%Hxypos(recidx,1) - refl_var%xs(isrc)
      r = sqrt(x**2 + y**2)

      beta = atan2(y,x)
      cosbeta = cos(beta)
      sinbeta = sin(beta)

      !check case r=0 only once...
      r_is_zero: if(r .eq. 0._real64) then

        !quick & dirty: skip the point if receiver is right at source point
        if(sz_eq_zr) then
          if(refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle receiver right at source point yet!'
          cycle
        endif

        ID1TEvmd = 0._real64
      else !r is not zero
        smallr: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .TRUE.
          call prepare_refcoef(refl_var,r,vmd,aniso)

          ibesord = 1
          ID1TEvmd = compute_1val(funcD1TEvmd,r,sz,zr,ibesord,wellbehaved(2))
        else !r larger than threshold radius for spline interpolation
          ID1TEvmd = splinterp_1val(refl_var,2,r)
        endif smallr
      endif r_is_zero

      Hr = JMz / j_om_mu * ID1TEvmd

      !Hx, Hy in original x-y coordinate system
      Hx(recidx) = Hx(recidx) + cosbeta*Hr
      Hy(recidx) = Hy(recidx) + sinbeta*Hr
    enddo !irec

    else

      have_hx: if(bgdat%nHx .gt. 0) then

    do irec=refl_var%irecstart,refl_var%irecend
      recidx = refl_var%irecperzHxy(irec)

      y = bgdat%Hxpos(recidx,2) - refl_var%ys(isrc)
      x = bgdat%Hxpos(recidx,1) - refl_var%xs(isrc)
      r = sqrt(x**2 + y**2)

      beta = atan2(y,x)
      cosbeta = cos(beta)

      !check case r=0 only once...
      r_is_zerohx: if(r .eq. 0._real64) then

        !quick & dirty: skip the point if receiver is right at source point
        if(sz_eq_zr) then
          if(refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle receiver right at source point yet!'
          cycle
        endif

        ID1TEvmd = 0._real64
      else !r is not zero
        smallrhx: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .TRUE.
          call prepare_refcoef(refl_var,r,vmd,aniso)

          ibesord = 1
          ID1TEvmd = compute_1val(funcD1TEvmd,r,sz,zr,ibesord,wellbehaved(2))
        else !r larger than threshold radius for spline interpolation
          ID1TEvmd = splinterp_1val(refl_var,2,r)
        endif smallrhx
      endif r_is_zerohx

      Hr = JMz / j_om_mu * ID1TEvmd

      !Hx, Hy in original x-y coordinate system
      Hx(recidx) = Hx(recidx) + cosbeta*Hr
    enddo !irec

      endif have_hx

      have_hy: if(bgdat%nHy .gt. 0) then

    do irec=refl_var%irecstart,refl_var%irecend
      recidx = refl_var%irecperzHxy(irec)

      y = bgdat%Hypos(recidx,2) - refl_var%ys(isrc)
      x = bgdat%Hypos(recidx,1) - refl_var%xs(isrc)
      r = sqrt(x**2 + y**2)

      beta = atan2(y,x)
      sinbeta = sin(beta)

      !check case r=0 only once...
      r_is_zerohy: if(r .eq. 0._real64) then

        !quick & dirty: skip the point if receiver is right at source point
        if(sz_eq_zr) then
          if(refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle receiver right at source point yet!'
          cycle
        endif

        ID1TEvmd = 0._real64
      else !r is not zero
        smallrhy: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .TRUE.
          call prepare_refcoef(refl_var,r,vmd,aniso)

          ibesord = 1
          ID1TEvmd = compute_1val(funcD1TEvmd,r,sz,zr,ibesord,wellbehaved(2))
        else !r larger than threshold radius for spline interpolation
          ID1TEvmd = splinterp_1val(refl_var,2,r)
        endif smallrhy
      endif r_is_zerohy

      Hr = JMz / j_om_mu * ID1TEvmd

      !Hx, Hy in original x-y coordinate system
      Hy(recidx) = Hy(recidx) + sinbeta*Hr
    enddo !irec

      endif have_hy
    endif hxy_equalpos

  enddo  !source elements

endsubroutine interp_intvals_vmd_Hxy


!------------------------------------------------------------
!  1D EM subroutine interp_intvals_vmd_Hz, Hz only
!
!  get interpolated field values at receiver locations, VMD source
!  NO change required for VTI since VMD fields do not depend on epsv
!  Rita Streich 2009-2011
!------------------------------------------------------------
subroutine interp_intvals_vmd_Hz(refl_var,src,ifreq,sz,zr,bgdat,Hz,j_om_mu,ommusq, funcA0TEvmd,ilay)

  implicit none

  !external variables
  type(refl_struct)               :: refl_var   !everything needed throughout 1D computations
  type(sorec),intent(in)          :: src        !source definition (we need the currents here)
  integer(kind=int32),intent(in)  :: ifreq      !frequency index
  real(kind=real64),intent(in)    :: sz,zr      !source and receiver depth
  type(backgrounddata) :: bgdat      !coordinate vectors and output EM fields
  complex(kind=real64),dimension(:)   :: Hz  !electric and magnetic field vectors
  complex(kind=real64),intent(in) :: j_om_mu    !j * omega * mu0
  real(kind=real64),intent(in)    :: ommusq     ! (omega * mu0) **2
  complex(kind=real64),external   :: funcA0TEvmd
  integer(kind=int32),intent(in)  :: ilay       !layer index for derivatives - leave at zero for forward modeling

  !internal variables
  integer(kind=int32)   :: isrc    !source element counter
  real(kind=real64)     :: x,y,r   !temp source-receiver distances
  integer(kind=int32)   :: irec    !receiver counter
  real(kind=real64)     :: beta    !temp angle

  complex(kind=real64)          :: IA0TEvmd !interpolated integral values

  logical,dimension(nintVMDdvti)    :: wellbehaved  !indicates if Hankel integration can be used
  logical                       :: sz_eq_zr     !indicates if source and reveicer are at the same depth
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")

  real(kind=real64)     :: cosbeta,sinbeta      !cos(beta) and sin(beta), precompute for efficiency
  complex(kind=real64)  :: cur      !temp source current
  complex(kind=real64)  :: JMz      !source current times constants
  integer(kind=int32)   :: idx      !source element index
  integer(kind=int32)   :: recidx   !receiver index

  !indicators for fast Hankel transform or adaptive integration
  wellbehaved = .TRUE.
  sz_eq_zr = .false.
  if(sz.eq.zr) then
    wellbehaved(1:2) = .false.
    sz_eq_zr = .TRUE.
  endif


  do isrc = refl_var%isrcstart,refl_var%isrcend

    !get source current times constant factors
    idx = refl_var%isrcperz(isrc)
    cur = conjg(src%cur(idx,ifreq))
    JMz = - cur * src%akz(idx) * j_om_mu / dfourpi

    do irec=refl_var%irecstart,refl_var%irecend
      recidx = refl_var%irecperzHz(irec)

      y = bgdat%Hzpos(recidx,2) - refl_var%ys(isrc)
      x = bgdat%Hzpos(recidx,1) - refl_var%xs(isrc)
      r = sqrt(x**2 + y**2)

      beta = atan2(y,x)
      cosbeta = cos(beta)
      sinbeta = sin(beta)

      !check case r=0 only once...
      r_is_zero: if(r .eq. 0._real64) then

        !quick & dirty: skip the point if receiver is right at source point
        if(sz_eq_zr) then
          if(refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle receiver right at source point yet!'
          !still add special contribution right at source point
          if(ilay.eq.0) Hz(recidx) = Hz(recidx) - JMz / j_om_mu
          !skip numerical integral evaluations - they would most likely fail anyway, need to look at this...
          cycle
        endif

        IA0TEvmd = compute_1valr0(funcA0TEvmd)
      else !r is not zero
        smallr: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .TRUE.
          call prepare_refcoef(refl_var,r,vmd,aniso)

          ibesord = 0
          IA0TEvmd = compute_1val(funcA0TEvmd,r,sz,zr,ibesord,wellbehaved(3))
        else !r larger than threshold radius for spline interpolation
          IA0TEvmd = splinterp_1val(refl_var,3,r)
        endif smallr
      endif r_is_zero

      !Hz in original x-y coordinate system
      Hz(recidx) = Hz(recidx) + (JMz/ommusq) * IA0TEvmd
    enddo !irec
  enddo  !source elements

endsubroutine interp_intvals_vmd_Hz


