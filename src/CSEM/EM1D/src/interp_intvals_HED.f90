!------------------------------------------------------------
!  1D EM subroutine interp_intvals_hed_allcomp
!
!  get interpolated field values at iReceiver locations, HED source,
!    all field components at the same coordinates
!
!  Rita Streich 2009-2011
!------------------------------------------------------------
subroutine interp_intvals_hed_allcomp(refl_var,src,ifreq,sz,zr,bgdat,Ex,Ey,Ez,Hx,Hy,Hz,omeps_recv,ommu, &
  funcA0TE,funcA0TM,funcA1TE,funcA1TM,funcDz1TM,funcD0TE,funcD0TM,funcD1TE,funcD1TM,funcAz1TE,ilay,funcDz1TMfwd, &
  funcA0TMv,funcA1TMv,funcDz1TMv,funcD0TMv,funcD1TMv,Exv,Eyv,Ezv,Hxv,Hyv,Hzv)

  implicit none

  !external variables
  type(refl_struct)               :: refl_var   !everything needed throughout 1D computations
  type(sorec),intent(in)          :: src        !source definition (we need the currents here)
  integer(kind=int32),intent(in)  :: ifreq      !frequency index
  real(kind=real64),intent(in)    :: sz,zr      !source and iReceiver depth
  type(backgrounddata) :: bgdat      !coordinate vectors and output EM fields
    !input separate fields because the same routine is used both for fields and derivatives!
  complex(kind=real64),dimension(:)   :: Ex,Ey,Hx,Hy,Hz    !EM fields or dervatives (horizontal/isotropic)
  complex(kind=real64),dimension(:),target   :: Ez    !EM fields or dervatives (horizontal/isotropic)
  complex(kind=real64),intent(in) :: omeps_recv !omega * epsilon in iReceiver layer
  real(kind=real64),intent(in)    :: ommu       !omega * mu0
  complex(kind=real64),external   :: funcA0TE,funcA0TM,funcA1TE,funcA1TM,funcDz1TM,funcD0TE,funcD0TM,funcD1TE,funcD1TM,funcAz1TE
  integer(kind=int32),intent(in)  :: ilay       !layer index for derivatives - leave at zero for forward modeling
  complex(kind=real64),external,optional   :: funcDz1TMfwd  !function only needed for derivatives in iReceiver layer
  complex(kind=real64),external,optional   :: funcA0TMv,funcA1TMv,funcDz1TMv,funcD0TMv,funcD1TMv !functions for TM epsv deriv.
  complex(kind=real64),dimension(:),optional   :: Exv,Eyv,Hxv,Hyv,Hzv !EM field for epsv derivatives
  complex(kind=real64),dimension(:),target,optional   :: Ezv !EM field for epsv derivatives

  !internal variables
  integer(kind=int32)   :: isrc    !source element counter
  real(kind=real64)     :: x,y,r   !temp source-iReceiver distances
  integer(kind=int32)   :: irec    !iReceiver counter
  real(kind=real64)     :: beta,betarot      !temp angles

  complex(kind=real64)          :: IA0TE,IA0TM,IA1TE,IA1TM,IDz1TM !interpolated integral values
  complex(kind=real64)          :: ID0TE,ID0TM,ID1TE,ID1TM,IAz1TE

  logical,dimension(nintHEDdvti)    :: wellbehaved  !indicates if Hankel integration can be used
  logical                       :: sz_eq_zr     !indicates if source and reveicer are at the same depth
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")

  complex(kind=real64)  :: Er,Ebeta,Hr,Hbeta    !temp field values in cylindrical coordinates
  real(kind=real64)     :: cosbeta,sinbeta      !cos(beta) and sin(beta), precompute for efficiency
  real(kind=real64)     :: cosbetarot,sinbetarot      !cos(betarot) and sin(betarot), precompute for efficiency
  complex(kind=real64)  :: cur      !temp source current
  complex(kind=real64)  :: JEh      !source current times constants
  integer(kind=int32)   :: idx      !source element index
  integer(kind=int32)   :: recidx   !iReceiver index
  !flag for computing epsv derivatives, not needed for forward computation, so "aniso" value cannot be used here
  logical               :: with_dvert
  complex(kind=real64),dimension(:),pointer   :: Erec  !points to Ez for isotropic and Ezv for VTI case

  !indicators for fast Hankel transform or adaptive integration
  wellbehaved = .true.
  sz_eq_zr = .false.
  if (sz.eq.zr) then
    wellbehaved(2) = .false.
    wellbehaved(5:7) = .false.
    sz_eq_zr = .true.
  endif
  
  if (present(funcD1TMv)) then
    with_dvert = .true.
    if (sz_eq_zr) then
      wellbehaved(13) = .false.
      wellbehaved(15:16) = .false.
    endif
  else
    with_dvert = .false.
  endif


  do isrc = refl_var%isrcstart,refl_var%isrcend

    !get source current times constant factors
    !complex comjugate of currents was taken when reading in source currents, so no need to do that again here
    idx = refl_var%isrcperz(isrc)
    cur = conjg(src%cur(idx,ifreq))
    JEh = cur * sqrt(src%ljx(idx)**2 + src%ljy(idx)**2) / dfourpi

    do irec=refl_var%irecstart,refl_var%irecend
      recidx = refl_var%irecperzExy(irec)

      y = bgdat%Exypos(recidx,2) - refl_var%ys(isrc)
      x = bgdat%Exypos(recidx,1) - refl_var%xs(isrc)

      r = sqrt(x**2 + y**2)

      beta = atan2(y,x)
      cosbeta = cos(beta)
      sinbeta = sin(beta)
      betarot = beta - refl_var%betasrc(isrc)
      cosbetarot = cos(betarot)
      sinbetarot = sin(betarot)


      r_is_zero: if (r.eq.0._real64) then
        !quick & dirty: skip the point if iReceiver is right at source point
        if (sz_eq_zr) then
          if (refl_var%infolevel.ge.output_more) then
            write(*,'(a)') 'WARNING: cannot handle iReceiver right at source point yet!'
          endif
          cycle
        endif

        !Bessel function J0(r=0) = 1, J1(r=0) = 0
        !--> need to evaluate zero-order integrals only
        !use adaptive integration
        IA0TE = compute_1valr0(funcA0TE)
        IA0TM = compute_1valr0(funcA0TM)
        ID0TE = compute_1valr0(funcD0TE)
        ID0TM = compute_1valr0(funcD0TM)

        !limits from (1/R) * integral(bessel1) become bessel0 ...
        Er = - JEh * cosbetarot * (0.5_real64 * (IA0TE + IA0TM))
        Ebeta = JEh * sinbetarot * (0.5_real64 * (IA0TE + IA0TM))
        Hr = - JEh * sinbetarot * (0.5_real64 * (ID0TE + ID0TM))
        Hbeta = - JEh * cosbetarot * (0.5_real64 * (ID0TE + ID0TM))
        
      else !r is not zero
        smallr: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .true.
          call prepare_refcoef(refl_var,r,hed,aniso)

          !evaluate all integrals just for this radius
          ibesord = 0
          IA0TE = compute_1val(funcA0TE,r,sz,zr,ibesord,wellbehaved(1))
          IA0TM = compute_1val(funcA0TM,r,sz,zr,ibesord,wellbehaved(2))
          ibesord = 1
          IA1TE = compute_1val(funcA1TE,r,sz,zr,ibesord,wellbehaved(3))
          IA1TM = compute_1val(funcA1TM,r,sz,zr,ibesord,wellbehaved(4))

          ibesord = 1
          IDz1TM = compute_1val(funcDz1TM,r,sz,zr,ibesord,wellbehaved(5))

          ibesord = 0
          ID0TE = compute_1val(funcD0TE,r,sz,zr,ibesord,wellbehaved(6))
          ID0TM = compute_1val(funcD0TM,r,sz,zr,ibesord,wellbehaved(7))
          ibesord = 1
          ID1TE = compute_1val(funcD1TE,r,sz,zr,ibesord,wellbehaved(8))
          ID1TM = compute_1val(funcD1TM,r,sz,zr,ibesord,wellbehaved(9))

          ibesord = 1
          IAz1TE = compute_1val(funcAz1TE,r,sz,zr,ibesord,wellbehaved(10))
          
        else !r larger than threshold radius for spline interpolation

          !get integral values for this radius by interpolation from precomputed integral values
          IA0TE = splinterp_1val(refl_var,1,r)
          IA0TM = splinterp_1val(refl_var,2,r)
          IA1TE = splinterp_1val(refl_var,3,r)
          IA1TM = splinterp_1val(refl_var,4,r)
          IDz1TM = splinterp_1val(refl_var,5,r)

          ID0TE = splinterp_1val(refl_var,6,r)
          ID0TM = splinterp_1val(refl_var,7,r)
          ID1TE = splinterp_1val(refl_var,8,r)
          ID1TM = splinterp_1val(refl_var,9,r)
          IAz1TE = splinterp_1val(refl_var,10,r)
          
        endif smallr

        !we get Ex and Ey from these
        Er = -JEh * cosbetarot * (IA0TM + (IA1TE - IA1TM)/r)
        Ebeta = -JEh * sinbetarot * (-IA0TE + (IA1TE - IA1TM)/r)

        !Ez is zero for zero radius
        Ez(recidx) = Ez(recidx) + JEh * dci * cosbetarot * IDz1TM / omeps_recv

        Hr = JEh * sinbetarot * (-ID0TE + (ID1TE - ID1TM)/r)
        Hbeta = -JEh * cosbetarot * (ID0TM + (ID1TE - ID1TM)/r)

        !Hz in original x-y coordinate system - no special case needed for radius zero
        Hz(recidx) = Hz(recidx) + JEh * dci * sinbetarot * IAz1TE / ommu

      endif r_is_zero

      !field values rotated back to original x-y coordinate system
      !don't take complex conjugate here, but at the very end!
      Ex(recidx) = Ex(recidx) + cosbeta*Er - sinbeta*Ebeta
      Ey(recidx) = Ey(recidx) + sinbeta*Er + cosbeta*Ebeta
      Hx(recidx) = Hx(recidx) + cosbeta*Hr - sinbeta*Hbeta
      Hy(recidx) = Hy(recidx) + sinbeta*Hr + cosbeta*Hbeta

      dvert: if (with_dvert) then
        r_is_zerov: if (r.eq.0._real64) then
          !not necessary to check if source and iReceiver coincide exactly
          !  since in that case the above check already effects going to next iReceiver
          !if (sz_eq_zr) cycle

          !Bessel function J0(r=0) = 1, J1(r=0) = 0
          !--> need to evaluate zero-order integrals only
          !use adaptive integration
          IA0TM = compute_1valr0(funcA0TMv)
          ID0TM = compute_1valr0(funcD0TMv)

          !limits from (1/R) * integral(bessel1) become bessel0 ...
          Er = - JEh * cosbetarot * 0.5_real64 * IA0TM
          Ebeta = JEh * sinbetarot * 0.5_real64 * IA0TM
          Hr = - JEh * sinbetarot * 0.5_real64 * ID0TM
          Hbeta = - JEh * cosbetarot * 0.5_real64 * ID0TM
        
        else !r is not zero
          smallrv: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,hed,aniso)

            !evaluate all integrals just for this radius
            ibesord = 0
            IA0TM = compute_1val(funcA0TMv,r,sz,zr,ibesord,wellbehaved(13))
            ibesord = 1
            IA1TM = compute_1val(funcA1TMv,r,sz,zr,ibesord,wellbehaved(12))
            IDz1TM = compute_1val(funcDz1TMv,r,sz,zr,ibesord,wellbehaved(16))

            ibesord = 0
            ID0TM = compute_1val(funcD0TMv,r,sz,zr,ibesord,wellbehaved(15))
            ibesord = 1
            ID1TM = compute_1val(funcD1TMv,r,sz,zr,ibesord,wellbehaved(14))

          else !r larger than threshold radius for spline interpolation

            !get integral values for this radius by interpolation from precomputed integral values
            IA0TM = splinterp_1val(refl_var,13,r)
            IA1TM = splinterp_1val(refl_var,12,r)
            IDz1TM = splinterp_1val(refl_var,16,r)

            ID0TM = splinterp_1val(refl_var,15,r)
            ID1TM = splinterp_1val(refl_var,14,r)
          endif smallrv

          !we get Ex and Ey from these
          Er = -JEh * cosbetarot * (IA0TM - IA1TM/r)
          Ebeta = -JEh * sinbetarot * (- IA1TM/r)

          !Ez is zero for zero radius
          Ezv(recidx) = Ezv(recidx) + JEh * dci * cosbetarot * IDz1TM / omeps_recv

          Hr = JEh * sinbetarot * (- ID1TM/r)
          Hbeta = -JEh * cosbetarot * (ID0TM - ID1TM/r)
          
          !deriv. dHz/depsv is zero since Hz does not contain epsv
        endif r_is_zerov

        !field values rotated back to original x-y coordinate system
        !don't take complex conjugate here, but at the very end!
        Exv(recidx) = Exv(recidx) + cosbeta*Er - sinbeta*Ebeta
        Eyv(recidx) = Eyv(recidx) + sinbeta*Er + cosbeta*Ebeta
        Hxv(recidx) = Hxv(recidx) + cosbeta*Hr - sinbeta*Hbeta
        Hyv(recidx) = Hyv(recidx) + sinbeta*Hr + cosbeta*Hbeta
      endif dvert

    enddo !irec
  enddo  !source elements
  
  
  !special term is needed for derivative of Ez in iReceiver layer
  !isotropic: add this to deriv. for overall complex permittivity --> E
  !VTI: this term does not exist for epsh deriv, but exists for epsv deriv --> nothing added to E, add to Ev
  deriv_ilayrec: if (ilay .eq. ilayrec) then
  
    if (.not. present(funcDz1TMfwd)) then
      write(*,'(a)') 'ERROR: function for Ez in iReceiver layer not given, cannot compute correct Ez derivative!'
      return
    endif
  
    !there is a factor epsv in the term before the Ez integral
    !-> for isotropic case, epsv becomes eps and contribution is added to E
    !-> for VTI case, derivative of this term only exists for epsv -> add to Ev, nothing added to E
    if (with_dvert) then
      Erec => Ezv
    else
      Erec => Ez
    endif
  
    !indicators for fast Hankel transform or adaptive integration
    wellbehaved(11) = wellbehaved(5)  !assume that derivative and fwd integrals behave equally well or badly - not tested!

    do isrc = refl_var%isrcstart,refl_var%isrcend

      !get source current times constant factors
      !complex comjugate of currents was taken when reading in source currents, so no need to do that again here
      idx = refl_var%isrcperz(isrc)
      cur = conjg(src%cur(idx,ifreq))
      JEh = cur * sqrt(src%ljx(idx)**2 + src%ljy(idx)**2) / dfourpi

      do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzExy(irec)

        y = bgdat%Exypos(recidx,2) - refl_var%ys(isrc)
        x = bgdat%Exypos(recidx,1) - refl_var%xs(isrc)

        r = sqrt(x**2 + y**2)

        if (r.eq.0._real64) then
          !Ez is zero for zero radius
          cycle
        else

          beta = atan2(y,x)
          betarot = beta - refl_var%betasrc(isrc)
          cosbetarot = cos(betarot)

          if (r.lt.rsplmin) then
            !reflection coeff. for this radius
            call prepare_refcoef(refl_var,r,ved,aniso) !use "ved" to compute TM refl. coeff. only
            ibesord = 1
            IDz1TM = compute_1val(funcDz1TMfwd,r,sz,zr,ibesord,wellbehaved(11))
          else
            IDz1TM = splinterp_1val(refl_var,11,r)
          endif

          Erec(recidx) = Erec(recidx) - JEh * dci * cosbetarot * IDz1TM / (omeps_recv * epsv(ilayrec))

        endif

      enddo !irec
    enddo  !source elements
  endif deriv_ilayrec

endsubroutine interp_intvals_hed_allcomp


!------------------------------------------------------------
!  1D EM subroutine interp_intvals_hed_Exy
!
!  get interpolated field values at iReceiver locations, HED source,
!    Ex / Ey only
!
!  Rita Streich 2009-2011
!------------------------------------------------------------
subroutine interp_intvals_hed_Exy(refl_var,src,ifreq,sz,zr,bgdat,Ex,Ey, funcA0TE,funcA0TM,funcA1TE,funcA1TM,ilay, &
  funcA0TMv,funcA1TMv,Exv,Eyv)

  implicit none

  !external variables
  type(refl_struct)               :: refl_var   !everything needed throughout 1D computations
  type(sorec),intent(in)          :: src        !source definition (we need the currents here)
  integer(kind=int32),intent(in)  :: ifreq      !frequency index
  real(kind=real64),intent(in)    :: sz,zr      !source and iReceiver depth
  type(backgrounddata) :: bgdat      !coordinate vectors and output EM fields
    !input separate fields because the same routine is used both for fields and derivatives!
  complex(kind=real64),dimension(:)   :: Ex,Ey    !EM fields or dervatives (horizontal/isotropic)
  complex(kind=real64),external   :: funcA0TE,funcA0TM,funcA1TE,funcA1TM
  integer(kind=int32),intent(in)  :: ilay       !layer index for derivatives - leave at zero for forward modeling
  complex(kind=real64),external,optional   :: funcA0TMv,funcA1TMv !functions for TM epsv deriv.
  complex(kind=real64),dimension(:),optional   :: Exv,Eyv !EM field for epsv derivatives

  !internal variables
  integer(kind=int32)   :: isrc    !source element counter
  real(kind=real64)     :: x,y,r   !temp source-iReceiver distances
  integer(kind=int32)   :: irec    !iReceiver counter
  real(kind=real64)     :: beta,betarot      !temp angles

  complex(kind=real64)          :: IA0TE,IA0TM,IA1TE,IA1TM !interpolated integral values

  logical,dimension(nintHEDdvti)    :: wellbehaved  !indicates if Hankel integration can be used
  logical                       :: sz_eq_zr     !indicates if source and reveicer are at the same depth
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")

  complex(kind=real64)  :: Er,Ebeta    !temp field values in cylindrical coordinates
  real(kind=real64)     :: cosbeta,sinbeta      !cos(beta) and sin(beta), precompute for efficiency
  real(kind=real64)     :: cosbetarot,sinbetarot      !cos(betarot) and sin(betarot), precompute for efficiency
  complex(kind=real64)  :: cur      !temp source current
  complex(kind=real64)  :: JEh      !source current times constants
  integer(kind=int32)   :: idx      !source element index
  integer(kind=int32)   :: recidx   !iReceiver index
  !flag for computing epsv derivatives, not needed for forward computation, so "aniso" value cannot be used here
  logical               :: with_dvert

  !indicators for fast Hankel transform or adaptive integration
  wellbehaved = .true.
  sz_eq_zr = .false.
  if (sz.eq.zr) then
    wellbehaved(2) = .false.
    wellbehaved(5:7) = .false.
    sz_eq_zr = .true.
  endif
  
  if (present(funcA0TMv)) then
    with_dvert = .true.
    if (sz_eq_zr) then
      wellbehaved(13) = .false.
      wellbehaved(15:16) = .false.
    endif
  else
    with_dvert = .false.
  endif


  do isrc = refl_var%isrcstart,refl_var%isrcend

    !get source current times constant factors
    !complex comjugate of currents was taken when reading in source currents, so no need to do that again here
    idx = refl_var%isrcperz(isrc)
    cur = conjg(src%cur(idx,ifreq))
    JEh = cur * sqrt(src%ljx(idx)**2 + src%ljy(idx)**2) / dfourpi

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
        betarot = beta - refl_var%betasrc(isrc)
        cosbetarot = cos(betarot)
        sinbetarot = sin(betarot)

        r_is_zero: if (r.eq.0._real64) then
          !quick & dirty: skip the point if iReceiver is right at source point
          if (sz_eq_zr) then
            if (refl_var%infolevel.ge.output_more) then
              write(*,'(a)') 'WARNING: cannot handle iReceiver right at source point yet!'
            endif
            cycle
          endif

          !Bessel function J0(r=0) = 1, J1(r=0) = 0
          !--> need to evaluate zero-order integrals only
          !use adaptive integration
          IA0TE = compute_1valr0(funcA0TE)
          IA0TM = compute_1valr0(funcA0TM)

          !limits from (1/R) * integral(bessel1) become bessel0 ...
          Er = - JEh * cosbetarot * (0.5_real64 * (IA0TE + IA0TM))
          Ebeta = JEh * sinbetarot * (0.5_real64 * (IA0TE + IA0TM))

        else !r is not zero
          smallr: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,hed,aniso)

            !evaluate all integrals just for this radius
            ibesord = 0
            IA0TE = compute_1val(funcA0TE,r,sz,zr,ibesord,wellbehaved(1))
            IA0TM = compute_1val(funcA0TM,r,sz,zr,ibesord,wellbehaved(2))
            ibesord = 1
            IA1TE = compute_1val(funcA1TE,r,sz,zr,ibesord,wellbehaved(3))
            IA1TM = compute_1val(funcA1TM,r,sz,zr,ibesord,wellbehaved(4))

          else !r larger than threshold radius for spline interpolation
            !get integral values for this radius by interpolation from precomputed integral values
            IA0TE = splinterp_1val(refl_var,1,r)
            IA0TM = splinterp_1val(refl_var,2,r)
            IA1TE = splinterp_1val(refl_var,3,r)
            IA1TM = splinterp_1val(refl_var,4,r)
          endif smallr

          !we get Ex and Ey from these
          Er = -JEh * cosbetarot * (IA0TM + (IA1TE - IA1TM)/r)
          Ebeta = -JEh * sinbetarot * (-IA0TE + (IA1TE - IA1TM)/r)
        endif r_is_zero

        !field values rotated back to original x-y coordinate system
        !don't take complex conjugate here, but at the very end!
        Ex(recidx) = Ex(recidx) + cosbeta*Er - sinbeta*Ebeta
        Ey(recidx) = Ey(recidx) + sinbeta*Er + cosbeta*Ebeta

        dvert: if (with_dvert) then
          r_is_zerov: if (r.eq.0._real64) then
            !not necessary to check if source and iReceiver coincide exactly
            !  since in that case the above check already effects going to next iReceiver
            !if (sz_eq_zr) cycle

            !Bessel function J0(r=0) = 1, J1(r=0) = 0
            !--> need to evaluate zero-order integrals only
            !use adaptive integration
            IA0TM = compute_1valr0(funcA0TMv)

            !limits from (1/R) * integral(bessel1) become bessel0 ...
            Er = - JEh * cosbetarot * 0.5_real64 * IA0TM
            Ebeta = JEh * sinbetarot * 0.5_real64 * IA0TM

          else !r is not zero
            smallrv: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

              !reflection coeff. for this radius
              refl_var%refcoef_changed = .true.
              call prepare_refcoef(refl_var,r,hed,aniso)

              !evaluate all integrals just for this radius
              ibesord = 0
              IA0TM = compute_1val(funcA0TMv,r,sz,zr,ibesord,wellbehaved(13))
              ibesord = 1
              IA1TM = compute_1val(funcA1TMv,r,sz,zr,ibesord,wellbehaved(12))
            else !r larger than threshold radius for spline interpolation

              !get integral values for this radius by interpolation from precomputed integral values
              IA0TM = splinterp_1val(refl_var,13,r)
              IA1TM = splinterp_1val(refl_var,12,r)
            endif smallrv

            !we get Ex and Ey from these
            Er = -JEh * cosbetarot * (IA0TM - IA1TM/r)
            Ebeta = -JEh * sinbetarot * (- IA1TM/r)
          endif r_is_zerov

          !field values rotated back to original x-y coordinate system
          !don't take complex conjugate here, but at the very end!
          Exv(recidx) = Exv(recidx) + cosbeta*Er - sinbeta*Ebeta
          Eyv(recidx) = Eyv(recidx) + sinbeta*Er + cosbeta*Ebeta
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
        betarot = beta - refl_var%betasrc(isrc)
        cosbetarot = cos(betarot)
        sinbetarot = sin(betarot)

        r_is_zeroex: if (r.eq.0._real64) then
          !quick & dirty: skip the point if iReceiver is right at source point
          if (sz_eq_zr) then
            if (refl_var%infolevel.ge.output_more) then
              write(*,'(a)') 'WARNING: cannot handle iReceiver right at source point yet!'
            endif
            cycle
          endif

          !Bessel function J0(r=0) = 1, J1(r=0) = 0
          !--> need to evaluate zero-order integrals only
          !use adaptive integration
          IA0TE = compute_1valr0(funcA0TE)
          IA0TM = compute_1valr0(funcA0TM)

          !limits from (1/R) * integral(bessel1) become bessel0 ...
          Er = - JEh * cosbetarot * (0.5_real64 * (IA0TE + IA0TM))
          Ebeta = JEh * sinbetarot * (0.5_real64 * (IA0TE + IA0TM))

        else !r is not zero
          smallrex: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,hed,aniso)

            !evaluate all integrals just for this radius
            ibesord = 0
            IA0TE = compute_1val(funcA0TE,r,sz,zr,ibesord,wellbehaved(1))
            IA0TM = compute_1val(funcA0TM,r,sz,zr,ibesord,wellbehaved(2))
            ibesord = 1
            IA1TE = compute_1val(funcA1TE,r,sz,zr,ibesord,wellbehaved(3))
            IA1TM = compute_1val(funcA1TM,r,sz,zr,ibesord,wellbehaved(4))

          else !r larger than threshold radius for spline interpolation
            !get integral values for this radius by interpolation from precomputed integral values
            IA0TE = splinterp_1val(refl_var,1,r)
            IA0TM = splinterp_1val(refl_var,2,r)
            IA1TE = splinterp_1val(refl_var,3,r)
            IA1TM = splinterp_1val(refl_var,4,r)
          endif smallrex

          !we get Ex and Ey from these
          Er = -JEh * cosbetarot * (IA0TM + (IA1TE - IA1TM)/r)
          Ebeta = -JEh * sinbetarot * (-IA0TE + (IA1TE - IA1TM)/r)
        endif r_is_zeroex

        !field values rotated back to original x-y coordinate system
        !don't take complex conjugate here, but at the very end!
        Ex(recidx) = Ex(recidx) + cosbeta*Er - sinbeta*Ebeta

        dvertex: if (with_dvert) then
          r_is_zerovex: if (r.eq.0._real64) then
            !not necessary to check if source and iReceiver coincide exactly
            !  since in that case the above check already effects going to next iReceiver
            !if (sz_eq_zr) cycle

            !Bessel function J0(r=0) = 1, J1(r=0) = 0
            !--> need to evaluate zero-order integrals only
            !use adaptive integration
            IA0TM = compute_1valr0(funcA0TMv)

            !limits from (1/R) * integral(bessel1) become bessel0 ...
            Er = - JEh * cosbetarot * 0.5_real64 * IA0TM
            Ebeta = JEh * sinbetarot * 0.5_real64 * IA0TM

          else !r is not zero
            smallrvex: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

              !reflection coeff. for this radius
              refl_var%refcoef_changed = .true.
              call prepare_refcoef(refl_var,r,hed,aniso)

              !evaluate all integrals just for this radius
              ibesord = 0
              IA0TM = compute_1val(funcA0TMv,r,sz,zr,ibesord,wellbehaved(13))
              ibesord = 1
              IA1TM = compute_1val(funcA1TMv,r,sz,zr,ibesord,wellbehaved(12))
            else !r larger than threshold radius for spline interpolation

              !get integral values for this radius by interpolation from precomputed integral values
              IA0TM = splinterp_1val(refl_var,13,r)
              IA1TM = splinterp_1val(refl_var,12,r)
            endif smallrvex

            !we get Ex and Ey from these
            Er = -JEh * cosbetarot * (IA0TM - IA1TM/r)
            Ebeta = -JEh * sinbetarot * (- IA1TM/r)
          endif r_is_zerovex

          !field values rotated back to original x-y coordinate system
          !don't take complex conjugate here, but at the very end!
          Exv(recidx) = Exv(recidx) + cosbeta*Er - sinbeta*Ebeta
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
        betarot = beta - refl_var%betasrc(isrc)
        cosbetarot = cos(betarot)
        sinbetarot = sin(betarot)

        r_is_zeroey: if (r.eq.0._real64) then
          !quick & dirty: skip the point if iReceiver is right at source point
          if (sz_eq_zr) then
            if (refl_var%infolevel.ge.output_more) then
              write(*,'(a)') 'WARNING: cannot handle iReceiver right at source point yet!'
            endif
            cycle
          endif

          !Bessel function J0(r=0) = 1, J1(r=0) = 0
          !--> need to evaluate zero-order integrals only
          !use adaptive integration
          IA0TE = compute_1valr0(funcA0TE)
          IA0TM = compute_1valr0(funcA0TM)

          !limits from (1/R) * integral(bessel1) become bessel0 ...
          Er = - JEh * cosbetarot * (0.5_real64 * (IA0TE + IA0TM))
          Ebeta = JEh * sinbetarot * (0.5_real64 * (IA0TE + IA0TM))

        else !r is not zero
          smallrey: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,hed,aniso)

            !evaluate all integrals just for this radius
            ibesord = 0
            IA0TE = compute_1val(funcA0TE,r,sz,zr,ibesord,wellbehaved(1))
            IA0TM = compute_1val(funcA0TM,r,sz,zr,ibesord,wellbehaved(2))
            ibesord = 1
            IA1TE = compute_1val(funcA1TE,r,sz,zr,ibesord,wellbehaved(3))
            IA1TM = compute_1val(funcA1TM,r,sz,zr,ibesord,wellbehaved(4))

          else !r larger than threshold radius for spline interpolation
            !get integral values for this radius by interpolation from precomputed integral values
            IA0TE = splinterp_1val(refl_var,1,r)
            IA0TM = splinterp_1val(refl_var,2,r)
            IA1TE = splinterp_1val(refl_var,3,r)
            IA1TM = splinterp_1val(refl_var,4,r)
          endif smallrey

          !we get Ex and Ey from these
          Er = -JEh * cosbetarot * (IA0TM + (IA1TE - IA1TM)/r)
          Ebeta = -JEh * sinbetarot * (-IA0TE + (IA1TE - IA1TM)/r)
        endif r_is_zeroey

        !field values rotated back to original x-y coordinate system
        !don't take complex conjugate here, but at the very end!
        Ey(recidx) = Ey(recidx) + sinbeta*Er + cosbeta*Ebeta

        dvertey: if (with_dvert) then
          r_is_zerovey: if (r.eq.0._real64) then
            !not necessary to check if source and iReceiver coincide exactly
            !  since in that case the above check already effects going to next iReceiver
            !if (sz_eq_zr) cycle

            !Bessel function J0(r=0) = 1, J1(r=0) = 0
            !--> need to evaluate zero-order integrals only
            !use adaptive integration
            IA0TM = compute_1valr0(funcA0TMv)

            !limits from (1/R) * integral(bessel1) become bessel0 ...
            Er = - JEh * cosbetarot * 0.5_real64 * IA0TM
            Ebeta = JEh * sinbetarot * 0.5_real64 * IA0TM

          else !r is not zero
            smallrvey: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

              !reflection coeff. for this radius
              refl_var%refcoef_changed = .true.
              call prepare_refcoef(refl_var,r,hed,aniso)

              !evaluate all integrals just for this radius
              ibesord = 0
              IA0TM = compute_1val(funcA0TMv,r,sz,zr,ibesord,wellbehaved(13))
              ibesord = 1
              IA1TM = compute_1val(funcA1TMv,r,sz,zr,ibesord,wellbehaved(12))
            else !r larger than threshold radius for spline interpolation

              !get integral values for this radius by interpolation from precomputed integral values
              IA0TM = splinterp_1val(refl_var,13,r)
              IA1TM = splinterp_1val(refl_var,12,r)
            endif smallrvey

            !we get Ex and Ey from these
            Er = -JEh * cosbetarot * (IA0TM - IA1TM/r)
            Ebeta = -JEh * sinbetarot * (- IA1TM/r)
          endif r_is_zerovey

          !field values rotated back to original x-y coordinate system
          !don't take complex conjugate here, but at the very end!
          Eyv(recidx) = Eyv(recidx) + sinbeta*Er + cosbeta*Ebeta
        endif dvertey

      enddo !irec
      endif have_ey
   endif exy_equalpos

  enddo  !source elements

endsubroutine interp_intvals_hed_Exy


!------------------------------------------------------------
!  1D EM subroutine interp_intvals_hed_Ez
!
!  get interpolated field values at iReceiver locations, HED source,
!    Ez only
!
!  Rita Streich 2009-2011
!------------------------------------------------------------
subroutine interp_intvals_hed_Ez(refl_var,src,ifreq,sz,zr,bgdat,Ez,omeps_recv, funcDz1TM,ilay,funcDz1TMfwd, funcDz1TMv,Ezv)

  implicit none

  !external variables
  type(refl_struct)               :: refl_var   !everything needed throughout 1D computations
  type(sorec),intent(in)          :: src        !source definition (we need the currents here)
  integer(kind=int32),intent(in)  :: ifreq      !frequency index
  real(kind=real64),intent(in)    :: sz,zr      !source and iReceiver depth
  type(backgrounddata) :: bgdat      !coordinate vectors and output EM fields
    !input separate fields because the same routine is used both for fields and derivatives!
  complex(kind=real64),dimension(:),target   :: Ez    !EM fields or dervatives (horizontal/isotropic)
  complex(kind=real64),intent(in) :: omeps_recv !omega * epsilon in iReceiver layer
  complex(kind=real64),external   :: funcDz1TM
  integer(kind=int32),intent(in)  :: ilay       !layer index for derivatives - leave at zero for forward modeling
  complex(kind=real64),external,optional   :: funcDz1TMfwd  !function only needed for derivatives in iReceiver layer
  complex(kind=real64),external,optional   :: funcDz1TMv !functions for TM epsv deriv.
  complex(kind=real64),dimension(:),target,optional   :: Ezv !EM field for epsv derivatives

  !internal variables
  integer(kind=int32)   :: isrc    !source element counter
  real(kind=real64)     :: x,y,r   !temp source-iReceiver distances
  integer(kind=int32)   :: irec    !iReceiver counter
  real(kind=real64)     :: beta,betarot      !temp angles

  complex(kind=real64)          :: IDz1TM !interpolated integral values

  logical,dimension(nintHEDdvti)    :: wellbehaved  !indicates if Hankel integration can be used
  logical                       :: sz_eq_zr     !indicates if source and reveicer are at the same depth
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")

  real(kind=real64)     :: cosbeta,sinbeta      !cos(beta) and sin(beta), precompute for efficiency
  real(kind=real64)     :: cosbetarot,sinbetarot      !cos(betarot) and sin(betarot), precompute for efficiency
  complex(kind=real64)  :: cur      !temp source current
  complex(kind=real64)  :: JEh      !source current times constants
  integer(kind=int32)   :: idx      !source element index
  integer(kind=int32)   :: recidx   !iReceiver index
  !flag for computing epsv derivatives, not needed for forward computation, so "aniso" value cannot be used here
  logical               :: with_dvert
  complex(kind=real64),dimension(:),pointer   :: Erec  !points to Ez for isotropic and Ezv for VTI case

  !indicators for fast Hankel transform or adaptive integration
  wellbehaved = .true.
  sz_eq_zr = .false.
  if (sz.eq.zr) then
    wellbehaved(2) = .false.
    wellbehaved(5:7) = .false.
    sz_eq_zr = .true.
  endif
  
  if (present(funcDz1TMv)) then
    with_dvert = .true.
    if (sz_eq_zr) then
      wellbehaved(13) = .false.
      wellbehaved(15:16) = .false.
    endif
  else
    with_dvert = .false.
  endif


  do isrc = refl_var%isrcstart,refl_var%isrcend

    !get source current times constant factors
    !complex comjugate of currents was taken when reading in source currents, so no need to do that again here
    idx = refl_var%isrcperz(isrc)
    cur = conjg(src%cur(idx,ifreq))
    JEh = cur * sqrt(src%ljx(idx)**2 + src%ljy(idx)**2) / dfourpi

    do irec=refl_var%irecstart,refl_var%irecend
      recidx = refl_var%irecperzEz(irec)

      y = bgdat%Ezpos(recidx,2) - refl_var%ys(isrc)
      x = bgdat%Ezpos(recidx,1) - refl_var%xs(isrc)

      r = sqrt(x**2 + y**2)

      beta = atan2(y,x)
      cosbeta = cos(beta)
      sinbeta = sin(beta)
      betarot = beta - refl_var%betasrc(isrc)
      cosbetarot = cos(betarot)
      sinbetarot = sin(betarot)


      r_is_zero: if (r.eq.0._real64) then
        !quick & dirty: skip the point if iReceiver is right at source point
        if (sz_eq_zr) then
          if (refl_var%infolevel.ge.output_more) then
            write(*,'(a)') 'WARNING: cannot handle iReceiver right at source point yet!'
          endif
          cycle
        endif

      else !r is not zero
        smallr: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .true.
          call prepare_refcoef(refl_var,r,hed,aniso)

          !evaluate all integrals just for this radius
          ibesord = 1
          IDz1TM = compute_1val(funcDz1TM,r,sz,zr,ibesord,wellbehaved(5))

        else !r larger than threshold radius for spline interpolation
          !get integral values for this radius by interpolation from precomputed integral values
          IDz1TM = splinterp_1val(refl_var,5,r)
        endif smallr

        !Ez is zero for zero radius
        Ez(recidx) = Ez(recidx) + JEh * dci * cosbetarot * IDz1TM / omeps_recv
      endif r_is_zero

      dvert: if (with_dvert) then
        r_is_zerov: if (r.eq.0._real64) then
          !not necessary to check if source and iReceiver coincide exactly
          !  since in that case the above check already effects going to next iReceiver
          !if (sz_eq_zr) cycle

        else !r is not zero
          smallrv: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,hed,aniso)

            !evaluate all integrals just for this radius
            ibesord = 1
            IDz1TM = compute_1val(funcDz1TMv,r,sz,zr,ibesord,wellbehaved(16))

          else !r larger than threshold radius for spline interpolation
            !get integral values for this radius by interpolation from precomputed integral values
            IDz1TM = splinterp_1val(refl_var,16,r)
          endif smallrv

          !Ez is zero for zero radius
          Ezv(recidx) = Ezv(recidx) + JEh * dci * cosbetarot * IDz1TM / omeps_recv
        endif r_is_zerov
      endif dvert

    enddo !irec
  enddo  !source elements
  
  
  !special term is needed for derivative of Ez in iReceiver layer
  !isotropic: add this to deriv. for overall complex permittivity --> E
  !VTI: this term does not exist for epsh deriv, but exists for epsv deriv --> nothing added to E, add to Ev
  deriv_ilayrec: if (ilay .eq. ilayrec) then
  
    if (.not. present(funcDz1TMfwd)) then
      write(*,'(a)') 'ERROR: function for Ez in iReceiver layer not given, cannot compute correct Ez derivative!'
      return
    endif
  
    !there is a factor epsv in the term before the Ez integral
    !-> for isotropic case, epsv becomes eps and contribution is added to E
    !-> for VTI case, derivative of this term only exists for epsv -> add to Ev, nothing added to E
    if (with_dvert) then
      Erec => Ezv
    else
      Erec => Ez
    endif
  
    !indicators for fast Hankel transform or adaptive integration
    wellbehaved(11) = wellbehaved(5)  !assume that derivative and fwd integrals behave equally well or badly - not tested!

    do isrc = refl_var%isrcstart,refl_var%isrcend

      !get source current times constant factors
      !complex comjugate of currents was taken when reading in source currents, so no need to do that again here
      idx = refl_var%isrcperz(isrc)
      cur = conjg(src%cur(idx,ifreq))
      JEh = cur * sqrt(src%ljx(idx)**2 + src%ljy(idx)**2) / dfourpi

      do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzEz(irec)

        y = bgdat%Ezpos(recidx,2) - refl_var%ys(isrc)
        x = bgdat%Ezpos(recidx,1) - refl_var%xs(isrc)

        r = sqrt(x**2 + y**2)

        if (r.eq.0._real64) then
          !Ez is zero for zero radius
          cycle
        else

          beta = atan2(y,x)
          betarot = beta - refl_var%betasrc(isrc)
          cosbetarot = cos(betarot)

          if (r.lt.rsplmin) then
            !reflection coeff. for this radius
            call prepare_refcoef(refl_var,r,ved,aniso) !use "ved" to compute TM refl. coeff. only
            ibesord = 1
            IDz1TM = compute_1val(funcDz1TMfwd,r,sz,zr,ibesord,wellbehaved(11))
          else
            IDz1TM = splinterp_1val(refl_var,11,r)
          endif

          Erec(recidx) = Erec(recidx) - JEh * dci * cosbetarot * IDz1TM / (omeps_recv * epsv(ilayrec))
        endif
      enddo !irec
    enddo  !source elements
  endif deriv_ilayrec

endsubroutine interp_intvals_hed_Ez


!------------------------------------------------------------
!  1D EM subroutine interp_intvals_hed_Hxy
!
!  get interpolated field values at iReceiver locations, HED source,
!    Hx and / or Hy only
!
!  Rita Streich 2009-2011
!------------------------------------------------------------
subroutine interp_intvals_hed_Hxy(refl_var,src,ifreq,sz,zr,bgdat,Hx,Hy, &
  funcD0TE,funcD0TM,funcD1TE,funcD1TM,ilay, funcD0TMv,funcD1TMv,Hxv,Hyv)

  implicit none

  !external variables
  type(refl_struct)               :: refl_var   !everything needed throughout 1D computations
  type(sorec),intent(in)          :: src        !source definition (we need the currents here)
  integer(kind=int32),intent(in)  :: ifreq      !frequency index
  real(kind=real64),intent(in)    :: sz,zr      !source and iReceiver depth
  type(backgrounddata) :: bgdat      !coordinate vectors and output EM fields
    !input separate fields because the same routine is used both for fields and derivatives!
  complex(kind=real64),dimension(:)   :: Hx,Hy    !EM fields or dervatives (horizontal/isotropic)
  complex(kind=real64),external   :: funcD0TE,funcD0TM,funcD1TE,funcD1TM
  integer(kind=int32),intent(in)  :: ilay       !layer index for derivatives - leave at zero for forward modeling
  complex(kind=real64),external,optional   :: funcD0TMv,funcD1TMv !functions for TM epsv deriv.
  complex(kind=real64),dimension(:),optional   :: Hxv,Hyv !EM field for epsv derivatives

  !internal variables
  integer(kind=int32)   :: isrc    !source element counter
  real(kind=real64)     :: x,y,r   !temp source-iReceiver distances
  integer(kind=int32)   :: irec    !iReceiver counter
  real(kind=real64)     :: beta,betarot      !temp angles

  complex(kind=real64)          :: ID0TE,ID0TM,ID1TE,ID1TM !interpolated integral values

  logical,dimension(nintHEDdvti)    :: wellbehaved  !indicates if Hankel integration can be used
  logical                       :: sz_eq_zr     !indicates if source and reveicer are at the same depth
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")

  complex(kind=real64)  :: Hr,Hbeta             !temp field values in cylindrical coordinates
  real(kind=real64)     :: cosbeta,sinbeta      !cos(beta) and sin(beta), precompute for efficiency
  real(kind=real64)     :: cosbetarot,sinbetarot      !cos(betarot) and sin(betarot), precompute for efficiency
  complex(kind=real64)  :: cur      !temp source current
  complex(kind=real64)  :: JEh      !source current times constants
  integer(kind=int32)   :: idx      !source element index
  integer(kind=int32)   :: recidx   !iReceiver index
  !flag for computing epsv derivatives, not needed for forward computation, so "aniso" value cannot be used here
  logical               :: with_dvert

  !indicators for fast Hankel transform or adaptive integration
  wellbehaved = .true.
  sz_eq_zr = .false.
  if (sz.eq.zr) then
    wellbehaved(2) = .false.
    wellbehaved(5:7) = .false.
    sz_eq_zr = .true.
  endif
  
  if (present(funcD1TMv)) then
    with_dvert = .true.
    if (sz_eq_zr) then
      wellbehaved(13) = .false.
      wellbehaved(15:16) = .false.
    endif
  else
    with_dvert = .false.
  endif


  do isrc = refl_var%isrcstart,refl_var%isrcend

    !get source current times constant factors
    !complex comjugate of currents was taken when reading in source currents, so no need to do that again here
    idx = refl_var%isrcperz(isrc)
    cur = conjg(src%cur(idx,ifreq))
    JEh = cur * sqrt(src%ljx(idx)**2 + src%ljy(idx)**2) / dfourpi

    equal_Hxypos: if (bgdat%nHxy.gt.0) then

    do irec=refl_var%irecstart,refl_var%irecend
      recidx = refl_var%irecperzHxy(irec)

      y = bgdat%Hxypos(recidx,2) - refl_var%ys(isrc)
      x = bgdat%Hxypos(recidx,1) - refl_var%xs(isrc)

      r = sqrt(x**2 + y**2)

      beta = atan2(y,x)
      cosbeta = cos(beta)
      sinbeta = sin(beta)
      betarot = beta - refl_var%betasrc(isrc)
      cosbetarot = cos(betarot)
      sinbetarot = sin(betarot)


      r_is_zero: if (r.eq.0._real64) then
        !quick & dirty: skip the point if iReceiver is right at source point
        if (sz_eq_zr) then
          if (refl_var%infolevel.ge.output_more) then
            write(*,'(a)') 'WARNING: cannot handle iReceiver right at source point yet!'
          endif
          cycle
        endif

        !Bessel function J0(r=0) = 1, J1(r=0) = 0
        !--> need to evaluate zero-order integrals only
        !use adaptive integration
        ID0TE = compute_1valr0(funcD0TE)
        ID0TM = compute_1valr0(funcD0TM)

        !limits from (1/R) * integral(bessel1) become bessel0 ...
        Hr = - JEh * sinbetarot * (0.5_real64 * (ID0TE + ID0TM))
        Hbeta = - JEh * cosbetarot * (0.5_real64 * (ID0TE + ID0TM))
        
      else !r is not zero
        smallr: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .true.
          call prepare_refcoef(refl_var,r,hed,aniso)

          !evaluate all integrals just for this radius
          ibesord = 0
          ID0TE = compute_1val(funcD0TE,r,sz,zr,ibesord,wellbehaved(6))
          ID0TM = compute_1val(funcD0TM,r,sz,zr,ibesord,wellbehaved(7))
          ibesord = 1
          ID1TE = compute_1val(funcD1TE,r,sz,zr,ibesord,wellbehaved(8))
          ID1TM = compute_1val(funcD1TM,r,sz,zr,ibesord,wellbehaved(9))
        else !r larger than threshold radius for spline interpolation

          !get integral values for this radius by interpolation from precomputed integral values
          ID0TE = splinterp_1val(refl_var,6,r)
          ID0TM = splinterp_1val(refl_var,7,r)
          ID1TE = splinterp_1val(refl_var,8,r)
          ID1TM = splinterp_1val(refl_var,9,r)
        endif smallr

        Hr = JEh * sinbetarot * (-ID0TE + (ID1TE - ID1TM)/r)
        Hbeta = -JEh * cosbetarot * (ID0TM + (ID1TE - ID1TM)/r)
      endif r_is_zero

      !field values rotated back to original x-y coordinate system
      !don't take complex conjugate here, but at the very end!
      Hx(recidx) = Hx(recidx) + cosbeta*Hr - sinbeta*Hbeta
      Hy(recidx) = Hy(recidx) + sinbeta*Hr + cosbeta*Hbeta

      dvert: if (with_dvert) then
        r_is_zerov: if (r.eq.0._real64) then
          !not necessary to check if source and iReceiver coincide exactly
          !  since in that case the above check already effects going to next iReceiver
          !if (sz_eq_zr) cycle

          !Bessel function J0(r=0) = 1, J1(r=0) = 0
          !--> need to evaluate zero-order integrals only
          !use adaptive integration
          ID0TM = compute_1valr0(funcD0TMv)

          !limits from (1/R) * integral(bessel1) become bessel0 ...
          Hr = - JEh * sinbetarot * 0.5_real64 * ID0TM
          Hbeta = - JEh * cosbetarot * 0.5_real64 * ID0TM
        
        else !r is not zero
          smallrv: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,hed,aniso)

            !evaluate all integrals just for this radius
            ibesord = 0
            ID0TM = compute_1val(funcD0TMv,r,sz,zr,ibesord,wellbehaved(15))
            ibesord = 1
            ID1TM = compute_1val(funcD1TMv,r,sz,zr,ibesord,wellbehaved(14))

          else !r larger than threshold radius for spline interpolation

            !get integral values for this radius by interpolation from precomputed integral values
            ID0TM = splinterp_1val(refl_var,15,r)
            ID1TM = splinterp_1val(refl_var,14,r)
          endif smallrv

          Hr = JEh * sinbetarot * (- ID1TM/r)
          Hbeta = -JEh * cosbetarot * (ID0TM - ID1TM/r)
        endif r_is_zerov

        !field values rotated back to original x-y coordinate system
        !don't take complex conjugate here, but at the very end!
        Hxv(recidx) = Hxv(recidx) + cosbeta*Hr - sinbeta*Hbeta
        Hyv(recidx) = Hyv(recidx) + sinbeta*Hr + cosbeta*Hbeta
      endif dvert
    enddo !irec

    else
      have_Hx: if (bgdat%nHx.gt.0) then

    do irec=refl_var%irecstart,refl_var%irecend
      recidx = refl_var%irecperzHxy(irec)

      y = bgdat%Hxpos(recidx,2) - refl_var%ys(isrc)
      x = bgdat%Hxpos(recidx,1) - refl_var%xs(isrc)

      r = sqrt(x**2 + y**2)

      beta = atan2(y,x)
      cosbeta = cos(beta)
      sinbeta = sin(beta)
      betarot = beta - refl_var%betasrc(isrc)
      cosbetarot = cos(betarot)
      sinbetarot = sin(betarot)


      r_is_zerohx: if (r.eq.0._real64) then
        !quick & dirty: skip the point if iReceiver is right at source point
        if (sz_eq_zr) then
          if (refl_var%infolevel.ge.output_more) then
            write(*,'(a)') 'WARNING: cannot handle iReceiver right at source point yet!'
          endif
          cycle
        endif

        !Bessel function J0(r=0) = 1, J1(r=0) = 0
        !--> need to evaluate zero-order integrals only
        !use adaptive integration
        ID0TE = compute_1valr0(funcD0TE)
        ID0TM = compute_1valr0(funcD0TM)

        !limits from (1/R) * integral(bessel1) become bessel0 ...
        Hr = - JEh * sinbetarot * (0.5_real64 * (ID0TE + ID0TM))
        Hbeta = - JEh * cosbetarot * (0.5_real64 * (ID0TE + ID0TM))
        
      else !r is not zero
        smallrhx: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .true.
          call prepare_refcoef(refl_var,r,hed,aniso)

          !evaluate all integrals just for this radius
          ibesord = 0
          ID0TE = compute_1val(funcD0TE,r,sz,zr,ibesord,wellbehaved(6))
          ID0TM = compute_1val(funcD0TM,r,sz,zr,ibesord,wellbehaved(7))
          ibesord = 1
          ID1TE = compute_1val(funcD1TE,r,sz,zr,ibesord,wellbehaved(8))
          ID1TM = compute_1val(funcD1TM,r,sz,zr,ibesord,wellbehaved(9))
        else !r larger than threshold radius for spline interpolation

          !get integral values for this radius by interpolation from precomputed integral values
          ID0TE = splinterp_1val(refl_var,6,r)
          ID0TM = splinterp_1val(refl_var,7,r)
          ID1TE = splinterp_1val(refl_var,8,r)
          ID1TM = splinterp_1val(refl_var,9,r)
        endif smallrhx

        Hr = JEh * sinbetarot * (-ID0TE + (ID1TE - ID1TM)/r)
        Hbeta = -JEh * cosbetarot * (ID0TM + (ID1TE - ID1TM)/r)
      endif r_is_zerohx

      !field values rotated back to original x-y coordinate system
      !don't take complex conjugate here, but at the very end!
      Hx(recidx) = Hx(recidx) + cosbeta*Hr - sinbeta*Hbeta

      dverthx: if (with_dvert) then
        r_is_zerovhx: if (r.eq.0._real64) then
          !not necessary to check if source and iReceiver coincide exactly
          !  since in that case the above check already effects going to next iReceiver
          !if (sz_eq_zr) cycle

          !Bessel function J0(r=0) = 1, J1(r=0) = 0
          !--> need to evaluate zero-order integrals only
          !use adaptive integration
          ID0TM = compute_1valr0(funcD0TMv)

          !limits from (1/R) * integral(bessel1) become bessel0 ...
          Hr = - JEh * sinbetarot * 0.5_real64 * ID0TM
          Hbeta = - JEh * cosbetarot * 0.5_real64 * ID0TM
        
        else !r is not zero
          smallrvhx: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,hed,aniso)

            !evaluate all integrals just for this radius
            ibesord = 0
            ID0TM = compute_1val(funcD0TMv,r,sz,zr,ibesord,wellbehaved(15))
            ibesord = 1
            ID1TM = compute_1val(funcD1TMv,r,sz,zr,ibesord,wellbehaved(14))

          else !r larger than threshold radius for spline interpolation

            !get integral values for this radius by interpolation from precomputed integral values
            ID0TM = splinterp_1val(refl_var,15,r)
            ID1TM = splinterp_1val(refl_var,14,r)
          endif smallrvhx

          Hr = JEh * sinbetarot * (- ID1TM/r)
          Hbeta = -JEh * cosbetarot * (ID0TM - ID1TM/r)
        endif r_is_zerovhx

        !field values rotated back to original x-y coordinate system
        !don't take complex conjugate here, but at the very end!
        Hxv(recidx) = Hxv(recidx) + cosbeta*Hr - sinbeta*Hbeta
      endif dverthx
    enddo !irec

      endif have_Hx

      have_Hy: if (bgdat%nHy.gt.0) then

    do irec=refl_var%irecstart,refl_var%irecend
      recidx = refl_var%irecperzHxy(irec)

      y = bgdat%Hypos(recidx,2) - refl_var%ys(isrc)
      x = bgdat%Hypos(recidx,1) - refl_var%xs(isrc)

      r = sqrt(x**2 + y**2)

      beta = atan2(y,x)
      cosbeta = cos(beta)
      sinbeta = sin(beta)
      betarot = beta - refl_var%betasrc(isrc)
      cosbetarot = cos(betarot)
      sinbetarot = sin(betarot)

      r_is_zerohy: if (r.eq.0._real64) then
        !quick & dirty: skip the point if iReceiver is right at source point
        if (sz_eq_zr) then
          if (refl_var%infolevel.ge.output_more) then
            write(*,'(a)') 'WARNING: cannot handle iReceiver right at source point yet!'
          endif
          cycle
        endif

        !Bessel function J0(r=0) = 1, J1(r=0) = 0
        !--> need to evaluate zero-order integrals only
        !use adaptive integration
        ID0TE = compute_1valr0(funcD0TE)
        ID0TM = compute_1valr0(funcD0TM)

        !limits from (1/R) * integral(bessel1) become bessel0 ...
        Hr = - JEh * sinbetarot * (0.5_real64 * (ID0TE + ID0TM))
        Hbeta = - JEh * cosbetarot * (0.5_real64 * (ID0TE + ID0TM))
        
      else !r is not zero
        smallrhy: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .true.
          call prepare_refcoef(refl_var,r,hed,aniso)

          !evaluate all integrals just for this radius
          ibesord = 0
          ID0TE = compute_1val(funcD0TE,r,sz,zr,ibesord,wellbehaved(6))
          ID0TM = compute_1val(funcD0TM,r,sz,zr,ibesord,wellbehaved(7))
          ibesord = 1
          ID1TE = compute_1val(funcD1TE,r,sz,zr,ibesord,wellbehaved(8))
          ID1TM = compute_1val(funcD1TM,r,sz,zr,ibesord,wellbehaved(9))
        else !r larger than threshold radius for spline interpolation

          !get integral values for this radius by interpolation from precomputed integral values
          ID0TE = splinterp_1val(refl_var,6,r)
          ID0TM = splinterp_1val(refl_var,7,r)
          ID1TE = splinterp_1val(refl_var,8,r)
          ID1TM = splinterp_1val(refl_var,9,r)
        endif smallrhy

        Hr = JEh * sinbetarot * (-ID0TE + (ID1TE - ID1TM)/r)
        Hbeta = -JEh * cosbetarot * (ID0TM + (ID1TE - ID1TM)/r)
      endif r_is_zerohy

      !field values rotated back to original x-y coordinate system
      !don't take complex conjugate here, but at the very end!
      Hy(recidx) = Hy(recidx) + sinbeta*Hr + cosbeta*Hbeta

      dverthy: if (with_dvert) then
        r_is_zerovhy: if (r.eq.0._real64) then
          !not necessary to check if source and iReceiver coincide exactly
          !  since in that case the above check already effects going to next iReceiver
          !if (sz_eq_zr) cycle

          !Bessel function J0(r=0) = 1, J1(r=0) = 0
          !--> need to evaluate zero-order integrals only
          !use adaptive integration
          ID0TM = compute_1valr0(funcD0TMv)

          !limits from (1/R) * integral(bessel1) become bessel0 ...
          Hr = - JEh * sinbetarot * 0.5_real64 * ID0TM
          Hbeta = - JEh * cosbetarot * 0.5_real64 * ID0TM
        
        else !r is not zero
          smallrvhy: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,hed,aniso)

            !evaluate all integrals just for this radius
            ibesord = 0
            ID0TM = compute_1val(funcD0TMv,r,sz,zr,ibesord,wellbehaved(15))
            ibesord = 1
            ID1TM = compute_1val(funcD1TMv,r,sz,zr,ibesord,wellbehaved(14))

          else !r larger than threshold radius for spline interpolation

            !get integral values for this radius by interpolation from precomputed integral values
            ID0TM = splinterp_1val(refl_var,15,r)
            ID1TM = splinterp_1val(refl_var,14,r)
          endif smallrvhy

          Hr = JEh * sinbetarot * (- ID1TM/r)
          Hbeta = -JEh * cosbetarot * (ID0TM - ID1TM/r)
        endif r_is_zerovhy

        !field values rotated back to original x-y coordinate system
        !don't take complex conjugate here, but at the very end!
        Hyv(recidx) = Hyv(recidx) + sinbeta*Hr + cosbeta*Hbeta
      endif dverthy
    enddo !irec

      endif have_Hy
    endif equal_Hxypos

  enddo  !source elements
  
endsubroutine interp_intvals_hed_Hxy


!------------------------------------------------------------
!  1D EM subroutine interp_intvals_hed_Hz
!
!  get interpolated field values at iReceiver locations, HED source,
!    Hz only
!
!  Rita Streich 2009-2011
!------------------------------------------------------------
subroutine interp_intvals_hed_Hz(refl_var,src,ifreq,sz,zr,bgdat,Hz,ommu, funcAz1TE,ilay)

  implicit none

  !external variables
  type(refl_struct)               :: refl_var   !everything needed throughout 1D computations
  type(sorec),intent(in)          :: src        !source definition (we need the currents here)
  integer(kind=int32),intent(in)  :: ifreq      !frequency index
  real(kind=real64),intent(in)    :: sz,zr      !source and iReceiver depth
  type(backgrounddata) :: bgdat      !coordinate vectors and output EM fields
    !input separate fields because the same routine is used both for fields and derivatives!
  complex(kind=real64),dimension(:)   :: Hz    !EM fields or dervatives (horizontal/isotropic)
  real(kind=real64),intent(in)    :: ommu       !omega * mu0
  complex(kind=real64),external   :: funcAz1TE
  integer(kind=int32),intent(in)  :: ilay       !layer index for derivatives - leave at zero for forward modeling

  !internal variables
  integer(kind=int32)   :: isrc    !source element counter
  real(kind=real64)     :: x,y,r   !temp source-iReceiver distances
  integer(kind=int32)   :: irec    !iReceiver counter
  real(kind=real64)     :: beta,betarot      !temp angles

  complex(kind=real64)          :: IAz1TE !interpolated integral values

  logical,dimension(nintHEDdvti)    :: wellbehaved  !indicates if Hankel integration can be used
  logical                       :: sz_eq_zr     !indicates if source and reveicer are at the same depth
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")

  real(kind=real64)     :: sinbetarot      !cos(betarot) and sin(betarot), precompute for efficiency
  complex(kind=real64)  :: cur      !temp source current
  complex(kind=real64)  :: JEh      !source current times constants
  integer(kind=int32)   :: idx      !source element index
  integer(kind=int32)   :: recidx   !iReceiver index

  !indicators for fast Hankel transform or adaptive integration
  wellbehaved = .true. !only need element (10), no deriv. with respect to vertical medium properties
  sz_eq_zr = .false.
  if (sz.eq.zr) then
    sz_eq_zr = .true.
  endif

  do isrc = refl_var%isrcstart,refl_var%isrcend

    !get source current times constant factors
    !complex comjugate of currents was taken when reading in source currents, so no need to do that again here
    idx = refl_var%isrcperz(isrc)
    cur = conjg(src%cur(idx,ifreq))
    JEh = cur * sqrt(src%ljx(idx)**2 + src%ljy(idx)**2) / dfourpi

    do irec=refl_var%irecstart,refl_var%irecend
      recidx = refl_var%irecperzHz(irec)

      y = bgdat%Hzpos(recidx,2) - refl_var%ys(isrc)
      x = bgdat%Hzpos(recidx,1) - refl_var%xs(isrc)

      r = sqrt(x**2 + y**2)

      beta = atan2(y,x)
      betarot = beta - refl_var%betasrc(isrc)
      sinbetarot = sin(betarot)


      r_is_zero: if (r.eq.0._real64) then
        !quick & dirty: skip the point if iReceiver is right at source point
        if (sz_eq_zr) then
          if (refl_var%infolevel.ge.output_more) then
            write(*,'(a)') 'WARNING: cannot handle iReceiver right at source point yet!'
          endif
          cycle
        endif

      else !r is not zero
        smallr: if (r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .true.
          call prepare_refcoef(refl_var,r,hed,aniso)

          !evaluate all integrals just for this radius
          ibesord = 1
          IAz1TE = compute_1val(funcAz1TE,r,sz,zr,ibesord,wellbehaved(10))
        else !r larger than threshold radius for spline interpolation

          !get integral values for this radius by interpolation from precomputed integral values
          IAz1TE = splinterp_1val(refl_var,10,r)
        endif smallr

        !Hz in original x-y coordinate system - no special case needed for radius zero
        Hz(recidx) = Hz(recidx) + JEh * dci * sinbetarot * IAz1TE / ommu
      endif r_is_zero

    enddo !irec
  enddo  !source elements
  
endsubroutine interp_intvals_hed_Hz



