!------------------------------------------------------------
!  1D EM subroutine interp_intvals_hmd_allcomp
!
!  get interpolated field values at receiver locations, HMD source,
!  all field components at the same locations
!
!  Rita Streich 2009-2011
!------------------------------------------------------------
subroutine interp_intvals_hmd_allcomp(refl_var,src,ifreq,sz,zr,bgdat,Ex,Ey,Ez,Hx,Hy,Hz,omeps_recv,j_om_mu, &
  funcB0TE,funcB0TM,funcB1TE,funcB1TM,funcCz1TM,funcC0TE,funcC0TM,funcC1TE,funcC1TM,funcBz1TE,ilay,funcCz1TMfwd, &
  funcB0TMv,funcB1TMv,funcCz1TMv,funcC0TMv,funcC1TMv,Exv,Eyv,Ezv,Hxv,Hyv,Hzv)

  implicit none

  !external variables
  type(refl_struct)               :: refl_var   !everything needed throughout 1D computations
  type(sorec),intent(in)          :: src        !source definition (we need the currents here)
  integer(kind=int32),intent(in)  :: ifreq      !frequency index
  real(kind=real64),intent(in)    :: sz,zr      !source and receiver depth
  type(backgrounddata) :: bgdat      !coordinate vectors and output EM fields
    !input separate fields because the same routine is used both for fields and derivatives!
  complex(kind=real64),dimension(:)   :: Ex,Ey,Hx,Hy,Hz  !electric and magnetic fields: nr of receivers
  complex(kind=real64),dimension(:),target   :: Ez  !electric fields: nr of receivers
  complex(kind=real64),intent(in) :: omeps_recv  !omega * epsilon in receiver layer
  complex(kind=real64),intent(in) :: j_om_mu    !j * omega * mu0
  complex(kind=real64),external   :: funcB0TE,funcB0TM,funcB1TE,funcB1TM,funcCz1TM,funcC0TE,funcC0TM,funcC1TE,funcC1TM,funcBz1TE
  integer(kind=int32),intent(in)  :: ilay       !layer index for derivatives - leave at zero for forward modeling
  complex(kind=real64),external,optional   :: funcCz1TMfwd  !function only needed for derivatives in receiver layer
  complex(kind=real64),external,optional   :: funcB0TMv,funcB1TMv,funcCz1TMv,funcC0TMv,funcC1TMv !derivative integrals for epsv
  complex(kind=real64),dimension(:),optional   :: Exv,Eyv,Hxv,Hyv,Hzv  !electric and magnetic fields: nr of receivers
  complex(kind=real64),dimension(:),target,optional   :: Ezv  !electric fields: nr of receivers

  !internal variables
  integer(kind=int32)   :: isrc    !source element counter
  real(kind=real64)     :: x,y,r   !temp source-receiver distances
  integer(kind=int32)   :: irec    !receiver counter
  real(kind=real64)     :: beta,betarot      !temp angles

  complex(kind=real64)          :: IB0TE,IB0TM,IB1TE,IB1TM,IBz1TE !interpolated integral values
  complex(kind=real64)          :: IC0TE,IC0TM,IC1TE,IC1TM,ICz1TM

  logical,dimension(nintHMDdvti)    :: wellbehaved  !indicates if Hankel integration can be used
  logical                       :: sz_eq_zr     !indicates if source and reveicer are at the same depth
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")

  complex(kind=real64)  :: Er,Ebeta,Hr,Hbeta    !temp field values in cylindrical coordinates
  real(kind=real64)     :: cosbeta,sinbeta      !cos(beta) and sin(beta), precompute for efficiency
  real(kind=real64)     :: cosbetarot,sinbetarot      !cos(betarot) and sin(betarot), precompute for efficiency
  complex(kind=real64)  :: cur      !temp source current
  complex(kind=real64)  :: JMh      !source current times constants
  integer(kind=int32)   :: idx      !source element index
  integer(kind=int32)   :: recidx   !receiver index
  !flag for computing epsv derivatives, not needed for forward computation, so "aniso" value cannot be used here
  logical               :: with_dvert
  complex(kind=real64),dimension(:),pointer   :: Erec  !points to either Ez (isotropic) or Ezv (VTI)

  !indicators for fast Hankel transform or adaptive integration
  wellbehaved = .true.
  sz_eq_zr = .false.
  if(sz.eq.zr) then
    wellbehaved(1:2) = .false.
    wellbehaved(6) = .false.
    wellbehaved(10) = .false.
    sz_eq_zr = .true.
  endif

  if(present(funcC1TMv)) then
    with_dvert = .true.
    if(sz_eq_zr) then
      wellbehaved(13) = .false.
    endif
  else
    with_dvert = .false.
  endif


  do isrc = refl_var%isrcstart,refl_var%isrcend

    !get source current times constant factors
    idx = refl_var%isrcperz(isrc)
    cur = conjg(src%cur(idx,ifreq))
    JMh = - cur * sqrt(src%akx(idx)**2 + src%aky(idx)**2) * j_om_mu / dfourpi

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

 
      r_is_zero: if(r.eq.0._real64) then

        !quick & dirty: skip the point if receiver is right at source point
        if(sz_eq_zr) then
          if(refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle receiver right at source point yet!'
          !skip numerical integral evaluations - they would most likely fail, need to look at this...
          cycle
        endif

        !Bessel function J0(r=0) = 1, J1(r=0) = 0
        !--> need to evaluate zero-order integrals only
        !use adaptive integration
        IB0TE = compute_1valr0(funcB0TE)
        IB0TM = compute_1valr0(funcB0TM)
        IC0TE = compute_1valr0(funcC0TE)
        IC0TM = compute_1valr0(funcC0TM)

        Er = - JMh * sinbetarot * (0.5_real64 * (IB0TE + IB0TM))
        Ebeta = - JMh * cosbetarot * (0.5_real64 * (IB0TE + IB0TM))
        Hr = JMh * cosbetarot * (0.5_real64 * (IC0TE + IC0TM))
        Hbeta = - JMh * sinbetarot * (0.5_real64 * (IC0TE + IC0TM))

      else !r is not zero
        smallr: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .true.
          call prepare_refcoef(refl_var,r,hmd,aniso)

          !evaluate all integrals just for this radius
          ibesord = 0
          IB0TE = compute_1val(funcB0TE,r,sz,zr,ibesord,wellbehaved(1))
          IB0TM = compute_1val(funcB0TM,r,sz,zr,ibesord,wellbehaved(2))
          ibesord = 1
          IB1TE = compute_1val(funcB1TE,r,sz,zr,ibesord,wellbehaved(3))
          IB1TM = compute_1val(funcB1TM,r,sz,zr,ibesord,wellbehaved(4))

          ibesord = 1
          ICz1TM = compute_1val(funcCz1TM,r,sz,zr,ibesord,wellbehaved(5))

          ibesord = 0
          IC0TE = compute_1val(funcC0TE,r,sz,zr,ibesord,wellbehaved(6))
          IC0TM = compute_1val(funcC0TM,r,sz,zr,ibesord,wellbehaved(7))
          ibesord = 1
          IC1TE = compute_1val(funcC1TE,r,sz,zr,ibesord,wellbehaved(8))
          IC1TM = compute_1val(funcC1TM,r,sz,zr,ibesord,wellbehaved(9))

          ibesord = 1
          IBz1TE = compute_1val(funcBz1TE,r,sz,zr,ibesord,wellbehaved(10))

        else !r larger than threshold radius for spline interpolation

          !get integral values for this radius by interpolation from precomputed integral values
          IB0TE = splinterp_1val(refl_var,1,r)
          IB0TM = splinterp_1val(refl_var,2,r)
          IB1TE = splinterp_1val(refl_var,3,r)
          IB1TM = splinterp_1val(refl_var,4,r)
          ICz1TM = splinterp_1val(refl_var,5,r)

          IC0TE = splinterp_1val(refl_var,6,r)
          IC0TM = splinterp_1val(refl_var,7,r)
          IC1TE = splinterp_1val(refl_var,8,r)
          IC1TM = splinterp_1val(refl_var,9,r)
          IBz1TE = splinterp_1val(refl_var,10,r)

        endif smallr

        !Er and Ebeta relative to horizontal source orientation (factor 4pi and dipole length absorbed into JEh)
        Er = - JMh * sinbetarot * (IB0TM + (IB1TE - IB1TM)/r)
        Ebeta = JMh * cosbetarot * (-IB0TE + (IB1TE - IB1TM)/r)

        !Ez in original x-y coordinate system, Ez is zero for r=0
        Ez(recidx) = Ez(recidx) - (JMh/(dci*omeps_recv)) * sinbetarot * ICz1TM

        Hr = - JMh * cosbetarot * (-IC0TE + (IC1TE - IC1TM)/r)
        Hbeta = - JMh * sinbetarot * (IC0TM + (IC1TE - IC1TM)/r)

        !Hz in original x-y coordinate system
        Hz(recidx) = Hz(recidx) + (JMh / j_om_mu) * cosbetarot * IBz1TE

      endif r_is_zero

      !Ex, Ey in original x-y coordinate system
      Ex(recidx) = Ex(recidx) + cosbeta*Er - sinbeta*Ebeta
      Ey(recidx) = Ey(recidx) + sinbeta*Er + cosbeta*Ebeta

      !Hx, Hy in original x-y coordinate system
      Hx(recidx) = Hx(recidx) + cosbeta*Hr - sinbeta*Hbeta
      Hy(recidx) = Hy(recidx) + sinbeta*Hr + cosbeta*Hbeta


      !VTI: deriv. for epsv
      dvert: if(with_dvert) then
        r_is_zerov: if(r.eq.0._real64) then

          !already cycled if receiver is exactly at source point
          !if(sz_eq_zr) cycle

          !Bessel function J0(r=0) = 1, J1(r=0) = 0
          !--> need to evaluate zero-order integrals only
          !use adaptive integration
          IB0TM = compute_1valr0(funcB0TMv)
          IC0TM = compute_1valr0(funcC0TMv)

          Er = - JMh * sinbetarot * 0.5_real64 * IB0TM
          Ebeta = - JMh * cosbetarot * 0.5_real64 * IB0TM
          Hr = JMh * cosbetarot * 0.5_real64 * IC0TM
          Hbeta = - JMh * sinbetarot * 0.5_real64 * IC0TM

        else !r is not zero
          smallrv: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,hmd,aniso)

            !evaluate all integrals just for this radius
            ibesord = 0
            IB0TM = compute_1val(funcB0TMv,r,sz,zr,ibesord,wellbehaved(13))
            ibesord = 1
            IB1TM = compute_1val(funcB1TMv,r,sz,zr,ibesord,wellbehaved(12))

            ibesord = 1
            ICz1TM = compute_1val(funcCz1TMv,r,sz,zr,ibesord,wellbehaved(16))

            ibesord = 0
            IC0TM = compute_1val(funcC0TMv,r,sz,zr,ibesord,wellbehaved(15))
            ibesord = 1
            IC1TM = compute_1val(funcC1TMv,r,sz,zr,ibesord,wellbehaved(14))

          else !r larger than threshold radius for spline interpolation

            !get integral values for this radius by interpolation from precomputed integral values
            IB0TM = splinterp_1val(refl_var,13,r)
            IB1TM = splinterp_1val(refl_var,12,r)
            ICz1TM = splinterp_1val(refl_var,16,r)

            IC0TM = splinterp_1val(refl_var,15,r)
            IC1TM = splinterp_1val(refl_var,14,r)
          endif smallrv

          !Er and Ebeta relative to horizontal source orientation (factor 4pi and dipole length absorbed into JEh)
          Er = - JMh * sinbetarot * (IB0TM - IB1TM/r)
          Ebeta = JMh * cosbetarot * (- IB1TM/r)

          !Ez in original x-y coordinate system, Ez is zero for r=0
          Ezv(recidx) = Ezv(recidx) - (JMh/(dci*omeps_recv)) * sinbetarot * ICz1TM

          Hr = - JMh * cosbetarot * (- IC1TM/r)
          Hbeta = - JMh * sinbetarot * (IC0TM - IC1TM/r)

          !deriv. dHz / dpesv is zero
        endif r_is_zerov

        !Ex, Ey in original x-y coordinate system
        Exv(recidx) = Exv(recidx) + cosbeta*Er - sinbeta*Ebeta
        Eyv(recidx) = Eyv(recidx) + sinbeta*Er + cosbeta*Ebeta

        !Hx, Hy in original x-y coordinate system
        Hxv(recidx) = Hxv(recidx) + cosbeta*Hr - sinbeta*Hbeta
        Hyv(recidx) = Hyv(recidx) + sinbeta*Hr + cosbeta*Hbeta
      endif dvert

    enddo !irec
  enddo  !source elements


  !special term for Ez for derivatives in receiver layer
  deriv_ilayrec: if(ilay .eq. ilayrec) then
  
    if(.not. present(funcCz1TMfwd)) then
      write(*,'(a)') 'ERROR: function for Ez in receiver layer not given, cannot compute correct Ez derivative!'
      return
    endif

    !there is a factor epsv in the term before the Ez integral
    !-> for isotropic case, epsv becomes eps and contribution is added to E
    !-> for VTI case, derivative of this term only exists for epsv -> add to Ev, nothing added to E
    if(with_dvert) then
      Erec => Ezv
    else
      Erec => Ez
    endif
  
    wellbehaved(11) = wellbehaved(5)

    do isrc = refl_var%isrcstart,refl_var%isrcend

      !get source current times constant factors
      idx = refl_var%isrcperz(isrc)
      cur = conjg(src%cur(idx,ifreq))
      JMh = - cur * sqrt(src%akx(idx)**2 + src%aky(idx)**2) * j_om_mu / dfourpi

      do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzExy(irec)

        y = bgdat%Exypos(recidx,2) - refl_var%ys(isrc)
        x = bgdat%Exypos(recidx,1) - refl_var%xs(isrc)
        r = sqrt(x**2 + y**2)

        !Ez is zero for r=0
        if(r.eq.0._real64) then
          cycle
        else
          beta = atan2(y,x)
          betarot = beta - refl_var%betasrc(isrc)
          sinbetarot = sin(betarot)

          if(r.lt.rsplmin) then
            !reflection coeff. for this radius
            call prepare_refcoef(refl_var,r,ved,aniso) !use "ved" to compute TM refl. coeff. only
            ibesord = 1
            ICz1TM = compute_1val(funcCz1TMfwd,r,sz,zr,ibesord,wellbehaved(11))
          else
            ICz1TM = splinterp_1val(refl_var,11,r)
          endif

          !Ez in original x-y coordinate system
          Erec(recidx) = Erec(recidx) + (JMh/(dci*omeps_recv*epsv(ilayrec))) * sinbetarot * ICz1TM
        endif

      enddo !irec
    enddo  !source elements

  endif deriv_ilayrec

endsubroutine interp_intvals_hmd_allcomp


!------------------------------------------------------------
!  1D EM subroutine interp_intvals_hmd_Exy
!
!  get interpolated field values at receiver locations, HMD source,
!  Ex and / or Ey only
!
!  Rita Streich 2009-2011
!------------------------------------------------------------
subroutine interp_intvals_hmd_Exy(refl_var,src,ifreq,sz,zr,bgdat,Ex,Ey,j_om_mu, funcB0TE,funcB0TM,funcB1TE,funcB1TM,ilay, &
  funcB0TMv,funcB1TMv,Exv,Eyv)

  implicit none

  !external variables
  type(refl_struct)               :: refl_var   !everything needed throughout 1D computations
  type(sorec),intent(in)          :: src        !source definition (we need the currents here)
  integer(kind=int32),intent(in)  :: ifreq      !frequency index
  real(kind=real64),intent(in)    :: sz,zr      !source and receiver depth
  type(backgrounddata) :: bgdat      !coordinate vectors and output EM fields
    !input separate fields because the same routine is used both for fields and derivatives!
  complex(kind=real64),dimension(:)   :: Ex,Ey  !electric and magnetic fields: nr of receivers
  complex(kind=real64),intent(in) :: j_om_mu    !j * omega * mu0
  complex(kind=real64),external   :: funcB0TE,funcB0TM,funcB1TE,funcB1TM
  integer(kind=int32),intent(in)  :: ilay       !layer index for derivatives - leave at zero for forward modeling
  complex(kind=real64),external,optional   :: funcB0TMv,funcB1TMv !derivative integrals for epsv
  complex(kind=real64),dimension(:),optional   :: Exv,Eyv  !electric and magnetic fields: nr of receivers

  !internal variables
  integer(kind=int32)   :: isrc    !source element counter
  real(kind=real64)     :: x,y,r   !temp source-receiver distances
  integer(kind=int32)   :: irec    !receiver counter
  real(kind=real64)     :: beta,betarot      !temp angles

  complex(kind=real64)          :: IB0TE,IB0TM,IB1TE,IB1TM !interpolated integral values

  logical,dimension(nintHMDdvti)    :: wellbehaved  !indicates if Hankel integration can be used
  logical                       :: sz_eq_zr     !indicates if source and reveicer are at the same depth
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")

  complex(kind=real64)  :: Er,Ebeta             !temp field values in cylindrical coordinates
  real(kind=real64)     :: cosbeta,sinbeta      !cos(beta) and sin(beta), precompute for efficiency
  real(kind=real64)     :: cosbetarot,sinbetarot      !cos(betarot) and sin(betarot), precompute for efficiency
  complex(kind=real64)  :: cur      !temp source current
  complex(kind=real64)  :: JMh      !source current times constants
  integer(kind=int32)   :: idx      !source element index
  integer(kind=int32)   :: recidx   !receiver index
  !flag for computing epsv derivatives, not needed for forward computation, so "aniso" value cannot be used here
  logical               :: with_dvert

  !indicators for fast Hankel transform or adaptive integration
  wellbehaved = .true.
  sz_eq_zr = .false.
  if(sz.eq.zr) then
    wellbehaved(1:2) = .false.
    wellbehaved(6) = .false.
    wellbehaved(10) = .false.
    sz_eq_zr = .true.
  endif

  if(present(funcB1TMv)) then
    with_dvert = .true.
    if(sz_eq_zr) then
      wellbehaved(13) = .false.
    endif
  else
    with_dvert = .false.
  endif

  do isrc = refl_var%isrcstart,refl_var%isrcend

    !get source current times constant factors
    idx = refl_var%isrcperz(isrc)
    cur = conjg(src%cur(idx,ifreq))
    JMh = - cur * sqrt(src%akx(idx)**2 + src%aky(idx)**2) * j_om_mu / dfourpi


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
      betarot = beta - refl_var%betasrc(isrc)
      cosbetarot = cos(betarot)
      sinbetarot = sin(betarot)

       r_is_zero: if(r.eq.0._real64) then

        !quick & dirty: skip the point if receiver is right at source point
        if(sz_eq_zr) then
          if(refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle receiver right at source point yet!'
          !skip numerical integral evaluations - they would most likely fail, need to look at this...
          cycle
        endif

        !Bessel function J0(r=0) = 1, J1(r=0) = 0
        !--> need to evaluate zero-order integrals only
        !use adaptive integration
        IB0TE = compute_1valr0(funcB0TE)
        IB0TM = compute_1valr0(funcB0TM)

        Er = - JMh * sinbetarot * (0.5_real64 * (IB0TE + IB0TM))
        Ebeta = - JMh * cosbetarot * (0.5_real64 * (IB0TE + IB0TM))

      else !r is not zero
        smallr: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .true.
          call prepare_refcoef(refl_var,r,hmd,aniso)

          !evaluate all integrals just for this radius
          ibesord = 0
          IB0TE = compute_1val(funcB0TE,r,sz,zr,ibesord,wellbehaved(1))
          IB0TM = compute_1val(funcB0TM,r,sz,zr,ibesord,wellbehaved(2))
          ibesord = 1
          IB1TE = compute_1val(funcB1TE,r,sz,zr,ibesord,wellbehaved(3))
          IB1TM = compute_1val(funcB1TM,r,sz,zr,ibesord,wellbehaved(4))

        else !r larger than threshold radius for spline interpolation

          !get integral values for this radius by interpolation from precomputed integral values
          IB0TE = splinterp_1val(refl_var,1,r)
          IB0TM = splinterp_1val(refl_var,2,r)
          IB1TE = splinterp_1val(refl_var,3,r)
          IB1TM = splinterp_1val(refl_var,4,r)
        endif smallr

        !Er and Ebeta relative to horizontal source orientation (factor 4pi and dipole length absorbed into JEh)
        Er = - JMh * sinbetarot * (IB0TM + (IB1TE - IB1TM)/r)
        Ebeta = JMh * cosbetarot * (-IB0TE + (IB1TE - IB1TM)/r)
      endif r_is_zero

      !Ex, Ey in original x-y coordinate system
      Ex(recidx) = Ex(recidx) + cosbeta*Er - sinbeta*Ebeta
      Ey(recidx) = Ey(recidx) + sinbeta*Er + cosbeta*Ebeta

      !VTI: deriv. for epsv
      dvert: if(with_dvert) then
        r_is_zerov: if(r.eq.0._real64) then
          !already cycled if receiver is exactly at source point
          !if(sz_eq_zr) cycle

          !Bessel function J0(r=0) = 1, J1(r=0) = 0
          !--> need to evaluate zero-order integrals only
          !use adaptive integration
          IB0TM = compute_1valr0(funcB0TMv)
          Er = - JMh * sinbetarot * 0.5_real64 * IB0TM
          Ebeta = - JMh * cosbetarot * 0.5_real64 * IB0TM
        else !r is not zero
          smallrv: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,hmd,aniso)

            !evaluate all integrals just for this radius
            ibesord = 0
            IB0TM = compute_1val(funcB0TMv,r,sz,zr,ibesord,wellbehaved(13))
            ibesord = 1
            IB1TM = compute_1val(funcB1TMv,r,sz,zr,ibesord,wellbehaved(12))
          else !r larger than threshold radius for spline interpolation
            !get integral values for this radius by interpolation from precomputed integral values
            IB0TM = splinterp_1val(refl_var,13,r)
            IB1TM = splinterp_1val(refl_var,12,r)
          endif smallrv

          !Er and Ebeta relative to horizontal source orientation (factor 4pi and dipole length absorbed into JEh)
          Er = - JMh * sinbetarot * (IB0TM - IB1TM/r)
          Ebeta = JMh * cosbetarot * (- IB1TM/r)
        endif r_is_zerov

        !Ex, Ey in original x-y coordinate system
        Exv(recidx) = Exv(recidx) + cosbeta*Er - sinbeta*Ebeta
        Eyv(recidx) = Eyv(recidx) + sinbeta*Er + cosbeta*Ebeta
      endif dvert

    enddo !irec

    else

      have_ex: if(bgdat%nEx .gt. 0) then

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

       r_is_zeroex: if(r.eq.0._real64) then

        !quick & dirty: skip the point if receiver is right at source point
        if(sz_eq_zr) then
          if(refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle receiver right at source point yet!'
          !skip numerical integral evaluations - they would most likely fail, need to look at this...
          cycle
        endif

        !Bessel function J0(r=0) = 1, J1(r=0) = 0
        !--> need to evaluate zero-order integrals only
        !use adaptive integration
        IB0TE = compute_1valr0(funcB0TE)
        IB0TM = compute_1valr0(funcB0TM)

        Er = - JMh * sinbetarot * (0.5_real64 * (IB0TE + IB0TM))
        Ebeta = - JMh * cosbetarot * (0.5_real64 * (IB0TE + IB0TM))

      else !r is not zero
        smallrex: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .true.
          call prepare_refcoef(refl_var,r,hmd,aniso)

          !evaluate all integrals just for this radius
          ibesord = 0
          IB0TE = compute_1val(funcB0TE,r,sz,zr,ibesord,wellbehaved(1))
          IB0TM = compute_1val(funcB0TM,r,sz,zr,ibesord,wellbehaved(2))
          ibesord = 1
          IB1TE = compute_1val(funcB1TE,r,sz,zr,ibesord,wellbehaved(3))
          IB1TM = compute_1val(funcB1TM,r,sz,zr,ibesord,wellbehaved(4))

        else !r larger than threshold radius for spline interpolation

          !get integral values for this radius by interpolation from precomputed integral values
          IB0TE = splinterp_1val(refl_var,1,r)
          IB0TM = splinterp_1val(refl_var,2,r)
          IB1TE = splinterp_1val(refl_var,3,r)
          IB1TM = splinterp_1val(refl_var,4,r)
        endif smallrex

        !Er and Ebeta relative to horizontal source orientation (factor 4pi and dipole length absorbed into JEh)
        Er = - JMh * sinbetarot * (IB0TM + (IB1TE - IB1TM)/r)
        Ebeta = JMh * cosbetarot * (-IB0TE + (IB1TE - IB1TM)/r)
      endif r_is_zeroex

      !Ex, Ey in original x-y coordinate system
      Ex(recidx) = Ex(recidx) + cosbeta*Er - sinbeta*Ebeta

      !VTI: deriv. for epsv
      dvertex: if(with_dvert) then
        r_is_zerovex: if(r.eq.0._real64) then
          !already cycled if receiver is exactly at source point
          !if(sz_eq_zr) cycle

          !Bessel function J0(r=0) = 1, J1(r=0) = 0
          !--> need to evaluate zero-order integrals only
          !use adaptive integration
          IB0TM = compute_1valr0(funcB0TMv)
          Er = - JMh * sinbetarot * 0.5_real64 * IB0TM
          Ebeta = - JMh * cosbetarot * 0.5_real64 * IB0TM
        else !r is not zero
          smallrvex: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,hmd,aniso)

            !evaluate all integrals just for this radius
            ibesord = 0
            IB0TM = compute_1val(funcB0TMv,r,sz,zr,ibesord,wellbehaved(13))
            ibesord = 1
            IB1TM = compute_1val(funcB1TMv,r,sz,zr,ibesord,wellbehaved(12))
          else !r larger than threshold radius for spline interpolation
            !get integral values for this radius by interpolation from precomputed integral values
            IB0TM = splinterp_1val(refl_var,13,r)
            IB1TM = splinterp_1val(refl_var,12,r)
          endif smallrvex

          !Er and Ebeta relative to horizontal source orientation (factor 4pi and dipole length absorbed into JEh)
          Er = - JMh * sinbetarot * (IB0TM - IB1TM/r)
          Ebeta = JMh * cosbetarot * (- IB1TM/r)
        endif r_is_zerovex

        !Ex, Ey in original x-y coordinate system
        Exv(recidx) = Exv(recidx) + cosbeta*Er - sinbeta*Ebeta
      endif dvertex
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
      sinbeta = sin(beta)
      betarot = beta - refl_var%betasrc(isrc)
      cosbetarot = cos(betarot)
      sinbetarot = sin(betarot)

       r_is_zeroey: if(r.eq.0._real64) then

        !quick & dirty: skip the point if receiver is right at source point
        if(sz_eq_zr) then
          if(refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle receiver right at source point yet!'
          !skip numerical integral evaluations - they would most likely fail, need to look at this...
          cycle
        endif

        !Bessel function J0(r=0) = 1, J1(r=0) = 0
        !--> need to evaluate zero-order integrals only
        !use adaptive integration
        IB0TE = compute_1valr0(funcB0TE)
        IB0TM = compute_1valr0(funcB0TM)

        Er = - JMh * sinbetarot * (0.5_real64 * (IB0TE + IB0TM))
        Ebeta = - JMh * cosbetarot * (0.5_real64 * (IB0TE + IB0TM))

      else !r is not zero
        smallrey: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .true.
          call prepare_refcoef(refl_var,r,hmd,aniso)

          !evaluate all integrals just for this radius
          ibesord = 0
          IB0TE = compute_1val(funcB0TE,r,sz,zr,ibesord,wellbehaved(1))
          IB0TM = compute_1val(funcB0TM,r,sz,zr,ibesord,wellbehaved(2))
          ibesord = 1
          IB1TE = compute_1val(funcB1TE,r,sz,zr,ibesord,wellbehaved(3))
          IB1TM = compute_1val(funcB1TM,r,sz,zr,ibesord,wellbehaved(4))

        else !r larger than threshold radius for spline interpolation

          !get integral values for this radius by interpolation from precomputed integral values
          IB0TE = splinterp_1val(refl_var,1,r)
          IB0TM = splinterp_1val(refl_var,2,r)
          IB1TE = splinterp_1val(refl_var,3,r)
          IB1TM = splinterp_1val(refl_var,4,r)
        endif smallrey

        !Er and Ebeta relative to horizontal source orientation (factor 4pi and dipole length absorbed into JEh)
        Er = - JMh * sinbetarot * (IB0TM + (IB1TE - IB1TM)/r)
        Ebeta = JMh * cosbetarot * (-IB0TE + (IB1TE - IB1TM)/r)
      endif r_is_zeroey

      !Ex, Ey in original x-y coordinate system
      Ey(recidx) = Ey(recidx) + sinbeta*Er + cosbeta*Ebeta

      !VTI: deriv. for epsv
      dvertey: if(with_dvert) then
        r_is_zerovey: if(r.eq.0._real64) then
          !already cycled if receiver is exactly at source point
          !if(sz_eq_zr) cycle

          !Bessel function J0(r=0) = 1, J1(r=0) = 0
          !--> need to evaluate zero-order integrals only
          !use adaptive integration
          IB0TM = compute_1valr0(funcB0TMv)
          Er = - JMh * sinbetarot * 0.5_real64 * IB0TM
          Ebeta = - JMh * cosbetarot * 0.5_real64 * IB0TM
        else !r is not zero
          smallrvey: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,hmd,aniso)

            !evaluate all integrals just for this radius
            ibesord = 0
            IB0TM = compute_1val(funcB0TMv,r,sz,zr,ibesord,wellbehaved(13))
            ibesord = 1
            IB1TM = compute_1val(funcB1TMv,r,sz,zr,ibesord,wellbehaved(12))
          else !r larger than threshold radius for spline interpolation
            !get integral values for this radius by interpolation from precomputed integral values
            IB0TM = splinterp_1val(refl_var,13,r)
            IB1TM = splinterp_1val(refl_var,12,r)
          endif smallrvey

          !Er and Ebeta relative to horizontal source orientation (factor 4pi and dipole length absorbed into JEh)
          Er = - JMh * sinbetarot * (IB0TM - IB1TM/r)
          Ebeta = JMh * cosbetarot * (- IB1TM/r)
        endif r_is_zerovey

        !Ex, Ey in original x-y coordinate system
        Eyv(recidx) = Eyv(recidx) + sinbeta*Er + cosbeta*Ebeta
      endif dvertey

    enddo !irec

      endif have_ey
    endif exy_equalpos

  enddo  !source elements

endsubroutine interp_intvals_hmd_Exy


!------------------------------------------------------------
!  1D EM subroutine interp_intvals_hmd_Ez
!
!  get interpolated field values at receiver locations, HMD source,
!  Ez only
!
!  Rita Streich 2009-2011
!------------------------------------------------------------
subroutine interp_intvals_hmd_Ez(refl_var,src,ifreq,sz,zr,bgdat,Ez,omeps_recv,j_om_mu, funcCz1TM,ilay,funcCz1TMfwd, funcCz1TMv,Ezv)

  implicit none

  !external variables
  type(refl_struct)               :: refl_var   !everything needed throughout 1D computations
  type(sorec),intent(in)          :: src        !source definition (we need the currents here)
  integer(kind=int32),intent(in)  :: ifreq      !frequency index
  real(kind=real64),intent(in)    :: sz,zr      !source and receiver depth
  type(backgrounddata) :: bgdat      !coordinate vectors and output EM fields
    !input separate fields because the same routine is used both for fields and derivatives!
  complex(kind=real64),dimension(:),target   :: Ez  !electric fields: nr of receivers
  complex(kind=real64),intent(in) :: omeps_recv  !omega * epsilon in receiver layer
  complex(kind=real64),intent(in) :: j_om_mu    !j * omega * mu0
  complex(kind=real64),external   :: funcCz1TM
  integer(kind=int32),intent(in)  :: ilay       !layer index for derivatives - leave at zero for forward modeling
  complex(kind=real64),external,optional   :: funcCz1TMfwd  !function only needed for derivatives in receiver layer
  complex(kind=real64),external,optional   :: funcCz1TMv !derivative integrals for epsv
  complex(kind=real64),dimension(:),target,optional   :: Ezv  !electric fields: nr of receivers

  !internal variables
  integer(kind=int32)   :: isrc    !source element counter
  real(kind=real64)     :: x,y,r   !temp source-receiver distances
  integer(kind=int32)   :: irec    !receiver counter
  real(kind=real64)     :: beta,betarot      !temp angles

  complex(kind=real64)          :: ICz1TM !interpolated integral values

  logical,dimension(nintHMDdvti)    :: wellbehaved  !indicates if Hankel integration can be used
  logical                       :: sz_eq_zr     !indicates if source and reveicer are at the same depth
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")

  real(kind=real64)     :: cosbeta,sinbeta      !cos(beta) and sin(beta), precompute for efficiency
  real(kind=real64)     :: cosbetarot,sinbetarot      !cos(betarot) and sin(betarot), precompute for efficiency
  complex(kind=real64)  :: cur      !temp source current
  complex(kind=real64)  :: JMh      !source current times constants
  integer(kind=int32)   :: idx      !source element index
  integer(kind=int32)   :: recidx   !receiver index
  !flag for computing epsv derivatives, not needed for forward computation, so "aniso" value cannot be used here
  logical               :: with_dvert
  complex(kind=real64),dimension(:),pointer   :: Erec  !points to either Ez (isotropic) or Ezv (VTI)

  !indicators for fast Hankel transform or adaptive integration
  wellbehaved = .true.
  sz_eq_zr = .false.
  if(sz.eq.zr) then
    wellbehaved(1:2) = .false.
    wellbehaved(6) = .false.
    wellbehaved(10) = .false.
    sz_eq_zr = .true.
  endif

  if(present(funcCz1TMv)) then
    with_dvert = .true.
    if(sz_eq_zr) then
      wellbehaved(13) = .false.
    endif
  else
    with_dvert = .false.
  endif


  do isrc = refl_var%isrcstart,refl_var%isrcend

    !get source current times constant factors
    idx = refl_var%isrcperz(isrc)
    cur = conjg(src%cur(idx,ifreq))
    JMh = - cur * sqrt(src%akx(idx)**2 + src%aky(idx)**2) * j_om_mu / dfourpi

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
 
      r_is_zero: if(r.eq.0._real64) then

        !quick & dirty: skip the point if receiver is right at source point
        if(sz_eq_zr) then
          if(refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle receiver right at source point yet!'
          cycle
        endif

      else !r is not zero
        smallr: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .true.
          call prepare_refcoef(refl_var,r,hmd,aniso)

          !evaluate all integrals just for this radius
          ibesord = 1
          ICz1TM = compute_1val(funcCz1TM,r,sz,zr,ibesord,wellbehaved(5))
        else !r larger than threshold radius for spline interpolation
          !get integral values for this radius by interpolation from precomputed integral values
          ICz1TM = splinterp_1val(refl_var,5,r)
        endif smallr

        !Ez in original x-y coordinate system, Ez is zero for r=0
        Ez(recidx) = Ez(recidx) - (JMh/(dci*omeps_recv)) * sinbetarot * ICz1TM
      endif r_is_zero

      !VTI: deriv. for epsv
      dvert: if(with_dvert) then
        r_is_zerov: if(r.ne.0._real64) then
          smallrv: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,hmd,aniso)

            ibesord = 1
            ICz1TM = compute_1val(funcCz1TMv,r,sz,zr,ibesord,wellbehaved(16))
          else !r larger than threshold radius for spline interpolation
            ICz1TM = splinterp_1val(refl_var,16,r)
          endif smallrv

          !Ez in original x-y coordinate system, Ez is zero for r=0
          Ezv(recidx) = Ezv(recidx) - (JMh/(dci*omeps_recv)) * sinbetarot * ICz1TM
        endif r_is_zerov
      endif dvert

    enddo !irec
  enddo  !source elements


  !special term for Ez for derivatives in receiver layer
  deriv_ilayrec: if(ilay .eq. ilayrec) then
  
    if(.not. present(funcCz1TMfwd)) then
      write(*,'(a)') 'ERROR: function for Ez in receiver layer not given, cannot compute correct Ez derivative!'
      return
    endif

    !there is a factor epsv in the term before the Ez integral
    !-> for isotropic case, epsv becomes eps and contribution is added to E
    !-> for VTI case, derivative of this term only exists for epsv -> add to Ev, nothing added to E
    if(with_dvert) then
      Erec => Ezv
    else
      Erec => Ez
    endif
  
    wellbehaved(11) = wellbehaved(5)

    do isrc = refl_var%isrcstart,refl_var%isrcend

      !get source current times constant factors
      idx = refl_var%isrcperz(isrc)
      cur = conjg(src%cur(idx,ifreq))
      JMh = - cur * sqrt(src%akx(idx)**2 + src%aky(idx)**2) * j_om_mu / dfourpi

      do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzEz(irec)

        y = bgdat%Ezpos(recidx,2) - refl_var%ys(isrc)
        x = bgdat%Ezpos(recidx,1) - refl_var%xs(isrc)
        r = sqrt(x**2 + y**2)

        !Ez is zero for r=0
        if(r.eq.0._real64) then
          cycle
        else
          beta = atan2(y,x)
          betarot = beta - refl_var%betasrc(isrc)
          sinbetarot = sin(betarot)

          if(r.lt.rsplmin) then
            !reflection coeff. for this radius
            call prepare_refcoef(refl_var,r,ved,aniso) !use "ved" to compute TM refl. coeff. only
            ibesord = 1
            ICz1TM = compute_1val(funcCz1TMfwd,r,sz,zr,ibesord,wellbehaved(11))
          else
            ICz1TM = splinterp_1val(refl_var,11,r)
          endif

          !Ez in original x-y coordinate system
          Erec(recidx) = Erec(recidx) + (JMh/(dci*omeps_recv*epsv(ilayrec))) * sinbetarot * ICz1TM
        endif

      enddo !irec
    enddo  !source elements

  endif deriv_ilayrec

endsubroutine interp_intvals_hmd_Ez


!------------------------------------------------------------
!  1D EM subroutine interp_intvals_hmd_Hxy
!
!  get interpolated field values at receiver locations, HMD source,
!  Hx and / or Hy only
!
!  Rita Streich 2009-2011
!------------------------------------------------------------
subroutine interp_intvals_hmd_Hxy(refl_var,src,ifreq,sz,zr,bgdat,Hx,Hy,j_om_mu, funcC0TE,funcC0TM,funcC1TE,funcC1TM,ilay, &
  funcC0TMv,funcC1TMv,Hxv,Hyv)

  implicit none

  !external variables
  type(refl_struct)               :: refl_var   !everything needed throughout 1D computations
  type(sorec),intent(in)          :: src        !source definition (we need the currents here)
  integer(kind=int32),intent(in)  :: ifreq      !frequency index
  real(kind=real64),intent(in)    :: sz,zr      !source and receiver depth
  type(backgrounddata) :: bgdat      !coordinate vectors and output EM fields
    !input separate fields because the same routine is used both for fields and derivatives!
  complex(kind=real64),dimension(:)   :: Hx,Hy  !electric and magnetic fields: nr of receivers
  complex(kind=real64),intent(in) :: j_om_mu    !j * omega * mu0
  complex(kind=real64),external   :: funcC0TE,funcC0TM,funcC1TE,funcC1TM
  integer(kind=int32),intent(in)  :: ilay       !layer index for derivatives - leave at zero for forward modeling
  complex(kind=real64),external,optional   :: funcC0TMv,funcC1TMv !derivative integrals for epsv
  complex(kind=real64),dimension(:),optional   :: Hxv,Hyv  !electric and magnetic fields: nr of receivers

  !internal variables
  integer(kind=int32)   :: isrc    !source element counter
  real(kind=real64)     :: x,y,r   !temp source-receiver distances
  integer(kind=int32)   :: irec    !receiver counter
  real(kind=real64)     :: beta,betarot      !temp angles

  complex(kind=real64)          :: IC0TE,IC0TM,IC1TE,IC1TM !interpolated integral values

  logical,dimension(nintHMDdvti)    :: wellbehaved  !indicates if Hankel integration can be used
  logical                       :: sz_eq_zr     !indicates if source and reveicer are at the same depth
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")

  complex(kind=real64)  :: Hr,Hbeta             !temp field values in cylindrical coordinates
  real(kind=real64)     :: cosbeta,sinbeta      !cos(beta) and sin(beta), precompute for efficiency
  real(kind=real64)     :: cosbetarot,sinbetarot      !cos(betarot) and sin(betarot), precompute for efficiency
  complex(kind=real64)  :: cur      !temp source current
  complex(kind=real64)  :: JMh      !source current times constants
  integer(kind=int32)   :: idx      !source element index
  integer(kind=int32)   :: recidx   !receiver index
  !flag for computing epsv derivatives, not needed for forward computation, so "aniso" value cannot be used here
  logical               :: with_dvert

  !indicators for fast Hankel transform or adaptive integration
  wellbehaved = .true.
  sz_eq_zr = .false.
  if(sz.eq.zr) then
    wellbehaved(1:2) = .false.
    wellbehaved(6) = .false.
    wellbehaved(10) = .false.
    sz_eq_zr = .true.
  endif

  if(present(funcC1TMv)) then
    with_dvert = .true.
    if(sz_eq_zr) then
      wellbehaved(13) = .false.
    endif
  else
    with_dvert = .false.
  endif


  do isrc = refl_var%isrcstart,refl_var%isrcend

    !get source current times constant factors
    idx = refl_var%isrcperz(isrc)
    cur = conjg(src%cur(idx,ifreq))
    JMh = - cur * sqrt(src%akx(idx)**2 + src%aky(idx)**2) * j_om_mu / dfourpi


    !same positions for Hx and Hy
    hxy_equalpos: if(bgdat%nHxy.gt.0) then

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

 
      r_is_zero: if(r.eq.0._real64) then

        !quick & dirty: skip the point if receiver is right at source point
        if(sz_eq_zr) then
          if(refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle receiver right at source point yet!'
          cycle
        endif

        !Bessel function J0(r=0) = 1, J1(r=0) = 0
        !--> need to evaluate zero-order integrals only
        !use adaptive integration
        IC0TE = compute_1valr0(funcC0TE)
        IC0TM = compute_1valr0(funcC0TM)

        Hr = JMh * cosbetarot * (0.5_real64 * (IC0TE + IC0TM))
        Hbeta = - JMh * sinbetarot * (0.5_real64 * (IC0TE + IC0TM))
      else !r is not zero
        smallr: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .true.
          call prepare_refcoef(refl_var,r,hmd,aniso)

          !evaluate all integrals just for this radius
          ibesord = 0
          IC0TE = compute_1val(funcC0TE,r,sz,zr,ibesord,wellbehaved(6))
          IC0TM = compute_1val(funcC0TM,r,sz,zr,ibesord,wellbehaved(7))
          ibesord = 1
          IC1TE = compute_1val(funcC1TE,r,sz,zr,ibesord,wellbehaved(8))
          IC1TM = compute_1val(funcC1TM,r,sz,zr,ibesord,wellbehaved(9))
        else !r larger than threshold radius for spline interpolation
          !get integral values for this radius by interpolation from precomputed integral values
          IC0TE = splinterp_1val(refl_var,6,r)
          IC0TM = splinterp_1val(refl_var,7,r)
          IC1TE = splinterp_1val(refl_var,8,r)
          IC1TM = splinterp_1val(refl_var,9,r)
        endif smallr

        Hr = - JMh * cosbetarot * (-IC0TE + (IC1TE - IC1TM)/r)
        Hbeta = - JMh * sinbetarot * (IC0TM + (IC1TE - IC1TM)/r)
      endif r_is_zero

      !Hx, Hy in original x-y coordinate system
      Hx(recidx) = Hx(recidx) + cosbeta*Hr - sinbeta*Hbeta
      Hy(recidx) = Hy(recidx) + sinbeta*Hr + cosbeta*Hbeta

      !VTI: deriv. for epsv
      dvert: if(with_dvert) then
        r_is_zerov: if(r.eq.0._real64) then
          !already cycled if receiver is exactly at source point
          !if(sz_eq_zr) cycle
          IC0TM = compute_1valr0(funcC0TMv)

          Hr = JMh * cosbetarot * 0.5_real64 * IC0TM
          Hbeta = - JMh * sinbetarot * 0.5_real64 * IC0TM
        else !r is not zero
          smallrv: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,hmd,aniso)

            ibesord = 0
            IC0TM = compute_1val(funcC0TMv,r,sz,zr,ibesord,wellbehaved(15))
            ibesord = 1
            IC1TM = compute_1val(funcC1TMv,r,sz,zr,ibesord,wellbehaved(14))
          else !r larger than threshold radius for spline interpolation
            IC0TM = splinterp_1val(refl_var,15,r)
            IC1TM = splinterp_1val(refl_var,14,r)
          endif smallrv

          Hr = - JMh * cosbetarot * (- IC1TM/r)
          Hbeta = - JMh * sinbetarot * (IC0TM - IC1TM/r)
        endif r_is_zerov

        !Hx, Hy in original x-y coordinate system
        Hxv(recidx) = Hxv(recidx) + cosbeta*Hr - sinbeta*Hbeta
        Hyv(recidx) = Hyv(recidx) + sinbeta*Hr + cosbeta*Hbeta
      endif dvert

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
      sinbeta = sin(beta)
      betarot = beta - refl_var%betasrc(isrc)
      cosbetarot = cos(betarot)
      sinbetarot = sin(betarot)
 
      r_is_zerohx: if(r.eq.0._real64) then

        !quick & dirty: skip the point if receiver is right at source point
        if(sz_eq_zr) then
          if(refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle receiver right at source point yet!'
          cycle
        endif

        !Bessel function J0(r=0) = 1, J1(r=0) = 0
        !--> need to evaluate zero-order integrals only
        !use adaptive integration
        IC0TE = compute_1valr0(funcC0TE)
        IC0TM = compute_1valr0(funcC0TM)

        Hr = JMh * cosbetarot * (0.5_real64 * (IC0TE + IC0TM))
        Hbeta = - JMh * sinbetarot * (0.5_real64 * (IC0TE + IC0TM))
      else !r is not zero
        smallrhx: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .true.
          call prepare_refcoef(refl_var,r,hmd,aniso)

          !evaluate all integrals just for this radius
          ibesord = 0
          IC0TE = compute_1val(funcC0TE,r,sz,zr,ibesord,wellbehaved(6))
          IC0TM = compute_1val(funcC0TM,r,sz,zr,ibesord,wellbehaved(7))
          ibesord = 1
          IC1TE = compute_1val(funcC1TE,r,sz,zr,ibesord,wellbehaved(8))
          IC1TM = compute_1val(funcC1TM,r,sz,zr,ibesord,wellbehaved(9))
        else !r larger than threshold radius for spline interpolation
          !get integral values for this radius by interpolation from precomputed integral values
          IC0TE = splinterp_1val(refl_var,6,r)
          IC0TM = splinterp_1val(refl_var,7,r)
          IC1TE = splinterp_1val(refl_var,8,r)
          IC1TM = splinterp_1val(refl_var,9,r)
        endif smallrhx

        Hr = - JMh * cosbetarot * (-IC0TE + (IC1TE - IC1TM)/r)
        Hbeta = - JMh * sinbetarot * (IC0TM + (IC1TE - IC1TM)/r)
      endif r_is_zerohx

      !Hx, Hy in original x-y coordinate system
      Hx(recidx) = Hx(recidx) + cosbeta*Hr - sinbeta*Hbeta

      !VTI: deriv. for epsv
      dverthx: if(with_dvert) then
        r_is_zerovhx: if(r.eq.0._real64) then
          !already cycled if receiver is exactly at source point
          !if(sz_eq_zr) cycle
          IC0TM = compute_1valr0(funcC0TMv)

          Hr = JMh * cosbetarot * 0.5_real64 * IC0TM
          Hbeta = - JMh * sinbetarot * 0.5_real64 * IC0TM
        else !r is not zero
          smallrvhx: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,hmd,aniso)

            ibesord = 0
            IC0TM = compute_1val(funcC0TMv,r,sz,zr,ibesord,wellbehaved(15))
            ibesord = 1
            IC1TM = compute_1val(funcC1TMv,r,sz,zr,ibesord,wellbehaved(14))
          else !r larger than threshold radius for spline interpolation
            IC0TM = splinterp_1val(refl_var,15,r)
            IC1TM = splinterp_1val(refl_var,14,r)
          endif smallrvhx

          Hr = - JMh * cosbetarot * (- IC1TM/r)
          Hbeta = - JMh * sinbetarot * (IC0TM - IC1TM/r)
        endif r_is_zerovhx

        !Hx, Hy in original x-y coordinate system
        Hxv(recidx) = Hxv(recidx) + cosbeta*Hr - sinbeta*Hbeta
      endif dverthx
    enddo !irec

      endif have_hx

      have_hy: if(bgdat%nHy .gt. 0) then

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

       r_is_zerohy: if(r.eq.0._real64) then

        !quick & dirty: skip the point if receiver is right at source point
        if(sz_eq_zr) then
          if(refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle receiver right at source point yet!'
          cycle
        endif

        !Bessel function J0(r=0) = 1, J1(r=0) = 0
        !--> need to evaluate zero-order integrals only
        !use adaptive integration
        IC0TE = compute_1valr0(funcC0TE)
        IC0TM = compute_1valr0(funcC0TM)

        Hr = JMh * cosbetarot * (0.5_real64 * (IC0TE + IC0TM))
        Hbeta = - JMh * sinbetarot * (0.5_real64 * (IC0TE + IC0TM))
      else !r is not zero
        smallrhy: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .true.
          call prepare_refcoef(refl_var,r,hmd,aniso)

          !evaluate all integrals just for this radius
          ibesord = 0
          IC0TE = compute_1val(funcC0TE,r,sz,zr,ibesord,wellbehaved(6))
          IC0TM = compute_1val(funcC0TM,r,sz,zr,ibesord,wellbehaved(7))
          ibesord = 1
          IC1TE = compute_1val(funcC1TE,r,sz,zr,ibesord,wellbehaved(8))
          IC1TM = compute_1val(funcC1TM,r,sz,zr,ibesord,wellbehaved(9))
        else !r larger than threshold radius for spline interpolation
          !get integral values for this radius by interpolation from precomputed integral values
          IC0TE = splinterp_1val(refl_var,6,r)
          IC0TM = splinterp_1val(refl_var,7,r)
          IC1TE = splinterp_1val(refl_var,8,r)
          IC1TM = splinterp_1val(refl_var,9,r)
        endif smallrhy

        Hr = - JMh * cosbetarot * (-IC0TE + (IC1TE - IC1TM)/r)
        Hbeta = - JMh * sinbetarot * (IC0TM + (IC1TE - IC1TM)/r)
      endif r_is_zerohy

      !Hx, Hy in original x-y coordinate system
      Hy(recidx) = Hy(recidx) + sinbeta*Hr + cosbeta*Hbeta

      !VTI: deriv. for epsv
      dverthy: if(with_dvert) then
        r_is_zerovhy: if(r.eq.0._real64) then
          !already cycled if receiver is exactly at source point
          !if(sz_eq_zr) cycle
          IC0TM = compute_1valr0(funcC0TMv)

          Hr = JMh * cosbetarot * 0.5_real64 * IC0TM
          Hbeta = - JMh * sinbetarot * 0.5_real64 * IC0TM
        else !r is not zero
          smallrvhy: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .true.
            call prepare_refcoef(refl_var,r,hmd,aniso)

            ibesord = 0
            IC0TM = compute_1val(funcC0TMv,r,sz,zr,ibesord,wellbehaved(15))
            ibesord = 1
            IC1TM = compute_1val(funcC1TMv,r,sz,zr,ibesord,wellbehaved(14))
          else !r larger than threshold radius for spline interpolation
            IC0TM = splinterp_1val(refl_var,15,r)
            IC1TM = splinterp_1val(refl_var,14,r)
          endif smallrvhy

          Hr = - JMh * cosbetarot * (- IC1TM/r)
          Hbeta = - JMh * sinbetarot * (IC0TM - IC1TM/r)
        endif r_is_zerovhy

        !Hx, Hy in original x-y coordinate system
        Hyv(recidx) = Hyv(recidx) + sinbeta*Hr + cosbeta*Hbeta
      endif dverthy
    enddo !irec

      endif have_hy
    endif hxy_equalpos

  enddo  !source elements

endsubroutine interp_intvals_hmd_Hxy


!------------------------------------------------------------
!  1D EM subroutine interp_intvals_hmd_Hz
!
!  get interpolated field values at receiver locations, HMD source,
!  Hz only
!
!  Rita Streich 2009-2011
!------------------------------------------------------------
subroutine interp_intvals_hmd_Hz(refl_var,src,ifreq,sz,zr,bgdat,Hz,j_om_mu, funcBz1TE,ilay)

  implicit none

  !external variables
  type(refl_struct)               :: refl_var   !everything needed throughout 1D computations
  type(sorec),intent(in)          :: src        !source definition (we need the currents here)
  integer(kind=int32),intent(in)  :: ifreq      !frequency index
  real(kind=real64),intent(in)    :: sz,zr      !source and receiver depth
  type(backgrounddata) :: bgdat      !coordinate vectors and output EM fields
    !input separate fields because the same routine is used both for fields and derivatives!
  complex(kind=real64),dimension(:)   :: Hz  !electric and magnetic fields: nr of receivers
  complex(kind=real64),intent(in) :: j_om_mu    !j * omega * mu0
  complex(kind=real64),external   :: funcBz1TE
  integer(kind=int32),intent(in)  :: ilay       !layer index for derivatives - leave at zero for forward modeling

  !internal variables
  integer(kind=int32)   :: isrc    !source element counter
  real(kind=real64)     :: x,y,r   !temp source-receiver distances
  integer(kind=int32)   :: irec    !receiver counter
  real(kind=real64)     :: beta,betarot      !temp angles

  complex(kind=real64)          :: IBz1TE !interpolated integral values

  logical,dimension(nintHMDdvti)    :: wellbehaved  !indicates if Hankel integration can be used
  logical                       :: sz_eq_zr     !indicates if source and reveicer are at the same depth
  integer(kind=int32),dimension(NREL)  :: ibesord   !bessel function order (integer "array")

  real(kind=real64)     :: cosbetarot      !cos(betarot) and sin(betarot), precompute for efficiency
  complex(kind=real64)  :: cur      !temp source current
  complex(kind=real64)  :: JMh      !source current times constants
  integer(kind=int32)   :: idx      !source element index
  integer(kind=int32)   :: recidx   !receiver index

  !indicators for fast Hankel transform or adaptive integration
  wellbehaved = .true.
  sz_eq_zr = .false.
  if(sz.eq.zr) then
    wellbehaved(1:2) = .false.
    wellbehaved(6) = .false.
    wellbehaved(10) = .false.
    sz_eq_zr = .true.
  endif

  do isrc = refl_var%isrcstart,refl_var%isrcend

    !get source current times constant factors
    idx = refl_var%isrcperz(isrc)
    cur = conjg(src%cur(idx,ifreq))
    JMh = - cur * sqrt(src%akx(idx)**2 + src%aky(idx)**2) * j_om_mu / dfourpi

    do irec=refl_var%irecstart,refl_var%irecend
      recidx = refl_var%irecperzHz(irec)

      y = bgdat%Hzpos(recidx,2) - refl_var%ys(isrc)
      x = bgdat%Hzpos(recidx,1) - refl_var%xs(isrc)
      r = sqrt(x**2 + y**2)

      beta = atan2(y,x)
      betarot = beta - refl_var%betasrc(isrc)
      cosbetarot = cos(betarot)
 
      r_is_zero: if(r.eq.0._real64) then

        !quick & dirty: skip the point if receiver is right at source point
        if(sz_eq_zr) then
          if(refl_var%infolevel.ge.output_more) &
            write(*,'(a)') 'WARNING: cannot handle receiver right at source point yet!'
          cycle
        endif

      else !r is not zero
        smallr: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

          !reflection coeff. for this radius
          refl_var%refcoef_changed = .true.
          call prepare_refcoef(refl_var,r,hmd,aniso)

          ibesord = 1
          IBz1TE = compute_1val(funcBz1TE,r,sz,zr,ibesord,wellbehaved(10))
        else !r larger than threshold radius for spline interpolation
          IBz1TE = splinterp_1val(refl_var,10,r)
        endif smallr

        !Hz in original x-y coordinate system
        Hz(recidx) = Hz(recidx) + (JMh / j_om_mu) * cosbetarot * IBz1TE
      endif r_is_zero

    enddo !irec
  enddo  !source elements

endsubroutine interp_intvals_hmd_Hz


