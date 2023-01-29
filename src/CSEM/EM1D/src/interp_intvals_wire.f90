!------------------------------------------------------------
!>  1D EM subroutine interp_intvals_wire_allcomp
!
!>  get interpolated field values at receiver locations, horizontal wire source,
!>    all field components at the same cordinates
!
!>  Rita Streich 2009-2011
!------------------------------------------------------------
subroutine interp_intvals_wire_allcomp(refl_var,bgdat,fld,src,sz,zr,omeps_recv,ommu, &
  func1Exwire,funcD0TE,funcAz1TE,func2Exwire,funcD0TM,funcdHxwire,ilay,funcD0TMfwd, &
  func2Exwirev,funcD0TMv,funcdHxwirev,fldv)

  implicit none

  !external variables
  type(refl_struct) :: refl_var   !everything needed throughout 1D computations
  type(backgrounddata) :: bgdat      !coordinate vectors and final output EM fields
  type(receiverdata),dimension(:) :: fld        !field or derivative vectors for each wire, all components
  type(sorec),intent( in ) :: src        !source definition (we need the currents here)
  real(kind=real64),intent( in ) :: sz,zr      !source and receiver depth
  complex(kind=real64),intent( in ) :: omeps_recv  !omega * eps in receiver layer
  real(kind=real64),intent( in ) :: ommu       !omega * mu0
  complex(kind=real64),external :: func1Exwire,funcD0TE,funcAz1TE,func2Exwire,funcD0TM,funcdHxwire
  integer(kind=int32),intent( in ) :: ilay       !layer index for derivatives - leave at zero for forward modeling
  complex(kind=real64),external,optional :: funcD0TMfwd  !function only needed for derivatives in receiver layer
  complex(kind=real64),external,optional :: func2Exwirev,funcD0TMv,funcdHxwirev !integral derivatives for epsv
  type(receiverdata),dimension(:),optional :: fldv       !field or derivative vectors for each wire, all components

  !internal variables
  integer(kind=int32) :: isrc    !source element counter
  real(kind=real64) :: x,y,r   !temp source-receiver distances
  integer(kind=int32) :: irec    !receiver counter
  real(kind=real64) :: beta,betarot      !temp angles

  complex(kind=real64) :: IExy,IHyx,IHz          !interpolated integral values for integrals along wire
  complex(kind=real64) :: IendExy,IendEz,IendHxy !interpolated integral values for end points

  logical,dimension(nintWiredvti) :: wellbehaved   !indicates if Hankel integration can be used
  logical :: sz_eq_zr      !indicates if source and reveicer are at the same depth
  integer(kind=int32),dimension(NREL) :: ibesord   !bessel function order (integer "array")

  complex(kind=real64) :: Exrot,Eyrot,Hxrot,Hyrot,Hzrot  !temp field values in cylindrical coordinates
  real(kind=real64) :: cosbetasrc,sinbetasrc !cos(betasrc) and sin(betasrc), precompute for efficiency
  real(kind=real64) :: cosbetarot,sinbetarot !cos(betarot) and sin(betarot), precompute for efficiency
  real(kind=real64) :: const    !geometric constant
  complex(kind=real64) :: constHz  !geometric constant for Hz
  complex(kind=real64) :: constEz  !geometric constant for Ez
  integer(kind=int32) :: idx      !source element index
  integer(kind=int32) :: recidx   !receiver index
  integer(kind=int32) :: ielemall,ielem       !wire element counters
  integer(kind=int32) :: iep      !wire end point counter
  real(kind=real64) :: xs,ys    !wire end point coordinates
  complex(kind=real64) :: term1    !temp result for end point field values
  !signs for wire end points (predefine for efficiency, need just elements 1 and 4)
  real(kind=real64),dimension(1:2),parameter :: signs = (/-1._real64,1._real64/)
  real(kind=real64) :: sign     !sign for end point field values
  !flag for computing epsv derivatives, not needed for forward computation, so "aniso" value cannot be used here
  logical :: with_dvert
  type(receiverdata),dimension(:),pointer :: Ewirerec  !points to Ez for isotropic and Ezv for VTI case
  integer(kind=int32) :: ierr     !error index
  integer(kind=int32) :: iwire    !wire counter


  !indicators for fast Hankel transform or adaptive integration
  wellbehaved = .TRUE.
  sz_eq_zr = .FALSE.
  if(sz.EQ.zr) then
    wellbehaved(1:2) = .FALSE.
    wellbehaved(4:5) = .FALSE.
    sz_eq_zr = .TRUE.
  endif

  if(present(funcdHxwirev)) then
    with_dvert = .TRUE.
    if(sz_eq_zr) then
      wellbehaved(8:9) = .FALSE.
    endif
  else
    with_dvert = .FALSE.
  endif

  !constant for Ez
  constEz = dci / (omeps_recv*dfourpi)


  !loop over wires
  !TE only --> no contribution to epsv derivatives
  ielemall = 0
  wires: do isrc = refl_var%isrcstart,refl_var%isrcend

    !get constant geometry factor: wire element length
    idx = refl_var%isrcperz(isrc)
    const = src%wire(idx)%dlw / dfourpi
    constHz = const * dci / ommu

    cosbetasrc = cos(refl_var%betasrc(isrc))
    sinbetasrc = sin(refl_var%betasrc(isrc))

    !loop over wire elements
    wire_elements: do ielem = 1,src%nelem(idx)
      ielemall = ielemall + 1

      receivers: do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzExy(irec)

        y = bgdat%Exypos(recidx,2) - refl_var%ys(ielemall)
        x = bgdat%Exypos(recidx,1) - refl_var%xs(ielemall)
        r = sqrt(x**2 + y**2)

        beta = atan2(y,x)
        betarot = beta - refl_var%betasrc(isrc)


        r_is_zero: if(r.EQ.0._real64) then

          !quick & dirty: skip the point if receiver is right at source point
          if(sz_eq_zr) then
            if(refl_var%infolevel.ge.output_more) &
              write(*,'(a)') 'WARNING: cannot handle a receiver at exactly the same location as a wire element yet, '// &
                'ignoring this contribution!'
            cycle
          endif

          !--> evaluate zero-order integrals only
          ibesord = 0
          IExy = compute_1valr0(func1Exwire)
          IHyx = compute_1valr0(funcD0TE)
          IHz = 0._real64

        else !r is not zero
          smallr: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .TRUE.
            call prepare_refcoef(refl_var,r,0,aniso)

            !evaluate all integrals just for this radius
            ibesord = 0
            IExy = compute_1val(func1Exwire,r,sz,zr,ibesord,wellbehaved(1))
            IHyx = compute_1val(funcD0TE,r,sz,zr,ibesord,wellbehaved(2))
            ibesord = 1
            IHz = compute_1val(funcAz1TE,r,sz,zr,ibesord,wellbehaved(3))

          else !r larger than threshold radius for spline interpolation

            !get integral values for this radius by interpolation from precomputed integral values
            IExy = splinterp_1val(refl_var,1,r)
            IHyx = splinterp_1val(refl_var,2,r)
            IHz = splinterp_1val(refl_var,3,r)
          endif smallr

        endif r_is_zero

        !field values in rotated coordinate system
        Exrot = IExy * const
        Hyrot = - IHyx * const
        Hzrot = IHz * constHz
 
        !field values rotated back to original x-y coordinate system
        !don't take complex conjugate here, but at the very end!
        fld(idx)%Ex(recidx) = fld(idx)%Ex(recidx) + Exrot*cosbetasrc
        fld(idx)%Ey(recidx) = fld(idx)%Ey(recidx) + Exrot*sinbetasrc

        fld(idx)%Hx(recidx) = fld(idx)%Hx(recidx) - Hyrot*sinbetasrc
        fld(idx)%Hy(recidx) = fld(idx)%Hy(recidx) + Hyrot*cosbetasrc
        fld(idx)%Hz(recidx) = fld(idx)%Hz(recidx) + Hzrot*sin(betarot)

      enddo receivers !iy
    enddo wire_elements


    !loop over wire end points
    wire_endpoints: do iep = 1,2
      xs = src%wire(idx)%endpos(1,iep)
      ys = src%wire(idx)%endpos(2,iep)

      sign = signs(iep)

      receivers1: do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzExy(irec)

        y = bgdat%Exypos(recidx,2) - ys
        x = bgdat%Exypos(recidx,1) - xs
        r = sqrt(x**2 + y**2)

        !distinguish case r=0 since here we need only 1 integral
        chkr0: if(r .EQ. 0._real64) then

          !quick & dirty: skip the point if receiver is right at source point
          if(sz_eq_zr) then
            if(refl_var%infolevel.ge.output_more) &
              write(*,'(a)') 'WARNING: cannot handle a receiver at exactly the same location as a grounding point yet, '// &
                'ignoring this contribution!'
            cycle
          endif

          !only need Bessel function of order 0, since Bessel function of order 1 is zero for r=0
          ibesord = 0
          IendEz = compute_1valr0(funcD0TM)

          term1 = sign * constEz * IendEz
          fld(idx)%Ez(recidx) = fld(idx)%Ez(recidx) + term1

        else !r is not zero

          beta = atan2(y,x)
          betarot = beta - refl_var%betasrc(isrc)
          cosbetarot = cos(betarot)
          sinbetarot = sin(betarot)

          smallr1: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .TRUE.
            call prepare_refcoef(refl_var,r,0,aniso)

            !evaluate all integrals just for this radius
            ibesord = 1
            IendExy = compute_1val(func2Exwire,r,sz,zr,ibesord,wellbehaved(4))
            ibesord = 0
            IendEz = compute_1val(funcD0TM,r,sz,zr,ibesord,wellbehaved(5))
            ibesord = 1
            IendHxy = compute_1val(funcdHxwire,r,sz,zr,ibesord,wellbehaved(6))

          else !r larger than threshold radius for spline interpolation

            !get integral values for this radius by interpolation from precomputed integral values
            IendExy = splinterp_1val(refl_var,4,r)
            IendEz = splinterp_1val(refl_var,5,r)
            IendHxy = splinterp_1val(refl_var,6,r)
          endif smallr1

          !include rotation here
          !> Exrot etc. are field components in coordiante system wih axes parallel to source wire
          !> E(irec,#) are field values in coord. system with original x,y axes
          term1 = -sign * IendExy / dfourpi
          Exrot = term1 * cosbetarot
          Eyrot = term1 * sinbetarot
          fld(idx)%Ex(recidx) = fld(idx)%Ex(recidx) + Exrot*cosbetasrc - Eyrot*sinbetasrc
          fld(idx)%Ey(recidx) = fld(idx)%Ey(recidx) + Exrot*sinbetasrc + Eyrot*cosbetasrc

          term1 = sign * constEz * IendEz
          fld(idx)%Ez(recidx) = fld(idx)%Ez(recidx) + term1

          term1 = sign * IendHxy / dfourpi
          Hxrot = term1 * sinbetarot
          Hyrot = -term1 * cosbetarot
          fld(idx)%Hx(recidx) = fld(idx)%Hx(recidx) + Hxrot*cosbetasrc - Hyrot*sinbetasrc
          fld(idx)%Hy(recidx) = fld(idx)%Hy(recidx) + Hxrot*sinbetasrc + Hyrot*cosbetasrc

        endif chkr0

        !VTI: epsv derivatives
        dvert: if(with_dvert) then
          chkr0v: if(r .EQ. 0._real64) then

            !cycling was already done above, no need to query again here is src and rec are at exactly the same position
            !if(sz_eq_zr) cycle

            !only need Bessel function of order 0, since Bessel function of order 1 is zero for r=0
            ibesord = 0
            IendEz = compute_1valr0(funcD0TMv)

            term1 = sign * constEz * IendEz
            fldv(idx)%Ez(recidx) = fldv(idx)%Ez(recidx) + term1

          else !r is not zero

            smallr1v: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

              !reflection coeff. for this radius
              refl_var%refcoef_changed = .TRUE.
              call prepare_refcoef(refl_var,r,0,aniso)

              !evaluate all integrals just for this radius
              ibesord = 1
              IendExy = compute_1val(func2Exwirev,r,sz,zr,ibesord,wellbehaved(8))
              ibesord = 0
              IendEz = compute_1val(funcD0TMv,r,sz,zr,ibesord,wellbehaved(9))
              ibesord = 1
              IendHxy = compute_1val(funcdHxwirev,r,sz,zr,ibesord,wellbehaved(10))

            else !r larger than threshold radius for spline interpolation

              !get integral values for this radius by interpolation from precomputed integral values
              IendExy = splinterp_1val(refl_var,8,r)
              IendEz = splinterp_1val(refl_var,9,r)
              IendHxy = splinterp_1val(refl_var,10,r)
            endif smallr1v

            !include rotation here
            !> Exrot etc. are field components in coordiante system wih axes parallel to source wire
            !> E(irec,#) are field values in coord. system with original x,y axes
            term1 = -sign * IendExy / dfourpi
            Exrot = term1 * cosbetarot
            Eyrot = term1 * sinbetarot
            fldv(idx)%Ex(recidx) = fldv(idx)%Ex(recidx) + Exrot*cosbetasrc - Eyrot*sinbetasrc
            fldv(idx)%Ey(recidx) = fldv(idx)%Ey(recidx) + Exrot*sinbetasrc + Eyrot*cosbetasrc

            term1 = sign * constEz * IendEz
            fldv(Idx)%Ez(recidx) = fldv(idx)%Ez(recidx) + term1

            term1 = sign * IendHxy / dfourpi
            Hxrot = term1 * sinbetarot
            Hyrot = -term1 * cosbetarot
            fldv(idx)%Hx(recidx) = fldv(idx)%Hx(recidx) + Hxrot*cosbetasrc - Hyrot*sinbetasrc
            fldv(idx)%Hy(recidx) = fldv(idx)%Hy(recidx) + Hxrot*sinbetasrc + Hyrot*cosbetasrc

          endif chkr0v
        endif dvert

      enddo receivers1 !receivers
    enddo wire_endpoints
  enddo wires  !wires


  !special contribution to derivative of Ez in receiver layer
  deriv_ilayrec: if(ilay .EQ. ilayrec) then
    allocate(Ewirerec(src%nwire),stat=ierr)
    if(ierr.NE. 0) call alloc_error(pid,'interp_intvals_wire','WEwirerec',ierr)

    !there is a factor epsv in the term before the Ez integral
    !-> for isotropic case, epsv becomes eps and contribution is added to E
    !-> for VTI case, derivative of this term only exists for epsv -> add to Ev, nothing added to E
    if(with_dvert) then
      do iwire=1,src%nwire
        Ewirerec(iwire)%Ez => fldv(iwire)%Ez
      enddo
    else
      do iwire=1,src%nwire
        Ewirerec(iwire)%Ez => fld(iwire)%Ez
      enddo
    endif

    !indicators for fast Hankel transform or adaptive integration
    wellbehaved(7) = wellbehaved(5)

    !constant for Ez
    constEz = - dci / (omeps_recv*dfourpi*epsv(ilayrec))

    !loop over wires
    wiresilayrec: do isrc = refl_var%isrcstart,refl_var%isrcend

      !get constant geometry factor: wire element length
      idx = refl_var%isrcperz(isrc)

      !loop over wire end points
      wire_endpointsilayrec: do iep = 1,2
        xs = src%wire(idx)%endpos(1,iep)
        ys = src%wire(idx)%endpos(2,iep)

        sign = signs(iep)

        receivers1ilayrec: do irec=refl_var%irecstart,refl_var%irecend
          recidx = refl_var%irecperzExy(irec)

          y = bgdat%Exypos(recidx,2) - ys
          x = bgdat%Exypos(recidx,1) - xs
          r = sqrt(x**2 + y**2)

          if(r.EQ.0._real64) then
            !skip integration right at source point - check this!
            if(sz_eq_zr) cycle
            IendEz = compute_1valr0(funcD0TMfwd)
          else
            if(r.lt.rsplmin) then
              !reflection coeff. for this radius
              call prepare_refcoef(refl_var,r,ved,aniso) !use "ved" to compute TM refl. coeff. only
              ibesord = 0
              IendEz = compute_1val(funcD0TMfwd,r,sz,zr,ibesord,wellbehaved(7))
            else
              IendEz = splinterp_1val(refl_var,7,r)
            endif
          endif

          term1 = sign * constEz * IendEz
          Ewirerec(idx)%Ez(recidx) = Ewirerec(idx)%Ez(recidx) + term1

        enddo receivers1ilayrec !receivers
      enddo wire_endpointsilayrec
    enddo wiresilayrec  !wires

    do iwire=1,src%nwire
      nullify(Ewirerec(iwire)%Ez)
    enddo
    deallocate(Ewirerec, stat=ierr)
  endif deriv_ilayrec

endsubroutine interp_intvals_wire_allcomp


!------------------------------------------------------------
!>  1D EM subroutine interp_intvals_wire_Exy
!
!>  get interpolated field values at receiver locations, horizontal wire source
!>    Ex and/or Ey only
!
!>  Rita Streich 2009-2011
!------------------------------------------------------------
subroutine interp_intvals_wire_Exy(refl_var,bgdat,fld,src,sz,zr, func1Exwire,func2Exwire,ilay, func2Exwirev,fldv)

  implicit none

  !external variables
  type(refl_struct) :: refl_var   !everything needed throughout 1D computations
  type(backgrounddata) :: bgdat      !coordinate vectors and final output EM fields
  type(receiverdata),dimension(:) :: fld        !field or derivative vectors for each wire, all components
  type(sorec),intent( in ) :: src        !source definition (we need the currents here)
  real(kind=real64),intent( in ) :: sz,zr      !source and receiver depth
  complex(kind=real64),external :: func1Exwire,func2Exwire
  integer(kind=int32),intent( in ) :: ilay       !layer index for derivatives - leave at zero for forward modeling
  complex(kind=real64),external,optional :: func2Exwirev !integral derivatives for epsv
  type(receiverdata),dimension(:),optional :: fldv         !field or derivative vectors for each wire, all components

  !internal variables
  integer(kind=int32) :: isrc    !source element counter
  real(kind=real64) :: x,y,r   !temp source-receiver distances
  integer(kind=int32) :: irec    !receiver counter
  real(kind=real64) :: beta,betarot      !temp angles

  complex(kind=real64) :: IExy          !interpolated integral values for integrals along wire
  complex(kind=real64) :: IendExy       !interpolated integral values for end points

  logical,dimension(nintWiredvti) :: wellbehaved   !indicates if Hankel integration can be used
  logical :: sz_eq_zr      !indicates if source and reveicer are at the same depth
  integer(kind=int32),dimension(NREL) :: ibesord   !bessel function order (integer "array")

  complex(kind=real64) :: Exrot,Eyrot  !temp field values in cylindrical coordinates
  real(kind=real64) :: cosbetasrc,sinbetasrc !cos(betasrc) and sin(betasrc), precompute for efficiency
  real(kind=real64) :: cosbetarot,sinbetarot !cos(betarot) and sin(betarot), precompute for efficiency
  real(kind=real64) :: const    !geometric constant
  integer(kind=int32) :: idx      !source element index
  integer(kind=int32) :: recidx   !receiver index
  integer(kind=int32) :: ielemall,ielem       !wire element counters
  integer(kind=int32) :: iep      !wire end point counter
  real(kind=real64) :: xs,ys    !wire end point coordinates
  complex(kind=real64) :: term1    !temp result for end point field values
  !signs for wire end points (predefine for efficiency, need just elements 1 and 4)
  real(kind=real64),dimension(1:2),parameter :: signs = (/-1._real64,1._real64/)
  real(kind=real64) :: sign     !sign for end point field values
  !flag for computing epsv derivatives, not needed for forward computation, so "aniso" value cannot be used here
  logical :: with_dvert


  !indicators for fast Hankel transform or adaptive integration
  wellbehaved = .TRUE.
  sz_eq_zr = .FALSE.
  if(sz.EQ.zr) then
    wellbehaved(1:2) = .FALSE.
    wellbehaved(4:5) = .FALSE.
    sz_eq_zr = .TRUE.
  endif

  if(present(func2Exwirev)) then
    with_dvert = .TRUE.
    if(sz_eq_zr) then
      wellbehaved(8:9) = .FALSE.
    endif
  else
    with_dvert = .FALSE.
  endif


  !loop over wires
  !TE only --> no contribution to epsv derivatives
  ielemall = 0
  wires: do isrc = refl_var%isrcstart,refl_var%isrcend

    !get constant geometry factor: wire element length
    idx = refl_var%isrcperz(isrc)
    const = src%wire(idx)%dlw / dfourpi

    cosbetasrc = cos(refl_var%betasrc(isrc))
    sinbetasrc = sin(refl_var%betasrc(isrc))

    !loop over wire elements
    wire_elements: do ielem = 1,src%nelem(idx)
      ielemall = ielemall + 1

      !same positions for Ex and Ey
      exy_equalpos: if(bgdat%nExy.gt.0) then

      receivers: do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzExy(irec)

        y = bgdat%Exypos(recidx,2) - refl_var%ys(ielemall)
        x = bgdat%Exypos(recidx,1) - refl_var%xs(ielemall)
        r = sqrt(x**2 + y**2)

        beta = atan2(y,x)
        betarot = beta - refl_var%betasrc(isrc)

        r_is_zero: if(r.EQ.0._real64) then

          !quick & dirty: skip the point if receiver is right at source point
          if(sz_eq_zr) then
            if(refl_var%infolevel.ge.output_more) &
              write(*,'(a)') 'WARNING: cannot handle a receiver at exactly the same location as a wire element yet, '// &
                'ignoring this contribution!'
            cycle
          endif

          !--> evaluate zero-order integrals only
          ibesord = 0
          IExy = compute_1valr0(func1Exwire)

        else !r is not zero
          smallr: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .TRUE.
            call prepare_refcoef(refl_var,r,0,aniso)

            !evaluate all integrals just for this radius
            ibesord = 0
            IExy = compute_1val(func1Exwire,r,sz,zr,ibesord,wellbehaved(1))
          else !r larger than threshold radius for spline interpolation
            !get integral values for this radius by interpolation from precomputed integral values
            IExy = splinterp_1val(refl_var,1,r)
          endif smallr
        endif r_is_zero

        !field values in rotated coordinate system
        Exrot = IExy * const
 
        !field values rotated back to original x-y coordinate system
        !don't take complex conjugate here, but at the very end!
        fld(idx)%Ex(recidx) = fld(idx)%Ex(recidx) + Exrot*cosbetasrc
        fld(idx)%Ey(recidx) = fld(idx)%Ey(recidx) + Exrot*sinbetasrc
      enddo receivers

      else

        have_ex: if(bgdat%nEx .gt. 0) then

      receiversex: do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzExy(irec)

        y = bgdat%Expos(recidx,2) - refl_var%ys(ielemall)
        x = bgdat%Expos(recidx,1) - refl_var%xs(ielemall)
        r = sqrt(x**2 + y**2)

        beta = atan2(y,x)
        betarot = beta - refl_var%betasrc(isrc)

        r_is_zeroex: if(r.EQ.0._real64) then

          !quick & dirty: skip the point if receiver is right at source point
          if(sz_eq_zr) then
            if(refl_var%infolevel.ge.output_more) &
              write(*,'(a)') 'WARNING: cannot handle a receiver at exactly the same location as a wire element yet, '// &
                'ignoring this contribution!'
            cycle
          endif

          !--> evaluate zero-order integrals only
          ibesord = 0
          IExy = compute_1valr0(func1Exwire)

        else !r is not zero
          smallrex: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .TRUE.
            call prepare_refcoef(refl_var,r,0,aniso)

            !evaluate all integrals just for this radius
            ibesord = 0
            IExy = compute_1val(func1Exwire,r,sz,zr,ibesord,wellbehaved(1))
          else !r larger than threshold radius for spline interpolation
            !get integral values for this radius by interpolation from precomputed integral values
            IExy = splinterp_1val(refl_var,1,r)
          endif smallrex
        endif r_is_zeroex

        !field values in rotated coordinate system
        Exrot = IExy * const
 
        !field values rotated back to original x-y coordinate system
        !don't take complex conjugate here, but at the very end!
        fld(idx)%Ex(recidx) = fld(idx)%Ex(recidx) + Exrot*cosbetasrc
      enddo receiversex

        endif have_ex

        have_ey: if(bgdat%nEy .gt. 0) then

      receiversey: do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzExy(irec)

        y = bgdat%Eypos(recidx,2) - refl_var%ys(ielemall)
        x = bgdat%Eypos(recidx,1) - refl_var%xs(ielemall)
        r = sqrt(x**2 + y**2)

        beta = atan2(y,x)
        betarot = beta - refl_var%betasrc(isrc)

        r_is_zeroey: if(r.EQ.0._real64) then

          !quick & dirty: skip the point if receiver is right at source point
          if(sz_eq_zr) then
            if(refl_var%infolevel.ge.output_more) &
              write(*,'(a)') 'WARNING: cannot handle a receiver at exactly the same location as a wire element yet, '// &
                'ignoring this contribution!'
            cycle
          endif

          !--> evaluate zero-order integrals only
          ibesord = 0
          IExy = compute_1valr0(func1Exwire)

        else !r is not zero
          smallrey: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .TRUE.
            call prepare_refcoef(refl_var,r,0,aniso)

            !evaluate all integrals just for this radius
            ibesord = 0
            IExy = compute_1val(func1Exwire,r,sz,zr,ibesord,wellbehaved(1))
          else !r larger than threshold radius for spline interpolation
            !get integral values for this radius by interpolation from precomputed integral values
            IExy = splinterp_1val(refl_var,1,r)
          endif smallrey
        endif r_is_zeroey

        !field values in rotated coordinate system
        Exrot = IExy * const
 
        !field values rotated back to original x-y coordinate system
        !don't take complex conjugate here, but at the very end!
        fld(idx)%Ey(recidx) = fld(idx)%Ey(recidx) + Exrot*sinbetasrc
      enddo receiversey !iy

        endif have_ey
      endif exy_equalpos
    enddo wire_elements


    !loop over wire end points
    wire_endpoints: do iep = 1,2
      xs = src%wire(idx)%endpos(1,iep)
      ys = src%wire(idx)%endpos(2,iep)

      sign = signs(iep)

      !same positions for Ex and Ey
      exy_equalpos_ep: if(bgdat%nExy.gt.0) then

      receivers1: do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzExy(irec)

        y = bgdat%Exypos(recidx,2) - ys
        x = bgdat%Exypos(recidx,1) - xs
        r = sqrt(x**2 + y**2)

        !distinguish case r=0 since here we need only 1 integral
        chkr0: if(r .EQ. 0._real64) then

          !quick & dirty: skip the point if receiver is right at source point
          if(sz_eq_zr) then
            if(refl_var%infolevel.ge.output_more) &
              write(*,'(a)') 'WARNING: cannot handle a receiver at exactly the same location as a grounding point yet, '// &
                'ignoring this contribution!'
            cycle
          endif
        else !r is not zero

          beta = atan2(y,x)
          betarot = beta - refl_var%betasrc(isrc)
          cosbetarot = cos(betarot)
          sinbetarot = sin(betarot)

          smallr1: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .TRUE.
            call prepare_refcoef(refl_var,r,0,aniso)

            !evaluate all integrals just for this radius
            ibesord = 1
            IendExy = compute_1val(func2Exwire,r,sz,zr,ibesord,wellbehaved(4))
          else !r larger than threshold radius for spline interpolation
            !get integral values for this radius by interpolation from precomputed integral values
            IendExy = splinterp_1val(refl_var,4,r)
          endif smallr1

          !include rotation here
          !> Exrot etc. are field components in coordiante system wih axes parallel to source wire
          !> E(irec,#) are field values in coord. system with original x,y axes
          term1 = -sign * IendExy / dfourpi
          Exrot = term1 * cosbetarot
          Eyrot = term1 * sinbetarot
          fld(idx)%Ex(recidx) = fld(idx)%Ex(recidx) + Exrot*cosbetasrc - Eyrot*sinbetasrc
          fld(idx)%Ey(recidx) = fld(idx)%Ey(recidx) + Exrot*sinbetasrc + Eyrot*cosbetasrc
        endif chkr0

        !VTI: epsv derivatives
        dvert: if(with_dvert) then
          chkr0v: if(r .NE. 0._real64) then
            smallr1v: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

              !reflection coeff. for this radius
              refl_var%refcoef_changed = .TRUE.
              call prepare_refcoef(refl_var,r,0,aniso)

              !evaluate all integrals just for this radius
              ibesord = 1
              IendExy = compute_1val(func2Exwirev,r,sz,zr,ibesord,wellbehaved(8))
            else !r larger than threshold radius for spline interpolation

              !get integral values for this radius by interpolation from precomputed integral values
              IendExy = splinterp_1val(refl_var,8,r)
            endif smallr1v

            !include rotation here
            !> Exrot etc. are field components in coordiante system wih axes parallel to source wire
            !> E(irec,#) are field values in coord. system with original x,y axes
            term1 = -sign * IendExy / dfourpi
            Exrot = term1 * cosbetarot
            Eyrot = term1 * sinbetarot
            fldv(idx)%Ex(recidx) = fldv(idx)%Ex(recidx) + Exrot*cosbetasrc - Eyrot*sinbetasrc
            fldv(idx)%Ey(recidx) = fldv(idx)%Ey(recidx) + Exrot*sinbetasrc + Eyrot*cosbetasrc
          endif chkr0v
        endif dvert

      enddo receivers1 !receivers

      else

        have_ex_ep: if(bgdat%nEx .gt. 0) then

      receivers1ex: do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzExy(irec)

        y = bgdat%Expos(recidx,2) - ys
        x = bgdat%Expos(recidx,1) - xs
        r = sqrt(x**2 + y**2)

        !distinguish case r=0 since here we need only 1 integral
        chkr0ex: if(r .EQ. 0._real64) then

          !quick & dirty: skip the point if receiver is right at source point
          if(sz_eq_zr) then
            if(refl_var%infolevel.ge.output_more) &
              write(*,'(a)') 'WARNING: cannot handle a receiver at exactly the same location as a grounding point yet, '// &
                'ignoring this contribution!'
            cycle
          endif
        else !r is not zero

          beta = atan2(y,x)
          betarot = beta - refl_var%betasrc(isrc)
          cosbetarot = cos(betarot)
          sinbetarot = sin(betarot)

          smallr1ex: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .TRUE.
            call prepare_refcoef(refl_var,r,0,aniso)

            !evaluate all integrals just for this radius
            ibesord = 1
            IendExy = compute_1val(func2Exwire,r,sz,zr,ibesord,wellbehaved(4))
          else !r larger than threshold radius for spline interpolation
            !get integral values for this radius by interpolation from precomputed integral values
            IendExy = splinterp_1val(refl_var,4,r)
          endif smallr1ex

          !include rotation here
          !> Exrot etc. are field components in coordiante system wih axes parallel to source wire
          !> E(irec,#) are field values in coord. system with original x,y axes
          term1 = -sign * IendExy / dfourpi
          Exrot = term1 * cosbetarot
          Eyrot = term1 * sinbetarot
          fld(idx)%Ex(recidx) = fld(idx)%Ex(recidx) + Exrot*cosbetasrc - Eyrot*sinbetasrc
        endif chkr0ex

        !VTI: epsv derivatives
        dvertex: if(with_dvert) then
          chkr0vex: if(r .NE. 0._real64) then
            smallr1vex: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

              !reflection coeff. for this radius
              refl_var%refcoef_changed = .TRUE.
              call prepare_refcoef(refl_var,r,0,aniso)

              !evaluate all integrals just for this radius
              ibesord = 1
              IendExy = compute_1val(func2Exwirev,r,sz,zr,ibesord,wellbehaved(8))
            else !r larger than threshold radius for spline interpolation

              !get integral values for this radius by interpolation from precomputed integral values
              IendExy = splinterp_1val(refl_var,8,r)
            endif smallr1vex

            !include rotation here
            !> Exrot etc. are field components in coordiante system wih axes parallel to source wire
            !> E(irec,#) are field values in coord. system with original x,y axes
            term1 = -sign * IendExy / dfourpi
            Exrot = term1 * cosbetarot
            Eyrot = term1 * sinbetarot
            fldv(idx)%Ex(recidx) = fldv(idx)%Ex(recidx) + Exrot*cosbetasrc - Eyrot*sinbetasrc
          endif chkr0vex
        endif dvertex
      enddo receivers1ex !receivers

        endif have_ex_ep

        have_ey_ep: if(bgdat%nEy .gt. 0) then

      receivers1ey: do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzExy(irec)

        y = bgdat%Eypos(recidx,2) - ys
        x = bgdat%Eypos(recidx,1) - xs
        r = sqrt(x**2 + y**2)

        !distinguish case r=0 since here we need only 1 integral
        chkr0ey: if(r .EQ. 0._real64) then

          !quick & dirty: skip the point if receiver is right at source point
          if(sz_eq_zr) then
            if(refl_var%infolevel.ge.output_more) &
              write(*,'(a)') 'WARNING: cannot handle a receiver at exactly the same location as a grounding point yet, '// &
                'ignoring this contribution!'
            cycle
          endif
        else !r is not zero

          beta = atan2(y,x)
          betarot = beta - refl_var%betasrc(isrc)
          cosbetarot = cos(betarot)
          sinbetarot = sin(betarot)

          smallr1ey: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .TRUE.
            call prepare_refcoef(refl_var,r,0,aniso)

            !evaluate all integrals just for this radius
            ibesord = 1
            IendExy = compute_1val(func2Exwire,r,sz,zr,ibesord,wellbehaved(4))
          else !r larger than threshold radius for spline interpolation
            !get integral values for this radius by interpolation from precomputed integral values
            IendExy = splinterp_1val(refl_var,4,r)
          endif smallr1ey

          !include rotation here
          !> Exrot etc. are field components in coordiante system wih axes parallel to source wire
          !> E(irec,#) are field values in coord. system with original x,y axes
          term1 = -sign * IendExy / dfourpi
          Exrot = term1 * cosbetarot
          Eyrot = term1 * sinbetarot
          fld(idx)%Ey(recidx) = fld(idx)%Ey(recidx) + Exrot*sinbetasrc + Eyrot*cosbetasrc
        endif chkr0ey

        !VTI: epsv derivatives
        dvertey: if(with_dvert) then
          chkr0vey: if(r .NE. 0._real64) then
            smallr1vey: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

              !reflection coeff. for this radius
              refl_var%refcoef_changed = .TRUE.
              call prepare_refcoef(refl_var,r,0,aniso)

              !evaluate all integrals just for this radius
              ibesord = 1
              IendExy = compute_1val(func2Exwirev,r,sz,zr,ibesord,wellbehaved(8))
            else !r larger than threshold radius for spline interpolation

              !get integral values for this radius by interpolation from precomputed integral values
              IendExy = splinterp_1val(refl_var,8,r)
            endif smallr1vey

            !include rotation here
            !> Exrot etc. are field components in coordiante system wih axes parallel to source wire
            !> E(irec,#) are field values in coord. system with original x,y axes
            term1 = -sign * IendExy / dfourpi
            Exrot = term1 * cosbetarot
            Eyrot = term1 * sinbetarot
            fldv(idx)%Ey(recidx) = fldv(idx)%Ey(recidx) + Exrot*sinbetasrc + Eyrot*cosbetasrc
          endif chkr0vey
        endif dvertey
      enddo receivers1ey !receivers

        endif have_ey_ep
      endif exy_equalpos_ep

    enddo wire_endpoints
  enddo wires  !wires

endsubroutine interp_intvals_wire_Exy


!------------------------------------------------------------
!>  1D EM subroutine interp_intvals_wire_Ez
!
!>  get interpolated field values at receiver locations, horizontal wire source,
!>    Ez only
!
!>  Rita Streich 2009-2011
!------------------------------------------------------------
subroutine interp_intvals_wire_Ez(refl_var,bgdat,fld,src,sz,zr,omeps_recv, funcD0TM,ilay,funcD0TMfwd, funcD0TMv,fldv)

  implicit none

  !external variables
  type(refl_struct) :: refl_var   !everything needed throughout 1D computations
  type(backgrounddata) :: bgdat      !coordinate vectors and final output EM fields
  type(receiverdata),dimension(:) :: fld        !field or derivative vectors for each wire, all components
  type(sorec),intent( in ) :: src        !source definition (we need the currents here)
  real(kind=real64),intent( in ) :: sz,zr      !source and receiver depth
  complex(kind=real64),intent( in ) :: omeps_recv  !omega * eps in receiver layer
  complex(kind=real64),external :: funcD0TM
  integer(kind=int32),intent( in ) :: ilay       !layer index for derivatives - leave at zero for forward modeling
  complex(kind=real64),external,optional :: funcD0TMfwd  !function only needed for derivatives in receiver layer
  complex(kind=real64),external,optional :: funcD0TMv    !integral derivatives for epsv
  type(receiverdata),dimension(:),optional :: fldv       !field or derivative vectors for each wire, all components

  !internal variables
  integer(kind=int32) :: isrc    !source element counter
  real(kind=real64) :: x,y,r   !temp source-receiver distances
  integer(kind=int32) :: irec    !receiver counter
  real(kind=real64) :: beta,betarot      !temp angles

  complex(kind=real64) :: IendEz !interpolated integral values for end points

  logical,dimension(nintWiredvti) :: wellbehaved   !indicates if Hankel integration can be used
  logical :: sz_eq_zr      !indicates if source and reveicer are at the same depth
  integer(kind=int32),dimension(NREL) :: ibesord   !bessel function order (integer "array")

  real(kind=real64) :: cosbetasrc,sinbetasrc !cos(betasrc) and sin(betasrc), precompute for efficiency
  real(kind=real64) :: cosbetarot,sinbetarot !cos(betarot) and sin(betarot), precompute for efficiency
  real(kind=real64) :: const    !geometric constant
  complex(kind=real64) :: constEz  !geometric constant for Ez
  integer(kind=int32) :: idx      !source element index
  integer(kind=int32) :: recidx   !receiver index
  integer(kind=int32) :: ielemall  !wire element counters
  integer(kind=int32) :: iep      !wire end point counter
  real(kind=real64) :: xs,ys    !wire end point coordinates
  complex(kind=real64) :: term1    !temp result for end point field values
  !signs for wire end points (predefine for efficiency, need just elements 1 and 4)
  real(kind=real64),dimension(1:2),parameter :: signs = (/-1._real64,1._real64/)
  real(kind=real64) :: sign     !sign for end point field values
  !flag for computing epsv derivatives, not needed for forward computation, so "aniso" value cannot be used here
  logical :: with_dvert
  type(receiverdata),dimension(:),pointer :: Ewirerec  !points to Ez for isotropic and Ezv for VTI case
  integer(kind=int32) :: ierr     !error index
  integer(kind=int32) :: iwire    !wire counter


  !indicators for fast Hankel transform or adaptive integration
  wellbehaved = .TRUE.
  sz_eq_zr = .FALSE.
  if(sz.EQ.zr) then
    wellbehaved(1:2) = .FALSE.
    wellbehaved(4:5) = .FALSE.
    sz_eq_zr = .TRUE.
  endif

  if(present(funcD0TMv)) then
    with_dvert = .TRUE.
    if(sz_eq_zr) then
      wellbehaved(8:9) = .FALSE.
    endif
  else
    with_dvert = .FALSE.
  endif

  !constant for Ez
  constEz = dci / (omeps_recv*dfourpi)


  !loop over wires
  !TE only --> no contribution to epsv derivatives
  ielemall = 0
  wires: do isrc = refl_var%isrcstart,refl_var%isrcend

    !get constant geometry factor: wire element length
    idx = refl_var%isrcperz(isrc)
    const = src%wire(idx)%dlw / dfourpi

    cosbetasrc = cos(refl_var%betasrc(isrc))
    sinbetasrc = sin(refl_var%betasrc(isrc))

    !loop over wire end points
    wire_endpoints: do iep = 1,2
      xs = src%wire(idx)%endpos(1,iep)
      ys = src%wire(idx)%endpos(2,iep)

      sign = signs(iep)

      receivers1: do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzEz(irec)

        y = bgdat%Ezpos(recidx,2) - ys
        x = bgdat%Ezpos(recidx,1) - xs
        r = sqrt(x**2 + y**2)

        !distinguish case r=0 since here we need only 1 integral
        chkr0: if(r .EQ. 0._real64) then

          !quick & dirty: skip the point if receiver is right at source point
          if(sz_eq_zr) then
            if(refl_var%infolevel.ge.output_more) &
              write(*,'(a)') 'WARNING: cannot handle a receiver at exactly the same location as a grounding point yet, '// &
                'ignoring this contribution!'
            cycle
          endif

          !only need Bessel function of order 0, since Bessel function of order 1 is zero for r=0
          ibesord = 0
          IendEz = compute_1valr0(funcD0TM)

          term1 = sign * constEz * IendEz
          fld(idx)%Ez(recidx) = fld(idx)%Ez(recidx) + term1

        else !r is not zero

          beta = atan2(y,x)
          betarot = beta - refl_var%betasrc(isrc)
          cosbetarot = cos(betarot)
          sinbetarot = sin(betarot)

          smallr1: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .TRUE.
            call prepare_refcoef(refl_var,r,0,aniso)

            !evaluate all integrals just for this radius
            ibesord = 0
            IendEz = compute_1val(funcD0TM,r,sz,zr,ibesord,wellbehaved(5))
          else !r larger than threshold radius for spline interpolation
            !get integral values for this radius by interpolation from precomputed integral values
            IendEz = splinterp_1val(refl_var,5,r)
          endif smallr1

          !include rotation here
          !> E(irec,#) are field values in coord. system with original x,y axes
          term1 = sign * constEz * IendEz
          fld(idx)%Ez(recidx) = fld(idx)%Ez(recidx) + term1
        endif chkr0

        !VTI: epsv derivatives
        dvert: if(with_dvert) then
          chkr0v: if(r .EQ. 0._real64) then
            !cycling was already done above, no need to query again here is src and rec are at exactly the same position
            !if(sz_eq_zr) cycle

            !only need Bessel function of order 0, since Bessel function of order 1 is zero for r=0
            ibesord = 0
            IendEz = compute_1valr0(funcD0TMv)

            term1 = sign * constEz * IendEz
            fldv(idx)%Ez(recidx) = fldv(idx)%Ez(recidx) + term1
          else !r is not zero
            smallr1v: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

              !reflection coeff. for this radius
              refl_var%refcoef_changed = .TRUE.
              call prepare_refcoef(refl_var,r,0,aniso)

              !evaluate all integrals just for this radius
              ibesord = 0
              IendEz = compute_1val(funcD0TMv,r,sz,zr,ibesord,wellbehaved(9))
            else !r larger than threshold radius for spline interpolation
              !get integral values for this radius by interpolation from precomputed integral values
              IendEz = splinterp_1val(refl_var,9,r)
            endif smallr1v

            !include rotation here
            !> E(irec,#) are field values in coord. system with original x,y axes
            term1 = sign * constEz * IendEz
            fldv(Idx)%Ez(recidx) = fldv(idx)%Ez(recidx) + term1
          endif chkr0v
        endif dvert

      enddo receivers1 !receivers
    enddo wire_endpoints
  enddo wires  !wires


  !special contribution to derivative of Ez in receiver layer
  deriv_ilayrec: if(ilay .EQ. ilayrec) then
    allocate(Ewirerec(src%nwire),stat=ierr)
    if(ierr.NE. 0) call alloc_error(pid,'interp_intvals_wire','WEwirerec',ierr)

    !there is a factor epsv in the term before the Ez integral
    !-> for isotropic case, epsv becomes eps and contribution is added to E
    !-> for VTI case, derivative of this term only exists for epsv -> add to Ev, nothing added to E
    if(with_dvert) then
      do iwire=1,src%nwire
        Ewirerec(iwire)%Ez => fldv(iwire)%Ez
      enddo
    else
      do iwire=1,src%nwire
        Ewirerec(iwire)%Ez => fld(iwire)%Ez
      enddo
    endif

    !indicators for fast Hankel transform or adaptive integration
    wellbehaved(7) = wellbehaved(5)

    !constant for Ez
    constEz = - dci / (omeps_recv*dfourpi*epsv(ilayrec))

    !loop over wires
    wiresilayrec: do isrc = refl_var%isrcstart,refl_var%isrcend

      !get constant geometry factor: wire element length
      idx = refl_var%isrcperz(isrc)

      !loop over wire end points
      wire_endpointsilayrec: do iep = 1,2
        xs = src%wire(idx)%endpos(1,iep)
        ys = src%wire(idx)%endpos(2,iep)

        sign = signs(iep)

        receivers1ilayrec: do irec=refl_var%irecstart,refl_var%irecend
          recidx = refl_var%irecperzEz(irec)

          y = bgdat%Ezpos(recidx,2) - ys
          x = bgdat%Ezpos(recidx,1) - xs
          r = sqrt(x**2 + y**2)

          if(r.EQ.0._real64) then
            !skip integration right at source point - check this!
            if(sz_eq_zr) cycle
            IendEz = compute_1valr0(funcD0TMfwd)
          else
            if(r.lt.rsplmin) then
              !reflection coeff. for this radius
              call prepare_refcoef(refl_var,r,ved,aniso) !use "ved" to compute TM refl. coeff. only
              ibesord = 0
              IendEz = compute_1val(funcD0TMfwd,r,sz,zr,ibesord,wellbehaved(7))
            else
              IendEz = splinterp_1val(refl_var,7,r)
            endif
          endif

          term1 = sign * constEz * IendEz
          Ewirerec(idx)%Ez(recidx) = Ewirerec(idx)%Ez(recidx) + term1

        enddo receivers1ilayrec !receivers
      enddo wire_endpointsilayrec
    enddo wiresilayrec  !wires

    do iwire=1,src%nwire
      nullify(Ewirerec(iwire)%Ez)
    enddo
    deallocate(Ewirerec, stat=ierr)
  endif deriv_ilayrec

endsubroutine interp_intvals_wire_Ez


!------------------------------------------------------------
!>  1D EM subroutine interp_intvals_wire_Hxy
!
!>  get interpolated field values at receiver locations, horizontal wire source,
!>    Hx and/or Hy only
!
!>  Rita Streich 2009-2011
!------------------------------------------------------------
subroutine interp_intvals_wire_Hxy(refl_var,bgdat,fld,src,sz,zr, funcD0TE,funcdHxwire,ilay, funcdHxwirev,fldv)

  implicit none

  !external variables
  type(refl_struct) :: refl_var   !everything needed throughout 1D computations
  type(backgrounddata) :: bgdat      !coordinate vectors and final output EM fields
  type(receiverdata),dimension(:) :: fld        !field or derivative vectors for each wire, all components
  type(sorec),intent( in ) :: src        !source definition (we need the currents here)
  real(kind=real64),intent( in ) :: sz,zr      !source and receiver depth
  complex(kind=real64),external :: funcD0TE,funcdHxwire
  integer(kind=int32),intent( in ) :: ilay       !layer index for derivatives - leave at zero for forward modeling
  complex(kind=real64),external,optional :: funcdHxwirev !integral derivatives for epsv
  type(receiverdata),dimension(:),optional :: fldv         !field or derivative vectors for each wire, all components

  !internal variables
  integer(kind=int32) :: isrc    !source element counter
  real(kind=real64) :: x,y,r   !temp source-receiver distances
  integer(kind=int32) :: irec    !receiver counter
  real(kind=real64) :: beta,betarot      !temp angles

  complex(kind=real64) :: IHyx          !interpolated integral values for integrals along wire
  complex(kind=real64) :: IendHxy       !interpolated integral values for end points

  logical,dimension(nintWiredvti) :: wellbehaved   !indicates if Hankel integration can be used
  logical :: sz_eq_zr      !indicates if source and reveicer are at the same depth
  integer(kind=int32),dimension(NREL) :: ibesord   !bessel function order (integer "array")

  complex(kind=real64) :: Hxrot,Hyrot           !temp field values in cylindrical coordinates
  real(kind=real64) :: cosbetasrc,sinbetasrc !cos(betasrc) and sin(betasrc), precompute for efficiency
  real(kind=real64) :: cosbetarot,sinbetarot !cos(betarot) and sin(betarot), precompute for efficiency
  real(kind=real64) :: const    !geometric constant
  integer(kind=int32) :: idx      !source element index
  integer(kind=int32) :: recidx   !receiver index
  integer(kind=int32) :: ielemall,ielem       !wire element counters
  integer(kind=int32) :: iep      !wire end point counter
  real(kind=real64) :: xs,ys    !wire end point coordinates
  complex(kind=real64) :: term1    !temp result for end point field values
  !signs for wire end points (predefine for efficiency, need just elements 1 and 4)
  real(kind=real64),dimension(1:2),parameter :: signs = (/-1._real64,1._real64/)
  real(kind=real64) :: sign     !sign for end point field values
  !flag for computing epsv derivatives, not needed for forward computation, so "aniso" value cannot be used here
  logical :: with_dvert


  !indicators for fast Hankel transform or adaptive integration
  wellbehaved = .TRUE.
  sz_eq_zr = .FALSE.
  if(sz.EQ.zr) then
    wellbehaved(1:2) = .FALSE.
    wellbehaved(4:5) = .FALSE.
    sz_eq_zr = .TRUE.
  endif

  if(present(funcdHxwirev)) then
    with_dvert = .TRUE.
    if(sz_eq_zr) then
      wellbehaved(8:9) = .FALSE.
    endif
  else
    with_dvert = .FALSE.
  endif

  !loop over wires
  !TE only --> no contribution to epsv derivatives
  ielemall = 0
  wires: do isrc = refl_var%isrcstart,refl_var%isrcend

    !get constant geometry factor: wire element length
    idx = refl_var%isrcperz(isrc)
    const = src%wire(idx)%dlw / dfourpi

    cosbetasrc = cos(refl_var%betasrc(isrc))
    sinbetasrc = sin(refl_var%betasrc(isrc))

    !loop over wire elements
    wire_elements: do ielem = 1,src%nelem(idx)
      ielemall = ielemall + 1

      !same positions for Hx and Hy
      hxy_equalpos: if(bgdat%nExy.gt.0) then

      receivers: do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzHxy(irec)

        y = bgdat%Hxypos(recidx,2) - refl_var%ys(ielemall)
        x = bgdat%Hxypos(recidx,1) - refl_var%xs(ielemall)
        r = sqrt(x**2 + y**2)

        beta = atan2(y,x)
        betarot = beta - refl_var%betasrc(isrc)

        r_is_zero: if(r.EQ.0._real64) then
          !quick & dirty: skip the point if receiver is right at source point
          if(sz_eq_zr) then
            if(refl_var%infolevel.ge.output_more) &
              write(*,'(a)') 'WARNING: cannot handle a receiver at exactly the same location as a wire element yet, '// &
                'ignoring this contribution!'
            cycle
          endif

          !--> evaluate zero-order integrals only
          ibesord = 0
          IHyx = compute_1valr0(funcD0TE)
        else !r is not zero
          smallr: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .TRUE.
            call prepare_refcoef(refl_var,r,0,aniso)

            !evaluate all integrals just for this radius
            ibesord = 0
            IHyx = compute_1val(funcD0TE,r,sz,zr,ibesord,wellbehaved(2))
          else !r larger than threshold radius for spline interpolation
            !get integral values for this radius by interpolation from precomputed integral values
            IHyx = splinterp_1val(refl_var,2,r)
          endif smallr
        endif r_is_zero

        !field values in rotated coordinate system
        Hyrot = - IHyx * const
 
        !field values rotated back to original x-y coordinate system
        fld(idx)%Hx(recidx) = fld(idx)%Hx(recidx) - Hyrot*sinbetasrc
        fld(idx)%Hy(recidx) = fld(idx)%Hy(recidx) + Hyrot*cosbetasrc
      enddo receivers !iy

      else
        have_hx: if(bgdat%nHx .gt. 0) then

      receivershx: do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzHxy(irec)

        y = bgdat%Hxpos(recidx,2) - refl_var%ys(ielemall)
        x = bgdat%Hxpos(recidx,1) - refl_var%xs(ielemall)
        r = sqrt(x**2 + y**2)

        beta = atan2(y,x)
        betarot = beta - refl_var%betasrc(isrc)

        r_is_zerohx: if(r.EQ.0._real64) then
          !quick & dirty: skip the point if receiver is right at source point
          if(sz_eq_zr) then
            if(refl_var%infolevel.ge.output_more) &
              write(*,'(a)') 'WARNING: cannot handle a receiver at exactly the same location as a wire element yet, '// &
                'ignoring this contribution!'
            cycle
          endif

          !--> evaluate zero-order integrals only
          ibesord = 0
          IHyx = compute_1valr0(funcD0TE)
        else !r is not zero
          smallrhx: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .TRUE.
            call prepare_refcoef(refl_var,r,0,aniso)

            !evaluate all integrals just for this radius
            ibesord = 0
            IHyx = compute_1val(funcD0TE,r,sz,zr,ibesord,wellbehaved(2))
          else !r larger than threshold radius for spline interpolation
            !get integral values for this radius by interpolation from precomputed integral values
            IHyx = splinterp_1val(refl_var,2,r)
          endif smallrhx
        endif r_is_zerohx

        !field values in rotated coordinate system
        Hyrot = - IHyx * const
 
        !field values rotated back to original x-y coordinate system
        fld(idx)%Hx(recidx) = fld(idx)%Hx(recidx) - Hyrot*sinbetasrc
      enddo receivershx !iy

        endif have_hx

        have_hy: if(bgdat%nHy .gt. 0) then

      receivershy: do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzHxy(irec)

        y = bgdat%Hypos(recidx,2) - refl_var%ys(ielemall)
        x = bgdat%Hypos(recidx,1) - refl_var%xs(ielemall)
        r = sqrt(x**2 + y**2)

        beta = atan2(y,x)
        betarot = beta - refl_var%betasrc(isrc)

        r_is_zerohy: if(r.EQ.0._real64) then
          !quick & dirty: skip the point if receiver is right at source point
          if(sz_eq_zr) then
            if(refl_var%infolevel.ge.output_more) &
              write(*,'(a)') 'WARNING: cannot handle a receiver at exactly the same location as a wire element yet, '// &
                'ignoring this contribution!'
            cycle
          endif

          !--> evaluate zero-order integrals only
          ibesord = 0
          IHyx = compute_1valr0(funcD0TE)
        else !r is not zero
          smallrhy: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .TRUE.
            call prepare_refcoef(refl_var,r,0,aniso)

            !evaluate all integrals just for this radius
            ibesord = 0
            IHyx = compute_1val(funcD0TE,r,sz,zr,ibesord,wellbehaved(2))
          else !r larger than threshold radius for spline interpolation
            !get integral values for this radius by interpolation from precomputed integral values
            IHyx = splinterp_1val(refl_var,2,r)
          endif smallrhy
        endif r_is_zerohy

        !field values in rotated coordinate system
        Hyrot = - IHyx * const
 
        !field values rotated back to original x-y coordinate system
        fld(idx)%Hy(recidx) = fld(idx)%Hy(recidx) + Hyrot*cosbetasrc
      enddo receivershy !iy

        endif have_hy
      endif hxy_equalpos
    enddo wire_elements


    !loop over wire end points
    wire_endpoints: do iep = 1,2
      xs = src%wire(idx)%endpos(1,iep)
      ys = src%wire(idx)%endpos(2,iep)

      sign = signs(iep)

      !same positions for Hx and Hy
      hxy_equalpos_ep: if(bgdat%nHxy.gt.0) then

      receivers1: do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzHxy(irec)

        y = bgdat%Hxypos(recidx,2) - ys
        x = bgdat%Hxypos(recidx,1) - xs
        r = sqrt(x**2 + y**2)

        !distinguish case r=0 since here we need only 1 integral
        chkr0: if(r .EQ. 0._real64) then

          !quick & dirty: skip the point if receiver is right at source point
          if(sz_eq_zr) then
            if(refl_var%infolevel.ge.output_more) &
              write(*,'(a)') 'WARNING: cannot handle a receiver at exactly the same location as a grounding point yet, '// &
                'ignoring this contribution!'
            cycle
          endif
        else !r is not zero

          beta = atan2(y,x)
          betarot = beta - refl_var%betasrc(isrc)
          cosbetarot = cos(betarot)
          sinbetarot = sin(betarot)

          smallr1: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .TRUE.
            call prepare_refcoef(refl_var,r,0,aniso)

            !evaluate all integrals just for this radius
            ibesord = 1
            IendHxy = compute_1val(funcdHxwire,r,sz,zr,ibesord,wellbehaved(6))
          else !r larger than threshold radius for spline interpolation

            !get integral values for this radius by interpolation from precomputed integral values
            IendHxy = splinterp_1val(refl_var,6,r)
          endif smallr1

          !include rotation here
          !> Exrot etc. are field components in coordiante system wih axes parallel to source wire
          !> E(irec,#) are field values in coord. system with original x,y axes
          term1 = sign * IendHxy / dfourpi
          Hxrot = term1 * sinbetarot
          Hyrot = -term1 * cosbetarot
          fld(idx)%Hx(recidx) = fld(Idx)%Hx(recidx) + Hxrot*cosbetasrc - Hyrot*sinbetasrc
          fld(idx)%Hy(recidx) = fld(idx)%Hy(recidx) + Hxrot*sinbetasrc + Hyrot*cosbetasrc
        endif chkr0

        !VTI: epsv derivatives
        dvert: if(with_dvert) then
          chkr0v: if(r .NE. 0._real64) then
            smallr1v: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

              !reflection coeff. for this radius
              refl_var%refcoef_changed = .TRUE.
              call prepare_refcoef(refl_var,r,0,aniso)

              !evaluate all integrals just for this radius
              ibesord = 1
              IendHxy = compute_1val(funcdHxwirev,r,sz,zr,ibesord,wellbehaved(10))
            else !r larger than threshold radius for spline interpolation
              !get integral values for this radius by interpolation from precomputed integral values
              IendHxy = splinterp_1val(refl_var,10,r)
            endif smallr1v

            !include rotation here
            !> Exrot etc. are field components in coordiante system wih axes parallel to source wire
            !> E(irec,#) are field values in coord. system with original x,y axes
            term1 = sign * IendHxy / dfourpi
            Hxrot = term1 * sinbetarot
            Hyrot = -term1 * cosbetarot
            fldv(idx)%Hx(recidx) = fldv(idx)%Hx(recidx) + Hxrot*cosbetasrc - Hyrot*sinbetasrc
            fldv(idx)%Hy(recidx) = fldv(idx)%Hy(recidx) + Hxrot*sinbetasrc + Hyrot*cosbetasrc
          endif chkr0v
        endif dvert
      enddo receivers1 !receivers

      else
        have_hx_ep: if(bgdat%nHx .gt. 0) then

      receivers1hx: do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzHxy(irec)

        y = bgdat%Hxpos(recidx,2) - ys
        x = bgdat%Hxpos(recidx,1) - xs
        r = sqrt(x**2 + y**2)

        !distinguish case r=0 since here we need only 1 integral
        chkr0hx: if(r .EQ. 0._real64) then

          !quick & dirty: skip the point if receiver is right at source point
          if(sz_eq_zr) then
            if(refl_var%infolevel.ge.output_more) &
              write(*,'(a)') 'WARNING: cannot handle a receiver at exactly the same location as a grounding point yet, '// &
                'ignoring this contribution!'
            cycle
          endif
        else !r is not zero

          beta = atan2(y,x)
          betarot = beta - refl_var%betasrc(isrc)
          cosbetarot = cos(betarot)
          sinbetarot = sin(betarot)

          smallr1hx: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .TRUE.
            call prepare_refcoef(refl_var,r,0,aniso)

            !evaluate all integrals just for this radius
            ibesord = 1
            IendHxy = compute_1val(funcdHxwire,r,sz,zr,ibesord,wellbehaved(6))
          else !r larger than threshold radius for spline interpolation

            !get integral values for this radius by interpolation from precomputed integral values
            IendHxy = splinterp_1val(refl_var,6,r)
          endif smallr1hx

          !include rotation here
          !> Exrot etc. are field components in coordiante system wih axes parallel to source wire
          !> E(irec,#) are field values in coord. system with original x,y axes
          term1 = sign * IendHxy / dfourpi
          Hxrot = term1 * sinbetarot
          Hyrot = -term1 * cosbetarot
          fld(idx)%Hx(recidx) = fld(Idx)%Hx(recidx) + Hxrot*cosbetasrc - Hyrot*sinbetasrc
        endif chkr0hx

        !VTI: epsv derivatives
        dverthx: if(with_dvert) then
          chkr0vhx: if(r .NE. 0._real64) then
            smallr1vhx: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

              !reflection coeff. for this radius
              refl_var%refcoef_changed = .TRUE.
              call prepare_refcoef(refl_var,r,0,aniso)

              !evaluate all integrals just for this radius
              ibesord = 1
              IendHxy = compute_1val(funcdHxwirev,r,sz,zr,ibesord,wellbehaved(10))
            else !r larger than threshold radius for spline interpolation
              !get integral values for this radius by interpolation from precomputed integral values
              IendHxy = splinterp_1val(refl_var,10,r)
            endif smallr1vhx

            !include rotation here
            !> Exrot etc. are field components in coordiante system wih axes parallel to source wire
            !> E(irec,#) are field values in coord. system with original x,y axes
            term1 = sign * IendHxy / dfourpi
            Hxrot = term1 * sinbetarot
            Hyrot = -term1 * cosbetarot
            fldv(idx)%Hx(recidx) = fldv(idx)%Hx(recidx) + Hxrot*cosbetasrc - Hyrot*sinbetasrc
          endif chkr0vhx
        endif dverthx
      enddo receivers1hx !receivers

        endif have_hx_ep

        have_hy_ep: if(bgdat%nHy .gt. 0) then

      receivers1hy: do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzHxy(irec)

        y = bgdat%Hypos(recidx,2) - ys
        x = bgdat%Hypos(recidx,1) - xs
        r = sqrt(x**2 + y**2)

        !distinguish case r=0 since here we need only 1 integral
        chkr0hy: if(r .EQ. 0._real64) then

          !quick & dirty: skip the point if receiver is right at source point
          if(sz_eq_zr) then
            if(refl_var%infolevel.ge.output_more) &
              write(*,'(a)') 'WARNING: cannot handle a receiver at exactly the same location as a grounding point yet, '// &
                'ignoring this contribution!'
            cycle
          endif
        else !r is not zero

          beta = atan2(y,x)
          betarot = beta - refl_var%betasrc(isrc)
          cosbetarot = cos(betarot)
          sinbetarot = sin(betarot)

          smallr1hy: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .TRUE.
            call prepare_refcoef(refl_var,r,0,aniso)

            !evaluate all integrals just for this radius
            ibesord = 1
            IendHxy = compute_1val(funcdHxwire,r,sz,zr,ibesord,wellbehaved(6))
          else !r larger than threshold radius for spline interpolation

            !get integral values for this radius by interpolation from precomputed integral values
            IendHxy = splinterp_1val(refl_var,6,r)
          endif smallr1hy

          !include rotation here
          !> Exrot etc. are field components in coordiante system wih axes parallel to source wire
          !> E(irec,#) are field values in coord. system with original x,y axes
          term1 = sign * IendHxy / dfourpi
          Hxrot = term1 * sinbetarot
          Hyrot = -term1 * cosbetarot
          fld(idx)%Hy(recidx) = fld(idx)%Hy(recidx) + Hxrot*sinbetasrc + Hyrot*cosbetasrc
        endif chkr0hy

        !VTI: epsv derivatives
        dverthy: if(with_dvert) then
          chkr0vhy: if(r .NE. 0._real64) then
            smallr1vhy: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

              !reflection coeff. for this radius
              refl_var%refcoef_changed = .TRUE.
              call prepare_refcoef(refl_var,r,0,aniso)

              !evaluate all integrals just for this radius
              ibesord = 1
              IendHxy = compute_1val(funcdHxwirev,r,sz,zr,ibesord,wellbehaved(10))
            else !r larger than threshold radius for spline interpolation
              !get integral values for this radius by interpolation from precomputed integral values
              IendHxy = splinterp_1val(refl_var,10,r)
            endif smallr1vhy

            !include rotation here
            !> Exrot etc. are field components in coordiante system wih axes parallel to source wire
            !> E(irec,#) are field values in coord. system with original x,y axes
            term1 = sign * IendHxy / dfourpi
            Hxrot = term1 * sinbetarot
            Hyrot = -term1 * cosbetarot
            fldv(idx)%Hy(recidx) = fldv(idx)%Hy(recidx) + Hxrot*sinbetasrc + Hyrot*cosbetasrc
          endif chkr0vhy
        endif dverthy
      enddo receivers1hy !receivers

        endif have_hy_ep
      endif hxy_equalpos_ep

    enddo wire_endpoints
  enddo wires  !wires


endsubroutine interp_intvals_wire_Hxy


!------------------------------------------------------------
!>  1D EM subroutine interp_intvals_wire_Hz
!
!>  get interpolated field values at receiver locations, horizontal wire source,
!>    Hz only
!
!>  Rita Streich 2009-2011
!------------------------------------------------------------
subroutine interp_intvals_wire_Hz(refl_var,bgdat,fld,src,sz,zr,ommu, funcAz1TE,ilay)

  implicit none

  !external variables
  type(refl_struct) :: refl_var   !everything needed throughout 1D computations
  type(backgrounddata) :: bgdat      !coordinate vectors and final output EM fields
  type(receiverdata),dimension(:) :: fld        !field or derivative vectors for each wire, all components
  type(sorec),intent( in ) :: src        !source definition (we need the currents here)
  real(kind=real64),intent( in ) :: sz,zr      !source and receiver depth
  real(kind=real64),intent( in ) :: ommu       !omega * mu0
  complex(kind=real64),external :: funcAz1TE
  integer(kind=int32),intent( in ) :: ilay       !layer index for derivatives - leave at zero for forward modeling

  !internal variables
  integer(kind=int32) :: isrc    !source element counter
  real(kind=real64) :: x,y,r   !temp source-receiver distances
  integer(kind=int32) :: irec    !receiver counter
  real(kind=real64) :: beta,betarot      !temp angles

  complex(kind=real64) :: IHz          !interpolated integral values for integrals along wire

  logical,dimension(nintWiredvti) :: wellbehaved   !indicates if Hankel integration can be used
  logical :: sz_eq_zr      !indicates if source and reveicer are at the same depth
  integer(kind=int32),dimension(NREL) :: ibesord   !bessel function order (integer "array")

  complex(kind=real64) :: Hzrot    !temp field values in cylindrical coordinates
  real(kind=real64) :: cosbetasrc,sinbetasrc !cos(betasrc) and sin(betasrc), precompute for efficiency
  real(kind=real64) :: const    !geometric constant
  complex(kind=real64) :: constHz  !geometric constant for Hz
  integer(kind=int32) :: idx      !source element index
  integer(kind=int32) :: recidx   !receiver index
  integer(kind=int32) :: ielemall,ielem       !wire element counters

  !indicators for fast Hankel transform or adaptive integration
  wellbehaved = .TRUE.
  sz_eq_zr = .FALSE.
  if(sz.EQ.zr) then
    wellbehaved(1:2) = .FALSE.
    wellbehaved(4:5) = .FALSE.
    sz_eq_zr = .TRUE.
  endif

  !loop over wires
  !TE only --> no contribution to epsv derivatives
  ielemall = 0
  wires: do isrc = refl_var%isrcstart,refl_var%isrcend

    !get constant geometry factor: wire element length
    idx = refl_var%isrcperz(isrc)
    const = src%wire(idx)%dlw / dfourpi
    constHz = const * dci / ommu

    cosbetasrc = cos(refl_var%betasrc(isrc))
    sinbetasrc = sin(refl_var%betasrc(isrc))

    !loop over wire elements
    wire_elements: do ielem = 1,src%nelem(idx)
      ielemall = ielemall + 1

      receivers: do irec=refl_var%irecstart,refl_var%irecend
        recidx = refl_var%irecperzHz(irec)

        y = bgdat%Hzpos(recidx,2) - refl_var%ys(ielemall)
        x = bgdat%Hzpos(recidx,1) - refl_var%xs(ielemall)
        r = sqrt(x**2 + y**2)

        beta = atan2(y,x)
        betarot = beta - refl_var%betasrc(isrc)

        r_is_zero: if(r.EQ.0._real64) then

          !quick & dirty: skip the point if receiver is right at source point
          if(sz_eq_zr) then
            if(refl_var%infolevel.ge.output_more) &
              write(*,'(a)') 'WARNING: cannot handle a receiver at exactly the same location as a wire element yet, '// &
                'ignoring this contribution!'
            cycle
          endif
          IHz = 0._real64
        else !r is not zero
          smallr: if(r.lt.rsplmin) then !r is but smaller than threshold radius for spline interpolation

            !reflection coeff. for this radius
            refl_var%refcoef_changed = .TRUE.
            call prepare_refcoef(refl_var,r,0,aniso)

            !evaluate all integrals just for this radius
            ibesord = 1
            IHz = compute_1val(funcAz1TE,r,sz,zr,ibesord,wellbehaved(3))
          else !r larger than threshold radius for spline interpolation
            !get integral values for this radius by interpolation from precomputed integral values
            IHz = splinterp_1val(refl_var,3,r)
          endif smallr
        endif r_is_zero

        !field values in rotated coordinate system
        Hzrot = IHz * constHz
        !field values rotated back to original x-y coordinate system
        fld(idx)%Hz(recidx) = fld(idx)%Hz(recidx) + Hzrot*sin(betarot)
      enddo receivers !iy
    enddo wire_elements

  enddo wires  !wires

endsubroutine interp_intvals_wire_Hz


