!---------------------------------------------------------
! EM1D subroutine addwirefields
!
! add fields for the separate wire segments of a multi-segment wire source
!
! Rita Streich 2011
!---------------------------------------------------------
subroutine addwirefields(bgdat,refl_var,src,icur,ifreq)

  implicit none

  !external variables
  type(backgrounddata) :: bgdat      !coordinate vectors and output EM fields
  type(refl_struct)    :: refl_var   !all variables that have to be remembered while computing 1D fields
  type(sorec)          :: src        !a single source
  integer(kind=int32)  :: icur       !source current counter
  integer(kind=int32)  :: ifreq      !index of frequency component

  !internal variables
  integer(kind=int32)            :: icurstart     !temp index for wire currents
  integer(kind=int32)            :: iwire         !wire counter
  complex(kind=real64)           :: cur           !temp source current
  integer(kind=int32)            :: irec          !iReceiver counter
  integer(kind=int32)            :: ilay          !layer counter for derivatives
  integer(kind=int32)            :: nrecEx,nrecEy,nrecEz,nrecHx,nrecHy,nrecHz   !nr of receivers for each field component


  nrecEx = size(bgdat%Ex,1)
  nrecEy = size(bgdat%Ey,1)
  nrecEz = size(bgdat%Ez,1)
  nrecHx = size(bgdat%Hx,1)
  nrecHy = size(bgdat%Hy,1)
  nrecHz = size(bgdat%Hz,1)


  !add up fields from separate wires
  icurstart = (icur-1)*src%nwire

  if ((bgdat%dowhat.eq.fwdmodel) .or. (bgdat%dowhat.eq.fwd_deriv)) then
    addfields: do iwire=1,src%nwire
      !take complex conjugate of current here to match Loeseth's sign convention
      !careful: some frequency components of the currents can be zero
      !--> do NOT overwrite values in EHwire, but multiply the current at the very end!
      cur = conjg(src%cur(icurstart+iwire,ifreq))
      do irec = 1,nrecEx
        bgdat%Ex(irec) = bgdat%Ex(irec) + refl_var%EHwire(iwire)%Ex(irec) * cur
      enddo
      do irec = 1,nrecEy
        bgdat%Ey(irec) = bgdat%Ey(irec) + refl_var%EHwire(iwire)%Ey(irec) * cur
      enddo
      do irec = 1,nrecEz
        bgdat%Ez(irec) = bgdat%Ez(irec) + refl_var%EHwire(iwire)%Ez(irec) * cur
      enddo
      do irec = 1,nrecHx
        bgdat%Hx(irec) = bgdat%Hx(irec) + refl_var%EHwire(iwire)%Hx(irec) * cur
      enddo
      do irec = 1,nrecHy
        bgdat%Hy(irec) = bgdat%Hy(irec) + refl_var%EHwire(iwire)%Hy(irec) * cur
      enddo
      do irec = 1,nrecHz
        bgdat%Hz(irec) = bgdat%Hz(irec) + refl_var%EHwire(iwire)%Hz(irec) * cur
      enddo
    enddo addfields
  endif

  if (bgdat%dowhat.ge.deriv) then
    addfieldsderiv: do iwire=1,src%nwire
      !careful: some frequency components of the currents can be zero
      !--> do NOT overwrite values in Ewire, but multiply the current at the very end!
      cur = conjg(src%cur(icurstart+iwire,ifreq))
      do ilay = 1,nlay
        do irec=1,nrecEx
          bgdat%dExdm(irec,ilay) = bgdat%dExdm(irec,ilay) + refl_var%EHwirederiv(iwire,ilay)%Ex(irec) * cur
        enddo
        do irec=1,nrecEy
          bgdat%dEydm(irec,ilay) = bgdat%dEydm(irec,ilay) + refl_var%EHwirederiv(iwire,ilay)%Ey(irec) * cur
        enddo
        do irec=1,nrecEz
          bgdat%dEzdm(irec,ilay) = bgdat%dEzdm(irec,ilay) + refl_var%EHwirederiv(iwire,ilay)%Ez(irec) * cur
        enddo
        do irec=1,nrecHx
          bgdat%dHxdm(irec,ilay) = bgdat%dHxdm(irec,ilay) + refl_var%EHwirederiv(iwire,ilay)%Hx(irec) * cur
        enddo
        do irec=1,nrecHy
          bgdat%dHydm(irec,ilay) = bgdat%dHydm(irec,ilay) + refl_var%EHwirederiv(iwire,ilay)%Hy(irec) * cur
        enddo
        do irec=1,nrecHz
          bgdat%dHzdm(irec,ilay) = bgdat%dHzdm(irec,ilay) + refl_var%EHwirederiv(iwire,ilay)%Hz(irec) * cur
        enddo
      enddo
    enddo addfieldsderiv
    if (bgdat%aniso .eq. vti) then
      addfieldsderivv: do iwire=1,src%nwire
        cur = conjg(src%cur(icurstart+iwire,ifreq))
        do ilay = 1,nlay
          do irec = 1,nrecEx
            bgdat%dExdmv(irec,ilay) = bgdat%dExdmv(irec,ilay) + refl_var%EHwirederivv(iwire,ilay)%Ex(irec) * cur
          enddo
          do irec = 1,nrecEy
            bgdat%dEydmv(irec,ilay) = bgdat%dEydmv(irec,ilay) + refl_var%EHwirederivv(iwire,ilay)%Ey(irec) * cur
          enddo
          do irec = 1,nrecEz
            bgdat%dEzdmv(irec,ilay) = bgdat%dEzdmv(irec,ilay) + refl_var%EHwirederivv(iwire,ilay)%Ez(irec) * cur
          enddo
          do irec = 1,nrecHx
            bgdat%dHxdmv(irec,ilay) = bgdat%dHxdmv(irec,ilay) + refl_var%EHwirederivv(iwire,ilay)%Hx(irec) * cur
          enddo
          do irec = 1,nrecHy
            bgdat%dHydmv(irec,ilay) = bgdat%dHydmv(irec,ilay) + refl_var%EHwirederivv(iwire,ilay)%Hy(irec) * cur
          enddo
          do irec = 1,nrecHz
            bgdat%dHzdmv(irec,ilay) = bgdat%dHzdmv(irec,ilay) + refl_var%EHwirederivv(iwire,ilay)%Hz(irec) * cur
          enddo
        enddo
      enddo addfieldsderivv
    endif
  endif !derivatives requested

endsubroutine addwirefields
