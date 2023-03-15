!---------------------------------------------------------
! EM1D subroutine conjugate_fields
!
! take complex conjugate of field / sensitivity values to match FFT convention of 
!   3D and 2.5D codes, which is opposite of the one used in Loeseth's formulas
! for derivatives: apply factor to take deriv. with respect to ln(sigma)
!
! Rita Streich 2011
!---------------------------------------------------------
subroutine conjugate_fields(bgdat,refl_var,src)

  implicit none

  !external variables
  type(backgrounddata)     :: bgdat              !coordinate vectors and output EM fields
  type(refl_struct)        :: refl_var   !all variables that have to be remembered while computing 1D fields
  type(sorec),intent(in)   :: src                !source description

  !internal variables
  complex(kind=real64)     :: logfact            !factor for getting derivatives with respect to log(sigma)
  integer(kind=int32)      :: ilay               !layer counter
  integer(kind=int32)      :: iw                 !wire counter


  if ((bgdat%dowhat .eq. fwdmodel) .or. (bgdat%dowhat.eq.fwd_deriv)) then
    if (src%type.eq.dipole) then
      bgdat%Ex = conjg(bgdat%Ex)
      bgdat%Ey = conjg(bgdat%Ey)
      bgdat%Ez = conjg(bgdat%Ez)
	 if ( bgdat%nHx .gt. 0) then
      bgdat%Hx = conjg(bgdat%Hx)
      bgdat%Hy = conjg(bgdat%Hy)
      bgdat%Hz = conjg(bgdat%Hz)
	 end if
    else
      do iw=1,src%nwire
        refl_var%EHwire(iw)%Ex = conjg(refl_var%EHwire(iw)%Ex)
        refl_var%EHwire(iw)%Ey = conjg(refl_var%EHwire(iw)%Ey)
        refl_var%EHwire(iw)%Ez = conjg(refl_var%EHwire(iw)%Ez)
		if ( bgdat%nHx .gt. 0) then
        refl_var%EHwire(iw)%Hx = conjg(refl_var%EHwire(iw)%Hx)
        refl_var%EHwire(iw)%Hy = conjg(refl_var%EHwire(iw)%Hy)
        refl_var%EHwire(iw)%Hz = conjg(refl_var%EHwire(iw)%Hz)
		end if
      enddo
    endif
  endif

  if (bgdat%dowhat.ge.deriv) then
    if (src%type.eq.dipole) then
      !apply factor for getting derivatives with respect to ln(sigma) instead of complex epsilon
      bgdat%dExdm = conjg(bgdat%dExdm)
      bgdat%dEydm = conjg(bgdat%dEydm)
      bgdat%dEzdm = conjg(bgdat%dEzdm)
	  if ( bgdat%nHx .gt. 0) then
      bgdat%dHxdm = conjg(bgdat%dHxdm)
      bgdat%dHydm = conjg(bgdat%dHydm)
      bgdat%dHzdm = conjg(bgdat%dHzdm)
	  end if
      do ilay=1,nlay
        logfact = conjg(dci * bgdat%sigh(ilay) / bgdat%omega)  !use sigma only, not complex conductivity!
        bgdat%dExdm(:,ilay) = bgdat%dExdm(:,ilay) * logfact
        bgdat%dEydm(:,ilay) = bgdat%dEydm(:,ilay) * logfact
        bgdat%dEzdm(:,ilay) = bgdat%dEzdm(:,ilay) * logfact
		if ( bgdat%nHx .gt. 0) then
        bgdat%dHxdm(:,ilay) = bgdat%dHxdm(:,ilay) * logfact
        bgdat%dHydm(:,ilay) = bgdat%dHydm(:,ilay) * logfact
        bgdat%dHzdm(:,ilay) = bgdat%dHzdm(:,ilay) * logfact
		end if
      enddo

      if (bgdat%aniso .eq. vti) then
        bgdat%dExdmv = conjg(bgdat%dExdmv)
        bgdat%dEydmv = conjg(bgdat%dEydmv)
        bgdat%dEzdmv = conjg(bgdat%dEzdmv)
		if ( bgdat%nHx .gt. 0) then
        bgdat%dHxdmv = conjg(bgdat%dHxdmv)
        bgdat%dHydmv = conjg(bgdat%dHydmv)
        bgdat%dHzdmv = conjg(bgdat%dHzdmv)
		end if
        do ilay=1,nlay
          logfact = conjg(dci * bgdat%sigv(ilay) / bgdat%omega)
          bgdat%dExdmv(:,ilay) = bgdat%dExdmv(:,ilay) * logfact
          bgdat%dEydmv(:,ilay) = bgdat%dEydmv(:,ilay) * logfact
          bgdat%dEzdmv(:,ilay) = bgdat%dEzdmv(:,ilay) * logfact
		  if ( bgdat%nHx .gt. 0) then
          bgdat%dHxdmv(:,ilay) = bgdat%dHxdmv(:,ilay) * logfact
          bgdat%dHydmv(:,ilay) = bgdat%dHydmv(:,ilay) * logfact
          bgdat%dHzdmv(:,ilay) = bgdat%dHzdmv(:,ilay) * logfact
		  end if
        enddo
      endif

    else !wire source
      do iw=1,src%nwire
        do ilay=1,nlay
          !apply factor for getting derivatives with respect to ln(sigma) instead of complex epsilon
          refl_var%EHwirederiv(iw,ilay)%Ex = conjg(refl_var%EHwirederiv(iw,ilay)%Ex)
          refl_var%EHwirederiv(iw,ilay)%Ey = conjg(refl_var%EHwirederiv(iw,ilay)%Ey)
          refl_var%EHwirederiv(iw,ilay)%Ez = conjg(refl_var%EHwirederiv(iw,ilay)%Ez)
		  if ( bgdat%nHx .gt. 0) then
          refl_var%EHwirederiv(iw,ilay)%Hx = conjg(refl_var%EHwirederiv(iw,ilay)%Hx)
          refl_var%EHwirederiv(iw,ilay)%Hy = conjg(refl_var%EHwirederiv(iw,ilay)%Hy)
          refl_var%EHwirederiv(iw,ilay)%Hz = conjg(refl_var%EHwirederiv(iw,ilay)%Hz)
		  end if
          logfact = conjg(dci * bgdat%sigh(ilay) / bgdat%omega)  !use sigma only, not complex conductivity!
          refl_var%EHwirederiv(iw,ilay)%Ex = refl_var%EHwirederiv(iw,ilay)%Ex * logfact
          refl_var%EHwirederiv(iw,ilay)%Ey = refl_var%EHwirederiv(iw,ilay)%Ey * logfact
          refl_var%EHwirederiv(iw,ilay)%Ez = refl_var%EHwirederiv(iw,ilay)%Ez * logfact
		  if ( bgdat%nHx .gt. 0) then
          refl_var%EHwirederiv(iw,ilay)%Hx = refl_var%EHwirederiv(iw,ilay)%Hx * logfact
          refl_var%EHwirederiv(iw,ilay)%Hy = refl_var%EHwirederiv(iw,ilay)%Hy * logfact
          refl_var%EHwirederiv(iw,ilay)%Hz = refl_var%EHwirederiv(iw,ilay)%Hz * logfact
		  end if
        enddo

        if (bgdat%aniso .eq. vti) then
          do ilay=1,nlay
            refl_var%EHwirederivv(iw,ilay)%Ex = conjg(refl_var%EHwirederivv(iw,ilay)%Ex)
            refl_var%EHwirederivv(iw,ilay)%Ey = conjg(refl_var%EHwirederivv(iw,ilay)%Ey)
            refl_var%EHwirederivv(iw,ilay)%Ez = conjg(refl_var%EHwirederivv(iw,ilay)%Ez)
			if ( bgdat%nHx .gt. 0) then
            refl_var%EHwirederivv(iw,ilay)%Hx = conjg(refl_var%EHwirederivv(iw,ilay)%Hx)
            refl_var%EHwirederivv(iw,ilay)%Hy = conjg(refl_var%EHwirederivv(iw,ilay)%Hy)
            refl_var%EHwirederivv(iw,ilay)%Hz = conjg(refl_var%EHwirederivv(iw,ilay)%Hz)
			end if
            logfact = conjg(dci * bgdat%sigv(ilay) / bgdat%omega)
            refl_var%EHwirederivv(iw,ilay)%Ex = refl_var%EHwirederivv(iw,ilay)%Ex * logfact
            refl_var%EHwirederivv(iw,ilay)%Ey = refl_var%EHwirederivv(iw,ilay)%Ey * logfact
            refl_var%EHwirederivv(iw,ilay)%Ez = refl_var%EHwirederivv(iw,ilay)%Ez * logfact
			if ( bgdat%nHx .gt. 0) then
            refl_var%EHwirederivv(iw,ilay)%Hx = refl_var%EHwirederivv(iw,ilay)%Hx * logfact
            refl_var%EHwirederivv(iw,ilay)%Hy = refl_var%EHwirederivv(iw,ilay)%Hy * logfact
            refl_var%EHwirederivv(iw,ilay)%Hz = refl_var%EHwirederivv(iw,ilay)%Hz * logfact
			end if
          enddo
        endif
      enddo
    endif
  endif

endsubroutine conjugate_fields
