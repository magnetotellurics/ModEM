! *****************************************************************************
module responses
	! Module responses contains the variables and the basic subroutines to do
	! with the data responses and the penalty functional and their derivatives

  use SolnSpace
  use receivers
  implicit none


Contains

  ! ***************************************************************************
  ! * refavg computes an analogue of Dst index: weighted sum of H_\theta for
  ! * reference mid-latitude observatories. This average value is used to scale
  ! * the raw magnetic field data at observatory locations to obtain responses.
  ! * Uses the reference observatory list refObsList stored in the receivers module.

  function refavg(H) result(Dst)

    type (cvector), intent(in)                      :: H
    complex(8)                                      :: Dst

    type (receiver_t)                               :: obs
    complex(8)                                      :: Hy
    integer                                         :: i,nobs

    Dst = C_ZERO
    nobs = 0

    do i = 1,refObsList%n

        obs = refObsList%info(i)

        if (.not.obs%located) then
            call LocateReceiver(H%grid,obs)
        end if

        if (.not.obs%defined) then
            write(0,*) 'Warning: Reference observatory ',trim(obs%code),' is not defined at refavg; skip it'
            cycle
        end if

        Hy = dotProd_noConj(obs%Ly,H) ! dotProd_noConj_scvector_f...could do real here

        Dst = Dst + Hy
        nobs = nobs + 1

     end do

     if (nobs > 0) then
        Dst = Dst / nobs
     else
        Dst = 1
     end if

  end function refavg ! refavg

  ! ***************************************************************************
  ! * drefavg computes the derivative of the refavg: 1/n is n > 0 and obs is
  ! * a reference observatory, otherwise zero; used in Lrows
  ! * assume that all reference observatories are located by now... done that
  ! * upon initialization, then again in refavg

  function drefavg(obs) result(davg)

    type (receiver_t), intent(in)                   :: obs
    real(8)                                         :: davg

    integer                                         :: i,nobs
    logical                                         :: refobs

    davg = R_ZERO
    refobs = .false.
    if (getObs(refObsList,obs%code) > 0) then
        refobs = .true.
    else
        return
    end if

    ! if observatory is in the list, count all the observatories used for refavg
    nobs = 0
    do i = 1,refObsList%n
        if (.not.refObsList%info(i)%defined) then
            cycle
        end if
        nobs = nobs + 1
    end do

    if (nobs > 0) then
       davg = 1. / nobs
    else
       davg = R_ZERO
    end if

  end function drefavg ! drefavg


  ! ***************************************************************************
  ! * fieldValue computes the raw magnetic fields at a given observatory location
  ! * using interpolation procedures through LocateReceiver->ComputeInterpWeights.

  function fieldValue_rx(obs,H) result(Resp)

    complex(8), external                            :: func
    type (receiver_t), intent(inout)                :: obs
    type (cvector), intent(in)                      :: H
    complex(8)                                      :: Resp(3)

    if (.not.obs%located) then
      call LocateReceiver(H%grid,obs)
    end if

    if (.not.obs%defined) then
      write(0,*) 'Observatory ',trim(obs%code),' is not defined at fieldValue; hence ignore'
      return
    end if

    Resp(1) = dotProd_noConj(obs%Lx,H) ! dotProd_noConj_scvector_f...could do real here
    Resp(2) = dotProd_noConj(obs%Ly,H)
    Resp(3) = dotProd_noConj(obs%Lz,H)

  end function fieldValue_rx ! fieldValue

  ! ***************************************************************************
  ! * dataResp performs the calculation that returns a response of a
  ! * chosen kind at a chosen observatory location, given H defined on grid
  ! * We are likely to need the responses at observatories (use this subroutine)
  ! * and at grid nodes on the surface. The latter is very easy to compute.
  ! * Use e.g. C(H%x(i,j,nzAir+1),H%y(i,j,nzAir+1),H%z(i,j,nzAir+1),grid%th(j)).
  ! * For responses at locations other than the observatory, construct an
  ! * observatory type (receiver_t) at the required location and call this routine.
  ! * Interpolation is automatic through LocateReceiver->ComputeInterpWeights.

  function dataResp_rx(func,obs,H) result(Resp)

	complex(8), external							:: func
	type (receiver_t), intent(inout)					:: obs
	type (cvector), intent(in)						:: H
	complex(8)										:: Resp

	type (sparsevecc)								:: Lx,Ly,Lz
	complex(8)										:: Hx,Hy,Hz


	if (.not.obs%located) then
	  call LocateReceiver(H%grid,obs)
	end if

	if (.not.obs%defined) then
	  !write(0,*) 'Observatory ',trim(obs%code),' is not defined at dataResp; hence ignore'
	  return
	end if

	Lx = obs%Lx
	Ly = obs%Ly
	Lz = obs%Lz

	Hx = dotProd_noConj(Lx,H) ! dotProd_noConj_scvector_f...could do real here
	Hy = dotProd_noConj(Ly,H)
	Hz = dotProd_noConj(Lz,H)

	! Calculate the response at this location
	Resp = func(Hx,Hy,Hz,obs%colat*d2r)

	!write(*,'(a40,2i6,5g15.7)') 'i,j,lon,colat,Hx,Hy,Hz = ',&
	!	obs%i,obs%j,obs%lon,obs%colat,Hx,Hy,Hz

	call deall_sparsevecc(Lx)
	call deall_sparsevecc(Ly)
	call deall_sparsevecc(Lz)

  end function dataResp_rx ! response


  ! ***************************************************************************
  ! * dataResp_ijk performs a calculation such that the output is
  ! * the response at a cell (defined at the midpoint of the interval from the
  ! * (i,j,k)'th to the (i,j+1,k)'th node), given H defined on grid
  ! * With the current interpolation routine, k is assumed to be the index
  ! * of the Earth's surface

  function dataResp_ijk(func,i,j,k,H) result(Resp)

	complex(8), external							:: func
	integer, intent(in)								:: i,j,k
	type (cvector), intent(in)						:: H
	complex(8)										:: Resp
	type (receiver_t)									:: obs
	real(8)											:: rad,lon,colat
	character(80)									:: name

	type (sparsevecc)								:: Lx,Ly,Lz
	complex(8)										:: Hx,Hy,Hz

	rad = H%grid%r(k)
	lon = H%grid%ph(i)*r2d
	colat = (H%grid%th(j)+H%grid%th(j+1))*r2d/2

	write(name,'(2i5)') i,j
	name = trim(name)

	! Build an observatory at this cell with code = name
	call CreateReceiver(obs,rad,lon,colat,name)

	! Compute interpolation parameters
	call LocateReceiver(H%grid,obs)

	if (.not.obs%defined) then
	  return
	end if

	Lx = obs%Lx
	Ly = obs%Ly
	Lz = obs%Lz

	Hx = dotProd_noConj(Lx,H) ! dotProd_noConj_scvector_f...could do real here
	Hy = dotProd_noConj(Ly,H)
	Hz = dotProd_noConj(Lz,H)

	! Calculate the response at this location
	Resp = func(Hx,Hy,Hz,obs%colat*d2r)

	!write(*,'(a40,2i6,5g15.7)') 'i,j,lon,colat,Hx,Hy,Hz = ',&
	!	obs%i,obs%j,obs%lon,obs%colat,Hx,Hy,Hz
	call deall_sparsevecc(Lx)
	call deall_sparsevecc(Ly)
	call deall_sparsevecc(Lz)

  end function dataResp_ijk ! dataResp_ijk


  ! ***************************************************************************
  ! * fieldValue_ij computes a set of magnetic field solutions
  ! * at the cells (defined at the midpoint of the interval from the
  ! * (i,j,k)'th to the (i,j+1,k)'th node), or at nodes, given H defined on grid, at
  ! * the specified grid radius

  subroutine fieldValue_ij(rad,H,Hij,type)

	real(8), intent(in)								:: rad
	type (cvector), intent(in)						:: H
	type (solution_t), intent(out)					:: Hij
	character(*), intent(in)                        :: type
	real(8)											:: lon,colat
	character(80)									:: name
	integer											:: i,j,istat
	integer											:: nx,ny

    if (trim(type) .eq. trim(CENTER)) then
        nx = H%grid%nx
        ny = H%grid%ny-2
    else if (trim(type) .eq. trim(CORNER)) then
        nx = H%grid%nx
        ny = H%grid%ny-1
    else
        call warning('Unknown solution type in fieldValue_ij; using '//trim(CENTER))
        nx = H%grid%nx
        ny = H%grid%ny-2
    end if

	allocate(Hij%x(nx,ny),Hij%y(nx,ny),Hij%z(nx,ny),Hij%o(nx,ny),STAT=istat)

	Hij%x = C_ZERO
	Hij%y = C_ZERO
	Hij%z = C_ZERO

	do j = 1,ny
	  do i = 1,nx

        if (trim(type) .eq. trim(CORNER)) then
		    lon = H%grid%ph(i)*r2d
            colat = H%grid%th(j+1)*r2d
        else
            !Use if obs is at midpoint between (i,j,k)'th to the (i,j+1,k)'th nodes
            !lon = grid%ph(i)*r2d
            !Use if obs is at the center of the (i,j,k)'th cell
            if (i == nx) then
              lon = (H%grid%ph(i)+2*PI)*r2d/2
            else
              lon = (H%grid%ph(i)+H%grid%ph(i+1))*r2d/2
            end if
            colat = (H%grid%th(j+1)+H%grid%th(j+2))*r2d/2
		end if

		write(name,'(2i4)') i,j+1
		name = trim(name)

		! Build an observatory at this cell with code = name
		call CreateReceiver(Hij%o(i,j),rad,lon,colat,name)

		! NB: The receiver is shifted down to the nearest grid radius!!!
		call LocateReceiver(H%grid,Hij%o(i,j))

		if (Hij%o(i,j)%located) then
		  Hij%x(i,j) = dotProd_noConj(Hij%o(i,j)%Lx,H) ! dotProd_noConj_scvector_f
		  Hij%y(i,j) = dotProd_noConj(Hij%o(i,j)%Ly,H)
		  Hij%z(i,j) = dotProd_noConj(Hij%o(i,j)%Lz,H)
		end if

	  end do
	end do

  end subroutine fieldValue_ij ! fieldValue_ij


end module responses
