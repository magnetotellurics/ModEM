! *****************************************************************************
module DataFunc
  ! Module containing the subroutines for a-posteriori analysis of the output
  ! vector h=lH from the routine SolveMaxwells

  use model_operators
  use SolnSpace
  use dataTypes
  use transmitters
  use receivers
  use functionals
  use responses
  implicit none

  public                        :: dataResp, Lrows, Qrows

Contains

!******************************************************************************
  subroutine dataResp(H,m0,iDt,iRx,Resp)
  ! given magnetic field solution and indices into data types and receiver
  ! dictionaries, compute the single complex response function; store output
  ! in a vector or real values

  implicit none
  type (solnVector_t), intent(in)   :: H
  type (modelParam_t), intent(in)   :: m0 ! currently not used
  integer, intent(in)           :: iDt
  integer, intent(in)           :: iRx
  real(kind=prec), intent(inout)    :: Resp(:)

  !  local variables
  complex(8)            :: res(3)
  type (functional_t), pointer   :: dataType
  type (receiver_t), pointer     :: obs
  integer               :: iComp,nComp

  dataType => TFList%info(iDt)
  obs => obsList%info(iRx)

  ! compute the complex data response (C/D ratio or the scaled H fields)

  select case (dataType%name)
     case('C')
        res(1) = dataResp_rx(C_ratio,obs,H%vec)
        ! save real output
        Resp(1) = dreal(res(1))
        Resp(2) = dimag(res(1))
     case('D')
        res(1) = dataResp_rx(D_ratio,obs,H%vec)
        ! save real output
        Resp(1) = dreal(res(1))
        Resp(2) = dimag(res(1))
     case('T') ! generalised magnetic field transfer functions
        res = fieldValue_rx(obs,H%vec)
        res = res / refavg(H%vec)
        ! save real output
        Resp(1) = dreal(res(1))
        Resp(2) = dimag(res(1))
        Resp(3) = dreal(res(2))
        Resp(4) = dimag(res(2))
        Resp(5) = dreal(res(3))
        Resp(6) = dimag(res(3))
     case('H') ! raw magnetic fields
        res = fieldValue_rx(obs,H%vec)
        ! save real output
        Resp(1) = dreal(res(1))
        Resp(2) = dimag(res(1))
        Resp(3) = dreal(res(2))
        Resp(4) = dimag(res(2))
        Resp(5) = dreal(res(3))
        Resp(6) = dimag(res(3))
     case default
        call errStop('Unknown data response specified in dataResp')
  end select


  end subroutine dataResp

  ! ***************************************************************************
  ! * Lrows is a subroutine to output a full vector g_j defined on edges,
  ! * such that for a single frequency and a single observatory,
  ! * $\pd{psi_{\omega}^j}{veca} = g_j^* \pd{vecH}{veca}$
  ! *
  ! * In my complex global code, a similar routine computed
  ! * g_sparse = conjg(gc_sparse)
  ! * i.e., g rather than g* (conjugated).
  ! * However, for Lrows as defined in the modular system, we want to keep it g*.

  subroutine Lrows(H,m0,iDt,iRx,L)

	! uses: grid
    type (solnVector_t), intent(in)                 :: H
    type (modelParam_t), intent(in)                 :: m0   ! not needed?
	integer, intent(in)					            :: iDt,iRx
	type (sparseVector_t), intent(inout)            :: L(:) ! nFunc = 1 or 3
	type (sparsevecc)								:: Lx,Ly,Lz
	complex(8)										:: Hx,Hy,Hz
	complex(8), dimension(3)					    :: pd_Hx,pd_Hy,pd_Hz
	type (sparsevecc)								:: gc_sparse  ! g*
	type (sparsevecc)                  				:: g_sparse	! g
    type (functional_t), pointer                    :: dataType
    type (receiver_t), pointer                      :: obs
    complex(8)                                      :: avg
	real(8)											:: davg,EARTH_R

	dataType => TFList%info(iDt)
	obs => obsList%info(iRx)

	if (.not.obs%located) then
	  write(0,*) 'Error: Observatory ',trim(obs%code),' not yet located in Lrows'
	  stop
	  !call LocateReceiver(grid,obs)
	end if

	if (.not.obs%defined) then
	  write(0,*) 'Error: Observatory ',trim(obs%code),' is not defined at Lrows'
	  stop
	end if


	Lx = obs%Lx
	Ly = obs%Ly
	Lz = obs%Lz

	Hx = dotProd_noConj(Lx,H%vec) ! dotProd_noConj_scvector_f...could do real here
	Hy = dotProd_noConj(Ly,H%vec)
	Hz = dotProd_noConj(Lz,H%vec)

	EARTH_R = H%grid%z(H%grid%nzAir+1) * m2km

	if (dataType%name == 'C') then

!	  pd_Hx(1) = C_ZERO
!	  pd_Hy(1) = - km2m * (EARTH_R/2) * dtan(obs%colat*d2r) * Hz/(Hy*Hy)
!	  pd_Hz(1) =   km2m * (EARTH_R/2) * dtan(obs%colat*d2r) * 1/Hy
	  pd_Hx(1) = C_ZERO
	  pd_Hy(1) = - km2m * (EARTH_R/2) * Hz/(Hy*Hy)
	  pd_Hz(1) =   km2m * (EARTH_R/2) * 1/Hy

	  call linComb_sparsevecc(Ly,pd_Hy(1),Lz,pd_Hz(1),gc_sparse)
      call create_sparseVector(H%grid,H%tx,L(1))
      L(1)%L = gc_sparse ! NOT g_sparse = conjg(gc_sparse)

	else if (dataType%name == 'D') then

!	  pd_Hx(1) =   km2m * (EARTH_R/2) * dsin(obs%colat*d2r) * 1/Hy
!	  pd_Hy(1) = - km2m * (EARTH_R/2) * dsin(obs%colat*d2r) * Hx/(Hy*Hy)
!	  pd_Hz(1) = C_ZERO
	  pd_Hx(1) =   km2m * (EARTH_R/2) * 1/Hy
	  pd_Hy(1) = - km2m * (EARTH_R/2) * Hx/(Hy*Hy)
	  pd_Hz(1) = C_ZERO

	  call linComb_sparsevecc(Ly,pd_Hy(1),Lx,pd_Hx(1),gc_sparse)
      call create_sparseVector(H%grid,H%tx,L(1))
      L(1)%L = gc_sparse ! NOT g_sparse = conjg(gc_sparse)

    else if (dataType%name == 'T') then

      avg = refavg(H%vec)
      davg = drefavg(obs)

      ! derivative of scaled Hx
      pd_Hx(1) =   1/avg
      pd_Hy(1) = - Hx/(avg**2) * davg
      pd_Hz(1) = C_ZERO

      call linComb_sparsevecc(Ly,pd_Hy(1),Lx,pd_Hx(1),gc_sparse)
      call create_sparseVector(H%grid,H%tx,L(1))
      L(1)%L = gc_sparse ! NOT g_sparse = conjg(gc_sparse)

      ! derivative of scaled Hy
      pd_Hx(2) = C_ZERO
      pd_Hy(2) = (1/avg) - Hy/(avg**2) * davg
      pd_Hz(2) = C_ZERO

      call linComb_sparsevecc(Ly,pd_Hy(2),Ly,C_ZERO,gc_sparse)
      call create_sparseVector(H%grid,H%tx,L(2))
      L(2)%L = gc_sparse ! NOT g_sparse = conjg(gc_sparse)

      ! derivative of scaled Hz
      pd_Hx(3) = C_ZERO
      pd_Hy(3) = - Hz/(avg**2) * davg
      pd_Hz(3) =   1/avg

      call linComb_sparsevecc(Ly,pd_Hy(3),Lz,pd_Hz(3),gc_sparse)
      call create_sparseVector(H%grid,H%tx,L(3))
      L(3)%L = gc_sparse ! NOT g_sparse = conjg(gc_sparse)

	else if (dataType%name == 'H') then

	  ! derivative of raw Hx
      pd_Hx(1) = C_ONE
      pd_Hy(1) = C_ZERO
      pd_Hz(1) = C_ZERO

      call linComb_sparsevecc(Lx,pd_Hx(1),Lx,C_ZERO,gc_sparse)
      call create_sparseVector(H%grid,H%tx,L(1))
      L(1)%L = gc_sparse ! NOT g_sparse = conjg(gc_sparse)

      ! derivative of raw Hy
      pd_Hx(2) = C_ZERO
      pd_Hy(2) = C_ONE
      pd_Hz(2) = C_ZERO

      call linComb_sparsevecc(Ly,pd_Hy(2),Ly,C_ZERO,gc_sparse)
      call create_sparseVector(H%grid,H%tx,L(2))
      L(2)%L = gc_sparse ! NOT g_sparse = conjg(gc_sparse)

      ! derivative of raw Hz
      pd_Hx(3) = C_ZERO
      pd_Hy(3) = C_ONE
      pd_Hz(3) = C_ZERO

      call linComb_sparsevecc(Lz,pd_Hz(3),Lz,C_ZERO,gc_sparse)
      call create_sparseVector(H%grid,H%tx,L(3))
      L(3)%L = gc_sparse ! NOT g_sparse = conjg(gc_sparse)

	else

	  write(0, *) 'Error: (Lrows) unknown data functional', dataType%name
	  stop

	end if

	call deall_sparsevecc(Lx)
	call deall_sparsevecc(Ly)
	call deall_sparsevecc(Lz)
	call deall_sparsevecc(gc_sparse)
    call deall_sparsevecc(g_sparse)

  end subroutine Lrows	! Lrows

!
!****************************************************************************
  subroutine Qrows(e0,m0,iDt,iRx,zeroValued,Qreal,Qimag)
  !  given input background solution vector (e0) and model parameter (Sigma0)
  !  and indices into data type and receiver dictionaries
  !  compute derivative of data functional with respect to model parameters
  !  for all components of the data type ...
  !             (ZERO VECTORS FOR EARTH!!!!)

  type (solnVector_t), intent(in)       :: e0
  type (modelParam_t), intent(in)       :: m0
  integer, intent(in)                   :: iDt, iRx
  logical, intent(out)                  :: zeroValued
  !   NOTE: Qreal and Qimag have to be declared as arrays for
  !     consistency with calling program (in general the
  !     number nFunc of complex data functionals that will
  !     be returned by data functional routines could be > 1)
  !   NOTE: Qreal and Qimag both exist regardless of whether the data
  !     are real or complex, since Q itself is complex
  type(modelParam_t), intent(inout)     :: Qreal(:), Qimag(:)

  ! local variables
  type(functional_t), pointer :: dataType
  type(receiver_t), pointer   :: obs
  integer       :: istat,nComp,iComp
  logical       :: isComplex

  dataType => TFList%info(iDt)
  obs => obsList%info(iRx)

  ! set the rows of Q to zero
  if((size(Qreal) .ne. dataType%nComp) .or. (size(Qimag) .ne. dataType%nComp)) then
    call errStop('incorrect output size in Qrows')
  endif

  ! for efficiency, if vectors are zero just set the logical to true and exit
  zeroValued = .true.

  !do iComp = 1, dataType%nComp
  !  Qreal(iComp) = m0
  !  call zero(Qreal(iComp))
  !  Qimag(iComp) = m0
  !  call zero(Qimag(iComp))
  !enddo

  end subroutine Qrows

end module DataFunc
