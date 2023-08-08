!-----------------------------------------------------------------
! FD EM subroutine get_dzsrc:
!   find layer thicknesses above and below source
!
! Rita Streich 2009
!-----------------------------------------------------------------
subroutine get_dzsrc(sz,zbound,ilaysrc,nlay)

  implicit none

  !external variables
  real(kind=real64),intent(in)                :: sz      !source depth
  real(kind=real64),dimension(:),intent(in)   :: zbound  !depths of layer boundaries
  integer(kind=int32),intent(in)              :: ilaysrc !source layer
  integer(kind=int32),intent(in)              :: nlay    !total number of layers

  !internal variables
  integer(kind=int32)    :: ierr  !error index
  integer(kind=int32)    :: ilay  !layer counter


  !transmission through intermediate layers without top & bottom
  !transmission through layers above source is only relevant if source is not in the top (infinite halfspace) layer
  !transm. through layers below source only relevent if src not in bottom infinite halfspace layer

  !these vectors are defined globally in refl. module
  allocate(trans_above_src(2:ilaysrc),trans_below_src(ilaysrc:nlay-1), &
           dz_above_src(2:ilaysrc), dz_below_src(ilaysrc:nlay-1), stat=ierr)
  if(ierr.ne.0) call alloc_error(pid,'get_dzsrc','homogeneous region transmission coefficients',ierr)


  !precompute layer thicknesses - source

  !above/below source
  do ilay=2,ilaysrc-1
    dz_above_src(ilay) = zbound(ilay) - zbound(ilay-1)
  enddo
  if(ilaysrc.gt.1) dz_above_src(ilaysrc) = sz - zbound(ilaysrc-1)
  if(ilaysrc.lt.nlay) dz_below_src(ilaysrc) = zbound(ilaysrc) - sz
  do ilay=ilaysrc+1,nlay-1
    dz_below_src(ilay) = zbound(ilay) - zbound(ilay-1)
  enddo

endsubroutine get_dzsrc


!-----------------------------------------------------------------
! FD EM subroutine get_dzrec:
!   find layer thicknesses above and below receiver
!
! Rita Streich 2009
!-----------------------------------------------------------------
subroutine get_dzrec(sz,zr,zbound,ilayrec,ilaysrc,nlay)

  implicit none

  !external variables
  real(kind=real64),intent(in)                   :: sz,zr  !source and receiver depths
  real(kind=real64),dimension(:),intent(in)      :: zbound !depths of layer boundaries
  integer(kind=int32),intent(in)                 :: ilayrec,ilaysrc !receiver and source layers
  integer(kind=int32),intent(in)                 :: nlay  !total number of layers

  !internal variables
  integer(kind=int32)    :: ierr  !error index
  integer(kind=int32)    :: ilay  !layer counter



  if(sz .ge. zr) then
  !if(sz .gt. zr) then

    !receivers above source: trans_below_rec is from source to receivers
    allocate(trans_below_rec(ilayrec:ilaysrc), trans_above_rec(2:ilayrec), &
             dz_below_rec(ilayrec:ilaysrc), dz_above_rec(2:ilayrec), stat=ierr)
    if(ierr.ne.0) call alloc_error(pid,'get_dzrec','homogeneous region transmission coefficients rec',ierr)


    !layers above receivers - independent of source
    do ilay=2,ilayrec-1
      dz_above_rec(ilay) = zbound(ilay) - zbound(ilay-1)
    enddo
    if(ilayrec .ge. 2) dz_above_rec(ilayrec) = zr - zbound(ilayrec-1)

    !layers below receivers - dependent on source depth
    if(ilayrec .eq. ilaysrc) then
      dz_below_rec(ilayrec) = sz - zr
    else
      dz_below_rec(ilayrec) = zbound(ilayrec) - zr
      do ilay=ilayrec+1,ilaysrc-1
        dz_below_rec(ilay) = zbound(ilay) - zbound(ilay-1)
      enddo
      dz_below_rec(ilaysrc) = sz - zbound(ilaysrc-1)
    endif

  else

    !receivers below source: trans_above_rec is from receivers to source
    allocate(trans_above_rec(ilaysrc:ilayrec), trans_below_rec(ilayrec:nlay-1), &
             dz_above_rec(ilaysrc:ilayrec), dz_below_rec(ilayrec:nlay-1), stat=ierr)
    if(ierr.ne.0) call alloc_error(pid,'get_dzrec','homogeneous region transmission coefficients rec',ierr)


    !layers above receivers - dependent on source depth
    if(ilayrec .eq. ilaysrc) then
      dz_above_rec(ilayrec) = zr - sz
    else
      dz_above_rec(ilaysrc) = zbound(ilaysrc) - sz
      do ilay=ilaysrc+1,ilayrec-1
        dz_above_rec(ilay) = zbound(ilay) - zbound(ilay-1)
      enddo
      dz_above_rec(ilayrec) = zr - zbound(ilayrec-1)
    endif

    !layers below receivers - independent of source
    if(ilayrec .lt. nlay) dz_below_rec(ilayrec) = zbound(ilayrec) - zr
    do ilay=ilayrec+1,nlay-1
      dz_below_rec(ilay) = zbound(ilay) - zbound(ilay-1)
    enddo

  endif

endsubroutine get_dzrec


!-------------------------------------------------------
! FD EM function getnpieces:
!   find number of intervals for numerical integration
!   in a quick & dirty way depending on distance from source
!
! Rita Streich 2009
!-------------------------------------------------------
function getnpieces(r,sz,zr) result(npieces)

  implicit none

  !external variables
  integer(kind=int32)   :: npieces
  real(kind=real64)     :: r      !horizontal radius
  real(kind=real64)     :: sz,zr  !source and receiver depths


  if(r .lt. 1._real64) then
    if(abs(zr-sz).gt.1._real64) then
      npieces = ((1._real64/r) * (log(abs(zr-sz))))
    else
      npieces = int((1._real64/r))
    endif
  else
    npieces = 1
  endif

endfunction getnpieces



!-------------------------------------------------------
! FD EM subroutine prepare_srcdepth
!   some initializations for one particular source depth
!
! Rita Streich 2011
!-------------------------------------------------------
subroutine prepare_srcdepth(sz,omeps_srcv,refl_var,src,izsrc,zbound,omega)

  implicit none

  !external variables
  real(kind=real64),intent(out)            :: sz          !source depth
  complex(kind=real64),intent(out)         :: omeps_srcv   !vertical omeps in source layer
  type(refl_struct),intent(inout)          :: refl_var    !all variables that have to be remembered while computing 1D fields
  type(sorec),intent(in)           :: src         ! source definition
  integer(kind=int32),intent(in)           :: izsrc       !counter for source depths
  real(kind=real64),dimension(:),pointer,intent(in)     :: zbound  !depths of layer boundaries
  real(kind=real64),intent(in)             :: omega       !angular frequency


  sz = refl_var%zsrc(izsrc)
  !find the layer where source depth is in

  ilaysrc = getilay(sz,zbound)

  ! prepare reflection and transmission matrices and get layer thicknesses, using info on source layer
  call get_dzsrc(sz,zbound,ilaysrc,nlay)

  !get source element coordinates for this source element depth
  call extract_srccoord_general(refl_var,src,izsrc)

  !medium property of source layer
  omeps_srcv = omega * epsv(ilaysrc)

endsubroutine prepare_srcdepth


!-------------------------------------------------------
! FD EM subroutine prepare_recdepth
!   some initializations for one particular receiver depth
!
! Rita Streich 2011
!-------------------------------------------------------
subroutine prepare_recdepth(zr,omeps_recv,refl_var,  zr_in,nrecperz,izrec,sz,zbound,aniso,omega,srctype)

  implicit none

  !external variables
  real(kind=real64),intent(out)            :: zr          !receiver depth, checked for near-equality to source depth
  complex(kind=real64),intent(out)         :: omeps_recv   !vertical omeps in receiver layer
  type(refl_struct),intent(inout)          :: refl_var    !container for stuff needed throughout 1D computations
  real(kind=real64),dimension(:),intent(in)             :: zr_in       !receiver depth
    integer(kind=int32),dimension(:),pointer :: nrecperz  !number of receivers at each receiver depth
  integer(kind=int32),intent(in)           :: izrec       !counter for receiver depths
  real(kind=real64),intent(in)             :: sz          !source depth
  real(kind=real64),dimension(:),pointer,intent(in)     :: zbound  !depths of layer boundaries
  integer(kind=int32),intent(in)           :: aniso       !anisotropy index
  real(kind=real64),intent(in)             :: omega       !angular frequency
  integer(kind=int32),intent(in)           :: srctype     !source type: HED etc.


  !internal variables
  real(kind=real64)              :: zd            !difference between source and receiver depths
  real(kind=real64),parameter    :: zdtiny = 1.e-10  !very small difference between source and receiver depth -> set sz = zr


  zr = zr_in(izrec)
  zd = abs((zr-sz)/sz)
  if(zd.lt.zdtiny) zr = sz

  refl_var%irecstart = sum(nrecperz(1:izrec-1)) + 1
  refl_var%irecend = refl_var%irecstart - 1 + nrecperz(izrec)


  ! find the layer where receivers are in
  ilayrec =  getilay(zr,zbound)

  ! get layer thicknesses and prepare reflection and transmission matrices - receivers:
  call get_dzrec(sz,zr,zbound,ilayrec,ilaysrc,nlay)

  !medium property of receiver layer
  omeps_recv = omega * epsv(ilayrec)

  !precompute interface refl. and transm. coeff.
   refl_var%refcoef_changed = .TRUE.
   call prepare_refcoef(refl_var,refl_var%rmax,srctype,aniso)
   refl_var%refcoef_changed = .false.

endsubroutine prepare_recdepth

