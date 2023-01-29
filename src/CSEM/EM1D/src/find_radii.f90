!******************************************************************************
!>  FD EM subroutine find_radii_global
!>    find overallmin/max radii for all source/receiver depths, all field components
!>    general routine for all instances of calling 1d bg computation
!
!>  Rita Streich 2011
!******************************************************************************
subroutine find_radii_global(refl_var,src,bgdat,comm)

  implicit none

  !external variables
  type(refl_struct) :: refl_var   !container for stuff needed throughout 1D computations
  type(sorec),intent( in ) :: src        !source specifications
  type(backgrounddata),intent( in ) :: bgdat      !coordinates and arrays for background fields
  integer(kind=int32),intent( in ) :: comm       !MPI communicator

  !internal variables
  integer(kind=int32) :: izrmax          !index for maximum receiver depth to check
  real(kind=real64) :: rmin,rmax       !temp radii


  !starting radii
  rmax = -1._real64
  rmin = 1.e35_real64

  !2x2 matrix of different cases:
  !same coordinates for all field components
  if(bgdat%allcomp_samecoord) then

    !same horizontal coordinates for all receiver (grid cell) depths?
    if(bgdat%allzrec_samecoord) then
      izrmax = 1
    else
      izrmax = refl_var%nzrecExy
    endif
    call update_maxminrad(src,bgdat%Exypos,refl_var%nrecperzExy,refl_var%irecperzExy,izrmax,rmax,rmin)

  else

    !radii for Ex and Ey
    if(bgdat%allzrec_samecoord) then
      izrmax = 1
    else
      izrmax = refl_var%nzrecExy
    endif
    have_commonExy: if(bgdat%nExy.gt.0) then
      call update_maxminrad(src,bgdat%Exypos,refl_var%nrecperzExy,refl_var%irecperzExy,izrmax,rmax,rmin)
    else

      if(bgdat%nEx.gt.0) then
        call update_maxminrad(src,bgdat%Expos,refl_var%nrecperzExy,refl_var%irecperzExy,izrmax,rmax,rmin)
      endif
      if(bgdat%nEy.gt.0) then
        call update_maxminrad(src,bgdat%Eypos,refl_var%nrecperzExy,refl_var%irecperzExy,izrmax,rmax,rmin)
      endif
    endif have_commonExy

    !radii for Ez
    if(bgdat%allzrec_samecoord) then
      izrmax = 1
    else
      izrmax = refl_var%nzrecEz
    endif
    have_Ez: if(bgdat%nEz.gt.0) then
      call update_maxminrad(src,bgdat%Ezpos,refl_var%nrecperzEz,refl_var%irecperzEz,izrmax,rmax,rmin)
    endif have_Ez

    !radii for Hx and Hy
    if(bgdat%allzrec_samecoord) then
      izrmax = 1
    else
      izrmax = refl_var%nzrecHxy
    endif
    have_commonHxy: if(bgdat%nHxy.gt.0) then
      call update_maxminrad(src,bgdat%Hxypos,refl_var%nrecperzHxy,refl_var%irecperzHxy,izrmax,rmax,rmin)
    else

      if(bgdat%nHx.gt.0) then
        call update_maxminrad(src,bgdat%Hxpos,refl_var%nrecperzHxy,refl_var%irecperzHxy,izrmax,rmax,rmin)
      endif
      if(bgdat%nHy.gt.0) then
        call update_maxminrad(src,bgdat%Hypos,refl_var%nrecperzHxy,refl_var%irecperzHxy,izrmax,rmax,rmin)
      endif
    endif have_commonHxy

    !radii for Hz
    if(bgdat%allzrec_samecoord) then
      izrmax = 1
    else
      izrmax = refl_var%nzrecHz
    endif
    have_Hz: if(bgdat%nHz.gt.0) then
      call update_maxminrad(src,bgdat%Hzpos,refl_var%nrecperzHz,refl_var%irecperzHz,izrmax,rmax,rmin)
    endif have_Hz

  endif

  !remember min and max radii for later Fast Hankel Transform
  !not every receiver depth will be present on all processes --> this routine is not executed by all procs simultaneously
  !--> quick & dirty: do not exchange min/max radius info
  !BUT: this causes little jumps of field values across domain boundaries
  call store_rmaxmin(refl_var,rmax,rmin, comm)

endsubroutine find_radii_global


!******************************************************************************
!>  FD EM subroutine update_maxminrad
!>    find maximum / minimum radii for 1d integral evaluations
!>    use original source definitions without pre-sorting
!>      (no benefit from sorting since we have to look at radii for all source elements anynway)
!>    receiver coordinates are pre-sorted by depth for higher efficiency
!
!>  Rita Streich 2011
!******************************************************************************
subroutine update_maxminrad(src,pos,nrecperz,recidx,izrmax,rmax,rmin)

  implicit none

  !external variables
  type(sorec),intent( in ) :: src     !source specifications
  real(kind=real64),dimension(:,:),intent( in ) :: pos     !positions where both Ex and Ey are computed 
  integer(kind=int32),dimension(:),intent( in ) :: nrecperz  !nr of points for each receiver depth
  integer(kind=int32),dimension(:),intent( in ) :: recidx       !indices for receivers at each depth
  integer(kind=int32) :: izrmax     !index for max. receiver depth to consider
  real(kind=real64),intent(inout) :: rmax,rmin  !maximum and minimum radii

  !internal variables
  integer(kind=int32) :: izr        !receiver depth index
  integer(kind=int32) :: ielem      !sorce element counter
  real(kind=real64) :: x,y                !source element - receiver distance
  integer(kind=int32) :: irecstart,irecend,irec  !receiver indices
  real(kind=real64) :: rtmp               !temp radius
  integer(kind=int32) :: iw         !wire counter


  !have to go through all source elements:
  !--> for dipole source: all dipole elements, independent of source type
  !--> for wires: all wire elements and wire end points
  select case(src%type)
  case(dipole)

    dipole_elems: do ielem = 1,src%nelem(1)
      recdepths: do izr = 1,izrmax

        irecstart = sum(nrecperz(1:izr-1)) + 1
        irecend = irecstart - 1 + nrecperz(izr)

        recs: do irec=irecstart,irecend

          y = pos(recidx(irec),2) - src%pos(2,ielem)
          x = pos(recidx(irec),1) - src%pos(1,ielem)

          rtmp = sqrt(x**2 + y**2)

          if(rtmp .NE. 0._real64) then
            if(rtmp .lt. rmin) rmin = rtmp
            if(rtmp .gt. rmax) rmax = rtmp
          endif

        enddo recs
      enddo recdepths
    enddo dipole_elems

  case(wire)

    wires: do iw=1,src%nwire
      wire_elems: do ielem = 1,src%wire(iw)%nelem
        recdepthsw: do izr = 1,izrmax

          irecstart = sum(nrecperz(1:izr-1)) + 1
          irecend = irecstart - 1 + nrecperz(izr)

          recsw: do irec=irecstart,irecend

            y = pos(recidx(irec),2) - src%wire(iw)%elempos(ielem,2)
            x = pos(recidx(irec),1) - src%wire(iw)%elempos(ielem,1)

            rtmp = sqrt(x**2 + y**2)

            if(rtmp .NE. 0._real64) then
              if(rtmp .lt. rmin) rmin = rtmp
              if(rtmp .gt. rmax) rmax = rtmp
            endif

          enddo recsw
        enddo recdepthsw
      enddo wire_elems

      wire_endpoints: do ielem=1,2
        recdepthswe: do izr = 1,izrmax

          irecstart = sum(nrecperz(1:izr-1)) + 1
          irecend = irecstart - 1 + nrecperz(izr)

          recswe: do irec=irecstart,irecend

            y = pos(recidx(irec),2) - src%wire(iw)%endpos(2,ielem)
            x = pos(recidx(irec),1) - src%wire(iw)%endpos(1,ielem)

            rtmp = sqrt(x**2 + y**2)

            if(rtmp .NE. 0._real64) then
              if(rtmp .lt. rmin) rmin = rtmp
              if(rtmp .gt. rmax) rmax = rtmp
            endif

          enddo recswe
        enddo recdepthswe
      enddo wire_endpoints

    enddo wires

  end select

endsubroutine update_maxminrad


!******************************************************************************
!>  FD EM subroutine extract_srccoord_general
!>    get source coordinates for one source (element) depth into refl_var structure
!
!>  Rita Streich 2009-2011
!******************************************************************************
subroutine extract_srccoord_general(refl_var,src,izsrc)

  implicit none

  !external variables
  type(refl_struct) :: refl_var   !stuff for 1D computations
  type(sorec),intent( in ) :: src        !source specifications
  integer(kind=int32),intent( in ) :: izsrc   !index for source depths

  !internal variables
  integer(kind=int32) :: isrcstart,isrcend  !start and end index of source elements
  integer(kind=int32) :: nelemall           !number of dipole or wire elements, for wire: sum over all wires of equal depth
  integer(kind=int32) :: ierr       !error index
  integer(kind=int32) :: isrc       !source element counter
  integer(kind=int32) :: idx        !wire index
  integer(kind=int32) :: nelem      !number of wire elements on one wire
  integer(kind=int32) :: ielem,ielemall  !counter for wire elements
  integer(kind=int32) :: nsrc       !number of dipole elements/wires of equal depth


  !dipole element or wire range within vector of all dipole elements or wires
  nsrc = refl_var%nsrcperz(izsrc)
  isrcstart = sum(refl_var%nsrcperz(1:izsrc-1)) + 1
  isrcend = isrcstart - 1 + refl_var%nsrcperz(izsrc)

  !remember these values
  refl_var%isrcstart = isrcstart
  refl_var%isrcend = isrcend


  select case(src%type)
  case(dipole)

    allocate(refl_var%xs(isrcstart:isrcend),refl_var%ys(isrcstart:isrcend), stat=ierr)
    if(ierr.NE.0) call alloc_error(pid,'extract_srccoord_general','xs,ys for dipole',ierr)

    !dipole sources: copy (and re-sort) the relevant positions from src structure
    do isrc = isrcstart,isrcend
      refl_var%xs(isrc) = src%pos(1,refl_var%isrcperz(isrc))
      refl_var%ys(isrc) = src%pos(2,refl_var%isrcperz(isrc))
    enddo

  case(wire)

    !wire sources: compute positions of wire elements

    !get total number of wire elements for this wire depth (can be from more than one wire)
    nelemall = 0
    do isrc = isrcstart,isrcend
      idx = refl_var%isrcperz(isrc)
      nelemall = nelemall + src%nelem(idx)
    enddo

    !vectors for all wire elements
    allocate(refl_var%xs(nelemall),refl_var%ys(nelemall), stat=ierr)
    if(ierr.NE.0) call alloc_error(pid,'extract_srccoord_general','xs,ys for wire',ierr)


    ielemall = 0
    do isrc=isrcstart,isrcend

      idx = refl_var%isrcperz(isrc)
      nelem = src%nelem(idx)

      do ielem = 1,nelem
        ielemall = ielemall + 1

        refl_var%xs(ielemall) = src%wire(idx)%elempos(ielem,1)
        refl_var%ys(ielemall) = src%wire(idx)%elempos(ielem,2)
      enddo

    enddo

    !re-define this for output
    isrcstart = 1
    isrcend = nelemall

  end select !source types: dipole or wire

endsubroutine extract_srccoord_general


!******************************************************************************
!>  FD EM subroutine store_rmaxmin
!>    put min and max radii into refl_var structure and allocate radius vectors
!
!>  Rita Streich 2009
!******************************************************************************
subroutine store_rmaxmin(refl_var,rmax,rmin, comm)

  implicit none

  !external variables
  type(refl_struct) :: refl_var  !stuff for 1D computations
  real(kind=real64),intent( in ) :: rmax,rmin   !min and max radius
  integer(kind=int32),intent( in ) :: comm     !MPI communicator

  !internal variable
  integer(kind=int32) :: ierr     !error index
  real(kind=real64) :: fact     !multiplication factor for radii
  integer(kind=int32) :: nrad     !number of radius values
  real(kind=real64) :: rmaxall  !max radius over all domains
  real(kind=real64) :: rminall  !min radius over all domains


#ifdef USE_MPI
  !get max radius from all domains
  call MPI_Allreduce(rmax,rmaxall,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
  if(ierr.NE.MPI_SUCCESS) call error_mpi(pid,'store_rmaxmin','MPI_Allreduce',ierr)
  call MPI_Allreduce(rmin,rminall,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
  if(ierr.NE.MPI_SUCCESS) call error_mpi(pid,'store_rmaxmin','MPI_Allreduce #2',ierr)
#else
  rmaxall = rmax
  rminall = rmin
#endif

  !fix if there is no positive non-zero radius
  if(rminall.gt.rmaxall) then
    rminall = 1.
    rmaxall = 1.
  endif


  fact = exp(logspace)        !logspace is global in reflectivity_mod
  refl_var%rmax = rmaxall * fact !take one more to prevent edge effects in later spline interpolation (need this???)
  refl_var%rmin = rminall

  
  !number of radius values at logarithmic spacing
  !use minimum non-zero radius here
  !rmax was initiated to -1 and stayed at this value if we only had a zero radius during radius search
  if(rmaxall .gt. 0._real64) then
    nrad = 1 + int(ceiling(10._real64*log(refl_var%rmax/rminall)))  !log is log10
  else
    nrad = 1
  endif

  !remember number of radius values
  refl_var%nrad = nrad

  !radius vectors for integration radii and temp integral values
  !radii are (re-)populated each time during fast Hankel integration
  allocate(refl_var%radlog(nrad),refl_var%intvaltmp(nrad,NRELmax), stat=ierr)
  if(ierr.NE.0) call alloc_error(pid,'find_radii_grid','refl_var%radlog etc.',ierr)

endsubroutine store_rmaxmin

