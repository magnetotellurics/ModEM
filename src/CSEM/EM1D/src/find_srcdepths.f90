!**************************************************************
!>  FD EM subroutine find_srcdepths
!>  Purpose:  get and sort unique source element depths for one 
!>    wire source or composite dipole source
!
!>  Rita Streich 2009-2011
!**************************************************************
subroutine find_srcdepths(src,refl_var,dipoletype)

  implicit none

  !external variables
  type(sorec),intent( in ) :: src           !source definitions for one source
  type(refl_struct) :: refl_var      !all variables that have to be remembered while computing 1D fields
  integer(kind=int32),intent( in ) :: dipoletype    !the type of dipole elements to look for

  !internal variables
  integer(kind=int32) :: ierr       !error index
  integer(kind=Int32) :: nelem      !temp number of source elements
  integer(kind=int32) :: nelemtot   !total number of dipole elements within source, can be the sum over different dipole types
  integer(kind=int32) :: ielem,nzs  !counters
  real(kind=real64),dimension(:),allocatable :: zstmp     !different depths
  integer(kind=int32),dimension(:),allocatable :: nsrcperz  !how many receivers / source elements at each depth
  integer(kind=int32),dimension(:),allocatable :: srcidx    !source element indices
  real(kind=real64),dimension(:),allocatable :: zstmpall  !depths of all source elements considered, can have duplicates
  integer(kind=int32) :: idx   !temp index for remapping source points
  real(kind=real64) :: lx,ly !x,y lengths of wires


  !------------------------------------------------------
  !extract elements for one particular source type
  !------------------------------------------------------

  select case(src%type)
  case(dipole)

    !check input
    if((dipoletype.lt.hed) .OR. (dipoletype.gt.vmd)) &
      call invalid_error(pid,'find_srcdepths','','dipoletype',intnum=dipoletype)

    nelemtot = src%nelem(1)
    allocate(zstmpall(nelemtot),srcidx(nelemtot), stat=ierr)
    if(ierr.NE.0) call alloc_error(pid,'find_srcdepths','zstmpall',ierr)


    nelem = 0

    !composite dipole sources may contain different dipole types
    select case(dipoletype)
    case(hed) !horizontal electric dipole

      do ielem=1,nelemtot
        if((src%ljx(ielem).NE.0._real64) .OR. (src%ljy(ielem).NE.0._real64)) then
          nelem = nelem + 1
          srcidx(nelem) = ielem
          zstmpall(nelem) = src%pos(3,ielem)
        endif
      enddo

      if(nelem .gt. 0) then
        allocate(refl_var%betasrc(1:nelem), stat=ierr)
        if(ierr.NE.0) call alloc_error(pid,'find_srcdepths','betasrc',ierr)
      endif

    case(ved) !vertical electric dipole
      do ielem=1,nelemtot
        if(src%ljz(ielem).NE.0._real64) then
          nelem = nelem + 1
          srcidx(nelem) = ielem
          zstmpall(nelem) = src%pos(3,ielem)
        endif
      enddo

    case(hmd) !horizontal magnetic dipole
      do ielem=1,nelemtot
        if((src%akx(ielem).NE.0._real64) .OR. (src%aky(ielem).NE.0._real64)) then
          nelem = nelem + 1
          srcidx(nelem) = ielem
          zstmpall(nelem) = src%pos(3,ielem)
        endif
      enddo

      if(nelem .gt. 0) then
        allocate(refl_var%betasrc(1:nelem), stat=ierr)
        if(ierr.NE.0) call alloc_error(pid,'find_srcdepths','betasrc',ierr)
      endif

    case(vmd) !vertical magnetic dipole
      do ielem=1,nelemtot
        if(src%akz(ielem).NE.0.d0) then
          nelem = nelem + 1
          srcidx(nelem) = ielem
          zstmpall(nelem) = src%pos(3,ielem)
        endif
      enddo

    end select


  !horizontal wire sources (so far!) contain horizontal electrical source elements only
  !no need to scan through wires (for now!) - should add option for vertical wire components!
  case(wire)

    nelem = src%nwire

    allocate(refl_var%betasrc(1:nelem), stat=ierr)
    if(ierr.NE.0) call alloc_error(pid,'find_srcdepths','betasrc',ierr)

    allocate(zstmpall(nelem), stat=ierr)
    if(ierr.NE.0) call alloc_error(pid,'find_srcdepths','zstmpall',ierr)

    do ielem=1,nelem
      zstmpall(ielem) = src%wire(ielem)%endpos(3,1)
    enddo

  end select


  !------------------------------------------------------
  !if source elements for the given dipole type exist (or we have a wire source),
  !>  then find unique depths of those source elements
  !------------------------------------------------------
  if(nelem .gt. 0) then

    allocate(nsrcperz(nelem),zstmp(nelem), stat=ierr)
    if(ierr.NE.0) call alloc_error(pid,'find_srcdepths','nsrcperz',ierr)

    allocate(refl_var%isrcperz(nelem), stat=ierr)
    if(ierr.NE.0) call alloc_error(pid,'find_srcdepths','refl_var%isrcperz',ierr)


    !get sort indices for source depths
    nzs = 1
  
    !get sort indices for the source elements selected above
    call indexx(nelem,zstmpall,refl_var%isrcperz)


    select case(src%type)
    case(dipole)

      !replace indices with those from original source vector
      !this will only have an effect if there are multiple dipole element types within the source
      do ielem = 1,nelem
        idx = refl_var%isrcperz(ielem)
        refl_var%isrcperz(ielem) = srcidx(idx)
      enddo


      !find different source depths
      zstmp = 0._real64
      zstmp(1) = src%pos(3,refl_var%isrcperz(1))
      nsrcperz(1) = 1

      do ielem=2,nelem

        if(src%pos(3,refl_var%isrcperz(ielem)) .EQ. src%pos(3,refl_var%isrcperz(ielem-1))) then
          nsrcperz(nzs) = nsrcperz(nzs) + 1
        else
          nzs = nzs + 1
          zstmp(nzs) = src%pos(3,refl_var%isrcperz(ielem))
          nsrcperz(nzs) = 1
        endif
      enddo


      !angles for HED sources
      if(dipoletype.EQ.hed) then
        do ielem=1,nelem
          idx = refl_var%isrcperz(ielem)
          refl_var%betasrc(ielem) = atan2(src%ljy(idx),src%ljx(idx))
        enddo

      !angles for HMD sources
      elseif(dipoletype.EQ.hmd) then
        do ielem=1,nelem
          idx = refl_var%isrcperz(ielem)
          refl_var%betasrc(ielem) = atan2(src%aky(idx),src%akx(idx))
        enddo
      endif


    !wire sources
    case(wire)
      zstmp = 0._real64
      zstmp(1) = src%wire(refl_var%isrcperz(1))%endpos(3,1)
      nsrcperz(1) = 1

      do ielem=2,nelem

        if(src%wire(refl_var%isrcperz(ielem))%endpos(3,1) .EQ. src%wire(refl_var%isrcperz(ielem-1))%endpos(3,1)) then
          nsrcperz(nzs) = nsrcperz(nzs) + 1
        else
          nzs = nzs + 1
          zstmp(nzs) = src%wire(refl_var%isrcperz(ielem))%endpos(3,1)
          nsrcperz(nzs) = 1
        endif
      enddo


      !angles between wires and x axis
      do ielem=1,nelem
        idx = refl_var%isrcperz(ielem)
        lx = src%wire(idx)%endpos(1,2) - src%wire(idx)%endpos(1,1)
        ly = src%wire(idx)%endpos(2,2) - src%wire(idx)%endpos(2,1)

        refl_var%betasrc(ielem) = atan2(ly,lx)
      enddo

    end select !source types: dipole or wire


    !copy results into refl_var structure
    refl_var%nzsrc = nzs
    allocate(refl_var%zsrc(nzs),refl_var%nsrcperz(nzs), stat=ierr)
    if(ierr.NE.0) call alloc_error(pid,'find_srcdepths','nzs',ierr)
    refl_var%zsrc = zstmp(1:nzs)
    refl_var%nsrcperz = nsrcperz(1:nzs)


    deallocate(zstmp,nsrcperz, stat=ierr)

  else !no source element for the given source dipole type
    refl_var%nzsrc = 0 !return info that no source elements for the given dipole type are present
  endif

  deallocate(zstmpall, stat=ierr)
  if(src%type .EQ. dipole) deallocate(srcidx, stat=ierr)

endsubroutine find_srcdepths


