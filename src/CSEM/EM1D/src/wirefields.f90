!**********************************************************************
!  FD EM subroutine init_wirefields
!
!  Purpose: allocate and initialize EM field for the separate wires of
!    a wire source, general routine for all instances of calling 1d background
!
!  Rita Streich 2011
!**********************************************************************
subroutine init_wirefields(src,refl_var,bgdat)

  implicit none

  !external variables
  type(sorec),intent(in)   :: src       !source specification
  type(refl_struct)        :: refl_var  !variable collection for 1D field computation
  type(backgrounddata),intent(in)  :: bgdat  !background field and position arrays

  !internal variables
  integer(kind=int32)      :: ierr  !error index
  integer(kind=int32)      :: iwire !wire counter
  integer(kind=int32)      :: ilay  !layer counter


  if ((bgdat%dowhat .eq. fwdmodel) .or. (bgdat%dowhat.eq.fwd_deriv)) then

    allocate(refl_var%EHwire(src%nwire), stat=ierr)
    if (ierr.ne.0) call alloc_error(pid,'init_wirefields','EHwire',ierr)

    do iwire=1,src%nwire
      call init_1wirefield(refl_var%EHwire(iwire),bgdat)
    enddo
  endif

  if (bgdat%dowhat.ge.deriv) then
    allocate(refl_var%EHwirederiv(src%nwire,bgdat%nlay), stat=ierr)
    if (ierr.ne.0) call alloc_error(pid,'init_wirefields','EHwirederiv',ierr)

    do ilay = 1,bgdat%nlay
      do iwire=1,src%nwire
        call init_1wirefield(refl_var%EHwirederiv(iwire,ilay),bgdat)
      enddo
    enddo

    if (bgdat%aniso.eq.vti) then
      allocate(refl_var%EHwirederivv(src%nwire,bgdat%nlay), stat=ierr)
      if (ierr.ne.0) call alloc_error(pid,'init_wirefields','EHwirederivv',ierr)

      do ilay = 1,bgdat%nlay
        do iwire=1,src%nwire
          call init_1wirefield(refl_var%EHwirederivv(iwire,ilay),bgdat)
        enddo
      enddo

    endif
  endif

endsubroutine init_wirefields


!**********************************************************************
!  FD EM subroutine init_1wirefield
!
!  Purpose: allocate and initialize EM field for one separate wires of
!    a wire source, or for derivatives for one wire and one layer
!
!  Rita Streich 2011
!**********************************************************************
subroutine init_1wirefield(EHwire,bgdat)

  implicit none

  !external variables
  type(receiverdata)  :: EHwire  !fields for one wire *or derivatives for one wire and one layer)
  type(backgrounddata),intent(in)  :: bgdat  !background field and position arrays

  !internal variables
  integer(kind=int32)   :: ierr  !error index


      if (bgdat%nExy.gt.0) then
        allocate(EHwire%Ex(bgdat%nExy),EHwire%Ey(bgdat%nExy), stat=ierr)
      else
        allocate(EHwire%Ex(bgdat%nEx),EHwire%Ey(bgdat%nEy), stat=ierr)
      endif
      if (ierr.ne.0) call alloc_error(pid,'init_1wirefield','Ex,Ey',ierr)

      if (bgdat%nHxy.gt.0) then
        allocate(EHwire%Hx(bgdat%nHxy),EHwire%Hy(bgdat%nHxy), stat=ierr)
      else
        allocate(EHwire%Hx(bgdat%nHx),EHwire%Hy(bgdat%nHy), stat=ierr)
      endif
      if (ierr.ne.0) call alloc_error(pid,'init_1wirefield','Hx,Hy',ierr)

      allocate(EHwire%Ez(bgdat%nEz),EHwire%Hz(bgdat%nHz), stat=ierr)
      if (ierr.ne.0) call alloc_error(pid,'init_1wirefield','Ez,Hz',ierr)

      EHwire%Ex = 0._real64
      EHwire%Ey = 0._real64
      EHwire%Ez = 0._real64
      EHwire%Hx = 0._real64
      EHwire%Hy = 0._real64
      EHwire%Hz = 0._real64

endsubroutine init_1wirefield


!**********************************************************************
!  FD EM subroutine clean_wirefields
!
!  Purpose: clean up EM field vectors of iReceiver locations for the separate wires of
!    a wire source
!
!  Rita Streich 2011
!**********************************************************************
subroutine clean_wirefields(src,refl_var,bgdat)

  implicit none

  !external variables
  type(sorec),intent(in)   :: src       !source specification
  type(refl_struct)        :: refl_var  !variable collection for 1D field computation
  type(backgrounddata),intent(in)  :: bgdat  !background field and position arrays

  !internal variables
  integer(kind=int32)      :: ierr  !error index
  integer(kind=int32)      :: iwire !wire counter
  type(receiverdata),pointer  :: fields  
  integer(kind=int32)      :: ilay  !layer counter


  if ((bgdat%dowhat.eq.fwdmodel) .or. (bgdat%dowhat.eq.fwd_deriv)) then
    do iwire=1,src%nwire
      fields => refl_var%EHwire(iwire)
      deallocate(fields%Ex,fields%Ey,fields%Ez,fields%Hx,fields%Hy,fields%Hz, stat=ierr)
    enddo
    deallocate(refl_var%EHwire, stat=ierr)
  endif

  if (bgdat%dowhat.ge.deriv) then
    do ilay=1,bgdat%nlay
      do iwire=1,src%nwire
        fields => refl_var%EHwirederiv(iwire,ilay)
        deallocate(fields%Ex,fields%Ey,fields%Ez,fields%Hx,fields%Hy,fields%Hz, stat=ierr)
      enddo
    enddo
    deallocate(refl_var%EHwirederiv, stat=ierr)

    if (bgdat%aniso.eq.vti) then
      do ilay=1,bgdat%nlay
        do iwire=1,src%nwire
          fields => refl_var%EHwirederivv(iwire,ilay)
          deallocate(fields%Ex,fields%Ey,fields%Ez,fields%Hx,fields%Hy,fields%Hz, stat=ierr)
        enddo
      enddo
      deallocate(refl_var%EHwirederivv, stat=ierr)
    endif
  endif
  nullify(fields)

endsubroutine clean_wirefields

