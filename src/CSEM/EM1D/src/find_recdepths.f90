!**************************************************************
!  FD EM subroutine find_recdepths_general
!  Purpose:  get unique iReceiver depths and sort receivers (or grid points) by depth
!    there may be different positions for different field components
!
!  Rita Streich 2011
!**************************************************************
subroutine find_recdepths_general(bgdat,refl_var)

  implicit none

  !external variables
  type(backgrounddata),intent(in) :: bgdat         !positions where fields should be computed
  type(refl_struct)               :: refl_var      !all variables that have to be remembered while computing 1D fields

  !internal variables


  if (bgdat%allcomp_samecoord) then
    !same coordinates for all field components: search depths for one component only
    call find_recdepths_1comp(bgdat%nExy,bgdat%Exypos,refl_var%irecperzExy,refl_var%nzrecExy,refl_var%zrecExy,refl_var%nrecperzExy)
    refl_var%irecperzEz => refl_var%irecperzExy
    refl_var%irecperzHxy => refl_var%irecperzExy
    refl_var%irecperzHz => refl_var%irecperzExy
    refl_var%nzrecEz = refl_var%nzrecExy
    refl_var%nzrecHxy = refl_var%nzrecExy
    refl_var%nzrecHz = refl_var%nzrecExy
    refl_var%zrecEz => refl_var%zrecExy
    refl_var%zrecHxy => refl_var%zrecExy
    refl_var%zrecHz => refl_var%zrecExy
    refl_var%nrecperzEz => refl_var%nrecperzExy
    refl_var%nrecperzHxy => refl_var%nrecperzExy
    refl_var%nrecperzHz => refl_var%nrecperzExy
  else
    !coordinates differ for different field components: search for depths for each comp. individually

    !number of receivers: common Exy points or Ex points or Ey points
    !permit only ONE set of depths of Ex and Ey
    !for computations on staggered grid, x,y differs for Ex and Ey but z is sthe same
    !for receivers, x,y,z are the same for Ex and Ey
    !at half-interpolated points, we run bg computations separately for each component --> we have EITHER Ex OR Ey
    !--> a mix of different depths for Ex and Ey is not supposed to occur!
    if (bgdat%nExy .gt. 0) then
      call find_recdepths_1comp(bgdat%nExy,bgdat%Exypos,refl_var%irecperzExy,refl_var%nzrecExy,refl_var%zrecExy,refl_var%nrecperzExy)
    elseif (bgdat%nEx .gt. 0) then
      call find_recdepths_1comp(bgdat%nEx,bgdat%Expos,refl_var%irecperzExy,refl_var%nzrecExy,refl_var%zrecExy,refl_var%nrecperzExy)
    else
      call find_recdepths_1comp(bgdat%nEy,bgdat%Eypos,refl_var%irecperzExy,refl_var%nzrecExy,refl_var%zrecExy,refl_var%nrecperzExy)
    endif
	
    if (bgdat%nEz .gt. 0) then
      call find_recdepths_1comp(bgdat%nEz,bgdat%Ezpos,refl_var%irecperzEz,refl_var%nzrecEz,refl_var%zrecEz,refl_var%nrecperzEz)
    end if
	
    if (bgdat%nHxy .gt. 0) then
      call find_recdepths_1comp(bgdat%nHxy,bgdat%Hxypos,refl_var%irecperzHxy,refl_var%nzrecHxy,refl_var%zrecHxy,refl_var%nrecperzHxy)
    elseif (bgdat%nHx .gt. 0) then
      call find_recdepths_1comp(bgdat%nHx,bgdat%Hxpos,refl_var%irecperzHxy,refl_var%nzrecHxy,refl_var%zrecHxy,refl_var%nrecperzHxy)
    elseif (bgdat%nHy .gt. 0) then
      call find_recdepths_1comp(bgdat%nHy,bgdat%Hypos,refl_var%irecperzHxy,refl_var%nzrecHxy,refl_var%zrecHxy,refl_var%nrecperzHxy)
    endif
	
    if (bgdat%nHz .gt. 0) then
     call find_recdepths_1comp(bgdat%nHz,bgdat%Hzpos,refl_var%irecperzHz,refl_var%nzrecHz,refl_var%zrecHz,refl_var%nrecperzHz)
	end if
	
  endif
endsubroutine find_recdepths_general


!**************************************************************
!  FD EM subroutine find_recdepths_1comp
!  Purpose:  get unique iReceiver depths and sort receivers (or grid points) by depth
!    for one field component
!
!  Rita Streich 2011
!**************************************************************
subroutine find_recdepths_1comp(nrec,pos,irecperz,nzrec,zrec,nrecperz_out)

  implicit none

  !external variables
  integer(kind=int32),intent(in)   :: nrec   !number of points at which fields have to be computed
  real(kind=real64),dimension(:,:), intent(in)  :: pos   !positions (nrec,3)
  integer(kind=int32),dimension(:),pointer,intent(out)     :: irecperz   !indices for receivers at each depth
  integer(kind=int32),intent(out)              :: nzrec  !number of unique iReceiver depths
  real(kind=real64),dimension(:),pointer       :: zrec   !unique iReceiver depths
  integer(kind=int32),dimension(:),pointer     :: nrecperz_out  !number of receivers at each iReceiver depth

  !internal variables
  integer(kind=int32)   :: ierr       !error index
  integer(kind=int32)   :: irec,nzr   !counters
  real(kind=real64),dimension(:),allocatable    :: zrtmp     !different depths
  integer(kind=int32),dimension(:),allocatable  :: nrecperz  !how many receivers at each depth


  allocate(irecperz(nrec), stat=ierr)
  if (ierr.ne.0) call alloc_error(pid,'find_recdepths_1comp','irecperz',ierr)

  !temp depth vector
  allocate(zrtmp(nrec),nrecperz(nrec), stat=ierr)
  if (ierr.ne.0) call alloc_error(pid,'find_recdepths_1comp','zrtmp,nrecperz',ierr)

  !get sort indices for iReceiver depths
  call indexx(nrec,pos(:,3),irecperz)

  if (nrec .eq. 0) then
    nzr = 0
  else

    !find different iReceiver depths
    zrtmp = 0._real64
    nzr = 1
    zrtmp(1) = pos(irecperz(1),3)
    nrecperz(1) = 1
    do irec=2,nrec
      if (pos(irecperz(irec),3) .eq. pos(irecperz(irec-1),3)) then
        nrecperz(nzr) = nrecperz(nzr) + 1
      else
        nzr = nzr + 1
        zrtmp(nzr) = pos(irecperz(irec),3)
        nrecperz(nzr) = 1
      endif
    enddo

  endif

  nzrec = nzr
  allocate(zrec(nzrec),nrecperz_out(nzrec), stat=ierr)
  if (ierr.ne.0) call alloc_error(pid,'find_recdepths_1comp','zrec,nzrec',ierr)
  zrec = zrtmp(1:nzrec)
  nrecperz_out = nrecperz(1:nzrec)

  deallocate(zrtmp,nrecperz, stat=ierr)

endsubroutine find_recdepths_1comp


!**************************************************************
!  FD EM subroutine clean_zrec
!  Purpose:  get unique iReceiver depths and sort receivers (or grid points) by depth
!    for one field componentdeallocate vectors for iReceiver depth sorting
!
!  Rita Streich 2011
!**************************************************************
subroutine clean_zrec(refl_var,bgdat)

  implicit none

  !external variables
  type(refl_struct)               :: refl_var      !all variables that have to be remembered while computing 1D fields
  type(backgrounddata),intent(in) :: bgdat         !positions where fields should be computed

  !internal variables
  integer(kind=int32)    :: ierr  !error index


  
  if (bgdat%allcomp_samecoord) then
    nullify(refl_var%irecperzEz,refl_var%zrecEz,refl_var%nrecperzEz)
    nullify(refl_var%irecperzHxy,refl_var%zrecHxy,refl_var%nrecperzHxy)
    nullify(refl_var%irecperzHz,refl_var%zrecHz,refl_var%nrecperzHz)
    deallocate(refl_var%irecperzExy,refl_var%zrecExy,refl_var%nrecperzExy, stat=ierr)
  else
    deallocate(refl_var%irecperzExy,refl_var%zrecExy,refl_var%nrecperzExy, stat=ierr)
    deallocate(refl_var%irecperzEz,refl_var%zrecEz,refl_var%nrecperzEz, stat=ierr)
    !deallocate(refl_var%irecperzHxy,refl_var%zrecHxy,refl_var%nrecperzHxy, stat=ierr)
    !deallocate(refl_var%irecperzHz,refl_var%zrecHz,refl_var%nrecperzHz, stat=ierr)
  endif

endsubroutine clean_zrec

