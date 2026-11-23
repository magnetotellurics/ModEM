module DistortionParam
  
  ! defines derived type that holds distortion matrices for each site in a 
  ! receiver dictionary, with basic operations on them (dot product, linear
  ! combination, etc.) analogue to the model parameter type. 
  ! the distortion matrice is a 2x2 tensor C for each site, where:
  ! C = g * T * S * A, where g is the scalar gain, T is the twist matrix, 
  ! S is the shear matrix, and A is the anisotropy matrix.
  use math_constants
  use utilities

  implicit none

  public :: distortionParam_t, &
            create_distortionParam, deall_distortionParam, &
            zero_distortionParam, identity_distortionParam, &
            dotProd_distortionParam, linComb_distortionParam, &
            scMult_distortionParam, copy_distortionParam, &
            read_distortionParam, write_distortionParam

  type :: distortionParam_t
     integer                                :: nSites = 0
     real(kind=prec), allocatable           :: C(:,:,:)   ! (2,2,nSites)
     integer, allocatable                   :: siteIndex(:) ! rxDict indices
  end type distortionParam_t

Contains

  subroutine create_distortionParam(nSites, D)
    integer, intent(in) :: nSites
    type(distortionParam_t), intent(out) :: D
    integer :: istat
    D%nSites = nSites
    allocate(D%C(2,2,nSites), D%siteIndex(nSites), STAT=istat)
    if (istat /= 0) call errStop('allocation error in create_distortionParam')
    call identity_distortionParam(D)
    D%siteIndex = 0
  end subroutine create_distortionParam

  subroutine deall_distortionParam(D)
    type(distortionParam_t), intent(inout) :: D
    integer :: istat
    if (allocated(D%C)) deallocate(D%C, STAT=istat)
    if (allocated(D%siteIndex)) deallocate(D%siteIndex, STAT=istat)
    D%nSites = 0
  end subroutine deall_distortionParam

  subroutine zero_distortionParam(D)
    type(distortionParam_t), intent(inout) :: D
    integer :: i
    if (.not. allocated(D%C)) return
    do i = 1, D%nSites
       D%C(:,:,i) = R_ZERO
    end do
  end subroutine zero_distortionParam

  subroutine identity_distortionParam(D)
    type(distortionParam_t), intent(inout) :: D
    integer :: i
    if (.not. allocated(D%C)) return
    do i = 1, D%nSites
       D%C(1,1,i) = ONE
       D%C(1,2,i) = R_ZERO
       D%C(2,1,i) = R_ZERO
       D%C(2,2,i) = ONE
    end do
  end subroutine identity_distortionParam

  function dotProd_distortionParam(D1, D2) result(dp)
    type(distortionParam_t), intent(in) :: D1, D2
    real(kind=prec) :: dp
    integer :: i, j
    dp = R_ZERO
    if (D1%nSites /= D2%nSites) call errStop('size mismatch in dotProd_distortionParam')
    do i = 1, D1%nSites
       do j = 1, 2
          dp = dp + D1%C(1,j,i)*D2%C(1,j,i) + D1%C(2,j,i)*D2%C(2,j,i)
       end do
    end do
  end function dotProd_distortionParam

  subroutine linComb_distortionParam(a, D1, b, D2, Dout)
    real(kind=prec), intent(in) :: a, b
    type(distortionParam_t), intent(in) :: D1, D2
    type(distortionParam_t), intent(inout) :: Dout
    integer :: i, j
    if (D1%nSites /= D2%nSites .or. D1%nSites /= Dout%nSites) &
         call errStop('size mismatch in linComb_distortionParam')
    do i = 1, D1%nSites
       do j = 1, 2
          Dout%C(1,j,i) = a*D1%C(1,j,i) + b*D2%C(1,j,i)
          Dout%C(2,j,i) = a*D1%C(2,j,i) + b*D2%C(2,j,i)
       end do
    end do
    Dout%siteIndex = D1%siteIndex
  end subroutine linComb_distortionParam

  subroutine scMult_distortionParam(a, D_in, D_out)
    real(kind=prec), intent(in) :: a
    type(distortionParam_t), intent(in) :: D_in
    type(distortionParam_t), intent(inout) :: D_out
    integer :: i, j
    if (D_in%nSites /= D_out%nSites) &
         call errStop('size mismatch in scMult_distortionParam')
    do i = 1, D_in%nSites
       do j = 1, 2
          D_out%C(1,j,i) = a * D_in%C(1,j,i)
          D_out%C(2,j,i) = a * D_in%C(2,j,i)
       end do
    end do
    D_out%siteIndex = D_in%siteIndex
  end subroutine scMult_distortionParam

  subroutine copy_distortionParam(D_in, D_out)
    type(distortionParam_t), intent(in) :: D_in
    type(distortionParam_t), intent(inout) :: D_out
    integer :: i, istat
    if (D_out%nSites /= D_in%nSites) then
       call deall_distortionParam(D_out)
       call create_distortionParam(D_in%nSites, D_out)
    end if
    do i = 1, D_in%nSites
       D_out%C(:,:,i) = D_in%C(:,:,i)
    end do
    D_out%siteIndex = D_in%siteIndex
  end subroutine copy_distortionParam

  subroutine read_distortionParam(fname, D)
    ! Read distortion matrices into an existing D (D%siteIndex set by the
    ! caller).  Entries are matched to D by receiver index (site_id), so the
    ! file may list the sites in any order; sites of D missing from the file
    ! keep their current value.
    character(*), intent(in) :: fname
    type(distortionParam_t), intent(inout) :: D
    integer :: i, k, nSites, iostat, iSite, iRec, idxDummy, nMatched
    real(kind=prec) :: c11, c12, c21, c22
    character(256) :: line
    character(80) :: str
    open(unit=77, file=fname, status='old', iostat=iostat)
    if (iostat /= 0) call errStop('Cannot open distortion file: '//trim(fname))
    read(77, *, iostat=iostat) str, nSites
    if (iostat /= 0) call errStop('Error reading distortion file header: '//trim(fname))
    nMatched = 0
    do i = 1, nSites
       read(77, '(A)', iostat=iostat) line
       if (iostat /= 0) call errStop('Error reading distortion line: '//trim(fname))
       ! Format: idx, site_id, Cxx, Cxy, Cyx, Cyy
       read(line, *, iostat=iRec) idxDummy, iSite, c11, c12, c21, c22
       if (iRec /= 0) then
          ! Backward-compatible format: site_id, Cxx, Cxy, Cyx, Cyy
          read(line, *, iostat=iRec) iSite, c11, c12, c21, c22
          if (iRec /= 0) then
             call errStop('Invalid distortion line format in: '//trim(fname))
          end if
       end if
       do k = 1, D%nSites
          if (D%siteIndex(k) == iSite) then
             D%C(1,1,k) = c11
             D%C(1,2,k) = c12
             D%C(2,1,k) = c21
             D%C(2,2,k) = c22
             nMatched = nMatched + 1
             exit
          end if
       end do
    end do
    close(77)
    if (nMatched < D%nSites) then
       write(0,'(a,i6,a,i6,a)') ' Warning: distortion file matched ', nMatched, &
            ' of ', D%nSites, ' sites; the others keep their current value'
    end if
  end subroutine read_distortionParam

  subroutine write_distortionParam(fname, D)
    character(*), intent(in) :: fname
    type(distortionParam_t), intent(in) :: D
    integer :: i, iostat
    open(unit=77, file=fname, status='replace', action='write', iostat=iostat)
    if (iostat /= 0) call errStop('Cannot open distortion file for writing: '//trim(fname))
    write(77, *) 'NSITE', D%nSites

    do i = 1, D%nSites
       write(77, '(i5,1x,i7,4es15.7)') i, D%siteIndex(i), D%C(1,1,i), D%C(1,2,i), &
            D%C(2,1,i), D%C(2,2,i)
    end do
    close(77)
  end subroutine write_distortionParam

end module DistortionParam
