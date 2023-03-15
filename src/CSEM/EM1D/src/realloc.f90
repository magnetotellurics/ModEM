!**************************************************************
!
!  FD EM functions
!  reallocate_int_1d
!  Rita Streich 2009
!
!**************************************************************

function reallocate_int_1d(trace, newdim)
    implicit none
    integer(kind=int32), pointer, dimension(:) :: trace, reallocate_int_1d
    integer(kind=int32), intent(in) :: newdim
    integer(kind=int32)             :: nold, ierr, i

    allocate(reallocate_int_1d(1:newdim), stat=ierr)
    if (ierr .ne. 0) call alloc_error(pid,'reallocate_int_1d','new vector',ierr)

    if(.not. associated(trace)) return
    nold = min(size(trace), newdim)
    !careful: prevent stack overflow
    do i=1,nold
      reallocate_int_1d(i) = trace(i)
    enddo
    deallocate(trace)
endfunction reallocate_int_1d


!**************************************************************
!
!  FD EM functions
!  reallocate_dcplx_1d
!  Rita Streich 2009
!
!**************************************************************

function reallocate_dcplx_1d(trace, newdim)
    implicit none
    complex(kind=real64), pointer, dimension(:) :: trace, reallocate_dcplx_1d
    integer(kind=int32), intent(in) :: newdim
    integer(kind=int32)             :: nold, ierr, i

    allocate(reallocate_dcplx_1d(1:newdim), stat=ierr)
    if (ierr .ne. 0) call alloc_error(pid,'reallocate_dcplx_1d','new vector',ierr)

    if(.not. associated(trace)) return
    nold = min(size(trace), newdim)
    !careful: prevent stack overflow
    do i=1,nold
      reallocate_dcplx_1d(i) = trace(i)
    enddo
    deallocate(trace)
endfunction reallocate_dcplx_1d

