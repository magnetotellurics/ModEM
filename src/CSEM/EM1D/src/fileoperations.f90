!------------------------------------------------------------------
! FD EM subroutine checkfileexist
!
! Purpose:    check file exstence and call crash routine if file does not exist
!             display a warning if an optional file does not exist
!
! Rita Streich 2009
!------------------------------------------------------------------
subroutine checkfileexist(filename,required,lexist)

  implicit none

  !external variable
  character(len=namlen),intent(inout) :: filename
  logical,intent(in)                  :: required  !is file definitely needed (true) or optional (false)?
  logical,intent(out)                 :: lexist    !return info if file exists

  lexist = .false.

  if (required) then
    !hard check, crash if file is not present
    inquire(file=trim(adjustl(filename)), exist=lexist)
    if (.not. lexist) call exist_error(pid,trim(adjustl(filename)))
  else
    !soft check, only if filename was not set to 'none'
    if (trim(adjustl(filename)).ne.'none') then
      inquire(file=trim(adjustl(filename)), exist=lexist)
      !set filename to none if file is not present
      if (.not.lexist) then
        if (pid .eq. 0) then
          write(*,fmt='(a26,a)') 'WARNING: cannot find file ',trim(adjustl(filename))
          write(*,fmt='(a)')     'File will not be used.'
        endif
        write(unit=filename,fmt='(a4)') 'none'
      endif
    endif
  endif

endsubroutine checkfileexist

!------------------------------------------------------------------
! INV3D subroutine deletefile
!
! Purpose: deletes a file
!
! Alexander Grayver 2011
!------------------------------------------------------------------
subroutine deletefile(fname)

  implicit none

  !external variables
  character(len=*)      :: fname  !file name

  !internal variables
  integer(kind=int32)   :: ierr   !error index
  
#ifdef USE_MPI
  call MPI_FILE_DELETE(fname, MPI_INFO_NULL,ierr)
#else
#endif

endsubroutine deletefile

!**********************************************************************
!  EM subroutine getfilesize
!
!  Purpose:  get size of a file
!
!  Rita Streich 2009
!  Alexander Grayver 2011
!**********************************************************************
integer(kind=int64) function getfilesize(fname)

  implicit none

  !external variable
  character(len=*),intent(in)  :: fname  !file name
  
  !internal variables
  integer(kind=int64)          :: fsize

  fsize = 0
  inquire(FILE=adjustl(fname),SIZE=fsize)
  
  getfilesize = fsize
  return
      
endfunction getfilesize
