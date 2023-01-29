!---------------------------------------------------------------------
!> FD EM  subroutine
!
!> alloc_error: display message if allocation failed and abort
!
!> Rita Streich 2009
!---------------------------------------------------------------------



subroutine alloc_error(pid,routine,array,err)
 use, intrinsic :: iso_fortran_env
  implicit none
  
  !external variables
  integer(kind=int32),intent( in ) :: pid     !process ID
  character(len=*),intent( in ) :: routine !name of routine where error occurred
  character(len=*),intent( in ) :: array   !name of array for which error occurred
  integer(kind=int32),intent( in ) :: err     !error index returned by allocate command
  
  !internal variable
  character(len=5) :: pidstr    !string to write process id to
#ifdef USE_MPI
  integer(kind=int32) :: mpierr  !error index for call to MPI_Abort
#endif

  !message from all processes where something went wrong
  write(unit=pidstr,fmt='(i5)') pid
  write(*,'(a)') 'Process '//trim(adjustl(pidstr))//': ERROR in routine '//trim(adjustl(routine))//':'
  write(*,'(a)') 'Process '//trim(adjustl(pidstr))//': Cannot allocate array '//trim(adjustl(array))//'!'
  write(*,'(a27,i6)') 'Process '//trim(adjustl(pidstr))//': Error code: ',err
#ifdef USE_MPI
  call MPI_Abort(MPI_COMM_WORLD,3,mpierr)
#endif
  stop 3

endsubroutine alloc_error


!------------------------------------------------------------------------------
!> read/write/open/close error
!------------------------------------------------------------------------------
subroutine io_error(pid,filename,action,err)
    USE, INTRINSIC :: ISO_FORTRAN_ENV
  implicit none

  !external variables
  integer(kind=int32),intent( in ) :: pid       !process ID
  character(len=*),intent( in ) :: filename  !name of file
  character(len=*),intent( in ) :: action    !description of what went wrong
  integer(kind=int32),intent( in ) :: err       !error index returned by I/O routine

  !internal variable
#ifdef USE_MPI
  integer(kind=int32) :: mpierr    !error index for call to MPI_Abort
#endif
  character(len=5) :: pidstr    !string to write process id to

  write(unit=pidstr,fmt='(i5)') pid
  write(*,'(a)') 'pid '//trim(adjustl(pidstr))//': Error '//trim(action)//' file '//trim(adjustl(filename))//'!'
  write( *, * ) 'Error code: ',err
#ifdef USE_MPI
  call MPI_Abort(MPI_COMM_WORLD,-2,mpierr)
#endif
  stop 2
  
endsubroutine io_error


!---------------------------------------------------------------------
!
!> FD EM subroutine
!
!> invalid_error: crash for invalid input
!
!> Rita Streich 2009
!---------------------------------------------------------------------
subroutine invalid_error(pid,routine,filename,var,intnum,realnum)
    USE, INTRINSIC :: ISO_FORTRAN_ENV
  implicit none

  !external variables
  integer(kind=int32),intent( in ) :: pid       !process ID
  character(len=*),intent( in ) :: routine   !name of routine where invalid input was detected
  character(len=*),intent( in ) :: filename  !name of file invalid input was read from
  character(len=*),intent( in ) :: var       !name of variable containing invalid input
  integer(kind=int32),intent( in ),optional :: intnum    !invalid value
  real(kind=real64),intent( in ),optional :: realnum   !invalid value

  !internal variable
  character(len=5) :: pidstr    !string to write process id to
#ifdef USE_MPI
  integer(kind=int32) :: mpierr    !error index for call to MPI_Abort
#endif

  write(unit=pidstr,fmt='(i5)') pid
  write(*,'(a)') 'Process '//trim(adjustl(pidstr))//': ERROR in '//trim(adjustl(routine))//':'
  write(*,'(a)') 'Process '//trim(adjustl(pidstr))//': Found invalid value for '//trim(adjustl(var))//' in file ' &
    //trim(adjustl(filename))//':'
  if(present(intnum)) write(*,'(a24,i10)') 'Value on process '//trim(adjustl(pidstr))//': ',intnum
  if(present(realnum)) write(*,'(a24,g20.12)') 'Value on process '//trim(adjustl(pidstr))//': ', realnum
  write(*,'(a)') 'Exiting!'
#ifdef USE_MPI
  call MPI_Abort(MPI_COMM_WORLD,-4,mpierr)
#endif
  stop 4
  
endsubroutine invalid_error


!---------------------------------------------------------------------
!
!> FD EM subroutine error_mpi
!
!> error_mpi:  display message if file open failed and abort
!>   Cannot call this routine "mpi_error" since that is already an MPI function!
!
!> Rita Streich 2009
!---------------------------------------------------------------------
#ifdef USE_MPI
subroutine error_mpi(pid,routine,mpifunc,ierr)

  implicit none

  !external
  integer(kind=int32),intent( in ) :: pid       !process ID
  character(len=*),intent( in ) :: routine   !name of routine where an MPI call failed
  character(len=*),intent( in ) :: mpifunc   !name of MPI function that failed
  integer(kind=int32),intent( in ) :: ierr      !error index returned from MPI call
  

  !internal
  character(len=5) :: pidstr    !string to write process id to
  character(len=MPI_MAX_ERROR_STRING) :: errstr    !string to store text of the error
  integer(kind=int32) :: reslen
  integer(kind=int32) :: mpierr    !error index for internal MPI call

  write(unit=pidstr,fmt='(i5)') pid
  write(*,'(a)') 'Process '//trim(adjustl(pidstr))//': ERROR in '//trim(adjustl(routine))//&
    ' during execution of MPI routine '//trim(mpifunc)
  call MPI_Error_string(ierr, errstr, reslen, mpierr)
  write(*,'(a12,a)') 'Error text: ', errstr
  
  call MPI_Abort(MPI_COMM_WORLD,-5,mpierr)
  stop

endsubroutine error_mpi
#endif


!---------------------------------------------------------------------
!
!> FD EM  subroutine
!
!> readwrite_error: display message if I/O failed and abort
!
!> Rita Streich 2009
!---------------------------------------------------------------------
subroutine readwrite_error(pid,routine,unitname,action,ierr,irec)
    USE, INTRINSIC :: ISO_FORTRAN_ENV
  implicit none

  !external
  integer(kind=int32),intent( in ) :: pid       !process ID
  character(len=*),intent( in ) :: routine   !name of subroutine where error occurred
  character(len=*),intent( in ) :: unitname  !name of I/O unit (file) where error occurred
  character,intent( in ) :: action    !error for reading (r) or writing (w)?
  integer(kind=int32),intent( in ) :: ierr      !error index from read/write operation
  integer(kind=int64),intent( in ),optional :: irec  !number of record at which I/O failed

  !internal
#ifdef USE_MPI
  integer(kind=int32) :: mpierr    !error index for internal MPI call
#endif
  character(len=5) :: pidstr    !string to write process id to
  character(len=24) :: actionstr !output message string

  if(action.EQ.'w') then
    actionstr = 'write to'
  elseif(action.EQ.'r') then
    actionstr = 'read from'
  else
    actionstr = 'do some I/O on'
  endif

  write(unit=pidstr,fmt='(i5)') pid

  write(*,'(a)') 'Process '//trim(adjustl(pidstr))//': ERROR in '//trim(adjustl(routine))//&
    ' when trying to '//trim(actionstr)//' '//trim(adjustl(unitname))
  write(*,'(a12,i5)') 'Error code: ',ierr
  if(present(irec)) write(*,'(a15,i10)') 'Record number: ',irec
#ifdef USE_MPI
  call MPI_Abort(MPI_COMM_WORLD,-42,mpierr)
#endif
  stop 42

endsubroutine readwrite_error


!---------------------------------------------------------------------
!
!> FD EM subroutine
!
!> open_error: display message if file open failed and abort
!
!> Rita Streich 2009
!---------------------------------------------------------------------
subroutine open_error(pid,routine,filename,ierr)
    USE, INTRINSIC :: ISO_FORTRAN_ENV
  implicit none

  !external
  integer(kind=int32),intent( in ) :: pid       !process ID
  character(len=*),intent( in ) :: routine   !name of subroutine where file opening failed
  character(len=*),intent( in ) :: filename  !name of file for which error occurred
  integer(kind=int32),intent( in ) :: ierr      !error index from open operation

  !internal
#ifdef USE_MPI
  integer(kind=int32) :: mpierr    !error index for internal MPI call
#endif
  character(len=5) :: pidstr    !string to write process id to

  write(unit=pidstr,fmt='(i5)') pid
  write(*,'(a)') 'Process '//trim(adjustl(pidstr))//': ERROR in '//trim(adjustl(routine))//&
    ' when trying to open file'//trim(filename)
  write(*,'(a12,i5)') 'Error code: ',ierr
#ifdef USE_MPI
  call MPI_Abort(MPI_COMM_WORLD,-43,mpierr)
#endif
  stop 43

endsubroutine open_error


!---------------------------------------------------------------------
!
!> FD EM subroutine
!
!> close_error: display message if close failed, but continue
!
!> Rita Streich 2009
!---------------------------------------------------------------------
subroutine close_error(pid,routine,filename,ierr)
    USE, INTRINSIC :: ISO_FORTRAN_ENV
  implicit none

  !external
  integer(kind=int32),intent( in ) :: pid       !process ID
  character(len=*),intent( in ) :: routine   !name of subroutine where file closing failed
  character(len=*),intent( in ) :: filename  !name of file for which error occurred
  integer(kind=int32),intent( in ) :: ierr      !error index from close operation

  !internal
  character(len=5) :: pidstr    !string to write process id to

  write(unit=pidstr,fmt='(i5)') pid
  write(*,'(a)') 'Process '//trim(adjustl(pidstr))//': WARNING!> Routine '//trim(adjustl(routine))//&
    ' returned an error when trying to close file:'
  write(*,'(a)') trim(filename)
  write(*,'(a12,i5)') 'Error code: ',ierr

endsubroutine close_error


!---------------------------------------------------------------------
!
!> FD EM subroutine
!
!> exist_error: crash for non-existing file
!
!> Rita Streich 2009
!---------------------------------------------------------------------
subroutine exist_error(pid,filename)
    USE, INTRINSIC :: ISO_FORTRAN_ENV
  implicit none

  !external
  integer(kind=int32),intent( in ) :: pid       !process ID
  character(len=*),intent( in ) :: filename  !name of non-existing file

  !internal
#ifdef USE_MPI
  integer(kind=int32) :: mpierr    !error index for internal MPI call
#endif
  character(len=5) :: pidstr    !string to write process id to


  write(unit=pidstr,fmt='(i5)') pid
  write(*,'(a)') 'Process '//trim(adjustl(pidstr))//' ERROR:'
  write(*,'(a)') 'File '//trim(adjustl(filename))//' does not exist.'
#ifdef USE_MPI
  call MPI_Abort(MPI_COMM_WORLD,-44,mpierr)
#endif
  stop 44

endsubroutine exist_error


