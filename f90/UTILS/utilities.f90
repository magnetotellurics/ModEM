module utilities


      use math_constants
      implicit none

  ! Variables required for storing the date and time in SECONDS. If used
  ! throughout the program, these make the routine profiling easier
      type  :: timer_t

        private
	real :: rtime = 0.0 ! run time
	real :: stime, etime ! start and end times

      end type timer_t

  !character(80)  :: msg

Contains

  ! **************************************************************************
  subroutine errStop(msg)

    character(*), intent(in)  :: msg
    write(0,'(a9)',advance='no') 'Error: '
    write(0,*) trim(msg)

    call ModEM_abort()

  end subroutine errStop

  subroutine ModEM_abort()

    use, intrinsic :: iso_fortran_env, only: stderr=>error_unit
#ifdef MPI
    use mpi
#endif

    implicit none

#ifdef MPI
    integer :: ierr, error_code
#endif

#ifdef MPI
    call MPI_Abort(MPI_COMM_WORLD, error_code, ierr)
#endif
    stop 1

  end subroutine

  ! **************************************************************************
  subroutine warning(msg)

    character(*), intent(in)  :: msg
    write(0,'(a9)',advance='no') 'Warning: '
    write(0,*) trim(msg)

  end subroutine warning

  ! **************************************************************************
  ! timer utilities: set timer
  subroutine reset_time(timer)

    type(timer_t), intent(inout) :: timer
    ! utility variable
    integer, dimension(8)	     :: tarray

 	! Restart the (portable) clock
	call date_and_time(values=tarray)
    timer%stime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)

  end subroutine reset_time

  ! **************************************************************************
  ! timer utilities: compute elapsed run time in seconds
  function elapsed_time(timer) result (rtime)

    type(timer_t), intent(inout) :: timer
    real                         :: rtime ! run time
    ! utility variable
    integer, dimension(8)	     :: tarray

	call date_and_time(values=tarray)
    timer%etime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)
    rtime = timer%etime - timer%stime

    ! Update total run time
    timer%rtime = timer%rtime + rtime

  end function elapsed_time

  ! **************************************************************************
  ! timer utilities: sum up the elapsed run times in seconds
  function saved_time(timer) result (rtime)

    type(timer_t), intent(inout) :: timer
    real                         :: rtime

    rtime = timer%rtime

  end function saved_time

  ! **************************************************************************
  ! timer utilities: clear saved timing information
  subroutine clear_time(timer)

    type(timer_t), intent(inout) :: timer
    ! utility variable
    integer, dimension(8)	     :: tarray

 	! Restart the (portable) clock
	call date_and_time(values=tarray)
    timer%stime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)

    ! Clear saved run time
    timer%rtime = 0.0

  end subroutine clear_time

  ! **************************************************************************
  function clean(x)
    ! This is a utility routine that provides an expression used to battle
	! against machine error problems. It returns the same real or real(8)
	! as the input, but without the extra digits at the end that are often
	! a cause of wrong comparisons in the if statements. ALWAYS use clean(x)
	! instead of x in an inequality!!!
	! LARGE_REAL is defined in the module math_constants
	! A.K.
    implicit none
    real (kind=prec), intent(in)                   :: x
    real (kind=prec)						       :: clean

	clean = dnint(x*LARGE_REAL)/LARGE_REAL

  end function clean

  ! **************************************************************************
  function nearest_meter(x) result(clean)
    ! This is a utility routine that provides an expression used to battle
	! against machine error problems. Both input and output are values in km.
	! The function rounds the value to the nearest meter. This is useful to
	! ensure that the grid read from a file does not depend on system precision.
	! A.K.
    implicit none
    real (kind=prec), intent(in)                   :: x
    real (kind=prec)						       :: clean

	clean = dnint(x*KM2M)/KM2M

  end function nearest_meter


  ! **************************************************************************
  subroutine find_index(array, xmin, xmax, imin, imax)
    ! For an increasing or decreasing array,
    ! find the minimum and maximum indices
    ! such that xmin <= array(i) < xmax
    ! author: A. Kelbert

    implicit none
    real (kind=prec), dimension(:), intent(in)     :: array
    real (kind=prec), intent(in)                   :: xmin,xmax
    integer, intent(out)                           :: imin,imax
    ! local
    logical                                     :: incr
    integer                                     :: i,n

    ! quick check to see what kind of array this is
    n = size(array)
    if (array(1) <= array(n)) then
        incr = .true.
    else
        incr = .false.
    endif

    if (incr) then
        ! for an increasing array...

        imin = 0
        do i = 1,n
           if(clean(array(i)) .ge. clean(xmin)) then
              imin = i
              exit
           endif
        enddo

        imax = 0
        do i = imin,n
           if(clean(array(i)) .lt. clean(xmax)) then
              imax = i
           endif
        enddo

    else
        ! for a decreasing array...

        imax = 0
        do i = n,1,-1
           if(clean(array(i)) .ge. clean(xmin)) then
              imax = i
              exit
           endif
        enddo

        imin = 0
        do i = imax,1,-1
           if(clean(array(i)) .lt. clean(xmax)) then
              imin = i
           endif
        enddo

    endif

  end subroutine find_index

  ! **************************************************************************
  function minNode(x, xNode) result(ix)
    !  This is a utility routine, used by several data functional
    !  set up routines, and for other interpolation functions
    !  Returns index ix such that  xNode(ix) <= x < xNode(ix+1)
    !  If x is out of range:
    !  x < xNode(1) returns 0; if x> xNode(nx) returns nx
    !  Assumes xNode is strictly increasing; does not check this
    !  NOTE: as presently coded, when xNode is called with center
    !  (face) node positions, this routine will return zero for
    !  the coordinates in the outer half cell nearest the boundary
    !  If evaluation over the complete model domain is to be allowed
    !  a more general interpolation rule will be required.
    !  A.K.: modified to allow input of any size, nx = size(xNode).

    implicit none
    real (kind=prec), intent(in)                   :: x
    real (kind=prec), dimension(:), intent(in)     :: xNode

    integer                                     :: ix
    integer                                     :: i

    ix = size(xNode)
    do i = 1,size(xNode)
       if(clean(xNode(i)) .gt. clean(x)) then
          ix = i-1
          exit
       endif
    enddo

  end function minNode


  ! **************************************************************************
  function maxNode(x, xNode) result(ix)
    !  This is a utility routine, used by several data functional
    !  set up routines, and for other interpolation functions
    !  Returns index ix such that  xNode(ix) <= x < xNode(ix+1)
    !  If x is out of range:
    !  x > xNode(1) returns 0; if x< xNode(nx) returns nx
    !  Assumes xNode is strictly decreasing; does not check this
    !  NOTE: as presently coded, when xNode is called with center
    !  (face) node positions, this routine will return zero for
    !  the coordinates in the outer half cell nearest the boundary
    !  If evaluation over the complete model domain is to be allowed
    !  a more general interpolation rule will be required.
    !  A.K.: modified to allow input of any size, nx = size(xNode).

    implicit none
    real (kind=prec), intent(in)                   :: x
    real (kind=prec), dimension(:), intent(in)     :: xNode

    integer                                     :: ix
    integer                                     :: i

    ix = size(xNode)
    do i = 1,size(xNode)
       if(clean(xNode(i)) .lt. clean(x)) then
          ix = i-1
          exit
       endif
    enddo

  end function maxNode

  ! **************************************************************************
  logical function ismember(n,Nvec)
    ! replicates the corresponding function in Matlab: for an integer array,
    ! outputs true if our integer is in the array, otherwise false.
    ! author: A. Kelbert

    integer                 :: n
    integer, dimension(:)   :: Nvec
    ! local
    integer                 :: i

    ismember = .false.
    do i = 1,size(Nvec)
        if (Nvec(i) == n) then
            ismember = .true.
            return
        end if
    end do

  end function

! *****************************************************************************

      integer function findstr(str1,str2)
      character*(*) str1, str2
!     returns the position of str2 in str1.  Ignores case.
!     returns 0 if str2 not found in str1

      integer i, j, capdif
      logical same

      capdif= ichar('a')-ichar('A')

      do 20 i= 1, len(str1)-len(str2)+1
         do 10 j=1,len(str2)

	      same= str1(i+j-1:i+j-1) .eq. str2(j:j) .or.  &
            'A'.le.str2(j:j) .and. str2(j:j).le.'Z' .and.  &
	     ichar(str1(i+j-1:i+j-1)) .eq. ichar(str2(j:j))+capdif .or.  &
            'a'.le.str2(j:j) .and. str2(j:j).le.'z' .and.  &
	     ichar(str1(i+j-1:i+j-1)) .eq. ichar(str2(j:j)) - capdif

	     if( .not.same) go to 20
10       continue
         findstr=i
         return
20    continue

      findstr=0
      return
      end function findstr



! *****************************************************************************

      integer function begwrd(string,iwrd)
      integer iwrd
      character*(*) string

!     Returns the index of the first non-blank character in the iwrd'th
!     non-blank word (word are seperated by spaces, tabs or commas).
!     Returns len if iwrd'th word is not found. integer i, nword

      logical wasblk
      intrinsic len
      integer  i,nword

      wasblk=.true.

      nword= 0
      do i=1,len(string)
         if( string(i:i).eq.' ' .or.string(i:i).eq.',' .or.  &
	     string(i:i).eq.'  '    )then

!           /* current character is blank */
             wasblk=.true.
	 else
	     if (wasblk) then
		nword= nword + 1
	     endif
	     wasblk= .false.
	     if(nword.eq.iwrd)then
	        begwrd= i
		return
	     end if
	 end if
      enddo

      begwrd = len(string)
      return

      end function begwrd

! *****************************************************************************


      integer function endwrd(string,iwrd)
      integer iwrd
      character*(*) string
!     Returns the index of the last non-blank character in the iwrd'th
!     non-blank word (word are seperated by spaces, tabs or commas).
!     Returns len if iwrd'th word is not found.
      integer i, nword
      logical wasblk
      intrinsic len

      wasblk=.true.
      nword= 0
      do 100 i=1,len(string)
	if( string(i:i).eq.' ' .or.  &
	    string(i:i).eq.',' .or.  &
	    string(i:i).eq.'  '    )then

!          /* current character is blank */
           wasblk=.true.
           if(nword.eq.iwrd) RETURN

	else
           if(wasblk) nword= nword + 1
           wasblk= .false.
           if(nword.eq.iwrd) endwrd= i
	end if
100   continue

      endwrd= len(string)

      return

      end function endwrd


! *****************************************************************************

      SUBROUTINE Lenb(string,length)
      IMPLICIT NONE
      CHARACTER*(*) string
      INTEGER nstr,istr,length

      nstr = len(string)
      DO istr=nstr,1,-1
         IF (string(istr:istr).ne.' ') THEN
	    length = istr
            RETURN
	 ENDIF
      ENDDO
      length = 0


      RETURN

      END Subroutine lenb

  ! **************************************************************************
  ! Naser Meqbel included this function: apparently, it is not supported by
  ! all compilers as an intrinsic
  logical function isnan(a)

        real (kind=prec), intent(in) ::a

        if (a .ne. a) then
        	isnan = .true.
        else
        	isnan = .false.
        end if

  end function isnan

  ! Some Fortran Character Utilities:
  ! See http://gbenthien.net/strings/Strings.pdf  for more information and addtional subroutines
  !#############################################################################################
  subroutine compact(str)

	! Converts multiple spaces and tabs to single spaces; deletes control characters;
	! removes initial spaces.
	character(len=*):: str
	character(len=1):: ch
	character(len=len_trim(str)):: outstr
	integer                     :: lenstr,isp,ich,k,i
	str=adjustl(str)
	lenstr=len_trim(str)
	outstr=' '
	isp=0
	k=0

	do i=1,lenstr
	  ch=str(i:i)
	  ich=iachar(ch)

	  select case(ich)

	    case(9,32)     ! space or tab character
	      if(isp==0) then
	        k=k+1
	        outstr(k:k)=' '
	      end if
	      isp=1

	    case(33:)      ! not a space, quote, or control character
	      k=k+1
	      outstr(k:k)=ch
	      isp=0

	  end select

	end do

	str=adjustl(outstr)

end subroutine compact

subroutine parse(str,delims,args,nargs)

! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
! the delimiters contained in the string 'delims'. Preceding a delimiter in
! 'str' by a backslash (\) makes this particular instance not a delimiter.
! The integer output variable nargs contains the number of arguments found.

	character(len=*) :: str,delims
	character(len=len_trim(str)) :: strsav
	character(len=*),dimension(:) :: args
	integer                     :: lenstr,isp,ich,k,i,nargs,na
	strsav=str
	call compact(str)
	na=size(args)
	do i=1,na
	  args(i)=' '
	end do
	nargs=0
	lenstr=len_trim(str)
	if(lenstr==0) return
	k=0

	do
	   if(len_trim(str) == 0) exit
	   nargs=nargs+1
	   call split(str,delims,args(nargs))
	   call removebksl(args(nargs))
	end do
	str=strsav

end subroutine parse
!**********************************************************************

subroutine split(str,delims,before,sep)

! Routine finds the first instance of a character from 'delims' in the
! the string 'str'. The characters before the found delimiter are
! output in 'before'. The characters after the found delimiter are
! output in 'str'. The optional output character 'sep' contains the
! found delimiter. A delimiter in 'str' is treated like an ordinary
! character if it is preceded by a backslash (\). If the backslash
! character is desired in 'str', then precede it with another backslash.

	character(len=*) :: str,delims,before
	character,optional :: sep
	logical :: pres
	character(1) :: ch
	character    :: cha
	integer                     :: lenstr,isp,ich,k,i,nargs,na,ibsl,iposa,ipos
	pres=present(sep)
	lenstr=len_trim(str)
	if(lenstr == 0) return        ! string str is empty
	k=0
	ibsl=0                        ! backslash initially inactive
	before=' '
	do i=1,lenstr
	   ch=str(i:i)
	   if(ibsl == 1) then          ! backslash active
	      k=k+1
	      before(k:k)=ch
	      ibsl=0
	      cycle
	   end if
	   if(ch == '\\') then          ! backslash with backslash inactive
	      k=k+1
	      before(k:k)=ch
	      ibsl=1
	      cycle
	   end if
	   ipos=index(delims,ch)
	   if(ipos == 0) then          ! character is not a delimiter
	      k=k+1
	      before(k:k)=ch
	      cycle
	   end if
	   if(ch /= ' ') then          ! character is a delimiter that is not a space
	      str=str(i+1:)
	      if(pres) sep=ch
	      exit
	   end if
	   cha=str(i+1:i+1)            ! character is a space delimiter
	   iposa=index(delims,cha)
	   if(iposa > 0) then          ! next character is a delimiter
	      str=str(i+2:)
	      if(pres) sep=cha
	      exit
	   else
	      str=str(i+1:)
	      if(pres) sep=ch
	      exit
	   end if
	end do
	if(i >= lenstr) str=''
	str=adjustl(str)              ! remove initial spaces
	return

end subroutine split

!**********************************************************************

subroutine removebksl(str)

! Removes backslash (\) characters. Double backslashes (\\) are replaced
! by a single backslash.

	character(len=*):: str
	character(len=1):: ch
	character(len=len_trim(str))::outstr
	integer                     :: lenstr,isp,ich,k,i,nargs,na,ibsl,iposa,ipos
	str=adjustl(str)
	lenstr=len_trim(str)
	outstr=' '
	k=0
	ibsl=0                        ! backslash initially inactive

	do i=1,lenstr
	  ch=str(i:i)
	  if(ibsl == 1) then          ! backslash active
	   k=k+1
	   outstr(k:k)=ch
	   ibsl=0
	   cycle
	  end if
	  if(ch == '\\') then          ! backslash with backslash inactive
	   ibsl=1
	   cycle
	  end if
	  k=k+1
	  outstr(k:k)=ch              ! non-backslash with backslash inactive
	end do

	str=adjustl(outstr)

end subroutine removebksl
!**********************************************************************
function is_letter(ch) result(res)

	! Returns .true. if ch is a letter and .false. otherwise
	
	character :: ch
	logical :: res
	
	select case(ch)
	case('A':'Z','a':'z')
	  res=.true.
	case default
	  res=.false.
	end select
	return

end function is_letter
!**********************************************************************
function is_digit(ch) result(res)

	! Returns .true. if ch is a digit (0,1,...,9) and .false. otherwise
	
	character :: ch
	logical :: res
	
	select case(ch)
	case('0':'9')
	  res=.true.
	case default
	  res=.false.
	end select
	return

end function is_digit

!**********************************************************************
recursive subroutine QSort(a, ia, i0, i1)
! a simple recursive quick sort routine using a middle pivot
! the average time complexity is O(nlog(n)), with the worst case of 
! O(n.^2)
  implicit none
  real(kind=prec),intent(inout),dimension(:) :: a
  real(kind=prec)                            :: pivot, t, random
  integer,intent(inout),dimension(:)         :: ia
  integer,intent(in),optional                :: i0, i1
  integer                                    :: first,last,i, j, it, nA
  if (size(a).ne.size(ia)) then
      write(6,*) 'error, array and array index is not of same size in QSort!'
      stop
  endif
  if (.not.present(i0)) then
      first = 1
      last = size(a)
  elseif (.not.present(i1)) then
      first = i0
      last = size(a)
  else 
      first = i0
      last = i1
  endif
  if(first.gt.last) then 
     write(6,*) 'error, first index is larger than the last in QSort!'
     stop
  elseif(first.eq.last) then !only one element, no need to sort now
!     write(6,*) 'no need to sort, only one element left'
     return
  endif
  !Na = j-i+1
  !call random_number(random)
  !write(6,*) Na,int(random*real(Na-1))
  !pivot = a(int(random*real(Na-1))+1)! using a random pivot
  pivot = a((first+last)/2) ! using midpoint
!  write(6,*)  'taking a pivot value of ',int(pivot)
  i = first
  j = last
  do
     do while (a(i).lt.pivot) !left half
        i=i+1
     end do
     do while (pivot.lt.a(j)) !right half
        j=j-1
     end do
     if (i.lt.j) then 
     ! swap array and index
        t = a(i)
        a(i) = a(j)
        a(j) = t
        it = ia(i)
        ia(i) = ia(j)
        ia(j) = it
        i=i+1
        j=j-1
!        write(6,*) 'swapping ai and aj: '
!        write(6,*) int(a)
     else
        exit
     endif
  end do
  if (first.lt.i-1) then 
!      write(6,*)  'taking care of the left part'
      call QSort(a, ia, first, i-1)
  endif
  if (j+1.lt.last) then 
!      write(6,*)  'taking care of the right part'
      call QSort(a, ia, j+1, last)
  endif
  return
end subroutine QSort

!**********************************************************************
end module utilities
