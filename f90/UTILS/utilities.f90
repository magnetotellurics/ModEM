module utilities

	use math_constants
	implicit none

  ! Variables required for storing the date and time in SECONDS. If used
  ! throughout the program, these make the routine profiling easier
  type  :: timer_t

    private
	real					:: rtime = 0.0 ! run time
	real					:: stime, etime ! start and end times

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

subroutine strcount(str,delims,nargs)

  ! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
  ! the delimiters contained in the string 'delims'. Preceding a delimiter in
  ! 'str' by a backslash (\) makes this particular instance not a delimiter.
  ! The integer output variable nargs contains the number of arguments found.

  ! Liu Zhongyin, 2019.08.27, copied from "parse", used to get str count
    implicit none
    character(len=*) :: str,delims
    character(len=len_trim(str)) :: strsav
    character(len=len_trim(str)) :: args
    integer                     :: lenstr,isp,ich,k,i,nargs
    strsav=str
    call compact(str)
    args = ' '
    nargs=0
    lenstr=len_trim(str)
    if(lenstr==0) return
    k=0
  
    do
       if(len_trim(str) == 0) exit
       nargs=nargs+1
       call split(str,delims,args)
       call removebksl(args)
    end do
    str=strsav
  
  end subroutine strcount

!**********************************************************************
! useful for grid operations & needed for I/O with HDF5 [Spencer Wilbur]
pure function cumsum(a) result (r)
real(kind=8), intent(in) :: a(:)
real(kind=8) :: r(size(a))
integer :: i
r(:) = [(sum(a(1:i)),i=1,size(a))]
end function cumsum

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
      call errStop('error, array and array index is not of same size in QSort!')
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
     call errStop('error, first index is larger than the last in QSort!')
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

! The following License applies to expand_string :
!
! Copyright (c) 2013-2019,  Los Alamos National Security, LLC (LANS) (Ocean: LA-CC-13-047;
! Land Ice: LA-CC-13-117) and the University Corporation for Atmospheric Research (UCAR).
!
! All rights reserved.
!
! LANS is the operator of the Los Alamos National Laboratory under Contract No.
! DE-AC52-06NA25396 with the U.S. Department of Energy.  UCAR manages the National
! Center for Atmospheric Research under Cooperative Agreement ATM-0753581 with the
! National Science Foundation.  The U.S. Government has rights to use, reproduce,
! and distribute this software.  NO WARRANTY, EXPRESS OR IMPLIED IS OFFERED BY
! LANS, UCAR OR THE GOVERNMENT AND NONE OF THEM ASSUME ANY LIABILITY FOR THE USE
! OF THIS SOFTWARE.  If software is modified to produce derivative works, such
! modified software should be clearly marked, so as not to confuse it with the
! version available from LANS and UCAR.
!
! Additionally, redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1) Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
!
! 2) Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation and/or
! other materials provided with the distribution.
!
! 3) None of the names of LANS, UCAR or the names of its contributors, if any, may
! be used to endorse or promote products derived from this software without
! specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
! ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
! ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!-----------------------------------------------------------------------
!  routine expand_string - (Renamed from mpas_log.F log_expand_string to expand_string)
!
!> \brief This is a utility routine that inserts formatted variables into a string.
!> \author Matt Hoffman
!> \date   02/20/2017
!> \details This routine inserts formatted variables into a string.
!>   The variables to be expanded are represented with a '$' symbol followed
!>   by one of these indicators:
!>   $i -> integer, formatted to be length of integer
!>   $l -> logical, fomatted as 'T' or 'F'
!>   $r -> real,  formatted as 9 digits of precision for SP mode, 17 for DP mode
!>                Floats are formatted using 'G' format which is smart about
!>                displaying as a decimal or scientific notation based on the value.
!>   The variable values to expand are supplied as optional arguments to this
!>   routine.  The substitution indicators are expanded as they are encountered.
!>   Thus, extra variable values will be ignored.  If the supplied variable values
!>   run out before the $ expansion indicators are all replaced, the remaining
!>   expansions will be filled with a fill value ('**').  The fill value is also
!>   used if the expansion indicator is of an unknown type, where the valid types
!>   are $i, $l, $r.
!>   If the user prefers more specific formatting, they have to do it external
!>   to this routine in a local string variable.  Similarly, character variables
!>   can be handled by the string concatenation command (//).
!>   This routine is based off of mpas_expand_string.
!-----------------------------------------------------------------------
subroutine expand_string(inString, outString, intArgs, logicArgs, realArgs)

  implicit none

  !-----------------------------------------------------------------
  ! input variables
  !-----------------------------------------------------------------
  character (len=*), intent(in) :: inString  !< Input: message to be expanded

  integer, dimension(:), intent(in), optional :: intArgs
     !< Input, Optional: array of integer variable values to be used in expansion
  logical, dimension(:), intent(in), optional :: logicArgs
     !< Input, Optional: array of logical variable values to be used in expansion
  real(kind=prec), dimension(:), intent(in), optional :: realArgs
     !< Input, Optional: array of real variable values to be used in expansion

  !-----------------------------------------------------------------
  ! input/output variables
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! output variables
  !-----------------------------------------------------------------
  character (len=512), intent(out) :: outString  !< Output: expanded version of input message after expansion

  !-----------------------------------------------------------------
  ! local variables
  !-----------------------------------------------------------------
  integer :: i, curLen
  integer :: nInts, nLogicals, nReals, nExps  !< the length of the variable arrays passed in
  integer :: iInt, iLogical, iReal !< Counter for the current index into each variable array
  character (len=64) :: realFormat  !< Format string to create to use for writing real variables to log file
  integer :: realPrecision !< precision of a real variable

  character (len=64) :: varPart
  character (len=64) :: errVarPart ! string to use if variable expansion fails
  logical :: charExpand

  ! Initialize the current index for each variable array to 1
  iInt = 1
  iLogical = 1
  iReal = 1

  ! For each variable array, get the size.  Size is 0 if not present.
  if (present(intArgs)) then
     nInts = size(intArgs)
  else
     nInts = 0
  endif

  if (present(logicArgs)) then
     nLogicals = size(logicArgs)
  else
     nLogicals = 0
  endif

  if (present(realArgs)) then
     nReals = size(realArgs)
  else
     nReals = 0
  endif

  ! Initialize strings
  write(outString,*) ''
  write(varPart,*) ''
  errVarPart = '**'  ! string to use if variable expansion fails

  !Initialize char info
  curLen = 0
  charExpand = .false.

  ! Loop over character positions in inString
  do i = 1, len_trim(inString)
     if (inString(i:i) == '$' .and. (.not. charExpand)) then
         charExpand = .true.
     else
         if (charExpand) then
            select case (inString(i:i))
               case ('i')
                  ! make the format large enough to include a large integer (up to 17 digits for 8-byte int)
                  ! it will be trimmed below
                  if (iInt <= nInts) then
                     write(varPart,'(i17)') intArgs(iInt)
                     iInt = iInt + 1
                  else
                     varPart = errVarPart
                  endif
               case ('l')
                  if (iLogical <= nLogicals) then
                     if (logicArgs(iLogical)) then
                        write(varPart ,'(a)') 'T'
                     else
                        write(varPart ,'(a)') 'F'
                     endif
                     iLogical = iLogical + 1
                  else
                     varPart = errVarPart
                  endif
               case ('r')
                  if (iReal <= nReals) then
                     realPrecision = precision(realArgs(iReal))
                     ! Note 7 additional characters may be needed beyond the precision
                     ! e.g.: -1234567.89012345   G13.6   ->   -0.123457E+07
                     write(realFormat, '(a, i2.2, a, i2.2, a)') '(G', realPrecision+7,'.', realPrecision, ')'
                     write(varPart, trim(realFormat)) realArgs(iReal)
                     iReal = iReal + 1
                  else
                     varPart = errVarPart
                  endif
               case ('$')
                     varPart = '$'
               case default
                  varPart = errVarPart
            end select
            outString = outString(1:curLen) // trim(adjustl(varPart))

            curLen = curLen + len_trim(adjustl(varPart))
            charExpand = .false.
         else
            outString(curLen+1:curLen+1) = inString(i:i)
            curLen = curLen+1
         end if
     end if
  end do

!--------------------------------------------------------------------
end subroutine expand_string
end module utilities
