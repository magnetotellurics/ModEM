module Utilities
    !
    use Constants
    use FileUnits
    !
    integer :: warning_counter = 0
    !
    character( len=200 ) :: str_msg
    !
    character( len=15 ) :: run_tag
    !
    character( len=50 ) :: outdir_name
    !
    public :: clean, minNode, maxNode
    public :: errStop, warning, writeWarning
    public :: QSort
    !
    !> Variables required for storing the date and time in SECONDS. If used
    !> throughout the program, these make the routine profiling easier
    !
    type :: timer_t
        !
        private
        real :: rtime = 0.0 ! run time
        real :: stime, etime ! start and end times
        !
    end type timer_t
    !
contains
    !
    !> No Subroutine Briefing
    !
    subroutine errStop( msg )
        implicit none
        !
        character(*), intent( in ) :: msg
        !
        !> System Bip!
        write( *, * ) char(7)
        !
        write( *, * ) achar(27)//"[31m# Error:"//achar(27)//"[0m "//trim( msg )
        !
#ifdef MPI
        !
        call MPI_Abort( main_comm, error_code, ierr )
        !
        call MPI_Finalize( ierr )
        !
#endif
        !
        stop
        !
    end subroutine errStop
    !
    !> No Subroutine Briefing
    !
    subroutine warning( msg )
        implicit none
        !
        character(*), intent( in ) :: msg
        !
        !> System Bip!
        write( *, * ) char(7)
        !
        write( *, * ) achar(27)//"[91m# Warning:"//achar(27)//"[0m "//trim( msg )
        !
        warning_counter = warning_counter + 1
        !
        call writeWarning( msg )
        !
    end subroutine warning
    !
    subroutine writeWarning( msg )
        implicit none
        !
        character(*), intent( in ) :: msg
        !
        integer :: ios
        !
        open( unit = ioWarning, &
        file = "warnings_"//run_tag//".log", &
        status = "unknown", position = "append", iostat = ios )
        !
        if( ios == 0 ) then
            !
            write( ioWarning, * ) warning_counter, ":", trim( msg )
            !
            close( ioWarning )
            !
        else
            call errStop( "writeWarning > cant open [warnings_"//run_tag//".log]" )
        endif
        !
    end subroutine writeWarning
    !
    !> Return the number of lines at warning file
    !
    subroutine printWarningBrief()
        implicit none
        !
        integer :: ios, counter
        logical :: exist_warnings
        character(100) :: line_text
        !
        counter = 0
        exist_warnings = .FALSE.
        !
        inquire( file = "warnings_"//run_tag//".log", exist = exist_warnings )
        !
        if( exist_warnings ) then
            !
            open( unit = ioWarning, file = "warnings_"//run_tag//".log", iostat = ios )
            !
            if( ios == 0 ) then
                !
                do
                    !
                    read( ioWarning, "(A)", iostat = ios ) line_text
                    !
                    if ( ios /= 0 ) exit
                    !
                    counter = counter + 1
                    !
                end do
                !
                close( ioWarning )
                !
                write( str_msg, "( I4, A49 )" ) counter, " entries listed in [warnings_"//run_tag//".log]"
                !
                write( *, * ) ""
                !
                write( *, * ) achar(27)//"[91m# Warning:"//achar(27)//"[0m "//str_msg
                !
            else
                call errStop( "printWarningBrief > cant open [warnings_"//run_tag//".log]." )
            endif
            !
        else
            !
            write( *, * )
            !
        endif
        !
    end subroutine printWarningBrief
    !
    !> Timer utilities: set timer
    !
    subroutine reset_time( timer )
        implicit none
        !
        type(timer_t), intent( inout ) :: timer
        ! utility variable
        integer, dimension(8) :: tarray
        !
        ! Restart the(portable) clock
        call date_and_time(values=tarray)
        timer%stime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)

    end subroutine reset_time
    !
    !> Timer utilities: compute elapsed run time in seconds
    !
    function elapsed_time( timer ) result( rtime )
        implicit none
        !
        type( timer_t ), intent( inout ) :: timer
        real :: rtime ! run time
        ! utility variable
        integer, dimension(8) :: tarray
        !
        call date_and_time(values=tarray)
        timer%etime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)
        rtime = timer%etime - timer%stime
        !
        ! Update total run time
        timer%rtime = timer%rtime + rtime
        !
    end function elapsed_time
    !
    !> Timer utilities: sum up the elapsed run times in seconds
    !
    function saved_time( timer ) result( rtime )
        implicit none
        !
        type( timer_t ), intent( inout ) :: timer
        real :: rtime
        !
        rtime = timer%rtime
        !
    end function saved_time
    !
    !> Timer utilities: clear saved timing information
    !
    subroutine clear_time( timer )
        implicit none
        !
        type( timer_t ), intent( inout ) :: timer
        ! utility variable
        integer, dimension(8) :: tarray
        !
        ! Restart the(portable) clock
        call date_and_time(values=tarray)
        timer%stime = tarray(5)*3600 + tarray(6)*60 + tarray(7) + 0.001*tarray(8)
        !
        ! Clear saved run time
        timer%rtime = 0.0
        !
    end subroutine clear_time
    !
    !> This is a utility routine that provides an expression used to battle
    !> against machine error problems. Both input and output are values in km.
    !> The function rounds the value to the nearest meter. This is useful to
    !> ensure that the grid read from a file does not depend on system precision.
    !> A.K.
    !
    function nearest_meter( x ) result( clean )
        implicit none
        !
        real( kind=prec ), intent( in ) :: x
        real( kind=prec ) :: clean
        !
        clean = dnint( x * KM2M ) / KM2M
        !
    end function nearest_meter
    !
    !> For an increasing or decreasing array, 
    !> find the minimum and maximum indices
    !> such that xmin <= array(i) < xmax
    !> author: A. Kelbert
    !
    subroutine find_index( array, xmin, xmax, imin, imax )
        implicit none
        !
        real( kind=prec ), dimension(:), intent( in ) :: array
        real( kind=prec ), intent( in ) :: xmin, xmax
        integer, intent( out ) :: imin, imax
        !
        logical :: incr
        integer :: i, n
        !
        ! quick check to see what kind of array this is
        n = size(array)
        !
        if(array(1) <= array(n)) then
            incr = .TRUE.
        else
            incr = .FALSE.
        endif

        if(incr) then
            !
            ! for an increasing array...
            imin = 0
            !
            do i = 1, n
                if(clean(array(i)) .ge. clean(xmin)) then
                    imin = i
                    exit
                endif
            enddo
            !
            imax = 0
            !
            do i = imin, n
                if(clean(array(i)) .LT. clean(xmax)) then
                    imax = i
                endif
            enddo
            !
        else
            !
            ! for a decreasing array...
            imax = 0
            do i = n, 1, -1
                if(clean(array(i)) .ge. clean(xmin)) then
                    imax = i
                    exit
                endif
            enddo
            !
            imin = 0
            do i = imax, 1, -1
                if(clean(array(i)) .LT. clean(xmax)) then
                    imin = i
                endif
            enddo
            !
        endif
        !
    end subroutine find_index
    !
    !> Replicates the corresponding function in Matlab: for an integer array, 
    !> outputs true if our integer is in the array, otherwise false.
    !> author: A. Kelbert
    !
    logical function ismember( n, Nvec )
        implicit none
        !
        integer :: n
        integer, dimension(:) :: Nvec
        !
        integer :: i
        !
        ismember = .FALSE.
        do i = 1, size(Nvec)
            if(Nvec(i) == n) then
                ismember = .TRUE.
                return
            endif
        enddo
        !
    end function
    !
    !> Return the position of str2 in str1.  Ignores case.
    !> Return 0 if str2 not found in str1
    !
    integer function findstr( str1, str2 )
        implicit none
        !
        character*(*) str1, str2
        integer i, j, capdif
        logical same
        !
        capdif= ichar("a")-ichar("A")
        !
        do 20 i= 1, len(str1)-len(str2)+1
        do 10 j=1, len(str2)
        !
        same= str1(i+j-1:i+j-1) .EQ. str2(j:j) .OR.  &
        "A".le.str2(j:j) .AND. str2(j:j).le."Z" .AND.  &
        ichar(str1(i+j-1:i+j-1)) .EQ. ichar(str2(j:j))+capdif .OR.  &
        "a".le.str2(j:j) .AND. str2(j:j).le."z" .AND.  &
        ichar(str1(i+j-1:i+j-1)) .EQ. ichar(str2(j:j)) - capdif
        !
        if( .NOT.same) go to 20
        10       continue
        findstr=i
        return
        20    continue
        !
        findstr=0
        !
    end function findstr
    !
    !> Return the index of the first non-blank character in the iwrd"th
    !> non-blank word(word are seperated by spaces, tabs or commas).
    !> Return len if iwrd"th word is not found. integer i, nword
    !
    integer function begwrd(string, iwrd)
        implicit none
        !
        integer iwrd
        character*(*) string
        !
        logical wasblk
        intrinsic len
        integer  i, nword
        !
        wasblk = .TRUE.
        !
        nword= 0
        do i=1, len(string)
            if( string(i:i).EQ." " .OR. string(i:i) .EQ. ", " .OR. string(i:i).EQ."  "    )then
                !
                !           /* current character is blank */
                wasblk=.TRUE.
            else
                if(wasblk) then
                    nword= nword + 1
                endif
                wasblk= .FALSE.
                if(nword.EQ.iwrd)then
                    begwrd= i
                    return
                endif
            endif
        enddo
        !
        begwrd = len(string)
        !
    end function begwrd
    !
    !> Return the index of the last non-blank character in the iwrd"th
    !> non-blank word(word are seperated by spaces, tabs or commas).
    !> Return len if iwrd"th word is not found.
    !
    integer function endwrd(string, iwrd)
        implicit none
        !
        integer iwrd
        character*(*) string
        integer i, nword
        logical wasblk
        intrinsic len
        !
        wasblk=.TRUE.
        nword= 0
        do 100 i=1, len(string)
        if( string(i:i).EQ." " .OR.  &
        string(i:i).EQ.", " .OR.  &
        string(i:i).EQ."  "    )then
        !
        !          /* current character is blank */
        wasblk=.TRUE.
        if(nword.EQ.iwrd) return
        !
        else
        if(wasblk) nword= nword + 1
        wasblk= .FALSE.
        if(nword.EQ.iwrd) endwrd= i
        endif
        100   continue
        !
        endwrd= len(string)
        !
    end function endwrd
    !
    !>  No subroutine briefing
    !
    subroutine lenb( string, length )
        implicit none
        !
        character*(*) string
        integer nstr, istr, length
        !
        nstr = len(string)
        do istr=nstr, 1, -1
            if(string(istr:istr) .NE. " ") then
                length = istr
                return
            endif
        enddo
        length = 0
        !
    end subroutine lenb
    !
    !> Naser Meqbel included this function: apparently, it is not supported by
    !> all compilers as an intrinsic
    !
    logical function isnan( a )
        implicit none
        !
        real( kind=prec ), intent( in ) :: a
        !
        if(a .NE. a) then
        isnan = .TRUE.
        else
        isnan = .FALSE.
        endif
        !
    end function isnan
    !
    !> Routine finds the first instance of a character from "delims" in the
    !> the string "str". The characters before the found delimiter are
    !> output in "before". The characters after the found delimiter are
    !> output in "str". The optional output character "sep" contains the
    !> found delimiter. A delimiter in "str" is treated like an ordinary
    !> character if it is preceded by a backslash(\). If the backslash
    !> character is desired in "str", then precede it with another backslash.
    !
    subroutine split( str, delims, before, sep )
        implicit none
        !
        character(len=*) :: str, delims, before
        character, optional :: sep
        logical :: pres
        character(1) :: ch
        character :: cha
        integer :: lenstr, isp, ich, k, i, nargs, na, ibsl, iposa, ipos
        pres=present(sep)
        lenstr=len_trim(str)
        if(lenstr == 0) return      ! string str is empty
        k=0
        ibsl=0                      ! backslash initially inactive
        before=" "
        do i=1, lenstr
        ch=str(i:i)
        if(ibsl == 1) then          ! backslash active
        k=k+1
        before(k:k)=ch
        ibsl=0
        cycle
        endif
        if(ch == "\\") then         ! backslash with backslash inactive
        k=k+1
        before(k:k)=ch
        ibsl=1
        cycle
        endif
        ipos=index(delims, ch)
        if(ipos == 0) then          ! character is not a delimiter
        k=k+1
        before(k:k)=ch
        cycle
        endif
        if(ch /= " ") then          ! character is a delimiter that is not a space
        str=str(i+1:)
        if(pres) sep=ch
        exit
        endif
        cha=str(i+1:i+1)            ! character is a space delimiter
        iposa=index(delims, cha)
        if(iposa > 0) then          ! next character is a delimiter
        str=str(i+2:)
        if(pres) sep=cha
        exit
        else
        str=str(i+1:)
        if(pres) sep=ch
        exit
        endif
        enddo
        if(i >= lenstr) str=""
        str=adjustl(str)            ! remove initial spaces
        !
    end subroutine split
    !
    !> Remove backslash(\) characters. Double backslashes(\\) are replaced
    !> by a single backslash.
    !
    subroutine removebksl( str )
        implicit none
        !
        character(len=*):: str
        character(len=1):: ch
        character(len=len_trim(str))::outstr
        integer :: lenstr, isp, ich, k, i, nargs, na, ibsl, iposa, ipos
        str=adjustl(str)
        lenstr=len_trim(str)
        outstr=" "
        k=0
        ibsl=0                    ! backslash initially inactive

        do i=1, lenstr
        ch=str(i:i)
        if(ibsl == 1) then        ! backslash active
        k=k+1
        outstr(k:k)=ch
        ibsl=0
        cycle
        endif
        if(ch == "\\") then        ! backslash with backslash inactive
        ibsl=1
        cycle
        endif
        k=k+1
        outstr(k:k)=ch            ! non-backslash with backslash inactive
        enddo

        str=adjustl(outstr)
        !
    end subroutine removebksl
    !
    !> Return .TRUE. if ch is a letter and .FALSE. otherwise
    !
    function is_letter( ch ) result( res )
        implicit none
        !
        character :: ch
        logical :: res

        select case(ch)
        case("A":"Z", "a":"z")
        res=.TRUE.
        case default
        res=.FALSE.
        end select
        !
    end function is_letter
    !
    !> Return .TRUE. if ch is a digit(0, 1, ..., 9) and .FALSE. otherwise
    !
    function is_digit( ch ) result( res )
        implicit none
        !
        character :: ch
        logical :: res
        !
        select case(ch)
            case("0":"9")
                res=.TRUE.
            case default
                res=.FALSE.
        end select
        !
    end function is_digit
    !
    !> This is a utility routine that provides an expression used to battle
    !> against machine error problems. It returns the same real or real(8)
    !> as the input, but without the extra digits at the end that are often
    !> a cause of wrong comparisons in the if statements. ALWAYS use clean(x)
    !> instead of x in an inequality!!!
    !> R_LARGE is defined in the module math_constants
    !> A.K.
    !
    function clean( x )
        implicit none
        !
        real( kind=prec ), intent( in ) :: x
        real( kind=prec ) :: clean
        !
        clean = dnint( x * R_LARGE) / R_LARGE
        !
    end function clean
    !
    !> This is a utility routine, used by several data functional
    !> set up routines, and for other interpolation functions
    !> Return index ix such that  xNode(ix) <= x < xNode(ix+1)
    !> If x is out of range:
    !> x < xNode(1) Return 0; if x> xNode(nx) Return nx
    !> Assumes xNode is strictly increasing; does not check this
    !> NOTE: as presently coded, when xNode is called with center
    !>(face) node positions, this routine will return zero for
    !> the coordinates in the outer half cell nearest the boundary
    !> If evaluation over the complete model domain is to be allowed
    !> a more general interpolation rule will be required.
    !> A.K.: modified to allow input of any size, nx = size(xNode).
    !
    function minNode( x, xNode ) result( ix )
        implicit none
        !
        real( kind=prec ), intent( in ) :: x
        real( kind=prec ), dimension(:), intent( in ) :: xNode
        !
        integer :: ix, i
        !
        do i = 1, size( xNode )
            if( clean( xNode(i) ) .GT. clean(x) ) then
                ix = i-1
                exit
            endif
        enddo
        !
    end function minNode
    !
    !>    This is a utility routine, used by several data functional
    !>    set up routines, and for other interpolation functions
    !>    Returns index ix such that    xNode(ix) <= x < xNode(ix+1)
    !>    If x is out of range:
    !>    x > xNode(1) returns 0; if x< xNode(nx) returns nx
    !>    Assumes xNode is strictly decreasing; does not check this
    !>    NOTE: as presently coded, when xNode is called with center
    !>    (face) node positions, this routine will return zero for
    !>    the coordinates in the outer half cell nearest the boundary
    !>    If evaluation over the complete model domain is to be allowed
    !>    a more general interpolation rule will be required.
    !>    A.K.: modified to allow input of any size, nx = size(xNode).
    !
    function maxNode( x, xNode ) result( ix )
        implicit none
        !
        real( kind=prec ), intent( in ) :: x
        real( kind=prec ), dimension(:), intent( in ) :: xNode
        !
        integer :: ix, i
        !
        do i = 1, size(xNode)
           if( clean( xNode(i) ) .LT. clean(x) ) then
                ix = i-1
                exit
           endif
        enddo
        !
    end function maxNode
    !
    !> A simple recursive quick sort routine using a middle pivot
    !> the average time complexity is O(nlog(n)), with the worst case of 
    !> O(n.^2)
    !
    recursive subroutine QSort( a, ia, i0, i1 )
        implicit none
        !
        real( kind=prec ), intent( inout ), dimension(:) :: a
        real( kind=prec ) :: pivot, t, random
        integer, intent( inout ), dimension(:) :: ia
        integer, intent( in ), optional :: i0, i1
        !
        integer :: first, last, i, j, it, nA
        !
        if( size(a) .NE. size(ia) ) then
            stop "Error: QSort > array and array index is not of same size in QSort!"
        endif
        !
        if( .NOT.present(i0) ) then
            first = 1
            last = size(a)
        elseif(.NOT.present(i1)) then
            first = i0
            last = size(a)
        else 
            first = i0
            last = i1
        endif
        !
        if( first .GT. last ) then 
            stop "Error: QSort > first index is larger than the last in QSort!"
        elseif( first .EQ. last ) then !only one element, no need to sort now
            !     write(6, *) "no need to sort, only one element left"
            return
        endif
        !Na = j-i+1
        !call random_number(random)
        !write(6, *) Na, int(random*real(Na-1))
        !pivot = a(int(random*real(Na-1))+1)! using a random pivot
        pivot = a((first+last)/2) ! using midpoint
        !  write(6, *)  "taking a pivot value of ", int(pivot)
        i = first
        j = last
        do
            !
            do while(a(i).LT.pivot) !left half
                i=i+1
            enddo
            !
            do while(pivot.LT.a(j)) !right half
                j=j-1
            enddo
            !
            if( i .LT. j ) then 
                ! swap array and index
                t = a(i)
                a(i) = a(j)
                a(j) = t
                it = ia(i)
                ia(i) = ia(j)
                ia(j) = it
                i=i+1
                j=j-1
                !        write(6, *) "swapping ai and aj: "
                !        write(6, *) int(a)
            else
                exit
            endif
            !
        enddo
        !
        if( first .LT. i-1 ) then 
            !      write(6, *)  "taking care of the left part"
            call QSort( a, ia, first, i-1 )
        endif
        !
        if( j+1 .LT. last ) then 
            !      write(6, *)  "taking care of the right part"
            call QSort( a, ia, j+1, last )
        endif
        !
    end subroutine QSort
    !
end module Utilities
