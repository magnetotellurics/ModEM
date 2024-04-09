!
!> SUMMARY
!>
!> Some Fortran Character Utilities:
!> See http://gbenthien.net/strings/Strings.pdf  for more information
!> and additional subroutines.
!> 
!
module String
    !
    type, public :: String_t
        character(:), allocatable :: str
    end type
    !
contains
    !
    !> COMPACT
    !
    !> Converts multiple spaces and tabs to single spaces; 
    !> Deletes control characters;
    !> Removes initial spaces.
    !
    subroutine compact (str)
        implicit none 
        !
        character(len = *) :: str
        !
        character(len = 1) :: ch
        character(len = len_trim (str)) :: outstr
        integer :: lenstr, isp, ich, k, i
        !
        str = adjustl (str)
        lenstr = len_trim (str)
        outstr = ' '
        isp = 0
        k = 0
        !
        do i = 1, lenstr
            ch = str(i:i)
            ich = iachar (ch)
            !
            select case( ich )
                !
                case(9,32)      ! space or tab character
                    if(isp == 0) then
                        k = k + 1
                        outstr(k:k) = ' '
                    endif
                    isp = 1
                    !
                case(33:)        ! not a space, quote, or control character
                    k = k + 1
                    outstr(k:k) = ch
                    isp = 0
            end select
        enddo
        !
        str = adjustl (outstr)
        !
    end subroutine compact
    !
    !> PARSE
    !>
    !> Parses the string 'str' into arguments args(1), ..., args(nargs) based on
    !> the delimiters contained in the string 'delims'. Preceding a delimiter in
    !> 'str' by a backslash (\) makes this particular instance not a delimiter.
    !> The integer output variable nargs contains the number of arguments found.
    !
    subroutine Parse (str, delims, args, p_nargs)
        implicit none
        !
        character(len = *) :: str, delims
        character(len = len_trim (str)) :: strsav
        character(len = *) :: args(:)
        integer :: p_nargs
        !
        integer :: lenstr, isp, ich, k, i, na
        !
        strsav = str
        call compact (str)
        na = size (args)
        do i = 1, na
            args(i) = ' '
        enddo
        p_nargs = 0
        lenstr = len_trim (str)
        if(lenstr == 0) return
        k = 0
        !
        do
            if(len_trim (str) == 0) exit
            p_nargs = p_nargs + 1
            call Split (str, delims, args(p_nargs))
            call Removebksl (args(p_nargs))
        enddo
        !
        str = strsav
        !
    end subroutine Parse
    !
    !> SPLIT
    !>
    !> Routine finds the first instance of a character from 'delims' in the
    !> the string 'str'. The characters before the found delimiter are
    !> output in 'before'. The characters after the found delimiter are
    !> output in 'str'. The optional output character 'sep' contains the
    !> found delimiter. A delimiter in 'str' is treated like an ordinary
    !> character if it is preceded by a backslash (\). If the backslash
    !> character is desired in 'str', then precede it with another backslash.
    !
    subroutine Split( str, delims, before, sep )
        implicit none
        !
        character( len=* ) :: str, delims, before
        character, optional :: sep
        !
        logical :: pres
        character(1) :: ch
        character :: cha
        integer :: lenstr, isp, ich, k, i, p_nargs, na, ibsl, iposa, ipos
        !
        pres = present (sep)
        lenstr = len_trim (str)
        if(lenstr == 0) return    ! string str is empty
        k = 0
        ibsl = 0                  ! backslash initially inactive
        before = ' '
        do i = 1, lenstr
            ch = str(i:i)
            if(ibsl == 1) then    ! backslash active
                k = k + 1
                before(k:k) = ch
                ibsl = 0
                cycle
            endif
            if(ch == '\\') then   ! backslash with backslash inactive
                k = k + 1
                before(k:k) = ch
                ibsl = 1
                cycle
            endif
            ipos = index (delims, ch)
            if(ipos == 0) then    ! character is not a delimiter
                k = k + 1
                before(k:k) = ch
                cycle
            endif
            if(ch /= ' ') then          ! character is a delimiter that is not a space
                str = str(i+1:)
                if(pres) sep = ch
                exit
            endif
            cha = str(i+1:i+1)          ! character is a space delimiter
            iposa = index (delims, cha)
            if(iposa > 0) then          ! next character is a delimiter
                str = str(i+2:)
                if(pres) sep=cha
                exit
            else
                str = str(i+1:)
                if(pres) sep=ch
                exit
            endif
        enddo
        !
        if(i >= lenstr) str = ''
        str = adjustl (str)            ! remove initial spaces
        !
        return
        !
    end subroutine Split
    !
    !> REMOVEBKSl
    !>
    !> Removes backslash (\) characters. Double backslashes (\\) are 
    !> replaced by a single backslash.
    subroutine Removebksl (str)
        implicit none
        !
        character(len = *):: str
        !
        character(len = 1):: ch
        character(len = len_trim (str)) :: outstr
        integer :: lenstr, isp, ich, k, i, p_nargs, na, ibsl, iposa, ipos
        !
        str = adjustl (str)
        lenstr = len_trim (str)
        outstr = ' '
        k = 0
        ibsl = 0                        ! backslash initially inactive
        !
        do i = 1, lenstr
            ch = str(i:i)
            if(ibsl == 1) then            ! backslash active
                k = k + 1
                outstr(k:k) = ch
                ibsl = 0
                cycle
            endif
            if(ch == '\\') then        ! backslash with backslash inactive
                ibsl = 1
                cycle
            endif
            k = k + 1
            outstr(k:k) = ch            ! non-backslash with backslash inactive
        enddo
        !
        str = adjustl (outstr)
        !
    end subroutine Removebksl
    !
    !> IS_LETTER
    !>
    !> Returns .TRUE. if ch is a letter and .FALSE. otherwise.
    function Is_letter (ch) result (res)
        implicit none
        !
        character :: ch
        !
        logical :: res
        !
        select case(ch)
        case ('A':'Z','a':'z')
            res = .TRUE.
        case default
            res = .FALSE.
        end select
        !
        return
        !
    end function Is_letter
    !
    !> IS_DIGIT
    !>
    !> Returns .TRUE. if ch is a digit (0,1,...,9) and .FALSE. otherwise.
    function Is_digit (ch) result (res)
        implicit none
        !
        character :: ch
        !
        logical :: res
        !
        select case(ch)
        case ('0':'9')
            res = .TRUE.
        case default
            res = .FALSE.
        end select
        !
        return
        !
    end function Is_digit
    !>
    !> UPPER
    !> 
    !> Convert string to all uppercase.
    !> From: http://fortranwiki.org/fortran/show/String_Functions
    function Upper (s1) result (s2)
        implicit none
        !
        character(*) :: s1
        !
        character(len (s1)) :: s2
        character :: ch
        integer, parameter :: DUC = ICHAR ('A') - ICHAR ('a')
        integer :: i
        !
        do i = 1, len (s1)
            ch = s1(i:i)
            if(ch >= 'a'.AND.ch <= 'z') ch = CHAR (ICHAR (ch) + DUC)
            s2(i:i) = ch
        enddo
    end function Upper
    !
end module String
!