!>*************
!>
!> Class to read a control file
!> And set Inversion parameters
!>
!>*************
!>
module InversionControlFile
    !
    use Constants
    use String
    !
    character(:), allocatable :: inversion_type
    !
    type :: InversionControlFile_t
        !
        !> Inversion parameters
        character(:), allocatable :: inversion_type
        character(:), allocatable :: max_inv_iters, max_grad_iters
        character(:), allocatable :: tolerance_error, tolerance_rms
        character(:), allocatable :: lambda
        !
        contains
            !
            final :: InversionControlFile_dtor
            !
    end type InversionControlFile_t
    !
    interface InversionControlFile_t
        module procedure InversionControlFile_ctor
    end interface InversionControlFile_t
!
contains
    !
    !> Procedure InversionControlFile_ctor
    !> Read line by line of the data file, create Data Entry objects(MT, MT_REF or CSEM)
    function InversionControlFile_ctor( funit, fname ) result( self )
        implicit none
        !
        integer, intent( in ) :: funit
        character(:), allocatable, intent( in ) :: fname
        !
        type( InversionControlFile_t ) :: self
        !
        character(1000) :: full_line_text
        character(len=200), dimension(20) :: args
        character(:), allocatable :: line_text
        integer :: line_counter, io_stat, p_nargs
        !
        !write( *,* ) "Constructor InversionControlFile_t"
        !
        call Compact( fname )
        !
        open( unit = funit, file = fname, iostat = io_stat, status = "old" )
        !
        if( io_stat == 0 ) then
            !
            do
                read( funit, "(a)", END = 10 ) full_line_text
                line_text = adjustl( full_line_text )
                line_text = trim( line_text )
                !
                call Parse( line_text, ":", args, p_nargs )
                !
                if( index( line_text, "#" ) == 0 .AND. index( line_text, ">" ) == 0 ) then
                    !
                    if( index( line_text, "inversion_type" ) > 0 ) then
                        self%inversion_type = trim( args(2) )
                    else if( index( line_text, "max_inv_iters" ) > 0 ) then
                        self%max_inv_iters = trim( args(2) )
                    else if( index( line_text, "max_grad_iters" ) > 0 ) then
                        self%max_grad_iters = trim( args(2) )
                    else if( index( line_text, "tolerance_error" ) > 0 ) then
                        self%tolerance_error = trim( args(2) )
                    else if( index( line_text, "tolerance_rms" ) > 0 ) then
                        self%tolerance_rms = trim( args(2) )
                    else if( index( line_text, "lambda" ) > 0 ) then
                        self%lambda = trim( args(2) )
                    else
                        write( *, * ) "Error: Unsupported Inversion parameter: ["//trim(line_text)//"]"
                        stop 
                    endif
                    !
                endif
                !
            enddo
            !
10          close( unit = funit )
            !
            ! Inversion type
            if( allocated( self%inversion_type ) ) then
                !
                select case( self%inversion_type )
                    !
                    case( "DCG" )
                        inversion_type = DCG
                    case( "NLCG" )
                        inversion_type = NLCG
                    case default
                        inversion_type = ""
                    stop "Error: Wrong inversion_type control, use [DCG|NLCG]"
                    !
                end select
                !
            endif
            !
        else
            write( *, * ) "Error opening [", fname, "] in InversionControlFile_ctor"
            stop
        endif
        !
    end function InversionControlFile_ctor
    !
    !> Deconstructor routine:
    !>     Deallocates inherent properties of this class.
    subroutine InversionControlFile_dtor( self )
        implicit none
        !
        type( InversionControlFile_t ), intent( inout ) :: self
        !
        !write( *,* ) "Destructor InversionControlFile_t"
        !
        if( allocated( self%inversion_type ) ) deallocate( self%inversion_type )
        if( allocated( self%max_inv_iters ) ) deallocate( self%max_inv_iters )
        if( allocated( self%max_grad_iters ) ) deallocate( self%max_grad_iters )
        if( allocated( self%tolerance_error ) ) deallocate( self%tolerance_error )
        if( allocated( self%tolerance_rms ) ) deallocate( self%tolerance_rms )
        if( allocated( self%lambda ) ) deallocate( self%lambda )
        !
    end subroutine InversionControlFile_dtor
    !
end module InversionControlFile
