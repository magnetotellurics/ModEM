!
!> Abstract Base class to hold all the information from a valid data file line
!
module DataEntry
    !
    use Constants
    !
    type, abstract :: DataEntry_t
        !
        integer :: id
        character(:), allocatable :: type, code, component
        real( kind=prec ) :: period, location(3)
        real( kind=prec ) :: real, imaginary, error
        !
    contains
        !
        procedure( interface_write ), deferred, public :: write
        procedure( interface_get_copy_data_entry ), deferred, public :: getCopy
        !
        procedure, public :: isEqual => isEqualDe
        procedure, public :: isComplex => isComplexDe
        !
    end type DataEntry_t
    !
    abstract interface
        !
        !> No interface subroutine briefing
        subroutine interface_write( self )
            import :: DataEntry_t
            class( DataEntry_t ), intent( in ) :: self
        end subroutine interface_write
        !
        !> No interface function briefing
        function interface_get_copy_data_entry( self ) result ( copy )
            import :: DataEntry_t
            !
            class( DataEntry_t ), intent( in ) :: self
            class( DataEntry_t ), allocatable :: copy
        end function interface_get_copy_data_entry
        !
    end interface
    !
contains
    !
    !> No function briefing
    function isEqualDe( self, other ) result ( equal )
        class( DataEntry_t ), intent( in ) :: self, other
        logical :: equal
        !
        equal = .FALSE.
        !
        if( self%type .Eq. other%type .AND.                          &
            ABS( self%period - other%period ) < TOL6 .AND.           & 
            ABS( self%location(1) - other%location(1) ) < TOL6 .AND. &
            ABS( self%location(2) - other%location(2) ) < TOL6 .AND. &
            ABS( self%location(3) - other%location(3) ) < TOL6 .AND. &
            self%code .Eq. other%code .AND.                          &
            self%component .Eq. other%component .AND.                &
            ABS( self%real - other%real ) < TOL6 .AND.               &
            ABS( self%imaginary - other%imaginary ) < TOL6 .AND.     &
            ABS( self%error - other%error ) < TOL6 ) then
                equal = .TRUE.
        endif
        !
    end function isEqualDe
    !
    !> No function briefing
    function isComplexDe( self ) result ( complex )
        class( DataEntry_t ), intent( in ) :: self
        logical :: complex
        !
        complex = .TRUE.
        !
        if( ( index( self%type, "Off_Diagonal_Rho_Phase" ) /= 0 ) .OR.   &
            ( index( self%type, "Phase_Tensor" ) /= 0 ) ) then
        !
            complex = .FALSE.
        !
        endif
        !
    end function isComplexDe
    !
end module DataEntry
