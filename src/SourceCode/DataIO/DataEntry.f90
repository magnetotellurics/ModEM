!
!> Abstract Base class to hold all the information from a valid data file line
!
module DataEntry
    !
    use Constants
    !
    type, abstract :: DataEntry_t
        !
        integer :: i_de
        character(:), allocatable :: dtype, code, component
        real( kind=prec ) :: period, location(3)
        real( kind=prec ) :: rvalue, imaginary, error, azimuth
        !
    contains
        !
        procedure( interface_write ), deferred, public :: write
        !
        procedure, public :: baseInit => initializeDataEntry
        !
        procedure, public :: baseDealloc => deallocateDataEntry
        !
        procedure, public :: isEqual => isEqualDe
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
    end interface
    !
contains
    !
    !> Initialize transmitter base variables (Avoid initialization on declaration).
    !
    subroutine initializeDataEntry( self )
        implicit none
        !
        class( DataEntry_t ), intent( inout ) :: self
        !
        self%i_de = 0
        self%period = R_ZERO
        self%location = R_ZERO
        self%rvalue = R_ZERO
        self%imaginary = R_ZERO
        self%error = R_ZERO
        self%azimuth = R_ZERO
        !
    end subroutine initializeDataEntry
    !
    !> Free the memory used by all allocatable variables belonging to this transmitter.
    !> Called before anything in the destructor of derived classes.
    !
    subroutine deallocateDataEntry( self )
        implicit none
        !
        class( DataEntry_t ), intent( inout ) :: self
        !
        if( allocated( self%dtype ) ) deallocate( self%dtype )
        !
        if( allocated( self%code ) ) deallocate( self%code )
        !
        if( allocated( self%component ) ) deallocate( self%component )
        !
    end subroutine deallocateDataEntry
    !
    !> No subroutine briefing
    !
    function isEqualDe( self, other ) result ( equal )
        class( DataEntry_t ), intent( in ) :: self, other
        logical :: equal
        !
        equal = .FALSE.
        !
        if( self%dtype .EQ. other%dtype .AND.                        &
            ABS( self%period - other%period ) < TOL6 .AND.           & 
            ABS( self%location(1) - other%location(1) ) < TOL6 .AND. &
            ABS( self%location(2) - other%location(2) ) < TOL6 .AND. &
            ABS( self%location(3) - other%location(3) ) < TOL6 .AND. &
            self%code .EQ. other%code .AND.                          &
            self%component .EQ. other%component .AND.                &
            ABS( self%rvalue - other%rvalue ) < TOL6 .AND.           &
            ABS( self%imaginary - other%imaginary ) < TOL6 .AND.     &
            ABS( self%error - other%error ) < TOL6 ) then
                equal = .TRUE.
        endif
        !
    end function isEqualDe
    !
end module DataEntry
