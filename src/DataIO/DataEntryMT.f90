!
!> Derived class to hold all the MT information from a valid data file line
!
module DataEntryMT
    !
    use DataEntry
    !
    type, extends( DataEntry_t ) :: DataEntryMT_t
        !
        real( kind=prec ) :: latitude, longitude
        !
    contains
        !
        final :: DataEntryMT_dtor
        !
        procedure, public :: write => writeDataEntryMT
        !
    end type DataEntryMT_t
    !
    interface DataEntryMT_t
        module procedure DataEntryMT_ctor
    end interface DataEntryMT_t
    !
contains
    !
    !> Parametrized constructor
    function DataEntryMT_ctor( i_de, dtype, period, code, latitude, longitude, location, component, rvalue, imaginary, error ) result ( self )
        implicit none
        !
        type( DataEntryMT_t ) :: self
        !
        integer, intent( in ) :: i_de
        character(:), allocatable, intent( in ) :: dtype, code, component
        real( kind=prec ), intent( in ) :: period, latitude, longitude, location(3)
        real( kind=prec ), intent( in ) :: rvalue, imaginary, error
        !
        !write( *, * ) "Constructor DataEntryMT_t"
        !
        call self%init
        !
        self%i_de = i_de
        self%dtype = dtype
        self%period = period
        self%code = code
        self%latitude = latitude
        self%longitude = longitude
        self%location = location
        self%component = component
        self%rvalue = rvalue
        self%imaginary = imaginary
        self%error = error
        !
    end function DataEntryMT_ctor
    !
    !> Deconstructor routine:
    !>     Calls the base routine dealloc().
    !
    subroutine DataEntryMT_dtor( self )
        implicit none
        !
        type( DataEntryMT_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor DataEntryMT_t:", self%id
        !
        call self%dealloc
        !
    end subroutine DataEntryMT_dtor
    !
    !> No subroutine briefing
    subroutine writeDataEntryMT( self )
        class( DataEntryMT_t ), intent( in ) :: self
        !
        write( *, * ) "Write DataEntryMT_t: ", self%i_de, self%dtype, self%period, self%code,    &
        self%latitude, self%longitude, self%location, self%component, self%rvalue, self%imaginary, self%error
        !
    end subroutine writeDataEntryMT
    !
end module DataEntryMT
