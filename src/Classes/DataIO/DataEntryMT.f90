!*************
!
! Derived class to hold all the MT information from a valid data file line
!
! Last modified at 05/2021 by Paulo Werdt
!
!*************
!
module DataEntryMT
    !
    use DataEntry
    !
    type, extends( DataEntry_t ) :: DataEntryMT_t
        !
        real( kind=prec ) :: latitude, longitude

    contains
        !
        procedure, public :: write => writeDataEntryMT
        procedure, public :: getCopy => getCopyDataEntryMT
        !
    end type DataEntryMT_t
    !
    interface DataEntryMT_t
        module procedure DataEntryMT_ctor
    end interface DataEntryMT_t
    !
contains
    !
    ! Parametrized constructor
    function DataEntryMT_ctor( id, type, period, code, latitude, longitude, location, component, real, imaginary, error ) result ( self )
        implicit none
        !
        type( DataEntryMT_t ) :: self
        !
        integer, intent( in )                   :: id
        character(:), allocatable, intent( in ) :: type, code, component
        real( kind=prec ), intent( in )         :: period, latitude, longitude, location(3)
        real( kind=prec ), intent( in )         :: real, imaginary, error
        !
        !write(*,*) "Constructor DataEntryMT_t"
        !
        self%id = id
        self%type = type
        self%period = period
        self%code = code
        self%latitude = latitude
        self%longitude = longitude
        self%location = location
        self%component = component
        self%real = real
        self%imaginary = imaginary
        self%error = error
        !
    end function DataEntryMT_ctor
    !
    function getCopyDataEntryMT( self ) result ( copy )
        implicit none
        !
        class( DataEntryMT_t ), intent( in ) :: self
        class( DataEntry_t ), allocatable     :: copy
        !
        allocate( copy, source = DataEntryMT_t( self%id, self%type, self%period, self%code, &
                  self%latitude, self%longitude, self%location, &
                  self%component, self%real, self%imaginary, self%error ) )
        !
    end function getCopyDataEntryMT
    !
    subroutine writeDataEntryMT( self )
        class( DataEntryMT_t ), intent( in )    :: self
        !
        write(*,*) "Write DataEntryMT_t: ", self%id, self%type, self%period, self%code,    &
        self%latitude, self%longitude, self%location, self%component, self%real, self%imaginary, self%error
        !
    end subroutine writeDataEntryMT
    !
end module DataEntryMT
