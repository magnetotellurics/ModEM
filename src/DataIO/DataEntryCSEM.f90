!
!> Derived class to hold all the CSEM information from a valid data file line
!
module DataEntryCSEM
    !
    use DataEntry
    !
    type, extends( DataEntry_t ) :: DataEntryCSEM_t
        !
        character(:), allocatable :: dipole
        real( kind=prec ) :: moment, tx_azimuth, dip, tx_location(3)
        !
        contains
            !
            final :: DataEntryCSEM_dtor
            !
            procedure, public :: write => writeDataEntryCSEM
            !
    end type DataEntryCSEM_t
    !
    interface DataEntryCSEM_t
        module procedure DataEntryCSEM_ctor
    end interface DataEntryCSEM_t
    !
contains
    !
    !> Parametrized constructor
    function DataEntryCSEM_ctor( i_de, dtype, dipole, period, moment, tx_azimuth, &
        dip, tx_location, code, location, component, rvalue, imaginary, error, azimuth ) result( self )
        implicit none
        !
        type( DataEntryCSEM_t ) :: self
        integer, intent( in ) :: i_de
        character(:), allocatable, intent( in ) :: dtype, dipole, code, component
        real( kind=prec ), intent( in ) :: period, moment, tx_azimuth, dip, location(3), tx_location(3)
        real( kind=prec ), intent( in ) :: rvalue, imaginary, error
        real( kind=prec ), intent( in ), optional :: azimuth
        !
        !write( *, * ) "Constructor DataEntryCSEM_t", dipole
        !
        call self%init
        !
        self%i_de = i_de
        self%dtype = trim( dtype )
        self%dipole = trim( dipole )
        self%period = period
        self%moment = moment
        self%tx_azimuth = tx_azimuth
        self%dip = dip
        self%tx_location = tx_location
        self%code = trim( code )
        self%location = location
        self%component = trim( component )
        self%rvalue = rvalue
        self%imaginary = imaginary
        self%error = error
        !
        if( present( azimuth ) ) then
            !
            self%azimuth = azimuth
        else
            self%azimuth = R_ZERO
        endif
        !
    end function DataEntryCSEM_ctor
    !
    !> Deconstructor routine:
    !>     Calls the base routine dealloc().
    !
    subroutine DataEntryCSEM_dtor( self )
        implicit none
        !
        type( DataEntryCSEM_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor DataEntryCSEM_t:", self%id
        !
        call self%dealloc
        !
        if( allocated( self%dipole ) ) deallocate( self%dipole )
        !
    end subroutine DataEntryCSEM_dtor
    !
    !> No subroutine briefing
    subroutine writeDataEntryCSEM( self )
        class( DataEntryCSEM_t ), intent( in ) :: self
        !
        write( *, * ) "Write DataEntryCSEM_t: ", self%i_de, self%dtype, self%dipole, self%period,    &
        self%moment, self%tx_azimuth, self%dip, self%tx_location, self%code, self%location, self%component,    &
        self%rvalue, self%imaginary, self%error
        !
    end subroutine writeDataEntryCSEM
    !
end module DataEntryCSEM
