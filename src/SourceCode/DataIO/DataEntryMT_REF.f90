!
!> Derived class to hold all the MT_REF information from a valid data file line
!
module DataEntryMT_REF
    !
    use DataEntryMT
    !
    implicit none
    !
    type, extends( DataEntryMT_t ), public :: DataEntryMT_REF_t
        !
        character(:), allocatable :: code_ref
        !
        real( kind=prec ) :: latitude_ref, longitude_ref, location_ref(3)
        !
    contains
        !
        final :: DataEntryMT_REF_dtor
        !
        procedure, public :: write => writeDataEntryMT_REF
        !
    end type DataEntryMT_REF_t
    !
    interface DataEntryMT_REF_t
        module procedure DataEntryMT_REF_ctor
    end interface DataEntryMT_REF_t
    !
contains
    !
    !> Parametrized constructor
    function DataEntryMT_REF_ctor( i_de, dtype,    &
        period, code, latitude, longitude, location, code_ref,    &
        latitude_ref, longitude_ref, location_ref, component, rvalue, imaginary, error ) result ( self )
        implicit none
        !
        type( DataEntryMT_REF_t ) :: self
        !
        integer, intent( in ) :: i_de
        character(:), allocatable, intent( in ) :: dtype, code, code_ref, component
        real( kind=prec ), intent( in ) :: period, latitude, longitude, &
                                           latitude_ref, longitude_ref, &
                                           location(3), location_ref(3)
        real( kind=prec ), intent( in ) :: rvalue, imaginary, error
        !
        !write( *, * ) "Constructor DataEntryMT_REF_t"
        !
        call self%baseInit
        !
        self%i_de = i_de
        self%dtype = trim( dtype )
        self%period = period
        self%code = trim( code )
        self%latitude = latitude
        self%longitude = longitude
        self%location = location
        self%code_ref = trim( code_ref )
        self%latitude_ref = latitude_ref
        self%longitude_ref = longitude_ref
        self%location_ref = location_ref
        self%component = trim( component )
        self%rvalue = rvalue
        self%imaginary = imaginary
        self%error = error
        !
    end function DataEntryMT_REF_ctor
    !
    !> Deconstructor routine:
    !>     Calls the base routine baseDealloc().
    !
    subroutine DataEntryMT_REF_dtor( self )
        implicit none
        !
        type( DataEntryMT_REF_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor DataEntryMT_REF_t:", self%id
        !
        call self%baseDealloc
        !
        if( allocated( self%code_ref ) ) deallocate( self%code_ref )
        !
    end subroutine DataEntryMT_REF_dtor
    !
    !> No subroutine briefing
    subroutine writeDataEntryMT_REF( self )
        class( DataEntryMT_REF_t ), intent( in ) :: self
        !
        write( *, * ) "Write DataEntryMT_REF_t: ", self%i_de, self%dtype, self%period, self%code,    &
        self%latitude, self%longitude, self%location, self%code_ref, self%latitude_ref, self%longitude_ref,    &
        self%location_ref, self%component, self%rvalue, self%imaginary, self%error
        !
    end subroutine writeDataEntryMT_REF
    !
end module DataEntryMT_REF
