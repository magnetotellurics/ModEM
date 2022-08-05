!*************
!
! Derived class to hold all the MT_REF information from a valid data file line
!
! Last modified at 05/2021 by Paulo Werdt
!
!*************
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
        real( kind=prec )            :: latitude_ref, longitude_ref, location_ref(3)
        !
    contains
        !
        procedure, public    :: write => writeDataEntryMT_REF
        procedure, public    :: getCopy => getCopyDataEntryMT_REF
        !
    end type DataEntryMT_REF_t
    !
    interface DataEntryMT_REF_t
        module procedure DataEntryMT_REF_ctor
    end interface DataEntryMT_REF_t
    !
contains
    !
    ! Parametrized constructor
    function DataEntryMT_REF_ctor( id, type,    &
            period, code, latitude, longitude, location, code_ref,    &
            latitude_ref, longitude_ref, location_ref, component, real, imaginary, error ) result ( self )
        implicit none
        !
        type( DataEntryMT_REF_t ) :: self
        !
        integer, intent( in )                         :: id
        character(:), allocatable, intent( in ) :: type, code, code_ref, component
        real( kind=prec ), intent( in )            :: period, latitude, longitude, &
                                                                latitude_ref, longitude_ref, &
                                                                 location(3), location_ref(3)
        real( kind=prec ), intent( in )            :: real, imaginary, error
        !
        !write(*,*) "Constructor DataEntryMT_REF_t"
        !
        self%id = id
        self%type = type
        self%period = period
        self%code = code
        self%latitude = latitude
        self%longitude = longitude
        self%location = location
        self%code_ref = code_ref
        self%latitude_ref = latitude_ref
        self%longitude_ref = longitude_ref
        self%location_ref = location_ref
        self%component = component
        self%real = real
        self%imaginary = imaginary
        self%error = error
        !
    end function DataEntryMT_REF_ctor
    !
    function getCopyDataEntryMT_REF( self ) result ( copy )
        implicit none
        !
        class( DataEntryMT_REF_t ), intent( in ) :: self
      class( DataEntry_t ), allocatable          :: copy
        !
        allocate( copy, source = DataEntryMT_REF_t( self%id, self%type,    &
                     self%period, self%code, self%latitude, self%longitude, self%location, self%code_ref,    &
                     self%latitude_ref, self%longitude_ref, self%location_ref, self%component, &
                self%real, self%imaginary, self%error ) )
        !
    end function getCopyDataEntryMT_REF
    !
    subroutine writeDataEntryMT_REF( self )
        class( DataEntryMT_REF_t ), intent(in)    :: self
        !
        write(*,*) "Write DataEntryMT_REF_t: ", self%id, self%type, self%period, self%code,    &
        self%latitude, self%longitude, self%location, self%code_ref, self%latitude_ref, self%longitude_ref,    &
        self%location_ref, self%component, self%real, self%imaginary, self%error
        !
    end subroutine writeDataEntryMT_REF
    !
end module DataEntryMT_REF
