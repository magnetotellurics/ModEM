!
!
!
module DataHandleCSEM
    !
    use DataHandle
    !
    type, extends( DataHandle_t ), public :: DataHandleCSEM_t
        !
        real( kind=prec )  :: tx_location(3), azimuth, dip, moment
        character(:), allocatable :: dipole
        !
    end type DataHandleCSEM_t
    !
    interface DataHandleCSEM_t
        module procedure DataHandleCSEM_ctor
    end interface DataHandleCSEM_t
    !
contains
    !
    function DataHandleCSEM_ctor( rx_type, code, component, period, tx_location, azimuth, dip, moment, dipole, rx_location, real, imaginary ) result( self )
        implicit none
        !
        integer, intent( in )                   :: rx_type
        character(:), allocatable, intent( in ) :: code, component
        real( kind=prec ), intent( in )         :: tx_location(3), azimuth, dip, moment
        real( kind=prec ), intent( in )         :: period, rx_location(3), real, imaginary
        character(:), allocatable :: dipole
        !
        !
        type( DataHandleCSEM_t ) :: self
        !
        !write(*,*) "Constructor DataHandleCSEM"
        !
        self%rx_type     = rx_type
        self%code        = code
        self%component   = component
        self%period      = period
        self%tx_location = tx_location
        self%azimuth     = azimuth
        self%dip         = dip
        self%moment      = moment
        self%dipole      = dipole
        self%rx_location = rx_location
        self%real        = real
        self%imaginary   = imaginary
        !
    end function DataHandleCSEM_ctor
    !
end module DataHandleCSEM
