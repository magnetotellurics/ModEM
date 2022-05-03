!
!
!
module DataHandleMT
    !
    use DataHandle
    !
    type, extends( DataHandle_t ), public :: DataHandleMT_t
        !
    end type DataHandleMT_t
    !
    interface DataHandleMT_t
        module procedure DataHandleMT_ctor
    end interface DataHandleMT_t
    !
contains
    !
    function DataHandleMT_ctor( rx_type, code, component, period, rx_location, real, imaginary ) result( self )
        implicit none
        !
        integer, intent( in )                   :: rx_type
        character(:), allocatable, intent( in ) :: code, component
        real( kind=prec ), intent( in )         :: period, rx_location(3), real, imaginary
        !
        type( DataHandleMT_t ) :: self
        !
        !write(*,*) "Constructor DataHandleMT"
        !
        self%rx_type     = rx_type
        self%code        = code
        self%component   = component
        self%period      = period
        self%rx_location = rx_location
        self%real        = real
        self%imaginary   = imaginary
        !
    end function DataHandleMT_ctor
	!
end module DataHandleMT
