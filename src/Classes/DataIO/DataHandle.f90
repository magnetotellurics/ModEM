!
!
!
module DataHandle
    !
    use Constants
    !
    ! Global file name for predicted data file
    character(:), allocatable :: predicted_data_file_name
    !
    type, abstract :: DataHandle_t
        !
        integer                   :: rx_type
        character(:), allocatable :: code, component
        real( kind=prec )         :: period, rx_location(3)
        real( kind=prec )         :: rvalue, imaginary
        !
    end type DataHandle_t
    !
end module DataHandle
