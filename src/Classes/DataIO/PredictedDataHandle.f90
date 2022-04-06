!
!
!
module PredictedDataHandle
    !
    use Constants
    !
    use DataEntryMT
    use DataEntryMT_REF
    use DataEntryCSEM
    !
    type :: PredictedDataHandle_t
        !
        character(:), allocatable :: code, component
        real( kind=prec )         :: period, xyz(3)
        real( kind=prec )         :: real, imaginary
        !
    end type PredictedDataHandle_t
    !
    public :: buildPredictedDataHandle
    !
contains
    !
    function buildPredictedDataHandle( code, component, period, xyz, real, imaginary ) result( data_entry_handle )
        implicit none
        !
        character(:), allocatable, intent( in ) :: code, component
        real( kind=prec ), intent( in )         :: period, xyz(3), real, imaginary
        !
        type( PredictedDataHandle_t ) :: data_entry_handle
        !
        data_entry_handle%code        = code
        data_entry_handle%component   = component
        data_entry_handle%period      = period
        data_entry_handle%xyz         = xyz
        data_entry_handle%real        = real
        data_entry_handle%imaginary   = imaginary
        !
    end function buildPredictedDataHandle
    !
end module PredictedDataHandle
