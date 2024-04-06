!
!> Abstract Base class to define a ModelReader
!
module ModelReader
    !
    use Grid
    use ModelParameter
    !
    type, abstract, public :: ModelReader_t
        !
        !> No base properties
        !
        contains
            !
            procedure( interface_read_model_reader ), deferred :: read
            !
    end type ModelReader_t
    !
    abstract interface
        !
        !> No interface subroutine briefing
        subroutine interface_read_model_reader( self, file_name, grid, model, param_grid )
            import :: ModelReader_t, Grid_t, ModelParameter_t
            !
            class( ModelReader_t ), intent( in ) :: self
            character(*), intent( in ) :: file_name
            class( Grid_t ), allocatable, intent( out ) :: grid
            class( ModelParameter_t ), allocatable, intent( out ) :: model
            class( Grid_t ), allocatable, intent( out ), optional :: param_grid
            !
        end subroutine interface_read_model_reader
        !
    end interface
    !
end module ModelReader
!