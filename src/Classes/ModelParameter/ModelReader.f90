module ModelReader
   !
   use Grid
   use ModelParameter
   !
   type, abstract, public :: ModelReader_t
    contains
       procedure( interface_read_model_reader ) , deferred :: Read
       procedure( interface_write_model_reader ), deferred :: Write
   end type ModelReader_t
   !
   abstract interface
       subroutine interface_read_model_reader(self, fileName, grid, model)
          import :: ModelReader_t, Grid_t, ModelParameter_t
          ! Arguments
          class( ModelReader_t ), intent( in )      :: self
          character(*), intent( in )                :: fileName
          class( Grid_t ), allocatable, intent(out) :: grid
          class( ModelParameter_t ), allocatable, intent(out) :: model
       end subroutine interface_read_model_reader

       subroutine interface_write_model_reader(self, fileName, grid, model)
          import :: ModelReader_t, Grid_t, ModelParameter_t
          ! Arguments
          class( ModelReader_t ), intent(in)    :: self
          character(*), intent(in)              :: fileName
          class( Grid_t), intent(in)            :: grid
          class( ModelParameter_t ), intent(in) :: model
       end subroutine interface_write_model_reader
   
   end interface
end module ModelReader
