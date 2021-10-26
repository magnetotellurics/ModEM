module ModelReader
  use Grid
  use ModelParameter
  
  implicit none

  private

  type, abstract, public :: ModelReader_t
   contains
     procedure(iface_Read) , deferred :: Read
     procedure(iface_Write), deferred :: Write
  end type ModelReader_t

  abstract interface
     subroutine iface_Read(self, fileName, grid, model)
       import :: ModelReader_t, Grid_t, ModelParameter_t
       ! Arguments
       class(ModelReader_t)   , intent(in) :: self
       character(*)           , intent(in) :: fileName
       class(Grid_t)          , pointer, intent(out) :: grid
       class(ModelParameter_t), pointer, intent(out) :: model
     end subroutine iface_Read

     subroutine iface_Write(self, fileName, grid, model)
       import :: ModelReader_t, Grid_t, ModelParameter_t
       ! Arguments
       class(ModelReader_t)   , intent(in) :: self
       character(*)           , intent(in) :: fileName
       class(Grid_t)          , intent(in) :: grid
       class(ModelParameter_t), intent(in) :: model
     end subroutine iface_Write
  
  end interface
end module ModelReader
