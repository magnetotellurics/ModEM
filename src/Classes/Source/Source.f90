module Source
  !
  use cVector
  use ModelOperator
  !
  type, abstract :: Source_t
     !
     class( ModelOperator_t ), pointer :: model_operator
     class( cVector_t ), pointer       :: rhs, bdry, E, e0
     !
     logical                           :: non_zero_source = .false., adjt = .false.
     !
   contains
     !
     procedure, public :: init    => initializeSource
     procedure, public :: dealloc => deallocateSource
     !
     procedure( interface_set_rhs ), deferred, public :: setRHS
     
  end type Source_t
  !
  abstract interface
     !
     subroutine interface_set_rhs(self)
       !
       import :: Source_t
       class( Source_t ), intent( inout ) :: self
       !
     end subroutine interface_set_rhs
     !
  end interface
  !
   contains
   !
   subroutine initializeSource( self )
      class( Source_t ), intent( inout ) :: self
      !
      self%rhs  => null()
      self%bdry => null()
      self%E    => null()
      self%e0   => null()
      !
   end subroutine initializeSource
   !
   subroutine deallocateSource( self )
      class( Source_t ), intent( inout ) :: self
      !
      deallocate( self%rhs )
      deallocate( self%bdry )
      deallocate( self%E )
      deallocate( self%e0 )
      !
   end subroutine deallocateSource
   !
end module Source
