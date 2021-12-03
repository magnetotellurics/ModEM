! *************
! 
! Base class to define a Source
! 
! Last modified at 10/11/2021 by Paulo Werdt
! 
! *************
! 
module Source
  !
  use cVector
  use ModelOperator
  use ModelParameter
  !
  type, abstract :: Source_t
     !
     real( kind=prec ) :: omega
     !
     class( ModelOperator_t ), pointer  :: model_operator
     class( ModelParameter_t ), pointer :: model_parameter
     class( cVector_t ), allocatable    :: rhs, bdry, e0, E
     !
     logical                            :: non_zero_source, adjt
     !
     contains
        !
        procedure, public :: init    => initializeSource
        procedure, public :: dealloc => deallocateSource
        !
        procedure( interface_set_rhs ), deferred, public :: setRHS
        procedure( interface_set_e0 ), deferred, public  :: setE0
        procedure( interface_set_e ), deferred, public   :: setE
        !
  end type Source_t
  !
  abstract interface
     !
     subroutine interface_set_rhs( self )
       !
       import :: Source_t
       class( Source_t ), intent( inout ) :: self
       !
     end subroutine interface_set_rhs
     !
     subroutine interface_set_e0( self )
       !
       import :: Source_t
       class( Source_t ), intent( inout ) :: self
       !
     end subroutine interface_set_e0
     !
     subroutine interface_set_e( self, omega, polarization )
       !
       import :: Source_t, prec
       class( Source_t ), intent( inout ) :: self
       real( kind=prec ), intent( in )    :: omega
       integer, intent( in )              :: polarization
       !
     end subroutine interface_set_e
     !
  end interface
  !
   contains
   !
   subroutine initializeSource( self )
      implicit none
      !
      class( Source_t ), intent( inout ) :: self
      !
      self%model_operator => null()
      self%model_parameter => null()
      !
      self%non_zero_source = .false.
      self%adjt            = .false.
      !
   end subroutine initializeSource
   !
   subroutine deallocateSource( self )
      implicit none
      !
      class( Source_t ), intent( inout ) :: self
      !
      if( allocated( self%rhs ) ) deallocate( self%rhs )
      if( allocated( self%bdry ) ) deallocate( self%bdry )
      if( allocated( self%e0 ) ) deallocate( self%e0 )
      if( allocated( self%E ) ) deallocate( self%E )
      !
   end subroutine deallocateSource
   !
end module Source
