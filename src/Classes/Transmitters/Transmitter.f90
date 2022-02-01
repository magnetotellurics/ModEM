! *************
! 
! Base class to define a Transmitter
! 
! Last modified at 10/11/2021 by Paulo Werdt
! 
! *************
! 
module Transmitter
   ! 
   use Constants
   use Source
   use ForwardSolver
   use IntegerArray
   use VectorArray
   !
   type, abstract :: Transmitter_t
      !
      integer           :: id, n_pol, fwd_key(8)
      !
      real( kind=prec ) :: period
      !
      class( ForwardSolver_t ), pointer :: forward_solver 
      class( Source_t ), pointer        :: source
      !
      class( cVector_t ), allocatable   :: e_all(:)
      class( IntegerArray_t ), pointer  :: receiver_indexes
      !
      character(:), allocatable :: DATA_TITLE
      !
   contains
      !
      procedure, public :: init    => initializeTx
      procedure, public :: dealloc => deallocateTx
      !
	  procedure, public :: setSource => setSourceTx
	  procedure, public :: setForwardSolver => setForwardSolverTx
	  !
      procedure, public :: updateFwdKey
      !
      procedure, public :: has    => hasReceiverTx
      procedure, public :: add    => addReceiverTx
      procedure, public :: get    => getReceiverTx
      procedure, public :: getNRx => getNumberOfReceivers
      !
      procedure( interface_solve_fwd_tx ), deferred, public  :: solveFWD
      !
      procedure( interface_get_type_tx ), deferred, public   :: getType
      procedure( interface_is_equal_tx ), deferred, public   :: isEqual
      procedure( interface_write_tx ), deferred, public      :: write
      !
   end type Transmitter_t
   !
   abstract interface
      !
      subroutine interface_solve_fwd_tx( self )
         import :: Transmitter_t
         class( Transmitter_t ), intent( inout ) :: self
      end subroutine interface_solve_fwd_tx
      !
      function interface_get_type_tx( self ) result( type )
         import :: Transmitter_t
         class( Transmitter_t ), intent( in ) :: self
         character(:), allocatable            :: type
      end function interface_get_type_tx
      !
      function interface_is_equal_tx( self, other ) result( equal )
         import :: Transmitter_t
         class( Transmitter_t ), intent( in ) :: self, other
         logical                              :: equal
      end function interface_is_equal_tx
      !
      subroutine interface_write_tx( self )
         import :: Transmitter_t
         class( Transmitter_t ), intent(in) :: self
      end subroutine interface_write_tx
      !
   end interface
   !
   contains
   !
   subroutine initializeTx( self )
      implicit none
      !
      class( Transmitter_t ), intent( inout ) :: self
      !
      call self%updateFwdKey()
      !
      self%forward_solver => null()
      self%source         => null()
      !
      allocate( self%receiver_indexes, source = IntegerArray_t() )
      !
   end subroutine initializeTx
   !
   subroutine deallocateTx( self )
      implicit none
      !
      class( Transmitter_t ), intent( inout ) :: self
      !
      !if( associated( self%forward_solver ) ) deallocate( self%forward_solver )
      !if( associated( self%source ) ) deallocate( self%source )
      !
      !if( associated( self%e_all ) ) deallocate( self%e_all )
      !deallocate( self%receiver_indexes )
      !
   end subroutine deallocateTx
   !
   subroutine updateFwdKey( self )
      implicit none
      !
      class( Transmitter_t ), intent( inout ) :: self
      !
      call date_and_time( values=self%fwd_key )
      !
   end subroutine updateFwdKey
   !
   function hasReceiverTx( self, receiver_index ) result( found )
      implicit none
      !
      class( Transmitter_t ), intent( in ) :: self
      integer, intent( in )                :: receiver_index
      !
      integer :: i_rx, n_rx
      logical :: found
      !
      found = .FALSE.
      !
      n_rx = self%receiver_indexes%size()
      !
      do i_rx = 1, n_rx
         !
         if( receiver_index == self%receiver_indexes%get( i_rx ) ) then
            found = .TRUE.
         end if
      end do
      !
   end function hasReceiverTx
   !
   subroutine addReceiverTx( self, receiver_index )
      implicit none
      !
      class( Transmitter_t ), intent( inout ) :: self
      integer, intent( in )                   :: receiver_index
      !
      call self%receiver_indexes%add( receiver_index )
      !
   end subroutine addReceiverTx
   !
   function getReceiverTx( self, index ) result( receiver_index )
      implicit none
      !
      class( Transmitter_t ), intent( in ) :: self
      integer                              :: index, receiver_index
      !
      receiver_index = self%receiver_indexes%Get( index )
      !
   end function getReceiverTx
   !
   function getNumberOfReceivers( self ) result( counter )
      implicit none
      !
      class( Transmitter_t ), intent( in ) :: self
      integer                              :: counter
      !
      counter = self%receiver_indexes%size()
      !
   end function getNumberOfReceivers
   !
   subroutine setForwardSolverTx( self, forward_solver )
      !
      class( Transmitter_t ), intent( inout ) :: self
      class( ForwardSolver_t ), target, intent( in )   :: forward_solver
      !
      self%forward_solver => forward_solver
      !
   end subroutine setForwardSolverTx
   !
   subroutine setSourceTx( self, source )
      !
      class( Transmitter_t ), intent( inout ) :: self
      class( Source_t ), target, intent( in ) :: source
      !
      self%source => source
      !
   end subroutine setSourceTx
   !
end module Transmitter
