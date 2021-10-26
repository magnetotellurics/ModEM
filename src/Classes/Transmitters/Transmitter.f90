! *************
! 
! Base class to define a Transmitter
! 
! Last modified at 08/06/2021 by Paulo Werdt
! 
! *************
! 
module Transmitter
   ! 
   use Constants
   use Source
   use ForwardSolverFromFile
   use IntegerArray
   use VectorArray
   use ModelOperator
   !
   type, abstract :: Transmitter_t
      !
      integer           :: id, n_pol, fwd_key(8)
      !
      real( kind=prec ) :: period
      !
      class( ForwardSolver_t ), pointer :: forward_solver 
      class( VectorArray_t ), pointer   :: e_solution
      class( Source_t ), allocatable    :: source
      !
      class( IntegerArray_t ), pointer  :: receiver_indexes
      !
	  character(:), allocatable :: DATA_TITLE
   contains
      !
      procedure, public :: init    => initializeTx
      procedure, public :: dealloc => deallocateTx
      !
      procedure, public :: updateFwdKey
      !
      procedure, public :: has    => hasReceiverTx
      procedure, public :: add    => addReceiverTx
      procedure, public :: get    => getReceiverTx
      procedure, public :: getNRx => getNumberOfReceivers
      !
      procedure( interface_solve_fwd_tx ), deferred, public  :: solveFWD
      procedure( interface_get_source_tx ), deferred, public :: getSource
      !
	  procedure( interface_get_type_tx ), deferred, public   :: getType
      procedure( interface_is_equal_tx ), deferred, public   :: isEqual
      procedure( interface_write_tx ), deferred, public      :: write
      !
   end type Transmitter_t
   !
   abstract interface
      !
      subroutine interface_solve_fwd_tx( self, model_operator )
         import :: Transmitter_t, ModelOperator_t
         class( Transmitter_t ), intent( inout )             :: self
         class( ModelOperator_t ), allocatable, intent( in ) :: model_operator
      end subroutine interface_solve_fwd_tx
      !
      subroutine interface_get_source_tx( self )
         import :: Transmitter_t
         class( Transmitter_t ), intent( in ) :: self
      end subroutine interface_get_source_tx
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
      class( Transmitter_t ), intent( inout ) :: self
      !
      call self%updateFwdKey()
      !
      self%e_solution => VectorArray_t()
      !
      self%receiver_indexes => IntegerArray_t()
      !
   end subroutine initializeTx
   !
   subroutine deallocateTx( self )
      class( Transmitter_t ), intent( inout ) :: self
      !
      deallocate( self%forward_solver )
      deallocate( self%e_solution )
      deallocate( self%source )
      !
      deallocate( self%receiver_indexes )
      !
   end subroutine deallocateTx
   !
   subroutine updateFwdKey( self )
      class( Transmitter_t ), intent( inout ) :: self
      !
      call date_and_time( values=self%fwd_key )
      !
   end subroutine updateFwdKey
   !
   function hasReceiverTx( self, receiver_index ) result( found )
      class( Transmitter_t ), intent( in ) :: self
      integer, intent( in )                :: receiver_index
      !
      logical :: found
      !
      found = .FALSE.
      !
      nRx = self%receiver_indexes%size()
      !
      do iRx = 1, nRx
         !
         if( receiver_index == self%receiver_indexes%get( iRx ) ) then
            found = .TRUE.
         end if
      end do
      !
   end function hasReceiverTx
   !
   subroutine addReceiverTx( self, receiver_index )
      class( Transmitter_t ), intent( inout )      :: self
      integer, pointer, intent( in )            :: receiver_index
      !
      call self%receiver_indexes%add( receiver_index )
      !
   end subroutine addReceiverTx
   !
   function getReceiverTx( self, index ) result( receiver_index )
      class( Transmitter_t ), intent( in )   :: self
      integer                            :: index, receiver_index
      !
      receiver_index = self%receiver_indexes%Get( index )
      !
   end function getReceiverTx
   !
   function getNumberOfReceivers( self ) result( counter )
      class( Transmitter_t ), intent( in )   :: self
      integer                            :: counter
      !
      counter = self%receiver_indexes%size()
      !
   end function getNumberOfReceivers
   !
end module Transmitter
