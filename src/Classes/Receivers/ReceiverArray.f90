!*************
!
! Class to provirx a dynamic and polymorphic array of Receivers
!
! Last modified at 05/2021 by Paulo Werdt
!
!*************
!
module ReceiverArray
   !
   use Receiver
   !
   implicit none
   !
   type, private :: Element_t
      !
      class( Receiver_t ), pointer :: receiver
      type( Element_t ), pointer   :: next => null()
      type( Element_t ), pointer   :: prev => null()
      !
   end type Element_t
   !
   type, public :: ReceiverArray_t
      !
      private
      !
      type( Element_t ), pointer   :: first
      type( Element_t ), pointer   :: last
      !
   contains
      !
      final :: ReceiverArray_dtor
      !
      procedure, public :: size   => getSizeReceiverArray
      procedure, public :: has   => hasReceiver
      procedure, public :: add   => addReceiver
      procedure, public :: get   => getReceiver
      !
   end type ReceiverArray_t
   !
   interface ReceiverArray_t
      module procedure ReceiverArray_ctor
   end interface ReceiverArray_t
   !
contains
   !
   function ReceiverArray_ctor() result( self )
      implicit none
      !
      class( ReceiverArray_t ), pointer   :: self
      !
      ! write(*,*) "Constructor ReceiverArray_t"
      !
      allocate( ReceiverArray_t :: self )
      !
      self%first   => null()
      self%last   => null()
      !
   end function ReceiverArray_ctor
   !
   subroutine ReceiverArray_dtor( self )
      implicit none
      !
      type( ReceiverArray_t ), intent( in out )   :: self
      !
      type( Element_t ), pointer            :: element
      !
      ! write(*,*) "Destructor ReceiverArray_t"
      !
      element => self%first
      do while( associated( element ) )
         deallocate( element%receiver )
         element => element%next
      end do
      !
   end subroutine ReceiverArray_dtor
   !
   function getSizeReceiverArray( self ) result( counter )
      implicit none
      !
      class( ReceiverArray_t ), intent( in )   :: self
      type( Element_t ), pointer            :: element
      integer                           :: counter
      !
      counter = 0
      element => self%first
      do while( associated( element ) )
         counter = counter + 1
         element => element%next
      end do
      !
   end function getSizeReceiverArray
   !
   function hasReceiver( self, receiver ) result( exist )
      class( ReceiverArray_t ), intent( in )   :: self
      class( Receiver_t ), intent( in )   :: receiver
      logical                        :: exist
      integer                        :: iRx, nRx
      !
      exist = .FALSE.
      !
      nRx = self%size()
      !
      do iRx = 1, nRx
         if( receiver%isEqual( self%get( iRx ) ) ) then
            exist = .TRUE.
            exit
         end if
      end do
      !
   end function hasReceiver
   !
   subroutine addReceiver( self, receiver )
      implicit none
      !
      class( ReceiverArray_t )   , intent(in out):: self
      class( Receiver_t ), pointer, intent( in )   :: receiver
      !
      type( Element_t ), pointer            :: element
      !
      allocate( element )
      !
      element%receiver => receiver
      element%next => null()
      !
      if( .not.associated( self%first ) ) then  
         element%prev => null()
         !
         self%first => element
         self%last => element
      else
         self%last%next => element  
         element%prev => self%last
         !
         self%last => element
      end if
      !
   end subroutine addReceiver
   !
   function getReceiver( self, index ) result( receiver )
      implicit none
      !
      class( ReceiverArray_t ), intent( in )   :: self
      integer, intent( in )                  :: index
      !
      type( Element_t ), pointer         :: element
      class( Receiver_t ), pointer         :: receiver
      integer                           :: counter
      !
      element => self%first
      counter = 1
      do while(associated( element ) )
         if( counter == index ) exit
       counter = counter + 1
         element => element%next
      end do
      !
      if(.not.associated(element)) then
         STOP 'ReceiverArray.f08: Receiver index not found.'
      end if
      !
      receiver => element%receiver
      !
   end function getReceiver
   !
end module ReceiverArray
