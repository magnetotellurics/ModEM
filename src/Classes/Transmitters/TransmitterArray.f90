!*************
!
! Class to provitx a dynamic and polymorphic array of Transmitters
!
! Last modified at 05/2021 by Paulo Werdt
!
!*************
!
module TransmitterArray
   !
   use Transmitter
   !
   implicit none
   !
   type, private :: Element_t
   !
   class( Transmitter_t ), pointer   :: transmitter
   type( Element_t ), pointer   :: next => null()
   type( Element_t ), pointer   :: prev => null()
   !
   end type Element_t
   !
   type, public :: TransmitterArray_t
   !
   private
   !
   type( Element_t ), pointer   :: first
   type( Element_t ), pointer   :: last
   !
   contains
      !
      final :: TransmitterArray_dtor
      !
      procedure, public :: size   => getSizeTransmitterArray
      procedure, public :: has   => hasTransmitter
      procedure, public :: add   => addTransmitter
      procedure, public :: get   => getTransmitter
      !
   end type TransmitterArray_t
   !
   interface TransmitterArray_t
      module procedure TransmitterArray_ctor
   end interface TransmitterArray_t
   !
contains
   !
   function TransmitterArray_ctor() result( self )
      implicit none
      !
      ! Local variables
      class( TransmitterArray_t ), pointer :: self
      !
      ! write(*,*) "Constructor TransmitterArray_t"
      !
      allocate( TransmitterArray_t :: self )
      !
      self%first   => null()
      self%last   => null()
      !
   end function TransmitterArray_ctor
   !
   subroutine TransmitterArray_dtor( self )
      implicit none
      !
      type( TransmitterArray_t ), intent( in out ) :: self
      !
      type( Element_t ), pointer   :: element
      !
      !write(*,*) "Destructor TransmitterArray_t"
      !
      element => self%first
      do while( associated( element ) )
         deallocate( element%transmitter )
         element => element%next
      end do
      !
   end subroutine TransmitterArray_dtor
   !
   function getSizeTransmitterArray( self ) result( counter )
      implicit none
      ! Arguments
      class( TransmitterArray_t ), intent( in )   :: self
      ! Local variables
      type( Element_t ), pointer      :: element
      integer                     :: counter
      !
      counter = 0
      element => self%first
      do while( associated( element ) )
         counter = counter + 1
         element => element%next
      end do
      !
   end function getSizeTransmitterArray
   !
   function hasTransmitter( self, transmitter ) result( exist )
      class( TransmitterArray_t ), intent( in )   :: self
      class( Transmitter_t ), intent( in )   :: transmitter
      !
      logical                        :: exist
      integer                        :: iTx, nTx
      !
      exist = .FALSE.
      !
      nTx = self%size()
      !
      do iTx = 1, nTx
         if( transmitter%isEqual( self%get( iTx ) ) ) then
            exist = .TRUE.
            exit
         end if
      end do
      !
   end function hasTransmitter
   !
   subroutine addTransmitter( self, transmitter )
      implicit none
      !
      class( TransmitterArray_t )   , intent( in out )   :: self
      class( Transmitter_t ), pointer, intent( in )   :: transmitter
      !
      type( Element_t ), pointer      :: element
      !
      allocate( element )
      !
      element%transmitter => transmitter
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
   end subroutine addTransmitter
   !
   function getTransmitter( self, index ) result( transmitter )
      implicit none
      !
      class( TransmitterArray_t ), intent( in )   :: self
      integer, intent( in )                  :: index
      !
      type( Element_t ), pointer            :: element
      class( Transmitter_t ), pointer         :: transmitter
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
         STOP 'TransmitterArray.f08: Transmitter index not found.'
      end if
      !
      transmitter => element%transmitter
      !
   end function getTransmitter
   !
end module TransmitterArray
