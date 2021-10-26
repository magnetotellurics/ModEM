!*************
!
! Class to provide a dynamic array of Integers
!
! Last modified at 06/2021 by Paulo Werdt
!
!*************
!
module IntegerArray
   !
   implicit none
   !
   type, private :: Element_t
   !
   integer                  :: int_value
   type( Element_t ), pointer   :: next => null()
   type( Element_t ), pointer   :: prev => null()
   !
   end type Element_t
   !
   type, public :: IntegerArray_t
   !
   private
   !
   type( Element_t ), pointer   :: first
   type( Element_t ), pointer   :: last
   !
   contains
      !
      final :: IntegerArray_dtor
      !
      procedure, public :: size   => getSizeIntegerArray
      procedure, public :: has   => hasInteger
      procedure, public :: add   => addInteger
      procedure, public :: get   => getInteger
      !
   end type IntegerArray_t
   !
   interface IntegerArray_t
      module procedure IntegerArray_ctor
   end interface IntegerArray_t
   !
contains
   !
   function IntegerArray_ctor() result( self )
      implicit none
      !
      ! Local variables
      class( IntegerArray_t ), pointer :: self
      !
      ! write(*,*) "Constructor IntegerArray_t"
      !
      allocate( IntegerArray_t :: self )
      !
      self%first   => null()
      self%last   => null()
      !
   end function IntegerArray_ctor
   !
   subroutine IntegerArray_dtor( self )
      implicit none
      !
      type( IntegerArray_t ), intent( in out ) :: self
      !
      type( Element_t ), pointer   :: element
      !
      ! write(*,*) "Destructor IntegerArray_t"
      !
   end subroutine IntegerArray_dtor
   !
   function getSizeIntegerArray( self ) result( counter )
      implicit none
      ! Arguments
      class( IntegerArray_t ), intent( in )   :: self
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
   end function getSizeIntegerArray
   !
   function hasInteger( self, int_value ) result( exist )
      class( IntegerArray_t ), intent( in )   :: self
      integer, intent( in )               :: int_value
      !
      logical                        :: exist
      integer                        :: index, size
      !
      exist = .FALSE.
      !
      size = self%size()
      !
      do index = 1, size
         if( int_value == self%get( index ) ) then
            exist = .TRUE.
            exit
         end if
      end do
      !
   end function hasInteger
   !
   subroutine addInteger( self, int_value )
      implicit none
      !
      class( IntegerArray_t )   , intent( in out )   :: self
      integer, intent( in )                  :: int_value
      !
      type( Element_t ), pointer      :: element
      !
      allocate( element )
      !
      element%int_value = int_value
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
   end subroutine addInteger
   !
   function getInteger( self, index ) result( int_value )
      implicit none
      !
      class( IntegerArray_t ), intent( in )   :: self
      integer, intent( in )               :: index
      !
      type( Element_t ), pointer   :: element
      integer                  :: int_value, counter
      !
      element => self%first
      counter = 1
      do while( associated( element ) )
         if( counter == index ) exit
       counter = counter + 1
         element => element%next
      end do
      !
      if(.not.associated(element)) then
         STOP 'IntegerArray.f08: index not found.'
      end if
      !
      int_value = element%int_value
      !
   end function getInteger
   !
end module IntegerArray
