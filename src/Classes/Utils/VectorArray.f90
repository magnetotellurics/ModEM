!*************
!
! Class to provide a dynamic and polymorphic array of cVector_t
!
! Last modified at 05/2021 by Paulo Werdt
!
!*************
!
module VectorArray
   !
   use cVector
   !
   implicit none
   !
   type, private :: Element_t
   !
   class( cVector_t ), pointer   :: vector
   type( Element_t ), pointer   :: next => null()
   type( Element_t ), pointer   :: prev => null()
   !
   end type Element_t
   !
   type, public :: VectorArray_t
   !
   private
   !
   type( Element_t ), pointer   :: first
   type( Element_t ), pointer   :: last
   !
   contains
      !
      final :: VectorArray_dtor
      !
      procedure, public :: size   => getSizeVectorArray
      procedure, public :: add   => addVector
      procedure, public :: get   => getVector
      !
   end type VectorArray_t
   !
   interface VectorArray_t
      module procedure VectorArray_ctor
   end interface VectorArray_t
   !
contains
   !
   function VectorArray_ctor() result( self )
      implicit none
      !
      ! Local variables
      class( VectorArray_t ), pointer :: self
      !
      ! write(*,*) "Constructor VectorArray_t"
      !
      allocate( VectorArray_t :: self )
      !
      self%first   => null()
      self%last   => null()
      !
   end function VectorArray_ctor
   !
   subroutine VectorArray_dtor( self )
      implicit none
      !
      type( VectorArray_t ), intent( in out ) :: self
      !
      type( Element_t ), pointer   :: element
      !
      ! write(*,*) "Destructor VectorArray_t"
      !
      element => self%first
      do while( associated( element ) )
         deallocate( element%Vector )
         element => element%next
      end do
      !
   end subroutine VectorArray_dtor
   !
   function getSizeVectorArray( self ) result( counter )
      implicit none
      ! Arguments
      class( VectorArray_t ), intent( in )   :: self
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
   end function getSizeVectorArray
   !
   subroutine addVector( self, vector )
      implicit none
      !
      class( VectorArray_t )   , intent( in out )   :: self
      class( cVector_t ), allocatable, intent( in )   :: vector
      !
      type( Element_t ), pointer      :: element
      !
      allocate( element )
      !
      allocate( element%vector, source = vector )
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
   end subroutine addVector
   !
   function getVector( self, index ) result( vector )
      implicit none
      !
      class( VectorArray_t ), intent( in ) :: self
      integer, intent( in )                :: index
      !
      type( Element_t ), pointer      :: element
      class( cVector_t ), allocatable :: vector
      integer                         :: counter
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
         STOP 'VectorArray.f08: Vector index not found.'
      end if
      !
      allocate( vector, source = element%vector )
      !
   end function getVector
   !
end module VectorArray
