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
   type, private :: Element_t
      !
      class( cVector_t ), allocatable :: vector
      class( Element_t ), pointer     :: next => null()
      class( Element_t ), pointer     :: prev => null()
      !
   end type Element_t
   !
   type, public :: VectorArray_t
	   !
	   private
	   !
	   class( Element_t ), pointer :: first
	   class( Element_t ), pointer :: last
	   !
	   contains
		  !
		  final :: VectorArray_dtor
		  !
		  procedure, public :: size => getSizeVectorArray
		  procedure, public :: add  => addVector
		  procedure, public :: get  => getVector
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
      type( VectorArray_t ) :: self
      !
      self%first => null()
      self%last  => null()
      !
   end function VectorArray_ctor
   !
   subroutine VectorArray_dtor( self )
      implicit none
      !
      type( VectorArray_t ), intent( inout ) :: self
      !
      class( Element_t ), pointer :: element
      !
      element => self%first
      do while( associated( element ) )
         deallocate( element%vector )
         element => element%next
      end do
      !
   end subroutine VectorArray_dtor
   !
   function getSizeVectorArray( self ) result( counter )
      implicit none
      !
      class( VectorArray_t ), intent( in ) :: self
      !
      class( Element_t ), pointer :: element
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
      class( VectorArray_t ), intent( inout ) :: self
      class( cVector_t ), intent( in )        :: vector
      !
      class( Element_t ), pointer :: element
      !
      allocate( element )
      !
      allocate( element%vector, source = vector )
	  !
      element%next => null()
      !
      if( .not. associated( self%first ) ) then  
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
      class( cVector_t ), pointer          :: vector
	  !
	  class( Element_t ), pointer :: element
      integer                     :: counter
      !
      element => self%first
      counter = 1
      do while(associated( element ) )
         if( counter == index ) exit
       counter = counter + 1
         element => element%next
      end do
      !
      if( .not. associated( element ) ) then
         STOP 'VectorArray.f08: Vector index not found.'
      end if
      !
      allocate( vector, source = element%vector )
      !
   end function getVector
   !
end module VectorArray
