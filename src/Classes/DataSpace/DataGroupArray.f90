!*************
!
! Class to providg a dynamic and polymorphic array of Data Groups
!
! Last modified at 05/2021 by Paulo Werdt
!
!*************
!
module DataGroupArray
   !
   use DataGroup
   !
   implicit none
   !
   type, private :: Element_t
      !
      class( DataGroup_t ), allocatable :: data_group
      type( Element_t ), pointer :: next => null()
      type( Element_t ), pointer :: prev => null()
      !
   end type Element_t
   !
   type, public :: DataGroupArray_t
      !
      type( Element_t ), pointer   :: first
      type( Element_t ), pointer   :: last
      !
      contains
         !
         final :: DataGroupArray_dtor
         !
         procedure, public :: size
         procedure, public :: has => hasDataGroup
         procedure, public :: add => addDataGroup
         procedure, public :: get => getDataGroup
         !
   end type DataGroupArray_t
   !
   interface DataGroupArray_t
      module procedure DataGroupArray_ctor
   end interface DataGroupArray_t
   !
contains
   !
   function DataGroupArray_ctor() result( self )
      implicit none
      !
      type( DataGroupArray_t ) :: self
      !
      self%first   => null()
      self%last   => null()
      !
   end function DataGroupArray_ctor
   !
   subroutine DataGroupArray_dtor( self )
      implicit none
      !
      type( DataGroupArray_t ), intent( inout ) :: self
      !
      type( Element_t ), pointer                :: element
      !
      write(*,*) "Destructor DataGroupArray_t"
      !
      element => self%first
      do while( associated( element ) )
         deallocate( element%data_group )
         element => element%next
      end do
      !
   end subroutine DataGroupArray_dtor
   !
   function size( self ) result( counter )
      implicit none
      !
      class( DataGroupArray_t ), intent( in ) :: self
	  integer                                 :: counter
      !
      type( Element_t ), pointer :: element
      !
      counter = 0
      element => self%first
      do while( associated( element ) )
         counter = counter + 1
         element => element%next
      end do
      !
   end function size
   !
   function hasDataGroup( self, data_group ) result( exist )
      class( DataGroupArray_t ), intent( in ) :: self
      class( DataGroup_t ), intent( in )      :: data_group
      logical                                 :: exist
	  !
      integer :: iDg, nDg
      !
      exist = .FALSE.
      !
      nDg = self%size()
      !
      do iDg = 1, nDg
         if( data_group%isEqual( self%get( iDg ) ) ) then
            exist = .TRUE.
            exit
         end if
      end do
      !
   end function hasDataGroup
   !
   subroutine addDataGroup( self, data_group )
      implicit none
      !
      class( DataGroupArray_t )   , intent( in out ) :: self
      class( DataGroup_t ), intent( in )             :: data_group
      !
      type( Element_t ), pointer :: element
      !
      allocate( element )
      !
      element%data_group = data_group
	  !
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
   end subroutine addDataGroup
   !
   function getDataGroup( self, index ) result( data_group )
      implicit none
      !
      class( DataGroupArray_t ), intent( in ) :: self
      integer, intent( in )                   :: index
	  class( DataGroup_t ), allocatable       :: data_group
      !
      type( Element_t ), pointer :: element
      integer                    :: counter
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
         STOP 'DataGroupArray.f08: DataGroup index not found.'
      end if
      !
      data_group = element%data_group
      !
   end function getDataGroup
   !
end module DataGroupArray
