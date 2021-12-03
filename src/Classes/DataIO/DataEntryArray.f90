!*************
!
! Class to provide a dynamic and polymorphic array of Data entries
!
! Last modified at 05/2021 by Paulo Werdt
!
!*************
!
module DataEntryArray
   !
   use DataEntry
   !
   implicit none
   !
   type, private :: Element_t
      !
      class( DataEntry_t ), allocatable :: data_entry
      type( Element_t ), pointer :: next => null()
      type( Element_t ), pointer :: prev => null()
      !
   end type Element_t
   !
   type, public :: DataEntryArray_t
      !
      type( Element_t ), pointer :: first
      type( Element_t ), pointer :: last
      !
      contains
        !
        final :: DataEntryArray_dtor
        !
        procedure, public :: size  => getSizeDataEntryArray
        procedure, public :: has   => hasDataEntry
        procedure, public :: add   => addDataEntry
        procedure, public :: get   => getDataEntry
        !
   end type DataEntryArray_t
   !
   interface DataEntryArray_t
      module procedure DataEntryArray_ctor
   end interface DataEntryArray_t
   !
contains
   !
   function DataEntryArray_ctor() result( self )
      implicit none
      !
      type( DataEntryArray_t ) :: self
      !
      !write(*,*) "Constructor DataEntryArray_t"
      !
      self%first => null()
      self%last  => null()
      !
   end function DataEntryArray_ctor
   !
   subroutine DataEntryArray_dtor( self )
      implicit none
      !
      type( DataEntryArray_t ), intent( in out ) :: self
      !
      type( Element_t ), pointer :: element
      !
      write(*,*) "Destructor DataEntryArray_t"
      !
      element => self%first
      do while( associated( element ) )
         deallocate( element%data_entry )
         element => element%next
      end do
      !
      self%first => null()
      self%last  => null()
      !
   end subroutine DataEntryArray_dtor
   !
   function getSizeDataEntryArray( self ) result( counter )
      implicit none
      !
      class( DataEntryArray_t ), intent( in ) :: self
	  integer                                 :: counter
      !
      type( Element_t ), pointer  :: element
      !
      counter = 0
      element => self%first
      do while( associated( element ) )
         counter = counter + 1
         element => element%next
      end do
      !
   end function getSizeDataEntryArray
   !
   function hasDataEntry( self, data_entry ) result( exist )
      class( DataEntryArray_t ), intent( in ) :: self
      class( DataEntry_t ), intent( in )      :: data_entry
      !
      logical                        :: exist
      integer                        :: iDe, nDe
      !
      exist = .FALSE.
      !
      nDe = self%size()
      do iDe = 1, nDe
         if( data_entry%isEqual( self%get( iDe ) ) ) then
            exist = .TRUE.
            exit
         end if
      end do
      !
   end function hasDataEntry
   !
   subroutine addDataEntry( self, data_entry )
      implicit none
      !
      class( DataEntryArray_t ), intent( inout ) :: self
      class( DataEntry_t ), intent( in )         :: data_entry
      !
      type( Element_t ), pointer :: element
      !
      allocate( element )
      !
      element%data_entry = data_entry
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
   end subroutine addDataEntry
   !
   function getDataEntry( self, de_index ) result( data_entry )
      implicit none
      !
      class( DataEntryArray_t ), intent( in ) :: self
      integer, intent( in )                   :: de_index
	  class( DataEntry_t ), pointer           :: data_entry
      !
      type( Element_t ), pointer :: element
      integer                    :: counter
      !
      element => self%first
      counter = 1
      do while( associated( element ) )
	     if( counter == de_index ) exit
         counter = counter + 1
         element => element%next
      end do
      !
      if( .not. associated( element ) ) then
         STOP 'DataEntryArray.f08: DataEntry index not found.'
      end if
      !
      allocate( data_entry, source = element%data_entry )
      !
   end function getDataEntry
   !
end module DataEntryArray
