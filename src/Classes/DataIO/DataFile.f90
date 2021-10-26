!*************
!
! Base class to read a data file
!
! Last modified at 05/2021 by Paulo Werdt
!
!*************
!
module DataFile
   !
   use Constants
   use String
   !
   use DataEntryArray
   !
   type, abstract :: DataFile_t
      !
      class( DataEntryArray_t ), pointer :: data_entries
      !
   contains
      !
      procedure, public :: init    => initializeDataFile
      procedure, public :: dealloc => deallocateDataFile
      !
      procedure, public :: get => getDataEntryDataFile
      !
   end type DataFile_t
   !
contains
   !
   function getDataEntryDataFile( self, iDe ) result( data_entry )
      class( DataFile_t ), intent( in out )   :: self
      integer, intent( in )               :: iDe
      class( DataEntry_t ), pointer         :: data_entry
      !
      data_entry => self%data_entries%Get( iDe )
      !
   end function getDataEntryDataFile
   !
   subroutine initializeDataFile( self )
      class( DataFile_t ), intent( inout ) :: self
      !
      self%data_entries => DataEntryArray_t()
      !
   end subroutine initializeDataFile
   !
   subroutine deallocateDataFile( self )
      class( DataFile_t ), intent( inout ) :: self
      !
      deallocate( self%data_entries )
      !
   end subroutine deallocateDataFile
   !
end module DataFile
