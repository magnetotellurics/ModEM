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
   end type DataFile_t
   !
contains
   !
   subroutine initializeDataFile( self )
      class( DataFile_t ), intent( inout ) :: self
      !
      allocate( self%data_entries, source = DataEntryArray_t() )
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
