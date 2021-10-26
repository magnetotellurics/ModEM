!*************
!
! Class to read a data file and create an array with all data entries(lines)
!
! Last modified at 08/2021 by Paulo Werdt
!
!*************
!
module ModEMControlFile
   !
   use String
   !
   type :: ModEMControlFile_t
      !
	  character(:), allocatable :: grid_type
	  character(:), allocatable :: grid_reader
	  !
   contains
      !
      final :: ModEMControlFile_dtor
      !
   end type ModEMControlFile_t
   !
   interface ModEMControlFile_t
      module procedure ModEMControlFile_ctor
   end interface ModEMControlFile_t
   !
contains
   !
   ! Read line by line of the data file, create Data Entry objects (MT, MT_REF or CSEM)
   function ModEMControlFile_ctor( funit, fname ) result( self )
      implicit none
      !
      integer, intent( in )                   :: funit
      character(:), allocatable, intent( in ) :: fname
	  !
	  class( ModEMControlFile_t ), pointer    :: self
      !
      character(1000)                   :: full_line_text
      character(len=200), dimension(20) :: args
      character(:), allocatable         :: line_text
	  integer                           :: line_counter, io_stat, p_nargs
      !
      !write(*,*) "Constructor ModEMControlFile_t"
      !
      allocate( ModEMControlFile_t :: self )
      !
      call Compact( fname )
      !
      open( unit = funit, file = fname, iostat = io_stat, status = 'old' )
      !
      if( io_stat /= 0 ) then
         write(*,*) 'Unable to open [', fname, '], Stat: ', io_stat
      else
         line_counter = getLineNumber( funit )
         !
         rewind( funit )
         !
         do
            read( funit, '(a)', END = 10 ) full_line_text
            line_text = adjustl( full_line_text )
            line_text = trim( line_text )
            !
            call Parse( line_text, ":", args, p_nargs )
            !
            if( index( line_text, "#" ) == 0 .and. index( line_text, ">" ) == 0 ) then
                !
				if( index( line_text, "Grid Type" ) > 0 ) then
				self%grid_type = trim( args(2) )
				end if
                !
				if( index( line_text, "Grid Reader" ) > 0 ) then
				self%grid_reader = trim( args(2) )
				end if
                !
            end if
         end do
         !
10       close( unit = funit )
         !
         write(*,*) 'Finish read file [', fname, ']:'
         if ( allocated( self%grid_type ) ) write(*,*) '     Set grid_type = ', self%grid_type
         if ( allocated( self%grid_reader ) ) write(*,*) '     Set grid_reader = ', self%grid_reader
         !
      end if
      !
   end function ModEMControlFile_ctor
   !
   subroutine ModEMControlFile_dtor( self )
      implicit none
      !
      type( ModEMControlFile_t ), intent( in out ) :: self
      !
      ! write(*,*) "Destructor ModEMControlFile_t"
      !
   end subroutine ModEMControlFile_dtor
   !
   ! Return the number of lines of a given file
   function getLineNumber( funit ) result( line_counter )
      implicit none
      !
      integer, intent( in ) :: funit
      !
      integer       :: line_counter
      character(80) :: line_text
      !
      line_counter = 0
      !
      rewind( funit )
      !
      do
         read( funit, "(a)", END = 10 ) line_text
         line_text = adjustl( line_text )
         if( index( line_text, "#" ).eq.0 ) then
            !
            line_counter = line_counter + 1
         end if
      end do
      !
10      return
      !
   end function getLineNumber
   !
end module ModEMControlFile
