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
   use Constants
   use String
   !
   type :: ModEMControlFile_t
      !
      character(:), allocatable :: grid_reader
      character(:), allocatable :: grid_type
      character(:), allocatable :: forward_solver
      character(:), allocatable :: source
      character(:), allocatable :: model_method, model_n_air_layer, model_max_height
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
      type( ModEMControlFile_t ) :: self
      !
      character(1000)                   :: full_line_text
      character(len=200), dimension(20) :: args
      character(:), allocatable         :: line_text
      integer                           :: line_counter, io_stat, p_nargs
      !
      !write( *,* ) "Constructor ModEMControlFile_t"
      !
      call Compact( fname )
      !
      open( unit = funit, file = fname, iostat = io_stat, status = "old" )
      !
      if( io_stat /= 0 ) then
         write( *,* ) "Unable to open [", fname, "], Stat: ", io_stat
      else
         line_counter = getLineNumber( funit )
         !
         rewind( funit )
         !
         do
            read( funit, "(a)", END = 10 ) full_line_text
            line_text = adjustl( full_line_text )
            line_text = trim( line_text )
            !
            call Parse( line_text, ":", args, p_nargs )
            !
            if( index( line_text, "#" ) == 0 .and. index( line_text, ">" ) == 0 ) then
                !
                if( index( line_text, "grid_header" ) > 0 ) then
                   self%grid_type = trim( args(2) )
                end if
                !
                if( index( line_text, "grid_type" ) > 0 ) then
                   self%grid_reader = trim( args(2) )
                end if
                !
                if( index( line_text, "forward_solver" ) > 0 ) then
                   self%forward_solver = trim( args(2) )
                end if
                !
                if( index( line_text, "source" ) > 0 ) then
                   self%source = trim( args(2) )
                end if
                !
                if( index( line_text, "model_method" ) > 0 ) then
                   self%model_method = trim( args(2) )
                end if
                !
                if( index( line_text, "model_n_air_layer" ) > 0 ) then
                   self%model_n_air_layer = trim( args(2) )
                end if
                !
                if( index( line_text, "model_max_height" ) > 0 ) then
                   self%model_max_height = trim( args(2) )
                end if
                !
            end if
         end do
         !
10       close( unit = funit )
         !
         ! GRID TYPE
         !
         if ( allocated( self%grid_type ) ) then
            write( *, "(A20, A10)" ) "      grid_type =", self%grid_type
            !
            select case ( self%grid_type )
               case( "SG" )
                  grid_type = GRID_SG
               case( "MR" )
                  grid_type = GRID_MR
               case default
                  grid_type = ""
                  STOP "Wrong grid_type control, use [SG|MR]"
            end select
            !
         endif
         !
         ! GRID READER
         !
         if ( allocated( self%grid_reader ) ) then
            write( *, "(A20, A10)" ) "      grid_reader =", self%grid_reader
         endif
         !
         ! FOWARD SOLVER
         !
         if ( allocated( self%forward_solver ) ) then
            write( *, "(A20, A10)" ) "      fwd_solver =", self%forward_solver
            !
            select case ( self%forward_solver )
               case( "FILE" )
                  forward_solver_type = FWD_FILE
                case( "IT" )
                  forward_solver_type = FWD_IT
               case( "IT_DC" )
                  forward_solver_type = FWD_IT_DC
               case default
                  forward_solver_type = ""
                  STOP "Wrong forward_solver control, use [FILE|IT|IT_DC]"
            end select
            !
         endif
         !
         ! SOURCE
         !
         if ( allocated( self%source ) ) then
            write( *, "(A20, A10)" ) "      source =", self%source
            !
            select case ( self%source )
               case( "1D" )
                  source_type = SRC_MT_1D
               case( "2D" )
                  source_type = SRC_MT_2D
               case default
                  source_type = ""
                  STOP "Wrong source control, use [1D|2D]"
            end select
            !
         endif
         !
         ! MODEL METHOD
         !
         if ( allocated( self%model_method ) ) then
            write( *, "(A20, A10)" ) "      model_method =", self%model_method
            !
            select case ( self%model_method )
               case( "fixed height" )
                  model_method = MM_METHOD_FIXED_H
               case default
                  model_method = ""
                  STOP "Wrong model_method control, use [none|fixed height]"
            end select
            !
         endif
         !
         ! MODEL NZAIR
         !
         if ( allocated( self%model_n_air_layer ) ) then
            write( *, "(A20, A30)" ) "      model_n_air_layer =", self%model_n_air_layer
            !
			read( self%model_n_air_layer, '(i5)' )   model_n_air_layer
            !
         endif
         !
         ! MODEL NZAIR
         !
         if ( allocated( self%model_max_height ) ) then
            write( *, "(A20, A30)" ) "      model_max_height =", self%model_max_height
            !
			read( self%model_max_height, '(f15.5)' )   model_max_height
            !
         endif
         !
      end if
      !
   end function ModEMControlFile_ctor
   !
   subroutine ModEMControlFile_dtor( self )
      implicit none
      !
      type( ModEMControlFile_t ), intent( inout ) :: self
      !
      !write( *,* ) "Destructor ModEMControlFile_t"
      !
      if( allocated( self%grid_reader ) ) deallocate( self%grid_reader )
      if( allocated( self%grid_type ) ) deallocate( self%grid_type )
      if( allocated( self%forward_solver ) ) deallocate( self%forward_solver )
      if( allocated( self%source ) ) deallocate( self%source )
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
