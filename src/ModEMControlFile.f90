!*************
!
! Class to read a data file and create an array with all data entries(lines)
!
!*************
!
module ModEMControlFile
    !
    use Constants
    use String
    use Grid
    use ForwardSolver
    use Solver
    use Source
    !
    type :: ModEMControlFile_t
        !
        character(:), allocatable :: grid_reader_type, grid_type, forward_solver_type, source_type
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
            if( io_stat == 0 ) then
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
                        if( index( line_text, "grid_reader" ) > 0 ) then
                            self%grid_type = trim( args(2) )
                        endif
                        !
                        if( index( line_text, "grid" ) > 0 ) then
                            self%grid_reader_type = trim( args(2) )
                        endif
                        !
                        if( index( line_text, "forward_solver" ) > 0 ) then
                            self%forward_solver_type = trim( args(2) )
                        endif
                        !
                        if( index( line_text, "source" ) > 0 ) then
                            self%source_type = trim( args(2) )
                        endif
                        !
                        if( index( line_text, "model_method" ) > 0 ) then
                            self%model_method = trim( args(2) )
                        endif
                        !
                        if( index( line_text, "model_n_air_layer" ) > 0 ) then
                            self%model_n_air_layer = trim( args(2) )
                        endif
                        !
                        if( index( line_text, "model_max_height" ) > 0 ) then
                            self%model_max_height = trim( args(2) )
                        endif
                        !
                    endif
                    !
                end do
                !
10              close( unit = funit )
                !
                ! Grid type
                if ( allocated( self%grid_type ) ) then
                    !
                    write( *, * ) "          grid = ", self%grid_type
                    !
                    select case ( self%grid_type )
                        case( "SG" )
                            grid_type = GRID_SG
                        case( "MR" )
                            grid_type = GRID_MR
                        case default
                            grid_type = ""
                            stop "Error: Wrong grid_type control, use [SG|MR]"
                    end select
                    !
                endif
                !
                ! Grid reader
                if ( allocated( self%grid_reader_type ) ) then
                    write( *, * ) "          grid_reader = ", self%grid_reader_type
                endif
                !
                ! Forward solver
                if ( allocated( self%forward_solver_type ) ) then
                    !
                    write( *, * ) "          fwd_solver = ", self%forward_solver_type
                    !
                    select case ( self%forward_solver_type )
                        !
                        case( "FILE" )
                            forward_solver_type = FWD_FILE
                        case( "IT" )
                            forward_solver_type = FWD_IT
                        case( "IT_DC" )
                            forward_solver_type = FWD_IT_DC
                        case default
                            forward_solver_type = ""
                            stop "Error: Wrong forward_solver control, use [FILE|IT|IT_DC]"
                        !
                    end select
                    !
                endif
                !
                ! Source_type
                if ( allocated( self%source_type ) ) then
                    write( *, * ) "          source = ", self%source_type
                    !
                    select case ( self%source_type )
                        !
                        case( "1D" )
                            source_type = SRC_MT_1D
                        case( "2D" )
                            source_type = SRC_MT_2D
                        case default
                            source_type = ""
                        stop "Error: Wrong source control, use [1D|2D]"
                        !
                    end select
                    !
                endif
                !
                ! Model method
                if ( allocated( self%model_method ) ) then
                    write( *, * ) "          model_method = ", self%model_method
                    !
                    select case ( self%model_method )
                        !
                        case( "fixed height" )
                            model_method = MM_METHOD_FIXED_H
                        case( "mirror" )
                            model_method = MM_METHOD_MIRROR
                        case default
                            model_method = ""
                        stop "Error: Wrong model_method control, use [mirror|fixed height]"
                        !
                    end select
                    !
                endif
                !
                ! Model nzAir
                if ( allocated( self%model_n_air_layer ) ) then
                    !
                    write( *, * ) "          model_n_air_layer =", self%model_n_air_layer
                    !
                    read( self%model_n_air_layer, "(i5)" ) model_n_air_layer
                    !
                endif
                !
                ! Model max height
                if ( allocated( self%model_max_height ) ) then
                    !
                    write( *, * ) "          model_max_height =", self%model_max_height
                    !
                    read( self%model_max_height, "(f15.5)" ) model_max_height
                    !
                endif
                    !
            else
                write( *, * ) "Error opening [", fname, "] in ModEMControlFile_ctor"
                stop
            endif
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
            if( allocated( self%grid_reader_type ) ) deallocate( self%grid_reader_type )
            if( allocated( self%grid_type ) ) deallocate( self%grid_type )
            if( allocated( self%forward_solver_type ) ) deallocate( self%forward_solver_type )
            if( allocated( self%source_type ) ) deallocate( self%source_type )
            !
        end subroutine ModEMControlFile_dtor
        !
end module ModEMControlFile
