!
!> Class to read a control file
!> And set Forward Modeling parameters
!
module ForwardControlFile
    !
    use Utilities
    use String
    use Grid3D_MR
    use ForwardSolver
    use Solver
    !
    integer, allocatable, dimension(:) :: grid_layers
    !
    type :: ForwardControlFile_t
        !
        !> FWD Components parameters
        character(:), allocatable :: model_operator_type
        !
        character(:), allocatable :: grid_reader_type, grid_format, solver_type, forward_solver_type
        character(:), allocatable :: model_method, model_n_air_layer, model_max_height
        !
        !> Solver parameters
        character(:), allocatable :: max_solver_iters, max_solver_calls, max_divcor_iters
        character(:), allocatable :: tolerance_divcor, tolerance_solver
        !
        contains
            !
            final :: ForwardControlFile_dtor
            !
    end type ForwardControlFile_t
    !
    !> Public Global ForwardControlFile object
    type( ForwardControlFile_t ), allocatable :: fwd_control_file
    !
    interface ForwardControlFile_t
        module procedure ForwardControlFile_ctor
    end interface ForwardControlFile_t
    !
contains
    !
    !> Procedure ForwardControlFile_ctor
    !> Read line by line of the data file, create Data Entry objects(MT, MT_REF or CSEM)
    !
    function ForwardControlFile_ctor( funit, fname ) result( self )
        implicit none
        !
        integer, intent( in ) :: funit
        character(:), allocatable, intent( in ) :: fname
        !
        type( ForwardControlFile_t ) :: self
        !
        character(1000) :: full_line_text
        character( len=200 ), dimension(20) :: args
        character(:), allocatable :: line_text
        integer :: i, line_counter, io_stat, p_nargs
        !
        !write( *,* ) "Constructor ForwardControlFile_t"
        !
        call compact( fname )
        !
        open( unit = funit, file = fname, iostat = io_stat, status = "old" )
        !
        if( io_stat /= 0 ) then
            !
            call errStop( "ForwardControlFile_ctor > Unable to open file ["//fname//"]" )
            !
        else
            !
            do
                !
                read( funit, "(a)", END = 10 ) full_line_text
                line_text = adjustl( full_line_text )
                line_text = trim( line_text )
                !
                call Parse( line_text, ":", args, p_nargs )
                !
                if( index( line_text, "#" ) == 0 .AND. index( line_text, ">" ) == 0 ) then
                    !
                    if( index( line_text, "grid_reader" ) > 0 ) then
                        self%grid_reader_type = trim( args(2) )
                    elseif( index( line_text, "grid_format" ) > 0 ) then
                        self%grid_format = trim( args(2) )
                    elseif( index( line_text, "model_operator_type" ) > 0 ) then
                        self%model_operator_type = trim( args(2) )
                    elseif( index( line_text, "forward_solver_type" ) > 0 ) then
                        self%forward_solver_type = trim( args(2) )
                    elseif( index( line_text, "solver_type" ) > 0 ) then
                        self%solver_type = trim( args(2) )
                    elseif( index( line_text, "model_method" ) > 0 ) then
                        self%model_method = trim( args(2) )
                    elseif( index( line_text, "model_n_air_layer" ) > 0 ) then
                        self%model_n_air_layer = trim( args(2) )
                    elseif( index( line_text, "model_max_height" ) > 0 ) then
                        self%model_max_height = trim( args(2) )
                    elseif( index( line_text, "max_solver_iters" ) > 0 ) then
                        self%max_solver_iters = trim( args(2) )
                    elseif( index( line_text, "max_solver_calls" ) > 0 ) then
                        self%max_solver_calls = trim( args(2) )
                    elseif( index( line_text, "max_divcor_iters" ) > 0 ) then
                        self%max_divcor_iters = trim( args(2) )
                    elseif( index( line_text, "tolerance_divcor" ) > 0 ) then
                        self%tolerance_divcor = trim( args(2) )
                    elseif( index( line_text, "tolerance_solver" ) > 0 ) then
                        self%tolerance_solver = trim( args(2) )
                    else
                        call errStop( "Unsupported Forward Modeling parameter: ["//trim(line_text)//"]" )
                    endif
                    !
                endif
                !
            enddo
            !
10          close( unit = funit )
            !
            ! Field type
            if( allocated( self%model_operator_type ) ) then
                !
                select case( self%model_operator_type )
                    case( "MF" )
                        model_operator_type = MODELOP_MF
                    case( "SP" )
                        model_operator_type = MODELOP_SP
                    case( "SP2" )
                        model_operator_type = MODELOP_SP2
                    case default
                        !
                        call errStop( "ForwardControlFile_ctor > Wrong model_operator_type, use [MF|SP|SP2]" )
                end select
                !
                write( *, "( A35, A20 )" ) "model_operator_type = ", model_operator_type
                !
            endif
            !
            ! Grid type
            if( allocated( self%grid_format ) ) then
                !
                line_text = trim( self%grid_format )
                !
                call Parse( line_text, ",", args, p_nargs )
                !
                allocate( grid_layers( p_nargs ) )
                !
                do i = 1, p_nargs
                    !
                    read( args(i), "( I8 )" ) grid_layers(i)
                    !
                enddo
                !
                grid_format = GRID_MR
                !
            else
                grid_format = GRID_SG
            endif
            !
            write( *, "( A35, A20 )" ) "Grid Format = ", grid_format
            !
            ! Grid reader
            if( allocated( self%grid_reader_type ) ) then
                !
                ! TO BE IMPLEMENTED
                write( *, "( A35, A20 )" ) "Grid Reader = ", self%grid_reader_type
                !
            endif
            !
            ! Forward solver
            if( allocated( self%solver_type ) ) then
                !
                select case( self%solver_type )
                    !
                    case( "QMR" )
                        solver_type = SLV_QMR
                    case( "BICG" )
                        solver_type = SLV_BICG
                    case default
                        !
                        call errStop( "ForwardControlFile_ctor > Wrong solver_type, use [QMR|BICG]" )
                    !
                end select
                !
            endif
            !
            ! Forward solver
            if( allocated( self%forward_solver_type ) ) then
                !
                select case( self%forward_solver_type )
                    !
                    case( "IT" )
                        forward_solver_type = FWD_IT
                    case( "IT_DC" )
                        forward_solver_type = FWD_IT_DC
                    case default
                        !
                        call errStop( "ForwardControlFile_ctor > Wrong forward_solver_type, use [IT|IT_DC]" )
                    !
                end select
                !
            endif
            !
            ! Model method
            if( allocated( self%model_method ) ) then
                !
                select case( self%model_method )
                    !
                    case( "fixed height" )
                        model_method = MM_METHOD_FIXED_H
                    case( "mirror" )
                        model_method = MM_METHOD_MIRROR
                    case default
                        !
                        call errStop( "ForwardControlFile_ctor > Wrong model_method, use [mirror|fixed height]" )
                    !
                end select
                !
            endif
            !
            ! Model nzAir
            if( allocated( self%model_n_air_layer ) ) then
                !
                read( self%model_n_air_layer, "(I8)" ) model_n_air_layer
                !
            endif
            !
            ! Model max height
            if( allocated( self%model_max_height ) ) then
                !
                read( self%model_max_height, "(f15.6)" ) model_max_height
                !
            endif
            !
            ! Solver max_solver_iters
            if( allocated( self%max_solver_iters ) ) then
                !
                read( self%max_solver_iters, "(I8)" ) max_solver_iters
                !
                write( *, "( A26, A9, I20)" ) solver_type, "Iters = ", max_solver_iters
                !
            endif
            !
            ! Solver max_solver_calls
            if( allocated( self%max_solver_calls ) ) then
                !
                read( self%max_solver_calls, "(I8)" ) max_solver_calls
                !
                write( *, "( A22, A4, A9, I20)" ) "Max ", solver_type, " Calls = ", max_solver_calls
                !
            endif
            !
            ! Solver tolerance_solver
            if( allocated( self%tolerance_solver ) ) then
                !
                read( self%tolerance_solver, * ) tolerance_solver
                !
                write( *, "( A22, A13, es20.2)" ) solver_type, " Tolerance = ", tolerance_solver
                !
            endif
            !
            ! Solver max_divcor_iters
            if( allocated( self%max_divcor_iters ) ) then
                !
                read( self%max_divcor_iters, "(I8)" ) max_divcor_iters
                !
                write( *, "( A35, I20)" ) "DivCorr Iters = ", max_divcor_iters
                !
            endif
            !
            ! Solver tolerance_divcor
            if( allocated( self%tolerance_divcor ) ) then
                !
                read( self%tolerance_divcor, * ) tolerance_divcor
                !
                write( *, "( A35, es20.2)" ) "DivCorr Tolerance = ", tolerance_divcor
                !
            endif
            !
            !> Performance tips for the users
            select case( forward_solver_type )
                !
                case( FWD_IT )
                    !
                    if( model_operator_type .EQ. MODELOP_MF .OR. model_operator_type .EQ. MODELOP_SP ) then
                        call warning( "Better to use forward_solver_type: IT_DC!" )
                    endif
                    !
                case( FWD_IT_DC )
                    !
                    if( model_operator_type .EQ. MODELOP_SP2 ) then
                        call warning( "Better to use forward_solver_type: IT!" )
                    endif
                    !
            end select
            !
            write( *, * ) ""
            !
        endif
        !
    end function ForwardControlFile_ctor
    !
    !> Deconstructor routine:
    !>     Deallocates inherent properties of this class.
    subroutine ForwardControlFile_dtor( self )
        implicit none
        !
        type( ForwardControlFile_t ), intent( inout ) :: self
        !
        !write( *,* ) "Destructor ForwardControlFile_t"
        !
        if( allocated( self%model_operator_type ) ) deallocate( self%model_operator_type )
        !
        if( allocated( self%grid_reader_type ) ) deallocate( self%grid_reader_type )
        if( allocated( self%grid_format ) ) deallocate( self%grid_format )
        !
        if( allocated( grid_layers ) ) deallocate( grid_layers )
        !
        if( allocated( self%forward_solver_type ) ) deallocate( self%forward_solver_type )
        !
        if( allocated( self%model_method ) ) deallocate( self%model_method )
        if( allocated( self%model_n_air_layer ) ) deallocate( self%model_n_air_layer )
        if( allocated( self%model_max_height ) ) deallocate( self%model_max_height )
        !
        if( allocated( self%max_solver_iters ) ) deallocate( self%max_solver_iters )
        !
        if( allocated( self%max_solver_calls ) ) deallocate( self%max_solver_calls )
        if( allocated( self%max_divcor_iters ) ) deallocate( self%max_divcor_iters )
        if( allocated( self%tolerance_divcor ) ) deallocate( self%tolerance_divcor )
        if( allocated( self%tolerance_solver ) ) deallocate( self%tolerance_solver )
        !
    end subroutine ForwardControlFile_dtor
    !
end module ForwardControlFile
!