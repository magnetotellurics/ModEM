!
!> Class to read a control file
!> And set Forward Modeling parameters
!
module ForwardControlFile
    !
    use Constants
    use String
    use Grid
    use ForwardSolver
    use Solver
    use Source
    !
    type :: ForwardControlFile_t
        !
        !> FWD Components parameters
        character(:), allocatable :: model_operator_type
        !
        character(:), allocatable :: grid_reader_type, grid_type, forward_solver_type
        character(:), allocatable :: source_type_mt, source_type_csem, get_1d_from
        character(:), allocatable :: model_method, model_n_air_layer, model_max_height
        !
        !> Solver parameters
        character(:), allocatable :: max_solver_iters, max_divcor_calls, max_divcor_iters
        character(:), allocatable :: tolerance_divcor, tolerance_solver
        !
        contains
            !
            final :: ForwardControlFile_dtor
            !
    end type ForwardControlFile_t
    !
    interface ForwardControlFile_t
        module procedure ForwardControlFile_ctor
    end interface ForwardControlFile_t
!
contains
    !
    !> Procedure ForwardControlFile_ctor
    !> Read line by line of the data file, create Data Entry objects(MT, MT_REF or CSEM)
    function ForwardControlFile_ctor( funit, fname ) result( self )
        implicit none
        !
        integer, intent( in ) :: funit
        character(:), allocatable, intent( in ) :: fname
        !
        type( ForwardControlFile_t ) :: self
        !
        character(1000) :: full_line_text
        character(len=200), dimension(20) :: args
        character(:), allocatable :: line_text
        integer :: line_counter, io_stat, p_nargs
        !
        !write( *,* ) "Constructor ForwardControlFile_t"
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
                if( index( line_text, "#" ) == 0 .AND. index( line_text, ">" ) == 0 ) then
                    !
                    if( index( line_text, "grid_reader" ) > 0 ) then
                        self%grid_reader_type = trim( args(2) )
                    elseif( index( line_text, "grid_type" ) > 0 ) then
                        self%grid_type = trim( args(2) )
                    elseif( index( line_text, "model_operator_type" ) > 0 ) then
                        self%model_operator_type = trim( args(2) )
                    elseif( index( line_text, "forward_solver_type" ) > 0 ) then
                        self%forward_solver_type = trim( args(2) )
                    elseif( index( line_text, "source_type_mt" ) > 0 ) then
                        self%source_type_mt = trim( args(2) )
                    elseif( index( line_text, "source_type_csem" ) > 0 ) then
                        self%source_type_csem = trim( args(2) )
                    elseif( index( line_text, "get_1d_from" ) > 0 ) then
                        self%get_1d_from = trim( args(2) )
                    elseif( index( line_text, "model_method" ) > 0 ) then
                        self%model_method = trim( args(2) )
                    elseif( index( line_text, "model_n_air_layer" ) > 0 ) then
                        self%model_n_air_layer = trim( args(2) )
                    elseif( index( line_text, "model_max_height" ) > 0 ) then
                        self%model_max_height = trim( args(2) )
                    elseif( index( line_text, "max_solver_iters" ) > 0 ) then
                        self%max_solver_iters = trim( args(2) )
                    elseif( index( line_text, "max_divcor_calls" ) > 0 ) then
                        self%max_divcor_calls = trim( args(2) )
                    elseif( index( line_text, "max_divcor_iters" ) > 0 ) then
                        self%max_divcor_iters = trim( args(2) )
                    elseif( index( line_text, "tolerance_divcor" ) > 0 ) then
                        self%tolerance_divcor = trim( args(2) )
                    elseif( index( line_text, "tolerance_solver" ) > 0 ) then
                        self%tolerance_solver = trim( args(2) )
                    else
                        write( *, * ) "     "//achar(27)//"[31m# Error:"//achar(27)//"[0m Unsupported Forward Modeling parameter: ["//trim(line_text)//"]"
                        stop 
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
                        model_operator_type = ""
                        stop "Error: Wrong model_operator_type control, use [MF|SP|SP2]"
                end select
                !
                write( *, "( A30, A20)" ) "          model_operator_type = ", model_operator_type
                !
            endif
            !
            ! Grid type
            if( allocated( self%grid_type ) ) then
                !
                select case( self%grid_type )
                    case( "SG" )
                        grid_type = GRID_SG
                    case( "MR" )
                        grid_type = GRID_MR
                    case default
                        grid_type = ""
                        stop "Error: Wrong grid_type control, use [SG|MR]"
                end select
                !
                write( *, "( A30, A20)" ) "          Grid Type = ", grid_type
                !
            endif
            !
            ! Grid reader
            if( allocated( self%grid_reader_type ) ) then
                !
                ! TO BE IMPLEMENTED
                write( *, "( A30, A20)" ) "          Grid Reader = ", self%grid_reader_type
                !
            endif
            !
            ! Forward solver
            if( allocated( self%forward_solver_type ) ) then
                !
                select case( self%forward_solver_type )
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
                write( *, "( A30, A20)" ) "          FWD Solver = ", forward_solver_type
                !
            endif
            !
            ! MT Source_type
            if( allocated( self%source_type_mt ) ) then
                !
                select case( self%source_type_mt )
                    !
                    case( "1D" )
                        source_type_mt = SRC_MT_1D
                    case( "2D" )
                        source_type_mt = SRC_MT_2D
                    case default
                        source_type_mt = ""
                        stop "Error: Wrong MT Source control, use [1D|2D]"
                        !
                end select
                !
                write( *, "( A30, A20)" ) "          MT Source = ", source_type_mt
                !
            endif
            !
            ! CSEM Source_type
            if( allocated( self%source_type_csem ) ) then
                !
                select case( self%source_type_csem )
                    !
                    case( "EM1D" )
                        source_type_csem = SRC_CSEM_EM1D
                    case( "Dipole1D" )
                        source_type_csem = SRC_CSEM_DIPOLE1D
                    case default
                        source_type_csem = ""
                        stop "Error: Wrong CSEM Source control, use [EM1D|Dipole1D]"
                        !
                end select
                !
                write( *, "( A30, A20)" ) "          CSEM Source = ", source_type_csem
                !
            endif
            !
            ! CSEM Source_type
            if( allocated( self%get_1d_from ) ) then
                !
                select case( self%get_1d_from )
                    !
                    case( "Fixed" )
                        get_1d_from = FROM_FIXED_VALUE
                    case( "Geometric_mean" )
                        get_1d_from = FROM_GEO_MEAN
                    case( "Mean_around_Tx" )
                        get_1d_from = FROM_TX_GEO_MEAN
                    case( "Tx_Position" )
                        get_1d_from = FROM_TX_LOCATION
                    case default
                        get_1d_from = ""
                        write( *, * ) "Error: Wrong get_1d_from, use [Fixed|Geometric_mean|Mean_around_Tx|Tx_Position]"
                        !
                end select
                !
                write( *, "( A30, A20)" ) "          Get 1D from = ", get_1d_from
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
                        model_method = ""
                    stop "Error: Wrong model_method control, use [mirror|fixed height]"
                    !
                end select
                write( *, "( A30, A20)" ) "          Model Method = ", model_method
                !
            endif
            !
            ! Model nzAir
            if( allocated( self%model_n_air_layer ) ) then
                !
                read( self%model_n_air_layer, "(I8)" ) model_n_air_layer
                !
                write( *, "( A30, I20)" ) "          N Air Layers = ", model_n_air_layer
                !
            endif
            !
            ! Model max height
            if( allocated( self%model_max_height ) ) then
                !
                read( self%model_max_height, "(f15.6)" ) model_max_height
                !
                write( *, "( A30, f20.2)" ) "          Model Max Height = ", model_max_height
                !
            endif
            !
            ! Solver max_solver_iters
            if( allocated( self%max_solver_iters ) ) then
                !
                read( self%max_solver_iters, "(I8)" ) max_solver_iters
                !
                write( *, "( A30, I20)" ) "          PCG Iters = ", max_solver_iters
                !
            endif
            !
            ! Solver max_divcor_calls
            if( allocated( self%max_divcor_calls ) ) then
                !
                read( self%max_divcor_calls, "(I8)" ) max_divcor_calls
                !
                write( *, "( A30, I20)" ) "          Max Divcor Calls = ", max_divcor_calls
                !
            endif
            !
            ! Solver max_divcor_iters
            if( allocated( self%max_divcor_iters ) ) then
                !
                read( self%max_divcor_iters, "(I8)" ) max_divcor_iters
                !
                write( *, "( A30, I20)" ) "          Max Divcor Iters = ", max_divcor_iters
                !
            endif
            !
            ! Solver tolerance_divcor
            if( allocated( self%tolerance_divcor ) ) then
                !
                read( self%tolerance_divcor, * ) tolerance_divcor
                !
                write( *, "( A30, es20.2)" ) "          Divcor Tolerance = ", tolerance_divcor
                !
            endif
            !
            ! Solver tolerance_solver
            if( allocated( self%tolerance_solver ) ) then
                !
                read( self%tolerance_solver, * ) tolerance_solver
                !
                write( *, "( A30, es20.2)" ) "          Solver Tolerance = ", tolerance_solver
                !
            endif
            !
        else
            write( *, * ) "Error opening [", fname, "] in ForwardControlFile_ctor"
            stop
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
        if( allocated( self%grid_type ) ) deallocate( self%grid_type )
        !
        if( allocated( self%forward_solver_type ) ) deallocate( self%forward_solver_type )
        !
        if( allocated( self%source_type_mt ) ) deallocate( self%source_type_mt )
        if( allocated( self%source_type_csem ) ) deallocate( self%source_type_csem )
        !
        if( allocated( self%model_method ) ) deallocate( self%model_method )
        if( allocated( self%model_n_air_layer ) ) deallocate( self%model_n_air_layer )
        if( allocated( self%model_max_height ) ) deallocate( self%model_max_height )
        !
        if( allocated( self%max_solver_iters ) ) deallocate( self%max_solver_iters )
        !
        if( allocated( self%max_divcor_calls ) ) deallocate( self%max_divcor_calls )
        if( allocated( self%max_divcor_iters ) ) deallocate( self%max_divcor_iters )
        if( allocated( self%tolerance_divcor ) ) deallocate( self%tolerance_divcor )
        if( allocated( self%tolerance_solver ) ) deallocate( self%tolerance_solver )
        !
    end subroutine ForwardControlFile_dtor
    !
end module ForwardControlFile
