!
!> Abstract Base class to define MetricElements
!
module MetricElements
    !
    use iScalar3D_SG
    use cScalar3D_SG
    use rScalar3D_MR
    use cVector3D_SG
    use rVector3D_MR
    use cVector3D_MR
    !
    type, abstract :: MetricElements_t
        !
        class( Grid_t ), pointer :: grid
        !
        class( Vector_t ), allocatable :: edge_length, dual_edge_length
        !
        class( Vector_t ), allocatable :: face_area, dual_face_area
        !
        class( Vector_t ), allocatable :: v_edge
        !
        class( Scalar_t ), allocatable :: v_node, v_cell
        !
     contains
        !
        procedure( interface_set_edge_length_metric_elements ), deferred, public :: setEdgeLength
        procedure( interface_set_face_area_metric_elements ), deferred, public :: setFaceArea
        procedure( interface_set_dual_edge_length_metric_elements ), deferred, public :: setDualEdgeLength
        procedure( interface_set_dual_face_area_metric_elements ), deferred, public :: setDualFaceArea
        procedure( interface_set_cell_volume_metric_elements ), deferred, public :: setCellVolume
        procedure( interface_set_edge_volume_metric_elements ), deferred, public :: setEdgeVolume
        procedure( interface_set_node_volume_metric_elements ), deferred, public :: setNodeVolume
        !
        procedure( interface_set_index_arrays_metric_elements ), deferred, public :: setIndexArrays
        !
        procedure( interface_set_all_index_arrays_metric_elements ), deferred, public :: setAllIndexArrays
        !
        !procedure( interface_boundary_index_metric_elements ), deferred, public :: boundaryIndex
        !
        procedure, public :: baseDealloc => deallocate_MetricElements
        !
        procedure, public :: alloc => allocate_MetricElements
        !
        procedure, public :: setup => setup_MetricElements
        !
        procedure, public :: createScalar, createVector
    !
    end type MetricElements_t
    !
    abstract interface
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_edge_length_metric_elements( self )
            import :: MetricElements_t
            !
            class( MetricElements_t ), intent( inout ) :: self
            !
        end subroutine interface_set_edge_length_metric_elements
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_face_area_metric_elements( self )
            import :: MetricElements_t
            !
            class( MetricElements_t ), intent( inout ) :: self
            !
        end subroutine interface_set_face_area_metric_elements
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_dual_edge_length_metric_elements( self )
            import :: MetricElements_t
            !
            class( MetricElements_t ), intent( inout ) :: self
            !
        end subroutine interface_set_dual_edge_length_metric_elements
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_dual_face_area_metric_elements( self )
            import :: MetricElements_t
            !
            class( MetricElements_t ), intent( inout ) :: self
            !
        end subroutine interface_set_dual_face_area_metric_elements
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_cell_volume_metric_elements( self )
            import :: MetricElements_t
            !
            class( MetricElements_t ), intent( inout ) :: self
            !
        end subroutine interface_set_cell_volume_metric_elements
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_edge_volume_metric_elements( self )
            import :: MetricElements_t
            !
            class( MetricElements_t ), intent( inout ) :: self
            !
        end subroutine interface_set_edge_volume_metric_elements
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_node_volume_metric_elements( self )
            import :: MetricElements_t
            !
            class( MetricElements_t ), intent( inout ) :: self
            !
        end subroutine interface_set_node_volume_metric_elements
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_index_arrays_metric_elements( self, grid_type, INDb, INDi, INDa )
            import :: MetricElements_t
            !
            class( MetricElements_t ), intent( in ) :: self
            character(*), intent( in ) :: grid_type
            integer, allocatable, dimension(:), intent( out ) :: INDb, INDi
            integer, dimension(:), allocatable, intent( out ), optional :: INDa
            !
        end subroutine interface_set_index_arrays_metric_elements
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_all_index_arrays_metric_elements( self )
            import :: MetricElements_t
            !
            class( MetricElements_t ), intent( in ) :: self
            !
        end subroutine interface_set_all_index_arrays_metric_elements
        ! !
        ! subroutine interface_boundary_index_metric_elements( self, grid_type, INDb, INDi )
            ! import :: MetricElements_t
            ! !
            ! class( MetricElements_t ), intent( in ) :: self
            ! character(*), intent( in ) :: grid_type
            ! integer, allocatable, dimension(:), intent( inout ) :: INDb, INDi
            ! !
        ! end subroutine interface_boundary_index_metric_elements
        ! !
    end interface
    !
contains
    !
    !> No subroutine briefing
    !
    subroutine setup_MetricElements( self )
        implicit none
        !
        class( MetricElements_t ), intent( inout ) :: self
        !
        call self%setEdgeLength
        call self%setFaceArea
        call self%setDualEdgeLength
        call self%setDualFaceArea
        call self%setEdgeVolume
        !
        call self%setCellVolume
        call self%setNodeVolume
        !
    end subroutine setup_MetricElements
    !
    !> No subroutine briefing
    !
    subroutine allocate_MetricElements( self )
        implicit none
        !
        class( MetricElements_t ), intent( inout ) :: self
        !
        call self%createVector( real_t, EDGE, self%edge_length )
        call self%createVector( real_t, FACE, self%dual_edge_length )
        call self%createVector( real_t, FACE, self%face_area )
        call self%createVector( real_t, EDGE, self%dual_face_area )
        !
        call self%createVector( real_t, EDGE, self%v_edge )
        !
        call self%createScalar( real_t, NODE, self%v_node )
        call self%createScalar( real_t, CELL, self%v_cell )
        !
    end subroutine allocate_MetricElements
    !
    !> No subroutine briefing
    !
    subroutine deallocate_MetricElements( self )
        implicit none
        !
        class( MetricElements_t ), intent( inout ) :: self
        !
        if( allocated( self%edge_length ) ) deallocate( self%edge_length )
        if( allocated( self%face_area ) ) deallocate( self%face_area )
        if( allocated( self%dual_face_area ) ) deallocate( self%dual_face_area )
        if( allocated( self%dual_edge_length ) ) deallocate( self%dual_edge_length )
        if( allocated( self%v_edge ) ) deallocate( self%v_edge )
        !
        if( allocated( self%v_node ) ) deallocate( self%v_node )
        if( allocated( self%v_cell ) ) deallocate( self%v_cell )
        !
    end subroutine deallocate_MetricElements
    !
    !> Create proper scalar from the Grid
    !
    subroutine createScalar( self, scalar_type, grid_type, scalar )
        implicit none
        !
        class( MetricElements_t ), intent( in ) :: self
        integer, intent( in ) :: scalar_type
        character( len=4 ), intent( in ) :: grid_type
        class( Scalar_t ), allocatable, intent( out ) :: scalar
        !
        if( grid_type /= NODE .AND. grid_type /= CELL .AND. grid_type /= CELL_EARTH ) then
            call errStop( "createScalar > grid_type must be NODE, CELL or CELL_EARTH" )
        else
            !
            select type( grid => self%grid )
                !
                class is( Grid3D_SG_t )
                    !
                    if( scalar_type == real_t ) then
                        allocate( scalar, source = rScalar3D_SG_t( grid, grid_type ) )
                    elseif( scalar_type == complex_t ) then
                        allocate( scalar, source = cScalar3D_SG_t( grid, grid_type ) )
                    elseif( scalar_type == integer_t ) then
                        allocate( scalar, source = iScalar3D_SG_t( grid, grid_type ) )
                    else
                        call errStop( "createScalar > choose SG real_t, complex_t or integer_t" )
                    endif
                    !
                class is( Grid3D_MR_t )
                    !
                    if( scalar_type == real_t ) then
                        allocate( scalar, source = rScalar3D_MR_t( grid, grid_type ) )
                    elseif( scalar_type == complex_t ) then
                        !allocate( scalar, source = cScalar3D_MR_t( grid, grid_type ) )
                        !
                        call errStop( "createScalar > MR complex_t to be implemented" )
                    elseif( scalar_type == integer_t ) then
                        !allocate( scalar, source = iScalar3D_MR_t( grid, grid_type ) )
                        !
                        call errStop( "createScalar > MR integer_t to be implemented" )
                    else
                        call errStop( "createScalar > choose MR real_t, complex_t or integer_t" )
                    endif
                    !
                class default
                   call errStop( "createScalar > Unclassified grid" )
                !
            end select
            !
        endif
        !
    end subroutine createScalar
    !
    !> Create proper vector from the Grid
    !
    subroutine createVector( self, vector_type, grid_type, vector )
        implicit none
        !
        class( MetricElements_t ), intent( in ) :: self
        integer, intent( in ) :: vector_type
        character( len=4 ), intent( in ) :: grid_type
        class( Vector_t ), allocatable, intent( out ) :: vector
        !
        if( grid_type /= EDGE .AND. grid_type /= FACE ) then
            call errStop( "createVector > grid_type must be EDGE or FACE" )
        else
            !
            select type( grid => self%grid )
                !
                class is( Grid3D_SG_t )
                    !
                    if( vector_type == real_t ) then
                        allocate( vector, source = rVector3D_SG_t( grid, grid_type ) )
                    elseif( vector_type == complex_t ) then
                        allocate( vector, source = cVector3D_SG_t( grid, grid_type ) )
                    elseif( vector_type == integer_t ) then
                        !allocate( vector, source = iVector3D_SG_t( grid, grid_type ) )
                        !
                        call errStop( "createVector > SG integer_t to be implemented" )
                    else
                        call errStop( "createVector > choose SG: real_t, complex_t or integer_t" )
                    endif
                    !
                class is( Grid3D_MR_t )
                    !
                    if( vector_type == real_t ) then
                        !
                        allocate( vector, source = rVector3D_MR_t( grid, grid_type ) )
                        !
                    elseif( vector_type == complex_t ) then
                        !
                        allocate( vector, source = cVector3D_MR_t( grid, grid_type ) )
                        !
                    elseif( vector_type == integer_t ) then
                        !allocate( vector, source = iVector3D_MR_t( grid, grid_type ) )
                        !
                        call errStop( "createVector > MR integer_t to be implemented" )
                    else
                        call errStop( "createVector > choose MR: real_t, complex_t or integer_t" )
                    endif
                    !
                class default
                   call errStop( "createVector > Unclassified grid" )
                !
            end select
            !
        endif
        !
    end subroutine createVector
    !
end module MetricElements
