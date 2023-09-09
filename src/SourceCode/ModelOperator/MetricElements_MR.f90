!
!> Derived class to define a ModelOperator_MF_SG
!>
!> This computes and stores Metric Elements for finite volume calculations
!> Based on Matlab development, code is take from ModEM module GridCalc.
!> In this Cartesian staggered grid_mr (CMR) version, metric elements are stored a
!> as TVector objects -- can be used to generate an SP version, and to 
!> generalize to MR; doing the same starting from GridCalcS can be used to
!> Create spherical versions
!> NOTE: there are other grid_mr mapping routines in GridCalc that are NOT to
!>            be included here -- I think these might be better as methods in the
!>            Vector/Scalar classes.
!
!> Variables that will be defined in base class    ... could add more as
!> in GridCalc -- but let's see if these are really useful.
!
module MetricElements_MR
    !
    use MetricElements_SG
    !
    type, extends( MetricElements_t ) :: MetricElements_MR_t
        !
        type( Grid3D_MR_t ), pointer :: grid_mr
        !
     contains
        !
        final :: MetricElements_MR_dtor
        !
        procedure, public :: setEdgeLength => setEdgeLength_MetricElements_MR
        !
        procedure, public :: setDualEdgeLength => setDualEdgeLength_MetricElements_MR
        !
        procedure, public :: setFaceArea => setFaceArea_MetricElements_MR
        !
        procedure, public :: setDualFaceArea => setDualFaceArea_MetricElements_MR
        !
        procedure, public :: setEdgeVolume => setEdgeVolume_MetricElements_MR
        !
        procedure, public :: setNodeVolume => setNodeVolume_MetricElements_MR
        !
        procedure, public :: setCellVolume => setCellVolume_MetricElements_MR
        !
        procedure, public :: setGridIndexArrays => setGridIndexArrays_MetricElements_MR
        !
        !procedure, public :: boundaryIndex => boundaryIndex_MetricElements_MR
        !
        procedure, public :: createScalar => createScalar_MetricElements_MR
        procedure, public :: createVector => createVector_MetricElements_MR
        !
    end type MetricElements_MR_t
    !
    interface MetricElements_MR_t
        module procedure MetricElements_MR_ctor
    end interface MetricElements_MR_t
    !
contains
    !
    !> No subroutine briefing
    !
    function MetricElements_MR_ctor( grid_mr ) result( self )
        implicit none
        !
        type( Grid3D_MR_t ), target :: grid_mr
        !
        type( MetricElements_MR_t ) :: self
        !
        !write( *, * ) "Constructor MetricElements_MR_t"
        !
        self%grid_mr => grid_mr
        !
        call self%alloc
        !
        !>    if were going to allocate storage for all, just set all now!
        call self%setup
        !
        call self%setGridIndexArrays( self%grid_mr )
        !
        call self%grid_mr%setXYZ
        !
        self%grid => grid_mr
        !
    end function MetricElements_MR_Ctor
    !
    !> No subroutine briefing
    !
    subroutine MetricElements_MR_dtor( self )
        implicit none
        !
        type( MetricElements_MR_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor MetricElements_MR"
        !
        call self%baseDealloc
        !
    end subroutine MetricElements_MR_dtor
    !
    !> No subroutine briefing
    !
    subroutine setEdgeLength_MetricElements_MR( self )
        implicit none
        !
        class( MetricElements_MR_t ), intent( inout ) :: self
        !
        integer :: i
        type( MetricElements_SG_t ) :: metric_sg
        !
        select type( edge_length => self%edge_length )
            !
            class is( rVector3D_MR_t )
                !
                do i = 1, self%grid_mr%n_grids
                    !
                    metric_sg = MetricElements_SG_t( self%grid_mr%sub_grid(i) )
                    !
                    edge_length%sub_vector(i) = metric_sg%edge_length
                    !
                enddo
                !
            class default
                call errStop( "setEdgeLength_MetricElements_MR > Unclassified edge_length" )
            !
        end select
        !
    end subroutine setEdgeLength_MetricElements_MR
    !
    !> No subroutine briefing
    !
    subroutine setDualEdgeLength_MetricElements_MR( self )
        implicit none
        !
        class( MetricElements_MR_t ), intent( inout ) :: self
        !
        integer :: i
        type( MetricElements_SG_t ) :: metric_sg
        !
        select type( dual_edge_length => self%dual_edge_length )
            !
            class is( rVector3D_MR_t )
                !
                do i = 1, self%grid_mr%n_grids
                    !
                    metric_sg = MetricElements_SG_t( self%grid_mr%sub_grid(i) )
                    !
                    dual_edge_length%sub_vector(i) = metric_sg%dual_edge_length
                    !
                enddo
                !
            class default
                call errStop( "setDualEdgeLength_MetricElements_MR > Unclassified dual_edge_length" )
            !
        end select
        !
    end subroutine setDualEdgeLength_MetricElements_MR
    !
    !> No subroutine briefing
    !
    subroutine setFaceArea_MetricElements_MR( self )
        implicit none
        !
        class( MetricElements_MR_t ), intent( inout ) :: self
        !
        integer :: i
        type( MetricElements_SG_t ) :: metric_sg
        !
        select type( face_area => self%face_area )
            !
            class is( rVector3D_MR_t )
                !
                do i = 1, self%grid_mr%n_grids
                    !
                    metric_sg = MetricElements_SG_t( self%grid_mr%sub_grid(i) )
                    !
                    face_area%sub_vector(i) = metric_sg%face_area
                    !
                enddo
                !
            class default
                call errStop( "setFaceArea_MetricElements_MR > Unclassified face_area" )
            !
        end select
        !
    end subroutine setFaceArea_MetricElements_MR
    !
    !> No subroutine briefing
    !
    subroutine setDualFaceArea_MetricElements_MR( self )
        implicit none
        !
        class( MetricElements_MR_t ), intent( inout ) :: self
        !
        integer :: i
        type( MetricElements_SG_t ) :: metric_sg
        !
        select type( dual_face_area => self%dual_face_area )
            !
            class is( rVector3D_MR_t )
                !
                do i = 1, self%grid_mr%n_grids
                    !
                    metric_sg = MetricElements_SG_t( self%grid_mr%sub_grid(i) )
                    !
                    dual_face_area%sub_vector(i) = metric_sg%dual_face_area
                    !
                enddo
                !
            class default
                call errStop( "setDualFaceArea_MetricElements_MR > Unclassified v_cell" )
            !
        end select
        !
    end subroutine setDualFaceArea_MetricElements_MR
    !
    !> No subroutine briefing
    !
    subroutine setEdgeVolume_MetricElements_MR( self )
        implicit none
        !
        class( MetricElements_MR_t ), intent( inout ) :: self
        !
        integer :: i
        type( MetricElements_SG_t ) :: metric_sg
        !
        select type( v_edge => self%v_edge )
            !
            class is( rVector3D_MR_t )
                !
                do i = 1, self%grid_mr%n_grids
                    !
                    metric_sg = MetricElements_SG_t( self%grid_mr%sub_grid(i) )
                    !
                    v_edge%sub_vector(i) = metric_sg%v_edge
                    !
                enddo
                !
            class default
                call errStop( "setEdgeVolume_MetricElements_MR > Unclassified grid_mr" )
            !
        end select
        !
    end subroutine setEdgeVolume_MetricElements_MR
    !
    !> No subroutine briefing
    !
    subroutine setNodeVolume_MetricElements_MR( self )
        implicit none
        !
        class( MetricElements_MR_t ), intent( inout ) :: self
        !
        integer :: i
        type( MetricElements_SG_t ) :: metric_sg
        !
        select type( v_node => self%v_node )
            !
            class is( rScalar3D_MR_t )
                !
                do i = 1, self%grid_mr%n_grids
                    !
                    metric_sg = MetricElements_SG_t( self%grid_mr%sub_grid(i) )
                    !
                    v_node%sub_scalar(i) = metric_sg%v_node
                    !
                enddo
                !
            class default
                call errStop( "setNodeVolume_MetricElements_MR > Unclassified v_node" )
            !
        end select
        !
    end subroutine setNodeVolume_MetricElements_MR
    !
    !> No subroutine briefing
    !
    subroutine setCellVolume_MetricElements_MR( self )
        implicit none
        !
        class( MetricElements_MR_t ), intent( inout ) :: self
        !
        integer :: i
        type( MetricElements_SG_t ) :: metric_sg
        !
        select type( v_cell => self%v_cell )
            !
            class is( rScalar3D_MR_t )
                !
                do i = 1, self%grid_mr%n_grids
                    !
                    metric_sg = MetricElements_SG_t( self%grid_mr%sub_grid(i) )
                    !
                    v_cell%sub_scalar(i) = metric_sg%v_cell
                    !
                enddo
                !
            class default
                call errStop( "setCellVolume_MetricElements_MR > Unclassified v_cell" )
            !
        end select
        !
    end subroutine setCellVolume_MetricElements_MR
    !
    !> For a given type find indexes for boundary and interior nodes
    !
    subroutine setGridIndexArrays_MetricElements_MR( self, grid )
        implicit none
        !
        class( MetricElements_MR_t ), intent( in ) :: self
        class( Grid_t ), intent( inout ) :: grid
        !
        integer :: i_grid
        class( Field_t ), allocatable :: temp_field
        !
        select type( grid )
            !
            class is( Grid3D_MR_t )
                !
                allocate( temp_field, source = rVector3D_MR_t( grid, EDGE ) )
                call temp_field%setIndexArrays( grid%EDGEb, grid%EDGEi, grid%EDGEa )
                deallocate( temp_field )
                !
                allocate( temp_field, source = rVector3D_MR_t( grid, FACE ) )
                call temp_field%setIndexArrays( grid%FACEb, grid%FACEi, grid%FACEa )
                deallocate( temp_field )
                !
                allocate( temp_field, source = rScalar3D_MR_t( grid, NODE ) )
                call temp_field%setIndexArrays( grid%NODEb, grid%NODEi, grid%NODEa )
                deallocate( temp_field )
                !
                do i_grid = 1, grid%n_grids
                    !
                    call self%setGridIndexArrays( grid%sub_grid(i_grid) )
                    !
                enddo
                !
            class is( Grid3D_SG_t )
                !
                allocate( temp_field, source = rVector3D_SG_t( grid, EDGE ) )
                call temp_field%setIndexArrays( grid%EDGEb, grid%EDGEi )
                deallocate( temp_field )
                !
                allocate( temp_field, source = rVector3D_SG_t( grid, FACE ) )
                call temp_field%setIndexArrays( grid%FACEb, grid%FACEi )
                deallocate( temp_field )
                !
                allocate( temp_field, source = rScalar3D_SG_t( grid, NODE ) )
                call temp_field%setIndexArrays( grid%NODEb, grid%NODEi )
                deallocate( temp_field )
                !
            class default
                call errStop( "setGridIndexArrays_MetricElements_MR > Unclassified v_cell" )
            !
        end select
        !
    end subroutine setGridIndexArrays_MetricElements_MR
    !
    !> For a given type find indexes for boundary and interior nodes
    ! !
    ! subroutine boundaryIndex_MetricElements_MR( self, grid_type, INDb, INDi )
        ! implicit none
        ! !
        ! class( MetricElements_MR_t ), intent( in ) :: self
        ! character(*), intent( in ) :: grid_type
        ! integer, allocatable, dimension(:), intent( inout ) :: INDb, INDi
        ! !
        ! !call errStop( "boundaryIndex_MetricElements_MR > Under implementation!" )
        ! !
        ! integer :: nVec(3), nVecT, nBdry, n, nb, ni, i
        ! type( rVector3D_MR_t ) :: temp_vector
        ! type( rScalar3D_MR_t ) :: temp_scalar
        ! complex( kind=prec ), allocatable, dimension(:) :: array
        ! !
        ! n = self%grid_mr%getNGrids()
        ! !
        ! selectcase( grid_type )
            ! !
            ! case( EDGE )
                ! !
                ! temp_vector = rVector3D_MR_t( self%grid_mr, EDGE )
                ! !
                ! nVec(1) = 0
                ! nVec(2) = 0
                ! nVec(3) = 0
                ! !
                ! do i = 1, n
                    ! !
                    ! nVec(1) = nVec(1) + size( temp_vector%sub_vector(i)%x )
                    ! nVec(2) = nVec(2) + size( temp_vector%sub_vector(i)%y )
                    ! nVec(3) = nVec(3) + size( temp_vector%sub_vector(i)%z )
                    ! !
                    ! !> Top sub_vector: Ignore bottom boundaries
                    ! if( i == 1 ) then
                        ! !
                        ! temp_vector%sub_vector(i)%x(:, 1, :) = 1
                        ! temp_vector%sub_vector(i)%x(:, temp_vector%ny+1, :) = 1
                        ! temp_vector%sub_vector(i)%x(:, :, 1) = 1
                        ! !temp_vector%sub_vector(i)%x(:, :, temp_vector%nz+1) = 1
                        ! temp_vector%sub_vector(i)%y(1, :, :) = 1
                        ! temp_vector%sub_vector(i)%y(temp_vector%nx+1, :, :) = 1
                        ! temp_vector%sub_vector(i)%y(:, :, 1) = 1
                        ! !temp_vector%sub_vector(i)%y(:, :, temp_vector%nz+1) = 1
                    ! !
                    ! !> Bottom sub_vector: Ignore top boundaries
                    ! elseif( i == n ) then
                        ! !
                        ! temp_vector%sub_vector(i)%x(:, 1, :) = 1
                        ! temp_vector%sub_vector(i)%x(:, temp_vector%ny+1, :) = 1
                        ! !temp_vector%sub_vector(i)%x(:, :, 1) = 1
                        ! temp_vector%sub_vector(i)%x(:, :, temp_vector%nz+1) = 1
                        ! temp_vector%sub_vector(i)%y(1, :, :) = 1
                        ! temp_vector%sub_vector(i)%y(temp_vector%nx+1, :, :) = 1
                        ! !temp_vector%sub_vector(i)%y(:, :, 1) = 1
                        ! temp_vector%sub_vector(i)%y(:, :, temp_vector%nz+1) = 1
                    ! !
                    ! !> Middle sub_vector: Ignore top and bottom boundaries
                    ! else
                        ! !
                        ! temp_vector%sub_vector(i)%x(:, 1, :) = 1
                        ! temp_vector%sub_vector(i)%x(:, temp_vector%ny+1, :) = 1
                        ! !temp_vector%sub_vector(i)%x(:, :, 1) = 1
                        ! !temp_vector%sub_vector(i)%x(:, :, temp_vector%nz+1) = 1
                        ! temp_vector%sub_vector(i)%y(1, :, :) = 1
                        ! temp_vector%sub_vector(i)%y(temp_vector%nx+1, :, :) = 1
                        ! !temp_vector%sub_vector(i)%y(:, :, 1) = 1
                        ! !temp_vector%sub_vector(i)%y(:, :, temp_vector%nz+1) = 1
                        ! !
                    ! endif
                    ! !
                    ! temp_vector%sub_vector(i)%z(1, :, :) = 1
                    ! temp_vector%sub_vector(i)%z(temp_vector%nx+1, :, :) = 1
                    ! temp_vector%sub_vector(i)%z(:, 1, :) = 1
                    ! temp_vector%sub_vector(i)%z(:, temp_vector%ny+1, :) = 1
                    ! !
                ! enddo
                ! !
                ! nVecT = nVec(1) + nVec(2) + nVec(3)
                ! !
                ! array = temp_vector%getArray()
                ! !
            ! case( FACE )
                ! !
                ! temp_vector = rVector3D_MR_t( self%grid_mr, FACE )
                ! !
                ! nVec(1) = 0
                ! nVec(2) = 0
                ! nVec(3) = 0
                ! !
                ! do i = 1, n
                    ! !
                    ! nVec(1) = nVec(1) + size( temp_vector%sub_vector(i)%x )
                    ! nVec(2) = nVec(2) + size( temp_vector%sub_vector(i)%y )
                    ! nVec(3) = nVec(3) + size( temp_vector%sub_vector(i)%z )
                    ! !
                    ! temp_vector%sub_vector(i)%x(1, :, :) = 1
                    ! temp_vector%sub_vector(i)%x(temp_vector%nx+1, :, :) = 1
                    ! temp_vector%sub_vector(i)%y(:, 1, :) = 1
                    ! temp_vector%sub_vector(i)%y(:, temp_vector%ny+1, :) = 1
                    ! !
                    ! !> Top sub_vector: Ignore bottom boundaries
                    ! if( i == 1 ) then
                        ! !
                        ! temp_vector%sub_vector(i)%z(:, :, 1) = 1
                        ! !temp_vector%sub_vector(i)%z(:, :, temp_vector%nz+1) = 1
                        ! !
                    ! !> Bottom sub_vector: Ignore top boundaries
                    ! elseif( i == n ) then
                        ! !
                        ! !temp_vector%sub_vector(i)%z(:, :, 1) = 1
                        ! temp_vector%sub_vector(i)%z(:, :, temp_vector%nz+1) = 1
                        ! !
                    ! endif
                    ! !
                ! enddo
                ! !
                ! nVecT = nVec(1) + nVec(2) + nVec(3)
                ! !
                ! array = temp_vector%getArray()
                ! !
            ! case( NODE )
                ! !
                ! temp_scalar = rScalar3D_MR_t( self%grid_mr, NODE )
                ! !
                ! nVecT = 0
                ! !
                ! do i = 1, n
                    ! !
                    ! nVecT = nVecT + size( temp_scalar%sub_scalar(i)%v )
                    ! !
                    ! temp_scalar%sub_scalar(i)%v(1, :, :) = 1
                    ! temp_scalar%sub_scalar(i)%v(temp_scalar%nx+1, :, :) = 1
                    ! temp_scalar%sub_scalar(i)%v(:, 1, :) = 1
                    ! temp_scalar%sub_scalar(i)%v(:, temp_scalar%ny+1, :) = 1
                    ! !
                    ! !> Top sub_vector: Ignore bottom boundaries
                    ! if( i == 1 ) then
                        ! !
                        ! temp_scalar%sub_scalar(i)%v(:, :, 1) = 1
                        ! !temp_scalar%sub_vector(i)%v(:, :, temp_scalar%nz+1) = 1
                        ! !
                    ! !> Bottom sub_vector: Ignore top boundaries
                    ! elseif( i == n ) then
                        ! !
                        ! !temp_scalar%sub_vector(i)%v(:, :, 1) = 1
                        ! temp_scalar%sub_scalar(i)%v(:, :, temp_scalar%nz+1) = 1
                        ! !
                    ! endif
                    ! !
                ! enddo
                ! !
                ! array = temp_scalar%getArray()
                ! !
            ! case default
                ! call errStop( "boundaryIndex_MetricElements_MR > Invalid grid_mr type ["//grid_type//"]" )
                ! !
        ! end select 
        ! !
        ! nBdry = 0
        ! do i = 1, nVecT
            ! nBdry = nBdry + nint( real( array(i), kind=prec ) )
        ! enddo
        ! !
        ! if( allocated( INDi ) ) then
            ! deallocate( INDi )
        ! endif
        ! !
        ! allocate( INDi( nVecT - nBdry ) )
        ! !
        ! if( allocated( INDb ) ) then
            ! deallocate( INDb )
        ! endif
        ! !
        ! allocate( INDb( nBdry ) )
        ! !
        ! nb = 0
        ! ni = 0
        ! !
        ! do i = 1, nVecT
            ! !
            ! if( nint( real( array(i), kind=prec ) ) .EQ. 1 ) then
                ! nb = nb+1
                ! INDb(nb) = i
            ! else
                ! ni = ni+1
                ! INDi(ni) = i
            ! endif
            ! !
        ! enddo
        ! !
    ! end subroutine boundaryIndex_MetricElements_MR
    !
    !> Create proper scalar from the Grid
    !
    subroutine createScalar_MetricElements_MR( self, scalar_type, grid_type, scalar )
        implicit none
        !
        class( MetricElements_MR_t ), intent( in ) :: self
        integer, intent( in ) :: scalar_type
        character( len=4 ), intent( in ) :: grid_type
        class( Scalar_t ), allocatable, intent( out ) :: scalar
        !
        if( grid_type /= NODE .AND. grid_type /= CELL .AND. grid_type /= CELL_EARTH ) then
            call errStop( "createScalar_MetricElements_MR > grid_type must be NODE, CELL or CELL_EARTH" )
        else
            !
            if( scalar_type == real_t ) then
                allocate( scalar, source = rScalar3D_MR_t( self%grid_mr, grid_type ) )
            elseif( scalar_type == complex_t ) then
                allocate( scalar, source = cScalar3D_MR_t( self%grid_mr, grid_type ) )
            elseif( scalar_type == integer_t ) then
                !allocate( scalar, source = iScalar3D_MR_t( self%grid_mr, grid_type ) )
                !
                call errStop( "createScalar_MetricElements_MR > MR integer_t to be implemented" )
            else
                call errStop( "createScalar_MetricElements_MR > choose MR real_t, complex_t or integer_t" )
            endif
            !
        endif
        !
    end subroutine createScalar_MetricElements_MR
    !
    !> Create proper vector from the Grid
    !
    subroutine createVector_MetricElements_MR( self, vector_type, grid_type, vector )
        implicit none
        !
        class( MetricElements_MR_t ), intent( in ) :: self
        integer, intent( in ) :: vector_type
        character( len=4 ), intent( in ) :: grid_type
        class( Vector_t ), allocatable, intent( out ) :: vector
        !
        if( grid_type /= EDGE .AND. grid_type /= FACE ) then
            call errStop( "createVector_MetricElements_MR > grid_type must be EDGE or FACE" )
        else
            !
            if( vector_type == real_t ) then
                !
                allocate( vector, source = rVector3D_MR_t( self%grid_mr, grid_type ) )
                !
            elseif( vector_type == complex_t ) then
                !
                allocate( vector, source = cVector3D_MR_t( self%grid_mr, grid_type ) )
                !
            elseif( vector_type == integer_t ) then
                !allocate( vector, source = iVector3D_MR_t( self%grid_mr, grid_type ) )
                !
                call errStop( "createVector_MetricElements_MR > MR integer_t to be implemented" )
            else
                call errStop( "createVector_MetricElements_MR > choose MR: real_t, complex_t or integer_t" )
            endif
            !
        endif
        !
    end subroutine createVector_MetricElements_MR
    !
end Module MetricElements_MR
