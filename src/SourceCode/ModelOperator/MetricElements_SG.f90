!
!> Derived class to define a ModelOperator_MF_SG
!>
!> This computes and stores Metric Elements for finite volume calculations
!> Based on Matlab development, code is take from ModEM module GridCalc.
!> In this Cartesian staggered grid (CSG) version, metric elements are stored a
!> as TVector objects -- can be used to generate an SP version, and to 
!> generalize to MR; doing the same starting from GridCalcS can be used to
!> Create spherical versions
!> NOTE: there are other grid mapping routines in GridCalc that are NOT to
!>            be included here -- I think these might be better as methods in the
!>            Vector/Scalar classes.
!
!> Variables that will be defined in base class    ... could add more as
!> in GridCalc -- but let's see if these are really useful.
!
module MetricElements_SG
    !
    use MetricElements
    !
    type, extends( MetricElements_t ) :: MetricElements_SG_t
        !
        !> No derived properties
        !
     contains
        !
        final :: MetricElements_SG_dtor
        !
        procedure, public :: setEdgeLength => setEdgeLength_MetricElements_SG
        !
        procedure, public :: setDualEdgeLength => setDualEdgeLength_MetricElements_SG
        !
        procedure, public :: setFaceArea => setFaceArea_MetricElements_SG
        !
        procedure, public :: setDualFaceArea => setDualFaceArea_MetricElements_SG
        !
        procedure, public :: setEdgeVolume => setEdgeVolume_MetricElements_SG
        !
        procedure, public :: setNodeVolume => setNodeVolume_MetricElements_SG
        !
        procedure, public :: setCellVolume => setCellVolume_MetricElements_SG
        !
        procedure, public :: setIndexArrays => setIndexArrays_MetricElements_SG
        !
        procedure, public :: setAllIndexArrays => setAllIndexArrays_MetricElements_SG
        !
        !procedure, public :: boundaryIndex => boundaryIndex_MetricElements_SG
        !
    end type MetricElements_SG_t
    !
    interface MetricElements_SG_t
        module procedure MetricElements_SG_ctor
    end interface MetricElements_SG_t
    !
contains
    !
    !> No subroutine briefing
    !
    function MetricElements_SG_ctor( grid ) result( self )
        implicit none
        !
        class( Grid_t ), target, intent( in ) :: grid
        type( MetricElements_SG_t ) :: self
        !
        !write( *, * ) "Constructor MetricElements_SG_t"
        !
        self%grid => grid
        !
        call self%alloc
        !
        !> if were going to allocate storage for all, just set all now!
        call self%setup
        !
    end function MetricElements_SG_Ctor
    !
    !> No subroutine briefing
    !
    subroutine MetricElements_SG_dtor( self )
        implicit none
        !
        type( MetricElements_SG_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor MetricElements_SG"
        !
        call self%baseDealloc
        !
    end subroutine MetricElements_SG_dtor
    !
    !> setEdgeLength_MetricElements_SG
    !> Creates line elements defined on edges of the primary grid.
    !> Edge length elements are defined on interior and boundary edges.
    !
    subroutine setEdgeLength_MetricElements_SG( self )
        implicit none
        !
        class( MetricElements_SG_t ), intent( inout ) :: self
        !
        integer :: ix, iy, iz
        !
        !> Instantiate the ModelOperator object
        select type( edge_length => self%edge_length )
            !
            class is( rVector3D_SG_t )
                !
                !> x-component edge length elements
                do ix = 1, self%grid%nx
                    edge_length%x(ix, :, :) = self%grid%dx(ix)
                enddo
                !
                !> y-component edge length elements
                do iy = 1, self%grid%ny
                    edge_length%y(:, iy, :) = self%grid%dy(iy)
                enddo
                !
                !> z-component edge length elements
                do iz = 1, self%grid%nz
                    edge_length%z(:, :, iz) = self%grid%dz(iz)
                enddo
                !
            class default
                call errStop( "setEdgeLength_MetricElements_SG > Unclassified target" )
            !
        end select
        !
    end subroutine setEdgeLength_MetricElements_SG
    !
    !> setDualEdgeLength_MetricElements_SG
    !> Creates line elements defined on edges of the dual grid.
    !> Edge length elements are defined on interior and boundary faces.
    !> Note that dual edge lengths are already defined in grid.
    !
    subroutine setDualEdgeLength_MetricElements_SG( self )
        implicit none
        !
        class( MetricElements_SG_t ), intent( inout ) :: self
        !
        integer :: ix, iy, iz
        !
        !> Instantiate the ModelOperator object
        select type( dual_edge_length => self%dual_edge_length )
            !
            class is( rVector3D_SG_t )
                !
                !> x-component edge length elements
                do ix = 1, self%grid%nx+1
                    dual_edge_length%x(ix, :, :) = self%grid%del_x(ix)
                enddo
                !
                !> y-component edge length elements
                do iy = 1, self%grid%ny+1
                    dual_edge_length%y(:, iy, :) = self%grid%del_y(iy)
                enddo
                !
                !> z-component edge length elements
                do iz = 1, self%grid%nz+1
                    dual_edge_length%z(:, :, iz) = self%grid%del_z(iz)
                enddo
                !
            class default
                call errStop( "setDualEdgeLength_MetricElements_SG > Unclassified target" )
            !
        end select
        !
    end subroutine setDualEdgeLength_MetricElements_SG
    !
    !> setFaceArea_MetricElements_SG
    !> Computes surface area elements on faces of the primary grid.
    !> Face surface area elements are defined on interior and boundary faces.
    !
    subroutine setFaceArea_MetricElements_SG( self )
        implicit none
        !
        class( MetricElements_SG_t ), intent( inout ) :: self
        !
        integer :: ix, iy, iz
        !
        !> Instantiate the ModelOperator object
        select type( face_area => self%face_area )
            !
            class is( rVector3D_SG_t )
                !
                !> x-components
                do ix = 1, self%grid%nx + 1
                    do iy = 1, self%grid%ny
                        do iz = 1, self%grid%nz
                            face_area%x(ix, iy, iz) = self%grid%dy(iy) * self%grid%dz(iz)
                        enddo
                    enddo
                enddo
                !
                !> y-components
                do ix = 1, self%grid%nx
                    do iy = 1, self%grid%ny + 1
                        do iz = 1,self%grid%nz
                            face_area%y(ix, iy, iz) = self%grid%dx(ix) * self%grid%dz(iz)
                        enddo
                    enddo
                enddo
                !
                !> z-components
                do ix = 1, self%grid%nx
                    do iy = 1, self%grid%ny
                        do iz = 1, self%grid%nz + 1
                            face_area%z(ix, iy, iz) = self%grid%dx(ix) * self%grid%dy(iy)
                        enddo
                    enddo
                enddo
                !
            class default
                call errStop( "setFaceArea_MetricElements_SG > Unclassified target" )
            !
        end select
        !
    end subroutine setFaceArea_MetricElements_SG
    !
    !> setDualFaceArea_MetricElements_SG
    !> Computes surface area elements on faces of the dual grid.
    !> Dual Face surface area elements are defined
    !> on interior and boundary edges.
    !> Note: dual edge lengths are already defined
    !> in grid, use these to compute dual-grid face areas.
    !
    subroutine setDualFaceArea_MetricElements_SG( self )
        implicit none
        !
        class( MetricElements_SG_t ), intent( inout ) :: self
        !
        integer :: ix, iy, iz
        !
        !> Instantiate the ModelOperator object
        select type( dual_face_area => self%dual_face_area )
            !
            class is( rVector3D_SG_t )
                !
                !> x-components
                do ix = 1, self%grid%nx
                    do iy = 1, self%grid%ny+1
                        do iz = 1, self%grid%nz+1
                            dual_face_area%x(ix, iy, iz) = self%grid%del_y(iy) * self%grid%del_z(iz)
                        enddo
                    enddo
                enddo
                !
                !> y-components
                do ix = 1,self%grid%nx + 1
                    do iy = 1,self%grid%ny
                        do iz = 1,self%grid%nz + 1
                            dual_face_area%y(ix, iy, iz) = self%grid%del_x(ix) * self%grid%del_z(iz)
                        enddo
                    enddo
                enddo
                !
                !> z-components
                do ix = 1, self%grid%nx + 1
                    do iy = 1, self%grid%ny + 1
                        do iz = 1, self%grid%nz
                            dual_face_area%z(ix, iy, iz) = self%grid%del_x(ix) * self%grid%del_y(iy)
                        enddo
                    enddo
                enddo
                !

            class default
                call errStop( "setDualFaceArea_MetricElements_SG > Unclassified target" )
            !
        end select
        !
    end subroutine setDualFaceArea_MetricElements_SG
    !
    !> setEdgeVolume_MetricElements_SG
    !> Creates volume elements centered around the edges of
    !> the grid, and stores them as real vectors with
    !> grid_type = EDGE.
    !
    !> v_edge = edge_length*dual_face_area --- simplest implementation
    !> is just to create these (but they might already be created--
    !>     let's assume they are--the way metric element objects are created
    !>      this is always true
    !call self%setEdgeLength_MetricElements_SG()
    !call self%setDualFaceArea_MetricElements_SG()
    !
    subroutine setEdgeVolume_MetricElements_SG( self )
        implicit none
        !
        class( MetricElements_SG_t ), intent( inout ) :: self
        !
        self%v_edge = self%edge_length
        !
        call self%v_edge%mult( self%dual_face_area )
        !
    end subroutine setEdgeVolume_MetricElements_SG
    !
    !> setNodeVolume_MetricElements_SG
    !
    subroutine setNodeVolume_MetricElements_SG( self )
        implicit none
        !
        class( MetricElements_SG_t ), intent( inout ) :: self
        !
        integer :: i, j, k
        !
        !> Instantiate the ModelOperator object
        select type( v_node => self%v_node )
            !
            class is( rScalar3D_SG_t )
                !
                do i = 1, self%grid%nx + 1
                    do j = 1, self%grid%ny + 1
                        do k = 1, self%grid%nz + 1
                            !> note that we are multiplying
                            !> using the distances with corner of a cell as a center
                            v_node%v(i, j, k) = self%grid%del_x(i) * self%grid%del_y(j) * self%grid%del_z(k)
                            !
                        enddo
                    enddo
                enddo
                !
            class default
                call errStop( "setNodeVolume_MetricElements_SG > Unclassified target" )
            !
        end select
        !
    end subroutine setNodeVolume_MetricElements_SG
    !
    !> CellVolume
    !> Creates volume elements for grid cells
    !> and stores them as real scalars with grid_type=CELL.
    !
    subroutine setCellVolume_MetricElements_SG( self )
        implicit none
        !
        class( MetricElements_SG_t ), intent( inout ) :: self
        !
        integer :: i, j, k
        !
        !> Instantiate the ModelOperator object
        select type( v_cell => self%v_cell )
            !
            class is( rScalar3D_SG_t )
                !
                do i = 1, self%grid%nx 
                    do j = 1, self%grid%ny
                        do k = 1, self%grid%nz
                            v_cell%v(i, j, k) = self%grid%dx(i) * self%grid%dy(j) * self%grid%dz(k)
                        enddo
                    enddo
                enddo
                !
            class default
                call errStop( "setCellVolume_MetricElements_SG > Unclassified target" )
            !
        end select
        !
    end subroutine setCellVolume_MetricElements_SG
    !
    !> For a given type find indexes for boundary and interior nodes
    !
    subroutine setIndexArrays_MetricElements_SG( self, grid_type, INDb, INDi, INDa )
        implicit none
        !
        class( MetricElements_SG_t ), intent( in ) :: self
        character(*), intent( in ) :: grid_type
        integer, allocatable, dimension(:), intent( out ) :: INDb, INDi
        integer, dimension(:), allocatable, intent( out ), optional :: INDa
        !
        class( Field_t ), allocatable :: temp_field
        !
        selectcase( grid_type )
            !
            case( EDGE, FACE )
                !
                allocate( temp_field, source = rVector3D_SG_t( self%grid, grid_type ) )
                !
            case( NODE )
                !
                allocate( temp_field, source = rScalar3D_SG_t( self%grid, grid_type ) )
                !
            case default
                call errStop( "setIndexArrays_MetricElements_SG > Invalid grid type ["//grid_type//"]" )
        end select 
        !
        call temp_field%setIndexArrays( INDb, INDi )
        !
    end subroutine setIndexArrays_MetricElements_SG
    !
    !> For a given type find indexes for boundary and interior nodes
    !
    subroutine setAllIndexArrays_MetricElements_SG( self )
        implicit none
        !
        class( MetricElements_SG_t ), intent( in ) :: self
        !
        call self%setIndexArrays( EDGE, self%grid%EDGEb, self%grid%EDGEi )
        write( *, "( a17, i8, a8, i8 )" ) "EDGEb=", size( self%grid%EDGEb ), ", EDGEi=", size( self%grid%EDGEi )
        !
        call self%setIndexArrays( FACE, self%grid%FACEb, self%grid%FACEi )
        write( *, "( a17, i8, a8, i8 )" ) "FACEb=", size( self%grid%FACEb ), ", FACEi=", size( self%grid%FACEi )
        !
        call self%setIndexArrays( NODE, self%grid%NODEb, self%grid%NODEi )
        write( *, "( a17, i8, a8, i8 )" ) "NODEb=", size( self%grid%NODEb ), ", NODEi=", size( self%grid%NODEi )
        !
    end subroutine setAllIndexArrays_MetricElements_SG
    !
    !> For a given type find indexes for boundary and interior nodes
    ! !
    ! subroutine boundaryIndex_MetricElements_SG( self, grid_type, INDb, INDi )
        ! implicit none
        ! !
        ! class( MetricElements_SG_t ), intent( in ) :: self
        ! character(*), intent( in ) :: grid_type
        ! integer, allocatable, dimension(:), intent( inout ) :: INDb, INDi
        ! !
        ! integer :: nVec(3), nVecT, nBdry, nb, ni, i
        ! class( Field_t ), allocatable :: temp_field
        ! complex( kind=prec ), allocatable, dimension(:) :: array
        ! !
        ! selectcase( grid_type )
            ! !
            ! case( EDGE )
                ! !
                ! allocate( temp_field, source = rVector3D_SG_t( self%grid, EDGE ) )
                ! !
                ! select type( temp_field )
                    ! !
                    ! class is( rVector3D_SG_t )
                        ! !
                        ! nVec(1) = size( temp_field%x )
                        ! nVec(2) = size( temp_field%y )
                        ! nVec(3) = size( temp_field%z )
                        ! !
                        ! nVecT = nVec(1)+nVec(2)+nVec(3)
                        ! !
                        ! temp_field%x(:, 1, :) = 1
                        ! temp_field%x(:, temp_field%ny+1, :) = 1
                        ! temp_field%x(:, :, 1) = 1
                        ! temp_field%x(:, :, temp_field%nz+1) = 1
                        ! temp_field%y(1, :, :) = 1
                        ! temp_field%y(temp_field%nx+1, :, :) = 1
                        ! temp_field%y(:, :, 1) = 1
                        ! temp_field%y(:, :, temp_field%nz+1) = 1
                        ! temp_field%z(1, :, :) = 1
                        ! temp_field%z(temp_field%nx+1, :, :) = 1
                        ! temp_field%z(:, 1, :) = 1
                        ! temp_field%z(:, temp_field%ny+1, :) = 1
                        ! !
                    ! class default
                        ! call errStop( "boundaryIndex_MetricElements_SG > Undefined EDGE" )
                        ! !
                ! end select
                ! !
            ! case( FACE )
                ! !
                ! allocate( temp_field, source = rVector3D_SG_t( self%grid, FACE ) )
                ! !
                ! select type( temp_field )
                    ! !
                    ! class is( rVector3D_SG_t )
                        ! !
                        ! nVec(1) = size( temp_field%x )
                        ! nVec(2) = size( temp_field%y )
                        ! nVec(3) = size( temp_field%z )
                        ! !
                        ! nVecT = nVec(1)+nVec(2)+nVec(3)
                        ! !
                        ! temp_field%x(1, :, :) = 1
                        ! temp_field%x(temp_field%nx+1, :, :) = 1
                        ! temp_field%y(:, 1, :) = 1
                        ! temp_field%y(:, temp_field%ny+1, :) = 1
                        ! temp_field%z(:, :, 1) = 1
                        ! temp_field%z(:, :, temp_field%nz+1) = 1
                        ! !
                    ! class default
                        ! call errStop( "boundaryIndex_MetricElements_SG > Undefined FACE" )
                        ! !
                ! end select
                ! !
            ! case( NODE )
                ! !
                ! allocate( temp_field, source = rScalar3D_SG_t( self%grid, NODE ) )
                ! !
                ! select type( temp_field )
                    ! !
                    ! class is( rScalar3D_SG_t )
                        ! !
                        ! nVecT = size( temp_field%v )
                        ! !
                        ! temp_field%v(1, :, :) = 1
                        ! temp_field%v(temp_field%nx+1, :, :) = 1
                        ! temp_field%v(:, 1, :) = 1
                        ! temp_field%v(:, temp_field%ny+1, :) = 1
                        ! temp_field%v(:, :, 1) = 1
                        ! temp_field%v(:, :, temp_field%nz+1) = 1
                        ! !
                    ! class default
                        ! call errStop( "boundaryIndex_MetricElements_SG > Undefined NODE" )
                        ! !
                ! end select
                ! !
            ! case default
                ! call errStop( "boundaryIndex_MetricElements_SG > Invalid grid type ["//grid_type//"]" )
        ! end select 
        ! !
        ! array = temp_field%getArray()
        ! !
        ! deallocate( temp_field )
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
    ! end subroutine boundaryIndex_MetricElements_SG
    ! !
end Module MetricElements_SG
