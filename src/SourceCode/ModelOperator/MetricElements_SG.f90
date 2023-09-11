!
!> Derived class to define a ModelOperator_MF_SG
!>
!> This computes and stores Metric Elements for finite volume calculations
!> Based on Matlab development, code is take from ModEM module GridCalc.
!> In this Cartesian staggered grid_sg (CSG) version, metric elements are stored a
!> as TVector objects -- can be used to generate an SP version, and to 
!> generalize to MR; doing the same starting from GridCalcS can be used to
!> Create spherical versions
!> NOTE: there are other grid_sg mapping routines in GridCalc that are NOT to
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
        type( Grid3D_SG_t ), pointer :: grid_sg
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
        procedure, public :: setGridIndexArrays => setGridIndexArrays_MetricElements_SG
        !
        procedure, public :: createScalar => createScalar_MetricElements_SG
        procedure, public :: createVector => createVector_MetricElements_SG
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
    function MetricElements_SG_ctor( grid_sg ) result( self )
        implicit none
        !
        type( Grid3D_SG_t ), target :: grid_sg
        !
        type( MetricElements_SG_t ) :: self
        !
        !write( *, * ) "Constructor MetricElements_SG_t"
        !
        self%grid_sg => grid_sg
        !
        call self%alloc
        !
        !> if were going to allocate storage for all, just set all now!
        call self%setup
        !
        call self%setGridIndexArrays( self%grid_sg )
        !
        self%grid => grid_sg
        !
        write( *, * ) "EDGE: ", self%grid%EDGEf, size( self%grid%EDGEb ), size( self%grid%EDGEi )
        write( *, * ) "FACE: ", self%grid%FACEf, size( self%grid%FACEb ), size( self%grid%FACEi )
        write( *, * ) "NODE: ", self%grid%NODEf, size( self%grid%NODEb ), size( self%grid%NODEi )
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
    !> Creates line elements defined on edges of the primary grid_sg.
    !> Edge length elements are defined on interior and boundary edges.
    !
    subroutine setEdgeLength_MetricElements_SG( self )
        implicit none
        !
        class( MetricElements_SG_t ), intent( inout ) :: self
        !
        integer :: ix, iy, iz
        !
        select type( edge_length => self%edge_length )
            !
            class is( rVector3D_SG_t )
                !
                !> x-component edge length elements
                do ix = 1, self%grid_sg%nx
                    edge_length%x(ix, :, :) = self%grid_sg%dx(ix)
                enddo
                !
                !> y-component edge length elements
                do iy = 1, self%grid_sg%ny
                    edge_length%y(:, iy, :) = self%grid_sg%dy(iy)
                enddo
                !
                !> z-component edge length elements
                do iz = 1, self%grid_sg%nz
                    edge_length%z(:, :, iz) = self%grid_sg%dz(iz)
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
    !> Creates line elements defined on edges of the dual grid_sg.
    !> Edge length elements are defined on interior and boundary faces.
    !> Note that dual edge lengths are already defined in grid_sg.
    !
    subroutine setDualEdgeLength_MetricElements_SG( self )
        implicit none
        !
        class( MetricElements_SG_t ), intent( inout ) :: self
        !
        integer :: ix, iy, iz
        !
        select type( dual_edge_length => self%dual_edge_length )
            !
            class is( rVector3D_SG_t )
                !
                !> x-component edge length elements
                do ix = 1, self%grid_sg%nx+1
                    dual_edge_length%x(ix, :, :) = self%grid_sg%del_x(ix)
                enddo
                !
                !> y-component edge length elements
                do iy = 1, self%grid_sg%ny+1
                    dual_edge_length%y(:, iy, :) = self%grid_sg%del_y(iy)
                enddo
                !
                !> z-component edge length elements
                do iz = 1, self%grid_sg%nz+1
                    dual_edge_length%z(:, :, iz) = self%grid_sg%del_z(iz)
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
    !> Computes surface area elements on faces of the primary grid_sg.
    !> Face surface area elements are defined on interior and boundary faces.
    !
    subroutine setFaceArea_MetricElements_SG( self )
        implicit none
        !
        class( MetricElements_SG_t ), intent( inout ) :: self
        !
        integer :: ix, iy, iz
        !
        select type( face_area => self%face_area )
            !
            class is( rVector3D_SG_t )
                !
                !> x-components
                do ix = 1, self%grid_sg%nx + 1
                    do iy = 1, self%grid_sg%ny
                        do iz = 1, self%grid_sg%nz
                            face_area%x(ix, iy, iz) = self%grid_sg%dy(iy) * self%grid_sg%dz(iz)
                        enddo
                    enddo
                enddo
                !
                !> y-components
                do ix = 1, self%grid_sg%nx
                    do iy = 1, self%grid_sg%ny + 1
                        do iz = 1,self%grid_sg%nz
                            face_area%y(ix, iy, iz) = self%grid_sg%dx(ix) * self%grid_sg%dz(iz)
                        enddo
                    enddo
                enddo
                !
                !> z-components
                do ix = 1, self%grid_sg%nx
                    do iy = 1, self%grid_sg%ny
                        do iz = 1, self%grid_sg%nz + 1
                            face_area%z(ix, iy, iz) = self%grid_sg%dx(ix) * self%grid_sg%dy(iy)
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
    !> Computes surface area elements on faces of the dual grid_sg.
    !> Dual Face surface area elements are defined
    !> on interior and boundary edges.
    !> Note: dual edge lengths are already defined
    !> in grid_sg, use these to compute dual-grid_sg face areas.
    !
    subroutine setDualFaceArea_MetricElements_SG( self )
        implicit none
        !
        class( MetricElements_SG_t ), intent( inout ) :: self
        !
        integer :: ix, iy, iz
        !
        select type( dual_face_area => self%dual_face_area )
            !
            class is( rVector3D_SG_t )
                !
                !> x-components
                do ix = 1, self%grid_sg%nx
                    do iy = 1, self%grid_sg%ny+1
                        do iz = 1, self%grid_sg%nz+1
                            dual_face_area%x(ix, iy, iz) = self%grid_sg%del_y(iy) * self%grid_sg%del_z(iz)
                        enddo
                    enddo
                enddo
                !
                !> y-components
                do ix = 1,self%grid_sg%nx + 1
                    do iy = 1,self%grid_sg%ny
                        do iz = 1,self%grid_sg%nz + 1
                            dual_face_area%y(ix, iy, iz) = self%grid_sg%del_x(ix) * self%grid_sg%del_z(iz)
                        enddo
                    enddo
                enddo
                !
                !> z-components
                do ix = 1, self%grid_sg%nx + 1
                    do iy = 1, self%grid_sg%ny + 1
                        do iz = 1, self%grid_sg%nz
                            dual_face_area%z(ix, iy, iz) = self%grid_sg%del_x(ix) * self%grid_sg%del_y(iy)
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
    !> the grid_sg, and stores them as real vectors with
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
        select type( v_node => self%v_node )
            !
            class is( rScalar3D_SG_t )
                !
                do i = 1, self%grid_sg%nx + 1
                    do j = 1, self%grid_sg%ny + 1
                        do k = 1, self%grid_sg%nz + 1
                            !> note that we are multiplying
                            !> using the distances with corner of a cell as a center
                            v_node%v(i, j, k) = self%grid_sg%del_x(i) * self%grid_sg%del_y(j) * self%grid_sg%del_z(k)
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
    !> Creates volume elements for grid_sg cells
    !> and stores them as real scalars with grid_type=CELL.
    !
    subroutine setCellVolume_MetricElements_SG( self )
        implicit none
        !
        class( MetricElements_SG_t ), intent( inout ) :: self
        !
        integer :: i, j, k
        !
        select type( v_cell => self%v_cell )
            !
            class is( rScalar3D_SG_t )
                !
                do i = 1, self%grid_sg%nx 
                    do j = 1, self%grid_sg%ny
                        do k = 1, self%grid_sg%nz
                            v_cell%v(i, j, k) = self%grid_sg%dx(i) * self%grid_sg%dy(j) * self%grid_sg%dz(k)
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
    subroutine setGridIndexArrays_MetricElements_SG( self, grid )
        implicit none
        !
        class( MetricElements_SG_t ), intent( in ) :: self
        class( Grid_t ), intent( inout ) :: grid
        class( Field_t ), allocatable :: temp_field
        !
        allocate( temp_field, source = rVector3D_SG_t( grid, EDGE ) )
        call temp_field%setIndexArrays( grid%EDGEf, grid%EDGEb, grid%EDGEi )
        deallocate( temp_field )
        !
        allocate( temp_field, source = rVector3D_SG_t( grid, FACE ) )
        call temp_field%setIndexArrays( grid%FACEf, grid%FACEb, grid%FACEi )
        deallocate( temp_field )
        !
        allocate( temp_field, source = rScalar3D_SG_t( grid, NODE ) )
        call temp_field%setIndexArrays( grid%NODEf, grid%NODEb, grid%NODEi )
        deallocate( temp_field )
        !
    end subroutine setGridIndexArrays_MetricElements_SG
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
                ! allocate( temp_field, source = rVector3D_SG_t( self%grid_sg, EDGE ) )
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
                ! allocate( temp_field, source = rVector3D_SG_t( self%grid_sg, FACE ) )
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
                ! allocate( temp_field, source = rScalar3D_SG_t( self%grid_sg, NODE ) )
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
                ! call errStop( "boundaryIndex_MetricElements_SG > Invalid grid_sg type ["//grid_type//"]" )
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
    !
    !> Create proper scalar from the Grid
    !
    subroutine createScalar_MetricElements_SG( self, scalar_type, grid_type, scalar )
        implicit none
        !
        class( MetricElements_SG_t ), intent( in ) :: self
        integer, intent( in ) :: scalar_type
        character( len=4 ), intent( in ) :: grid_type
        class( Scalar_t ), allocatable, intent( out ) :: scalar
        !
        if( grid_type /= NODE .AND. grid_type /= CELL .AND. grid_type /= CELL_EARTH ) then
            call errStop( "createScalar_MetricElements_SG > grid_type must be NODE, CELL or CELL_EARTH" )
        else
            !
            if( scalar_type == real_t ) then
                allocate( scalar, source = rScalar3D_SG_t( self%grid_sg, grid_type ) )
            elseif( scalar_type == complex_t ) then
                allocate( scalar, source = cScalar3D_SG_t( self%grid_sg, grid_type ) )
            elseif( scalar_type == integer_t ) then
                allocate( scalar, source = iScalar3D_SG_t( self%grid_sg, grid_type ) )
            else
                call errStop( "createScalar > choose SG real_t, complex_t or integer_t" )
            endif
            !
        endif
        !
    end subroutine createScalar_MetricElements_SG
    !
    !> Create proper vector from the Grid
    !
    subroutine createVector_MetricElements_SG( self, vector_type, grid_type, vector )
        implicit none
        !
        class( MetricElements_SG_t ), intent( in ) :: self
        integer, intent( in ) :: vector_type
        character( len=4 ), intent( in ) :: grid_type
        class( Vector_t ), allocatable, intent( out ) :: vector
        !
        if( grid_type /= EDGE .AND. grid_type /= FACE ) then
            call errStop( "createVector_MetricElements_SG > grid_type must be EDGE or FACE" )
        else
            !
            if( vector_type == real_t ) then
                allocate( vector, source = rVector3D_SG_t( self%grid_sg, grid_type ) )
            elseif( vector_type == complex_t ) then
                allocate( vector, source = cVector3D_SG_t( self%grid_sg, grid_type ) )
            elseif( vector_type == integer_t ) then
                !allocate( vector, source = iVector3D_SG_t( self%grid_sg, grid_type ) )
                !
                call errStop( "createVector_MetricElements_SG > SG integer_t to be implemented" )
            else
                call errStop( "createVector_MetricElements_SG > choose SG: real_t, complex_t or integer_t" )
            endif
            !
        endif
        !
    end subroutine createVector_MetricElements_SG
    !
end Module MetricElements_SG
