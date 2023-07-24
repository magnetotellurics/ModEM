!
!> Derived class to define a ModelOperator_MF
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
module MetricElements_CSG
    !
    use MetricElements
    use cVector3D_SG
    !
    type, extends( MetricElements_t ) :: MetricElements_CSG_t
        !
        !> No derived properties
        !
     contains
        !
        final :: MetricElements_CSG_dtor
        !
        procedure, public :: setEdgeLength
        procedure, public :: setFaceArea
        procedure, public :: setDualEdgeLength
        procedure, public :: setDualFaceArea
        procedure, public :: setCellVolume
        procedure, public :: setEdgeVolume
        procedure, public :: setNodeVolume
        !
        procedure, public :: alloc => allocateMetricElements_CSG
        !
    end type MetricElements_CSG_t
    !
    interface MetricElements_CSG_t
        module procedure MetricElements_CSG_ctor
    end interface MetricElements_CSG_t
    !
contains
    !
    !> No subroutine briefing
    !
    function MetricElements_CSG_ctor( grid ) result( self )
        implicit none
        !
        class( Grid_t ), target, intent( in ) :: grid
        type( MetricElements_CSG_t ) :: self
        !
        !write( *, * ) "Constructor MetricElements_CSG_t"
        !
        self%grid => grid
        !
        call self%alloc
        !
        !>    if were going to allocate storage for all, just set all now!
        call self%setMetricElements
        !
    end function MetricElements_CSG_Ctor
    !
    !> No subroutine briefing
    subroutine MetricElements_CSG_dtor( self )
        implicit none
        !
        type( MetricElements_CSG_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor MetricElements_CSG"
        !
        call self%baseDealloc
        !
    end subroutine MetricElements_CSG_dtor
    !
    !> No subroutine briefing
    subroutine allocateMetricElements_CSG( self )
        implicit none
        !
        class( MetricElements_CSG_t ), intent( inout ) :: self
        !
        allocate( self%edge_length, source = rVector3D_SG_t( self%grid, EDGE ) )
        allocate( self%face_area, source = rVector3D_SG_t( self%grid, FACE ) )
        allocate( self%dual_edge_length, source = rVector3D_SG_t( self%grid, FACE ) )
        allocate( self%dual_face_area, source = rVector3D_SG_t( self%grid, EDGE ) )
        allocate( self%v_edge, source = rVector3D_SG_t( self%grid, EDGE ) )
        !
        allocate( self%v_node, source = rScalar3D_SG_t( self%grid, NODE ) )
        allocate( self%v_cell, source = rScalar3D_SG_t( self%grid, CELL ) )
        !
    end subroutine allocateMetricElements_CSG
    !
    !> setEdgeLength
    !> Creates line elements defined on edges of the primary grid.
    !> Edge length elements are defined on interior and boundary edges.
    subroutine setEdgeLength( self )
        implicit none
        !
        class( MetricElements_CSG_t ), intent( inout ) :: self
        !
        integer :: ix, iy, iz
        !
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
            !
        end select
        !
    end subroutine setEdgeLength
    !
    !> setDualEdgeLength
    !> Creates line elements defined on edges of the dual grid.
    !> Edge length elements are defined on interior and boundary faces.
    !> Note that dual edge lengths are already defined in grid.
    subroutine setDualEdgeLength( self )
        implicit none
        !
        class( MetricElements_CSG_t ), intent( inout ) :: self
        !
        integer :: ix, iy, iz
        !
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
        end select
        !
    end subroutine setDualEdgeLength
    !
    !> setFaceArea
    !> Computes surface area elements on faces of the primary grid.
    !> Face surface area elements are defined on interior and boundary faces.
    subroutine setFaceArea( self )
        implicit none
        !
        class( MetricElements_CSG_t ), intent( inout ) :: self
        !
        integer :: ix, iy, iz
        !
        select type( face_area => self%face_area )
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
        end select
        !
    end subroutine setFaceArea
    !
    !> setDualFaceArea
    !> Computes surface area elements on faces of the dual grid.
    !> Dual Face surface area elements are defined
    !> on interior and boundary edges.
    !> Note: dual edge lengths are already defined
    !> in grid, use these to compute dual-grid face areas.
    subroutine setDualFaceArea( self )
        implicit none
      !
        class( MetricElements_CSG_t ), intent( inout ) :: self
        !
        integer :: ix, iy, iz
        !
        select type( dual_face_area => self%dual_face_area )
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
        end select
        !
    end subroutine setDualFaceArea
    !
    !> setNodeVolume
    subroutine setNodeVolume( self )
        implicit none
        !
        class( MetricElements_CSG_t ), intent( inout ) :: self
        !
        integer :: i, j, k
        !
        select type( v_node => self%v_node )
            class is( rScalar3D_SG_t )
               !
               do i = 1, self%grid%nx + 1
                   do j = 1, self%grid%ny + 1
                      do k = 1, self%grid%nz + 1
                          !> note that we are multiplying
                          !> using the distances with corner of a cell as a center
                          v_node%v(i, j, k) = self%grid%del_x(i) * self%grid%del_y(j) * self%grid%del_z(k)
                      enddo
                   enddo
               enddo
               !
        end select
        !
    end subroutine setNodeVolume
    !
    !> CellVolume
    !> Creates volume elements for grid cells
    !> and stores them as real scalars with grid_type=CELL.
    subroutine setCellVolume( self )
        implicit none
        !
        class( MetricElements_CSG_t ), intent( inout ) :: self
        !
        integer :: i, j, k
        !
        select type( v_cell => self%v_cell )
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
        end select
        !
    end subroutine setCellVolume
    !
    !> setEdgeVolume
    !> Creates volume elements centered around the edges of
    !> the grid, and stores them as real vectors with
    !> grid_type = EDGE.
    subroutine setEdgeVolume( self )
        implicit none
        !
        class( MetricElements_CSG_t ), intent( inout ) :: self

        !> v_edge = edge_length*dual_face_area --- simplest implementation
        !> is just to create these (but they might already be created--
        !>     let's assume they are--the way metric element objects are created
        !>      this is always true
        !call self%setEdgeLength()
        !call self%setDualFaceArea()
        !
        self%v_edge = self%edge_length
        call self%v_edge%mult( self%dual_face_area )
        !
    end subroutine setEdgeVolume
    !
end Module MetricElements_CSG
