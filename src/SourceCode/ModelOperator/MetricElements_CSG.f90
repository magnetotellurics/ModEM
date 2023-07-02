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
        procedure, public :: SetEdgelength
        procedure, public :: SetFaceArea
        procedure, public :: SetDualEdgelength
        procedure, public :: SetDualFaceArea
        procedure, public :: SetCellVolume
        procedure, public :: SetEdgeVolume
        procedure, public :: SetNodeVolume
        !
        procedure, public :: allocate => allocateMetricElements_CSG
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
        class( Grid3D_SG_t ), target, intent( in ) :: grid
        type( MetricElements_CSG_t ) :: self
        !
        !write( *, * ) "Constructor MetricElements_CSG_t"
        !
        self%grid => grid
        !
        call self%allocate()
        !
        !>    if were going to allocate storage for all, just set all now!
        call self%setMetricElements()
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
        call self%dealloc
        !
    end subroutine MetricElements_CSG_dtor
    !
    !> No subroutine briefing
    subroutine allocateMetricElements_CSG( self )
        implicit none
        !
        class( MetricElements_CSG_t ), intent( inout ) :: self
        !
        allocate( self%Edgelength, source = rVector3D_SG_t( self%grid, EDGE ) )
        allocate( self%FaceArea, source = rVector3D_SG_t( self%grid, FACE ) )
        allocate( self%DualEdgelength, source = rVector3D_SG_t( self%grid, FACE ) )
        allocate( self%DualFaceArea, source = rVector3D_SG_t( self%grid, EDGE ) )
        allocate( self%Vedge, source = rVector3D_SG_t( self%grid, EDGE ) )
        !
        allocate( self%VNode, source = rScalar3D_SG_t( self%grid, NODE ) )
        allocate( self%Vcell, source = rScalar3D_SG_t( self%grid, CELL ) )
        !
    end subroutine allocateMetricElements_CSG
    !
    !> SetEdgelength
    !> Creates line elements defined on edges of the primary grid.
    !> Edge length elements are defined on interior and boundary edges.
    subroutine SetEdgelength( self )
        implicit none
        !
        class( MetricElements_CSG_t ), intent( inout ) :: self
        !
        integer :: ix, iy, iz
        !
        select type( edge_length => self%Edgelength )
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
    end subroutine SetEdgelength
    !
    !> SetDualEdgelength
    !> Creates line elements defined on edges of the dual grid.
    !> Edge length elements are defined on interior and boundary faces.
    !> Note that dual edge lengths are already defined in grid.
    subroutine SetDualEdgelength( self )
        implicit none
        !
        class( MetricElements_CSG_t ), intent( inout ) :: self
        !
        integer :: ix, iy, iz
        !
        select type( dual_edge_length => self%DualEdgelength )
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
    end subroutine SetDualEdgelength
    !
    !> SetFaceArea
    !> Computes surface area elements on faces of the primary grid.
    !> Face surface area elements are defined on interior and boundary faces.
    subroutine SetFaceArea( self )
        implicit none
        !
        class( MetricElements_CSG_t ), intent( inout ) :: self
        !
        integer :: ix, iy, iz
        !
        select type( face_area => self%FaceArea )
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
    end subroutine SetFaceArea
    !
    !> SetDualFaceArea
    !> Computes surface area elements on faces of the dual grid.
    !> Dual Face surface area elements are defined
    !> on interior and boundary edges.
    !> Note: dual edge lengths are already defined
    !> in grid, use these to compute dual-grid face areas.
    subroutine SetDualFaceArea( self )
        implicit none
      !
        class( MetricElements_CSG_t ), intent( inout ) :: self
        !
        integer :: ix, iy, iz
        !
        select type( dual_face_area => self%DualFaceArea )
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
    end subroutine SetDualFaceArea
    !
    !> SetNodeVolume
    subroutine SetNodeVolume( self )
        implicit none
        !
        class( MetricElements_CSG_t ), intent( inout ) :: self
        !
        integer :: i, j, k
        !
        select type( vnode => self%Vnode )
            class is( rScalar3D_SG_t )
               !
               do i = 1, self%grid%nx + 1
                   do j = 1, self%grid%ny + 1
                      do k = 1, self%grid%nz + 1
                          !> note that we are multiplying
                          !> using the distances with corner of a cell as a center
                          vnode%v(i, j, k) = self%grid%del_x(i) * self%grid%del_y(j) * self%grid%del_z(k)
                      enddo
                   enddo
               enddo
               !
        end select
        !
    end subroutine SetNodeVolume
    !
    !> CellVolume
    !> Creates volume elements for grid cells
    !> and stores them as real scalars with gridType=CELL.
    subroutine SetCellVolume( self )
        implicit none
        !
        class( MetricElements_CSG_t ), intent( inout ) :: self
        !
        integer :: i, j, k
        !
        select type( vcell => self%Vcell )
            class is( rScalar3D_SG_t )
               !
               do i = 1, self%grid%nx 
                   do j = 1, self%grid%ny
                      do k = 1, self%grid%nz
                          vcell%v(i, j, k) = self%grid%dx(i) * self%grid%dy(j) * self%grid%dz(k)
                      enddo
                   enddo
               enddo
               !
        end select
        !
    end subroutine SetCellVolume
    !
    !> SetEdgeVolume
    !> Creates volume elements centered around the edges of
    !> the grid, and stores them as real vectors with
    !> gridType = EDGE.
    subroutine SetEdgeVolume( self )
        implicit none
        !
        class( MetricElements_CSG_t ), intent( inout ) :: self

        !> Vedge = Edgelength*DualFaceArea --- simplest implementation
        !> is just to create these (but they might already be created--
        !>     let's assume they are--the way metric element objects are created
        !>      this is always true
        !call self%setEdgelength()
        !call self%setDualFaceArea()
        !
        self%Vedge = self%Edgelength
        call self%Vedge%mult( self%DualFaceArea )
        !
    end subroutine setEdgeVolume
    !
end Module MetricElements_CSG
