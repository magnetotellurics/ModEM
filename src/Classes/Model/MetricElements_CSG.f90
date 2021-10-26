!**
! This computes and stores Metric Elements for finite volume calculations
! Based on Matlab development, code is take from ModEM module GridCalc.
! In this Cartesian staggered grid (CSG) version, metric elements are stored a
! as TVector objects -- can be used to generate an SP version, and to 
! generalize to MR; doing the same starting from GridCalcS can be used to
! Create spherical versions
! NOTE: there are other grid mapping routines in GridCalc that are NOT to
!       be included here -- I think these might be better as methods in the
!       Vector/Scalar classes.
!
! Variables that will be defined in base class  ... could add more as
! in GridCalc -- but let's see if these are really useful.
!
!*
module MetricElements_CSG
  use Constants
  use Grid
  use Grid3D_SG
  use MetricElements
  use rVector3D_SG
  use rScalar3D_SG
  
  type, extends(MetricElements_t) :: MetricElements_CSG_t
     class(Grid3D_SG_t), pointer :: grid
     
     type(rVector3D_SG_t) :: EdgeLength
     type(rVector3D_SG_t) :: FaceArea
     type(rVector3D_SG_t) :: DualFaceArea
     type(rVector3D_SG_t) :: DualEdgeLength
     type(rScalar3D_SG_t) :: Vnode
     type(rScalar3D_SG_t) :: Vcell
     type(rVector3D_SG_t) :: Vedge
     
     integer :: nx = 0
     integer :: ny = 0
     integer :: nz = 0
     
   contains
     
     procedure, public :: SetEdgeLength
     procedure, public :: SetFaceArea
     procedure, public :: SetDualEdgeLength
     procedure, public :: SetDualFaceArea
     procedure, public :: SetCellVolume
     procedure, public :: SetEdgeVolume
     procedure, public :: SetNodeVolume

     procedure, public :: create     => createMetricElements_CSG
     procedure, public :: allocate   => allocateMetricElements_CSG
     procedure, public :: deallocate => deallocateMetricElements_CSG
     
  end type MetricElements_CSG_t
  
  interface MetricElements_CSG_t
     module procedure MetricElements_CSG_ctor
  end interface MetricElements_CSG_t
  
contains

  !**
  ! MetricElements_CSG constructor
  !*
  function MetricElements_CSG_ctor(inGrid) result(Metric)
    implicit none
    class(Grid3D_SG_t), intent(in) :: inGrid
    type(MetricElements_CSG_t) :: Metric
    !
    call Metric%create( inGrid )
    
  end function MetricElements_CSG_Ctor
  
  !**
  ! createMetricElements_CSG
  !*
  subroutine createMetricElements_CSG(self, inGrid)
    implicit none
    class(MetricElements_CSG_t), intent(inout) :: self
    class(Grid3D_SG_t)         , intent(in), target :: inGrid
    
    self%Grid => inGrid
    self%nx = inGrid%nx
    self%ny = inGrid%ny
    self%nz = inGrid%nz
    
    call self%allocate()
    
    call self%SetFaceArea()    
    call self%SetEdgeLength()
    
  end subroutine createMetricElements_CSG
  
  !**
  ! allocateMetricElements_CSG
  !*
  subroutine allocateMetricElements_CSG(self)
    implicit none
    class(MetricElements_CSG_t), intent(inout) :: self

    self%EdgeLength     = rVector3D_SG_t(self%grid, EDGE)
    self%FaceArea       = rVector3D_SG_t(self%grid, FACE)
    self%DualEdgeLength = rVector3D_SG_t(self%grid, FACE)
    self%DualFaceArea   = rVector3D_SG_t(self%grid, EDGE)
    self%Vedge          = rVector3D_SG_t(self%grid, EDGE)
    self%VNode          = rScalar3D_SG_t(self%grid, NODE)

  end subroutine allocateMetricElements_CSG

  !**
  ! deallocateMetricElements_CSG
  !*
  subroutine deallocateMetricElements_CSG(self)
    implicit none
    class(MetricElements_CSG_t), intent(inout) :: self
    
  end subroutine deallocateMetricElements_CSG
  
  !**
  ! SetEdgeLength
  ! Creates line elements defined on edges of the primary grid.
  ! Edge length elements are defined on interior and boundary edges.
  !*
  subroutine SetEdgeLength(self)
    implicit none
    class(MetricElements_CSG_t), intent(inout) :: self
    ! Local variables
    integer :: ix, iy, iz
    
    ! x-component edge length elements
    do ix = 1, self%nx
       self%EdgeLength%x(ix, :, :) = self%grid%dx(ix)
    end do
    
    ! y-component edge length elements
    do iy = 1, self%ny
       self%EdgeLength%y(:, iy, :) = self%grid%dy(iy)
    end do
    
    ! z-component edge length elements
    do iz = 1, self%nz
       self%EdgeLength%z(:, :, iz) = self%grid%dz(iz)
    end do
  end subroutine SetEdgeLength
  
  !**
  ! SetDualEdgeLength
  ! Creates line elements defined on edges of the dual grid.
  ! Edge length elements are defined on interior and boundary faces.
  ! Note that dual edge lengths are already defined in grid.
  !*
  subroutine SetDualEdgeLength(self)
    implicit none
    class(MetricElements_CSG_t), intent(inout) :: self
    ! Local variables
    integer :: ix, iy, iz
    
    ! x-component edge length elements
    do ix = 1, self%nx
       self%DualEdgeLength%x(ix, :, :) = self%grid%delX(ix)
    end do
    
    ! y-component edge length elements
    do iy = 1, self%ny
       self%DualEdgeLength%y(:, iy, :) = self%grid%delY(iy)
    end do
    
    ! z-component edge length elements
    do iz = 1, self%nz
       self%DualEdgeLength%z(:, :, iz) = self%grid%delZ(iz)
    end do
    
  end subroutine SetDualEdgeLength
  
  !**
  ! SetFaceArea
  ! Computes surface area elements on faces of the primary grid.
  ! Face surface area elements are defined on interior and boundary faces.
  !*
  subroutine SetFaceArea(self)
    implicit none
    class(MetricElements_CSG_t), intent(inout) :: self
    ! Local variables
    integer :: ix, iy, iz
    
    ! x-components
    do ix = 1, self%nx + 1
       do iy = 1, self%ny
          do iz = 1, self%nz
             self%FaceArea%x(ix, iy, iz) = self%grid%dy(iy) * self%grid%dz(iz)
          end do
       end do
    end do
    
    ! y-components
    do ix = 1, self%nx
       do iy = 1, self%ny + 1
          do iz = 1,self%nz
             self%FaceArea%y(ix, iy, iz) = self%grid%dx(ix) * self%grid%dz(iz)
          end do
       end do
    end do
    
    ! z-components
    do ix = 1, self%nx
       do iy = 1, self%ny
          do iz = 1, self%nz + 1
             self%FaceArea%z(ix, iy, iz) = self%grid%dx(ix) * self%grid%dy(iy)
          end do
       end do
    end do
    
  end subroutine SetFaceArea
  
  !**
  ! SetDualFaceArea
  ! Computes surface area elements on faces of the dual grid.
  ! Dual Face surface area elements are defined
  ! on interior and boundary edges.
  ! Note: dual edge lengths are already defined
  ! in grid, use these to compute dual-grid face areas.
  !*
  subroutine SetDualFaceArea(self)
    implicit none
    class(MetricElements_CSG_t), intent(inout) :: self
    ! Local variables
    integer :: ix, iy, iz
    
    ! x-components
    do ix = 1, self%nx + 1
       do iy = 1, self%ny
          do iz = 1, self%nz
             self%DualFaceArea%x(ix, iy, iz) = self%grid%delY(iy) * self%grid%delZ(iz)
          end do
       end do
    end do
    
    ! y-components
    do ix = 1,self%nx
       do iy = 1,self%ny + 1
          do iz = 1,self%nz
             self%DualFaceArea%y(ix, iy, iz) = self%grid%delX(ix) * self%grid%delZ(iz)
          end do
       end do
    end do
    
    ! z-components
    do ix = 1, self%nx
       do iy = 1, self%ny
          do iz = 1, self%nz + 1
             self%DualFaceArea%z(ix, iy, iz) = self%grid%delX(ix) * self%grid%delY(iy)
          end do
       end do
    end do
  end subroutine SetDualFaceArea
  
  !**
  ! SetNodeVolume
  !*
  subroutine SetNodeVolume(self)
    implicit none
    class(MetricElements_CSG_t), intent(inout) :: self
    ! Local variables
    integer :: i, j, k
    
    do i = 1, self%nx + 1
       do j = 1, self%ny + 1
          do k = 1, self%nz + 1
             ! note that we are multiplying
             ! using the distances with corner of a cell as a center
             self%Vnode%v(i, j, k) = self%grid%delX(i) * self%grid%delY(j) * self%grid%delZ(k)
          end do
       end do
    end do
  end subroutine SetNodeVolume
  
  !**
  ! CellVolume
  ! Creates volume elements for grid cells
  ! and stores them as real scalars with gridType=CELL.
  !*
  subroutine SetCellVolume(self)
    implicit none
    class(MetricElements_CSG_t), intent(inout) :: self
    ! Local variables
    integer :: i, j, k
    
    do i = 1, self%nx + 1
       do j = 1, self%ny + 1
          do k = 1, self%nz + 1
             self%Vcell%v(i, j, k) = self%grid%dx(i) * self%grid%dy(j) * self%grid%dz(k)
          end do
       end do
    end do
  end subroutine SetCellVolume
  
  !**
  ! SetEdgeVolume
  ! Creates volume elements centered around the edges of
  ! the grid, and stores them as real vectors with
  ! gridType = EDGE.
  !*
  subroutine SetEdgeVolume(self)
    implicit none
    class(MetricElements_CSG_t), intent(inout) :: self
    
    ! Vedge = EdgeLength*DualFaceArea --- simplest implementation
    ! is just to create these (but they might already be created--
    call self%SetEdgeLength()
    call self%SetDualFaceArea()
    
    self%Vedge = self%EdgeLength * self%DualFaceArea
  end subroutine setEdgeVolume
  
end Module MetricElements_CSG
