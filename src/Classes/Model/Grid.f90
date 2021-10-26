!**
! Definition of the base class for grids. 
! Also defines some auxiliary constants and data types.
! 
!*
module Grid
  use Constants

  implicit none
  
  private

  public :: Grid_t
  
  public :: NODE, FACE, EDGE, CELL, CENTER, CORNER, CELL_EARTH
  public :: XFACE, XEDGE, YFACE, YEDGE, ZFACE, ZEDGE

  !**
  ! Possible grid types for EMfield, storing
  ! the intention of use for types such as vectors
  ! and scalars.
  !*
  character (len = 4), parameter  :: NODE   = 'NODE'
  character (len = 4), parameter  :: FACE   = 'FACE'
  character (len = 4), parameter  :: EDGE   = 'EDGE'
  character (len = 6), parameter  :: CENTER = 'CELL'
  character (len = 6), parameter  :: CORNER = 'NODE'
  character (len = 4), parameter  :: CELL = 'CELL'  
  character (len = 10), parameter :: CELL_EARTH = 'CELL EARTH'

  !**
  ! Possible node types:
  !*  
  character(len = 5), parameter :: XFACE = 'XFACE'
  character(len = 5), parameter :: XEDGE = 'XEDGE'
  character(len = 5), parameter :: YFACE = 'YFACE'
  character(len = 5), parameter :: YEDGE = 'YEDGE'
  character(len = 5), parameter :: ZFACE = 'ZFACE'
  character(len = 5), parameter :: ZEDGE = 'ZEDGE'

  type, abstract :: Grid_t
     
     !**
     ! Grid geometry:
     ! regional or global grid geometry; important - used in EMfield
     ! This only refers to full sphere vs cube (region). The coordinates
     ! (either cartesian or spherical) are given by a global variable
     ! gridCoords, defined in GridCalc module.
     !*
     character(len = 80) :: geometry = REGION

     !**
     ! the origin of the model, by default set to zero
     !*
     real(kind = prec) :: ox = 0.0, oy = 0.0, oz = 0.0

     !**
     ! The rotation angle in degrees, by default set to zero
     !*
     real(kind = prec) :: rotDeg = 0.0 

     !**
     ! Grid Dimensions:
     ! nx is grid dimension (number of cells) in the x-direction
     ! ny is grid dimension (number of cells) in the y-direction
     ! nzEarth is number of earth layers used in the grid modeling
     ! nz is grid dimension (number of cells) in the z-direction:
     ! nz = nzAir + nzEarth
     !*
     integer :: nx = 0, ny = 0, nz = 0  ! Number of grid cells in x,y,z
                            ! directions

     integer :: nzAir = 0       ! Number of air layers
     integer :: nzEarth = 0     ! Number of earth layers

          !**
     ! Grid geometry:
     ! Dx, DDx are arrays of grid spacings in x-direction
     ! Dx denotes spacings betgween cell edges: dimension: Dx(nx)
     ! DDx denotes spacings between cell centers: dimension: DDx(nx+1)
     ! why dimensions: DDx(nx+1) (DDx(2) is distance between centers of cells
     ! 2 and 1)
     !*
     real(kind = prec), allocatable, dimension(:) :: dx, dy, dz
     real(kind = prec), allocatable, dimension(:) :: dxInv, dyInv, dzInv
     
     real(kind = prec), allocatable, dimension(:) :: delX, delY, delZ
     real(kind = prec), allocatable, dimension(:) :: delXInv, delYInv, delZInv

     !**
     ! Book-keeping on cumulative distances
     ! xEdge is the array for cumulative distance of the
     ! edge faces from the coordinate axis with dimensions nx+1
     ! xCenter is the array for cumulative distance of the
     ! edge center from the coordinate axis with dimensions nx
     ! yEdge, yCenter, zEdge, zCenter are analagous
     ! arrays for other directions.
     !*
     real(kind = prec), allocatable, dimension(:) :: xEdge
     real(kind = prec), allocatable, dimension(:) :: yEdge
     real(kind = prec), allocatable, dimension(:) :: zEdge
     real(kind = prec), allocatable, dimension(:) :: xCenter
     real(kind = prec), allocatable, dimension(:) :: yCenter
     real(kind = prec), allocatable, dimension(:) :: zCenter

     !**
     ! Total thickness of the air above
     !*
     real(kind = prec) :: zAirThick

     logical :: allocated = .False.

   contains
     private
     
     procedure(iface_NumberOfEdges)  , deferred, public :: NumberOfEdges
     procedure(iface_NumberOfFaces)  , deferred, public :: NumberOfFaces
     procedure(iface_NumberOfNodes)  , deferred, public :: NumberOfNodes
     procedure(iface_GridIndex)      , deferred, public :: GridIndex
     procedure(iface_VectorIndex)    , deferred, public :: VectorIndex
     procedure(iface_Limits)         , deferred, public :: Limits
     procedure(iface_IsAllocated)    , deferred, public :: IsAllocated

     procedure(iface_Copy_from), deferred, public :: Copy_from

     procedure(iface_SetCellSizes), deferred, public :: SetCellSizes
     
     procedure, public :: SetOrigin
     procedure, public :: GetOrigin

     procedure, public :: SetGridRotation
     procedure, public :: GetGridRotation     

     procedure, public :: SetGridGeometry
     procedure, public :: GetGridGeometry
          
  end type Grid_t


  abstract interface
     
     !**
     ! NumberOfEdges
     !
     !*
     subroutine iface_NumberOfEdges(self, nXedge, nYedge, nZedge)
       import :: Grid_t
       class(Grid_t), intent(in)  :: self
       integer      , intent(out) :: nXedge, nYedge, nZedge
     end subroutine iface_NumberOfEdges
     
     !**
     ! NumberOfFaces
     !
     !*
     subroutine iface_NumberOfFaces(self, nXface, nYface, nZface) 
       import :: Grid_t
       class(Grid_t), intent(in)  :: self
       integer      , intent(out) :: nXface, nYface, nZface
     end subroutine iface_NumberOfFaces
     
     !**
     ! NumberOfNodes
     !
     !*
     function iface_NumberOfNodes(self) result(n)
       import :: Grid_t
       class(Grid_t), intent(in)  :: self
       integer :: n
     end function iface_NumberOfNodes

     !**
     ! GridIndex
     !
     ! Based on matlab method of same name in class Grid_t
     ! IndVec is the index within the list of nodes of a fixed type
     ! e.g., among the list of y-Faces.   An offset needs to be
     ! added to get index in list of all faces (for example).
     !*
     subroutine iface_GridIndex(self, nodeType, indVec, i, j, k)
       import :: Grid_t
       class(Grid_t)         , intent(in)  :: self
       character(*)          , intent(in)  :: nodeType
       integer, dimension (:), intent(in)  :: indVec
       integer, dimension (:), intent(out) :: i, j, k
     end subroutine iface_GridIndex
     
     !**
     ! VectorIndex
     !
     ! Based on matlab method of same name in class Grid_t
     ! returned array IndVec gives numbering of nodes within
     ! the list for nodeType; need to add an offset for position
     ! in full list of all faces or edges (not nodes and cells).
     !*
     subroutine iface_VectorIndex(self, nodeType, i, j, k, indVec)
       import :: Grid_t
       class(Grid_t)         , intent(in)  :: self
       character(*)          , intent(in)  :: nodeType
       integer, dimension (:), intent(in)  :: i, j, k
       integer, dimension (:), intent(out) :: indVec
     end subroutine iface_VectorIndex
     
     !**
     ! Limits
     !
     !*
     subroutine iface_Limits(self, nodeType, nx, ny, nz)
       import :: Grid_t
       class(Grid_t) , intent(in)  :: self
       character(*)  , intent(in)  :: nodeType
       integer       , intent(out) :: nx, ny, nz
     end subroutine iface_Limits

     !**
     !
     !*
     function iface_IsAllocated(self) result(f)
       import :: Grid_t
       class(Grid_t), intent(in) :: self
       logical :: f       
     end function iface_IsAllocated
     
     subroutine iface_Copy_from(self, g)
       import :: Grid_t
       class(Grid_t), intent(inout) :: self
       class(Grid_t), intent(in)    :: g
     end subroutine iface_Copy_from

     subroutine iface_SetCellSizes(self, dx, dy, dz)
       import :: Grid_t, prec
       class(Grid_t), intent(inout) :: self
       real(kind = prec) , dimension(:), intent(in) :: dx, dy, dz
     end subroutine iface_SetCellSizes
     
  end interface
  
contains
  
  subroutine SetOrigin(self, ox, oy, oz)
    ! Arguments
    class(Grid_t)    , intent(inout) :: self
    real(kind = prec), intent(in) :: ox, oy, oz
    
    self%ox = ox
    self%oy = oy
    self%oz = oz
  end subroutine SetOrigin

  subroutine GetOrigin(self, ox, oy, oz)
    ! Arguments
    class(Grid_t)    , intent(in)  :: self
    real(kind = prec), intent(out) :: ox, oy, oz
    
    ox = self%ox
    oy = self%oy
    oz = self%oz
  end subroutine GetOrigin
    
  subroutine SetGridRotation(self, rotDeg)
    ! Arguments
    class(Grid_t)    , intent(inout) :: self
    real(kind = prec), intent(in)    :: rotDeg
    
    self%rotDeg = rotDeg
  end subroutine SetGridRotation

  function GetGridRotation(self) result(rotDeg)
    ! Arguments
    class(Grid_t), intent(in) :: self
    ! Local variables
    real(kind = prec) :: rotDeg
    
    rotDeg = self%rotDeg
  end function GetGridRotation
  
  subroutine SetGridGeometry(self, s)
    ! Arguments
    class(Grid_t), intent(inout) :: self
    character(*) , intent(in)    :: s

    self%geometry = s
  end subroutine SetGridGeometry

  function GetGridGeometry(self) result(s)
    ! Arguments
    class(Grid_t), intent(in) :: self
    ! Local variables
    character(80) :: s

    s = self%geometry 
  end function GetGridGeometry
  
end module Grid
