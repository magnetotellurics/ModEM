!**
! Stripped down grid for 1D -- most of what was in the base 3D grid is not needed
!
module Grid1D
  use Constants
  
  implicit none
  
  private
  
  !    I am going to assume that Air Layers are provided -- not created
  !    with a method in this object (as for 2D)
  public  :: Grid1D_t 

  type :: Grid1D_t
     !**
     ! Grid Dimensions:
     integer :: nz = 0  ! Number of grid cells in z direction
     integer :: nzAir = 0       ! Number of air layers
     integer :: nzEarth = 0     ! Number of earth layers

    !  Grid variables
     real(kind = prec), allocatable, dimension(:) :: dz
     real(kind = prec), allocatable, dimension(:) :: dzInv
     real(kind = prec), allocatable, dimension(:) :: delZ
     real(kind = prec), allocatable, dimension(:) :: delZInv
     
     ! Book-keeping on cumulative distances
     real(kind = prec), allocatable, dimension(:) :: zEdge
     real(kind = prec), allocatable, dimension(:) :: zCenter

     !**
     ! Total thickness of the air above
     !*
     real(kind = prec) :: zAirThick

     logical :: allocated = .False.

   contains

     !   destructor
     final :: Grid1D_dtor
     
     procedure, public :: IsAllocated
     procedure, public :: Create
     procedure, public :: Allocate
     procedure, public :: deallocate
     procedure, public :: Setup
     procedure, public :: GetDimensions
     procedure, public :: SetCellSizes
     procedure, public :: GetCellSizes

  end type Grid1D_t

  interface Grid1D_t
     module procedure Grid1D_ctor
  end interface Grid1D_t
  
contains

  !**
  ! Class constructor for simple 1D FD grid
  ! Dz is cell dimensions for z direction
  ! Nza is number of air layers to allow (included in Dz)
  !*
  function Grid1D_ctor(nzAir, nzEarth, dz) result(grid)
    ! Arguments
    integer, intent(in) :: nzAir, nzEarth
    real(kind = prec) , dimension(:), intent(in) :: dz
    ! Local variables
    type(Grid1D_t) :: grid
    
    call grid%Create(nzAir, nzEarth)
    call grid%SetCellSizes(dz)
    call grid%Setup()
    
  end function Grid1D_ctor

  subroutine Create(self, nzAir, nzEarth)
    ! Arguments
    class(Grid1D_t), intent(inout) :: self
    integer           , intent(in) :: nzAir, nzEarth
    ! Local variables
    integer :: nz
    
    nz = nzEarth + nzAir
    
    self%nzAir = nzAir
    self%nzEarth = nzEarth
    
    self%nz = nz
    
    call self%Allocate()

  end subroutine Create
  
  subroutine Allocate(self)
    ! Arguments
    class(Grid1D_t), intent(inout) :: self
    ! Local variables
    integer :: nz
    
    if (self%Allocated) call self%deallocate()

    nz = self%nz
    allocate(self%dz(nz))
    allocate(self%dzInv(nz))
    allocate(self%delZ(nz + 1))
    allocate(self%delZInv(nz + 1))
    allocate(self%zEdge(nz + 1))
    allocate(self%zCenter(nz))
    
    self%allocated = .true.
    
  end subroutine Allocate

  subroutine deallocate(self)
    ! Arguments
    class(Grid1D_t), intent(inout) :: self
    
    if (.not.self%allocated) return

    deallocate(self%dz)
    deallocate(self%dzInv)
    deallocate(self%delZ)
    deallocate(self%delZInv)
    deallocate(self%zEdge)
    deallocate(self%zCenter)
    
    self%nz = 0
    self%nzAir = 0
    self%nzEarth = 0
    self%zAirThick = R_ZERO
  end subroutine deallocate
    
  subroutine Grid1D_dtor(self)
    type(Grid1D_t)  :: self
    call self%deallocate()
  end subroutine Grid1D_dtor

  !**
  ! Setup does calculations for grid geometry, which cannot be done
  ! until dz is set
  !
  !*
  subroutine Setup(self)
    implicit none
    ! Arguments
    class(Grid1D_t), intent(inout)        :: self
    ! Local variables
    integer :: iz
    integer :: nzAir
    real(kind = prec) :: zCum
    
    self%dzInv = 1/self%dz

    zCum = 0.0
    do iz = 1, self%nz
       zCum = zCum + self%dz(iz)
       self%zEdge(iz + 1) = zCum
    end do

    nzAir = self%nzAir
    self%zAirThick = self%zEdge(nzAir + 1)

    ! Distance between center of the selfs
    self%delZ(1) = self%dz(1)
    do iz = 2, self%nz
       self%delZ(iz) = self%dz(iz - 1) + self%dz(iz)
    end do
    self%delZ(self%nz + 1) = self%dz(self%nz)
    self%delZ = self%delZ/2.0
    
    self%delZInv = 1/self%delZ

    ! Cumulative distance between the centers, adjusted to model origin
    zCum = 0.0
    do iz = 1, self%nz
       zCum = zCum + self%delZ(iz)
       self%zCenter(iz) = zCum
    end do
    
    do iz = 1, self%nz
       self%zCenter(iz) = self%zCenter(iz) - self%zAirThick
       self%zEdge(iz) = self%zEdge(iz) - self%zAirThick
    end do
    self%zEdge(self%nz + 1) = self%zEdge(self%nz + 1) - &
         self%zAirThick
    
  end subroutine Setup

  subroutine GetDimensions(self, nz, nzAir)
    ! Arguments
    class(Grid1D_t), intent(in)  :: self
    integer           , intent(out) :: nz, nzAir

    nz = self%nz
    nzAir = self%nzAir
  end subroutine GetDimensions

  subroutine SetCellSizes(self, dz)
    ! Arguments
    class(Grid1D_t), intent(inout) :: self
    real(kind = prec) , dimension(:), intent(in) :: dz

    if (.not.self%IsAllocated()) then
       write(*, *) 'ERROR:Grid1D_t:SetCellSizes:'
       write(*, *) '  Grid not allocated.'

       STOP
    end if

    ! Check dimensions
    if ( (size(dz).ne.size(self%dz))) then
       write(*, *) 'ERROR:Grid1D_t:SetCellSizes:'
       write(*, *) '  Incompatible sizes for cell arrays.'

       STOP
    end if

    self%dz = dz

  end subroutine SetCellSizes
  
  subroutine GetCellSizes(self, dz)
    ! Arguments
    class(Grid1D_t), intent(in) :: self
    real(kind = prec) , intent(out) :: dz(:)

    if (.not.self%IsAllocated()) then
       write(*, *) 'ERROR:Grid1D_t:SetCellSizes:'
       write(*, *) '  Grid not allocated.'
       STOP
    end if

    ! Check dimensions
    if ( (size(dz).ne.size(self%dz))) then
       write(*, *) 'ERROR:Grid1D_t:SetCellSizes:'
       write(*, *) '  Incompatible sizes for cell arrays.'
       STOP
    end if
    
    dz = self%dz    
    
  end subroutine GetCellSizes

  function IsAllocated(self) result(f)
    ! Arguments
    class(Grid1D_t), intent(in) :: self
    ! Local variables
    logical :: f

    f = self%allocated
  end function IsAllocated

end module Grid1D
