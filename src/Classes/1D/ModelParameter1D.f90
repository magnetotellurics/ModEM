!**
! Basic model parameter class for 1D
! cell conductivities defined on numerical grid
!*
module ModelParameter1D
  use Constants
  use Grid1D

  implicit none
  
  private

  public :: ModelParameter1D_t, SetLayeredModel
  
  type :: ModelParameter1D_t
     ! This will be slighty modified from numerical grid --
     ! NzAir = 0 -- more generally this might be a
     ! completely different grid.
     real(kind=prec), allocatable, dimension(:) :: cellCond
     integer             :: mKey(8)
     character(len = 80) :: paramType   = ''
     real(kind = prec)   :: airCond     = 1E-10
     !   model grid
     type(Grid1D_t) :: grid
     !   parameter grid
     class(Grid1D_t), allocatable :: paramGrid
     logical      :: isAllocated = .false.
     logical      :: zeroValued = .false.
     
   contains
     !   do not include all of the vector space operations
     procedure, public :: Zeros
     procedure, public :: SetLayeredModel
     procedure, public :: SetConductivity
!     procedure, public :: Copy   !  omit copy for now`

  end type ModelParameter1D_t

  interface ModelParameter1D_t
     module procedure ModelParameter1D_ctor
  end interface ModelParameter1D_t
  
contains

  !**
  ! Class constructor
  !*
  function ModelParameter1D_ctor(grid) result(m)
    !    since this will (for now) be used for boundary and initial conditions
    !    the creator will not set the conductivity, just create arrays 
    type(Grid1D_t),  intent(in) :: grid
    ! Local variables
    integer :: nz, nzAir, status
    type(ModelParameter1D_t) :: m
    
    !   set pointer to external grid
    m%grid = grid

    nz = grid%nz - grid%nzAir
    nzAir = 0
    m%paramGrid = Grid1D_t( nzAir, nz, grid%dz(grid%nzAir+1:grid%nz) )
    
    allocate(m%cellCond(nz),STAT = status)
    
    !m%mKey = datetime
    m%isAllocated = .true.
    
    !    set parameter array to zeros ...
    call m%Zeros()

  end function ModelParameter1D_ctor
  
  subroutine Zeros(self)
    implicit none
    ! Arguments
    class(ModelParameter1D_t), intent(inout) :: self
    
    self%cellCond = 0
    self%zeroValued = .true.

  end subroutine Zeros
  !
  !**************
  !
  subroutine SetLayeredModel(self, nLayers, h, sigma)
    !    model parameter creator sets the parameter grid to the model
    !      grid, with Nza set to zero
    !    This allows setting for a layered model with nLayers, of
    !     thickness h, and conductivity (always linear) sigma 

    implicit none
    ! Arguments
    class(ModelParameter1D_t), intent(inout) :: self
    integer, intent(in) :: nLayers
    real(kind=prec),intent(in), dimension(nLayers) :: h
    real(kind=prec),intent(in), dimension(nLayers) :: sigma

    !  local variables
    integer   :: NzAir
    
    NzAir = 0
    self%ParamGrid = Grid1D_t(NzAir, nLayers, h)
    self%CellCond = sigma
    self%paramTYpe = 'LINEAR'

  end subroutine SetLayeredModel
  !
  !*********
  !
  subroutine SetConductivity(self, CellCond, AirCond, paramType, mKey)
    !   sets conductivity array, assuming the ModelParameter1D_t
    !    object is already created

    implicit none
    ! Arguments
    class(ModelParameter1D_t), intent(inout) :: self
    character(len=80),intent(in) :: paramType
    real(kind=prec),intent(in), dimension(:) :: CellCond
    real(kind=prec), intent(in) :: AirCond
    integer, intent(in)   :: mKey(8)
    !  local variables
    integer :: nz
    
    if (.not. self%isAllocated) then
       write(*, *) 'ERROR: ModelParameter1D (SetConductivity)'
       write(*, *) '  input object not allocated'
       STOP
    endif
    
    nz = size(CellCond)
	!
    if(nz .ne. self%ParamGrid%nz) then
       write(*, *) 'ERROR: ModelParameter1D (SetConductivity)'
       write(*, *) '  input condutivity not consistent with grid'
       STOP
    endif

    self%cellCond  = cellCond
    self%AirCond = AirCond
    self%paramType = trim(paramType)
    self%mKey      = mKey

  end subroutine SetConductivity
  !**
  ! Add
  ! Adds two model parameters.
  !*
!  function Add(self, rhs) result(Eout)
!  end function Add
  
  !**
  ! Sub
  ! Subtracts two model parameters.
  !*
!  function Sub(self, rhs) result(Eout)
!  end function Sub
  
  !**
  ! Mult
  ! Multiply model parameter by real scalar c.
  !*
!  function Mult(c, self) result(Eout)
!  end function Mult
  
  !**
  ! DotProd
  ! Dot product of two model parameters.
  !*
!  function DotProd(self, rhs) result(r)
!  end function DotProd
  
  !**
  ! Zeros
  ! Zero model parameter
  !*

  !**
  ! SetType
  !*
  
end Module ModelParameter1D
