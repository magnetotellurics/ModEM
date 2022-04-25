!**
! Basic model parameter class --
! cell conductivities defined on numerical grid
! in Earth cells only -- as coded this is specific
!    to usual CSG.     Again, the abstract base class is too specific to 3D
!        (this is wrong!     ModelParamter base class should not be explicit at
!         all with regard to structure of the model parameter!!!!)
!        but I wil make this independent of any base class for now
!*
module ModelParameter2D
    use Constants
!    use ModelParameter     !     we use this to get access to "SigMap"
!        don't have a simple enough working version of ModelParameter abstract
!         class to "use"
    use Grid2D
    use Esoln2DTM                !     need this for PDEmapping
    use Grid1D                    !     need 1D Grid, ModelParameter for 1D "slice" routines
    use ModelParameter1D

    implicit none
    
    private

    public :: ModelParameter2D_t
    
    type :: ModelParameter2D_t
         ! This will be slightly modified from numerical grid --
         ! NzAir = 0 -- more generally this might be a
         ! completely different grid.
         real(kind=prec), allocatable, dimension(:,:) :: cellCond
         integer                         :: mKey(8) = 0
         character(:), allocatable :: paramType
         real(kind = prec)     :: airCond         = 1E-10
         class(Grid2D_t), allocatable :: grid 
         class(Grid2D_t), allocatable :: ParamGrid
         logical            :: isAllocated = .false.
         logical            :: zeroValued = .false.
         
     contains
         !     do not include all of the vector space operations
         procedure, public :: Length
         procedure, public :: Zeros
!         procedure, public :: Copy     !     no copy unless we need it
         procedure, public :: SetConductivity
         procedure, public :: PDEmapping
         procedure, public :: Slice1D => Slice1DModelParameter2D

    end type ModelParameter2D_t

    interface ModelParameter2D_t
         module procedure ModelParameter2D_ctor
    end interface ModelParameter2D_t
    
contains

    !**
    ! Class constructor
    !*
    function ModelParameter2D_ctor(grid) result(m)
        !        since this will (for now) be used for boundary and initial conditions
        !        the creator will not set the conductivity, just create arrays 
        type(Grid2D_t), intent(in) :: grid
        ! Local variables
        integer :: ny, nz, nzAir, status
        type(ModelParameter2D_t) :: m
        
        ny = grid%ny
        nz = grid%nz - grid%nzAir
        !     model grid
        m%grid = grid
        nzAir = 0
        
        !    Create ParamGrid using model grid parameters
        m%ParamGrid = Grid2D_t(ny,nzAir,nz,grid%dy, &
            grid%dz(grid%nzAir+1:grid%nz))

        allocate(m%cellCond(ny,nz),STAT = status)

        !m%mKey = datetime
        m%isAllocated = .true.
        
        !        set parameter array to zeros ...
        call m%Zeros()

    end function ModelParameter2D_ctor
    
    function Length(self) result(nParam)
        implicit none
        class(ModelParameter2D_t), intent(in) :: self
        integer :: nParam
        
        nParam = self%grid%Length()

    end function Length

    subroutine Zeros(self)
        implicit none
        ! Arguments
        class(ModelParameter2D_t), intent(inout) :: self
        
        self%cellCond = 0
        self%zeroValued = .true.

    end subroutine Zeros

    !**
    ! Copy rhs to self.
    !*
!    subroutine Copy(self, rhs)
        !     copies RHS to self, creating object if necessary

!        implicit none
        ! Arguments
!        class(ModelParameter2D_t), intent(inout) :: rhs
!        class(ModelParameter2D_t), intent(inout) :: self
        
!        if (.not. rhs%isAllocated) then
!             write(*, *) 'ERROR: ModelParameter2D (copy)'
!             write(*, *) '    RHS not allocated'
!             STOP
!        endif
!
!        if (.not. self%isAllocated) then
!             write(*, *) 'ERROR: ModelParameter2D (copy)'
!             write(*, *) '    self not allocated'
!             STOP
!        endif
!
!        self%cellCond    = rhs%cellCond
!        self%airCond     = rhs%airCond
!        self%paramType = rhs%paramType
!        self%mKey            = rhs%mKey
!
!    end subroutine Copy
    
    subroutine SetConductivity(self, CellCond, AirCond, paramType, mKey)
        !     sets conductivity array, assuming the ModelParameter2D_t
        !        object is already created

        implicit none
        ! Arguments
        class(ModelParameter2D_t), intent(inout)        :: self
        character(:), allocatable, intent(in)             :: paramType
        real(kind=prec), intent(in), dimension(:,:) :: CellCond
        real(kind=prec), intent(in)                                 :: AirCond
        integer, intent(in)                                                 :: mKey(8)
        !    local variables
        integer    :: nyz(2)
        
        if (.not. self%isAllocated) then
             write(*, *) 'ERROR: ModelParameter2D (SetConductivity)'
             write(*, *) '    input object not allocated'
             STOP
        endif

        nyz = shape( CellCond )
    !
        if((nyz(1) .ne. self%ParamGrid%ny) .or. &
            (nyz(2) .ne. self%ParamGrid%nz)) then
         write(*, *) 'nyz()                     ', nyz(1), nyz(2)
             
         write(*, *) 'self%ParamGrid%n', self%ParamGrid%ny, self%ParamGrid%nz
             write(*, *) 'ERROR: ModelParameter2D (SetConductivity)'
             write(*, *) '    input condutivity not consistent with grid'
             STOP
        endif

        self%cellCond    = cellCond
        self%AirCond = AirCond
        self%paramType = trim(paramType)
        self%mKey            = mKey

    end subroutine SetConductivity
    !
    !******
    !
    function PDEmapping(self) result(sigmaEdge)
        !    Output E is a real 1D array -- but easiest to do this using
        ! preliminary    mapping onto an Esoln2DTM object -- which has complex fields
        implicit none
        ! Arguments
        class(ModelParameter2D_t), intent(inout)    :: self
        real(kind=prec), allocatable,dimension(:) ::    sigmaEdge     !     this    allocated before input 
        !
        !     local variables
        type(Esoln2DTM_t) :: E2D        !     use this to store Ey, Ez on edges
        complex(kind=prec), allocatable, dimension(:) :: x     !    temporary complex storage
                                                     !     for edge conductivities
        integer :: ny, nz, nzEarth, nzAir, j, k, k1
        real(kind=prec)    :: w,AirCond
         
    if( trim( self%paramType ) .NE. LINEAR ) then
            print*,'ERROR SO FAR ONLY CODED FOR LINEAR paramTYpe:', trim( self%paramType )
            stop
        endif                                 
        E2D = Esoln2DTM_t(self%grid)         !     create 2D soln object form numerical grid

        ny = self%grid%ny
        nz = self%grid%nz
        nzAir = self%grid%nzAir
    nzEarth = nz-nzAir
        AirCond = self%AirCond
        !     set interior air edges to AirCond     -- Ey is comlpex and AirCond real
        E2D%Ey(:,2:NzAir) = AirCond    !     leave boundary set to zero
        E2D%Ez(2:Ny,1:NzAir) = AirCond     ! leave bondary set to zero
        !    air-earth interface
        w = self%grid%dz(NzAir+1)/(self%grid%dz(NzAir+1)+    &
                 self%grid%dz(NzAir))
        E2D%Ey(:,NzAir+1) = (AirCond*(1.0-w)) + self%CellCond(:,1)*w
        !     interior earth y-edges
        k1 = 0
        do k = NzAir+2,Nz
             k1 = k1+1
             w = self%grid%dz(k-1)/(self%grid%dz(k-1)+self%grid%dz(k))
             E2D%Ey(:,k) = w*self%CellCond(:,k1)+(1.0-w)*self%CellCond(:,k1+1)
        enddo
        !     interior earth z-edges
        do j = 2,Ny
             k1 = NzAir+1
             w = self%grid%dy(j-1)/(self%grid%dy(j-1)+self%grid%dy(j))
             E2D%Ez(j,k1:nz) = w*self%CellCond(j-1,1:nzEarth)+(1.0-w)* &
                 self%CellCond(j,1:nzEarth)
        enddo
        x = E2D%GetArray()
        sigmaEdge = real(x,kind=prec)

    end function    PDEmapping
    !
    !        **********
    !
    function Slice1DModelParameter2D(self,j) result(m1D)
        !     extracts slice corresponding to column j of model parameter
        implicit none
         ! Arguments
         class(ModelParameter2D_t), intent(in) :: self
         integer, intent(in)    :: j
         type(ModelParameter1D_t) ::    m1D 
         !     local variables
         type(Grid1D_t)     :: grid1
         real(kind=prec), allocatable, dimension(:)     :: CondSlice

         !     extract 1D grid
         grid1 = self%grid%Slice1D()
         !     create 1D model parameter
         m1D = ModelParameter1D_t(grid1)
         !     comnductivity slice
         allocate(CondSlice(grid1%nzEarth))
         CondSlice = self%CellCond(j,:)
         call m1d%SetConductivity(CondSlice, self%AirCond, &
             self%paramType, self%mKey)

    end function Slice1DModelParameter2D
    !**
    ! Add
    ! Adds two model parameters.
    !*
!    function Add(self, rhs) result(Eout)
!    end function Add
    
    !**
    ! Sub
    ! Subtracts two model parameters.
    !*
!    function Sub(self, rhs) result(Eout)
!    end function Sub
    
    !**
    ! Mult
    ! Multiply model parameter by real scalar c.
    !*
!    function Mult(c, self) result(Eout)
!    end function Mult
    
    !**
    ! DotProd
    ! Dot product of two model parameters.
    !*
!    function DotProd(self, rhs) result(r)
!    end function DotProd
    
    !**
    ! Zeros
    ! Zero model parameter
    !*

    !**
    ! SetType
    !*
    
end Module ModelParameter2D
