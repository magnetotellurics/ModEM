!
!> Basic model parameter class --
!> cell conductivities defined on numerical grid
!> in Earth cells only -- as coded this is specific
!>    to usual CSG.     Again, the abstract base class is too specific to 3D
!>        (this is wrong!>     ModelParamter base class should not be explicit at
!>         all with regard to structure of the model parameter!!!!)
!>        but I wil make this independent of any base class for now
!
module ModelParameter2D
    !
    use Constants
    !
    use Grid2D
    use Esoln2DTM    !> need Esoln2DTM for PDEmapping
    use Grid1D        !> Grid1D, ModelParameter1D for 1D "slice" routines
    use ModelParameter1D
    !
    type :: ModelParameter2D_t
         !
         real( kind=prec ), allocatable, dimension(:,:) :: cellCond
         !
         integer :: mKey(8)
         !
         character(:), allocatable :: paramType
         !
         real( kind=prec ) :: airCond
         !
         type(Grid2D_t) :: grid, ParamGrid
         !
         logical :: is_allocated, zeroValued
         !
         contains
            !
            procedure, public :: length => lengthModelParameter2D
            procedure, public :: zeros => zerosParameter2D
            procedure, public :: setConductivity => setConductivityModelParameter2D
            procedure, public :: PDEmapping => PDEmappingParameter2D
            procedure, public :: slice1D => slice1DModelParameter2D
            !
    end type ModelParameter2D_t
    !
    interface ModelParameter2D_t
         module procedure ModelParameter2D_ctor
    end interface ModelParameter2D_t
    !
contains
    !
    !> No subroutine briefing
    !
    function ModelParameter2D_ctor(grid) result( self )
        implicit none
        !
        type( Grid2D_t ), intent( in ) :: grid
        !
        integer :: ny, nz, nzAir
        type( ModelParameter2D_t ) :: self
        !
        call date_and_time( values=self%mKey )
        !
        self%paramType = ""
        self%airCond = 1E-10
        self%is_allocated = .FALSE.
        self%zeroValued = .FALSE.
        !
        ny = grid%ny
        nz = grid%nz - grid%nzAir
        !>     model grid
        self%grid = grid
        nzAir = 0
        
        !>    Create ParamGrid using model grid parameters
        self%ParamGrid = Grid2D_t(ny,nzAir,nz,grid%dy, &
            grid%dz(grid%nzAir+1:grid%nz))

        allocate(self%cellCond(ny,nz))
        self%is_allocated = .TRUE.
        
        !>        set parameter array to zeros ...
        call self%zeros

    end function ModelParameter2D_ctor
    !
    !> No subroutine briefing
    !
    function lengthModelParameter2D( self ) result(nParam)
        implicit none
        class(ModelParameter2D_t), intent( in ) :: self
        integer :: nParam
        
        nParam = self%grid%Length()

    end function lengthModelParameter2D
    !
    !> No subroutine briefing
    subroutine zerosParameter2D( self )
        implicit none
        !
        class(ModelParameter2D_t), intent(inout) :: self
        
        self%cellCond = 0
        self%zeroValued = .TRUE.

    end subroutine zerosParameter2D
    !
    !> Procedure setConductivityModelParameter2D
    subroutine setConductivityModelParameter2D(self, CellCond, AirCond, paramType, mKey)
        !>     sets conductivity array, assuming the ModelParameter2D_t
        !>        object is already created
        !
        implicit none
        !
        class(ModelParameter2D_t), intent(inout) :: self
        character(:), allocatable, intent( in ) :: paramType
        real( kind=prec ), intent( in ), dimension(:,:) :: CellCond
        real( kind=prec ), intent( in ) :: AirCond
        integer, intent( in ) :: mKey(8)
        !>    local variables
        integer :: nyz(2)
        
        if(.NOT. self%is_allocated) then
             write(*, *) 'ERROR: ModelParameter2D (SetConductivity)'
             write(*, *) '    input object not allocated'
             STOP
        endif

        nyz = shape( CellCond )
    !
        if((nyz(1) .NE. self%ParamGrid%ny) .OR. &
            (nyz(2) .NE. self%ParamGrid%nz)) then
         write(*, *) 'nyz()                     ', nyz(1), nyz(2)
             
         write(*, *) 'self%ParamGrid%n', self%ParamGrid%ny, self%ParamGrid%nz
             write(*, *) 'ERROR: ModelParameter2D (SetConductivity)'
             write(*, *) '    input condutivity not consistent with grid'
             STOP
        endif

        self%cellCond = cellCond
        self%AirCond = AirCond
        self%paramType = trim(paramType)
        self%mKey = mKey

    end subroutine setConductivityModelParameter2D
    !
    !> No subroutine briefing
    !
    function PDEmappingParameter2D( self ) result(sigmaEdge)
        !>    Output E is a real 1D array -- but easiest to do this using
        !> preliminary    mapping onto an Esoln2DTM object -- which has complex fields
        implicit none
        !
        class(ModelParameter2D_t), intent(inout) :: self
        real( kind=prec ), allocatable,dimension(:) ::    sigmaEdge     !>     this    allocated before input 
        !
        !>     local variables
        type(Esoln2DTM_t) :: E2D        !>     use this to store Ey, Ez on edges
        complex( kind=prec ), allocatable, dimension(:) :: x     !>    temporary complex storage
                                                     !>     for edge conductivities
        integer :: ny, nz, nzEarth, nzAir, j, k, k1
        real( kind=prec ) :: w,AirCond
         
        if( trim( self%paramType ) .NE. LINEAR ) then
            print*,'ERROR SO FAR ONLY CODED FOR LINEAR paramTYpe:', trim( self%paramType )
            stop
        endif                                 
        E2D = Esoln2DTM_t(self%grid)         !>     create 2D soln object form numerical grid

        ny = self%grid%ny
        nz = self%grid%nz
        nzAir = self%grid%nzAir
        nzEarth = nz-nzAir
        AirCond = self%AirCond
        !>     set interior air edges to AirCond     -- Ey is comlpex and AirCond real
        E2D%Ey(:,2:NzAir) = AirCond    !>     leave boundary set to zero
        E2D%Ez(2:Ny,1:NzAir) = AirCond     !> leave bondary set to zero
        !>    air-earth interface
        w = self%grid%dz(NzAir+1)/(self%grid%dz(NzAir+1)+    &
                 self%grid%dz(NzAir))
        E2D%Ey(:,NzAir+1) = (AirCond*(1.0-w)) + self%CellCond(:,1)*w
        !>     interior earth y-edges
        k1 = 0
        do k = NzAir+2,Nz
             k1 = k1+1
             w = self%grid%dz(k-1)/(self%grid%dz(k-1)+self%grid%dz(k))
             E2D%Ey(:,k) = w*self%CellCond(:,k1)+(1.0-w)*self%CellCond(:,k1+1)
        enddo
        !>     interior earth z-edges
        do j = 2,Ny
             k1 = NzAir+1
             w = self%grid%dy(j-1)/(self%grid%dy(j-1)+self%grid%dy(j))
             E2D%Ez(j,k1:nz) = w*self%CellCond(j-1,1:nzEarth)+(1.0-w)* &
                 self%CellCond(j,1:nzEarth)
        enddo
        x = E2D%GetArray()
        sigmaEdge = real(x,kind=prec)

    end function PDEmappingParameter2D
    !
    !> No subroutine briefing
    !
    function slice1DModelParameter2D(self,j) result(m1D)
        !>     extracts slice corresponding to column j of model parameter
        implicit none
         !
         class(ModelParameter2D_t), intent( in ) :: self
         integer, intent( in ) :: j
         type(ModelParameter1D_t) ::    m1D 
         !>     local variables
         type(Grid1D_t) :: grid1
         real( kind=prec ), allocatable, dimension(:) :: CondSlice

         !>     extract 1D grid
         grid1 = self%grid%Slice1D()
         !>     create 1D model parameter
         m1D = ModelParameter1D_t(grid1)
         !>     comnductivity slice
         allocate(CondSlice(grid1%nzEarth))
         CondSlice = self%CellCond(j,:)
         call m1d%SetConductivity(CondSlice, self%AirCond, &
             self%paramType, self%mKey)

    end function slice1DModelParameter2D
    !
end Module ModelParameter2D
