!
!> Basic model parameter class for 1D
!> cell conductivities defined on numerical grid
!
module ModelParameter1D
    !
    use Constants
    use Grid1D
    !
    type :: ModelParameter1D_t
        !
        real( kind=prec ), allocatable, dimension(:) :: cellCond
        !
        integer :: mKey(8)
        !
        character( len=80 ) :: paramType
        !
        real( kind=prec ) :: airCond
        !
        type( Grid1D_t ) :: grid, paramGrid
        !
        logical :: is_allocated, zero_valued
        !
        contains
            !
            final :: ModelParameter1D_dtor
            !
            procedure, public :: zeros => zerosModelParameter1D
            procedure, public :: setLayeredModel => setLayeredModelModelParameter1D
            procedure, public :: setConductivity => setConductivityModelParameter1D
            !
    end type ModelParameter1D_t
    !
    interface ModelParameter1D_t
        module procedure ModelParameter1D_ctor
    end interface ModelParameter1D_t
    !
contains
    !
    !> No subroutine briefing
    !
    function ModelParameter1D_ctor( grid ) result( self )
        implicit none
        !
        type( Grid1D_t ), intent( in ) :: grid
        type( ModelParameter1D_t ) :: self
        !
        integer :: nz, nzAir
        !
        !write( *, * ) "Constructor ModelParameter1D_t"
        !
        self%paramType = ""
        self%airCond = 1E-10
        self%is_allocated = .FALSE.
        self%zero_valued = .FALSE.
        !
        self%grid = grid
        !
        nz = grid%nz - grid%nzAir
        nzAir = 0
        self%paramGrid = Grid1D_t( nzAir, nz, grid%dz(grid%nzAir+1:grid%nz) )
        !
        allocate( self%cellCond( nz ) )
        !
        call date_and_time( values=self%mKey )
        !
        self%is_allocated = .TRUE.
        !
        call self%zeros
        !
    end function ModelParameter1D_ctor
    !
    !> Deconstructor routine:
    !>     Calls the base routine baseDealloc().
    !>     Deallocates inherent properties of this class.
    subroutine ModelParameter1D_dtor( self )
        implicit none
        !
        type( ModelParameter1D_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor ModelParameter1D_t"
        !
        if(allocated(self%cellCond)) deallocate( self%cellCond )
        !
    end subroutine ModelParameter1D_dtor
    !
    !> No subroutine briefing
    subroutine zerosModelParameter1D( self )
        implicit none
        !
        class( ModelParameter1D_t ), intent( inout ) :: self
        !
        self%cellCond = 0
        self%zero_valued = .TRUE.
        !
    end subroutine zerosModelParameter1D
    !
    !> No subroutine briefing
    subroutine setLayeredModelModelParameter1D( self, nLayers, h, sigma )
        implicit none
        !
        class( ModelParameter1D_t ), intent( inout ) :: self
        integer, intent( in ) :: nLayers
        real( kind=prec ),intent( in ), dimension( nLayers ) :: h
        real( kind=prec ),intent( in ), dimension( nLayers ) :: sigma
        !
        integer :: NzAir
        !
        NzAir = 0
        self%ParamGrid = Grid1D_t( NzAir, nLayers, h )
        self%CellCond = sigma
        self%paramTYpe = "LINEAR"
        !
    end subroutine setLayeredModelModelParameter1D
    !
    !> No subroutine briefing
    subroutine setConductivityModelParameter1D( self, CellCond, AirCond, paramType, mKey )
        implicit none
        !
        class( ModelParameter1D_t ), intent( inout ) :: self
        real( kind=prec ),intent( in ), dimension(:) :: CellCond
        real( kind=prec ), intent( in ) :: AirCond
        character(:), allocatable, intent( in ) :: paramType
        integer, intent( in ) :: mKey(8)
        !
        integer :: nz
        !
        if( .NOT. self%is_allocated) then
            write( *, * ) "     "//achar(27)//"[31m# Error:"//achar(27)//"[0m ModelParameter1D (SetConductivity)"
            stop "     input object not allocated"
        endif
        !
        nz = size(CellCond)
        !
        if( nz .NE. self%ParamGrid%nz ) then
            write( *, * ) "     "//achar(27)//"[31m# Error:"//achar(27)//"[0m ModelParameter1D (SetConductivity)"
            stop "     input condutivity not consistent with grid"
        endif
        !
        self%cellCond = cellCond
        self%AirCond = AirCond
        self%paramType = trim(paramType)
        self%mKey = mKey
        !
    end subroutine setConductivityModelParameter1D
    !
end Module ModelParameter1D
