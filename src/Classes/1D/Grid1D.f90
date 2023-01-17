!
!> Stripped down grid for 1D -- most of what was in the base 3D grid is not needed
!
module Grid1D
    !
    use Constants
    !
    type :: Grid1D_t
        !
        !> Grid Dimensions:
        integer :: nz, nzAir, nzEarth        !> Number of earth layers
        !
        !>    Grid variables
        real( kind=prec ), allocatable, dimension(:) :: dz
        real( kind=prec ), allocatable, dimension(:) :: dzInv
        real( kind=prec ), allocatable, dimension(:) :: delZ
        real( kind=prec ), allocatable, dimension(:) :: delZInv
        !
        !> Book-keeping on cumulative distances
        real( kind=prec ), allocatable, dimension(:) :: zEdge
        real( kind=prec ), allocatable, dimension(:) :: zCenter
        !
        !> Total thickness of the air above
        real( kind=prec ) :: zAirThick
        !
        logical :: is_allocated
        !
        contains
            !
            final :: Grid1D_dtor
            !
            procedure, public :: alloc => allocateGrid1D
            procedure, public :: dealloc => deallocateGrid1D
            procedure, public :: setup => setupGrid1D
            procedure, public :: getDimensions => getDimensionsGrid1D
            procedure, public :: setCellSizes => setCellSizesGrid1D
            procedure, public :: getCellSizes => getCellSizesGrid1D

    end type Grid1D_t

    interface Grid1D_t
        module procedure Grid1D_ctor
    end interface Grid1D_t
    
contains
    !
    !> Class constructor for simple 1D FD grid
    !> Dz is cell dimensions for z direction
    !> Nza is number of air layers to allow (included in Dz)
    function Grid1D_ctor( nzAir, nzEarth, dz ) result( self )
        !
        integer, intent( in ) :: nzAir, nzEarth
        real( kind=prec ), dimension(:), intent( in ) :: dz
        !
        type( Grid1D_t ) :: self
        !
        !write( *, * ) "Constructor Grid1D"
        !
        self%nzAir = nzAir
        self%nzEarth = nzEarth
        !
        self%zAirThick = R_ZERO
        !
        self%is_allocated = .FALSE.
        !
        call self%alloc()
        !
        call self%setCellSizes( dz )
        call self%setup()
        !
    end function Grid1D_ctor
    !
    !> Deconstructor routine:
    !>     Calls the base routine dealloc().
    subroutine Grid1D_dtor( self )
        implicit none
        !
        type( Grid1D_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor Grid1D"
        !
        call self%dealloc()
        !
    end subroutine Grid1D_dtor
    !
    !> No subroutine briefing
    subroutine AllocateGrid1D( self )
        implicit none
        !
        class( Grid1D_t ), intent( inout ) :: self
        !
        self%nz = self%nzEarth + self%nzAir
        !
        allocate( self%dz( self%nz ) )
        !
        self%dz = R_ZERO
        !
        allocate( self%dzInv( self%nz ) )
        !
        self%dzInv = R_ZERO
        !
        allocate( self%delZ( self%nz + 1 ) )
        !
        self%delZ = R_ZERO
        !
        allocate( self%delZInv( self%nz + 1 ) )
        !
        self%delZInv = R_ZERO
        !
        allocate( self%zEdge( self%nz + 1 ) )
        !
        self%zEdge = R_ZERO
        !
        allocate( self%zCenter( self%nz ) )
        !
        self%zCenter = R_ZERO
        !
        self%is_allocated = .TRUE.
        !
    end subroutine AllocateGrid1D
    !
    !> No subroutine briefing
    subroutine deallocateGrid1D( self )
        implicit none
        !
        class( Grid1D_t ), intent( inout ) :: self
        !
        if( allocated(self%dz) ) deallocate(self%dz)
        if( allocated(self%dzInv) ) deallocate(self%dzInv)
        if( allocated(self%delZ) ) deallocate(self%delZ)
        if( allocated(self%delZInv) ) deallocate(self%delZInv)
        if( allocated(self%zEdge) ) deallocate(self%zEdge)
        if( allocated(self%zCenter) ) deallocate(self%zCenter)
        !
        self%nz = 0
        self%nzAir = 0
        self%nzEarth = 0
        self%zAirThick = R_ZERO
        !
        self%is_allocated = .FALSE.
        !
    end subroutine deallocateGrid1D
    !
    !> Procedure setupGrid1D
    !> Setup does calculations for grid geometry, which cannot be done
    !> until dz is set
    subroutine setupGrid1D( self )
        implicit none
        !
        class( Grid1D_t ), intent( inout ) :: self
        !
        integer :: iz, nzAir
        real( kind=prec ) :: zCum
        !
        self%dzInv = 1/self%dz
        !
        zCum = R_ZERO
        do iz = 1, self%nz
            zCum = zCum + self%dz(iz)
            self%zEdge(iz + 1) = zCum
        enddo
        !
        nzAir = self%nzAir
        self%zAirThick = self%zEdge(nzAir + 1)
        !
        !> Distance between center of the selfs
        self%delZ(1) = self%dz(1)
        do iz = 2, self%nz
            self%delZ(iz) = self%dz(iz - 1) + self%dz(iz)
        enddo
        self%delZ(self%nz + 1) = self%dz(self%nz)
        self%delZ = self%delZ/2.0
        !
        self%delZInv = 1/self%delZ
        !
        !> Cumulative distance between the centers, adjusted to model origin
        zCum = R_ZERO
        do iz = 1, self%nz
            zCum = zCum + self%delZ(iz)
            self%zCenter(iz) = zCum
        enddo
        !
        do iz = 1, self%nz
            self%zCenter(iz) = self%zCenter(iz) - self%zAirThick
            self%zEdge(iz) = self%zEdge(iz) - self%zAirThick
        enddo
        !
        self%zEdge(self%nz + 1) = self%zEdge(self%nz + 1) - self%zAirThick
        !
    end subroutine setupGrid1D
    !
    !> No subroutine briefing
    subroutine getDimensionsGrid1D( self, nz, nzAir )
        implicit none
        !
        class( Grid1D_t ), intent( in ) :: self
        integer, intent( out ) :: nz, nzAir
        !
        nz = self%nz
        nzAir = self%nzAir
        !
    end subroutine getDimensionsGrid1D
    !
    !> No subroutine briefing
    subroutine setCellSizesGrid1D( self, dz )
        implicit none
        !
        class( Grid1D_t ), intent( inout ) :: self
        real( kind=prec ), dimension(:), intent( in ) :: dz
        !
        if( .NOT. self%is_allocated ) then
            write( *, * ) "ERROR:Grid1D_t:SetCellSizes:"
            stop "    Grid not is_allocated."
        endif
        !
        !> Check dimensions
        if( (size(dz).NE.size(self%dz))) then
            write( *, * ) "ERROR:Grid1D_t:SetCellSizes:"
            stop "    Incompatible sizes for cell arrays."
        endif
        !
        self%dz = dz
        !
    end subroutine setCellSizesGrid1D
    !
    !> No subroutine briefing
    subroutine getCellSizesGrid1D( self, dz )
        implicit none
        !
        class( Grid1D_t ), intent( in ) :: self
        real( kind=prec ), intent( out ) :: dz(:)

        if( .NOT. self%is_allocated ) then
            write( *, * ) "ERROR:Grid1D_t:SetCellSizes:"
            stop "    Grid not is_allocated."
        endif

        !> Check dimensions
        if( ( size( dz ) .NE. size( self%dz ) ) ) then
            write( *, * ) "ERROR:Grid1D_t:SetCellSizes:"
            stop "    Incompatible sizes for cell arrays."
        endif
        !
        dz = self%dz
        !
    end subroutine getCellSizesGrid1D
    !
end module Grid1D
