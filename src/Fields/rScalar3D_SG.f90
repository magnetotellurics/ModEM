!
!> Derived class to define a rScalar3D_SG 
!
module rScalar3D_SG
    !
    use Scalar
    use Grid3D_SG
    !
    type, extends( Scalar_t ) :: rScalar3D_SG_t
        !
        real( kind=prec ), allocatable, dimension(:, :, :) :: v
        !
        real( kind=prec ), allocatable, dimension(:) :: sv
        !
        contains
            !
            !> Destructor
            final :: rScalar3D_SG_dtor
            !
            !> Boundary operations
            procedure, public :: setAllBoundary => setAllBoundaryRScalar3D_SG
            procedure, public :: setOneBoundary => setOneBoundaryRScalar3D_SG
            procedure, public :: setAllInterior => setAllInteriorRScalar3D_SG
            procedure, public :: intBdryIndices => intBdryIndicesRScalar3D_SG
            !
            !> Dimensioning operations
            procedure, public :: length => lengthRScalar3D_SG
            procedure, public :: setVecComponents => setVecComponentsRScalar3D_SG
            !
            !> Arithmetic/algebraic unary operations
            procedure, public :: zeros => zerosRScalar3D_SG
            procedure, public :: sumEdges => sumEdgesRScalar3D_SG
            procedure, public :: avgCells => avgCellsRScalar3D_SG
            procedure, public :: conjugate => conjugateRScalar3D_SG
            !
            !> Arithmetic/algebraic binary operations
            procedure, public :: add => addRScalar3D_SG
            !
            procedure, public :: linComb => linCombRScalar3D_SG
            !
            procedure, public :: subValue => subValueRScalar3D_SG
            procedure, public :: subField => subFieldRScalar3D_SG
            !
            procedure, public :: multByReal => multByRealRScalar3D_SG
            procedure, public :: multByComplex => multByComplexRScalar3D_SG
            procedure, public :: multByField => multByFieldRScalar3D_SG
            !
            procedure, public :: multAdd => multAddRScalar3D_SG
            !
            procedure, public :: dotProd => dotProdRScalar3D_SG
            !
            procedure, public :: divByField => divByFieldRScalar3D_SG
            procedure, public :: divByValue => divByValueRScalar3D_SG
            !
            !> Miscellaneous
            procedure, public :: getReal => getRealRScalar3D_SG
            procedure, public :: getArray => getArrayRScalar3D_SG
            procedure, public :: setArray => setArrayRScalar3D_SG
            procedure, public :: switchStoreState => switchStoreStateRScalar3D_SG
            procedure, public :: copyFrom => copyFromRScalar3D_SG
            !
            !> I/O operations
            procedure, public :: read => readRScalar3D_SG
            procedure, public :: write => writeRScalar3D_SG
            procedure, public :: print => printRScalar3D_SG
            !
    end type rScalar3D_SG_t
    !
    interface rScalar3D_SG_t
        module procedure rScalar3D_SG_ctor
    end interface rScalar3D_SG_t
    !
contains
    !
    !> No subroutine briefing
    !
    function rScalar3D_SG_ctor( grid, grid_type ) result ( self )
        implicit none
        !
        class( Grid_t ), target, intent( in ) :: grid
        character( len=4 ), intent( in ) :: grid_type
        !
        type( rScalar3D_SG_t ) :: self
        !
        integer :: nx, ny, nz, nzAir, nz_earth
        integer :: status
        !
        !write( *, * ) "Constructor rScalar3D_SG"
        !
        call self%init
        !
        self%grid => grid
        self%grid_type = grid_type
        !
        !> Grid dimensions
        call grid%getDimensions( nx, ny, nz, nzAir )
        nz_earth = nz - nzAir
        !
        self%nx = nx
        self%ny = ny
        self%nz = nz
        !
        !> allocate memory for x,y,z ;
        !> self%allocated will be true if all allocations succeed
        self%is_allocated = .TRUE.
        !
        if( grid_type == CORNER ) then
             allocate(self%v(nx + 1, ny + 1, nz + 1), STAT = status)    
             self%NdV = (/self%nx + 1, self%ny + 1, self%nz + 1/)
        !
        else if( grid_type == CENTER ) then
             allocate(self%v(nx, ny, nz), STAT = status) 
             self%NdV = (/self%nx, self%ny, self%nz/)
        !
        else if( grid_type == CELL_EARTH) then
             self%nz = nz_earth
             allocate(self%v(nx, ny, nz_earth), STAT = status)
             self%NdV = (/nx, ny, nz_earth/)
        !
        else
            write( *, * ) "Error: rScalar3D_SG_ctor > unrecognized grid type: [", grid_type, "]"
            stop
        endif
        !
        self%is_allocated = self%is_allocated .AND. (status .EQ. 0)
        !
        if( self%is_allocated ) then
             self%v = R_ZERO
        else
             stop "Error: rScalar3D_SG_ctor > Unable to allocate rScalar - invalid grid supplied"
        endif
        !
        self%Nxyz = product( self%NdV )
        !
    end function rScalar3D_SG_ctor
    !
    !> No subroutine briefing
    !
    subroutine rScalar3D_SG_dtor( self )
        implicit none
        !
        type( rScalar3D_SG_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor rScalar3D_SG"
        !
        if( allocated( self%v ) ) deallocate( self%v )
        if( allocated( self%sv ) ) deallocate( self%sv )
        !
        self%nx = 0
        self%ny = 0
        self%nz = 0
        !
        self%grid_type = ""
        self%is_allocated = .FALSE.
        !
    end subroutine rScalar3D_SG_dtor
    !
    !> No subroutine briefing
    !
    subroutine setAllBoundaryRScalar3D_SG( self, cvalue )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        if( self%store_state /= compound ) then
             call self%switchStoreState
        endif
        !
        select case( self%grid_type )
            case( CORNER ) 
                 self%v((/1, self%NdV(1)/), :, :) = real( cvalue, kind=prec )
                 self%v(:, (/1, self%NdV(2)/), :) = real( cvalue, kind=prec )
                 self%v(:, :, (/1, self%NdV(3)/)) = real( cvalue, kind=prec )
                 !
            case default
                 stop "Error: setAllBoundaryRScalar3D_SG > Grid type not recognized."
        end select
        !
    end subroutine setAllBoundaryRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setOneBoundaryRScalar3D_SG( self, bdry, cvalue, int_only )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        character(*), intent( in ) :: bdry
        complex( kind=prec ), intent( in ) :: cvalue
        logical, intent( in ), optional :: int_only
        !
        logical :: int_only_p
        !
        if( self%store_state /= compound ) then
             call self%switchStoreState
        endif
        !
        if( .NOT. present (int_only)) then
             int_only_p = .FALSE.
        else 
             int_only_p = int_only
        endif
        !
        select case( self%grid_type )
        case(CORNER)
             if( int_only_p) then
                select case(bdry)
                case("x1")
                     self%v(1, 2:self%NdV(2)-1, 2:self%NdV(3)-1) = real( cvalue, kind=prec ) 
                case("x2")
                     self%v(self%NdV(1), 2:self%NdV(2)-1, 2:self%NdV(3)-1) = real( cvalue, kind=prec )
                case("y1")
                     self%v(2:self%NdV(1)-1, 1, 2:self%NdV(3)-1) = real( cvalue, kind=prec )
                case("y2")
                     self%v(2:self%NdV(1)-1, self%NdV(2), 2:self%NdV(3)-1) = real( cvalue, kind=prec )
                case("z1")
                     self%v(2:self%NdV(1)-1, 2:self%NdV(2)-1, 1) = real( cvalue, kind=prec )
                case("z2")
                     self%v(2:self%NdV(1)-1, 2:self%NdV(2)-1, self%NdV(3)) = real( cvalue, kind=prec )
                end select
             else
                select case(bdry)
                case("x1")
                     self%v(1, :, :) = real( cvalue, kind=prec )
                case("x2")
                     self%v(self%NdV(1), :, :) = real( cvalue, kind=prec )
                case("y1")
                     self%v(:, 1, :) = real( cvalue, kind=prec )
                case("y2")
                     self%v(:, self%NdV(2), :) = real( cvalue, kind=prec )
                case("z1")
                     self%v(:, :, 1) = real( cvalue, kind=prec )
                case("z2")
                     self%v(:, :, self%NdV(3)) = real( cvalue, kind=prec )
                end select
             endif
             !
        case(FACE)
             select case(bdry)
                 case("x1")
                    self%v(1, :, :) = real( cvalue, kind=prec )
                 case("x2")
                    self%v(self%NdV(1), :, :) = real( cvalue, kind=prec )
                 case("y1")
                    self%v(:, 1, :) = real( cvalue, kind=prec )
                 case("y2")
                    self%v(:, self%NdV(2), :) = real( cvalue, kind=prec )
                 case("z1")
                    self%v(:, :, 1) = real( cvalue, kind=prec )
                 case("z2")
                    self%v(:, :, self%NdV(3)) = real( cvalue, kind=prec )
             end select
             !
        case default
             stop "Error: setOneBoundaryRScalar3D_SG > Invalid grid type"
        end select
        !
    end subroutine setOneBoundaryRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setAllInteriorRScalar3D_SG( self, cvalue )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        if( self%store_state /= compound ) then
             call self%switchStoreState
        endif
        !
        write( *, * ) "Error: setAllInteriorRScalar3D_SG to be implement: ", cvalue
        stop
        !
    end subroutine setAllInteriorRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine intBdryIndicesRScalar3D_SG( self, ind_i, ind_b )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        integer, allocatable, intent( out ) :: ind_i(:), ind_b(:)
        !
        integer :: nVecT, nBdry, nb, ni, i
        complex( kind=prec ), allocatable :: temp(:)
        type( rScalar3D_SG_t ) :: phi
        !
        if( self%is_allocated ) then
            !
            phi = rScalar3D_SG_t( self%grid, self%grid_type )
            !
        else
            stop "Error: intBdryIndicesRScalar3D_SG > Not allocated. Exiting."
        endif
        !
        if( self%store_state /= compound ) then
             call self%switchStoreState
        endif
        !
        select case( self%grid_type )
            !
            case(CORNER)
                 !
                 phi%v(1,:,:) = 1
                 phi%v(phi%nx+1,:,:) = 1
                 phi%v(:,1,:) = 1
                 phi%v(:,phi%ny+1,:) = 1
                 phi%v(:,:,1) = 1
                 phi%v(:,:,phi%nz+1) = 1
                 !
                 temp = phi%getArray()
                 !
            case default
                 stop "Error: intBdryIndicesRScalar3D_SG: Unknown self%grid_type"
                 !
        end select
        !
        nVecT = size( phi%v )
        nBdry = 0
        do i = 1, nVecT
             nBdry = nBdry + nint( real( temp(i), kind=prec ) )
        enddo
        !
        if( allocated(ind_i)) deallocate(ind_i)
        allocate(ind_i(nVecT - nBdry))
        !
        if( allocated(ind_b)) deallocate(ind_b)
        allocate(ind_b(nBdry))
        !
        nb = 0
        ni = 0
        do i = 1, nVecT
             if( nint( real( temp(i), kind=prec ) ) .EQ. 1) then
                nb = nb+1
                ind_b(nb) = i
             else
                ni = ni+1
                ind_i(ni) = i
             endif
        enddo
        !
        deallocate( temp )
        !
    end subroutine intBdryIndicesRScalar3D_SG
    !
    !> No subroutine briefing
    !
    function lengthRScalar3D_SG( self ) result( field_length )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( in ) :: self
        !
        integer :: field_length
        !
        field_length = self%Nxyz
        !
    end function lengthRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setVecComponentsRScalar3D_SG( self, xyz, &
                                             xmin, xstep, xmax, &
                                             ymin, ystep, ymax, &
                                             zmin, zstep, zmax, rvalue )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        character, intent( in ) :: xyz
        integer, intent( in ) :: xmin, xstep, xmax
        integer, intent( in ) :: ymin, ystep, ymax
        integer, intent( in ) :: zmin, zstep, zmax
        real( kind=prec ), intent( in ) :: rvalue
        !
        integer :: x1, x2
        integer :: y1, y2
        integer :: z1, z2
        !
        if( self%store_state /= compound ) then
             call self%switchStoreState
        endif
        !
        x1 = xmin; x2 = xmax
        y1 = ymin; y2 = ymax
        z1 = zmin; z2 = zmax
        !
        if( xmin == 0) x1 = self%NdV(1)
        if( xmax <= 0) x2 = self%NdV(1) + xmax
        !
        if( ymin == 0) y1 = self%NdV(2)
        if( ymax <= 0) y2 = self%NdV(2) + ymax
        !
        if( zmin == 0) z1 = self%NdV(3)
        if( zmax <= 0) z2 = self%NdV(3) + zmax
        !
        self%v(x1:x2:xstep, y1:y2:ystep, z1:z2:zstep) = rvalue
        !
    end subroutine setVecComponentsRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine zerosRScalar3D_SG( self )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        !
        if( .NOT. self%is_allocated) then
             stop "Error: zerosRScalar3D_SG > self not allocated."
        endif
        !
        if( self%store_state .EQ. compound ) then
            !
            self%v = R_ZERO
            !
        else if( self%store_state .EQ. singleton ) then
            !
            self%sv = R_ZERO
            !
        else
            stop "Error: zerosRScalar3D_SG > Unknown store_state!"
        endif
        !
    end subroutine zerosRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine sumEdgesRScalar3D_SG( self, cell_obj, interior_only )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), allocatable, intent( inout ) :: cell_obj
        logical, optional, intent( in ) :: interior_only
        !
        stop "Error: sumEdgesRScalar3D_SG not implemented yet"
        !
    end subroutine sumEdgesRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine avgCellsRScalar3D_SG( self, E_in, ptype )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: E_in
        character(*), intent( in ), optional :: ptype
        !
        stop "Error: avgCellsRScalar3D_SG not implemented yet"
        !
    end subroutine avgCellsRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine conjugateRScalar3D_SG( self )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        !
        stop "Error: conjugateRScalar3D_SG: Do not try to conjugate a real scalar!"
        !
    end subroutine conjugateRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine addRScalar3D_SG( self, rhs )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( .NOT. rhs%is_allocated) then
             stop "Error: addRScalar3D_SG > rhs not allocated."
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            if( self%store_state /= rhs%store_state ) call self%switchStoreState
            !
            select type( rhs )
                !
                class is( rScalar3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%v = self%v + rhs%v
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%sv = self%sv + rhs%sv
                        !
                    else
                        stop "Error: addRScalar3D_SG > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: addRScalar3D_SG: undefined rhs"
                !
            end select
            !
        else
            stop "Error: addRScalar3D_SG > Incompatible inputs."
        endif
        !
    end subroutine addRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine linCombRScalar3D_SG( self, rhs, c1, c2 )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        complex( kind=prec ), intent( in ) :: c1, c2
        !
        !>  linear combination, in place: self = c1*self+c2*rhs
        if( self%isCompatible(rhs)) then
            !
            if( self%store_state /= rhs%store_state ) call self%switchStoreState
            !
            select type( rhs )
                !
                class is( rScalar3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%v = c1 * self%v + c2 * rhs%v
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%sv = c1 * self%sv + c2 * rhs%sv
                        !
                    else
                        stop "Error: linCombRScalar3D_SG > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: linCombRScalar3D_SG: undefined rhs"
                !
            end select
        else
            stop "Error: linCombRScalar3D_SG > Incompatible rhs"
        endif
        !
    end subroutine linCombRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine subValueRScalar3D_SG( self, cvalue )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        if( self%store_state .EQ. compound ) then
            !
            self%v = self%v - cvalue
            !
        else if( self%store_state .EQ. singleton ) then
            !
            self%sv = self%sv - cvalue
            !
        else
            stop "Error: subValueRScalar3D_SG > Unknown self store_state!"
        endif
        !
    end subroutine subValueRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine subFieldRScalar3D_SG( self, rhs )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( self%isCompatible( rhs ) ) then
            !
            if( self%store_state /= rhs%store_state ) call self%switchStoreState
            !
            select type( rhs )
                !
                class is( rScalar3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%v = self%v - rhs%v
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%sv = self%sv - rhs%sv
                        !
                    else
                        stop "Error: subFieldRScalar3D_SG > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: subFieldRScalar3D_SG: undefined rhs"
                !
            end select
            !
        else
            stop "Error: subFieldRScalar3D_SG > Incompatible inputs."
        endif
        !
    end subroutine subFieldRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multByRealRScalar3D_SG( self, rvalue )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: rvalue
        !
        if( self%store_state .EQ. compound ) then
            !
            self%v = self%v * rvalue
            !
        else if( self%store_state .EQ. singleton ) then
            !
            self%sv = self%sv * rvalue
            !
        else
            stop "Error: multByRealRScalar3D_SG > Unknown self store_state!"
        endif
        !
    end subroutine multByRealRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multByComplexRScalar3D_SG( self, cvalue )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        if( self%store_state .EQ. compound ) then
            !
            self%v = self%v * cvalue
            !
        else if( self%store_state .EQ. singleton ) then
            !
            self%sv = self%sv * cvalue
            !
        else
            stop "Error: multByComplexRScalar3D_SG > Unknown self store_state!"
        endif
        !
    end subroutine multByComplexRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multByFieldRScalar3D_SG( self, rhs )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( self%isCompatible( rhs ) ) then
            !
            if( self%store_state /= rhs%store_state ) call self%switchStoreState
            !
            select type( rhs )
                !
                class is( rScalar3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%v = self%v * rhs%v
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%sv = self%sv * rhs%sv
                        !
                    else
                        stop "Error: multByFieldRScalar3D_SG > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: multByFieldRScalar3D_SG: undefined rhs"
                !
            end select
            !
        else
            stop "Error: multByFieldRScalar3D_SG: incompatible rhs"
        endif
        !
    end subroutine multByFieldRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multAddRScalar3D_SG( self, cvalue, rhs )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        class( Field_t ), intent( in ) :: rhs
        !
        if( self%isCompatible( rhs ) ) then
            !
            if( self%store_state /= rhs%store_state ) call self%switchStoreState
            !
            select type( rhs )
                !
                class is( rScalar3D_SG_t ) 
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%v = self%v + cvalue * rhs%v
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%sv = self%sv + cvalue * rhs%sv
                        !
                    else
                        stop "Error: multAddRScalar3D_SG > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: multAddRScalar3D_SG > rhs undefined."
                    !
            end select
            !
            !
        else
            stop "Error: multAddRScalar3D_SG >Incompatible inputs."
        endif
        !
    end subroutine multAddRScalar3D_SG
    !
    !> No subroutine briefing
    !
    function dotProdRScalar3D_SG( self, rhs ) result( cvalue )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( in ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        complex( kind=prec ) :: cvalue
        !
        if( self%isCompatible( rhs ) ) then
            !
            if( self%store_state == rhs%store_state ) then
                !
                select type( rhs )
                    !
                    class is( rScalar3D_SG_t )
                        !
                        if( rhs%store_state .EQ. compound ) then
                            !
                            cvalue = cmplx( sum( self%v * rhs%v ), 0.0, kind=prec )
                            !
                        else if( rhs%store_state .EQ. singleton ) then
                            !
                            cvalue = cmplx( sum( self%sv * rhs%sv ), 0.0, kind=prec )
                            !
                        else
                            stop "Error: dotProdRVector3D_SG > Unknown rhs store_state!"
                        endif
                        !
                    class default
                        stop "Error: dotProdRScalar3D_SG > undefined rhs"
                    !
                end select
                !
            else
                stop "Error: dotProdRScalar3D_SG > Incompatible store_state"
            endif
            !
        else
            stop "Error: dotProdRScalar3D_SG > Incompatible rhs"
        endif
        !
    end function dotProdRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine divByValueRScalar3D_SG( self, cvalue )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        if( self%store_state .EQ. compound ) then
            !
            self%v = self%v / cvalue
            !
        else if( self%store_state .EQ. singleton ) then
            !
            self%sv = self%sv / cvalue
            !
        else
            stop "Error: divByValueRScalar3D_SG > Unknown self store_state!"
        endif
        !
    end subroutine divByValueRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine divByFieldRScalar3D_SG( self, rhs )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( self%isCompatible( rhs ) ) then
            !
            if( self%store_state /= rhs%store_state ) call self%switchStoreState
            !
            select type( rhs )
                !
                class is( rScalar3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%v = self%v / rhs%v
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%sv = self%sv / rhs%sv
                        !
                    else
                        stop "Error: divByFieldRScalar3D_SG > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: divByFieldRScalar3D_SG: undefined rhs"
                !
            end select
            !
        else
            stop "Error: divByFieldRScalar3D_SG: incompatible rhs"
        endif
        !
    end subroutine divByFieldRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine getRealRScalar3D_SG( self, r_field )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( in ) :: self
        class( Field_t ), allocatable, intent( out ) :: r_field
        !
        allocate( r_field, source = rScalar3D_SG_t( self%grid, self%grid_type ) )
        !
        call r_field%copyFrom( self )
        !
    end subroutine getRealRScalar3D_SG
    !
    !> No subroutine briefing
    !
    function getArrayRScalar3D_SG( self ) result( array )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( in ) :: self
        complex( kind=prec ), allocatable, dimension(:) :: array
        !
        if( self%store_state .EQ. compound ) then
            !
            allocate( array( self%length() ) )
            array = (/reshape( cmplx( self%v, 0.0, kind=prec ), (/self%Nxyz, 1/))/)
            !
        else if( self%store_state .EQ. singleton ) then
            !
            array = self%sv
            !
        else
            stop "Error: getArrayRScalar3D_SG > Unknown store_state!"
        endif
        !
    end function getArrayRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setArrayRScalar3D_SG( self, array )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), dimension(:), intent( in ) :: array
        !
        if( self%store_state .EQ. compound ) then
            !
            self%v = reshape( real( array, kind=prec ), (/self%NdV(1), self%NdV(2), self%NdV(3)/) )
            !
        else if( self%store_state .EQ. singleton ) then
            !
            self%sv = array
            !
        else
            stop "Error: setArrayRScalar3D_SG > Unknown store_state!"
        endif
        !
    end subroutine setArrayRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine switchStoreStateRScalar3D_SG( self )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        !
        integer :: nzAir
        !
        select case( self%store_state )
            !
            case( compound )
                !
                allocate( self%sv( self%length() ) )
                !
                self%sv = (/reshape( self%v, (/self%Nxyz, 1/))/)
                !
                deallocate( self%v )
                !
                self%store_state = singleton
                !
            case( singleton )
                !
                if( self%grid_type == CORNER ) then
                    !
                    allocate( self%v( self%nx + 1, self%ny + 1, self%nz + 1 ) )
                    !
                else if( self%grid_type == CENTER ) then
                    !
                    allocate( self%v( self%nx, self%ny, self%nz ) )
                    !
                else if( self%grid_type == CELL_EARTH ) then
                    !
                    call self%grid%getDimensions( self%nx, self%ny, self%nz, nzAir )
                    !
                    allocate( self%v( self%nx, self%ny, self%nz - nzAir ) )
                    !
                else
                     write( *, * ) "Error: switchStoreStateRScalar3D_SG > unrecognized grid type: [", self%grid_type, "]"
                     stop
                endif
                !
                self%v = reshape( self%sv, (/self%NdV(1), self%NdV(2), self%NdV(3)/) )
                !
                deallocate( self%sv )
                !
                self%store_state = compound
                !
            case default
                write( *, * ) "Error: switchStoreStateRScalar3D_SG > Unknown store_state :[", self%store_state, "]"
                stop
            !
        end select
        !
    end subroutine switchStoreStateRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine copyFromRScalar3D_SG( self, rhs )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( .NOT. rhs%is_allocated ) then
            stop "Error: copyFromRScalar3D_SG > rhs not allocated"
        endif
        !
        self%grid => rhs%grid
        self%grid_type = rhs%grid_type
        self%nx = rhs%nx
        self%ny = rhs%ny
        self%nz = rhs%nz
        self%store_state = rhs%store_state
        !
        select type( rhs )
            !
            class is( rScalar3D_SG_t )
                !
                self%NdV = rhs%NdV
                self%Nxyz = rhs%Nxyz
                !
                if( rhs%store_state .EQ. compound ) then
                    !
                    self%v = rhs%v
                    !
                else if( rhs%store_state .EQ. singleton ) then
                    !
                    self%sv = rhs%sv
                    !
                else
                    stop "Error: copyFromRScalar3D_SG > Unknown store_state!"
                endif
                !
            class default
                stop "Error: copyFromRScalar3D_SG > Unclassified rhs"
            !
        end select
        !
        self%is_allocated = .TRUE.
        !
    end subroutine copyFromRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine readRScalar3D_SG( self, funit, ftype )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        integer, intent( in ) :: funit
        character(:), allocatable, intent( in ), optional :: ftype
        !
        integer :: Nx, Ny, Nz
        character(4) :: grid_type
        integer :: i, j, k, k1, k2, istat
        real( kind=prec ), allocatable :: temp(:)
        logical :: ok, hasname, binary
        character(:), allocatable :: fname, isbinary
        !
        if( self%store_state /= compound ) then
             call self%switchStoreState
        endif
        !
        if( .NOT. present( ftype ) ) then
             binary = .FALSE.
        else if( index( ftype, "b" ) > 0) then
             binary = .TRUE.
        else
             binary = .FALSE.
        endif
        !
        inquire( funit, opened = ok, named = hasname, name = fname, unformatted = isbinary )
        !
        if( ok ) then
            !
            !> check that the file is unformatted if binary, formatted if ascii
            if( (index(isbinary, "yes") > 0 .OR. index(isbinary, "YES") > 0) &
                     .AND.  .NOT. binary) then             
                 write( *, * ) "Error: rScalar3D_SG_t::readRScalar3D_SG: "
                 write( *, * ) "            Unable to read scalar from unformatted file ", &
                            trim(fname), ".Exiting."
                 stop
            else if( (index(isbinary, "no") > 0 .OR. index(isbinary, "NO") > 0) &
                     .AND.binary) then
                 write( *, * ) "Error: rScalar3D_SG_t::readRScalar3D_SG: "
                 write( *, * ) "            Unable to read scalar from formatted file ", &
                            trim(fname), "."
                 stop
            endif
            !
            if( binary) then
                 !> read binary from unformatted files
                 read(funit) self%Nx, self%Ny, self%Nz, grid_type
                 read(funit) self%v
            endif
            !
            Nx = size(self%v, 1)
            Ny = size(self%v, 2)
            Nz = size(self%v, 3)
            !
            allocate(temp(Ny), STAT = istat)
            !
            i = 1
            do
                 read(funit, *, iostat = istat) k1, k2
                 if( istat /= 0) exit
                 !
                 if( (k1 < 0) .OR. (k2 > Nz)) then
                        write( *, * ) "Error: rScalar3D_SG::readRScalar3D_SG: "
                        write( *, * ) "      While reading the ", i, "th block."
                        stop
                 else if( k1 > k2) then
                        write( *, * ) "Warning: rScalar3D_SG::readRScalar3D_SG: "
                        write( *, * ) "                Block ", i, " will be ignored."
                 endif
                 !
                 do j = Nx, 1, -1
                        read(funit, *, iostat = istat) temp
                        
                        if( istat /= 0) then
                             write( *, * ) "Error: rScalar3D_SG::readRScalar3D_SG: "
                             write( *, * ) "            While reading the ", j, "th row in ", i,"th block."
                             stop
                        endif
                        
                        do k = k1, k2
                             self%v(j, :, k) = temp
                        enddo
                 enddo
                 !
                 if( k == Nz) exit
                 !
                 i = i + 1
                 !
            enddo
            !
            deallocate( temp )
            !
        else
            stop "Error: readRScalar3D_SG: unable to open file"
        endif
        !
    end subroutine readRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine writeRScalar3D_SG( self, funit, ftype )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        integer, intent( in ) :: funit
        character(:), allocatable, intent( in ), optional :: ftype
        !
        integer :: Nx, Ny, Nz
        integer :: i, j, k, k1, k2, istat
        real( kind=prec ), allocatable, dimension(:, :) :: temp
        logical :: ok, hasname, binary
        character(80) :: fname, isbinary
        !
        if(  .NOT. self%is_allocated) then
             stop "Error: writeRScalar3D_SG > Not allocated"
        endif
        !
        !> Make sure the store_state is compound
        if( self%store_state /= compound ) then
             call self%switchStoreState
        endif
        !
        if(  .NOT. present( ftype ) ) then
             binary = .FALSE.
        else if( index( ftype, "b" ) > 0) then
             binary = .TRUE.
        else
             binary = .FALSE.
        endif
        !
        inquire( funit, opened = ok, named = hasname, name = fname, unformatted = isbinary )
        !
        if( ok ) then
            !
            if( (index(isbinary, "yes") > 0 .OR. index(isbinary, "YES") > 0) &
                     .AND. .NOT. binary) then             
                 write( *, * ) "Error: writeRScalar3D_SG > Unable to write vector to unformatted file ", &
                            trim(fname), "."
                 !
                 stop
            else if( (index(isbinary,"no") > 0 .OR. index(isbinary,"NO") > 0) &
                     .AND.binary) then
                 write( *, * ) "Error: writeRScalar3D_SG > Unable to write vector to formatted file ", &
                            trim(fname), "."
                 !
                 stop
            endif
            !
            if( binary) then
                 write(funit) self%nx, self%ny, self%nz, self%grid_type
                 write(funit) self%v             
                 return
            endif
            !
            !
            !> ASCII format
            !
            write(funit, "(3i5,a10)", iostat = istat) self%nx, self%ny, self%nz, trim(self%grid_type)
            !
            Nx = size(self%v, 1)
            Ny = size(self%v, 2)
            Nz = size(self%v, 3)
            !
            allocate(temp(Nx, Ny), STAT = istat)
            !
            k1 = 1
            do
                k2 = Nz
                do k = k1, Nz - 1
                    temp = abs( self%v(:, :, k + 1) - self%v(:, :, k) )
                    if( maxval( real( temp ) ) > TOL6 ) then
                        k2 = k
                        exit
                    endif
                enddo
                !
                write( funit, "(2i5)", iostat = istat ) k1, k2
                !
                if( istat /= 0) then
                    stop "Error: writeRScalar3D_SG > Failed while writing to file."
                endif
                !
                temp = self%v(:, :, k1)
                !
                do i = Nx, 1, -1
                    do j = 1, Ny
                         write(funit, "(es13.5)", iostat = istat, &
                                    advance = "no") self%v(i, j, k1)
                    enddo
                    write(funit, *)
                enddo
                !
                k1 = k2 + 1
                !
                if( k1 > Nz) exit
            enddo
            !
            deallocate( temp )
            !
        else
            stop "Error: readRVector3D_SG: unable to open file"
        endif
        !
    end subroutine writeRScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine printRScalar3D_SG( self, io_unit, title, append )
        implicit none
        !
        class( rScalar3D_SG_t ), intent( inout ) :: self
        integer, intent( in ), optional :: io_unit
        character(*), intent( in ), optional :: title
        logical, intent( in ), optional :: append
        !
        integer :: ix, iy, iz,funit
        !
        if( self%store_state /= compound ) then
             call self%switchStoreState
        endif
        !
        if( present( io_unit ) ) then
            funit = io_unit
        else
            funit = 0    !>    usually this will work to write to standard output
        endif
        if(present(title)) then
          write(funit,*) title
        endif
        !
        write( funit, * ) self%nx, self%ny, self%nz
        !
        write(funit,*) "scalar field"
        do ix = 1, self%nx
             do iy = 1, self%ny
                  do iz = 1, self%nz
                        if( self%v( ix, iy, iz ) /= 0 ) then
                            write(funit,*) ix,iy,iz, ":[", self%v( ix, iy, iz ), "]"
                        endif
                  enddo
             enddo
        enddo
        !
    end subroutine printRScalar3D_SG
    !
end module rScalar3D_SG
