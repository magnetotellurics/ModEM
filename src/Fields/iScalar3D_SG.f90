!
!> Derived class to define a iScalar3D_SG 
!
module iScalar3D_SG
    !
    use Scalar
    use Grid3D_SG
    !
    type, extends( Scalar_t ) :: iScalar3D_SG_t
        !
        integer( kind=prec ), allocatable, dimension(:, :, :) :: v
        !
        integer( kind=prec ), allocatable, dimension(:) :: sv
        !
        contains
            !
            !> Destructor
            final :: iScalar3D_SG_dtor
            !
            !> Boundary operations
            procedure, public :: setAllBoundary => setAllBoundaryiScalar3D_SG
            procedure, public :: setOneBoundary => setOneBoundaryiScalar3D_SG
            procedure, public :: setAllInterior => setAllInterioriScalar3D_SG
            procedure, public :: intBdryIndices => intBdryIndicesiScalar3D_SG
            !
            !> Dimensioning operations
            procedure, public :: length => lengthiScalar3D_SG
            procedure, public :: setVecComponents => setVecComponentsiScalar3D_SG
            !
            !> Arithmetic/algebraic unary operations
            procedure, public :: zeros => zerosiScalar3D_SG
            procedure, public :: sumEdges => sumEdgesiScalar3D_SG
            procedure, public :: avgCells => avgCellsiScalar3D_SG
            procedure, public :: conjugate => conjugateiScalar3D_SG
            !
            !> Arithmetic/algebraic binary operations
            procedure, public :: add => addiScalar3D_SG
            !
            procedure, public :: linComb => linCombiScalar3D_SG
            !
            procedure, public :: subValue => subValueiScalar3D_SG
            procedure, public :: subField => subFieldiScalar3D_SG
            !
            procedure, public :: multByReal => multByRealiScalar3D_SG
            procedure, public :: multByComplex => multByComplexiScalar3D_SG
            procedure, public :: multByField => multByFieldiScalar3D_SG
            !
            procedure, public :: multAdd => multAddiScalar3D_SG
            !
            procedure, public :: dotProd => dotProdiScalar3D_SG
            !
            procedure, public :: divByField => divByFieldiScalar3D_SG
            procedure, public :: divByValue => divByValueiScalar3D_SG
            !
            !> Miscellaneous
            procedure, public :: getReal => getRealiScalar3D_SG
            procedure, public :: getArray => getArrayiScalar3D_SG
            procedure, public :: setArray => setArrayiScalar3D_SG
            procedure, public :: switchStoreState => switchStoreStateiScalar3D_SG
            procedure, public :: copyFrom => copyFromiScalar3D_SG
            !
            !> I/O operations
            procedure, public :: read => readiScalar3D_SG
            procedure, public :: write => writeiScalar3D_SG
            procedure, public :: print => printiScalar3D_SG
            !
    end type iScalar3D_SG_t
    !
    interface iScalar3D_SG_t
        module procedure iScalar3D_SG_ctor
    end interface iScalar3D_SG_t
    !
contains
    !
    !> No subroutine briefing
    !
    function iScalar3D_SG_ctor( grid, grid_type ) result ( self )
        implicit none
        !
        class( Grid_t ), target, intent( in ) :: grid
        character( len=4 ), intent( in ) :: grid_type
        !
        type( iScalar3D_SG_t ) :: self
        !
        integer :: nx, ny, nz, nzAir, nz_earth
        integer :: status
        !
        !write( *, * ) "Constructor iScalar3D_SG"
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
        if( grid_type == CORNER) then
             allocate(self%v(nx + 1, ny + 1, nz + 1), STAT = status)    
             self%NdV = (/self%nx + 1, self%ny + 1, self%nz + 1/)
             
        elseif( grid_type == CENTER) then             
             allocate(self%v(nx, ny, nz), STAT = status) 
             self%NdV = (/self%nx, self%ny, self%nz/)
             
        elseif( grid_type == CELL_EARTH) then
             self%nz = nz_earth
             allocate(self%v(nx, ny, nz_earth), STAT = status)
             self%NdV = (/nx, ny, nz_earth/)
             
        else
            write( *, * ) "     "//achar(27)//"[31m# Error:"//achar(27)//"[0m iScalar3D_SG_ctor > unrecognized grid type: [", grid_type, "]"
            stop
        endif
        !
        self%is_allocated = self%is_allocated.AND.(status .EQ. 0)
        !
        if( self%is_allocated ) then
             self%v = R_ZERO
        else
             stop "Error: iScalar3D_SG_ctor > Unable to allocate rScalar - invalid grid supplied"
        endif
        !
        self%Nxyz = product( self%NdV )
        !
    end function iScalar3D_SG_ctor
    !
    !> No subroutine briefing
    !
    subroutine iScalar3D_SG_dtor( self )
        implicit none
        !
        type( iScalar3D_SG_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor iScalar3D_SG"
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
    end subroutine iScalar3D_SG_dtor
    !
    !> No subroutine briefing
    !
    subroutine setAllBoundaryiScalar3D_SG( self, cvalue )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
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
                 stop "Error: setAllBoundaryiScalar3D_SG > Grid type not recognized."
        end select
        !
    end subroutine setAllBoundaryiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setOneBoundaryiScalar3D_SG( self, bdry, cvalue, int_only )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        character(:), allocatable, intent( in ) :: bdry
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
             stop "Error: setOneBoundaryiScalar3D_SG > Invalid grid type"
        end select
        !
    end subroutine setOneBoundaryiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setAllInterioriScalar3D_SG( self, cvalue )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        if( self%store_state /= compound ) then
             call self%switchStoreState
        endif
        !
        write( *, * ) "     "//achar(27)//"[31m# Error:"//achar(27)//"[0m setAllInterioriScalar3D_SG to be implement: ", cvalue
        stop
        !
    end subroutine setAllInterioriScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine intBdryIndicesiScalar3D_SG( self, ind_i, ind_b )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        integer, allocatable, intent( out ) :: ind_i(:), ind_b(:)
        !
        integer :: nVecT, nBdry, nb, ni, i
        complex( kind=prec ), allocatable :: temp(:)
        type( iScalar3D_SG_t ) :: phi
        !
        if( self%is_allocated ) then
            !
            phi = iScalar3D_SG_t( self%grid, self%grid_type )
            !
        else
            stop "Error: intBdryIndicesiScalar3D_SG > Not allocated. Exiting."
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
                 stop "Error: intBdryIndicesiScalar3D_SG: Unknown self%grid_type"
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
    end subroutine intBdryIndicesiScalar3D_SG
    !
    !> No subroutine briefing
    !
    function lengthiScalar3D_SG( self ) result( field_length )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( in ) :: self
        !
        integer :: field_length
        !
        field_length = self%Nxyz
        !
    end function lengthiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setVecComponentsiScalar3D_SG( self, xyz, &
                                             xmin, xstep, xmax, &
                                             ymin, ystep, ymax, &
                                             zmin, zstep, zmax, rvalue )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
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
    end subroutine setVecComponentsiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine zerosiScalar3D_SG( self )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        !
        if( .NOT. self%is_allocated) then
             stop "Error: zerosiScalar3D_SG > self not allocated."
        endif
        !
        if( self%store_state .EQ. compound ) then
            !
            self%v = R_ZERO
            !
        elseif( self%store_state .EQ. singleton ) then
            !
            self%sv = R_ZERO
            !
        else
            stop "Error: zerosiScalar3D_SG > Unknown store_state!"
        endif
        !
    end subroutine zerosiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine sumEdgesiScalar3D_SG( self, cell_obj, interior_only )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), allocatable, intent( inout ) :: cell_obj
        logical, optional, intent( in ) :: interior_only
        !
        stop "Error: sumEdgesiScalar3D_SG not implemented yet"
        !
    end subroutine sumEdgesiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine avgCellsiScalar3D_SG( self, E_in, ptype )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: E_in
        character(*), intent( in ), optional :: ptype
        !
        stop "Error: avgCellsiScalar3D_SG not implemented yet"
        !
    end subroutine avgCellsiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine conjugateiScalar3D_SG( self )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        !
        stop "Error: conjugateiScalar3D_SG: Do not try to conjugate a real scalar!"
        !
    end subroutine conjugateiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine addiScalar3D_SG( self, rhs )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( .NOT. rhs%is_allocated) then
             stop "Error: addiScalar3D_SG > rhs not allocated."
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            if( self%store_state /= rhs%store_state ) call self%switchStoreState
            !
            select type( rhs )
                !
                class is( iScalar3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%v = self%v + rhs%v
                        !
                    elseif( rhs%store_state .EQ. singleton ) then
                        !
                        self%sv = self%sv + rhs%sv
                        !
                    else
                        stop "Error: addiScalar3D_SG > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: addiScalar3D_SG: undefined rhs"
                !
            end select
            !
        else
            stop "Error: addiScalar3D_SG > Incompatible inputs."
        endif
        !
    end subroutine addiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine linCombiScalar3D_SG( self, rhs, c1, c2 )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
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
                class is( iScalar3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%v = c1 * self%v + c2 * rhs%v
                        !
                    elseif( rhs%store_state .EQ. singleton ) then
                        !
                        self%sv = c1 * self%sv + c2 * rhs%sv
                        !
                    else
                        stop "Error: linCombiScalar3D_SG > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: linCombiScalar3D_SG: undefined rhs"
                !
            end select
        else
            stop "Error: linCombiScalar3D_SG > Incompatible rhs"
        endif
        !
    end subroutine linCombiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine subValueiScalar3D_SG( self, cvalue )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        if( self%store_state .EQ. compound ) then
            !
            self%v = self%v - cvalue
            !
        elseif( self%store_state .EQ. singleton ) then
            !
            self%sv = self%sv - cvalue
            !
        else
            stop "Error: subValueiScalar3D_SG > Unknown self store_state!"
        endif
        !
    end subroutine subValueiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine subFieldiScalar3D_SG( self, rhs )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( self%isCompatible( rhs ) ) then
            !
            if( self%store_state /= rhs%store_state ) call self%switchStoreState
            !
            select type( rhs )
                !
                class is( iScalar3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%v = self%v - rhs%v
                        !
                    elseif( rhs%store_state .EQ. singleton ) then
                        !
                        self%sv = self%sv - rhs%sv
                        !
                    else
                        stop "Error: subFieldiScalar3D_SG > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: subFieldiScalar3D_SG: undefined rhs"
                !
            end select
            !
        else
            stop "Error: subFieldiScalar3D_SG > Incompatible inputs."
        endif
        !
    end subroutine subFieldiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multByRealiScalar3D_SG( self, rvalue )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: rvalue
        !
        if( self%store_state .EQ. compound ) then
            !
            self%v = self%v * rvalue
            !
        elseif( self%store_state .EQ. singleton ) then
            !
            self%sv = self%sv * rvalue
            !
        else
            stop "Error: multByRealiScalar3D_SG > Unknown self store_state!"
        endif
        !
    end subroutine multByRealiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multByComplexiScalar3D_SG( self, cvalue )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        if( self%store_state .EQ. compound ) then
            !
            self%v = self%v * cvalue
            !
        elseif( self%store_state .EQ. singleton ) then
            !
            self%sv = self%sv * cvalue
            !
        else
            stop "Error: multByComplexiScalar3D_SG > Unknown self store_state!"
        endif
        !
    end subroutine multByComplexiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multByFieldiScalar3D_SG( self, rhs )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( self%isCompatible( rhs ) ) then
            !
            if( self%store_state /= rhs%store_state ) call self%switchStoreState
            !
            select type( rhs )
                !
                class is( iScalar3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%v = self%v * rhs%v
                        !
                    elseif( rhs%store_state .EQ. singleton ) then
                        !
                        self%sv = self%sv * rhs%sv
                        !
                    else
                        stop "Error: multByFieldiScalar3D_SG > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: multByFieldiScalar3D_SG: undefined rhs"
                !
            end select
            !
        else
            stop "Error: multByFieldiScalar3D_SG: incompatible rhs"
        endif
        !
    end subroutine multByFieldiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multAddiScalar3D_SG( self, cvalue, rhs )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        class( Field_t ), intent( in ) :: rhs
        !
        if( self%isCompatible( rhs ) ) then
            !
            if( self%store_state /= rhs%store_state ) call self%switchStoreState
            !
            select type( rhs )
                !
                class is( iScalar3D_SG_t ) 
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%v = self%v + cvalue * rhs%v
                        !
                    elseif( rhs%store_state .EQ. singleton ) then
                        !
                        self%sv = self%sv + cvalue * rhs%sv
                        !
                    else
                        stop "Error: multAddiScalar3D_SG > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: multAddiScalar3D_SG > rhs undefined."
                    !
            end select
            !
            !
        else
            stop "Error: multAddiScalar3D_SG >Incompatible inputs."
        endif
        !
    end subroutine multAddiScalar3D_SG
    !
    !> No subroutine briefing
    !
    function dotProdiScalar3D_SG( self, rhs ) result( cvalue )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( in ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        complex( kind=prec ) :: cvalue
        !
        type( iScalar3D_SG_t ) :: aux_vec
        !
        if( self%isCompatible( rhs ) ) then
            !
            aux_vec = self
            if( aux_vec%store_state /= rhs%store_state ) call aux_vec%switchStoreState
            !
            select type( rhs )
                !
                class is( iScalar3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        cvalue = cmplx( sum( aux_vec%v * rhs%v ), 0.0, kind=prec )
                        !
                    elseif( rhs%store_state .EQ. singleton ) then
                        !
                        cvalue = cmplx( sum( aux_vec%sv * rhs%sv ), 0.0, kind=prec )
                        !
                    else
                        stop "Error: dotProdRVector3D_SG > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: dotProdiScalar3D_SG > undefined rhs"
                !
            end select
            !
        else
            stop "Error: dotProdiScalar3D_SG > Incompatible rhs"
        endif
        !
    end function dotProdiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine divByValueiScalar3D_SG( self, cvalue )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        if( self%store_state .EQ. compound ) then
            !
            self%v = self%v / cvalue
            !
        elseif( self%store_state .EQ. singleton ) then
            !
            self%sv = self%sv / cvalue
            !
        else
            stop "Error: divByValueiScalar3D_SG > Unknown self store_state!"
        endif
        !
    end subroutine divByValueiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine divByFieldiScalar3D_SG( self, rhs )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( self%isCompatible( rhs ) ) then
            !
            if( self%store_state /= rhs%store_state ) call self%switchStoreState
            !
            select type( rhs )
                !
                class is( iScalar3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%v = self%v / rhs%v
                        !
                    elseif( rhs%store_state .EQ. singleton ) then
                        !
                        self%sv = self%sv / rhs%sv
                        !
                    else
                        stop "Error: divByFieldiScalar3D_SG > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: divByFieldiScalar3D_SG: undefined rhs"
                !
            end select
            !
        else
            stop "Error: divByFieldiScalar3D_SG: incompatible rhs"
        endif
        !
    end subroutine divByFieldiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine getRealiScalar3D_SG( self, r_field )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( in ) :: self
        class( Field_t ), allocatable, intent( out ) :: r_field
        !
        allocate( r_field, source = iScalar3D_SG_t( self%grid, self%grid_type ) )
        !
        call r_field%copyFrom( self )
        !
    end subroutine getRealiScalar3D_SG
    !
    !> No subroutine briefing
    !
    function getArrayiScalar3D_SG( self ) result( array )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( in ) :: self
        complex( kind=prec ), allocatable, dimension(:) :: array
        !
        if( self%store_state .EQ. compound ) then
            !
            allocate( array( self%length() ) )
            array = (/reshape( cmplx( self%v, 0.0, kind=prec ), (/self%Nxyz, 1/))/)
            !
        elseif( self%store_state .EQ. singleton ) then
            !
            array = self%sv
            !
        else
            stop "Error: getArrayiScalar3D_SG > Unknown store_state!"
        endif
        !
    end function getArrayiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setArrayiScalar3D_SG( self, array )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), dimension(:), intent( in ) :: array
        !
        if( self%store_state .EQ. compound ) then
            !
            self%v = reshape( real( array, kind=prec ), (/self%NdV(1), self%NdV(2), self%NdV(3)/) )
            !
        elseif( self%store_state .EQ. singleton ) then
            !
            self%sv = array
            !
        else
            stop "Error: setArrayiScalar3D_SG > Unknown store_state!"
        endif
        !
    end subroutine setArrayiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine switchStoreStateiScalar3D_SG( self )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
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
                elseif( self%grid_type == CENTER ) then
                    !
                    allocate( self%v( self%nx, self%ny, self%nz ) )
                    !
                elseif( self%grid_type == CELL_EARTH ) then
                    !
                    call self%grid%getDimensions( self%nx, self%ny, self%nz, nzAir )
                    !
                    allocate( self%v( self%nx, self%ny, self%nz - nzAir ) )
                    !
                else
                     write( *, * ) "     "//achar(27)//"[31m# Error:"//achar(27)//"[0m switchStoreStateCScalar3D_SG > unrecognized grid type: [", self%grid_type, "]"
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
                write( *, * ) "     "//achar(27)//"[31m# Error:"//achar(27)//"[0m switchStoreStateiScalar3D_SG > Unknown store_state :[", self%store_state, "]"
                stop
            !
        end select
        !
    end subroutine switchStoreStateiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine copyFromiScalar3D_SG( self, rhs )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( .NOT. rhs%is_allocated ) then
            stop "Error: copyFromiScalar3D_SG > rhs not allocated"
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
            class is( iScalar3D_SG_t )
                !
                self%NdV = rhs%NdV
                self%Nxyz = rhs%Nxyz
                !
                if( rhs%store_state .EQ. compound ) then
                    !
                    self%v = rhs%v
                    !
                elseif( rhs%store_state .EQ. singleton ) then
                    !
                    self%sv = rhs%sv
                    !
                else
                    stop "Error: copyFromiScalar3D_SG > Unknown store_state!"
                endif
                !
            class default
                stop "Error: copyFromiScalar3D_SG > Unclassified rhs"
            !
        end select
        !
        self%is_allocated = .TRUE.
        !
    end subroutine copyFromiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine readiScalar3D_SG( self, funit, ftype )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
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
        elseif( index( ftype, "b" ) > 0) then
             binary = .TRUE.
        else
             binary = .FALSE.
        endif
        !
        !inquire( funit, opened = ok, named = hasname, name = fname, unformatted = isbinary )
        !
        !if( ok ) then
            !
            !> check that the file is unformatted if binary, formatted if ascii
            !
            if( binary ) then
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
                        write( *, * ) "     "//achar(27)//"[31m# Error:"//achar(27)//"[0m iScalar3D_SG::readiScalar3D_SG: "
                        write( *, * ) "      While reading the ", i, "th block."
                        stop
                 elseif( k1 > k2) then
                        write( *, * ) "Warning: iScalar3D_SG::readiScalar3D_SG: "
                        write( *, * ) "                Block ", i, " will be ignored."
                 endif
                 !
                 do j = Nx, 1, -1
                        read(funit, *, iostat = istat) temp
                        
                        if( istat /= 0) then
                             write( *, * ) "     "//achar(27)//"[31m# Error:"//achar(27)//"[0m iScalar3D_SG::readiScalar3D_SG: "
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
        !else
            !stop "Error: readiScalar3D_SG: unable to open file"
        !endif
        !
    end subroutine readiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine writeiScalar3D_SG( self, funit, ftype )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
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
             stop "Error: writeiScalar3D_SG > Not allocated"
        endif
        !
        !> Make sure the store_state is compound
        if( self%store_state /= compound ) then
             call self%switchStoreState
        endif
        !
        if(  .NOT. present( ftype ) ) then
             binary = .FALSE.
        elseif( index( ftype, "b" ) > 0) then
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
                 write( *, * ) "     "//achar(27)//"[31m# Error:"//achar(27)//"[0m writeiScalar3D_SG > Unable to write vector to unformatted file ", &
                            trim(fname), "."
                 !
                 stop
            elseif( (index(isbinary,"no") > 0 .OR. index(isbinary,"NO") > 0) &
                     .AND.binary) then
                 write( *, * ) "     "//achar(27)//"[31m# Error:"//achar(27)//"[0m writeiScalar3D_SG > Unable to write vector to formatted file ", &
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
                    stop "Error: writeiScalar3D_SG > Failed while writing to file."
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
    end subroutine writeiScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine printiScalar3D_SG( self, io_unit, title, append )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
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
    end subroutine printiScalar3D_SG
    !
end module iScalar3D_SG
