!
!> Derived class to define a iScalar3D_SG 
!
module iScalar3D_SG
    !
    use Scalar
    !
    type, extends( Scalar_t ) :: iScalar3D_SG_t
        !
        integer( kind=prec ), allocatable, dimension(:, :, :) :: v
        !
        integer( kind=prec ), allocatable, dimension(:) :: s_v
        !
        contains
            !
            !> Destructor
            final :: iScalar3D_SG_dtor
            !
            !> Boundary operations
            procedure, public :: setAllBoundary => setAllBoundary_iScalar3D_SG
            procedure, public :: setOneBoundary => setOneBoundary_iScalar3D_SG
            procedure, public :: intBdryIndices => intBdryIndices_iScalar3D_SG
            !
            !> Dimensioning operations
            procedure, public :: setVecComponents => setVecComponents_iScalar3D_SG
            !
            !> Arithmetic/algebraic unary operations
            procedure, public :: zeros => zeros_iScalar3D_SG
            procedure, public :: conjugate => conjugate_iScalar3D_SG
            !
            !> Arithmetic/algebraic binary operations
            procedure, public :: add => add_iScalar3D_SG
            !
            procedure, public :: linComb => linComb_iScalar3D_SG
            !
            procedure, public :: subValue => subValue_iScalar3D_SG
            procedure, public :: subField => subField_iScalar3D_SG
            !
            procedure, public :: multByReal => multByReal_iScalar3D_SG
            procedure, public :: multByComplex => multByComplex_iScalar3D_SG
            procedure, public :: multByField => multByField_iScalar3D_SG
            !
            procedure, public :: multAdd => multAdd_iScalar3D_SG
            !
            procedure, public :: dotProd => dotProd_iScalar3D_SG
            !
            procedure, public :: divByField => divByField_iScalar3D_SG
            procedure, public :: divByValue => divByValue_iScalar3D_SG
            !
            !> Getters & Setters
            procedure, public :: getV => getV_iScalar3D_SG
            procedure, public :: setV => setV_iScalar3D_SG
            !
            procedure, public :: getSV => getSV_iScalar3D_SG
            procedure, public :: setSV => setSV_iScalar3D_SG
            !
            procedure, public :: getArray => getArray_iScalar3D_SG
            procedure, public :: setArray => setArray_iScalar3D_SG
            !
            !> Miscellaneous
            procedure, public :: copyFrom => copyFrom_iScalar3D_SG
            !
            !> I/O operations
            procedure, public :: read => read_iScalar3D_SG
            procedure, public :: write => write_iScalar3D_SG
            procedure, public :: print => print_iScalar3D_SG
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
        call self%baseInit
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
        !
        self%is_allocated = .TRUE.
        !
        if( grid_type == NODE ) then
             !
             allocate( self%v(nx + 1, ny + 1, nz + 1), STAT = status )
             self%NdV = (/self%nx + 1, self%ny + 1, self%nz + 1/)
             !
        elseif( grid_type == CELL ) then
             !
             allocate(self%v(nx, ny, nz), STAT = status) 
             self%NdV = (/self%nx, self%ny, self%nz/)
             !
        elseif( grid_type == CELL_EARTH ) then
             !
             self%nz = nz_earth
             allocate(self%v(nx, ny, nz_earth), STAT = status)
             self%NdV = (/nx, ny, nz_earth/)
             !
        else
            call errStop( "iScalar3D_SG_ctor > unrecognized grid type: ["//grid_type//"]" )
        endif
        !
        self%is_allocated = self%is_allocated .AND. ( status .EQ. 0 )
        !
        if( self%is_allocated ) then
             self%v = R_ZERO
        else
            call errStop( "iScalar3D_SG_ctor > Unable to allocate rScalar - invalid grid supplied" )
        endif
        !
        self%Nxyz = product( self%NdV )
        !
        call self%setIndexArrays
        call self%zeros
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
        if( .NOT. self%is_allocated ) then
            call errStop( "iScalar3D_SG_dtor > self not allocated." )
        endif
        !
        call self%baseDealloc
        !
        if( allocated( self%v ) ) deallocate( self%v )
        if( allocated( self%s_v ) ) deallocate( self%s_v )
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
    subroutine setAllBoundary_iScalar3D_SG( self, cvalue )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "setAllBoundary_iScalar3D_SG > self not allocated." )
        endif
        !
        call self%switchStoreState( compound )
        !
        select case( self%grid_type )
            !
            case( NODE, CELL, CELL_EARTH ) 
                !
                self%v((/1, self%NdV(1)/), :, :) = cvalue
                self%v(:, (/1, self%NdV(2)/), :) = cvalue
                self%v(:, :, (/1, self%NdV(3)/)) = cvalue
                !
            case default
                call errStop( "setAllBoundary_iScalar3D_SG > grid_type ["//self%grid_type//"] not recognized." )
            !
        end select
        !
    end subroutine setAllBoundary_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setOneBoundary_iScalar3D_SG( self, bdry, cvalue, int_only )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        character(*), intent( in ) :: bdry
        complex( kind=prec ), intent( in ) :: cvalue
        logical, intent( in ), optional :: int_only
        !
        logical :: int_only_p
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "setOneBoundary_iScalar3D_SG > self not allocated." )
        endif
        !
        call self%switchStoreState( compound )
        !
        if( .NOT. present( int_only ) ) then
             int_only_p = .FALSE.
        else 
             int_only_p = int_only
        endif
        !
        select case( self%grid_type )
            !
            case( NODE )
                if( int_only_p ) then
                    select case(bdry)
                        case("x1")
                            self%v(1, 2:self%NdV(2)-1, 2:self%NdV(3)-1) = cvalue 
                        case("x2")
                            self%v(self%NdV(1), 2:self%NdV(2)-1, 2:self%NdV(3)-1) = cvalue
                        case("y1")
                            self%v(2:self%NdV(1)-1, 1, 2:self%NdV(3)-1) = cvalue
                        case("y2")
                            self%v(2:self%NdV(1)-1, self%NdV(2), 2:self%NdV(3)-1) = cvalue
                        case("z1")
                            self%v(2:self%NdV(1)-1, 2:self%NdV(2)-1, 1) = cvalue
                        case("z2")
                            self%v(2:self%NdV(1)-1, 2:self%NdV(2)-1, self%NdV(3)) = cvalue
                    end select
                else
                    select case(bdry)
                        case("x1")
                            self%v(1, :, :) = cvalue
                        case("x2")
                            self%v(self%NdV(1), :, :) = cvalue
                        case("y1")
                            self%v(:, 1, :) = cvalue
                        case("y2")
                            self%v(:, self%NdV(2), :) = cvalue
                        case("z1")
                            self%v(:, :, 1) = cvalue
                        case("z2")
                            self%v(:, :, self%NdV(3)) = cvalue
                    end select
                endif
                !
            case( FACE )
                select case(bdry)
                    case("x1")
                        self%v(1, :, :) = cvalue
                    case("x2")
                        self%v(self%NdV(1), :, :) = cvalue
                    case("y1")
                        self%v(:, 1, :) = cvalue
                    case("y2")
                        self%v(:, self%NdV(2), :) = cvalue
                    case("z1")
                        self%v(:, :, 1) = cvalue
                    case("z2")
                        self%v(:, :, self%NdV(3)) = cvalue
                end select
                !
            case default
                call errStop( "setOneBoundary_iScalar3D_SG > Invalid grid type" )
            !
        end select
        !
    end subroutine setOneBoundary_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine intBdryIndices_iScalar3D_SG( self, ind_i, ind_b )
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
            call errStop( "intBdryIndices_iScalar3D_SG > Not allocated. Exiting." )
        endif
        !
        select case( self%grid_type )
            !
            case( NODE )
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
                call errStop( "intBdryIndices_iScalar3D_SG > Unknown self%grid_type" )
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
    end subroutine intBdryIndices_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setVecComponents_iScalar3D_SG( self, xyz, &
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
        if( .NOT. self%is_allocated ) then
            call errStop( "setVecComponents_iScalar3D_SG > self not allocated." )
        endif
        !
        call self%switchStoreState( compound )
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
    end subroutine setVecComponents_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine zeros_iScalar3D_SG( self )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "zeros_iScalar3D_SG > self not allocated." )
        endif
        !
        if( self%store_state .EQ. compound ) then
            !
            self%v = R_ZERO
            !
        elseif( self%store_state .EQ. singleton ) then
            !
            self%s_v = R_ZERO
            !
        else
            call errStop( "zeros_iScalar3D_SG > Unknown store_state!" )
        endif
        !
    end subroutine zeros_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine conjugate_iScalar3D_SG( self )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        !
        call errStop( "conjugate_iScalar3D_SG: do not try to conjugate a integer scalar!" )
        !
    end subroutine conjugate_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine add_iScalar3D_SG( self, rhs )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "add_iScalar3D_SG > self not allocated." )
        endif
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "add_iScalar3D_SG > rhs not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            call self%switchStoreState( rhs%store_state )
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
                        self%s_v = self%s_v + rhs%s_v
                        !
                    else
                        call errStop( "add_iScalar3D_SG > Unknown rhs store_state!" )
                    endif
                    !
                class default
                    call errStop( "add_iScalar3D_SG: undefined rhs" )
                !
            end select
            !
        else
            call errStop( "add_iScalar3D_SG > Incompatible inputs." )
        endif
        !
    end subroutine add_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine linComb_iScalar3D_SG( self, rhs, c1, c2 )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        complex( kind=prec ), intent( in ) :: c1, c2
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "linComb_iScalar3D_SG > self not allocated." )
        endif
        !
        !>  linear combination, in place: self = c1*self+c2*rhs
        if( self%isCompatible(rhs)) then
            !
            call self%switchStoreState( rhs%store_state )
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
                        self%s_v = c1 * self%s_v + c2 * rhs%s_v
                        !
                    else
                        call errStop( "linComb_iScalar3D_SG > Unknown rhs store_state!" )
                    endif
                    !
                class default
                    call errStop( "linComb_iScalar3D_SG: undefined rhs" )
                !
            end select
        else
            call errStop( "linComb_iScalar3D_SG > Incompatible rhs" )
        endif
        !
    end subroutine linComb_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine subValue_iScalar3D_SG( self, cvalue )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "subValue_iScalar3D_SG > self not allocated." )
        endif
        !
        if( self%store_state .EQ. compound ) then
            !
            self%v = self%v - cvalue
            !
        elseif( self%store_state .EQ. singleton ) then
            !
            self%s_v = self%s_v - cvalue
            !
        else
            call errStop( "subValue_iScalar3D_SG > Unknown self store_state!" )
        endif
        !
    end subroutine subValue_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine subField_iScalar3D_SG( self, rhs )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "subField_iScalar3D_SG > self not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            call self%switchStoreState( rhs%store_state )
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
                        self%s_v = self%s_v - rhs%s_v
                        !
                    else
                        call errStop( "subField_iScalar3D_SG > Unknown rhs store_state!" )
                    endif
                    !
                class default
                    call errStop( "subField_iScalar3D_SG: undefined rhs" )
                !
            end select
            !
        else
            call errStop( "subField_iScalar3D_SG > Incompatible inputs." )
        endif
        !
    end subroutine subField_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multByReal_iScalar3D_SG( self, rvalue )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: rvalue
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "multByReal_iScalar3D_SG > self not allocated." )
        endif
        !
        if( self%store_state .EQ. compound ) then
            !
            self%v = self%v * rvalue
            !
        elseif( self%store_state .EQ. singleton ) then
            !
            self%s_v = self%s_v * rvalue
            !
        else
            call errStop( "multByReal_iScalar3D_SG > Unknown self store_state!" )
        endif
        !
    end subroutine multByReal_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multByComplex_iScalar3D_SG( self, cvalue )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "multByComplex_iScalar3D_SG > self not allocated." )
        endif
        !
        if( self%store_state .EQ. compound ) then
            !
            self%v = self%v * cvalue
            !
        elseif( self%store_state .EQ. singleton ) then
            !
            self%s_v = self%s_v * cvalue
            !
        else
            call errStop( "multByComplex_iScalar3D_SG > Unknown self store_state!" )
        endif
        !
    end subroutine multByComplex_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multByField_iScalar3D_SG( self, rhs )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "multByField_iScalar3D_SG > self not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            call self%switchStoreState( rhs%store_state )
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
                        self%s_v = self%s_v * rhs%s_v
                        !
                    else
                        call errStop( "multByField_iScalar3D_SG > Unknown rhs store_state!" )
                    endif
                    !
                class default
                    call errStop( "multByField_iScalar3D_SG > undefined rhs" )
                !
            end select
            !
        else
            call errStop( "multByField_iScalar3D_SG > incompatible rhs" )
        endif
        !
    end subroutine multByField_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multAdd_iScalar3D_SG( self, cvalue, rhs )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        class( Field_t ), intent( in ) :: rhs
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "multAdd_iScalar3D_SG > self not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            call self%switchStoreState( rhs%store_state )
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
                        self%s_v = self%s_v + cvalue * rhs%s_v
                        !
                    else
                        call errStop( "multAdd_iScalar3D_SG > Unknown rhs store_state!" )
                    endif
                    !
                class default
                    call errStop( "multAdd_iScalar3D_SG > rhs undefined." )
                    !
            end select
            !
        else
            call errStop( "multAdd_iScalar3D_SG > Incompatible inputs." )
        endif
        !
    end subroutine multAdd_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    function dotProd_iScalar3D_SG( self, rhs ) result( cvalue )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( in ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        complex( kind=prec ) :: cvalue
        !
        type( iScalar3D_SG_t ) :: copy
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "dotProd_iScalar3D_SG > self not allocated." )
        endif
        !
        copy = self
        !
        if( copy%isCompatible( rhs ) ) then
            !
            call copy%switchStoreState( rhs%store_state )
            !
            select type( rhs )
                !
                class is( iScalar3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        cvalue = cmplx( sum( copy%v * rhs%v ), 0.0, kind=prec )
                        !
                    elseif( rhs%store_state .EQ. singleton ) then
                        !
                        cvalue = cmplx( sum( copy%s_v * rhs%s_v ), 0.0, kind=prec )
                        !
                    else
                        call errStop( "dotProd_iScalar3D_SG > Unknown rhs store_state!" )
                    endif
                    !
                class default
                    call errStop( "dotProd_iScalar3D_SG > undefined rhs" )
                !
            end select
            !
        else
            call errStop( "dotProd_iScalar3D_SG > Incompatible rhs" )
        endif
        !
    end function dotProd_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine divByValue_iScalar3D_SG( self, cvalue )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "divByValue_iScalar3D_SG > self not allocated." )
        endif
        !
        if( self%store_state .EQ. compound ) then
            !
            self%v = self%v / cvalue
            !
        elseif( self%store_state .EQ. singleton ) then
            !
            self%s_v = self%s_v / cvalue
            !
        else
            call errStop( "divByValue_iScalar3D_SG > Unknown self store_state!" )
        endif
        !
    end subroutine divByValue_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine divByField_iScalar3D_SG( self, rhs )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "divByField_iScalar3D_SG > self not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            call self%switchStoreState( rhs%store_state )
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
                        self%s_v = self%s_v / rhs%s_v
                        !
                    else
                        call errStop( "divByField_iScalar3D_SG > Unknown rhs store_state!" )
                    endif
                    !
                class default
                    call errStop( "divByField_iScalar3D_SG > undefined rhs" )
                !
            end select
            !
        else
            call errStop( "divByField_iScalar3D_SG: incompatible rhs" )
        endif
        !
    end subroutine divByField_iScalar3D_SG
    !
    !> No function briefing
    !
    function getV_iScalar3D_SG( self ) result( v )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( in ) :: self
        !
        complex( kind=prec ), allocatable :: v(:, :, :)
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "getV_iScalar3D_SG > self not allocated." )
        endif
        !
        if( .NOT. allocated( v ) ) then
            call errStop( "getV_iScalar3D_SG > v not allocated." )
        else
            !
            v = cmplx( self%v, 0.0, kind=prec )
            !
        endif
        !
    end function getV_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setV_iScalar3D_SG( self, v )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), allocatable, intent( in ) :: v(:, :, :)
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "setV_iScalar3D_SG > self not allocated." )
        endif
        !
        if( .NOT. allocated( v ) ) then
            call errStop( "setV_iScalar3D_SG > v not allocated." )
        endif
        !
        call self%switchStoreState( compound )
        !
        if( allocated( self%s_v ) ) deallocate( self%s_v )
        !
        self%v = real( v, kind=prec )
        !
    end subroutine setV_iScalar3D_SG
    !
    !> No function briefing
    !
    function getSV_iScalar3D_SG( self ) result( s_v )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( in ) :: self
        !
        complex( kind=prec ), allocatable :: s_v(:)
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "getSV_iScalar3D_SG > self not allocated." )
        endif
        !
        if( .NOT. allocated( self%s_v ) ) then
            call errStop( "getSV_iScalar3D_SG > self%s_v not allocated." )
        else
            !
            s_v = cmplx( self%s_v, 0.0, kind=prec )
            !
        endif
        !
    end function getSV_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setSV_iScalar3D_SG( self, s_v )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), allocatable, intent( in ) :: s_v(:)
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "setSV_iScalar3D_SG > self not allocated." )
        endif
        !
        if( .NOT. allocated( s_v ) ) then
            call errStop( "setSV_iScalar3D_SG > s_v not allocated." )
        endif
        !
        call self%switchStoreState( singleton )
        !
        self%s_v = real( s_v, kind=prec )
        !
    end subroutine setSV_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    function getArray_iScalar3D_SG( self ) result( array )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( in ) :: self
        complex( kind=prec ), allocatable, dimension(:) :: array
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "getArray_iScalar3D_SG > self not allocated." )
        endif
        !
        if( self%store_state .EQ. compound ) then
            !
            allocate( array( self%length() ) )
            array = (/reshape( cmplx( self%v, 0.0, kind=prec ), (/self%Nxyz, 1/))/)
            !
        elseif( self%store_state .EQ. singleton ) then
            !
            array = self%s_v
            !
        else
            call errStop( "getArray_iScalar3D_SG > Unknown store_state!" )
        endif
        !
    end function getArray_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setArray_iScalar3D_SG( self, array )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), dimension(:), intent( in ) :: array
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "setArray_iScalar3D_SG > self not allocated." )
        endif
        !
        if( self%store_state .EQ. compound ) then
            !
            if( allocated( self%s_v ) ) deallocate( self%s_v )
            !
            self%v = reshape( real( array, kind=prec ), (/self%NdV(1), self%NdV(2), self%NdV(3)/) )
            !
        elseif( self%store_state .EQ. singleton ) then
            !
            if( allocated( self%v ) ) deallocate( self%v )
            !
            self%s_v = array
            !
        else
            call errStop( "setArray_iScalar3D_SG > Unknown store_state!" )
        endif
        !
    end subroutine setArray_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine copyFrom_iScalar3D_SG( self, rhs )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        ! !
        ! if( .NOT. self%is_allocated ) then
            ! call errStop( "copyFrom_iScalar3D_SG > self not allocated." )
        ! endif
        ! !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "copyFrom_iScalar3D_SG > rhs not allocated" )
        endif
        !
        self%grid => rhs%grid
        self%grid_type = rhs%grid_type
        self%nx = rhs%nx
        self%ny = rhs%ny
        self%nz = rhs%nz
        self%store_state = rhs%store_state
        !
        if( allocated( rhs%ind_interior ) ) &
        self%ind_interior = rhs%ind_interior
        !
        if( allocated( rhs%ind_boundaries ) ) &
        self%ind_boundaries = rhs%ind_boundaries
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
                    self%s_v = rhs%s_v
                    !
                else
                    call errStop( "copyFrom_iScalar3D_SG > Unknown store_state!" )
                endif
                !
            class default
                    call errStop( "copyFrom_iScalar3D_SG > Unclassified rhs" )
            !
        end select
        !
        self%is_allocated = .TRUE.
        !
    end subroutine copyFrom_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine read_iScalar3D_SG( self, funit, ftype )
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
        if( .NOT. self%is_allocated ) then
            call errStop( "read_iScalar3D_SG > self not allocated." )
        endif
        !
        call self%switchStoreState( compound )
        !
        if( .NOT. present( ftype ) ) then
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
                        write( *, * ) "     "//achar(27)//"[31m# Error:"//achar(27)//"[0m iScalar3D_SG::read_iScalar3D_SG: "
                        write( *, * ) "      While reading the ", i, "th block."
                        stop
                 elseif( k1 > k2) then
                        write( *, * ) "Warning: iScalar3D_SG::read_iScalar3D_SG: "
                        write( *, * ) "                Block ", i, " will be ignored."
                 endif
                 !
                 do j = Nx, 1, -1
                        read(funit, *, iostat = istat) temp
                        
                        if( istat /= 0) then
                             write( *, * ) "     "//achar(27)//"[31m# Error:"//achar(27)//"[0m iScalar3D_SG::read_iScalar3D_SG: "
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
            stop "Error: read_iScalar3D_SG: unable to open file"
        endif
        !
    end subroutine read_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine write_iScalar3D_SG( self, funit, ftype )
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
        if( .NOT. self%is_allocated ) then
            call errStop( "write_iScalar3D_SG > self not allocated." )
        endif
        !
        call self%switchStoreState( compound )
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
                 write( *, * ) "     "//achar(27)//"[31m# Error:"//achar(27)//"[0m write_iScalar3D_SG > Unable to write vector to unformatted file ", &
                            trim(fname), "."
                 !
                 stop
            elseif( (index(isbinary,"no") > 0 .OR. index(isbinary,"NO") > 0) &
                     .AND.binary) then
                 write( *, * ) "     "//achar(27)//"[31m# Error:"//achar(27)//"[0m write_iScalar3D_SG > Unable to write vector to formatted file ", &
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
                    call errStop( "write_iScalar3D_SG > Failed while writing to file." )
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
            call errStop( "readRVector3D_SG: unable to open file" )
        endif
        !
    end subroutine write_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine print_iScalar3D_SG( self, io_unit, title, append )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        integer, intent( in ), optional :: io_unit
        character(*), intent( in ), optional :: title
        logical, intent( in ), optional :: append
        !
        integer :: ix, iy, iz,funit
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "print_iScalar3D_SG > self not allocated." )
        endif
        !
        call self%switchStoreState( compound )
        !
        if( present( io_unit ) ) then
            funit = io_unit
        else
            funit = 0    !>    usually this will work to write to standard output
        endif
        !
        if( present(title) ) then
            write(funit,*) title
        endif
        !
        write( funit, * ) self%nx, self%ny, self%nz
        !
        write(funit,*) "iScalar3D_SG"
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
    end subroutine print_iScalar3D_SG
    !
end module iScalar3D_SG
