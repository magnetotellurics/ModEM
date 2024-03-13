!
!> Derived class to define a iScalar3D_SG 
!
module iScalar3D_SG
    !
    use Scalar
    !
    type, extends( Scalar_t ) :: iScalar3D_SG_t
        !
        integer, dimension(3) :: NdV
        !
        integer :: Nxyz
        !
        integer( kind=prec ), allocatable, dimension(:,:,:) :: v
        !
        contains
            !
            !> Destructor
            final :: iScalar3D_SG_dtor
            !
            !> Boundary operations
            procedure, public :: setAllBoundary => setAllBoundary_iScalar3D_SG
            procedure, public :: setOneBoundary => setOneBoundary_iScalar3D_SG
            !
            !> Dimensioning operations
            procedure, public :: length => length_iScalar3D_SG
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
            procedure, public :: sumToNode => sumToNode_iScalar3D_SG
            !
            !> Getters & Setters
            !
            procedure, public :: getArray => getArray_iScalar3D_SG
            procedure, public :: setArray => setArray_iScalar3D_SG
            !
            !> Miscellaneous
            procedure, public :: copyFrom => copyFrom_iScalar3D_SG
            !
            procedure, public :: getReal => getReal_iScalar3D_SG
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
        integer :: nx, ny, nz, nzAir
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
             allocate( self%v(nx + 1, ny + 1, nz + 1), stat=status )
             self%NdV = (/self%nx + 1, self%ny + 1, self%nz + 1/)
             !
        elseif( grid_type == CELL ) then
             !
             allocate(self%v(nx, ny, nz), stat=status) 
             self%NdV = (/self%nx, self%ny, self%nz/)
             !
        else
            call errStop( "iScalar3D_SG_ctor > unrecognized grid type: ["//grid_type//"]" )
        endif
        !
        self%is_allocated = self%is_allocated .AND. ( status .EQ. 0 )
        !
        if( self%is_allocated ) then
            !
            self%v = R_ZERO
            !
            self%Nxyz = product( self%NdV )
            !
        else
            call errStop( "iScalar3D_SG_ctor > Unable to allocate rScalar - invalid grid supplied" )
        endif
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
        if( allocated( self%v ) ) deallocate( self%v )
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
        select case( self%grid_type )
            !
            case( NODE, CELL ) 
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
    function length_iScalar3D_SG( self ) result( field_length )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( in ) :: self
        !
        integer :: field_length
        !
        field_length = self%Nxyz
        !
    end function length_iScalar3D_SG
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
        self%v( x1:x2:xstep, y1:y2:ystep, z1:z2:zstep ) = rvalue
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
        self%v = R_ZERO
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
            select type( rhs )
                !
                class is( iScalar3D_SG_t )
                    !
                    self%v = self%v + rhs%v
                    !
                class default
                    call errStop( "add_iScalar3D_SG > rhs must be Scalar (try vec%scl)!" )
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
        if( self%isCompatible( rhs ) ) then
            !
            select type( rhs )
                !
                class is( iScalar3D_SG_t )
                    !
                    self%v = c1 * self%v + c2 * rhs%v
                    !
                class default
                    call errStop( "linComb_iScalar3D_SG > rhs must be Scalar (try vec%scl)!" )
                !
            end select
            !
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
        self%v = self%v - cvalue
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
            select type( rhs )
                !
                class is( iScalar3D_SG_t )
                    !
                    self%v = self%v - rhs%v
                    !
                class default
                    call errStop( "subField_iScalar3D_SG > rhs must be Scalar (try vec%scl)!" )
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
        self%v = self%v * rvalue
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
        self%v = self%v * cvalue
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
            select type( rhs )
                !
                class is( iScalar3D_SG_t )
                    !
                    self%v = self%v * rhs%v
                    !
                class default
                    call errStop( "multByField_iScalar3D_SG > rhs must be Scalar (try vec%scl)!" )
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
            select type( rhs )
                !
                class is( iScalar3D_SG_t )
                    !
                    self%v = self%v + cvalue * rhs%v
                    !
                class default
                    call errStop( "multAdd_iScalar3D_SG > rhs must be Scalar (try vec%scl)!" )
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
        if( .NOT. self%is_allocated ) then
            call errStop( "dotProd_iScalar3D_SG > self not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type( rhs )
                !
                class is( iScalar3D_SG_t )
                    !
                    cvalue = sum( self%v * rhs%v )
                    !
                class default
                    call errStop( "dotProd_iScalar3D_SG > rhs must be Scalar (try vec%scl)!" )
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
        self%v = self%v / cvalue
        !
    end subroutine divByValue_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine sumToNode_iScalar3D_SG( self, node_scalar, interior_only )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( inout ) :: self
        class( Scalar_t ), intent( inout ) :: node_scalar
        logical, intent( in ), optional :: interior_only
        !
        type( iScalar3D_SG_t ) :: temp_node
        integer :: v_xend, v_yend, v_zend
        logical :: is_interior_only
        !
        if( .NOT. self%is_allocated ) then
             call errStop( "sumToNode_iScalar3D_SG > self not allocated." )
        endif
        !
        if( .NOT. node_scalar%is_allocated ) then
             call errStop( "sumToNode_iScalar3D_SG > node_scalar not allocated." )
        endif
        !
        is_interior_only = .FALSE.
        !
        if( present( interior_only ) ) is_interior_only = interior_only
        !
        if( is_interior_only ) then
            call self%setAllBoundary( C_ZERO )
        endif
        !
        temp_node = iScalar3D_SG_t( self%grid, NODE )
        !
        select case( self%grid_type )
            !
            case( CELL )
                !
                v_xend = size( self%v, 1 )
                v_yend = size( self%v, 2 )
                v_zend = size( self%v, 3 )
                !
                !> Interior
                temp_node%v( 2:v_xend-1, 2:v_yend-1, 2:v_zend-1 ) = &
                self%v( 1:v_xend-1, 1:v_yend-1, 1:v_zend-1 ) + &
                self%v( 2:v_xend  , 1:v_yend-1, 1:v_zend-1 ) + &
                self%v( 1:v_xend-1, 2:v_yend  , 1:v_zend-1 ) + &
                self%v( 1:v_xend-1, 1:v_yend-1, 2:v_zend   ) + &
                self%v( 2:v_xend  , 2:v_yend  , 1:v_zend-1 ) + &
                self%v( 2:v_xend  , 1:v_yend-1, 2:v_zend   ) + &
                self%v( 1:v_xend-1, 2:v_yend  , 2:v_zend   ) + &
                self%v( 2:v_xend  , 2:v_yend  , 2:v_zend   )
                !
                node_scalar = temp_node
                !
                !call node_scalar%mult( cmplx( 0.125_prec, 0.0, kind=prec ) )
                !
            case default
                call errStop( "sumToNode_iScalar3D_SG: undefined self%grid_type" )
        end select
        !
    end subroutine sumToNode_iScalar3D_SG
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
            select type( rhs )
                !
                class is( iScalar3D_SG_t )
                    !
                    self%v = self%v / rhs%v
                    !
                class default
                    call errStop( "divByField_iScalar3D_SG > rhs must be Scalar (try vec%scl)!" )
                !
            end select
            !
        else
            call errStop( "divByField_iScalar3D_SG: incompatible rhs" )
        endif
        !
    end subroutine divByField_iScalar3D_SG
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
        allocate( array( self%length() ) )
        !
        array = (/reshape( self%v, (/self%Nxyz, 1/))/)
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
        complex( kind=prec ), allocatable, dimension(:,:,:) :: v
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "setArray_iScalar3D_SG > self not allocated." )
        endif
        !
        v = reshape( array, (/self%NdV(1), self%NdV(2), self%NdV(3)/) )
        !
        self%v = v
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
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "copyFrom_iScalar3D_SG > rhs not allocated" )
        endif
        !
        self%grid => rhs%grid
        self%grid_type = rhs%grid_type
        self%nx = rhs%nx
        self%ny = rhs%ny
        self%nz = rhs%nz
        !
        select type( rhs )
            !
            class is( iScalar3D_SG_t )
                !
                self%NdV = rhs%NdV
                self%Nxyz = rhs%Nxyz
                !
                self%v = rhs%v
                !
                self%is_allocated = .TRUE.
                !
            class default
                call errStop( "copyFrom_iScalar3D_SG > Unclassified rhs" )
            !
        end select
        !
    end subroutine copyFrom_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine getReal_iScalar3D_SG( self, r_field )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( in ) :: self
        class( Field_t ), allocatable, intent( out ) :: r_field
		!
		call errStop( "getReal_iScalar3D_SG > not implemented." )
		!
    end subroutine getReal_iScalar3D_SG
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
        real( kind=prec ), allocatable, dimension(:) :: temp
        logical :: ok, hasname, binary
        character(:), allocatable :: fname, isbinary
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "read_iScalar3D_SG > self not allocated." )
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
                        write( *, * ) "While reading the ", i, "th block."
                        call errStop( "read_iScalar3D_SG." )
                 elseif( k1 > k2) then
                        write( *, * ) "Block ", i, " will be ignored."
                        call errStop( "read_iScalar3D_SG." )
                 endif
                 !
                 do j = Nx, 1, -1
                        read(funit, *, iostat = istat) temp
                        
                        if( istat /= 0) then
                             write( *, * ) "While reading the ", j, "th row in ", i,"th block."
                             call errStop( "read_iScalar3D_SG." )
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
            call errStop( "read_iScalar3D_SG: unable to open file" )
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
                 call errStop( "write_iScalar3D_SG > Unable to write vector to unformatted file ["//trim(fname)//"]." )
            elseif( (index(isbinary,"no") > 0 .OR. index(isbinary,"NO") > 0) &
                     .AND.binary) then
                 call errStop( "write_iScalar3D_SG > Unable to write vector to formatted file ["//trim(fname)//"]." )
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
            call errStop( "write_iScalar3D_SG: unable to open file" )
        endif
        !
    end subroutine write_iScalar3D_SG
    !
    !> No subroutine briefing
    !
    subroutine print_iScalar3D_SG( self, io_unit, title, append )
        implicit none
        !
        class( iScalar3D_SG_t ), intent( in ) :: self
        integer, intent( in ), optional :: io_unit
        character(*), intent( in ), optional :: title
        logical, intent( in ), optional :: append
        !
        integer :: ix, iy, iz, funit
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "print_iScalar3D_SG > self not allocated." )
        endif
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
