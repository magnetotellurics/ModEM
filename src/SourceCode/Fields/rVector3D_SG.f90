!
!> Derived class to define a rVector3D_SG
!
module rVector3D_SG
    !
    use MatUtils
    use Vector
    use cScalar3D_SG
    !
    type, extends( Vector_t ) :: rVector3D_SG_t
        !
        real( kind=prec ), allocatable, dimension(:, :, :) :: x, y, z
        !
        real( kind=prec ), allocatable, dimension(:) :: s_v
        !
        contains
            !
            final :: rVector3D_SG_dtor
            !
            !> Boundary operations
            procedure, public :: setAllBoundary => setAllBoundary_rVector3D_SG
            procedure, public :: setOneBoundary => setOneBoundary_rVector3D_SG
            procedure, public :: intBdryIndices => intBdryIndices_rVector3D_SG
            !
            !> Dimensioning operations
            procedure, public :: setVecComponents => setVecComponents_rVector3D_SG
            !
            !> Arithmetic/algebraic unary operations
            procedure, public :: zeros => zeros_rVector3D_SG
            !
            procedure, public :: sumEdge => sumEdge_rVector3D_SG
            procedure, public :: sumEdgeVTI => sumEdgeVTI_rVector3D_SG
            !
            procedure, public :: avgCell => avgCell_rVector3D_SG
            procedure, public :: avgCellVTI => avgCellVTI_rVector3D_SG
            !
            procedure, public :: conjugate => conjugate_rVector3D_SG
            !
            !> Arithmetic/algebraic binary operations
            procedure, public :: add => add_rVector3D_SG
            !
            procedure, public :: linComb => linComb_rVector3D_SG
            !
            procedure, public :: subValue => subValue_rVector3D_SG
            procedure, public :: subField => subField_rVector3D_SG
            !
            procedure, public :: multByReal => multByReal_rVector3D_SG
            procedure, public :: multByComplex => multByComplex_rVector3D_SG
            procedure, public :: multByField => multByField_rVector3D_SG
            !
            procedure, public :: diagMult => diagMult_rVector3D_SG
            !
            procedure, public :: multAdd => multAdd_rVector3D_SG
            !
            procedure, public :: dotProd => dotProd_rVector3D_SG
            !
            procedure, public :: divByField => divByField_rVector3D_SG
            procedure, public :: divByValue => divByValue_rVector3D_SG
            !
            procedure, public :: interpFunc => interpFunc_rVector3D_SG
            !
            !> Miscellaneous
            procedure, public :: getAxis => getAxis_rVector3D_SG
            !
            procedure, public :: getReal => getReal_rVector3D_SG
            !
            procedure, public :: getX => getX_rVector3D_SG
            procedure, public :: setX => setX_rVector3D_SG
            procedure, public :: getY => getY_rVector3D_SG
            procedure, public :: setY => setY_rVector3D_SG
            procedure, public :: getZ => getZ_rVector3D_SG
            procedure, public :: setZ => setZ_rVector3D_SG
            !
            procedure, public :: getSV => getSV_rVector3D_SG
            procedure, public :: setSV => setSV_rVector3D_SG
            !
            procedure, public :: getArray => getArray_rVector3D_SG
            procedure, public :: setArray => setArray_rVector3D_SG
            procedure, public :: copyFrom => copyFrom_rVector3D_SG
            !
            !> I/O operations
            procedure, public :: read => read_rVector3D_SG
            procedure, public :: write => write_rVector3D_SG
            procedure, public :: print => print_rVector3D_SG
            !
    end type rVector3D_SG_t
    !
    public :: getRvector, setRvector, edgeLength
    !
    interface rVector3D_SG_t
        module procedure rVector3D_SG_ctor
    end interface rVector3D_SG_t
    !
contains
    !
    !> No subroutine briefing
    !
    function rVector3D_SG_ctor( grid, grid_type ) result ( self )
        implicit none
        !
        class( Grid_t ), target, intent( in ) :: grid
        character( len=4 ), intent( in ) :: grid_type
        !
        type( rVector3D_SG_t ) :: self
        !
        integer :: status
        !
        !write( *, * ) "Constructor rVector3D_SG"
        !
        call self%baseInit
        !
        self%grid => grid
        !
        self%nx = self%grid%nx
        self%ny = self%grid%ny
        self%nz = self%grid%nz
        !
        self%grid_type = trim( grid_type )
        self%is_allocated = .FALSE.
        !
        if( self%grid_type == EDGE ) then
            !
            allocate(self%x(self%nx, self%ny + 1, self%nz + 1), STAT = status)
            self%is_allocated = status.EQ.0
            !
            allocate(self%y(self%nx + 1, self%ny, self%nz + 1), STAT = status)
            self%is_allocated = self%is_allocated.AND.(status.EQ.0)
            !
            allocate(self%z(self%nx + 1, self%ny + 1, self%nz), STAT = status)
            self%is_allocated = self%is_allocated.AND.(status.EQ.0)
            !
            self%NdX = (/self%nx, self%ny + 1, self%nz + 1/)
            self%NdY = (/self%nx + 1, self%ny, self%nz + 1/)
            self%NdZ = (/self%nx + 1, self%ny + 1, self%nz/)
            !
        else if(self%grid_type == FACE) then
            !
            allocate(self%x(self%nx + 1, self%ny, self%nz), STAT = status)
            self%is_allocated = status.EQ.0
            !
            allocate(self%y(self%nx, self%ny + 1, self%nz), STAT = status)
            self%is_allocated = self%is_allocated.AND.(status.EQ.0)
            !
            allocate(self%z(self%nx, self%ny, self%nz + 1), STAT = status)
            self%is_allocated = self%is_allocated.AND.(status.EQ.0)
            !
            self%NdX = (/self%nx + 1, self%ny, self%nz/)
            self%NdY = (/self%nx, self%ny + 1, self%nz/)
            self%NdZ = (/self%nx, self%ny, self%nz + 1/)
            !
        else
            call errStop( "rVector3D_SG_ctor > Only EDGE or FACE types allowed." )
        endif
        !
        if(self%is_allocated) then
            self%x = C_ZERO
            self%y = C_ZERO
            self%z = C_ZERO
        else
            call errStop( "rVector3D_SG_ctor > Unable to allocate vector." )
        endif
        !
        self%Nxyz = (/product(self%NdX), product(self%NdY), product(self%NdZ)/)
        !
        call self%setIndexArrays
        call self%zeros
        !
    end function rVector3D_SG_ctor
    !
    !> No subroutine briefing
    !
    subroutine rVector3D_SG_dtor( self )
        implicit none
        !
        type( rVector3D_SG_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor rVector3D_SG"
        !
        call self%baseDealloc
        !
        if( allocated( self%x ) ) deallocate( self%x )
        if( allocated( self%y ) ) deallocate( self%y )
        if( allocated( self%z ) ) deallocate( self%z )
        !
        if( allocated( self%s_v ) ) deallocate( self%s_v )
        !
        self%nx = 0
        self%ny = 0
        self%nz = 0
        !
        self%grid_type = ""
        self%is_allocated = .FALSE.
        !
    end subroutine rVector3D_SG_dtor
    !
    !> No subroutine briefing
    !
    subroutine setAllBoundary_rVector3D_SG( self, cvalue )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        call self%switchStoreState( compound )
        !
        select case( self%grid_type )
            !
            case( EDGE )
                !
                self%x(:, (/1, self%NdX(2)/), :) = cvalue
                self%x(:, :, (/1, self%NdX(3)/)) = cvalue
                self%y((/1, self%NdY(1)/), :, :) = cvalue
                self%y(:, :, (/1, self%NdY(3)/)) = cvalue
                self%z(:, (/1, self%NdZ(2)/), :) = cvalue
                self%z((/1, self%NdZ(1)/), :, :) = cvalue
                !
            case( FACE )
                !
                self%x((/1, self%NdX(1)/), :, :) = cvalue
                self%y(:, (/1, self%NdY(2)/), :) = cvalue
                self%z(:, :, (/1, self%NdZ(3)/)) = cvalue
                !
            case default
                call errStop( "setAllBoundary_rVector3D_SG > Invalid grid type." )
            !
        end select
        !
    end subroutine setAllBoundary_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setOneBoundary_rVector3D_SG( self, bdry, cvalue, int_only )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        character(*), intent( in ) :: bdry
        complex( kind=prec ), intent( in ) :: cvalue
        logical, intent( in ), optional :: int_only
        !
        logical :: int_only_p
        !
        call self%switchStoreState( compound )
        !
        if( .NOT. present( int_only )) then
            int_only_p = .FALSE.
        else
            int_only_p = int_only
        endif
        !
        select case( self%grid_type )
            !
            case( EDGE )
                !
                if( int_only_p ) then
                    !
                    select case( bdry )
                        !
                        case("x1")
                            self%z(1, 2:self%NdZ(2)-1, :) = cvalue
                            self%y(1, :, 2:self%NdY(3)-1) = cvalue
                        case("x2")
                            self%z(self%NdZ(1), 2:self%NdZ(2)-1, :) = cvalue
                            self%y(self%NdY(1), :, 2:self%NdY(3)-1) = cvalue
                        case("y1")
                            self%z(2:self%NdZ(1)-1, 1, :) = cvalue
                            self%x(:, 1, 2:self%NdX(3)-1) = cvalue
                        case("y2")
                            self%z(2:self%NdZ(1)-1, self%NdZ(2), :) = cvalue
                            self%x(:, self%NdX(2), 2:self%NdX(3)-1) = cvalue
                        case("z1")
                            self%x(:, 2:self%NdX(2)-1, 1) = cvalue
                            self%y(2:self%NdY(1)-1, :, 1) = cvalue
                        case("z2")
                            self%x(:, 2:self%NdX(2)-1, self%NdX(3)) = cvalue
                            self%y(2:self%NdY(1)-1, :, self%NdY(3)) = cvalue
                        case("z1_x")
                            self%x(:, 2:self%NdX(2)-1, 1) = cvalue
                        case("z2_x")
                            self%x(:, 2:self%NdX(2)-1, self%NdX(3)) = cvalue
                        case("z1_y")
                            self%y(2:self%NdY(1)-1, :, 1) = cvalue
                        case("z2_y")
                            self%y(2:self%NdY(1)-1, :, self%NdY(3)) = cvalue
                        case default
                            call errStop( "setOneBoundary_rVector3D_SG > Invalid int_only_p bdry." )
                        !
                    end select
                    !
                else
                    !
                    select case( bdry )
                        !
                        case("x1")
                            self%z(1, :, :) = cvalue
                            self%y(1, :, :) = cvalue
                        case("x2")
                            self%z(self%NdZ(1), :, :) = cvalue
                            self%y(self%NdY(1), :, :) = cvalue
                        case("y1")
                            self%z(:, 1, :) = cvalue
                            self%x(:, 1, :) = cvalue
                        case("y2")
                            self%z(:, self%NdZ(2), :) = cvalue
                            self%x(:, self%NdX(2), :) = cvalue
                        case("z1")
                            self%x(:, :, 1) = cvalue
                            self%y(:, :, 1) = cvalue
                        case("z2")
                            self%x(:, :, self%NdX(3)) = cvalue
                            self%y(:, :, self%NdY(3)) = cvalue
                        case("z1_x")
                            self%x(:, :, 1) = cvalue
                        case("z2_x")
                            self%x(:, :, self%NdX(3)) = cvalue
                        case("z1_y")
                            self%y(:, :, 1) = cvalue
                        case("z2_y")
                            self%y(:, :, self%NdY(3)) = cvalue
                        case default
                            call errStop( "setOneBoundary_rVector3D_SG > Invalid bdry." )
                        !
                    end select
                    !
                endif
                !
            case( FACE )
                !
                select case( bdry )
                    !
                    case("x1")
                        self%x(1, :, :) = cvalue
                    case("x2")
                        self%x(self%NdX(1), :, :) = cvalue
                    case("y1")
                        self%y(:, 1, :) = cvalue
                    case("y2")
                        self%y(:, self%NdY(2), :) = cvalue
                    case("z1")
                        self%z(:, :, 1) = cvalue
                    case("z2")
                        self%z(:, :, self%NdZ(3)) = cvalue
                    case default
                        call errStop( "setOneBoundary_rVector3D_SG > Invalid FACE bdry." )
                    !
                end select
                !
            case default
                call errStop( "setOneBoundary_rVector3D_SG > Invalid grid type." )
        end select
        !
    end subroutine setOneBoundary_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine intBdryIndices_rVector3D_SG( self, ind_i, ind_b )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        integer, allocatable, intent( out ) :: ind_i(:), ind_b(:)
        !
        integer :: nVecT, nBdry, nb, ni, i
        real( kind=prec ), dimension(:), allocatable :: temp
        type( rVector3D_SG_t ) :: E
        !
        if( self%is_allocated ) then
            !
            E = rVector3D_SG_t( self%grid, self%grid_type )
            !
        else
            call errStop( "intBdryIndices_rVector3D_SG > Not allocated. Exiting." )
        endif
        !
        select case( self%grid_type )
            !
            case( EDGE )
                !
                E%x(:, 1, :) = 1
                E%x(:, E%ny + 1, :) = 1
                E%x(:, :, 1) = 1
                E%x(:, :, E%nz + 1) = 1
                E%y(1, :, :) = 1
                E%y(E%nx + 1, :, :) = 1
                E%y(:, :, 1) = 1
                E%y(:, :, E%nz + 1) = 1
                E%z(1, :, :) = 1
                E%z(E%nx + 1, :, :) = 1
                E%z(:, 1, :) = 1
                E%z(:, E%ny + 1, :) = 1
                !
            case( FACE )
                !
                E%x(1, :, :) = 1
                E%x(E%nx + 1, :, :) = 1
                E%y(:, 1, :) = 1
                E%y(:, E%ny + 1, :) = 1
                E%z(:, :, 1) = 1
                E%z(:, :, E%nz + 1) = 1
                !
            case default
                call errStop( "intBdryIndices_rVector3D_SG > Undefined self%grid_type" )
                !
        end select
        !
        temp = E%getArray()
        !
        nVecT = size( E%x ) + size( E%y ) + size( E%z )
        nBdry = 0
        do i = 1, nVecT
            nBdry = nBdry + nint( real( temp(i) ) )
        enddo
        !
        if( allocated( ind_i ) ) deallocate( ind_i )
        allocate( ind_i( nVecT - nBdry ) )
        !
        if( allocated( ind_b ) ) deallocate( ind_b )
        allocate( ind_b( nBdry ) )
        !
        nb = 0
        ni = 0
        do i = 1, nVecT
            if( nint( real( temp(i) ) ) .EQ. 1 ) then
                nb = nb + 1
                ind_b(nb) = i
            else
                ni = ni + 1
                ind_i(ni) = i
            endif
        enddo
        !
        deallocate( temp )
        !
    end subroutine intBdryIndices_rVector3D_SG
    !
    subroutine setVecComponents_rVector3D_SG( self, xyz, &
            &                                 xmin, xstep, xmax, &
            &                                 ymin, ystep, ymax, &
            &                                 zmin, zstep, zmax, cvalue )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        character, intent( in ) :: xyz
        integer, intent( in ) :: xmin, xstep, xmax
        integer, intent( in ) :: ymin, ystep, ymax
        integer, intent( in ) :: zmin, zstep, zmax
        complex( kind=prec ), intent ( in ) :: cvalue
        !
        integer :: x1, x2
        integer :: y1, y2
        integer :: z1, z2
        !
        call self%switchStoreState( compound )
        !
        x1 = xmin; x2 = xmax
        y1 = ymin; y2 = ymax
        z1 = zmin; z2 = zmax
        !
        select case( xyz )
            !
            case( "x" )
                !
                if(xmin == 0) x1 = self%NdX(1)
                if(xmax <= 0) x2 = self%NdX(1) + xmax
                !
                if(ymin == 0) y1 = self%NdX(2)
                if(ymax <= 0) y2 = self%NdX(2) + ymax
                !
                if(zmin == 0) z1 = self%NdX(3)
                if(zmax <= 0) z2 = self%NdX(3) + zmax
                !
                self%x(x1:x2:xstep, y1:y2:ystep, z1:z2:zstep) = cvalue
                !
            case( "y" )
                !
                if(xmin == 0) x1 = self%NdY(1)
                if(xmax <= 0) x2 = self%NdY(1) + xmax
                !
                if(ymin == 0) y1 = self%NdY(2)
                if(ymax <= 0) y2 = self%NdY(2) + ymax
                !
                if(zmin == 0) z1 = self%NdY(3)
                if(zmax <= 0) z2 = self%NdY(3) + zmax
                !
                self%y(x1:x2:xstep, y1:y2:ystep, z1:z2:zstep) = cvalue
                !
            case( "z" )
                !
                if(xmin == 0) x1 = self%NdZ(1)
                if(xmax <= 0) x2 = self%NdZ(1) + xmax
                !
                if(ymin == 0) y1 = self%NdZ(2)
                if(ymax <= 0) y2 = self%NdZ(2) + ymax
                !
                if(zmin == 0) z1 = self%NdZ(3)
                if(zmax <= 0) z2 = self%NdZ(3) + zmax
                !
                self%z(x1:x2:xstep, y1:y2:ystep, z1:z2:zstep) = cvalue
                !
            case default
                call errStop( "setVecComponents_rVector3D_SG > Invalid xyz argument." )
        end select
        !
    end subroutine setVecComponents_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine zeros_rVector3D_SG( self )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        !
        if( .NOT. self%is_allocated ) then
             call errStop( "zeros_rVector3D_SG > self not allocated." )
        endif
        !
        if( self%store_state .EQ. compound ) then
            !
            self%x = R_ZERO
            self%y = R_ZERO
            self%z = R_ZERO
            !
        else if( self%store_state .EQ. singleton ) then
            !
            self%s_v = R_ZERO
            !
        else
            call errStop( "zeros_rVector3D_SG > Unknown store_state!" )
        endif
        !
    end subroutine zeros_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine sumEdge_rVector3D_SG( self, cell_out, interior_only )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Scalar_t ), allocatable, intent( out ) :: cell_out
        logical, intent( in ), optional :: interior_only
        !
        integer :: x_xend, x_yend, x_zend
        integer :: y_xend, y_yend, y_zend
        integer :: z_xend, z_yend, z_zend
        logical :: is_interior_only
        !
        if( .NOT. self%is_allocated ) then
             call errStop( "sumEdge_rVector3D_SG > self not allocated." )
        endif
        !
        call self%switchStoreState( compound )
        !
        is_interior_only = .FALSE.
        !
        if( present( interior_only ) ) is_interior_only = interior_only
        !
        if( is_interior_only ) then
            call self%setAllBoundary( C_ZERO )
        endif
        !
        allocate( cell_out, source = rScalar3D_SG_t( self%grid, CELL ) )
        !
        select type( cell_out )
            !
            class is( rScalar3D_SG_t )
                !
                select case( self%grid_type )
                    !
                    case( EDGE )
                        !
                        x_xend = size(self%x, 1)
                        x_yend = size(self%x, 2)
                        x_zend = size(self%x, 3)
                        !
                        y_xend = size(self%y, 1)
                        y_yend = size(self%y, 2)
                        y_zend = size(self%y, 3)
                        !
                        z_xend = size(self%z, 1)
                        z_yend = size(self%z, 2)
                        z_zend = size(self%z, 3)
                        !
                        cell_out%v = self%x(:,1:x_yend-1,1:x_zend-1) + &
                        self%x(:,2:x_yend,1:x_zend-1)       + &
                        self%x(:,1:x_yend-1,2:x_zend)       + &
                        self%x(:,2:x_yend,2:x_zend)         + &
                        self%y(1:y_xend-1,:,1:y_zend-1)     + &
                        self%y(2:y_xend,:,1:y_zend-1)       + &
                        self%y(1:y_xend-1,:,2:y_zend)       + &
                        self%y(2:y_xend,:,2:y_zend)         + &
                        self%z(1:z_xend-1,1:z_yend-1,:)     + &
                        self%z(2:z_xend,1:z_yend-1,:)       + &
                        self%z(1:z_xend-1,2:z_yend,:)       + &
                        self%z(2:z_xend,2:z_yend,:)
                        !
                    case( FACE )
                        !
                        x_xend = size(self%x, 1)
                        y_xend = size(self%y, 1)
                        z_xend = size(self%z, 1)
                        !
                        cell_out%v = self%x(1:x_xend-1,:,:) + self%x(2:x_xend,:,:) + &
                        self%y(:,1:y_yend-1,:) + self%y(:,2:y_yend,:) + &
                        self%z(:,:,1:z_zend-1) + self%z(:,:,2:z_zend)
                        !
                    case default
                        call errStop( "sumEdge_rVector3D_SG: undefined self%grid_type" )
                end select
                !
            class default
                call errStop( "sumEdge_rVector3D_SG: Unclassified cell_out" )
            !
        end select
        !
    end subroutine sumEdge_rVector3D_SG
    !
    subroutine sumEdgeVTI_rVector3D_SG( self, cell_h_out, cell_v_out, interior_only )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Scalar_t ), allocatable, intent( out ) :: cell_h_out, cell_v_out
        logical, optional, intent( in ) :: interior_only
        !
        complex( kind=prec ), allocatable :: sigma_v(:, :, :)
        integer :: x_xend, x_yend, x_zend
        integer :: y_xend, y_yend, y_zend
        integer :: z_xend, z_yend, z_zend
        logical :: is_interior_only
        !
        if( .NOT. self%is_allocated ) then
             call errStop( "sumEdgeVTI_rVector3D_SG > self not allocated." )
        endif
        !
        call self%switchStoreState( compound )
        !
        is_interior_only = .FALSE.
        !
        if( present( interior_only ) ) is_interior_only = interior_only
        !
        if( is_interior_only ) then
            call self%setAllBoundary( C_ZERO )
        endif
        !
        allocate( cell_h_out, source = rScalar3D_SG_t( self%grid, CELL ) )
        !
        allocate( cell_v_out, source = rScalar3D_SG_t( self%grid, CELL ) )
        !
        select type( cell_h_out )
            !
            class is( cScalar3D_SG_t )
                !
                select type( cell_v_out )
                    !
                    class is( cScalar3D_SG_t )
                        !
                        select case( self%grid_type )
                            !
                            case( EDGE )
                                !
                                x_xend = size( self%x, 1 )
                                x_yend = size( self%x, 2 )
                                x_zend = size( self%x, 3 )
                                !
                                y_xend = size( self%y, 1 )
                                y_yend = size( self%y, 2 )
                                y_zend = size( self%y, 3 )
                                !
                                z_xend = size( self%z, 1 )
                                z_yend = size( self%z, 2 )
                                z_zend = size( self%z, 3 )
                                !
                                sigma_v = self%x(:,1:x_yend-1,1:x_zend-1) + &
                                self%x(:,2:x_yend,1:x_zend-1)       + &
                                self%x(:,1:x_yend-1,2:x_zend)       + &
                                self%x(:,2:x_yend,2:x_zend)         + &
                                self%y(1:y_xend-1,:,1:y_zend-1)     + &
                                self%y(2:y_xend,:,1:y_zend-1)       + &
                                self%y(1:y_xend-1,:,2:y_zend)       + &
                                self%y(2:y_xend,:,2:y_zend)
                                !
                                call cell_h_out%setV( sigma_v )
                                !
                                sigma_v = self%z(1:z_xend-1,1:z_yend-1,:)     + &
                                self%z(2:z_xend,1:z_yend-1,:)       + &
                                self%z(1:z_xend-1,2:z_yend,:)       + &
                                self%z(2:z_xend,2:z_yend,:)
                                !
                                call cell_v_out%setV( sigma_v )
                                !
                            case( FACE )
                                !
                                x_xend = size( self%x, 1 )
                                y_xend = size( self%y, 1 )
                                z_xend = size( self%z, 1 )
                                !
                                sigma_v = self%x(1:x_xend-1,:,:) + self%x(2:x_xend,:,:) + &
                                          self%y(:,1:y_yend-1,:) + self%y(:,2:y_yend,:)
                                !
                                call cell_h_out%setV( sigma_v )
                                !
                                sigma_v = self%z(:,:,1:z_zend-1) + self%z(:,:,2:z_zend)
                                !
                                call cell_v_out%setV( sigma_v )
                                !
                            case default
                                call errStop( "sumEdgeVTI_rVector3D_SG: undefined self%grid_type" )
                            !
                        end select
                        !
                    class default
                        call errStop( "sumEdgeVTI_rVector3D_SG: Unclassified cell_v_out" )
                end select
                !
            class default
                call errStop( "sumEdgeVTI_rVector3D_SG: Unclassified cell_h_out" )
            !
        end select
        !
    end subroutine sumEdgeVTI_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine avgCell_rVector3D_SG( self, cell_in, ptype )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Scalar_t ), intent( inout ) :: cell_in
        character(*), intent( in ), optional :: ptype
        !
        complex( kind=prec ), allocatable :: cell_in_v(:, :, :)
        character(10) :: grid_type
        integer :: xend, yend, zend
        integer :: v_xend, v_yend, v_zend
        integer :: ix, iy, iz
        !
        if( .NOT. self%is_allocated ) then
             call errStop( "avgCell_rVector3D_SG > self not allocated." )
        endif
        !
        if( .NOT. cell_in%is_allocated ) then
             call errStop( "avgCell_rVector3D_SG > cell_in not allocated." )
        endif
        !
        if( index( self%grid_type, CELL ) > 0 ) then
            call errStop( "avgCell_rVector3D_SG > Only CELL type supported." )
        endif
        !
        if( .NOT. present( ptype ) ) then
            grid_type = EDGE
        else
            grid_type = ptype
        endif
        !
        cell_in_v = cell_in%getV()
        !
        call self%switchStoreState( compound )
        !
        v_xend = size( cell_in_v, 1 )
        v_yend = size( cell_in_v, 2 )
        v_zend = size( cell_in_v, 3 )
        !
        select case( grid_type )
            !
            case( EDGE )
                !
                !> for x-components inside the domain
                do ix = 1, self%grid%nx
                    do iy = 2, self%grid%ny
                        do iz = 2, self%grid%nz
                            self%x(ix, iy, iz) = (cell_in_v(ix, iy-1, iz-1) + cell_in_v(ix, iy, iz-1) + &
                            cell_in_v(ix, iy-1, iz) + cell_in_v(ix, iy, iz))/4.0d0
                        enddo
                    enddo
                enddo
                !
                !> for y-components inside the domain
                do ix = 2, self%grid%nx
                    do iy = 1, self%grid%ny
                        do iz = 2, self%grid%nz
                            self%y(ix, iy, iz) = (cell_in_v(ix-1, iy, iz-1) + cell_in_v(ix, iy, iz-1) + &
                            cell_in_v(ix-1, iy, iz) + cell_in_v(ix, iy, iz))/4.0d0
                        enddo
                    enddo
                enddo
                !
                !> for z-components inside the domain
                do ix = 2, self%grid%nx
                    do iy = 2, self%grid%ny
                        do iz = 1, self%grid%nz
                            self%z(ix, iy, iz) = (cell_in_v(ix-1, iy-1, iz) + cell_in_v(ix-1, iy, iz) + &
                            cell_in_v(ix, iy-1, iz) + cell_in_v(ix, iy, iz))/4.0d0
                        enddo
                    enddo
                enddo
                !
            case( FACE )
                !
                xend = size(self%x, 1)
                self%x(2:xend-1,:,:) = cell_in_v(1:v_xend-1,:,:) + cell_in_v(2:v_xend,:,:)
                !
                yend = size(self%y, 1)
                self%y(:, 2:yend-1, :) = cell_in_v(:, 1:v_yend-1, :) + cell_in_v(:, 2:v_yend, :)
                !
                zend = size(self%z, 1) 
                self%z(:, :, 2:zend-1) = cell_in_v(:, :, 1:v_zend-1) + cell_in_v(:, :, 2:v_zend)
                !
            case default
                call errStop( "avgCell_rVector3D_SG: Unknown type" )
            !
        end select !type
        !
    end subroutine avgCell_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine avgCellVTI_rVector3D_SG( self, cell_h_in, cell_v_in, ptype )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Scalar_t ), intent( inout ) :: cell_h_in, cell_v_in
        character(*), intent( in ), optional :: ptype
        !
        complex( kind=prec ), allocatable :: v_h(:, :, :), v_v(:, :, :)
        character(10) :: grid_type
        integer :: xend, yend, zend
        integer :: v_xend, v_yend, v_zend
        integer :: ix, iy, iz
        !
        if( .NOT. self%is_allocated ) then
             call errStop( "avgCellVTI_rVector3D_SG > self not allocated." )
        endif
        !
        if( .NOT. cell_h_in%is_allocated ) then
             call errStop( "avgCellVTI_rVector3D_SG > cell_h_in not allocated." )
        endif
        !
        if( .NOT. cell_v_in%is_allocated ) then
             call errStop( "avgCellVTI_rVector3D_SG > cell_v_in not allocated." )
        endif
        !
        if( index( self%grid_type, CELL ) > 0 ) then
            call errStop( "avgCellVTI_rVector3D_SG > Only CELL type supported." )
        endif
        !
        if( .NOT. present( ptype ) ) then
            grid_type = EDGE
        else
            grid_type = ptype
        endif
        !
        call self%switchStoreState( compound )
        !
        v_h = cell_h_in%getV()
        v_v = cell_v_in%getV()
        !
        v_xend = size( v_h, 1 )
        v_yend = size( v_h, 2 )
        v_zend = size( v_v, 3 )
        !
        select case( grid_type )
            !
            case( EDGE )
                !
                !> for x-components inside the domain
                do ix = 1, self%grid%nx
                    do iy = 2, self%grid%ny
                        do iz = 2, self%grid%nz
                            self%x(ix, iy, iz) = ( v_h(ix, iy-1, iz-1) + v_h(ix, iy, iz-1) + &
                            v_h(ix, iy-1, iz) + v_h(ix, iy, iz) ) / 4.0d0
                        enddo
                    enddo
                enddo
                !
                !> for y-components inside the domain
                do ix = 2, self%grid%nx
                    do iy = 1, self%grid%ny
                        do iz = 2, self%grid%nz
                            self%y(ix, iy, iz) = ( v_h(ix-1, iy, iz-1) + v_h(ix, iy, iz-1) + &
                            v_h(ix-1, iy, iz) + v_h(ix, iy, iz) ) / 4.0d0
                        enddo
                    enddo
                enddo
                !
                !> for z-components inside the domain
                do ix = 2, self%grid%nx
                    do iy = 2, self%grid%ny
                        do iz = 1, self%grid%nz
                            self%z(ix, iy, iz) = ( v_v(ix-1, iy-1, iz) + v_v(ix-1, iy, iz) + &
                            v_v(ix, iy-1, iz) + v_v(ix, iy, iz) ) / 4.0d0
                        enddo
                    enddo
                enddo
                !
            case( FACE )
                !
                xend = size(self%x, 1)
                self%x(2:xend-1,:,:) = v_h(1:v_xend-1,:,:) + v_h(2:v_xend,:,:)
                !
                yend = size(self%y, 1)
                self%y(:, 2:yend-1, :) = v_h(:, 1:v_yend-1, :) + v_h(:, 2:v_yend, :)
                !
                zend = size(self%z, 1) 
                self%z(:, :, 2:zend-1) = v_v(:, :, 1:v_zend-1) + v_v(:, :, 2:v_zend)
                !
            case default
                call errStop( "avgCellVTI_rVector3D_SG: Unknown type" )
            !
        end select !grid_type
        !
    end subroutine avgCellVTI_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine conjugate_rVector3D_SG( self )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        !
        call errStop( "conjugate_rVector3D_SG: Do not try to conjugate a real vector!" )
        !
    end subroutine conjugate_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine add_rVector3D_SG( self, rhs )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( ( .NOT. self%is_allocated ) .OR. ( .NOT. rhs%is_allocated ) ) then
            call errStop( "add_rVector3D_SG > Input vectors not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            call self%switchStoreState( rhs%store_state )
            !
            select type( rhs )
                !
                class is( rVector3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%x = self%x + rhs%x
                        self%y = self%y + rhs%y
                        self%z = self%z + rhs%z
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%s_v = self%s_v + rhs%s_v
                        !
                    else
                        call errStop( "add_rVector3D_SG > Unknown rVector3D_SG_t rhs store_state!" )
                    endif
                    !
                class is( cScalar3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%x = self%x + rhs%v
                        self%y = self%y + rhs%v
                        self%z = self%z + rhs%v
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%s_v = self%s_v + rhs%s_v
                        !
                    else
                        call errStop( "add_rVector3D_SG > Unknown cScalar3D_SG_t rhs store_state!" )
                    endif
                    !
                class is( rScalar3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%x = self%x + rhs%v
                        self%y = self%y + rhs%v
                        self%z = self%z + rhs%v
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%s_v = self%s_v + rhs%s_v
                        !
                    else
                        call errStop( "add_rVector3D_SG > Unknown rScalar3D_SG_t rhs store_state!" )
                    endif
                    !
                class default
                    call errStop( "add_rVector3D_SG > Undefined rhs" )
                !
            end select
            !
        else
            call errStop( "add_rVector3D_SG > Incompatible inputs." )
        endif
        !
    end subroutine add_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine linComb_rVector3D_SG( self, rhs, c1, c2 )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        complex( kind=prec ), intent( in ) :: c1, c2
        !
        if( ( .NOT. self%is_allocated ) .OR. ( .NOT. rhs%is_allocated ) ) then
            call errStop( "linComb_rVector3D_SG > Input vectors not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            call self%switchStoreState( rhs%store_state )
            !
            select type( rhs )
                !
                class is( rVector3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%x = c1 * self%x + c2 * rhs%x
                        self%y = c1 * self%y + c2 * rhs%y
                        self%z = c1 * self%z + c2 * rhs%z
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%s_v = c1 * self%s_v + c2 * rhs%s_v
                        !
                    else
                        call errStop( "linComb_rVector3D_SG > Unknown rhs store_state!" )
                    endif
                    !
                class default
                    call errStop( "linComb_rVector3D_SG > rhs undefined." )
            end select
            !
        else
            call errStop( "linComb_rVector3D_SG > Incompatible inputs." )
        endif
        !
    end subroutine linComb_rVector3D_SG
    !
    !> No subroutine briefing
    subroutine subValue_rVector3D_SG( self, cvalue )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        if( ( .NOT. self%is_allocated ) ) then
            call errStop( "subValue_rVector3D_SG > Self not allocated." )
        endif
        !
        if( self%store_state .EQ. compound ) then
            !
            self%x = self%x - cvalue
            self%y = self%y - cvalue
            self%z = self%z - cvalue
            !
        else if( self%store_state .EQ. singleton ) then
            !
            self%s_v = self%s_v - cvalue
            !
        else
            call errStop( "subValue_rVector3D_SG > Unknown self store_state!" )
        endif
        !
    end subroutine subValue_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine subField_rVector3D_SG( self, rhs )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( ( .NOT. self%is_allocated ) .OR. ( .NOT. rhs%is_allocated ) ) then
            call errStop( "subField_rVector3D_SG > Input vectors not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            call self%switchStoreState( rhs%store_state )
            !
            select type( rhs )
                !
                class is( rVector3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%x = self%x - rhs%x
                        self%y = self%y - rhs%y
                        self%z = self%z - rhs%z
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%s_v = self%s_v - rhs%s_v
                        !
                    else
                        stop "Error: subField_rVector3D_SG > Unknown rhs store_state!"
                    endif
                    !
                class is( cScalar3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%x = self%x - rhs%v
                        self%y = self%y - rhs%v
                        self%z = self%z - rhs%v
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%s_v = self%s_v - rhs%s_v
                        !
                    else
                        call errStop( "subField_rVector3D_SG > Unknown cScalar3D_SG_t rhs store_state!" )
                    endif
                    !
                class is( rScalar3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%x = self%x - rhs%v
                        self%y = self%y - rhs%v
                        self%z = self%z - rhs%v
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%s_v = self%s_v - rhs%s_v
                        !
                    else
                        call errStop( "subField_rVector3D_SG > Unknown rhs store_state!" )
                    endif
                    !
                class default
                    call errStop( "subField_rVector3D_SG > Undefined rhs" )
                !
            end select
            !
        else
            call errStop( "subField_rVector3D_SG > Incompatible inputs." )
        endif
        !
    end subroutine subField_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multByReal_rVector3D_SG( self, rvalue )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: rvalue
        !
        if( self%store_state .EQ. compound ) then
            !
            self%x = self%x * rvalue
            self%y = self%y * rvalue
            self%z = self%z * rvalue
            !
        else if( self%store_state .EQ. singleton ) then
            !
            self%s_v = self%s_v * rvalue
            !
        else
            call errStop( "multByReal_rVector3D_SG > Unknown self store_state!" )
        endif
        !
    end subroutine multByReal_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multByComplex_rVector3D_SG( self, cvalue )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        if( self%store_state .EQ. compound ) then
            !
            self%x = self%x * cvalue
            self%y = self%y * cvalue
            self%z = self%z * cvalue
            !
        else if( self%store_state .EQ. singleton ) then
            !
            self%s_v = self%s_v * cvalue
            !
        else
            call errStop( "multByComplex_rVector3D_SG > Unknown self store_state!" )
        endif
        !
    end subroutine multByComplex_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multByField_rVector3D_SG( self, rhs )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( ( .NOT. self%is_allocated ) .OR. ( .NOT. rhs%is_allocated ) ) then
            call errStop( "multByField_rVector3D_SG > Input vectors not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            call self%switchStoreState( rhs%store_state )
            !
            select type( rhs )
                !
                class is( rVector3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%x = self%x * rhs%x
                        self%y = self%y * rhs%y
                        self%z = self%z * rhs%z
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%s_v = self%s_v * rhs%s_v
                        !
                    else
                        call errStop( "multByField_rVector3D_SG > Unknown rVector3D_SG_t rhs store_state!" )
                    endif
                    !
                class is( rScalar3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%x = self%x * rhs%v
                        self%y = self%y * rhs%v
                        self%z = self%z * rhs%v
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%s_v = self%s_v * rhs%s_v
                        !
                    else
                        call errStop( "multByField_rVector3D_SG > Unknown rScalar3D_SG_t rhs store_state!" )
                    endif
                    !
                class is( cScalar3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%x = self%x * rhs%v
                        self%y = self%y * rhs%v
                        self%z = self%z * rhs%v
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%s_v = self%s_v * rhs%s_v
                        !
                    else
                        call errStop( "multByField_rVector3D_SG > Unknown cScalar3D_SG_t rhs store_state!" )
                    endif
                    !
                class default
                    call errStop( "multByField_rVector3D_SG: undefined rhs" )
                !
            end select
            !
        else
            call errStop( "multByField_rVector3D_SG: incompatible rhs" )
        endif
        !
    end subroutine multByField_rVector3D_SG
    !
    !> No subroutine briefing
    !
    function diagMult_rVector3D_SG( self, rhs ) result( diag_mult )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Vector_t ), intent( in ) :: rhs
        !
        class( Vector_t ), allocatable :: diag_mult
        !
        if( ( .NOT. self%is_allocated ) .OR. ( .NOT. rhs%is_allocated ) ) then
            call errStop( "diagMult_rVector3D_SG > Input vectors not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            allocate( diag_mult, source = rVector3D_SG_t( self%grid, self%grid_type ) )
            !
            call diag_mult%switchStoreState( rhs%store_state )
            !
            call self%switchStoreState( rhs%store_state )
            !
            select type( diag_mult )
                !
                class is( rVector3D_SG_t )
                    !
                    select type( rhs )
                        !
                        class is( rVector3D_SG_t )
                            !
                            if( rhs%store_state .EQ. compound ) then
                                !
                                diag_mult%x = self%x * rhs%x
                                diag_mult%y = self%y * rhs%y
                                diag_mult%z = self%z * rhs%z
                                !
                            else if( rhs%store_state .EQ. singleton ) then
                                !
                                diag_mult%s_v = self%s_v * rhs%s_v
                                !
                            else
                                call errStop( "diagMult_rVector3D_SG > Unknown rVector3D_SG_t rhs store_state!" )
                            endif
                            !
                        class default
                            call errStop( "diagMult_rVector3D_SG > Undefined rhs" )
                        !
                    end select
                !
                class default
                    call errStop( "diagMult_rVector3D_SG > Undefined diag_mult" )
                !
            end select
        else
            call errStop( "diagMult_rVector3D_SG > Incompatible inputs." )
        endif
        !
    end function diagMult_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multAdd_rVector3D_SG( self, cvalue, rhs )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        class( Field_t ), intent( in ) :: rhs
        !
        if( ( .NOT. self%is_allocated ) .OR. ( .NOT. rhs%is_allocated ) ) then
            call errStop( "multAdd_rVector3D_SG > Input vectors not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            call self%switchStoreState( rhs%store_state )
            !
            select type( rhs )
                !
                class is( rVector3D_SG_t ) 
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%x = self%x + cvalue * rhs%x
                        self%y = self%y + cvalue * rhs%y
                        self%z = self%z + cvalue * rhs%z
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%s_v = self%s_v + cvalue * rhs%s_v
                        !
                    else
                        call errStop( "multAdd_rVector3D_SG > Unknown rVector3D_SG_t rhs store_state!" )
                    endif
                    !
                class default
                    call errStop( "multAdd_rVector3D_SG > rhs undefined." )
            end select
            !
        else
            call errStop( "multAdd_rVector3D_SG > Incompatible inputs." )
        endif
        !
    end subroutine multAdd_rVector3D_SG
    !
    !> No subroutine briefing
    !
    function dotProd_rVector3D_SG( self, rhs ) result( cvalue )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( inout ) :: rhs
        !
        complex( kind=prec ) :: cvalue
        !
        cvalue = C_ZERO
        !
        if( ( .NOT. self%is_allocated ) .OR. ( .NOT. rhs%is_allocated ) ) then
            call errStop( "dotProd_rVector3D_SG > Input vectors not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            call self%switchStoreState( rhs%store_state )
            !
            select type( rhs )
                !
                class is( rVector3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        cvalue = cvalue + sum( self%x * rhs%x )
                        cvalue = cvalue + sum( self%y * rhs%y )
                        cvalue = cvalue + sum( self%z * rhs%z )
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        cvalue = cvalue + sum( self%s_v * rhs%s_v )
                        !
                    else
                        call errStop( "dotProd_rVector3D_SG > Unknown rhs store_state!" )
                    endif
                    !
                class default
                    call errStop( "dotProd_rVector3D_SG: undefined rhs" )
                !
            end select
            !
        else
            call errStop( "dotProd_rVector3D_SG > Incompatible rhs" )
        endif
        !
    end function dotProd_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine divByValue_rVector3D_SG( self, cvalue )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        if( ( .NOT. self%is_allocated ) ) then
            call errStop( "divByValue_rVector3D_SG > Self not allocated." )
        endif
        !
        if( self%store_state .EQ. compound ) then
            !
            self%x = self%x / cvalue
            self%y = self%y / cvalue
            self%z = self%z / cvalue
            !
        else if( self%store_state .EQ. singleton ) then
            !
            self%s_v = self%s_v / cvalue
            !
        else
            call errStop( "divByValue_rVector3D_SG > Unknown self store_state!" )
        endif
        !
    end subroutine divByValue_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine divByField_rVector3D_SG( self, rhs )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( ( .NOT. self%is_allocated ) .OR. ( .NOT. rhs%is_allocated ) ) then
            call errStop( "divByField_rVector3D_SG > Input vectors not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            call self%switchStoreState( rhs%store_state )
            !
            select type( rhs )
                !
                class is( rVector3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%x = self%x / rhs%x
                        self%y = self%y / rhs%y
                        self%z = self%z / rhs%z
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%s_v = self%s_v / rhs%s_v
                        !
                    else
                        call errStop( "divByField_rVector3D_SG > Unknown rVector3D_SG_t rhs store_state!" )
                    endif
                    !
                class is( rScalar3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%x = self%x / rhs%v
                        self%y = self%y / rhs%v
                        self%z = self%z / rhs%v
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%s_v = self%s_v / rhs%s_v
                        !
                    else
                        call errStop( "divByField_rVector3D_SG > Unknown rScalar3D_SG_t rhs store_state!" )
                    endif
                    !
                class is( cScalar3D_SG_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%x = self%x / rhs%v
                        self%y = self%y / rhs%v
                        self%z = self%z / rhs%v
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%s_v = self%s_v / rhs%s_v
                        !
                    else
                        call errStop( "divByField_rVector3D_SG > Unknown cScalar3D_SG_t rhs store_state!" )
                    endif
                    !
                class default
                    call errStop( "divByField_rVector3D_SG: undefined rhs" )
                !
            end select
            !
        else
            call errStop( "divByField_rVector3D_SG: incompatible rhs" )
        endif
        !
    end subroutine divByField_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine interpFunc_rVector3D_SG( self, location, xyz, interp )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        real( kind=prec ), intent( in ) :: location(3)
        character, intent( in ) :: xyz
        class( Vector_t ), allocatable, intent( inout ) :: interp
        !
        real( kind=prec ), allocatable, dimension(:) :: xC, yC, zC
        integer :: ix, iy, iz, i
        real( kind=prec ) :: wx, wy, wz
        logical, dimension(:), allocatable :: tmp
        !
        if( ( .NOT. self%is_allocated ) ) then
            call errStop( "interpFunc_rVector3D_SG > Self not allocated." )
        endif
        !
        select case( self%grid_type )
            !
            case( EDGE )
                !
                allocate( interp, source = rVector3D_SG_t( self%grid, EDGE ) )
                !
                select case( xyz )
                    !
                    case("x")
                        !
                        allocate(xC(size(self%grid%del_x)))
                        allocate(yC(size(self%grid%dy + 1)))
                        allocate(zC(size(self%grid%dz + 1)))
                        !
                        xC = CumSum(self%grid%del_x)
                        yC = CumSum([0._prec, self%grid%dy])
                        zC = CumSum([0._prec, self%grid%dz])
                        !
                    case("y")
                        !
                        allocate(xC(size(self%grid%dx + 1)))
                        allocate(yC(size(self%grid%del_y)))
                        allocate(zC(size(self%grid%dz)))
                        
                        xC = CumSum([0._prec, self%grid%dx])
                        yC = CumSum([self%grid%del_y])
                        zC = CumSum([0._prec, self%grid%dz])
                        !
                    case("z")
                        !
                        allocate(xC(size(self%grid%dx + 1)))
                        allocate(yC(size(self%grid%dy + 1)))
                        allocate(zC(size(self%grid%del_z)))
                        !
                        xC = CumSum([0._prec, self%grid%dx])
                        yC = CumSum([0._prec, self%grid%dy])
                        zC = CumSum([self%grid%del_z])
                        !
                    case default
                        call errStop( "interpFunc_rVector3D_SG: Unknown xyz" )
                    !
                end select
                !
            case( FACE )
                !
                allocate( interp, source = rVector3D_SG_t( self%grid, FACE ) )
                !
                select case( xyz )
                    !
                    case( "x" )
                        !
                        allocate(xC(size(self%grid%dx + 1)))
                        allocate(yC(size(self%grid%del_y)))
                        allocate(zC(size(self%grid%del_z)))
                        !
                        xC = CumSum([0._prec, self%grid%dx])
                        yC = CumSum([self%grid%del_y])
                        zC = CumSum([self%grid%del_z])
                        !
                    case( "y" )
                        !
                        allocate(xC(size(self%grid%del_x)))
                        allocate(yC(size(self%grid%dy + 1)))
                        allocate(zC(size(self%grid%del_z)))
                        !
                        xC = CumSum([self%grid%del_x])
                        yC = CumSum([0._prec, self%grid%dy])
                        zC = CumSum([self%grid%del_z])
                        !
                    case( "z" )
                        !
                        allocate(xC(size(self%grid%del_x)))
                        allocate(yC(size(self%grid%del_y)))
                        allocate(zC(size(self%grid%dz + 1)))
                        !
                        xC = CumSum([self%grid%del_x])
                        yC = CumSum([self%grid%del_y])
                        zC = CumSum([0._prec, self%grid%dz])
                        !
                    case default
                        call errStop( "interpFunc_rVector3D_SG: Unknown xyz" )
                    !
                end select
                !
            case default
                call errStop( "interpFunc_rVector3D_SG: Unknown grid_type" )
            !
        end select
        !
        xC = xC + self%grid%ox
        yC = yC + self%grid%oy
        zC = zC - sum(self%grid%dz(1:self%grid%nzAir)) - self%grid%oz
        !
        tmp = location(1) > xC
        !
        ix = size( tmp )
        !
        do i = size( tmp ), 1, -1 
            if(tmp(i)) then
                ix = i
                exit
            endif
        enddo
        !
        tmp = location(2) > yC
        !
        iy = size( tmp )
        !
        do i = size( tmp ), 1, -1 
            if(tmp(i)) then
                iy = i
                exit
            endif
        enddo
        !
        tmp = location(3) > zC
        !
        iz = size( tmp )
        !
        do i = size( tmp ), 1, -1 
            if(tmp(i)) then
                iz = i
                exit
            endif
        enddo
        !
        deallocate( tmp )
        !
        !> ????
        !ix = findloc( location(1) > xC, .TRUE., back = .TRUE., dim = 1 )
        !iy = findloc( location(2) > yC, .TRUE., back = .TRUE., dim = 1 )
        !iz = findloc( location(3) > zC, .TRUE., back = .TRUE., dim = 1 )
        !
        ! Find weights
        wx = (xC(ix + 1) - location(1))/(xC(ix + 1) - xC(ix))
        !
        deallocate( xC )
        !
        wy = (yC(iy + 1) - location(2))/(yC(iy + 1) - yC(iy))
        !
        deallocate( yC )
        !
        wz = (zC(iz + 1) - location(3))/(zC(iz + 1) - zC(iz))
        !
        deallocate( zC )
        !
        select type( interp )
            !
            class is( rVector3D_SG_t )
                !
                select case( xyz )
                    !
                    case("x")
                        !
                        interp%x(ix,iy,iz) = wx*wy*wz
                        interp%x(ix+1,iy,iz) = (1-wx)*wy*wz
                        interp%x(ix,iy+1,iz) = wx*(1-wy)*wz
                        interp%x(ix,iy,iz+1) = wx*wy*(1-wz)
                        interp%x(ix,iy+1,iz+1) = wx*(1-wy)*(1-wz)
                        interp%x(ix+1,iy,iz+1) = (1-wx)*wy*(1-wz)
                        interp%x(ix+1,iy+1,iz) = (1-wx)*(1-wy)*wz
                        interp%x(ix+1,iy+1,iz+1) = (1-wx)*(1-wy)*(1-wz)
                        !
                    case("y")
                        !
                        interp%y(ix,iy,iz) = wx*wy*wz
                        interp%y(ix+1,iy,iz) = (1-wx)*wy*wz
                        interp%y(ix,iy+1,iz) = wx*(1-wy)*wz
                        interp%y(ix,iy,iz+1) = wx*wy*(1-wz)
                        interp%y(ix,iy+1,iz+1) = wx*(1-wy)*(1-wz)
                        interp%y(ix+1,iy,iz+1) = (1-wx)*wy*(1-wz)
                        interp%y(ix+1,iy+1,iz) = (1-wx)*(1-wy)*wz
                        interp%y(ix+1,iy+1,iz+1) = (1-wx)*(1-wy)*(1-wz)
                        !
                    case("z")
                        !
                        interp%z(ix,iy,iz) = wx*wy*wz
                        interp%z(ix+1,iy,iz) = (1-wx)*wy*wz
                        interp%z(ix,iy+1,iz) = wx*(1-wy)*wz
                        interp%z(ix,iy,iz+1) = wx*wy*(1-wz)
                        interp%z(ix,iy+1,iz+1) = wx*(1-wy)*(1-wz)
                        interp%z(ix+1,iy,iz+1) = (1-wx)*wy*(1-wz)
                        interp%z(ix+1,iy+1,iz) = (1-wx)*(1-wy)*wz
                        interp%z(ix+1,iy+1,iz+1) = (1-wx)*(1-wy)*(1-wz)
                        !
                    case default
                        call errStop( "interpFunc_rVector3D_SG: Unknown xyz" )
                    !
                end select !XYZ
            !
            class default
                call errStop( "interpFunc_rVector3D_SG: undefined interp" )
        !
        end select !XYZ
            !
    end subroutine interpFunc_rVector3D_SG
    !
    !> No function briefing
    !
    function getAxis_rVector3D_SG( self, comp_lbl ) result( comp )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        character, intent( in ) :: comp_lbl
        !
        complex( kind=prec ), allocatable :: comp(:, :, :)
        !
        if( ( .NOT. self%is_allocated ) ) then
            call errStop( "interpFunc_rVector3D_SG > Self not allocated." )
        endif
        !
        call self%switchStoreState( compound )
        !
        if( comp_lbl == "x" .OR. comp_lbl == "X" ) then
            comp = self%x
        elseif( comp_lbl == "y" .OR. comp_lbl == "Y" ) then
            comp = self%y
        elseif( comp_lbl == "z" .OR. comp_lbl == "Z" ) then
            comp = self%z
        else
            call errStop( "getAxis_rVector3D_SG > wrong component label" )
        endif
        !
    end function getAxis_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine getReal_rVector3D_SG( self, r_vector )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Vector_t ), allocatable, intent( out ) :: r_vector
        !
        allocate( r_vector, source = rVector3D_SG_t( self%grid, self%grid_type ) )
        !
        call r_vector%copyFrom( self )
        !
    end subroutine getReal_rVector3D_SG
    !
    !> No interface function briefing
    !
    function getX_rVector3D_SG( self ) result( x )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        !
        complex( kind=prec ), allocatable :: x(:, :, :)
        !
        call self%switchStoreState( compound )
        !
        x = self%x
        !
    end function getX_rVector3D_SG
    !
    !> No interface subroutine briefing
    !
    subroutine setX_rVector3D_SG( self, x )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), allocatable, intent( in ) :: x(:, :, :)
        !
        self%store_state = compound
        !
        if( allocated( self%s_v ) ) deallocate( self%s_v )
        !
        self%x = x
        !
    end subroutine setX_rVector3D_SG
    !
    !> No interface function briefing
    !
    function getY_rVector3D_SG( self ) result( y )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        !
        complex( kind=prec ), allocatable :: y(:, :, :)
        !
        call self%switchStoreState( compound )
        !
        y = self%y
        !
    end function getY_rVector3D_SG
    !
    !> No interface subroutine briefing
    !
    subroutine setY_rVector3D_SG( self, y )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), allocatable, intent( in ) :: y(:, :, :)
        !
        self%store_state = compound
        !
        if( allocated( self%s_v ) ) deallocate( self%s_v )
        !
        self%y = y
        !
    end subroutine setY_rVector3D_SG
    !
    !> No interface function briefing
    !
    function getZ_rVector3D_SG( self ) result( z )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        !
        complex( kind=prec ), allocatable :: z(:, :, :)
        !
        call self%switchStoreState( compound )
        !
        z = self%z
        !
    end function getZ_rVector3D_SG
    !
    !> No interface subroutine briefing
    !
    subroutine setZ_rVector3D_SG( self, z )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), allocatable, intent( in ) :: z(:, :, :)
        !
        self%store_state = compound
        !
        if( allocated( self%s_v ) ) deallocate( self%s_v )
        !
        self%z = z
        !
    end subroutine setZ_rVector3D_SG
    !
    !> No function briefing
    !
    function getSV_rVector3D_SG( self ) result( s_v )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        !
        complex( kind=prec ), allocatable :: s_v(:)
        !
        call self%switchStoreState( singleton )
        !
        s_v = self%s_v
        !
    end function getSV_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setSV_rVector3D_SG( self, s_v )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), allocatable, intent( in ) :: s_v(:)
        !
        self%store_state = singleton
        !
        if( allocated( self%x ) ) deallocate( self%x )
        if( allocated( self%y ) ) deallocate( self%y )
        if( allocated( self%z ) ) deallocate( self%z )
        !
        self%s_v = s_v
        !
    end subroutine setSV_rVector3D_SG
    !
    !> No function briefing
    !
    function getArray_rVector3D_SG( self ) result( array )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        !
        complex( kind=prec ), allocatable, dimension(:) :: array
        !
        if( ( .NOT. self%is_allocated ) ) then
            call errStop( "getArray_rVector3D_SG > Self not allocated." )
        endif
        !
        if( self%store_state .EQ. compound ) then
            !
            allocate( array( self%length() ) )
            !
            array = (/reshape(self%x, (/self%Nxyz(1), 1/)), &
                      reshape(self%y, (/self%Nxyz(2), 1/)), &
                      reshape(self%z, (/self%Nxyz(3), 1/))/)
            !
        else if( self%store_state .EQ. singleton ) then
            !
            array = self%s_v
            !
        else
            call errStop( "getArray_rVector3D_SG > Unknown store_state!" )
        endif
        !
    end function getArray_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setArray_rVector3D_SG( self, array )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), dimension(:), intent( in ) :: array
        !
        integer :: i1, i2
        !
        if( ( .NOT. self%is_allocated ) ) then
            call errStop( "setArray_rVector3D_SG > Self not allocated." )
        endif
        !
        if( self%store_state .EQ. compound ) then
            !
            !> Ex
            i1 = 1; i2 = self%Nxyz(1)
            !
            self%x = reshape( array(i1:i2), self%NdX )
            !> Ey
            i1 = i2 + 1; i2 = i2 + self%Nxyz(2)
            !
            self%y = reshape( array(i1:i2), self%NdY )
            !> Ez
            i1 = i2 + 1; i2 = i2 + self%Nxyz(3)
            !
            self%z = reshape(array(i1:i2), self%NdZ)
            !
        else if( self%store_state .EQ. singleton ) then
            !
            self%s_v = array
            !
        else
            call errStop( "setArray_rVector3D_SG > Unknown store_state!" )
        endif
        !
    end subroutine setArray_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine copyFrom_rVector3D_SG( self, rhs )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( .NOT. rhs%is_allocated) then
            call errStop( "copyFrom_rVector3D_SG > rhs not allocated" )
        endif
        !
        self%grid => rhs%grid
        self%grid_type = rhs%grid_type
        self%nx = rhs%nx
        self%ny = rhs%ny
        self%nz = rhs%nz
        self%store_state = rhs%store_state
        !
        if( allocated( rhs%ind_interior ) ) then
            self%ind_interior = rhs%ind_interior
        else
            call errStop( "copyFrom_rVector3D_SG > rhs self%ind_interior not allocated" )
        endif
        !
        if( allocated( rhs%ind_boundaries ) ) then
            self%ind_boundaries = rhs%ind_boundaries
        else
            call errStop( "copyFrom_rVector3D_SG > rhs self%ind_interior not allocated" )
        endif
        !
        select type( rhs )
            !
            class is( rVector3D_SG_t )
                !
                self%NdX = rhs%NdX
                self%NdY = rhs%NdY
                self%NdZ = rhs%NdZ
                self%Nxyz = rhs%Nxyz
                !
                if( rhs%store_state .EQ. compound ) then
                    !
                    if( allocated( self%x ) ) deallocate( self%x )
                    allocate( self%x, source = rhs%x )
                    !
                    if( allocated( self%y ) ) deallocate( self%y )
                    allocate( self%y, source = rhs%y )
                    !
                    if( allocated( self%z ) ) deallocate( self%z )
                    allocate( self%z, source = rhs%z )
                    !
                else if( rhs%store_state .EQ. singleton ) then
                    !
                    if( allocated( self%s_v ) ) deallocate( self%s_v )
                    allocate( self%s_v, source = rhs%s_v )
                    !
                else
                    call errStop( "copyFrom_rVector3D_SG > Unknown store_state!" )
                endif
                !
            class default
                call errStop( "copyFrom_rVector3D_SG > Different type of rhs" )
            !
        end select
        !
        self%is_allocated = .TRUE.
        !
    end subroutine copyFrom_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine read_rVector3D_SG( self, funit, ftype )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        integer, intent( in ) :: funit
        character(:), allocatable, intent( in ), optional :: ftype
        !
        integer :: Nx, Ny, Nz
        character(4) :: grid_type
        logical :: ok, hasname, binary
        character(80) :: fname, isbinary
        !
        call self%switchStoreState( compound )
        !
        binary = .TRUE.
        !
        inquire( funit, opened = ok, named = hasname, name = fname, unformatted = isbinary )
        !
        if( ok ) then
            !
            !> Check that the file is unformatted if binary, formatted if ascii.
            if((index(isbinary, "yes") > 0 .OR.index(isbinary, "YES") > 0) &
                  .AND.   .NOT. binary ) then
                write( *, * ) "Error: read_rVector3D_SG > Unable to read_rVector3D_SG vector from unformatted file. ", &
                        trim(fname), "."
                stop
            else if((index(isbinary, "no") > 0 .OR.index(isbinary, "NO") > 0) &
                  .AND.binary) then
                write( *, * ) "Error: read_rVector3D_SG > Unable to read_rVector3D_SG vector from formatted file ", &
                        trim(fname), "."
                stop
            endif
            !
            read(funit) Nx, Ny, Nz, grid_type
            !
            if(  .NOT. self%is_allocated) then
                write( *, * ) "Error: read_rVector3D_SG > Vector must be allocated before read_rVector3D_SGing from ", &
                        trim(fname), "."
                stop
            else if(self%grid_type.NE.grid_type) then
                write( *, * ) "Error: read_rVector3D_SG > Vector must be of type ", grid_type, &
                        &            "           before read_rVector3D_SGing from ", trim (fname), "."
                stop
            else if((self%nx.NE.Nx).OR. &
                  (self%ny.NE.Ny).OR.(self%nz.NE.Nz)) then
                write( *, * ) "Error: read_rVector3D_SG > Wrong size of vector on input from ", trim (fname), "."
                stop
            endif
            !
            read(funit) self%x
            read(funit) self%y
            read(funit) self%z
            !
        else
            stop "Error: read_rVector3D_SG: unable to open file"
        endif
        !
    end subroutine read_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine write_rVector3D_SG( self, funit, ftype )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        integer, intent( in ) :: funit
        character(:), allocatable, intent( in ), optional :: ftype
        !
        logical :: ok, hasname, binary
        character(80) :: fname, isbinary
        !
        if( .NOT. self%is_allocated) then
            stop "Error: write_rVector3D_SG > Not allocated."
        endif
        !
        call self%switchStoreState( compound )
        !
        binary = .TRUE.
        !
        inquire( funit, opened = ok, named = hasname, name = fname, unformatted = isbinary )
        !
        if( ok ) then
            !
            !> Check that the file is unformatted if binary, formatted if ascii.
            if((index(isbinary, "yes") > 0.OR.index(isbinary, "YES") > 0) &
                  .AND.   .NOT. binary) then
                write( *, * ) "Error: write_rVector3D_SG > Unable to write_rVector3D_SG vector to unformatted file. ", &
                        trim(fname), "."
                stop
            else if((index(isbinary,"no") > 0.OR.index(isbinary,"NO") > 0) &
                  .AND.binary) then
                write( *, * ) "Error: write_rVector3D_SG > Unable to write_rVector3D_SG vector to formatted file. ", &
                        trim(fname), "."
                stop
            endif
            !
            write(funit) self%nx, self%ny, self%nz, self%grid_type
            write(funit) self%x
            write(funit) self%y
            write(funit) self%z
            !
        else
            stop "Error: write_rVector3D_SG > unable to open file"
        endif
        !
    end subroutine write_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine print_rVector3D_SG( self, io_unit, title, append )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        integer, intent( in ), optional :: io_unit
        character(*), intent( in ), optional :: title
        logical, intent( in ), optional :: append
        !
        integer :: ix, iy, iz, funit
        !
        call self%switchStoreState( compound )
        !
        if( present( io_unit ) ) then
            funit = io_unit
        else
            funit = 0
        endif
        !
        if( present( title ) ) write( funit, * ) title
        !
        write( funit, * ) self%nx, self%ny, self%nz
        write(funit, * ) "x-component",self%NdX
        do ix = 1, self%NdX(1)
             do iy = 1, self%NdX(2)
                do iz = 1, self%NdX(3)
                     if( self%x( ix, iy, iz ) /= 0 ) then
                        write(funit,*) ix,iy,iz, ":[", self%x( ix, iy, iz ), "]"
                     endif
                enddo
             enddo
        enddo
        !
        write(funit,*) "y-component",self%NdY
        do ix = 1, self%NdY(1)
             do iy = 1, self%NdY(2)
                do iz = 1, self%NdY(3)
                     if( self%y( ix, iy, iz ) /= 0 ) then
                        write(funit,*) ix,iy,iz, ":[", self%y( ix, iy, iz ), "]"
                     endif
                enddo
             enddo
        enddo
        !
        write(funit,*) "z-component",self%NdZ
        do ix = 1, self%NdZ(1)
             do iy = 1, self%NdZ(2)
                do iz = 1, self%NdZ(3)
                     if( self%z( ix, iy, iz ) /= 0 ) then
                        write(funit,*) ix,iy,iz, ":[", self%z( ix, iy, iz ), "]"
                     endif
                enddo
             enddo
        enddo
        !
    end subroutine print_rVector3D_SG
    !
    !> Convert an input Rvector E to a 1-D real array,
    !> following standard staggered grid ordering
    !
    subroutine getRvector( E, v )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: E
        real( kind=prec ), dimension(:), pointer, intent( inout ) :: v
        !
        integer :: nVec(3), nVecT, id(1), i1, i2
        !
        if( .NOT. E%is_allocated ) then
            stop "Error: getRvector > E not allocated"
        endif
        !
        nVec(1) = size(E%x)
        nVec(2) = size(E%y)
        nVec(3) = size(E%z)
        nVecT = nVec(1)+nVec(2)+nVec(3)
        !
        if(associated(v)) then
            if(nVect.ne.size(v)) then
                deallocate(v)
            endif
        endif
        if(.not.associated(v)) then
            allocate(v(nVecT))
        endif
        !   now that we know v is allocated, an of proper size
        !     just copy contents of E into v
        id(1) = nVec(1)
        i1 = 1
        i2 = nVec(1)
        v(i1:i2) = reshape(E%x,id)
        id(1) = nVec(2)
        i1 = i2+1
        i2 = i2+nVec(2)
        v(i1:i2) = reshape(E%y,id)
        id(1) = nVec(3)
        i1 = i2+1
        i2 = i2+nVec(3)
        v(i1:i2) = reshape(E%z,id)
        !
    end subroutine getRvector
    !
    !> Copy contents of v into an already created and allocated Rvector
    !
    subroutine setRvector( v, E )
        implicit none
        !
        real( kind=prec ), dimension(:), intent( in ) :: v
        class( rVector3D_SG_t ), intent( inout ) :: E
        !
        integer :: nVec(3,3), nVecT, id(3), i, i1, i2
        !
        if( .NOT. E%is_allocated ) then
            stop "Error: setRvector > E not allocated"
        endif
        !
        do i = 1,3 
            nVec(1,i) = size(E%x,i)
            nVec(2,i) = size(E%y,i)
            nVec(3,i) = size(E%z,i)
        enddo
        !
        nVect = 0
        !
        do i =1,3
            nVecT = nVecT + nVec(i,1)*nVec(i,2)*nVec(i,3)
        enddo
        !
        if( nVecT .NE. size(v) ) then
            stop "Error: setRvector > Input vector of incorrect size"
        endif
        !     copy contents of v into E
        i1 = 1
        i2 = nVec(1,1)*nVec(1,2)*nVec(1,3)
        id = nVec(1,:)
        E%x = reshape(v(i1:i2),id)
        i1 = i2+1
        i2 = i2+nVec(2,1)*nVec(2,2)*nVec(2,3)
        id = nVec(2,:)
        E%y = reshape(v(i1:i2),id)
        i1 = i2+1
        i2 = i2+nVec(3,1)*nVec(3,2)*nVec(3,3)
        id = nVec(3,:)
        E%z = reshape(v(i1:i2),id)
        !
    end subroutine setRvector
    !
    subroutine edgeLength( grid, l_e )
        implicit none
        !
        class( Grid_t ), intent( in ) :: grid
        type( rVector3D_SG_t ), intent( inout )  :: l_e
        !
        integer :: ix, iy, iz
        !
        l_e = rVector3D_SG_t( grid, EDGE )
        !
        ! x-component edge length elements
        do ix = 1,grid%nx
            l_e%x(ix, :, :) = grid%dx(ix)
        enddo
        !
        ! y-component edge length elements
        do iy = 1,grid%ny
            l_e%y(:, iy, :) = grid%dy(iy)
        enddo
        !
        ! z-component edge length elements
        do iz = 1,grid%nz
            l_e%z(:, :, iz) = grid%dz(iz)
        enddo
        !
    end subroutine edgeLength
    !
end module rVector3D_SG
