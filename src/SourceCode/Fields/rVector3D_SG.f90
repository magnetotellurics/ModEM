!
!> Derived class to define a rVector3D_SG
!
module rVector3D_SG
    !
    use MatUtils
    use Vector
    use Grid3D_SG
    use rScalar3D_SG
    !
    type, extends( Vector_t ) :: rVector3D_SG_t
        !
        integer, dimension(3) :: NdX, NdY, NdZ, Nxyz
        !
        real( kind=prec ), allocatable, dimension(:,:,:) :: x, y, z
        !
        contains
            !
            final :: rVector3D_SG_dtor
            !
            !> Boundary operations
            procedure, public :: setAllBoundary => setAllBoundary_rVector3D_SG
            procedure, public :: setOneBoundary => setOneBoundary_rVector3D_SG
            !
            !> Dimensioning operations
            procedure, public :: length => length_rVector3D_SG
            procedure, public :: setVecComponents => setVecComponents_rVector3D_SG
            !
            !> Arithmetic/algebraic unary operations
            procedure, public :: zeros => zeros_rVector3D_SG
            !
            procedure, public :: sumEdge => sumEdge_rVector3D_SG
            procedure, public :: sumEdgeVTI => sumEdgeVTI_rVector3D_SG
            !
            procedure, public :: sumCell => sumCell_rVector3D_SG
            procedure, public :: sumCellVTI => sumCellVTI_rVector3D_SG
            !
            procedure, public :: conjugate => conjugate_rVector3D_SG
            !
            !> Arithmetic/algebraic binary operations
            procedure, public :: add => add_rVector3D_SG
            !
            procedure, public :: subValue => subValue_rVector3D_SG
            procedure, public :: subField => subField_rVector3D_SG
            !
            procedure, public :: linComb => linComb_rVector3D_SG
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
            procedure, public :: getArray => getArray_rVector3D_SG
            procedure, public :: setArray => setArray_rVector3D_SG
            !
            procedure, public :: copyFrom => copyFrom_rVector3D_SG
            !
            procedure, public :: getReal => getReal_rVector3D_SG
            !
            procedure, public :: edgeLength => edgeLength_rVector3D_SG
            !
            !> I/O operations
            procedure, public :: read => read_rVector3D_SG
            procedure, public :: write => write_rVector3D_SG
            procedure, public :: print => print_rVector3D_SG
            !
    end type rVector3D_SG_t
    !
    public :: getRvector, setRvector, Edgelength
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
        integer :: alloc_stat
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
        !
        if( self%grid_type == EDGE ) then
            !
            allocate( self%x(self%nx, self%ny + 1, self%nz + 1), stat=alloc_stat )
            self%is_allocated = alloc_stat .EQ. 0
            !
            allocate( self%y(self%nx + 1, self%ny, self%nz + 1), stat=alloc_stat )
            self%is_allocated = self%is_allocated .AND. ( alloc_stat .EQ. 0 )
            !
            allocate( self%z(self%nx + 1, self%ny + 1, self%nz), stat=alloc_stat )
            self%is_allocated = self%is_allocated .AND. ( alloc_stat .EQ. 0 )
            !
            self%NdX = (/self%nx, self%ny + 1, self%nz + 1/)
            self%NdY = (/self%nx + 1, self%ny, self%nz + 1/)
            self%NdZ = (/self%nx + 1, self%ny + 1, self%nz/)
            !
        elseif( self%grid_type == FACE ) then
            !
            allocate( self%x(self%nx + 1, self%ny, self%nz), stat=alloc_stat )
            self%is_allocated = alloc_stat .EQ. 0
            !
            allocate( self%y(self%nx, self%ny + 1, self%nz), stat=alloc_stat )
            self%is_allocated = self%is_allocated .AND. ( alloc_stat .EQ. 0 )
            !
            allocate( self%z(self%nx, self%ny, self%nz + 1), stat=alloc_stat)
            self%is_allocated = self%is_allocated .AND. ( alloc_stat .EQ. 0 )
            !
            self%NdX = (/self%nx + 1, self%ny, self%nz/)
            self%NdY = (/self%nx, self%ny + 1, self%nz/)
            self%NdZ = (/self%nx, self%ny, self%nz + 1/)
            !
        else
            call errStop( "rVector3D_SG_ctor > Only EDGE or FACE types allowed." )
        endif
        !
        if( self%is_allocated ) then
            !
            self%x = C_ZERO
            self%y = C_ZERO
            self%z = C_ZERO
            !
            self%Nxyz = (/product(self%NdX), product(self%NdY), product(self%NdZ)/)
            !
        else
            call errStop( "rVector3D_SG_ctor > Unable to allocate vector." )
        endif
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
        if( .NOT. present( int_only ) ) then
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
    function length_rVector3D_SG( self ) result( field_length )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        !
        integer :: field_length
        !
        field_length = self%Nxyz(1) + self%Nxyz(2) + self%Nxyz(3)
        !
    end function length_rVector3D_SG
    !
    !> No subroutine briefing
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
                self%z( x1:x2:xstep, y1:y2:ystep, z1:z2:zstep ) = cvalue
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
        self%x = R_ZERO
        self%y = R_ZERO
        self%z = R_ZERO
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
                                    self%x(:,2:x_yend,1:x_zend-1)   + &
                                    self%x(:,1:x_yend-1,2:x_zend)   + &
                                    self%x(:,2:x_yend,2:x_zend)     + &
                                    self%y(1:y_xend-1,:,1:y_zend-1) + &
                                    self%y(2:y_xend,:,1:y_zend-1)   + &
                                    self%y(1:y_xend-1,:,2:y_zend)   + &
                                    self%y(2:y_xend,:,2:y_zend)     + &
                                    self%z(1:z_xend-1,1:z_yend-1,:) + &
                                    self%z(2:z_xend,1:z_yend-1,:)   + &
                                    self%z(1:z_xend-1,2:z_yend,:)   + &
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
                call errStop( "sumEdge_rVector3D_SG > Unclassified cell_out" )
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
        integer :: x_xend, x_yend, x_zend
        integer :: y_xend, y_yend, y_zend
        integer :: z_xend, z_yend, z_zend
        logical :: is_interior_only
        !
        if( .NOT. self%is_allocated ) then
             call errStop( "sumEdgeVTI_rVector3D_SG > self not allocated." )
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
        allocate( cell_h_out, source = rScalar3D_SG_t( self%grid, CELL ) )
        !
        allocate( cell_v_out, source = rScalar3D_SG_t( self%grid, CELL ) )
        !
        select type( cell_h_out )
            !
            class is( rScalar3D_SG_t )
                !
                select type( cell_v_out )
                    !
                    class is( rScalar3D_SG_t )
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
                                cell_h_out%v = self%x(:,1:x_yend-1,1:x_zend-1) + &
                                            self%x(:,2:x_yend,1:x_zend-1)   + &
                                            self%x(:,1:x_yend-1,2:x_zend)   + &
                                            self%x(:,2:x_yend,2:x_zend)     + &
                                            self%y(1:y_xend-1,:,1:y_zend-1) + &
                                            self%y(2:y_xend,:,1:y_zend-1)   + &
                                            self%y(1:y_xend-1,:,2:y_zend)   + &
                                            self%y(2:y_xend,:,2:y_zend)
                                !
                                cell_v_out%v = self%z(1:z_xend-1,1:z_yend-1,:) + &
                                            self%z(2:z_xend,1:z_yend-1,:)   + &
                                            self%z(1:z_xend-1,2:z_yend,:)   + &
                                            self%z(2:z_xend,2:z_yend,:)
                                !
                            case( FACE )
                                !
                                x_xend = size( self%x, 1 )
                                y_xend = size( self%y, 1 )
                                z_xend = size( self%z, 1 )
                                !
                                cell_h_out%v = self%x(1:x_xend-1,:,:) + self%x(2:x_xend,:,:) + &
                                                    self%y(:,1:y_yend-1,:) + self%y(:,2:y_yend,:)
                                !
                                cell_v_out%v = self%z(:,:,1:z_zend-1) + self%z(:,:,2:z_zend)
                                !
                            case default
                                call errStop( "sumEdgeVTI_rVector3D_SG: undefined self%grid_type" )
                            !
                        end select
                        !
                    class default
                        call errStop( "sumEdgeVTI_rVector3D_SG > Unclassified cell_v_out" )
                    !
                end select
                !
            class default
                call errStop( "sumEdgeVTI_rVector3D_SG > Unclassified cell_h_out" )
            !
        end select
        !
    end subroutine sumEdgeVTI_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine sumCell_rVector3D_SG( self, cell_in, ptype )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Scalar_t ), intent( in ) :: cell_in
        character(*), intent( in ), optional :: ptype
        !
        character( len=4 ) :: grid_type
        integer :: xend, yend, zend
        integer :: v_xend, v_yend, v_zend
        integer :: ix, iy, iz
        !
        if( .NOT. self%is_allocated ) then
             call errStop( "sumCell_rVector3D_SG > self not allocated." )
        endif
        !
        if( .NOT. cell_in%is_allocated ) then
             call errStop( "sumCell_rVector3D_SG > cell_in not allocated." )
        endif
        !
        if( .NOT. cell_in%grid_type == CELL ) then
            call errStop( "sumCell_rVector3D_SG > Only CELL type supported." )
        endif
        !
        if( .NOT. present( ptype ) ) then
            grid_type = EDGE
        else
            grid_type = ptype
        endif
        !
        select type( cell_in )
            !
            class is( rScalar3D_SG_t )
                !
                v_xend = size( cell_in%v, 1 )
                v_yend = size( cell_in%v, 2 )
                v_zend = size( cell_in%v, 3 )
                !
                select case( grid_type )
                    !
                    case( EDGE )
                        !
                        !> for x-components inside the domain
                        do ix = 1, self%grid%nx
                            do iy = 2, self%grid%ny
                                do iz = 2, self%grid%nz
                                    self%x(ix, iy, iz) = ( cell_in%v(ix, iy-1, iz-1) + cell_in%v(ix, iy, iz-1) + &
                                    cell_in%v(ix, iy-1, iz) + cell_in%v(ix, iy, iz) )
                                enddo
                                !
                                self%x(ix, iy, 1) =  cell_in%v(ix, iy-1, 1) + cell_in%v(ix, iy, 1 )
                                self%x(ix, iy, self%grid%nz+1) =  cell_in%v(ix, iy-1, self%grid%nz) + cell_in%v(ix, iy, self%grid%nz )
                                !
                            enddo
                        enddo
                        !
                        !> for y-components inside the domain
                        do ix = 2, self%grid%nx
                            do iy = 1, self%grid%ny
                                do iz = 2, self%grid%nz
                                    self%y(ix, iy, iz) = ( cell_in%v(ix-1, iy, iz-1) + cell_in%v(ix, iy, iz-1) + &
                                    cell_in%v(ix-1, iy, iz) + cell_in%v(ix, iy, iz) )
                                enddo
                                !
                                self%y(ix, iy, 1) =  cell_in%v(ix-1, iy, 1) + cell_in%v(ix, iy, 1)
                                self%y(ix, iy, self%grid%nz+1) =  cell_in%v(ix-1, iy, self%grid%nz) + cell_in%v(ix, iy, self%grid%nz )
                                !
                            enddo
                        enddo
                        !
                        !> for z-components inside the domain
                        do ix = 2, self%grid%nx
                            do iy = 2, self%grid%ny
                                do iz = 1, self%grid%nz
                                    self%z(ix, iy, iz) = ( cell_in%v(ix-1, iy-1, iz) + cell_in%v(ix-1, iy, iz) + &
                                    cell_in%v(ix, iy-1, iz) + cell_in%v(ix, iy, iz) )
                                enddo
                            enddo
                        enddo
                        !
                    case( FACE )
                        !
                        xend = size(self%x, 1)
                        self%x(2:xend-1,:,:) = cell_in%v(1:v_xend-1,:,:) + cell_in%v(2:v_xend,:,:)
                        !
                        yend = size(self%y, 1)
                        self%y(:, 2:yend-1, :) = cell_in%v(:, 1:v_yend-1, :) + cell_in%v(:, 2:v_yend, :)
                        !
                        zend = size(self%z, 1) 
                        self%z(:, :, 2:zend-1) = cell_in%v(:, :, 1:v_zend-1) + cell_in%v(:, :, 2:v_zend)
                        !
                    case default
                        call errStop( "sumCell_rVector3D_SG: Unknown type" )
                    !
                end select !type
                !
            class default
                call errStop( "sumCell_rVector3D_SG > Unclassified cell_in" )
            !
        end select
        !
    end subroutine sumCell_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine sumCellVTI_rVector3D_SG( self, cell_h_in, cell_v_in, ptype )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Scalar_t ), intent( in ) :: cell_h_in, cell_v_in
        character(*), intent( in ), optional :: ptype
        !
        character( len=4 ) :: grid_type
        integer :: xend, yend, zend
        integer :: v_xend, v_yend, v_zend
        integer :: ix, iy, iz
        !
        if( .NOT. self%is_allocated ) then
             call errStop( "sumCellVTI_rVector3D_SG > self not allocated." )
        endif
        !
        if( .NOT. cell_h_in%is_allocated ) then
             call errStop( "sumCellVTI_rVector3D_SG > cell_h_in not allocated." )
        endif
        !
        if( .NOT. cell_v_in%is_allocated ) then
             call errStop( "sumCellVTI_rVector3D_SG > cell_v_in not allocated." )
        endif
        !
        if( .NOT. self%grid_type == CELL ) then
            call errStop( "sumCellVTI_rVector3D_SG > Only CELL type supported." )
        endif
        !
        if( .NOT. present( ptype ) ) then
            grid_type = EDGE
        else
            grid_type = ptype
        endif
        !
        select type( cell_h_in )
            !
            class is( rScalar3D_SG_t )
                !
                select type( cell_v_in )
                    !
                    class is( rScalar3D_SG_t )
                        !
                        v_xend = size( cell_h_in%v, 1 )
                        v_yend = size( cell_h_in%v, 2 )
                        v_zend = size( cell_v_in%v, 3 )
                        !
                        select case( grid_type )
                            !
                            case( EDGE )
                                !
                                !> for x-components inside the domain
                                do ix = 1, self%grid%nx
                                    do iy = 2, self%grid%ny
                                        do iz = 2, self%grid%nz
                                            self%x(ix, iy, iz) = ( cell_h_in%v(ix, iy-1, iz-1) + cell_h_in%v(ix, iy, iz-1) + &
                                            cell_h_in%v(ix, iy-1, iz) + cell_h_in%v(ix, iy, iz) )! / 4.0d0
                                        enddo
                                        !
                                        self%x(ix, iy, 1) =  cell_h_in%v(ix, iy-1, 1) + cell_h_in%v(ix, iy, 1 )
                                        self%x(ix, iy, self%grid%nz+1) =  cell_h_in%v(ix, iy-1, self%grid%nz) + cell_h_in%v(ix, iy, self%grid%nz )
                                        !
                                    enddo
                                enddo
                                !
                                !> for y-components inside the domain
                                do ix = 2, self%grid%nx
                                    do iy = 1, self%grid%ny
                                        do iz = 2, self%grid%nz
                                            self%y(ix, iy, iz) = ( cell_h_in%v(ix-1, iy, iz-1) + cell_h_in%v(ix, iy, iz-1) + &
                                            cell_h_in%v(ix-1, iy, iz) + cell_h_in%v(ix, iy, iz) )! / 4.0d0
                                        enddo
                                        !
                                        self%y(ix, iy, 1) =  cell_h_in%v(ix-1, iy, 1) + cell_h_in%v(ix, iy, 1)
                                        self%y(ix, iy, self%grid%nz+1) =  cell_h_in%v(ix-1, iy, self%grid%nz) + cell_h_in%v(ix, iy, self%grid%nz )
                                        !
                                    enddo
                                enddo
                                !
                                !> for z-components inside the domain
                                do ix = 2, self%grid%nx
                                    do iy = 2, self%grid%ny
                                        do iz = 1, self%grid%nz
                                            self%z(ix, iy, iz) = ( cell_v_in%v(ix-1, iy-1, iz) + cell_v_in%v(ix-1, iy, iz) + &
                                            cell_v_in%v(ix, iy-1, iz) + cell_v_in%v(ix, iy, iz) )! / 4.0d0
                                        enddo
                                    enddo
                                enddo
                                !
                            case( FACE )
                                !
                                xend = size(self%x, 1)
                                self%x(2:xend-1,:,:) = cell_h_in%v(1:v_xend-1,:,:) + cell_h_in%v(2:v_xend,:,:)
                                !
                                yend = size(self%y, 1)
                                self%y(:, 2:yend-1, :) = cell_h_in%v(:, 1:v_yend-1, :) + cell_h_in%v(:, 2:v_yend, :)
                                !
                                zend = size(self%z, 1) 
                                self%z(:, :, 2:zend-1) = cell_v_in%v(:, :, 1:v_zend-1) + cell_v_in%v(:, :, 2:v_zend)
                                !
                            case default
                                call errStop( "sumCellVTI_rVector3D_SG: Unknown type" )
                            !
                        end select !grid_type
                        !
                    class default
                        call errStop( "sumCellVTI_rVector3D_SG > Unclassified cell_v_in" )
                    !
                end select
                !
            class default
                call errStop( "sumCellVTI_rVector3D_SG > Unclassified cell_h_in" )
            !
        end select
        !
    end subroutine sumCellVTI_rVector3D_SG
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
            select type( rhs )
                !
                class is( rVector3D_SG_t )
                    !
                    self%x = self%x + rhs%x
                    self%y = self%y + rhs%y
                    self%z = self%z + rhs%z
                    !
                class default
                    call errStop( "add_rVector3D_SG > Undefined compound rhs" )
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
    subroutine subValue_rVector3D_SG( self, cvalue )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        if( ( .NOT. self%is_allocated ) ) then
            call errStop( "subValue_rVector3D_SG > self not allocated." )
        endif
        !
        self%x = self%x - cvalue
        self%y = self%y - cvalue
        self%z = self%z - cvalue
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
        if( .NOT. self%is_allocated ) then
            call errStop( "subField_rVector3D_SG > self not allocated." )
        endif
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "subField_rVector3D_SG > rhs not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type( rhs )
                !
                class is( rVector3D_SG_t )
                    !
                    self%x = self%x - rhs%x
                    self%y = self%y - rhs%y
                    self%z = self%z - rhs%z
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
    subroutine linComb_rVector3D_SG( self, rhs, c1, c2 )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        complex( kind=prec ), intent( in ) :: c1, c2
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "linComb_rVector3D_SG > self not allocated." )
        endif
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "linComb_rVector3D_SG > rhs not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type( rhs )
                !
                class is( rVector3D_SG_t )
                    !
                    self%x = c1 * self%x + c2 * rhs%x
                    self%y = c1 * self%y + c2 * rhs%y
                    self%z = c1 * self%z + c2 * rhs%z
                    !
                class default
                    call errStop( "linComb_rVector3D_SG > Undefined rhs" )
                    !
            end select
            !
        else
            call errStop( "linComb_rVector3D_SG > Incompatible inputs." )
        endif
        !
    end subroutine linComb_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multByReal_rVector3D_SG( self, rvalue )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: rvalue
        !
        self%x = self%x * rvalue
        self%y = self%y * rvalue
        self%z = self%z * rvalue
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
        self%x = self%x * cvalue
        self%y = self%y * cvalue
        self%z = self%z * cvalue
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
        if( .NOT. self%is_allocated ) then
            call errStop( "multByField_rVector3D_SG > self not allocated." )
        endif
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "multByField_rVector3D_SG > rhs not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type( rhs )
                !
                class is( rVector3D_SG_t )
                    !
                    self%x = self%x * rhs%x
                    self%y = self%y * rhs%y
                    self%z = self%z * rhs%z
                    !
                class default
                    call errStop( "multByField_rVector3D_SG > Undefined rhs" )
                    !
            end select
            !
        else
            call errStop( "multByField_rVector3D_SG > Incompatible inputs." )
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
        class( Field_t ), intent( in ) :: rhs
        !
        class( Vector_t ), allocatable :: diag_mult
        !
        type( rVector3D_SG_t ) :: diag_mult_temp
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "diagMult_rVector3D_SG > self not allocated." )
        endif
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "diagMult_rVector3D_SG > rhs not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            diag_mult_temp = rVector3D_SG_t( self%grid, self%grid_type )
            !
            select type( rhs )
                !
                class is( rVector3D_SG_t )
                    !
                    diag_mult_temp%x = self%x * rhs%x
                    diag_mult_temp%y = self%y * rhs%y
                    diag_mult_temp%z = self%z * rhs%z
                    !
                class default
                    call errStop( "diagMult_rVector3D_SG > Undefined rhs" )
                    !
            end select
            !
            allocate( diag_mult, source = diag_mult_temp )
            !
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
        if( .NOT. self%is_allocated ) then
            call errStop( "multAdd_rVector3D_SG > self not allocated." )
        endif
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "multAdd_rVector3D_SG > rhs not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type( rhs )
                !
                class is( rVector3D_SG_t )
                    !
                    self%x = self%x + cvalue * rhs%x
                    self%y = self%y + cvalue * rhs%y
                    self%z = self%z + cvalue * rhs%z
                    !
                class default
                    call errStop( "multAdd_rVector3D_SG > Undefined rhs" )
                    !
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
        class( rVector3D_SG_t ), intent( in ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        complex( kind=prec ) :: cvalue
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "dotProd_rVector3D_SG > self not allocated." )
        endif
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "dotProd_rVector3D_SG > rhs not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type( rhs )
                !
                class is( rVector3D_SG_t )
                    !
                    cvalue = sum( self%x * rhs%x )
                    cvalue = cvalue + sum( self%y * rhs%y )
                    cvalue = cvalue + sum( self%z * rhs%z )
                    !
                class default
                    call errStop( "dotProd_rVector3D_SG > Undefined rhs" )
                    !
            end select
            !
        else
            call errStop( "dotProd_rVector3D_SG > Incompatible inputs." )
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
            call errStop( "divByValue_rVector3D_SG > self not allocated." )
        endif
        !
        self%x = self%x / cvalue
        self%y = self%y / cvalue
        self%z = self%z / cvalue
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
        if( .NOT. self%is_allocated ) then
            call errStop( "divByField_rVector3D_SG > self not allocated." )
        endif
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "divByField_rVector3D_SG > rhs not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type( rhs )
                !
                class is( rVector3D_SG_t )
                    !
                    self%x = self%x / rhs%x
                    self%y = self%y / rhs%y
                    self%z = self%z / rhs%z
                    !
                class default
                    call errStop( "divByField_rVector3D_SG > Undefined rhs" )
                    !
            end select
            !
        else
            call errStop( "divByField_rVector3D_SG > Incompatible inputs." )
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
        class( Vector_t ), intent( inout ) :: interp
        !
        integer :: ix, iy, iz, i
        real( kind=prec ) :: wx, wy, wz
        logical, allocatable, dimension(:) :: tmp
        real( kind=prec ), allocatable, dimension(:) :: xC, yC, zC
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "interpFunc_rVector3D_SG > self not allocated." )
        endif
        !
        if( .NOT. interp%is_allocated ) then
            call errStop( "interpFunc_rVector3D_SG > interp not allocated." )
        endif
        !
		write( *, * ) "interpFunc_rVector3D_SG"
		!
        select type( grid => self%grid )
            !
            class is( Grid3D_SG_t )
                !
                select case( self%grid_type )
                    !
                    case( EDGE )
                        !
                        select case( xyz )
                            !
                            case("x")
                                !
                                allocate(xC(size(grid%del_x)))
                                allocate(yC(size(grid%dy) + 1))
                                allocate(zC(size(grid%dz) + 1))
                                !
                                xC = CumSum(grid%del_x)          ! grid%xCenter
                                yC = CumSum([0._prec, grid%dy])  ! grid%yEdge
                                zC = CumSum([0._prec, grid%dz])  ! grid%zEdge
                                !
                            case("y")
                                !
                                allocate(xC(size(grid%dx) + 1))   !   etc.
                                allocate(yC(size(grid%del_y)))
                                allocate(zC(size(grid%dz)))
                                
                                xC = CumSum([0._prec, grid%dx])
                                yC = CumSum([grid%del_y])
                                zC = CumSum([0._prec, grid%dz])
                                !
                            case("z")
                                !
                                allocate(xC(size(grid%dx) + 1))
                                allocate(yC(size(grid%dy) + 1))
                                allocate(zC(size(grid%del_z)))
                                !
                                xC = CumSum([0._prec, grid%dx])
                                yC = CumSum([0._prec, grid%dy])
                                zC = CumSum([grid%del_z])
                                !
                            case default
                                call errStop( "interpFunc_rVector3D_SG: Unknown xyz" )
                            !
                        end select
                        !
                    case( FACE )
                        !
                        select case( xyz )
                            !
                            case( "x" )
                                !
                                allocate(xC(size(grid%dx) + 1))
                                allocate(yC(size(grid%del_y)))
                                allocate(zC(size(grid%del_z)))
                                !
                                xC = CumSum([0._prec, grid%dx])
                                yC = CumSum([grid%del_y])
                                zC = CumSum([grid%del_z])
                                !
                            case( "y" )
                                !
                                allocate(xC(size(grid%del_x)))
                                allocate(yC(size(grid%dy) + 1))
                                allocate(zC(size(grid%del_z)))
                                !
                                xC = CumSum([grid%del_x])
                                yC = CumSum([0._prec, grid%dy])
                                zC = CumSum([grid%del_z])
                                !
                            case( "z" )
                                !
                                allocate(xC(size(grid%del_x)))
                                allocate(yC(size(grid%del_y)))
                                allocate(zC(size(grid%dz) + 1))
                                !
                                xC = CumSum([grid%del_x])
                                yC = CumSum([grid%del_y])
                                zC = CumSum([0._prec, grid%dz])
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
                xC = xC + grid%ox
                yC = yC + grid%oy
                zC = zC - sum( grid%dz(1:grid%nzAir) ) + grid%oz
                !
            class default
                call errStop( "interpFunc_rVector3D_SG > Undefined grid" )
        end select
        !
        tmp = location(1) >= xC
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
        tmp = location(2) >= yC
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
        tmp = location(3) >= zC
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
                call errStop( "interpFunc_rVector3D_SG: Unknown interp" )
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
        complex( kind=prec ), allocatable :: comp(:,:,:)
        !
        if( ( .NOT. self%is_allocated ) ) then
            call errStop( "interpFunc_rVector3D_SG > self not allocated." )
        endif
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
    function getArray_rVector3D_SG( self ) result( array )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        !
        complex( kind=prec ), allocatable, dimension(:) :: array
        !
        if( ( .NOT. self%is_allocated ) ) then
            call errStop( "getArray_rVector3D_SG > self not allocated." )
        endif
        !
        allocate( array( self%length() ) )
        !
        array = (/reshape( cmplx( self%x, 0.0, kind=prec ), (/self%Nxyz(1), 1/)), &
        reshape( cmplx( self%y, 0.0, kind=prec ), (/self%Nxyz(2), 1/)), &
        reshape( cmplx( self%z, 0.0, kind=prec ), (/self%Nxyz(3), 1/))/)
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
            call errStop( "setArray_rVector3D_SG > self not allocated." )
        endif
        !
        !> Ex
        i1 = 1; i2 = self%Nxyz(1)
        !
        self%x = reshape( real( array(i1:i2), kind=prec ), self%NdX )
        !
        !> Ey
        i1 = i2 + 1; i2 = i2 + self%Nxyz(2)
        !
        self%y = reshape( real( array(i1:i2), kind=prec ), self%NdY )
        !
        !> Ez
        i1 = i2 + 1; i2 = i2 + self%Nxyz(3)
        !
        self%z = reshape( real( array(i1:i2), kind=prec ), self%NdZ )
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
                self%x = rhs%x
                self%y = rhs%y
                self%z = rhs%z
                !
                self%is_allocated = .TRUE.
                !
            class default
                call errStop( "copyFrom_rVector3D_SG > Different type of rhs" )
            !
        end select
        !
    end subroutine copyFrom_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine getReal_rVector3D_SG( self, r_vector )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        class( Vector_t ), allocatable, intent( out ) :: r_vector
        !
        allocate( r_vector, source = self )
        !
        call warning( "getReal_rVector3D_SG > Getting Real Field from already Real Field" )
        !
    end subroutine getReal_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine edgeLength_rVector3D_SG( self, edge_length )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        type( rVector3D_SG_t ), intent( inout ) :: edge_length
        !
        integer :: ix, iy, iz
        !
        edge_length = self
        !
        select type( grid => self%grid )
            !
            class is( Grid3D_SG_t )
                !
                ! x-component edge length elements
                do ix = 1, grid%nx
                    edge_length%x(ix, :, :) = grid%dx(ix)
                enddo
                !
                ! y-component edge length elements
                do iy = 1, grid%ny
                    edge_length%y(:, iy, :) = grid%dy(iy)
                enddo
                !
                ! z-component edge length elements
                do iz = 1, grid%nz
                    edge_length%z(:, :, iz) = grid%dz(iz)
                enddo
                !
            class default
                call errStop( "edgeLength_rVector3D_SG > Undefined grid" )
        end select
        !
    end subroutine edgeLength_rVector3D_SG
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
        binary = .TRUE.
        !
        inquire( funit, opened = ok, named = hasname, name = fname, unformatted = isbinary )
        !
        if( ok ) then
            !
            !> Check that the file is unformatted if binary, formatted if ascii.
            if((index(isbinary, "yes") > 0 .OR.index(isbinary, "YES") > 0) &
                  .AND.   .NOT. binary ) then
                call errStop( "read_rVector3D_SG > Unable to read from unformatted file ["//trim(fname)//"]." )
            elseif((index(isbinary, "no") > 0 .OR.index(isbinary, "NO") > 0) &
                  .AND.binary) then
                call errStop( "read_rVector3D_SG > Unable to read from formatted file ["//trim(fname)//"]." )
            endif
            !
            read(funit) Nx, Ny, Nz, grid_type
            !
            if(  .NOT. self%is_allocated ) then
                call errStop( "read_rVector3D_SG > self not allocated." )
            elseif( self%grid_type .NE. grid_type ) then
                call errStop( "read_rVector3D_SG > Incompatible grid_type." )
            elseif( ( self%nx .NE. Nx ).OR. &
                    ( self%ny .NE. Ny ).OR.( self%nz .NE. Nz ) ) then
                call errStop( "read_rVector3D_SG > Wrong size on input from ["//trim(fname)//"]." )
            endif
            !
            read(funit) self%x
            read(funit) self%y
            read(funit) self%z
            !
        else
            call errStop( "read_rVector3D_SG: unable to open file" )
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
        if( .NOT. self%is_allocated ) then
            call errStop( "write_rVector3D_SG > self not allocated." )
        endif
        !
        binary = .TRUE.
        !
        inquire( funit, opened = ok, named = hasname, name = fname, unformatted = isbinary )
        !
        if( ok ) then
            !
            !> Check that the file is unformatted if binary, formatted if ascii.
            if( ( index( isbinary, "yes" ) > 0 .OR. index(isbinary, "YES") > 0 ) &
                  .AND.   .NOT. binary) then
                call errStop( "write_rVector3D_SG > Unable to write to unformatted file ["//trim(fname)//"]." )
            elseif( ( index( isbinary,"no" ) > 0 .OR. index( isbinary, "NO" ) > 0 ) &
                  .AND.binary) then
                call errStop( "write_rVector3D_SG > Unable to write to formatted file ["//trim(fname)//"]." )
            endif
            !
            write(funit) self%nx, self%ny, self%nz, self%grid_type
            write(funit) self%x
            write(funit) self%y
            write(funit) self%z
            !
        else
            call errStop( "write_rVector3D_SG > unable to open file" )
        endif
        !
    end subroutine write_rVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine print_rVector3D_SG( self, io_unit, title, append )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        integer, intent( in ), optional :: io_unit
        character(*), intent( in ), optional :: title
        logical, intent( in ), optional :: append
        !
        integer :: ix, iy, iz, funit
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
end module rVector3D_SG
!