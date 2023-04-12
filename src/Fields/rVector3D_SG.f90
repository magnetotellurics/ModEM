!
!> Derived class to define a rVector3D_SG_t 
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
        real( kind=prec ), allocatable, dimension(:) :: sv
        !
        contains
            !
            final :: rVector3D_SG_dtor
            !
            !> Boundary operations
            procedure, public :: setAllBoundary => setAllBoundaryRVector3D_SG
            procedure, public :: setOneBoundary => setOneBoundaryRVector3D_SG
            procedure, public :: intBdryIndices => intBdryIndicesRVector3D_SG
            !
            !> Dimensioning operations
            procedure, public :: length => lengthRVector3D_SG
            procedure, public :: setVecComponents => setVecComponentsRVector3D_SG
            !
            !> Arithmetic/algebraic unary operations
            procedure, public :: zeros => zerosRVector3D_SG
            procedure, public :: sumEdges => sumEdgesRVector3D_SG
            procedure, public :: avgCells => avgCellsRVector3D_SG
            procedure, public :: conjugate => conjugateRVector3D_SG
            !
            !> Arithmetic/algebraic binary operations
            procedure, public :: add => addRVector3D_SG
            !
            procedure, public :: linComb => linCombRVector3D_SG
            !
            procedure, public :: subValue => subValueRVector3D_SG
            procedure, public :: subField => subFieldRVector3D_SG
            !
            procedure, public :: multByReal => multByRealRVector3D_SG
            procedure, public :: multByComplex => multByComplexRVector3D_SG
            procedure, public :: multByField => multByFieldRVector3D_SG
            !
            procedure, public :: diagMult => diagMultRVector3D_SG
            !
            procedure, public :: multAdd => multAddRVector3D_SG
            !
            procedure, public :: dotProd => dotProdRVector3D_SG
            !
            procedure, public :: divByValue => divByValueRVector3D_SG
            procedure, public :: divByField => divByFieldRVector3D_SG
            !
            procedure, public :: interpFunc => interpFuncRVector3D_SG
            !
            !> Miscellaneous
            procedure, public :: getReal => getRealRVector3D_SG
            procedure, public :: setArray => setArrayRVector3D_SG
            procedure, public :: getArray => getArrayRVector3D_SG
            procedure, public :: switchStoreState => switchStoreStateRVector3D_SG
            procedure, public :: copyFrom => copyFromRVector3D_SG
            !
            !> I/O operations
            procedure, public :: print => printRVector3D_SG
            procedure, public :: read => readRVector3D_SG
            procedure, public :: write => writeRVector3D_SG
            !
    end type rVector3D_SG_t
    !
    public :: getRvector, setRvector, EdgeLength
    !
    interface rVector3D_SG_t
        module procedure rVector3D_SG_ctor
    end interface rVector3D_SG_t
    !
contains
    !
    !> No subroutine briefing
    !
    function rVector3D_SG_ctor( igrid, grid_type ) result( self )
        implicit none
        !
        class( Grid_t ), target, intent( in ) :: igrid
        character( len=4 ), intent( in ) :: grid_type
        !
        type( rVector3D_SG_t ) :: self
        !
        integer :: nx, ny, nz, nzAir, nz_earth, status
        !
        !write( *, * ) "Constructor rVector3D_SG"
        !
        call self%init
        !
        self%grid => igrid
        !
        !> Grid dimensions
        call igrid%getDimensions(nx, ny, nz, nzAir)
        nz_earth = nz - nzAir
        !
        self%nx = nx
        self%ny = ny
        self%nz = nz
        !
        self%grid_type = trim( grid_type )
        self%is_allocated = .FALSE.
        !
        if(self%grid_type == EDGE) then
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
            stop "Error: rVector3D_SG_ctor > Only EDGE or FACE types allowed."
        endif
        !
        if(self%is_allocated) then
            self%x = R_ZERO
            self%y = R_ZERO
            self%z = R_ZERO
        else
            stop "Error: rVector3D_SG_ctor > Unable to allocate vector."
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
        call self%dealloc
        !
        if( allocated( self%x ) ) deallocate( self%x )
        if( allocated( self%y ) ) deallocate( self%y )
        if( allocated( self%z ) ) deallocate( self%z )
        if( allocated( self%sv ) ) deallocate( self%sv )
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
    subroutine setAllBoundaryRVector3D_SG( self, cvalue )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        if( self%store_state /= compound ) then
             call self%switchStoreState
        endif
        !
        select case(self%grid_type)
            case(EDGE)
                self%x(:, (/1, self%NdX(2)/), :) = real( cvalue, kind=prec )
                self%x(:, :, (/1, self%NdX(3)/)) = real( cvalue, kind=prec )
                self%y((/1, self%NdY(1)/), :, :) = real( cvalue, kind=prec )
                self%y(:, :, (/1, self%NdY(3)/)) = real( cvalue, kind=prec )
                self%z(:, (/1, self%NdZ(2)/), :) = real( cvalue, kind=prec )
                self%z((/1, self%NdZ(1)/), :, :) = real( cvalue, kind=prec )
                !
            case(FACE)
                self%x((/1, self%NdX(1)/), :, :) = real( cvalue, kind=prec )
                self%y(:, (/1, self%NdY(2)/), :) = real( cvalue, kind=prec )
                self%z(:, :, (/1, self%NdZ(3)/)) = real( cvalue, kind=prec )
                !
            case default
                stop "Error: setAllBoundaryRVector3D_SG > Invalid grid type."
        end select
        !
    end subroutine setAllBoundaryRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setOneBoundaryRVector3D_SG( self, bdry, cvalue, int_only )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
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
        if( .NOT. present(int_only)) then
            int_only_p = .FALSE.
        else
            int_only_p = int_only
        endif
        !
        select case(self%grid_type)
            case(EDGE)
                if(int_only_p) then
                  select case(bdry)
                      case("x1")
                          self%z(1, 2:self%NdZ(2)-1, :) = real( cvalue, kind=prec )
                          self%y(1, :, 2:self%NdY(3)-1) = real( cvalue, kind=prec )
                      case("x2")
                          self%z(self%NdZ(1), 2:self%NdZ(2)-1, :) = real( cvalue, kind=prec )
                          self%y(self%NdY(1), :, 2:self%NdY(3)-1) = real( cvalue, kind=prec )
                      case("y1")
                          self%z(2:self%NdZ(1)-1, 1, :) = real( cvalue, kind=prec )
                          self%x(:, 1, 2:self%NdX(3)-1) = real( cvalue, kind=prec )
                      case("y2")
                          self%z(2:self%NdZ(1)-1, self%NdZ(2), :) = real( cvalue, kind=prec )
                          self%x(:, self%NdX(2), 2:self%NdX(3)-1) = real( cvalue, kind=prec )
                      case("z1")
                          self%x(:, 2:self%NdX(2)-1, 1) = real( cvalue, kind=prec )
                          self%y(2:self%NdY(1)-1, :, 1) = real( cvalue, kind=prec )
                      case("z2")
                          self%x(:, 2:self%NdX(2)-1, self%NdX(3)) = real( cvalue, kind=prec )
                          self%y(2:self%NdY(1)-1, :, self%NdY(3)) = real( cvalue, kind=prec )
                      case("z1_x")
                          self%x(:, 2:self%NdX(2)-1, 1) = real( cvalue, kind=prec )
                      case("z2_x")
                          self%x(:, 2:self%NdX(2)-1, self%NdX(3)) = real( cvalue, kind=prec )
                      case("z1_y")
                          self%y(2:self%NdY(1)-1, :, 1) = real( cvalue, kind=prec )
                      case("z2_y")
                          self%y(2:self%NdY(1)-1, :, self%NdY(3)) = real( cvalue, kind=prec )
                  end select
                else
                  select case(bdry)
                      case("x1")
                          self%z(1, :, :) = real( cvalue, kind=prec )
                          self%y(1, :, :) = real( cvalue, kind=prec )
                      case("x2")
                          self%z(self%NdZ(1), :, :) = real( cvalue, kind=prec )
                          self%y(self%NdY(1), :, :) = real( cvalue, kind=prec )
                      case("y1")
                          self%z(:, 1, :) = real( cvalue, kind=prec )
                          self%x(:, 1, :) = real( cvalue, kind=prec )
                      case("y2")
                          self%z(:, self%NdZ(2), :) = real( cvalue, kind=prec )
                          self%x(:, self%NdX(2), :) = real( cvalue, kind=prec )
                      case("z1")
                          self%x(:, :, 1) = real( cvalue, kind=prec )
                          self%y(:, :, 1) = real( cvalue, kind=prec )
                      case("z2")
                          self%x(:, :, self%NdX(3)) = real( cvalue, kind=prec )
                          self%y(:, :, self%NdY(3)) = real( cvalue, kind=prec )
                      case("z1_x")
                          self%x(:, :, 1) = real( cvalue, kind=prec )
                      case("z2_x")
                          self%x(:, :, self%NdX(3)) = real( cvalue, kind=prec )
                      case("z1_y")
                          self%y(:, :, 1) = real( cvalue, kind=prec )
                      case("z2_y")
                          self%y(:, :, self%NdY(3)) = real( cvalue, kind=prec )
                  end select
                endif
                !
            case(FACE)
                select case(bdry)
                  case("x1")
                      self%x(1, :, :) = real( cvalue, kind=prec )
                  case("x2")
                      self%x(self%NdX(1), :, :) = real( cvalue, kind=prec )
                  case("y1")
                      self%y(:, 1, :) = real( cvalue, kind=prec )
                  case("y2")
                      self%y(:, self%NdY(2), :) = real( cvalue, kind=prec )
                  case("z1")
                      self%z(:, :, 1) = real( cvalue, kind=prec )
                  case("z2")
                      self%z(:, :, self%NdZ(3)) = real( cvalue, kind=prec )
                end select
            case default
                stop "Error: setOneBoundaryRVector3D_SG > Invalid grid type."
        end select
        !
    end subroutine setOneBoundaryRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine intBdryIndicesRVector3D_SG( self, ind_i, ind_b )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        integer, allocatable, intent( out ) :: ind_i(:), ind_b(:)
        !
        integer :: nVecT, nBdry, nb, ni, i
        complex( kind=prec ), allocatable :: temp(:)
        type( rVector3D_SG_t ) :: E
        !
        if( self%is_allocated ) then
            !
            E = rVector3D_SG_t( self%grid, self%grid_type )
            !
        else
            stop "Error: intBdryIndicesRVector3D_SG > Not allocated. Exiting."
        endif
        !
        if( self%store_state /= compound ) then
             call self%switchStoreState
        endif
        !
        select case(self%grid_type)
            case(EDGE)
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
            case(FACE)
                !
                E%x(1, :, :) = 1
                E%x(E%nx + 1, :, :) = 1
                E%y(:, 1, :) = 1
                E%y(:, E%ny + 1, :) = 1
                E%z(:, :, 1) = 1
                E%z(:, :, E%nz + 1) = 1
                !
        end select
        !
        temp = E%getArray()
        nVecT = size( E%x ) + size( E%y ) + size( E%z )
        nBdry = 0
        do i = 1, nVecT
            nBdry = nBdry + nint( real( temp(i), kind=prec ) )
        enddo
        !
        if( allocated(ind_i) ) deallocate(ind_i)
        allocate(ind_i(nVecT - nBdry))
        !
        if( allocated(ind_b) ) deallocate(ind_b)
        allocate(ind_b(nBdry))
        !
        nb = 0
        ni = 0
        do i = 1, nVecT
            if( nint( real( temp(i), kind=prec ) ) .EQ. 1 ) then
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
    end subroutine intBdryIndicesRVector3D_SG
    !
    !> No subroutine briefing
    !
    function lengthRVector3D_SG( self ) result( n )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        !
        integer :: n
        !
        n = self%Nxyz(1) + self%Nxyz(2) + self%Nxyz(3)
        !
    end function lengthRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setVecComponentsRVector3D_SG( self, xyz, &
            &                              xmin, xstep, xmax, &
            &                              ymin, ystep, ymax, &
            &                              zmin, zstep, zmax, rvalue )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        character, intent( in ) :: xyz
        integer, intent( in ) :: xmin, xstep, xmax
        integer, intent( in ) :: ymin, ystep, ymax
        integer, intent( in ) :: zmin, zstep, zmax
        real( kind=prec ), intent ( in ) :: rvalue
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
        select case(xyz)
            case("x")
                if(xmin == 0) x1 = self%NdX(1)
                if(xmax <= 0) x2 = self%NdX(1) + xmax
                !
                if(ymin == 0) y1 = self%NdX(2)
                if(ymax <= 0) y2 = self%NdX(2) + ymax
                !
                if(zmin == 0) z1 = self%NdX(3)
                if(zmax <= 0) z2 = self%NdX(3) + zmax
                !
                self%x(x1:x2:xstep, y1:y2:ystep, z1:z2:zstep) = rvalue
                !
            case("y")
                if(xmin == 0) x1 = self%NdY(1)
                if(xmax <= 0) x2 = self%NdY(1) + xmax
                !
                if(ymin == 0) y1 = self%NdY(2)
                if(ymax <= 0) y2 = self%NdY(2) + ymax
                !
                if(zmin == 0) z1 = self%NdY(3)
                if(zmax <= 0) z2 = self%NdY(3) + zmax
                !
                self%y(x1:x2:xstep, y1:y2:ystep, z1:z2:zstep) = rvalue
                !
            case("z")
                if(xmin == 0) x1 = self%NdZ(1)
                if(xmax <= 0) x2 = self%NdZ(1) + xmax
                !
                if(ymin == 0) y1 = self%NdZ(2)
                if(ymax <= 0) y2 = self%NdZ(2) + ymax
                !
                if(zmin == 0) z1 = self%NdZ(3)
                if(zmax <= 0) z2 = self%NdZ(3) + zmax
                !
                self%z(x1:x2:xstep, y1:y2:ystep, z1:z2:zstep) = rvalue
                !
            case default
                stop "Error: setVecComponentsRVector3D_SG > Invalid xyz argument."
        end select
        !
    end subroutine setVecComponentsRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine zerosRVector3D_SG( self )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        !
        if( .NOT. self%is_allocated) then
             stop "Error: zerosRVector3D_SG > Not allocated."
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
            self%sv = R_ZERO
            !
        else
            stop "Error: zerosRVector3D_SG > Unknown store_state!"
        endif
        !
    end subroutine zerosRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine sumEdgesRVector3D_SG( self, cell_obj, interior_only )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Field_t ), allocatable, intent( inout ) :: cell_obj
        logical, optional, intent( in ) :: interior_only
        !
        integer :: x_xend, x_yend, x_zend
        integer :: y_xend, y_yend, y_zend
        integer :: z_xend, z_yend, z_zend
        type( rVector3D_SG_t ) :: E_tmp
        logical :: is_interior_only
        !
        if( self%store_state /= compound ) then
             call self%switchStoreState
        endif
        !
        E_tmp = self
        !
        is_interior_only = .FALSE.
        !
        if( present( interior_only ) ) is_interior_only = interior_only
        !
        if( is_interior_only ) then
            call E_tmp%setAllBoundary( C_ZERO )
        endif
        !
        if( allocated( cell_obj ) ) then
            cell_obj = rScalar3D_SG_t( self%grid, CELL )
        else
            allocate( cell_obj, source = rScalar3D_SG_t( self%grid, CELL ) )
        endif
        !
        select type( cell_obj )
            !
            class is( rScalar3D_SG_t )
                !
                select case( E_tmp%grid_type )
                    !
                    case(EDGE)
                        x_xend = size(E_tmp%x, 1)
                        x_yend = size(E_tmp%x, 2)
                        x_zend = size(E_tmp%x, 3)
                        !
                        y_xend = size(E_tmp%y, 1)
                        y_yend = size(E_tmp%y, 2)
                        y_zend = size(E_tmp%y, 3)
                        !
                        z_xend = size(E_tmp%z, 1)
                        z_yend = size(E_tmp%z, 2)
                        z_zend = size(E_tmp%z, 3)
                        !
                        cell_obj%v = E_tmp%x(:,1:x_yend-1,1:x_zend-1) + &
                        E_tmp%x(:,2:x_yend,1:x_zend-1)       + &
                        E_tmp%x(:,1:x_yend-1,2:x_zend)       + &
                        E_tmp%x(:,2:x_yend,2:x_zend)         + &
                        E_tmp%y(1:y_xend-1,:,1:y_zend-1)     + &
                        E_tmp%y(2:y_xend,:,1:y_zend-1)       + &
                        E_tmp%y(1:y_xend-1,:,2:y_zend)       + &
                        E_tmp%y(2:y_xend,:,2:y_zend)         + &
                        E_tmp%z(1:z_xend-1,1:z_yend-1,:)     + &
                        E_tmp%z(2:z_xend,1:z_yend-1,:)       + &
                        E_tmp%z(1:z_xend-1,2:z_yend,:)       + &
                        E_tmp%z(2:z_xend,2:z_yend,:)
                        !
                    case(FACE)
                        x_xend = size(E_tmp%x, 1)
                        y_xend = size(E_tmp%y, 1)
                        z_xend = size(E_tmp%z, 1)
                        !
                        cell_obj%v =    E_tmp%x(1:x_xend-1,:,:) + E_tmp%x(2:x_xend,:,:) + &
                                        E_tmp%y(:,1:y_yend-1,:) + E_tmp%y(:,2:y_yend,:) + &
                                        E_tmp%z(:,:,1:z_zend-1) + E_tmp%z(:,:,2:z_zend)
                    !
                    case default
                        stop "Error: sumEdgesRVector3D_SG: undefined E_tmp%grid_type"
                    !
                  end select
                  !
            class default
                stop "Error: sumEdgesRVector3D_SG: undefined grid"
                !
        end select
        !
    end subroutine sumEdgesRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine avgCellsRVector3D_SG( self, E_in, ptype )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: E_in
        character(*), intent( in ), optional :: ptype
        !
        character(10) :: type
        integer :: xend, yend, zend
        integer :: v_xend, v_yend, v_zend
        integer :: ix, iy, iz
        !
        if( index( self%grid_type, CELL ) > 0 ) then
            stop "Error: avgCellsRVector3D_SG > Only CELL type supported."
        endif
        !
        if( .NOT. present(ptype)) then
            type = EDGE
        else
            type = ptype
        endif
        !
        if( self%store_state /= compound ) then
             call self%switchStoreState
        endif
        !
        select type( E_in )
            class is( rScalar3D_SG_t )
                !
                v_xend = size( E_in%v, 1 )
                v_yend = size( E_in%v, 2 )
                v_zend = size( E_in%v, 3 )
                !
                select case(type)
                    case(EDGE)

                        !> for x-components inside the domain
                        do ix = 1, self%grid%nx
                           do iy = 2, self%grid%ny
                              do iz = 2, self%grid%nz
                                 self%x(ix, iy, iz) = (E_in%v(ix, iy-1, iz-1) + E_in%v(ix, iy, iz-1) + &
                                      E_in%v(ix, iy-1, iz) + E_in%v(ix, iy, iz)) / 4.0d0
                              enddo
                           enddo
                        enddo
                        
                        !> for y-components inside the domain
                        do ix = 2, self%grid%nx
                           do iy = 1, self%grid%ny
                              do iz = 2, self%grid%nz
                                 self%y(ix, iy, iz) = (E_in%v(ix-1, iy, iz-1) + E_in%v(ix, iy, iz-1) + &
                                      E_in%v(ix-1, iy, iz) + E_in%v(ix, iy, iz)) / 4.0d0
                              enddo
                           enddo
                        enddo
                        
                        do ix = 2, self%grid%nx
                              do iy = 2, self%grid%ny
                                 do iz = 1, self%grid%nz
                                    self%z(ix, iy, iz) = (E_in%v(ix-1, iy-1, iz) + E_in%v(ix-1, iy, iz) + &
                                         E_in%v(ix, iy-1, iz) + E_in%v(ix, iy, iz)) / 4.0d0
                                 enddo
                              enddo
                           enddo
                        !
                    case(FACE)
                        xend = size(self%x, 1)
                        self%x(2:xend-1,:,:) = E_in%v(1:v_xend-1,:,:) + E_in%v(2:v_xend,:,:)
                        !
                        yend = size(self%y, 1)
                        self%y(:, 2:yend-1, :) = E_in%v(:, 1:v_yend-1, :) + E_in%v(:, 2:v_yend, :)
                        !
                        zend = size(self%z, 1) 
                        self%z(:, :, 2:zend-1) = E_in%v(:, :, 1:v_zend-1) + E_in%v(:, :, 2:v_zend)
                end select
            class default
                stop "Error: avgCellsRVector3D_SG > Incompatible input Scalar_t."
        end select
        !
    end subroutine avgCellsRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine conjugateRVector3D_SG( self )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        !
        stop "Error: conjugateRVector3D_SG: Do not try to conjugate a real vector!"
        !
    end subroutine conjugateRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine addRVector3D_SG( self, rhs )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( .NOT. rhs%is_allocated) then
             stop "Error: addRVector3D_SG > rhs not allocated."
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            if( self%store_state /= rhs%store_state ) call self%switchStoreState
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
                        self%sv = self%sv + rhs%sv
                        !
                    else
                        stop "Error: addRVector3D_SG > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: addRVector3D_SG > Undefined rhs"
                !
            end select
            !
        else
            stop "Error: addRVector3D_SG > Incompatible inputs."
        endif
        !
    end subroutine addRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine linCombRVector3D_SG( self, rhs, c1, c2 )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        complex( kind=prec ), intent( in ) :: c1, c2
        !
        if( self%isCompatible( rhs ) ) then
            !
            if( self%store_state /= rhs%store_state ) call self%switchStoreState
            !
            select type(rhs)
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
                        self%sv = c1 * self%sv + c2 * rhs%sv
                        !
                    else
                        stop "Error: linCombRVector3D_SG > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: linCombRVector3D_SG > rhs undefined."
            end select
            !
        else
            stop "Error: linCombRVector3D_SG > Incompatible inputs."
        endif
        !
    end subroutine linCombRVector3D_SG
    !
    !> No subroutine briefing
    subroutine subValueRVector3D_SG( self, cvalue )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        if( self%store_state .EQ. compound ) then
            !
            self%x = self%x - cvalue
            self%y = self%y - cvalue
            self%z = self%z - cvalue
            !
        else if( self%store_state .EQ. singleton ) then
            !
            self%sv = self%sv - cvalue
            !
        else
            stop "Error: subValueRVector3D_SG > Unknown self store_state!"
        endif
        !
    end subroutine subValueRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine subFieldRVector3D_SG( self, rhs )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( self%isCompatible( rhs ) ) then
            !
            if( self%store_state /= rhs%store_state ) call self%switchStoreState
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
                        self%sv = self%sv - rhs%sv
                        !
                    else
                        stop "Error: subFieldRVector3D_SG > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: subFieldRVector3D_SG > Undefined rhs"
                !
            end select
            !
        else
            stop "Error: subFieldRVector3D_SG > Incompatible inputs."
        endif
        !
    end subroutine subFieldRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multByRealRVector3D_SG( self, rvalue )
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
            self%sv = self%sv * rvalue
            !
        else
            stop "Error: multByRealRVector3D_SG > Unknown self store_state!"
        endif
        !
    end subroutine multByRealRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multByComplexRVector3D_SG( self, cvalue )
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
            self%sv = self%sv * cvalue
            !
        else
            stop "Error: multByComplexRVector3D_SG > Unknown self store_state!"
        endif
        !
    end subroutine multByComplexRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multByFieldRVector3D_SG( self, rhs )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( self%isCompatible( rhs ) ) then
            !
            if( self%store_state /= rhs%store_state ) call self%switchStoreState
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
                        self%sv = self%sv * rhs%sv
                        !
                    else
                        stop "Error: multByFieldRVector3D_SG > Unknown rhs store_state!"
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
                        self%sv = self%sv * rhs%sv
                        !
                    else
                        stop "Error: multByFieldRVector3D_SG > Unknown rhs store_state!"
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
                        self%sv = self%sv * rhs%sv
                        !
                    else
                        stop "Error: multByFieldRVector3D_SG > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: multByFieldRVector3D_SG: undefined rhs"
                !
            end select
            !
        else
            stop "Error: multByFieldRVector3D_SG: incompatible rhs"
        endif
        !
    end subroutine multByFieldRVector3D_SG
    !
    !> No subroutine briefing
    !
    function diagMultRVector3D_SG( self, rhs ) result( diag_mult )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Vector_t ), intent( in ) :: rhs
        !
        class( Vector_t ), allocatable :: diag_mult
        !
        if( self%isCompatible( rhs ) ) then
            !
            allocate( diag_mult, source = rVector3D_SG_t( self%grid, self%grid_type ) )
            !
            if( diag_mult%store_state /= rhs%store_state ) call diag_mult%switchStoreState
            !
            if( self%store_state /= rhs%store_state ) call self%switchStoreState
            !
            select type( diag_mult )
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
                                diag_mult%sv = self%sv * rhs%sv
                                !
                            else
                                stop "Error: diagMultRVector3D_SG > Unknown rhs store_state!"
                            endif
                            !
                        class default
                            stop "Error: diagMultRVector3D_SG > Undefined rhs"
                        !
                    end select
                !
                class default
                    stop "Error: diagMultRVector3D_SG > Undefined diag_mult"
                !
            end select
            !
        else
            stop "Error: diagMultRVector3D_SG > Incompatible inputs."
        endif
        !
    end function diagMultRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multAddRVector3D_SG( self, cvalue, rhs )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        class( Field_t ), intent( in ) :: rhs
        !
        if( self%isCompatible( rhs ) ) then
            !
            if( self%store_state /= rhs%store_state ) call self%switchStoreState
            !
            select type(rhs)
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
                        self%sv = self%sv + cvalue * rhs%sv
                        !
                    else
                        stop "Error: multAddRVector3D_SG > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: multAddRVector3D_SG > rhs undefined."
                !
            end select
            !
        else
            stop "Error: multAddRVector3D_SG >Incompatible inputs."
        endif
        !
    end subroutine multAddRVector3D_SG
    !
    !> No subroutine briefing
    !
    function dotProdRVector3D_SG( self, rhs ) result( cvalue )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        complex( kind=prec ) :: cvalue
        !
        cvalue = C_ZERO
        !
        if(( .NOT. self%is_allocated ) .OR. ( .NOT. rhs%is_allocated )) then
            stop "Error: dotProdRVector3D_SG > Input vectors not allocated."
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            if( self%store_state == rhs%store_state ) then
                !
                select type( rhs )
                    !
                    class is( rVector3D_SG_t )
                        !
                        if( rhs%store_state .EQ. compound ) then
                            !
                            cvalue = cvalue + cmplx( sum( self%x * rhs%x ), 0.0, kind=prec )
                            cvalue = cvalue + cmplx( sum( self%y * rhs%y ), 0.0, kind=prec )
                            cvalue = cvalue + cmplx( sum( self%z * rhs%z ), 0.0, kind=prec )
                            !
                        else if( rhs%store_state .EQ. singleton ) then
                            !
                            cvalue = cvalue + cmplx( sum( self%sv * rhs%sv ), 0.0, kind=prec )
                            !
                        else
                            stop "Error: dotProdRVector3D_SG > Unknown rhs store_state!"
                        endif
                        !
                    class default
                        stop "Error: dotProdRVector3D_SG: undefined rhs"
                    !
                end select
                !
            else
                stop "Error: dotProdRVector3D_SG > Incompatible store_state"
            endif
            !
        else
            stop "Error: dotProdRVector3D_SG > Incompatible rhs"
        endif
        !
    end function dotProdRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine divByValueRVector3D_SG( self, cvalue )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        if( self%store_state .EQ. compound ) then
            !
            self%x = self%x / cvalue
            self%y = self%y / cvalue
            self%z = self%z / cvalue
            !
        else if( self%store_state .EQ. singleton ) then
            !
            self%sv = self%sv / cvalue
            !
        else
            stop "Error: divByValueRVector3D_SG > Unknown self store_state!"
        endif
        !
    end subroutine divByValueRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine divByFieldRVector3D_SG( self, rhs )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( self%isCompatible( rhs ) ) then
            !
            if( self%store_state /= rhs%store_state ) call self%switchStoreState
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
                        self%sv = self%sv / rhs%sv
                        !
                    else
                        stop "Error: divByFieldRVector3D_SG > Unknown rhs store_state!"
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
                        self%sv = self%sv / rhs%sv
                        !
                    else
                        stop "Error: divByFieldRVector3D_SG > Unknown rhs store_state!"
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
                        self%sv = self%sv / rhs%sv
                        !
                    else
                        stop "Error: divByFieldRVector3D_SG > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: divByFieldRVector3D_SG: undefined rhs"
                !
            end select
            !
        else
            stop "Error: divByFieldRVector3D_SG: incompatible rhs"
        endif
        !
    end subroutine divByFieldRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine interpFuncRVector3D_SG( self, location, xyz, interp )
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
        select type( grid => self%grid )
            class is( Grid3D_SG_t )
                !
                select case( self%grid_type )
                    !
                    case( EDGE )
                        !
                        interp = rVector3D_SG_t( grid, EDGE )
                        !
                        select case( xyz )
                            case("x")
                                allocate(xC(size(grid%del_x)))
                                allocate(yC(size(grid%dy + 1)))
                                allocate(zC(size(grid%dz + 1)))
                                !
                                xC = CumSum(grid%del_x)
                                yC = CumSum([0._prec, grid%dy])
                                zC = CumSum([0._prec, grid%dz])
                            case("y")
                                allocate(xC(size(grid%dx + 1)))
                                allocate(yC(size(grid%del_y)))
                                allocate(zC(size(grid%dz)))
                                !
                                xC = CumSum([0._prec, grid%dx])
                                yC = CumSum([grid%del_y])
                                zC = CumSum([0._prec, grid%dz])
                            case("z")
                                allocate(xC(size(grid%dx + 1)))
                                allocate(yC(size(grid%dy + 1)))
                                allocate(zC(size(grid%del_z)))
                                !
                                xC = CumSum([0._prec, grid%dx])
                                yC = CumSum([0._prec, grid%dy])
                                zC = CumSum([grid%del_z])
                        end select !xyz
                        !
                    case( FACE )
                        !
                        interp = rVector3D_SG_t( grid, FACE )
                        !
                        select case(xyz)
                            case("x")
                                allocate(xC(size(grid%dx + 1)))
                                allocate(yC(size(grid%del_y)))
                                allocate(zC(size(grid%del_z)))
                                !
                                xC = CumSum([0._prec, grid%dx])
                                yC = CumSum([grid%del_y])
                                zC = CumSum([grid%del_z])
                            case("y")
                                allocate(xC(size(grid%del_x)))
                                allocate(yC(size(grid%dy + 1)))
                                allocate(zC(size(grid%del_z)))
                                !
                                xC = CumSum([grid%del_x])
                                yC = CumSum([0._prec, grid%dy])
                                zC = CumSum([grid%del_z])
                            case("z")
                                allocate(xC(size(grid%del_x)))
                                allocate(yC(size(grid%del_y)))
                                allocate(zC(size(grid%dz + 1)))
                                !
                                xC = CumSum([grid%del_x])
                                yC = CumSum([grid%del_y])
                                zC = CumSum([0._prec, grid%dz])
                        end select !xyz
                        !
                    case default
                        stop "Error: multAddRVector3D_SG > self%grid_type undefined."
                    !
                end select !GRID TYPE
                !
            class default
                stop "Error: interpFuncRVector3D_SG: undefined grid"
                !
        end select !GRID
        !
        xC = xC - self%grid%ox
        yC = yC - self%grid%oy
        zC = zC - sum(self%grid%dz(1:self%grid%nzAir)) - self%grid%oz
        !
        tmp = location(1) > xC
        do i = size(tmp), 1, -1 
            if(tmp(i)) then
                ix = i
                exit
            endif
        enddo
        tmp = location(2) > yC
        do i = size(tmp), 1, -1 
            if(tmp(i)) then
                iy = i
                exit
            endif
        enddo
        tmp = location(3) > zC
        do i = size(tmp), 1, -1 
            if(tmp(i)) then
                iz = i
                exit
            endif
        enddo
        !
        deallocate( tmp )
        !
        !ix = findloc(location(1) > xC, .TRUE., back = .TRUE., dim = 1)
        !iy = findloc(location(2) > yC, .TRUE., back = .TRUE., dim = 1)
        !iz = findloc(location(3) > zC, .TRUE., back = .TRUE., dim = 1)
        !> Find weights
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
            class is( rVector3D_SG_t )
                select case( xyz )
                    case("x")
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
                        interp%z(ix,iy,iz) = wx*wy*wz
                        interp%z(ix+1,iy,iz) = (1-wx)*wy*wz
                        interp%z(ix,iy+1,iz) = wx*(1-wy)*wz
                        interp%z(ix,iy,iz+1) = wx*wy*(1-wz)
                        interp%z(ix,iy+1,iz+1) = wx*(1-wy)*(1-wz)
                        interp%z(ix+1,iy,iz+1) = (1-wx)*wy*(1-wz)
                        interp%z(ix+1,iy+1,iz) = (1-wx)*(1-wy)*wz
                        interp%z(ix+1,iy+1,iz+1) = (1-wx)*(1-wy)*(1-wz)
                end select
                !
            class default
                stop "Error: interpFuncRVector3D_SG: undefined interp"
                !
        end select
        !
    end subroutine interpFuncRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine getRealRVector3D_SG( self, r_field )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        class( Field_t ), allocatable, intent( out ) :: r_field
        !
        allocate( r_field, source = rVector3D_SG_t( self%grid, self%grid_type ) )
        !
        call r_field%copyFrom( self )
        !
    end subroutine getRealRVector3D_SG
    !
    !> No subroutine briefing
    !
    function getArrayRVector3D_SG( self ) result( array )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        !
        complex( kind=prec ), allocatable, dimension(:) :: array
        !
        if( self%store_state .EQ. compound ) then
            !
            allocate( array( self%length() ) )
            !
            array = (/reshape( cmplx( self%x, 0.0, kind=prec ), (/self%Nxyz(1), 1/) ), &
                    reshape( cmplx( self%y, 0.0, kind=prec ), (/self%Nxyz(2), 1/) ), &
                    reshape( cmplx( self%z, 0.0, kind=prec ), (/self%Nxyz(3), 1/) )/)
            !
        else if( self%store_state .EQ. singleton ) then
            !
            array = cmplx( self%sv, 0.0, kind=prec )
            !
        else
            stop "Error: getArrayRVector3D_SG > Unknown store_state!"
        endif
        !
		write( *, * ) "getArrayRVector3D_SG: ", size( array )
		!
    end function getArrayRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setArrayRVector3D_SG( self, array )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), dimension(:), intent( in ) :: array
        !
        integer :: i1, i2
        !
        if( self%store_state .EQ. compound ) then
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
        else if( self%store_state .EQ. singleton ) then
            !
            self%sv = array
            !
        else
            stop "Error: setArrayRVector3D_SG > Unknown store_state!"
        endif
        !
    end subroutine setArrayRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine switchStoreStateRVector3D_SG( self )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        !
        integer i1, i2
        !
        select case( self%store_state )
            !
            case( compound )
                !
                allocate( self%sv( self%length() ) )
                !
                self%sv = &
                (/reshape( self%x, (/self%Nxyz(1), 1/) ), &
                  reshape( self%y, (/self%Nxyz(2), 1/) ), &
                  reshape( self%z, (/self%Nxyz(3), 1/) )/)
                !
                deallocate( self%x )
                deallocate( self%y )
                deallocate( self%z )
                !
                self%store_state = singleton
                !
            case( singleton )
                !
                if( self%grid_type == EDGE ) then
                    !
                    allocate( self%x( self%nx, self%ny + 1, self%nz + 1 ) )
                    allocate( self%y( self%nx + 1, self%ny, self%nz + 1 ) )
                    allocate( self%z( self%nx + 1, self%ny + 1, self%nz ) )
                    !
                else if( self%grid_type == FACE ) then
                    !
                    allocate( self%x( self%nx + 1, self%ny, self%nz ) )
                    allocate( self%y( self%nx, self%ny + 1, self%nz ) )
                    allocate( self%z( self%nx, self%ny, self%nz + 1 ) )
                    !
                else
                    stop "Error: switchStoreStateRVector3D_SG > Only EDGE or FACE types allowed."
                endif
                !
                ! Ex
                i1 = 1; i2 = self%Nxyz(1)
                self%x = reshape( self%sv( i1 : i2 ), self%NdX )
                !
                ! Ey
                i1 = i2 + 1; i2 = i2 + self%Nxyz(2)
                self%y = reshape( self%sv( i1 : i2 ), self%NdY )
                !
                ! Ez
                i1 = i2 + 1; i2 = i2 + self%Nxyz(3)
                self%z = reshape( self%sv( i1 : i2 ), self%NdZ )
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
    end subroutine switchStoreStateRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine copyFromRVector3D_SG( self, rhs )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( .NOT. rhs%is_allocated ) then
            stop "Error: copyFromRVector3D_SG > rhs not allocated"
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
            class is( rVector3D_SG_t )
                !
                self%NdX = rhs%NdX
                self%NdY = rhs%NdY
                self%NdZ = rhs%NdZ
                self%Nxyz = rhs%Nxyz
                !
                if( rhs%store_state .EQ. compound ) then
                    !
                    self%x = rhs%x
                    self%y = rhs%y
                    self%z = rhs%z
                    !
                else if( rhs%store_state .EQ. singleton ) then
                    !
                    self%sv = rhs%sv
                    !
                else
                    stop "Error: copyFromRVector3D_SG > Unknown store_state!"
                endif
                !
            class default
                stop "Error: copyFromRVector3D_SG > Undefined rhs"
        end select
        !
        self%is_allocated = .TRUE.
        !
    end subroutine copyFromRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine readRVector3D_SG( self, funit, ftype )
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
        !> Make sure the store_state is compound
        if( self%store_state /= compound ) then
             call self%switchStoreState
        endif
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
                write( *, * ) "Error: readRVector3D_SG > Unable to readRVector3D_SG vector from unformatted file. ", &
                        trim(fname), "."
                stop
            else if((index(isbinary, "no") > 0 .OR.index(isbinary, "NO") > 0) &
                  .AND.binary) then
                write( *, * ) "Error: readRVector3D_SG > Unable to readRVector3D_SG vector from formatted file ", &
                        trim(fname), "."
                stop
            endif
            !
            read(funit) Nx, Ny, Nz, grid_type
            !
            if(  .NOT. self%is_allocated) then
                write( *, * ) "Error: readRVector3D_SG > Vector must be allocated before readRVector3D_SGing from ", &
                        trim(fname), "."
                stop
            else if(self%grid_type.NE.grid_type) then
                write( *, * ) "Error: readRVector3D_SG > Vector must be of type ", grid_type, &
                        &            "           before readRVector3D_SGing from ", trim (fname), "."
                stop
            else if((self%nx.NE.Nx).OR. &
                  (self%ny.NE.Ny).OR.(self%nz.NE.Nz)) then
                write( *, * ) "Error: readRVector3D_SG > Wrong size of vector on input from ", trim (fname), "."
                stop
            endif
            !
            read(funit) self%x
            read(funit) self%y
            read(funit) self%z
            !
        else
            stop "Error: readRVector3D_SG: unable to open file"
        endif
        !
    end subroutine readRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine writeRVector3D_SG( self, funit, ftype )
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
            stop "Error: writeRVector3D_SG > Not allocated."
        endif
        !
        !> Make sure the store_state is compound
        if( self%store_state /= compound ) then
             call self%switchStoreState
        endif
        !
        binary = .TRUE.
        !
        inquire( funit, opened = ok, named = hasname, name = fname, unformatted = isbinary )
        !
        if( ok ) then
            !
            !> Check that the file is unformatted if binary, formatted if ascii.
            if((index(isbinary, "yes") > 0 .OR. index(isbinary, "YES") > 0) &
                  .AND.   .NOT. binary) then
                write( *, * ) "Error: writeRVector3D_SG > Unable to writeRVector3D_SG vector to unformatted file. ", &
                        trim(fname), "."
                stop
            else if((index(isbinary,"no") > 0.OR.index(isbinary,"NO") > 0) &
                  .AND.binary) then
                write( *, * ) "Error: writeRVector3D_SG > Unable to writeRVector3D_SG vector to formatted file. ", &
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
            stop "Error: writeRVector3D_SG > unable to open file"
        endif
        !
    end subroutine writeRVector3D_SG
    !
    !> No subroutine briefing
    !
    subroutine printRVector3D_SG( self, io_unit, title, append )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
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
            funit = 0
        endif
        !
        write( funit, * ) "ModEM-OO RVector3D_SG"
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
    end subroutine printRVector3D_SG
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
    end subroutine
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
    subroutine EdgeLength( grid, l_E )
        implicit none
        !
        class( Grid_t ), intent( in ) :: grid
        type( rVector3D_SG_t ), intent( inout )  :: l_E
        !
        integer :: ix, iy, iz
        !
        l_E = rVector3D_SG_t( grid, EDGE )
        !
        ! x-component edge length elements
        do ix = 1,grid%nx
            l_E%x(ix, :, :) = grid%dx(ix)
        enddo
        !
        ! y-component edge length elements
        do iy = 1,grid%ny
            l_E%y(:, iy, :) = grid%dy(iy)
        enddo
        !
        ! z-component edge length elements
        do iz = 1,grid%nz
            l_E%z(:, :, iz) = grid%dz(iz)
        enddo
        !
    end subroutine EdgeLength
    !
end module rVector3D_SG
