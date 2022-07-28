!
module rVector3D_SG
    !
    use MatUtils
    use Grid3D_SG
    use Vector
    use rScalar3D_SG
    use cScalar3D_SG
    !
    type, extends( Vector_t ) :: rVector3D_SG_t
        !
        real( kind=prec ), allocatable, dimension(:, :, :) :: x, y, z
        !
    contains
        !
        final :: rVector3D_SG_dtor
        !
        procedure, public :: read  => readRVector3D_SG
        procedure, public :: write => writeRVector3D_SG
        !
        procedure, public :: setAllBoundary => setAllBoundaryRVector3D_SG
        procedure, public :: setOneBoundary => setOneBoundaryRVector3D_SG
        procedure, public :: setAllInterior => setAllInteriorRVector3D_SG
        procedure, public :: intBdryIndices => intBdryIndicesRVector3D_SG
        procedure, public :: boundary => boundaryRVector3D_SG
        procedure, public :: interior => interiorRVector3D_SG
        !
        procedure, public :: length => lengthRVector3D_SG
        !
        procedure, public :: getRealArray    => getRealArrayRVector3D_SG
        procedure, public :: getComplexArray => getComplexArrayRVector3D_SG
        !
        procedure, public :: setRealArray    => setRealArrayRVector3D_SG
        procedure, public :: setComplexArray => setComplexArrayRVector3D_SG
        !
        procedure, public :: setVecComponents => setVecComponentsRVector3D_SG
        !
        procedure, public :: zeros => zerosRVector3D_SG
        procedure, public :: add => addRVector3D_SG
        procedure, public :: sub => subRVector3D_SG
        procedure, public :: multByField => multByFieldRVector3D_SG
        procedure, public :: multByValue => multByValueRVector3D_SG
        procedure, public :: divByField => divByFieldRVector3D_SG
        procedure, public :: divByValue => divByValueRVector3D_SG
        procedure, public :: dotProd => dotProdRVector3D_SG
        procedure, public :: diagMult => diagMultRVector3D_SG
        !
        procedure, public :: copyFrom => copyFromRVector3D_SG
        !
        procedure, public :: sumEdges => sumEdgesRVector3D_SG
        procedure, public :: sumCells => sumCellsRVector3D_SG
        !
        procedure, public :: linCombS   => linCombSRVector3D_SG
        procedure, public :: scMultAddS => scMultAddSRVector3D_SG
        procedure, public :: interpFunc => interpFuncRVector3D_SG
        !
        procedure, public :: print => printRVector3D_SG
        !
    end type rVector3D_SG_t
    !
    interface rVector3D_SG_t
        module procedure rVector3D_SG_ctor
    end interface rVector3D_SG_t
    !
contains
    !
    function rVector3D_SG_ctor( igrid, grid_type ) result ( self )
        implicit none
        !
        class( Grid3D_SG_t ), target, intent( in ) :: igrid
        character( len=4 ), intent( in )           :: grid_type
        !
        type( rVector3D_SG_t ) :: self
        !
        integer :: nx, ny, nz, nzAir, nz_earth, status
        !
        !write( *, * ) "Constructor rVector3D_SG"
        !
        call self%init()
        !
        self%grid => igrid
        !
        ! Grid dimensions
        call igrid%GetDimensions(nx, ny, nz, nzAir)
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
            self%is_allocated = status.eq.0
            !
            allocate(self%y(self%nx + 1, self%ny, self%nz + 1), STAT = status)
            self%is_allocated = self%is_allocated.AND.(status.eq.0)
            !
            allocate(self%z(self%nx + 1, self%ny + 1, self%nz), STAT = status)
            self%is_allocated = self%is_allocated.AND.(status.eq.0)
            !
            self%NdX = (/self%nx, self%ny + 1, self%nz + 1/)
            self%NdY = (/self%nx + 1, self%ny, self%nz + 1/)
            self%NdZ = (/self%nx + 1, self%ny + 1, self%nz/)
            !
        else if(self%grid_type == FACE) then
            allocate(self%x(self%nx + 1, self%ny, self%nz), STAT = status)
            self%is_allocated = status.eq.0
            !
            allocate(self%y(self%nx, self%ny + 1, self%nz), STAT = status)
            self%is_allocated = self%is_allocated.AND.(status.eq.0)
            !
            allocate(self%z(self%nx, self%ny, self%nz + 1), STAT = status)
            self%is_allocated = self%is_allocated.AND.(status.eq.0)
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
    end function rVector3D_SG_ctor
    !
    subroutine rVector3D_SG_dtor( self )
        implicit none
        !
        type( rVector3D_SG_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor rVector3D_SG"
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
    subroutine readRVector3D_SG( self, funit, ftype )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        integer, intent( in )                    :: funit
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
            ! Check that the file is unformatted if binary, formatted if ascii.
            if((index(isbinary, "yes") > 0 .OR.index(isbinary, "YES") > 0) &
                  .AND.  .NOT.binary ) then
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
            if( .NOT.self%is_allocated) then
                write( *, * ) "Error: readRVector3D_SG > Vector must be allocated before readRVector3D_SGing from ", &
                        trim(fname), "."
                stop
            else if(self%grid_type.ne.grid_type) then
                write( *, * ) "Error: readRVector3D_SG > Vector must be of type ", grid_type, &
                        &            "           before readRVector3D_SGing from ", trim (fname), "."
                stop
            else if((self%nx.ne.Nx).OR. &
                  (self%ny.ne.Ny).OR.(self%nz.ne.Nz)) then
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
    subroutine writeRVector3D_SG( self, funit, ftype )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        integer, intent( in )                 :: funit
        character(:), allocatable, intent( in ), optional :: ftype
        !
        logical :: ok, hasname, binary
        character(80) :: fname, isbinary
        !
        if( .NOT. self%is_allocated) then
            stop "Error: writeRVector3D_SG > Not allocated."
        endif
        !
        binary = .TRUE.
        !
        inquire( funit, opened = ok, named = hasname, name = fname, unformatted = isbinary )
        !
        if( ok ) then
            !
            ! Check that the file is unformatted if binary, formatted if ascii.
            if((index(isbinary, "yes") > 0.OR.index(isbinary, "YES") > 0) &
                  .AND.  .NOT.binary) then
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
    subroutine setAllBoundaryRVector3D_SG( self, cvalue )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in )       :: cvalue
        !
        select case(self%grid_type)
            case(EDGE)
                self%x(:, (/1, self%NdX(2)/), :) = cvalue
                self%x(:, :, (/1, self%NdX(3)/)) = cvalue
                self%y((/1, self%NdY(1)/), :, :) = cvalue
                self%y(:, :, (/1, self%NdY(3)/)) = cvalue
                self%z(:, (/1, self%NdZ(2)/), :) = cvalue
                self%z((/1, self%NdZ(1)/), :, :) = cvalue
                !
            case(FACE)
                self%x((/1, self%NdX(1)/), :, :) = cvalue
                self%y(:, (/1, self%NdY(2)/), :) = cvalue
                self%z(:, :, (/1, self%NdZ(3)/)) = cvalue
                !
            case default
                stop "Error: setAllBoundaryRVector3D_SG > Invalid grid type."
        end select
        !
    end subroutine setAllBoundaryRVector3D_SG
    !
    subroutine setAllInteriorRVector3D_SG( self, cvalue )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in )       :: cvalue
        !
        select case(self%grid_type)
            case(EDGE)
                self%x(:, 2:self%NdX(2)-1, :) = cvalue
                self%x(:, :, 2:self%NdX(3)-1) = cvalue
                self%y(2:self%NdY(1)-1, :, :) = cvalue
                self%y(:, :, 2:self%NdY(3)-1) = cvalue
                self%z(:, 2:self%NdZ(2)-1, :) = cvalue
                self%z(2:self%NdZ(1)-1, :, :) = cvalue
            case(FACE)
                self%x(2:self%NdX(1)-1, :, :) = cvalue
                self%y(:, 2:self%NdY(2)-1, :) = cvalue
                self%z(:, :, 2:self%NdZ(3)-1) = cvalue
            case default
                stop "Error: setAllInteriorRVector3D_SG > Invalid grid type."
        end select
        !
    end subroutine setAllInteriorRVector3D_SG
    !
    subroutine setOneBoundaryRVector3D_SG( self, bdry, cvalue, int_only )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        character(:), allocatable, intent( in )  :: bdry
        complex( kind=prec ), intent( in )       :: cvalue
        logical, intent( in ), optional          :: int_only
        !
        logical :: int_only_p
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
                  select case (bdry)
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
                  end select
                else
                  select case (bdry)
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
                  end select
                endif
                !
            case(FACE)
                select case (bdry)
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
                end select
            case default
                stop "Error: setOneBoundaryRVector3D_SG > Invalid grid type."
        end select
        !
    end subroutine setOneBoundaryRVector3D_SG
    !
    subroutine intBdryIndicesRVector3D_SG( self, ind_i, ind_b )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        integer, allocatable, intent( out )   :: ind_i(:), ind_b(:)
        !
        integer :: nVec(3), nVecT, nBdry, nb, ni, i
        real( kind=prec ), allocatable :: temp(:)
        type( rVector3D_SG_t ) :: E
        !
        if(self%is_allocated) then
            select type(grid => self%grid)
                class is( Grid3D_SG_t )
                  E = rVector3D_SG_t( grid, self%grid_type )
                class default
                    stop "Error: intBdryIndicesRVector3D_SG: undefined grid"
            end select
            !
        else
            stop "Error: intBdryIndicesRVector3D_SG > Not allocated."
        endif
        !
        select case(self%grid_type)
            case(EDGE)
                nVec(1) = size(E%x)
                nVec(2) = size(E%y)
                nVec(3) = size(E%z)
                nVecT = nVec(1) + nVec(2) + nVec(3)
                !
                allocate(temp(nVecT))
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
                call E%getArray(temp)
                !
            case(FACE)
                nVec(1) = size(E%x)
                nVec(2) = size(E%y)
                nVec(3) = size(E%z)
                nVecT = nVec(1) + nVec(2) + nVec(3)
                !
                allocate(temp(nVecT))
                !
                E%x(1, :, :) = 1
                E%x(E%nx + 1, :, :) = 1
                E%y(:, 1, :) = 1
                E%y(:, E%ny + 1, :) = 1
                E%z(:, :, 1) = 1
                E%z(:, :, E%nz + 1) = 1
                !
                call E%getArray(temp) 
                !
        end select
        !
        nBdry = 0
        do i = 1, nVecT
            nBdry = nBdry + nint(temp(i))
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
            if(nint(temp(i)).eq.1) then
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
    subroutine boundaryRVector3D_SG( self, boundary )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        !
        class( Vector_t ), allocatable :: boundary
        !
        allocate( boundary, source = self )
        !
        call boundary%setAllInterior( C_ZERO )
       !
    end subroutine boundaryRVector3D_SG
    !
    subroutine interiorRVector3D_SG( self, interior )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        !
        class( Vector_t ), allocatable :: interior
        !
        allocate( interior, source = self )
        !
        call interior%setAllboundary( C_ZERO )
        !
    end subroutine interiorRVector3D_SG
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
    subroutine getRealArrayRVector3D_SG( self, array )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in )         :: self
        real( kind=prec ), allocatable, intent( out ) :: array(:)
        allocate (array(self%length()))
        array = (/reshape(self%x, (/self%Nxyz(1), 1/)), &
                reshape(self%y, (/self%Nxyz(2), 1/)), &
                reshape(self%z, (/self%Nxyz(3), 1/))/)
        !
    end subroutine getRealArrayRVector3D_SG
    !
    subroutine getComplexArrayRVector3D_SG( self, array )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in )            :: self
        complex( kind=prec ), allocatable, intent( out ) :: array(:)
        !
        allocate (array(self%length()))
        array = (/reshape(self%x, (/self%Nxyz(1), 1/)), &
                reshape(self%y, (/self%Nxyz(2), 1/)), &
                reshape(self%z, (/self%Nxyz(3), 1/))/)
        !
    end subroutine getComplexArrayRVector3D_SG
    !
    subroutine setRealArrayRVector3D_SG( self, array )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        real( kind=prec ), intent( in )          :: array(:)
        !
        integer :: i1, i2
        ! Ex
        i1 = 1; i2 = self%Nxyz(1)
        self%x = reshape(array(i1:i2), self%NdX)
        ! Ey
        i1 = i2 + 1; i2 = i2 + self%Nxyz(2)
        self%y = reshape(array(i1:i2), self%NdY)
        ! Ez
        i1 = i2 + 1; i2 = i2 + self%Nxyz(3)
        self%z = reshape(array(i1:i2), self%NdZ)
        !
    end subroutine setRealArrayRVector3D_SG
    !
    subroutine setComplexArrayRVector3D_SG( self, array )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in )       :: array(:)
        !
        integer :: i1, i2
        ! Ex
        i1 = 1; i2 = self%Nxyz(1)
        self%x = reshape(array(i1:i2), self%NdX)
        ! Ey
        i1 = i2 + 1; i2 = i2 + self%Nxyz(2)
        self%y = reshape(array(i1:i2), self%NdY)
        ! Ez
        i1 = i2 + 1; i2 = i2 + self%Nxyz(3)
        self%z = reshape(array(i1:i2), self%NdZ)
        !
    end subroutine setComplexArrayRVector3D_SG
    !
    subroutine setVecComponentsRVector3D_SG( self, xyz, &
            &                              xmin, xstep, xmax, &
            &                              ymin, ystep, ymax, &
            &                              zmin, zstep, zmax, c )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        character, intent( in )                :: xyz
        integer, intent( in )                  :: xmin, xstep, xmax
        integer, intent( in )                  :: ymin, ystep, ymax
        integer, intent( in )                  :: zmin, zstep, zmax
        real( kind=prec ), intent ( in )       :: c
        !
        integer :: x1, x2
        integer :: y1, y2
        integer :: z1, z2
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
                self%x(x1:x2:xstep, y1:y2:ystep, z1:z2:zstep) = c
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
                self%y(x1:x2:xstep, y1:y2:ystep, z1:z2:zstep) = c
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
                self%z(x1:x2:xstep, y1:y2:ystep, z1:z2:zstep) = c
                !
            case default
                stop "Error: setVecComponentsRVector3D_SG > Invalid xyz argument."
        end select
        !
    end subroutine setVecComponentsRVector3D_SG
    !
    subroutine zerosRVector3D_SG( self )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        !
        self%x = R_ZERO
        self%y = R_ZERO
        self%z = R_ZERO
        !
    end subroutine zerosRVector3D_SG
    !
    subroutine addRVector3D_SG( self, rhs )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in )           :: rhs
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type( rhs )
                !
                class is( rVector3D_SG_t )
                    self%x = self%x + rhs%x
                    self%y = self%y + rhs%y
                    self%z = self%z + rhs%z
                class default
                    stop "Error: addRVector3D_SG > Undefined rhs"
                !
            end select
            !
        else
            stop "Error: addRVector3D_SG > Incompatible inputs."
        end if
        !
    end subroutine addRVector3D_SG
    !
    subroutine subRVector3D_SG( self, rhs )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in )           :: rhs
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type( rhs )
                !
                class is( rVector3D_SG_t )
                    self%x = self%x - rhs%x
                    self%y = self%y - rhs%y
                    self%z = self%z - rhs%z
                class default
                    stop "Error: subRVector3D_SG > Undefined rhs"
                !
            end select
            !
        else
            stop "Error: subRVector3D_SG > Incompatible inputs."
        endif
        !
    end subroutine subRVector3D_SG
    !
    subroutine multByFieldRVector3D_SG( self, rhs )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in )           :: rhs
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type( rhs )
                !
                class is( rVector3D_SG_t )
                    self%x = self%x * rhs%x
                    self%y = self%y * rhs%y
                    self%z = self%z * rhs%z
                !
                class is( rScalar3D_SG_t )
                    self%x = self%x * rhs%v
                    self%y = self%y * rhs%v
                    self%z = self%z * rhs%v
                !
                class is( cScalar3D_SG_t )
                    self%x = self%x * rhs%v
                    self%y = self%y * rhs%v
                    self%z = self%z * rhs%v
                    !
                class default
                    stop "Error: multByFieldRVector3D_SG: undefined rhs"
                !
            end select
            !
        else
            stop "Error: multByFieldRVector3D_SG: incompatible rhs"
        end if
        !
    end subroutine multByFieldRVector3D_SG
    !
    subroutine multByValueRVector3D_SG( self, cvalue )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in )       :: cvalue
        !
        self%x = cvalue * self%x
        self%y = cvalue * self%y
        self%z = cvalue * self%z
        !
    end subroutine multByValueRVector3D_SG
    !
    subroutine divByFieldRVector3D_SG( self, rhs )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in )           :: rhs
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type( rhs )
                !
                class is( rVector3D_SG_t )
                    self%x = self%x / rhs%x
                    self%y = self%y / rhs%y
                    self%z = self%z / rhs%z
                !
                class is( rScalar3D_SG_t )
                    self%x = self%x / rhs%v
                    self%y = self%y / rhs%v
                    self%z = self%z / rhs%v
                !
                class is( cScalar3D_SG_t )
                    self%x = self%x / rhs%v
                    self%y = self%y / rhs%v
                    self%z = self%z / rhs%v
                    !
                class default
                    stop "Error: divByFieldRVector3D_SG: undefined rhs"
                !
            end select
            !
        else
            stop "Error: divByFieldRVector3D_SG: incompatible rhs"
        end if
        !
    end subroutine divByFieldRVector3D_SG
    !
    subroutine divByValueRVector3D_SG( self, cvalue )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in )       :: cvalue
        !
        self%x = self%x / cvalue
        self%y = self%y / cvalue
        self%z = self%z / cvalue
        !
    end subroutine divByValueRVector3D_SG
    !
    function dotProdRVector3D_SG( self, rhs ) result( cvalue )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        class( Vector_t ), intent( in )       :: rhs
        complex( kind=prec )                  :: cvalue
        !
        cvalue = C_ZERO
        !
        if(( .NOT. self%is_allocated) .OR. ( .NOT. rhs%is_allocated)) then
            stop "Error: dotProdRVector3D_SG > Input vectors not allocated."
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type( rhs )
                class is( rVector3D_SG_t )
                    cvalue = cvalue + sum( self%x * rhs%x )
                    cvalue = cvalue + sum( self%y * rhs%y )
                    cvalue = cvalue + sum( self%z * rhs%z )
                class default
                    stop "Error: dotProdRVector3D_SG: undefined rhs"
            end select
            !
        else
            stop "Error: dotProdRVector3D_SG > Incompatible rhs"
        end if
        !
    end function dotProdRVector3D_SG
    !
    function diagMultRVector3D_SG( self, rhs ) result ( diag_mult )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        class( Vector_t ), intent( in )      :: rhs
        class( Vector_t ), allocatable       :: diag_mult
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type( grid => self%grid )
                class is( Grid3D_SG_t )
                    !
                    diag_mult = rVector3D_SG_t( grid, self%grid_type )
                    !
                    select type( diag_mult )
                        class is( rVector3D_SG_t )
                            !
                            select type( rhs )
                                !
                                class is( rVector3D_SG_t )
                                    diag_mult%x = self%x * rhs%x
                                    diag_mult%y = self%y * rhs%y
                                    diag_mult%z = self%z * rhs%z
                                class default
                                    stop "Error: diagMultRVector3D_SG > Undefined rhs"
                                !
                            end select
                            !
                    end select
                class default
                    stop "Error: diagMultRVector3D_SG: undefined grid"
                    !
            end select
        else
            stop "Error: diagMultRVector3D_SG > Incompatible inputs."
        endif
        !
    end function diagMultRVector3D_SG
    !
    function sumEdgesRVector3D_SG( self, interior_only ) result( cell_obj )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        logical, optional, intent( in )       :: interior_only
        class( Scalar_t ), allocatable        :: cell_obj
        !
        integer :: x_xend, x_yend, x_zend
        integer :: y_xend, y_yend, y_zend
        integer :: z_xend, z_yend, z_zend
        !
        type( rVector3D_SG_t ) :: E_tmp
        !
        E_tmp = self
        !
        if( interior_only ) then
            call E_tmp%setAllBoundary( C_ZERO )
        endif
        !
        select type( grid => self%grid )
            class is( Grid3D_SG_t )
                !
                cell_obj = rScalar3D_SG_t( grid, CELL )
                !
            class default
                stop "Error: sumEdgesRVector3D_SG: undefined grid"
                !
        end select
        !
        select type( cell_obj )
            class is(rScalar3D_SG_t)
                select case(E_tmp%grid_type)
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
                  case(FACE)
                      x_xend = size(E_tmp%x, 1)
                      y_xend = size(E_tmp%y, 1)
                      z_xend = size(E_tmp%z, 1)
                      !
                      cell_obj%v = E_tmp%x(1:x_xend-1,:,:) + E_tmp%x(2:x_xend,:,:) + &
                                E_tmp%y(:,1:y_yend-1,:) + E_tmp%y(:,2:y_yend,:) + &
                                E_tmp%z(:,:,1:z_zend-1) + E_tmp%z(:,:,2:z_zend)
                end select
        end select
        !
    end function sumEdgesRVector3D_SG
    !
    subroutine sumCellsRVector3D_SG( self, E_in, ptype )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Scalar_t ), intent( in )          :: E_in
        character(*), intent( in ), optional     :: ptype
        !
        character(10) :: type
        integer :: xend, yend, zend
        integer :: v_xend, v_yend, v_zend
        integer :: ix, iy, iz
        !
        if(index(self%grid_type, CELL) > 0) then
            stop "Error: sumCellsRVector3D_SG > Only CELL type supported."
        endif
        !
        if( .NOT. present(ptype)) then
            type = EDGE
        else
            type = ptype
        endif
        !
        select type( E_in )
            class is( rScalar3D_SG_t )
                !
                v_xend = size(E_in%v, 1)
                v_yend = size(E_in%v, 2)
                v_zend = size(E_in%v, 3)
                !
                select case(type)
                    case(EDGE)

                        ! for x-components inside the domain
                        do ix = 1, self%grid%nx
                           do iy = 2, self%grid%ny
                              do iz = 2, self%grid%nz
                                 self%x(ix, iy, iz) = (E_in%v(ix, iy-1, iz-1) + E_in%v(ix, iy, iz-1) + &
                                      E_in%v(ix, iy-1, iz) + E_in%v(ix, iy, iz))/4.0d0
                              end do
                           end do
                        end do
                        
                        ! for y-components inside the domain
                        do ix = 2, self%grid%nx
                           do iy = 1, self%grid%ny
                              do iz = 2, self%grid%nz
                                 self%y(ix, iy, iz) = (E_in%v(ix-1, iy, iz-1) + E_in%v(ix, iy, iz-1) + &
                                      E_in%v(ix-1, iy, iz) + E_in%v(ix, iy, iz))/4.0d0
                              end do
                           end do
                        end do
                        
                        do ix = 2, self%grid%nx
                              do iy = 2, self%grid%ny
                                 do iz = 1, self%grid%nz
                                    self%z(ix, iy, iz) = (E_in%v(ix-1, iy-1, iz) + E_in%v(ix-1, iy, iz) + &
                                         E_in%v(ix, iy-1, iz) + E_in%v(ix, iy, iz))/4.0d0
                                 end do
                              end do
                           end do
                          ! upper boundary
                        iz = 1
                        do iy = 1, self%grid%ny
                           do ix = 2, self%grid%nx
                              self%y(ix, iy, iz) = (E_in%v(ix-1, iy, iz) + E_in%v(ix, iy, iz))/2.0d0
                           end do
                        end do
                        do ix = 1, self%grid%nx
                           do iy = 2,self%grid%ny
                              self%x(ix, iy, iz) = (E_in%v(ix, iy-1, iz) + E_in%v(ix, iy, iz))/2.0d0
                           end do
                        end do
                        
                        ! lower boundary
                        iz = self%grid%nz + 1
                        do iy = 1, self%grid%ny
                           do ix = 2, self%grid%nx
                              self%y(ix, iy, iz) = (E_in%v(ix-1, iy, iz-1) + E_in%v(ix, iy, iz-1))/2.0d0
                           end do
                        end do
                        do ix = 1, self%grid%nx
                           do iy = 2, self%grid%ny
                              self%x(ix, iy, iz) = (E_in%v(ix, iy-1, iz-1) + E_in%v(ix, iy, iz-1))/2.0d0
                           end do
                        end do  
                        
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
                stop "Error: sumCellsRVector3D_SG > Incompatible input Scalar_t."
        end select
        !
    end subroutine sumCellsRVector3D_SG
    !
    subroutine linCombSRVector3D_SG( self, rhs, c1, c2 )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Vector_t ), intent( in )          :: rhs
        complex( kind=prec ), intent( in )       :: c1, c2
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type(rhs)
                !
                class is( rVector3D_SG_t )
                    self%x = c1 * self%x + c2 * rhs%x
                    self%y = c1 * self%y + c2 * rhs%y
                    self%z = c1 * self%z + c2 * rhs%z
                class default
                    stop "Error: linCombSRVector3D_SG > rhs undefined."
            end select
            !
        else
            stop "Error: linCombSRVector3D_SG > Incompatible inputs."
        end if
        !
    end subroutine linCombSRVector3D_SG
    !
    subroutine scMultAddSRVector3D_SG( self, rhs, cvalue )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        class( Vector_t ), intent( inout )   :: rhs
        complex( kind=prec ), intent( in )    :: cvalue
        !
        if ( self%isCompatible( rhs ) ) then
            !
            select type(rhs)
                class is( rVector3D_SG_t ) 
                    rhs%x = rhs%x + cvalue * self%x
                    rhs%y = rhs%y + cvalue * self%y
                    rhs%z = rhs%z + cvalue * self%z
            end select
            !
        else
            stop "Error: scMultAddSRVector3D_SG >Incompatible inputs."
        end if
        !
    end subroutine scMultAddSRVector3D_SG
    !
    subroutine interpFuncRVector3D_SG( self, location, xyz, interp )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        real( kind=prec ), intent( in )       :: location(3)
        character, intent( in )               :: xyz
        class( Vector_t ), intent( inout )    :: interp
        !
        real( kind=prec ), allocatable, dimension(:) :: xC, yC, zC
        integer :: ix, iy, iz, i
        real( kind=prec ) :: wx, wy, wz
        logical, dimension(:), allocatable :: tmp
        !
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
                                allocate(xC(size(grid%delX)))
                                allocate(yC(size(grid%dy + 1)))
                                allocate(zC(size(grid%dz + 1)))
                                !
                                xC = CumSum(grid%delX)
                                yC = CumSum([0._prec, grid%dy])
                                zC = CumSum([0._prec, grid%dz])
                            case("y")
                                allocate(xC(size(grid%dx + 1)))
                                allocate(yC(size(grid%delY)))
                                allocate(zC(size(grid%dz)))
                                !
                                xC = CumSum([0._prec, grid%dx])
                                yC = CumSum([grid%delY])
                                zC = CumSum([0._prec, grid%dz])
                            case("z")
                                allocate(xC(size(grid%dx + 1)))
                                allocate(yC(size(grid%dy + 1)))
                                allocate(zC(size(grid%delZ)))
                                !
                                xC = CumSum([0._prec, grid%dx])
                                yC = CumSum([0._prec, grid%dy])
                                zC = CumSum([grid%delZ])
                        end select !xyz
                        !
                    case( FACE )
                        !
                        interp = rVector3D_SG_t( grid, FACE )
                        !
                        select case(xyz)
                            case("x")
                                allocate(xC(size(grid%dx + 1)))
                                allocate(yC(size(grid%delY)))
                                allocate(zC(size(grid%delZ)))
                                !
                                xC = CumSum([0._prec, grid%dx])
                                yC = CumSum([grid%delY])
                                zC = CumSum([grid%delZ])
                            case("y")
                                allocate(xC(size(grid%delX)))
                                allocate(yC(size(grid%dy + 1)))
                                allocate(zC(size(grid%delZ)))
                                !
                                xC = CumSum([grid%delX])
                                yC = CumSum([0._prec, grid%dy])
                                zC = CumSum([grid%delZ])
                            case("z")
                                allocate(xC(size(grid%delX)))
                                allocate(yC(size(grid%delY)))
                                allocate(zC(size(grid%dz + 1)))
                                !
                                xC = CumSum([grid%delX])
                                yC = CumSum([grid%delY])
                                zC = CumSum([0._prec, grid%dz])
                        end select !xyz
                        !
                end select !GRID TYPE
                !
                xC = xC - grid%ox
                yC = yC - grid%oy
                zC = zC - sum(grid%dz(1:grid%nzAir)) - grid%oz
                !
            class default
                stop "Error: interpFuncRVector3D_SG: undefined grid"
                !
        end select !GRID
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
        end select
        !
    end subroutine interpFuncRVector3D_SG
    !
    subroutine copyFromRVector3D_SG( self, rhs )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in )           :: rhs
        !
        if( .NOT. rhs%is_allocated) then
            stop "Error: copyFromRVector3D_SG > rhs not allocated"
        endif
        !
        select type( rhs )
            class is( rVector3D_SG_t )
                !
                self%grid => rhs%grid
                self%grid_type = rhs%grid_type
                self%nx = rhs%nx
                self%ny = rhs%ny
                self%nz = rhs%nz
                self%NdX = rhs%NdX
                self%NdY = rhs%NdY
                self%NdZ = rhs%NdZ
                self%Nxyz = rhs%Nxyz
                self%x = rhs%x
                self%y = rhs%y
                self%z = rhs%z
                self%is_allocated = .TRUE.
                !
            class default
                stop "Error: copyFromRVector3D_SG > Incompatible rhs"
        end select
        !
    end subroutine copyFromRVector3D_SG
    !
    subroutine printRVector3D_SG( self, io_unit, title, append )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in )             :: self
        integer, intent( in ), optional                   :: io_unit
        character(:), allocatable, intent( in ), optional :: title
        logical, intent( in ), optional                   :: append
        !
        integer :: ix, iy, iz,funit
        !
        if( present( io_unit ) ) then
            funit = io_unit
        else
            funit = 0
        endif
        !
        !
        write( funit, * ) title
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
end module rVector3D_SG
