!**
! SUMMARY
!
! Standard cartesian grid vectors.
!*
module rVector3D_SG
    use Constants
    use MatUtils
    use Grid
    use Grid3D_SG
    use rScalar
    use rScalar3D_SG
    use rVector
    !**
    !*
    type, extends( rVector_t ) :: rVector3D_SG_t
        !**
        ! pointer to parent grid
        class( Grid3D_SG_t ), pointer :: grid
        !**
        ! Store the intention of the use in a character
        ! string defined in GridDef as a parameter: EDGE
        ! or FACE are two possibilities.
        !*
        character( len=4 ) :: gridType
        !**
        ! Grid Dimensions:
        ! nx is grid dimension (number of cells) in the x-direction
        ! ny is grid dimension (number of cells) in the y-direction
        ! nz is grid dimension (number of cells) in the z-direction:
        !*
        integer :: nx, ny, nz
        !
        integer, dimension(3) :: NdX, NdY, NdZ
        integer, dimension(3) :: Nxyz
        !
        !**
        ! Typical usage:    electrical fields on cell edges of
        ! staggered grid
        ! For example, in an EDGE, the dimensions would be
        ! x: edge nodes in x-direction: dimension Nx, Ny+1, Nz+1
        ! y: edge nodes in y-direction: dimension Nx+1, Ny, Nz+1
        ! z: edge nodes in z-direction: dimension Nx+1, Ny+1, Nz
        ! Note that the arrays are defined through dynamic
        ! memory allocation
        !*
        real( kind=prec ), allocatable, dimension(:, :, :) :: x, y, z
        !
    contains
        !**
        ! Initialization and finalization
        !*
        final :: rVector3D_SG_dtor
        !**
        ! Input/Output
        !*
        procedure, public :: read    => readRVector3D_SG
        procedure, public :: write => writeRVector3D_SG
        !**
        ! Boundary operations
        !*
        procedure, public :: setAllBoundary => setAllBoundaryRVector3D_SG
        procedure, public :: setOneBoundary => setOneBoundaryRVector3D_SG
        procedure, public :: setAllInterior => setAllInteriorRVector3D_SG
        procedure, public :: intBdryIndices => intBdryIndicesRVector3D_SG
        procedure, public :: boundary => boundaryRVector3D_SG
        procedure, public :: interior => interiorRVector3D_SG
        !**
        ! Data access
        !*
        procedure, public :: length => lengthRVector3D_SG
        procedure, public :: getArray => getArrayRVector3D_SG
        procedure, public :: setArray => setArrayRVector3D_SG
        procedure, public :: setVecComponents => setVecComponentsRVector3D_SG
        !**
        ! Arithmetic operations
        !*
        procedure, public :: zeros => zerosRVector3D_SG
        procedure, public :: add1 => add1RVector3D_SG
        procedure, public :: sub1 => sub1RVector3D_SG
        procedure, public :: mult1 => mult1RVector3D_SG
        procedure, public, pass( self ) :: mult2 => mult2RVector3D_SG
        procedure, public :: div1 => div1RVector3D_SG
        procedure, public :: dotProd => dotProdRVector3D_SG
        procedure, public :: diagMult => diagMultRVector3D_SG
        !**
        !     Subroutine versions -- first argument is overwritten
        !
        procedure, public :: divS1 => divS1RVector3D_SG
        procedure, public :: multS1 => multS1RVector3D_SG
        procedure, public :: multS2 => multS2RVector3D_SG
        !**
        ! Miscellaneous
        !*
        procedure, public :: isCompatible => isCompatibleRVector3D_SG
        procedure, public :: copyFrom => copyFromRVector3D_SG
        !
        procedure, public :: sumEdges => sumEdgesRVector3D_SG
        procedure, public :: SumCells => SumCellsRVector3D_SG
        !
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
    !**
    ! Constructor for REAL vectors.
    !
    ! Arguments
    !     igrid           Underlying grid.
    !     gridtype    Definied in GridDef.f90
    !
    !*
    function rVector3D_SG_ctor( igrid, gridType ) result ( self )
        implicit none
        !
        class( Grid3D_SG_t ), target, intent( in ) :: igrid
        character(*), intent( in )                 :: gridType
        !
        integer :: status
        type( rVector3D_SG_t ) :: self
        !
        !write( *, * ) "Constructor rVector3D_SG"
        !
        call self%init()
        !
        self%grid => igrid
        !
        ! Grid dimensions
        self%nx = igrid%nx
        self%ny = igrid%ny
        self%nz = igrid%nz
        !
        self%gridType = trim( gridType )
        self%is_allocated = .FALSE.
        !
        if(self%gridType == EDGE) then
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
        else if(self%gridType == FACE) then
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
            write( *, * ) "ERROR:rVector3D_SG_ctor:"
            stop "           Only EDGE or FACE types allowed. Exiting."
        endif
        !
        if(self%is_allocated) then
            self%x = R_ZERO
            self%y = R_ZERO
            self%z = R_ZERO
        else
            write( *, * ) "ERROR:rVector3D_SG_ctor"
            stop "           Unable to allocate vector. Exiting."
        endif
        !
        self%Nxyz = (/product(self%NdX), &
                product(self%NdY), &
                product(self%NdZ)/)
        !
    end function rVector3D_SG_ctor
    
    subroutine rVector3D_SG_dtor( self )
        implicit none
        !
        type( rVector3D_SG_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor rVector3D_SG"
        !
        !if( self%is_allocated ) then
            !
            if( allocated( self%x ) ) deallocate( self%x )
            if( allocated( self%y ) ) deallocate( self%y )
            if( allocated( self%z ) ) deallocate( self%z )
            !
            self%nx = 0
            self%ny = 0
            self%nz = 0
            !
            self%gridType = ""
            self%is_allocated = .FALSE.
            !
        !endif
        !
    end subroutine rVector3D_SG_dtor
    !
    !************************************************
    ! Input/Output
    !************************************************
    !
    !**
    ! readRVector3D_SG
    !*
    subroutine readRVector3D_SG( self, fid )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        integer, intent( in )                  :: fid
        !
        integer :: Nx, Ny, Nz
        character(4) :: gridType
        logical :: ok, hasname, binary
        character(80) :: fname, isbinary
        !
        binary = .TRUE.
        !
        inquire( fid, opened = ok, named = hasname, name = fname, unformatted = isbinary )
        !
        if( ok ) then
            !
            ! Check that the file is unformatted if binary, formatted if ascii.
            if((index(isbinary, "yes") > 0 .OR.index(isbinary, "YES") > 0) &
                  .AND.  .NOT.binary ) then
                write( *, * ) "ERROR:rVector3D_SG::readRVector3D_SG: "
                write( *, * ) "           Unable to readRVector3D_SG vector from unformatted file. ", &
                        trim(fname), ". Exiting."
                stop
            else if((index(isbinary, "no") > 0 .OR.index(isbinary, "NO") > 0) &
                  .AND.binary) then
                write( *, * ) "ERROR:rVector3D_SG::readRVector3D_SG: "
                write( *, * ) "           Unable to readRVector3D_SG vector from formatted file ", &
                        trim(fname), ". Exiting."
                stop
            endif
            !
            read(fid) Nx, Ny, Nz, gridType
            !
            if( .NOT.self%is_allocated) then
                write( *, * ) "ERROR:rVector3D_SG::readRVector3D_SG: "
                write( *, * ) "           Vector must be allocated before readRVector3D_SGing from ", &
                        trim(fname), ". Exiting."
                stop
            else if(self%gridType.ne.gridType) then
                write( *, * ) "ERROR:rVector3D_SG::readRVector3D_SG: "
                write( *, * ) "           Vector must be of type ", gridType, &
                        &            "           before readRVector3D_SGing from ", trim (fname), ". Exiting."
                stop
            else if((self%nx.ne.Nx).OR. &
                  (self%ny.ne.Ny).OR.(self%nz.ne.Nz)) then
                write( *, * ) "ERROR:rVector3D_SG_t::readRVector3D_SG: "
                write( *, * ) "           Wrong size of vector on input from ", trim (fname), ". Exiting."
                stop
            endif
            !
            read(fid) self%x
            read(fid) self%y
            read(fid) self%z
            !
        else
            stop "readRVector3D_SG: unable to open file"
        endif
        !
    end subroutine readRVector3D_SG
    
    !**
    ! writeRVector3D_SG
    !*
    subroutine writeRVector3D_SG( self, fid )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        integer, intent( in )                 :: fid
        !
        logical :: ok, hasname, binary
        character(80) :: fname, isbinary
        !
        if( .NOT. self%is_allocated) then
            write( *, * ) "ERROR:rVector3D_SG::writeRVector3D_SG: "
            stop "           Not allocated. Exiting."
        endif
        !
        binary = .TRUE.
        !
        inquire( fid, opened = ok, named = hasname, name = fname, unformatted = isbinary )
        !
        if( ok ) then
            !
            ! Check that the file is unformatted if binary, formatted if ascii.
            if((index(isbinary, "yes") > 0.OR.index(isbinary, "YES") > 0) &
                  .AND.  .NOT.binary) then
                write( *, * ) "ERROR:rVector3D_SG::writeRVector3D_SG: "
                write( *, * ) "           Unable to writeRVector3D_SG vector to unformatted file. ", &
                        trim(fname), ". Exiting."
                stop
            else if((index(isbinary,"no") > 0.OR.index(isbinary,"NO") > 0) &
                  .AND.binary) then
                write( *, * ) "ERROR:rVector3D_SG::writeRVector3D_SG: "
                write( *, * ) "           Unable to writeRVector3D_SG vector to formatted file. ", &
                        trim(fname), ". Exiting."
                stop
            endif
            !
            write(fid) self%nx, self%ny, self%nz, self%gridType
            write(fid) self%x
            write(fid) self%y
            write(fid) self%z
            !
        else
            stop "writeRVector3D_SG: unable to open file"
        endif
        !
    end subroutine writeRVector3D_SG
    !
    !************************************************
    ! Boundary operations
    !************************************************
    !
    !** 
    ! setAllBoundaryRVector3D_SG
    !
    !*
    subroutine setAllBoundaryRVector3D_SG( self, c_in )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        real( kind=prec ), intent( in )          :: c_in
        !
        select case(self%gridType)
            case(EDGE)
                self%x(:, (/1, self%NdX(2)/), :) = c_in
                self%x(:, :, (/1, self%NdX(3)/)) = c_in
                self%y((/1, self%NdY(1)/), :, :) = c_in
                self%y(:, :, (/1, self%NdY(3)/)) = c_in
                self%z(:, (/1, self%NdZ(2)/), :) = c_in
                self%z((/1, self%NdZ(1)/), :, :) = c_in
                !
            case(FACE)
                self%x((/1, self%NdX(1)/), :, :) = c_in
                self%y(:, (/1, self%NdY(2)/), :) = c_in
                self%z(:, :, (/1, self%NdZ(3)/)) = c_in
                !
            case default
                write( *, * ) "ERROR:rVector3D_SG::setAllBoundaryRVector3D_SG: "
                stop "           Invalid grid type. Exiting."
        end select
        !
    end subroutine setAllBoundaryRVector3D_SG
    !** 
    ! setAllInteriorRVector3D_SG
    !
    !*
    subroutine setAllInteriorRVector3D_SG( self, c_in )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        real( kind=prec ), intent( in )          :: c_in
        !
        select case(self%gridType)
            case(EDGE)
                self%x(:, 2:self%NdX(2)-1, :) = c_in
                self%x(:, :, 2:self%NdX(3)-1) = c_in
                self%y(2:self%NdY(1)-1, :, :) = c_in
                self%y(:, :, 2:self%NdY(3)-1) = c_in
                self%z(:, 2:self%NdZ(2)-1, :) = c_in
                self%z(2:self%NdZ(1)-1, :, :) = c_in
            case(FACE)
                self%x(2:self%NdX(1)-1, :, :) = c_in
                self%y(:, 2:self%NdY(2)-1, :) = c_in
                self%z(:, :, 2:self%NdZ(3)-1) = c_in
            case default
                write( *, * ) "ERROR:rVector3D_SG::setAllInteriorRVector3D_SG:"
                stop "           Invalid grid type. Exiting."
        end select
        !
    end subroutine setAllInteriorRVector3D_SG
    !**
    ! setOneBoundaryRVector3D_SG
    !*
    subroutine setOneBoundaryRVector3D_SG( self, bdry, c, int_only )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        character(*), intent( in )               :: bdry
        real( kind=prec ), intent( in )          :: c
        logical, optional, intent( in )          :: int_only
        !
        logical :: int_only_p
        !
        if( .NOT. present(int_only)) then
            int_only_p = .FALSE.
        else
            int_only_p = int_only
        endif
        !
        select case(self%gridType)
            case(EDGE)
                if(int_only_p) then
                  select case (bdry)
                      case("x1")
                          self%z(1, 2:self%NdZ(2)-1, :) = c
                          self%y(1, :, 2:self%NdY(3)-1) = c
                      case("x2")
                          self%z(self%NdZ(1), 2:self%NdZ(2)-1, :) = c
                          self%y(self%NdY(1), :, 2:self%NdY(3)-1) = c
                      case("y1")
                          self%z(2:self%NdZ(1)-1, 1, :) = c
                          self%x(:, 1, 2:self%NdX(3)-1) = c
                      case("y2")
                          self%z(2:self%NdZ(1)-1, self%NdZ(2), :) = c
                          self%x(:, self%NdX(2), 2:self%NdX(3)-1) = c
                      case("z1")
                          self%x(:, 2:self%NdX(2)-1, 1) = c
                          self%y(2:self%NdY(1)-1, :, 1) = c
                      case("z2")
                          self%x(:, 2:self%NdX(2)-1, self%NdX(3)) = c
                          self%y(2:self%NdY(1)-1, :, self%NdY(3)) = c
                      case("z1_x")
                          self%x(:, 2:self%NdX(2)-1, 1) = c
                      case("z2_x")
                          self%x(:, 2:self%NdX(2)-1, self%NdX(3)) = c
                      case("z1_y")
                          self%y(2:self%NdY(1)-1, :, 1) = c
                      case("z2_y")
                          self%y(2:self%NdY(1)-1, :, self%NdY(3)) = c
                  end select
                else
                  select case (bdry)
                      case("x1")
                          self%z(1, :, :) = c
                          self%y(1, :, :) = c
                      case("x2")
                          self%z(self%NdZ(1), :, :) = c
                          self%y(self%NdY(1), :, :) = c
                      case("y1")
                          self%z(:, 1, :) = c
                          self%x(:, 1, :) = c
                      case("y2")
                          self%z(:, self%NdZ(2), :) = c
                          self%x(:, self%NdX(2), :) = c
                      case("z1")
                          self%x(:, :, 1) = c
                          self%y(:, :, 1) = c
                      case("z2")
                          self%x(:, :, self%NdX(3)) = c
                          self%y(:, :, self%NdY(3)) = c
                      case("z1_x")
                          self%x(:, :, 1) = c
                      case("z2_x")
                          self%x(:, :, self%NdX(3)) = c
                      case("z1_y")
                          self%y(:, :, 1) = c
                      case("z2_y")
                          self%y(:, :, self%NdY(3)) = c
                  end select
                endif
                !
            case(FACE)
                select case (bdry)
                  case("x1")
                      self%x(1, :, :) = c
                  case("x2")
                      self%x(self%NdX(1), :, :) = c
                  case("y1")
                      self%y(:, 1, :) = c
                  case("y2")
                      self%y(:, self%NdY(2), :) = c
                  case("z1")
                      self%z(:, :, 1) = c
                  case("z2")
                      self%z(:, :, self%NdZ(3)) = c
                end select
            case default
                write( *, * ) "ERROR:rVector3D_SG::setOneBoundaryRVector3D_SG:"
                stop "           Invalid grid type. Exiting."
        end select
        !
    end subroutine setOneBoundaryRVector3D_SG

    !**
    ! intBdryIndicesRVector3D_SG
    !*
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
                  E = rVector3D_SG_t( grid, self%gridType )
            end select
        else
            write( *, * ) "ERROR:rVector3D_SG:intBdryIndicesRVector3D_SG:"
            stop "           Not allocated. Exiting."
        endif
        !
        select case(self%gridType)
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
    !**
    ! boundaryRVector3D_SG
    ! Returns a copy of this vector with all interiorRVector3D_SG elements ser to zero.
    !*
    function boundaryRVector3D_SG( self ) result( E )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        !
        class( rVector_t ), allocatable :: E
        !
        allocate( E, source = self )
        !
        call E%setAllInterior(R_ZERO)
       !
    end function boundaryRVector3D_SG
    !**
    ! interiorCVector3D_SG
    ! Returns a copy of this vector with all boundaryCVector3D_SG elements ser to zero.
    !*
    function interiorRVector3D_SG( self ) result( E )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        !
        class( rVector_t ), allocatable :: E
        !
        allocate( E, source = self )
        !
        call E%setAllboundary(R_ZERO)
        !
    end function interiorRVector3D_SG
    !
    !************************************************
    ! Data access
    !************************************************
    !
    !**
    ! lengthRVector3D_SG
    !
    !*
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
    !**
    ! getArrayRVector3D_SG
    !*
    subroutine getArrayRVector3D_SG( self, v )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in )         :: self
        real( kind=prec ), allocatable, intent( out ) :: v(:)
        allocate (v(self%length()))
        v = (/reshape(self%x, (/self%Nxyz(1), 1/)), &
                reshape(self%y, (/self%Nxyz(2), 1/)), &
                reshape(self%z, (/self%Nxyz(3), 1/))/)
        !
    end subroutine getArrayRVector3D_SG
    !**
    ! setArrayRVector3D_SG
    !*
    subroutine setArrayRVector3D_SG( self, v )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        real( kind=prec ), intent( in )          :: v(:)
        !
        integer :: i1, i2
        ! Ex
        i1 = 1; i2 = self%Nxyz(1)
        self%x = reshape(v(i1:i2), self%NdX)
        ! Ey
        i1 = i2 + 1; i2 = i2 + self%Nxyz(2)
        self%y = reshape(v(i1:i2), self%NdY)
        ! Ez
        i1 = i2 + 1; i2 = i2 + self%Nxyz(3)
        self%z = reshape(v(i1:i2), self%NdZ)
        !
    end subroutine setArrayRVector3D_SG
    !**
    ! setVecComponentsRVector3D_SG
    !
    !*
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
        real( kind=prec ), intent ( in )         :: c
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
                write( *, * ) "ERROR:rVector3D_SG::setVecComponentsRVector3D_SG:"
                stop "           Invalid xyz argument. Exiting."
        end select
        !
    end subroutine setVecComponentsRVector3D_SG
    !
    !************************************************
    ! Arithmetic/algebraic operations
    !************************************************
    !
    !**
    ! zerosRVector3D_SG
    !*
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
    !**
    ! add1RVector3D_SG
    !*
    function add1RVector3D_SG( lhs, rhs ) result( Eout )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: lhs
        class( rVector_t ), intent( in )      :: rhs
        class( rVector_t ), allocatable       :: Eout
        !
        if(lhs%isCompatible(rhs)) then
            !
            allocate(Eout, source = rVector3D_SG_t(lhs%grid, lhs%gridType))
            !
            select type( Eout )
                class is( rVector3D_SG_t )
                select type(rhs)
                  class is( rVector3D_SG_t )
                      Eout%x = lhs%x + rhs%x
                      Eout%y = lhs%y + rhs%y
                      Eout%z = lhs%z + rhs%z
                end select
            end select
        else
            write( *, * ) "ERROR:rVector3D_SG::add1RVector3D_SG"
            stop "    Incompatible inputs. Exiting."
        endif
        !
    end function add1RVector3D_SG
    !**
    ! sub1RVector3D_SG
    !*
    subroutine sub1RVector3D_SG( lhs, rhs )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: lhs
        class( rVector_t ), intent( in )         :: rhs
        !
        if( lhs%isCompatible( rhs ) ) then
            !
            select type( rhs )
                class is( rVector3D_SG_t )
                    lhs%x = lhs%x - rhs%x
                    lhs%y = lhs%y - rhs%y
                    lhs%z = lhs%z - rhs%z
            end select
            !
        else
            write( *, * ) "ERROR:rVector3D_SG::sub1RVector3D_SG"
            stop "    Incompatible inputs. Exiting."
        endif
        !
    end subroutine sub1RVector3D_SG
    !**
    ! mult1RVector3D_SG
    !*
    subroutine mult1RVector3D_SG( self, rhs )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( rVector_t ), intent( in )         :: rhs
        !
        if( self%isCompatible(rhs) ) then
            !
            select type(rhs)
                class is( rVector3D_SG_t )
                    self%x = self%x * rhs%x
                    self%y = self%y * rhs%y
                    self%z = self%z * rhs%z
            end select
            !
        else
            write( *, * ) "ERROR:rVector3D_SG::mult_1"
            stop "    Incompatible inputs. Exiting."
        endif
        !
    end subroutine mult1RVector3D_SG
    !**
    ! mult2RVector3D_SG
    !*
    subroutine mult2RVector3D_SG( self, c )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        real( kind=prec ), intent( in )          :: c
        !
        self%x = c * self%x
        self%y = c * self%y
        self%z = c * self%z
        !
    end subroutine mult2RVector3D_SG
    !**
    ! div1RVector3D_SG
    !*
    function div1RVector3D_SG( lhs, rhs ) result( Eout )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: lhs
        class( rVector_t ), intent( in )      :: rhs
        class( rVector_t ), allocatable       :: Eout
        !
        if(lhs%isCompatible(rhs)) then
            !
            allocate(Eout, source = rVector3D_SG_t(lhs%grid, lhs%gridType))
            !
            select type( Eout )
                class is( rVector3D_SG_t )
                  select type(rhs)
                      class is( rVector3D_SG_t )
                          Eout%x = lhs%x / rhs%x
                          Eout%y = lhs%y / rhs%y
                          Eout%z = lhs%z / rhs%z
                  end select
            end select
        else
            write( *, * ) "ERROR:rVector3D_SG::div1RVector3D_SG"
            stop "    Incompatible inputs. Exiting."
        endif
    end function div1RVector3D_SG
    !
    !    Arithmetic operations -- subroutine versions, first argument overwritten
    subroutine divS1RVector3D_SG( lhs, rhs )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: lhs
        class( rVector_t ), intent( in )         :: rhs
        !
        if(lhs%isCompatible(rhs)) then
            select type(rhs)
                class is( rVector3D_SG_t )
                  lhs%x = lhs%x / rhs%x
                  lhs%y = lhs%y / rhs%y
                  lhs%z = lhs%z / rhs%z
            end select
        else
             write( *, * ) "ERROR:rVector3D_SG::divS1RVector3D_SG"
             stop "    Incompatible inputs. Exiting."
        endif
        !
    end subroutine divS1RVector3D_SG
    !
    subroutine multS1RVector3D_SG( lhs, rhs )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: lhs
        class( rVector_t ), intent( in )         :: rhs
        !
        if(lhs%isCompatible(rhs)) then
            select type(rhs)
               class is( rVector3D_SG_t )
                 lhs%x = lhs%x * rhs%x
                 lhs%y = lhs%y * rhs%y
                 lhs%z = lhs%z * rhs%z
            end select
        else
            write( *, * ) "ERROR:rVector3D_SG::multS1RVector3D_SG"
            stop "    Incompatible inputs. Exiting."
        endif
    end subroutine multS1RVector3D_SG
    !
    subroutine multS2RVector3D_SG( lhs, r )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: lhs
        real(kind=prec) ,      intent( in )      :: r
        !
        lhs%x = lhs%x * r
        lhs%y = lhs%y * r
        lhs%z = lhs%z * r
        !
    end subroutine multS2RVector3D_SG

    !    end subroutine versions
    
    function dotProdRVector3D_SG( lhs, rhs ) result( r )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: lhs
        class( rVector_t ), intent( in )      :: rhs
        !
        real( kind=prec ) :: r
        !
        r = R_ZERO
        !
        if(( .NOT. lhs%is_allocated) .OR. ( .NOT. rhs%is_allocated)) then
            write( *, * ) "ERROR:rVector3D_SG::dotProdRVector3D_SG: "
            stop "           Input vectors not allocated. Exiting."
        endif
        !
        if(lhs%isCompatible(rhs)) then
            select type(rhs)
                class is( rVector3D_SG_t )
                    r = r + sum(lhs%x * rhs%x)
                    r = r + sum(lhs%y * rhs%y)
                    r = r + sum(lhs%z * rhs%z)
            end select
        else
            write( *, * ) "ERROR:rVector3D_SG:dotProdRVector3D_SG:"
            stop "           Incompatible input. Exiting."
        endif
        !
    end function dotProdRVector3D_SG
    !
    !****************************************************************************
    ! diagMultRVector3D_SG_rvector multiplies two vectors self, rhs stored as derived data
    ! type rvector pointwise; subroutine version
    ! mult can overwriteRVector3D_SG self or rhs
    function diagMultRVector3D_SG( self, rhs ) result ( diag_mult )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        class( rVector_t ), intent( in )      :: rhs
        class( rVector_t ), allocatable       :: diag_mult
        !
        if(self%isCompatible(rhs)) then
            !
            allocate(diag_mult, source = rVector3D_SG_t(self%grid, self%gridType))
            !
            select type(diag_mult)
                class is( rVector3D_SG_t )
                  select type(rhs)
                      class is( rVector3D_SG_t )
                          diag_mult%x = self%x * rhs%x
                          diag_mult%y = self%y * rhs%y
                          diag_mult%z = self%z * rhs%z
                end select
            end select
            else
                write( *, * ) "ERROR:rVector3D_SG::diagMultRVector3D_SG"
                stop "    Incompatible inputs. Exiting."
        endif
        !
    end function diagMultRVector3D_SG ! diagMultRVector3D_SG_rvector
    !**
    ! sumEdgesRVector3D_SG
    ! Sum all edges (or faces) that bound a cell, return as
    ! rScalar object.
    ! If interior only is true, only sum over interior edges (faces).
    !*
    function sumEdgesRVector3D_SG( self, interiorOnly ) result( cellObj )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        logical, optional, intent( in )       :: interiorOnly
        !
        class( rScalar_t ), allocatable       :: cellObj
        !
        integer :: x_xend, x_yend, x_zend
        integer :: y_xend, y_yend, y_zend
        integer :: z_xend, z_yend, z_zend
        !
        type( rVector3D_SG_t ) :: E_tmp
        !
        E_tmp = self
        !
        if(interiorOnly) then
            call E_tmp%setAllBoundary(R_ZERO)
        endif
        !
        allocate( cellObj, source = rScalar3D_SG_t(E_tmp%grid, CELL) )
        !
        select type(cellObj)
            class is(rScalar3D_SG_t)
                select case(E_tmp%gridType)
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
                      cellObj%v = E_tmp%x(:,1:x_yend-1,1:x_zend-1) + &
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
                      cellObj%v = E_tmp%x(1:x_xend-1,:,:) + E_tmp%x(2:x_xend,:,:) + &
                                E_tmp%y(:,1:y_yend-1,:) + E_tmp%y(:,2:y_yend,:) + &
                                E_tmp%z(:,:,1:z_zend-1) + E_tmp%z(:,:,2:z_zend)
                end select
        end select
        !
        deallocate( cellObj )
        !
    end function sumEdgesRVector3D_SG

    !**
    ! SumCellsRVector3D_SG
    ! Sum scalar object (cell type only) onto rVector object;
    ! default is to sum onto all bounding interior edges;
    ! boundary edges are presently zero;
    ! faces case is coded now.
    !*
    subroutine SumCellsRVector3D_SG( self, E_in, ptype )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( rScalar_t ), intent( in )         :: E_in
        character(*), intent( in ), optional     :: ptype
        !
        character(10) :: type
        integer :: xend, yend, zend
        integer :: v_xend, v_yend, v_zend
        integer :: ix, iy, iz
        !
        if(index(self%gridType, CELL) > 0) then
            write( *, * ) "ERROR:rVector3D_SG_t::SumCellsRVector3D_SG:"
            stop "           Only CELL type supported. Exiting."
        endif
        !
        if( .NOT. present(ptype)) then
            type = EDGE
        else
            type = ptype
        endif
        !
        select type(E_in)
            class is(rScalar3D_SG_t)
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
                write( *, * ) "ERROR:rVector3D_SG::SumCellsRVector3D_SG:"
                stop "           Incompatible input rScalar_t. Exiting."
        end select
        !
    end subroutine SumCellsRVector3D_SG
    !**
    ! InterpFunc
    ! Creates a Vector object containing weights needed for
    ! interpolation of xyz component of obj1 to location.
    !*
    function interpFuncRVector3D_SG( self, location, xyz ) result( E )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        real( kind=prec ), intent( in )       :: location(3)
        character, intent( in )               :: xyz
        !
        class( rVector_t ), allocatable :: E
        !
        real( kind=prec ), allocatable, dimension(:) :: xC, yC, zC
        integer :: ix, iy, iz, i
        real( kind=prec ) :: wx, wy, wz
        logical, dimension(:), allocatable :: tmp
        !
        select case(self%gridType)
        case(EDGE)
            !
            allocate(E, source = rVector3D_SG_t(self%grid, EDGE))
            !
            select case(xyz)
                case("x")
                    allocate(xC(size(self%grid%delX)))
                    allocate(yC(size(self%grid%dy + 1)))
                    allocate(zC(size(self%grid%dz + 1)))
                    !
                    xC = CumSum(self%grid%delX)
                    yC = CumSum([0._prec, self%grid%dy])
                    zC = CumSum([0._prec, self%grid%dz])
                case("y")
                    allocate(xC(size(self%grid%dx + 1)))
                    allocate(yC(size(self%grid%delY)))
                    allocate(zC(size(self%grid%dz)))
                    !
                    xC = CumSum([0._prec, self%grid%dx])
                    yC = CumSum([self%grid%delY])
                    zC = CumSum([0._prec, self%grid%dz])
                case("z")
                    allocate(xC(size(self%grid%dx + 1)))
                    allocate(yC(size(self%grid%dy + 1)))
                    allocate(zC(size(self%grid%delZ)))
                    !
                    xC = CumSum([0._prec, self%grid%dx])
                    yC = CumSum([0._prec, self%grid%dy])
                    zC = CumSum([self%grid%delZ])
            end select
        case(FACE)
            !
            allocate(E, source = rVector3D_SG_t(self%grid, FACE))
            !
            select case(xyz)
                case("x")
                    allocate(xC(size(self%grid%dx + 1)))
                    allocate(yC(size(self%grid%delY)))
                    allocate(zC(size(self%grid%delZ)))
                    !
                    xC = CumSum([0._prec, self%grid%dx])
                    yC = CumSum([self%grid%delY])
                    zC = CumSum([self%grid%delZ])
                case("y")
                    allocate(xC(size(self%grid%delX)))
                    allocate(yC(size(self%grid%dy + 1)))
                    allocate(zC(size(self%grid%delZ)))
                    !
                    xC = CumSum([self%grid%delX])
                    yC = CumSum([0._prec, self%grid%dy])
                    zC = CumSum([self%grid%delZ])
                case("z")
                    allocate(xC(size(self%grid%delX)))
                    allocate(yC(size(self%grid%delY)))
                    allocate(zC(size(self%grid%dz + 1)))
                    !
                    xC = CumSum([self%grid%delX])
                    yC = CumSum([self%grid%delY])
                    zC = CumSum([0._prec, self%grid%dz])
            end select
        end select
        !
        xC = xC - self%grid%ox
        yC = yC - self%grid%oy
        zC = zC - sum(self%grid%dz(1:self%grid%nzAir)) - self%grid%oz
        ! Find indicies of "minimal corner"
        ! Code below should be replaced by intrinsic
        ! "findloc" if it is avaiable.
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
        select type(E)
            class is( rVector3D_SG_t )
                select case(xyz)
                    case("x")
                        E%x(ix,iy,iz) = wx*wy*wz
                        E%x(ix+1,iy,iz) = (1-wx)*wy*wz
                        E%x(ix,iy+1,iz) = wx*(1-wy)*wz
                        E%x(ix,iy,iz+1) = wx*wy*(1-wz)
                        E%x(ix,iy+1,iz+1) = wx*(1-wy)*(1-wz)
                        E%x(ix+1,iy,iz+1) = (1-wx)*wy*(1-wz)
                        E%x(ix+1,iy+1,iz) = (1-wx)*(1-wy)*wz
                        E%x(ix+1,iy+1,iz+1) = (1-wx)*(1-wy)*(1-wz)
                        !
                    case("y")
                        E%y(ix,iy,iz) = wx*wy*wz
                        E%y(ix+1,iy,iz) = (1-wx)*wy*wz
                        E%y(ix,iy+1,iz) = wx*(1-wy)*wz
                        E%y(ix,iy,iz+1) = wx*wy*(1-wz)
                        E%y(ix,iy+1,iz+1) = wx*(1-wy)*(1-wz)
                        E%y(ix+1,iy,iz+1) = (1-wx)*wy*(1-wz)
                        E%y(ix+1,iy+1,iz) = (1-wx)*(1-wy)*wz
                        E%y(ix+1,iy+1,iz+1) = (1-wx)*(1-wy)*(1-wz)
                        !
                    case("z")
                        E%z(ix,iy,iz) = wx*wy*wz
                        E%z(ix+1,iy,iz) = (1-wx)*wy*wz
                        E%z(ix,iy+1,iz) = wx*(1-wy)*wz
                        E%z(ix,iy,iz+1) = wx*wy*(1-wz)
                        E%z(ix,iy+1,iz+1) = wx*(1-wy)*(1-wz)
                        E%z(ix+1,iy,iz+1) = (1-wx)*wy*(1-wz)
                        E%z(ix+1,iy+1,iz) = (1-wx)*(1-wy)*wz
                        E%z(ix+1,iy+1,iz+1) = (1-wx)*(1-wy)*(1-wz)
                end select
        end select
        !
    end function interpFuncRVector3D_SG
    !
    function isCompatibleRVector3D_SG( self, rhs ) result( status )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        class( rVector_t ), intent( in )      :: rhs
        logical :: status
        !
        status = .FALSE.
        !
        !if(same_type_as(self, rhs)) then
            select type(rhs)
                class is( rVector3D_SG_t )
                    if( self%nx == rhs%nx .AND.self%ny == rhs%ny .AND. self%nz == rhs%nz .AND.   &
                        self%gridType == rhs%gridType ) then
                        status = .TRUE.
                    end if
            end select
        !endif
    end function isCompatibleRVector3D_SG
    !
    subroutine copyFromRVector3D_SG( self, rhs )
        implicit none
        !
        class( rVector3D_SG_t ), intent( inout ) :: self
        class( rVector_t ), intent( in )         :: rhs
        !
        if( .NOT. rhs%is_allocated) then
            write( *, * ) "ERROR:rVector3D_SG::copyFromRVector3D_SG"
            stop "    Input not allocated. Exiting."
        endif
        !
        select type( rhs )
            class is( rVector3D_SG_t )
                !
                self%grid => rhs%grid
                self%gridType = rhs%gridType
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
                write( *, * ) "ERROR:rVector3D_SG_t::copyFromRVector3D_SG:"
                stop "            Incompatible input type. Exiting."
        end select
        !
    end subroutine copyFromRVector3D_SG
    !
    subroutine printRVector3D_SG( self, io_unit, title )
        implicit none
        !
        class( rVector3D_SG_t ), intent( in ) :: self
        integer, intent( in ), optional       :: io_unit
        character(*), intent( in ), optional  :: title
        !
        integer :: ix, iy, iz,fid
        !
        if( present( io_unit ) ) then
            fid = io_unit
        else
            fid = 0
        endif
        !
        !
        write( fid, * ) title
        write( fid, * ) self%nx, self%ny, self%nz
        write(fid, * ) "x-component",self%NdX
        do ix = 1, self%NdX(1)
             do iy = 1, self%NdX(2)
                do iz = 1, self%NdX(3)
                     if( self%x( ix, iy, iz ) /= 0 ) then
                        write(fid,*) ix,iy,iz, ":[", self%x( ix, iy, iz ), "]"
                     endif
                enddo
             enddo
        enddo
        !
        write(fid,*) "y-component",self%NdY
        do ix = 1, self%NdY(1)
             do iy = 1, self%NdY(2)
                do iz = 1, self%NdY(3)
                     if( self%y( ix, iy, iz ) /= 0 ) then
                        write(fid,*) ix,iy,iz, ":[", self%y( ix, iy, iz ), "]"
                     endif
                enddo
             enddo
        enddo
        !
        write(fid,*) "z-component",self%NdZ
        do ix = 1, self%NdZ(1)
             do iy = 1, self%NdZ(2)
                do iz = 1, self%NdZ(3)
                     if( self%z( ix, iy, iz ) /= 0 ) then
                        write(fid,*) ix,iy,iz, ":[", self%z( ix, iy, iz ), "]"
                     endif
                enddo
             enddo
        enddo
        !
    end subroutine printRVector3D_SG
    !
end module rVector3D_SG
