!**
! SUMMARY
!
! Standard cartesian grid vectors.
!*
module cVector3D_SG
    !
    use Constants
    use MatUtils
    use Grid
    use Grid3D_SG
    use rScalar
    use rScalar3D_SG
    use rVector
    use rVector3D_SG
    use cVector
    !
    type, extends( cVector_t ) :: cVector3D_SG_t
        !
        class( Grid3D_SG_t ), pointer :: grid
        !
        character( len=4 ) :: gridType
        !
        integer :: nx, ny, nz
        !
        integer, dimension(3) :: NdX, NdY, NdZ, Nxyz
        !
        complex( kind=prec ), allocatable, dimension(:, :, :) :: x, y, z
        !
     contains
        !**
        ! Initialization and finalization
        !*
        final :: cVector3D_SG_dtor
        !**
        ! Input/Output
        !*
        procedure, public :: read  => readCVector3D_SG
        procedure, public :: write => writeCVector3D_SG
        !**
        ! boundaryCVector3D_SG operations
        !*
        procedure, public :: setAllboundary => setAllboundaryCVector3D_SG
        procedure, public :: setOneboundary => setOneboundaryCVector3D_SG
        procedure, public :: setAllinterior => setAllinteriorCVector3D_SG
        procedure, public :: intBdryIndices => intBdryIndicesCVector3D_SG
        procedure, public :: boundary       => boundaryCVector3D_SG
        procedure, public :: interior       => interiorCVector3D_SG
        !**
        ! Data access
        !*
        procedure, public :: length => lengthCVector3D_SG
        procedure, public :: getArray => getArrayCVector3D_SG
        procedure, public :: setArray => setArrayCVector3D_SG
        procedure, public :: setVecComponentsArray => setVecComponentsArrayCVector3D_SG
        !**
        ! Arithmetic operations
        !*
        procedure, public :: zeros    => zerosCVector3D_SG
        procedure, public :: add1     => add1CVector3D_SG
        procedure, public :: sub1     => sub1CVector3D_SG
        procedure, public :: mult1    => mult1CVector3D_SG
        procedure, public, pass(self) :: mult2 => mult2CVector3D_SG
        procedure, public :: mult3 => mult3CVector3D_SG
        procedure, public :: mults1 => mults1CVector3D_SG
        procedure, public :: mults3 => mults3CVector3D_SG
        procedure, public :: div1 => div1CVector3D_SG
        procedure, public :: div2 => div2CVector3D_SG        
        procedure, public :: divs2 => divs2CVector3D_SG        
        procedure, public :: dotProd => dotProdCVector3D_SG
        procedure, public :: linCombS => linCombSCVector3D_SG        
        procedure, public :: scMultAddS => scMultAddSCVector3D_SG        
        !**
        ! Miscellaneous
        !*
        procedure, public :: isCompatible1 => isCompatible1CVector3D_SG
        procedure, public :: isCompatible2 => isCompatible2CVector3D_SG
        generic :: isCompatible => isCompatible1, isCompatible2
        !
        procedure, public :: copyFrom     => copyFromCVector3D_SG
        procedure, public :: interpFunc  => interpFuncCVector3D_SG
        !
        procedure, public :: print => printCVector3D_SG
    end type cVector3D_SG_t
    
    interface cVector3D_SG_t
        module procedure cVector3D_SG_ctor
    end interface cVector3D_SG_t
  
contains
    !
    function cVector3D_SG_ctor( igrid, gridType ) result ( self )
        implicit none
        !
        class( Grid3D_SG_t ), target, intent( in ) :: igrid
        character(*), intent( in )                 :: gridType
        !
        type( cVector3D_SG_t ) :: self
        !
        integer :: status
        !
        !write(*,*) "Constructor cVector3D_SG"
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
        self%gridType = trim(gridType)
        self%is_allocated = .FALSE.
        !
        if (self%gridType == EDGE) then
            allocate(self%x(self%nx, self%ny + 1, self%nz + 1), STAT = status)
            self%is_allocated = status.eq.0
            !
            allocate(self%y(self%nx + 1, self%ny, self%nz + 1), STAT = status)
            self%is_allocated = self%is_allocated .AND.(status.eq.0)
            !
            allocate(self%z(self%nx + 1, self%ny + 1, self%nz), STAT = status)
            self%is_allocated = self%is_allocated .AND.(status.eq.0)
            !
            self%NdX = (/self%nx, self%ny + 1, self%nz + 1/)
            self%NdY = (/self%nx + 1, self%ny, self%nz + 1/)
            self%NdZ = (/self%nx + 1, self%ny + 1, self%nz/)
            !
        else if (self%gridType == FACE) then
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
            write( *, * ) "ERROR:cVector3D_SG_ctor:"
            stop "           Only EDGE or FACE types allowed. Exiting."
        end if
        !
        if ( self%is_allocated ) then
            self%x = R_ZERO
            self%y = R_ZERO
            self%z = R_ZERO
        else
            write( *, * ) "ERROR:cVector3D_SG_ctor:"
            stop "           Unable to allocate vector. Exiting."
        !
        end if
        !
        self%Nxyz = (/product(self%NdX), &
                product(self%NdY), &
                product(self%NdZ)/)
        !
    end function cVector3D_SG_ctor
    
    subroutine cVector3D_SG_dtor( self )
        implicit none
        !
        type( cVector3D_SG_t ), intent( in out ) :: self
        !
        !write(*,*) "Destructor cVector3D_SG"
        !
        deallocate( self%x )
        deallocate( self%y )
        deallocate( self%z )
        !
        self%nx = 0
        self%ny = 0
        self%nz = 0
        
        self%gridType = ""
        self%is_allocated = .FALSE.
        !
    end subroutine cVector3D_SG_dtor
    
    !
    !************************************************
    ! Input/Output
    !************************************************
    !

    !**
    ! readCVector3D_SG
    !*
    subroutine readCVector3D_SG( self, fid )
        implicit none
        !
        class( cVector3D_SG_t ), intent( inout ) :: self
        integer, intent( in ) :: fid
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
            if ( (index( isbinary, "yes" ) > 0 .OR. index( isbinary, "YES" ) > 0 ) &
                  .AND. .NOT. binary ) then
                write( *, * ) "ERROR:cVector3D_SG_t::readCVector3D_SG: "
                write( *, * ) "           Unable to readCVector3D_SG vector from unformatted file. ", &
                       trim(fname), ". Exiting."
                stop
            else if ((index(isbinary, "no") > 0 .OR.index(isbinary, "NO") > 0) &
                  .AND.binary) then
                write( *, * ) "ERROR:cVector3D_SG_t::readCVector3D_SG: "
                write( *, * ) "           Unable to readCVector3D_SG vector from formatted file ", &
                       trim(fname), ". Exiting."
                stop
            end if
            !
            read( fid ) Nx, Ny, Nz, gridType
            write( *, * ) "Nx, Ny, Nz", Nx, Ny, Nz
            write( *, * ) "gridType", gridType
            !
            if ( .NOT. self%is_allocated ) then
                write( *, * ) "ERROR:cVector3D_SG_t::readCVector3D_SG: "
                write( *, * ) "           Vector must be allocated before readCVector3D_SGing from ", &
                       trim(fname), ". Exiting."
                stop
            else if (self%gridType.ne.gridType) then
                write( *, * ) "ERROR:cVector3D_SG_t::readCVector3D_SG: "
                write( *, * ) "           Vector must be of type ", gridType, &
                       &            "           before readCVector3D_SGing from ", trim (fname), ". Exiting."
                stop
            else if ((self%nx.ne.Nx).OR. &
                  (self%ny.ne.Ny).OR.(self%nz.ne.Nz)) then
                write( *, * ) "ERROR:cVector3D_SG_t::readCVector3D_SG: "
                write( *, * ) "           Wrong size of vector on input from ", &
                       trim (fname), ". Exiting."
                stop
            end if
            !
            read( fid ) self%x
            read( fid ) self%y
            read( fid ) self%z
            !
        else
            stop "readCVector3D_SG: unable to open file"
        endif
        !
    end subroutine readCVector3D_SG
    
    !**
    ! writeCVector3D_SG
    !*
    subroutine writeCVector3D_SG( self, fid )
        implicit none
        !
        class( cVector3D_SG_t ), intent( in ) :: self
        integer, intent( in )                 :: fid
        !
        integer :: Nx, Ny, Nz
        logical :: ok, hasname, binary
        character(80) :: fname, isbinary
        !
        binary = .TRUE.
        !
        if ( .NOT. self%is_allocated ) then
            write( *, * ) "ERROR:cVector3D_SG::writeCVector3D_SG: "
            stop "           Vector not allocated. Exiting."
        end if
        !
        inquire( fid, opened = ok, named = hasname, name = fname, unformatted = isbinary )
        !
        if( ok ) then
            !
            ! Check that the file is unformatted if binary, formatted if ascii.
            if ( ( index( isbinary, "yes" ) > 0 .OR. index( isbinary, "YES" ) > 0 ) &
                  .AND. .NOT. binary) then
                write( *, * ) "ERROR:cVector3D_SG::writeCVector3D_SG: "
                write( *, * ) "           Unable to writeCVector3D_SG vector to unformatted file. ", &
                       trim(fname), ". Exiting."
                stop
            else if ((index(isbinary,"no") > 0 .OR. index(isbinary,"NO") > 0) &
                  .AND. binary) then
                write( *, * ) "ERROR:cVector3D_SG::writeCVector3D_SG: "
                write( *, * ) "           Unable to writeCVector3D_SG vector to formatted file. ", &
                       trim(fname), ". Exiting."
                stop
            end if
            !
            write( fid ) self%nx, self%ny, self%nz, self%gridType
            write( fid ) self%x
            write( fid ) self%y
            write( fid ) self%z
        else
            stop "writeCVector3D_SG: unable to open file"
        endif
        !
    end subroutine writeCVector3D_SG
    !
    !************************************************
    ! boundaryCVector3D_SG operations
    !************************************************
    !
    !** 
    ! setAllboundaryCVector3D_SG
    !
    !*
    subroutine setAllboundaryCVector3D_SG( self, c_in )
        implicit none
        ! Arguments
        class( cVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in )       :: c_in
        !
        select case( self%gridType )
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
                write( *, * ) "ERROR:cVector3D_SG::setAllboundaryCVector3D_SG: "
                stop "           Invalid grid type. Exiting."
        end select
        !
    end subroutine setAllboundaryCVector3D_SG
    !** 
    ! setAllinteriorCVector3D_SG
    !
    !*
    subroutine setAllinteriorCVector3D_SG( self, c_in )
        implicit none
        ! Arguments
        class( cVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ) , intent( in )      :: c_in
        !
        select case(self%gridType)
            case(EDGE)
                self%x(:, 2:self%NdX(2)-1, 2:self%NdX(3)-1) = c_in
                self%y(2:self%NdY(1)-1, :, 2:self%NdY(3)-1) = c_in
                self%z(2:self%NdZ(1), 2:self%NdZ(2)-1, :) = c_in
            case(FACE)
                self%x(2:self%NdX(1)-1, :, :) = c_in
                self%y(:, 2:self%NdY(2)-1, :) = c_in
                self%z(:, :, 2:self%NdZ(3)-1) = c_in
                
            case default
                write( *, * ) "ERROR:cVector3D_SG::setAllinteriorCVector3D_SG: "
                stop "           Invalid grid type. Exiting."
        end select
        !
    end subroutine setAllinteriorCVector3D_SG
    !**
    ! setOneboundaryCVector3D_SG
    !
    !*
    subroutine setOneboundaryCVector3D_SG(self, bdry, c, int_only)
        implicit none
        !
        class( cVector3D_SG_t ), intent( inout ) :: self
        character(*), intent( in )               :: bdry
        complex( kind=prec ) , intent( in )      :: c
        logical, optional        , intent( in )  :: int_only
        !
        logical :: int_only_p
        
        if (.NOT.present(int_only)) then
            int_only_p = .FALSE.
        else
            int_only_p = int_only
        end if
        
        select case(self%gridType)
        case(EDGE)
            if (int_only_p) then
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
            end if
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
            write( *, * ) "ERROR:cVector3D_SG::setOneboundaryCVector3D_SG: "
            stop "           Invalid grid type. Exiting."
        end select
        !
    end subroutine setOneboundaryCVector3D_SG
    !**
    ! intBdryIndicesCVector3D_SG
    !*
    subroutine intBdryIndicesCVector3D_SG( self, ind_i, ind_b )
        implicit none
        !
        class( cVector3D_SG_t ), intent( in ) :: self
        integer, allocatable, intent( out ) :: ind_i(:), ind_b(:)
        !
        integer :: nVec(3), nVecT, nBdry, nb, ni, i
        complex( kind=prec ), allocatable   :: temp(:)
        type( cVector3D_SG_t ), allocatable :: E
        
        if( self%is_allocated ) then
            select type( grid => self%grid )
                class is(Grid3D_SG_t)
                  E = cVector3D_SG_t( grid, self%gridType )
            end select
        else
            write( *, * ) "ERROR:cVector3D_SG_t:intBdryIndicesCVector3D_SG"
            stop "           Not allocated. Exiting."
        end if
        
        select case(self%gridType)
            case( EDGE )
                !
                nVec(1) = size(E%x)
                nVec(2) = size(E%y)
                nVec(3) = size(E%z)
                nVecT = nVec(1) + nVec(2) + nVec(3)
                !
                allocate( temp( nVecT ) )
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
            case( FACE )
                nVec(1) = size(E%x)
                nVec(2) = size(E%y)
                nVec(3) = size(E%z)
                nVecT = nVec(1) + nVec(2) + nVec(3)
                !
                allocate( temp( nVecT ) )
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
        deallocate( E )
        !
        nBdry = 0
        do i = 1, nVecT
            nBdry = nBdry + nint(real(temp(i)))
        end do
        !
        if ( allocated( ind_i ) ) deallocate(ind_i)
        allocate( ind_i( nVecT - nBdry ) )
        !
        if( allocated( ind_b ) ) deallocate(ind_b)
        allocate( ind_b( nBdry ) )
        !
        nb = 0
        ni = 0
        do i = 1, nVecT
            if (nint(real(temp(i))).eq.1) then
                nb = nb + 1
                ind_b(nb) = i
            else
                ni = ni + 1
                ind_i(ni) = i
            end if
        end do
        !
        deallocate( temp )
        !
    end subroutine intBdryIndicesCVector3D_SG
    !**
    ! boundaryCVector3D_SG
    ! Returns a copy of this vector with all interiorCVector3D_SG elements ser to zero.
    !*
    function boundaryCVector3D_SG( self ) result( E )
        implicit none
        !
        class( cVector3D_SG_t ), intent( in ) :: self
        !
        class( cVector_t ), allocatable :: E
        !
        allocate( E, source = self )
        !
        call E%setAllInterior( C_ZERO )
        !
    end function boundaryCVector3D_SG
    !**
    ! interiorCVector3D_SG
    ! Returns a copy of this vector with all boundaryCVector3D_SG elements ser to zero.
    !*
    function interiorCVector3D_SG(self) result(E)
        implicit none
        !
        class( cVector3D_SG_t ), intent( in ) :: self
        !
        class( cVector_t ), allocatable :: E
        !
        allocate( E, source = self )
        !
        call E%setAllboundary( C_ZERO )
        !
    end function interiorCVector3D_SG
    !
    !************************************************
    ! Data access
    !************************************************
    !
    !**
    ! lengthCVector3D_SG
    !
    !*
    function lengthCVector3D_SG( self ) result( n )
        implicit none
        !
        class( cVector3D_SG_t ), intent( in ) :: self
        !
        integer :: n
        !
        n = self%Nxyz(1) + self%Nxyz(2) + self%Nxyz(3)
        !
    end function lengthCVector3D_SG
    !**
    ! getArrayCVector3D_SG
    !
    !*
    subroutine getArrayCVector3D_SG( self, v )
        implicit none
        !
        class( cVector3D_SG_t ), intent( in ) :: self
        complex( kind=prec ), allocatable, intent( out ) :: v(:)
        allocate ( v( self%length() ) )
        v = (/reshape(self%x, (/self%Nxyz(1), 1/)), &
                reshape(self%y, (/self%Nxyz(2), 1/)), &
                reshape(self%z, (/self%Nxyz(3), 1/))/)
        
    end subroutine getArrayCVector3D_SG
    !**
    ! setArrayCVector3D_SG
    !
    !*
    subroutine setArrayCVector3D_SG(self, v)
        implicit none
        ! Arguments
        class( cVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in )       :: v(:)
        ! Local variables
        integer :: i1, i2
        !
        ! Ex
        i1 = 1
        i2 = self%Nxyz(1)
        self%x = reshape(v(i1:i2), self%NdX)
        !
        ! Ey
        i1 = i2 + 1
        i2 = i2 + self%Nxyz(2)
        self%y = reshape(v(i1:i2), self%NdY)
        !
        ! Ez
        i1 = i2 + 1
        i2 = i2 + self%Nxyz(3)
        self%z = reshape(v(i1:i2), self%NdZ)
        !
    end subroutine setArrayCVector3D_SG
    !**
    ! setVecComponentsArrayCVector3D_SG
    !
    !*
    subroutine setVecComponentsArrayCVector3D_SG( self, xyz, &
            &                                 xmin, xstep, xmax, &
            &                                 ymin, ystep, ymax, &
            &                                 zmin, zstep, zmax, c )
        implicit none
        !
        class( cVector3D_SG_t ), intent( inout ) :: self
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
                if (xmin == 0) x1 = self%NdX(1)
                if (xmax <= 0) x2 = self%NdX(1) + xmax
                !
                if (ymin == 0) y1 = self%NdX(2)
                if (ymax <= 0) y2 = self%NdX(2) + ymax
                !
                if (zmin == 0) z1 = self%NdX(3)
                if (zmax <= 0) z2 = self%NdX(3) + zmax
                !
                self%x(x1:x2:xstep, y1:y2:ystep, z1:z2:zstep) = c
                !
            case("y")
                if (xmin == 0) x1 = self%NdY(1)
                if (xmax <= 0) x2 = self%NdY(1) + xmax
                !
                if (ymin == 0) y1 = self%NdY(2)
                if (ymax <= 0) y2 = self%NdY(2) + ymax
                !
                if (zmin == 0) z1 = self%NdY(3)
                if (zmax <= 0) z2 = self%NdY(3) + zmax
                !
                self%y(x1:x2:xstep, y1:y2:ystep, z1:z2:zstep) = c
                !
            case("z")
                if (xmin == 0) x1 = self%NdZ(1)
                if (xmax <= 0) x2 = self%NdZ(1) + xmax

                if (ymin == 0) y1 = self%NdZ(2)
                if (ymax <= 0) y2 = self%NdZ(2) + ymax
                !
                if (zmin == 0) z1 = self%NdZ(3)
                if (zmax <= 0) z2 = self%NdZ(3) + zmax
                !
                self%z(x1:x2:xstep, y1:y2:ystep, z1:z2:zstep) = c
                !
            case default
                write( *, * ) "ERROR:cVector3D_SG::setVecComponentsArrayCVector3D_SG: "
                stop "           Invalid xyz argument. Exiting."
        end select
        
    end subroutine setVecComponentsArrayCVector3D_SG
    !
    !************************************************
    ! Arithmetic/algebraic operations
    !************************************************
    !
    !**
    ! zerosCVector3D_SG
    !*
    subroutine zerosCVector3D_SG(self)
        implicit none
        !
        class( cVector3D_SG_t ), intent( inout ) :: self
        !
        self%x = R_ZERO
        self%y = R_ZERO
        self%z = R_ZERO
        
    end subroutine zerosCVector3D_SG
    !**
    ! add1CVector3D_SG
    !*
    function add1CVector3D_SG( lhs, rhs ) result( Eout )
        implicit none
        !
        class( cVector3D_SG_t ), intent( in ) :: lhs
        class( cVector_t ), intent( in )      :: rhs
        class( cVector_t ), allocatable       :: Eout
        !
        if (lhs%isCompatible(rhs)) then
            allocate(Eout, source = cVector3D_SG_t(lhs%grid, lhs%gridType))
         !
         select type(Eout)
            class is( cVector3D_SG_t )
                !
               select type(rhs)
                class is( cVector3D_SG_t )
                 !
                Eout%x = lhs%x + rhs%x
                    Eout%y = lhs%y + rhs%y
                    Eout%z = lhs%z + rhs%z
                end select
            end select
        else
            write( *, * ) "ERROR:cVector3D_SG::add1CVector3D_SG"
            stop "    Incompatible inputs. Exiting."
        end if
        !
    end function add1CVector3D_SG
    !**
    ! sub1CVector3D_SG
    !*
    function sub1CVector3D_SG(lhs, rhs) result(Eout)
        implicit none
        !
        class( cVector3D_SG_t ), intent( in ) :: lhs
        class( cVector_t ), intent( in ) :: rhs
        class( cVector_t ), allocatable :: Eout
        !
        if (lhs%isCompatible(rhs)) then
            !
            allocate( Eout, source = cVector3D_SG_t( lhs%grid, lhs%gridType ) )
            !
            select type(Eout)
                class is( cVector3D_SG_t )
                    select type(rhs)
                        class is( cVector3D_SG_t )
                            Eout%x = lhs%x - rhs%x
                            Eout%y = lhs%y - rhs%y
                            Eout%z = lhs%z - rhs%z
                    end select
            end select
            !
        else
            write( *, * ) "ERROR:cVector3D_SG::sub1CVector3D_SG"
            stop "    Incompatible inputs. Exiting."
        end if
    end function sub1CVector3D_SG
    !**
    ! mult1CVector3D_SG
    !*
    function mult1CVector3D_SG(lhs, rhs) result(Eout)
        implicit none
        !
        class( cVector3D_SG_t ), intent( in ) :: lhs
        class( cVector_t ), intent( in )      :: rhs
        class( cVector_t ), allocatable       :: Eout
        !
        if (lhs%isCompatible(rhs)) then
            !
            allocate(Eout, source = cVector3D_SG_t(lhs%grid, lhs%gridType))
            !
            select type(Eout)
                class is( cVector3D_SG_t )
                    select type(rhs)
                        class is( cVector3D_SG_t )
                            Eout%x = lhs%x * rhs%x
                            Eout%y = lhs%y * rhs%y
                            Eout%z = lhs%z * rhs%z
                    end select
            end select
            !
        else
            write( *, * ) "ERROR:cVector3D_SG::mult_1"
            stop "    Incompatible inputs. Exiting."
        end if
    end function mult1CVector3D_SG

    !**
    ! mult2CVector3D_SG
    !*
    function mult2CVector3D_SG( self, c ) result( Eout )
        implicit none
        !
        class( cVector3D_SG_t ), intent( in ) :: self
        complex( kind=prec ) , intent( in )   :: c
        !
        class( cVector_t ), allocatable       :: Eout
        !
        allocate( Eout, source = cVector3D_SG_t( self%grid, self%gridType ) )
        !
        select type(Eout)
            class is( cVector3D_SG_t )
                Eout%x = c * self%x
                Eout%y = c * self%y
                Eout%z = c * self%z
        end select
        !
    end function mult2CVector3D_SG
    !**
    ! mult3CVector3D_SG
    !*
    function mult3CVector3D_SG(lhs, rhs) result(Eout)
        implicit none
        !
        class( cVector3D_SG_t ), intent( in )  :: lhs
        class( rVector_t )        , intent( in ) :: rhs
        class( cVector_t ), allocatable        :: Eout
        !
        if (lhs%isCompatible(rhs)) then
            !
            allocate(Eout, source = cVector3D_SG_t(lhs%grid, lhs%gridType))
            !
            select type(Eout)
                class is( cVector3D_SG_t )
                    select type(rhs)
                        class is(rVector3D_SG_t)
                            Eout%x = lhs%x * rhs%x
                            Eout%y = lhs%y * rhs%y
                            Eout%z = lhs%z * rhs%z
                    end select
            end select
            !
        else
            write( *, * ) "ERROR:cVector3D_SG::mult3"
            stop "    Incompatible inputs. Exiting."
        end if
    end function mult3CVector3D_SG
    !**
    ! mults3CVector3D_SG
    !*
    !    subroutine version of mult3CVector3D_SG -- overwrites input lhs
    subroutine mults3CVector3D_SG(lhs, rhs)
        implicit none
        !
        class( cVector3D_SG_t ), intent( inout ) :: lhs
        class( rVector_t ), intent( in )         :: rhs
        !
        if (lhs%isCompatible(rhs)) then
            select type(rhs)
                class is(rVector3D_SG_t)
                    lhs%x = lhs%x * rhs%x
                    lhs%y = lhs%y * rhs%y
                    lhs%z = lhs%z * rhs%z
            end select
        else
            write( *, * ) "ERROR:cVector3D_SG::multS3"
            stop "    Incompatible inputs. Exiting."
        end if
        !
    end subroutine mults3CVector3D_SG
    !**
    ! mults1CVector3D_SG
    ! Subroutine version of mult2CVector3D_SG
    !*
    subroutine mults1CVector3D_SG(self, c)
        implicit none
        !
        class( cVector3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ) , intent( in )      :: c
        !
        self%x = c * self%x
        self%y = c * self%y
        self%z = c * self%z
        !
    end subroutine mults1CVector3D_SG
    !**
    ! linCombSCVector3D_SG
    !*
    subroutine linCombSCVector3D_SG(lhs, rhs, c1, c2)
        implicit none
        !
        class( cVector3D_SG_t ), intent( inout ) :: lhs
        class( cVector_t )        , intent( in ) :: rhs
        complex( kind=prec ), intent( in )       :: c1, c2
        !
        if (lhs%isCompatible(rhs)) then
            select type(rhs)
                class is( cVector3D_SG_t )
                    lhs%x = c1*lhs%x + c2*rhs%x
                    lhs%y = c1*lhs%y + c2*rhs%y
                    lhs%z = c1*lhs%z + c2*rhs%z
            end select
        else
            write( *, * ) "ERROR:cVector3D_SG::linCombS"
            stop "    Incompatible inputs. Exiting."
        end if
    end subroutine linCombSCVector3D_SG
    !*
    ! scMultAddSCScalar3D_SG
    !*
    subroutine scMultAddSCVector3D_SG( lhs, rhs, c )
        implicit none
        !
        class( cVector3D_SG_t ), intent( in ) :: lhs
        class( cVector_t ), intent( inout )   :: rhs
        complex( kind=prec ), intent( in )    :: c
        !
        if (lhs%isCompatible(rhs)) then
            select type(rhs)
                class is( cVector3D_SG_t ) 
                    rhs%x = rhs%x + c*lhs%x
                    rhs%y = rhs%y + c*lhs%y
                    rhs%z = rhs%z + c*lhs%z
            end select
        else
            write( *, * ) "ERROR:cVector3D_SG::scMultAdd"
            stop "    Incompatible inputs. Exiting."
        end if
    end subroutine scMultAddSCVector3D_SG
    !**
    ! div1CVector3D_SG
    !*
    function div1CVector3D_SG(lhs, rhs) result(Eout)
        implicit none
        !
        class( cVector3D_SG_t ), intent( in )    :: lhs
        class( cVector_t )        , intent( in ) :: rhs
        class( cVector_t ), allocatable          :: Eout
        !
        if (lhs%isCompatible(rhs)) then
            !
            allocate( Eout, source = cVector3D_SG_t( lhs%grid, lhs%gridType ) )
            !
            select type(Eout)
            class is( cVector3D_SG_t )
                select type(rhs)
                class is( cVector3D_SG_t )
                    Eout%x = lhs%x / rhs%x
                    Eout%y = lhs%y / rhs%y
                    Eout%z = lhs%z / rhs%z
                end select
            end select
            !
        else
            write( *, * ) "ERROR:cVector3D_SG::div1CVector3D_SG"
            stop "    Incompatible inputs. Exiting."
        end if
        !
    end function div1CVector3D_SG
    !**
    ! div2CVector3D_SG
    !*
    function div2CVector3D_SG(lhs, rhs) result(Eout)
        implicit none
        !
        class( cVector3D_SG_t ), intent( in )  :: lhs
        class( rVector_t )        , intent( in ) :: rhs
        class( cVector_t ), allocatable        :: Eout
        !
        if (lhs%isCompatible(rhs)) then
            !
            allocate(Eout, source = cVector3D_SG_t(lhs%grid, lhs%gridType))
            !
            select type(Eout)
                class is( cVector3D_SG_t )
                    select type(rhs)
                        class is(rVector3D_SG_t)
                            Eout%x = lhs%x / rhs%x
                            Eout%y = lhs%y / rhs%y
                            Eout%z = lhs%z / rhs%z
                end select
            end select
            !
        else
            write( *, * ) "ERROR:cVector3D_SG::div2CVector3D_SG"
            stop "    Incompatible inputs. Exiting."
        end if
    end function div2CVector3D_SG
    !**
    ! divs2CVector3D_SG
    !*
    !     Subroutine version of div2CVector3D_SG -- overwrites lhs
    subroutine divs2CVector3D_SG(lhs, rhs)
        implicit none
        !
        class( cVector3D_SG_t ), intent( inout ) :: lhs
        class( rVector_t ), intent( in ) :: rhs
        !
        if (lhs%isCompatible(rhs)) then
           select type(rhs)
               class is(rVector3D_SG_t)
                lhs%x = lhs%x / rhs%x
                lhs%y = lhs%y / rhs%y
                lhs%z = lhs%z / rhs%z
           end select
        else
            write( *, * ) "ERROR:cVector3D_SG::div2CVector3D_SG"
            stop "    Incompatible inputs. Exiting."
        end if
        !
    end subroutine divs2CVector3D_SG
    !
    function dotProdCVector3D_SG(lhs, rhs) result( c )
        implicit none
        !
        class( cVector3D_SG_t ), intent( in ) :: lhs
        class( cVector_t ), intent( in )      :: rhs
        !
        complex( kind=prec ) :: c
        !
        c = C_ZERO
        !
        if ((.NOT.lhs%is_allocated).OR.(.NOT.rhs%is_allocated)) then
            write( *, * ) "ERROR:cVector3D_SG::dotProdCVector3D_SG: "
            stop "           Input vectors not allocated. Exiting."
        end if
        !
        if (lhs%isCompatible(rhs)) then
            select type(rhs)
                class is( cVector3D_SG_t )
                    c = c + sum(conjg(lhs%x) * rhs%x)
                    c = c + sum(conjg(lhs%y) * rhs%y)
                    c = c + sum(conjg(lhs%z) * rhs%z)
            end select
        else
            write( *, * ) "ERROR:cVector3D_SG:dotProdCVector3D_SG:"
            stop         "           Incompatible input. Exiting"
        end if
    end function dotProdCVector3D_SG
    !**
    ! InterpFunc
    ! Creates a Vector object containing weights needed for
    ! interpolation of xyz component of obj1 to location.
    !*
    subroutine interpFuncCVector3D_SG( self, location, xyz, E )
        implicit none
        !
        class( cVector3D_SG_t ), intent( in )          :: self
        real( kind=prec ), intent( in )                :: location(3)
        character, intent( in )                        :: xyz
        class( cVector_t ), intent( out ), allocatable :: E
        !
        real( kind=prec ), allocatable, dimension(:) :: xC, yC, zC
        integer :: ix, iy, iz, i
        real( kind=prec ) :: wx, wy, wz
        logical, dimension(:), allocatable :: tmp
        !
        select case(self%gridType)
        case(EDGE)
            allocate(E, source = cVector3D_SG_t(self%grid, EDGE))
            select case(xyz)
            case("x")
                allocate(xC(size(self%grid%delX)))
                allocate(yC(size(self%grid%dy + 1)))
                allocate(zC(size(self%grid%dz + 1)))
                
                xC = CumSum(self%grid%delX)
                yC = CumSum([0._prec, self%grid%dy])
                zC = CumSum([0._prec, self%grid%dz])
            case("y")
                allocate(xC(size(self%grid%dx + 1)))
                allocate(yC(size(self%grid%delY)))
                allocate(zC(size(self%grid%dz)))
                
                xC = CumSum([0._prec, self%grid%dx])
                yC = CumSum([self%grid%delY])
                zC = CumSum([0._prec, self%grid%dz])
            case("z")
                allocate(xC(size(self%grid%dx + 1)))
                allocate(yC(size(self%grid%dy + 1)))
                allocate(zC(size(self%grid%delZ)))
                
                xC = CumSum([0._prec, self%grid%dx])
                yC = CumSum([0._prec, self%grid%dy])
                zC = CumSum([self%grid%delZ])
            end select
                      
        case(FACE)
            allocate(E, source = cVector3D_SG_t(self%grid, FACE))
            select case(xyz)
            case("x")
                allocate(xC(size(self%grid%dx + 1)))
                allocate(yC(size(self%grid%delY)))
                allocate(zC(size(self%grid%delZ)))
                
                xC = CumSum([0._prec, self%grid%dx])
                yC = CumSum([self%grid%delY])
                zC = CumSum([self%grid%delZ])
            case("y")
                allocate(xC(size(self%grid%delX)))
                allocate(yC(size(self%grid%dy + 1)))
                allocate(zC(size(self%grid%delZ)))
                
                xC = CumSum([self%grid%delX])
                yC = CumSum([0._prec, self%grid%dy])
                zC = CumSum([self%grid%delZ])
            case("z")
                allocate(xC(size(self%grid%delX)))
                allocate(yC(size(self%grid%delY)))
                allocate(zC(size(self%grid%dz + 1)))
                
                xC = CumSum([self%grid%delX])
                yC = CumSum([self%grid%delY])
                zC = CumSum([0._prec, self%grid%dz])
            end select
        end select

        xC = xC + self%grid%ox
        yC = yC + self%grid%oy
        zC = zC - sum(self%grid%dz(1:self%grid%nzAir)) - self%grid%oz
        !
        tmp = location(1) > xC
        do i = size(tmp), 1, -1 
            if (tmp(i)) then
                ix = i
                exit
            end if
        end do
        tmp = location(2) > yC
        do i = size(tmp), 1, -1 
            if (tmp(i)) then
                iy = i
                exit
            end if
        end do
        tmp = location(3) > zC
        do i = size(tmp), 1, -1 
            if (tmp(i)) then
                iz = i
                exit
            end if
        end do
        !
        deallocate( tmp )
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
        select type(E)
        class is( cVector3D_SG_t )
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
                
            case("y")
                E%y(ix,iy,iz) = wx*wy*wz
                E%y(ix+1,iy,iz) = (1-wx)*wy*wz
                E%y(ix,iy+1,iz) = wx*(1-wy)*wz
                E%y(ix,iy,iz+1) = wx*wy*(1-wz)
                E%y(ix,iy+1,iz+1) = wx*(1-wy)*(1-wz)
                E%y(ix+1,iy,iz+1) = (1-wx)*wy*(1-wz)
                E%y(ix+1,iy+1,iz) = (1-wx)*(1-wy)*wz
                E%y(ix+1,iy+1,iz+1) = (1-wx)*(1-wy)*(1-wz)
                
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
    end subroutine interpFuncCVector3D_SG

    subroutine copyFromCVector3D_SG(self, rhs)
        implicit none
        !
        class( cVector3D_SG_t ), intent( inout ) :: self
        class( cVector_t )        , intent( in ) :: rhs

        if (.NOT.rhs%is_allocated) then
            write( *, * ) "ERROR:cVector3D_SG::copyFromCVector3D_SG"
            stop "    Inputs not allocated. Exiting."
        end if
        !
        select type( rhs )
            class is( cVector3D_SG_t )
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
                self%is_allocated = rhs%is_allocated
                !
            class default
                write( *, * ) "ERROR:cVector3D_SG_t::copyFromCVector3D_SG:"
                stop "            Incompatible input type. Exiting."
        end select
        !
    end subroutine copyFromCVector3D_SG
    !
    function isCompatible1CVector3D_SG( self, rhs ) result( status )
        implicit none
        !
        class( cVector3D_SG_t ), intent( in ) :: self
        class( cVector_t ), intent( in )      :: rhs
        logical :: status
        !
        status = .FALSE.
        !
        select type(rhs)
            class is( cVector3D_SG_t )
                if( self%nx == rhs%nx .AND.self%ny == rhs%ny .AND. self%nz == rhs%nz .AND.   &
                    self%gridType == rhs%gridType ) then
                    status = .TRUE.
                end if
        end select
        !
    end function isCompatible1CVector3D_SG
    
    function isCompatible2CVector3D_SG( self, rhs) result( status )
        implicit none
        !
        class( cVector3D_SG_t ), intent( in ) :: self
        class( rVector_t ), intent( in )      :: rhs
        logical :: status
        !
        status = .FALSE.
        !
        select type(rhs)
            class is(rVector3D_SG_t)
                if( self%nx == rhs%nx .AND.self%ny == rhs%ny .AND. self%nz == rhs%nz .AND.   &
                    self%gridType == rhs%gridType ) then
                    status = .TRUE.
                end if
        end select
        !
    end function isCompatible2CVector3D_SG
    !
    subroutine printCVector3D_SG( self, io_unit, title, append )
        implicit none
        !
        class( cVector3D_SG_t ), intent( in ) :: self
        integer, intent( in ), optional       :: io_unit
        character(*), intent( in ), optional  :: title
        logical, intent( in ), optional       :: append
        !
        integer :: ix, iy, iz,fid
        !
        if( present( io_unit ) ) then
           fid = io_unit
        else
           fid = 0    !    usually this will work to write to standard output
        endif
        !
        if(present(title)) then
         write( fid, * ) title
        end if
        !
        !if ( present( append ) .AND. append ) open( fid, position="append", action="write" )
        !
        write( fid, * ) self%nx, self%ny, self%nz
        !
        !
        write( fid, * ) "x-component",self%NdX
        do iz = 1, self%NdX(3)
            do iy = 1, self%NdX(2)
                do ix = 1, self%NdX(1)
                    if( self%x( ix, iy, iz ) /= 0 ) then
                      write(fid,*) ix,iy,iz, ":[", self%x( ix, iy, iz ), "]"
                    endif
                enddo
            enddo
        enddo
        !
        write(fid,*) "y-component",self%NdY
        do iz = 1, self%NdY(3)
            do iy = 1, self%NdY(2)
                do ix = 1, self%NdY(1)
                    if( self%y( ix, iy, iz ) /= 0 ) then
                      write(fid,*) ix,iy,iz, ":[", self%y( ix, iy, iz ), "]"
                    endif
                enddo
            enddo
        enddo
        !
        write(fid,*) "z-component",self%NdZ
        do iz = 1, self%NdZ(3)
            do iy = 1, self%NdZ(2)
                do ix = 1, self%NdZ(1)
                    if( self%z( ix, iy, iz ) /= 0 ) then
                      write(fid,*) ix,iy,iz, ":[", self%z( ix, iy, iz ), "]"
                    endif
                enddo
            enddo
        enddo
        !
    end subroutine printCVector3D_SG
    !
end module cVector3D_SG
