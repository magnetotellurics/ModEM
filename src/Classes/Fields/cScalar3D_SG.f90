!
module cScalar3D_SG
    !
    use Constants
    use Grid3D_SG
    use rScalar3D_SG
    !
    type, extends( Scalar_t ) :: cScalar3D_SG_t
        !
        complex( kind=prec ), allocatable :: v(:, :, :)
        !
    contains
        !
        ! Destructor
        final :: cScalar3D_SG_dtor
        !
        ! I/O operation
        procedure, public :: read  => readCScalar3D_SG
        procedure, public :: write => writeCScalar3D_SG
        !
        ! Boundary operations
        procedure, public :: setAllBoundary => setAllBoundaryCScalar3D_SG
        procedure, public :: setOneBoundary => setOneBoundaryCScalar3D_SG
        procedure, public :: setAllInterior => setAllInteriorCScalar3D_SG
        procedure, public :: intBdryIndices => intBdryIndicesCScalar3D_SG
        !
        ! Dimensioning operations
        procedure, public :: length => lengthCScalar3D_SG
        !
        procedure, public :: getRealArray    => getRealArrayCScalar3D_SG
        procedure, public :: getComplexArray => getComplexArrayCScalar3D_SG
        !
        procedure, public :: setRealArray    => setRealArrayCScalar3D_SG
        procedure, public :: setComplexArray => setComplexArrayCScalar3D_SG
        !
        procedure, public :: setVecComponents => setVecComponentsCScalar3D_SG
        !
        ! Arithmetic/algebraic operations
        procedure, public :: zeros => zerosCScalar3D_SG
        procedure, public :: add   => addCScalar3D_SG
        procedure, public :: sub   => subCScalar3D_SG
        !
        procedure, public :: multByField => multByFieldCScalar3D_SG
        procedure, public :: multByValue => multByValueCScalar3D_SG
        !
        procedure, public :: divByField => divByFieldCScalar3D_SG
        procedure, public :: divByValue => divByValueCScalar3D_SG
        !
        procedure, public :: dotProd => dotProdCScalar3D_SG
        !
        ! Miscellaneous
        procedure, public :: linCombS   => linCombSCScalar3D_SG
        procedure, public :: scMultAddS => scMultAddSCScalar3D_SG
        !
        procedure, public :: copyFrom => copyFromCScalar3D_SG
        !
        procedure, public :: print => printCScalar3D_SG
        !
    end type cScalar3D_SG_t
    !
    interface cScalar3D_SG_t
        module procedure cScalar3D_SG_ctor
    end interface cScalar3D_SG_t
    !
contains
    !
    function cScalar3D_SG_ctor( grid, grid_type ) result ( self )
        implicit none
        !
        class( Grid3D_SG_t ), target, intent( in ) :: grid
        character( len=4 ), intent( in )           :: grid_type
        !
        type( cScalar3D_SG_t ) :: self
        !
        !
        integer :: nx, ny, nz, nzAir, nz_earth, status
        !
        !write( *, * ) "Constructor cScalar3D_SG"
        !
        call self%init()
        !
        self%grid => grid
        self%grid_type = grid_type
        !
        ! Grid dimensions
        call grid%GetDimensions(nx, ny, nz, nzAir)
        nz_earth = nz - nzAir
        !
        self%nx = nx
        self%ny = ny
        self%nz = nz
        !
        ! allocate memory for x,y,z ;
        ! self%allocated will be true if all allocations succeed
        self%is_allocated = .TRUE.
        !
        if( grid_type == CORNER) then
             allocate(self%v(nx + 1, ny + 1, nz + 1), STAT = status)    
             self%NdV = (/self%nx + 1, self%ny + 1, self%nz + 1/)
             
        else if( grid_type == CENTER) then             
             allocate(self%v(nx, ny, nz), STAT = status) 
             self%NdV = (/self%nx, self%ny, self%nz/)
             
        else if( grid_type == CELL_EARTH) then
             self%nz = nz_earth
             allocate(self%v(nx, ny, nz_earth), STAT = status)
             self%NdV = (/nx, ny, nz_earth/)
             
        else
             write( *, * ) "Error: cScalar3D_SG_ctor > unrecognized grid type: [", grid_type, "]"
             stop
        end if
        !
        self%is_allocated = self%is_allocated.AND.(status .EQ. 0)
        if( self%is_allocated) then
             self%v = R_ZERO
        else
             stop "Error: cScalar3D_SG_ctor > Unable to allocate cScalar - invalid grid supplied"
        end if
        !
        self%Nxyz = product( self%NdV )
        !
    end function cScalar3D_SG_ctor
    !
    subroutine cScalar3D_SG_dtor( self )
        implicit none
        !
        type( cScalar3D_SG_t ), intent( in out ) :: self
        !
        !write( *, * ) "Destructor cScalar3D_SG"
        !
        deallocate( self%v )
        !
        self%nx = 0
        self%ny = 0
        self%nz = 0
        !
        self%grid_type = ""
        self%is_allocated = .FALSE.
        !
    end subroutine cScalar3D_SG_dtor
    !
    subroutine readCScalar3D_SG( self, funit, ftype )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( inout ) :: self
        integer, intent( in )                    :: funit
        character(:), allocatable, intent( in ), optional :: ftype
        !
        integer :: Nx, Ny, Nz
        character(4) :: grid_type
        integer :: i, j, k, k1, k2, istat
        complex( kind=prec ), allocatable :: temp(:)
        logical :: ok, hasname, binary
        character(:), allocatable :: fname, isbinary
        !
        if( .NOT. present (ftype)) then
             binary = .FALSE.
        else if( index (ftype, "b") > 0) then
             binary = .TRUE.
        else
             binary = .FALSE.
        end if
        !
        inquire( funit, opened = ok, named = hasname, name = fname, unformatted = isbinary )
        !
        if( ok ) then
            !
            ! check that the file is unformatted if binary, formatted if ascii
            if( (index(isbinary, "yes") > 0 .or. index(isbinary, "YES") > 0) &
                     .AND. .NOT.binary) then             
                 write( *, * ) "Error: cScalar3D_SG_t::readCScalar3D_SG: "
                 write( *, * ) "            Unable to read scalar from unformatted file ", &
                            trim(fname), ".Exiting."
                 stop
            else if( (index(isbinary, "no") > 0 .or. index(isbinary, "NO") > 0) &
                     .AND.binary) then
                 write( *, * ) "Error: cScalar3D_SG_t::readCScalar3D_SG: "
                 write( *, * ) "            Unable to read scalar from formatted file ", &
                            trim(fname), ". Exiting."
                 stop
            end if
            !
            if( binary) then
                 ! read binary from unformatted files
                 read(funit) self%Nx, self%Ny, self%Nz, grid_type
                 read(funit) self%v
            end if
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
                 if( (k1 < 0) .or. (k2 > Nz)) then
                        write( *, * ) "Error: cScalar3D_SG::readCScalar3D_SG: "
                        write( *, * ) "      While reading the ", i, "th block. Exiting."
                        stop
                 else if( k1 > k2) then
                        write( *, * ) "Warning: cScalar3D_SG::readCScalar3D_SG: "
                        write( *, * ) "                Block ", i, " will be ignored."
                 end if
                 !
                 do j = Nx, 1, -1
                        read(funit, *, iostat = istat) temp
                        
                        if( istat /= 0) then
                             write( *, * ) "Error: cScalar3D_SG::readCScalar3D_SG: "
                             write( *, * ) "            While reading the ", j, "th row in ", i,"th block. Exiting."
                             stop
                        end if
                        
                        do k = k1, k2
                             self%v(j, :, k) = temp
                        end do
                 end do
                 !
                 if( k == Nz) exit
                 !
                 i = i + 1
                 !
            end do
            !
            deallocate( temp )
            !
        else
            stop "Error: readCScalar3D_SG: unable to open file"
        endif
        !
    end subroutine readCScalar3D_SG
    !
    subroutine writeCScalar3D_SG( self, funit, ftype )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( in ) :: self
        integer, intent( in )                 :: funit
        character(:), allocatable, intent( in ), optional :: ftype
        !
        integer :: Nx, Ny, Nz
        integer :: i, j, k, k1, k2, istat
        complex( kind=prec ), allocatable, dimension(:, :) :: temp
        logical :: ok, hasname, binary
        character(80) :: fname, isbinary
        !
        if( .NOT. self%is_allocated) then
             stop "Error: cScalar3D_SG::writeCScalar3D_SG > Not allocated"
        end if
        !
        if( .NOT.present(ftype)) then
             binary = .FALSE.
        else if( index(ftype, "b") > 0) then
             binary = .TRUE.
        else
             binary = .FALSE.
        end if
        !
        inquire( funit, opened = ok, named = hasname, name = fname, unformatted = isbinary )
        !
        if( ok ) then
            !
            if( (index(isbinary, "yes") > 0.or.index(isbinary, "YES") > 0) &
                     .AND..NOT.binary) then             
                 write( *, * ) "Error: cScalar3D_SG::writeCScalar3D_SG: "
                 write( *, * ) "            Unable to write vector to unformatted file ", &
                            trim(fname), ". Exiting."
                 !
                 stop
            else if( (index(isbinary,"no") > 0.or.index(isbinary,"NO") > 0) &
                     .AND.binary) then
                 write( *, * ) "Error: cScalar3D_SG::writeCScalar3D_SG: "
                 write( *, * ) " Unable to write vector to formatted file ", &
                            trim(fname), ". Exiting."
                 !
                 stop
            end if
            !
            if( binary) then
                 write(funit) self%nx, self%ny, self%nz, self%grid_type
                 write(funit) self%v             
                 return
            end if
            !
            !**
            ! ASCII format
            !*
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
                        temp = abs(self%v(:, :, k + 1) - self%v(:, :, k))
                        if( maxval(real(temp)) > TOL6) then
                             k2 = k
                             exit
                        end if
                 end do
                 !
                 write(funit, "(2i5)", iostat = istat) k1, k2
                 !
                 if( istat /= 0) then
                        write( *, * ) "Error: cScalar3D_SG::writeCScalar3D_SG: "
                        write( *, * ) "            Failed while writing to file. Exiting."
                        
                        stop
                 end if
                 !
                 temp = self%v(:, :, k1)
                 !
                 do i = Nx, 1, -1
                        do j = 1, Ny
                             write(funit, "(es13.5)", iostat = istat, &
                                        advance = "no") self%v(i, j, k1)
                        end do
                        write(funit, *)
                 end do
                 !
                 k1 = k2 + 1
                 !
                 if( k1 > Nz) exit
            end do
            !
            deallocate( temp )
            !
        else
            stop "Error: readRVector3D_SG: unable to open file"
        endif
        !
    end subroutine writeCScalar3D_SG
    !
    subroutine setAllBoundaryCScalar3D_SG( self, cvalue )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in )       :: cvalue
        !
        select case( self%grid_type )
            case (CORNER) 
                 self%v((/1, self%NdV(1)/), :, :) = cvalue
                 self%v(:, (/1, self%NdV(2)/), :) = cvalue
                 self%v(:, :, (/1, self%NdV(3)/)) = cvalue
                 !
            case default
                 stop "Error: setAllBoundaryCScalar3D_SG > Grid type not recognized. Exiting."
        end select
        !
    end subroutine setAllBoundaryCScalar3D_SG
    !
    subroutine setOneBoundaryCScalar3D_SG( self, bdry, cvalue, int_only )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( inout ) :: self
        character(:), allocatable, intent( in )  :: bdry
        complex( kind=prec ), intent( in )       :: cvalue
        logical, intent( in ), optional          :: int_only
        !
        logical :: int_only_p
        !
        if( .NOT. present (int_only)) then
             int_only_p = .FALSE.
        else 
             int_only_p = int_only
        end if
        !
        select case( self%grid_type )
        case (CORNER)
             if( int_only_p) then
                select case (bdry)
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
             end if
             !
        case(FACE)
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
             stop "Error: setOneBoundaryCScalar3D_SG > Invalid grid type"
        end select
        !
    end subroutine setOneBoundaryCScalar3D_SG
    !
    subroutine setAllInteriorCScalar3D_SG( self, cvalue )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in )       :: cvalue
        !
        stop "Error: setAllInteriorCScalar3D_SG to be implemented!"
        !
    end subroutine setAllInteriorCScalar3D_SG
    !
    subroutine intBdryIndicesCScalar3D_SG( self, ind_i, ind_b )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( in ) :: self
        integer, allocatable, intent( out )   :: ind_i(:), ind_b(:)
        !
        integer :: nVec(3), nVecT, nBdry, nb, ni, i
        complex( kind=prec ), allocatable :: temp(:)
        type( cScalar3D_SG_t ) :: phi
        !
        !
        if( self%is_allocated) then
            select type(grid => self%grid)
                class is( Grid3D_SG_t )
                    phi = cScalar3D_SG_t( grid, self%grid_type )
                class default
                    stop "Error: intBdryIndicesCScalar3D_SG: undefined grid"
            end select
        else
            stop "Error: intBdryIndicesCScalar3D_SG > Not allocated."
        end if
        !
        select case( self%grid_type )
        case(CORNER)
            nVecT = size(phi%v)
            !
            allocate(temp(nVecT))
            !
            phi%v(1,:,:) = 1
            phi%v(phi%nx+1,:,:) = 1
            phi%v(:,1,:) = 1
            phi%v(:,phi%ny+1,:) = 1
            phi%v(:,:,1) = 1
            phi%v(:,:,phi%nz+1) = 1
            !
            call phi%getArray(temp)
            !
            case default
                stop "Error: intBdryIndicesCScalar3D_SG: Unknown self%grid_type"
        end select
        !
        nBdry = 0
        do i = 1, nVecT
             nBdry = nBdry + nint(real(temp(i)))
        end do
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
             if( nint(real(temp(i))).eq.1) then
                nb = nb+1
                ind_b(nb) = i
             else
                ni = ni+1
                ind_i(ni) = i
             end if
        end do
        !
        deallocate( temp )
        !
    end subroutine intBdryIndicesCScalar3D_SG
    !
    function lengthCScalar3D_SG( self ) result( field_length )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( in ) :: self
        !
        integer :: field_length
        !
        field_length = self%Nxyz
        !
    end function lengthCScalar3D_SG
    !
    subroutine getRealArrayCScalar3D_SG( self, array )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( in )         :: self
        real( kind=prec ), allocatable, intent( out ) :: array(:)
        !
        allocate(array(self%length()))
        array = (/reshape(real( self%v%re, kind=prec ), (/self%Nxyz, 1/))/)
        !
    end subroutine getRealArrayCScalar3D_SG
    !
    subroutine getComplexArrayCScalar3D_SG( self, array )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( in )            :: self
        complex( kind=prec ), allocatable, intent( out ) :: array(:)
        !
        allocate(array(self%length()))
        array = (/reshape(self%v, (/self%Nxyz, 1/))/)
        !
    end subroutine getComplexArrayCScalar3D_SG
    !
    subroutine setRealArrayCScalar3D_SG( self, array )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( inout ) :: self
        real( kind=prec ), intent( in )          :: array(:)
        !
        self%v = reshape( cmplx( array, 0.0, kind=prec ), (/self%NdV(1), self%NdV(2), self%NdV(3)/))
        !
    end subroutine setRealArrayCScalar3D_SG
    !
    subroutine setComplexArrayCScalar3D_SG( self, array )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in )       :: array(:)
        !
        self%v = reshape(array, (/self%NdV(1), self%NdV(2), self%NdV(3)/))
        !
    end subroutine setComplexArrayCScalar3D_SG
    !
    subroutine setVecComponentsCScalar3D_SG( self, xyz, &
                                             xmin, xstep, xmax, &
                                             ymin, ystep, ymax, &
                                             zmin, zstep, zmax, cvalue )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( inout ) :: self
        character, intent( in )                  :: xyz
        integer, intent( in )                    :: xmin, xstep, xmax
        integer, intent( in )                    :: ymin, ystep, ymax
        integer, intent( in )                    :: zmin, zstep, zmax
        complex( kind=prec ), intent( in ) :: cvalue
        !
        integer :: x1, x2
        integer :: y1, y2
        integer :: z1, z2
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
        self%v(x1:x2:xstep, y1:y2:ystep, z1:z2:zstep) = cvalue
        !
    end subroutine setVecComponentsCScalar3D_SG
    !
    subroutine zerosCScalar3D_SG( self )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( inout ) :: self
        !
        if( .NOT. self%is_allocated) then
             stop "Error: zerosCScalar3D_SG > Not allocated."
        end if
        !
        self%v = C_ZERO
        !
    end subroutine zerosCScalar3D_SG
    !
    subroutine addCScalar3D_SG( self, rhs )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in )           :: rhs
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type( rhs )
                class is( cScalar3D_SG_t )
                    self%v = self%v + rhs%v
                class is( rScalar3D_SG_t )
                    self%v = self%v + cmplx( rhs%v, 0.0, kind=prec )
                class default
                    stop "Error: addCScalar3D_SG: undefined rhs"
            end select
            !
        else
            stop "Error: addCScalar3D_SG > Incompatible inputs."
        end if
        !
    end subroutine addCScalar3D_SG
    !
    subroutine subCScalar3D_SG( self, rhs )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in )           :: rhs
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type( rhs )
                class is( cScalar3D_SG_t )
                    self%v = self%v - rhs%v
                class is( rScalar3D_SG_t )
                    self%v = self%v - cmplx( rhs%v, 0.0, kind=prec )
                class default
                    stop "Error: subCScalar3D_SG: undefined rhs"
            end select
            !
        else
            stop "Error: subCScalar3D_SG > Incompatible inputs."
        end if
        !
    end subroutine subCScalar3D_SG
    !
    subroutine multByFieldCScalar3D_SG( self, rhs )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in )           :: rhs
        !
        if( self%isCompatible(rhs)) then
            !
            select type( rhs )
                class is( cScalar3D_SG_t )
                    self%v = self%v * rhs%v
                class is( rScalar3D_SG_t )
                    self%v = self%v * cmplx( rhs%v, 0.0, kind=prec )
                class default
                    stop "Error: multByFieldCScalar3D_SG: undefined rhs"
            end select
            !
        else
            stop "Error: multByFieldCScalar3D_SG: incompatible rhs"
        end if
        !
    end subroutine multByFieldCScalar3D_SG
    !
    subroutine multByValueCScalar3D_SG( self, cvalue )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in )       :: cvalue
        !
        self%v = self%v * cvalue
        !
    end subroutine multByValueCScalar3D_SG
    !
    subroutine divByFieldCScalar3D_SG( self, rhs )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in )           :: rhs
        !
        if( self%isCompatible(rhs)) then
            !
            select type( rhs )
                class is( cScalar3D_SG_t )
                    self%v = self%v / rhs%v
                class is( rScalar3D_SG_t )
                    self%v = self%v / cmplx( rhs%v, 0.0, kind=prec )
                class default
                    stop "Error: divByFieldCScalar3D_SG: undefined rhs"
            end select
            !
        else
            stop "Error: divByFieldCScalar3D_SG: incompatible rhs"
        end if
        !
    end subroutine divByFieldCScalar3D_SG
    !
    subroutine divByValueCScalar3D_SG( self, cvalue )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in )       :: cvalue
        !
        self%v = self%v / cvalue
        !
    end subroutine divByValueCScalar3D_SG
    !
    function dotProdCScalar3D_SG( self, rhs ) result( cvalue )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( in ) :: self
        class( Scalar_t ), intent( in )       :: rhs
        complex( kind=prec )                  :: cvalue
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type( rhs )
                class is( cScalar3D_SG_t )
                    cvalue = sum( conjg( self%v ) * rhs%v )
                class default
                    stop "Error: dotProdCScalar3D_SG > undefined rhs"
            end select
            !
        else
            stop "Error: dotProdCScalar3D_SG > Incompatible rhs"
        end if
        !
    end function dotProdCScalar3D_SG
    !
    subroutine linCombSCScalar3D_SG( self, rhs, c1, c2 )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( inout ) :: self
        class( Scalar_t ), intent( in )          :: rhs
        complex( kind=prec ), intent( in )       :: c1, c2
        !
        !  linear combination, in place: self = c1*self+c2*rhs
        if( self%isCompatible(rhs)) then
            !
            select type( rhs )
                class is( cScalar3D_SG_t )
                    self%v = c1 * self%v + c2 * rhs%v
                class default
                    stop "Error: linCombSCScalar3D_SG: undefined rhs"
            !
            end select
        else
            stop "Error: linCombSCScalar3D_SG > Incompatible rhs"
        end if
        !
    end subroutine linCombSCScalar3D_SG
    !
    subroutine scMultAddSCScalar3D_SG( self, rhs, cvalue )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( in ) :: self
        class( Scalar_t ), intent( inout )    :: rhs
        complex( kind=prec ), intent( in )    :: cvalue
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type( rhs )
                class is( cScalar3D_SG_t )
                    rhs%v = rhs%v + cvalue * self%v
                class default
                    stop "Error: scMultAddSCScalar3D_SG: undefined rhs"
            end select
            !
        else
            stop "Error: scMultAddSCScalar3D_SG > Incompatible rhs"
        end if
        !
    end subroutine scMultAddSCScalar3D_SG
    !
    subroutine copyFromCScalar3D_SG( self, rhs )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in )           :: rhs
        !
        if( .NOT. rhs%is_allocated ) then
            stop "Error: copyFromCScalar3D_SG > rhs not allocated"
        end if
        !
        self%grid => rhs%grid
        self%grid_type = rhs%grid_type
        self%nx = rhs%nx
        self%ny = rhs%ny
        self%nz = rhs%nz
        !
        select type( rhs )
            class is( cScalar3D_SG_t )
                !
                self%NdV = rhs%NdV
                self%Nxyz = rhs%Nxyz
                !
                self%v = rhs%v
                !
            class default
                stop "Error: copyFromCScalar3D_SG > Incompatible rhs"
        end select
        !
        self%is_allocated = .TRUE.
        !
    end subroutine copyFromCScalar3D_SG

    subroutine printCScalar3D_SG( self, io_unit, title, append )
        implicit none
        !
        class( cScalar3D_SG_t ), intent( in )             :: self
        integer, intent( in ), optional                   :: io_unit
        character(:), allocatable, intent( in ), optional :: title
        logical, intent( in ), optional                   :: append
        !
        integer :: ix, iy, iz,funit
        !
        if( present( io_unit ) ) then
            funit = io_unit
        else
            funit = 0    !    usually this will work to write to standard output
        endif
        if(present(title)) then
          write(funit,*) title
        end if
        !
        write( funit, * ) self%nx, self%ny, self%nz
        !
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
    end subroutine printCScalar3D_SG
    !
end module cScalar3D_SG
