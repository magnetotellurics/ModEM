!
!> Derived class to define a rScalar3D_MR 
!
module rScalar3D_MR
    !
    use rScalar3D_SG
    use Grid3D_MR
    !
    type, extends( rScalar3D_SG_t ) :: rScalar3D_MR_t
        !
        type( rScalar3D_SG_t ), allocatable :: sub_scalars(:)
        !
        integer, dimension(:), allocatable :: ind_active
        !
        contains
            !
            !> Destructor
            final :: rScalar3D_MR_dtor
            !
            !> MR Routines
            procedure, public :: initializeSub => initializeSubRScalar3D_MR
            !
            procedure, public :: setIndexArrays => setIndexArraysRScalar3D_MR
            !
            procedure, public :: setFull => setFullRScalar3D_MR
            procedure, public :: getFull => getFullRScalar3D_MR
            !
            procedure, public :: lengthFull => lengthFullRScalar3D_MR
            procedure, public :: findFull => findFullRScalar3D_MR
            !
            procedure, public :: findValue => findValueRScalar3D_MR
            !
            procedure, public :: sgToMR => sgToRScalar3D_MR
            !
            !> Boundary operations
            procedure, public :: setAllBoundary => setAllBoundaryRScalar3D_MR
            procedure, public :: setOneBoundary => setOneBoundaryRScalar3D_MR
            procedure, public :: intBdryIndices => intBdryIndicesRScalar3D_MR
            !
            !> Dimensioning operations
            procedure, public :: length => lengthRScalar3D_MR
            procedure, public :: setVecComponents => setVecComponentsRScalar3D_MR
            !
            !> Arithmetic/algebraic unary operations
            procedure, public :: zeros => zerosRScalar3D_MR
            procedure, public :: sumEdges => sumEdgesRScalar3D_MR
            procedure, public :: avgCells => avgCellsRScalar3D_MR
            procedure, public :: conjugate => conjugateRScalar3D_MR
            !
            !> Arithmetic/algebraic binary operations
            procedure, public :: add => addRScalar3D_MR
            !
            procedure, public :: linComb => linCombRScalar3D_MR
            !
            procedure, public :: subValue => subValueRScalar3D_MR
            procedure, public :: subField => subFieldRScalar3D_MR
            !
            procedure, public :: multByReal => multByRealRScalar3D_MR
            procedure, public :: multByComplex => multByComplexRScalar3D_MR
            procedure, public :: multByField => multByFieldRScalar3D_MR
            !
            procedure, public :: multAdd => multAddRScalar3D_MR
            !
            procedure, public :: dotProd => dotProdRScalar3D_MR
            !
            procedure, public :: divByField => divByFieldRScalar3D_MR
            procedure, public :: divByValue => divByValueRScalar3D_MR
            !
            !> Miscellaneous
            procedure, public :: getReal => getRealRScalar3D_MR
            procedure, public :: getArray => getArrayRScalar3D_MR
            procedure, public :: setArray => setArrayRScalar3D_MR
            procedure, public :: switchStoreState => switchStoreStateRScalar3D_MR
            procedure, public :: copyFrom => copyFromRScalar3D_MR
            !
            !> I/O operations
            procedure, public :: read => readRScalar3D_MR
            procedure, public :: write => writeRScalar3D_MR
            procedure, public :: print => printRScalar3D_MR
            !
    end type rScalar3D_MR_t
    !
    interface rScalar3D_MR_t
        module procedure rScalar3D_MR_ctor_copy
        module procedure rScalar3D_MR_ctor_default
    end interface rScalar3D_MR_t
    !
contains
    !
    !> No function briefing
    !
    function rScalar3D_MR_ctor_copy( E_in ) result ( self )
        implicit none
        !
        type( rScalar3D_MR_t ), intent( in ) :: E_in
        !
        type( rScalar3D_MR_t ) :: self
        !
        integer :: i
        !
        self%grid => E_in%grid
        self%grid_type = E_in%grid_type
        !
        call self%initializeSub
        !
        select type( grid => E_in%grid )
            !
            class is( Grid3D_MR_t )
                !
                do i = 1, grid%n_grids
                    self%sub_scalars(i) = E_in%sub_scalars(i)
                end do
                !
            class default
                stop "Error: rScalar3D_MR_ctor_copy > Unclassified grid"
            !
        end select
        !
    end function rScalar3D_MR_ctor_copy
    !
    !> No function briefing
    !
    function rScalar3D_MR_ctor_default( grid, grid_type ) result ( self )
        implicit none
        !
        class( Grid_t ), target, intent( in ) :: grid
        character( len=4 ), intent( in ) :: grid_type
        !
        type( rScalar3D_MR_t ) :: self
        !
        self%grid => grid
        self%grid_type = grid_type
        !
        call self%initializeSub
        !
    end function rScalar3D_MR_ctor_default
    !
    !> No subroutine briefing
    !
    subroutine initializeSubRScalar3D_MR( self ) 
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        !
        integer :: i, status
        !
        select type( grid => self%grid )
            !
            class is( Grid3D_MR_t )
                !
                self%is_allocated = .TRUE.
                allocate( self%sub_scalars( grid%n_grids ), STAT = status )
                self%is_allocated = self%is_allocated .AND. ( status .EQ. 0 )
                !
                do i = 1, grid%n_grids
                    self%sub_scalars(i) = rScalar3D_MR_t( grid%sub_grids(i), grid_type )
                end do
                !
                call self%setIndexArrays
                call self%zeros
                !
            class default
                stop "Error: initializeSubRScalar3D_MR > Unclassified grid"
            !
        end select
        !
    end subroutine initializeSubRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setIndexArraysRScalar3D_MR( self, xy_in ) 
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        logical, intent( in ), optional :: xy_in
        !
        logical :: xy, int_only
        integer :: i, k
        integer :: n_full, n_active, n_interior, n_boundaries
        real( kind=prec ), dimension(:), allocatable :: v_1, v_2
        !
        if ( .NOT. present( xy_in ) ) then
            xy = .FALSE.
        else
            xy = xy_in
        endif
        !
        select type( grid => self%grid )
            !
            class is( Grid3D_MR_t )
                !
                ! Loop over sub-grids, setting boundary edges to one,
                ! interior to  zero
                do k = 1, grid%n_grids
                    call self%sub_scalars(k)%setAllBoundary( cmplx( 1._prec, 0.0, kind=prec ) )
                end do
                !
                ! Loop over interfaces: set redundant interface edges to 2
                select case( self%grid_type )
                    !
                    case( EDGE )
                        int_only = .TRUE.
                    case( FACE )
                        int_only = .FALSE.
                    case( NODE )
                        int_only = .TRUE.
                    case default
                        !
                        stop "Error: setIndexArraysRScalar3D_MR > Invalid grid type option!"
                    !
                end select
                !
                do k = 2, grid%n_grids
                    !
                    if( grid%Coarseness(k - 1, 1) < grid%Coarseness(k, 1) ) then
                        ! upper grid is finer: grid k-1 interface nodes are
                        ! not active; also reset interior part of interface
                        ! edges to 0
                        if( xy ) then
                            call self%sub_scalars(k-1)%setOneBoundary( "z2_x", cmplx( -1.0_prec, 0.0, kind=prec ) )
                            call self%sub_scalars(k-1)%setOneBoundary( "z2_y", cmplx( -10.0_prec, 0.0, kind=prec ) )
                        else
                            call self%sub_scalars(k-1)%setOneBoundary( "z2", cmplx( -1.0_prec, 0.0, kind=prec ) )
                        endif
                        !
                        call self%sub_scalars(k)%setOneBoundary( "z1", cmplx( 0._prec, 0.0, kind=prec ), int_only )
                    else
                        if( xy ) then
                            call self%sub_scalars(k)%setOneBoundary( "z1_x", cmplx( -1.0_prec, 0.0, kind=prec ) )
                            call self%sub_scalars(k)%setOneBoundary( "z1_y", cmplx( -10.0_prec, 0.0, kind=prec ) )
                        else
                            call self%sub_scalars(k)%setOneBoundary( "z1", cmplx( -1.0_prec, 0.0, kind=prec ) )
                        endif
                        !
                        call self%sub_scalars(k-1)%setOneBoundary( "z2", cmplx( 0._prec, 0.0, kind=prec ), int_only )
                        !
                    endif
                    !
                end do
                !
            class default
                stop "Error: setIndexArraysRScalar3D_MR > Unclassified grid"
            !
        end select
        !
        ! Set active, interior, and boundary edges. ***
        !
        call self%getFull( v_1 )
        !
        n_full = size( v_1 )
        !
        n_active = 0
        do k = 1, n_full
            if (v_1(k) >= 0) then
                n_active = n_active + 1
            endif
        end do
        !
        if (allocated (self%ind_active)) then
            deallocate (self%ind_active)
        endif
        !
        allocate (self%ind_active(n_active))
        !
        i = 0
        do k = 1, n_full
            if (v_1(k) >= 0) then
                i = i + 1
                self%ind_active(i) = k
            endif
        end do
        !
        n_interior = 0
        do k = 1, n_full
            if (v_1(k) == 0) then
                n_interior = n_interior + 1
            endif
        end do
        !
        allocate (v_2(n_active))
        v_2 = v_1(self%ind_active)
        !
        if( allocated( self%ind_interior ) ) then
            deallocate( self%ind_interior )
        endif
        !
        allocate( self%ind_interior( n_interior ) )
        !
        i = 0
        do k = 1, n_active
            if( v_2(k) == 0 ) then
                i = i + 1
                self%ind_interior(i) = k
            endif
        end do
        !!
        n_boundaries = 0
        do k = 1, n_active
            if (v_2(k) == 1) then
                n_boundaries = n_boundaries + 1
            endif
        end do
        !
        if (allocated (self%ind_boundaries)) then
            deallocate (self%ind_boundaries)
        endif
        !
        allocate (self%ind_boundaries(n_boundaries)) 
        !
        i = 0
        do k = 1, n_active
            if (v_2(k) == 1) then
                i = i + 1
                self%ind_boundaries(i) = k
            endif
        end do
        !
    end subroutine setIndexArraysRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setFullRScalar3D_MR( self, v )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        real( kind=prec ), intent (in) :: v(:)
        !
        integer :: i1, i2, k, n
        !
        select type( grid => self%grid )
            !
            class is( Grid3D_MR_t )
                !
                i1 = 1; i2 = 0;
                do k = 1, grid%n_grids
                    n = self%sub_scalars(k)%length ()
                    i2 = i2 + n
                    call self%sub_scalars(k)%setArray( cmplx( v(i1:i2), 0.0, kind=prec ) )
                    i1 = i1 + n
                end do
                !
            class default
                stop "Error: setFullRScalar3D_MR > Unclassified grid"
            !
        end select
        !
    end subroutine setFullRScalar3D_MR
    !
    !> Creates standard (1-D array) for all sub-scalars,
    !> INCLUDING redundant interface nodes.
    !
    subroutine getFullRScalar3D_MR( self, v )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( in ) :: self
        real( kind=prec ), allocatable, intent( out ) :: v(:)
        !
        real( kind=prec ), allocatable :: v_temp(:)
        integer :: n, i1, i2, k
        !
        n = self%lengthFull()
        allocate( v(n) )
        !
        v = R_ZERO
        !
        i1 = 1
        i2 = 0
        !
        select type( grid => self%grid )
            !
            class is( Grid3D_MR_t )
                !
                do k = 1, grid%n_grids
                    n = self%sub_scalars(k)%length()
                    i2 = i2 + n
                    v_temp = self%sub_scalars(k)%getArray()
                    v(i1:i2) = v_temp
                    i1 = i1 + n
                end do
                !
            class default
                stop "Error: getFullRScalar3D_MR > Unclassified grid"
            !
        end select
        !
    end subroutine getFullRScalar3D_MR
    !
    !> No function briefing
    !
    function lengthFullRScalar3D_MR( self ) result( n )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( in ) :: self
        !
        integer :: k
        integer :: n
        !
        select type( grid => self%grid )
            !
            class is( Grid3D_MR_t )
                !
                n = 0
                do k = 1, grid%n_grids
                    n = n + self%sub_scalars(k)%length()
                end do
                !
            class default
                stop "Error: lengthFullRScalar3D_MR > Unclassified grid"
            !
        end select
        !
    end function lengthFullRScalar3D_MR
    !
    !> No function briefing
    !
    function findFullRScalar3D_MR( self, c ) result( I )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( in ) :: self
        real( kind=prec ), intent (in) :: c
        !
        integer, dimension(:), allocatable :: I
        real( kind=prec ), dimension(:), allocatable :: v
        integer :: n, n_I, k
        !
        n = self%lengthFull ()
        call self%getFull(v)
        !
        n_I = 0
        do k = 1, n
            if (v(k) == c) n_I = n_I + 1
        end do
        !
        allocate (I(n_I))
        !
        n_I = 0
        do k = 1, n
            if (v(k) == c) then
                n_I = n_I + 1
                I(n_I) = k
            endif
        end do
        !
    end function findFullRScalar3D_MR
    !
    !> No function briefing
    !
    function findValueRScalar3D_MR(self, c) result (I)
        implicit none
        !
        class( rScalar3D_MR_t ), intent( in ) :: self
        real( kind=prec ), intent( in ) :: c
        !
        integer, dimension(:), allocatable :: I
        real( kind=prec ), dimension(:), allocatable :: v
        integer :: n, n_I, k
        !
        n = self%length()
        allocate (v(n))
        v = self%getArray()
        !
        n_I = 0
        do k = 1, n
            if (v(k) == c) n_I = n_I + 1
        end do
        !
        allocate (I(n_I))
        !
        n_I = 0
        do k = 1, n
            if (v(k) == c) then
                n_I = n_I + 1
                I(n_I) = k
            endif
        end do
        !
    end function findValueRScalar3D_MR
    !
    subroutine sgToRScalar3D_MR( self, sg_v )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        type ( rScalar3D_SG_t ), intent( in ) :: sg_v
        !
        class( Grid_t ), pointer :: grid
        !
        integer :: x_nx, x_ny, x_nz
        integer :: last, Cs, i1, i2, i, k
        !
        grid => sg_v%grid
        !
        x_nx = size(sg_v%v, 1)
        x_ny = size(sg_v%v, 2)
        x_nz = size(sg_v%v, 3)
        !
        select type( grid => self%grid )
            !
            class is( Grid3D_MR_t )
                !
                do k = 1, grid%n_grids
                    !
                    Cs = 2**grid%Coarseness(k, 1)
                    i1 = grid%Coarseness(k, 3)
                    i2 = grid%Coarseness(k, 4)
                    !
                    do i = 1, Cs
                        !
                        last = size(grid%Dx)
                        self%sub_scalars(k)%v = self%sub_scalars(k)%v + &
                        sg_v%v(i:x_nx:Cs, 1:x_ny:Cs, i1:i2+1) *    &
                        repMat(grid%Dx(i:last:Cs), 1, &
                        grid%sub_grids(k)%Ny + 1, &
                        grid%sub_grids(k)%Nz + 1, .FALSE. )
                        !
                    end do
                    !
                    self%sub_scalars(k)%v = self%sub_scalars(k)%v / &
                    repMat(grid%sub_grids(k)%Dx, &
                    1, &
                    grid%sub_grids(k)%Ny + 1, &
                    grid%sub_grids(k)%Nz + 1, .FALSE. )
                    !
                end do
                !
            class default
                stop "Error: lengthFullRScalar3D_MR > Unclassified grid"
            !
        end select
        !
    end subroutine sgToRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine rScalar3D_MR_dtor( self )
        implicit none
        !
        type( rScalar3D_MR_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor rScalar3D_MR"
        !
        call self%dealloc
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
    end subroutine rScalar3D_MR_dtor
    !
    !> No subroutine briefing
    !
    subroutine setAllBoundaryRScalar3D_MR( self, cvalue )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        stop "Error: setAllBoundaryRScalar3D_MR not implemented!"
        !
    end subroutine setAllBoundaryRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setOneBoundaryRScalar3D_MR( self, bdry, cvalue, int_only )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
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
        if( .NOT. present (int_only)) then
             int_only_p = .FALSE.
        else 
             int_only_p = int_only
        endif
        !
        select case( self%grid_type )
            case( NODE )
                !
                if( int_only_p ) then
                    !
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
                    !
                 else
                    !
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
            case( FACE )
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
                 stop "Error: setOneBoundaryRScalar3D_MR > Invalid grid type"
        end select
        !
    end subroutine setOneBoundaryRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine intBdryIndicesRScalar3D_MR( self, ind_i, ind_b )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        integer, allocatable, intent( out ) :: ind_i(:), ind_b(:)
        !
        integer :: m, n
        !
        m = size( self%ind_interior )
        n = size( self%ind_boundaries )
        !
        allocate( ind_i(m) )
        allocate( ind_b(n) )
        !
        ind_i = self%ind_interior
        ind_b = self%ind_boundaries
        !
    end subroutine intBdryIndicesRScalar3D_MR
    !
    !> No subroutine briefing
    !
    function lengthRScalar3D_MR( self ) result( field_length )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( in ) :: self
        !
        integer :: field_length
        !
        field_length = size( self%ind_active )
        !
    end function lengthRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setVecComponentsRScalar3D_MR( self, xyz, &
                                             xmin, xstep, xmax, &
                                             ymin, ystep, ymax, &
                                             zmin, zstep, zmax, rvalue )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
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
    end subroutine setVecComponentsRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine zerosRScalar3D_MR( self )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        !
        integer :: i
        !
        select type( grid => self%grid )
            !
            class is( Grid3D_MR_t )
                !
                do i = 1, grid%n_grids
                    call self%sub_scalars(i)%zeros()
                end do
                !
            class default
                stop "Error: zerosRScalar3D_MR > Unclassified grid"
            !
        end select
        !
    end subroutine zerosRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine sumEdgesRScalar3D_MR( self, cell_obj, interior_only )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        class( Field_t ), allocatable, intent( inout ) :: cell_obj
        logical, optional, intent( in ) :: interior_only
        !
        stop "Error: sumEdgesRScalar3D_MR not implemented yet"
        !
    end subroutine sumEdgesRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine avgCellsRScalar3D_MR( self, E_in, ptype )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: E_in
        character(*), intent( in ), optional :: ptype
        !
        stop "Error: avgCellsRScalar3D_MR not implemented yet"
        !
    end subroutine avgCellsRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine conjugateRScalar3D_MR( self )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        !
        stop "Error: conjugateRScalar3D_MR: Do not try to conjugate a real scalar!"
        !
    end subroutine conjugateRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine addRScalar3D_MR( self, rhs )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        stop "Error: addRScalar3D_MR not implemented!"
        !
    end subroutine addRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine linCombRScalar3D_MR( self, rhs, c1, c2 )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        complex( kind=prec ), intent( in ) :: c1, c2
        !
        stop "Error: linCombRScalar3D_MR not implemented!"
        !
    end subroutine linCombRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine subValueRScalar3D_MR( self, cvalue )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        stop "Error: subValueRScalar3D_MR not implemented!"
        !
    end subroutine subValueRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine subFieldRScalar3D_MR( self, rhs )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        stop "Error: subFieldRScalar3D_MR not implemented!"
        !
    end subroutine subFieldRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine multByRealRScalar3D_MR( self, rvalue )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: rvalue
        !
        stop "Error: multByRealRScalar3D_MR not implemented!"
        !
    end subroutine multByRealRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine multByComplexRScalar3D_MR( self, cvalue )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        stop "Error: multByComplexRScalar3D_MR not implemented!"
        !
    end subroutine multByComplexRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine multByFieldRScalar3D_MR( self, rhs )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        stop "Error: multByFieldRScalar3D_MR not implemented!"
        !
    end subroutine multByFieldRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine multAddRScalar3D_MR( self, cvalue, rhs )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        class( Field_t ), intent( in ) :: rhs
        !
        stop "Error: multAddRScalar3D_MR not implemented!"
        !
    end subroutine multAddRScalar3D_MR
    !
    !> No subroutine briefing
    !
    function dotProdRScalar3D_MR( self, rhs ) result( cvalue )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( in ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        complex( kind=prec ) :: cvalue
        !
        stop "Error: dotProdRScalar3D_MR not implemented!"
        !
    end function dotProdRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine divByValueRScalar3D_MR( self, cvalue )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        stop "Error: divByValueRScalar3D_MR not implemented!"
        !
    end subroutine divByValueRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine divByFieldRScalar3D_MR( self, rhs )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        stop "Error: divByFieldRScalar3D_MR not implemented!"
        !
    end subroutine divByFieldRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine getRealRScalar3D_MR( self, r_field )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( in ) :: self
        class( Field_t ), allocatable, intent( out ) :: r_field
        !
        allocate( r_field, source = rScalar3D_MR_t( self%grid, self%grid_type ) )
        !
        call r_field%copyFrom( self )
        !
    end subroutine getRealRScalar3D_MR
    !
    !> No subroutine briefing
    !
    function getArrayRScalar3D_MR( self ) result( array )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( in ) :: self
        complex( kind=prec ), allocatable, dimension(:) :: array
        !
        real( kind=prec ), allocatable :: v_full(:)
        !
        call self%getFull( v_full )
        !
        allocate( array( size( self%ind_active ) ) )
        !
        array = v_full( self%ind_active )
        !
        deallocate( v_full )
        !
    end function getArrayRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setArrayRScalar3D_MR( self, array )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), dimension(:), intent( in ) :: array
        !
        real( kind=prec ), allocatable :: vFull(:)
        !
        allocate( vFull( self%lengthFull() ) )
        !
        vFull( self%ind_active ) = array
        !
        call self%setFull( vFull )
        !
        deallocate( vFull )
        !
    end subroutine setArrayRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine switchStoreStateRScalar3D_MR( self )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
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
                if( self%grid_type == NODE ) then
                    !
                    allocate( self%v( self%nx + 1, self%ny + 1, self%nz + 1 ) )
                    !
                else if( self%grid_type == CELL ) then
                    !
                    allocate( self%v( self%nx, self%ny, self%nz ) )
                    !
                else if( self%grid_type == CELL_EARTH ) then
                    !
                    call self%grid%getDimensions( self%nx, self%ny, self%nz, nzAir )
                    !
                    allocate( self%v( self%nx, self%ny, self%nz - nzAir ) )
                    !
                else
                     write( *, * ) "Error: switchStoreStateRScalar3D_MR > unrecognized grid type: [", self%grid_type, "]"
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
                write( *, * ) "Error: switchStoreStateRScalar3D_MR > Unknown store_state :[", self%store_state, "]"
                stop
            !
        end select
        !
    end subroutine switchStoreStateRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine copyFromRScalar3D_MR( self, rhs )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        integer :: i
        !
        if( .NOT. rhs%is_allocated ) then
            stop "Error: copyFromRScalar3D_MR > rhs not allocated"
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
        if( allocated( rhs%ind_active ) ) &
        self%ind_active = rhs%ind_active
        !
        select type( rhs )
            !
            class is( rScalar3D_MR_t )
                !
                self%NdV = rhs%NdV
                self%Nxyz = rhs%Nxyz
                !
                if( rhs%store_state .EQ. compound ) then
                    !
                    self%v = rhs%v
                    !
                else if( rhs%store_state .EQ. singleton ) then
                    !
                    self%sv = rhs%sv
                    !
                else
                    stop "Error: copyFromRScalar3D_MR > Unknown store_state!"
                endif
                !
                do i = 1, size( self%sub_scalars )
                    self%sub_scalars(i) = rhs%sub_scalars(i)
                end do
                !
                self%is_allocated = .TRUE.
                !
            class default
                stop "Error: copyFromRScalar3D_MR > Unclassified rhs"
            !
        end select
        !
    end subroutine copyFromRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine readRScalar3D_MR( self, funit, ftype )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        integer, intent( in ) :: funit
        character(:), allocatable, intent( in ), optional :: ftype
        !
        stop "Error: readRScalar3D_MR not implemented!"
        !
    end subroutine readRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine writeRScalar3D_MR( self, funit, ftype )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        integer, intent( in ) :: funit
        character(:), allocatable, intent( in ), optional :: ftype
        !
        stop "Error: writeRScalar3D_MR not implemented!"
        !
    end subroutine writeRScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine printRScalar3D_MR( self, io_unit, title, append )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
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
    end subroutine printRScalar3D_MR
    !
end module rScalar3D_MR
