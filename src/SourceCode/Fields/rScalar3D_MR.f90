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
            procedure, public :: initializeSub => initializeSub_rScalar3D_MR
            !
            procedure, public :: setIndexArrays => setIndexArrays_rScalar3D_MR
            !
            procedure, public :: setFull => setFull_rScalar3D_MR
            procedure, public :: getFull => getFull_rScalar3D_MR
            !
            procedure, public :: lengthFull => lengthFull_rScalar3D_MR
            procedure, public :: findFull => findFull_rScalar3D_MR
            !
            procedure, public :: findValue => findValue_rScalar3D_MR
            !
            procedure, public :: sgToMR => sgTo_rScalar3D_MR
            !
            !> Boundary operations
            procedure, public :: setAllBoundary => setAllBoundary_rScalar3D_MR
            procedure, public :: setOneBoundary => setOneBoundary_rScalar3D_MR
            procedure, public :: intBdryIndices => intBdryIndices_rScalar3D_MR
            !
            !> Dimensioning operations
            procedure, public :: length => length_rScalar3D_MR
            procedure, public :: setVecComponents => setVecComponents_rScalar3D_MR
            !
            !> Arithmetic/algebraic unary operations
            procedure, public :: zeros => zeros_rScalar3D_MR
            procedure, public :: conjugate => conjugate_rScalar3D_MR
            !
            !> Arithmetic/algebraic binary operations
            procedure, public :: add => add_rScalar3D_MR
            !
            procedure, public :: linComb => linComb_rScalar3D_MR
            !
            procedure, public :: subValue => subValue_rScalar3D_MR
            procedure, public :: subField => subField_rScalar3D_MR
            !
            procedure, public :: multByReal => multByReal_rScalar3D_MR
            procedure, public :: multByComplex => multByComplex_rScalar3D_MR
            procedure, public :: multByField => multByField_rScalar3D_MR
            !
            procedure, public :: multAdd => multAdd_rScalar3D_MR
            !
            procedure, public :: dotProd => dotProd_rScalar3D_MR
            !
            procedure, public :: divByField => divByField_rScalar3D_MR
            procedure, public :: divByValue => divByValue_rScalar3D_MR
            !
            !> Miscellaneous
            procedure, public :: getV => getV_rScalar3D_MR
            procedure, public :: setV => setV_rScalar3D_MR
            !
            procedure, public :: getSV => getSV_rScalar3D_MR
            procedure, public :: setSV => setSV_rScalar3D_MR
            !
            procedure, public :: getArray => getArray_rScalar3D_MR
            procedure, public :: setArray => setArray_rScalar3D_MR
            !
            procedure, public :: copyFrom => copyFrom_rScalar3D_MR
            !
            !> I/O operations
            procedure, public :: read => read_rScalar3D_MR
            procedure, public :: write => write_rScalar3D_MR
            procedure, public :: print => print_rScalar3D_MR
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
    subroutine initializeSub_rScalar3D_MR( self ) 
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
                    self%sub_scalars(i) = rScalar3D_MR_t( grid%sub_grids(i), self%grid_type )
                end do
                !
                call self%setIndexArrays
                call self%zeros
                !
            class default
                stop "Error: initializeSub_rScalar3D_MR > Unclassified grid"
            !
        end select
        !
    end subroutine initializeSub_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setIndexArrays_rScalar3D_MR( self, xy_in ) 
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
                        stop "Error: setIndexArrays_rScalar3D_MR > Invalid grid type option!"
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
                stop "Error: setIndexArrays_rScalar3D_MR > Unclassified grid"
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
    end subroutine setIndexArrays_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setFull_rScalar3D_MR( self, v )
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
                stop "Error: setFull_rScalar3D_MR > Unclassified grid"
            !
        end select
        !
    end subroutine setFull_rScalar3D_MR
    !
    !> Creates standard (1-D array) for all sub-scalars,
    !> INCLUDING redundant interface nodes.
    !
    subroutine getFull_rScalar3D_MR( self, v )
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
                stop "Error: getFull_rScalar3D_MR > Unclassified grid"
            !
        end select
        !
    end subroutine getFull_rScalar3D_MR
    !
    !> No function briefing
    !
    function lengthFull_rScalar3D_MR( self ) result( n )
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
                stop "Error: lengthFull_rScalar3D_MR > Unclassified grid"
            !
        end select
        !
    end function lengthFull_rScalar3D_MR
    !
    !> No function briefing
    !
    function findFull_rScalar3D_MR( self, c ) result( I )
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
    end function findFull_rScalar3D_MR
    !
    !> No function briefing
    !
    function findValue_rScalar3D_MR(self, c) result (I)
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
    end function findValue_rScalar3D_MR
    !
    subroutine sgTo_rScalar3D_MR( self, sg_v )
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
                stop "Error: lengthFull_rScalar3D_MR > Unclassified grid"
            !
        end select
        !
    end subroutine sgTo_rScalar3D_MR
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
    end subroutine rScalar3D_MR_dtor
    !
    !> No subroutine briefing
    !
    subroutine setAllBoundary_rScalar3D_MR( self, cvalue )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        stop "Error: setAllBoundary_rScalar3D_MR not implemented!"
        !
    end subroutine setAllBoundary_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setOneBoundary_rScalar3D_MR( self, bdry, cvalue, int_only )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        character(*), intent( in ) :: bdry
        complex( kind=prec ), intent( in ) :: cvalue
        logical, intent( in ), optional :: int_only
        !
        logical :: int_only_p
        !
        call self%switchStoreState( compound )
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
                 stop "Error: setOneBoundary_rScalar3D_MR > Invalid grid type"
        end select
        !
    end subroutine setOneBoundary_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine intBdryIndices_rScalar3D_MR( self, ind_i, ind_b )
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
    end subroutine intBdryIndices_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    function length_rScalar3D_MR( self ) result( field_length )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( in ) :: self
        !
        integer :: field_length
        !
        field_length = size( self%ind_active )
        !
    end function length_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setVecComponents_rScalar3D_MR( self, xyz, &
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
    end subroutine setVecComponents_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine zeros_rScalar3D_MR( self )
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
                stop "Error: zeros_rScalar3D_MR > Unclassified grid"
            !
        end select
        !
    end subroutine zeros_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine conjugate_rScalar3D_MR( self )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        !
        stop "Error: conjugate_rScalar3D_MR: Do not try to conjugate a real scalar!"
        !
    end subroutine conjugate_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine add_rScalar3D_MR( self, rhs )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        stop "Error: add_rScalar3D_MR not implemented!"
        !
    end subroutine add_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine linComb_rScalar3D_MR( self, rhs, c1, c2 )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        complex( kind=prec ), intent( in ) :: c1, c2
        !
        stop "Error: linComb_rScalar3D_MR not implemented!"
        !
    end subroutine linComb_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine subValue_rScalar3D_MR( self, cvalue )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        stop "Error: subValue_rScalar3D_MR not implemented!"
        !
    end subroutine subValue_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine subField_rScalar3D_MR( self, rhs )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        stop "Error: subField_rScalar3D_MR not implemented!"
        !
    end subroutine subField_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine multByReal_rScalar3D_MR( self, rvalue )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: rvalue
        !
        stop "Error: multByReal_rScalar3D_MR not implemented!"
        !
    end subroutine multByReal_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine multByComplex_rScalar3D_MR( self, cvalue )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        stop "Error: multByComplex_rScalar3D_MR not implemented!"
        !
    end subroutine multByComplex_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine multByField_rScalar3D_MR( self, rhs )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        stop "Error: multByField_rScalar3D_MR not implemented!"
        !
    end subroutine multByField_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine multAdd_rScalar3D_MR( self, cvalue, rhs )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        class( Field_t ), intent( in ) :: rhs
        !
        stop "Error: multAdd_rScalar3D_MR not implemented!"
        !
    end subroutine multAdd_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    function dotProd_rScalar3D_MR( self, rhs ) result( cvalue )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( in ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        complex( kind=prec ) :: cvalue
        !
        call errStop( "dotProd_rScalar3D_MR still not implemented" )
        !
    end function dotProd_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine divByValue_rScalar3D_MR( self, cvalue )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        stop "Error: divByValue_rScalar3D_MR not implemented!"
        !
    end subroutine divByValue_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine divByField_rScalar3D_MR( self, rhs )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        call errStop( "divByField_rScalar3D_MR not implemented!" )
        !
    end subroutine divByField_rScalar3D_MR
    !
    !> No function briefing
    !
    function getV_rScalar3D_MR( self ) result( v )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( in ) :: self
        !
        complex( kind=prec ), allocatable :: v(:,:,:)
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "getV_rScalar3D_MR > self not allocated." )
        endif
        !
        if( .NOT. allocated( v ) ) then
            call errStop( "getV_rScalar3D_MR > v not allocated." )
        else
            !
            v = cmplx( self%v, 0.0, kind=prec )
            !
        endif
        !
    end function getV_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setV_rScalar3D_MR( self, v )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), allocatable, intent( in ) :: v(:,:,:)
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "setV_rScalar3D_MR > self not allocated." )
        endif
        !
        if( .NOT. allocated( v ) ) then
            call errStop( "setV_rScalar3D_MR > v not allocated." )
        endif
        !
        call self%switchStoreState( compound )
        !
        if( allocated( self%s_v ) ) deallocate( self%s_v )
        !
        self%v = real( v, kind=prec )
        !
    end subroutine setV_rScalar3D_MR
    !
    !> No function briefing
    !
    function getSV_rScalar3D_MR( self ) result( s_v )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( in ) :: self
        !
        complex( kind=prec ), allocatable :: s_v(:)
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "getSV_rScalar3D_MR > self not allocated." )
        endif
        !
        if( .NOT. allocated( self%s_v ) ) then
            call errStop( "getSV_rScalar3D_MR > self%s_v not allocated." )
        else
            !
            s_v = cmplx( self%s_v, 0.0, kind=prec )
            !
        endif
        !
    end function getSV_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setSV_rScalar3D_MR( self, s_v )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), allocatable, intent( in ) :: s_v(:)
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "setSV_rScalar3D_MR > self not allocated." )
        endif
        !
        if( .NOT. allocated( s_v ) ) then
            call errStop( "setSV_rScalar3D_MR > s_v not allocated." )
        endif
        !
        call self%switchStoreState( singleton )
        !
        if( allocated( self%v ) ) deallocate( self%v )
        !
        self%s_v = s_v
        !
    end subroutine setSV_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    function getArray_rScalar3D_MR( self ) result( array )
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
    end function getArray_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setArray_rScalar3D_MR( self, array )
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
    end subroutine setArray_rScalar3D_MR
    !
    subroutine copyFrom_rScalar3D_MR( self, rhs )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        integer :: i
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "copyFrom_rScalar3D_MR > rhs not allocated" )
        endif
        !
        self%grid => rhs%grid
        self%grid_type = rhs%grid_type
        self%nx = rhs%nx
        self%ny = rhs%ny
        self%nz = rhs%nz
        self%store_state = rhs%store_state
        !
        select type( rhs )
            !
            class is( rScalar3D_MR_t )
                !
                if( allocated( rhs%ind_active ) ) then
                    self%ind_active = rhs%ind_active
                else
                    call errStop( "copyFrom_rScalar3D_MR > rhs%ind_active not allocated" )
                endif
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
                    self%s_v = rhs%s_v
                    !
                else
                    stop "Error: copyFrom_rScalar3D_MR > Unknown store_state!"
                endif
                !
                do i = 1, size( self%sub_scalars )
                    self%sub_scalars(i) = rhs%sub_scalars(i)
                end do
                !
                self%is_allocated = .TRUE.
                !
                call self%setIndexArrays
                !
            class default
                stop "Error: copyFrom_rScalar3D_MR > Unclassified rhs"
            !
        end select
        !
    end subroutine copyFrom_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine read_rScalar3D_MR( self, funit, ftype )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        integer, intent( in ) :: funit
        character(:), allocatable, intent( in ), optional :: ftype
        !
        stop "Error: read_rScalar3D_MR not implemented!"
        !
    end subroutine read_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine write_rScalar3D_MR( self, funit, ftype )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        integer, intent( in ) :: funit
        character(:), allocatable, intent( in ), optional :: ftype
        !
        stop "Error: write_rScalar3D_MR not implemented!"
        !
    end subroutine write_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine print_rScalar3D_MR( self, io_unit, title, append )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        integer, intent( in ), optional :: io_unit
        character(*), intent( in ), optional :: title
        logical, intent( in ), optional :: append
        !
        integer :: ix, iy, iz,funit
        !
        call self%switchStoreState( compound )
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
    end subroutine print_rScalar3D_MR
    !
end module rScalar3D_MR
