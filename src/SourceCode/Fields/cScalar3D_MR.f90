!
!> Derived class to define a cScalar3D_MR 
!
module cScalar3D_MR
    !
    use rScalar3D_MR
    use cScalar3D_SG
    !
    type, extends( Scalar_t ) :: cScalar3D_MR_t
        !
        type( cScalar3D_SG_t ), allocatable, dimension(:) :: sub_scalar
        !
        contains
            !
            !> Destructor
            final :: cScalar3D_MR_dtor
            !
            !> MR Routines
            procedure, public :: initializeSub => initializeSub_cScalar3D_MR
            !
            procedure, public :: setIndexArrays => setIndexArrays_cScalar3D_MR
            !
            procedure, public :: getArray => getArray_cScalar3D_MR
            procedure, public :: setArray => setArray_cScalar3D_MR
            !
            procedure, public :: getFullArray => getFullArray_cScalar3D_MR
            procedure, public :: setFullArray => setFullArray_cScalar3D_MR
            !
            procedure, public :: lengthFull => lengthFull_cScalar3D_MR
            procedure, public :: findFull => findFull_cScalar3D_MR
            !
            procedure, public :: divFine => divFine_cScalar3D_MR
            !
            procedure, public :: toSG => toSG_cScalar3D_MR
            procedure, public :: fromSG => fromSG_cScalar3D_MR
            !
            !> Boundary operations
            procedure, public :: setAllBoundary => setAllBoundary_cScalar3D_MR
            procedure, public :: setOneBoundary => setOneBoundary_cScalar3D_MR
            !
            !> Dimensioning operations
            procedure, public :: length => length_cScalar3D_MR
            !
            !> Arithmetic/algebraic unary operations
            procedure, public :: zeros => zeros_cScalar3D_MR
            procedure, public :: conjugate => conjugate_cScalar3D_MR
            !
            !> Arithmetic/algebraic binary operations
            procedure, public :: add => add_cScalar3D_MR
            !
            procedure, public :: linComb => linComb_cScalar3D_MR
            !
            procedure, public :: subValue => subValue_cScalar3D_MR
            procedure, public :: subField => subField_cScalar3D_MR
            !
            procedure, public :: multByReal => multByReal_cScalar3D_MR
            procedure, public :: multByComplex => multByComplex_cScalar3D_MR
            procedure, public :: multByField => multByField_cScalar3D_MR
            !
            procedure, public :: multAdd => multAdd_cScalar3D_MR
            !
            procedure, public :: dotProd => dotProd_cScalar3D_MR
            !
            procedure, public :: divByField => divByField_cScalar3D_MR
            procedure, public :: divByValue => divByValue_cScalar3D_MR
            !
            procedure, public :: sumToNode => sumToNode_cScalar3D_MR
            !
            procedure, public :: copyFrom => copyFrom_cScalar3D_MR
            !
            !> I/O operations
            procedure, public :: read => read_cScalar3D_MR
            procedure, public :: write => write_cScalar3D_MR
            procedure, public :: print => print_cScalar3D_MR
            !
    end type cScalar3D_MR_t
    !
    interface cScalar3D_MR_t
        module procedure cScalar3D_MR_ctor
    end interface cScalar3D_MR_t
    !
contains
    !
    !> No function briefing
    !
    function cScalar3D_MR_ctor( grid, grid_type ) result( self )
        implicit none
        !
        class( Grid_t ), target, intent( in ) :: grid
        character( len=4 ), intent( in ) :: grid_type
        !
        type( cScalar3D_MR_t ) :: self
        !
        integer :: nzAir
        !
        call self%baseInit
        !
        self%grid => grid
        self%grid_type = grid_type
        !
        !> Grid dimensions
        call self%grid%getDimensions( self%nx, self%ny, self%nz, nzAir )
        !
        call self%initializeSub
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "cScalar3D_MR_ctor > Unable to allocate self" )
        endif
        !
    end function cScalar3D_MR_ctor
    !
    !> No subroutine briefing
    !
    subroutine initializeSub_cScalar3D_MR( self ) 
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        !
        integer :: i, alloc_stat
        !
        select type( grid => self%grid )
            !
            class is( Grid3D_MR_t )
                !
                self%is_allocated = .TRUE.
                allocate( self%sub_scalar( grid%n_grids ), stat = alloc_stat )
                self%is_allocated = self%is_allocated .AND.( alloc_stat .EQ. 0 )
                !
                do i = 1, grid%n_grids
                    !
                    self%sub_scalar(i) = rScalar3D_SG_t( grid%sub_grid(i), self%grid_type )
                    !
                    !write( *, * ) "cSubScalar", i, "-nx=", self%sub_scalar(i)%nx, ", ny=", self%sub_scalar(i)%ny, "nz=", self%sub_scalar(i)%nz
                    !
                enddo
                !
                !write( *, * ) "cMainScalar-nx=", self%nx, ", ny=", self%ny, "nz=", self%nz, self%nx*self%ny*self%nz
                !
            class default
                call errStop( "initializeSub_cScalar3D_MR > Unclassified grid" )
            !
        end select
        !
    end subroutine initializeSub_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setIndexArrays_cScalar3D_MR( self, n_full, ind_boundary, ind_interior, ind_active, xy_in ) 
        implicit none
        !
        class( cScalar3D_MR_t ), intent( in ) :: self
        integer, intent( inout ) :: n_full
        integer, allocatable, dimension(:), intent( out ) :: ind_boundary, ind_interior
        integer, allocatable, dimension(:), intent( out ), optional :: ind_active
        logical, intent( in ), optional :: xy_in
        !
        type( cScalar3D_MR_t ) :: temp_scalar
        logical :: xy, int_only
        integer :: i, k
        integer :: n_active, n_interior, n_boundaries
        real( kind=prec ), allocatable, dimension(:) :: v_1, v_2
        !
        if( .NOT. present( xy_in ) ) then
            xy = .FALSE.
        else
            xy = xy_in
        endif
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "setIndexArrays_cScalar3D_MR > self not allocated." )
        endif
        !
        temp_scalar = self
        !
        select type( grid => temp_scalar%grid )
            !
            class is( Grid3D_MR_t )
                !
                ! Loop over sub-grids, setting boundary edges to one,
                ! interior to  zero
                do k = 1, grid%n_grids
                    call temp_scalar%sub_scalar(k)%setAllBoundary( cmplx( 1._prec, 0.0, kind=prec ) )
                enddo
                !
                ! Loop over interfaces: set redundant interface edges to 2
                select case( temp_scalar%grid_type )
                    !
                    case( CELL )
                        int_only = .FALSE.
                    case( NODE )
                        int_only = .TRUE.
                    case default
                        !
                        call errStop( "setIndexArrays_cScalar3D_MR > Invalid grid type option!" )
                    !
                end select
                !
                do k = 2, grid%n_grids
                    !
                    if( grid%coarseness(k - 1, 1) < grid%coarseness(k, 1) ) then
                        !
                        ! upper grid is finer: grid k-1 interface nodes are
                        ! not active; also reset interior part of interface
                        ! edges to 0
                        if( xy ) then
                            call temp_scalar%sub_scalar(k-1)%setOneBoundary( "z2_x", cmplx( -1.0_prec, 0.0, kind=prec ) )
                            call temp_scalar%sub_scalar(k-1)%setOneBoundary( "z2_y", cmplx( -10.0_prec, 0.0, kind=prec ) )
                        else
                            call temp_scalar%sub_scalar(k-1)%setOneBoundary( "z2", cmplx( -1.0_prec, 0.0, kind=prec ) )
                        endif
                        !
                        call temp_scalar%sub_scalar(k)%setOneBoundary( "z1", cmplx( 0._prec, 0.0, kind=prec ), int_only )
                    else
                        if( xy ) then
                            call temp_scalar%sub_scalar(k)%setOneBoundary( "z1_x", cmplx( -1.0_prec, 0.0, kind=prec ) )
                            call temp_scalar%sub_scalar(k)%setOneBoundary( "z1_y", cmplx( -10.0_prec, 0.0, kind=prec ) )
                        else
                            call temp_scalar%sub_scalar(k)%setOneBoundary( "z1", cmplx( -1.0_prec, 0.0, kind=prec ) )
                        endif
                        !
                        call temp_scalar%sub_scalar(k-1)%setOneBoundary( "z2", cmplx( 0._prec, 0.0, kind=prec ), int_only )
                        !
                    endif
                    !
                enddo
                !
            class default
                call errStop( "setIndexArrays_cScalar3D_MR > Unclassified grid" )
            !
        end select
        !
        ! Set active, interior, and boundary edges. ***
        !
        v_1 = temp_scalar%getFullArray()
        !
        n_full = size( v_1 )
        !
        n_active = 0
        do k = 1, n_full
            if( v_1(k) >= 0 ) then
                n_active = n_active + 1
            endif
        enddo
        !
        allocate( ind_active( n_active ) )
        !
        i = 0
        do k = 1, n_full
            if( v_1(k) >= 0 ) then
                i = i + 1
                ind_active(i) = k
            endif
        enddo
        !
        n_interior = 0
        do k = 1, n_full
            if( v_1(k) == 0 ) then
                n_interior = n_interior + 1
            endif
        enddo
        !
        allocate( v_2( n_active ) )
        v_2 = v_1( ind_active )
        !
        allocate( ind_interior( n_interior ) )
        !
        i = 0
        do k = 1, n_active
            if( v_2(k) == 0 ) then
                i = i + 1
                ind_interior(i) = k
            endif
        enddo
        !
        n_boundaries = 0
        do k = 1, n_active
            if( v_2(k) == 1 ) then
                n_boundaries = n_boundaries + 1
            endif
        enddo
        !
        allocate( ind_boundary( n_boundaries ) ) 
        !
        i = 0
        do k = 1, n_active
            if( v_2(k) == 1 ) then
                i = i + 1
                ind_boundary(i) = k
            endif
        enddo
        !
    end subroutine setIndexArrays_cScalar3D_MR
    !
    !> Creates standard(1-D array) for all sub_scalar,
    !> INCLUDING redundant interface nodes.
    !
    function getFullArray_cScalar3D_MR( self ) result( array )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( in ) :: self
        !
        complex( kind=prec ), allocatable, dimension(:) :: array
        !
        integer :: n, i, i1, i2
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "getFullArray_cScalar3D_MR > self not allocated." )
        endif
        !
        n = self%lengthFull()
        !
        allocate( array(n) )
        !
        array = C_ZERO
        !
        i1 = 1
        i2 = 0
        !
        do i = 1, self%grid%n_grids
            !
            n = self%sub_scalar(i)%length()
            !
            i2 = i2 + n
            array(i1:i2) = self%sub_scalar(i)%getArray()
            i1 = i1 + n
            !
        enddo
        !
    end function getFullArray_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setFullArray_cScalar3D_MR( self, array )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), dimension(:), intent( in ) :: array
        !
        integer :: i, i1, i2, n
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "setFullArray_cScalar3D_MR > self not allocated." )
        endif
        !
        i1 = 1
        i2 = 0
        !
        do i = 1, self%grid%n_grids
            !
            n = self%sub_scalar(i)%length()
            i2 = i2 + n
            call self%sub_scalar(i)%setArray( array(i1:i2) )
            i1 = i1 + n
            !
        enddo
        !
    end subroutine setFullArray_cScalar3D_MR
    !
    !> No function briefing
    !
    function lengthFull_cScalar3D_MR( self ) result( n )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( in ) :: self
        !
        integer :: i, n
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "lengthFull_cScalar3D_MR > self not allocated." )
        endif
        !
        n = 0
        !
        do i = 1, self%grid%n_grids
            !
            n = n + self%sub_scalar(i)%length()
            !
        enddo
        !
    end function lengthFull_cScalar3D_MR
    !
    !> No function briefing
    !
    function findFull_cScalar3D_MR( self, c ) result( I )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( in ) :: self
        real( kind=prec ), intent(in) :: c
        !
        integer, allocatable, dimension(:) :: I
        real( kind=prec ), dimension(:), allocatable :: v
        integer :: n, n_I, k
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "findFull_cScalar3D_MR > self not allocated." )
        endif
        !
        n = self%lengthFull()
        v = self%getFullArray()
        !
        n_I = 0
        do k = 1, n
            if(v(k) == c) n_I = n_I + 1
        enddo
        !
        allocate(I(n_I))
        !
        n_I = 0
        do k = 1, n
            if(v(k) == c) then
                n_I = n_I + 1
                I(n_I) = k
            endif
        enddo
        !
    end function findFull_cScalar3D_MR
    !
    !> toSG
    !
    !> input self is of class rScalar3D_MR , output SGscalar isi of class rScalar3D_SG
    !> this just copies contents of an MR cell into all subdividing fine grid cells
    !
    subroutine toSG_cScalar3D_MR( self, scalar_sg )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( in ) :: self
        type( cScalar3D_SG_t ), intent( out ) :: scalar_sg
        !
        integer :: i_grid, i, j, k, z, cs
        integer :: i1, i2, j1, j2, k1, k2
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "toSG_cScalar3D_MR > self not allocated." )
        endif
        !
        !> Using a temporary Grid SG with AirLayers, for instantiate the scalar_sg output
        scalar_sg = cScalar3D_SG_t( self%grid, self%grid_type )
        !
        select type( grid => self%grid )
            !
            class is( Grid3D_MR_t )
                !
                do i_grid = 1, grid%n_grids
                    !
                    !> vertical layers in fine grid
                    k1 = grid%coarseness( i_grid, 3 )
                    k2 = grid%coarseness( i_grid, 4 )
                    !
                    cs = 2 ** grid%coarseness( i_grid, 1 )
                    !
                    z = 1
                    !
                    do k = k1, k2
                        !
                        do i = 1, self%sub_scalar( i_grid )%nx
                            !
                            i1 = (i-1) * cs + 1
                            i2 = i * cs
                            !
                            do j = 1, self%sub_scalar( i_grid )%ny
                                !
                                j1 = (j-1)*cs+1
                                j2 = j * cs
                                !
                                scalar_sg%v( i1:i2, j1:j2, k ) = self%sub_scalar( i_grid )%v(i,j,z)
                                !
                            enddo
                        enddo
                        !
                        z = z + 1
                        !
                    enddo
                enddo
                !
            class default
                call errStop( "toSG_cScalar3D_MR > Unclassified grid" )
            !
        end select
        !
    end subroutine toSG_cScalar3D_MR
    !
    !> divFine
    !
    !> Rescaling each sub_scalar by dividing by dividing each
    !> coarsened cell by the number of fine-grid cells it contains.  This combined with MR2SG_T
    !> can be used to compute average (of fine grid scalar) on the MR grid cells.   It will also
    !> be needed in the PDEmapping routines (including transposes needed for inversion)
    !
    subroutine divFine_cScalar3D_MR( self )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        !
        integer :: i_grid
        real( kind=prec ) :: c
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "divFine_cScalar3D_MR > self not allocated." )
        endif
        !
        select type( grid => self%grid )
            !
            class is( Grid3D_MR_t )
                !
                do i_grid = 1, grid%n_grids
                    !
                    !> c*c is total number of fine grid cells in each cell in this subgrid
                    c = 2 ** grid%coarseness( i_grid, 1 )
                    !
                    c = 1. / ( c * c )
                    !
                    call self%sub_scalar( i_grid )%mult( c )
                    !
                enddo
                !
            class default
                call errStop( "divFine_cScalar3D_MR > Unclassified grid" )
            !
        end select
        !
    end subroutine divFine_cScalar3D_MR
    !
    !> fromSG
    !> Gary's implementation
    !
    !> this is adjoint (transpose) of MR2SG
    !
    !> self is of class cScalar3D_MR (output/modified), SGscalar (input, not modified)
    !> is of class cScalar3D_SG -- should be compatible
    !
    subroutine fromSG_cScalar3D_MR( self, scalar_sg )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        type( cScalar3D_SG_t ), intent( in ) :: scalar_sg
        !
        integer :: i_grid, i, j, k, cs, z
        integer :: i1, i2, j1, j2, k1, k2
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "fromSG_cScalar3D_MR > self not allocated." )
        endif
        !
        if( .NOT. scalar_sg%is_allocated ) then
            call errStop( "fromSG_cScalar3D_MR > scalar_sg not allocated" )
        endif
        !
        select case( self%grid_type )
            !
            case( CELL )
                !
                select type( grid => self%grid )
                    !
                    class is( Grid3D_MR_t )
                        !
                        do i_grid = 1, grid%n_grids
                            !
                            !> vertical layers in fine grid
                            k1 = grid%coarseness( i_grid, 3 )
                            k2 = grid%coarseness( i_grid, 4 )
                            !
                            cs = 2 ** grid%coarseness( i_grid, 1 )
                            !
                            z = 1
                            !
                            do k = k1, k2
                                !
                                do i = 1, self%sub_scalar( i_grid )%nx
                                    !
                                    i1 = (i-1) * cs + 1
                                    i2 = i * cs
                                    !
                                    do j = 1, self%sub_scalar( i_grid )%ny
                                        !
                                        j1 = (j-1) * cs + 1
                                        j2 = j * cs
                                        !
                                        self%sub_scalar( i_grid )%v(i,j,z) = sum( scalar_sg%v( i1:i2, j1:j2, k ) )
                                        !
                                    enddo
                                    !
                                enddo
                                !
                                z = z + 1
                                !
                            enddo
                        enddo
                        !
                    class default
                        call errStop( "fromSG_cScalar3D_MR > Unclassified grid" )
                    !
                end select
                !
            case default
                call errStop( "fromSG_cScalar3D_MR > Implemented just for type CELL" )
            !
        end select
        !
    end subroutine fromSG_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine cScalar3D_MR_dtor( self )
        implicit none
        !
        type( cScalar3D_MR_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor cScalar3D_MR"
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "cScalar3D_MR_dtor > self not allocated." )
        endif
        !
        call self%baseDealloc
        !
        if( allocated( self%sub_scalar ) ) deallocate( self%sub_scalar )
        !
        self%nx = 0
        self%ny = 0
        self%nz = 0
        !
        self%grid_type = ""
        self%is_allocated = .FALSE.
        !
    end subroutine cScalar3D_MR_dtor
    !
    !> No subroutine briefing
    !
    subroutine setAllBoundary_cScalar3D_MR( self, cvalue )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        call errStop( "setAllBoundary_cScalar3D_MR just implemented for SG!" )
        !
    end subroutine setAllBoundary_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setOneBoundary_cScalar3D_MR( self, bdry, cvalue, int_only )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        character(*), intent( in ) :: bdry
        complex( kind=prec ), intent( in ) :: cvalue
        logical, intent( in ), optional :: int_only
        !
        integer :: i_grid
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "setAllBoundary_cScalar3D_MR > self not allocated." )
        endif
        !
        do i_grid = 1, self%grid%n_grids
            !
            call self%sub_scalar( i_grid )%setOneBoundary( bdry, cvalue, int_only )
            !
        enddo
        !
    end subroutine setOneBoundary_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    function length_cScalar3D_MR( self ) result( field_length )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( in ) :: self
        !
        integer :: field_length
        !
        field_length = size( self%indActive() )
        !
    end function length_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine zeros_cScalar3D_MR( self )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        !
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "zeros_cScalar3D_MR > self not allocated." )
        endif
        !
        select type( grid => self%grid )
            !
            class is( Grid3D_MR_t )
                !
                do i = 1, grid%n_grids
                    call self%sub_scalar(i)%zeros
                enddo
                !
            class default
                call errStop( "zeros_cScalar3D_MR > Unclassified grid" )
            !
        end select
        !
    end subroutine zeros_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine conjugate_cScalar3D_MR( self )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        !
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "conjugate_cScalar3D_MR > self not allocated." )
        endif
        !
        do i = 1, self%grid%n_grids
            !
            call self%sub_scalar(i)%conjugate
            !
        enddo
        !
    end subroutine conjugate_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine add_cScalar3D_MR( self, rhs )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "add_cScalar3D_MR > self not allocated." )
        endif
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "add_cScalar3D_MR > rhs not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type( rhs )
                !
                class is( cScalar3D_MR_t )
                    !
                    do i = 1, self%grid%n_grids
                        !
                        self%sub_scalar(i)%v = self%sub_scalar(i)%v + rhs%sub_scalar(i)%v
                        !
                    enddo
                !
                class is( rScalar3D_MR_t )
                    !
                    do i = 1, self%grid%n_grids
                        !
                        self%sub_scalar(i)%v = self%sub_scalar(i)%v + rhs%sub_scalar(i)%v
                        !
                    enddo
                    !
                class default
                    call errStop( "add_cScalar3D_MR > rhs must be Scalar (try vec%scl)!" )
                    !
            end select
            !
        else
            call errStop( "add_cScalar3D_MR > Incompatible inputs." )
        endif
        !
    end subroutine add_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine linComb_cScalar3D_MR( self, rhs, c1, c2 )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        complex( kind=prec ), intent( in ) :: c1, c2
        !
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "linComb_cScalar3D_MR > self not allocated." )
        endif
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "linComb_cScalar3D_MR > rhs not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type( rhs )
                !
                class is( cScalar3D_MR_t )
                    !
                    do i = 1, self%grid%n_grids
                        !
                        self%sub_scalar(i)%v = c1 * self%sub_scalar(i)%v + c2 * rhs%sub_scalar(i)%v
                        !
                    enddo
                    !
                class is( rScalar3D_MR_t )
                    !
                    do i = 1, self%grid%n_grids
                        !
                        self%sub_scalar(i)%v = c1 * self%sub_scalar(i)%v + c2 * rhs%sub_scalar(i)%v
                        !
                    enddo
                    !
                class default
                    call errStop( "linComb_cScalar3D_MR > rhs must be Scalar (try vec%scl)!" )
                    !
            end select
            !
        else
            call errStop( "linComb_cScalar3D_MR > Incompatible inputs." )
        endif
        !
    end subroutine linComb_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine subValue_cScalar3D_MR( self, cvalue )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "subValue_cScalar3D_MR > self not allocated." )
        endif
        !
        do i = 1, self%grid%n_grids
            !
            self%sub_scalar(i)%v = self%sub_scalar(i)%v - cvalue
            !
        enddo
        !
    end subroutine subValue_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine subField_cScalar3D_MR( self, rhs )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "subField_cScalar3D_MR > self not allocated." )
        endif
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "subField_cScalar3D_MR > rhs not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type( rhs )
                !
                class is( cScalar3D_MR_t )
                    !
                    do i = 1, self%grid%n_grids
                        !
                        self%sub_scalar(i)%v = self%sub_scalar(i)%v - rhs%sub_scalar(i)%v
                        !
                    enddo
                    !
                class is( rScalar3D_MR_t )
                    !
                    do i = 1, self%grid%n_grids
                        !
                        self%sub_scalar(i)%v = self%sub_scalar(i)%v - rhs%sub_scalar(i)%v
                        !
                    enddo
                    !
                class default
                    call errStop( "subField_cScalar3D_MR > rhs must be Scalar (try vec%scl)!" )
                !
            end select
            !
        else
            call errStop( "subField_cScalar3D_MR > Incompatible inputs." )
        endif
        !
    end subroutine subField_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine multByReal_cScalar3D_MR( self, rvalue )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: rvalue
        !
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "multByReal_cScalar3D_MR > self not allocated." )
        endif
        !
        do i = 1, self%grid%n_grids
            !
            self%sub_scalar(i)%v = self%sub_scalar(i)%v * rvalue
            !
        enddo
        !
    end subroutine multByReal_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine multByComplex_cScalar3D_MR( self, cvalue )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "multByComplex_cScalar3D_MR > self not allocated." )
        endif
        !
        do i = 1, self%grid%n_grids
            !
            self%sub_scalar(i)%v = self%sub_scalar(i)%v * cvalue
            !
        enddo
        !
    end subroutine multByComplex_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine multByField_cScalar3D_MR( self, rhs )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "multByField_cScalar3D_MR > self not allocated." )
        endif
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "multByField_cScalar3D_MR > rhs not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type( rhs )
                !
                class is( cScalar3D_MR_t )
                    !
                    do i = 1, self%grid%n_grids
                        !
                        self%sub_scalar(i)%v = self%sub_scalar(i)%v * rhs%sub_scalar(i)%v
                        !
                    enddo
                    !
                class is( rScalar3D_MR_t )
                    !
                    do i = 1, self%grid%n_grids
                        !
                        self%sub_scalar(i)%v = self%sub_scalar(i)%v * rhs%sub_scalar(i)%v
                        !
                    enddo
                    !
                    ! MR FWD ENTER HERE !!!!
                    !call warning( "multByField_cScalar3D_MR > rhs is rScalar3D_MR_t" )
                    !
                class default
                    call errStop( "multByField_cScalar3D_MR > rhs must be Scalar (try vec%scl)!" )
                    !
            end select
            !
        else
            call errStop( "multByField_cScalar3D_MR > Incompatible inputs." )
        endif
        !
    end subroutine multByField_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine multAdd_cScalar3D_MR( self, cvalue, rhs )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        class( Field_t ), intent( in ) :: rhs
        !
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "multAdd_cScalar3D_MR > self not allocated." )
        endif
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "multAdd_cScalar3D_MR > rhs not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type( rhs )
                !
                class is( cScalar3D_MR_t )
                    !
                    do i = 1, self%grid%n_grids
                        !
                        self%sub_scalar(i)%v = self%sub_scalar(i)%v + cvalue * rhs%sub_scalar(i)%v
                        !
                    enddo
                !
                class is( rScalar3D_MR_t )
                !
                    do i = 1, self%grid%n_grids
                        !
                        self%sub_scalar(i)%v = self%sub_scalar(i)%v + cvalue * rhs%sub_scalar(i)%v
                        !
                    enddo
                    !
                class default
                    call errStop( "multAdd_cScalar3D_MR > rhs must be Scalar [try vec%mult(scl)]." )
                    !
            end select
            !
        else
            call errStop( "multAdd_cScalar3D_MR > Incompatible inputs." )
        endif
        !
    end subroutine multAdd_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    function dotProd_cScalar3D_MR( self, rhs ) result( cvalue )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( in ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        complex( kind=prec ) :: cvalue
        !
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "dotProd_cScalar3D_MR > self not allocated." )
        endif
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "dotProd_cScalar3D_MR > rhs not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            cvalue = C_ZERO
            !
            select type( rhs )
                !
                class is( cScalar3D_MR_t )
                    !
                    do i = 1, self%grid%n_grids
                        !
                        cvalue = cvalue + self%sub_scalar(i)%dotProd( rhs%sub_scalar(i) )
                        !
                    enddo
                    !
                class is( rScalar3D_MR_t )
                    !
                    do i = 1, self%grid%n_grids
                        !
                        cvalue = cvalue + self%sub_scalar(i)%dotProd( rhs%sub_scalar(i) )
                        !
                    enddo
                    !
                class default
                    call errStop( "dotProd_cScalar3D_MR > rhs must be Scalar (try vec%scl)!" )
                    !
            end select
            !
        else
            call errStop( "dotProd_cScalar3D_MR > Incompatible rhs" )
        endif
        !
    end function dotProd_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine divByValue_cScalar3D_MR( self, cvalue )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "divByValue_cScalar3D_MR > self not allocated." )
        endif
        !
        do i = 1, self%grid%n_grids
            !
            self%sub_scalar(i)%v = self%sub_scalar(i)%v / cvalue
            !
        enddo
        !
    end subroutine divByValue_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine divByField_cScalar3D_MR( self, rhs )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "divByField_cScalar3D_MR > self not allocated." )
        endif
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "divByField_cScalar3D_MR > rhs not allocated." )
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            select type( rhs )
                !
                class is( cScalar3D_MR_t )
                    !
                    do i = 1, self%grid%n_grids
                        !
                        self%sub_scalar(i)%v = self%sub_scalar(i)%v / rhs%sub_scalar(i)%v
                        !
                    enddo
                    !
                class is( rScalar3D_MR_t )
                    !
                    do i = 1, self%grid%n_grids
                        !
                        self%sub_scalar(i)%v = self%sub_scalar(i)%v / rhs%sub_scalar(i)%v
                        !
                    enddo
                    !
                class default
                    call errStop( "divByField_cScalar3D_MR > rhs must be Scalar (try vec%scl)!" )
                    !
            end select
            !
        else
            call errStop( "divByField_cScalar3D_MR > Incompatible inputs." )
        endif
        !
    end subroutine divByField_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine sumToNode_cScalar3D_MR( self, node_scalar, interior_only )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        class( Scalar_t ), intent( inout ) :: node_scalar
        logical, intent( in ), optional :: interior_only
        !
        integer :: v_xend, v_yend, v_zend
        logical :: is_interior_only
        integer :: i, nxF, nyF, nzF, nxC, nyC, nzC
        !
        if( .NOT. self%is_allocated ) then
             call errStop( "sumToNode_cScalar3D_MR > self not allocated." )
        endif
        !
        if( .NOT. node_scalar%is_allocated ) then
             call errStop( "sumToNode_cScalar3D_MR > node_scalar not allocated." )
        endif
        !
        select type( node_scalar )
            !
            class is( cScalar3D_MR_t )
                !
                select case( self%grid_type )
                    !
                    case( CELL )
                        !
                        !> set nodes for interior of all sub-scalars
                        do i = 1, self%grid%n_grids
                            !
                            call self%sub_scalar(i)%sumToNode( node_scalar%sub_scalar(i) )
                            !
                        enddo
                        !
                        !>set coarse grid nodes on interfaces
                        do i = 1, self%grid%n_grids-1
                            !
                            if( self%sub_scalar(i)%grid%nx .LT. self%sub_scalar(i+1)%grid%nx ) then
                                !
                                !> upper layer is coarser  -- fill in bottom level nodes
                                nxC = self%sub_scalar(i)%grid%nx
                                nyC = self%sub_scalar(i)%grid%ny
                                nzC = self%sub_scalar(i)%grid%nz
                                !
                                nxF = self%sub_scalar(i+1)%grid%nx
                                nyF = self%sub_scalar(i+1)%grid%ny
                                nzF = self%sub_scalar(i+1)%grid%nz
                                !
                                node_scalar%sub_scalar(i)%v( 2:nxC, 2:nyC, nzC+1 ) = &
                                  self%sub_scalar(i)%v( 1:nxC-1,   1:nyC-1,   nzC ) + &
                                  self%sub_scalar(i)%v( 2:nxC,     1:nyC-1,   nzC ) + &
                                  self%sub_scalar(i)%v( 1:nxC-1,   2:nyC,     nzC ) + &
                                  self%sub_scalar(i)%v( 2:nxC,     2:nyC,     nzC ) + &
                                self%sub_scalar(i+1)%v( 2:2:nxF-2, 2:2:nyF-2, 1   ) + &
                                self%sub_scalar(i+1)%v( 3:2:nxF-1, 2:2:nyF-2, 1   ) + &
                                self%sub_scalar(i+1)%v( 2:2:nxF-2, 3:2:nyF-1, 1   ) + &
                                self%sub_scalar(i+1)%v( 3:2:nxF-1, 3:2:nyF-1, 1   )
                                !
                            else
                                !
                                nxF = self%sub_scalar(i)%grid%nx
                                nyF = self%sub_scalar(i)%grid%ny
                                nzF = self%sub_scalar(i)%grid%nz
                                !
                                nxC = self%sub_scalar(i+1)%grid%nx
                                nyC = self%sub_scalar(i+1)%grid%ny
                                nzC = self%sub_scalar(i+1)%grid%nz
                                !
                                node_scalar%sub_scalar(i+1)%v( 2:nxC, 2:nyC, 1 ) = &
                                self%sub_scalar(i+1)%v( 1:nxC-1,   1:nyC-1,   1   ) + &
                                self%sub_scalar(i+1)%v( 2:nxC,     1:nyC-1,   1   ) + &
                                self%sub_scalar(i+1)%v( 1:nxC-1,   2:nyC,     1   ) + &
                                self%sub_scalar(i+1)%v( 2:nxC,     2:nyC,   1     ) + &
                                  self%sub_scalar(i)%v( 2:2:nxF-2, 2:2:nyF-2, nzF ) + &
                                  self%sub_scalar(i)%v( 3:2:nxF-1, 2:2:nyF-2, nzF ) + &
                                  self%sub_scalar(i)%v( 2:2:nxF-2, 3:2:nyF-1, nzF ) + &
                                  self%sub_scalar(i)%v( 3:2:nxF-1, 3:2:nyF-1, nzF )
                                !
                            endif
                            !
                        enddo
                        !
                    case default
                        call errStop( "sumToNode_cScalar3D_MR just for CELL type" )
                    !
                end select
                !
            class default
                call errStop( "sumToNode_cScalar3D_MR > Unclassified node_scalar" )
            !
        end select
        !
        ! select type( node_scalar )
            ! !
            ! class is( cScalar3D_MR_t )
                ! !
                ! select case( self%grid_type )
                    ! !
                    ! case( CELL )
                        ! !
                        ! do i = 1, self%grid%n_grids
                            ! !
                            ! is_interior_only = .FALSE.
                            ! !
                            ! if( present( interior_only ) ) is_interior_only = interior_only
                            ! !
                            ! if( is_interior_only ) then
                                ! call self%sub_scalar(i)%setAllBoundary( C_ZERO )
                            ! endif
                            ! !
                            ! v_xend = size( self%sub_scalar(i)%v, 1 )
                            ! v_yend = size( self%sub_scalar(i)%v, 2 )
                            ! v_zend = size( self%sub_scalar(i)%v, 3 )
                            ! !
                            ! !> Interior
                            ! node_scalar%sub_scalar(i)%v( 2:v_xend, 2:v_yend, 2:v_zend ) = &
                            ! self%sub_scalar(i)%v( 1:v_xend-1, 1:v_yend-1, 1:v_zend-1 ) + &
                            ! self%sub_scalar(i)%v( 2:v_xend  , 1:v_yend-1, 1:v_zend-1 ) + &
                            ! self%sub_scalar(i)%v( 1:v_xend-1, 2:v_yend  , 1:v_zend-1 ) + &
                            ! self%sub_scalar(i)%v( 1:v_xend-1, 1:v_yend-1, 2:v_zend   ) + &
                            ! self%sub_scalar(i)%v( 2:v_xend  , 2:v_yend  , 1:v_zend-1 ) + &
                            ! self%sub_scalar(i)%v( 2:v_xend  , 1:v_yend-1, 2:v_zend   ) + &
                            ! self%sub_scalar(i)%v( 1:v_xend-1, 2:v_yend  , 2:v_zend   ) + &
                            ! self%sub_scalar(i)%v( 2:v_xend  , 2:v_yend  , 2:v_zend   )
                            ! !
                        ! enddo
                        ! !
                    ! case default
                        ! call errStop( "sumToNode_cScalar3D_MR just for CELL type" )
                ! end select
                ! !
            ! class default
                ! call errStop( "sumToNode_cScalar3D_MR > Unclassified node_scalar" )
            ! !
        ! end select
        ! !
    end subroutine sumToNode_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setSV_cScalar3D_MR( self, s_v )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), dimension(:), intent( in ) :: s_v
        !
        call errStop( "setSV_cScalar3D_MR not implemented!" )
        ! !
        ! if( .NOT. self%is_allocated ) then
            ! call errStop( "setSV_cScalar3D_MR > self not allocated." )
        ! endif
        ! !
        ! call self%switchStoreState( singleton )
        ! !
        ! if( allocated( self%v ) ) deallocate( self%v )
        ! !
        ! self%s_v = s_v
        ! !
    end subroutine setSV_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    function getArray_cScalar3D_MR( self ) result( array )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( in ) :: self
        complex( kind=prec ), allocatable, dimension(:) :: array
        !
        complex( kind=prec ), allocatable, dimension(:) :: v_full
        !
        if( .NOT. self%is_allocated ) then
             call errStop( "getArray_cScalar3D_MR > self not allocated." )
        endif
        !
        v_full = self%getFullArray()
        !
        array = v_full( self%indActive() )
        !
        deallocate( v_full )
        !
    end function getArray_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setArray_cScalar3D_MR( self, array )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), dimension(:), intent( in ) :: array
        !
        complex( kind=prec ), allocatable, dimension(:) :: vFull
        !
        if( .NOT. self%is_allocated ) then
             call errStop( "setArray_cScalar3D_MR > self not allocated." )
        endif
        !
        allocate( vFull( self%lengthFull() ) )
        !
        vFull = C_ZERO
        !
        vFull( self%indActive() ) = array
        !
        call self%setFullArray( vFull )
        !
        deallocate( vFull )
        !
    end subroutine setArray_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine copyFrom_cScalar3D_MR( self, rhs )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "copyFrom_cScalar3D_MR > rhs not allocated" )
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
            class is( cScalar3D_MR_t )
                !
                if( allocated( rhs%sub_scalar ) ) then
                    !
                    self%sub_scalar = rhs%sub_scalar
                    !
                else
                    call errStop( "copyFrom_cScalar3D_MR > rhs%sub_scalar not allocated" )
                endif
                !
                self%is_allocated = .TRUE.
                !
            class default
                call errStop( "copyFrom_cScalar3D_MR > Unclassified rhs" )
            !
        end select
        !
    end subroutine copyFrom_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine read_cScalar3D_MR( self, funit, ftype )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        integer, intent( in ) :: funit
        character(:), allocatable, intent( in ), optional :: ftype
        !
        call errStop( "read_cScalar3D_MR not implemented!" )
        !
    end subroutine read_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine write_cScalar3D_MR( self, funit, ftype )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( inout ) :: self
        integer, intent( in ) :: funit
        character(:), allocatable, intent( in ), optional :: ftype
        !
        integer :: i_grid
        !
        !    header for base fine grid
        write( funit )  self%nx, self%ny, self%nz, self%grid%n_grids
        !
        do i_grid = 1, self%grid%n_grids
            !
            !   then write out each sub_vector
            call self%sub_scalar(i_grid)%write(funit)
            !
        enddo
        !
    end subroutine write_cScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine print_cScalar3D_MR( self, io_unit, title, append )
        implicit none
        !
        class( cScalar3D_MR_t ), intent( in ) :: self
        integer, intent( in ), optional :: io_unit
        character(*), intent( in ), optional :: title
        logical, intent( in ), optional :: append
        !
        integer :: i, ix, iy, iz, funit
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "print_cScalar3D_MR > self not allocated." )
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
        write(funit,*) "cScalar3D_MR field"
        do i = 1, self%grid%n_grids
            !
            write( funit, * ) self%sub_scalar(i)%nx, self%sub_scalar(i)%ny, self%sub_scalar(i)%nz
            !
            write(funit,*) "sub_scalar ", i
            do ix = 1, self%sub_scalar(i)%nx
                do iy = 1, self%sub_scalar(i)%ny
                    do iz = 1, self%sub_scalar(i)%nz
                        if( self%sub_scalar(i)%v( ix, iy, iz ) /= 0 ) then
                            write(funit,*) ix,iy,iz, ":[", self%sub_scalar(i)%v( ix, iy, iz ), "]"
                        endif
                    enddo
                enddo
            enddo
            !
        enddo
        !
    end subroutine print_cScalar3D_MR
    !
end module cScalar3D_MR
