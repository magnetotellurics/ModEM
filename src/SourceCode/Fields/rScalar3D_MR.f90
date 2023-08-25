!
!> Derived class to define a rScalar3D_MR 
!
module rScalar3D_MR
    !
    use rScalar3D_SG
    use Grid3D_MR
    !
    type, extends( Scalar_t ) :: rScalar3D_MR_t
        !
        type( rScalar3D_SG_t ), allocatable, dimension(:) :: sub_scalar
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
            procedure, public :: getArray => getArray_rScalar3D_MR
            procedure, public :: setArray => setArray_rScalar3D_MR
            !
            procedure, public :: getFullArray => getFullArray_rScalar3D_MR
            procedure, public :: setFullArray => setFullArray_rScalar3D_MR
            !
            procedure, public :: lengthFull => lengthFull_rScalar3D_MR
            procedure, public :: findFull => findFull_rScalar3D_MR
            !
            procedure, public :: findValue => findValue_rScalar3D_MR
            !
            procedure, public :: MRtoSG => MRtoSG_rScalar3D_MR
            procedure, public :: SGtoMR => SGtoMR_rScalar3D_MR
            procedure, public :: divFine => divFine_rScalar3D_MR
            !
            !> Boundary operations
            procedure, public :: setAllBoundary => setAllBoundary_rScalar3D_MR
            procedure, public :: setOneBoundary => setOneBoundary_rScalar3D_MR
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
            procedure, public :: toNode => toNode_rScalar3D_MR
            !
            !> Miscellaneous
            procedure, public :: getV => getV_rScalar3D_MR
            procedure, public :: setV => setV_rScalar3D_MR
            !
            procedure, public :: getSV => getSV_rScalar3D_MR
            procedure, public :: setSV => setSV_rScalar3D_MR
            !
            procedure, public :: deallOtherState => deallOtherState_rScalar3D_MR
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
    function rScalar3D_MR_ctor_copy( E_in ) result( self )
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
        select type( grid => self%grid )
            !
            class is( Grid3D_MR_t )
                !
                do i = 1, grid%n_grids
                    self%sub_scalar(i) = E_in%sub_scalar(i)
                enddo
                !
            class default
                call errStop( "rScalar3D_MR_ctor_copy > Unclassified grid" )
            !
        end select
        !
    end function rScalar3D_MR_ctor_copy
    !
    !> No function briefing
    !
    function rScalar3D_MR_ctor_default( grid, grid_type ) result( self )
        implicit none
        !
        class( Grid_t ), target, intent( in ) :: grid
        character( len=4 ), intent( in ) :: grid_type
        !
        type( rScalar3D_MR_t ) :: self
        !
        integer :: i, nx, ny, nz, nzAir, nz_earth, alloc_stat
        !
        self%grid => grid
        self%grid_type = grid_type
        !
        write(*,*) "grid_type: ", grid_type
        !
        !> Grid dimensions
        call self%grid%getDimensions( nx, ny, nz, nzAir )
        nz_earth = nz - nzAir
        !
        self%nx = nx
        self%ny = ny
        self%nz = nz
        !
        call self%initializeSub
        !
        if( self%is_allocated ) then
            !
            self%Nxyz = product( self%NdV )
            !
        else
            call errStop( "rScalar3D_MR_ctor_default > Unable to allocate scalar." )
        endif
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
        integer :: i, nx, ny, nz, nzAir, nz_earth, alloc_stat
        !
        !> Grid dimensions
        call self%grid%getDimensions( nx, ny, nz, nzAir )
        nz_earth = nz - nzAir
        !
        self%nx = nx
        self%ny = ny
        self%nz = nz
        !
        select type( grid => self%grid )
            !
            class is( Grid3D_MR_t )
                !
                self%is_allocated = .TRUE.
                allocate( self%sub_scalar( grid%n_grids ), stat = alloc_stat )
                self%is_allocated = self%is_allocated .AND.( alloc_stat .EQ. 0 )
                !
                do i = 1, self%grid%getNGrids()
                    !
                    self%sub_scalar(i) = rScalar3D_SG_t( grid%sub_grids(i), self%grid_type )
                    !
                    write( *, * ) "SubScalar", i, "-nx=", self%sub_scalar(i)%nx, ", ny=", self%sub_scalar(i)%ny, "nz=", self%sub_scalar(i)%nz
                    !
                enddo
                !
                write( *, * ) "MainScalar-nx=", self%nx, ", ny=", self%ny, "nz=", self%nz, self%nx*self%ny*self%nz
                !
            class default
                call errStop( "initializeSub_rScalar3D_MR > Unclassified grid" )
            !
        end select
        !
    end subroutine initializeSub_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setIndexArrays_rScalar3D_MR( self, ind_boundary, ind_interior, ind_active, xy_in ) 
        implicit none
        !
        class( rScalar3D_MR_t ), intent( in ) :: self
        integer, dimension(:), allocatable, intent( out ) :: ind_boundary, ind_interior
        integer, dimension(:), allocatable, intent( out ), optional :: ind_active
        logical, intent( in ), optional :: xy_in
        !
        type( rScalar3D_MR_t ) :: temp_scalar
        logical :: xy, int_only
        integer :: i, k
        integer :: n_full, n_active, n_interior, n_boundaries
        real( kind=prec ), dimension(:), allocatable :: v_1, v_2
        !
        if( .NOT. present( xy_in ) ) then
            xy = .FALSE.
        else
            xy = xy_in
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
                    case( CELL, CELL_EARTH )
                        int_only = .FALSE.
                    case( NODE )
                        int_only = .TRUE.
                    case default
                        !
                        call errStop( "setIndexArrays_rScalar3D_MR > Invalid grid type option!" )
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
                call errStop( "setIndexArrays_rScalar3D_MR > Unclassified grid" )
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
    end subroutine setIndexArrays_rScalar3D_MR
    !
    !> Creates standard(1-D array) for all sub_scalar,
    !> INCLUDING redundant interface nodes.
    !
    function getFullArray_rScalar3D_MR( self ) result( array )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( in ) :: self
        !
        complex( kind=prec ), allocatable, dimension(:) :: array
        !
        integer :: n, i, i1, i2
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
        do i = 1, self%grid%getNGrids()
            !
            n = self%sub_scalar(i)%length()
            !
            i2 = i2 + n
            array(i1:i2) = self%sub_scalar(i)%getArray()
            i1 = i1 + n
            !
        enddo
        !
    end function getFullArray_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setFullArray_rScalar3D_MR( self, array )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), dimension(:), intent( in ) :: array
        !
        integer :: i, i1, i2, n
        !
        i1 = 1
        i2 = 0
        !
        do i = 1, self%grid%getNGrids()
            !
            n = self%sub_scalar(i)%length()
            i2 = i2 + n
            call self%sub_scalar(i)%setArray( array(i1:i2) )
            i1 = i1 + n
            !
        enddo
        !
    end subroutine setFullArray_rScalar3D_MR
    !
    !> No function briefing
    !
    function lengthFull_rScalar3D_MR( self ) result( n )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( in ) :: self
        !
        integer :: i, n
        !
        n = 0
        !
        do i = 1, self%grid%getNGrids()
            !
            n = n + self%sub_scalar(i)%length()
            !
        enddo
        !
    end function lengthFull_rScalar3D_MR
    !
    !> No function briefing
    !
    function findFull_rScalar3D_MR( self, c ) result( I )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( in ) :: self
        real( kind=prec ), intent(in) :: c
        !
        integer, dimension(:), allocatable :: I
        real( kind=prec ), dimension(:), allocatable :: v
        integer :: n, n_I, k
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
    end function findFull_rScalar3D_MR
    !
    !> No function briefing
    !
    function findValue_rScalar3D_MR(self, c) result( I )
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
        allocate( v(n) )
        v = self%getArray()
        !
        n_I = 0
        do k = 1, n
            if( v(k) == c ) n_I = n_I + 1
        enddo
        !
        allocate( I(n_I) )
        !
        n_I = 0
        do k = 1, n
            if( v(k) == c ) then
                n_I = n_I + 1
                I(n_I) = k
            endif
        enddo
        !
    end function findValue_rScalar3D_MR
    !
    !> MRtoSG
    !
    !> input self is of class rScalar3D_MR , output SGscalar isi of class rScalar3D_SG
    !> this just copies contents of an MR cell into all subdividing fine grid cells
    !
    subroutine MRtoSG_rScalar3D_MR( self, scalar_sg )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( in ) :: self
        type( rScalar3D_SG_t ), intent( out ) :: scalar_sg
        !
        integer :: i_grid, i, j, k, cs
        integer :: i1, i2, j1, j2, k1, k2
        !
        scalar_sg = rScalar3D_SG_t( self%grid, NODE )
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
                    cs = 2 ** grid%coarseness( i_grid, 1 )
                    !
                    do k = k1, k2
                        !
                        do i = 1, self%sub_scalar( i_grid )%nx
                            !
                            i1 = (i-1) * cs + 1
                            i2 = i1 + cs
                            !
                            do j = 1, self%sub_scalar( i_grid )%ny
                                !
                                j1 = (j-1)*cs+1
                                j2 = j1+cs
                                !
                                scalar_sg%v( i1:i2, j1:j2, k ) = self%sub_scalar( i_grid )%v(i,j,k)
                                !
                            enddo
                        enddo
                    enddo
                enddo
                !
            class default
                call errStop( "MRtoSG_rScalar3D_MR > Unclassified grid" )
            !
        end select
        !
    end subroutine MRtoSG_rScalar3D_MR
    !
    !> divFine
    !
    !> Rescaling each sub_scalar by dividing by dividing each
    !> coarsened cell by the number of fine-grid cells it contains.  This combined with MR2SG_T
    !> can be used to compute average (of fine grid scalar) on the MR grid cells.   It will also
    !> be needed in the PDEmapping routines (including transposes needed for inversion)
    !
    subroutine divFine_rScalar3D_MR( self )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        !
        integer :: i_grid
        real( kind=prec ) :: c
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
                call errStop( "divFine_rScalar3D_MR > Unclassified grid" )
            !
        end select
        !
    end subroutine divFine_rScalar3D_MR
    !
    !> SGtoMR
    !> Gary's implementation
    !
    !> this is adjoint (transpose) of MR2SG
    !
    !> self is of class rScalar3D_MR (output/modified), SGscalar (input, not modified)
    !> is of class rScalar3D_SG -- should be compatible
    !
    subroutine SGtoMR_rScalar3D_MR( self, scalar_sg )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        type( rScalar3D_SG_t ), intent( in ) :: scalar_sg
        !
        integer :: i_grid, i, j, k, cs
        integer :: i1, i2, j1, j2, k1, k2
        !
        select type( grid => scalar_sg%grid )
            !
            class is( Grid3D_MR_t )
                !
                do i_grid = 1,self%grid%n_grids
                    !
                    !> vertical layers in fine grid
                    k1 = grid%coarseness( i_grid, 3 )
                    k2 = grid%coarseness( i_grid, 4 )
                    cs = 2 ** grid%coarseness( i_grid, 1 )
                    !
                    do k = k1, k2
                        !
                        do i = 1, self%sub_scalar( i_grid )%nx
                            !
                            i1 = (i-1) * cs + 1
                            i2 = i1 + cs
                            !
                            do j = 1, self%sub_scalar( i_grid )%ny
                                !
                                j1 = (j-1) * cs + 1
                                j2 = j1 + cs
                                !  not sure sum-sum works here, but guess it should
                                !    Note that only this line changes from MR2SG
                                self%sub_scalar( i_grid )%v(i,j,k) = sum( scalar_sg%v( i1:i2, j1:j2, k ) )
                                !
                            enddo
                        enddo
                    enddo
                enddo
                !
            class default
                call errStop( "SGtoMR_rScalar3D_MR > Unclassified grid" )
            !
        end select
        !
    end subroutine SGtoMR_rScalar3D_MR
    !
    !> SGtoMR
    !> Williams implementation
    ! !
    ! subroutine SGtoMR_rScalar3D_MR( self, scalar_sg )
        ! implicit none
        ! !
        ! class( rScalar3D_MR_t ), intent( inout ) :: self
        ! type( rScalar3D_SG_t ), intent( in ) :: scalar_sg
        ! !
        ! class( Grid_t ), pointer :: grid
        ! !
        ! integer :: x_nx, x_ny, x_nz
        ! integer :: last, Cs, i1, i2, i, k
        ! !
        ! select type( grid => scalar_sg%grid )
            ! !
            ! class is( Grid3D_MR_t )
                ! !
                ! x_nx = size( scalar_sg%v, 1 )
                ! x_ny = size( scalar_sg%v, 2 )
                ! x_nz = size( scalar_sg%v, 3 )
                ! !
                ! do k = 1, grid%n_grids
                    ! !
                    ! Cs = 2**grid%coarseness(k, 1)
                    ! i1 = grid%coarseness(k, 3)
                    ! i2 = grid%coarseness(k, 4)
                    ! !
                    ! do i = 1, Cs
                        ! !
                        ! last = size( grid%Dx )
                        ! !
                        ! self%sub_scalar(k)%v = self%sub_scalar(k)%v + &
                        ! scalar_sg%v( i:x_nx:Cs, 1:x_ny:Cs, i1:i2+1 ) * &
                        ! repMat( grid%Dx(i:last:Cs), 1, &
                        ! grid%sub_grids(k)%Ny + 1, &
                        ! grid%sub_grids(k)%Nz + 1, .FALSE. )
                        ! !
                    ! enddo
                    ! !
                    ! self%sub_scalar(k)%v = self%sub_scalar(k)%v / &
                    ! repMat(grid%sub_grids(k)%Dx, &
                    ! 1, &
                    ! grid%sub_grids(k)%Ny + 1, &
                    ! grid%sub_grids(k)%Nz + 1, .FALSE. )
                    ! !
                ! enddo
                ! !
            ! class default
                ! call errStop( "SGtoMR_rScalar3D_MR > Unclassified grid" )
            ! !
        ! end select
        ! !
    ! end subroutine SGtoMR_rScalar3D_MR
    ! !
    !> No subroutine briefing
    !
    subroutine rScalar3D_MR_dtor( self )
        implicit none
        !
        type( rScalar3D_MR_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor rScalar3D_MR"
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "rScalar3D_MR_dtor > self not allocated." )
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
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "setAllBoundary_rScalar3D_MR > self not allocated." )
        endif
        !
        call self%switchStoreState( compound )
        !
        do i = 1, self%grid%getNGrids()
            !
            select case( self%grid_type )
                !
                case( NODE, CELL, CELL_EARTH ) 
                    !
                    self%sub_scalar(i)%v((/1, self%NdV(1)/), :, :) = cvalue
                    self%sub_scalar(i)%v(:, (/1, self%NdV(2)/), :) = cvalue
                    self%sub_scalar(i)%v(:, :, (/1, self%NdV(3)/)) = cvalue
                    !
                case default
                    call errStop( "setAllBoundary_rScalar3D_MR > grid_type ["//self%grid_type//"] not recognized." )
                !
            end select
            !
        enddo
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
        integer :: i
        logical :: int_only_p
        !
        call self%switchStoreState( compound )
        !
        do i = 1, self%grid%getNGrids()
            !
            if( .NOT. present(int_only)) then
                 int_only_p = .FALSE.
            else 
                 int_only_p = int_only
            endif
            !
            select case( self%sub_scalar(i)%grid_type )
                !
                case( NODE )
                    !
                    if( int_only_p ) then
                        !
                        select case( bdry )
                            !
                            case("x1")
                                self%sub_scalar(i)%v(1, 2:self%sub_scalar(i)%NdV(2)-1, 2:self%sub_scalar(i)%NdV(3)-1) = real( cvalue, kind=prec ) 
                            case("x2")
                                self%sub_scalar(i)%v(self%sub_scalar(i)%NdV(1), 2:self%sub_scalar(i)%NdV(2)-1, 2:self%sub_scalar(i)%NdV(3)-1) = real( cvalue, kind=prec )
                            case("y1")
                                self%sub_scalar(i)%v(2:self%sub_scalar(i)%NdV(1)-1, 1, 2:self%sub_scalar(i)%NdV(3)-1) = real( cvalue, kind=prec )
                            case("y2")
                                self%sub_scalar(i)%v(2:self%sub_scalar(i)%NdV(1)-1, self%sub_scalar(i)%NdV(2), 2:self%sub_scalar(i)%NdV(3)-1) = real( cvalue, kind=prec )
                            case("z1")
                                self%sub_scalar(i)%v(2:self%sub_scalar(i)%NdV(1)-1, 2:self%sub_scalar(i)%NdV(2)-1, 1) = real( cvalue, kind=prec )
                            case("z2")
                                self%sub_scalar(i)%v(2:self%sub_scalar(i)%NdV(1)-1, 2:self%sub_scalar(i)%NdV(2)-1, self%sub_scalar(i)%NdV(3)) = real( cvalue, kind=prec )
                            !
                        end select
                        !
                    else
                        !
                        select case( bdry )
                            !
                            case("x1")
                                self%sub_scalar(i)%v(1, :, :) = real( cvalue, kind=prec )
                            case("x2")
                                self%sub_scalar(i)%v(self%sub_scalar(i)%NdV(1), :, :) = real( cvalue, kind=prec )
                            case("y1")
                                self%sub_scalar(i)%v(:, 1, :) = real( cvalue, kind=prec )
                            case("y2")
                                self%sub_scalar(i)%v(:, self%sub_scalar(i)%NdV(2), :) = real( cvalue, kind=prec )
                            case("z1")
                                self%sub_scalar(i)%v(:, :, 1) = real( cvalue, kind=prec )
                            case("z2")
                                self%sub_scalar(i)%v(:, :, self%sub_scalar(i)%NdV(3)) = real( cvalue, kind=prec )
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
                            self%sub_scalar(i)%v(1, :, :) = real( cvalue, kind=prec )
                        case("x2")
                            self%sub_scalar(i)%v(self%sub_scalar(i)%NdV(1), :, :) = real( cvalue, kind=prec )
                        case("y1")
                            self%sub_scalar(i)%v(:, 1, :) = real( cvalue, kind=prec )
                        case("y2")
                            self%sub_scalar(i)%v(:, self%sub_scalar(i)%NdV(2), :) = real( cvalue, kind=prec )
                        case("z1")
                            self%sub_scalar(i)%v(:, :, 1) = real( cvalue, kind=prec )
                        case("z2")
                            self%sub_scalar(i)%v(:, :, self%sub_scalar(i)%NdV(3)) = real( cvalue, kind=prec )
                        !
                    end select
                    !
                case default
                    call errStop( "setOneBoundary_rScalar3D_MR > Invalid grid type" )
                !
            end select
            !
        enddo
        !
    end subroutine setOneBoundary_rScalar3D_MR
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
        field_length = size( self%indActive() )
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
        integer :: i
        integer :: x1, x2
        integer :: y1, y2
        integer :: z1, z2
        !
        call self%switchStoreState( compound )
        !
        do i = 1, self%grid%getNGrids()
            !
            x1 = xmin; x2 = xmax
            y1 = ymin; y2 = ymax
            z1 = zmin; z2 = zmax
            !
            if( xmin == 0) x1 = self%sub_scalar(i)%NdV(1)
            if( xmax <= 0) x2 = self%sub_scalar(i)%NdV(1) + xmax
            !
            if( ymin == 0) y1 = self%sub_scalar(i)%NdV(2)
            if( ymax <= 0) y2 = self%sub_scalar(i)%NdV(2) + ymax
            !
            if( zmin == 0) z1 = self%sub_scalar(i)%NdV(3)
            if( zmax <= 0) z2 = self%sub_scalar(i)%NdV(3) + zmax
            !
            self%sub_scalar(i)%v( x1:x2:xstep, y1:y2:ystep, z1:z2:zstep ) = rvalue
            !
        enddo
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
                    call self%sub_scalar(i)%zeros()
                enddo
                !
            class default
                call errStop( "zeros_rScalar3D_MR > Unclassified grid" )
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
        call errStop( "conjugate_rScalar3D_MR: Do not try to conjugate a real scalar!" )
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
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "add_rScalar3D_MR > self not allocated." )
        endif
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "add_rScalar3D_MR > rhs not allocated." )
        endif
        !
        call self%switchStoreState( rhs%store_state )
        !
        do i = 1, self%grid%getNGrids()
            !
            if( self%isCompatible( rhs ) ) then
                !
                select type( rhs )
                    !
                    class is( rScalar3D_MR_t )
                        !
                        if( rhs%store_state .EQ. compound ) then
                            !
                            self%sub_scalar(i)%v = self%sub_scalar(i)%v + rhs%sub_scalar(i)%v
                            !
                        elseif( rhs%store_state .EQ. singleton ) then
                            !
                            self%sub_scalar(i)%s_v = self%sub_scalar(i)%s_v + rhs%sub_scalar(i)%s_v
                            !
                        else
                            call errStop( "add_rScalar3D_MR > Unknown rhs store_state!" )
                        endif
                        !
                    class default
                        call errStop( "add_rScalar3D_MR > rhs must be Scalar (try vec%scl)!" )
                    !
                end select
                !
            else
                call errStop( "add_rScalar3D_MR > Incompatible inputs." )
            endif
            !
        enddo
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
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "linComb_rScalar3D_MR > self not allocated." )
        endif
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "linComb_rScalar3D_MR > rhs not allocated." )
        endif
        !
        call self%switchStoreState( rhs%store_state )
        !
        do i = 1, self%grid%getNGrids()
            !
            if( self%isCompatible( rhs ) ) then
                !
                select type( rhs )
                    !
                    class is( rScalar3D_MR_t )
                        !
                        if( rhs%store_state .EQ. compound ) then
                            !
                            self%sub_scalar(i)%v = c1 * self%sub_scalar(i)%v + c2 * rhs%sub_scalar(i)%v
                            !
                        elseif( rhs%store_state .EQ. singleton ) then
                            !
                            self%sub_scalar(i)%s_v = c1 * self%sub_scalar(i)%s_v + c2 * rhs%sub_scalar(i)%s_v
                            !
                        else
                            call errStop( "linComb_rScalar3D_MR > Unknown rhs store_state!" )
                        endif
                        !
                    class default
                        call errStop( "linComb_rScalar3D_MR > rhs must be Scalar (try vec%scl)!" )
                    !
                end select
                !
            else
                call errStop( "linComb_rScalar3D_MR > Incompatible inputs." )
            endif
            !
        enddo
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
        integer :: i
        !
        do i = 1, self%grid%getNGrids()
            !
            if( self%store_state .EQ. compound ) then
                !
                self%sub_scalar(i)%v = self%sub_scalar(i)%v - cvalue
                !
            elseif( self%store_state .EQ. singleton ) then
                !
                self%sub_scalar(i)%s_v = self%sub_scalar(i)%s_v - cvalue
                !
            else
                call errStop( "subValue_rScalar3D_MR > Unknown store_state!" )
            endif
            !
        enddo
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
        integer :: i
        !
        call self%switchStoreState( rhs%store_state )
        !
        do i = 1, self%grid%getNGrids()
            !
            if( self%isCompatible( rhs ) ) then
                !
                select type( rhs )
                    !
                    class is( rScalar3D_MR_t )
                        !
                        if( rhs%store_state .EQ. compound ) then
                            !
                            self%sub_scalar(i)%v = self%sub_scalar(i)%v - rhs%sub_scalar(i)%v
                            !
                        elseif( rhs%store_state .EQ. singleton ) then
                            !
                            self%sub_scalar(i)%s_v = self%sub_scalar(i)%s_v - rhs%sub_scalar(i)%s_v
                            !
                        else
                            call errStop( "subField_rScalar3D_MR > Unknown rhs store_state!" )
                        endif
                        !
                    class default
                        call errStop( "subField_rScalar3D_MR > rhs must be Scalar (try vec%scl)!" )
                    !
                end select
                !
            else
                call errStop( "subField_rScalar3D_MR > Incompatible inputs." )
            endif
            !
        enddo
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
        integer :: i
        !
        do i = 1, self%grid%getNGrids()
            !
            if( self%store_state .EQ. compound ) then
                !
                self%sub_scalar(i)%v = self%sub_scalar(i)%v * rvalue
                !
            elseif( self%store_state .EQ. singleton ) then
                !
                self%sub_scalar(i)%s_v = self%sub_scalar(i)%s_v * rvalue
                !
            else
                call errStop( "multByReal_rScalar3D_MR > Unknown store_state!" )
            endif
            !
        enddo
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
        integer :: i
        !
        do i = 1, self%grid%getNGrids()
            !
            if( self%store_state .EQ. compound ) then
                !
                self%sub_scalar(i)%v = self%sub_scalar(i)%v * cvalue
                !
            elseif( self%store_state .EQ. singleton ) then
                !
                self%sub_scalar(i)%s_v = self%sub_scalar(i)%s_v * cvalue
                !
            else
                call errStop( "multByComplex_rScalar3D_MR > Unknown store_state!" )
            endif
            !
        enddo
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
        integer :: i
        !
        call self%switchStoreState( rhs%store_state )
        !
        do i = 1, self%grid%getNGrids()
            !
            if( self%isCompatible( rhs ) ) then
                !
                select type( rhs )
                    !
                    class is( rScalar3D_MR_t )
                        !
                        if( rhs%store_state .EQ. compound ) then
                            !
                            self%sub_scalar(i)%v = self%sub_scalar(i)%v * rhs%sub_scalar(i)%v
                            !
                        elseif( rhs%store_state .EQ. singleton ) then
                            !
                            self%sub_scalar(i)%s_v = self%sub_scalar(i)%s_v * rhs%sub_scalar(i)%s_v
                            !
                        else
                            call errStop( "multByField_rScalar3D_MR > Unknown rhs store_state!" )
                        endif
                        !
                    class default
                        call errStop( "multByField_rScalar3D_MR > rhs must be Scalar (try vec%scl)!" )
                    !
                end select
                !
            else
                call errStop( "multByField_rScalar3D_MR > Incompatible inputs." )
            endif
            !
        enddo
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
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "multAdd_rScalar3D_MR > self not allocated." )
        endif
        !
        call self%switchStoreState( rhs%store_state )
        !
        do i = 1, self%grid%getNGrids()
            !
            if( self%isCompatible( rhs ) ) then
                !
                select type( rhs )
                    !
                    class is( rScalar3D_MR_t )
                        !
                        if( rhs%sub_scalar(i)%store_state .EQ. compound ) then
                            !
                            self%sub_scalar(i)%v = self%sub_scalar(i)%v + cvalue * rhs%sub_scalar(i)%getV()
                            !
                        elseif( rhs%sub_scalar(i)%store_state .EQ. singleton ) then
                            !
                            self%sub_scalar(i)%s_v = self%sub_scalar(i)%s_v + cvalue * rhs%sub_scalar(i)%getSV()
                            !
                        else
                            call errStop( "multAdd_rScalar3D_MR > Unknown rhs store_state!" )
                        endif
                        !
                    class default
                        call errStop( "multAdd_rScalar3D_MR > rhs must be Scalar [try vec%mult(scl)]." )
                    !
                end select
                !
            else
                call errStop( "multAdd_rScalar3D_MR > Incompatible inputs." )
            endif
            !
        enddo
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
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "dotProd_rScalar3D_MR > self not allocated." )
        endif
        !
        cvalue = C_ZERO
        !
        do i = 1, self%grid%getNGrids()
            !
            if( self%isCompatible( rhs ) ) then
                !
                select type( rhs )
                    !
                    class is( rScalar3D_MR_t )
                        !
                        cvalue = cvalue + self%sub_scalar(i)%dotProd( rhs%sub_scalar(i) )
                        !
                    class default
                        call errStop( "dotProd_rScalar3D_MR > rhs must be Scalar (try vec%scl)!" )
                    !
                end select
                !
            else
                call errStop( "dotProd_rScalar3D_MR > Incompatible rhs" )
            endif
            !
        enddo
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
        integer :: i
        !
        do i = 1, self%grid%getNGrids()
            !
            if( self%store_state .EQ. compound ) then
                !
                self%sub_scalar(i)%v = self%sub_scalar(i)%v / cvalue
                !
            elseif( self%store_state .EQ. singleton ) then
                !
                self%sub_scalar(i)%s_v = self%sub_scalar(i)%s_v / cvalue
                !
            else
                call errStop( "divByValue_rScalar3D_MR > Unknown store_state!" )
            endif
            !
        enddo
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
        integer :: i
        !
        call self%switchStoreState( rhs%store_state )
        !
        do i = 1, self%grid%getNGrids()
            !
            if( self%isCompatible( rhs ) ) then
                !
                select type( rhs )
                    !
                    class is( rScalar3D_MR_t )
                        !
                        if( rhs%store_state .EQ. compound ) then
                            !
                            self%sub_scalar(i)%v = self%sub_scalar(i)%v / rhs%sub_scalar(i)%v
                            !
                        elseif( rhs%store_state .EQ. singleton ) then
                            !
                            self%sub_scalar(i)%s_v = self%sub_scalar(i)%s_v / rhs%sub_scalar(i)%s_v
                            !
                        else
                            call errStop( "divByField_rScalar3D_MR > Unknown rhs store_state!" )
                        endif
                        !
                    class default
                        call errStop( "divByField_rScalar3D_MR > rhs must be Scalar (try vec%scl)!" )
                    !
                end select
                !
            else
                call errStop( "divByField_rScalar3D_MR > Incompatible inputs." )
            endif
            !
        enddo
        !
    end subroutine divByField_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine toNode_rScalar3D_MR( self, node_scalar, interior_only )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        class( Scalar_t ), intent( inout ) :: node_scalar
        logical, intent( in ), optional :: interior_only
        !
        call errStop( "toNode_rScalar3D_MR > Under implementation" )
        ! !
        ! type( rScalar3D_MR_t ) :: node_out_temp
        ! integer :: v_xend, v_yend, v_zend
        ! logical :: is_interior_only
        ! integer :: i
        ! !
        ! if( .NOT. self%is_allocated ) then
             ! call errStop( "toNode_rScalar3D_SG > self not allocated." )
        ! endif
        ! !
        ! call self%switchStoreState( compound )
        ! !
        ! node_out_temp = rScalar3D_MR_t( self%grid, NODE )
        ! !
        ! do i = 1, self%grid%getNGrids()
            ! !
            ! is_interior_only = .FALSE.
            ! !
            ! if( present( interior_only ) ) is_interior_only = interior_only
            ! !
            ! if( is_interior_only ) then
                ! call self%sub_scalar(i)%setAllBoundary( C_ZERO )
            ! endif
            ! !
            ! select case( self%sub_scalar(i)%grid_type )
                ! !
                ! case( CELL )
                    ! !
                    ! v_xend = size( self%sub_scalar(i)%v, 1 )
                    ! v_yend = size( self%sub_scalar(i)%v, 2 )
                    ! v_zend = size( self%sub_scalar(i)%v, 3 )
                    ! !
                    ! !> Interior
                    ! node_out_temp%sub_scalar(i)%v( 2:v_xend-1, 2:v_yend-1, 2:v_zend-1 ) = &
                    ! self%sub_scalar(i)%v( 1:v_xend-1, 1:v_yend-1, 1:v_zend-1 ) + &
                    ! self%sub_scalar(i)%v( 2:v_xend  , 1:v_yend-1, 1:v_zend-1 ) + &
                    ! self%sub_scalar(i)%v( 1:v_xend-1, 2:v_yend  , 1:v_zend-1 ) + &
                    ! self%sub_scalar(i)%v( 1:v_xend-1, 1:v_yend-1, 2:v_zend   ) + &
                    ! self%sub_scalar(i)%v( 2:v_xend  , 2:v_yend  , 1:v_zend-1 ) + &
                    ! self%sub_scalar(i)%v( 2:v_xend  , 1:v_yend-1, 2:v_zend   ) + &
                    ! self%sub_scalar(i)%v( 1:v_xend-1, 2:v_yend  , 2:v_zend   ) + &
                    ! self%sub_scalar(i)%v( 2:v_xend  , 2:v_yend  , 2:v_zend   )
                    ! !
                ! case default
                    ! call errStop( "toNode_rScalar3D_SG: undefined self%grid_type" )
            ! end select
            ! !
            ! !node_out = node_out_temp
            ! !
        ! enddo
        ! !
    end subroutine toNode_rScalar3D_MR
    !
    !> No function briefing
    !
    function getV_rScalar3D_MR( self ) result( v )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( in ) :: self
        !
        complex( kind=prec ), allocatable, dimension(:,:,:) :: v
        !
        call errStop( "getV_rScalar3D_MR not implemented!" )
        !
    end function getV_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setV_rScalar3D_MR( self, v )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), dimension(:,:,:), intent( in ) :: v
        !
        call errStop( "setV_rScalar3D_MR not implemented!" )
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
        complex( kind=prec ), allocatable, dimension(:) :: s_v
        !
        call errStop( "getSV_rScalar3D_MR not implemented!" )
        ! !
        ! if( .NOT. self%is_allocated ) then
            ! call errStop( "getSV_rScalar3D_MR > self not allocated." )
        ! endif
        ! !
        ! if( .NOT. allocated( self%s_v ) ) then
            ! call errStop( "getSV_rScalar3D_MR > self%s_v not allocated." )
        ! else
            ! !
            ! s_v = cmplx( self%s_v, 0.0, kind=prec )
            ! !
        ! endif
        ! !
    end function getSV_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine deallOtherState_rScalar3D_MR( self )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        !
        call errStop( "deallOtherState_rScalar3D_MR not implemented!" )
        ! !
        ! if( ( .NOT. self%is_allocated ) ) then
            ! call errStop( "deallOtherState_rScalar3D_MR > Self not allocated." )
        ! endif
        ! !
        ! if( self%store_state .EQ. compound ) then
            ! !
            ! if( allocated( self%s_v ) ) deallocate( self%s_v )
            ! !
        ! elseif( self%store_state .EQ. singleton ) then
            ! !
            ! if( allocated( self%v ) ) deallocate( self%v )
            ! !
        ! else
            ! call errStop( "deallOtherState_rScalar3D_MR > Unknown store_state!" )
        ! endif
        ! !
    end subroutine deallOtherState_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setSV_rScalar3D_MR( self, s_v )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), dimension(:), intent( in ) :: s_v
        !
        call errStop( "setSV_rScalar3D_MR not implemented!" )
        ! !
        ! if( .NOT. self%is_allocated ) then
            ! call errStop( "setSV_rScalar3D_MR > self not allocated." )
        ! endif
        ! !
        ! call self%switchStoreState( singleton )
        ! !
        ! if( allocated( self%v ) ) deallocate( self%v )
        ! !
        ! self%s_v = s_v
        ! !
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
        v_full = self%getFullArray()
        !
        array = v_full( self%indActive() )
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
        complex( kind=prec ), allocatable, dimension(:) :: vFull
        !
        allocate( vFull( self%lengthFull() ) )
        !
        vFull( self%indActive() ) = array
        !
        call self%setFullArray( vFull )
        !
        deallocate( vFull )
        !
    end subroutine setArray_rScalar3D_MR
    !
    !> No subroutine briefing
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
                self%NdV = rhs%NdV
                self%Nxyz = rhs%Nxyz
                !
                if( allocated( rhs%sub_scalar ) ) then
                    !
                    if( allocated( self%sub_scalar ) ) deallocate( self%sub_scalar )
                    allocate( self%sub_scalar( size( rhs%sub_scalar ) ) )
                    !
                else
                    call errStop( "copyFrom_rScalar3D_MR > rhs%sub_scalar not allocated" )
                endif
                !
                do i = 1, size( self%sub_scalar )
                    self%sub_scalar(i) = rhs%sub_scalar(i)
                enddo
                !
                self%is_allocated = .TRUE.
                !
            class default
                call errStop( "copyFrom_rScalar3D_MR > Unclassified rhs" )
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
        call errStop( "read_rScalar3D_MR not implemented!" )
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
        call errStop( "write_rScalar3D_MR not implemented!" )
        !
    end subroutine write_rScalar3D_MR
    !
    !> No subroutine briefing
    !
    subroutine print_rScalar3D_MR( self, io_unit, title, append )
        implicit none
        !
        class( rScalar3D_MR_t ), intent( in ) :: self
        integer, intent( in ), optional :: io_unit
        character(*), intent( in ), optional :: title
        logical, intent( in ), optional :: append
        !
        type( rScalar3D_MR_t ) :: copy
        integer :: i, ix, iy, iz, funit
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "print_rScalar3D_MR > self not allocated." )
        endif
        !
        copy = self
        !
        call copy%switchStoreState( compound )
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
        write(funit,*) "rScalar3D_MR field"
        do i = 1, copy%grid%getNGrids()
            !
            write( funit, * ) copy%sub_scalar(i)%nx, copy%sub_scalar(i)%ny, copy%sub_scalar(i)%nz
            !
            write(funit,*) "sub_scalar ", i
            do ix = 1, copy%sub_scalar(i)%nx
                do iy = 1, copy%sub_scalar(i)%ny
                    do iz = 1, copy%sub_scalar(i)%nz
                        if( copy%sub_scalar(i)%v( ix, iy, iz ) /= 0 ) then
                            write(funit,*) ix,iy,iz, ":[", copy%sub_scalar(i)%v( ix, iy, iz ), "]"
                        endif
                    enddo
                enddo
            enddo
            !
        enddo
        !
    end subroutine print_rScalar3D_MR
    !
end module rScalar3D_MR
