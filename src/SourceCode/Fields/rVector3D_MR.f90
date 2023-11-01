!
!> This module specializes the abstract Vector3D_real class for
!> real vector fields on a multi-resolution staggered grid.
!
module rVector3D_MR
    !
    use MatUtils
    use rVector3D_SG
    use cScalar3D_MR
    !
    type, extends( Vector_t ) :: rVector3D_MR_t
        !
        type( rVector3D_SG_t ), allocatable, dimension(:) :: sub_vector
        !
        contains
            !
            final :: rVector3D_MR_dtor
            !
            !> MR Routines
            procedure, public :: initializeSub => initializeSub_rVector3D_MR
            !
            procedure, public :: setIndexArrays => setIndexArrays_rVector3D_MR
            !
            procedure, public :: getArray => getArray_rVector3D_MR
            procedure, public :: setArray => setArray_rVector3D_MR
            !
            procedure, public :: getFullArray => getFullArray_rVector3D_MR
            procedure, public :: setFullArray => setFullArray_rVector3D_MR
            !
            procedure, public :: lengthFull => lengthFull_rVector3D_MR
            procedure, public :: findFull => findFull_rVector3D_MR
            !
            procedure, public :: MRtoSG => MRtoSG_rVector3D_MR
            procedure, public :: SGtoMRE0 => SGtoMRE0_rVector3D_MR
            !
            procedure, public :: fromSG => fromSG_rVector3D_MR
            !
            !> Boundary operations
            procedure, public :: setAllBoundary => setAllBoundary_rVector3D_MR
            procedure, public :: setOneBoundary => setOneBoundary_rVector3D_MR
            !
            !> Dimensioning operations
            procedure, public :: length => length_rVector3D_MR
            !
            !> Arithmetic/algebraic unary operations
            procedure, public :: zeros => zeros_rVector3D_MR
            !
            procedure, public :: sumEdge => sumEdge_rVector3D_MR
            procedure, public :: sumEdgeVTI => sumEdgeVTI_rVector3D_MR
            !
            procedure, public :: sumCell => sumCell_rVector3D_MR
            procedure, public :: sumCellVTI => sumCellVTI_rVector3D_MR
            !
            procedure, public :: conjugate => conjugate_rVector3D_MR
            !
            !> Arithmetic/algebraic binary operations
            procedure, public :: add => add_rVector3D_MR
            !
            procedure, public :: linComb => linComb_rVector3D_MR
            !
            procedure, public :: subValue => subValue_rVector3D_MR
            procedure, public :: subField => subField_rVector3D_MR
            !
            procedure, public :: multByReal => multByReal_rVector3D_MR
            procedure, public :: multByComplex => multByComplex_rVector3D_MR
            procedure, public :: multByField => multByField_rVector3D_MR
            !
            procedure, public :: diagMult => diagMult_rVector3D_MR
            !
            procedure, public :: multAdd => multAdd_rVector3D_MR
            !
            procedure, public :: dotProd => dotProd_rVector3D_MR
            !
            procedure, public :: divByValue => divByValue_rVector3D_MR
            procedure, public :: divByField => divByField_rVector3D_MR
            !
            procedure, public :: interpFunc => interpFunc_rVector3D_MR
            !
            !> Miscellaneous
            procedure, public :: getAxis => getAxis_rVector3D_MR
            !
            procedure, public :: getReal => getReal_rVector3D_MR
            !
            procedure, public :: deallOtherState => deallOtherState_rVector3D_MR
            !
            procedure, public :: copyFrom => copyFrom_rVector3D_MR
            !
            !> I/O operations
            procedure, public :: print => print_rVector3D_MR
            procedure, public :: read => read_rVector3D_MR
            procedure, public :: write => write_rVector3D_MR
            !
    end type rVector3D_MR_t
    !
    interface rVector3D_MR_t
        module procedure rVector3D_MR_ctor
    end interface rVector3D_MR_t
    !
    public :: setInactiveEdge_rVector3D_MR
    public :: addCellFromAdjacentGrid_rVector3D_MR
    !
contains
    !
    !> No function briefing
    !
    function rVector3D_MR_ctor( grid, grid_type ) result(  self )
        implicit none
        !
        class( Grid_t ), target, intent( in ) :: grid
        character( len=4 ), intent( in ) :: grid_type
        !
        type( rVector3D_MR_t ) :: self
        !
        call self%baseInit
        !
        self%grid => grid
        self%grid_type = grid_type
        !
        self%nx = self%grid%nx
        self%ny = self%grid%ny
        self%nz = self%grid%nz
        !
        call self%initializeSub
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "rVector3D_MR_ctor > Unable to allocate vector." )
        endif
        !
    end function rVector3D_MR_ctor
    !
    !> No subroutine briefing
    !
    subroutine initializeSub_rVector3D_MR( self ) 
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        !
        integer :: i, alloc_stat
        !
        select type( grid => self%grid )
            !
            class is( Grid3D_MR_t )
                !
                self%is_allocated = .TRUE.
                allocate( self%sub_vector( grid%n_grids ), stat = alloc_stat )
                self%is_allocated = self%is_allocated .AND. (  alloc_stat .EQ. 0 )
                !
                do i = 1, grid%n_grids
                    !
                    self%sub_vector(i) = rVector3D_SG_t( grid%sub_grid(i), self%grid_type )
                    !
                    !write( *, * ) "rSubVector", i, "-nx=", self%sub_vector(i)%nx, ", ny=", self%sub_vector(i)%ny, "nz=", self%sub_vector(i)%nz
                    !
                enddo
                !
                !write( *, * ) "rMainVector-nx=", self%nx, ", ny=", self%ny, "nz=", self%nz, self%nx*self%ny*self%nz
                !
            class default
                call errStop( "initializeSub_rVector3D_MR > Unclassified grid" )
            !
        end select
        !
    end subroutine initializeSub_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setIndexArrays_rVector3D_MR( self, n_full, ind_boundary, ind_interior, ind_active, xy_in ) 
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
        integer, intent( inout ) :: n_full
        integer, dimension(:), allocatable, intent( out ) :: ind_boundary, ind_interior
        integer, dimension(:), allocatable, intent( out ), optional :: ind_active
        logical, intent( in ), optional :: xy_in
        !
        type( rVector3D_MR_t ) :: temp_vector
        logical :: xy, int_only
        integer :: i, k
        integer :: n_active, n_interior, n_boundaries
        real( kind=prec ), dimension(:), allocatable :: v_1, v_2
        !
        if( .NOT. present( xy_in ) ) then
            xy = .FALSE.
        else
            xy = xy_in
        endif
        !
        temp_vector = self
        !
        select type( grid => temp_vector%grid )
            !
            class is( Grid3D_MR_t )
                !
                ! Loop over sub-grids, setting boundary edges to one,
                ! interior to  zero
                do k = 1, grid%n_grids
                    call temp_vector%sub_vector(k)%setAllBoundary( cmplx( 1._prec, 0.0, kind=prec ) )
                enddo
                !
                ! Loop over interfaces: set redundant interface edges to 2
                select case( temp_vector%grid_type )
                    !
                    case( EDGE )
                        int_only = .TRUE.
                    case( FACE )
                        int_only = .FALSE.
                    case default
                        !
                        call errStop( "setIndexArrays_rVector3D_MR > Invalid grid type option!" )
                    !
                end select
                !
                do k = 2, grid%n_grids
                    !
                    if( grid%coarseness(k - 1, 1) < grid%coarseness(k, 1) ) then
                        ! upper grid is finer: grid k-1 interface nodes are
                        ! not active; also reset interior part of interface
                        ! edges to 0
                        if( xy ) then
                            call temp_vector%sub_vector(k-1)%setOneBoundary( "z2_x", cmplx( -1.0_prec, 0.0, kind=prec ) )
                            call temp_vector%sub_vector(k-1)%setOneBoundary( "z2_y", cmplx( -10.0_prec, 0.0, kind=prec ) )
                        else
                            call temp_vector%sub_vector(k-1)%setOneBoundary( "z2", cmplx( -1.0_prec, 0.0, kind=prec ) )
                        endif
                        !
                        call temp_vector%sub_vector(k)%setOneBoundary( "z1", cmplx( 0._prec, 0.0, kind=prec ), int_only )
                    else
                        if( xy ) then
                            call temp_vector%sub_vector(k)%setOneBoundary( "z1_x", cmplx( -1.0_prec, 0.0, kind=prec ) )
                            call temp_vector%sub_vector(k)%setOneBoundary( "z1_y", cmplx( -10.0_prec, 0.0, kind=prec ) )
                        else
                            call temp_vector%sub_vector(k)%setOneBoundary( "z1", cmplx( -1.0_prec, 0.0, kind=prec ) )
                        endif
                        !
                        call temp_vector%sub_vector(k-1)%setOneBoundary( "z2", cmplx( 0._prec, 0.0, kind=prec ), int_only )
                        !
                    endif
                    !
                enddo
                !
            class default
                call errStop( "setIndexArrays_rVector3D_MR > Unclassified grid" )
            !
        end select
        !
        ! Set active, interior, and boundary edges. ***
        !
        v_1 = temp_vector%getFullArray()
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
            if(v_1(k) >= 0) then
                i = i + 1
                ind_active(i) = k
            endif
        enddo
        !
        n_interior = 0
        do k = 1, n_full
            if(v_1(k) == 0) then
                n_interior = n_interior + 1
            endif
        enddo
        !
        allocate( v_2(n_active))
        v_2 = v_1(ind_active)
        !
        if(allocated( ind_interior)) then
            deallocate( ind_interior)
        endif
        allocate( ind_interior(n_interior))
        !
        i = 0
        do k = 1, n_active
            if(v_2(k) == 0) then
                i = i + 1
                ind_interior(i) = k
            endif
        enddo
        !!
        n_boundaries = 0
        do k = 1, n_active
            if(v_2(k) == 1) then
                n_boundaries = n_boundaries + 1
            endif
        enddo
        !
        if(allocated( ind_boundary)) then
            deallocate( ind_boundary)
        endif
        allocate( ind_boundary(n_boundaries)) 
        !
        i = 0
        do k = 1, n_active
            if(v_2(k) == 1) then
                i = i + 1
                ind_boundary(i) = k
            endif
        enddo
        !
    end subroutine setIndexArrays_rVector3D_MR
    !
    !> Creates standard( 1-D array) for all sub-scalars,
    !> INCLUDING redundant interface nodes.
    !
    function getFullArray_rVector3D_MR( self ) result( array )
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
        !
        complex( kind=prec ), allocatable, dimension(:) :: array
        !
        integer :: n, i, i1, i2
        !
        n = self%lengthFull()
        allocate( array(n) )
        !
        array = C_ZERO
        !
        i1 = 1
        i2 = 0
        !
        do i = 1, self%grid%getNGrids()
            !
            n = self%sub_vector(i)%length()
            i2 = i2 + n
            array(i1:i2) = self%sub_vector(i)%getArray()
            i1 = i1 + n
            !
        enddo
        !
    end function getFullArray_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setFullArray_rVector3D_MR( self, array )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), dimension(:), intent( in ) :: array
        !
        integer :: i1, i2, k, n
        !
        select type( grid => self%grid )
            !
            class is( Grid3D_MR_t )
                !
                i1 = 1
                i2 = 0
                !
                do k = 1, grid%n_grids
                    !
                    n = self%sub_vector(k)%length()
                    i2 = i2 + n
                    call self%sub_vector(k)%setArray( array(i1:i2) )
                    i1 = i1 + n
                    !
                enddo
                !
            class default
                call errStop( "setFullArray_rVector3D_MR > Unclassified grid" )
            !
        end select
        !
    end subroutine setFullArray_rVector3D_MR
    !
    !> No function briefing
    !
    function lengthFull_rVector3D_MR( self ) result( n )
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
        !
        integer :: i, n
        !
        n = 0
        !
        do i = 1, self%grid%getNGrids()
            !
            n = n + self%sub_vector(i)%length()
            !
        enddo
        !
    end function lengthFull_rVector3D_MR
    !
    !> No function briefing
    !
    function findFull_rVector3D_MR( self, c ) result( I )
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
        real( kind=prec ), intent( in) :: c
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
            if( v(k) == c ) n_I = n_I + 1
        enddo
        !
        allocate( I(n_I))
        !
        n_I = 0
        do k = 1, n
            if( v(k) == c ) then
                n_I = n_I + 1
                I(n_I) = k
            endif
        enddo
        !
    end function findFull_rVector3D_MR
    !
    !> Convert an MR TVector to a full SG, filling in the full fine grid.
    !> Converts MR Vector object to SG
    !> copying from variable resolution sub-grids to completely fill in the
    !> underlying fine grid.
    !
    subroutine MRtoSG_rVector3D_MR( self, sg_v )
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
        type( rVector3D_SG_t ), intent( inout ) :: sg_v
        !
        type( rVector3D_SG_t ) :: temp
        integer :: x_nx, x_ny, x_nz
        integer :: y_nx, y_ny, y_nz
        integer :: z_nx, z_ny, z_nz
        integer :: last, Cs, i1, i2, i, k
        real( kind=prec ) :: w1, w2
        !
        temp = rVector3D_SG_t( self%grid, self%grid_type )
        !
        temp%x = 0; temp%y = 0; temp%z = 0
        !
        x_nx = size(temp%x, 1)
        x_ny = size(temp%x, 2)
        x_nz = size(temp%x, 3)
        !
        y_nx = size(temp%y, 1)
        y_ny = size(temp%y, 2)
        y_nz = size(temp%y, 3)
        !
        z_nx = size(temp%z, 1)
        z_ny = size(temp%z, 2)
        z_nz = size(temp%z, 3)
        !
        sg_v%x = 0; sg_v%y = 0; sg_v%z = 0
        !
        select case( self%grid_type )
            !
            case( EDGE )
                !
                select type( grid => self%grid )
                    !
                    class is( Grid3D_MR_t )
                        !
                        do k = 1, grid%n_grids
                            !
                            Cs = 2**grid%coarseness(k, 1)
                            i1 = grid%coarseness(k, 3)
                            i2 = grid%coarseness(k, 4)
                            ! Copy  x and y components in x and y directions
                            ! edges that aligned with sub-grid edge.
                            do i = 1, Cs
                                !
                                temp%x(i:x_nx:Cs, 1:x_ny:Cs, i1:i2+1) = self%sub_vector(k)%x
                                temp%y(1:y_nx:Cs, i:y_ny:Cs, i1:i2+1) = self%sub_vector(k)%y
                                !
                                w1 = 1. -( i - 1.)/Cs
                                w2 = 1. - w1
                                !
                                if(i == 1) then
                                    temp%z(1:z_nx:Cs, 1:z_ny:Cs, i1:i2) = self%sub_vector(k)%z
                                else
                                    last = size(self%sub_vector(k)%z(:, 1, 1))
                                    temp%z(i:z_nx:Cs, 1:z_ny:Cs, i1:i2) = &
                                    self%sub_vector(k)%z(1:last-1, :, :) * &
                                    w1 + self%sub_vector(k)%z(2:last, :, :) * w2
                                endif
                                !
                            enddo
                            ! edges that subdivide the sub-grid
                            ! interpolate  in y and x directions
                            ! copy/interpolate x in y direction
                            ! copy x and y along x and y directions
                            ! respectively
                            do i = 2, Cs
                                w1 = 1. -( i - 1.)/Cs
                                w2 = 1. - w1

                                temp%x(:, i:x_ny:Cs, i1:i2+1) = temp%x(:, 1:x_ny-Cs:Cs, i1:i2+1)*w1 + &
                                temp%x(:, Cs+1:x_ny:Cs, i1:i2+1)*w2

                                temp%y(i:y_nx:Cs, :, i1:i2+1) = temp%y(1:y_nx-Cs:Cs, :, i1:i2+1)*w1 + &
                                temp%y(Cs+1:y_nx:Cs, :, i1:i2+1)*w2

                                ! added by zhhq, 2017
                                temp%z(:, i:z_ny:Cs, i1:i2) = temp%z(:, 1:z_ny-Cs:Cs, i1:i2)*w1 + &
                                temp%z(:, Cs+1:z_ny:Cs, i1:i2) * w2
                                ! temp.z(i:Cs:end,i:Cs:end,i1:i2) = temp.z(:,1:Cs:end-Cs,i1:i2)*w1+ ...
                                ! temp.z(:,Cs+1:Cs:end,i1:i2)*w2;
                                ! added by zhhq, 2017
                            enddo
                            !
                            sg_v%x(:, :, i1:i2+1) = sg_v%x(:, :, i1:i2+1) + temp%x(:, :, i1:i2+1)
                            sg_v%y(:, :, i1:i2+1) = sg_v%y(:, :, i1:i2+1) + temp%y(:, :, i1:i2+1)
                            sg_v%z(:, :, i1:i2)   = temp%z(:, :, i1:i2)
                            !
                        enddo
                        !
                    class default
                        call errStop( "MRtoSG_rVector3D_MR > Unclassified grid" )
                    !
                end select
                !
            case default
                call errStop( "MRtoSG_rVector3D_MR > Unrecognized grid type." )
            !
        end select
        !
    end subroutine MRtoSG_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine fromSG_rVector3D_MR( self, vector_sg )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        type( rVector3D_SG_t ), intent( in ) :: vector_sg
        !
        type( rVector3D_SG_t ) :: temp_edge_length_sg, sg_vector_l, sg_vector_el
        real( kind=prec ), allocatable, dimension(:) :: array_l, array_e, array_le
        real( kind=prec ), allocatable, dimension(:,:,:) :: x_length, y_length
        integer :: sx1, sx2, sx3, sy1, sy2, sy3, s1, s2
        integer :: i_grid, Cs, i, i1, i2
        !
        if( .NOT. vector_sg%is_allocated ) then
            call errStop( "fromSG_rVector3D_MR > vector_sg not allocated" )
        endif
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "fromSG_rVector3D_MR > self not allocated" )
        endif
        !
        !>
        write( *, * ) "fromSG_rVector3D_MR SG: ", vector_sg%grid%nx, vector_sg%grid%ny, vector_sg%grid%nz, vector_sg%grid%n_grids
        write( *, * ) "                    MR: ", self%grid%nx, self%grid%ny, self%grid%nz, self%grid%n_grids
        !
        call vector_sg%edgeLength( temp_edge_length_sg )
        !
        array_l = temp_edge_length_sg%getArray()
        !
        array_e = vector_sg%getArray()
        !
        allocate( array_le( size( array_l ) ) )
        !
        array_le = array_l * array_e
        !
        sg_vector_l = rVector3D_SG_t( vector_sg%grid, EDGE )
        !
        call sg_vector_l%setArray( cmplx( array_l, 0.0, kind=prec ) )
        !
        sg_vector_el = rVector3D_SG_t( vector_sg%grid, EDGE )
        !
        call sg_vector_el%setArray( cmplx( array_le, 0.0, kind=prec ) )
        !
        select type( grid => vector_sg%grid )
            !
            class is( Grid3D_MR_t )
                !
                do i_grid = 1, grid%n_grids
                    !
                    self%sub_vector( i_grid )%grid => grid%sub_grid( i_grid )
                    !
                    Cs = 2 ** grid%coarseness( i_grid, 1 )
                    i1 = grid%coarseness( i_grid, 3 )
                    i2 = grid%coarseness( i_grid, 4 )
                    !
                    sx1 = size( self%sub_vector(i_grid)%x, 1 )
                    sx2 = size( self%sub_vector(i_grid)%x, 2 )
                    sx3 = size( self%sub_vector(i_grid)%x, 3 )
                    allocate( x_length( sx1, sx2, sx3 ) )
                    x_length = 0.0
                    !
                    sy1 = size( self%sub_vector(i_grid)%y, 1 )
                    sy2 = size( self%sub_vector(i_grid)%y, 2 )
                    sy3 = size( self%sub_vector(i_grid)%y, 3 )
                    allocate( y_length( sy1, sy2, sy3 ) )
                    y_length = 0.0
                    !
                    do i = 1, Cs
                        !
                        s1 = size( sg_vector_el%x, 1 )
                        s2 = size( sg_vector_el%x, 2 )
                        self%sub_vector(i_grid)%x = self%sub_vector(i_grid)%x + &
                        sg_vector_el%x( i:s1:Cs, 1:s2:Cs, i1:i2+1 )
                        !
                        s1 = size( sg_vector_el%y, 1 )
                        s2 = size( sg_vector_el%y, 2 )
                        self%sub_vector(i_grid)%y = self%sub_vector(i_grid)%y + &
                        sg_vector_el%y( 1:s1:Cs, i:s2:Cs, i1:i2+1 )
                        !
                        s1 = size( sg_vector_l%x, 1 )
                        s2 = size( sg_vector_l%x, 2 )
                        x_length = x_length + sg_vector_l%x( i:s1:Cs, 1:s2:Cs, i1:i2+1 )
                        !
                        s1 = size( sg_vector_l%y, 1 )
                        s2 = size( sg_vector_l%y, 2 )
                        y_length = y_length + sg_vector_l%y( 1:s1:Cs, i:s2:Cs, i1:i2+1 )
                        !
                    enddo
                    !
                    self%sub_vector(i_grid)%x = self%sub_vector(i_grid)%x / x_length
                    self%sub_vector(i_grid)%y = self%sub_vector(i_grid)%y / y_length
                    !
                    s1 = size( vector_sg%z, 1 )
                    s2 = size( vector_sg%z, 2 )
                    self%sub_vector(i_grid)%z = vector_sg%z( 1:s1:Cs, 1:s2:Cs, i1:i2 )
                    !
                    deallocate( x_length, y_length )
                    !
                enddo
                !
            class default
                call errStop( "fromSG_rVector3D_MR > Unclassified grid" )
            !
        end select
        !
    end subroutine fromSG_rVector3D_MR
    !
    !> Converts SG Vector object to MR by averaging
    !> used only for e0
    !
    subroutine SGtoMRE0_rVector3D_MR(self, sg_v)
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        type( rVector3D_SG_t ), intent( in ) :: sg_v
        !
        class( Grid_t ), pointer :: grid
        !
        integer :: x_nx, x_ny, x_nz
        integer :: y_nx, y_ny, y_nz
        integer :: z_nx, z_ny, z_nz
        integer :: last, Cs, i1, i2, i, k
        !
        grid => sg_v%grid
        !
        x_nx = size(sg_v%x, 1)
        x_ny = size(sg_v%x, 2)
        x_nz = size(sg_v%x, 3)
        !
        y_nx = size(sg_v%y, 1)
        y_ny = size(sg_v%y, 2)
        y_nz = size(sg_v%y, 3)
        !
        z_nx = size(sg_v%z, 1)
        z_ny = size(sg_v%z, 2)
        z_nz = size(sg_v%z, 3)
        !
        select type( grid => self%grid )
            !
            class is( Grid3D_MR_t )
                !
                do k = 1, grid%n_grids
                    !
                    Cs = 2**grid%coarseness(k, 1)
                    i1 = grid%coarseness(k, 3)
                    i2 = grid%coarseness(k, 4)
                    !
                    do i = 1, Cs
                        !
                        last = size(grid%Dx)
                        self%sub_vector(k)%x = self%sub_vector(k)%x + &
                        sg_v%x(i:x_nx:Cs, 1:x_ny:Cs, i1:i2+1) *    &
                        repMat(grid%Dx(i:last:Cs), &
                        1, &
                        grid%sub_grid(k)%Ny + 1, &
                        grid%sub_grid(k)%Nz + 1, .false.)

                        last = size(grid%Dy)
                        self%sub_vector(k)%y = self%sub_vector(k)%y + &
                        sg_v%y(1:y_nx:Cs, i:y_ny:Cs, i1:i2+1) *  & 
                        repMat(grid%Dy(i:last:Cs), &
                        grid%sub_grid(k)%Nx + 1, &
                        1, &
                        grid%sub_grid(k)%Nz + 1, .TRUE.)
                        !
                    enddo
                    !
                    self%sub_vector(k)%x = self%sub_vector(k)%x / &
                    repMat(grid%sub_grid(k)%Dx, &
                    1, &
                    grid%sub_grid(k)%Ny + 1, &
                    grid%sub_grid(k)%Nz + 1, .false.)
                    !
                    self%sub_vector(k)%y = self%sub_vector(k)%y / &
                    repMat(grid%sub_grid(k)%Dy, &
                    grid%sub_grid(k)%Nx + 1, &
                    1, &
                    grid%sub_grid(k)%Nz + 1, .TRUE.)
                    !
                    self%sub_vector(k)%z = sg_v%z(1:z_nx:Cs, 1:z_ny:Cs, i1:i2)
                    !
                enddo
                !
            class default
                call errStop( "SGtoMRE0_rVector3D_MR > Unclassified grid" )
            !
        end select
        !
    end subroutine SGtoMRE0_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine rVector3D_MR_dtor( self )
        implicit none
        !
        type( rVector3D_MR_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor rVector3D_MR"
        !
        call self%baseDealloc
        !
        if( allocated( self%sub_vector ) ) deallocate( self%sub_vector )
        !
        self%nx = 0
        self%ny = 0
        self%nz = 0
        !
        self%grid_type = ""
        self%is_allocated = .FALSE.
        !
    end subroutine rVector3D_MR_dtor
    !
    !> No subroutine briefing
    !
    subroutine setAllBoundary_rVector3D_MR( self, cvalue )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        call errStop( "setAllBoundary_rVector3D_MR just implemented for SG!" )
        !
    end subroutine setAllBoundary_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setOneBoundary_rVector3D_MR( self, bdry, cvalue, int_only )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        character(*), intent( in ) :: bdry
        complex( kind=prec ), intent( in ) :: cvalue
        logical, intent( in ), optional :: int_only
        !
        logical :: int_only_p
        integer :: i
        !
        if( .NOT. present(int_only)) then
            int_only_p = .FALSE.
        else
            int_only_p = int_only
        endif
        !
        call self%switchStoreState( compound )
        !
        do i = 1, self%grid%getNGrids()
            !
            select case( self%sub_vector(i)%grid_type )
                !
                case( EDGE )
                    !
                    if( int_only_p ) then
                        !
                        select case( bdry )
                            !
                            case("x1")
                                self%sub_vector(i)%z(1, 2:self%sub_vector(i)%NdZ(2)-1, :) = real( cvalue, kind=prec )
                                self%sub_vector(i)%y(1, :, 2:self%sub_vector(i)%NdY(3)-1) = real( cvalue, kind=prec )
                            case("x2")
                                self%sub_vector(i)%z(self%sub_vector(i)%NdZ(1), 2:self%sub_vector(i)%NdZ(2)-1, :) = real( cvalue, kind=prec )
                                self%sub_vector(i)%y(self%sub_vector(i)%NdY(1), :, 2:self%sub_vector(i)%NdY(3)-1) = real( cvalue, kind=prec )
                            case("y1")
                                self%sub_vector(i)%z(2:self%sub_vector(i)%NdZ(1)-1, 1, :) = real( cvalue, kind=prec )
                                self%sub_vector(i)%x(:, 1, 2:self%sub_vector(i)%NdX(3)-1) = real( cvalue, kind=prec )
                            case("y2")
                                self%sub_vector(i)%z(2:self%sub_vector(i)%NdZ(1)-1, self%sub_vector(i)%NdZ(2), :) = real( cvalue, kind=prec )
                                self%sub_vector(i)%x(:, self%sub_vector(i)%NdX(2), 2:self%sub_vector(i)%NdX(3)-1) = real( cvalue, kind=prec )
                            case("z1")
                                self%sub_vector(i)%x(:, 2:self%sub_vector(i)%NdX(2)-1, 1) = real( cvalue, kind=prec )
                                self%sub_vector(i)%y(2:self%sub_vector(i)%NdY(1)-1, :, 1) = real( cvalue, kind=prec )
                            case("z2")
                                self%sub_vector(i)%x(:, 2:self%sub_vector(i)%NdX(2)-1, self%sub_vector(i)%NdX(3)) = real( cvalue, kind=prec )
                                self%sub_vector(i)%y(2:self%sub_vector(i)%NdY(1)-1, :, self%sub_vector(i)%NdY(3)) = real( cvalue, kind=prec )
                            case("z1_x")
                                self%sub_vector(i)%x(:, 2:self%sub_vector(i)%NdX(2)-1, 1) = real( cvalue, kind=prec )
                            case("z2_x")
                                self%sub_vector(i)%x(:, 2:self%sub_vector(i)%NdX(2)-1, self%sub_vector(i)%NdX(3)) = real( cvalue, kind=prec )
                            case("z1_y")
                                self%sub_vector(i)%y(2:self%sub_vector(i)%NdY(1)-1, :, 1) = real( cvalue, kind=prec )
                            case("z2_y")
                                self%sub_vector(i)%y(2:self%sub_vector(i)%NdY(1)-1, :, self%sub_vector(i)%NdY(3)) = real( cvalue, kind=prec )
                            !
                        end select
                        !
                    else
                        !
                        select case( bdry )
                            !
                            case("x1")
                                self%sub_vector(i)%z(1, :, :) = real( cvalue, kind=prec )
                                self%sub_vector(i)%y(1, :, :) = real( cvalue, kind=prec )
                            case("x2")
                                self%sub_vector(i)%z(self%sub_vector(i)%NdZ(1), :, :) = real( cvalue, kind=prec )
                                self%sub_vector(i)%y(self%sub_vector(i)%NdY(1), :, :) = real( cvalue, kind=prec )
                            case("y1")
                                self%sub_vector(i)%z(:, 1, :) = real( cvalue, kind=prec )
                                self%sub_vector(i)%x(:, 1, :) = real( cvalue, kind=prec )
                            case("y2")
                                self%sub_vector(i)%z(:, self%sub_vector(i)%NdZ(2), :) = real( cvalue, kind=prec )
                                self%sub_vector(i)%x(:, self%sub_vector(i)%NdX(2), :) = real( cvalue, kind=prec )
                            case("z1")
                                self%sub_vector(i)%x(:, :, 1) = real( cvalue, kind=prec )
                                self%sub_vector(i)%y(:, :, 1) = real( cvalue, kind=prec )
                            case("z2")
                                self%sub_vector(i)%x(:, :, self%sub_vector(i)%NdX(3)) = real( cvalue, kind=prec )
                                self%sub_vector(i)%y(:, :, self%sub_vector(i)%NdY(3)) = real( cvalue, kind=prec )
                            case("z1_x")
                                self%sub_vector(i)%x(:, :, 1) = real( cvalue, kind=prec )
                            case("z2_x")
                                self%sub_vector(i)%x(:, :, self%sub_vector(i)%NdX(3)) = real( cvalue, kind=prec )
                            case("z1_y")
                                self%sub_vector(i)%y(:, :, 1) = real( cvalue, kind=prec )
                            case("z2_y")
                                self%sub_vector(i)%y(:, :, self%sub_vector(i)%NdY(3)) = real( cvalue, kind=prec )
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
                            self%sub_vector(i)%x(1, :, :) = real( cvalue, kind=prec )
                        case("x2")
                            self%sub_vector(i)%x(self%sub_vector(i)%NdX(1), :, :) = real( cvalue, kind=prec )
                        case("y1")
                            self%sub_vector(i)%y(:, 1, :) = real( cvalue, kind=prec )
                        case("y2")
                            self%sub_vector(i)%y(:, self%sub_vector(i)%NdY(2), :) = real( cvalue, kind=prec )
                        case("z1")
                            self%sub_vector(i)%z(:, :, 1) = real( cvalue, kind=prec )
                        case("z2")
                            self%sub_vector(i)%z(:, :, self%sub_vector(i)%NdZ(3)) = real( cvalue, kind=prec )
                        !
                    end select
                    !
                case default
                    call errStop( "setOneBoundary_rVector3D_MR > Invalid grid type." )
            end select
            !
        enddo
        !
    end subroutine setOneBoundary_rVector3D_MR
    !
    !> No subroutine briefing
    !
    function length_rVector3D_MR( self ) result( field_length )
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
        !
        integer :: field_length
        !
        if( self%grid_type == EDGE ) then
            !
            field_length = size( self%grid%EDGEa )
            !
        elseif( self%grid_type == FACE ) then
            !
            field_length = size( self%grid%FACEa )
            !
        else
            call errStop( "length_rVector3D_MR > unrecognized grid type: ["//self%grid_type//"]" )
        endif
        !
    end function length_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine zeros_rVector3D_MR( self )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        !
        integer :: i
        !
        do i = 1, self%grid%getNGrids()
            call self%sub_vector(i)%zeros()
        enddo
        !
    end subroutine zeros_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine sumEdge_rVector3D_MR( self, cell_out, interior_only )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        class( Scalar_t ), allocatable, intent( out ) :: cell_out
        logical, intent( in ), optional :: interior_only
        !
        integer :: i
        logical :: is_interior_only
        class( Scalar_t ), allocatable :: aux_scalar
        !
        if( .NOT. self%is_allocated ) then
             call errStop( "sumEdge_rVector3D_MR > self not allocated." )
        endif
        !
        call self%switchStoreState( compound )
        !
        is_interior_only = .FALSE.
        !
        if( present( interior_only ) ) is_interior_only = interior_only
        !
        if( is_interior_only ) then
            call self%setAllBoundary( C_ZERO )
        endif
        !
        allocate( cell_out, source = rScalar3D_MR_t( self%grid, CELL ) )
        !
        select type( cell_out )
            !
            class is( rScalar3D_MR_t )
                !
                select case( self%grid_type )
                    !
                    case( EDGE )
                        !
                        !> loop over INTERFACES (one less than n_grids) and fill in inactive edges
                        do i = 1, self%grid%n_grids-1
                            call setInactiveEdge_rVector3D_MR( self%sub_vector(i), self%sub_vector(i+1) )
                        enddo
                        !
                        !> loop over sub-vectors and sum edges onto cells
                        do i = 1, self%grid%n_grids
                            !
                            call self%sub_vector(i)%sumEdges( aux_scalar )
                            !
                            cell_out%sub_scalar(i) = aux_scalar
                            !
                            deallocate( aux_scalar )
                            !
                        enddo
                        !
                    case( FACE )
                        !
                        call errStop( "sumEdge_rVector3D_MR: Unknown type FACE not implemented" )
                        !
                    case default
                        call errStop( "sumEdge_rVector3D_MR: Unknown type" )
                    !
                end select !type
                !
            class default
                call errStop( "sumEdge_rVector3D_MR > Unclassified cell_in" )
            !
        end select
        !
    end subroutine sumEdge_rVector3D_MR
    !
    subroutine sumEdgeVTI_rVector3D_MR( self, cell_h_out, cell_v_out, interior_only )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        class( Scalar_t ), allocatable, intent( out ) :: cell_h_out, cell_v_out
        logical, optional, intent( in ) :: interior_only
        !
        call errStop( "sumEdgeVTI_rVector3D_MR > Under implementation" )
        !
    end subroutine sumEdgeVTI_rVector3D_MR
    !
    !> need to sort out face type (not really used as far as I recall!)
    !> this should not be too hard, but I need to think about it for now
    !> could just code for ptype = EDGE
    !> NOTE: MAKE EVERYTHING sumCells -- get rid of the divide by 4 (etc.)
    !> 
    !> algorithm -- details (declarations, allocation, etc.)  need to be filled in
    !> 
    !> self is of class cScalar3D_MR (or rSscalar3D_MR); edge_obj (output)
    !> is Vector of corresponding type
    !> 
    !> This is adjoint of sumEdges -- so reverse order of operations i
    !> (e.g., for matrix multiplies (AB)' = B'A')
    !
    subroutine sumCell_rVector3D_MR( self, cell_in, ptype )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        class( Scalar_t ), intent( in ) :: cell_in
        character(*), intent( in ), optional :: ptype
        !
        integer :: i
        logical :: top_coarser
        character(4) :: grid_type
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "sumCell_rVector3D_MR > self not allocated." )
        endif
        !
        if( .NOT. cell_in%is_allocated ) then
            call errStop( "sumCell_rVector3D_MR > cell_in not allocated." )
        endif
        !
        if( index( self%grid_type, CELL ) > 0 ) then
            call errStop( "sumCell_rVector3D_MR > Only CELL type supported." )
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
            class is( rScalar3D_MR_t )
                !
                select case( grid_type )
                    !
                    case( EDGE )
                        !
                        !> loop over sub-sclars and sum cells onto edge for each
                        do i = 1, self%grid%n_grids
                            !
                            call self%sub_vector(i)%sumCells( cell_in%sub_scalar(i) )
                            !
                        enddo
                        !
                        !> then loop over INTERFACES (one less than n_grids) and fill in inactive edges
                        do i = 1, cell_in%grid%n_grids-1
                            !
                            !> 
                            top_coarser = cell_in%sub_scalar(i)%grid%nx .LT.  & 
                            cell_in%sub_scalar(i+1)%grid%nx
                            !
                            if( top_coarser ) then
                                call addCellFromAdjacentGrid_rVector3D_MR( self%sub_vector(i), &
                                cell_in%sub_scalar(i+1), top_coarser )
                            else
                                call addCellFromAdjacentGrid_rVector3D_MR( self%sub_vector(i+1), &
                                cell_in%sub_scalar(i), top_coarser )
                            endif
                            !
                        enddo
                        !
                    case( FACE )
                        !
                        call errStop( "sumCell_rVector3D_MR: Unknown type FACE not implemented" )
                        !
                    case default
                        !
                        call errStop( "sumCell_rVector3D_MR: Unknown type" )
                        !
                end select !type
                !
            class default
                call errStop( "sumCell_rVector3D_MR > Unclassified cell_in" )
            !
        end select
        !
    end subroutine sumCell_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine sumCellVTI_rVector3D_MR( self, cell_h_in, cell_v_in, ptype )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        class( Scalar_t ), intent( in ) :: cell_h_in, cell_v_in
        character(*), intent( in ), optional :: ptype
        !
        complex( kind=prec ), allocatable :: v_h(:,:,:), v_v(:,:,:)
        character(10) :: grid_type
        integer :: xend, yend, zend
        integer :: v_xend, v_yend, v_zend
        integer :: i, ix, iy, iz
        !
        call errStop( "sumCellVTI_rVector3D_MR > Under implementation" )
        !
    end subroutine sumCellVTI_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine conjugate_rVector3D_MR( self )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        !
        call errStop( "conjugate_rVector3D_MR: Do not try to conjugate a real vector!" )
        !
    end subroutine conjugate_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine add_rVector3D_MR( self, rhs )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "add_rVector3D_MR > self not allocated." )
        endif
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "add_rVector3D_MR > rhs not allocated." )
        endif
        !
        call self%switchStoreState( rhs%store_state )
        !
        if( self%isCompatible( rhs ) ) then
            !
            do i = 1, self%grid%getNGrids()
                !
                if( rhs%store_state .EQ. compound ) then
                    !
                    select type( rhs )
                        !
                        class is( rVector3D_MR_t )
                            !
                            self%sub_vector(i)%x = self%sub_vector(i)%x + rhs%sub_vector(i)%x
                            self%sub_vector(i)%y = self%sub_vector(i)%y + rhs%sub_vector(i)%y
                            self%sub_vector(i)%z = self%sub_vector(i)%z + rhs%sub_vector(i)%z
                            !
                        class is( rScalar3D_MR_t )
                            !
                            self%sub_vector(i)%x = self%sub_vector(i)%x + rhs%sub_scalar(i)%v
                            self%sub_vector(i)%y = self%sub_vector(i)%y + rhs%sub_scalar(i)%v
                            self%sub_vector(i)%z = self%sub_vector(i)%z + rhs%sub_scalar(i)%v
                            !
                        class default
                            call errStop( "add_rVector3D_MR > Undefined compound rhs" )
                            !
                    end select
                    !
                elseif( rhs%store_state .EQ. singleton ) then
                    !
                    select type( rhs )
                        !
                        class is( rVector3D_MR_t )
                            !
                            self%sub_vector(i)%s_v = self%sub_vector(i)%s_v + rhs%sub_vector(i)%s_v
                            !
                        class is( rScalar3D_MR_t )
                            !
                            self%sub_vector(i)%s_v = self%sub_vector(i)%s_v + rhs%sub_scalar(i)%s_v
                            !
                        class default
                            call errStop( "add_rVector3D_MR > Undefined singleton rhs" )
                            !
                    end select
                    !
                else
                    call errStop( "add_rVector3D_MR > Unknow store_state." )
                endif
                !
            enddo
            !
        else
            call errStop( "add_rVector3D_MR > Incompatible inputs." )
        endif
        !
    end subroutine add_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine linComb_rVector3D_MR( self, rhs, c1, c2 )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        complex( kind=prec ), intent( in ) :: c1, c2
        !
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "linComb_rVector3D_MR > self not allocated." )
        endif
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "linComb_rVector3D_MR > rhs not allocated." )
        endif
        !
        call self%switchStoreState( rhs%store_state )
        !
        if( self%isCompatible( rhs ) ) then
            !
            do i = 1, self%grid%getNGrids()
                !
                if( rhs%store_state .EQ. compound ) then
                    !
                    select type( rhs )
                        !
                        class is( rVector3D_MR_t )
                            !
                            self%sub_vector(i)%x = c1* self%sub_vector(i)%x + c2 * rhs%sub_vector(i)%x
                            self%sub_vector(i)%y = c1* self%sub_vector(i)%y + c2 * rhs%sub_vector(i)%y
                            self%sub_vector(i)%z = c1* self%sub_vector(i)%z + c2 * rhs%sub_vector(i)%z
                            !
                        class is( rScalar3D_MR_t )
                            !
                            self%sub_vector(i)%x = c1* self%sub_vector(i)%x + c2 * rhs%sub_scalar(i)%v
                            self%sub_vector(i)%y = c1* self%sub_vector(i)%y + c2 * rhs%sub_scalar(i)%v
                            self%sub_vector(i)%z = c1* self%sub_vector(i)%z + c2 * rhs%sub_scalar(i)%v
                            !
                        class default
                            call errStop( "linComb_rVector3D_MR > Undefined rhs" )
                            !
                    end select
                    !
                elseif( rhs%store_state .EQ. singleton ) then
                    !
                    select type( rhs )
                        !
                        class is( rVector3D_MR_t )
                            !
                            self%sub_vector(i)%s_v = c1* self%sub_vector(i)%s_v + c2 * rhs%sub_vector(i)%s_v
                            !
                        class is( rScalar3D_MR_t )
                            !
                            self%sub_vector(i)%s_v = c1* self%sub_vector(i)%s_v + c2 * rhs%sub_scalar(i)%s_v
                            !
                        class default
                            call errStop( "linComb_rVector3D_MR > Undefined singleton rhs" )
                            !
                    end select
                    !
                else
                    call errStop( "linComb_rVector3D_MR > Unknow store_state." )
                endif
                !
            enddo
            !
        else
            call errStop( "linComb_rVector3D_MR > Incompatible inputs." )
        endif
        !
    end subroutine linComb_rVector3D_MR
    !
    !> No subroutine briefing
    subroutine subValue_rVector3D_MR( self, cvalue )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        integer :: i
        !
        do i = 1, self%grid%getNGrids()
            !
            if( self%store_state .EQ. compound ) then
                !
                self%sub_vector(i)%x = self%sub_vector(i)%x - cvalue
                self%sub_vector(i)%y = self%sub_vector(i)%y - cvalue
                self%sub_vector(i)%z = self%sub_vector(i)%z - cvalue
                !
            elseif( self%store_state .EQ. singleton ) then
                !
                self%sub_vector(i)%s_v = self%sub_vector(i)%s_v - cvalue
                !
            else
                call errStop( "subValue_rVector3D_MR > Unknown store_state!" )
            endif
            !
        enddo
        !
    end subroutine subValue_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine subField_rVector3D_MR( self, rhs )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "subField_rVector3D_MR > self not allocated." )
        endif
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "subField_rVector3D_MR > rhs not allocated." )
        endif
        !
        call self%switchStoreState( rhs%store_state )
        !
        if( self%isCompatible( rhs ) ) then
            !
            do i = 1, self%grid%getNGrids()
                !
                if( rhs%store_state .EQ. compound ) then
                    !
                    select type( rhs )
                        !
                        class is( rVector3D_MR_t )
                            !
                            self%sub_vector(i)%x = self%sub_vector(i)%x - rhs%sub_vector(i)%x
                            self%sub_vector(i)%y = self%sub_vector(i)%y - rhs%sub_vector(i)%y
                            self%sub_vector(i)%z = self%sub_vector(i)%z - rhs%sub_vector(i)%z
                            !
                        class is( rScalar3D_MR_t )
                            !
                            self%sub_vector(i)%x = self%sub_vector(i)%x - rhs%sub_scalar(i)%v
                            self%sub_vector(i)%y = self%sub_vector(i)%y - rhs%sub_scalar(i)%v
                            self%sub_vector(i)%z = self%sub_vector(i)%z - rhs%sub_scalar(i)%v
                            !
                        class default
                            call errStop( "subField_rVector3D_MR > Undefined rhs" )
                            !
                    end select
                    !
                elseif( rhs%store_state .EQ. singleton ) then
                    !
                    select type( rhs )
                        !
                        class is( rVector3D_MR_t )
                            !
                            self%sub_vector(i)%s_v = self%sub_vector(i)%s_v - rhs%sub_vector(i)%s_v
                            !
                        class is( rScalar3D_MR_t )
                            !
                            self%sub_vector(i)%s_v = self%sub_vector(i)%s_v - rhs%sub_scalar(i)%s_v
                            !
                        class default
                            call errStop( "add_rVector3D_MR > Undefined singleton rhs" )
                            !
                    end select
                    !
                else
                    call errStop( "subField_rVector3D_MR > Unknow store_state." )
                endif
                !
            enddo
            !
        else
            call errStop( "subField_rVector3D_MR > Incompatible inputs." )
        endif
        !
    end subroutine subField_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine multByReal_rVector3D_MR( self, rvalue )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: rvalue
        !
        integer :: i
        !
        do i = 1, self%grid%getNGrids()
            !
            if( self%store_state .EQ. compound ) then
                !
                self%sub_vector(i)%x = self%sub_vector(i)%x * rvalue
                self%sub_vector(i)%y = self%sub_vector(i)%y * rvalue
                self%sub_vector(i)%z = self%sub_vector(i)%z * rvalue
                !
            elseif( self%store_state .EQ. singleton ) then
                !
                self%sub_vector(i)%s_v = self%sub_vector(i)%s_v * rvalue
                !
            else
                call errStop( "multByReal_rVector3D_MR > Unknown store_state!" )
            endif
            !
        enddo
        !
    end subroutine multByReal_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine multByComplex_rVector3D_MR( self, cvalue )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        integer :: i
        !
        do i = 1, self%grid%getNGrids()
            !
            if( self%store_state .EQ. compound ) then
                !
                self%sub_vector(i)%x = self%sub_vector(i)%x * cvalue
                self%sub_vector(i)%y = self%sub_vector(i)%y * cvalue
                self%sub_vector(i)%z = self%sub_vector(i)%z * cvalue
                !
            elseif( self%store_state .EQ. singleton ) then
                !
                self%sub_vector(i)%s_v = self%sub_vector(i)%s_v * cvalue
                !
            else
                call errStop( "multByComplex_rVector3D_MR > Unknown store_state!" )
            endif
            !
        enddo
        !
    end subroutine multByComplex_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine multByField_rVector3D_MR( self, rhs )
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "multByField_rVector3D_MR > self not allocated." )
        endif
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "multByField_rVector3D_MR > rhs not allocated." )
        endif
        !
        call self%switchStoreState( rhs%store_state )
        !
        if( self%isCompatible( rhs ) ) then
            !
            do i = 1, self%grid%getNGrids()
                !
                if( rhs%store_state .EQ. compound ) then
                    !
                    select type( rhs )
                        !
                        class is( rVector3D_MR_t )
                            !
                            self%sub_vector(i)%x = self%sub_vector(i)%x * rhs%sub_vector(i)%x
                            self%sub_vector(i)%y = self%sub_vector(i)%y * rhs%sub_vector(i)%y
                            self%sub_vector(i)%z = self%sub_vector(i)%z * rhs%sub_vector(i)%z
                            !
                        class is( rScalar3D_MR_t )
                            !
                            self%sub_vector(i)%x = self%sub_vector(i)%x * rhs%sub_scalar(i)%v
                            self%sub_vector(i)%y = self%sub_vector(i)%y * rhs%sub_scalar(i)%v
                            self%sub_vector(i)%z = self%sub_vector(i)%z * rhs%sub_scalar(i)%v
                            !
                        class default
                            call errStop( "multByField_rVector3D_MR > Undefined rhs" )
                            !
                    end select
                    !
                elseif( rhs%store_state .EQ. singleton ) then
                    !
                    select type( rhs )
                        !
                        class is( rVector3D_MR_t )
                            !
                            self%sub_vector(i)%s_v = self%sub_vector(i)%s_v * rhs%sub_vector(i)%s_v
                            !
                        class is( rScalar3D_MR_t )
                            !
                            self%sub_vector(i)%s_v = self%sub_vector(i)%s_v * rhs%sub_scalar(i)%s_v
                            !
                        class default
                            call errStop( "add_rVector3D_MR > Undefined singleton rhs" )
                            !
                    end select
                    !
                else
                    call errStop( "multByField_rVector3D_MR > Unknow store_state." )
                endif
                !
            enddo
            !
        else
            call errStop( "multByField_rVector3D_MR > Incompatible inputs." )
        endif
        !
    end subroutine multByField_rVector3D_MR
    !
    !> No subroutine briefing
    !
    function diagMult_rVector3D_MR( self, rhs ) result( diag_mult )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        class( Vector_t ), allocatable :: diag_mult
        !
        integer :: i
        type( rVector3D_MR_t ) :: diag_mult_temp
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "diagMult_rVector3D_MR > self not allocated." )
        endif
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "diagMult_rVector3D_MR > rhs not allocated." )
        endif
        !
        diag_mult_temp = rVector3D_MR_t( self%grid, self%grid_type )
        !
        call diag_mult_temp%switchStoreState( rhs%store_state )
        !
        call self%switchStoreState( rhs%store_state )
        !
        if( self%isCompatible( rhs ) ) then
            !
            do i = 1, self%grid%getNGrids()
                !
                if( rhs%store_state .EQ. compound ) then
                    !
                    select type( rhs )
                        !
                        class is( rVector3D_MR_t )
                            !
                            diag_mult_temp%sub_vector(i)%x = self%sub_vector(i)%x * rhs%sub_vector(i)%x
                            diag_mult_temp%sub_vector(i)%y = self%sub_vector(i)%y * rhs%sub_vector(i)%y
                            diag_mult_temp%sub_vector(i)%z = self%sub_vector(i)%z * rhs%sub_vector(i)%z
                            !
                        class is( rScalar3D_MR_t )
                            !
                            diag_mult_temp%sub_vector(i)%x = self%sub_vector(i)%x * rhs%sub_scalar(i)%v
                            diag_mult_temp%sub_vector(i)%y = self%sub_vector(i)%y * rhs%sub_scalar(i)%v
                            diag_mult_temp%sub_vector(i)%z = self%sub_vector(i)%z * rhs%sub_scalar(i)%v
                            !
                        class default
                            call errStop( "diagMult_rVector3D_MR > Undefined rhs" )
                            !
                    end select
                    !
                elseif( rhs%store_state .EQ. singleton ) then
                    !
                    select type( rhs )
                        !
                        class is( rVector3D_MR_t )
                            !
                            diag_mult_temp%sub_vector(i)%s_v = self%sub_vector(i)%s_v * rhs%sub_vector(i)%s_v
                            !
                        class is( rScalar3D_MR_t )
                            !
                            diag_mult_temp%sub_vector(i)%s_v = self%sub_vector(i)%s_v * rhs%sub_scalar(i)%s_v
                            !
                        class default
                            call errStop( "add_rVector3D_MR > Undefined singleton rhs" )
                            !
                    end select
                    !
                else
                    call errStop( "diagMult_rVector3D_MR > Unknow store_state." )
                endif
                !
            enddo
            !
            allocate( diag_mult, source = diag_mult_temp )
            !
        else
            call errStop( "diagMult_rVector3D_MR > Incompatible inputs." )
        endif
        !
    end function diagMult_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine multAdd_rVector3D_MR( self, cvalue, rhs )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        class( Field_t ), intent( in ) :: rhs
        !
        integer :: i
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "multAdd_rVector3D_MR > rhs not allocated." )
        endif
        !
        call self%switchStoreState( rhs%store_state )
        !
        if( self%isCompatible( rhs ) ) then
            !
            do i = 1, self%grid%getNGrids()
                !
                call self%switchStoreState( rhs%store_state )
                !
                select type( rhs )
                    !
                    class is( rVector3D_MR_t ) 
                        !
                        if( rhs%store_state .EQ. compound ) then
                            !
                            self%sub_vector(i)%x = self%sub_vector(i)%x + cvalue * rhs%sub_vector(i)%x
                            self%sub_vector(i)%y = self%sub_vector(i)%y + cvalue * rhs%sub_vector(i)%y
                            self%sub_vector(i)%z = self%sub_vector(i)%z + cvalue * rhs%sub_vector(i)%z
                            !
                        elseif( rhs%store_state .EQ. singleton ) then
                            !
                            self%sub_vector(i)%s_v = self%sub_vector(i)%s_v + cvalue * rhs%sub_vector(i)%s_v
                            !
                        else
                            call errStop( "multAdd_rVector3D_MR > Unknown rhs store_state!" )
                        endif
                        !
                    class default
                        call errStop( "multAdd_rVector3D_MR > rhs undefined." )
                    !
                end select
                !
            enddo
            !
        else
            call errStop( "multAdd_rVector3D_MR >Incompatible inputs." )
        endif
        !
    end subroutine multAdd_rVector3D_MR
    !
    !> No subroutine briefing
    !
    function dotProd_rVector3D_MR( self, rhs ) result( cvalue )
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
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
                    class is( rVector3D_MR_t )
                        !
                        cvalue = cvalue + self%sub_vector(i)%dotProd( rhs%sub_vector(i) )
                        !
                    class default
                        call errStop( "dotProd_rVector3D_MR > rhs must be rVector3D_MR!" )
                    !
                end select
                !
            else
                call errStop( "dotProd_rVector3D_MR > Incompatible rhs" )
            endif
            !
        enddo
        !
    end function dotProd_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine divByValue_rVector3D_MR( self, cvalue )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        integer :: i
        !
        do i = 1, self%grid%getNGrids()
            !
            if( self%store_state .EQ. compound ) then
                !
                self%sub_vector(i)%x = self%sub_vector(i)%x / cvalue
                self%sub_vector(i)%y = self%sub_vector(i)%y / cvalue
                self%sub_vector(i)%z = self%sub_vector(i)%z / cvalue
                !
            elseif( self%store_state .EQ. singleton ) then
                !
                self%sub_vector(i)%s_v = self%sub_vector(i)%s_v / cvalue
                !
            else
                call errStop( "divByValue_rVector3D_MR > Unknown self store_state!" )
            endif
            !
        enddo
        !
    end subroutine divByValue_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine divByField_rVector3D_MR( self, rhs )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        integer :: i
        !
        if( self%isCompatible( rhs ) ) then
            !
            call self%switchStoreState( rhs%store_state )
            !
            do i = 1, self%grid%getNGrids()
                !
                select type( rhs )
                    !
                    class is( rVector3D_MR_t )
                        !
                        if( rhs%store_state .EQ. compound ) then
                            !
                            self%sub_vector(i)%x = self%sub_vector(i)%x / rhs%sub_vector(i)%x
                            self%sub_vector(i)%y = self%sub_vector(i)%y / rhs%sub_vector(i)%y
                            self%sub_vector(i)%z = self%sub_vector(i)%z / rhs%sub_vector(i)%z
                            !
                        elseif( rhs%store_state .EQ. singleton ) then
                            !
                            self%sub_vector(i)%s_v = self%sub_vector(i)%s_v / rhs%sub_vector(i)%s_v
                            !
                        else
                            call errStop( "divByField_rVector3D_MR > Unknown rhs store_state!" )
                        endif
                        !
                    class is( rScalar3D_MR_t )
                        !
                        if( rhs%store_state .EQ. compound ) then
                            !
                            self%sub_vector(i)%x = self%sub_vector(i)%x / rhs%sub_scalar(i)%v
                            self%sub_vector(i)%y = self%sub_vector(i)%y / rhs%sub_scalar(i)%v
                            self%sub_vector(i)%z = self%sub_vector(i)%z / rhs%sub_scalar(i)%v
                            !
                        elseif( rhs%store_state .EQ. singleton ) then
                            !
                            self%sub_vector(i)%s_v = self%sub_vector(i)%s_v / rhs%sub_scalar(i)%s_v
                            !
                        else
                            call errStop( "divByField_rVector3D_MR > Unknown rhs store_state!" )
                        endif
                        !
                    class default
                        call errStop( "divByField_rVector3D_MR: undefined rhs" )
                    !
                end select
                !
            enddo
            !
        else
            call errStop( "divByField_rVector3D_MR: incompatible rhs" )
        endif
        !
    end subroutine divByField_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine interpFunc_rVector3D_MR( self, location, xyz, interp )
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
        real( kind=prec ), intent( in ) :: location(3)
        character, intent( in ) :: xyz
        class( Vector_t ), allocatable, intent( inout ) :: interp
        !
        call errStop( "interpFunc_rVector3D_MR still not implemented" )
        !
    end subroutine interpFunc_rVector3D_MR
    !
    !> No function briefing
    !
    function getAxis_rVector3D_MR( self, comp_lbl ) result( comp )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        character, intent( in ) :: comp_lbl
        !
        complex( kind=prec ), allocatable :: comp(:,:,:)
        !
        call errStop( "getAxis_rVector3D_MR still not implemented" )
        !
    end function getAxis_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine getReal_rVector3D_MR( self, r_vector )
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
        class( Vector_t ), allocatable, intent( out ) :: r_vector
        !
        allocate( r_vector, source = rVector3D_MR_t( self%grid, self%grid_type ) )
        !
        call r_vector%copyFrom( self )
        !
    end subroutine getReal_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine deallOtherState_rVector3D_MR( self )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        !
        call errStop( "deallOtherState_rVector3D_MR not implemented!" )
        ! !
        ! if( ( .NOT. self%is_allocated ) ) then
            ! call errStop( "deallOtherState_rVector3D_MR > Self not allocated." )
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
            ! call errStop( "deallOtherState_rVector3D_MR > Unknown store_state!" )
        ! endif
        ! !
    end subroutine deallOtherState_rVector3D_MR
    !
    !> No subroutine briefing
    !
    function getArray_rVector3D_MR( self ) result( array )
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
        complex( kind=prec ), allocatable, dimension(:) :: array
        !
        real( kind=prec ), allocatable, dimension(:) :: v_full
        !
        v_full = self%getFullArray()
        !
        array = v_full( self%indActive() )
        !
        deallocate( v_full )
        !
    end function getArray_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setArray_rVector3D_MR( self, array )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), dimension(:), intent( in ) :: array
        !
        complex( kind=prec ), allocatable, dimension(:) :: vFull
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
    end subroutine setArray_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine switchStoreState_rVector3D_MR( self, store_state )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        integer, intent( in ), optional :: store_state
        !
        call errStop( "switchStoreState_rVector3D_MR not implemented!" )
        !
    end subroutine switchStoreState_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine copyFrom_rVector3D_MR( self, rhs )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( .NOT. rhs%is_allocated) then
            call errStop( "copyFrom_rVector3D_MR > rhs not allocated" )
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
            class is( rVector3D_MR_t )
                !
                if( allocated( rhs%sub_vector ) ) then
                    !
                    self%sub_vector = rhs%sub_vector
                    !
                else
                    call errStop( "copyFrom_rVector3D_MR > rhs%sub_vector not allocated" )
                endif
                !
                self%is_allocated = .TRUE.
                !
            class default
                call errStop( "copyFrom_rVector3D_MR > Undefined rhs" )
        end select
        !
    end subroutine copyFrom_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine read_rVector3D_MR( self, funit, ftype )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        integer, intent( in ) :: funit
        character(:), allocatable, intent( in ), optional :: ftype
        !
        call errStop( "read_rVector3D_MR not implemented!" )
        !
    end subroutine read_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine write_rVector3D_MR( self, funit, ftype )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        integer, intent( in ) :: funit
        character(:), allocatable, intent( in ), optional :: ftype
        !
        integer :: i_grid
        !
        write( *, * ) "rVector3D_MR: ", self%nx, self%ny, self%nz
        !
        do i_grid = 1, size( self%sub_vector)
            !
            write( *, * ) "   ", i_grid, ":", self%sub_vector(i_grid)%nx, self%sub_vector(i_grid)%ny, self%sub_vector(i_grid)%nz
            !
        enddo
        !
    end subroutine write_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine print_rVector3D_MR( self, io_unit, title, append )
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
        integer, intent( in ), optional :: io_unit
        character(*), intent( in ), optional :: title
        logical, intent( in ), optional :: append
        !
        type( rVector3D_MR_t ) :: copy
        integer :: i, ix, iy, iz,funit
        !
        copy = self
        !
        call copy%switchStoreState( compound )
        !
        if( present( io_unit ) ) then
            funit = io_unit
        else
            funit = 0
        endif
        !
        write(funit,*) "rVector3D_MR field"
        !
        do i = 1, copy%grid%getNGrids()
            !
            if( present( title ) ) write( funit, * ) title
            !
            write( funit, * ) copy%sub_vector(i)%nx, copy%sub_vector(i)%ny, copy%sub_vector(i)%nz
            !
            write(funit, * ) "x-component",copy%sub_vector(i)%NdX
            do ix = 1, copy%sub_vector(i)%NdX(1)
                 do iy = 1, copy%sub_vector(i)%NdX(2)
                    do iz = 1, copy%sub_vector(i)%NdX(3)
                         if( copy%sub_vector(i)%x( ix, iy, iz ) /= 0 ) then
                            write(funit,*) ix,iy,iz, ":[", copy%sub_vector(i)%x( ix, iy, iz ), "]"
                         endif
                    enddo
                 enddo
            enddo
            !
            write(funit,*) "y-component",copy%sub_vector(i)%NdY
            do ix = 1, copy%sub_vector(i)%NdY(1)
                 do iy = 1, copy%sub_vector(i)%NdY(2)
                    do iz = 1, copy%sub_vector(i)%NdY(3)
                         if( copy%sub_vector(i)%y( ix, iy, iz ) /= 0 ) then
                            write(funit,*) ix,iy,iz, ":[", copy%sub_vector(i)%y( ix, iy, iz ), "]"
                         endif
                    enddo
                 enddo
            enddo
            !
            write(funit,*) "z-component",copy%sub_vector(i)%NdZ
            do ix = 1, copy%sub_vector(i)%NdZ(1)
                 do iy = 1, copy%sub_vector(i)%NdZ(2)
                    do iz = 1, copy%sub_vector(i)%NdZ(3)
                         if( copy%sub_vector(i)%z( ix, iy, iz ) /= 0 ) then
                            write(funit,*) ix,iy,iz, ":[", copy%sub_vector(i)%z( ix, iy, iz ), "]"
                         endif
                    enddo
                 enddo
            enddo
            !
        enddo
        !
    end subroutine print_rVector3D_MR
    !
    !> vec1 and vec2 are both SG cVectors (or rVectors)
    !> vec1 represents layers on top.   vec1 and vec2 have resolutions difffering
    !> by a factor of 2 -- the coarser resolution could be the upper layers or
    !> lower -- the routine can figure this out by comparing grid parameters for
    !> vec1%grid and vec2%grid.
    !> assume for now that vec1 is coarser grid -- then on the interface
    !> (bottom of vec1, top of vec2)  the x/y edges of vec1 are active those of vec2
    !> inactive.   The routine  loops over fine grid x and y edges on the interface:
    !>     1) generates x,y,z coordinates of centers for each interface edge.  In the case
    !>     we consider explicitly here, z will be 0 since these points are at top of the vec2 grid.
    !>     2) change z to z_bottom when the fine gris is on top, and the coarse grid below.
    !>     3) call vec1%interpFunc(x,y,z, ...) to get the (sparse vector)
    !>     interpolation functional, form dot product with vec1, and set
    !>     appropriate edge component of vec2
    !> 
    !> Afer call vec1 is unchanged, and vec_2 has interface x/y edges filled in
    !> by interpolation of  coarse grid x/y edeges from vec1
    !
    subroutine setInactiveEdge_rVector3D_MR( vec_1, vec_2 )
        implicit none
        !
        type( rVector3D_SG_t ), intent( inout ) :: vec_1, vec_2
        !
        integer :: i, j
        real( kind=prec ) :: location(3)
        class( Vector_t ), allocatable :: interp
        type( Grid3D_SG_t ) :: sub_grid_1, sub_grid_2
        !
        sub_grid_1 = vec_1%grid
        sub_grid_2 = vec_2%grid
        !
        !> vec1 is coarser (by factor of two)
        if( vec_2%nx .GT. vec_1%nx ) then
            !
            !> fill in inactive edges (top boundary) of vec_2 first x-edges
            do i = 1, vec_2%nx
                do j = 1, vec_2%ny+1
                    !
                    !> these are supposed to be coordinates of centers of x-edges
                    !> in vec_2, but at bottom of vec_1
                    location(1) = sub_grid_2%x_center(i)
                    location(2) = sub_grid_2%y_edge(j)
                    !
                    !> yes, use vec_1 to compute vertical position
                    location(3) = sub_grid_1%z_edge( vec_1%nz+1 )
                    !
                    !> interpolate to location (in vec_1 grid) and store in top x/y layer of vec_2
                    call vec_1%interpFunc( location, 'x', interp )
                    !
                    vec_2%x(i,j,1) = vec_1%dotProd( interp )
                    !
                enddo
            enddo
            !
            !  then y-edges
            do i = 1, vec_2%nx+1
                do j = 1, vec_2%ny
                    !
                    !> these are supposed to be coordinates of centers of y-edges
                    !> in vec_2, but at bottom of vec_1
                    location(1) = sub_grid_2%x_edge(i)
                    location(2) = sub_grid_2%y_center(j)
                    !
                    !> yes, use vec_1 to compute vertical position
                    location(3) = sub_grid_1%z_edge( vec_1%nz + 1 )
                    !
                    !> interpolate to location (in vec_1 grid) and store in top x/y layer of vec_2
                    call vec_1%interpFunc( location, 'y', interp )
                    !
                    vec_2%y(i,j,1) = vec_1%dotProd( interp )
                    !
                enddo
            enddo
        !
        !> vec_2 is coarser
        else
            !If vec_2 is coarser
            !--> switch vec_2, vec_2 
            location(3) = 0
            !last line of x-edge block is vec_2%x(i,j,vec_2%nz+1) = vec_2%dot(interp)
            !and analaously for y block, etc.
            !
            !> fill in inactive edges (top boundary) of vec_2 first x-edges
            do i = 1, vec_1%nx
                do j = 1, vec_1%ny+1
                    !
                    !> these are supposed to be coordinates of centers of x-edges
                    !> in vec_1, but at bottom of vec_2
                    location(1) = sub_grid_1%x_center(i)
                    location(2) = sub_grid_1%y_edge(j)
                    !
                    !> interpolate to location (in vec_2 grid) and store in top x/y layer of vec_1
                    call vec_2%interpFunc( location, 'x', interp )
                    !
                    vec_1%x( i, j, vec_1%nz+1 ) = vec_2%dotProd( interp )
                    !
                enddo
            enddo
            !
            !  then y-edges
            do i = 1, vec_1%nx+1
                do j = 1, vec_1%ny
                    !
                    !> these are supposed to be coordinates of centers of y-edges
                    !> in vec_1, but at bottom of vec_2
                    location(1) = sub_grid_1%x_edge(i)
                    location(2) = sub_grid_1%y_center(j)
                    !
                    !> interpolate to location (in vec_2 grid) and store in top x/y layer of vec_1
                    call vec_2%interpFunc( location, 'y', interp )
                    !
                    vec_1%x( i, j, vec_1%nz+1 ) = vec_2%dotProd( interp )
                    !
                enddo
            enddo
            !
        endif
        !
    end subroutine setInactiveEdge_rVector3D_MR
    !
    !> Now inputs are an SG vector (vec) and and SG scalar (scalar) and a logical "topCoarser"
    !> the vector is always on the coarser grid -- sitting above or below the finer grid
    !> depending on the value of topCoarser (obviously if .true., the vector is defined on
    !> the upper, coarser subgrid.  In this routine the vector will be modified, using values
    !> from the scalar (not modified)
    !> This is the routine needed for PDEmapping
    !
    subroutine addCellFromAdjacentGrid_rVector3D_MR( CoarseGridVector, FineGridScalar, topCoarser )
        implicit none
        !
        type( rVector3D_SG_t ), intent( inout ) :: CoarseGridVector
        type( rScalar3D_SG_t ), intent( in ) :: FineGridScalar
        logical, intent( in ) :: topCoarser
        !
        integer :: i, j, k, kFine, iFine(2), jFine(2), ii, jj
        !
        if( topCoarser ) then
            !   vertical layer index for fine scalar/vector
            kFine = 1
            k = CoarseGridVector%grid%nz+1
        else
            kFine = FineGridScalar%grid%nz
            k = 1
        endif

        !  modify interior x-edges of coarse grid vector on interface
        do i = 1, CoarseGridVector%grid%nx
            iFine(2) = 2*i
            iFine(1) = iFine(2)-1
            !   exclude boundary of MR grid
            do j = 2, CoarseGridvector%grid%ny
                jFine(1) = (j-1)*2
                jFine(2) = jFine(1)+1
                do ii = 1, 2
                    do jj = 1, 2
                        CoarseGridVector%x(i,j,k) = CoarseGridVector%x(i,j,k)+  &
                        FineGridScalar%v(iFine(ii),jFine(jj),kFine)
                    enddo
                enddo
            enddo
        enddo
        !
        !  modify interior y-edges of coarse grid vector on interface
        do j = 1,CoarseGridVector%grid%ny
            jFine(2) = 2*j
            jFine(1) = jFine(2)-1
            !   exclude boundary of MR grid
            do i = 2,CoarseGridVector%grid%nx
                iFine(1) = (i-1)*2
                iFine(2) = iFine(1)+1
                do ii = 1,2
                    do jj  = 1,2
                        CoarseGridVector%y(i,j,k) = CoarseGridVector%y(i,j,k)+   &
                        FineGridScalar%v(iFine(ii),jFine(jj),kFine)
                    enddo
                enddo
            enddo
        enddo
        !
    end subroutine addCellFromAdjacentGrid_rVector3D_MR
    !
end module rVector3D_MR
!