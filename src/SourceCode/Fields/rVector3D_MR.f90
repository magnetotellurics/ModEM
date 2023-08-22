!
!> This module specializes the abstract Vector3D_real class for
!> real vector fields on a multi-resolution staggered grid.
!
module rVector3D_MR
    !
    use MatUtils
    use rVector3D_SG
    use rScalar3D_MR
    !
    type, extends( Vector_t ) :: rVector3D_MR_t
        !
        type( rVector3D_SG_t ), allocatable :: sub_vectors(:)
        !
        integer, dimension(:), allocatable :: ind_active
        !
        contains
            !
            final :: rVector3D_MR_dtor
            !
            !> MR Routines
            procedure, public :: initializeSub => initializeSub_rVector3D_MR
            !
            procedure, public :: setIndexArrays => setIndexArraysVector3D_MR
            !
            procedure, public :: getArray => getArray_rVector3D_MR
            procedure, public :: setArray => setArray_rVector3D_MR
            !
            procedure, public :: lengthFull => lengthFull_rVector3D_MR
            procedure, public :: findFull => findFull_rVector3D_MR
            !
            procedure, public :: findValue => findValue_rVector3D_MR
            !
            procedure, public :: mrToSg => mrToSg_rVector3D_MR
            procedure, public :: sgToMr => sgToMr_rVector3D_MR
            procedure, public :: sgToMrE0 => sgToMrE0_rVector3D_MR
            !
            !> Boundary operations
            procedure, public :: setAllBoundary => setAllBoundary_rVector3D_MR
            procedure, public :: setOneBoundary => setOneBoundary_rVector3D_MR
            procedure, public :: intBdryIndices => intBdryIndices_rVector3D_MR
            !
            !> Dimensioning operations
            procedure, public :: length => length_rVector3D_MR
            procedure, public :: setVecComponents => setVecComponents_rVector3D_MR
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
            procedure, public :: getX => getX_rVector3D_MR
            procedure, public :: setX => setX_rVector3D_MR
            procedure, public :: getY => getY_rVector3D_MR
            procedure, public :: setY => setY_rVector3D_MR
            procedure, public :: getZ => getZ_rVector3D_MR
            procedure, public :: setZ => setZ_rVector3D_MR
            !
            procedure, public :: getSV => getSV_rVector3D_MR
            procedure, public :: setSV => setSV_rVector3D_MR
            !
            procedure, public :: deallOtherState => deallOtherState_rVector3D_MR
            !
            procedure, public :: getActive => getActive_rVector3D_MR
            procedure, public :: setActive => setActive_rVector3D_MR
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
        module procedure rVector3D_MR_ctor_copy
        module procedure rVector3D_MR_ctor_default
    end interface rVector3D_MR_t
    !
contains
    !
    !> No function briefing
    !
    function rVector3D_MR_ctor_copy( E_in ) result(  self )
        implicit none
        !
        type( rVector3D_MR_t ), intent( in ) :: E_in
        !
        type( rVector3D_MR_t ) :: self
        !
        integer :: i
        !
        self%grid => E_in%grid
        self%grid_type = E_in%grid_type
        !
        call self%initializeSub
        !
        self%ind_active = E_in%ind_active
        self%ind_interior = E_in%ind_interior
        self%ind_boundary = E_in%ind_boundary
        !
        select type( grid => E_in%grid )
            !
            class is( Grid3D_MR_t )
                !
                do i = 1, grid%n_grids
                    self%sub_vectors(i) = E_in%sub_vectors(i)
                enddo
                !
            class default
                call errStop( "rVector3D_MR_ctor_copy > Unclassified grid" )
            !
        end select
        !
    end function rVector3D_MR_ctor_copy
    !
    !> No function briefing
    !
    function rVector3D_MR_ctor_default( grid, grid_type ) result(  self )
        implicit none
        !
        class( Grid_t ), target, intent( in ) :: grid
        character( len=4 ), intent( in ) :: grid_type
        !
        type( rVector3D_MR_t ) :: self
        !
        self%grid => grid
        self%grid_type = grid_type
        !
        write(*,*) "grid_type: ", grid_type
        !
        self%nx = self%grid%nx
        self%ny = self%grid%ny
        self%nz = self%grid%nz
        !
        call self%initializeSub
        !
        if( self%is_allocated ) then
            !
            call self%setIndexArrays
            call self%zeros
            !
        else
            call errStop( "rVector3D_MR_ctor_default > Unable to allocate vector." )
        endif
        !
    end function rVector3D_MR_ctor_default
    !
    !> No subroutine briefing
    !
    subroutine initializeSub_rVector3D_MR( self ) 
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        !
        integer :: i, alloc_stat, ndx, ndy, ndz
        !
        select type( grid => self%grid )
            !
            class is( Grid3D_MR_t )
                !
                self%is_allocated = .TRUE.
                allocate( self%sub_vectors( grid%n_grids ), stat = alloc_stat )
                self%is_allocated = self%is_allocated .AND. (  alloc_stat .EQ. 0 )
                !
                ndx = 0
                ndy = 0
                ndz = 0
                !
                do i = 1, grid%n_grids
                    !
                    self%sub_vectors(i) = rVector3D_SG_t( grid%sub_grids(i), self%grid_type )
                    !
                    write( *, * ) "SubVector", i, "-nx=", self%sub_vectors(i)%nx, ", ny=", self%sub_vectors(i)%ny, "nz=", self%sub_vectors(i)%nz
                    !
                    ndx = ndx + self%sub_vectors(i)%nx / ( 2 ** grid%coarseness( i, 1 ) )
                    ndy = ndy + self%sub_vectors(i)%ny / ( 2 ** grid%coarseness( i, 1 ) )
                    ndz = ndz + self%sub_vectors(i)%nz
                    !
                enddo
                !
                ndx = ndx / self%grid%getNGrids()
                ndy = ndy / self%grid%getNGrids()
                ndz = ndz / self%grid%getNGrids()
                !
                write( *, * ) "MainVector-nx=", self%nx, ", ny=", self%ny, "nz=", self%nz, self%nx*self%ny*self%nz
                !
                if( self%grid_type == EDGE ) then
                    !
                    self%NdX = (/ndx, ndy + 1, ndz + 1/)
                    self%NdY = (/ndx + 1, ndy, ndz + 1/)
                    self%NdZ = (/ndx + 1, ndy + 1, ndz/)
                    !
                elseif( self%grid_type == FACE ) then
                    !
                    self%NdX = (/ndx + 1, ndy, ndz/)
                    self%NdY = (/ndx, ndy + 1, ndz/)
                    self%NdZ = (/ndx, ndy, ndz + 1/)
                    !
                else
                    call errStop( "initializeSub_rVector3D_MR > Only EDGE or FACE types allowed." )
                endif
                !
                write( *, * ) "       NdX=", self%NdX(1), ", ny=", self%NdX(2), "nz=", self%NdX(3), self%NdX(1)*self%NdX(2)*self%NdX(3)
                write( *, * ) "       NdY=", self%NdY(1), ", ny=", self%NdY(2), "nz=", self%NdY(3), self%NdY(1)*self%NdY(2)*self%NdY(3)
                write( *, * ) "       NdZ=", self%NdZ(1), ", ny=", self%Ndz(2), "nz=", self%Ndz(3), self%Ndz(1)*self%Ndz(2)*self%Ndz(3)
                !
                self%Nxyz = (/product(self%NdX), product(self%NdY), product(self%NdZ)/)
                !
                write( *, * ) "      Nxyz=", self%Nxyz(1), ", ny=", self%Nxyz(2), "nz=", self%Nxyz(3)
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
    subroutine setIndexArraysVector3D_MR( self, xy_in ) 
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        logical, intent( in ), optional :: xy_in
        !
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
        select type( grid => self%grid )
            !
            class is( Grid3D_MR_t )
                !
                ! Loop over sub-grids, setting boundary edges to one,
                ! interior to  zero
                do k = 1, grid%n_grids
                    call self%sub_vectors(k)%setAllBoundary( cmplx( 1._prec, 0.0, kind=prec ) )
                enddo
                !
                ! Loop over interfaces: set redundant interface edges to 2
                select case( self%grid_type )
                    !
                    case( EDGE )
                        int_only = .TRUE.
                    case( FACE )
                        int_only = .FALSE.
                    case default
                        !
                        call errStop( "setIndexArraysVector3D_MR > Invalid grid type option!" )
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
                            call self%sub_vectors(k-1)%setOneBoundary( "z2_x", cmplx( -1.0_prec, 0.0, kind=prec ) )
                            call self%sub_vectors(k-1)%setOneBoundary( "z2_y", cmplx( -10.0_prec, 0.0, kind=prec ) )
                        else
                            call self%sub_vectors(k-1)%setOneBoundary( "z2", cmplx( -1.0_prec, 0.0, kind=prec ) )
                        endif
                        !
                        call self%sub_vectors(k)%setOneBoundary( "z1", cmplx( 0._prec, 0.0, kind=prec ), int_only )
                    else
                        if( xy ) then
                            call self%sub_vectors(k)%setOneBoundary( "z1_x", cmplx( -1.0_prec, 0.0, kind=prec ) )
                            call self%sub_vectors(k)%setOneBoundary( "z1_y", cmplx( -10.0_prec, 0.0, kind=prec ) )
                        else
                            call self%sub_vectors(k)%setOneBoundary( "z1", cmplx( -1.0_prec, 0.0, kind=prec ) )
                        endif
                        !
                        call self%sub_vectors(k-1)%setOneBoundary( "z2", cmplx( 0._prec, 0.0, kind=prec ), int_only )
                        !
                    endif
                    !
                enddo
                !
            class default
                call errStop( "setIndexArraysVector3D_MR > Unclassified grid" )
            !
        end select
        !
        ! Set active, interior, and boundary edges. ***
        !
        v_1 = self%getArray()
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
        if( allocated( self%ind_active ) ) then
            deallocate( self%ind_active)
        endif
        !
        allocate( self%ind_active( n_active ) )
        !
        i = 0
        do k = 1, n_full
            if(v_1(k) >= 0) then
                i = i + 1
                self%ind_active(i) = k
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
        v_2 = v_1(self%ind_active)
        !
        if(allocated( self%ind_interior)) then
            deallocate( self%ind_interior)
        endif
        allocate( self%ind_interior(n_interior))
        !
        i = 0
        do k = 1, n_active
            if(v_2(k) == 0) then
                i = i + 1
                self%ind_interior(i) = k
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
        if(allocated( self%ind_boundary)) then
            deallocate( self%ind_boundary)
        endif
        allocate( self%ind_boundary(n_boundaries)) 
        !
        i = 0
        do k = 1, n_active
            if(v_2(k) == 1) then
                i = i + 1
                self%ind_boundary(i) = k
            endif
        enddo
        !
    end subroutine setIndexArraysVector3D_MR
    !
    !> Creates standard( 1-D array) for all sub-scalars,
    !> INCLUDING redundant interface nodes.
    !
    function getArray_rVector3D_MR( self ) result( array )
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
            n = self%sub_vectors(i)%length()
            i2 = i2 + n
            array(i1:i2) = self%sub_vectors(i)%getArray()
            i1 = i1 + n
            !
        enddo
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
                    n = self%sub_vectors(k)%length()
                    i2 = i2 + n
                    call self%sub_vectors(k)%setArray( array(i1:i2) )
                    i1 = i1 + n
                    !
                enddo
                !
            class default
                call errStop( "setArray_rVector3D_MR > Unclassified grid" )
            !
        end select
        !
    end subroutine setArray_rVector3D_MR
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
            n = n + self%sub_vectors(i)%length()
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
        v = self%getArray()
        !
        n_I = 0
        do k = 1, n
            if(v(k) == c) n_I = n_I + 1
        enddo
        !
        allocate( I(n_I))
        !
        n_I = 0
        do k = 1, n
            if(v(k) == c) then
                n_I = n_I + 1
                I(n_I) = k
            endif
        enddo
        !
    end function findFull_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine findValue_rVector3D_MR( self, I, c )
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
        integer, allocatable, intent( out ) :: I(:)
        real( kind=prec ), intent( in ) :: c
        !
        real( kind=prec ), allocatable :: v(:)
        integer :: n, n_I, k
        real( kind=prec ) :: TOL
        !
        TOL = 1E-5
        !
        n = self%length()
        allocate( v(n))
        v = self%getActive()
        !
        n_I = 0
        do k = 1, n
            if( abs( v(k) - c ) / abs(c) <= TOL ) n_I = n_I + 1
        enddo
        !
        allocate( I( n_I ) )
        !
        n_I = 0
        do k = 1, n
            if( abs( v(k) - c ) / abs(c) <= TOL ) then
                n_I = n_I + 1
                I(n_I) = k
            endif
        enddo
        !
    end subroutine findValue_rVector3D_MR
    !
    !> Convert an MR TVector to a full SG, filling in the full fine grid.
    !> Converts MR Vector object to SG
    !> copying from variable resolution sub-grids to completely fill in the
    !> underlying fine grid.
    !
    subroutine mrToSg_rVector3D_MR( self, sg_v, sg_grid )
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
        type( rVector3D_SG_t ), intent( inout ) :: sg_v
        class( Grid_t ), pointer, intent( in ) :: sg_grid
        !
        type( rVector3D_SG_t ) :: temp
        integer :: x_nx, x_ny, x_nz
        integer :: y_nx, y_ny, y_nz
        integer :: z_nx, z_ny, z_nz
        integer :: last, Cs, i1, i2, i, k
        real( kind=prec ) :: w1, w2
        !
        temp = rVector3D_SG_t( sg_grid, self%grid_type )
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
                                temp%x(i:x_nx:Cs, 1:x_ny:Cs, i1:i2+1) = self%sub_vectors(k)%x
                                temp%y(1:y_nx:Cs, i:y_ny:Cs, i1:i2+1) = self%sub_vectors(k)%y
                                !
                                w1 = 1. -( i - 1.)/Cs
                                w2 = 1. - w1
                                !
                                if(i == 1) then
                                    temp%z(1:z_nx:Cs, 1:z_ny:Cs, i1:i2) = self%sub_vectors(k)%z
                                else
                                    last = size(self%sub_vectors(k)%z(:, 1, 1))
                                    temp%z(i:z_nx:Cs, 1:z_ny:Cs, i1:i2) = &
                                    self%sub_vectors(k)%z(1:last-1, :, :) * &
                                    w1 + self%sub_vectors(k)%z(2:last, :, :) * w2
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
                            sg_v%z(:, :, i1:i2)   = sg_v%z(:, :, i1:i2)
                            !
                        enddo
                        !
                    class default
                        call errStop( "mrToSg_rVector3D_MR > Unclassified grid" )
                    !
                end select
                !
            case default
                call errStop( "mrToSg_rVector3D_MR > Unrecognized grid type." )
            !
        end select
        !
    end subroutine mrToSg_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine sgToMr_rVector3D_MR( self, sg_v )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        type( rVector3D_SG_t ), intent( in ) :: sg_v
        !
        class( Grid_t ), pointer :: grid
        real( kind=prec ), pointer, dimension(:) :: tempe
        type( rVector3D_SG_t ) :: templ_r
        real( kind=prec ), pointer, dimension(:) :: templ, temple
        real( kind=prec ), allocatable, dimension(:,:,:) :: lengthx, lengthy
        integer :: sx1, sx2, sx3, sy1, sy2, sy3, s1, s2
        integer :: Cs, i1, i2
        !
        type( rVector3D_SG_t ) :: temp_L, tempEL
        integer :: k, i
        !
        grid => sg_v%grid
        !
        call edgeLength( grid, templ_r )
        templ => null()
        call getRVector(templ_r, templ)
        !
        tempe => null()
        call getRVector(sg_v, tempe)
        !
        allocate(temple(size(templ)))
        temple = templ*tempe
        !
        temp_L = rVector3D_SG_t( grid, EDGE )
        !
        call setRVector(templ, temp_L)
        !
        tempEL = rVector3D_SG_t( grid, EDGE )
        !
        call setRVector(temple, tempEL)
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
                    sx1 = size(self%sub_vectors(k)%x, 1)
                    sx2 = size(self%sub_vectors(k)%x, 2)
                    sx3 = size(self%sub_vectors(k)%x, 3)
                    allocate( lengthx(sx1, sx2, sx3))
                    lengthx = 0.0
                    !
                    sy1 = size(self%sub_vectors(k)%y, 1)
                    sy2 = size(self%sub_vectors(k)%y, 2)
                    sy3 = size(self%sub_vectors(k)%y, 3)
                    allocate( lengthy(sy1, sy2, sy3))
                    lengthy = 0.0
                    !
                    do i = 1, Cs
                        !
                        s1 = size(tempEL%x, 1)
                        s2 = size(tempEL%x, 2)
                        self%sub_vectors(k)%x = self%sub_vectors(k)%x + &
                        tempEL%x(i:s1:Cs, 1:s2:Cs, i1:i2+1)
                        !
                        s1 = size(tempEL%y, 1)
                        s2 = size(tempEL%y, 2)
                        self%sub_vectors(k)%y = self%sub_vectors(k)%y + &
                        tempEL%y(1:s1:Cs,i:s2:Cs, i1:i2+1)
                        !
                        s1 = size(temp_L%x, 1)
                        s2 = size(temp_L%x, 2)
                        lengthx = lengthx + temp_L%x(i:s1:Cs, 1:s2:Cs, i1:i2+1)
                        !
                        s1 = size(temp_L%y, 1)
                        s2 = size(temp_L%y, 2)
                        lengthy = lengthy + temp_L%y(1:s1:Cs, i:s2:Cs, i1:i2+1)
                        !
                    enddo
                    !
                    self%sub_vectors(k)%x = self%sub_vectors(k)%x/lengthx
                    self%sub_vectors(k)%y = self%sub_vectors(k)%y/lengthy
                    !
                    s1 = size(sg_v%z, 1)
                    s2 = size(sg_v%z, 2)
                    self%sub_vectors(k)%z = sg_v%z(1:s1:Cs, 1:s2:Cs, i1:i2)
                    !
                    deallocate( lengthx, lengthy)
                    !
                enddo
                !
            class default
                call errStop( "sgToMr_rVector3D_MR > Unclassified grid" )
            !
        end select
        !
    end subroutine sgToMr_rVector3D_MR
    !
    !> Converts SG Vector object to MR by averaging
    !> used only for e0
    !
    subroutine sgToMrE0_rVector3D_MR(self, sg_v)
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
                        self%sub_vectors(k)%x = self%sub_vectors(k)%x + &
                        sg_v%x(i:x_nx:Cs, 1:x_ny:Cs, i1:i2+1) *    &
                        repMat(grid%Dx(i:last:Cs), &
                        1, &
                        grid%sub_grids(k)%Ny + 1, &
                        grid%sub_grids(k)%Nz + 1, .false.)

                        last = size(grid%Dy)
                        self%sub_vectors(k)%y = self%sub_vectors(k)%y + &
                        sg_v%y(1:y_nx:Cs, i:y_ny:Cs, i1:i2+1) *  & 
                        repMat(grid%Dy(i:last:Cs), &
                        grid%sub_grids(k)%Nx + 1, &
                        1, &
                        grid%sub_grids(k)%Nz + 1, .TRUE.)
                        !
                    enddo
                    !
                    self%sub_vectors(k)%x = self%sub_vectors(k)%x / &
                    repMat(grid%sub_grids(k)%Dx, &
                    1, &
                    grid%sub_grids(k)%Ny + 1, &
                    grid%sub_grids(k)%Nz + 1, .false.)
                    !
                    self%sub_vectors(k)%y = self%sub_vectors(k)%y / &
                    repMat(grid%sub_grids(k)%Dy, &
                    grid%sub_grids(k)%Nx + 1, &
                    1, &
                    grid%sub_grids(k)%Nz + 1, .TRUE.)
                    !
                    self%sub_vectors(k)%z = sg_v%z(1:z_nx:Cs, 1:z_ny:Cs, i1:i2)
                    !
                enddo
                !
            class default
                call errStop( "sgToMrE0_rVector3D_MR > Unclassified grid" )
            !
        end select
        !
    end subroutine sgToMrE0_rVector3D_MR
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
        if( allocated( self%sub_vectors ) ) deallocate( self%sub_vectors )
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
        integer :: i
        !
        call self%switchStoreState( compound )
        !
        do i = 1, self%grid%getNGrids()
            !
            select case( self%grid_type )
                !
                case( EDGE )
                    !
                    self%sub_vectors(i)%x(:,(/1, self%sub_vectors(i)%NdX(2)/), :) = real( cvalue, kind=prec )
                    self%sub_vectors(i)%x(:, :,(/1, self%sub_vectors(i)%NdX(3)/)) = real( cvalue, kind=prec )
                    self%sub_vectors(i)%y((/1, self%sub_vectors(i)%NdY(1)/), :, :) = real( cvalue, kind=prec )
                    self%sub_vectors(i)%y(:, :,(/1, self%sub_vectors(i)%NdY(3)/)) = real( cvalue, kind=prec )
                    self%sub_vectors(i)%z(:,(/1, self%sub_vectors(i)%NdZ(2)/), :) = real( cvalue, kind=prec )
                    self%sub_vectors(i)%z((/1, self%sub_vectors(i)%NdZ(1)/), :, :) = real( cvalue, kind=prec )
                    !
                case( FACE )
                    !
                    self%sub_vectors(i)%x((/1, self%sub_vectors(i)%NdX(1)/), :, :) = real( cvalue, kind=prec )
                    self%sub_vectors(i)%y(:,(/1, self%sub_vectors(i)%NdY(2)/), :) = real( cvalue, kind=prec )
                    self%sub_vectors(i)%z(:, :,(/1, self%sub_vectors(i)%NdZ(3)/)) = real( cvalue, kind=prec )
                    !
                case default
                    call errStop( "setAllBoundary_rVector3D_MR > Invalid grid type." )
            end select
            !
        enddo
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
            select case( self%sub_vectors(i)%grid_type )
                !
                case( EDGE )
                    !
                    if( int_only_p ) then
                        !
                        select case( bdry )
                            !
                            case("x1")
                                self%sub_vectors(i)%z(1, 2:self%sub_vectors(i)%NdZ(2)-1, :) = real( cvalue, kind=prec )
                                self%sub_vectors(i)%y(1, :, 2:self%sub_vectors(i)%NdY(3)-1) = real( cvalue, kind=prec )
                            case("x2")
                                self%sub_vectors(i)%z(self%sub_vectors(i)%NdZ(1), 2:self%sub_vectors(i)%NdZ(2)-1, :) = real( cvalue, kind=prec )
                                self%sub_vectors(i)%y(self%sub_vectors(i)%NdY(1), :, 2:self%sub_vectors(i)%NdY(3)-1) = real( cvalue, kind=prec )
                            case("y1")
                                self%sub_vectors(i)%z(2:self%sub_vectors(i)%NdZ(1)-1, 1, :) = real( cvalue, kind=prec )
                                self%sub_vectors(i)%x(:, 1, 2:self%sub_vectors(i)%NdX(3)-1) = real( cvalue, kind=prec )
                            case("y2")
                                self%sub_vectors(i)%z(2:self%sub_vectors(i)%NdZ(1)-1, self%sub_vectors(i)%NdZ(2), :) = real( cvalue, kind=prec )
                                self%sub_vectors(i)%x(:, self%sub_vectors(i)%NdX(2), 2:self%sub_vectors(i)%NdX(3)-1) = real( cvalue, kind=prec )
                            case("z1")
                                self%sub_vectors(i)%x(:, 2:self%sub_vectors(i)%NdX(2)-1, 1) = real( cvalue, kind=prec )
                                self%sub_vectors(i)%y(2:self%sub_vectors(i)%NdY(1)-1, :, 1) = real( cvalue, kind=prec )
                            case("z2")
                                self%sub_vectors(i)%x(:, 2:self%sub_vectors(i)%NdX(2)-1, self%sub_vectors(i)%NdX(3)) = real( cvalue, kind=prec )
                                self%sub_vectors(i)%y(2:self%sub_vectors(i)%NdY(1)-1, :, self%sub_vectors(i)%NdY(3)) = real( cvalue, kind=prec )
                            case("z1_x")
                                self%sub_vectors(i)%x(:, 2:self%sub_vectors(i)%NdX(2)-1, 1) = real( cvalue, kind=prec )
                            case("z2_x")
                                self%sub_vectors(i)%x(:, 2:self%sub_vectors(i)%NdX(2)-1, self%sub_vectors(i)%NdX(3)) = real( cvalue, kind=prec )
                            case("z1_y")
                                self%sub_vectors(i)%y(2:self%sub_vectors(i)%NdY(1)-1, :, 1) = real( cvalue, kind=prec )
                            case("z2_y")
                                self%sub_vectors(i)%y(2:self%sub_vectors(i)%NdY(1)-1, :, self%sub_vectors(i)%NdY(3)) = real( cvalue, kind=prec )
                            !
                        end select
                        !
                    else
                        !
                        select case( bdry )
                            !
                            case("x1")
                                self%sub_vectors(i)%z(1, :, :) = real( cvalue, kind=prec )
                                self%sub_vectors(i)%y(1, :, :) = real( cvalue, kind=prec )
                            case("x2")
                                self%sub_vectors(i)%z(self%sub_vectors(i)%NdZ(1), :, :) = real( cvalue, kind=prec )
                                self%sub_vectors(i)%y(self%sub_vectors(i)%NdY(1), :, :) = real( cvalue, kind=prec )
                            case("y1")
                                self%sub_vectors(i)%z(:, 1, :) = real( cvalue, kind=prec )
                                self%sub_vectors(i)%x(:, 1, :) = real( cvalue, kind=prec )
                            case("y2")
                                self%sub_vectors(i)%z(:, self%sub_vectors(i)%NdZ(2), :) = real( cvalue, kind=prec )
                                self%sub_vectors(i)%x(:, self%sub_vectors(i)%NdX(2), :) = real( cvalue, kind=prec )
                            case("z1")
                                self%sub_vectors(i)%x(:, :, 1) = real( cvalue, kind=prec )
                                self%sub_vectors(i)%y(:, :, 1) = real( cvalue, kind=prec )
                            case("z2")
                                self%sub_vectors(i)%x(:, :, self%sub_vectors(i)%NdX(3)) = real( cvalue, kind=prec )
                                self%sub_vectors(i)%y(:, :, self%sub_vectors(i)%NdY(3)) = real( cvalue, kind=prec )
                            case("z1_x")
                                self%sub_vectors(i)%x(:, :, 1) = real( cvalue, kind=prec )
                            case("z2_x")
                                self%sub_vectors(i)%x(:, :, self%sub_vectors(i)%NdX(3)) = real( cvalue, kind=prec )
                            case("z1_y")
                                self%sub_vectors(i)%y(:, :, 1) = real( cvalue, kind=prec )
                            case("z2_y")
                                self%sub_vectors(i)%y(:, :, self%sub_vectors(i)%NdY(3)) = real( cvalue, kind=prec )
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
                            self%sub_vectors(i)%x(1, :, :) = real( cvalue, kind=prec )
                        case("x2")
                            self%sub_vectors(i)%x(self%sub_vectors(i)%NdX(1), :, :) = real( cvalue, kind=prec )
                        case("y1")
                            self%sub_vectors(i)%y(:, 1, :) = real( cvalue, kind=prec )
                        case("y2")
                            self%sub_vectors(i)%y(:, self%sub_vectors(i)%NdY(2), :) = real( cvalue, kind=prec )
                        case("z1")
                            self%sub_vectors(i)%z(:, :, 1) = real( cvalue, kind=prec )
                        case("z2")
                            self%sub_vectors(i)%z(:, :, self%sub_vectors(i)%NdZ(3)) = real( cvalue, kind=prec )
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
    subroutine intBdryIndices_rVector3D_MR( self, ind_i, ind_b )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        integer, allocatable, intent( out ) :: ind_i(:), ind_b(:)
        !
        integer :: m, n
        !
        m = size( self%ind_interior )
        n = size( self%ind_boundary )
        !
        allocate( ind_i(m) )
        allocate( ind_b(n) )
        !
        ind_i = self%ind_interior
        ind_b = self%ind_boundary
        !
    end subroutine intBdryIndices_rVector3D_MR
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
        field_length = size( self%ind_active )
        !
    end function length_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setVecComponents_rVector3D_MR( self, xyz, &
                                              xmin, xstep, xmax, &
                                              ymin, ystep, ymax, &
                                              zmin, zstep, zmax, cvalue )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        character, intent( in ) :: xyz
        integer, intent( in ) :: xmin, xstep, xmax
        integer, intent( in ) :: ymin, ystep, ymax
        integer, intent( in ) :: zmin, zstep, zmax
        complex( kind=prec ), intent(  in ) :: cvalue
        !
        integer :: i
        integer :: x1, x2
        integer :: y1, y2
        integer :: z1, z2
        !
        x1 = xmin; x2 = xmax
        y1 = ymin; y2 = ymax
        z1 = zmin; z2 = zmax
        !
        call self%switchStoreState( compound )
        !
        do i = 1, self%grid%getNGrids()
            !
            select case( xyz )
                !
                case("x")
                    !
                    if(xmin == 0) x1 = self%NdX(1)
                    if(xmax <= 0) x2 = self%NdX(1) + xmax
                    !
                    if(ymin == 0) y1 = self%NdX(2)
                    if(ymax <= 0) y2 = self%NdX(2) + ymax
                    !
                    if(zmin == 0) z1 = self%NdX(3)
                    if(zmax <= 0) z2 = self%NdX(3) + zmax
                    !
                    self%sub_vectors(i)%x( x1:x2:xstep, y1:y2:ystep, z1:z2:zstep ) = cvalue
                    !
                case("y")
                    !
                    if(xmin == 0) x1 = self%NdY(1)
                    if(xmax <= 0) x2 = self%NdY(1) + xmax
                    !
                    if(ymin == 0) y1 = self%NdY(2)
                    if(ymax <= 0) y2 = self%NdY(2) + ymax
                    !
                    if(zmin == 0) z1 = self%NdY(3)
                    if(zmax <= 0) z2 = self%NdY(3) + zmax
                    !
                    self%sub_vectors(i)%y( x1:x2:xstep, y1:y2:ystep, z1:z2:zstep ) = cvalue
                    !
                case("z")
                    !
                    if(xmin == 0) x1 = self%NdZ(1)
                    if(xmax <= 0) x2 = self%NdZ(1) + xmax
                    !
                    if(ymin == 0) y1 = self%NdZ(2)
                    if(ymax <= 0) y2 = self%NdZ(2) + ymax
                    !
                    if(zmin == 0) z1 = self%NdZ(3)
                    if(zmax <= 0) z2 = self%NdZ(3) + zmax
                    !
                    self%sub_vectors(i)%z( x1:x2:xstep, y1:y2:ystep, z1:z2:zstep ) = cvalue
                    !
                case default
                    call errStop( "setVecComponents_rVector3D_MR > Invalid xyz argument." )
            end select
            !
        enddo
        !
    end subroutine setVecComponents_rVector3D_MR
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
            call self%sub_vectors(i)%zeros()
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
        call errStop( "sumEdge_rVector3D_MR > Under implementation" )
        ! !
        ! type( rScalar3D_SG_t ) :: cell_out
        ! integer :: i
        ! integer :: x_xend, x_yend, x_zend
        ! integer :: y_xend, y_yend, y_zend
        ! integer :: z_xend, z_yend, z_zend
        ! logical :: is_interior_only
        ! !
        ! if( .NOT. self%is_allocated ) then
             ! call errStop( "sumEdge_rVector3D_MR > self not allocated." )
        ! endif
        ! !
        ! call self%switchStoreState( compound )
        ! !
        ! is_interior_only = .FALSE.
        ! !
        ! if( present( interior_only ) ) is_interior_only = interior_only
        ! !
        ! if( is_interior_only ) then
            ! call self%setAllBoundary( C_ZERO )
        ! endif
        ! !
        ! select type( grid => self%grid )
            ! !
            ! class is( Grid3D_MR_t )
                ! !
                ! do i = 1, grid%n_grids
                    ! !
                    ! cell_out = rScalar3D_SG_t( grid%sub_grids(i), CELL )
                    ! !
                    ! select case( self%grid_type )
                        ! !
                        ! case( EDGE )
                            ! !
                            ! x_xend = size( self%sub_vectors(i)%x, 1 )
                            ! x_yend = size( self%sub_vectors(i)%x, 2 )
                            ! x_zend = size( self%sub_vectors(i)%x, 3 )
                            ! !
                            ! y_xend = size( self%sub_vectors(i)%y, 1 )
                            ! y_yend = size( self%sub_vectors(i)%y, 2 )
                            ! y_zend = size( self%sub_vectors(i)%y, 3 )
                            ! !
                            ! z_xend = size( self%sub_vectors(i)%z, 1 )
                            ! z_yend = size( self%sub_vectors(i)%z, 2 )
                            ! z_zend = size( self%sub_vectors(i)%z, 3 )
                            ! !
                            ! cell_out%v = self%sub_vectors(i)%x(:,1:x_yend-1,1:x_zend-1) + &
                                              ! self%sub_vectors(i)%x(:,2:x_yend,1:x_zend-1)   + &
                                              ! self%sub_vectors(i)%x(:,1:x_yend-1,2:x_zend)   + &
                                              ! self%sub_vectors(i)%x(:,2:x_yend,2:x_zend)     + &
                                              ! self%sub_vectors(i)%y(1:y_xend-1,:,1:y_zend-1) + &
                                              ! self%sub_vectors(i)%y(2:y_xend,:,1:y_zend-1)   + &
                                              ! self%sub_vectors(i)%y(1:y_xend-1,:,2:y_zend)   + &
                                              ! self%sub_vectors(i)%y(2:y_xend,:,2:y_zend)     + &
                                              ! self%sub_vectors(i)%z(1:z_xend-1,1:z_yend-1,:) + &
                                              ! self%sub_vectors(i)%z(2:z_xend,1:z_yend-1,:)   + &
                                              ! self%sub_vectors(i)%z(1:z_xend-1,2:z_yend,:)   + &
                                              ! self%sub_vectors(i)%z(2:z_xend,2:z_yend,:)
                            ! !
                        ! case( FACE )
                            ! !
                            ! x_xend = size( self%sub_vectors(i)%x, 1 )
                            ! y_xend = size( self%sub_vectors(i)%y, 1 )
                            ! z_xend = size( self%sub_vectors(i)%z, 1 )
                            ! !
                            ! cell_out%v = self%sub_vectors(i)%x(1:x_xend-1,:,:) + self%sub_vectors(i)%x(2:x_xend,:,:) + &
                                              ! self%sub_vectors(i)%y(:,1:y_yend-1,:) + self%sub_vectors(i)%y(:,2:y_yend,:) + &
                                              ! self%sub_vectors(i)%z(:,:,1:z_zend-1) + self%sub_vectors(i)%z(:,:,2:z_zend)
                            ! !
                        ! case default
                            ! call errStop( "sumEdge_rVector3D_MR: undefined grid_type" )
                    ! end select
                    ! !
                    ! allocate( cell_out, source = rScalar3D_MR_t( grid, CELL ) )
                    ! !
                    ! select type( cell_out )
                        ! !
                        ! class is( rScalar3D_MR_t )
                            ! !
                            ! cell_out%sub_scalars(i) = cell_out
                            ! !
                        ! class default
                            ! call errStop( "sumEdge_rVector3D_MR: undefined cell_out" )
                        ! !
                    ! end select
                    ! !
                ! enddo
                ! !
            ! class default
                ! call errStop( "setArray_rVector3D_MR > Unclassified grid" )
            ! !
        ! end select
        ! !
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
        ! !
        ! type( rScalar3D_SG_t ) :: cell_h_out_temp, cell_v_out_temp
        ! integer :: i
        ! integer :: x_xend, x_yend, x_zend
        ! integer :: y_xend, y_yend, y_zend
        ! integer :: z_xend, z_yend, z_zend
        ! logical :: is_interior_only
        ! !
        ! if( .NOT. self%is_allocated ) then
             ! call errStop( "sumEdgeVTI_rVector3D_MR > self not allocated." )
        ! endif
        ! !
        ! call self%switchStoreState( compound )
        ! !
        ! is_interior_only = .FALSE.
        ! !
        ! if( present( interior_only ) ) is_interior_only = interior_only
        ! !
        ! if( is_interior_only ) then
            ! call self%setAllBoundary( C_ZERO )
        ! endif
        ! !
        ! select type( grid => self%grid )
            ! !
            ! class is( Grid3D_MR_t )
                ! !
                ! do i = 1, grid%n_grids
                    ! !
                    ! cell_h_out_temp = rScalar3D_SG_t( grid%sub_grids(i), CELL )
                    ! !
                    ! cell_v_out_temp = rScalar3D_SG_t( grid%sub_grids(i), CELL )
                    ! !
                    ! select case( self%sub_vectors(i)%grid_type )
                        ! !
                        ! case( EDGE )
                            ! !
                            ! x_xend = size( self%sub_vectors(i)%x, 1 )
                            ! x_yend = size( self%sub_vectors(i)%x, 2 )
                            ! x_zend = size( self%sub_vectors(i)%x, 3 )
                            ! !
                            ! y_xend = size( self%sub_vectors(i)%y, 1 )
                            ! y_yend = size( self%sub_vectors(i)%y, 2 )
                            ! y_zend = size( self%sub_vectors(i)%y, 3 )
                            ! !
                            ! z_xend = size( self%sub_vectors(i)%z, 1 )
                            ! z_yend = size( self%sub_vectors(i)%z, 2 )
                            ! z_zend = size( self%sub_vectors(i)%z, 3 )
                            ! !
                            ! cell_h_out_temp%v = self%sub_vectors(i)%x(:,1:x_yend-1,1:x_zend-1) + &
                                                ! self%sub_vectors(i)%x(:,2:x_yend,1:x_zend-1)   + &
                                                ! self%sub_vectors(i)%x(:,1:x_yend-1,2:x_zend)   + &
                                                ! self%sub_vectors(i)%x(:,2:x_yend,2:x_zend)     + &
                                                ! self%sub_vectors(i)%y(1:y_xend-1,:,1:y_zend-1) + &
                                                ! self%sub_vectors(i)%y(2:y_xend,:,1:y_zend-1)   + &
                                                ! self%sub_vectors(i)%y(1:y_xend-1,:,2:y_zend)   + &
                                                ! self%sub_vectors(i)%y(2:y_xend,:,2:y_zend)
                            ! !
                            ! cell_v_out_temp%v = self%sub_vectors(i)%z(1:z_xend-1,1:z_yend-1,:) + &
                                                ! self%sub_vectors(i)%z(2:z_xend,1:z_yend-1,:)   + &
                                                ! self%sub_vectors(i)%z(1:z_xend-1,2:z_yend,:)   + &
                                                ! self%sub_vectors(i)%z(2:z_xend,2:z_yend,:)
                            ! !
                        ! case( FACE )
                            ! !
                            ! x_xend = size( self%sub_vectors(i)%x, 1 )
                            ! y_xend = size( self%sub_vectors(i)%y, 1 )
                            ! z_xend = size( self%sub_vectors(i)%z, 1 )
                            ! !
                            ! cell_h_out_temp%v = self%sub_vectors(i)%x(1:x_xend-1,:,:) + self%sub_vectors(i)%x(2:x_xend,:,:) + &
                                                ! self%sub_vectors(i)%y(:,1:y_yend-1,:) + self%sub_vectors(i)%y(:,2:y_yend,:)
                            ! !
                            ! cell_v_out_temp%v = self%sub_vectors(i)%z(:,:,1:z_zend-1) + self%sub_vectors(i)%z(:,:,2:z_zend)
                            ! !
                        ! case default
                            ! call errStop( "sumEdgeVTI_rVector3D_MR: undefined grid_type" )
                        ! !
                    ! end select
                    ! !
                    ! allocate( cell_h_out, source = rScalar3D_MR_t( grid, CELL ) )
                    ! !
                    ! select type( cell_h_out )
                        ! !
                        ! class is( rScalar3D_MR_t )
                            ! !
                            ! cell_h_out%sub_scalars(i) = cell_h_out_temp
                            ! !
                        ! class default
                            ! call errStop( "sumEdgeVTI_rVector3D_MR: undefined cell_h_out" )
                        ! !
                    ! end select
                    ! !
                    ! allocate( cell_v_out, source = rScalar3D_MR_t( grid, CELL ) )
                    ! !
                    ! select type( cell_v_out )
                        ! !
                        ! class is( rScalar3D_MR_t )
                            ! !
                            ! cell_v_out%sub_scalars(i) = cell_v_out_temp
                            ! !
                        ! class default
                            ! call errStop( "sumEdgeVTI_rVector3D_MR: undefined cell_v_out" )
                        ! !
                    ! end select
                    ! !
                ! enddo
                ! !
            ! class default
                ! call errStop( "setArray_rVector3D_MR > Unclassified grid" )
            ! !
        ! end select
        ! !
    end subroutine sumEdgeVTI_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine sumCell_rVector3D_MR( self, cell_in, ptype )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        class( Scalar_t ), intent( in ) :: cell_in
        character(*), intent( in ), optional :: ptype
        !
        complex( kind=prec ), allocatable :: cell_in_v(:,:,:)
        character(10) :: grid_type
        integer :: xend, yend, zend
        integer :: v_xend, v_yend, v_zend
        integer :: i, ix, iy, iz
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
                call self%switchStoreState( compound )
                !
                do i = 1, self%grid%getNgrids()
                    !
                    cell_in_v = cell_in%sub_scalars(i)%getV()
                    !
                    v_xend = size( cell_in_v, 1 )
                    v_yend = size( cell_in_v, 2 )
                    v_zend = size( cell_in_v, 3 )
                    !
                    select case( grid_type )
                        !
                        case( EDGE )
                            !
                            !> for x-components inside the domain
                            do ix = 1, self%sub_vectors(i)%grid%nx
                                do iy = 2, self%sub_vectors(i)%grid%ny
                                    do iz = 2, self%sub_vectors(i)%grid%nz
                                        self%sub_vectors(i)%x(ix, iy, iz) = ( cell_in_v(ix, iy-1, iz-1) + cell_in_v(ix, iy, iz-1) + &
                                        cell_in_v(ix, iy-1, iz) + cell_in_v(ix, iy, iz) ) / 4.0d0
                                    enddo
                                enddo
                            enddo
                            !
                            !> for y-components inside the domain
                            do ix = 2, self%sub_vectors(i)%grid%nx
                                do iy = 1, self%sub_vectors(i)%grid%ny
                                    do iz = 2, self%sub_vectors(i)%grid%nz
                                        self%sub_vectors(i)%y(ix, iy, iz) = ( cell_in_v(ix-1, iy, iz-1) + cell_in_v(ix, iy, iz-1) + &
                                        cell_in_v(ix-1, iy, iz) + cell_in_v(ix, iy, iz) ) / 4.0d0
                                    enddo
                                enddo
                            enddo
                            !
                            !> for z-components inside the domain
                            do ix = 2, self%sub_vectors(i)%grid%nx
                                do iy = 2, self%sub_vectors(i)%grid%ny
                                    do iz = 1, self%sub_vectors(i)%grid%nz
                                        self%sub_vectors(i)%z(ix, iy, iz) = ( cell_in_v(ix-1, iy-1, iz) + cell_in_v(ix-1, iy, iz) + &
                                        cell_in_v(ix, iy-1, iz) + cell_in_v(ix, iy, iz) ) / 4.0d0
                                    enddo
                                enddo
                            enddo
                            !
                        case( FACE )
                            !
                            xend = size(self%sub_vectors(i)%x, 1)
                            self%sub_vectors(i)%x(2:xend-1,:,:) = cell_in_v(1:v_xend-1,:,:) + cell_in_v(2:v_xend,:,:)
                            !
                            yend = size(self%sub_vectors(i)%y, 1)
                            self%sub_vectors(i)%y(:, 2:yend-1, :) = cell_in_v(:, 1:v_yend-1, :) + cell_in_v(:, 2:v_yend, :)
                            !
                            zend = size(self%sub_vectors(i)%z, 1) 
                            self%sub_vectors(i)%z(:, :, 2:zend-1) = cell_in_v(:, :, 1:v_zend-1) + cell_in_v(:, :, 2:v_zend)
                            !
                        case default
                            call errStop( "sumCell_rVector3D_MR: Unknown type" )
                        !
                    end select !type
                    !
                enddo
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
        if( .NOT. self%is_allocated ) then
             call errStop( "sumCellVTI_rVector3D_MR > self not allocated." )
        endif
        !
        if( .NOT. cell_h_in%is_allocated ) then
             call errStop( "sumCellVTI_rVector3D_MR > cell_h_in not allocated." )
        endif
        !
        if( .NOT. cell_v_in%is_allocated ) then
             call errStop( "sumCellVTI_rVector3D_MR > cell_v_in not allocated." )
        endif
        !
        if( index( self%grid_type, CELL ) > 0 ) then
            call errStop( "sumCellVTI_rVector3D_MR > Only CELL type supported." )
        endif
        !
        if( .NOT. present( ptype ) ) then
            grid_type = EDGE
        else
            grid_type = ptype
        endif
        !
        select type( cell_h_in )
            !
            class is( rScalar3D_MR_t )
                !
                select type( cell_v_in )
                    !
                    class is( rScalar3D_MR_t )
                        !
                        call self%switchStoreState( compound )
                        !
                        do i = 1, self%grid%getNgrids()
                            !
                            v_h = cell_h_in%getV()
                            v_v = cell_v_in%getV()
                            !
                            v_xend = size( v_h, 1 )
                            v_yend = size( v_h, 2 )
                            v_zend = size( v_v, 3 )
                            !
                            select case( grid_type )
                                !
                                case( EDGE )
                                    !
                                    !> for x-components inside the domain
                                    do ix = 1, self%sub_vectors(i)%grid%nx
                                        do iy = 2, self%sub_vectors(i)%grid%ny
                                            do iz = 2, self%sub_vectors(i)%grid%nz
                                                self%sub_vectors(i)%x(ix, iy, iz) = ( v_h(ix, iy-1, iz-1) + v_h(ix, iy, iz-1) + &
                                                v_h(ix, iy-1, iz) + v_h(ix, iy, iz) ) / 4.0d0
                                            enddo
                                        enddo
                                    enddo
                                    !
                                    !> for y-components inside the domain
                                    do ix = 2, self%sub_vectors(i)%grid%nx
                                        do iy = 1, self%sub_vectors(i)%grid%ny
                                            do iz = 2, self%sub_vectors(i)%grid%nz
                                                self%sub_vectors(i)%y(ix, iy, iz) = ( v_h(ix-1, iy, iz-1) + v_h(ix, iy, iz-1) + &
                                                v_h(ix-1, iy, iz) + v_h(ix, iy, iz) ) / 4.0d0
                                            enddo
                                        enddo
                                    enddo
                                    !
                                    !> for z-components inside the domain
                                    do ix = 2, self%sub_vectors(i)%grid%nx
                                        do iy = 2, self%sub_vectors(i)%grid%ny
                                            do iz = 1, self%sub_vectors(i)%grid%nz
                                                self%sub_vectors(i)%z(ix, iy, iz) = ( v_v(ix-1, iy-1, iz) + v_v(ix-1, iy, iz) + &
                                                v_v(ix, iy-1, iz) + v_v(ix, iy, iz) ) / 4.0d0
                                            enddo
                                        enddo
                                    enddo
                                    !
                                case( FACE )
                                    !
                                    xend = size(self%sub_vectors(i)%x, 1)
                                    self%sub_vectors(i)%x(2:xend-1,:,:) = v_h(1:v_xend-1,:,:) + v_h(2:v_xend,:,:)
                                    !
                                    yend = size(self%sub_vectors(i)%y, 1)
                                    self%sub_vectors(i)%y(:, 2:yend-1, :) = v_h(:, 1:v_yend-1, :) + v_h(:, 2:v_yend, :)
                                    !
                                    zend = size(self%sub_vectors(i)%z, 1) 
                                    self%sub_vectors(i)%z(:, :, 2:zend-1) = v_v(:, :, 1:v_zend-1) + v_v(:, :, 2:v_zend)
                                    !
                                case default
                                    call errStop( "sumCellVTI_rVector3D_MR: Unknown grid_type" )
                                    !
                            end select
                            !
                        enddo
                        !
                    class default
                        call errStop( "sumCellVTI_rVector3D_MR > Unclassified cell_in_v" )
                    !
                end select
                !
            class default
                call errStop( "sumCellVTI_rVector3D_MR > Unclassified cell_in_h" )
            !
        end select
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
                            self%sub_vectors(i)%x = self%sub_vectors(i)%x + rhs%sub_vectors(i)%x
                            self%sub_vectors(i)%y = self%sub_vectors(i)%y + rhs%sub_vectors(i)%y
                            self%sub_vectors(i)%z = self%sub_vectors(i)%z + rhs%sub_vectors(i)%z
                            !
                        class is( rScalar3D_MR_t )
                            !
                            self%sub_vectors(i)%x = self%sub_vectors(i)%x + rhs%sub_scalars(i)%v
                            self%sub_vectors(i)%y = self%sub_vectors(i)%y + rhs%sub_scalars(i)%v
                            self%sub_vectors(i)%z = self%sub_vectors(i)%z + rhs%sub_scalars(i)%v
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
                            self%sub_vectors(i)%s_v = self%sub_vectors(i)%s_v + rhs%sub_vectors(i)%s_v
                            !
                        class is( rScalar3D_MR_t )
                            !
                            self%sub_vectors(i)%s_v = self%sub_vectors(i)%s_v + rhs%sub_scalars(i)%s_v
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
                            self%sub_vectors(i)%x = c1* self%sub_vectors(i)%x + c2 * rhs%sub_vectors(i)%getX()
                            self%sub_vectors(i)%y = c1* self%sub_vectors(i)%y + c2 * rhs%sub_vectors(i)%getY()
                            self%sub_vectors(i)%z = c1* self%sub_vectors(i)%z + c2 * rhs%sub_vectors(i)%getZ()
                            !
                        class is( rScalar3D_MR_t )
                            !
                            self%sub_vectors(i)%x = c1* self%sub_vectors(i)%x + c2 * rhs%sub_scalars(i)%getV()
                            self%sub_vectors(i)%y = c1* self%sub_vectors(i)%y + c2 * rhs%sub_scalars(i)%getV()
                            self%sub_vectors(i)%z = c1* self%sub_vectors(i)%z + c2 * rhs%sub_scalars(i)%getV()
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
                            self%sub_vectors(i)%s_v = c1* self%sub_vectors(i)%s_v + c2 * rhs%sub_vectors(i)%getSV()
                            !
                        class is( rScalar3D_MR_t )
                            !
                            self%sub_vectors(i)%s_v = c1* self%sub_vectors(i)%s_v + c2 * rhs%sub_scalars(i)%getSV()
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
                self%sub_vectors(i)%x = self%sub_vectors(i)%x - cvalue
                self%sub_vectors(i)%y = self%sub_vectors(i)%y - cvalue
                self%sub_vectors(i)%z = self%sub_vectors(i)%z - cvalue
                !
            elseif( self%store_state .EQ. singleton ) then
                !
                self%sub_vectors(i)%s_v = self%sub_vectors(i)%s_v - cvalue
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
                            self%sub_vectors(i)%x = self%sub_vectors(i)%x - rhs%sub_vectors(i)%getX()
                            self%sub_vectors(i)%y = self%sub_vectors(i)%y - rhs%sub_vectors(i)%getY()
                            self%sub_vectors(i)%z = self%sub_vectors(i)%z - rhs%sub_vectors(i)%getZ()
                            !
                        class is( rScalar3D_MR_t )
                            !
                            self%sub_vectors(i)%x = self%sub_vectors(i)%x - rhs%sub_scalars(i)%getV()
                            self%sub_vectors(i)%y = self%sub_vectors(i)%y - rhs%sub_scalars(i)%getV()
                            self%sub_vectors(i)%z = self%sub_vectors(i)%z - rhs%sub_scalars(i)%getV()
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
                            self%sub_vectors(i)%s_v = self%sub_vectors(i)%s_v - rhs%sub_vectors(i)%s_v
                            !
                        class is( rScalar3D_MR_t )
                            !
                            self%sub_vectors(i)%s_v = self%sub_vectors(i)%s_v - rhs%sub_scalars(i)%s_v
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
                self%sub_vectors(i)%x = self%sub_vectors(i)%x * rvalue
                self%sub_vectors(i)%y = self%sub_vectors(i)%y * rvalue
                self%sub_vectors(i)%z = self%sub_vectors(i)%z * rvalue
                !
            elseif( self%store_state .EQ. singleton ) then
                !
                self%sub_vectors(i)%s_v = self%sub_vectors(i)%s_v * rvalue
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
                self%sub_vectors(i)%x = self%sub_vectors(i)%x * cvalue
                self%sub_vectors(i)%y = self%sub_vectors(i)%y * cvalue
                self%sub_vectors(i)%z = self%sub_vectors(i)%z * cvalue
                !
            elseif( self%store_state .EQ. singleton ) then
                !
                self%sub_vectors(i)%s_v = self%sub_vectors(i)%s_v * cvalue
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
                            self%sub_vectors(i)%x = self%sub_vectors(i)%x * rhs%sub_vectors(i)%getX()
                            self%sub_vectors(i)%y = self%sub_vectors(i)%y * rhs%sub_vectors(i)%getY()
                            self%sub_vectors(i)%z = self%sub_vectors(i)%z * rhs%sub_vectors(i)%getZ()
                            !
                        class is( rScalar3D_MR_t )
                            !
                            self%sub_vectors(i)%x = self%sub_vectors(i)%x * rhs%sub_scalars(i)%getV()
                            self%sub_vectors(i)%y = self%sub_vectors(i)%y * rhs%sub_scalars(i)%getV()
                            self%sub_vectors(i)%z = self%sub_vectors(i)%z * rhs%sub_scalars(i)%getV()
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
                            self%sub_vectors(i)%s_v = self%sub_vectors(i)%s_v * rhs%sub_vectors(i)%s_v
                            !
                        class is( rScalar3D_MR_t )
                            !
                            self%sub_vectors(i)%s_v = self%sub_vectors(i)%s_v * rhs%sub_scalars(i)%s_v
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
                            diag_mult_temp%sub_vectors(i)%x = self%sub_vectors(i)%x * rhs%sub_vectors(i)%x
                            diag_mult_temp%sub_vectors(i)%y = self%sub_vectors(i)%y * rhs%sub_vectors(i)%y
                            diag_mult_temp%sub_vectors(i)%z = self%sub_vectors(i)%z * rhs%sub_vectors(i)%z
                            !
                        class is( rScalar3D_MR_t )
                            !
                            diag_mult_temp%sub_vectors(i)%x = self%sub_vectors(i)%x * rhs%sub_scalars(i)%v
                            diag_mult_temp%sub_vectors(i)%y = self%sub_vectors(i)%y * rhs%sub_scalars(i)%v
                            diag_mult_temp%sub_vectors(i)%z = self%sub_vectors(i)%z * rhs%sub_scalars(i)%v
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
                            diag_mult_temp%sub_vectors(i)%s_v = self%sub_vectors(i)%s_v * rhs%sub_vectors(i)%s_v
                            !
                        class is( rScalar3D_MR_t )
                            !
                            diag_mult_temp%sub_vectors(i)%s_v = self%sub_vectors(i)%s_v * rhs%sub_scalars(i)%s_v
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
                            self%sub_vectors(i)%x = self%sub_vectors(i)%x + cvalue * rhs%sub_vectors(i)%x
                            self%sub_vectors(i)%y = self%sub_vectors(i)%y + cvalue * rhs%sub_vectors(i)%y
                            self%sub_vectors(i)%z = self%sub_vectors(i)%z + cvalue * rhs%sub_vectors(i)%z
                            !
                        elseif( rhs%store_state .EQ. singleton ) then
                            !
                            self%sub_vectors(i)%s_v = self%sub_vectors(i)%s_v + cvalue * rhs%sub_vectors(i)%s_v
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
        cvalue = C_ZERO
        !
        call errStop( "dotProd_rVector3D_MR still not implemented" )
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
                self%sub_vectors(i)%x = self%sub_vectors(i)%x / cvalue
                self%sub_vectors(i)%y = self%sub_vectors(i)%y / cvalue
                self%sub_vectors(i)%z = self%sub_vectors(i)%z / cvalue
                !
            elseif( self%store_state .EQ. singleton ) then
                !
                self%sub_vectors(i)%s_v = self%sub_vectors(i)%s_v / cvalue
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
            select type( rhs )
                !
                class is( rVector3D_MR_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%sub_vectors(i)%x = self%sub_vectors(i)%x / rhs%sub_vectors(i)%x
                        self%sub_vectors(i)%y = self%sub_vectors(i)%y / rhs%sub_vectors(i)%y
                        self%sub_vectors(i)%z = self%sub_vectors(i)%z / rhs%sub_vectors(i)%z
                        !
                    elseif( rhs%store_state .EQ. singleton ) then
                        !
                        self%sub_vectors(i)%s_v = self%sub_vectors(i)%s_v / rhs%sub_vectors(i)%s_v
                        !
                    else
                        call errStop( "divByField_rVector3D_MR > Unknown rhs store_state!" )
                    endif
                    !
                class is( rScalar3D_MR_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%sub_vectors(i)%x = self%sub_vectors(i)%x / rhs%sub_scalars(i)%v
                        self%sub_vectors(i)%y = self%sub_vectors(i)%y / rhs%sub_scalars(i)%v
                        self%sub_vectors(i)%z = self%sub_vectors(i)%z / rhs%sub_scalars(i)%v
                        !
                    elseif( rhs%store_state .EQ. singleton ) then
                        !
                        self%sub_vectors(i)%s_v = self%sub_vectors(i)%s_v / rhs%sub_scalars(i)%s_v
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
        type( rVector3D_SG_t ) :: interp_temp
        real( kind=prec ), allocatable, dimension(:) :: xC, yC, zC
        integer :: ix, iy, iz, i, ii
        real( kind=prec ) :: wx, wy, wz
        logical, dimension(:), allocatable :: tmp
        !
        if( ( .NOT. self%is_allocated ) ) then
            call errStop( "interpFunc_rVector3D_SG > Self not allocated." )
        endif
        !
        do ii = 1, self%grid%getNGrids()
            !
            select case( self%grid_type )
                !
                case( EDGE )
                    !
                    interp_temp = rVector3D_SG_t( self%grid, EDGE )
                    !
                    select case( xyz )
                        !
                        case("x")
                            !
                            allocate(xC(size(self%grid%del_x)))
                            allocate(yC(size(self%grid%dy + 1)))
                            allocate(zC(size(self%grid%dz + 1)))
                            !
                            xC = CumSum(self%grid%del_x)
                            yC = CumSum([0._prec, self%grid%dy])
                            zC = CumSum([0._prec, self%grid%dz])
                            !
                        case("y")
                            !
                            allocate(xC(size(self%grid%dx + 1)))
                            allocate(yC(size(self%grid%del_y)))
                            allocate(zC(size(self%grid%dz)))
                            
                            xC = CumSum([0._prec, self%grid%dx])
                            yC = CumSum([self%grid%del_y])
                            zC = CumSum([0._prec, self%grid%dz])
                            !
                        case("z")
                            !
                            allocate(xC(size(self%grid%dx + 1)))
                            allocate(yC(size(self%grid%dy + 1)))
                            allocate(zC(size(self%grid%del_z)))
                            !
                            xC = CumSum([0._prec, self%grid%dx])
                            yC = CumSum([0._prec, self%grid%dy])
                            zC = CumSum([self%grid%del_z])
                            !
                        case default
                            call errStop( "interpFunc_rVector3D_MR: Unknown xyz" )
                        !
                    end select
                    !
                case( FACE )
                    !
                    interp_temp = rVector3D_SG_t( self%grid, FACE )
                    !
                    select case( xyz )
                        !
                        case( "x" )
                            !
                            allocate(xC(size(self%grid%dx + 1)))
                            allocate(yC(size(self%grid%del_y)))
                            allocate(zC(size(self%grid%del_z)))
                            !
                            xC = CumSum([0._prec, self%grid%dx])
                            yC = CumSum([self%grid%del_y])
                            zC = CumSum([self%grid%del_z])
                            !
                        case( "y" )
                            !
                            allocate(xC(size(self%grid%del_x)))
                            allocate(yC(size(self%grid%dy + 1)))
                            allocate(zC(size(self%grid%del_z)))
                            !
                            xC = CumSum([self%grid%del_x])
                            yC = CumSum([0._prec, self%grid%dy])
                            zC = CumSum([self%grid%del_z])
                            !
                        case( "z" )
                            !
                            allocate(xC(size(self%grid%del_x)))
                            allocate(yC(size(self%grid%del_y)))
                            allocate(zC(size(self%grid%dz + 1)))
                            !
                            xC = CumSum([self%grid%del_x])
                            yC = CumSum([self%grid%del_y])
                            zC = CumSum([0._prec, self%grid%dz])
                            !
                        case default
                            call errStop( "interpFunc_rVector3D_MR: Unknown xyz" )
                        !
                    end select
                    !
                case default
                    call errStop( "interpFunc_rVector3D_MR: Unknown grid_type" )
                !
            end select
            !
            xC = xC + self%grid%ox
            yC = yC + self%grid%oy
            zC = zC - sum(self%grid%dz(1:self%grid%nzAir)) - self%grid%oz
            !
            tmp = location(1) > xC
            !
            ix = size( tmp )
            !
            do i = size( tmp ), 1, -1 
                if(tmp(i)) then
                    ix = i
                    exit
                endif
            enddo
            !
            tmp = location(2) > yC
            !
            iy = size( tmp )
            !
            do i = size( tmp ), 1, -1 
                if(tmp(i)) then
                    iy = i
                    exit
                endif
            enddo
            !
            tmp = location(3) > zC
            !
            iz = size( tmp )
            !
            do i = size( tmp ), 1, -1 
                if(tmp(i)) then
                    iz = i
                    exit
                endif
            enddo
            !
            deallocate( tmp )
            !
            !> ????
            !ix = findloc( location(1) > xC, .TRUE., back = .TRUE., dim = 1 )
            !iy = findloc( location(2) > yC, .TRUE., back = .TRUE., dim = 1 )
            !iz = findloc( location(3) > zC, .TRUE., back = .TRUE., dim = 1 )
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
            select case( xyz )
                !
                case("x")
                    !
                    interp_temp%x(ix,iy,iz) = wx*wy*wz
                    interp_temp%x(ix+1,iy,iz) = (1-wx)*wy*wz
                    interp_temp%x(ix,iy+1,iz) = wx*(1-wy)*wz
                    interp_temp%x(ix,iy,iz+1) = wx*wy*(1-wz)
                    interp_temp%x(ix,iy+1,iz+1) = wx*(1-wy)*(1-wz)
                    interp_temp%x(ix+1,iy,iz+1) = (1-wx)*wy*(1-wz)
                    interp_temp%x(ix+1,iy+1,iz) = (1-wx)*(1-wy)*wz
                    interp_temp%x(ix+1,iy+1,iz+1) = (1-wx)*(1-wy)*(1-wz)
                    !
                case("y")
                    !
                    interp_temp%y(ix,iy,iz) = wx*wy*wz
                    interp_temp%y(ix+1,iy,iz) = (1-wx)*wy*wz
                    interp_temp%y(ix,iy+1,iz) = wx*(1-wy)*wz
                    interp_temp%y(ix,iy,iz+1) = wx*wy*(1-wz)
                    interp_temp%y(ix,iy+1,iz+1) = wx*(1-wy)*(1-wz)
                    interp_temp%y(ix+1,iy,iz+1) = (1-wx)*wy*(1-wz)
                    interp_temp%y(ix+1,iy+1,iz) = (1-wx)*(1-wy)*wz
                    interp_temp%y(ix+1,iy+1,iz+1) = (1-wx)*(1-wy)*(1-wz)
                    !
                case("z")
                    !
                    interp_temp%z(ix,iy,iz) = wx*wy*wz
                    interp_temp%z(ix+1,iy,iz) = (1-wx)*wy*wz
                    interp_temp%z(ix,iy+1,iz) = wx*(1-wy)*wz
                    interp_temp%z(ix,iy,iz+1) = wx*wy*(1-wz)
                    interp_temp%z(ix,iy+1,iz+1) = wx*(1-wy)*(1-wz)
                    interp_temp%z(ix+1,iy,iz+1) = (1-wx)*wy*(1-wz)
                    interp_temp%z(ix+1,iy+1,iz) = (1-wx)*(1-wy)*wz
                    interp_temp%z(ix+1,iy+1,iz+1) = (1-wx)*(1-wy)*(1-wz)
                    !
                case default
                    call errStop( "interpFunc_rVector3D_MR: Unknown xyz" )
                !
            end select !XYZ
            !
            allocate( interp, source = self )
            !
            select type( interp )
                !
                class is( rVector3D_MR_t )
                    !
                    interp%sub_vectors(ii) = interp_temp
                    !
                class default
                    call errStop( "interpFunc_rVector3D_MR: undefined interp" )
                !
            end select
            !
        enddo
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
    !> No interface function briefing
    !
    function getX_rVector3D_MR( self ) result( x )
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
        !
        complex( kind=prec ), allocatable, dimension(:,:,:) :: x
        !
        type( rVector3D_MR_t ) :: temp
        complex( kind=prec ), allocatable, dimension(:) :: array
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "getX_rVector3D_MR > self not allocated." )
        endif
        !
        !> USE getActive() or getArray() ????
        array = self%getArray()
        !
        write( *, * ) "getX_rVector3D_MR: ", self%NdX(1), self%NdX(2), self%NdX(3), &
        self%NdX(1)*self%NdX(2)*self%NdX(3), size( array )
        !
        x = reshape( real( array, kind=prec ), (/self%NdX(1), self%NdX(2), self%NdX(3)/) )
        !
    end function getX_rVector3D_MR
    !
    !> No interface subroutine briefing
    !
    subroutine setX_rVector3D_MR( self, x )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), dimension(:,:,:), intent( in ) :: x
        !
        type( rVector3D_SG_t ) :: temp
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "setX_rVector3D_MR > self not allocated." )
        endif
        !
        temp = rVector3D_SG_t( self%grid, self%grid_type )
        !
        call temp%setX( x )
        !
        call self%sgToMR( temp )
        !
    end subroutine setX_rVector3D_MR
    !
    !> No interface function briefing
    !
    function getY_rVector3D_MR( self ) result( y )
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
        !
        complex( kind=prec ), allocatable, dimension(:,:,:) :: y
        !
        type( rVector3D_MR_t ) :: temp
        complex( kind=prec ), allocatable, dimension(:) :: array
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "getY_rVector3D_MR > self not allocated." )
        endif
        !
        !> USE getActive() or getArray() ????
        array = self%getArray()
        !
        write( *, * ) "getY_rVector3D_MR: ", self%NdY(1), self%NdY(2), self%NdY(3), &
        self%NdY(1)*self%NdY(2)*self%NdY(3), size( array )
        !
        y = reshape( real( array, kind=prec ), (/self%NdY(1), self%NdY(2), self%NdY(3)/) )
        !
    end function getY_rVector3D_MR
    !
    !> No interface subroutine briefing
    !
    subroutine setY_rVector3D_MR( self, y )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), dimension(:,:,:), intent( in ) :: y
        !
        type( rVector3D_SG_t ) :: temp
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "setY_rVector3D_MR > self not allocated." )
        endif
        !
        temp = rVector3D_SG_t( self%grid, self%grid_type )
        !
        call temp%setY( y )
        !
        call self%sgToMR( temp )
        !
    end subroutine setY_rVector3D_MR
    !
    !> No interface function briefing
    !
    function getZ_rVector3D_MR( self ) result( z )
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
        !
        complex( kind=prec ), allocatable, dimension(:,:,:) :: z
        !
        type( rVector3D_MR_t ) :: temp
        complex( kind=prec ), allocatable, dimension(:) :: array
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "getZ_rVector3D_MR > self not allocated." )
        endif
        !
        !> USE getActive() or getArray() ????
        array = self%getArray()
        !
        write( *, * ) "getZ_rVector3D_MR: ", self%NdZ(1), self%NdZ(2), self%NdZ(3), &
        self%NdZ(1)*self%NdZ(2)*self%NdZ(3), size( array )
        !
        z = reshape( real( array, kind=prec ), (/self%NdZ(1), self%NdZ(2), self%NdZ(3)/) )
        !
    end function getZ_rVector3D_MR
    !
    !> No interface subroutine briefing
    !
    subroutine setZ_rVector3D_MR( self, z )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), dimension(:,:,:), intent( in ) :: z
        !
        type( rVector3D_SG_t ) :: temp
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "setZ_rVector3D_MR > self not allocated." )
        endif
        !
        temp = rVector3D_SG_t( self%grid, self%grid_type )
        !
        call temp%setZ( z )
        !
        call self%sgToMR( temp )
        !
    end subroutine setZ_rVector3D_MR
    !
    !> No function briefing
    !
    function getSV_rVector3D_MR( self ) result( s_v )
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
        !
        complex( kind=prec ), allocatable, dimension(:) :: s_v
        !
        call errStop( "getSV_rVector3D_MR not implemented!" )
        !
    end function getSV_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setSV_rVector3D_MR( self, s_v )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), dimension(:), intent( in ) :: s_v
        !
        call errStop( "setSV_rVector3D_MR not implemented!" )
        !
    end subroutine setSV_rVector3D_MR
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
    function getActive_rVector3D_MR( self ) result( array )
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
        complex( kind=prec ), allocatable, dimension(:) :: array
        !
        real( kind=prec ), allocatable :: v_full(:)
        !
        v_full = self%getArray()
        !
        allocate( array( size( self%ind_active ) ) )
        !
        array = v_full( self%ind_active )
        !
        deallocate( v_full )
        !
    end function getActive_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setActive_rVector3D_MR( self, array )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), dimension(:), intent( in ) :: array
        !
        complex( kind=prec ), allocatable, dimension(:) :: vFull
        !
        allocate( vFull( self%lengthFull() ) )
        !
        vFull( self%ind_active ) = array
        !
        call self%setArray( vFull )
        !
        deallocate( vFull )
        !
    end subroutine setActive_rVector3D_MR
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
        integer :: i
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
        if( allocated( rhs%ind_interior ) ) then
            self%ind_interior = rhs%ind_interior
        else
            call errStop( "copyFrom_rVector3D_MR > rhs%ind_interior not allocated" )
        endif
        !
        if( allocated( rhs%ind_boundary ) ) then
            self%ind_boundary = rhs%ind_boundary
        else
            call errStop( "copyFrom_rVector3D_MR > rhs%ind_boundary not allocated" )
        endif
        !
        select type( rhs )
            !
            class is( rVector3D_MR_t )
                !
                if( allocated( rhs%sub_vectors ) ) then
                    !
                    if( allocated( self%sub_vectors ) ) deallocate( self%sub_vectors )
                    allocate( self%sub_vectors( size( rhs%sub_vectors ) ) )
                    !
                else
                    call errStop( "copyFrom_rVector3D_MR > rhs%sub_vectors not allocated" )
                endif
                !
                do i = 1, size( self%sub_vectors )
                    self%sub_vectors(i) = rhs%sub_vectors(i)
                enddo
                !
                if( allocated( rhs%ind_active ) ) then
                    self%ind_active = rhs%ind_active
                else
                    call errStop( "copyFrom_rVector3D_MR > rhs%ind_active not allocated" )
                endif
                !
                self%is_allocated = .TRUE.
                !
            class default
                call errStop( "copyFrom_rVector3D_MR > Undefined rhs" )
        end select
        !
        self%is_allocated = .TRUE.
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
        call errStop( "write_rVector3D_MR not implemented!" )
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
            write( funit, * ) copy%sub_vectors(i)%nx, copy%sub_vectors(i)%ny, copy%sub_vectors(i)%nz
            !
            write(funit, * ) "x-component",copy%sub_vectors(i)%NdX
            do ix = 1, copy%sub_vectors(i)%NdX(1)
                 do iy = 1, copy%sub_vectors(i)%NdX(2)
                    do iz = 1, copy%sub_vectors(i)%NdX(3)
                         if( copy%sub_vectors(i)%x( ix, iy, iz ) /= 0 ) then
                            write(funit,*) ix,iy,iz, ":[", copy%sub_vectors(i)%x( ix, iy, iz ), "]"
                         endif
                    enddo
                 enddo
            enddo
            !
            write(funit,*) "y-component",copy%sub_vectors(i)%NdY
            do ix = 1, copy%sub_vectors(i)%NdY(1)
                 do iy = 1, copy%sub_vectors(i)%NdY(2)
                    do iz = 1, copy%sub_vectors(i)%NdY(3)
                         if( copy%sub_vectors(i)%y( ix, iy, iz ) /= 0 ) then
                            write(funit,*) ix,iy,iz, ":[", copy%sub_vectors(i)%y( ix, iy, iz ), "]"
                         endif
                    enddo
                 enddo
            enddo
            !
            write(funit,*) "z-component",copy%sub_vectors(i)%NdZ
            do ix = 1, copy%sub_vectors(i)%NdZ(1)
                 do iy = 1, copy%sub_vectors(i)%NdZ(2)
                    do iz = 1, copy%sub_vectors(i)%NdZ(3)
                         if( copy%sub_vectors(i)%z( ix, iy, iz ) /= 0 ) then
                            write(funit,*) ix,iy,iz, ":[", copy%sub_vectors(i)%z( ix, iy, iz ), "]"
                         endif
                    enddo
                 enddo
            enddo
            !
        enddo
        !
    end subroutine print_rVector3D_MR
    !
end module rVector3D_MR
