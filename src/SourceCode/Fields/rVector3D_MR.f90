!
!> This file is part of the ModEM modeling and inversion package.
!>
!>LICENSING information
!>
!> Copyright (C) 2020 ModEM research group.
!> Contact: http://
!>
!> GNU General Public License Usage
!> This file may be used under the terms of the GNU
!> General Public License version 3.0 as published by the Free Software
!> Foundation and appearing in the file LICENSE.GPL included in the
!> packaging of this file.  Please review the following information to
!> ensure the GNU General Public License version 3.0 requirements will be
!> met: http://www.gnu.org/copyleft/gpl.html.
!>
!> SUMMARY
!>
!> This module specializes the abstract Vector3D_real class for
!> real vector fields on a multi-resolution staggered grid.
!>
!
module rVector3D_MR
    !
    use MatUtils
    use rVector3D_SG
    use rScalar3D_MR
    !
    type, extends( rVector3D_SG_t ) :: rVector3D_MR_t
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
            procedure, public :: setFull => setFull_rVector3D_MR
            procedure, public :: getFull => getFull_rVector3D_MR
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
            procedure, public :: getArray => getArray_rVector3D_MR
            procedure, public :: setArray => setArray_rVector3D_MR
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
    function rVector3D_MR_ctor_copy( E_in ) result ( self )
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
        self%ind_boundaries = E_in%ind_boundaries
        !
        select type( grid => E_in%grid )
            !
            class is( Grid3D_MR_t )
                !
                do i = 1, grid%n_grids
                    self%sub_vectors(i) = E_in%sub_vectors(i)
                end do
                !
            class default
                stop "Error: rVector3D_MR_ctor_copy > Unclassified grid"
            !
        end select
        !
    end function rVector3D_MR_ctor_copy
    !
    !> No function briefing
    !
    function rVector3D_MR_ctor_default( grid, grid_type ) result ( self )
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
        call self%initializeSub
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
        integer :: i, status
        !
        select type( grid => self%grid )
            !
            class is( Grid3D_MR_t )
                !
                self%is_allocated = .TRUE.
                allocate( self%sub_vectors( grid%n_grids ), STAT = status )
                self%is_allocated = self%is_allocated .AND. ( status .EQ. 0 )
                !
                do i = 1, grid%n_grids
                    self%sub_vectors(i) = rScalar3D_MR_t( grid%sub_grids(i), self%grid_type )
                end do
                !
                call self%setIndexArrays
                call self%zeros
                !
            class default
                stop "Error: initializeSub_rVector3D_MR > Unclassified grid"
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
        if ( .NOT. present( xy_in ) ) then
            xy = .FALSE.
        else
            xy = xy_in
        end if
        !
        select type( grid => self%grid )
            !
            class is( Grid3D_MR_t )
                !
                ! Loop over sub-grids, setting boundary edges to one,
                ! interior to  zero
                do k = 1, grid%n_grids
                    call self%sub_vectors(k)%setAllBoundary( cmplx( 1._prec, 0.0, kind=prec ) )
                end do
                !
                ! Loop over interfaces: set redundant interface edges to 2
                select case( self%grid_type )
                    !
                    case (EDGE)
                        int_only = .TRUE.
                    case (FACE)
                        int_only = .FALSE.
                    case (NODE)
                        int_only = .TRUE.
                    case default
                        !
                        stop "Error: setIndexArraysVector3D_MR > Invalid grid type option!"
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
                            call self%sub_vectors(k-1)%setOneBoundary ("z2_x", cmplx( -1.0_prec, 0.0, kind=prec ) )
                            call self%sub_vectors(k-1)%setOneBoundary ("z2_y", cmplx( -10.0_prec, 0.0, kind=prec ) )
                        else
                            call self%sub_vectors(k-1)%setOneBoundary ("z2", cmplx( -1.0_prec, 0.0, kind=prec ) )
                        endif
                        !
                        call self%sub_vectors(k)%setOneBoundary ("z1", cmplx( 0._prec, 0.0, kind=prec ), int_only )
                    else
                        if( xy ) then
                            call self%sub_vectors(k)%setOneBoundary ("z1_x", cmplx( -1.0_prec, 0.0, kind=prec ) )
                            call self%sub_vectors(k)%setOneBoundary ("z1_y", cmplx( -10.0_prec, 0.0, kind=prec ) )
                        else
                            call self%sub_vectors(k)%setOneBoundary ("z1", cmplx( -1.0_prec, 0.0, kind=prec ) )
                        endif
                        !
                        call self%sub_vectors(k-1)%setOneBoundary ("z2", cmplx( 0._prec, 0.0, kind=prec ), int_only )
                        !
                    end if
                    !
                end do
                !
            class default
                stop "Error: setIndexArraysVector3D_MR > Unclassified grid"
            !
        end select
        !
        ! Set active, interior, and boundary edges. ***
        !
        call self%getFull(v_1)
        !
        n_full = size (v_1)
        !
        n_active = 0
        do k = 1, n_full
            if (v_1(k) >= 0) then
                n_active = n_active + 1
            end if
        end do
        !
        if (allocated (self%ind_active)) then
            deallocate (self%ind_active)
        end if
        allocate (self%ind_active(n_active))
        !
        i = 0
        do k = 1, n_full
            if (v_1(k) >= 0) then
                i = i + 1
                self%ind_active(i) = k
            end if
        end do
        !
        n_interior = 0
        do k = 1, n_full
            if (v_1(k) == 0) then
                n_interior = n_interior + 1
            end if
        end do
        !
        allocate (v_2(n_active))
        v_2 = v_1(self%ind_active)
        !
        if (allocated (self%ind_interior)) then
            deallocate (self%ind_interior)
        end if
        allocate (self%ind_interior(n_interior))
        !
        i = 0
        do k = 1, n_active
            if (v_2(k) == 0) then
                i = i + 1
                self%ind_interior(i) = k
            end if
        end do
        !!
        n_boundaries = 0
        do k = 1, n_active
            if (v_2(k) == 1) then
                n_boundaries = n_boundaries + 1
            end if
        end do
        !
        if (allocated (self%ind_boundaries)) then
            deallocate (self%ind_boundaries)
        end if
        allocate (self%ind_boundaries(n_boundaries)) 
        !
        i = 0
        do k = 1, n_active
            if (v_2(k) == 1) then
                i = i + 1
                self%ind_boundaries(i) = k
            end if
        end do
        !
    end subroutine setIndexArraysVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine setFull_rVector3D_MR( self, v )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
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
                    n = self%sub_vectors(k)%length ()
                    i2 = i2 + n
                    call self%sub_vectors(k)%setArray( cmplx( v(i1:i2), 0.0, kind=prec ) )
                    i1 = i1 + n
                end do
                !
            class default
                stop "Error: setFull_rVector3D_MR > Unclassified grid"
            !
        end select
        !
    end subroutine setFull_rVector3D_MR
    !
    !> Creates standard (1-D array) for all sub-scalars,
    !> INCLUDING redundant interface nodes.
    !
    subroutine getFull_rVector3D_MR( self, v )
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
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
                    n = self%sub_vectors(k)%length()
                    i2 = i2 + n
                    v_temp = self%sub_vectors(k)%getArray()
                    v(i1:i2) = v_temp
                    i1 = i1 + n
                end do
                !
            class default
                stop "Error: getFull_rVector3D_MR > Unclassified grid"
            !
        end select
        !
    end subroutine getFull_rVector3D_MR
    !
    !> No function briefing
    !
    function lengthFull_rVector3D_MR( self ) result( n )
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
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
                    n = n + self%sub_vectors(k)%length()
                end do
                !
            class default
                stop "Error: lengthFull_rVector3D_MR > Unclassified grid"
            !
        end select
        !
    end function lengthFull_rVector3D_MR
    !
    !> No function briefing
    !
    function findFull_rVector3D_MR( self, c ) result( I )
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
        real( kind=prec ), intent (in) :: c
        !
        integer, dimension(:), allocatable :: I
        real( kind=prec ), dimension(:), allocatable :: v
        integer :: n, n_I, k
        !
        n = self%lengthFull()
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
            end if
        end do
        !
    end function findFull_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine findValue_rVector3D_MR(self, I, c)
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
        integer, allocatable, intent( out ) :: I(:)
        real( kind=prec ), intent( in ) :: c
        !
        real( kind=prec ), allocatable :: v(:)
        integer :: n, n_I, k
        real( kind=prec ), parameter :: TOL = 1E-5
        !
        n = self%length()
        allocate (v(n))
        v = self%getArray()
        !
        n_I = 0
        do k = 1, n
            if (abs (v(k) - c)/abs (c) <= TOL) n_I = n_I + 1
        end do
        !
        allocate (I(n_I))
        !
        n_I = 0
        do k = 1, n
            if (abs (v(k) - c)/abs (c) <= TOL) then
                n_I = n_I + 1
                I(n_I) = k
            end if
        end do
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
                            Cs = 2**grid%Coarseness(k, 1)
                            i1 = grid%Coarseness(k, 3)
                            i2 = grid%Coarseness(k, 4)
                            ! Copy  x and y components in x and y directions
                            ! edges that aligned with sub-grid edge.
                            do i = 1, Cs
                                !
                                temp%x(i:x_nx:Cs, 1:x_ny:Cs, i1:i2+1) = self%sub_vectors(k)%x
                                temp%y(1:y_nx:Cs, i:y_ny:Cs, i1:i2+1) = self%sub_vectors(k)%y
                                !
                                w1 = 1. - (i - 1.)/Cs
                                w2 = 1. - w1
                                !
                                if (i == 1) then
                                    temp%z(1:z_nx:Cs, 1:z_ny:Cs, i1:i2) = self%sub_vectors(k)%z
                                else
                                    last = size(self%sub_vectors(k)%z(:, 1, 1))
                                    temp%z(i:z_nx:Cs, 1:z_ny:Cs, i1:i2) = &
                                    self%sub_vectors(k)%z(1:last-1, :, :) * &
                                    w1 + self%sub_vectors(k)%z(2:last, :, :) * w2
                                end if
                                !
                            end do
                            ! edges that subdivide the sub-grid
                            ! interpolate  in y and x directions
                            ! copy/interpolate x in y direction
                            ! copy x and y along x and y directions
                            ! respectively
                            do i = 2, Cs
                                w1 = 1. - (i - 1.)/Cs
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
                            end do
                            !
                            sg_v%x(:, :, i1:i2+1) = sg_v%x(:, :, i1:i2+1) + temp%x(:, :, i1:i2+1)
                            sg_v%y(:, :, i1:i2+1) = sg_v%y(:, :, i1:i2+1) + temp%y(:, :, i1:i2+1)
                            sg_v%z(:, :, i1:i2)   = sg_v%z(:, :, i1:i2)
                            !
                        end do
                        !
                    class default
                        stop "Error: mrToSg_rVector3D_MR > Unclassified grid"
                    !
                end select
                !
            case default
                stop "Error: mrToSg_rVector3D_MR > Unrecognized grid type."
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
        real( kind=prec ), allocatable, dimension(:, :, :) :: lengthx, lengthy
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
                    Cs = 2**grid%Coarseness(k, 1)
                    i1 = grid%Coarseness(k, 3)
                    i2 = grid%Coarseness(k, 4)
                    !
                    sx1 = size(self%sub_vectors(k)%x, 1)
                    sx2 = size(self%sub_vectors(k)%x, 2)
                    sx3 = size(self%sub_vectors(k)%x, 3)
                    allocate (lengthx(sx1, sx2, sx3))
                    lengthx = 0.0
                    !
                    sy1 = size(self%sub_vectors(k)%y, 1)
                    sy2 = size(self%sub_vectors(k)%y, 2)
                    sy3 = size(self%sub_vectors(k)%y, 3)
                    allocate (lengthy(sy1, sy2, sy3))
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
                    end do
                    !
                    self%sub_vectors(k)%x = self%sub_vectors(k)%x/lengthx
                    self%sub_vectors(k)%y = self%sub_vectors(k)%y/lengthy
                    !
                    s1 = size(sg_v%z, 1)
                    s2 = size(sg_v%z, 2)
                    self%sub_vectors(k)%z = sg_v%z(1:s1:Cs, 1:s2:Cs, i1:i2)
                    !
                    deallocate (lengthx, lengthy)
                    !
                end do
                !
            class default
                stop "Error: sgToMr_rVector3D_MR > Unclassified grid"
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
                    Cs = 2**grid%Coarseness(k, 1)
                    i1 = grid%Coarseness(k, 3)
                    i2 = grid%Coarseness(k, 4)
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
                    end do
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
                end do
                !
            class default
                stop "Error: sgToMrE0_rVector3D_MR > Unclassified grid"
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
        if( allocated( self%x ) ) deallocate( self%x )
        if( allocated( self%y ) ) deallocate( self%y )
        if( allocated( self%z ) ) deallocate( self%z )
        if( allocated( self%s_v ) ) deallocate( self%s_v )
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
        call self%switchStoreState( compound )
        !
        select case(self%grid_type)
            !
            case(EDGE)
                !
                self%x(:, (/1, self%NdX(2)/), :) = real( cvalue, kind=prec )
                self%x(:, :, (/1, self%NdX(3)/)) = real( cvalue, kind=prec )
                self%y((/1, self%NdY(1)/), :, :) = real( cvalue, kind=prec )
                self%y(:, :, (/1, self%NdY(3)/)) = real( cvalue, kind=prec )
                self%z(:, (/1, self%NdZ(2)/), :) = real( cvalue, kind=prec )
                self%z((/1, self%NdZ(1)/), :, :) = real( cvalue, kind=prec )
                !
            case(FACE)
                !
                self%x((/1, self%NdX(1)/), :, :) = real( cvalue, kind=prec )
                self%y(:, (/1, self%NdY(2)/), :) = real( cvalue, kind=prec )
                self%z(:, :, (/1, self%NdZ(3)/)) = real( cvalue, kind=prec )
                !
            case default
                stop "Error: setAllBoundary_rVector3D_MR > Invalid grid type."
        end select
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
        !
        call self%switchStoreState( compound )
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
                  select case(bdry)
                      case("x1")
                          self%z(1, 2:self%NdZ(2)-1, :) = real( cvalue, kind=prec )
                          self%y(1, :, 2:self%NdY(3)-1) = real( cvalue, kind=prec )
                      case("x2")
                          self%z(self%NdZ(1), 2:self%NdZ(2)-1, :) = real( cvalue, kind=prec )
                          self%y(self%NdY(1), :, 2:self%NdY(3)-1) = real( cvalue, kind=prec )
                      case("y1")
                          self%z(2:self%NdZ(1)-1, 1, :) = real( cvalue, kind=prec )
                          self%x(:, 1, 2:self%NdX(3)-1) = real( cvalue, kind=prec )
                      case("y2")
                          self%z(2:self%NdZ(1)-1, self%NdZ(2), :) = real( cvalue, kind=prec )
                          self%x(:, self%NdX(2), 2:self%NdX(3)-1) = real( cvalue, kind=prec )
                      case("z1")
                          self%x(:, 2:self%NdX(2)-1, 1) = real( cvalue, kind=prec )
                          self%y(2:self%NdY(1)-1, :, 1) = real( cvalue, kind=prec )
                      case("z2")
                          self%x(:, 2:self%NdX(2)-1, self%NdX(3)) = real( cvalue, kind=prec )
                          self%y(2:self%NdY(1)-1, :, self%NdY(3)) = real( cvalue, kind=prec )
                      case("z1_x")
                          self%x(:, 2:self%NdX(2)-1, 1) = real( cvalue, kind=prec )
                      case("z2_x")
                          self%x(:, 2:self%NdX(2)-1, self%NdX(3)) = real( cvalue, kind=prec )
                      case("z1_y")
                          self%y(2:self%NdY(1)-1, :, 1) = real( cvalue, kind=prec )
                      case("z2_y")
                          self%y(2:self%NdY(1)-1, :, self%NdY(3)) = real( cvalue, kind=prec )
                  end select
                else
                  select case(bdry)
                      case("x1")
                          self%z(1, :, :) = real( cvalue, kind=prec )
                          self%y(1, :, :) = real( cvalue, kind=prec )
                      case("x2")
                          self%z(self%NdZ(1), :, :) = real( cvalue, kind=prec )
                          self%y(self%NdY(1), :, :) = real( cvalue, kind=prec )
                      case("y1")
                          self%z(:, 1, :) = real( cvalue, kind=prec )
                          self%x(:, 1, :) = real( cvalue, kind=prec )
                      case("y2")
                          self%z(:, self%NdZ(2), :) = real( cvalue, kind=prec )
                          self%x(:, self%NdX(2), :) = real( cvalue, kind=prec )
                      case("z1")
                          self%x(:, :, 1) = real( cvalue, kind=prec )
                          self%y(:, :, 1) = real( cvalue, kind=prec )
                      case("z2")
                          self%x(:, :, self%NdX(3)) = real( cvalue, kind=prec )
                          self%y(:, :, self%NdY(3)) = real( cvalue, kind=prec )
                      case("z1_x")
                          self%x(:, :, 1) = real( cvalue, kind=prec )
                      case("z2_x")
                          self%x(:, :, self%NdX(3)) = real( cvalue, kind=prec )
                      case("z1_y")
                          self%y(:, :, 1) = real( cvalue, kind=prec )
                      case("z2_y")
                          self%y(:, :, self%NdY(3)) = real( cvalue, kind=prec )
                  end select
                endif
                !
            case(FACE)
                select case(bdry)
                  case("x1")
                      self%x(1, :, :) = real( cvalue, kind=prec )
                  case("x2")
                      self%x(self%NdX(1), :, :) = real( cvalue, kind=prec )
                  case("y1")
                      self%y(:, 1, :) = real( cvalue, kind=prec )
                  case("y2")
                      self%y(:, self%NdY(2), :) = real( cvalue, kind=prec )
                  case("z1")
                      self%z(:, :, 1) = real( cvalue, kind=prec )
                  case("z2")
                      self%z(:, :, self%NdZ(3)) = real( cvalue, kind=prec )
                end select
            case default
                stop "Error: setOneBoundary_rVector3D_MR > Invalid grid type."
        end select
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
        n = size( self%ind_boundaries )
        !
        allocate( ind_i(m) )
        allocate( ind_b(n) )
        !
        ind_i = self%ind_interior
        ind_b = self%ind_boundaries
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
            &                              xmin, xstep, xmax, &
            &                              ymin, ystep, ymax, &
            &                              zmin, zstep, zmax, rvalue )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        character, intent( in ) :: xyz
        integer, intent( in ) :: xmin, xstep, xmax
        integer, intent( in ) :: ymin, ystep, ymax
        integer, intent( in ) :: zmin, zstep, zmax
        real( kind=prec ), intent ( in ) :: rvalue
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
                self%x(x1:x2:xstep, y1:y2:ystep, z1:z2:zstep) = rvalue
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
                self%y(x1:x2:xstep, y1:y2:ystep, z1:z2:zstep) = rvalue
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
                self%z(x1:x2:xstep, y1:y2:ystep, z1:z2:zstep) = rvalue
                !
            case default
                stop "Error: setVecComponents_rVector3D_MR > Invalid xyz argument."
        end select
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
        select type( grid => self%grid )
            !
            class is( Grid3D_MR_t )
                !
                do i = 1, grid%n_grids
                    call self%sub_vectors(i)%zeros()
                end do
                !
            class default
                stop "Error: zeros_rVector3D_MR > Unclassified grid"
            !
        end select
        !
    end subroutine zeros_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine conjugate_rVector3D_MR( self )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        !
        stop "Error: conjugate_rVector3D_MR: Do not try to conjugate a real vector!"
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
        if( .NOT. rhs%is_allocated) then
             stop "Error: add_rVector3D_MR > rhs not allocated."
        endif
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
                        self%x = self%x + rhs%x
                        self%y = self%y + rhs%y
                        self%z = self%z + rhs%z
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%s_v = self%s_v + rhs%s_v
                        !
                    else
                        stop "Error: add_rVector3D_MR > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: add_rVector3D_MR > Undefined rhs"
                !
            end select
            !
        else
            stop "Error: add_rVector3D_MR > Incompatible inputs."
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
        if( self%isCompatible( rhs ) ) then
            !
            call self%switchStoreState( rhs%store_state )
            !
            select type(rhs)
                !
                class is( rVector3D_MR_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%x = c1 * self%x + c2 * rhs%x
                        self%y = c1 * self%y + c2 * rhs%y
                        self%z = c1 * self%z + c2 * rhs%z
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%s_v = c1 * self%s_v + c2 * rhs%s_v
                        !
                    else
                        stop "Error: linComb_rVector3D_MR > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: linComb_rVector3D_MR > rhs undefined."
            end select
            !
        else
            stop "Error: linComb_rVector3D_MR > Incompatible inputs."
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
        if( self%store_state .EQ. compound ) then
            !
            self%x = self%x - cvalue
            self%y = self%y - cvalue
            self%z = self%z - cvalue
            !
        else if( self%store_state .EQ. singleton ) then
            !
            self%s_v = self%s_v - cvalue
            !
        else
            stop "Error: subValue_rVector3D_MR > Unknown self store_state!"
        endif
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
                        self%x = self%x - rhs%x
                        self%y = self%y - rhs%y
                        self%z = self%z - rhs%z
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%s_v = self%s_v - rhs%s_v
                        !
                    else
                        stop "Error: subField_rVector3D_MR > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: subField_rVector3D_MR > Undefined rhs"
                !
            end select
            !
        else
            stop "Error: subField_rVector3D_MR > Incompatible inputs."
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
        if( self%store_state .EQ. compound ) then
            !
            self%x = self%x * rvalue
            self%y = self%y * rvalue
            self%z = self%z * rvalue
            !
        else if( self%store_state .EQ. singleton ) then
            !
            self%s_v = self%s_v * rvalue
            !
        else
            stop "Error: multByReal_rVector3D_MR > Unknown self store_state!"
        endif
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
        if( self%store_state .EQ. compound ) then
            !
            self%x = self%x * cvalue
            self%y = self%y * cvalue
            self%z = self%z * cvalue
            !
        else if( self%store_state .EQ. singleton ) then
            !
            self%s_v = self%s_v * cvalue
            !
        else
            stop "Error: multByComplex_rVector3D_MR > Unknown self store_state!"
        endif
        !
    end subroutine multByComplex_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine multByField_rVector3D_MR( self, rhs )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
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
                        self%x = self%x * rhs%x
                        self%y = self%y * rhs%y
                        self%z = self%z * rhs%z
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%s_v = self%s_v * rhs%s_v
                        !
                    else
                        stop "Error: multByField_rVector3D_MR > Unknown rhs store_state!"
                    endif
                    !
                class is( rScalar3D_MR_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%x = self%x * rhs%v
                        self%y = self%y * rhs%v
                        self%z = self%z * rhs%v
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%s_v = self%s_v * rhs%s_v
                        !
                    else
                        stop "Error: multByField_rVector3D_MR > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: multByField_rVector3D_MR: undefined rhs"
                !
            end select
            !
        else
            stop "Error: multByField_rVector3D_MR: incompatible rhs"
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
        class( Vector_t ), intent( in ) :: rhs
        !
        class( Vector_t ), allocatable :: diag_mult
        !
        if( self%isCompatible( rhs ) ) then
            !
            allocate( diag_mult, source = rVector3D_MR_t( self%grid, self%grid_type ) )
            !
            call diag_mult%switchStoreState( rhs%store_state )
            !
            call self%switchStoreState( rhs%store_state )
            !
            select type( diag_mult )
                class is( rVector3D_MR_t )
                    !
                    select type( rhs )
                        !
                        class is( rVector3D_MR_t )
                            !
                            if( rhs%store_state .EQ. compound ) then
                                !
                                diag_mult%x = self%x * rhs%x
                                diag_mult%y = self%y * rhs%y
                                diag_mult%z = self%z * rhs%z
                                !
                            else if( rhs%store_state .EQ. singleton ) then
                                !
                                diag_mult%s_v = self%s_v * rhs%s_v
                                !
                            else
                                stop "Error: diagMult_rVector3D_MR > Unknown rhs store_state!"
                            endif
                            !
                        class default
                            stop "Error: diagMult_rVector3D_MR > Undefined rhs"
                        !
                    end select
                !
                class default
                    stop "Error: diagMult_rVector3D_MR > Undefined diag_mult"
                !
            end select
            !
        else
            stop "Error: diagMult_rVector3D_MR > Incompatible inputs."
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
        if( self%isCompatible( rhs ) ) then
            !
            call self%switchStoreState( rhs%store_state )
            !
            select type(rhs)
                !
                class is( rVector3D_MR_t ) 
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%x = self%x + cvalue * rhs%x
                        self%y = self%y + cvalue * rhs%y
                        self%z = self%z + cvalue * rhs%z
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%s_v = self%s_v + cvalue * rhs%s_v
                        !
                    else
                        stop "Error: multAdd_rVector3D_MR > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: multAdd_rVector3D_MR > rhs undefined."
                !
            end select
            !
        else
            stop "Error: multAdd_rVector3D_MR >Incompatible inputs."
        endif
        !
    end subroutine multAdd_rVector3D_MR
    !
    !> No subroutine briefing
    !
    function dotProd_rVector3D_MR( self, rhs ) result( cvalue )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        complex( kind=prec ) :: cvalue
        !
        cvalue = C_ZERO
        !
        if(( .NOT. self%is_allocated ) .OR. ( .NOT. rhs%is_allocated )) then
            stop "Error: dotProd_rVector3D_MR > Input vectors not allocated."
        endif
        !
        if( self%isCompatible( rhs ) ) then
            !
            if( self%store_state == rhs%store_state ) then
                !
                select type( rhs )
                    !
                    class is( rVector3D_MR_t )
                        !
                        if( rhs%store_state .EQ. compound ) then
                            !
                            cvalue = cvalue + cmplx( sum( self%x * rhs%x ), 0.0, kind=prec )
                            cvalue = cvalue + cmplx( sum( self%y * rhs%y ), 0.0, kind=prec )
                            cvalue = cvalue + cmplx( sum( self%z * rhs%z ), 0.0, kind=prec )
                            !
                        else if( rhs%store_state .EQ. singleton ) then
                            !
                            cvalue = cvalue + cmplx( sum( self%s_v * rhs%s_v ), 0.0, kind=prec )
                            !
                        else
                            stop "Error: dotProd_rVector3D_MR > Unknown rhs store_state!"
                        endif
                        !
                    class default
                        stop "Error: dotProd_rVector3D_MR: undefined rhs"
                    !
                end select
                !
            else
                stop "Error: dotProd_rVector3D_MR > Incompatible store_state"
            endif
            !
        else
            stop "Error: dotProd_rVector3D_MR > Incompatible rhs"
        endif
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
        if( self%store_state .EQ. compound ) then
            !
            self%x = self%x / cvalue
            self%y = self%y / cvalue
            self%z = self%z / cvalue
            !
        else if( self%store_state .EQ. singleton ) then
            !
            self%s_v = self%s_v / cvalue
            !
        else
            stop "Error: divByValue_rVector3D_MR > Unknown self store_state!"
        endif
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
                        self%x = self%x / rhs%x
                        self%y = self%y / rhs%y
                        self%z = self%z / rhs%z
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%s_v = self%s_v / rhs%s_v
                        !
                    else
                        stop "Error: divByField_rVector3D_MR > Unknown rhs store_state!"
                    endif
                    !
                class is( rScalar3D_MR_t )
                    !
                    if( rhs%store_state .EQ. compound ) then
                        !
                        self%x = self%x / rhs%v
                        self%y = self%y / rhs%v
                        self%z = self%z / rhs%v
                        !
                    else if( rhs%store_state .EQ. singleton ) then
                        !
                        self%s_v = self%s_v / rhs%s_v
                        !
                    else
                        stop "Error: divByField_rVector3D_MR > Unknown rhs store_state!"
                    endif
                    !
                class default
                    stop "Error: divByField_rVector3D_MR: undefined rhs"
                !
            end select
            !
        else
            stop "Error: divByField_rVector3D_MR: incompatible rhs"
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
        stop "Error: interpFunc_rVector3D_MR still not implemented"
        !
    end subroutine interpFunc_rVector3D_MR
    !
    !> No function briefing
    !
    function getAxis_rVector3D_MR( self, comp_lbl ) result( comp )
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
        character, intent( in ) :: comp_lbl
        !
        complex( kind=prec ), allocatable :: comp(:, :, :)
        !
        stop "Error: getAxis_rVector3D_MR still not implemented"
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
        class( rVector3D_MR_t ), intent( inout ) :: self
        !
        complex( kind=prec ), allocatable :: x(:, :, :)
        !
        x = self%x
        !
    end function getX_rVector3D_MR
    !
    !> No interface subroutine briefing
    !
    subroutine setX_rVector3D_MR( self, x )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), allocatable, intent( in ) :: x(:, :, :)
        !
        self%x = x
        !
    end subroutine setX_rVector3D_MR
    !
    !> No interface function briefing
    !
    function getY_rVector3D_MR( self ) result( y )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        !
        complex( kind=prec ), allocatable :: y(:, :, :)
        !
        y = self%y
        !
    end function getY_rVector3D_MR
    !
    !> No interface subroutine briefing
    !
    subroutine setY_rVector3D_MR( self, y )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), allocatable, intent( in ) :: y(:, :, :)
        !
        self%y = y
        !
    end subroutine setY_rVector3D_MR
    !
    !> No interface function briefing
    !
    function getZ_rVector3D_MR( self ) result( z )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        !
        complex( kind=prec ), allocatable :: z(:, :, :)
        !
        z = self%z
        !
    end function getZ_rVector3D_MR
    !
    !> No interface subroutine briefing
    !
    subroutine setZ_rVector3D_MR( self, z )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        complex( kind=prec ), allocatable, intent( in ) :: z(:, :, :)
        !
        self%z = z
        !
    end subroutine setZ_rVector3D_MR
    !
    !> No function briefing
    !
    function getSV_rVector3D_MR( self ) result( s_v )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
        !
        complex( kind=prec ), allocatable :: s_v(:)
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
        complex( kind=prec ), allocatable, intent( in ) :: s_v(:)
        !
        call errStop( "setSV_rVector3D_MR not implemented!" )
        !
    end subroutine setSV_rVector3D_MR
    !
    !> No subroutine briefing
    !
    function getArray_rVector3D_MR( self ) result( array )
        implicit none
        !
        class( rVector3D_MR_t ), intent( in ) :: self
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
        integer :: i
        !
        if( .NOT. rhs%is_allocated) then
            stop "Error: copyFrom_rVector3D_MR > rhs not allocated"
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
        if( allocated( rhs%ind_boundaries ) ) then
            self%ind_boundaries = rhs%ind_boundaries
        else
            call errStop( "copyFrom_rVector3D_MR > rhs%ind_boundaries not allocated" )
        endif
        !
        select type( rhs )
            !
            class is( rVector3D_MR_t )
                !
                if( allocated( rhs%ind_active ) ) then
                    self%ind_active = rhs%ind_active
                else
                    call errStop( "copyFrom_rVector3D_MR > rhs%ind_active not allocated" )
                endif
                !
                self%NdX = rhs%NdX
                self%NdY = rhs%NdY
                self%NdZ = rhs%NdZ
                self%Nxyz = rhs%Nxyz
                !
                if( rhs%store_state .EQ. compound ) then
                    !
                    self%x = rhs%x
                    self%y = rhs%y
                    self%z = rhs%z
                    !
                else if( rhs%store_state .EQ. singleton ) then
                    !
                    self%s_v = rhs%s_v
                    !
                else
                    stop "Error: copyFrom_rVector3D_MR > Unknown store_state!"
                endif
                !
                do i = 1, size( self%sub_vectors )
                    self%sub_vectors(i) = rhs%sub_vectors(i)
                end do
                !
                self%is_allocated = .TRUE.
                !
            class default
                stop "Error: copyFrom_rVector3D_MR > Undefined rhs"
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
        integer :: Nx, Ny, Nz
        character(4) :: grid_type
        logical :: ok, hasname, binary
        character(80) :: fname, isbinary
        !
        call self%switchStoreState( compound )
        !
        binary = .TRUE.
        !
        inquire( funit, opened = ok, named = hasname, name = fname, unformatted = isbinary )
        !
        if( ok ) then
            !
            !> Check that the file is unformatted if binary, formatted if ascii.
            if((index(isbinary, "yes") > 0 .OR.index(isbinary, "YES") > 0) &
                  .AND.   .NOT. binary ) then
                write( *, * ) "Error: read_rVector3D_MR > Unable to read_rVector3D_MR vector from unformatted file. ", &
                        trim(fname), "."
                stop
            else if((index(isbinary, "no") > 0 .OR.index(isbinary, "NO") > 0) &
                  .AND.binary) then
                write( *, * ) "Error: read_rVector3D_MR > Unable to read_rVector3D_MR vector from formatted file ", &
                        trim(fname), "."
                stop
            endif
            !
            read(funit) Nx, Ny, Nz, grid_type
            !
            if(  .NOT. self%is_allocated) then
                write( *, * ) "Error: read_rVector3D_MR > Vector must be allocated before read_rVector3D_MRing from ", &
                        trim(fname), "."
                stop
            else if(self%grid_type.NE.grid_type) then
                write( *, * ) "Error: read_rVector3D_MR > Vector must be of type ", grid_type, &
                        &            "           before read_rVector3D_MRing from ", trim (fname), "."
                stop
            else if((self%nx.NE.Nx).OR. &
                  (self%ny.NE.Ny).OR.(self%nz.NE.Nz)) then
                write( *, * ) "Error: read_rVector3D_MR > Wrong size of vector on input from ", trim (fname), "."
                stop
            endif
            !
            read(funit) self%x
            read(funit) self%y
            read(funit) self%z
            !
        else
            stop "Error: read_rVector3D_MR: unable to open file"
        endif
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
        logical :: ok, hasname, binary
        character(80) :: fname, isbinary
        !
        if( .NOT. self%is_allocated) then
            stop "Error: write_rVector3D_MR > Not allocated."
        endif
        !
        call self%switchStoreState( compound )
        !
        binary = .TRUE.
        !
        inquire( funit, opened = ok, named = hasname, name = fname, unformatted = isbinary )
        !
        if( ok ) then
            !
            !> Check that the file is unformatted if binary, formatted if ascii.
            if((index(isbinary, "yes") > 0 .OR. index(isbinary, "YES") > 0) &
                  .AND.   .NOT. binary) then
                write( *, * ) "Error: write_rVector3D_MR > Unable to write_rVector3D_MR vector to unformatted file. ", &
                        trim(fname), "."
                stop
            else if((index(isbinary,"no") > 0.OR.index(isbinary,"NO") > 0) &
                  .AND.binary) then
                write( *, * ) "Error: write_rVector3D_MR > Unable to write_rVector3D_MR vector to formatted file. ", &
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
            stop "Error: write_rVector3D_MR > unable to open file"
        endif
        !
    end subroutine write_rVector3D_MR
    !
    !> No subroutine briefing
    !
    subroutine print_rVector3D_MR( self, io_unit, title, append )
        implicit none
        !
        class( rVector3D_MR_t ), intent( inout ) :: self
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
            funit = 0
        endif
        !
        write( funit, * ) "ModEM-OO rVector3D_MR"
        if( present( title ) ) write( funit, * ) title
        !
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
    end subroutine print_rVector3D_MR
    !
end module rVector3D_MR
