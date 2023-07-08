!
!> Basic model parameter class --
!> cell conductivities defined on numerical grid
!> in Earth cells only -- as coded this is specific
!> to usual CSG.
!
module ModelParameterCell_SG_VTI
    !
    use ModelParameterCell_SG
    !
    use Constants
    use FileUnits
    use rScalar3D_SG
    use rVector3D_SG 
    use Grid3D_SG
    use ModelParameter1D
    use ModelParameter2D
    !
    type, extends( ModelParameterCell_SG_t ) :: ModelParameterCell_SG_VTI_t
        !
        ! No derived properties
        !
        contains
            !
            procedure, public :: dotProd => dotProd_ModelParameterCell_SG_VTI
            !
            procedure, public :: linComb => linComb_ModelParameterCell_SG_VTI
            !
            procedure, public :: PDEmapping => PDEmapping_ModelParameterCell_SG_VTI
            procedure, public :: dPDEmapping => dPDEmapping_ModelParameterCell_SG_VTI
            procedure, public :: dPDEmapping_T => dPDEmapping_T_ModelParameterCell_SG_VTI
            !
            procedure, public :: write => write_ModelParameterCell_SG_VTI
            !
    end type ModelParameterCell_SG_VTI_t
    !
    interface ModelParameterCell_SG_VTI_t
         module procedure ModelParameterCell_SG_VTI_ctor
    end interface ModelParameterCell_SG_VTI_t
    !
contains
    !
    !> No subroutine briefing
    !
    function ModelParameterCell_SG_VTI_ctor( grid, cell_cond, param_type ) result( self )
        implicit none
        !
        class( Grid_t ), intent( in ) :: grid
        class( Scalar_t ), intent( in ) :: cell_cond
        character(:), allocatable, optional, intent( in ) :: param_type
        !
        type( ModelParameterCell_SG_VTI_t ) :: self
        !
        integer :: nzAir
        !
        !write( *, * ) "Constructor ModelParameterCell_SG_VTI_t"
        !
        call self%init
        !
        if( .NOT. present( param_type ) ) then
            self%param_type = LOGE
        else
            self%param_type = trim( param_type )
        endif
        !if
        !
        nzAir = 0
        !
        allocate( self%param_grid, source = Grid3D_SG_t( grid%nx, grid%ny, nzAir, &
                    ( grid%nz - grid%nzAir ), grid%dx, grid%dy, &
                    grid%dz( grid%nzAir+1:grid%nz ) ) )
        !
        allocate( rScalar3D_SG_t :: self%cell_cond(2) )
        self%cell_cond(1) = cell_cond
        !
        if( present( param_type ) ) then
            !
            call self%setsigMap( param_type )
            !
        endif
        !
        self%is_allocated = .TRUE.
        !
    end function ModelParameterCell_SG_VTI_ctor
    !
    !>
    subroutine linComb_ModelParameterCell_SG_VTI( self, a1, a2, rhs )
        implicit none
        !
        class( ModelParameterCell_SG_VTI_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: a1, a2
        class( ModelParameter_t ), intent( in ) :: rhs
        !
        integer :: i
        complex( kind=prec ), allocatable :: v(:, :, :)
        !
        select type( rhs )
            !
            class is( ModelParameterCell_SG_VTI_t )
                !
                do i = 1, size( self%cell_cond )
                    !
                    if( self%cell_cond(i)%isCompatible( rhs%cell_cond(i) ) ) then
                        !
                        v = a1 * self%cell_cond(i)%getV() + a2 * rhs%cell_cond(i)%getV()
                        !
                        call self%cell_cond(i)%setV( v )
                        !
                    else
                        write( *, * ) "Error: linComb_ModelParameterCell_SG_VTI > Incompatible rhs cell_cond (", i, ")!"
                        stop
                    endif
                    !
                enddo
                !
            class default
                stop "Error: linComb_ModelParameterCell_SG_VTI > undefined rhs"
            !
        end select
        !
        !self%air_cond = rhs%air_cond
        !
    end subroutine linComb_ModelParameterCell_SG_VTI
    !
    !> No subroutine briefing
    !
    function dotProd_ModelParameterCell_SG_VTI( self, rhs ) result( rvalue )
        implicit none
        !
        class( ModelParameterCell_SG_VTI_t ), intent( in ) :: self
        class( ModelParameter_t ), intent( in ) :: rhs
        !
        real( kind=prec ) :: rvalue
        !
        integer :: i
        !
        rvalue = R_ZERO
        !
        select type( rhs )
            !
            class is( ModelParameterCell_SG_VTI_t )
                !
                do i = 1, size( self%cell_cond )
                    !
                    if( self%cell_cond(i)%isCompatible( rhs%cell_cond(i)) ) then
                        !
                        rvalue = rvalue + sum( self%cell_cond(i)%getV() * rhs%cell_cond(i)%getV() )
                        !
                    else
                        write( *, * ) "Error: dotProd_ModelParameterCell_SG_VTI > Incompatible rhs cell_cond (", i, ")!"
                        stop
                    endif
                    !
                enddo
                !
            class default
                stop "Error: dotProd_ModelParameterCell_SG_VTI > Unclassified rhs"
            !
        end select
        !
    end function dotProd_ModelParameterCell_SG_VTI
    !
    !> Map the entire model cells into a single edge Vector_t (eVec).
    !
    subroutine PDEmapping_ModelParameterCell_SG_VTI( self, eVec )
        implicit none
        !
        class( ModelParameterCell_SG_VTI_t ), intent( in ) :: self
        class( Vector_t ), allocatable, intent( inout ) :: eVec
        !
        type( rScalar3D_SG_t ) :: sigma_cell(2)
        integer :: i, k0, k1, k2
        !
        if( .NOT. allocated( eVec ) ) then
            allocate( eVec, source = rVector3D_SG_t( self%metric%grid, EDGE ) )
        else
            eVec = rVector3D_SG_t( self%metric%grid, EDGE )
        endif
        !
        k0 = self%metric%grid%nzAir
        k1 = k0 + 1
        k2 = self%metric%grid%Nz
        !
        do i = 1, size( self%cell_cond )
            !
            sigma_cell(i) = rScalar3D_SG_t( self%metric%grid, CELL )
            !
            sigma_cell(i)%v( :, :, 1:k0 ) = self%air_cond
            !
            sigma_cell(i)%v( :, :, k1:k2 ) = self%sigMap( real( self%cell_cond(i)%getV(), kind=prec ) )
            !
            call sigma_cell(i)%mult( self%metric%VCell )
            !
        enddo
        !> Call avgCells interface for both directions
        call eVec%avgCells( sigma_cell(1), sigma_cell(2) )
        !
        call eVec%div( self%metric%VEdge )
        !
    end subroutine PDEmapping_ModelParameterCell_SG_VTI
    !
    !> Map the perturbation between two models onto a single Vector_t (eVec).
    !
    subroutine dPDEmapping_ModelParameterCell_SG_VTI( self, dsigma, eVec )
        implicit none
        !
        class( ModelParameterCell_SG_VTI_t ), intent( in ) :: self
        class( ModelParameter_t ), intent( in ) :: dsigma
        class( Vector_t ), allocatable, intent( inout ) :: eVec
        !
        type( rScalar3D_SG_t ) :: sigma_cell(2)
        character( len=5 ), parameter :: JOB = "DERIV"
        integer :: i, k0, k1, k2
        !
        if( .NOT. allocated( eVec ) ) then
            allocate( eVec, source = rVector3D_SG_t( self%metric%grid, EDGE ) )
        else
            eVec = rVector3D_SG_t( self%metric%grid, EDGE )
        endif
        !
        call eVec%zeros
        !
        do i = 1, size( self%cell_cond )
            !
            sigma_cell(i) = rScalar3D_SG_t( self%metric%grid, CELL )
            !
            k0 = self%metric%grid%NzAir
            k1 = k0 + 1
            k2 = self%metric%grid%Nz
            !
            call sigma_cell(i)%zeros
            !
            sigma_cell(i)%v( :, :, k1:k2 ) = self%sigMap( real( self%cell_cond(i)%getV(), kind=prec ), JOB )
            !
            select type( dsigma )
                !
                class is( ModelParameterCell_SG_VTI_t )
                    !
                    sigma_cell(i)%v(:,:,k1:k2) = sigma_cell(i)%v(:,:,k1:k2) * dsigma%cell_cond(i)%getV()
                    !
                class default
                    stop "Error: dPDEmapping_ModelParameterCell_SG_VTI > Unclassified dsigma"
                !
            end select
            !
            call sigma_cell(i)%mult( self%metric%Vcell )
            !
        enddo
        !
        !> Call avgCells interface for both directions
        call eVec%avgCells( sigma_cell(1), sigma_cell(2) )
        !
        call eVec%div( self%metric%Vedge )
        !
    end subroutine dPDEmapping_ModelParameterCell_SG_VTI
    !
    !> Transpose the perturbation represented in a Vector_t (eVec), to a new dsigma model.
    !
    subroutine dPDEmapping_T_ModelParameterCell_SG_VTI( self, eVec, dsigma )
        implicit none
        !
        class( ModelParameterCell_SG_VTI_t ), intent( in ) :: self
        class( Vector_t ), intent( in ) :: eVec
        class( ModelParameter_t ), allocatable, intent( out ) :: dsigma
        !
        class( Vector_t ), allocatable :: evec_interior
        type( rScalar3D_SG_t ) :: sigma_h_cell, sigma_v_cell
        complex( kind=prec ), allocatable :: sigma_v(:, :, :)
        character( len=5 ), parameter :: JOB = "DERIV"
        integer :: i, k0, k1, k2
        !
        allocate( dsigma, source = ModelParameterCell_SG_VTI_t( self%param_grid, self%cell_cond(1), self%param_type ) )
        !
        select type( dsigma )
            !
            class is( ModelParameterCell_SG_VTI_t )
                !
                call dsigma%setCond( self%cell_cond(2), 2 )
                !
                call eVec%interior( evec_interior )
                !
                call evec_interior%div( self%metric%Vedge )
                !
                call evec_interior%mult( cmplx( 0.25_prec, 0.0, kind=prec ) )
                !
                k0 = self%metric%grid%NzAir
                k1 = k0 + 1
                k2 = self%metric%grid%Nz
                !
                sigma_h_cell = rScalar3D_SG_t( self%metric%grid, CELL )
                !
                sigma_v_cell = rScalar3D_SG_t( self%metric%grid, CELL )
                !
                call evec_interior%sumEdges( sigma_h_cell, sigma_v_cell, .TRUE. )
                !
                deallocate( evec_interior )
                !
                !> Horizontal
                sigma_v = self%sigMap( real( self%cell_cond(1)%getV(), kind=prec ), JOB )
                !
                call dsigma%cell_cond(1)%setV( sigma_v )
                !
                call sigma_h_cell%mult( self%metric%Vcell )
                !
                sigma_v = sigma_h_cell%getV()
                !
                sigma_v = sigma_v(:,:,k1:k2) * dsigma%cell_cond(1)%getV()
                !
                call dsigma%cell_cond(1)%setV( sigma_v )
                !
                !> Vertical
                sigma_v = self%sigMap( real( self%cell_cond(2)%getV(), kind=prec ), JOB )
                !
                call dsigma%cell_cond(2)%setV( sigma_v )
                !
                call sigma_v_cell%mult( self%metric%Vcell )
                !
                sigma_v = sigma_v_cell%getV()
                !
                sigma_v = sigma_v(:,:,k1:k2) * dsigma%cell_cond(2)%getV()
                !
                call dsigma%cell_cond(2)%setV( sigma_v )
                !
            class default
                stop "Error: dPDEmapping_T_ModelParameterCell_SG_VTI > Unclassified dsigma"
            !
        end select
        !
    end subroutine dPDEmapping_T_ModelParameterCell_SG_VTI
    !
    !> opens cfile on unit ioModelParam, writes out object of
    !> type modelParam in Weerachai Siripunvaraporn"s format,
    !> closes file.
    !
    subroutine write_ModelParameterCell_SG_VTI( self, file_name, comment )
        implicit none
        !
        class( ModelParameterCell_SG_VTI_t ), intent( in ) :: self
        character(*), intent( in ) :: file_name
        character(*), intent( in ), optional :: comment
        !
        type( rScalar3D_SG_t ) :: rho_v, rho_h, ccond_v
        integer :: Nx, Ny, NzEarth, i, j, k, ios
        !
        ! Verbose
        !write( *, * ) "     > Write Model to file: [", file_name, "]"
        !
        open( ioModelParam, file = file_name, action = "write", form = "formatted", iostat = ios )
        !
        if( ios == 0 ) then
            !
            if( present( comment ) ) then
                write( ioModelParam, * ) "# ", trim( comment )
            else
                write( ioModelParam, * ) "# 3D MT model written by ModEM-OO in WS format"
            endif
            !
            !> Write grid geometry definitions
            Nx = self%metric%grid%nx
            Ny = self%metric%grid%ny
            NzEarth = self%metric%grid%nz - self%metric%grid%nzAir
            !
            write( ioModelParam, "(4i5)", advance = "no" ) Nx, Ny, NzEarth, 0
            !
            write( ioModelParam, "(a10)", advance = "no" ) trim( self%param_type )
            write( ioModelParam, * ) " VTI"
            !
            !> Write self%metric%grid spacings
            do j = 1, self%metric%grid%nx
                write( ioModelParam, "(f12.3)", advance = "no" ) self%metric%grid%dx(j)
            enddo
            !
            write( ioModelParam, * )
            !
            do j = 1, self%metric%grid%ny
                write( ioModelParam, "(f12.3)", advance = "no" ) self%metric%grid%dy(j)
            enddo
            !
            write( ioModelParam, * )
            !
            do j = self%metric%grid%nzAir + 1, self%metric%grid%nz
                write( ioModelParam, "(f12.3)", advance = "no" ) self%metric%grid%dz(j)
            enddo
            !
            write( ioModelParam, * )
            !
            !> Convert (horizontal) conductivity to resistivity
            rho_h = self%cell_cond(1)
            if((index(self%param_type,"LOGE" ) > 0) .OR. (index(self%param_type,"LOG10" ) > 0)) then
                rho_h%v = -self%cell_cond(1)%getV()
            elseif(index(self%param_type,"LINEAR" ) > 0) then
                rho_h%v = ONE/self%cell_cond(1)%getV()
            endif
            !
            !> Write the (horizontal) resistivity
            !
            write( ioModelParam, * )
            !
            do k = 1, nzEarth
                do j = 1, Ny
                    do i = Nx, 1, -1
                        write( ioModelParam, "(es13.5)", iostat = ios, advance = "no" ) rho_h%v(i,j,k)
                    enddo
                    !
                    write( ioModelParam, * )
                    !
                enddo
                !
                write( ioModelParam, * )
                !
            enddo
            !
            !> Convert (vertical) conductivity to resistivity
            rho_v = self%cell_cond(2)
            if((index(self%param_type,"LOGE" ) > 0) .OR. (index(self%param_type,"LOG10" ) > 0)) then
                rho_v%v = -self%cell_cond(2)%getV()
            elseif(index(self%param_type,"LINEAR" ) > 0) then
                rho_v%v = ONE/self%cell_cond(2)%getV()
            endif
            !
            !> Write the (vertical) resistivity
            !
            write( ioModelParam, * )
            !
            do k = 1, nzEarth
                do j = 1, Ny
                    do i = Nx, 1, -1
                        write( ioModelParam, "(es13.5)", iostat = ios, advance = "no" ) rho_v%v(i,j,k)
                        enddo
                        !
                    write( ioModelParam, * )
                    !
                enddo
                !
                write( ioModelParam, * )
                !
            enddo
            !
            !> Note that our standard subroutine doesn"t work with Weerachai"s
            !> real value format. It is still better than either Mackie"s or WS"s...
            !> call write_rscalar(ioModelParam,rho)
            !> Also write the self%metric%grid origin (in metres!) and rotation (in degrees)...
            !
            write( ioModelParam, "(3f16.3)", iostat = ios) self%metric%grid%ox, self%metric%grid%oy, self%metric%grid%oz
            write( ioModelParam, "(f9.3)", iostat = ios)  self%metric%grid%rotdeg
            !
            close( ioModelParam )
            !
        else
            !
            write( *, * ) "Error opening file in write_ModelParameterCell_SG_VTI [", file_name, "]!"
            stop
            !
        endif
        !
    end subroutine write_ModelParameterCell_SG_VTI

end Module ModelParameterCell_SG_VTI
