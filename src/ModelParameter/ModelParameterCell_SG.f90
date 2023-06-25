!
!> Basic model parameter class --
!> cell conductivities defined on numerical grid
!> in Earth cells only -- as coded this is specific
!> to usual CSG.
!
module ModelParameterCell_SG
    !
    use Constants
    use FileUnits
    use rScalar3D_SG
    use rVector3D_SG 
    use ModelParameter
    use Grid3D_SG
    use ModelParameter1D
    use ModelParameter2D
    !
    type, extends( ModelParameter_t ) :: ModelParameterCell_SG_t
        !
        class( Grid_t ), allocatable :: param_grid
        !
        class( Scalar_t ), allocatable :: cell_cond
        !
        contains
            !
            final :: ModelParameterCell_SG_dtor
            !
            procedure, public :: getCond => getCondModel_ParameterCell_SG
            procedure, public :: setCond => setCondModel_ParameterCell_SG
            !
            procedure, public :: zeros => zeros_ModelParameterCell_SG
            !
            procedure, public :: copyFrom => copyFrom_ModelParameterCell_SG
            !
            procedure, public :: countModel => countModel_ModelParameterCell_SG
            !
            procedure, public :: dotProd => dotProd_ModelParameterCell_SG
            !
            procedure, public :: linComb => linComb_ModelParameterCell_SG
            !
            procedure, public :: PDEmapping => PDEmapping_ModelParameterCell_SG
            procedure, public :: dPDEmapping => dPDEmapping_ModelParameterCell_SG
            procedure, public :: dPDEmapping_T => dPDEmapping_T_ModelParameterCell_SG
            !
            procedure, public :: slice1D => slice1D_ModelParameterCell_SG
            procedure, public :: slice2D => slice2D_ModelParameterCell_SG
            !
            procedure, public :: avgModel1D => avgModel1D_ModelParameterCell_SG
            !
            procedure, public :: setType => setType_ModelParameterCell_SG
            !
            procedure, public :: write => write_ModelParameterCell_SG
            !
            procedure, public :: print => print_ModelParameterCell_SG
            !
    end type ModelParameterCell_SG_t
    !
    interface ModelParameterCell_SG_t
         module procedure ModelParameterCell_SG_ctor
    end interface ModelParameterCell_SG_t
    !
contains
    !
    !> No subroutine briefing
    !
    function ModelParameterCell_SG_ctor( grid, cell_cond, param_type ) result( self )
        implicit none
        !
        class( Grid_t ), intent( in ) :: grid
        class( Scalar_t ), intent( in ) :: cell_cond
        character(:), allocatable, optional, intent( in ) :: param_type
        !
        type( ModelParameterCell_SG_t ) :: self
        !
        integer :: nzAir
        !
        !write( *, * ) "Constructor ModelParameterCell_SG_t"
        !
        call self%init
        !
        if( .NOT. present( param_type ) ) then
            self%param_type = LOGE
        else
            self%param_type = trim( param_type )
        endif
        !
        nzAir = 0
        !
        allocate( self%param_grid, source = Grid3D_SG_t( grid%nx, grid%ny, nzAir, &
                    ( grid%nz - grid%nzAir ), grid%dx, grid%dy, &
                    grid%dz( grid%nzAir+1:grid%nz ) ) )
        !
        allocate( self%cell_cond, source = cell_cond )
        !
        self%cell_cond%store_state = compound
        !
        if( present( param_type ) ) then
            !
            call self%setsigMap( param_type )
        !
        endif
        !
        self%is_allocated = .TRUE.
        !
    end function ModelParameterCell_SG_ctor
    !
    !> No subroutine briefing
    subroutine ModelParameterCell_SG_dtor( self )
        implicit none
        !
        type( ModelParameterCell_SG_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor ModelParameterCell_SG"
        !
        if( allocated( self%param_grid ) ) deallocate( self%param_grid )
        !
        deallocate( self%cell_cond )
        !
    end subroutine ModelParameterCell_SG_dtor
    !
    !> No subroutine briefing
    !
    function slice1D_ModelParameterCell_SG( self, ix, iy ) result( model_param_1D )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( in ) :: self
        integer, intent( in ) :: ix, iy
        !
        type( ModelParameter1D_t ) ::  model_param_1D 
        !
        complex( kind=prec ), allocatable :: v(:, :, :)
        real( kind=prec ), allocatable, dimension(:) :: cond_slice
        !
        model_param_1D = ModelParameter1D_t( self%metric%grid%Slice1D() )
        !
        allocate( cond_slice( model_param_1D%grid%nz ) )
        !
        v = self%cell_cond%getV()
        cond_slice = self%sigMap( real( v( ix, iy, : ), kind=prec ) )
        !
        call model_param_1D%setConductivity( cond_slice, self%air_cond, self%param_type, self%mKey )
        !
        deallocate( cond_slice )
        !
    end function slice1D_ModelParameterCell_SG
    !
    !> No subroutine briefing
    !
    function avgModel1D_ModelParameterCell_SG( self ) result( model_param_1D )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( in ) :: self
        !
        type( ModelParameter1D_t ) ::  model_param_1D 
        !
        complex( kind=prec ), allocatable :: v(:, :, :)
        real( kind=prec ), allocatable, dimension(:) :: cond_slice
        real( kind=prec ) :: wt, temp_sigma_value
        integer :: i, j, k
        !
        v = self%cell_cond%getV()
        !
        model_param_1D = ModelParameter1D_t( self%metric%grid%Slice1D() )
        !
        allocate( cond_slice( self%metric%grid%nzEarth ) )
        !
        do k = 1, self%metric%grid%nzEarth
            !
            wt = R_ZERO
            temp_sigma_value = R_ZERO
            do i = 1, self%metric%grid%Nx
                do j = 1, self%metric%grid%Ny
                    wt = wt + self%metric%grid%dx(i) * self%metric%grid%dy(j)
                    temp_sigma_value = temp_sigma_value + v( i, j, k ) * &
                    self%metric%grid%dx(i) * self%metric%grid%dy(j)
                enddo
            enddo
            !
            cond_slice(k) = self%sigMap( temp_sigma_value / wt )
            !
        enddo
        !
        call model_param_1D%setConductivity( cond_slice, self%air_cond, self%param_type, self%mKey )
        !
        deallocate( cond_slice )
        !
    end function avgModel1D_ModelParameterCell_SG
    !
    !> No subroutine briefing
    !
    function slice2D_ModelParameterCell_SG( self, axis, j ) result( m2D )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( in ) :: self
        integer, intent( in ) :: axis, j
        !
        type( ModelParameter2D_t ) :: m2D 
        !
        character(:), allocatable :: param_type
        complex( kind=prec ), allocatable :: v(:, :, :)
        real( kind=prec ), allocatable, dimension(:,:) :: cond_slice
        !
        param_type = LINEAR
        !
        m2D = ModelParameter2D_t( self%metric%grid%Slice2D() )
        !
        allocate( cond_slice( self%metric%grid%ny, self%metric%grid%nzEarth ) )
        !
        v = self%cell_cond%getV()
        !
        if( axis == 1 ) then
            cond_slice = self%sigMap( real( v(j,:,:), kind=prec ) )
        elseif( axis == 2 ) then
            cond_slice = self%sigMap( real( v(:,j,:), kind=prec ) )
        elseif( axis == 3 ) then
            cond_slice = self%sigMap( real( v(:,:,j), kind=prec ) )
        else
            stop "Error: slice2D_ModelParameterCell_SG > wrong axis"
        endif
        !
        call m2D%setConductivity( cond_slice, self%air_cond, param_type, self%mKey )
        !
        deallocate( cond_slice )
        !
    end function slice2D_ModelParameterCell_SG
    !
    !> No interface subroutine briefing
    !
    subroutine getCondModel_ParameterCell_SG( self, ccond )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( in ) :: self
        class( Scalar_t ), allocatable, intent( inout ) :: ccond
        !
        if( allocated( ccond ) ) deallocate( ccond )
        allocate( ccond, source = self%cell_cond )
        !
    end subroutine getCondModel_ParameterCell_SG
    !
    !> No interface subroutine briefing
    !
    subroutine setCondModel_ParameterCell_SG( self, ccond, i_cond )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( inout ) :: self
        class( Scalar_t ), allocatable, intent( in ) :: ccond
        integer, intent( in ) :: i_cond
        !
        if( i_cond /= 1 ) then
            stop "Error: setCondModel_ParameterCell_SG > This class does not support anisotropy"
        endif
        !
        if( allocated( self%cell_cond ) ) deallocate( self%cell_cond )
        allocate( self%cell_cond, source = ccond )
        !
    end subroutine setCondModel_ParameterCell_SG
!
    !> No subroutine briefing
    subroutine zeros_ModelParameterCell_SG( self )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( inout ) :: self
        !
        call self%cell_cond%zeros
        !
    end subroutine zeros_ModelParameterCell_SG
    !
    !> No subroutine briefing
    subroutine copyFrom_ModelParameterCell_SG( self, rhs )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( inout ) :: self
        class( ModelParameter_t ), intent( in ) :: rhs
        !
        select type( rhs )
            !
            class is( ModelParameterCell_SG_t )
                !
                self%metric => rhs%metric
                !
                self%mKey = rhs%mKey
                !
                self%air_cond = rhs%air_cond
                !
                self%param_type = rhs%param_type
                !
                self%is_allocated = rhs%is_allocated
                !
                self%param_grid = rhs%param_grid
                !
                self%cell_cond = rhs%cell_cond
                !
                self%sigMap_ptr => rhs%sigMap_ptr
                !
            class default
               stop "Error: copyFrom_ModelParameterCell_SG > Unclassified rhs."
            !
        end select
        !
    end subroutine copyFrom_ModelParameterCell_SG
    !
    !> ????
    function countModel_ModelParameterCell_SG( self ) result( counter )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( in ) :: self
        !
        integer :: counter, nx, ny, nz, nzAir, nz_earth
        !
        if( .NOT. self%cell_cond%is_allocated ) then
            stop "Error: countModel_ModelParameterCell_SG > cell_cond not allocated!"
        endif
        !
        !
        !> Grid dimensions
        call self%cell_cond%grid%getDimensions( nx, ny, nz, nzAir )
        nz_earth = nz - nzAir
        !
        counter = self%cell_cond%Nx * self%cell_cond%Ny * nz_earth
        !
    end function countModel_ModelParameterCell_SG
    !
    !>
    subroutine linComb_ModelParameterCell_SG( self, a1, a2, rhs )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: a1, a2
        class( ModelParameter_t ), intent( in ) :: rhs
        !
        complex( kind=prec ), allocatable :: v(:, :, :)
        !
        select type( rhs )
            !
            class is( ModelParameterCell_SG_t )
                !
                if( self%cell_cond%isCompatible( rhs%cell_cond ) ) then
                    !
                    v = a1 * self%cell_cond%getV() + a2 * rhs%cell_cond%getV()
                    call self%cell_cond%setV( v )
                    !
                else
                    stop "Error: linComb_ModelParameterCell_SG > Incompatible rhs"
                endif
                !
            class default
                stop "Error: linComb_ModelParameterCell_SG > undefined rhs"
            !
        end select
        !
        !self%air_cond = rhs%air_cond
        !
    end subroutine linComb_ModelParameterCell_SG
    !
    !> No subroutine briefing
    !
    function dotProd_ModelParameterCell_SG( self, rhs ) result( rvalue )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( in ) :: self
        class( ModelParameter_t ), intent( in ) :: rhs
        real( kind=prec ) :: rvalue
        !
        select type( rhs )
            !
            class is( ModelParameterCell_SG_t )
                !
                if( self%cell_cond%isCompatible( rhs%cell_cond ) ) then
                    !
                    rvalue = sum( self%cell_cond%getV() * rhs%cell_cond%getV() )
                    !
                else
                    stop "Error: dotProd_ModelParameterCell_SG > Incompatible rhs"
                endif
                !
            class default
                stop "Error: dotProd_ModelParameterCell_SG > undefined rhs"
            !
        end select
        !
    end function dotProd_ModelParameterCell_SG
    !
    !> Map the entire model cells into a single edge Vector_t (eVec).
    !
    subroutine PDEmapping_ModelParameterCell_SG( self, eVec )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( in ) :: self
        class( Vector_t ), allocatable, intent( inout ) :: eVec
        !
        type( rScalar3D_SG_t ) :: sigma_cell
        integer :: k0, k1, k2
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
        sigma_cell = rScalar3D_SG_t( self%metric%grid, CELL )
        !
        sigma_cell%v( :, :, 1:k0 ) = self%air_cond
        !
        sigma_cell%v( :, :, k1:k2 ) = self%sigMap( real( self%cell_cond%getV(), kind=prec ) )
        !
        call sigma_cell%mult( self%metric%VCell )
        !
        call eVec%avgCells( sigma_cell )
        !
        call eVec%div( self%metric%VEdge )
        !
    end subroutine PDEmapping_ModelParameterCell_SG
    !
    !> Map the perturbation between two models onto a single Vector_t (eVec).
    !
    subroutine dPDEmapping_ModelParameterCell_SG( self, dsigma, eVec )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( in ) :: self
        class( ModelParameter_t ), intent( in ) :: dsigma
        class( Vector_t ), allocatable, intent( inout ) :: eVec
        !
        type( rScalar3D_SG_t ) :: sigma_cell
        character( len=5 ), parameter :: JOB = "DERIV"
        integer :: k0, k1, k2
        !
        if( .NOT. allocated( eVec ) ) then
            allocate( eVec, source = rVector3D_SG_t( self%metric%grid, EDGE ) )
        else
            eVec = rVector3D_SG_t( self%metric%grid, EDGE )
        endif
        !
        call eVec%zeros
        !
        sigma_cell = rScalar3D_SG_t( self%metric%grid, CELL )
        !
        k0 = self%metric%grid%NzAir
        k1 = k0 + 1
        k2 = self%metric%grid%Nz
        !
        !> Ensure values in air are zero.
        call sigma_cell%zeros
        !
        sigma_cell%v( :, :, k1:k2 ) = self%sigMap( real( self%cell_cond%getV(), kind=prec ), JOB )
        !
        !> Required to access the cell_cond attribute of ModelParameterCell_SG
        select type( dsigma )
            !
            class is( ModelParameterCell_SG_t )
                !
                sigma_cell%v(:,:,k1:k2) = sigma_cell%v(:,:,k1:k2) * dsigma%cell_cond%getV()
                !
            class default
                stop "Error: dPDEmapping_ModelParameterCell_SG > Unclassified dsigma"
            !
        end select
        !
        call sigma_cell%mult( self%metric%Vcell )
        !
        call eVec%avgCells( sigma_cell )
        !
        call eVec%div( self%metric%Vedge )
        !
    end subroutine dPDEmapping_ModelParameterCell_SG
    !
    !> Transpose the perturbation represented in a Vector_t (eVec), to a new dsigma model.
    !
    subroutine dPDEmapping_T_ModelParameterCell_SG( self, eVec, dsigma )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( in ) :: self
        class( Vector_t ), intent( in ) :: eVec
        class( ModelParameter_t ), allocatable, intent( out ) :: dsigma
        !
        class( Vector_t ), allocatable :: temp_interior
        class( Scalar_t ), allocatable :: sigma_cell
        complex( kind=prec ), allocatable :: v(:, :, :), s_v(:, :, :)
        character( len=5 ), parameter :: JOB = "DERIV"
        integer :: k0, k1, k2
        !
        allocate( dsigma, source = ModelParameterCell_SG_t( self%param_grid, self%cell_cond, self%param_type ) )
        !
        select type( dsigma )
            !
            class is( ModelParameterCell_SG_t )
                !
                call eVec%interior( temp_interior )
                !
                call temp_interior%div( self%metric%Vedge )
                !
                call temp_interior%mult( cmplx( 0.25_prec, 0.0, kind=prec ) )
                !
                call temp_interior%sumEdges( sigma_cell, .TRUE. )
                !
                deallocate( temp_interior )
                !
                call sigma_cell%mult( self%metric%Vcell )
                !
                v = self%sigMap( real( self%cell_cond%getV(), kind=prec ), JOB )
                call dsigma%cell_cond%setV( v )
                !
                k0 = self%metric%grid%NzAir
                k1 = k0 + 1
                k2 = self%metric%grid%Nz
                !
                s_v = sigma_cell%getV()
                !
                v = dsigma%cell_cond%getV() * s_v(:,:,k1:k2)
                !
                call dsigma%cell_cond%setV( v )
                !
            class default
                stop "Error: dPDEmapping_T_ModelParameterCell_SG > Incompatible input [eVec]."
        end select
        !
    end subroutine dPDEmapping_T_ModelParameterCell_SG
    !
    !> No subroutine briefing
    !
    subroutine setType_ModelParameterCell_SG( self, param_type )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( inout ) :: self
        character(:), allocatable, intent( in ) :: param_type
        !
        complex( kind=prec ), allocatable :: v(:, :, :)
        !
        if( .NOT. self%is_allocated ) then
            stop "Error: setType_ModelParameterCell_SG > Not allocated."
        endif
        !
        v = self%cell_cond%getV()
        if( trim( param_type ) .EQ. trim( self%param_type ) ) then
            ! Nothing to be done
        elseif( self%param_type == "" ) then
            self%param_type = trim( param_type )
        elseif( self%param_type == LINEAR ) then
            !
            if( param_type == LOGE ) then
                v = log( v )
                call self%cell_cond%setV( v )
            elseif( param_type == LOG_10) then
                v = log10( real( v, kind=prec ) )
                call self%cell_cond%setV( v )
            endif
            !
        elseif( param_type == LINEAR ) then
            !
            if( self%param_type == LOGE ) then
                v = exp( v )
                call self%cell_cond%setV( v )
            elseif( self%param_type == LOG_10 ) then
                v = exp( v * log(10.) )
                call self%cell_cond%setV( v )
            endif
            !
        elseif( ( self%param_type == LOGE ) .AND. ( param_type == LOG_10 ) ) then
            v = v / log(10.)
            call self%cell_cond%setV( v )
        elseif( ( self%param_type == LOG_10 ) .AND. ( param_type == LOGE ) ) then
            v = v * log(10.)
            call self%cell_cond%setV( v )
        else
            stop "Error: setType_ModelParameterCell_SG > Unknown param_type."
        endif
        !
        self%param_type = param_type 
        !
    end subroutine setType_ModelParameterCell_SG
    !
    !> No subroutine briefing
    !
    subroutine print_ModelParameterCell_SG( self )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( in ) :: self
        !
        write( *, * ) "ModelParameterCell_SG_t:", self%mKey, self%air_cond, self%param_type, &
        self%is_allocated, self%param_grid%nx, self%param_grid%ny, self%param_grid%nz, self%param_grid%nzAir
        !
        !call self%cell_cond%print
        !
    end subroutine print_ModelParameterCell_SG
    !
    !> opens cfile on unit ioModelParam, writes out object of
    !> type modelParam in Weerachai Siripunvaraporn"s format,
    !> closes file.
    !
    subroutine write_ModelParameterCell_SG( self, file_name, comment )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( in ) :: self
        character(*), intent( in ) :: file_name
        character(*), intent( in ), optional :: comment
        !
        type( rScalar3D_SG_t ) :: rho_v, rho_h, ccond_v
        integer :: Nx, Ny, NzEarth, i, j, k, ios
        !
        ! Verbose
        !write( *, * ) "     > Write Model to file: [", file_name, "]"
        !
        !> Convert modelParam to natural log or log10 for output
        !paramType = userParamType
        !
        !if( self%%is_vti ) then
        !    call getValue_modelParam(m, paramType, self%cell_cond, v_v=ccond_v)
        !else
        !    call getValue_modelParam(m, paramType, self%cell_cond)
        !endif
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
            write( ioModelParam, "(a10)", advance = "yes" ) trim( self%param_type )
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
            rho_h = self%cell_cond
            if((index(self%param_type,"LOGE" ) > 0) .OR. (index(self%param_type,"LOG10" ) > 0)) then
                rho_h%v = -self%cell_cond%getV()
            elseif(index(self%param_type,"LINEAR" ) > 0) then
                rho_h%v = ONE/self%cell_cond%getV()
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
            write( *, * ) "Error opening file in write_ModelParameterCell_SG [", file_name, "]!"
            stop
            !
        endif
        !
    end subroutine write_ModelParameterCell_SG

end Module ModelParameterCell_SG
