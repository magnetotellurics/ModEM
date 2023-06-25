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
    type, extends( ModelParameter_t ) :: ModelParameterCell_SG_VTI_t
        !
        class( Grid_t ), allocatable :: param_grid
        !
        class( Scalar_t ), allocatable :: cell_cond_h, cell_cond_v
        !
        contains
            !
            final :: ModelParameterCell_SG_VTI_dtor
            !
            procedure, public :: getCond => getCond_ModelParameterCell_SG_VTI
            procedure, public :: setCond => setCond_Model_ParameterCell_SG_VTI
            !
            procedure, public :: zeros => zeros_ModelParameterCell_SG_VTI
            !
            procedure, public :: copyFrom => copyFrom_ModelParameterCell_SG_VTI
            !
            procedure, public :: countModel => countModel_ModelParameterCell_SG_VTI
            !
            procedure, public :: dotProd => dotProd_ModelParameterCell_SG_VTI
            !
            procedure, public :: linComb => linComb_ModelParameterCell_SG_VTI
            !
            procedure, public :: PDEmapping => PDEmapping_ModelParameterCell_SG_VTI
            procedure, public :: dPDEmapping => dPDEmapping_ModelParameterCell_SG_VTI
            procedure, public :: dPDEmapping_T => dPDEmapping_T_ModelParameterCell_SG_VTI
            !
            procedure, public :: slice1D => slice1D_ModelParameterCell_SG_VTI
            procedure, public :: slice2D => slice2D_ModelParameterCell_SG_VTI
            !
            procedure, public :: avgModel1D => avgModel1D_ModelParameterCell_SG_VTI
            !
            procedure, public :: setType => setType_ModelParameterCell_SG_VTI
            !
            procedure, public :: write => write_ModelParameterCell_SG_VTI
            !
            procedure, public :: print => print_ModelParameterCell_SG_VTI
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
    function ModelParameterCell_SG_VTI_ctor( grid, cell_cond, param_type, anisotropic_level ) result( self )
        implicit none
        !
        class( Grid_t ), intent( in ) :: grid
        class( Scalar_t ), intent( in ) :: cell_cond
        character(:), allocatable, optional, intent( in ) :: param_type
        integer, optional, intent( in ) :: anisotropic_level
        !
        type( ModelParameterCell_SG_VTI_t ) :: self
        !
        integer :: n_cond, nzAir
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
        !
        if( present( anisotropic_level ) ) then
            n_cond = anisotropic_level
        else
            n_cond = 1
        endif
        !
        nzAir = 0
        !
        allocate( self%param_grid, source = Grid3D_SG_t( grid%nx, grid%ny, nzAir, &
                    ( grid%nz - grid%nzAir ), grid%dx, grid%dy, &
                    grid%dz( grid%nzAir+1:grid%nz ) ) )
        !
        allocate( self%cell_cond_h, source = cell_cond )
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
    !> No subroutine briefing
    subroutine ModelParameterCell_SG_VTI_dtor( self )
        implicit none
        !
        type( ModelParameterCell_SG_VTI_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor ModelParameterCell_SG_VTI"
        !
        deallocate( self%param_grid )
        !
        deallocate( self%cell_cond_h )
        deallocate( self%cell_cond_v )
        !
    end subroutine ModelParameterCell_SG_VTI_dtor
    !
    !> No subroutine briefing
    !
    function slice1D_ModelParameterCell_SG_VTI( self, ix, iy ) result( model_param_1D )
        implicit none
        !
        class( ModelParameterCell_SG_VTI_t ), intent( in ) :: self
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
        v = self%cell_cond_h%getV()
        cond_slice = self%sigMap( real( v( ix, iy, : ), kind=prec ) )
        !
        call model_param_1D%setConductivity( cond_slice, self%air_cond, self%param_type, self%mKey )
        !
        deallocate( cond_slice )
        !
    end function slice1D_ModelParameterCell_SG_VTI
    !
    !> No subroutine briefing
    !
    function avgModel1D_ModelParameterCell_SG_VTI( self ) result( model_param_1D )
        implicit none
        !
        class( ModelParameterCell_SG_VTI_t ), intent( in ) :: self
        !
        type( ModelParameter1D_t ) ::  model_param_1D 
        !
        complex( kind=prec ), allocatable :: v(:, :, :)
        real( kind=prec ), allocatable, dimension(:) :: cond_slice
        complex( kind=prec ) :: wt, temp_sigma_value
        integer :: i, j, k
        !
        v = self%cell_cond_h%getV()
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
            cond_slice(k) = self%sigMap( real( temp_sigma_value / wt, kind=prec ) )
            !
        enddo
        !
        call model_param_1D%setConductivity( cond_slice, self%air_cond, self%param_type, self%mKey )
        !
        deallocate( cond_slice )
        !
    end function avgModel1D_ModelParameterCell_SG_VTI
    !
    !> No subroutine briefing
    !
    function slice2D_ModelParameterCell_SG_VTI( self, axis, j ) result( m2D )
        implicit none
        !
        class( ModelParameterCell_SG_VTI_t ), intent( in ) :: self
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
        v = self%cell_cond_h%getV()
        !real( self%cell_cond_h%getV(), kind=prec )
        if( axis == 1 ) then
            cond_slice = self%sigMap( real( v(j,:,:), kind=prec ) )
        elseif( axis == 2 ) then
            cond_slice = self%sigMap( real( v(:,j,:), kind=prec ) )
        elseif( axis == 3 ) then
            cond_slice = self%sigMap( real( v(:,:,j), kind=prec ) )
        else
            stop "Error: slice2D_ModelParameterCell_SG_VTI > wrong axis"
        endif
        !
        call m2D%setConductivity( cond_slice, self%air_cond, param_type, self%mKey )
        !
        deallocate( cond_slice )
        !
    end function slice2D_ModelParameterCell_SG_VTI
    !
    !> No interface subroutine briefing
    !
    subroutine getCond_ModelParameterCell_SG_VTI( self, ccond )
        implicit none
        !
        class( ModelParameterCell_SG_VTI_t ), intent( in ) :: self
            class( Scalar_t ), allocatable, intent( inout ) :: ccond
        !
        if( allocated( ccond ) ) deallocate( ccond )
        allocate( ccond, source = self%cell_cond_h )
        !
    end subroutine getCond_ModelParameterCell_SG_VTI
    !
    !> No interface subroutine briefing
    !
    subroutine setCond_Model_ParameterCell_SG_VTI( self, ccond, i_cond )
        implicit none
        !
        class( ModelParameterCell_SG_VTI_t ), intent( inout ) :: self
        class( Scalar_t ), allocatable, intent( in ) :: ccond
        integer, intent( in ) :: i_cond
        !
        integer :: i
        !
        if( i_cond == 1 ) then
            !
            if( allocated( self%cell_cond_h ) ) deallocate( self%cell_cond_h )
            allocate( self%cell_cond_h, source = ccond )
            !
        elseif( i_cond == 2 ) then
            !
            if( allocated( self%cell_cond_v ) ) deallocate( self%cell_cond_v )
            allocate( self%cell_cond_v, source = ccond )
            !
        else
            !
            stop "Error: setCond_Model_ParameterCell_SG_VTI > VTI not support anisotropy level greater than 2"
            !
        endif
        !
    end subroutine setCond_Model_ParameterCell_SG_VTI
    !
    !> No subroutine briefing
    subroutine zeros_ModelParameterCell_SG_VTI( self )
        implicit none
        !
        class( ModelParameterCell_SG_VTI_t ), intent( inout ) :: self
        !
        call self%cell_cond_h%zeros
        call self%cell_cond_v%zeros
        !
    end subroutine zeros_ModelParameterCell_SG_VTI
    !
    !> No subroutine briefing
    subroutine copyFrom_ModelParameterCell_SG_VTI( self, rhs )
        implicit none
        !
        class( ModelParameterCell_SG_VTI_t ), intent( inout ) :: self
        class( ModelParameter_t ), intent( in ) :: rhs
        !
        select type( rhs )
            !
            class is( ModelParameterCell_SG_VTI_t )
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
                self%cell_cond_h = rhs%cell_cond_h
                !
                self%cell_cond_v = rhs%cell_cond_v
                !
                self%sigMap_ptr => rhs%sigMap_ptr
                !
            class default
               stop "Error: copyFrom_ModelParameterCell_SG_VTI > Incompatible input."
            !
        end select
        !
    end subroutine copyFrom_ModelParameterCell_SG_VTI
    !
    !> ????
    function countModel_ModelParameterCell_SG_VTI( self ) result( counter )
        implicit none
        !
        class( ModelParameterCell_SG_VTI_t ), intent( in ) :: self
        !
        integer :: counter, nx, ny, nz, nzAir, nz_earth
        !
        counter = 0
        !
        if( .NOT. self%cell_cond_h%is_allocated ) then
            stop "Error: countModel_ModelParameterCell_SG_VTI > cell_cond_h not allocated!"
        endif
        !
        !> Horizontal
        call self%cell_cond_h%grid%getDimensions( nx, ny, nz, nzAir )
        nz_earth = nz - nzAir
        !
        counter = counter + self%cell_cond_h%Nx * self%cell_cond_h%Ny * nz_earth
        !
        !> Vertical
        call self%cell_cond_v%grid%getDimensions( nx, ny, nz, nzAir )
        nz_earth = nz - nzAir
        !
        counter = counter + self%cell_cond_v%Nx * self%cell_cond_v%Ny * nz_earth
        !
    end function countModel_ModelParameterCell_SG_VTI
    !
    !>
    subroutine linComb_ModelParameterCell_SG_VTI( self, a1, a2, rhs )
        implicit none
        !
        class( ModelParameterCell_SG_VTI_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: a1, a2
        class( ModelParameter_t ), intent( in ) :: rhs
        !
        complex( kind=prec ), allocatable :: v(:, :, :)
        !
        select type( rhs )
            !
            class is( ModelParameterCell_SG_VTI_t )
                !
                !> Horizontal
                if( self%cell_cond_h%isCompatible( rhs%cell_cond_h ) ) then
                    !
                    v = a1 * self%cell_cond_h%getV() + a2 * rhs%cell_cond_h%getV()
                    !
                    call self%cell_cond_h%setV( v )
                    !
                else
                    stop "Error: linComb_ModelParameterCell_SG_VTI > Incompatible rhs cell_cond_h"
                endif
                !
                !> Vertical
                if( self%cell_cond_h%isCompatible( rhs%cell_cond_h ) ) then
                    !
                    v = a1 * self%cell_cond_v%getV() + a2 * rhs%cell_cond_v%getV()
                    !
                    call self%cell_cond_v%setV( v )
                    !
                else
                    stop "Error: linComb_ModelParameterCell_SG_VTI > Incompatible rhs cell_cond_v"
                endif
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
        real( kind=prec ) :: rvalue
        integer :: i
        !
        rvalue = R_ZERO
        !
        select type( rhs )
            !
            class is( ModelParameterCell_SG_VTI_t )
                !
                !> Horizontal
                if( self%cell_cond_h%isCompatible( rhs%cell_cond_h) ) then
                    !
                    rvalue = rvalue + sum( self%cell_cond_h%getV() * rhs%cell_cond_h%getV() )
                    !
                else
                    stop "Error: dotProd_ModelParameterCell_SG_VTI > Incompatible rhs"
                endif
                !
                !> Vertical
                if( self%cell_cond_v%isCompatible( rhs%cell_cond_v ) ) then
                    !
                    rvalue = rvalue + sum( self%cell_cond_v%getV() * rhs%cell_cond_v%getV() )
                    !
                else
                    stop "Error: dotProd_ModelParameterCell_SG_VTI > Incompatible rhs"
                endif
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
        type( rScalar3D_SG_t ) :: sigma_cell_h, sigma_cell_v
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
        !> Horizontal
        !
        sigma_cell_h = rScalar3D_SG_t( self%metric%grid, CELL )
        !
        sigma_cell_h%v( :, :, 1:k0 ) = self%air_cond
        !
        sigma_cell_h%v( :, :, k1:k2 ) = self%sigMap( real( self%cell_cond_h%getV(), kind=prec ) )
        !
        call sigma_cell_h%mult( self%metric%VCell )
        !
        !> Vertical
        !
        sigma_cell_v = rScalar3D_SG_t( self%metric%grid, CELL )
        !
        sigma_cell_v%v( :, :, 1:k0 ) = self%air_cond
        !
        sigma_cell_v%v( :, :, k1:k2 ) = self%sigMap( real( self%cell_cond_v%getV(), kind=prec ) )
        !
        call sigma_cell_v%mult( self%metric%VCell )
        !
        !> Call avgCells interface for both directions
        call eVec%avgCells( sigma_cell_h, sigma_cell_v )
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
        type( rScalar3D_SG_t ) :: sigma_cell_h, sigma_cell_v
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
        !> Horizontal
        sigma_cell_h = rScalar3D_SG_t( self%metric%grid, CELL )
        !
        k0 = self%metric%grid%NzAir
        k1 = k0 + 1
        k2 = self%metric%grid%Nz
        !
        call sigma_cell_h%zeros
        !
        sigma_cell_h%v( :, :, k1:k2 ) = self%sigMap( real( self%cell_cond_h%getV(), kind=prec ), JOB )
        !
        select type( dsigma )
            !
            class is( ModelParameterCell_SG_VTI_t )
                !
                sigma_cell_h%v(:,:,k1:k2) = sigma_cell_h%v(:,:,k1:k2) * dsigma%cell_cond_h%getV()
                !
            class default
                stop "Error: dPDEmapping_ModelParameterCell_SG_VTI > Unclassified dsigma"
            !
        end select
        !
        call sigma_cell_h%mult( self%metric%Vcell )
        !
        !> Vertical
        sigma_cell_v = rScalar3D_SG_t( self%metric%grid, CELL )
        !
        k0 = self%metric%grid%NzAir
        k1 = k0 + 1
        k2 = self%metric%grid%Nz
        !
        call sigma_cell_v%zeros
        !
        sigma_cell_v%v( :, :, k1:k2 ) = self%sigMap( real( self%cell_cond_v%getV(), kind=prec ), JOB )
        !
        select type( dsigma )
            !
            class is( ModelParameterCell_SG_VTI_t )
                !
                sigma_cell_v%v(:,:,k1:k2) = sigma_cell_v%v(:,:,k1:k2) * dsigma%cell_cond_v%getV()
                !
            class default
                stop "Error: dPDEmapping_ModelParameterCell_SG_VTI > Unclassified dsigma"
            !
        end select
        !
        call sigma_cell_v%mult( self%metric%Vcell )
        !
        !> Call avgCells interface for both directions
        call eVec%avgCells( sigma_cell_h, sigma_cell_v )
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
        class( Vector_t ), allocatable :: temp_interior
        class( Scalar_t ), allocatable :: sigma_cell_h, sigma_cell_v
        complex( kind=prec ), allocatable :: v(:, :, :), s_v(:, :, :)
        character( len=5 ), parameter :: JOB = "DERIV"
        integer :: k0, k1, k2
        !
        allocate( dsigma, source = ModelParameterCell_SG_VTI_t( self%param_grid, self%cell_cond_h, self%param_type ) )
        !
        call dsigma%setCond( self%cell_cond_v, 2 )
        !
        select type( dsigma )
            !
            class is( ModelParameterCell_SG_VTI_t )
                !
                call eVec%interior( temp_interior )
                !
                call temp_interior%div( self%metric%Vedge )
                !
                call temp_interior%mult( cmplx( 0.25_prec, 0.0, kind=prec ) )
                !
                !> Horizontal
                call temp_interior%sumEdges( sigma_cell_h, .TRUE. )
                !
                call sigma_cell_h%mult( self%metric%Vcell )
                !
                v = self%sigMap( real( self%cell_cond_h%getV(), kind=prec ), JOB )
                call dsigma%cell_cond_h%setV( v )
                !
                k0 = self%metric%grid%NzAir
                k1 = k0 + 1
                k2 = self%metric%grid%Nz
                !
                s_v = sigma_cell_h%getV()
                !
                v = dsigma%cell_cond_h%getV() * s_v(:,:,k1:k2)
                !
                call dsigma%cell_cond_h%setV( v )
                !
                !> Vertical
                call temp_interior%sumEdges( sigma_cell_v, .TRUE. )
                !
                deallocate( temp_interior )
                !
                call sigma_cell_v%mult( self%metric%Vcell )
                !
                v = self%sigMap( real( self%cell_cond_v%getV(), kind=prec ), JOB )
                call dsigma%cell_cond_v%setV( v )
                !
                k0 = self%metric%grid%NzAir
                k1 = k0 + 1
                k2 = self%metric%grid%Nz
                !
                s_v = sigma_cell_v%getV()
                !
                v = dsigma%cell_cond_v%getV() * s_v(:,:,k1:k2)
                !
                call dsigma%cell_cond_v%setV( v )
                !
            class default
                stop "Error: dPDEmapping_T_ModelParameterCell_SG_VTI > Unclassified dsigma"
            !
        end select
        !
    end subroutine dPDEmapping_T_ModelParameterCell_SG_VTI
    !
    !> No subroutine briefing
    !
    subroutine setType_ModelParameterCell_SG_VTI( self, param_type )
        implicit none
        !
        class( ModelParameterCell_SG_VTI_t ), intent( inout ) :: self
        character(:), allocatable, intent( in ) :: param_type
        !
        complex( kind=prec ), allocatable :: v_h(:, :, :), v_v(:, :, :)
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            stop "Error: setType_ModelParameterCell_SG_VTI > Self not allocated."
        endif
        !
        !> Horizontal & Vertical
        v_h = self%cell_cond_h%getV()
        v_v = self%cell_cond_v%getV()
        !
        if( trim( param_type ) .EQ. trim( self%param_type ) ) then
            ! Nothing to be done
        elseif( self%param_type == "" ) then
            self%param_type = trim( param_type )
        elseif( self%param_type == LINEAR ) then
            !
            if( param_type == LOGE ) then
                !
                v_h = log( v_h )
                call self%cell_cond_h%setV( v_h )
                !
                v_v = log( v_v )
                call self%cell_cond_v%setV( v_v )
                !
            elseif( param_type == LOG_10) then
                !
                v_h = log10( real( v_h, kind=prec ) )
                call self%cell_cond_h%setV( v_h )
                !
                v_v = log10( real( v_v, kind=prec ) )
                call self%cell_cond_v%setV( v_v )
                !
            endif
            !
        elseif( param_type == LINEAR ) then
            !
            if( self%param_type == LOGE ) then
                !
                v_h = exp( v_h )
                call self%cell_cond_h%setV( v_h )
                !
                v_v = exp( v_v )
                call self%cell_cond_v%setV( v_v )
                !
            elseif( self%param_type == LOG_10 ) then
                !
                v_h = exp( v_h * log(10.) )
                call self%cell_cond_h%setV( v_h )
                !
                v_v = exp( v_v * log(10.) )
                call self%cell_cond_v%setV( v_v )
                !
            endif
            !
        elseif( ( self%param_type == LOGE ) .AND. ( param_type == LOG_10 ) ) then
            !
            v_h = v_h / log(10.)
            call self%cell_cond_h%setV( v_h )
            !
            v_v = v_v / log(10.)
            call self%cell_cond_v%setV( v_v )
            !
        elseif( ( self%param_type == LOG_10 ) .AND. ( param_type == LOGE ) ) then
            !
            v_h = v_h * log(10.)
            call self%cell_cond_h%setV( v_h )
            !
            v_v = v_v * log(10.)
            call self%cell_cond_v%setV( v_v )
            !
        else
            stop "Error: setType_ModelParameterCell_SG_VTI > Unknown param_type."
        endif
        !
        self%param_type = param_type 
        !
    end subroutine setType_ModelParameterCell_SG_VTI
    !
    !> No subroutine briefing
    !
    subroutine print_ModelParameterCell_SG_VTI( self )
        implicit none
        !
        class( ModelParameterCell_SG_VTI_t ), intent( in ) :: self
        !
        write( *, * ) "ModelParameterCell_SG_VTI_t:", self%mKey, self%air_cond, self%param_type, &
        self%is_allocated, self%param_grid%nx, self%param_grid%ny, self%param_grid%nz, self%param_grid%nzAir
        !
        !call self%cell_cond%print
        !
    end subroutine print_ModelParameterCell_SG_VTI
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
        !> Convert modelParam to natural log or log10 for output
        !param_type = userparam_type
        !
        !if( self%%is_vti ) then
        !    call getValue_modelParam(m, param_type, self%cell_cond, v_v=ccond_v)
        !else
        !    call getValue_modelParam(m, param_type, self%cell_cond)
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
            rho_h = self%cell_cond_h
            if((index(self%param_type,"LOGE" ) > 0) .OR. (index(self%param_type,"LOG10" ) > 0)) then
                rho_h%v = -self%cell_cond_h%getV()
            elseif(index(self%param_type,"LINEAR" ) > 0) then
                rho_h%v = ONE/self%cell_cond_h%getV()
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
            rho_v = self%cell_cond_v
            if((index(self%param_type,"LOGE" ) > 0) .OR. (index(self%param_type,"LOG10" ) > 0)) then
                rho_v%v = -self%cell_cond_v%getV()
            elseif(index(self%param_type,"LINEAR" ) > 0) then
                rho_v%v = ONE/self%cell_cond_v%getV()
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
