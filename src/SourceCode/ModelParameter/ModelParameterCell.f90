!
!> Basic model parameter class --
!> cell conductivities defined on numerical grid
!> in Earth cells only -- as coded this is specific
!> to usual CSG.
!
module ModelParameterCell
    !
    use ModelParameter
    use FileUnits
    use Grid3D_SG
    !
    type, extends( ModelParameter_t ) :: ModelParameterCell_t
        !
        class( Grid_t ), allocatable :: param_grid
        !
        type( GenScalar_t ), allocatable, dimension(:) :: cell_cond
        !
        contains
            !
            final :: ModelParameterCell_dtor
            !
            procedure, public :: getOneCond => getOneCond_ModelParameterCell
            procedure, public :: getAllCond => getAllCond_ModelParameterCell
            !
            procedure, public :: setOneCond => setOneCond_ModelParameterCell
            procedure, public :: setAllCond => setAllCond_ModelParameterCell
            !
            procedure, public :: zeros => zeros_ModelParameterCell
            !
            procedure, public :: copyFrom => copyFrom_ModelParameterCell
            !
            procedure, public :: countModel => countModel_ModelParameterCell
            !
            procedure, public :: dotProd => dotProd_ModelParameterCell
            !
            procedure, public :: linComb => linComb_ModelParameterCell
            !
            procedure, public :: cellToNode => cellToNode_ModelParameterCell
            !
            procedure, public :: PDEmapping => PDEmapping_ModelParameterCell
            procedure, public :: dPDEmapping => dPDEmapping_ModelParameterCell
            procedure, public :: dPDEmapping_T => dPDEmapping_T_ModelParameterCell
            !
            procedure, public :: slice1D => slice1D_ModelParameterCell
            procedure, public :: slice2D => slice2D_ModelParameterCell
            !
            procedure, public :: avgModel1D => avgModel1D_ModelParameterCell
            !
            procedure, public :: setType => setType_ModelParameterCell
            !
            procedure, public :: write => write_ModelParameterCell
            !
            procedure, public :: print => print_ModelParameterCell
            !
    end type ModelParameterCell_t
    !
    interface ModelParameterCell_t
         module procedure ModelParameterCell_ctor_one_cond
         module procedure ModelParameterCell_ctor_all_conds
    end interface ModelParameterCell_t
    !
contains
    !
    !> No subroutine briefing
    !
    function ModelParameterCell_ctor_one_cond( grid, cell_cond, anisotropic_level, param_type ) result( self )
        implicit none
        !
        class( Grid_t ), intent( in ) :: grid
        class( Scalar_t ), allocatable, intent( in ) :: cell_cond
        integer, intent( in ) :: anisotropic_level
        character(:), allocatable, optional, intent( in ) :: param_type
        !
        type( ModelParameterCell_t ) :: self
        !
        integer :: nzAir
        !
        !write( *, * ) "Constructor ModelParameterCell_ctor_one_cond"
        !
        call self%baseInit
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
        self%anisotropic_level = anisotropic_level
        !
        allocate( self%cell_cond( anisotropic_level ) )
        !
        allocate( self%cell_cond(1)%s, source = cell_cond )
        !
        self%cell_cond(1)%s%store_state = compound
        !
        if( present( param_type ) ) then
            !
            call self%setsigMap( param_type )
        !
        endif
        !
        self%is_allocated = .TRUE.
        !
    end function ModelParameterCell_ctor_one_cond
    !
    !> No subroutine briefing
    !
    function ModelParameterCell_ctor_all_conds( grid, cell_cond, param_type ) result( self )
        implicit none
        !
        class( Grid_t ), intent( in ) :: grid
        type( GenScalar_t ), allocatable, dimension(:), intent( in ) :: cell_cond
        character(:), allocatable, optional, intent( in ) :: param_type
        !
        type( ModelParameterCell_t ) :: self
        !
        integer :: i, nzAir
        !
        !write( *, * ) "Constructor ModelParameterCell_ctor_all_conds"
        !
        call self%baseInit
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
                ( grid%nz - grid%nzAir ), grid%dx, grid%dy, grid%dz( grid%nzAir+1:grid%nz ) ) )
        !
        self%anisotropic_level = size( cell_cond )
        !
        self%cell_cond = cell_cond
        !
        !do i = 1, self%anisotropic_level
            !self%cell_cond(i)%s%store_state = compound
        !enddo
        !
        if( present( param_type ) ) then
            !
            call self%setsigMap( param_type )
        !
        endif
        !
        self%is_allocated = .TRUE.
        !
    end function ModelParameterCell_ctor_all_conds
    !
    !> No subroutine briefing
    subroutine ModelParameterCell_dtor( self )
        implicit none
        !
        type( ModelParameterCell_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor ModelParameterCell"
        !
        if( allocated( self%param_grid ) ) deallocate( self%param_grid )
        !
        !> CRASHING ????
        if( allocated( self%cell_cond ) ) deallocate( self%cell_cond )
        !
    end subroutine ModelParameterCell_dtor
    !
    !> No subroutine briefing
    !
    function slice1D_ModelParameterCell( self, ix, iy ) result( model_param_1D )
        implicit none
        !
        class( ModelParameterCell_t ), intent( inout ) :: self
        integer, intent( in ) :: ix, iy
        !
        type( ModelParameter1D_t ) ::  model_param_1D 
        !
        real( kind=prec ), allocatable, dimension(:) :: cond_slice
        complex( kind=prec ), allocatable, dimension(:,:,:) :: self_cond_v
        !
        model_param_1D = ModelParameter1D_t( self%metric%grid%slice1D() )
        !
        allocate( cond_slice( model_param_1D%grid%nz ) )
        !
        self_cond_v = self%cell_cond(1)%s%getV()
        !
        cond_slice = self%sigMap( real( self_cond_v( ix, iy, : ), kind=prec ) )
        !
        call model_param_1D%setConductivity( cond_slice, self%air_cond, self%param_type, self%mKey )
        !
        deallocate( cond_slice )
        !
    end function slice1D_ModelParameterCell
    !
    !> No subroutine briefing
    !
    function avgModel1D_ModelParameterCell( self ) result( model_param_1D )
        implicit none
        !
        class( ModelParameterCell_t ), intent( inout ) :: self
        !
        type( ModelParameter1D_t ) ::  model_param_1D 
        !
        real( kind=prec ), allocatable, dimension(:) :: cond_slice
        complex( kind=prec ), allocatable, dimension(:,:,:) :: self_cond_v
        real( kind=prec ) :: wt, temp_sigma_value
        integer :: i, j, k
        !
        model_param_1D = ModelParameter1D_t( self%metric%grid%slice1D() )
        !
        allocate( cond_slice( self%metric%grid%nzEarth ) )
        !
        self_cond_v = self%cell_cond(1)%s%getV()
        !
        do k = 1, self%metric%grid%nzEarth
            !
            wt = R_ZERO
            temp_sigma_value = R_ZERO
            do i = 1, self%metric%grid%Nx
                do j = 1, self%metric%grid%Ny
                    wt = wt + self%metric%grid%dx(i) * self%metric%grid%dy(j)
                    temp_sigma_value = temp_sigma_value + self_cond_v( i, j, k ) * &
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
    end function avgModel1D_ModelParameterCell
    !
    !> No subroutine briefing
    !
    function slice2D_ModelParameterCell( self, axis, j ) result( m2D )
        implicit none
        !
        class( ModelParameterCell_t ), intent( inout ) :: self
        integer, intent( in ) :: axis, j
        !
        type( ModelParameter2D_t ) :: m2D 
        !
        character(:), allocatable :: param_type
        real( kind=prec ), allocatable, dimension(:,:) :: cond_slice
        complex( kind=prec ), allocatable, dimension(:,:,:) :: self_cond_v
        !
        param_type = LINEAR
        !
        m2D = ModelParameter2D_t( self%metric%grid%slice2D() )
        !
        self_cond_v = self%cell_cond(1)%s%getV()
        !
        allocate( cond_slice( self%metric%grid%ny, self%metric%grid%nzEarth ) )
        !
        if( axis == 1 ) then
            cond_slice = self%sigMap( real( self_cond_v(j,:,:), kind=prec ) )
        elseif( axis == 2 ) then
            cond_slice = self%sigMap( real( self_cond_v(:,j,:), kind=prec ) )
        elseif( axis == 3 ) then
            cond_slice = self%sigMap( real( self_cond_v(:,:,j), kind=prec ) )
        else
            call errStop( "slice2D_ModelParameterCell > wrong axis" )
        endif
        !
        call m2D%setConductivity( cond_slice, self%air_cond, param_type, self%mKey )
        !
        deallocate( cond_slice )
        !
    end function slice2D_ModelParameterCell
    !
    !> No function briefing
    !
    function getOneCond_ModelParameterCell( self, i_cond ) result( cell_cond )
        implicit none
        !
        class( ModelParameterCell_t ), intent( in ) :: self
        integer, intent( in ) :: i_cond
        !
        class( Scalar_t ), allocatable :: cell_cond
        !
        if( i_cond .GT. self%anisotropic_level ) then
            !
            call errStop( "getOneCond_ModelParameterCell > conductivity level too high" )
            !
        endif
        !
        allocate( cell_cond, source = self%cell_cond( i_cond )%s )
        !
    end function getOneCond_ModelParameterCell
    !
    !> No function briefing
    !
    function getAllCond_ModelParameterCell( self ) result( cell_cond )
        implicit none
        !
        class( ModelParameterCell_t ), intent( in ) :: self
        !
        type( GenScalar_t ), allocatable, dimension(:) :: cell_cond
        !
        cell_cond = self%cell_cond
        !
    end function getAllCond_ModelParameterCell
    !
    !> No interface subroutine briefing
    !
    subroutine setOneCond_ModelParameterCell( self, cell_cond, i_cond )
        implicit none
        !
        class( ModelParameterCell_t ), intent( inout ) :: self
        class( Scalar_t ), intent( in ) :: cell_cond
        integer, intent( in ) :: i_cond
        !
        if( .NOT. cell_cond%is_allocated ) then
            call errStop( "setOneCond_ModelParameterCell > cell_cond not allocated" )
        endif
        !
        if( i_cond .LE. self%anisotropic_level ) then
            !
            if( allocated( self%cell_cond( i_cond )%s ) ) deallocate( self%cell_cond( i_cond )%s )
            allocate( self%cell_cond( i_cond )%s, source = cell_cond )
            !
        else
            !
            call errStop( "setOneCond_ModelParameterCell > Unsupport general anisotropy yet" )
            !
        endif
        !
    end subroutine setOneCond_ModelParameterCell
    !
    !> No interface subroutine briefing
    !
    subroutine setAllCond_ModelParameterCell( self, cell_cond )
        implicit none
        !
        class( ModelParameterCell_t ), intent( inout ) :: self
        type( GenScalar_t ), allocatable, dimension(:), intent( in ) :: cell_cond
        !
        integer :: i
        !
        do i = 1, self%anisotropic_level
            !
            if( .NOT. cell_cond(i)%s%is_allocated ) then
                !
                call errStop( "setAllCond_ModelParameterCell > cell_cond has no V" )
                !
            elseif( .NOT. self%cell_cond(i)%s%is_allocated )then
                !
                call errStop( "setAllCond_ModelParameterCell > self has no V" )
                !
            endif
            !
            call self%setCond( cell_cond(i)%s, i )
            !
        enddo
        !
    end subroutine setAllCond_ModelParameterCell
    !
    !> No subroutine briefing
    subroutine zeros_ModelParameterCell( self )
        implicit none
        !
        class( ModelParameterCell_t ), intent( inout ) :: self
        !
        integer :: i
        !
        do i = 1, self%anisotropic_level
            !
            call self%cell_cond(i)%s%zeros
            !
        enddo
        !
    end subroutine zeros_ModelParameterCell
    !
    !> No subroutine briefing
    subroutine copyFrom_ModelParameterCell( self, rhs )
        implicit none
        !
        class( ModelParameterCell_t ), intent( inout ) :: self
        class( ModelParameter_t ), intent( in ) :: rhs
        !
        if( .NOT. rhs%is_allocated ) then
            call errStop( "copyFrom_ModelParameterCell > rhs not allocated." )
        endif
        !
        select type( rhs )
            !
            class is( ModelParameterCell_t )
                !
                self%metric => rhs%metric
                !
                self%mKey = rhs%mKey
                !
                self%anisotropic_level = rhs%anisotropic_level
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
                call errStop( "copyFrom_ModelParameterCell > Unclassified rhs." )
            !
        end select
        !
    end subroutine copyFrom_ModelParameterCell
    !
    !> ????
    function countModel_ModelParameterCell( self ) result( counter )
        implicit none
        !
        class( ModelParameterCell_t ), intent( in ) :: self
        !
        integer :: i, counter, nx, ny, nz, nzAir, nz_earth
        !
        counter = 0
        !
        do i = 1, self%anisotropic_level
            !
            if( .NOT. self%cell_cond(i)%s%is_allocated ) then
                write( *, * ) "Error: countModel_ModelParameterCell > cell_cond (", i, ") not allocated!"
                stop
            endif
            !
            call self%cell_cond(i)%s%grid%getDimensions( nx, ny, nz, nzAir )
            nz_earth = nz - nzAir
            !
            counter = counter + self%cell_cond(i)%s%Nx * self%cell_cond(i)%s%Ny * nz_earth
            !
        enddo
        !
    end function countModel_ModelParameterCell
    !
    !>
    subroutine linComb_ModelParameterCell( self, a1, a2, rhs )
        implicit none
        !
        class( ModelParameterCell_t ), intent( inout ) :: self
        real( kind=prec ), intent( in ) :: a1, a2
        class( ModelParameter_t ), intent( inout ) :: rhs
        !
        integer :: i
        complex( kind=prec ), allocatable :: v(:,:,:)
        !
        select type( rhs )
            !
            class is( ModelParameterCell_t )
                !
                do i = 1, self%anisotropic_level
                    !
                    if( self%cell_cond(i)%s%isCompatible( rhs%cell_cond(i)%s ) ) then
                        !
                        v = a1 * self%cell_cond(i)%s%getV() + a2 * rhs%cell_cond(i)%s%getV()
                        !
                        call self%cell_cond(i)%s%setV( v )
                        !
                    else
                        write( *, * ) "Error: linComb_ModelParameterCell > Incompatible rhs cell_cond (", i, ")!"
                        stop
                    endif
                    !
                enddo
                !
            class default
                call errStop( "linComb_ModelParameterCell > undefined rhs" )
            !
        end select
        !
        !> ????
        !self%air_cond = rhs%air_cond
        !
    end subroutine linComb_ModelParameterCell
    !
    !> No subroutine briefing
    !
    function dotProd_ModelParameterCell( self, rhs ) result( rvalue )
        implicit none
        !
        class( ModelParameterCell_t ), intent( inout ) :: self
        class( ModelParameter_t ), intent( inout ) :: rhs
        !
        real( kind=prec ) :: rvalue
        !
        integer :: i
        !
        rvalue = R_ZERO
        !
        select type( rhs )
            !
            class is( ModelParameterCell_t )
                !
                do i = 1, self%anisotropic_level
                    !
                    if( self%cell_cond(i)%s%isCompatible( rhs%cell_cond(i)%s ) ) then
                        !
                        rvalue = rvalue + sum( self%cell_cond(i)%s%getV() * rhs%cell_cond(i)%s%getV() )
                        !
                    else
                        write( *, * ) "Error: dotProd_ModelParameterCell > Incompatible rhs cell_cond (", i, ")!"
                        stop
                    endif
                    !
                enddo
                !
            class default
                call errStop( "dotProd_ModelParameterCell > Unclassified rhs" )
            !
        end select
        !
    end function dotProd_ModelParameterCell
    !
    !> Map the entire model cells into a single edge Vector_t (e_vec).
    !> Need to implement for VTI ????
    !
    subroutine cellToNode_ModelParameterCell( self, sigma_node )
        implicit none
        !
        class( ModelParameterCell_t ), intent( in ) :: self
        class( Scalar_t ), allocatable, intent( inout ) :: sigma_node
        !
        integer :: k0, k1, k2
        class( Scalar_t ), allocatable:: sigma_cell, sigma_cond
        type( GenScalar_t ), allocatable, dimension(:) :: sigma_cells
        complex( kind=prec ), allocatable, dimension(:,:,:) :: sigma_cond_v, sigma_cell_v
        !
        if( .NOT. sigma_node%is_allocated ) then
            call errStop( "cellToNode_ModelParameterCell > sigma_node not allocated" )
        endif
        !
        sigma_cells = self%getCond()
        !
        !> Just get the proper conductivity into sigma_cond
        !> Directly use sigma_cells(1)%s%getV() causes compiler error ????
        allocate( sigma_cond, source = sigma_cells(1)%s )
        sigma_cond_v = sigma_cond%getV()
        deallocate( sigma_cond )
        !
        k0 = self%metric%grid%NzAir
        k1 = k0 + 1
        k2 = self%metric%grid%Nz
        !
        call self%metric%createScalar( real_t, CELL, sigma_cell )
        !
        sigma_cell_v = sigma_cell%getV()
        !
        sigma_cell_v( :, :, k1:k2 ) = self%sigMap( real( sigma_cond_v, kind=prec ) )
        !
        call sigma_cell%setV( sigma_cell_v )
        !
        call sigma_cell%mult( self%metric%v_cell )
        !
        call sigma_cell%sumCell( sigma_node, .TRUE. )
        !
        call sigma_node%mult( cmplx( 0.125_prec, 0.0, kind=prec ) )
        !
        !> This should be the right thing, but no ????
        !call sigma_node%div( self%metric%v_node )
        !
        deallocate( sigma_cell )
        !
    end subroutine cellToNode_ModelParameterCell
    !
    !> Map the entire model cells into a single edge Vector_t (e_vec).
    !
    subroutine PDEmapping_ModelParameterCell( self, e_vec )
        implicit none
        !
        class( ModelParameterCell_t ), intent( in ) :: self
        class( Vector_t ), allocatable, intent( inout ) :: e_vec
        !
        type( GenScalar_t ), allocatable, dimension(:) :: sigma_cells
        class( Scalar_t ), allocatable :: sigma_cell
        complex( kind=prec ), allocatable, dimension(:,:,:) :: sigma_cell_v
        integer :: i, k0, k1, k2
        !
        if( .NOT. e_vec%is_allocated ) then
            call errStop( "PDEmapping_ModelParameterCell > e_vec not allocated yet" )
        endif
        !
        allocate( sigma_cells( self%anisotropic_level ) )
        !
        k0 = self%metric%grid%nzAir
        k1 = k0 + 1
        k2 = self%metric%grid%Nz
        !
        do i = 1, self%anisotropic_level
            !
            !> Create and initialize e_vec with zeros
            call self%metric%createScalar( real_t, CELL, sigma_cell )
            !
            call sigma_cell%zeros
            !
            sigma_cell_v = sigma_cell%getV()
            !
            sigma_cell_v( :, :, 1:k0 ) = self%air_cond
            !
            sigma_cell_v( :, :, k1:k2 ) = self%sigMap( real( self%cell_cond(i)%s%getV(), kind=prec ) )
            !
            call sigma_cell%setV( sigma_cell_v )
            !
            call sigma_cell%mult( self%metric%v_cell )
            !
            allocate( sigma_cells(i)%s, source = sigma_cell )
            !
        enddo
        !
        !> Call due avgCells based on anisotropic_level
        if( self%anisotropic_level == 1 ) then
            !
            call e_vec%avgCells( sigma_cells(1)%s )
            !
        elseif( self%anisotropic_level == 2 ) then
            !
            call e_vec%avgCells( sigma_cells(1)%s, sigma_cells(2)%s )
            !
        else
            call errStop( "PDEmapping_ModelParameterCell > unsupported anisotropy level" )
        endif
        !
        call e_vec%div( self%metric%v_edge )
        !
    end subroutine PDEmapping_ModelParameterCell
    !
    !> Map the perturbation between two models onto a single Vector_t (e_vec).
    !
    subroutine dPDEmapping_ModelParameterCell( self, dsigma, e_vec )
        implicit none
        !
        class( ModelParameterCell_t ), intent( inout ) :: self
        class( ModelParameter_t ), intent( in ) :: dsigma
        class( Vector_t ), allocatable, intent( inout ) :: e_vec
        !
        type( GenScalar_t ), allocatable, dimension(:) :: sigma_cells
        class( Scalar_t ), allocatable :: sigma_cell, dsigma_cond
        complex( kind=prec ), allocatable, dimension(:,:,:) :: sigma_cell_v
        character( len=5 ), parameter :: JOB = "DERIV"
        integer :: i, k0, k1, k2
        !
        if( .NOT. dsigma%is_allocated ) then
            call errStop( "dPDEmapping_ModelParameterCell > dsigma not allocated" )
        endif
        !
        if( .NOT. e_vec%is_allocated ) then
            call errStop( "dPDEmapping_ModelParameterCell > e_vec not allocated" )
        endif
        !
        call e_vec%zeros
        !
        k0 = self%metric%grid%NzAir
        k1 = k0 + 1
        k2 = self%metric%grid%Nz
        !
        allocate( sigma_cells( self%anisotropic_level ) )
        !
        do i = 1, self%anisotropic_level
            !
            !> Create and initialize e_vec with zeros
            call self%metric%createScalar( real_t, CELL, sigma_cell )
            !
            call sigma_cell%zeros
            !
            sigma_cell_v = sigma_cell%getV()
            !
            sigma_cell_v( :, :, k1:k2 ) = self%sigMap( real( self%cell_cond(i)%s%getV(), kind=prec ), JOB )
            !
            allocate( dsigma_cond, source = dsigma%getCond(i) )
            !
            sigma_cell_v( :, :, k1:k2 ) = sigma_cell_v( :, :, k1:k2 ) * dsigma_cond%getV()
            !
            deallocate( dsigma_cond )
            !
            call sigma_cell%setV( sigma_cell_v )
            !
            call sigma_cell%mult( self%metric%v_cell )
            !
            allocate( sigma_cells(i)%s, source = sigma_cell )
            !
        enddo
        !
        !> Call specific avgCells based on anisotropic_level
        if( self%anisotropic_level == 1 ) then
            !
            call e_vec%avgCells( sigma_cells(1)%s )
            !
        elseif( self%anisotropic_level == 2 ) then
            !
            call e_vec%avgCells( sigma_cells(1)%s, sigma_cells(2)%s )
            !
        else
            call errStop( "dPDEmapping_ModelParameterCell > unsupported anisotropy level" )
        endif
        !
        call e_vec%div( self%metric%v_edge )
        !
    end subroutine dPDEmapping_ModelParameterCell
    !
    !> Transpose the perturbation represented in a Vector_t (e_vec), to a new dsigma model.
    !
    subroutine dPDEmapping_T_ModelParameterCell( self, e_vec, dsigma )
        implicit none
        !
        class( ModelParameterCell_t ), intent( inout ) :: self
        class( Vector_t ), intent( in ) :: e_vec
        class( ModelParameter_t ), allocatable, intent( out ) :: dsigma
        !
        class( Scalar_t ), allocatable :: sigma_cell, dsigma_cond
        class( Vector_t ), allocatable :: e_vec_interior
        type( GenScalar_t ), allocatable, dimension(:) :: sigma_cells
        complex( kind=prec ), allocatable, dimension(:,:,:) :: self_cond_v, dsigma_cond_v, sigma_cell_v
        character( len=5 ), parameter :: JOB = "DERIV"
        integer :: i, k0, k1, k2
        !
        if( .NOT. e_vec%is_allocated ) then
            call errStop( "dPDEmapping_T_ModelParameterCell > e_vec not allocated" )
        endif
        !
        !> e_vec
        call e_vec%interior( e_vec_interior )
        !
        call e_vec_interior%div( self%metric%v_edge )
        !
        call e_vec_interior%mult( cmplx( 0.25_prec, 0.0, kind=prec ) )
        !
        !> sigma_cells
        allocate( sigma_cells( self%anisotropic_level ) )
        !
        !> Call specific sumEdges based on anisotropic_level
        if( self%anisotropic_level == 1 ) then
            !
            call e_vec_interior%sumEdges( sigma_cells(1)%s, .TRUE. )
            !
        elseif( self%anisotropic_level == 2 ) then
            !
            call e_vec_interior%sumEdges( sigma_cells(1)%s, sigma_cells(2)%s, .TRUE. )
            !
        else
            call errStop( "dPDEmapping_T_ModelParameterCell > unsupported anisotropy level" )
        endif
        !
        deallocate( e_vec_interior )
        !
        k0 = self%metric%grid%NzAir
        k1 = k0 + 1
        k2 = self%metric%grid%Nz
        !
        allocate( dsigma, source = ModelParameterCell_t( self%param_grid, self%cell_cond, self%param_type ) )
        !
        do i = 1, self%anisotropic_level
            !
            allocate( dsigma_cond, source = dsigma%getCond(i) )
            !
            call dsigma_cond%zeros
            !
            dsigma_cond_v = dsigma_cond%getV()
            !
            self_cond_v = self%cell_cond(i)%s%getV()
            !
            dsigma_cond_v = self%sigMap( real( self_cond_v, kind=prec ), JOB )
            !
            call sigma_cells(i)%s%mult( self%metric%v_cell )
            !
            allocate( sigma_cell, source = sigma_cells(i)%s )
            !
            sigma_cell_v = sigma_cell%getV()
            !
            deallocate( sigma_cell )
            !
            dsigma_cond_v = dsigma_cond_v * sigma_cell_v( :, :, k1:k2 )
            !
            call dsigma_cond%setV( dsigma_cond_v )
            !
            call dsigma%setCond( dsigma_cond, i )
            !
            deallocate( dsigma_cond )
            !
        enddo
        !
    end subroutine dPDEmapping_T_ModelParameterCell
    !
    !> No subroutine briefing
    !
    subroutine setType_ModelParameterCell( self, param_type )
        implicit none
        !
        class( ModelParameterCell_t ), intent( inout ) :: self
        character(:), allocatable, intent( in ) :: param_type
        !
        integer :: i
        real( kind=prec ), allocatable, dimension(:,:,:) :: self_cond_v
        complex( kind=prec ), allocatable, dimension(:,:,:) :: self_cond_cv
        !
        if( .NOT. self%is_allocated ) then
                call errStop( "setType_ModelParameterCell > Self not allocated." )
        endif
        !
        do i = 1, self%anisotropic_level
            !
            self_cond_v = real( self%cell_cond(i)%s%getV(), kind=prec )
            !
            if( trim( param_type ) .EQ. trim( self%param_type ) ) then
                ! Nothing to be done
            elseif( self%param_type == "" ) then
                self%param_type = trim( param_type )
            elseif( self%param_type == LINEAR ) then
                !
                if( param_type == LOGE ) then
                    !
                    self_cond_v = log( self_cond_v )
                    !
                elseif( param_type == LOG_10) then
                    !
                    self_cond_v = log10( self_cond_v )
                    !
                endif
                !
            elseif( param_type == LINEAR ) then
                !
                if( self%param_type == LOGE ) then
                    !
                    self_cond_v = exp( self_cond_v )
                    !
                elseif( self%param_type == LOG_10 ) then
                    !
                    self_cond_v = exp( self_cond_v * log(10.) )
                    !
                endif
                !
            elseif( ( self%param_type == LOGE ) .AND. ( param_type == LOG_10 ) ) then
                !
                self_cond_v = self_cond_v / log(10.)
                !
            elseif( ( self%param_type == LOG_10 ) .AND. ( param_type == LOGE ) ) then
                !
                self_cond_v = self_cond_v * log(10.)
                !
            else
                call errStop( "setType_ModelParameterCell > Unknown param_type." )
            endif
            !
            self_cond_cv = cmplx( self_cond_v, 0.0, kind=prec )
            !
            call self%cell_cond(i)%s%setV( self_cond_cv )
            !
        enddo
        !
        self%param_type = param_type 
        !
    end subroutine setType_ModelParameterCell
    !
    !> No subroutine briefing
    !
    subroutine print_ModelParameterCell( self )
        implicit none
        !
        class( ModelParameterCell_t ), intent( inout ) :: self
        !
        integer :: i
        !
        write( *, * ) "ModelParameterCell_t:", self%mKey, self%air_cond, self%param_type, &
        self%is_allocated, self%param_grid%nx, self%param_grid%ny, self%param_grid%nz, self%param_grid%nzAir
        !
        do i = 1, self%anisotropic_level
            !
            call self%cell_cond(i)%s%print
            !
        enddo
        !
    end subroutine print_ModelParameterCell
    !
    !> opens cfile on unit ioModelParam, writes out object of
    !> type modelParam in Weerachai Siripunvaraporn"s format,
    !> closes file.
    !
    subroutine write_ModelParameterCell( self, file_name, comment )
        implicit none
        !
        class( ModelParameterCell_t ), intent( inout ) :: self
        character(*), intent( in ) :: file_name
        character(*), intent( in ), optional :: comment
        !
        integer :: Nx, Ny, NzEarth, ii, i, j, k, ios
        real( kind=prec ), allocatable, dimension(:,:,:) :: self_cond_v
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
            !
            if( self%anisotropic_level == 2 ) then
                !
                write( ioModelParam, * ) " VTI"
                !
            else
                !
                write( ioModelParam, * )
                !
            endif
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
            write( *, * ) self%anisotropic_level
            !
            do ii = 1, self%anisotropic_level
                !
                !> Convert (horizontal) conductivity to resistivity
                self_cond_v = self%cell_cond(ii)%s%getV()
                !
                if( index( self%param_type,"LOGE" ) > 0 .OR. index( self%param_type,"LOG10" ) > 0 ) then
                    self_cond_v = -self_cond_v
                elseif( index(self%param_type,"LINEAR" ) > 0 ) then
                    self_cond_v = ONE / self_cond_v
                endif
                !
                !> Write the (horizontal) resistivity
                !
                write( ioModelParam, * )
                !
                do k = 1, nzEarth
                    do j = 1, Ny
                        do i = Nx, 1, -1
                            write( ioModelParam, "(es13.5)", iostat = ios, advance = "no" ) self_cond_v(i,j,k)
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
            call errStop( "write_ModelParameterCell > Error opening file ["//file_name//"]!" )
        endif
        !
    end subroutine write_ModelParameterCell
    !
end Module ModelParameterCell
