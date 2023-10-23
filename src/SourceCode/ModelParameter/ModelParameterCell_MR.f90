!
!> Basic model parameter class --
!> cell conductivities defined on numerical grid
!> in Earth cells only -- as coded this is specific
!> to usual MR Grid.
!
module ModelParameterCell_MR
    !
    use MetricElements_SG
    use ModelParameterCell_SG
    use Grid3D_MR
    !
    type, extends( ModelParameterCell_t ) :: ModelParameterCell_MR_t
        !
        !> No derived properties
        !
        !> Base rScalar cell_cond treated initialezed
        !> as rScalar3D_MR array here (base class property)
        !
        contains
            !
            final :: ModelParameterCell_MR_dtor
            !
            procedure, public :: nodeCond => nodeCond_ModelParameterCell_MR
            !
            procedure, public :: PDEmapping => PDEmapping_ModelParameterCell_MR
            procedure, public :: dPDEmapping => dPDEmapping_ModelParameterCell_MR
            procedure, public :: dPDEmapping_T => dPDEmapping_T_ModelParameterCell_MR
            !
            !> Dimensioned operations
            procedure, public :: slice1D => slice1D_ModelParameterCell_MR
            procedure, public :: slice2D => slice2D_ModelParameterCell_MR
            !
            procedure, public :: avgModel1D => avgModel1D_ModelParameterCell_MR
            !
            procedure, public :: write => write_ModelParameterCell_MR
            !
    end type ModelParameterCell_MR_t
    !
    interface ModelParameterCell_MR_t
         module procedure ModelParameterCell_MR_ctor_one_cond
         module procedure ModelParameterCell_MR_ctor_all_conds
    end interface ModelParameterCell_MR_t
    !
contains
    !
    !> No subroutine briefing
    !
    function ModelParameterCell_MR_ctor_one_cond( cell_cond, anisotropic_level, param_type, layers ) result( self )
        implicit none
        !
        type( rScalar3D_SG_t ), intent( in ) :: cell_cond
        integer, intent( in ) :: anisotropic_level
        character(:), allocatable, optional, intent( in ) :: param_type
        integer, allocatable, dimension(:), intent( in ) :: layers
        !
        type( ModelParameterCell_MR_t ) :: self
        !
        !write( *, * ) "Constructor ModelParameterCell_SG_ctor_one_cond"
        !
        call self%baseInit
        !
        if( .NOT. present( param_type ) ) then
            self%param_type = LOGE
        else
            self%param_type = trim( param_type )
        endif
        !
        allocate( self%param_grid, source = Grid3D_SG_t( cell_cond%grid%nx, cell_cond%grid%ny, 0, &
        ( cell_cond%grid%nz - cell_cond%grid%nzAir ), cell_cond%grid%dx, cell_cond%grid%dy, &
        cell_cond%grid%dz( cell_cond%grid%nzAir+1 : cell_cond%grid%nz ) ) )
        !
        self%anisotropic_level = anisotropic_level
        !
        allocate( self%cell_cond( anisotropic_level ) )
        !
        self%cell_cond(1) = cell_cond
        !
        self%cell_cond(1)%grid => self%param_grid
        !
        if( present( param_type ) ) then
            !
            call self%setSigMap( param_type )
            !
        endif
        !
        self%is_allocated = .TRUE.
        !
    end function ModelParameterCell_MR_ctor_one_cond
    !
    !> No subroutine briefing
    !
    function ModelParameterCell_MR_ctor_all_conds( cell_cond, param_type, layers ) result( self )
        implicit none
        !
        type( rScalar3D_SG_t ), dimension(:), intent( in ) :: cell_cond
        character(:), allocatable, optional, intent( in ) :: param_type
        integer, allocatable, dimension(:), intent( in ) :: layers
        !
        type( ModelParameterCell_MR_t ) :: self
        !
        integer :: i
        !
        !write( *, * ) "Constructor ModelParameterCell_SG_ctor_all_conds"
        !
        call self%baseInit
        !
        if( .NOT. present( param_type ) ) then
            self%param_type = LOGE
        else
            self%param_type = trim( param_type )
        endif
        !
        allocate( self%param_grid, source = Grid3D_SG_t( cell_cond(1)%grid%nx, cell_cond(1)%grid%ny, 0, &
        ( cell_cond(1)%grid%nz - cell_cond(1)%grid%nzAir ), cell_cond(1)%grid%dx, cell_cond(1)%grid%dy, &
        cell_cond(1)%grid%dz( cell_cond(1)%grid%nzAir+1 : cell_cond(1)%grid%nz ) ) )
        !
        self%anisotropic_level = size( cell_cond )
        !
        self%cell_cond = cell_cond
        !
        do i = 1, self%anisotropic_level
            self%cell_cond(i)%grid => self%param_grid
        enddo
        !
        if( present( param_type ) ) then
            !
            call self%setSigMap( param_type )
        !
        endif
        !
        self%is_allocated = .TRUE.
        !
    end function ModelParameterCell_MR_ctor_all_conds
    !
    !> No subroutine briefing
    !
    subroutine ModelParameterCell_MR_dtor( self )
        implicit none
        !
        type( ModelParameterCell_MR_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor ModelParameterCell_MR"
        !
        call self%baseDealloc
        !
        call self%deallocCell
        !
    end subroutine ModelParameterCell_MR_dtor
    !
    !> Map the entire model cells into a single edge Vector_t(e_vec).
    !> Need to implement for VTI ????
    !
    subroutine nodeCond_ModelParameterCell_MR( self, sigma_node )
        implicit none
        !
        class( ModelParameterCell_MR_t ), intent( in ) :: self
        class( Scalar_t ), intent( inout ) :: sigma_node
        !
        integer :: i_grid, nz_air
        type( Grid3D_MR_t ) :: temp_grid_mr, temp_grid_al_mr
        type( rScalar3D_MR_t ) :: sigma_cell_mr, sigma_cell_al_mr
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "PDEmapping_ModelParameterCell_SG > self not allocated" )
        endif
        !
        if( .NOT. sigma_node%is_allocated ) then
            call errStop( "nodeCond_ModelParameterCell_SG > sigma_node not allocated" )
        endif
        !
        !> Grid MR with AirLayers
        temp_grid_al_mr = self%metric%grid
        !
        !> Grid MR without AirLayers
        temp_grid_mr = Grid3D_MR_t( self%param_grid%nx, self%param_grid%ny, &
        self%param_grid%nzAir, self%param_grid%nzEarth, self%param_grid%dx, &
        self%param_grid%dy, self%param_grid%dz, temp_grid_al_mr%cs )
        !
        call self%metric%setGridIndexArrays( temp_grid_mr )
        !
        !> cell cond as MR
        sigma_cell_mr = rScalar3D_MR_t( temp_grid_mr, CELL )
        !
        call sigma_cell_mr%fromSG( self%cell_cond(1) )
        !
        !> cell cond as MR with AirLayers
        sigma_cell_al_mr = rScalar3D_MR_t( self%metric%grid, CELL )
        !
        !> Considering AirLayers just for the first sub_grid
        nz_air = temp_grid_al_mr%NzAir
        !
        sigma_cell_al_mr%sub_scalar(1)%v( :, :, 1:nz_air ) = self%air_cond
        !
        sigma_cell_al_mr%sub_scalar(1)%v( :, :, nz_air+1:sigma_cell_mr%sub_scalar(1)%nz ) = self%sigMap( real( sigma_cell_mr%sub_scalar(1)%v( :, :, : ), kind=prec ) )
        !
        !> sigMapping for the next sub-grids
        do i_grid = 2, size( sigma_cell_mr%sub_scalar )
            !
            sigma_cell_al_mr%sub_scalar( i_grid )%v = self%sigMap( real( sigma_cell_mr%sub_scalar( i_grid )%v, kind=prec ) )
            !
        enddo
        !
        call sigma_cell_al_mr%mult( self%metric%v_cell )
        !
        call sigma_cell_al_mr%toNode( sigma_node, .TRUE. )
        !
    end subroutine nodeCond_ModelParameterCell_MR
    !
    !> Map the entire model cells into a single edge Vector_t(e_vec).
    !
    subroutine PDEmapping_ModelParameterCell_MR( self, e_vec )
        implicit none
        !
        class( ModelParameterCell_MR_t ), intent( in ) :: self
        class( Vector_t ), intent( inout ) :: e_vec
        !
        integer :: i_grid, nz_air
        type( Grid3D_MR_t ) :: temp_grid_mr
        type( rScalar3D_MR_t ) :: sigma_cell_mr, sigma_cell_al_mr
        class( Vector_t ), allocatable :: e_vol
        complex( kind=prec ), allocatable, dimension(:) :: e_vec_v, e_vol_v
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "PDEmapping_ModelParameterCell_SG > self not allocated" )
        endif
        !
        if( .NOT. e_vec%is_allocated ) then
            call errStop( "PDEmapping_ModelParameterCell_SG > e_vec not allocated" )
        endif
        !
        select type( grid => self%metric%grid )
            !
            class is( Grid3D_MR_t )
                !
                temp_grid_mr = Grid3D_MR_t( self%param_grid%nx, self%param_grid%ny, &
                self%param_grid%nzAir, self%param_grid%nzEarth, self%param_grid%dx, &
                self%param_grid%dy, self%param_grid%dz, grid%cs )
                !
            class default
                call errStop( "PDEmapping_ModelParameterCell_MR > Grid must be MR!" )
            !
        end select
        !
        !> cell cond as MR without AirLayers
        sigma_cell_mr = rScalar3D_MR_t( temp_grid_mr, CELL )
        !
        call sigma_cell_mr%fromSG( self%cell_cond(1) )
        !
        !> cell cond as MR with AirLayers
        sigma_cell_al_mr = rScalar3D_MR_t( self%metric%grid, CELL )
        !
        !> Considering AirLayers just for the first sub_grid
        nz_air = self%metric%grid%NzAir
        !
        sigma_cell_al_mr%sub_scalar(1)%v( :, :, 1:nz_air ) = self%air_cond
        !
        sigma_cell_al_mr%sub_scalar(1)%v( :, :, nz_air+1:sigma_cell_mr%sub_scalar(1)%nz ) = self%sigMap( real( sigma_cell_mr%sub_scalar(1)%v, kind=prec ) )
        !
        !> sigMapping for the next sub-grids
        do i_grid = 2, size( sigma_cell_mr%sub_scalar )
            !
            sigma_cell_al_mr%sub_scalar( i_grid )%v = self%sigMap( real( sigma_cell_mr%sub_scalar( i_grid )%v, kind=prec ) )
            !
        enddo
        !
        !> E_VEC: Boundaries set to zero.
        call sigma_cell_al_mr%mult( self%metric%v_cell )
        !
        call e_vec%sumCells( sigma_cell_al_mr )
        !
        e_vec_v = e_vec%getArray()
        !
        e_vec_v( e_vec%indBoundary() ) = C_ZERO
        !
        call e_vec%setArray( e_vec_v )
        !
        !> E_VOL: Boundaries set to one to avoid NaNs at the division in the end.
        call self%metric%createVector( real_t, EDGE, e_vol )
        !
        call e_vol%sumCells( self%metric%v_cell )
        !
        e_vol_v = e_vol%getArray()
        !
        e_vol_v( e_vol%indBoundary() ) = C_ONE
        !
        call e_vol%setArray( e_vol_v )
        !
        call e_vec%div( e_vol )
        !
        deallocate( e_vol )
        !
    end subroutine PDEmapping_ModelParameterCell_MR
    !
    !> Map the perturbation between two models onto a single Vector_t(e_vec).
    !
    subroutine dPDEmapping_ModelParameterCell_MR( self, dsigma, e_vec )
        implicit none
        !
        class( ModelParameterCell_MR_t ), intent( in ) :: self
        class( ModelParameter_t ), intent( in ) :: dsigma
        class( Vector_t ), allocatable, intent( inout ) :: e_vec
        !
        call errStop( "dPDEmapping_ModelParameterCell_MR > must be implemented!" )
        !
    end subroutine dPDEmapping_ModelParameterCell_MR
    !
    !> Transpose the perturbation represented in a Vector_t(e_vec), to a new dsigma model.
    !
    subroutine dPDEmapping_T_ModelParameterCell_MR( self, e_vec, dsigma )
        implicit none
        !
        class( ModelParameterCell_MR_t ), intent( in ) :: self
        class( Vector_t ), intent( in ) :: e_vec
        class( ModelParameter_t ), allocatable, intent( inout ) :: dsigma
        !
        call errStop( "dPDEmapping_T_ModelParameterCell_MR > must be implemented!" )
        !
    end subroutine dPDEmapping_T_ModelParameterCell_MR
    !
    !> No subroutine briefing
    !
    function slice1D_ModelParameterCell_MR( self, ix, iy ) result( model_param_1D )
        implicit none
        !
        class( ModelParameterCell_MR_t ), intent( in ) :: self
        integer, intent( in ) :: ix, iy
        !
        type( ModelParameter1D_t ) ::  model_param_1D 
        !
        real( kind=prec ), allocatable, dimension(:) :: cond_slice
        !
        model_param_1D = ModelParameter1D_t( self%metric%grid%slice1D() )
        !
        allocate( cond_slice( model_param_1D%grid%nz ) )
        !
        cond_slice = self%sigMap( real( self%cell_cond(1)%v( ix, iy, : ), kind=prec ) )
        !
        call model_param_1D%setConductivity( cond_slice, self%air_cond, self%param_type, self%mKey )
        !
        deallocate( cond_slice )
        !
    end function slice1D_ModelParameterCell_MR
    !
    !> No subroutine briefing
    !
    function avgModel1D_ModelParameterCell_MR( self ) result( model_param_1D )
        implicit none
        !
        class( ModelParameterCell_MR_t ), intent( in ) :: self
        !
        type( ModelParameter1D_t ) ::  model_param_1D 
        !
        call errStop( "avgModel1D_ModelParameterCell_MR > must be implemented!" )
        !
    end function avgModel1D_ModelParameterCell_MR
    !
    !> No subroutine briefing
    !
    function slice2D_ModelParameterCell_MR( self, axis, j ) result( m2D )
        implicit none
        !
        class( ModelParameterCell_MR_t ), intent( in ) :: self
        integer, intent( in ) :: axis, j
        !
        type( ModelParameter2D_t ) :: m2D 
        !
        call errStop( "slice2D_ModelParameterCell_MR > must be implemented!" )
        !
    end function slice2D_ModelParameterCell_MR
    !
    !> opens cfile on unit ioModelParam, writes out object of
    !> type modelParam in Weerachai Siripunvaraporn"s format,
    !> closes file.
    !
    subroutine write_ModelParameterCell_MR( self, file_name, comment )
        implicit none
        !
        class( ModelParameterCell_MR_t ), intent( in ) :: self
        character(*), intent( in ) :: file_name
        character(*), intent( in ), optional :: comment
        !
        call errStop( "write_ModelParameterCell_MR > must be implemented!" )
        !
    end subroutine write_ModelParameterCell_MR
    !
end Module ModelParameterCell_MR
