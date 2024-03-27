!
!> Basic model parameter class --
!> cell conductivities defined on numerical grid
!> in Earth cells only -- as coded this is specific
!> to usual MR Grid.
!
module ModelParameterCell_MR
    !
    use MetricElements_MR
    use ModelParameterCell_SG
    use Grid3D_MR
    !
    type, extends( ModelParameterCell_t ) :: ModelParameterCell_MR_t
        !
        !> No derived properties
        !
        !> Base rScalar cell_cond treated initialized
        !> as rScalar3D_MR array here (base class property)
        !
        contains
            !
            final :: ModelParameterCell_MR_dtor
            !
            procedure, public :: PDEmapping => PDEmapping_ModelParameterCell_MR
            procedure, public :: dPDEmapping => dPDEmapping_ModelParameterCell_MR
            procedure, public :: dPDEmapping_T => dPDEmapping_T_ModelParameterCell_MR
            !
            procedure, public :: nodeCond => nodeCond_ModelParameterCell_MR
            !
            procedure, public :: modelToCell => modelToCell_ModelParameterCell_MR
            procedure, public :: cellToModel => cellToModel_ModelParameterCell_MR
            !
            !> Dimensioned operations
            procedure, public :: slice1D => slice1D_ModelParameterCell_MR
            procedure, public :: slice2D => slice2D_ModelParameterCell_MR
            !
            procedure, public :: avgModel1D => avgModel1D_ModelParameterCell_MR
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
        self%param_grid => cell_cond%grid
        !
        self%anisotropic_level = anisotropic_level
        !
        allocate( self%cell_cond( anisotropic_level ) )
        !
        self%cell_cond(1) = cell_cond
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
        self%param_grid => cell_cond(1)%grid
        !
        self%anisotropic_level = size( cell_cond )
        !
        self%cell_cond = cell_cond
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
        call self%deallCell
        !
    end subroutine ModelParameterCell_MR_dtor
    !
    !> perhaps this, and similar, can now go into Base ModelParameter_Cell class?
    !
    subroutine PDEmapping_ModelParameterCell_MR( self, e_vec )
        implicit none
        !
        class( ModelParameterCell_MR_t ), intent( in ) :: self
        class( Vector_t ), intent( inout ) :: e_vec
        !
        type( rScalar3D_MR_t ) :: sigma_cell_al_mr
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "PDEmapping_ModelParameterCell_MR > self not allocated" )
        endif
        !
        if( .NOT. e_vec%is_allocated ) then
            call errStop( "PDEmapping_ModelParameterCell_MR > e_vec not allocated" )
        endif
        !
        !> cell cond as MR with AirLayers
        sigma_cell_al_mr = rScalar3D_MR_t( self%metric%grid, CELL )
        !
        call self%modelToCell( self%air_cond, sigma_cell_al_mr )
        !
        call e_vec%zeros
        !
        call self%cellToEdge( sigma_cell_al_mr, e_vec )
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
        type( rScalar3D_MR_t ) :: sigma_cell_al_mr
        !
        !> cell cond as MR with AirLayers
        sigma_cell_al_mr = rScalar3D_MR_t( self%metric%grid, CELL )
        !
        ! difference with PDEmapping: set air layers to zero, pass optional argument dsigma
        call self%modelToCell( R_ZERO, sigma_cell_al_mr, dsigma )
        !
        call e_vec%zeros
        !
        call self%cellToEdge( sigma_cell_al_mr, e_vec )
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
        class( Scalar_t ), allocatable :: sigma_cell_al
        type( rScalar3D_MR_t ) :: sigma_cell_al_mr
        integer :: i
        !
        if( .NOT. e_vec%is_allocated ) then
            call errStop( "dPDEmapping_T_ModelParameterCell_MR > e_vec not allocated" )
        endif
        !
        if( .NOT. dsigma%is_allocated ) then
            call errStop( "dPDEmapping_T_ModelParameterCell_MR > dsigma not allocated" )
        endif
        !
        !>
        call self%metric%createScalar( real_t, CELL, sigma_cell_al )    !sigma_cell_al SHOULD BE GENERIC FOR edgeToCell...
        !
        call self%edgeToCell( e_vec, sigma_cell_al )                    !... HERE!
        !
        sigma_cell_al_mr = rScalar3D_MR_t( sigma_cell_al%grid, sigma_cell_al%grid_type )
        !
        select type( sigma_cell_al )
            !
            class is( rScalar3D_MR_t )
                !
                do i = 1, size( sigma_cell_al%sub_scalar )
                    !
                    sigma_cell_al_mr%sub_scalar(i)%v = real( sigma_cell_al%sub_scalar(i)%v, kind=prec )
                    !
                enddo
                !
            class default
                call errStop( "dPDEmapping_T_ModelParameterCell_MR > Unclassified sigma_cell_al" )
            !
        end select
        !
        deallocate( sigma_cell_al )
        !
        ! difference with PDEmapping: set air layers to zero, pass optional argument dsigma
        call self%cellToModel( sigma_cell_al_mr, dsigma )
        !
    end subroutine dPDEmapping_T_ModelParameterCell_MR
    !
    !> No subroutine briefing
    !
    subroutine nodeCond_ModelParameterCell_MR( self, sigma_node )
        implicit none
        !
        class( ModelParameterCell_MR_t ), intent( in ) :: self
        class( Scalar_t ), intent( inout ) :: sigma_node
        !
        type( rScalar3D_MR_t ) :: sigma_cell_al_mr
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "nodeCond_ModelParameterCell_MR > self not allocated" )
        endif
        !
        if( .NOT. sigma_node%is_allocated ) then
            call errStop( "nodeCond_ModelParameterCell_MR > sigma_node not allocated" )
        endif
        !
        sigma_cell_al_mr = rScalar3D_MR_t( self%metric%grid, CELL )
        !
        call self%modelToCell( self%air_cond, sigma_cell_al_mr )
        !
        call sigma_node%zeros
        !
        call self%cellToNode( sigma_cell_al_mr, sigma_node )
        !
    end subroutine nodeCond_ModelParameterCell_MR
    !
    !> Idea: no matter what the grid is (SG or MR) a routine like this
    !> uses conductivity array(s) in self (always SG at present -- but might change this???)
    !> and returns the "full array" including air layers on the appropriate (SG or MR) grid
    !> this is where SG and MR differ - so might be able to make more routines generic
    !> (for example define key generic routines in base ModelParameter_Cell class)
    !> 
    !> COMPLICATION so that this can also be used for both PDE mapping and dPDE mapping:
    !> optional argument dsigma also a model parameter defined on SG grid needs to be converted
    !> to MR and multiplied by sigma (from self)
    !
    subroutine modelToCell_ModelParameterCell_MR( self, air_value, sigma_cell, dsigma )
        implicit none
        !
        class( ModelParameterCell_MR_t ), intent( in ) :: self
        real( kind=prec ), intent ( in ) :: air_value
        class( Scalar_t ), intent( inout )  :: sigma_cell
        class( ModelParameter_t ), intent( in ), optional :: dsigma
        !
        character(:), allocatable :: job
        integer :: i_grid, nz_air
        logical :: dPDE
        type( Grid3D_MR_t ) :: temp_grid_mr
        type( rScalar3D_MR_t ) :: sigma_cell_mr, dsigma_cell_mr, temp_sigma_cell_mr
        !
        dPDE = present( dsigma )
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
                call errStop( "modelToCell_ModelParameterCell_MR > Grid must be MR!" )
            !
        end select
        !
        !> cell cond as MR without AirLayers
        sigma_cell_mr = rScalar3D_MR_t( temp_grid_mr, CELL )
        !
        call sigma_cell_mr%fromSG( self%cell_cond(1) )
        !
        if ( dPDE ) then
            !
            !> convert dsigma to MR w/o airlayers
            job = DERIV
            !
            dsigma_cell_mr = rScalar3D_MR_t( temp_grid_mr, CELL )
            !
            call dsigma_cell_mr%fromSG( dsigma%getCond(1) )
            !
        else
            job = FORWARD
        endif
        !
        !> Temporary MR Scalar allowing access its specific properties
        temp_sigma_cell_mr = sigma_cell
        !
        !> Considering AirLayers just for the first sub_grid
        !    THIS NEEDS TO BE GENERALIZED
        nz_air = self%metric%grid%NzAir
        !
        temp_sigma_cell_mr%sub_scalar(1)%v( :, :, 1:nz_air ) = air_value 
        !
        temp_sigma_cell_mr%sub_scalar(1)%v( :, :, nz_air+1:sigma_cell_mr%sub_scalar(1)%nz ) = &
        self%sigMap( real( sigma_cell_mr%sub_scalar(1)%v, kind=prec ), job )
        !
        !> sigMapping for the next sub-grids
        do i_grid = 2, size( sigma_cell_mr%sub_scalar )
            !
            temp_sigma_cell_mr%sub_scalar( i_grid )%v = self%sigMap( real( sigma_cell_mr%sub_scalar( i_grid )%v, kind=prec ) )
            !
        enddo
        !
        if( dPDE ) then
            !
            !> multiply temp_sigma_cell_mr by dsigma_cell_mr
            temp_sigma_cell_mr%sub_scalar(1)%v( :, :, nz_air+1:sigma_cell_mr%sub_scalar(1)%nz ) = &
            temp_sigma_cell_mr%sub_scalar(1)%v( :, :, nz_air+1:sigma_cell_mr%sub_scalar(1)%nz ) * &
            dsigma_cell_mr%sub_scalar(1)%v
            !
            do i_grid = 2, size( sigma_cell_mr%sub_scalar )
                !
                temp_sigma_cell_mr%sub_scalar( i_grid )%v = temp_sigma_cell_mr%sub_scalar( i_grid )%v * &
                dsigma_cell_mr%sub_scalar(i_grid)%v
                !
            enddo
            !
        endif
        !
        sigma_cell = temp_sigma_cell_mr
        !
    end subroutine modelToCell_ModelParameterCell_MR
    !
    !> This is the adjoint of modelToCell_ModelParameterCell_MR -- only used for dPDEmapping
    !
    subroutine cellToModel_ModelParameterCell_MR( self, sigma_cell, dsigma )
        implicit none
        !
        class( ModelParameterCell_MR_t ), intent( in ) :: self
        class( Scalar_t ), intent( in ) :: sigma_cell
        class( ModelParameter_t ), intent( inout ) :: dsigma
        !
        integer :: i_grid, nz_air
        type( Grid3D_MR_t ) :: temp_grid_mr
        type( rScalar3D_MR_t ) :: temp_sigma_cell_mr, sigma_cell_mr, dsigma_cell_mr
        type( rScalar3D_SG_t ) :: dsigma_cell_sg
        !
        if( .NOT. sigma_cell%is_allocated ) then
            call errStop( "cellToModel_ModelParameterCell_MR > sigma_cell not allocated" )
        endif
        !
        if( .NOT. dsigma%is_allocated ) then
            call errStop( "cellToModel_ModelParameterCell_MR > dsigma not allocated" )
        endif
        !
        !  set MR grid -- need this as a preliminary step for fwd and adjt
        select type( grid => self%metric%grid )
            !
            class is( Grid3D_MR_t )
                !
                temp_grid_mr = Grid3D_MR_t( self%param_grid%nx, self%param_grid%ny, &
                self%param_grid%nzAir, self%param_grid%nzEarth, self%param_grid%dx, &
                self%param_grid%dy, self%param_grid%dz, grid%cs )
                !
            class default
                !
                call errStop( "cellToModel_ModelParameterCell_MR > Grid must be MR!" )
                !
        end select
        !
        !> cell cond as MR without AirLayers  -- will need this to compute sigma on MR grid cells
        !   again same preliminaries for fwd and adjt
        sigma_cell_mr = rScalar3D_MR_t( temp_grid_mr, CELL )
        !
        call sigma_cell_mr%fromSG( self%cell_cond(1) )
        !
        !> sigMapping for all sub-grids (applied to background model parameter)
        do i_grid = 1, size( sigma_cell_mr%sub_scalar )
            !
            sigma_cell_mr%sub_scalar( i_grid )%v =   & 
            self%sigMap( real( sigma_cell_mr%sub_scalar( i_grid )%v, kind=prec ), DERIV )
            !
        enddo
        !
        !> Temporary MR Scalar allowing access its specific properties
        temp_sigma_cell_mr = sigma_cell
        !
        !> Now the actual adjoint mapping -- reverse order from fwd
        !
        !> multiply temp_sigma_cell_mr by dsigma_cell_mr, considering AirLayers just for the first sub_grid2
        !> This needs to be generalized !!!!
        nz_air = self%metric%grid%NzAir
        !
        sigma_cell_mr%sub_scalar(1)%v = sigma_cell_mr%sub_scalar(1)%v * &
        temp_sigma_cell_mr%sub_scalar(1)%v( :,:,nz_air+1:temp_sigma_cell_mr%sub_scalar(1)%nz )
        !
        do i_grid = 2, size( sigma_cell_mr%sub_scalar )
            !
            sigma_cell_mr%sub_scalar( i_grid )%v =  sigma_cell_mr%sub_scalar( i_grid )%v * &
            temp_sigma_cell_mr%sub_scalar( i_grid )%v
            !
        enddo
        !
        call sigma_cell_mr%toSG( dsigma_cell_sg )
        !
        call dsigma%setCond( dsigma_cell_sg, 1 )
        !
    end subroutine cellToModel_ModelParameterCell_MR
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
end Module ModelParameterCell_MR
!
