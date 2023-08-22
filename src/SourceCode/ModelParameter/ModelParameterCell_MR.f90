!
!> Basic model parameter class --
!> cell conductivities defined on numerical grid
!> in Earth cells only -- as coded this is specific
!> to usual MR Grid.
!
module ModelParameterCell_MR
    !
    use ModelParameterCell_SG
    !
    type, extends( ModelParameterCell_t ) :: ModelParameterCell_MR_t
        !
        !> No derived properties
        !> cell_cond treated as SG field array here
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
    function ModelParameterCell_MR_ctor_one_cond( grid, cell_cond, anisotropic_level, param_type ) result( self )
        implicit none
        !
        class( Grid_t ), intent( in ) :: grid
        type( rScalar3D_SG_t ), intent( in ) :: cell_cond
        integer, intent( in ) :: anisotropic_level
        character(:), allocatable, optional, intent( in ) :: param_type
        !
        type( ModelParameterCell_MR_t ) :: self
        !
        integer :: nzAir
        !
        !write( *, * ) "Constructor ModelParameterCell_MR_ctor_one_cond"
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
        self%cell_cond(1) = cell_cond
        !
        self%cell_cond(1)%store_state = compound
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
    function ModelParameterCell_MR_ctor_all_conds( grid, cell_cond, param_type ) result( self )
        implicit none
        !
        class( Grid_t ), intent( in ) :: grid
        type( rScalar3D_SG_t ), dimension(:), intent( in ) :: cell_cond
        character(:), allocatable, optional, intent( in ) :: param_type
        !
        type( ModelParameterCell_MR_t ) :: self
        !
        integer :: i, nzAir
        !
        !write( *, * ) "Constructor ModelParameterCell_MR_ctor_all_conds"
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
        if( allocated( self%cell_cond ) ) deallocate( self%cell_cond )
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
        class( Scalar_t ), allocatable, intent( inout ) :: sigma_node
        !
        call errStop( "nodeCond_ModelParameterCell_MR > must be implemented!" )
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
        call errStop( "PDEmapping_ModelParameterCell_MR > under implementation!" )
        ! !
        ! type( Grid3D_SG_t ) :: temp_grid_sg, sub_grid
        ! type( ModelParameterCell_SG_t ) :: sub_model
        ! class( Vector_t ), allocatable :: temp_vector_sg
        ! type( rScalar3D_SG_t ) :: temp_scalar
        ! type( rVector3D_MR_t ) :: temp_e_vec
        ! integer :: kk1, kk2, k1, k2, s1, s2, Cs, Nz, kDown, kUp, i, j, k 
        ! real( kind=prec ), allocatable, dimension(:,:,:) :: temp_v
        ! !
        ! !>
        ! k1 = self%metric%grid%nzAir + 1
        ! k2 = self%metric%grid%nz
        ! !
        ! call self%metric%createScalar( real_t, temp_scalar, CELL )
        ! !
        ! !> USE SIGMAP !!!!
        ! select case( self%param_type )
            ! !
            ! case( LOGE )
                ! !
                ! temp_scalar%v = exp( self%air_cond )
                ! temp_scalar%v(:,:,k1:k2) = exp( self%cell_cond(1)%s%getV() )
                ! !
            ! case( LINEAR )
                ! !
                ! temp_scalar%v = self%air_cond
                ! temp_scalar%v(:, :, k1:k2) = self%cell_cond(1)%s%getV()
                ! !
            ! case default
                ! !
                ! call errStop( "PDEmapping_ModelParameterCell_MR > Unknow param_grid!" )
                ! !
        ! end select
        ! !
        ! select type( grid => self%metric%grid )
            ! !
            ! class is( Grid3D_MR_t )
                ! !
                ! k1 = 1
                ! k2 = 0
                ! !
                ! do k = 1, grid%n_grids
                    ! !
                    ! ! First extract relevant part of cellcond
                    ! !(including cells on both sides of any fine interface)
                    ! sub_grid = grid%sub_grids(k)
                    ! k1 = k2 + 1
                    ! k2 = k2 + sub_grid%nz
                    ! kUp = k1
                    ! kk1 = 1
                    ! !
                    ! if( k > 1 ) then
                        ! !
                        ! if( grid%Coarseness( k-1, 1 ) < grid%Coarseness( k, 1 ) ) then
                            ! !
                            ! ! switch polarity of if test, relative to MR case
                            ! kUp = kUp - 1
                            ! kk1 = 2
                            ! !
                        ! endif
                        ! !
                    ! endif
                    ! !
                    ! kk2 = kk1 + sub_grid%nz
                    ! kDown = k2
                    ! !
                    ! if( k < grid%n_grids ) then
                        ! !
                        ! if( grid%Coarseness( k, 1 ) > grid%Coarseness( k+1, 1 ) ) then
                            ! !
                            ! ! switch polarity of if test, relative to MR case
                            ! kDown = kDown + 1
                            ! !
                        ! end if
                        ! !
                    ! end if
                    ! !
                    ! Nz = kDown - kUp + 1
                    ! !
                    ! allocate( temp_v( sub_grid%nx, sub_grid%ny, Nz ) )
                    ! temp_v = R_ZERO
                    ! !
                    ! Cs = 2 ** grid%Coarseness(k, 1)
                    ! !
                    ! s1 = size( temp_scalar%v, 1 )
                    ! s2 = size( temp_scalar%v, 2 )
                    ! !
                    ! do i = 1, Cs
                        ! !
                        ! do j = 1, Cs
                            ! temp_v = temp_v + temp_scalar%v( i:s1:Cs, j:s2:Cs, kUp:kDown )
                        ! enddo
                        ! !
                    ! enddo
                    ! !
                    ! !> Create grid object for extended sub-grid
                    ! temp_grid_sg = Grid3D_SG_t( size(sub_grid%dx), size(sub_grid%dy), 0, &
                    ! size( self%metric%grid%dz( kUp:kDown ) ) )
                    ! !
                    ! temp_grid_sg%dx = sub_grid%dx
                    ! temp_grid_sg%dy = sub_grid%dy
                    ! temp_grid_sg%dz = self%metric%grid%dz( kUp:kDown )
                    ! !
                    ! call temp_grid_sg%setup
                    ! !
                    ! !> then scalar field object containing averaged
                    ! !> conductivities.
                    ! sub_model = ModelParameterCell_SG_t( temp_grid_sg, self%cell_cond(1)%s, LINEAR )
                    ! !
                    ! call sub_model%cell_cond(1)%s%setV( temp_v / (Cs*Cs) )
                    ! !
                    ! !> Finally average onto edges of extended subgrid
                    ! allocate( temp_vector_sg, source = rVector3D_SG_t( temp_grid_sg, EDGE ) )
                    ! !
                    ! call sub_model%PDEmapping( temp_vector_sg )
                    ! !
                    ! !> ... and then extract appropriate part into
                    ! !> appropriate subgrid of MR vector object.
                    ! !> (NOTE: redundant edges will not be set(correctly at least)
                    ! !> ... but these are not used in the equations!
                    ! temp_e_vec = e_vec
                    ! !
                    ! temp_e_vec%sub_vectors(k)%x = temp_vector_sg%x(:,:, kk1:kk2)
                    ! temp_e_vec%sub_vectors(k)%y = temp_vector_sg%y(:,:, kk1:kk2)
                    ! temp_e_vec%sub_vectors(k)%z = temp_vector_sg%z(:,:, kk1:kk2-1)
                    ! !
                    ! e_vec = temp_e_vec
                    ! !
                ! enddo
                ! !
                ! !mr_m = e_vec%get_array() !????
                ! !
            ! class default
                ! call errStop( "sgToMR_rScalar3D_MR > Unclassified grid" )
            ! !
        ! end select
        ! !
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
        call errStop( "slice1D_ModelParameterCell_MR > must be implemented!" )
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
