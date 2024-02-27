!
!> Basic model parameter class --
!> cell conductivities defined on numerical grid
!> in Earth cells only -- as coded this is specific
!> to usual SG Grid.
!
module ModelParameterCell_SG
    !
    use ModelParameterCell
    !
    type, extends( ModelParameterCell_t ) :: ModelParameterCell_SG_t
        !
        !> No derived properties
        !
        contains
            !
            final :: ModelParameterCell_SG_dtor
            !
            !> Mappings
            procedure, public :: nodeCond => nodeCond_ModelParameterCell_SG
            !
            procedure, public :: PDEmapping => PDEmapping_ModelParameterCell_SG
            procedure, public :: dPDEmapping => dPDEmapping_ModelParameterCell_SG
            procedure, public :: dPDEmapping_T => dPDEmapping_T_ModelParameterCell_SG
            !
            !> Dimensioned operations
            procedure, public :: slice1D => slice1D_ModelParameterCell_SG
            procedure, public :: slice2D => slice2D_ModelParameterCell_SG
            !
            procedure, public :: avgModel1D => avgModel1D_ModelParameterCell_SG
            !
            procedure, public :: write => write_ModelParameterCell_SG
            !
    end type ModelParameterCell_SG_t
    !
    interface ModelParameterCell_SG_t
         module procedure ModelParameterCell_SG_ctor_one_cond
         module procedure ModelParameterCell_SG_ctor_all_conds
    end interface ModelParameterCell_SG_t
    !
contains
    !
    !> No subroutine briefing
    !
    function ModelParameterCell_SG_ctor_one_cond( cell_cond, anisotropic_level, param_type ) result( self )
        implicit none
        !
        type( rScalar3D_SG_t ), intent( in ) :: cell_cond
        integer, intent( in ) :: anisotropic_level
        character(:), allocatable, optional, intent( in ) :: param_type
        !
        type( ModelParameterCell_SG_t ) :: self
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
    end function ModelParameterCell_SG_ctor_one_cond
    !
    !> No subroutine briefing
    !
    function ModelParameterCell_SG_ctor_all_conds( cell_cond, param_type ) result( self )
        implicit none
        !
        type( rScalar3D_SG_t ), dimension(:), intent( in ) :: cell_cond
        character(:), allocatable, optional, intent( in ) :: param_type
        !
        type( ModelParameterCell_SG_t ) :: self
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
    end function ModelParameterCell_SG_ctor_all_conds
    !
    !> No subroutine briefing
    !
    subroutine ModelParameterCell_SG_dtor( self )
        implicit none
        !
        type( ModelParameterCell_SG_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor ModelParameterCell_SG"
        !
        call self%baseDealloc
        !
        call self%deallocCell
        !
    end subroutine ModelParameterCell_SG_dtor
    !
    !> Map cell_cond to nodes
    !
    subroutine nodeCond_ModelParameterCell_SG( self, sigma_node )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( in ) :: self
        class( Scalar_t ), intent( inout ) :: sigma_node
        !
        integer :: k0, k1, k2
        type( rScalar3D_SG_t ) :: sigma_cell
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "nodeCond_ModelParameterCell_SG > self not allocated" )
        endif
        !
        if( .NOT. sigma_node%is_allocated ) then
            call errStop( "nodeCond_ModelParameterCell_SG > sigma_node not allocated" )
        endif
        !
        k0 = self%metric%grid%NzAir
        k1 = k0 + 1
        k2 = self%metric%grid%Nz
        !
        sigma_cell = rScalar3D_SG_t( sigma_node%grid, CELL )
        !
        sigma_cell%v( :, :, 1:k0 ) = self%air_cond
        !
        sigma_cell%v( :, :, k1:k2 ) = self%sigMap( real( self%cell_cond(1)%v, kind=prec ) )
        !
        call sigma_cell%mult( self%metric%v_cell )
        !
        call sigma_cell%sumToNode( sigma_node, .TRUE. )
        !
        !> Later fix for SP2 - 27/02/2024!!!!
        call sigma_node%div( self%metric%v_node )
        !
        call sigma_node%mult( cmplx( 0.125_prec, 0.0, kind=prec ) )
        !
    end subroutine nodeCond_ModelParameterCell_SG
    !
    !> Map the entire model cells into a single edge Vector_t (e_vec).
    !
    subroutine PDEmapping_ModelParameterCell_SG( self, e_vec )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( in ) :: self
        class( Vector_t ), intent( inout ) :: e_vec
        !
        type( rScalar3D_SG_t ), allocatable, dimension(:) :: sigma_cells
        integer :: i, k0, k1, k2
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "PDEmapping_ModelParameterCell_SG > self not allocated" )
        endif
        !
        if( .NOT. e_vec%is_allocated ) then
            call errStop( "PDEmapping_ModelParameterCell_SG > e_vec not allocated" )
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
            sigma_cells(i) = rScalar3D_SG_t( self%metric%grid, CELL )
            !
            sigma_cells(i)%v( :, :, 1:k0 ) = self%air_cond
            !
            sigma_cells(i)%v( :, :, k1:k2 ) = self%sigMap( real( self%cell_cond(i)%v, kind=prec ) )
            !
            call sigma_cells(i)%mult( self%metric%v_cell )
            !
        enddo
        !
        !> Call due sumCells based on anisotropic_level
        if( self%anisotropic_level == 1 ) then
            !
            call e_vec%sumCells( sigma_cells(1) )
            !
        elseif( self%anisotropic_level == 2 ) then
            !
            call e_vec%sumCells( sigma_cells(1), sigma_cells(2) )
            !
        else
            call errStop( "PDEmapping_ModelParameterCell_SG > unsupported anisotropy level" )
        endif
        !
        call e_vec%mult( cmplx( 0.25_prec, 0.0, kind=prec ) )
        !
        call e_vec%div( self%metric%v_edge )
        !
    end subroutine PDEmapping_ModelParameterCell_SG
    !
    !> Map the perturbation between two models onto a single Vector_t (e_vec).
    !
    subroutine dPDEmapping_ModelParameterCell_SG( self, dsigma, e_vec )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( in ) :: self
        class( ModelParameter_t ), intent( in ) :: dsigma
        class( Vector_t ), allocatable, intent( inout ) :: e_vec
        !
        type( rScalar3D_SG_t ), allocatable, dimension(:) :: sigma_cells
        type( rScalar3D_SG_t ) :: dsigma_cond
        integer :: i, k0, k1, k2
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "dPDEmapping_ModelParameterCell_SG > self not allocated" )
        endif
        !
        if( .NOT. dsigma%is_allocated ) then
            call errStop( "dPDEmapping_ModelParameterCell_SG > dsigma not allocated" )
        endif
        !
        if( .NOT. e_vec%is_allocated ) then
            call errStop( "dPDEmapping_ModelParameterCell_SG > e_vec not allocated" )
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
            !> Create and initialize sigma_cells with zeros
            sigma_cells(i) = rScalar3D_SG_t( self%metric%grid, CELL )
            !
            sigma_cells(i)%v( :, :, k1:k2 ) = self%sigMap( real( self%cell_cond(i)%v, kind=prec ), DERIV )
            !
            dsigma_cond = dsigma%getCond(i)
            !
            sigma_cells(i)%v( :, :, k1:k2 ) = sigma_cells(i)%v( :, :, k1:k2 ) * dsigma_cond%v
            !
            call sigma_cells(i)%mult( self%metric%v_cell )
            !
        enddo
        !
        !> Call specific sumCells based on anisotropic_level
        if( self%anisotropic_level == 1 ) then
            !
            call e_vec%sumCells( sigma_cells(1) )
            !
        elseif( self%anisotropic_level == 2 ) then
            !
            call e_vec%sumCells( sigma_cells(1), sigma_cells(2) )
            !
        else
            call errStop( "dPDEmapping_ModelParameterCell_SG > unsupported anisotropy level" )
        endif
        !
        call e_vec%mult( cmplx( 0.25_prec, 0.0, kind=prec ) )
        !
        call e_vec%div( self%metric%v_edge )
        !
    end subroutine dPDEmapping_ModelParameterCell_SG
    !
    !> Transpose the perturbation represented in a Vector_t (e_vec), to a new dsigma model.
    !
    subroutine dPDEmapping_T_ModelParameterCell_SG( self, e_vec, dsigma )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( in ) :: self
        class( Vector_t ), intent( in ) :: e_vec
        class( ModelParameter_t ), allocatable, intent( inout ) :: dsigma
        !
        class( Vector_t ), allocatable :: e_vec_interior
        type( rScalar3D_SG_t ) :: dsigma_cond, sigma_cell
        type( GenScalar_t ), allocatable, dimension(:) :: sigma_cells
        complex( kind=prec ), allocatable, dimension(:,:,:) :: sigma_cell_v
        integer :: i, k0, k1, k2
        !
        if( .NOT. self%is_allocated ) then
            call errStop( "dPDEmapping_T_ModelParameterCell_SG > self not allocated" )
        endif
        !
        if( .NOT. e_vec%is_allocated ) then
            call errStop( "dPDEmapping_T_ModelParameterCell_SG > e_vec not allocated" )
        endif
        !
        if( .NOT. dsigma%is_allocated ) then
            call errStop( "dPDEmapping_T_ModelParameterCell_SG > dsigma not allocated" )
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
            call errStop( "dPDEmapping_T_ModelParameterCell_SG > unsupported anisotropy level" )
        endif
        !
        deallocate( e_vec_interior )
        !
        k0 = self%metric%grid%NzAir
        k1 = k0 + 1
        k2 = self%metric%grid%Nz
        !
        !allocate( dsigma, source = self )
        !dsigma = self
        !
        do i = 1, self%anisotropic_level
            !
            dsigma_cond = dsigma%getCond(i)
            !
            call dsigma_cond%zeros
            !
            dsigma_cond%v = self%sigMap( real( self%cell_cond(i)%v, kind=prec ), DERIV )
            !
            call sigma_cells(i)%s%mult( self%metric%v_cell )
            !
            sigma_cell = sigma_cells(i)%s
            !
            sigma_cell_v = sigma_cell%v
            !
            dsigma_cond%v = dsigma_cond%v * sigma_cell_v( :, :, k1:k2 )
            !
            call dsigma%setCond( dsigma_cond, i )
            !
        enddo
        !
        deallocate( sigma_cells )
        !
    end subroutine dPDEmapping_T_ModelParameterCell_SG
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
        real( kind=prec ), allocatable, dimension(:) :: cond_slice
        real( kind=prec ) :: wt, temp_sigma_value
        integer :: i, j, k
        !
        model_param_1D = ModelParameter1D_t( self%metric%grid%slice1D() )
        !
        allocate( cond_slice( self%metric%grid%nzEarth ) )
        !
        do k = 1, self%metric%grid%nzEarth
            !
            wt = R_ZERO
            temp_sigma_value = R_ZERO
            do i = 1, self%metric%grid%Nx
                do j = 1, self%metric%grid%Ny
                    !
                    wt = wt + self%metric%grid%dx(i) * self%metric%grid%dy(j)
                    temp_sigma_value = temp_sigma_value + self%cell_cond(1)%v( i, j, k ) * &
                    self%metric%grid%dx(i) * self%metric%grid%dy(j)
                    !
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
        real( kind=prec ), allocatable, dimension(:,:) :: cond_slice
        !
        param_type = LINEAR
        !
        m2D = ModelParameter2D_t( self%metric%grid%slice2D() )
        !
        allocate( cond_slice( self%metric%grid%ny, self%metric%grid%nzEarth ) )
        !
        if( axis == 1 ) then
            cond_slice = self%sigMap( real( self%cell_cond(1)%v(j,:,:), kind=prec ) )
        elseif( axis == 2 ) then
            cond_slice = self%sigMap( real( self%cell_cond(1)%v(:,j,:), kind=prec ) )
        elseif( axis == 3 ) then
            cond_slice = self%sigMap( real( self%cell_cond(1)%v(:,:,j), kind=prec ) )
        else
            call errStop( "slice2D_ModelParameterCell_SG > wrong axis" )
        endif
        !
        call m2D%setConductivity( cond_slice, self%air_cond, param_type, self%mKey )
        !
        deallocate( cond_slice )
        !
    end function slice2D_ModelParameterCell_SG
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
        real( kind=prec ), allocatable, dimension(:,:,:) :: cond_v
        integer :: Nx, Ny, NzEarth, ii, i, j, k, ios
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
            do ii = 1, self%anisotropic_level
                !
                !> Convert (horizontal) conductivity to resistivity
                !
                cond_v = self%cell_cond(ii)%v
                !
                if( index( self%param_type, "LOGE" ) > 0 .OR. index( self%param_type, "LOG10" ) > 0 ) then
                    cond_v = -cond_v
                elseif( index(self%param_type, "LINEAR" ) > 0 ) then
                    cond_v = ONE / cond_v
                endif
                !
                !> Write the (horizontal) resistivity
                !
                write( ioModelParam, * )
                !
                do k = 1, nzEarth
                    do j = 1, Ny
                        do i = Nx, 1, -1
                            write( ioModelParam, "(es13.5)", iostat = ios, advance = "no" ) cond_v(i,j,k)
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
            call errStop( "write_ModelParameterCell_SG > Error opening file ["//file_name//"]!" )
        endif
        !
    end subroutine write_ModelParameterCell_SG
    !
end Module ModelParameterCell_SG
!