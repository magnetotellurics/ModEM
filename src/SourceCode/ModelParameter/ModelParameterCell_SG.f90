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
        type( rScalar3D_SG_t ), allocatable, dimension(:) :: cell_cond
        !
        contains
            !
            final :: ModelParameterCell_SG_dtor
            !
            procedure, public :: getOneCond => getOneCond_ModelParameterCell_SG
            procedure, public :: getAllCond => getAllCond_ModelParameterCell_SG
            !
            procedure, public :: setOneCond => setOneCond_ModelParameterCell_SG
            procedure, public :: setAllCond => setAllCond_ModelParameterCell_SG
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
         module procedure ModelParameterCell_SG_ctor_one_cond
         module procedure ModelParameterCell_SG_ctor_all_conds
    end interface ModelParameterCell_SG_t
    !
contains
    !
    !> No subroutine briefing
    !
    function ModelParameterCell_SG_ctor_one_cond( grid, cell_cond, anisotropic_level, param_type ) result( self )
        implicit none
        !
        class( Grid_t ), intent( in ) :: grid
        type( rScalar3D_SG_t ), intent( in ) :: cell_cond
        integer, intent( in ) :: anisotropic_level
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
            call self%setsigMap( param_type )
        !
        endif
        !
        self%is_allocated = .TRUE.
        !
    end function ModelParameterCell_SG_ctor_one_cond
    !
    !> No subroutine briefing
    !
    function ModelParameterCell_SG_ctor_all_conds( grid, cell_cond, param_type ) result( self )
        implicit none
        !
        class( Grid_t ), intent( in ) :: grid
        type( rScalar3D_SG_t ), dimension(:), intent( in ) :: cell_cond
        character(:), allocatable, optional, intent( in ) :: param_type
        !
        type( ModelParameterCell_SG_t ) :: self
        !
        integer :: i, nzAir
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
                ( grid%nz - grid%nzAir ), grid%dx, grid%dy, grid%dz( grid%nzAir+1:grid%nz ) ) )
        !
        self%anisotropic_level = size( cell_cond )
        !
        self%cell_cond = cell_cond
        !
        do i = 1, self%anisotropic_level
            self%cell_cond(i)%store_state = compound
        enddo
        !
        if( present( param_type ) ) then
            !
            call self%setsigMap( param_type )
        !
        endif
        !
        self%is_allocated = .TRUE.
        !
    end function ModelParameterCell_SG_ctor_all_conds
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
        !> CRASHING ????
        if( allocated( self%cell_cond ) ) deallocate( self%cell_cond )
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
                    wt = wt + self%metric%grid%dx(i) * self%metric%grid%dy(j)
                    temp_sigma_value = temp_sigma_value + self%cell_cond(1)%v( i, j, k ) * &
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
            stop "Error: slice2D_ModelParameterCell_SG > wrong axis"
        endif
        !
        call m2D%setConductivity( cond_slice, self%air_cond, param_type, self%mKey )
        !
        deallocate( cond_slice )
        !
    end function slice2D_ModelParameterCell_SG
    !
    !> No function briefing
    !
    function getOneCond_ModelParameterCell_SG( self, i_cond ) result( cell_cond )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( in ) :: self
        integer, intent( in ) :: i_cond
        !
        type( rScalar3D_SG_t ) :: cell_cond
        !
        if( i_cond .GT. self%anisotropic_level ) then
            !
            call errStop( "getOneCond_ModelParameterCell_SG > conductivity level too high" )
            !
        endif
        !
        cell_cond = self%cell_cond( i_cond )
        !
    end function getOneCond_ModelParameterCell_SG
    !
    !> No function briefing
    !
    function getAllCond_ModelParameterCell_SG( self ) result( cell_cond )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( in ) :: self
        !
        type( rScalar3D_SG_t ), allocatable, dimension(:) :: cell_cond
        !
        cell_cond = self%cell_cond
        !
    end function getAllCond_ModelParameterCell_SG
    !
    !> No interface subroutine briefing
    !
    subroutine setOneCond_ModelParameterCell_SG( self, cell_cond, i_cond )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( inout ) :: self
        class( Scalar_t ), intent( in ) :: cell_cond
        integer, intent( in ) :: i_cond
        !
        if( i_cond .LE. self%anisotropic_level ) then
            !
            self%cell_cond( i_cond ) = cell_cond
            !
        else
            !
            call errStop( "setOneCond_ModelParameterCell_SG > Unsupport general anisotropy yet" )
            !
        endif
        !
    end subroutine setOneCond_ModelParameterCell_SG
    !
    !> No interface subroutine briefing
    !
    subroutine setAllCond_ModelParameterCell_SG( self, cell_cond )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( inout ) :: self
        class( Scalar_t ), dimension(:), intent( in ) :: cell_cond
        !
        integer :: i
        !
        do i = 1, self%anisotropic_level
            !
            if( .NOT. cell_cond(i)%is_allocated ) then
                !
                call errStop( "setAllCond_ModelParameterCell_SG > cell_cond has no V" )
                !
            else if( .NOT. self%cell_cond(i)%is_allocated )then
                !
                call errStop( "setAllCond_ModelParameterCell_SG > self has no V" )
                !
            endif
            !
            self%cell_cond(i) = cell_cond(i)
            !
        enddo
        !
    end subroutine setAllCond_ModelParameterCell_SG
    !
    !> No subroutine briefing
    subroutine zeros_ModelParameterCell_SG( self )
        implicit none
        !
        class( ModelParameterCell_SG_t ), intent( inout ) :: self
        !
        integer :: i
        !
        do i = 1, self%anisotropic_level
            !
            call self%cell_cond(i)%zeros
            !
        enddo
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
        if( .NOT. rhs%is_allocated ) then
            call errStop( "copyFrom_ModelParameterCell_SG > rhs not allocated." )
        endif
        !
        select type( rhs )
            !
            class is( ModelParameterCell_SG_t )
                !
                write(*,*) "##### ENTER COPY #####"
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
        integer :: i, counter, nx, ny, nz, nzAir, nz_earth
        !
        counter = 0
        !
        do i = 1, self%anisotropic_level
            !
            if( .NOT. self%cell_cond(i)%is_allocated ) then
                write( *, * ) "Error: countModel_ModelParameterCell_SG > cell_cond (", i, ") not allocated!"
                stop
            endif
            !
            call self%cell_cond(i)%grid%getDimensions( nx, ny, nz, nzAir )
            nz_earth = nz - nzAir
            !
            counter = counter + self%cell_cond(i)%Nx * self%cell_cond(i)%Ny * nz_earth
            !
        enddo
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
        integer :: i
        complex( kind=prec ), allocatable :: v(:, :, :)
        !
        select type( rhs )
            !
            class is( ModelParameterCell_SG_t )
                !
                do i = 1, self%anisotropic_level
                    !
                    if( self%cell_cond(i)%isCompatible( rhs%cell_cond(i) ) ) then
                        !
                        v = a1 * self%cell_cond(i)%v + a2 * rhs%cell_cond(i)%v
                        !
                        call self%cell_cond(i)%setV( v )
                        !
                    else
                        write( *, * ) "Error: linComb_ModelParameterCell_SG > Incompatible rhs cell_cond (", i, ")!"
                        stop
                    endif
                    !
                enddo
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
        !
        real( kind=prec ) :: rvalue
        !
        integer :: i
        !
        rvalue = R_ZERO
        !
        select type( rhs )
            !
            class is( ModelParameterCell_SG_t )
                !
                do i = 1, self%anisotropic_level
                    !
                    if( self%cell_cond(i)%isCompatible( rhs%cell_cond(i) ) ) then
                        !
                        rvalue = rvalue + sum( self%cell_cond(i)%v * rhs%cell_cond(i)%v )
                        !
                    else
                        write( *, * ) "Error: dotProd_ModelParameterCell_SG > Incompatible rhs cell_cond (", i, ")!"
                        stop
                    endif
                    !
                enddo
                !
            class default
                call errStop( "dotProd_ModelParameterCell_SG > Unclassified rhs" )
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
        type( rScalar3D_SG_t ), allocatable, dimension(:) :: sigma_cell
        integer :: i, k0, k1, k2
        !
        if( .NOT. allocated( eVec ) ) then
            allocate( eVec, source = rVector3D_SG_t( self%metric%grid, EDGE ) )
        else
            eVec = rVector3D_SG_t( self%metric%grid, EDGE )
        endif
        !
        allocate( sigma_cell( self%anisotropic_level ) )
        !
        k0 = self%metric%grid%nzAir
        k1 = k0 + 1
        k2 = self%metric%grid%Nz
        !
        do i = 1, self%anisotropic_level
            !
            sigma_cell(i) = rScalar3D_SG_t( self%metric%grid, CELL )
            !
            sigma_cell(i)%v( :, :, 1:k0 ) = self%air_cond
            !
            sigma_cell(i)%v( :, :, k1:k2 ) = self%sigMap( self%cell_cond(i)%v )
            !
            call sigma_cell(i)%mult( self%metric%VCell )
            !
        enddo
        !
        !> Call due avgCells based on anisotropic_level
        if( self%anisotropic_level == 1 ) then
            !
            call eVec%avgCells( sigma_cell(1) )
            !
        elseif( self%anisotropic_level == 2 ) then
            !
            call eVec%avgCells( sigma_cell(1), sigma_cell(2) )
            !
        else
            call errStop( "PDEmapping_ModelParameterCell_SG > unsupported anisotropy level" )
        endif
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
        type( rScalar3D_SG_t ), allocatable, dimension(:) :: sigma_cell
        character( len=5 ), parameter :: JOB = "DERIV"
        integer :: i, k0, k1, k2
        !
        if( .NOT. allocated( eVec ) ) then
            allocate( eVec, source = rVector3D_SG_t( self%metric%grid, EDGE ) )
        else
            eVec = rVector3D_SG_t( self%metric%grid, EDGE )
        endif
        !
        allocate( sigma_cell( self%anisotropic_level ) )
        !
        call eVec%zeros
        !
        do i = 1, self%anisotropic_level
            !
            sigma_cell(i) = rScalar3D_SG_t( self%metric%grid, CELL )
            !
            k0 = self%metric%grid%NzAir
            k1 = k0 + 1
            k2 = self%metric%grid%Nz
            !
            call sigma_cell(i)%zeros
            !
            sigma_cell(i)%v( :, :, k1:k2 ) = self%sigMap( self%cell_cond(i)%v, JOB )
            !
            select type( dsigma )
                !
                class is( ModelParameterCell_SG_t )
                    !
                    sigma_cell(i)%v(:,:,k1:k2) = sigma_cell(i)%v(:,:,k1:k2) * dsigma%cell_cond(i)%v
                    !
                class default
                    call errStop( "dPDEmapping_ModelParameterCell_SG > Unclassified dsigma" )
                !
            end select
            !
            call sigma_cell(i)%mult( self%metric%Vcell )
            !
        enddo
        !
        !> Call specific avgCells based on anisotropic_level
        if( self%anisotropic_level == 1 ) then
            !
            call eVec%avgCells( sigma_cell(1) )
            !
        elseif( self%anisotropic_level == 2 ) then
            !
            call eVec%avgCells( sigma_cell(1), sigma_cell(2) )
            !
        else
            call errStop( "dPDEmapping_ModelParameterCell_SG > unsupported anisotropy level" )
        endif
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
        class( Vector_t ), allocatable :: evec_interior
        type( rScalar3D_SG_t ), allocatable, dimension(:) :: sigma_cell
        complex( kind=prec ), allocatable :: sigma_v(:, :, :)
        character( len=5 ), parameter :: JOB = "DERIV"
        integer :: i, k0, k1, k2
        !
        allocate( sigma_cell( self%anisotropic_level ) )
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
        do i = 1, self%anisotropic_level
			sigma_cell(i) = rScalar3D_SG_t( self%metric%grid, CELL )
        enddo
        !
        !> Call specific sumEdges based on anisotropic_level
        if( self%anisotropic_level == 1 ) then
			!
			call evec_interior%sumEdges( sigma_cell(1), .TRUE. )
			!
        elseif( self%anisotropic_level == 2 ) then
			!
			call evec_interior%sumEdges( sigma_cell(1), sigma_cell(2), .TRUE. )
			!
        else
		!
        call errStop( "dPDEmapping_T_ModelParameterCell_SG > unsupported anisotropy level" )
        endif
        !
        deallocate( evec_interior )
        !
        allocate( dsigma, source = ModelParameterCell_SG_t( self%param_grid, self%cell_cond, self%param_type ) )
        !
        select type( dsigma )
            !
            class is( ModelParameterCell_SG_t )
                !
                do i = 1, self%anisotropic_level
                    !
                    dsigma%cell_cond(i)%v = self%sigMap( self%cell_cond(i)%v, JOB )
                    !
                    call sigma_cell(i)%mult( self%metric%Vcell )
                    !
                    dsigma%cell_cond(i)%v = sigma_cell(i)%v(:,:,k1:k2) * dsigma%cell_cond(i)%v
                    !
                enddo
                !
            class default
                stop "Error: dPDEmapping_T_ModelParameterCell_SG > Unclassified dsigma"
            !
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
        integer :: i
        !
        if( .NOT. self%is_allocated ) then
            stop "Error: setType_ModelParameterCell_SG > Self not allocated."
        endif
        !
        do i = 1, self%anisotropic_level
            !
            if( trim( param_type ) .EQ. trim( self%param_type ) ) then
                ! Nothing to be done
            elseif( self%param_type == "" ) then
                self%param_type = trim( param_type )
            elseif( self%param_type == LINEAR ) then
                !
                if( param_type == LOGE ) then
                    !
                    self%cell_cond(i)%v = log( self%cell_cond(i)%v )
                    !
                elseif( param_type == LOG_10) then
                    !
                    self%cell_cond(i)%v = log10( self%cell_cond(i)%v )
                    !
                endif
                !
            elseif( param_type == LINEAR ) then
                !
                if( self%param_type == LOGE ) then
                    !
                    self%cell_cond(i)%v = exp( self%cell_cond(i)%v )
                    !
                elseif( self%param_type == LOG_10 ) then
                    !
                    self%cell_cond(i)%v = exp( self%cell_cond(i)%v * log(10.) )
                    !
                endif
                !
            elseif( ( self%param_type == LOGE ) .AND. ( param_type == LOG_10 ) ) then
                !
                self%cell_cond(i)%v = self%cell_cond(i)%v / log(10.)
                !
            elseif( ( self%param_type == LOG_10 ) .AND. ( param_type == LOGE ) ) then
                !
                self%cell_cond(i)%v = self%cell_cond(i)%v * log(10.)
                !
            else
                stop "Error: setType_ModelParameterCell_SG > Unknown param_type."
            endif
            !
        enddo
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
        class( ModelParameterCell_SG_t ), intent( inout ) :: self
        !
        integer :: i
        !
        write( *, * ) "ModelParameterCell_SG_t:", self%mKey, self%air_cond, self%param_type, &
        self%is_allocated, self%param_grid%nx, self%param_grid%ny, self%param_grid%nz, self%param_grid%nzAir
        !
        do i = 1, self%anisotropic_level
            !
            call self%cell_cond(i)%print
            !
        enddo
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
        integer :: Nx, Ny, NzEarth, ii, i, j, k, ios
        !
        ! Verbose
        !write( *, * ) "     > Write Model to file: [", file_name, "]"
        !
        !> Convert modelParam to natural log or log10 for output
        !paramType = userParamType
        !
        !if( self%%is ) then
        !    call getValue_modelParam(m, paramType, self%cell_cond(1), v_v=ccond_v)
        !else
        !    call getValue_modelParam(m, paramType, self%cell_cond(1))
        !endif
        !
        open( unit = ioModelParam, file = file_name, action = "write", form = "formatted", iostat = ios )
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
            do ii = 1, self%anisotropic_level
                !
                !> Convert (horizontal) conductivity to resistivity
                rho_h = self%cell_cond(ii)
                if( index( self%param_type,"LOGE" ) > 0 .OR. index( self%param_type,"LOG10" ) > 0 ) then
                    rho_h%v = -self%cell_cond(ii)%v
                elseif( index(self%param_type,"LINEAR" ) > 0) then
                    rho_h%v = ONE/self%cell_cond(ii)%v
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
