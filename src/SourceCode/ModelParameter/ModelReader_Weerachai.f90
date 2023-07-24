!
!> Derived class to define a ModelReader_Weerachai
!
module ModelReader_Weerachai
    !
    use Constants
    use String
    use Grid
    use Grid3D_SG
    use Grid3D_MR
    use rScalar3D_SG
    use ModelParameterCell_SG
    use ModelReader
    use ForwardControlFile
    !
    type, extends( ModelReader_t ), public :: ModelReader_Weerachai_t
        !
        !> No derived properties
        !
     contains
        !
        procedure, public :: read => readModelReaderWeerachai
        !
    end type ModelReader_Weerachai_t
    !
contains
    !
    !> No subroutine briefing
    subroutine readModelReaderWeerachai( self, file_name, grid, model ) 
        implicit none
        !
        class( ModelReader_Weerachai_t ), intent( in ) :: self
        character(*), intent( in ) :: file_name
        class( Grid_t ), allocatable, intent( out ) :: grid
        class( ModelParameter_t ), allocatable, intent( out ) :: model
        !
        character( len=80 ) :: someChar 
        character(:), allocatable :: paramType 
        integer :: nx, ny, nzEarth, nzAir, someIndex, i, ii, j, k, ioPrm, io_stat, p_nargs, anisotropic_level
        real( kind=prec ), dimension(:), allocatable :: dx, dy, dz
        real( kind=prec ) :: ox, oy, oz, rotDeg
        real( kind=prec ), dimension(:, :, :), allocatable :: rho
        type( rScalar3D_SG_t ) :: ccond
        real( kind=prec ) :: ALPHA
        character(len=200), dimension(20) :: args
        !
        integer, allocatable, dimension(:) :: layers, levels
        !
        layers = (/ 0, 3, 1, 4, 2, 4 /)
        !
        someChar = ""
        paramType = ""
        someIndex = 0
        ALPHA = 3.0
        !
        open( newunit = ioPrm, file = trim(file_name),status = "old", iostat = io_stat)
        !
        if( io_stat == 0 ) then
            !
            !> First read the comment line
            read(ioPrm, "(a80)") someChar
            !
            !> Now read the second line with the grid dimensions
            read( ioPrm, "(a80)" ) someChar
            read( someChar, * ) nx, ny, nzEarth, someIndex
            !
            !> Now read the second line with the grid dimensions
            nzAir = 0
            !
            allocate(dx(nx))
            allocate(dy(ny))
            allocate(dz(nzAir + nzEarth))
            !
            read(ioPrm, *) (dx(j), j = 1, nx)
            read(ioPrm, *) (dy(j), j = 1, ny)
            read(ioPrm, *) (dz(j), j = nzAir + 1, nzAir + nzEarth)
            !
            if( someIndex /= 0 ) then
                call errStop( "readModelReaderWeerachai > Mapping not supported." )
            endif
            !
            !> By default assume "LINEAR RHO" -
            !> Weerachai"s linear resistivity format
            if( index( someChar, "LOGE" ) > 0 ) then
                paramType = LOGE
            elseif( index( someChar, "LOG10" ) > 0) then
                paramType = LOG_10
            else
                paramType = LINEAR
            endif
            !
            !> The default method for creating air layers in the grid has been deleted
            !
            select case( grid_format )
                !
                case( GRID_SG )
                    allocate( grid, source = Grid3D_SG_t( nx, ny, nzAir, nzEarth, dx, dy, dz ) )
                case( GRID_MR )
                    allocate( grid, source = Grid3D_MR_t( nx, ny, nzAir, nzEarth, dx, dy, dz, layers ) )
                case default
                    call errStop( "readModelReaderWeerachai > Unknow grid_format ["//grid_format//"]." )
            end select
            !
            !> Consider isotope at first
            anisotropic_level = 1
            !
            !> Read conductivity values in a model parameter object.
            if( index( someChar, "VTI" ) > 0 ) then
                !
                anisotropic_level = 2
                !
            end if
            !
            !>
            do ii = 1, anisotropic_level
                !
                allocate( rho( nx, ny, nzEarth ) )
                !
                do k = 1, nzEarth
                    do j = 1, ny
                        read(ioPrm, *, iostat = io_stat) (rho(i, j, k), i = nx, 1, -1)
                    enddo
                enddo
                !
                ccond = rScalar3D_SG_t( grid, CELL_EARTH )
                !
                if( index( paramType, "LOGE" ) > 0 .OR. &
                    index( paramType, "LOG10" ) > 0 ) then
                    !
                    ccond%v = -rho
                    !
                elseif( index(paramType, "LINEAR") > 0 ) then
                    !
                    ccond%v = ONE/rho
                    !
                endif
                !
                deallocate( rho )
                !
                if( anisotropic_level == 1 ) then
                    !
                    allocate( model, source = ModelParameterCell_SG_t( grid, ccond, 1, paramType ) )
                    !
                else
                    !
                    if( allocated( model ) ) then
                        !
                        call model%setCond( ccond, ii )
                        !
                    else
                        !
                        allocate( model, source = ModelParameterCell_SG_t( grid, ccond, 2, paramType ) )
                        !
                    endif
                    !
                endif
                !
            end do
            !
            !> ALWAYS convert modelParam to natural log for computations ????
            paramType = LOGE
            call model%setType( paramType )
            !
            !> End - Reading cells conductivity values.
            !
            !> In case the grid origin is stored next (in metres!)...
            read(ioPrm, *, iostat = io_stat) ox, oy, oz
            !
            !> Defaults to the grid center at the Earth"s surface
            if(io_stat /= 0) then
                 ox = -sum(dx)/2.0
                 oy = -sum(dy)/2.0
                 oz = R_ZERO
            endif
            !
            call grid%setOrigin( ox, oy, oz )
            !
            read(ioPrm, *, iostat = io_stat) rotDeg
            if(io_stat /= 0) then
                 rotDeg = R_ZERO
            endif
            !
            call grid%setRotation( rotDeg )
            !
            close( ioPrm )
            !
        else
            call errStop( "Error opening ["//file_name//"] in readModelReaderWeerachai" )
        endif
        !
    end subroutine readModelReaderWeerachai
    !
end module ModelReader_Weerachai
