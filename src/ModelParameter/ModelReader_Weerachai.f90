!
!> Derived class to define a ModelReader_Weerachai
!
module ModelReader_Weerachai
    !
    use Constants
    use Grid
    use Grid3D_SG
    use rScalar3D_SG
    use ModelParameterCell_SG
    use ModelReader
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
        integer :: nx, ny, nzEarth, nzAir, someIndex, i, j, k, ioPrm, io_stat
        real( kind=prec ), dimension(:), allocatable :: dx, dy, dz
        real( kind=prec ) :: ox, oy, oz, rotDeg
        real( kind=prec ), dimension(:, :, :), allocatable :: rho
        type( rScalar3D_SG_t  ) :: ccond
        real( kind=prec ) :: ALPHA
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
            read(ioPrm, "(a80)") someChar
            read(someChar, *) nx, ny, nzEarth, someIndex
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
            if(someIndex /= 0) then
                stop "Error: readModelReaderWeerachai > Mapping not supported."
            endif
            !
            !> By default assume "LINEAR RHO" -
            !> Weerachai"s linear resistivity format
            if(index(someChar, "LOGE") > 0) then
                 paramType = LOGE
            elseif(index(someChar, "LOG10") > 0) then
                 paramType = LOG_10
            else
                 paramType = LINEAR
            endif
            !
            !> The default method for creating air layers in the grid has been deleted
            !
            !> Create the grid object with nzAir = 0 -- no air layers so far
            allocate( grid, source = Grid3D_SG_t( nx, ny, nzAir, nzEarth, dx, dy, dz ) )
            !
            !> Read conductivity values in a model parameter object.
            allocate( rho( nx, ny, nzEarth ) )
            do k = 1, nzEarth
                do j = 1, ny
                    read(ioPrm, *, iostat = io_stat) (rho(i, j, k), i = nx, 1, -1)
                enddo
            enddo
            !
            !> Convert from resistivity to conductivity
            select type( grid )
                !
                class is( Grid3D_SG_t )
                    !
                    ccond = rScalar3D_SG_t( grid, CELL_EARTH ) 
                    !
                    if((index(paramType, "LOGE") > 0) .OR. &
                    (index(paramType, "LOG10") > 0)) then
                        ccond%v = -rho
                    elseif(index(paramType, "LINEAR") > 0) then
                        ccond%v = ONE/rho
                    endif
                    !
                    deallocate( rho )
                    !
                    allocate( model, source = ModelParameterCell_SG_t( grid, ccond, paramType ) )
                    !
                class default
                    stop "Error: readModelReaderWeerachai > Unclassified grid"
                !
            end select
            !
            !> ALWAYS convert modelParam to natural log for computations ????
            paramType = LOGE
            call model%SetType( paramType )
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
            call grid%SetOrigin( ox, oy, oz )
            !
            read(ioPrm, *, iostat = io_stat) rotDeg
            if(io_stat /= 0) then
                 rotDeg = R_ZERO
            endif
            !
            call grid%SetGridRotation( rotDeg )
            !
            close( ioPrm )
            !
        else
            write( *, * ) "Error opening [", file_name, "] in readModelReaderWeerachai"
            stop
        endif
        !
    end subroutine readModelReaderWeerachai
    !
end module ModelReader_Weerachai
