module ModelReader_Weerachai
    !
    use Constants
    use Grid
    use Grid3D_SG
    use rScalar3D_SG
    use ModelParameter
    use ModelReader
    use ModelParameterCell_SG
    !
    type, extends( ModelReader_t ), public :: ModelReader_Weerachai_t
        !
        ! PROPERTIES HERE
        !
     contains
        !
        procedure, public :: Read
        procedure, public :: Write
        !
    end type ModelReader_Weerachai_t
    !
contains
    !
    subroutine Read( self, fileName, grid, model ) 
        implicit none
        !
        class( ModelReader_Weerachai_t ), intent( in )        :: self
        character(*), intent( in )                            :: fileName
        class( Grid_t ), allocatable, intent( out )           :: grid
        class( ModelParameter_t ), allocatable, intent( out ) :: model
        !
        character( len = 80 ) :: someChar 
        character(:), allocatable :: paramType 
        integer :: nx, ny, nzEarth, nzAir, someIndex, i, j, k, ioPrm, istat
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
        open(newunit = ioPrm, file = trim(fileName),status = "old", iostat = istat)
        !
        if (istat /= 0) then
             write(0, *) "ERROR:WeerachaiSG_GridReader_t:"//&
                        &            "Could not open input file."
             !grid => null()
             return
        end if
        !
        ! First read the comment line
        read(ioPrm, "(a80)") someChar
        !
        ! Now read the second line with the grid dimensions
        read(ioPrm, "(a80)") someChar
        read(someChar, *) nx, ny, nzEarth, someIndex
        !
        ! Now read the second line with the grid dimensions
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
        if (someIndex /= 0) then
            write(*, *) "ERROR:WeerachaiSG_GridReader_t:"
            write(*, *) "    Mapping not supported."
            !grid => null()
            STOP
        end if
        !
        ! By default assume "LINEAR RHO" -
        ! Weerachai"s linear resistivity format
        if (index(someChar, "LOGE") > 0) then
             paramType = LOGE
        else if (index(someChar, "LOG10") > 0) then
             paramType = LOG_10
        else
             paramType = LINEAR
        end if
        !
        ! The default method for creating air layers in the grid has been deleted
        !**
        ! Create the grid object with nzAir = 0 -- no air layers so far
        grid = Grid3D_SG_t( nx, ny, nzAir, nzEarth, dx, dy, dz )
        !**
        ! Read conductivity values in a model parameter object.
        allocate(rho(nx, ny, nzEarth))
        do k = 1, nzEarth
            do j = 1, ny
                read(ioPrm, *, iostat = istat) &
                        (rho(i, j, k), i = nx, 1, -1)
            end do
        end do
        !**
        ! Convert from resistivity to conductivity
        !*        
        select type(grid)
            class is(Grid3D_SG_t)
                !
                ccond = rScalar3D_SG_t( grid, CELL_EARTH ) 
                !
                if ((index(paramType, "LOGE") > 0).or.&
                (index(paramType, "LOG10") > 0)) then
                    ccond%v = -rho
                else if (index(paramType, "LINEAR") > 0) then
                    ccond%v = ONE/rho
                end if
                !
                allocate( model, source = ModelParameterCell_SG_t( grid, ccond, paramType ) )
        end select
        !
        ! ALWAYS convert modelParam to natural log for computations ????
        paramType = LOGE
        call model%SetType( paramType )
        !
        ! End - Reading cells conductivity values.
        
        !*
        ! In case the grid origin is stored next (in metres!)...
        read(ioPrm, *, iostat = istat) ox, oy, oz
        !
        ! Defaults to the grid center at the Earth"s surface
        if (istat /= 0) then
             ox = -sum(dx)/2.0
             oy = -sum(dy)/2.0
             oz = 0.0
        end if
        !
        call grid%SetOrigin(ox, oy, oz)
        !
        read(ioPrm, *, iostat = istat) rotDeg
        if (istat /= 0) then
             rotDeg = 0.0
        end if
        call grid%SetGridRotation(rotDeg)
        !**
        ! Clean up
        !*
        close(ioPrm)
        !
    end subroutine Read
    !
    subroutine Write( self, fileName, grid, model )
        implicit none
        !
        class( ModelReader_Weerachai_t ), intent( in ) :: self
        character(*), intent( in )                     :: fileName
        class( Grid_t ), intent( in )                  :: grid
        class( ModelParameter_t ), intent( in )        :: model
    end subroutine Write
    
end module ModelReader_Weerachai
