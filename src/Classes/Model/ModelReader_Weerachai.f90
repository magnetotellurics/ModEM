module ModelReader_Weerachai
  use Constants
  use Grid
  use Grid3D_SG
  use rScalar3D_SG
  use ModelParameter
  use ModelReader
  use ModelParameter
  use ModelParameterCell_SG

  implicit none
  
  private  
  
  type, extends(ModelReader_t), public :: ModelReader_Weerachai_t
     private
     
   contains
     private
     
     procedure, public :: Read
     procedure, public :: Write
     
  end type ModelReader_Weerachai_t
  
contains

  subroutine Read(self, fileName, grid, model) 
    ! Arguments
    class(ModelReader_Weerachai_t), intent(in)  :: self
    character(*)                  , intent(in)  :: fileName
    class(Grid_t)          , pointer, intent(out) :: grid
    class(ModelParameter_t), pointer, intent(out) :: model
    ! Local variables     
    character(80) :: someChar = '', paramType = ''
    integer :: nx, ny, nzEarth, nzAir
    integer :: someIndex = 0, i, j, k
    integer :: ioPrm, istat
    real(kind = prec), dimension(:), allocatable :: dx, dy, dz
    real(kind = prec) :: ox, oy, oz, rotDeg
    real(kind = prec), dimension(:, :, :), allocatable :: rho
    type(rScalar3D_SG_t) :: ccond
    real(kind = prec) :: ALPHA = 3.0
    
    open(newunit = ioPrm, file = trim(fileName), &
         status = 'old', iostat = istat)
    
    if (istat /= 0) then
       write(0, *) 'ERROR:WeerachaiSG_GridReader_t:'//&
            &      'Could not open input file.'       
       grid => null()
       return
    end if
    
    ! First read the comment line
    read(ioPrm, '(a80)') someChar

    ! Now read the second line with the grid dimensions
    read(ioPrm, '(a80)') someChar
    read(someChar, *) nx, ny, nzEarth, someIndex
    
    if (someIndex /= 0) then
       write(0, *) 'ERROR:WeerachaiSG_GridReader_t:'
       write(*, *) '  Mapping not supported.'
       grid => null()
       STOP
    end if
    
    ! By default assume 'LINEAR RHO' -
    ! Weerachai's linear resistivity format
    if (index(someChar, 'LOGE') > 0) then
       paramType = 'LOGE'
    else if (index(someChar, 'LOG10') > 0) then
       paramType = 'LOG10'
    else
       paramType = 'LINEAR'
    end if

    ! No information about the air layers in file. Hardcoded here.
    ! For backwards compatibility, keeping this logic for now...
    ! but it may be overwritten by forward solver configuration file.
    nzAir = 10
    
    allocate(dx(nx))
    allocate(dy(ny))
    allocate(dz(nzAir + nzEarth))
    
    read(ioPrm, *) (dx(j), j = 1, nx)
    read(ioPrm, *) (dy(j), j = 1, ny)
    read(ioPrm, *) (dz(j), j = nzAir + 1, nzAir + nzEarth)

    !**
    ! ** AirLayers setup **
    !
    ! Following is Kush's approach to setting air layers:
    ! mirror imaging the dz values in the air layer with respect to
    ! earth layer as far as we can using the following formulation
    ! air layer(bottom:top) = (alpha)^(j-1) * earth layer(top:bottom).
    i = nzAir + 1
    j = 0
    do k = nzAir, 1, -1
       j = j + 1
       dz(k) = ((ALPHA)**(j-1))*dz(i)
       i = i + 1
    end do

    ! The topmost air layer has to be at least 30 km.
    if (dz(1).lt.30000) then
       dz(1) = 30000
    end if
    !
    ! End - Setting air layers.
    !*
    
    !**
    ! Create the grid object
    allocate(grid, source = Grid3D_SG_t(nx, ny, nzAir, nzEarth, dx, dy, dz))

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
       ccond = rScalar3D_SG_t(grid, CELL_EARTH)       
       
       if ((index(paramType, 'LOGE') > 0).or.&
            (index(paramType, 'LOG10') > 0)) then
          ccond%v = -rho
       else if (index(paramType, 'LINEAR') > 0) then
          ccond%v = ONE/rho
       end if

       allocate(model, source = ModelParameterCell_SG_t(grid, ccond, paramType))
    end select
    
    ! ALWAYS convert modelParam to natural log for computations
    paramType = 'LOGE'
    call model%SetType(paramType)
    !
    ! End - Reading cells conductivity values.
    
    !*
    ! In case the grid origin is stored next (in metres!)...
    read(ioPrm, *, iostat = istat) ox, oy, oz
    
    ! Defaults to the grid centre at the Earth's surface
    if (istat /= 0) then
       ox = -sum(dx)/2.0
       oy = -sum(dy)/2.0
       oz = 0.0
    end if

    call grid%SetOrigin(ox, oy, oz)
    read(ioPrm, *, iostat = istat) rotDeg
    if (istat /= 0) then
       rotDeg = 0.0
    end if
    call grid%SetGridRotation(rotDeg)
    
    !**
    ! Clean up
    !*
    close(ioPrm)
    
  end subroutine Read

  subroutine Write(self, fileName, grid, model)
    ! Arguments
    class(ModelReader_Weerachai_t), intent(in) :: self
    character(*)                  , intent(in) :: fileName
    class(Grid_t)                 , intent(in) :: grid
    class(ModelParameter_t)       , intent(in) :: model
  end subroutine Write
  
end module ModelReader_Weerachai
