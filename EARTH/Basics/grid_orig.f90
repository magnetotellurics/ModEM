! *****************************************************************************
module grid_orig
	! This module contains the definition of the model grid only.

  implicit none

  integer               						:: nx, ny, nz, nzEarth, nzAir
	  ! storing the (spherical) grid in Randie Mackie's format
	  ! when allocated, dimensions will be x(nx), y(ny+1), z(nz+1)
  real(8), allocatable, dimension(:)			:: x,y,z

end module grid_orig