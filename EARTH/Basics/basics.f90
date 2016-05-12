! *****************************************************************************
module basics
	! This module will contain any constants used throughout in other scripts.
	! To ensure precision, and to follow the original script, standard math
	! constants are declared as variables rather than parameters for now. 

  use math_constants
  implicit none

! pi        - 4.0d0*datan(1.0d0) or 3.1415926536
! d2r       - degrees to radians conversion 180.d0/pi
! mu0       - the permeability of freespace 4.0d-07*pi
! rair	    - the resistivity of the air layer 1.0d10
! Rearth    - the equatorial Earth's radius in km 6378.164
 
! real(8)					:: pi,d2r,r2d,mu0,rair,Rearth

  public					:: init_const


Contains

  ! ***************************************************************************
  ! * init_const has to be called before any other subroutines using PI etc

	  subroutine init_const(pi,d2r,r2d,mu0,rair,Rearth)

	    real(8), intent(out)  :: pi,d2r,r2d,mu0,rair,Rearth

		pi=4.d0*datan(1.d0)
        d2r=pi/180.d0
        r2d=180.d0/pi
		mu0 = 4.0d-07 * pi
		rair = 1.0d10
		Rearth = 6371.0d0 !6378.164

	  end subroutine init_const	! init_const



end module basics !	basics