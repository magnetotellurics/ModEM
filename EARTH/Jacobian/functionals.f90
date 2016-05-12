! *****************************************************************************
module functionals
  ! Module containing the functions to compute the data functionals from the
  ! magnetic field on the surface

  use math_constants
  implicit none
    
  ! We would like to be able to compute responses not only
  ! at an observatory, but at any point on the Earth's surface
  ! (other radii are unlikely to be required right now).
  ! Try to keep the arguments consistent, as we would like to
  ! be able to pass the procedures as arguments to subroutines.
  ! Currently using them in a more straight-forward way.

	! NB: Converted C and D responses into D and C response ratios
	! by getting rid of the tan(theta) and sin(theta) terms, resp.
	! This is a necessity since C responses do not make sense around
	! the equator and ruin the interpolation & other operations.
	! From now on, everything is done in terms of the ratios.
	! Any C and D response data is converted into ratios at the
	! input and all the outputs now only include the fields and
	! the ratios. - Anna Kelbert, 18 Jan 2007 
  interface C_ratio
	 module procedure C_ratio
  end interface

  interface D_ratio
	 module procedure D_ratio
  end interface

Contains

	! EM response C = \dfrac{R\tan(theta)}{2} \times Hz/Hy
	! (theta is colatitude in radians)
    complex(8) function C_ratio(Hx,Hy,Hz,theta)
	  
	  real(8), intent(in)			:: theta
	  complex(8), intent(in)		:: Hz,Hy,Hx

	   C_ratio = km2m * (EARTH_R/2.0) * Hz/Hy
	   !C = km2m * (EARTH_R/2.0) * dtan(theta) * Hz/Hy

    end function C_ratio

	! EM response D = \dfrac{R\cos(theta)}{2} \times Hx/Hy
	! (theta is colatitude in radians)
	! THIS EXPRESSION WAS USED IN THE ORIGINAL FORWARD SOLVER
	! NEED TO MAKE SURE THIS IS THE SAME FORMULA AS THAT USED
	! BY IKUKO FUJII - THERE IS AN INCONSISTENCY IN THE PAPER
	! AS TO COS vs SIN and the sign of the expression
    complex(8) function D_ratio(Hx,Hy,Hz,theta)

	  real(8), intent(in)			:: theta
	  complex(8), intent(in)		:: Hx,Hy,Hz

	   D_ratio = km2m * (EARTH_R/2.0) * Hx/Hy
	   !D = km2m * (EARTH_R/2.0) * dsin(theta) * Hx/Hy
	   !D = km2m * (EARTH_R/2.0) * dcos(theta) * Hx/Hy

    end function D_ratio


end module functionals
