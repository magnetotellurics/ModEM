! *****************************************************************************
module modeldef
	! This module contains the type definitions for the structures that contain
	! the full information about the model parametrization

  use sg_scalar
  use griddef
  implicit none


  ! ***************************************************************************
  ! * type modelPoint_t contains information about a point in the Earth: coords
  ! * and relative location on the grid
  type :: modelPoint_t

	  real(8)								  :: phi,theta,r

  end type modelPoint_t


  ! ***************************************************************************
  ! * type modelShell_t contains the thin shell S-conductance on the uppermost
  ! * few km of the solid Earth - mainly Earth crust. A typical depth is 12650 m,
  ! * the values correspond to the given grid in GM coordinates so that
  ! * rho(i,j) = 1/(cond(i,j)/depth) = depth/cond(i,j) (where depth in in metres)
  type :: modelShell_t

	  real(8), dimension(:,:), pointer		  :: cond !(nx,ny) - conductance
	  logical                                 :: allocated

  end type modelShell_t	! modelShell_t


  ! ***************************************************************************
  ! * type modelFunc_t contains the information about a single function that is
  ! * multiplied by the respective coefficient in the model parametrization;
  ! * add more information entries if necessary; for this, we may want to use
  ! * an analogue of "union" structure available in C, once Fortran has it.
  type :: modelFunc_t

	  ! order number in functionals structure
	  integer								  :: num

	  ! which function is used?
	  character(80)							  :: name

	  ! degree and order, if F is a spherical harmonic term
	  integer								  :: l,m

	  ! range of values, if F is an identity function
	  ! real(8)								  :: xmin,xmax,ymin,ymax

  end type modelFunc_t	! modelFunc_t


  ! ***************************************************************************
  ! * type modelCoeff_t contains the information about a single coefficient of the
  ! * model parametrization
  type :: modelCoeff_t

     type (modelLayer_t), pointer			  :: L	! which layer
     type (modelFunc_t) , pointer			  :: F	! which function it corresponds to
     type (modelCoeff_t), pointer			  :: prev	! same coefficient prev layer
     type (modelCoeff_t), pointer			  :: next	! same coefficient next layer
	 character(80)				  :: code ! information to identify the coefficient
     real(8)					  :: value	! value of coeff.
     real(8)                      :: min,max ! range information for inversion
	 ! if coefficient is fixed/frozen, no derivative will be provided in the output
     logical					  :: frozen
	 ! it is possible to create a padding zero-valued coefficient that does not exist
	 logical					  :: exists

  end type modelCoeff_t	! modelCoeff_t


  ! ***************************************************************************
  ! * type modelLayer_t contains the information for a single layer in the layer
  ! * structure of the model parametrization, that is defined in modeldef.f90;
  ! * going from the first (uppermost) layer to the last that ends at CMB
  ! * Air layer is layer 0.
  ! * If used for the source parametrization, layers correspond to periods, so
  ! * add that variable to the layer type explicitly in the hope of avoiding confusion
  type :: modelLayer_t

	  ! order number in layer structure
	  integer											  :: num

	  ! upper layer boundary and the width of the layer in kilometers
	  real(8)											  :: lbound,ubound
	  real(8)											  :: width,depth

	  ! layer regularization parameters
	  real(8)											  :: alpha,beta

	  ! the indicator of whether the linear combination of F_i gives log_10(rho)
	  logical											  :: if_log

	  ! the filename of additional rho info if it exists; if so, add these
	  ! values to either rho or log_10(rho) (as specified by if_log) on the grid
	  character(80)										  :: fn=''

      ! period (in days) instead of the layer depth needed for the source structure
      real(8)                                             :: period

  end type modelLayer_t	! modelLayer_t



  ! ***************************************************************************
  ! * type modelParam_t combines together all information about the parametrization
  ! * in the most general case. When mapped onto grid, this information provides
  ! * a resistivity structure. It is also used to compute the gradients.
  ! * Since the allocatable attribute inside a derived data type definition
  ! * is not universally supported, we use pointers instead.

  type :: modelParam_t

      ! defines the type of the model parameter ('harmonic'/'grid' etc)
      character(80)                                       :: type=''

	  ! all information about the layer structure
	  integer											  :: nF
	  type (modelFunc_t),  dimension(:), pointer		  :: F

	  ! all information about the layer structure
	  integer											  :: nL
	  type (modelLayer_t), dimension(:), pointer		  :: L

	  ! all coefficient and function information
	  integer											  :: nc ! number of active coeffs
	  type (modelCoeff_t), dimension(:,:), pointer		  :: c	! structure of all coeffs

	  ! the near-surface shell is also part of model parameter
	  type (modelShell_t)                                 :: crust

      ! the full resistivity on a grid which may be defined instead of the harmonics
      type (rscalar)                                      :: rho

	  ! also contains a pointer to the grid, used for mappings
	  ! NOTE: this may need rethinking unless one global grid is used!
	  type (grid_t), pointer                              :: grid

	  logical											  :: allocated=.FALSE.

      !  this logical is set to true by zero_modelParam ONLY
      logical                                             :: zeroValued = .FALSE.

      ! this is needed for output: true only after CmSqrt has been applied to the model
      logical                                             :: smoothed=.FALSE.

	  ! temporary:  .true. for function outputs only; necessary to avoid memory leaks
	  ! (probably will not be needed in the future when compilers will support
	  ! ISO/IEC 15581 - the "allocatable array extension")
	  logical											  :: temporary=.FALSE.

  end type modelParam_t	! modelParam_t


  ! ***************************************************************************
  ! * type misfitDef_t contains information that specifies how to calculate the
  ! * misfit; normally this is user-defined.
  type :: misfitDef_t

	character(80)					:: name	  ! type of data functional
	!real(8)							:: alpha  ! horizontal smoothing
	!real(8)							:: beta	  ! vertical smoothing
	real(8)							:: mu  ! damping parameter

  end type misfitDef_t


  ! ***************************************************************************
  ! * type misfitInfo_t contains information about misfit and derivative
  ! * for a single frequency and functional type; normally we would create
  ! * an array type (misfitInfo_t) misfit(nfreq,nfunc). In future this will
  ! * replace the current way of storing misfit information.
  type :: misfitInfo_t

	real(8)							:: w  ! weight
	real(8)							:: v  ! non-averaged misfit value
	integer							:: n  ! number of data values (to average)
	type (modelParam_t)				:: d  ! derivative

  end type misfitInfo_t


  ! ***************************************************************************
  ! * type intArray_t contains an array of integer indices

  type :: intArray_t

	  integer, dimension(:), pointer					  :: num

  end type intArray_t


  ! ***************************************************************************
  ! * type realArray_t contains an array of real(8) indices

  type :: realArray_t

	  real(8), dimension(:), pointer					  :: value

  end type realArray_t


  ! ***************************************************************************
  ! * type grid_t contains the necessary information about user-defined grid
  ! This is now part of SG_Basics/stag_grid.f90
!  type :: grid_t
!
!	  ! storing the grid dimensions
!	  integer               					:: nx, ny, nz, nzEarth, nzAir
!	  ! storing the (spherical) grid in Randie Mackie's format
!	  ! when allocated, dimensions will be x(nx), y(ny+1), z(nz+1)
!	  real(8), pointer, dimension(:)			:: x,y,z
!	  ! storing the cell node coordinates, in radians (ph,th) and km (r)
!	  ! dimensions are ph(nx+1), th(ny+1), r(nz+1)
!	  real(8), pointer, dimension(:)  			:: ph,th,r
!
!  end type grid_t	! grid_t



  ! ***************************************************************************
  ! * type gridCell_t contains the values of phi, theta in radians and radius
  ! * in kilometers for each cell centre on the user-defined grid. Will be used
  ! * as a 3-D array of dimensions (nx,ny,nzEarth)
  type :: gridCell_t

	  integer									:: i,j,k
	  ! the co-ordinates of the cell centre, in radians (ph,th) and km (r)
	  real(8)									:: ph,th,r

  end type gridCell_t	! gridCell_t


Contains


  ! ***************************************************************************
  ! * These subroutines will be used as basic subroutines for various purposes


  ! Evaluate order number of Ylm in a layer

  integer function iY(l,m)

	implicit none

	integer, intent(in)	:: l,m
	integer				:: k

	iY=0
	if ((l<0).or.(m>l)) then
	  write(0,*) 'Error: (iY) wrong usage of spherical harmonics'
	  stop
	else if (l==0) then
	  iY=1
	else
	  do k=0,l-1
		iY=iY+(2*k+1)
	  end do
	  if (m==0) then
		iY=iY+1
	  else if (m<0) then
		iY=iY+2*iabs(m)+1
	  else if (m>0) then
		iY=iY+2*m
	  end if
	end if

  end function iY


  ! Compute the value of the function at given phi and theta values

  real(8) function F_at_point(f,pt)

	use paramfunc
	implicit none

	type (modelFunc_t), intent(in)		  :: f
	type (modelPoint_t),intent(in)		  :: pt

	select case ( f%name )
	case ('Y')
	  F_at_point  = SphHarm(f%l,f%m,pt%phi,pt%theta)
	case default
	  write(0,*) 'Error: (F_at_point) sorry, only spherical harmonics are currently implemented'
	  stop
	end select

  end function F_at_point


  ! Decide whether given radius is in the given parametrization layer

  pure logical function in_layer(rad,layer)

	type (modelLayer_t), intent(in)					:: layer
	real(8), intent(in)								:: rad

	! Consistent with <sigma> layer definitions
	if ((rad>layer%lbound).and.(rad<=layer%ubound)) then

	  in_layer = .TRUE.

	else

	  in_layer = .FALSE.

	end if

  end function in_layer


end module modeldef
