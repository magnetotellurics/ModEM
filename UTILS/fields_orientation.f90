! ***************************************************************************************
! fields_orientation.f90 defines general local azimuths and tilts of all field components

module fields_orientation
  use math_constants
  implicit none

  ! Defines a general real value associated with each local field component
  type :: fieldsInfo_t

      real (kind=prec)  :: Hx, Hy, Hz, Ex, Ey, Hx_ref, Hy_ref

  end type fieldsInfo_t


  ! Defines complete orientation in space for a local set of geomagnetic fields.
  ! Here, azimuth is always defined relative to geographic reference frame.
  ! The azimuth refers to the horizontal angle relative to north positive clockwise,
  ! and the tilt refers to the vertical angle with respect to the horizontal plane.
  ! In this reference frame, a tilt angle of 90 points downward, 0 is parallel
  ! with the surface, and -90 points upwards. All angles are in degrees.
  ! REF: Peacock et al (2022) https://doi.org/10.1016/j.cageo.2022.105102
  type :: orient_t

      type(fieldsInfo_t)  :: azimuth
      type(fieldsInfo_t)  :: tilt

  end type orient_t

  ! This is a global variable that needs to be set up only once
  type(orient_t), save, public   :: orient0

Contains

  ! **************************************************************************
  ! Initializes and sets up default orientation of all field components.
  ! This is the orientation that is used for internal computations as well.
  subroutine setup_default_orientation(orient0)

    type(orient_t), intent(inout)   :: orient0

    orient0%azimuth%Hx = R_ZERO
    orient0%azimuth%Hy = 90.0*ONE
    orient0%azimuth%Hz = R_ZERO
    orient0%azimuth%Ex = R_ZERO
    orient0%azimuth%Ey = 90.0*ONE
    orient0%azimuth%Hx_ref = R_ZERO
    orient0%azimuth%Hy_ref = 90.0*ONE

    orient0%tilt%Hx = R_ZERO
    orient0%tilt%Hy = R_ZERO
    orient0%tilt%Hz = 90.0*ONE
    orient0%tilt%Ex = R_ZERO
    orient0%tilt%Ey = R_ZERO
    orient0%tilt%Hx_ref = R_ZERO
    orient0%tilt%Hy_ref = R_ZERO

  end subroutine setup_default_orientation

  ! **************************************************************************
  ! Logical function to compare azimuths and tilts
  logical function compare_orientation(orient1,orient2)

      type(orient_t), intent(in)   :: orient1,orient2

      compare_orientation = .TRUE.

      if ((abs(orient1%azimuth%Hx - orient2%azimuth%Hx)>TOL4) .or. &
          (abs(orient1%azimuth%Hy - orient2%azimuth%Hy)>TOL4) .or. &
          (abs(orient1%azimuth%Hz - orient2%azimuth%Hz)>TOL4) .or. &
          (abs(orient1%azimuth%Ex - orient2%azimuth%Ex)>TOL4) .or. &
          (abs(orient1%azimuth%Ey - orient2%azimuth%Ey)>TOL4) .or. &
          (abs(orient1%azimuth%Hx_ref - orient2%azimuth%Hx_ref)>TOL4) .or. &
          (abs(orient1%azimuth%Hy_ref - orient2%azimuth%Hy_ref)>TOL4)) then
          compare_orientation = .FALSE.
      endif

      if ((abs(orient1%tilt%Hx - orient2%tilt%Hx)>TOL4) .or. &
          (abs(orient1%tilt%Hy - orient2%tilt%Hy)>TOL4) .or. &
          (abs(orient1%tilt%Hz - orient2%tilt%Hz)>TOL4) .or. &
          (abs(orient1%tilt%Ex - orient2%tilt%Ex)>TOL4) .or. &
          (abs(orient1%tilt%Ey - orient2%tilt%Ey)>TOL4) .or. &
          (abs(orient1%tilt%Hx_ref - orient2%tilt%Hx_ref)>TOL4) .or. &
          (abs(orient1%tilt%Hy_ref - orient2%tilt%Hy_ref)>TOL4)) then
          compare_orientation = .FALSE.
      endif

  end function compare_orientation

end module fields_orientation !fields_orientation
