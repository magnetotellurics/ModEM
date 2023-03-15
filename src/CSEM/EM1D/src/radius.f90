!**********************************************************************
!  FD EM function radius
!
!  Purpose           :  distance between points
!
!  Rita Streich 2009
!
!**********************************************************************
real(kind=real64) function radius(x1,y1,z1,x2,y2,z2)

  implicit none

  !external variables
  real(kind=real64)  :: x1,y1,z1,x2,y2,z2

  radius = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

end function radius


!**********************************************************************
!  FD EM function radius0
!
!  Purpose           :  distance of a point from origin
!
!  Rita Streich 2009
!
!**********************************************************************
real(kind=real64) function radius0(x,y,z)

  implicit none

  !external variables
  real(kind=real64)  :: x,y,z

  radius0 = sqrt(x**2 + y**2 + z**2)

end function radius0

