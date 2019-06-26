program test_sub01
  use Sub

  real :: x, y, z
  real, parameter :: tol = 1E-4
  real :: tval

  x = 3.0
  y = 4.0
  tval = 5.0

  z = sub01(x, y)

  if(abs(z - tval) <= tol) then
     STOP 0
  else
     STOP 1
  end if
  
end program test_sub01
