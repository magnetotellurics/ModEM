module Sub
  
contains
  function sub01(x, y) result(r)
    real :: x, y
    real :: r

    r = sqrt(x**2 + y**2)
    
  end function sub01
end module Sub
