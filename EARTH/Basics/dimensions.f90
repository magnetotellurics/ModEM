! *****************************************************************************
module dimensions
	! This module only contains the definitions of integer dimensions.

  implicit none

  ! ***************************************************************************
  ! * define the vector dimensions used throughout in the subroutines

  ! np1=nx*(ny+1)*(nz+1);np2=3*np1;np3=nx*(ny+1);np4=9*np1;np5=21*np1;np6=4*np1
  ! in practice, using np3 = nx(ny-1)+2 - number of surface nodes
  integer									    :: np1,np2,np3,np4,np5,np6

  ! dimensions for boundary values
  integer										:: bv1,bv2,bv3
  !---------------------------------------------------------
  ! np3=nx*(ny+1); use for variables like
  !c    complex exs(nx*ny),eys(nx*ny),ezs(nx*ny)
  !c    complex hxs(nx*ny),hys(nx*ny),hzs(nx*ny)
  !-----------> surface EM field
  !---------------------------------------------------------


end module dimensions