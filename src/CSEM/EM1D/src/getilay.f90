!------------------------------------------------------------------
! FD EM function getilay
!   find out which layer a depth point is in
!
! Rita Streich 2009
!------------------------------------------------------------------
function getilay(z,zbound) result(lay)

  implicit none

  !external variables
  integer(kind=int32)   :: lay !the layer whcih z is in
  real(kind=real64)     :: z   !the depth to test
  real(kind=real64),dimension(:)  :: zbound

  !internal variable
  integer(kind=int32)   :: ilay !layer counter
  integer(kind=int32)   :: nzb  !number of layer boundaries

  nzb = size(zbound)

  lay = 1
  do ilay=1,nzb
    !if(z .gt. zbound(ilay)) lay = lay + 1
    ! for source right on boundary, define source to be in lower medium
    ! --> this is more stable and more appropriate for grounded sources
    if(z .ge. zbound(ilay)) lay = lay + 1
  enddo

endfunction getilay
