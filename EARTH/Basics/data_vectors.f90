! *****************************************************************************
module data_vectors
	! Module data_vectors contains the variables and the basic subroutines to do
	! with the data vectors and the penalty functional and their derivatives

  use dataspace
  implicit none

  ! ***************************************************************************
  ! * storing data and transfer functions, dim(nfreq,nfunc,nobs)
  ! * Keep in mind that in general nobs(ifreq) <= nobs; use logical indicators.
  !type (dataValue_t), allocatable, dimension(:,:,:)	:: dat,psi,res,wres
  type (dataVectorMTX_t), save				:: dat,psi,res,wres

end module data_vectors
