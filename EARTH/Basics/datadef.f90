! *****************************************************************************
module datadef
	! This module contains the type definitions for the structures that contain
	! the full information about the data and data functionals, frequencies and
	! other data-related information

  use griddef
  use sg_scalar
  use sg_sparse_vector
  implicit none

  ! --------------------------------------
  ! Define in a globally accessible module
  ! integer										:: nobs
  ! integer										:: nfreq
  ! integer										:: nfunc
  ! type (dataValue_t), dimension(nfreq,nobs)	:: dat
  ! type (Obs_List)								:: obsList
  ! type (Freq_List)							:: freqList
  ! type (TF_List)								:: funcList
  ! Instances of data_func: datC, psiC, resC etc
  ! dat is data with errors
  ! psi are computed functionals
  ! res are residuals
  ! These are vectors dimension(nfreq,nobs)
  ! Use datC%resp%exists == .FALSE. to indicate
  ! that data is not present at a given
  ! frequency + observatory pair.
  ! --------------------------------------





end module datadef
