! *****************************************************************************
module UserData
	! Contains global structures used directly by upper level modules only,
	! including I/O modules and the main program.
	! These are initialized in module Main.

  use DataSpace
  use ModelDef
  use iotypes
  implicit none

  ! ***************************************************************************
  ! *  allData: Contains all input data, saved for computations and output
  ! * outFiles: Information about output file names; extensions are hard-coded
  ! * cUserDef: Vital character-based information specified by the user
  ! * fwdCtrls: User-specified information about the forward solver relaxations
  type (dataVectorMTX_t), save  :: allData
  type (output_info), save      :: outFiles
  type (userdef_control), save       :: cUserDef
  type (fwdCtrl_t), save        :: fwdCtrls,adjCtrls

  ! ***************************************************************************
  ! * misfitType: preconditioning and regularisation parameters
  type (misfitDef_t)                                :: misfitType

  ! ***************************************************************************
  ! * some common integers are stored here to simplify the code
  integer										:: nfreq
  integer               						:: nobs
  integer										:: nfunc
  integer										:: nvar
  integer										:: ncoeff


  ! Dimensions:
  ! nvar = total number of variable parameters (n)
  ! nobs = total number of observatories in the analysis
  !		  (compare to $m_\omega$ = total number of observations for freq. omega)
  ! 3*np1 = number of vector components on cell edges (|E|)
  ! nx*ny*nz = number of scalar values at cell centers (|G|)
  ! nfreq = number of frequency loops

end module UserData
