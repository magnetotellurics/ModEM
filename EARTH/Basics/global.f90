! *****************************************************************************
module global
	! Module global is used by all other modules, except for the input module,
	! which initializes the variables. There, to make the code clear, we pass
	! them as arguments. Otherwise, these variables are assumed known.
	! This is crucial, since some of the data contains pointers to values
	! stored here. These global variables are initialized in the module main.

  use griddef
  use modeldef
  use datadef
  use iotypes
  implicit none

  !save

  ! ***************************************************************************
  ! * dictionaries of receivers, transmitters and transfer functions
  !type (Obs_List), save								:: obsList
  !type (Freq_List), save								:: freqList
  !type (TF_List), save								:: TFList

  ! ***************************************************************************
  ! * user-defined list of radii in meters at which we output the solution
  !type (Rad_List), save								:: slices

  ! ***************************************************************************
  ! * rho: storing the model on the grid - resistivity defined in cell centres
  !real(8), allocatable, dimension(:,:,:), save		:: rho	!(nx,ny,nz)
  type(rscalar), save     :: rho

  ! ***************************************************************************
  ! * grid: Contains the information about the user-specified grid
  !type (grid_t), target, save							:: targetGrid

  ! ***************************************************************************
  ! * shell: Contains the information about thin-shell conductance
  ! - now part of the model parameter
  !type (modelShell_t), save							:: crust

  ! ***************************************************************************
  ! * param: Contains the information about the user-specified parametrization
  type (modelParam_t), save								:: p_input
  type (modelParam_t), save                             :: p0_background
  type (modelParam_t), save								:: p0_input
  type (modelParam_t), save								:: p_smooth
  type (modelParam_t), save								:: p_diff
  type (modelParam_t), save								:: param
  type (modelParam_t), save								:: param0
  type (modelParam_t), save                             :: p_delta

  ! ***************************************************************************
  ! *  allData: Contains all input data, saved for computations and output
  ! * outFiles: Information about output file names; extensions are hard-coded
  ! * cUserDef: Vital character-based information specified by the user
  ! * fwdCtrls: User-specified information about the forward solver relaxations
  !type (dataVectorMTX_t), save  :: allData
  !type (output_info), save      :: outFiles
  !type (userdef_control), save       :: cUserDef
  !type (fwdCtrl_t), save        :: fwdCtrls


  ! ***************************************************************************
  ! * some common integers are stored here to simplify the code
  !integer										:: nfreq
  !integer               						:: nobs
  !integer										:: nfunc
  !integer										:: nvar
  !integer										:: ncoeff


  ! Dimensions:
  ! nvar = total number of variable parameters (n)
  ! nobs = total number of observatories in the analysis
  !		  (compare to $m_\omega$ = total number of observations for freq. omega)
  ! 3*np1 = number of vector components on cell edges (|E|)
  ! nx*ny*nz = number of scalar values at cell centers (|G|)
  ! nfreq = number of frequency loops

end module global
