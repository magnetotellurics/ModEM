! *****************************************************************************
module Main
	! These subroutines are called from the main program only

  use ModelSpace
  use dataspace ! dataVectorMTX_t
  use datafunc ! to deallocate rxDict, typeDict
  use ForwardSolver ! txDict, solnVectorMTX
  use sensmatrix
  use userctrl
  use ioascii
  use dataio
  implicit none

      ! I/O units ... reuse generic read/write units if
     !   possible; for those kept open during program run,
     !   reserve a specific unit here
     integer (kind=4), save :: fidRead = 1
     integer (kind=4), save :: fidWrite = 2
     integer (kind=4), save :: fidError = 99

  ! ***************************************************************************
  ! * fwdCtrls: User-specified information about the forward solver relaxations
  !type (emsolve_control), save								:: fwdCtrls
  !type (inverse_control), save								:: invCtrls

  ! this is used to set up the numerical grid in SensMatrix
  type(grid_t), save	        :: grid

  ! impedance data structure
  type(dataVectorMTX_t), save		:: allData

  !  storage for the "background" conductivity parameter
  type(modelParam_t), save		:: sigma0
  !  storage for a perturbation to conductivity
  type(modelParam_t), save		:: dsigma
  !  storage for the inverse solution
  type(modelParam_t), save		:: sigma1
  !  currently only used for TEST_GRAD feature (otherwise, use allData)
  type(dataVectorMTX_t), save       :: predData
  !  also for TEST_GRAD feature...
  type(modelParam_t), save      :: sigmaGrad
  real(kind=prec), save         :: rms,mNorm,f1,f2,alpha
  !  storage for multi-Tx outputed from JT computation
  type(modelParam_t),pointer, dimension(:), save :: JT_multi_Tx_vec

  !  storage for the full sensitivity matrix (dimension nTx)
  type(sensMatrix_t), pointer, dimension(:), save	:: sens

  !  storage for EM solutions
  type(solnVectorMTX_t), save            :: eAll

  logical                   :: write_model, write_data, write_EMsoln, write_EMrhs



Contains

  ! ***************************************************************************
  ! * InitGlobalData is the routine to call to initialize all derived data types
  ! * and other variables defined in modules basics, modeldef, datadef and
  ! * in this module. These include:
  ! * 1) constants
  ! * 2) grid information: nx,ny,nz,x,y,z, cell centre coordinates
  ! * 3) model information: layers, coefficients, rho on the grid
  ! * 4) periods/freq info: nfreq, freq
  ! * 5) forward solver controls
  ! * 6) output file names

  subroutine initGlobalData(cUserDef)

	implicit none
	type (userdef_control), intent(in)          :: cUserDef
	integer										:: i,ios=0,istat=0
	logical                                     :: exists
	character(80)                               :: paramtype,header
	integer                                     :: nTx

	!--------------------------------------------------------------------------
	! Set global output level
	output_level = cUserDef%output_level

	!--------------------------------------------------------------------------
	! Check whether model parametrization file exists and read it, if exists
	inquire(FILE=cUserDef%rFile_Model,EXIST=exists)

	if (exists) then
	   ! Read background conductivity parameter and grid
       call read_modelParam(grid,sigma0,cUserDef%rFile_Model)

       !  set array size parameters in WS forward code module
       !   these stay fixed for all forward modeling with this grid
       call setWSparams(grid%Ny,grid%Nz,grid%Nza)

       ! Finish setting up the grid (if that is not done in the read subroutine)
       call setup_grid(grid)

	else
	  call warning('No input model parametrization')

	  ! set up an empty grid to avoid segmentation faults in sensitivity tests
          call create_grid(1,1,1,grid)
	end if

	!--------------------------------------------------------------------------
	!  Read in data file (only a template on input--periods/sites)
	inquire(FILE=cUserDef%rFile_Data,EXIST=exists)

	if (exists) then
       !  This also sets up dictionaries
	   call read_dataVectorMTX(allData,cUserDef%rFile_Data)
    else
       call warning('No input data file - unable to set up dictionaries')
    end if

	!--------------------------------------------------------------------------
	!  Initialize additional data as necessary
	select case (cUserDef%job)

     case (MULT_BY_J)
	   inquire(FILE=cUserDef%rFile_dModel,EXIST=exists)
	   if (exists) then
	      call deall_grid(grid)
	   	  call read_modelParam(grid,dsigma,cUserDef%rFile_dModel)
	   else
	      call warning('The input model perturbation file does not exist')
	   end if


     case (INVERSE)
       inquire(FILE=cUserDef%rFile_Cov,EXIST=exists)
       if (exists) then
          call create_CmSqrt(sigma0,cUserDef%rFile_Cov)
       else
          call create_CmSqrt(sigma0)
       end if
	   inquire(FILE=cUserDef%rFile_dModel,EXIST=exists)
	   if (exists) then
	      call deall_grid(grid)
	   	  call read_modelParam(grid,dsigma,cUserDef%rFile_dModel)
	      if (output_level > 0) then
	        write(*,*) 'Using the initial model perturbations from file ',trim(cUserDef%rFile_dModel)
	      endif
	   else
	      dsigma = sigma0
	      call zero_modelParam(dsigma)
	      if (output_level > 0) then
	        write(*,*) 'Starting search from the prior model ',trim(cUserDef%rFile_Model)
	      endif
	   end if
       select case (cUserDef%search)
       case ('NLCG')
       case ('DCG')
       case ('Hybrid')
       case default
          ! a placeholder for anything specific to a particular inversion algorithm;
          ! currently empty
       end select

     case (APPLY_COV)
       inquire(FILE=cUserDef%rFile_Cov,EXIST=exists)
       if (exists) then
          call create_CmSqrt(sigma0,cUserDef%rFile_Cov)
       else
          call create_CmSqrt(sigma0)
       end if
       dsigma = sigma0
       inquire(FILE=cUserDef%rFile_Prior,EXIST=exists)
       if (exists) then
           call deall_grid(grid)
           call read_modelParam(grid,sigma0,cUserDef%rFile_Prior)
       else
           call zero(sigma0)
       end if
       sigma1 = sigma0
       call zero(sigma1)

     case (TEST_GRAD, TEST_SENS)
         inquire(FILE=cUserDef%rFile_dModel,EXIST=exists)
         if (exists) then
             call deall_grid(grid)
             call read_modelParam(grid,dsigma,cUserDef%rFile_dModel)
         else
             call warning('The input model perturbation file does not exist')
         end if

     case (TEST_ADJ)
       select case (cUserDef%option)
           case('J','Q','P')
               inquire(FILE=cUserDef%rFile_dModel,EXIST=exists)
               if (exists) then
                  call deall_grid(grid)
                  call read_modelParam(grid,dsigma,cUserDef%rFile_dModel)
               else
                  call warning('The input model perturbation file does not exist')
               end if
           case default
       end select
       select case (cUserDef%option)
           case('L','P','e')
               inquire(FILE=cUserDef%rFile_EMsoln,EXIST=exists)
               if (exists) then
                  call read_solnVectorMTX(grid,eAll,cUserDef%rFile_EMsoln)
               else
                  call warning('The input EM solution file does not exist')
               end if
           case('S','b')
               inquire(FILE=cUserDef%rFile_EMrhs,EXIST=exists)
               if (exists) then
                  call read_rhsVectorMTX(grid,bAll,cUserDef%rFile_EMrhs)
               else
                  call warning('The input EM RHS file does not exist')
               end if
           case default
       end select

    end select

	!--------------------------------------------------------------------------
	!  Only write out these files if the output file names are specified
    write_model = .false.
    if (len_trim(cUserDef%wFile_Model)>1) then
       write_model = .true.
    end if
    write_data = .false.
    if (len_trim(cUserDef%wFile_Data)>1) then
       write_data = .true.
    end if
    write_EMsoln = .false.
    if (len_trim(cUserDef%wFile_EMsoln)>1) then
       write_EMsoln = .true.
    end if
    write_EMrhs = .false.
    if (len_trim(cUserDef%wFile_EMrhs)>1) then
       write_EMrhs = .true.
    end if

	return

  end subroutine initGlobalData	! initGlobalData


  ! ***************************************************************************
  ! * DeallGlobalData deallocates all allocatable data defined globally.

  subroutine deallGlobalData()

	integer	:: i, istat

    write(0,*) 'Cleaning up...'

	! Deallocate global variables that have been allocated by InitGlobalData()
	if (output_level > 3) then
	   write(0,*) 'Cleaning up grid...'
	endif
	call deall_grid(grid)

	if (output_level > 3) then
	   write(0,*) 'Cleaning up data...'
	endif
	call deall_dataVectorMTX(allData)

	if (output_level > 3) then
	   write(0,*) 'Cleaning up EM soln...'
	endif
	call deall_solnVectorMTX(eAll)

	if (output_level > 3) then
	   write(0,*) 'Cleaning up models...'
	endif
	call deall_modelParam(sigma0)
	call deall_modelParam(dsigma)
	call deall_modelParam(sigma1)

    if (output_level > 3) then
       write(0,*) 'Cleaning up dictionaries...'
    endif
	call deall_txDict() ! 2D_MT/DICT/transmitters.f90
	call deall_rxDict() ! 2D_MT/DICT/receivers.f90
	call deall_typeDict() ! 2D_MT/DICT/dataTypes.f90

	if (associated(sens)) then
		call deall_sensMatrixMTX(sens)
	end if

	call deall_CmSqrt()

    if (output_level > 3) then
       write(0,*) 'All done.'
    endif

  end subroutine deallGlobalData	! deallGlobalData


end module Main
