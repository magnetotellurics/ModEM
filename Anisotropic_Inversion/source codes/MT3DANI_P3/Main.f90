! *****************************************************************************
module Main
	! These subroutines are called from the main program only

  use ModelSpace
  use dataspace ! dataVectorMTX_t
  use datafunc ! to deallocate rxDict
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
     logical, public,  save :: Read_Efield_from_file

  ! ***************************************************************************
  ! * fwdCtrls: User-specified information about the forward solver relaxations
  !type (emsolve_control), save								:: fwdCtrls
  !type (inverse_control), save								:: invCtrls

  ! forward solver control defined in EMsolve3D
  type(EMsolve_control),save  :: solverParams

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
  !  storage for multi-Tx outputed from JT computation
  type(modelParam_t),pointer, dimension(:), save :: JT_multi_Tx_vec

  !  storage for the full sensitivity matrix (dimension nTx)
  type(sensMatrix_t), pointer, dimension(:), save	:: sens

  !  storage for EM solutions
  type(solnVectorMTX_t), save            :: eAll

  !  storage for EM rhs (only for S test)
  type(rhsVectorMTX_t), save            :: RHS

  logical                   :: write_model, write_data, write_EMsoln



Contains
  ! ***************************************************************************
   subroutine read_Efiled_from_file
    character (len=80)			                :: inFile
    character (len=20)                          :: fileVersion
    Integer                                     :: iTx,iMod,filePer
    integer			                            :: fileMode,ioNum
    type(solnVector_t)                          :: e_temp
    real (kind=prec)                        	:: Omega
    integer										:: i,ios=0,istat=0
    
        write(6,*) 'The Master reads E field from: ',trim(solverParams%E0fileName)
        inFile = trim(solverParams%E0fileName)
        call read_solnVectorMTX(Larg_Grid,eAll_larg,inFile)
        nTx_nPol=eAll_larg%nTx*eAll_larg%solns(1)%nPol   
    end subroutine read_Efiled_from_file    

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
! Variables used to read E solution from a file.
	integer                    					:: filePer
    integer			                            :: fileMode,ioNum
    character (len=20)                          :: fileVersion
    character (len=80)			                :: inFile
    Integer                                     :: iTx,iMod
    type(solnVector_t)                          :: e_temp
    real (kind=prec)                        	:: Omega
	!--------------------------------------------------------------------------
	! Set global output level stored with the file units
	output_level = cUserDef%output_level

	!--------------------------------------------------------------------------
	! Check whether model parametrization file exists and read it, if exists
	inquire(FILE=cUserDef%rFile_Model,EXIST=exists)

	if (exists) then
	   ! Read background conductivity parameter and grid
       call read_modelParam(grid,sigma0,cUserDef%rFile_Model)

       ! Finish setting up the grid (if that is not done in the read subroutine)
       call setup_grid(grid)

	else
	  call warning('No input model parametrization')
	end if



	!--------------------------------------------------------------------------
    !  Read forward solver control in EMsolve3D (or use defaults)
    call readEMsolveControl(solverParams,cUserDef%rFile_fwdCtrl,exists,cUserDef%eps)

	!--------------------------------------------------------------------------
	!  Read in data file (only a template on input--periods/sites)
	inquire(FILE=cUserDef%rFile_Data,EXIST=exists)

	if (exists) then
       !  This also sets up dictionaries
	   call read_dataVectorMTX(allData,cUserDef%rFile_Data)
    else
       call warning('No input data file - unable to set up dictionaries')
    end if


    ! Check if a larg grid file with E field is defined:
    ! NOTE: right now both grids share the same transmitters.
    ! This why, reading and setting the large grid and its E solution comes after setting the trasnmitters Dictionary.

     if (solverParams%read_E0_from_File) then
         Read_Efield_from_file = .true.
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
		  call setup_grid(grid) ! Added by Oat
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
           case('L','P')
               inquire(FILE=cUserDef%rFile_EMsoln,EXIST=exists)
               if (exists) then
                  call read_solnVectorMTX(grid,eAll,cUserDef%rFile_EMsoln)
               else
                  call warning('The input EM solution file does not exist')
               end if
           case default
       end select
       select case (cUserDef%option)
           case('S')
            call create_rhsVectorMTX(allData%ntx,RHS)
            do iTx = 1,allData%ntx
                RHS%combs(iTx)%nonzero_source = .true.
                call create_rhsVector(grid,iTx,RHS%combs(iTx))
            end do
            call random_rhsVectorMTX(RHS,cUserDef%eps)
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

	return

  end subroutine initGlobalData	! initGlobalData


  ! ***************************************************************************
  ! * DeallGlobalData deallocates all allocatable data defined globally.

  subroutine deallGlobalData()

	integer	:: i, istat

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
	call deall_txDict() ! 3D_MT/DICT/transmitters.f90
	call deall_rxDict() ! 3D_MT/DICT/receivers.f90
	call deall_typeDict() ! 3D_MT/DICT/dataTypes.f90

	if (associated(sens)) then
		call deall_sensMatrixMTX(sens)
	end if

	call deallEMsolveControl() ! 3D_MT/FWD/EMsolve3D.f90

	call deall_CmSqrt()

    if (output_level > 3) then
       write(0,*) 'All done.'
    endif

  end subroutine deallGlobalData	! deallGlobalData


end module Main
