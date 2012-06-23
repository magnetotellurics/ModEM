! *****************************************************************************
module main
	! These subroutines are called from the main program earth3d only

  use input
  use iospec
  use iotypes
  use math_constants
  use model_operators
  use modeldef
  use modelmap
  use global
  use grid_orig
  use field_vectors
  use sg_scalar
  use dataspace
  use dataio
  use userdata
  use dataTypes
  use transmitters
  use receivers
  use senscomp

  implicit none

  ! ***************************************************************************
  ! * targetGrid: Target for all grid pointers for use in the main program
  type (grid_t), target, save                   :: targetGrid

  ! ***************************************************************************
  ! * targetRho0: Target background resistivity for model mappings
  type (rscalar), target, save                  :: targetRho0

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
  ! * The information about the grid, frequency and resistivity map is saved
  ! * in the module global, and the grid is also saved in module griddef.
  ! * Frequencies are listed in ascending order (by value)
  ! * Observatories are listed in ascending order (alphabetic)

  subroutine InitGlobalData(ctrl,eps)

	implicit none
    type (userdef_control), intent(inout)			:: ctrl
	integer										:: i,ios=0,istat=0
	character(100)								:: label
	real(8),intent(in),optional	                :: eps
	type (modelParam_t)                         :: p0_grid
	type (modelCoeff_t)							:: coeff
	type (modelShell_t)                         :: crust

	!--------------------------------------------------------------------------
	! Initialize the math/physics constants defined in module basics
	!call init_const(pi,d2r,r2d,mu0,rair,Rearth)
	!--------------------------------------------------------------------------
	! Initialize the user-defined switches and file names
    !call readStartFile(fn_startup,cUserDef)
    !--------------------------------------------------------------------------
    ! Initialize the global user control and file names structure
    cUserDef = ctrl
	!--------------------------------------------------------------------------
	! Initialize preconditioning information
	call initFunctional(cUserDef,misfitType)
	!--------------------------------------------------------------------------
	! Read and compute grid information
	call initGrid(cUserDef,targetGrid)
    !--------------------------------------------------------------------------
    ! ... and set the grid for sensitivity computations
    call setGrid(targetGrid)
    !--------------------------------------------------------------------------
    ! Right after grid, read the background resistivity model
    call initModelParam(cUserDef,cUserDef%fn_param0,targetGrid,p0_background)
    !--------------------------------------------------------------------------
    ! set the grid for background model parameter
    call setGrid_modelParam(targetGrid,p0_background)
	!--------------------------------------------------------------------------
	! Read and save thin shell conductance information, if present
	call initCrust(cUserDef,targetGrid,crust)
	!--------------------------------------------------------------------------
	! Initialize thin shell conductance in the background resistivity model
	call setCrust_modelParam(crust,p0_background)
    !--------------------------------------------------------------------------
    ! ... and initialize background resistivity on the grid
    call mapToGrid(p0_background,targetRho0)
    !--------------------------------------------------------------------------
    ! Check whether variable parametrization exists
    inquire(FILE=cUserDef%fn_param,EXIST=exists)
	!--------------------------------------------------------------------------
	! Read variable parametrization; if undefined, use background model for
	! forward modelling and to define the parametrization structure
	if (exists) then
	    ! read it in and initialize thinsheet in the variable model parameter
	    call initModelParam(cUserDef,cUserDef%fn_param,targetGrid,p_input)
        call setCrust_modelParam(crust,p_input)
        call setGrid_modelParam(targetGrid,p_input)
        call setBackground_modelParam(targetRho0,p_input)
	else
	    ! use zero starting model parameter
	    p_input = p0_background
	    call zero(p_input)
        call setBackground_modelParam(targetRho0,p_input)
	end if
	!--------------------------------------------------------------------------
	! Adjust the layer boundaries to match the grid
	call adjustLayers_modelParam(p_input,targetGrid%r)
	!--------------------------------------------------------------------------
	! Create a skeleton prior for the inversion
	write(0,*) node_info,'Using zero prior model; background model contains prior information'
	p0_input = p_input
	call zero(p0_input)
	! Very important for model output: mark prior model as "smoothed"; then linear
	! combinations such as Cm^1/2 mhat + m0 are also smoothed.
	p0_input%smoothed = .TRUE.
	!--------------------------------------------------------------------------
	! Compute the correction (only needed if run for a test perturbation)
	if (present(eps)) then
	  p_delta = p_input
	  call random_modelParam(p_delta,eps)
	  call linComb(ONE,p_input,ONE,p_delta,p_input)
	  !call deall_modelParam(p_delta)
	  !param%p(:)%value = param%p(:)%value + da(:)
	end if
	!--------------------------------------------------------------------------
	! 'Smooth' the parametrization by applying inverted regularization operator
	param = multBy_CmSqrt(p_input)
	!--------------------------------------------------------------------------
	! Compute parametrization to use (for model norm we will still use p_input)
	! NO LONGER COMPUTE THIS LINEAR COMB. INSTEAD USE RHO0 AS BACKGROUND!
	! call linComb(ONE,param,ONE,p0_input,param)
    !--------------------------------------------------------------------------
    ! Test to make sure no grid is defined outside the layered region
    if (trim(param%type) .eq. 'harmonic') then
        if (targetGrid%r(targetGrid%nz+1) < param%L(param%nL)%lbound) then
            param%L(param%nL)%lbound = targetGrid%r(targetGrid%nz+1)
        end if
    end if
	!--------------------------------------------------------------------------
	! Compute model information everywhere in the domain, including air & crust
	! call initModel(param,rho,targetGrid,targetRho0)
	rho = targetRho0
    !--------------------------------------------------------------------------
    ! Check whether the optional interior source file exists; read it
    ! (NOT IMPLEMENTED YET)
    inquire(FILE=cUserDef%fn_intsource,EXIST=exists)
    !if (exists) then
    !    call read_sparsevecc(source,cUserDef%fn_intsource)
    !end if
	!--------------------------------------------------------------------------
	! Read the information about the frequencies
	call initFreq(cUserDef,freqList)
	!--------------------------------------------------------------------------
	! Initialize observatory locations required for output
	call initCoords(cUserDef,obsList)
	!--------------------------------------------------------------------------
	! Compute the observatory locations on the grid and interpolation weights
	call initObsList(targetGrid,obsList)
    !--------------------------------------------------------------------------
    ! Initialize reference observatory locations, if present
    call initRefCoords(cUserDef,refObsList)
    !--------------------------------------------------------------------------
    ! Compute the observatory locations on the grid and interpolation weights
    call initObsList(targetGrid,refObsList)
	!--------------------------------------------------------------------------
	! Initialize transfer function information
	call initTF(cUserDef,TFList)
	!--------------------------------------------------------------------------
	! Read radii at which we will output the full solution
	call initSlices(cUserDef,slices)
	!--------------------------------------------------------------------------
	! Read forward solver controls to store in fwdCtrls
	call initControls(fwdCtrls,cUserDef%fn_fwdctrl)
    !--------------------------------------------------------------------------
    ! Read forward solver controls to store in adjCtrls
    inquire(FILE=cUserDef%fn_adjctrl,EXIST=exists)
    if (exists) then
        call initControls(adjCtrls,cUserDef%fn_adjctrl)
    else
        write(0,*) node_info,'Warning: Using forward solver configuration for adjoint computations'
        adjCtrls = fwdCtrls
    end if
	!--------------------------------------------------------------------------
	! Initialize the file names to store in outFiles
    call initOutput(cUserDef,outFiles)
	!--------------------------------------------------------------------------
	nfreq=freqList%n
	nobs =obsList%n
	nfunc=TFList%n
	ncoeff=param%nc

	select case (cUserDef%verbose)
	case ('debug')
	  print *,node_info,'Output all information including debugging lines.'
	  output_level = 5
	case ('full')
	  print *,node_info,'Output full information to screen and to files.'
	  output_level = 4
	case ('regular')
	  print *,node_info,'Output information to files, and progress report to screen (default).'
	  output_level = 3
	case ('compact')
	  print *,node_info,'Output information to files, and compact summary to screen.'
	  output_level = 2
	case ('result')
	  print *,node_info,'Output information to files, and result to screen.'
	  output_level = 1
	case ('none')
	  print *,node_info,'Output nothing at all except result to screen and to files.'
	  output_level = 0
	case default
	  output_level = 3
	end select

	!--------------------------------------------------------------------------
	! Helpful output
#ifdef MPI
    if (taskid==0) then
        !call print_modelParam(p0_input,output_level-1,"Prior model m_0 = ")
        call print_modelParam(p_input,output_level-1,"Input model \hat{m} = ")
        call print_modelParam(param,output_level-1,"Final model m = ")
    end if
#else
	!call print_modelParam(p0_input,output_level-1,"Prior model m_0 = ")
	call print_modelParam(p_input,output_level-1,"Input model \hat{m} = ")
	call print_modelParam(param,output_level-1,"Final model m = ")
#endif
	!--------------------------------------------------------------------------

	! If this information is required, initialize data functionals
	if (cUserDef%calculate == 'original') then
	else
	end if

	! If this information is required, initialize data functionals
	if (cUserDef%calculate == 'original') then
	else
	  !call create_dataVectorMTX(nfreq,nfunc,nobs,dat)
	  !call create_dataVectorMTX(nfreq,nfunc,nobs,psi)
	  !call create_dataVectorMTX(nfreq,nfunc,nobs,res)
	  allocate(ndat(nfreq,nfunc),STAT=istat)
	  allocate(misfitValue(nfunc))
	  call initData(cUserDef,allData,obsList,freqList,TFList)
	  call initMisfit(misfitType,TFList,freqList,allData,misfit)
	  if (cUserDef%calculate /= 'responses') then
		allocate(dmisfitValue(nfunc,ncoeff))
		!call create_dataVectorMTX(nfreq,nfunc,nobs,wres) ! weighted residuals
		allocate(misfit%dRda(nfreq,nfunc,ncoeff),STAT=istat)
		call create_rscalar(targetGrid,sens%drho_real,CENTER)
		call create_rscalar(targetGrid,sens%drho_imag,CENTER)
		allocate(sens%da_real(nfreq,nfunc,nobs,ncoeff),STAT=istat)
		allocate(sens%da_imag(nfreq,nfunc,nobs,ncoeff),STAT=istat)
		allocate(sens%da(nfreq,nfunc,nobs,ncoeff),STAT=istat)
	  end if
	end if

    ! Output the frequencies
    write(0,*) node_info,'nfreq =',nfreq
    do i=1,nfreq
	  write(0,*)node_info,' freq(',trim(freqList%info(i)%code),')=',freqList%info(i)%value
	end do

	! Too complicated to rewrite input to all subroutines that use x,y,z,nx,ny,nz
	! in terms of the grid variable, but that would be the way to do it in the
	! future. Currently, use this patch instead. x,y,z stored in module griddef
	nx = targetGrid%nx; ny = targetGrid%ny; nz = targetGrid%nz
	nzEarth = targetGrid%nzEarth; nzAir = targetGrid%nzAir
	allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
	x = targetGrid%x; y = targetGrid%y; z = targetGrid%z

	call InitGlobalArrays()

	return

  end subroutine InitGlobalData	! InitGlobalData


  ! ***************************************************************************
  ! * DeleteGlobalData deallocates all allocatable data defined globally.
  subroutine DeleteGlobalData()

	integer	:: i,istat

	! Deallocate global variables that have been allocated by InitGlobalData()
	call deall_modelParam(param)
	call deall_modelParam(param0)
	call deall_modelParam(p_input)
	call deall_modelParam(p0_input)
    call deall_modelParam(p0_background)
	call deall_modelParam(p_smooth)
	call deall_modelParam(p_diff)
    call deall_rscalar(rho)
	deallocate(ndat,STAT=istat)
	deallocate(misfit%value,STAT=istat)
	deallocate(misfit%ndat,STAT=istat)
	deallocate(misfit%weight,STAT=istat)
	deallocate(misfit%dRda,STAT=istat)
	deallocate(sens%da_real,STAT=istat)
	deallocate(sens%da_imag,STAT=istat)
	deallocate(sens%da,STAT=istat)
	call deall_rscalar(sens%drho_real)
	call deall_rscalar(sens%drho_imag)
	call deall_dataVectorMTX(allData)

    ! Deallocate dictionaries
	call deall_obsList()
	call deall_freqList()
	call deall_TFList()

	! Deallocate the forward solver
	call cleanUp()

	if (allocated(misfitValue)) then
	  deallocate(misfitValue)
	end if
	if (allocated(dmisfitValue)) then
	  deallocate(dmisfitValue)
	end if

	call DeleteGlobalArrays()

  end subroutine DeleteGlobalData	! DeleteGlobalData


  ! ***************************************************************************
  ! * InitDivergenceCorrection allocates and initializes with zeros all global
  ! * field vectors that are required to keep the divergence correction valid
  subroutine initDivergenceCorrection()

	!use field_vectors
	use dimensions

	integer					  :: istat

	allocate(hx(np1),hy(np1),hz(np1),STAT=istat)
	hx = (0.0d0,0.0d0)
	hy = (0.0d0,0.0d0)
	hz = (0.0d0,0.0d0)
	allocate(sx(np1),sy(np1),sz(np1),STAT=istat)
	sx = (0.0d0,0.0d0)
	sy = (0.0d0,0.0d0)
	sz = (0.0d0,0.0d0)
	allocate(divr(np1),divi(np1),STAT=istat)
	divr = R_ZERO
	divi = R_ZERO
	allocate(divsr(np1),divsi(np1),STAT=istat)
	divsr = R_ZERO
	divsi = R_ZERO

  end subroutine initDivergenceCorrection


  ! ***************************************************************************
  ! * InitGlobalArrays allocates and initializes with zeros all global arrays
  ! * with an allocatable attribute defined in various modules. No allocatable
  ! * arrays should be defined in the main program. Use modules instead.
  subroutine InitGlobalArrays()

	use dimensions
	use boundaries
	use ringcurrent
	!use initFields
	integer	:: istat

	np1=   nx*(ny+1)*(nz+1)	! number of scalars on grid (with some liberal extras)
	np2=3 *nx*(ny+1)*(nz+1)	! number of vectors on grid (with some liberal extras)
	np3=   nx*(ny-1)+2		! number of surface nodes (exact)
  	np4=9 *nx*(ny+1)*(nz+1)
	np5=21*nx*(ny+1)*(nz+1)
	np6=4 *nx*(ny+1)*(nz+1)

	bv1= 2*(ny-1)*(nx-1) + (ny-1)
	bv2= 2*nx*(ny-1) + nx
	bv3= nx*(ny-1)+2

  ! Dynamic allocation of all global vectors
!	allocate(vectorh(np2),vectorb(np2),STAT=istat)

  ! Initialize all global vectors here
!	vectorh = C_ZERO
!	vectorb = C_ZERO

	call initDivergenceCorrection()

  end subroutine InitGlobalArrays ! InitGlobalArrays


  ! ***************************************************************************
  ! * DeleteGlobalArrays does exactly what it is meant to do: deallocates
  ! * all allocatable arrays. No checks are currently performed, since currently
  ! * all memory is deallocated in this subroutine. If memory is ever deallocated
  ! * otherwise, if(allocated(var)) checks are required.
  subroutine DeleteGlobalArrays()

	integer	:: istat

	! Deallocate fields and temporary variables
	deallocate(hx,hy,hz,STAT=istat)
	deallocate(sx,sy,sz,STAT=istat)
	deallocate(divr,divi,STAT=istat)
	deallocate(divsr,divsi,STAT=istat)
!	deallocate(vectorh,vectorb,STAT=istat)

  end subroutine DeleteGlobalArrays	! DeleteGlobalArrays


end module main
