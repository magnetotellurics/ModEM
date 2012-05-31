! *****************************************************************************
module input
  ! This module defines input file information and input routines

  use math_constants
  use file_units
  use griddef
  use utilities
  use modeldef
  use dataMisfit
  use DataSpace
  use SolnSpace
  use iotypes
  use UserData
  use dataTypes
  use transmitters
  use receivers
  implicit none

  !integer, parameter							:: ioStartup=101
  !integer, parameter							:: ioMdl=1
  !integer, parameter							:: ioGrd=2
  !integer, parameter							:: ioShell=3
  !integer, parameter							:: ioPrm=4
  !integer, parameter							:: ioPt=32
  !integer, parameter							:: ioPer=31
  !integer, parameter							:: ioCtrl=16
  !integer, parameter							:: ioCond=23
  !integer, parameter							:: ioDat=17
  !integer, parameter							:: ioObs=18
  !integer, parameter							:: ioFunc=19
  !integer, parameter							:: ioRad=15

Contains

  ! ***************************************************************************
  ! * readStartFile reads the filename from the screen, if it is not specified.
  ! * This file contains essential information for the program to run.

  subroutine readStartFile(fn_startup,cUserDef)

    character(*), intent(inout)  		:: fn_startup
	type (userdef_control), intent(out)		:: cUserDef
    integer								:: ios
	character(20)						:: string

    ! passed an empty string
    if (fn_startup == '') then
       ! Prompt the user to give the name for the startup file
       write(0, *) '**********************************************************'
       write(0, *) 'Please, ENTER THE START FILE:'
       read(*, '(a80)') fn_startup
    end if

    open (unit=ioStartup,file=fn_startup,status='old',iostat=ios)

    if(ios/=0) then
       write(0,*) 'Error opening file:', fn_startup
    endif

    ! This is the list of options specified in the startup file
    read (ioStartup,'(a17,a80)') string,cUserDef%paramname;
    read (ioStartup,'(a17,a80)') string,cUserDef%modelname;
    read (ioStartup,'(a17,a80)') string,cUserDef%verbose;
    read (ioStartup,'(a17,a80)') string,cUserDef%calculate;
    read (ioStartup,'(a17,a80)') string,cUserDef%secondary_field;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_thinsheet;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_grid;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_param0;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_param;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_period;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_extsource;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_intsource;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_fwdctrl;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_adjctrl;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_invctrl;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_slices;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_coords;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_func;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_cdata;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_ddata;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_hdata;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_href;

    ! if no inverse control file, use this default damping parameter
    cUserDef%damping = 10.

    ! old output files; probably won't be needed anymore
    cUserDef%fn_misfit = 'MISFIT'
    cUserDef%fn_gradient = 'GRADIENT'
    cUserDef%fn_point = 'POINT'

    close(ioStartup)

  end subroutine readStartFile  ! readStartFile


  ! ***************************************************************************
  ! * initGrid reads the modelfile fn_grid to store the grid inGrid only.
  ! * Traditionally, for global spherical grid, assume the following directions:
  ! * x -> phi (longitude, varies from 0 to 360)
  ! * y -> theta (co-latitude = 90 - latitude; varies from 0 to 180)
  ! * z -> -r (radius from Earth's centre, r=2871.0 at CMB,
  ! *									     6371.0 at Earth/air interface)
  ! *
  ! * For consistency with the original Randy Mackie, 06-27-85, 3-D model, and
  ! * also for consistency with the current forward solver subroutines, we are
  ! * keeping the following grid structure in this forward solver:
  ! *      line 1: dimensions of model (nx,ny,nzAir,nzCrust,nzEarth)
  ! *      line 2: x(*) in degrees (interval)
  ! *      line 3: y(*) in degrees (position from n-pole)
  ! *      line 4: z(*) in km (distance from center of the earth, decreasing)

  subroutine initGrid(cUserDef,mygrid)

    type (userdef_control), intent(in)				:: cUserDef
    type (grid_t) , intent(out)					    :: mygrid
    integer				                            :: i

    call read_grid(mygrid,cUserDef%fn_grid)

	inquire(FILE=cUserDef%fn_thinsheet,EXIST=exists)
	! If no thin sheet distribution present, include the thinsheet layers into the model domain
	if (.not.exists) then
	  write(0,*) node_info,'Warning: No thin shell conductance distribution specified; assume no crust'
      mygrid%nzEarth = mygrid%nzCrust + mygrid%nzEarth
      mygrid%nzCrust = 0
	end if

  end subroutine initGrid	! initGrid

  ! ***************************************************************************
  ! * initRho reads the resistivities on the grid from a simple file

  subroutine initRho(cUserDef,mygrid,myrho)

    type (userdef_control), intent(in)          :: cUserDef
    type (grid_t), intent(in)                   :: mygrid
    type (rscalar), intent(inout)               :: myrho  !(nx,ny,nz)
    ! local
    real(8)                                     :: lon,lat,depth
    integer                                     :: ios,i,j,k,nx,ny,nzEarth


    inquire(FILE=cUserDef%fn_param,EXIST=exists)
    if(.not.exists) then
      write(6,*) node_info,'Model resistivities will not be initialized: ',trim(cUserDef%fn_param)," not found"
      return
    end if

    open(ioMdl,file=cUserDef%fn_param,status='old',iostat=ios)

    write(6,*) node_info,'Reading the resistivities from file ',trim(cUserDef%fn_param)
    read(ioMdl,*) ! header line
    read(ioMdl,*) nx,ny,nzEarth

    if ((nx .ne. mygrid%nx) .or. (ny .ne. mygrid%ny) .or. (nzEarth .ne. (mygrid%nzCrust+mygrid%nzEarth))) then
      write(6,*) 'Warning: Model resistivities do not match grid size in ',trim(cUserDef%fn_param)
    end if

    call create_rscalar(mygrid,myrho,CENTER)
    myrho%v(:,:,:) = 1./SIGMA_AIR
    do k=mygrid%nzAir+1,mygrid%nz
      do i=1,mygrid%nx
        do j=1,mygrid%ny
          read(ioMdl,*) lon,lat,depth,myrho%v(i,j,k)
        end do
      end do
    end do

    close(ioMdl)

  end subroutine initRho  ! initRho

  ! ***************************************************************************
  ! * initField reads in the full field solution from fn_field.
  ! * This is used for reading in the primary (radial) magnetic fields.
  ! * If the file exists, it contains the fields for the prior model.

  subroutine initField(cUserDef,mygrid,H)

    type (userdef_control), intent(in)					:: cUserDef
    type (grid_t) , intent(in)						:: mygrid
    type (solnVectorMTX_t) , intent(inout)				:: H
    ! local
    integer				                            :: ios,istat,i

    write(6,*) node_info,'Reading from the EM solution file ',trim(cUserDef%fn_field)

	inquire(FILE=cUserDef%fn_field,EXIST=exists)
	if(.not.exists) then
      write(6,*) node_info,'Field solution will not be initialized: ',trim(cUserDef%fn_field)," not found"
	  call deall_solnVectorMTX(H)
	  return
	end if

	call read_solnVectorMTX(cUserDef%fn_field,H,mygrid)
	call write_solnVectorMTX('test.field',H)

  end subroutine initField	! initField

  ! ***************************************************************************
  ! * initCrust reads the modelfile fn_thinsheet to store S-conductance of Earth's
  ! * crust in GM coordinates. It also stores the depth of the crust in km.
  ! * This should be called after initializing the grid, since it uses the grid
  ! * dimensions to determine the number of entries in file.

  subroutine initCrust(cUserDef,mygrid,mycrust)

    type (userdef_control), intent(in)					:: cUserDef
    type (grid_t) , intent(in)					:: mygrid
    type (modelShell_t) , intent(out)					:: mycrust
    integer				                            :: ios,istat,i


	inquire(FILE=cUserDef%fn_thinsheet,EXIST=exists)
	if(.not.exists) then
	  mycrust%allocated = .false.
	  return
	end if

	allocate(mycrust%cond(mygrid%nx,mygrid%ny),STAT=istat)

	open(ioShell,file=cUserDef%fn_thinsheet,status='old',iostat=ios)

    write(6,*) node_info,'Reading from the thin sheet conductance file ',trim(cUserDef%fn_thinsheet)

	do i=1,mygrid%nx
	  read(ioShell,*,iostat=ios) mycrust%cond(i,1:mygrid%ny)
	end do

	close(ioShell)

    mycrust%variable = .false. ! default; may be changed later
	mycrust%allocated = .true.
	return

  end subroutine initCrust	! initCrust


  ! ***************************************************************************
  ! * initMisfit initializes the values required to compute the full penalty
  ! * functional
  ! * Uses: freqList,TFList
  subroutine initMisfit(misfitType,TFList,freqList,allData,misfit)

	type (misfitDef_t), intent(in)						:: misfitType
	type (dataVectorMTX_t), intent(in)						:: allData
	type (Freq_List), intent(in)						:: freqList
	type (TF_List), intent(in)							:: TFList
	type (misfit_t), intent(out)						:: misfit
	!integer, dimension(:,:), intent(out)				:: ndat
	integer												:: ifunc,ifreq,iobs
	integer												:: istat,i,j,iTx,iDt
	real(8)												:: total_weight

	allocate(misfit%value(freqList%n,TFList%n),STAT=istat)
	allocate(misfit%ndat(freqList%n,TFList%n),STAT=istat)
	allocate(misfit%weight(TFList%n),STAT=istat)

	misfit%name = misfitType%name
	misfit%damping = misfitType%mu

	misfit%ndat(:,:) = 0
	do i = 1,allData%nTx
	   do j = 1,allData%d(i)%nDt
	       iTx = allData%d(i)%tx
	       iDt = allData%d(i)%data(j)%dataType
	       misfit%ndat(iTx,iDt) = allData%d(i)%data(j)%nSite * allData%d(i)%data(j)%nComp / 2 !countData(allData%d(i)%data(j))
	   end do
	end do

	misfit%value = 0.0d0

	total_weight = dot_product(TFList%info%w,sum(misfit%ndat,1))
	do ifunc = 1,TFList%n
	  misfit%weight(ifunc) = TFList%info(ifunc)%w * sum(misfit%ndat(:,ifunc))/total_weight
	end do

  end subroutine initMisfit

  ! ***************************************************************************
  ! * initFunctional specifies basic information about the data functional

  subroutine initFunctional(cUserDef,misfitType)

	implicit none
	type (userdef_control), intent(in)					:: cUserDef
	type (misfitDef_t), intent(out)					:: misfitType

	misfitType%name = 'Mean Squared'  ! Default data functional

	misfitType%mu   = cUserDef%damping

  end subroutine initFunctional	! initFunctional




  ! ***************************************************************************
  ! * initModelParam reads the parametrization info to store in the derived data
  ! * type variables (modelCoeff_t, modelLayer_t) in the special case when the
  ! * parametrization is in terms of spherical harmonics of various degree/order
  ! * per layer, each of the parameters being either variable (keyword 'range') or
  ! * constant (keyword 'const'). We keep this information internally in such
  ! * a format, that generalizations of this case would be easy to implement.
  ! * In the parametrization type variables, each layer has the same (maximum)
  ! * degree and order of spherical harmonics, and hence identical numbers of
  ! * coefficients to store. Out of these coefficients, those that have been
  ! * defined in the script bear logical keyword 'exists'. Out of those, variable
  ! * coefficients have logical 'frozen==.FALSE.'
  ! * Obviously, the information on the range is not required for the forward
  ! * solver to operate. This is provided for the inversion, which will share
  ! * the same input format for now.

  subroutine initModelParam(cUserDef,mygrid,myparam,p0)

	use model_operators

    type (userdef_control), intent(in)					:: cUserDef
    type (grid_t), target, intent(in)                   :: mygrid
    type (modelParam_t), intent(inout)					:: myparam
	logical, intent(in), optional		:: p0


    if (trim(cUserDef%paramname) .eq. 'harmonic') then

        if(present(p0)) then
          if(p0) then
            call read_modelParam(myparam,cUserDef%fn_param0)
          end if
        else
          call read_modelParam(myparam,cUserDef%fn_param)
        end if
        myparam%grid => mygrid

    else if (trim(cUserDef%paramname) .eq. 'grid') then

        call initRho(cUserDef,mygrid,myparam%rho)
        myparam%allocated = .true.

    else
        write(0,*) node_info,'Warning: model parametrization ',trim(cUserDef%paramname),' not implemented yet'
        stop
    end if

    myparam%type = trim(cUserDef%paramname)


  end subroutine initModelParam	! initModelParam


  ! ***************************************************************************
  ! * initControls reads the file fn_fwdctrl and sets the values of main control
  ! * parameters for the forward solver, stored in fwdCtrls

  subroutine initControls(fwdCtrls,fname)

	implicit none
	character(len=*), intent(in)						:: fname
	type (fwdCtrl_t), intent(out)						:: fwdCtrls

	  open(ioFwdCtrl,file=fname,form='formatted',status='old')

      write(6,*) node_info,'Reading from the forward solver controls file ',trim(fname)
      read(ioFwdCtrl,*) fwdCtrls%ipotloopmax ! max number of divergence correction loops
      read(ioFwdCtrl,*) fwdCtrls%errend ! tolerance on the solution update (herr)
      read(ioFwdCtrl,*) fwdCtrls%nrelmax    ! max number of solution updates between div. corrections
      read(ioFwdCtrl,*) fwdCtrls%n_reldivh  ! number of divergence correction iterations
      read(ioFwdCtrl,*) fwdCtrls%ipot0,fwdCtrls%ipotint,fwdCtrls%ipot_max ! check and run div. correction
            ! if needed every ipotmax solution updates, ipotmax = min(ipot0+ndivcorr*ipotint,ipot_max)

      ! write them all out to save in the output
      write(6,*) node_info,fwdCtrls%ipotloopmax
      write(6,*) node_info,fwdCtrls%errend
      write(6,*) node_info,fwdCtrls%nrelmax
      write(6,*) node_info,fwdCtrls%n_reldivh
      write(6,*) node_info,fwdCtrls%ipot0,fwdCtrls%ipotint,fwdCtrls%ipot_max

      close(ioFwdCtrl)

  end subroutine initControls	! initControls


  ! ***************************************************************************
  ! * initOutput initialises all the output file names stored in outFiles
  ! * Assuming type (output_info) outFiles and type (userdef_control) cUserDef are
  ! * not available.
  subroutine initOutput(cUserDef,outFiles)

	implicit none
    type (userdef_control) ,intent(in)					:: cUserDef
    type (output_info),intent(out)					:: outFiles
	integer											:: i

    if (cUserDef%modelname == '') then
       write(0, *) node_info,'modelname not specified yet for FileInfoInit'
       stop
    end if

	i=index(cUserDef%modelname,' ')
	i=i-1

	outFiles%fn_hx	  =cUserDef%modelname(1:i)//'.hx';
	outFiles%fn_hy	  =cUserDef%modelname(1:i)//'.hy';
	outFiles%fn_hz	  =cUserDef%modelname(1:i)//'.hz';
	outFiles%fn_jxjyjz =cUserDef%modelname(1:i)//'.jxjyjz';
	outFiles%fn_hxhyhz =cUserDef%modelname(1:i)//'.hxhyhz';
	outFiles%fn_err	  =cUserDef%modelname(1:i)//'.err';
	outFiles%fn_cresp  =cUserDef%modelname(1:i)//'.cresp';
	outFiles%fn_dresp  =cUserDef%modelname(1:i)//'.dresp';
	outFiles%fn_cdat  =cUserDef%modelname(1:i)//'.cout';
	outFiles%fn_ddat  =cUserDef%modelname(1:i)//'.dout';
    outFiles%fn_hdat  =cUserDef%modelname(1:i)//'.hout';
	outFiles%fn_cjac  =cUserDef%modelname(1:i)//'.cj';
	outFiles%fn_djac  =cUserDef%modelname(1:i)//'.dj';
	outFiles%fn_model  =cUserDef%modelname(1:i)//'.rho';
	outFiles%fn_residuals  =cUserDef%modelname(1:i)//'.res';
	outFiles%fn_bv	  =cUserDef%modelname(1:i)//'.bv';

	outFiles%fn_avg_cresp='earth_response.cdat';
	outFiles%fn_avg_dresp='earth_response.ddat';


  end subroutine initOutput	! initOutput

end module input
