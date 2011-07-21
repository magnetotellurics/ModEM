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
    read (ioStartup,'(a17,g15.7)') string,cUserDef%damping;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_grid;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_shell;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_field;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_period;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_coords;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_func;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_ctrl;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_invctrl;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_slices;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_param0;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_param;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_source;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_hdata;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_cdata;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_ddata;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_misfit;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_gradient;
    read (ioStartup,'(a17,a80)') string,cUserDef%fn_point;

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
  ! *      line 1: dimensions of model (x,y,z)
  ! *      line 2: x(*) in degrees (interval)
  ! *      line 3: y(*) in degrees (position from n-pole)
  ! *      line 4: z(*) in km (distance from center of the earth, decreasing)

  subroutine initGrid(cUserDef,mygrid)

    type (userdef_control), intent(in)				:: cUserDef
    type (grid_t) , intent(out)					    :: mygrid
    integer				                            :: i

    call read_grid(mygrid,cUserDef%fn_grid)

	grid%nzCrust = grid%nzAir
	inquire(FILE=cUserDef%fn_shell,EXIST=exists)
	! If no thin shell information present, assume no crust for model computations
	if (.not.exists) then
	  write(0,*) node_info,'Warning: No thin shell conductance distribution specified; assume no crust'
	else
	  do i = grid%nzAir+1,grid%nz
		if (clean(grid%r(i)) > CRUST_R + EPS_GRID) then ! Note: all distances are in km
		  grid%nzCrust = grid%nzCrust + 1
		end if
	  end do
	end if

  end subroutine initGrid	! initGrid

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
  ! * initCrust reads the modelfile fn_shell to store S-conductance of Earth's
  ! * crust in GM coordinates. It also stores the depth of the crust in km.
  ! * This should be called after initializing the grid, since it uses the grid
  ! * dimensions to determine the number of entries in file.

  subroutine initCrust(cUserDef,mygrid,mycrust)

    type (userdef_control), intent(in)					:: cUserDef
    type (grid_t) , intent(in)					:: mygrid
    type (modelShell_t) , intent(out)					:: mycrust
    integer				                            :: ios,istat,i


	inquire(FILE=cUserDef%fn_shell,EXIST=exists)
	if(.not.exists) then
	  mycrust%allocated = .false.
	  return
	end if

	allocate(mycrust%cond(mygrid%nx,mygrid%ny),STAT=istat)

	open(ioShell,file=cUserDef%fn_shell,status='old',iostat=ios)

    write(6,*) node_info,'Reading from the thin sheet conductance file ',trim(cUserDef%fn_shell)

	do i=1,mygrid%nx
	  read(ioShell,*,iostat=ios) mycrust%cond(i,1:mygrid%ny)
	end do

	close(ioShell)

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
  ! * initData initializes the data according to the information stored in
  ! * the list of transfer functions TFList; then merges all data types into
  ! * a single data vector allData
  subroutine initData(cUserDef,allData,obsList,freqList,TFList)

	implicit none
    type (userdef_control), intent(in)							:: cUserDef
	type (dataVectorMTX_t), intent(inout)					:: allData
	type (Obs_List), target, intent(inout)					:: obsList
	type (Freq_List), target, intent(inout)					:: freqList
	type (TF_List), target, intent(in)						:: TFList
	! local
	type (dataVectorMTX_t)     :: cdata,ddata
	real(8)                    :: theta
	integer                    :: i,j,k,iTx,iDt,iRx,istat

	! Read the data from files into data blocks, and save in the data vector
    do iDt=1,TFList%n

      select case ( TFList%info(iDt)%name )

      case ('C')
        ! Read C responses into a data vector with a single data type (cdata)
        call initDataList(cdata,iDt,obsList,freqList,cUserDef%fn_cdata)
        if (.not. cdata%allocated) then
          write(0,*) node_info,'Warning: Data misfit will not be calculated for C responses: data file not available'
          call initDataList(cdata,iDt,obsList,freqList)
        end if
        do i=1,cdata%ntx
          do k=1,cdata%d(i)%data(1)%nSite
                iRx = cdata%d(i)%data(1)%rx(k)
                theta = obsList%info(iRx)%colat*d2r;
                if ((theta >= pi/2-EPS_GRID).and.(theta <= pi/2+EPS_GRID)) then
                  write(0,*) node_info,'Error: (initData) receiver ',trim(obsList%info(iRx)%code),' located too close to equator;'
                  write(0,*) node_info,'Error: (initData) please delete this receiver from the list and try again.'
                    exit
                end if
                ! convert from C responses to C ratios
                cdata%d(i)%data(1)%value(:,k) = cdata%d(i)%data(1)%value(:,k)/dtan(theta)
                cdata%d(i)%data(1)%error(:,k) = cdata%d(i)%data(1)%error(:,k)/abs(dtan(theta))
          end do
        end do

      case ('D')
        ! Read D responses into a data vector with a single data type (ddata)
        call initDataList(ddata,iDt,obsList,freqList,cUserDef%fn_ddata)
        if (.not. ddata%allocated) then
          write(0,*) node_info,'Warning: Data misfit will not be calculated for D responses: data file not available'
          call initDataList(ddata,iDt,obsList,freqList)
        end if
        do i=1,ddata%ntx
          do k=1,ddata%d(i)%data(1)%nSite
                iRx = ddata%d(i)%data(1)%rx(k)
                theta = obsList%info(iRx)%colat*d2r;
                if ((theta >= pi/2-EPS_GRID).and.(theta <= pi/2+EPS_GRID)) then
                  write(0,*) node_info,'Error: (initData) receiver ',trim(obsList%info(iRx)%code),' located too close to equator;'
                  write(0,*) node_info,'Error: (initData) please delete this receiver from the list and try again.'
                    exit
                end if
                ! convert from D responses to D ratios
                ddata%d(i)%data(1)%value(:,k) = ddata%d(i)%data(1)%value(:,k)/dsin(theta)
                ddata%d(i)%data(1)%error(:,k) = ddata%d(i)%data(1)%error(:,k)/abs(dsin(theta))
          end do
        end do

      case default

        write(0,*) node_info,'Please specify correct transfer functions in the file ',trim(cUserDef%fn_func)
        stop

      end select

    end do

    ! Now, merge the data types into a single data vector
    call merge_dataVectorMTX(cdata,ddata,allData)

    ! Clean up
    call deall_dataVectorMTX(cdata)
    call deall_dataVectorMTX(ddata)

  end subroutine initData ! initData

  ! ***************************************************************************
  ! * initData reads the file fn_data (in future could be several files) that
  ! * contain all the information about the available data: the values of the
  ! * responses, data errors corresponding to an observatory and a frequency.

  subroutine initDataList(mydat,itype,obsList,freqList,fname)

	implicit none
	type (dataVectorMTX_t), intent(inout)      :: mydat
	integer, intent(in)                        :: itype
    type (Obs_List), target, intent(in)        :: obsList
    type (Freq_List), target, intent(in)       :: freqList
	character(*), intent(in), optional         :: fname
	! local
    integer                 :: ifreq,iobs
    character(100)          :: label
    integer                 :: i,j,ios=0,istat=0
	real(8)                 :: const,large
	real(8)                 :: freq,days,prev
	character(3)            :: code
	real(8)                 :: lon,lat
	real(8)                 :: creal,cimag,cerr
	integer                 :: nfreq,nobs,ncomp,nsite
    integer                 :: countData,countFreq
	complex(8), allocatable :: value(:,:) ! (nfreq,nobs)
	real(8), allocatable    :: error(:,:) ! (nfreq,nobs)
    logical, allocatable    :: exist(:,:) ! (nfreq,nobs)
	logical                 :: new,isComplex,errorBar

	large = 2.0e15
	const = 1.0d-2
	prev = 0.0d0

	! Initialize the data functionals
	nfreq = freqList%n
	nobs  =  obsList%n
    allocate(value(nfreq,nobs),STAT=istat)
    allocate(error(nfreq,nobs),STAT=istat)
    allocate(exist(nfreq,nobs),STAT=istat)
	value(:,:) = dcmplx(0.0d0,0.0d0)
	error(:,:) = large
    exist(:,:) = .FALSE.

    ! If called without a data file, use dictionaries to create data vector
    if (.not. present(fname)) then
	    call create_dataVectorMTX(nfreq,mydat)
	    mydat%allocated = .TRUE.
	    do i = 1,nfreq
	       isComplex = .TRUE.
	       errorBar = .TRUE.
	       call create_dataVector(1,mydat%d(i))
	       mydat%d(i)%tx = i
	       mydat%d(i)%allocated = .TRUE.
	       call create_dataBlock(2,nobs,mydat%d(i)%data(1),isComplex,errorBar)
	       do j = 1,nobs
	           exist(i,j) = .TRUE.
	           if(exist(i,j)) then
	               mydat%d(i)%data(1)%value(1,j) = real(value(i,j))
	               mydat%d(i)%data(1)%value(2,j) = imag(value(i,j))
	               mydat%d(i)%data(1)%error(1,j) = error(i,j)
	               mydat%d(i)%data(1)%error(2,j) = error(i,j)
                   mydat%d(i)%data(1)%exist(1,j) = exist(i,j)
                   mydat%d(i)%data(1)%exist(2,j) = exist(i,j)
	               mydat%d(i)%data(1)%rx(j) = j
	           end if
	       end do
	       mydat%d(i)%data(1)%dataType = itype
	       mydat%d(i)%data(1)%tx = i
	       mydat%d(i)%data(1)%allocated = .TRUE.
	    end do
	    deallocate(value,STAT=istat)
	    deallocate(error,STAT=istat)
	    deallocate(exist,STAT=istat)
	    return
    end if

    ! Otherwise, if data file does not exist, exit
	inquire(FILE=trim(fname),EXIST=exists)
	if (.not.exists) then
      deallocate(value,STAT=istat)
      deallocate(error,STAT=istat)
      deallocate(exist,STAT=istat)
	  return
	end if

	! If data file is present, initialize data, functionals and residuals
    open(ioDat,file=trim(fname),status='old',form='formatted',iostat=ios)

    write(6,*) node_info,'Reading from the data file ',trim(fname)
    read(ioDat,'(a)') label
    write(6,*) label
    read(ioDat,'(a)') label

    countData = 0
    countFreq = 0


	READ_DATA: do

	  !read(ioDat,'(a12,1x)',advance='no',iostat=ios)  label !'# freq, obs: '
	  !read(ioDat,*,iostat=ios) freq, num

	  if (ios /= 0) exit

	  read(ioDat,*,iostat=ios) days,code,lon,lat,creal,cimag,cerr

	  ! Find the relevant frequency in the list
	  freq  = 1/(days*24*60*60)
	  new = .FALSE.
	  if (abs((prev-freq)/freq) > const) then
		new = .TRUE.
		prev = freq
	  end if
	  if (new) then
	    ifreq = getFreq(freqList,freq)
	    if (ifreq > 0) then
	       write(6,'(a12,a32,i6,a2,es12.6,a5)') node_info,'Reading data for the period ',ifreq,': ',days,' days'!freqList%info(i)%value
	       countFreq = countFreq + 1
		end if
	  end if
      if (ifreq == 0) then
        !print *, 'Frequency ',freq,' is not in the frequency list; ignore'
        !do j=1,num
          !read(ioDat,*,iostat=ios)
        !end do
        cycle READ_DATA
      end if
	  ! Find the relevant observatory in the list
	  iobs = getObs(obsList,code)
	  if (iobs > 0) then
	    countData = countData + 1
	  else
		!print *, 'Observatory ',trim(code),' in file ',trim(fname),' is not in our list'
		!print *, 'This is an error; exiting...'
		!stop
		cycle READ_DATA
	  end if

	  ! If everything is correct, save the data into the 2D arrays
	  value(ifreq,iobs) = km2m * dcmplx(creal,cimag)
	  error(ifreq,iobs) = km2m * cerr
	  exist(ifreq,iobs) = .TRUE.

	end do READ_DATA

	close(ioDat)

	if (countData==0) then
  	  write(0,*) node_info,'Warning: No data matches given frequencies and observatories'
	else
  	  write(0,*) node_info,'Number of complex data values in total: ',countData
	end if

	! Finally, store the data in the data vector
	call create_dataVectorMTX(countFreq,mydat)
	mydat%allocated = .TRUE.
	i = 0
	ifreq = 0
	SAVE_DATA: do ifreq = 1,nfreq

	   nComp = 2
	   nSite = count(exist(ifreq,:))
	   isComplex = .TRUE.
	   errorBar = .TRUE.
	   if(nSite == 0) then
	       ! no data for this frequency ... try next one
	       cycle SAVE_DATA
	   end if
       i = i + 1
	   call create_dataVector(1,mydat%d(i))
	   mydat%d(i)%tx = ifreq
       mydat%d(i)%allocated = .TRUE.
	   call create_dataBlock(nComp,nSite,mydat%d(i)%data(1),isComplex,errorBar)
	   j = 1
	   do iobs = 1,nobs
	       if(exist(ifreq,iobs)) then
	           mydat%d(i)%data(1)%value(1,j) = real(value(ifreq,iobs))
               mydat%d(i)%data(1)%value(2,j) = imag(value(ifreq,iobs))
               mydat%d(i)%data(1)%error(1,j) = error(ifreq,iobs)
               mydat%d(i)%data(1)%error(2,j) = error(ifreq,iobs)
               mydat%d(i)%data(1)%exist(1,j) = exist(ifreq,iobs)
               mydat%d(i)%data(1)%exist(2,j) = exist(ifreq,iobs)
               mydat%d(i)%data(1)%rx(j) = iobs
               j = j + 1
           end if
       end do
       mydat%d(i)%data(1)%dataType = itype
       mydat%d(i)%data(1)%tx = ifreq
       mydat%d(i)%data(1)%allocated = .TRUE.

    end do SAVE_DATA

    deallocate(value,STAT=istat)
    deallocate(error,STAT=istat)
    deallocate(exist,STAT=istat)

    return

  end subroutine initDataList !	initData


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

  subroutine initModelParam(cUserDef,myparam,p0)

	use model_operators

    type (userdef_control), intent(in)					:: cUserDef
    type (modelParam_t), intent(inout)					:: myparam
	logical, intent(in), optional		:: p0

	if(present(p0)) then
	  if(p0) then
	    call read_modelParam(myparam,cUserDef%fn_param0)
	  end if
	else
      call read_modelParam(myparam,cUserDef%fn_param)
	end if

  end subroutine initModelParam	! initModelParam


  ! ***************************************************************************
  ! * initControls reads the file fn_ctrl and sets the values of main control
  ! * parameters for the forward solver, stored in fwdCtrls

  subroutine initControls(cUserDef,fwdCtrls)

	implicit none
	type (userdef_control), intent(in)						:: cUserDef
	type (fwdCtrl_t), intent(out)						:: fwdCtrls

	  open(ioFwdCtrl,file=cUserDef%fn_ctrl,form='formatted',status='old')

      write(6,*) node_info,'Reading from the forward solver controls file ',trim(cUserDef%fn_ctrl)
      read(ioFwdCtrl,*) fwdCtrls%ipotloopmax
      read(ioFwdCtrl,*) fwdCtrls%errend
      read(ioFwdCtrl,*) fwdCtrls%nrelmax
      read(ioFwdCtrl,*) fwdCtrls%n_reldivh
      read(ioFwdCtrl,*) fwdCtrls%ipot0,fwdCtrls%ipotint,fwdCtrls%ipot_max

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
	outFiles%fn_cjac  =cUserDef%modelname(1:i)//'.cj';
	outFiles%fn_djac  =cUserDef%modelname(1:i)//'.dj';
	outFiles%fn_model  =cUserDef%modelname(1:i)//'.rho';
	outFiles%fn_residuals  =cUserDef%modelname(1:i)//'.res';
	outFiles%fn_bv	  =cUserDef%modelname(1:i)//'.bv';

	outFiles%fn_avg_cresp='earth_response.cdat';
	outFiles%fn_avg_dresp='earth_response.ddat';


  end subroutine initOutput	! initOutput

end module input
