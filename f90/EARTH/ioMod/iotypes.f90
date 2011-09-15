! *****************************************************************************
module iotypes
  ! This module defines the derived data type structure with all filenames
  ! and any other types necessary for streamlined output

  use math_constants
  use file_units
  implicit none

  public            :: initSlices, deall_Slices

  ! ***************************************************************************
  ! * storing the forward solver control parameters
  type :: fwdCtrl_t

      integer                                   :: ipotloopmax
      integer                                   :: nrelmax
      integer                                   :: n_reldivh
      integer                                   :: ipot0
      integer                                   :: ipotint
      integer                                   :: ipot_max
      real(8)                                   :: errend

  end type fwdCtrl_t

  ! ***************************************************************************
  ! * userdef_control contains the list of all essential input information currently
  ! * read in from fn_startup.
  type :: userdef_control

	character(80)				:: paramname  ! 'harmonic'/'mixed'/'grid' parametrization
	character(80)				:: modelname  ! specify how to call the output files
	character(80)				:: verbose  ! 'none'/'model'/'responses'/'all'/'debug'
	character(80)				:: calculate  ! 'responses'/'jacobian'/'derivative'
    character(80)               :: secondary_field  ! 'yes'/'no'
	character(80)				:: fn_thinsheet	! GM thin sheet conductance values
	character(80)				:: fn_grid	! grid information
	character(80)				:: fn_param0	! base model parametrization
	character(80)				:: fn_param	! information about the parametrization
	character(80)				:: fn_period  ! periods or frequencies
    character(80)               :: fn_extsource ! external source
    character(80)               :: fn_intsource ! additional interior source (not yet implemented)
	character(80)				:: fn_fwdctrl	! forward solver control
    character(80)               :: fn_adjctrl   ! adjoint solver control
	character(80)				:: fn_invctrl	! inverse solver control
	character(80)				:: fn_slices  ! grid radii at which we output the data
	character(80)				:: fn_coords  ! observatory coordinates file
	character(80)				:: fn_func  ! information about data functionals
	character(80)				:: fn_cdata	! name of C responses data file
	character(80)				:: fn_ddata	! name of D responses data file
    character(80)               :: fn_hdata ! name of magnetic field data file
    character(80)               :: fn_href ! reference observatories for magnetic field inversion

	! obsolete
	character(80)				:: fn_field	! input radial field solution (obsolete: computed internally)
	character(80)				:: fn_precond	! preconditioning parameters
	character(80)				:: fn_misfit  ! output file for data misfit
	character(80)				:: fn_gradient	! output file for derivative
	character(80)				:: fn_point	! output file for point parametrization
	real(8)             		:: damping ! the value of damping parameter mu
	real(8)             		:: step_size ! initial step size parameter for inversion

  end type userdef_control ! userdef_control


  ! ***************************************************************************
  ! * output_info contains the list of all output file names. Files such as the
  ! * boundary condition or starting solution, or diagnostic output are not
  ! * currently required. Will be added if necessary.
  type :: output_info

	! Output files:
	character(80)				:: fn_hx, fn_hy, fn_hz
	character(80)				:: fn_jxjyjz, fn_hxhyhz
	character(80)				:: fn_err
	character(80)				:: fn_cresp, fn_dresp
	character(80)				:: fn_cdat, fn_ddat, fn_hdat
	character(80)				:: fn_cjac, fn_djac
	character(80)				:: fn_model
	character(80)				:: fn_residuals
	character(80)				:: fn_avg_cresp, fn_avg_dresp
	! Optional input file:
    character(80)				:: fn_bv

  end type output_info  ! output_info

  logical                                       :: exists ! for I/O inquiries
  character(100)                                :: label  ! first line in files


  ! ***************************************************************************
  ! * list of radii ("slices") at which we output the full solution ...
  ! * ... this doesn't necessarily belong here, but good enough for now!
  type :: Rad_List

    character(80)                               :: type ! 'CELL'/'NODE'
    integer                                     :: n
    real(8), pointer, dimension(:)              :: r    !nrad

  end type Rad_List

  ! * user-defined list of radii in meters at which we output the solution
  type (Rad_List), save                             :: slices


Contains

  ! *****************************************************************************
  ! * initSlices reads the file fn_slices that contains the information about the
  ! * number and the values of grid radii at which we output full solution
  subroutine initSlices(cUserDef,slices)

    implicit none
    type (userdef_control), intent(in)                      :: cUserDef
    type (Rad_List), intent(out)                            :: slices
    character(80)                                           :: type
    integer                                                 :: num,i,ios

    inquire(FILE=cUserDef%fn_slices,EXIST=exists)
    ! If file is present, initialize the output radii ("slices")
    if (.not.exists) then
      write(0,*) 'Warning: No radii specified; we will only output surface solution at cell centers'
      allocate(slices%r(1))
      slices%type = 'CELL'
      slices%n = 1
      slices%r(1) = EARTH_R
    end if

    open(ioRad,file=cUserDef%fn_slices,status='old',form='formatted',iostat=ios)

    write(6,*) node_info,'Reading from the output radii file ',trim(cUserDef%fn_slices)
    read(ioRad,'(a)') label

    read(ioRad,*) type
    if ((trim(type) .ne. 'CELL') .and. (trim(type) .ne. 'NODE')) then
      write(6,*) 'Error: Unknown output solution type ',trim(type)
      stop
    end if

    read(ioRad,*) num
    allocate(slices%r(num))
    do i=1,num

      read(ioRad,*) slices%r(i)

    end do

    close(ioRad)

    slices%n = num
    slices%type = type

  end subroutine initSlices ! initSlices


  ! **************************************************************************
  ! Cleans up and deletes slices dictionary at end of program execution
  subroutine deall_Slices()

    integer     :: i,istat

    if (associated(slices%r)) then
        deallocate(slices%r,STAT=istat)
    end if

  end subroutine deall_Slices

end module iotypes
