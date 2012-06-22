! *****************************************************************************
module ForwardSolver
  !  High level interface/control module used by top level routines
  !   for initializing and using the solver.  The key public routines
  !   in this module have only generic (abstract data type) arguments
  !   and can thus be called from the top-level inversion routines.
  !
  !  All public routines must have the same generic names, parameter lists
  !   and abstract functionality; this implementation is for EARTH
  !
  !  Idea is to always call initSolver before calling fwdSolver,
  !   in particular before the first solution for each transmitter
  !   (frequency).  If called for the first time (in a program run,
  !   or after a call to exitSolver), full initialization
  !   (after deallocation/cleanup if required) is performed.

  use math_constants
  use SolnSpace
  use jacobian
  use transmitters
  use input
  use output
  use UserData
  use initFields
  use modelmap
  use field1d

  implicit none

  public :: initSolver, fwdSolve, sensSolve, exitSolver

  save
  private
  type(timer_t)     :: timer
  type(rscalar), target     :: rho0 ! fixed background resistivity on the grid
  logical           :: solverInitialized = .false.
  logical           :: secondaryField = .false.
  logical           :: newModelParam = .false.
  type(cvector)     :: b, e
  type(rvector)     :: drhoF
  type(rscalar)     :: rho1d, drho
  type(modelParam_t):: m1d
  type(rscalar)     :: rho ! resistivity on the grid
  type(cvector)     :: source ! user-specified interior source
  type(sparsevecc)  :: BC ! boundary conditions set from P10
  type(modelParam_t)    :: mBackground  ! use it to get the radial 1d model for SFF
  type(modelParam_t)    :: mPrev  ! store the previous model for efficiency
  type(solnVector_t)    :: hPrev  ! store the previous forward solver solution
  !type(rhsVector_t)     :: b   ! needed to store the RHS for forward solver
  type(solnVector_t) :: h1d, dh ! needed for secondary field formulation
  type(modelShell_t) :: crust
  ! debugging variables
  character(80)     :: cfile

Contains

  !**********************************************************************
  subroutine initSolver(iTx,m0,grid,h0,h,comb)
   !  Initializes both forward and sens solvers for transmitter iTx.

   integer, intent(in)                      :: iTx
   type(modelParam_t),intent(in)	        :: m0
   type(grid_t), intent(in), target         :: grid
   type(solnVector_t), intent(inout), optional :: h0
   type(solnVector_t), intent(inout), optional :: h
   type(rhsVector_t), intent(inout), optional :: comb
   !  local variables
   type(transmitter_t), pointer             :: freq
   type(modelParam_t)                       :: mgrid
   logical                                  :: initFwd, initForSens

   initFwd = present(h0)
   initForSens = present(comb)

   ! First of all, initialize a local copy of the GRID.
   ! All vectors will point to that grid, which should be
   ! local to each processor.
   !grid = igrid

   freq => freqList%info(iTx)
   write(*,'(a12,a46,es9.3,a5)') node_info, &
        'Initializing 3D SGFD global solver for period ',freq%period,' days'

   ! reset timer
   call reset_time(timer)

   ! If h0 is already allocated, do not reinitialize - use the previous
   ! forward solution as starting solution for this frequency ...
   if(initFwd) then
     if(.not. h0%allocated) then
        call create_solnVector(grid,iTx,h0)
        call initialize_fields(grid,h0%vec,BC)
     endif
   endif

   if(initForSens) then
        ! also initialize the optional outputs
      call create_solnVector(grid,iTx,h)
      comb%nonzero_source = .true.
      call create_rhsVector(grid,iTx,comb)
   endif

   secondaryField = freq%secondaryField

   if(.not. solverInitialized) then
      ! only do this when called for the first time; rho0 is background resistivity
      ! that is already initialized and saved in m0 ... usually 1D. However, recomputing
      ! rho1d from m1d allows rho0 to have more general form. Also note that setting
      ! rho1d = rho0 will not work since rho0 has the thinsheet in it already.
      rho0 = m0%rho0
      if(secondaryField) then
	    ! Make 1D parametrization out of full, and map to grid (primary cell centers)
	    write(0,*) node_info,'Reading the background resistivity model for SFF; assuming it is log10 and 1D.'
	    call read_modelParam(mBackground,cUserDef%fn_param0)
	    call setGrid_modelParam(grid,mBackground)
        call initCrust(cUserDef,grid,crust)
        call setCrust_modelParam(crust,mBackground)
	    call getRadial(m1d,mBackground)
        write(0,*) node_info,'Crust average conductance is ',m1d%crust%avg
	    call create_rscalar(grid,rho1d,CENTER)
	    call initModel(m1d,rho1d,grid) ! do NOT use background resistivity for this
	    !call mapToGrid(mgrid,m1d)
	    rho1d = mgrid%rho
        call create_solnVector(grid,iTx,h1d)
        !call fwdSolve1d(iTx,m1d,h1d) ! should be saved in solnVectorMTX if computed once
        !call outputModel('test1d.rho',grid,rho1d%v)
	  endif
      solverInitialized = .true.
   endif

   !write(0,*) 'Output background model...'
   !call outputModel('background.rho',grid,rho0%v)

   ! compute the 1D fields; could be computed once for all frequencies & saved in solnVectorMTX_t
   if(secondaryField) then
!      ! instead of reading the fields, compute them
!      !call read_solnVector(cUserDef%fn_field,grid,iTx,h1d)
!      !call read_solnVector('p10_5deg',grid,iTx,h1d)
      call create_solnVector(grid,iTx,h1d)
      call fwdSolve1d(iTx,m1d,h1d)
!
!      ! DEBUG: output forward solver forcing
!        write(cfile,'(a11,i3.3,a6)') 'compute_1D_',iTx,'.field'
!        write(*,*) 'Writing to file: ',cfile
!        open(ioWRITE,file=cfile,status='unknown',form='formatted')
!        write(ioWRITE,'(a45,f9.3,a6)') "# Full EM field solution output"
!        write(ioWRITE,'(i3)') 1
!        call write_cvector(ioWRITE,h1d%vec)
!        close(ioWRITE)
   endif

   ! this will only be true when new model is supplied (m0 /= mPrev)
   mPrev = m0
   newModelParam = .true.

   if(newModelParam) then
      ! compute the resistivity on the grid
      write(0,*) node_info,'Mapping new model parameter to grid'
      !call mapToGrid(mgrid,m0)
      !rho = mgrid%rho
      !call create_rscalar(grid,rho,CENTER)
      !write(0,*) 'Output background model...'
      !call outputModel('background.rho',grid,rho0%v)
      !if (iTx == 1) then
      !  call write_modelParam(m0,'newModelParam.prm')
      !end if
      call initModel(m0,rho,grid,rho0)

      if(secondaryField) then
	    ! Take the difference on the grid, to avoid the problem with zero resistivity
	    call create_rscalar(grid,drho,CENTER)
	    call linComb_rscalar(ONE,rho,MinusONE,rho1d,drho)
	    !call outputModel('testdelta.rho',grid,drho%v)

	    ! Map the resistivity vector to primary cell faces (dual edges)
	    call operatorL(drho,drhoF,grid)
	  endif

   endif

   ! Initialize interior source from file (currently, set to zero)
   if(.not. source%allocated) then
    call create_cvector(grid,source,EDGE)
   endif

   call deall_modelParam(mgrid)

  end subroutine initSolver

  !**********************************************************************
  subroutine fwdSolve1d(iTx,m1d,h1d)
!  1D Forward solver. External source is stored in the transmitter.
!  Output solution vector has to be pre-initialized with the grid.

   integer, intent(in)                          :: iTx
   type(modelParam_t), intent(in)               :: m1d
   type(solnVector_t), intent(inout)            :: h1d
   ! local variables
   type(conf1d_t)                               :: conf1d
   type(transmitter_t), pointer                 :: freq
   real(kind=prec)                              :: tsdepth,tssigma,tau,period ! secs
   real(kind=prec), allocatable, dimension(:)   :: depths,logrho
   integer                                      :: i,nL,lmax,istat

   ! IMPORTANT: FIRST update pointer to the transmitter in solnVector
   h1d%tx = iTx

   freq => freqList%info(iTx)

   ! run FWD solver
   write (*,'(a12,a12,a3,a20,i4,a2,es12.6,a5)') node_info, &
    'Solving the 1D ',FWD,' problem for period ',iTx,': ',1/freq%value,' secs'

    period  = 1./freq%value     ! period in secs

    nL = m1d%nL
    allocate(depths(nL),logrho(nL), STAT=istat)
    do i = 1,nL
        depths(i) = m1d%L(i)%depth
        call getCoeffValue_modelParam(m1d,i,0,0,logrho(i))
    end do

    ! source file should only have one layer
    lmax = freq%degree

    ! set earth radius and domain top radius (in meters)
    conf1d%r0  = 6371.0e3
    conf1d%rmax= 1.0e3 * h1d%grid%r(1)

    ! set tolerance on toroidal potential
    conf1d%tol = 1.e-9

    ! if thinsheet exists, replace it with an infinitely thin layer
    if (h1d%grid%nzCrust > 0) then
        tsdepth = h1d%grid%z(h1d%grid%nzAir+1) - h1d%grid%z(h1d%grid%nzAir+h1d%grid%nzCrust+1)
        tssigma = 10.0**( -logrho(1) )
    end if
    tau = m1d%crust%avg - tsdepth * tssigma

    ! ... and set surface conductance
     write(0,'(a12,a26,f12.3,a19)') node_info,'Using surface conductance ',1.0d0,' S for 1D modelling'
     conf1d%tau = 1.0d0
!    if (h1d%grid%nzCrust > 0) then
!        write(0,'(a12,a26,f12.3,a19)') node_info,'Using surface conductance ',tau,' S for 1D modelling'
!        conf1d%tau = tau
!    else
!        write(0,'(a12,a71)') node_info,'No thinsheet; using negligible surface conductance 1 S for 1D modelling'
!        conf1d%tau = 1.0d0
!    end if

    ! save model in 1D configuration structure: layers include the core
    allocate(conf1d%layer(nL+1),conf1d%sigma(nL+1), STAT=istat)
    conf1d%layer(1) = 0.0d0
    conf1d%layer(2:nL+1) = 1.0e3 * depths(1:nL)
    conf1d%sigma(1:nL) = 10.0**( -logrho(1:nL) )
    conf1d%sigma(nL+1) = 10.0**( 5.0 ) ! core conductivity

    write(*,*) 'Tops of model layers: ',conf1d%layer
    write(*,*) 'Conductivity values:  ',conf1d%sigma

    ! compute full magnetic field for the layered model
    call sourceField1d(conf1d,lmax,freq%jExt,period,h1d%grid,h1d%vec)

    deallocate(depths,logrho, STAT=istat)
    deallocate(conf1d%layer,conf1d%sigma, STAT=istat)

   if (output_level > 1) then
      write (*,*) node_info, ' time taken (mins) ', elapsed_time(timer)/60.0
   end if
  end subroutine fwdSolve1d

  !**********************************************************************
  subroutine fwdSolve(iTx,h)
!  Forward solver; uses b as the RHS. Similar to sensSolve, less general.
!  Temporary routine (will be merged with sensSolve).

   integer, intent(in)                        :: iTx
   type(solnVector_t), intent(inout)            :: h
   ! local variables
   type(transmitter_t), pointer                 :: freq
   real(kind=prec)                              :: omega
   logical                                      :: adjoint,sens

   ! IMPORTANT: FIRST update pointer to the transmitter in solnVector
   h%tx = iTx

   freq => freqList%info(iTx)

   ! run FWD/ADJ solver
   write (*,'(a12,a12,a3,a20,i4,a2,es12.6,a5)') node_info, &
    'Solving the ',FWD,' problem for period ',iTx,': ',1/freq%value,' secs'

   omega  = 2.0d0*pi*freq%value     ! angular frequency (radians/sec)

   if (secondaryField) then

      write(*,*) node_info, 'Using the secondary field formulation ...'

      ! Compute the RHS = - del x drho (del x H)
      h = h1d
      call operatorD_l_mult(h%vec,h%grid)
      call operatorC(h%vec,e,h%grid)
      call diagMult(drhoF,e,e)
      call operatorCt(e,b,h%grid)
      call operatorD_Si_divide(b,h%grid)
      call linComb(C_MinusONE,b,C_ZERO,b,b)

      ! solve S_m <h> = <b> for vector <h>
      dh = h
      call zero_solnVector(dh)
      adjoint = .false.
      sens = .true.
      call linComb_cvector(C_ONE,source,C_ONE,b,b)
      call operatorMii(dh%vec,b,omega,rho,h%grid,fwdCtrls,h%errflag,adjoint)

      ! Full solution for one frequency is the sum H1D + dH
      call linComb_solnVector(C_ONE,h1d,C_ONE,dh,h)

   else

      write(*,*) node_info, 'Using the forward solver ...'

      ! solve S_m <h> = <s> for vector <h>
      adjoint = .false.
      sens = .true.
      call operatorMii(h%vec,source,omega,rho,h%grid,fwdCtrls,h%errflag,adjoint,BC)

   end if

   ! compute and output fields & C and D responses at cells
   call outputSolution(freq,h%vec,slices,h%grid,cUserDef,rho%v,'h')

   ! output full H-field cvector
   if (output_level > 3) then
      call write_solnVector(cUserDef%modelname,h)
   end if

   if (output_level > 1) then
      write (*,*) node_info, ' time taken (mins) ', elapsed_time(timer)/60.0
   end if
  end subroutine fwdSolve

  !**********************************************************************
  subroutine sensSolve(iTx,FWDorADJ,h,comb)
!  Generic forward solver; if comb is specified, it is used as forcing;
!  otherwise, assume that the RHS is already stored in b.
!  Use the input e for starting solution
!  (normally, e will be zero on input for sensitivities).

   integer, intent(in)                        :: iTx
   character (len=3), intent(in)                :: FWDorADJ
   type(solnVector_t), intent(inout)            :: h
   type(rhsVector_t), intent(inout), optional    :: comb
   ! local variables
   type(transmitter_t), pointer                 :: freq
   real(kind=prec)                              :: omega
   logical                                      :: adjoint,sens

   ! update pointer to the transmitter in solnVector
   h%tx = iTx

   freq => freqList%info(iTx)

   ! run FWD/ADJ solver
   write (*,'(a12,a12,a3,a20,i4,a2,es12.6,a5)') node_info, &
    'Solving the ',FWDorADJ,' problem for period ',iTx,': ',1/freq%value,' secs'

   omega  = 2.0d0*pi*freq%value     ! angular frequency (radians/sec)

   ! solve S_m <h> = <s> for vector <h>
   if(.not. present(comb)) then

    ! assume interior forcing s + BC; starting solution already initialized
    adjoint = .false.
    sens = .false.
    call operatorMii(h%vec,source,omega,rho,h%grid,fwdCtrls,h%errflag,adjoint,BC)

   else

    ! use comb for forcing and assume zero BC; starting solution should be zero
    adjoint = (FWDorADJ .ne. FWD)
    sens = .true.
    call operatorMii(h%vec,comb%source,omega,rho,h%grid,adjCtrls,h%errflag,adjoint)

   endif

   if (output_level > 1) then
      write (*,*) node_info, ' time taken (mins) ', elapsed_time(timer)/60.0
   end if


  end subroutine sensSolve


  !**********************************************************************
  subroutine exitSolver(h0,h,comb)
   ! Cleans up after fwdSolver and deallocates all solver data. Call
   ! initSolver to update solver data; call this to exit completely.
   ! Optionally, deallocates e0,e,comb
   type(solnVector_t), intent(inout), optional  :: h0
   type(solnVector_t), intent(inout), optional  :: h
   type(rhsVector_t), intent(inout), optional   :: comb

   ! local variables
   logical          :: initForSens

   initForSens = present(comb)

   if(present(h0)) then
      call deall_solnVector(h0)
   endif

   if(initForSens) then
      call deall_rhsVector(comb)
      call deall_solnVector(h)
   endif

   if(solverInitialized) then
      ! cleanup/deallocation routines for model operators
      if(secondaryField) then
        call deall_solnVector(h1d)
        call deall_solnVector(dh)
        call deall_cvector(b)
        call deall_cvector(e)
        call deall_rvector(drhoF)
        call deall_rscalar(rho1d)
        call deall_rscalar(drho)
        call deall_modelParam(m1d)
      endif
      !call deall_rhsVector(b)
      call deall_solnVector(hPrev)
      call deall_modelParam(mPrev)
      call deall_modelParam(mBackground)
      call deall_cvector(source)
      call deall_sparsevecc(BC)
      call deall_rscalar(rho)
      call deall_rscalar(rho0)
      solverInitialized = .false.
      newModelParam = .false.
   endif

   ! restart the clock
   call reset_time(timer)

  end subroutine exitSolver


end module ForwardSolver
