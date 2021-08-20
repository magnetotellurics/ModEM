module ForwardSolver

!  High level interface/control module used by top level routines
!   for initializing and using the solver.  The key public routines
!   in this module have only generic (abstract data type) arguments
!   and can thus be called from the top-level inversion routines.
!  A similar solver interface will be required for any specific use
!   of the top-level inversion modules
!  All public routines must have the same generic names, parameter lists
!   and abstract functionality; this implementation is for 3D-MT
!
!  Main difference from 2D: no need to keep track of TE/TM modes

use math_constants
use datafunc
use dataspace
use solnspace
use emsolve3d
use transmitters
!=======================================================================
!Mar. 13, 2011================ Special Module added for CSEM calculation 
!=======================================================================
! Dipole1D module is module for calculating primary field.
! This module was created by Key 2009
   use Dipole1D
   use DC_GeoElec   

implicit none

! options for boundary conditions input or computation
logical, save, public   :: COMPUTE_BC = .true.
logical, save, public   :: NESTED_BC = .false.
logical, save, public   :: BC_FROM_RHS_FILE = .false.
logical, save, public   :: BC_FROM_E0_FILE = .false.

! option to read a primary solution from file for SFF
logical, save, public   :: PRIMARY_E_FROM_FILE = .false.


!=======================================================================
!Mar. 13, 2011============================== Use only for MT Calculation
!=======================================================================

!! Used for interpolating the BC from a larger grid.
  type(grid_t),save,public               	 ::	 Larg_Grid   
  type(solnVectorMTX_t),save,public          ::  eAll_larg
  integer              ,save,public          ::   nTx_nPol
  logical              ,save,public          ::  nestedEM_initialized

!=======================================================================
!May 15, 2018== AK == New general RHS for all transmitters now stored
!=======================================================================
  type(rhsVectorMTX_t),save,public           ::  bAll

!=======================================================================
!Aug 18, 2021== AK == Also store an array of primary fields (not good
!in terms of memory usage, will read each from file as needed ...)
!=======================================================================
  type(solnVectorMTX_t),save,public           ::  eAllPrimary
  type(modelParam_t),save,public              ::  sigmaPrimary

!  initialization routines (call Fwd version if no sensitivities are
!     are calculated).  Note that these routines are set up to
!    automatically manage memory and to figure out which initialization
!    (or reinitialization) steps are required (e.g., when the frequency
!    changes from the previous solver call, appropriate solver
!    coefficients are updated, matrices factored, etc.).  This
!    functionality needs to be maintained in implementations for new
!    problems!

!=======================================================================
!Mar. 13, 2011============== Special Variable added for CSEM calculation 
!=======================================================================
   type(rvector), save, private :: condAnomaly ! Anomalous conductivity (on edge)
   type(cvector), save, private :: E_p,V_p         ! Primary field
   type(rvector), save :: condNomaly ! Nomalous conductivity (on edge)
      
   type (rscalar),save,public   	:: phi 
 
!   real(kind=prec), save, private :: nlay1D ! [S/m]
!   real(kind=prec), save, private allocatable, dimension(:) :: zlay1D ! [m]
!   real(kind=prec), save, private allocatable, dimension(:) :: sig1D ! [S/m]


public initSolver

!  cleanup/deallocation routines
public exitSolver

! solver routines
public fwdSolve, sensSolve

logical, save, private		:: modelDataInitialized = .false.
logical, save, private		:: BC_from_file_Initialized = .false.
!  logical, save, private		:: sigmaNotCurrent = .true.

Contains
   !**********************************************************************
   ! edited by AK 24 May 2018 to store the output in bAll. Strongly
   ! suspect that BC_from_file can be completely replaced by using bAll
   ! but for now, trying to minimize changes to the extent possible...
   ! cleaning up for nested modeling is not complete yet.
   ! Note that as coded by NM, it looks like the BC from large E-solution
   ! are the same for ALL periods and modes. This is actually not true.
   ! Here, BC were merely input used to create BC_from_file in nestedEM.f90.
   ! I now initialize them there.
subroutine Interpolate_BC(grid)
  type(grid_t), intent(in)        :: grid

  ! local
  integer       :: iTx,ix,iy,iz
  

  !In case of interpolating the BC from eAll_larg   
  ! If eAll_ larg solution is already allocated, then use that to interpolate the BC from it
     
  if (eAll_larg%allocated) then
  	write(15,*) ' Start interpolating',grid%nx,grid%ny,grid%nz
	call Interpolate_BC_from_E_soln (eAll_larg,Larg_Grid,Grid)
        !Once we are ready from eAll_larg, deallocate it, and keep track, that BC_from_file are already Initialized.
        call deall(eAll_larg)
        call deall_grid(Larg_Grid)
        BC_from_file_Initialized=.true.
  end if
  write(15,*) ' End interpolating',BC_from_file(1)%yXMin(10,11)
  


        
      
end subroutine Interpolate_BC
!**********************************************************************
subroutine ini_BC_from_file(grid)
  type(grid_t), intent(in)        :: grid

  ! local
  type(cboundary)   :: BC
  integer           :: iTx,j,status
  

     call create_cboundary(grid,BC)
     
     
     allocate (BC_from_file(nTx_nPol), STAT=status)
     do j=1,nTx_nPol
       BC_from_file(j)=BC
     end do
     BC_from_file_Initialized=.true.
     

end subroutine ini_BC_from_file     

!**********************************************************************
subroutine copyE0fromFile()
   !   this is just another kluge --- eAll_larg is not available to SetBound,
   !         a routine in ModelOperator3D; copy to E0_from_file (in NestedEM) which is
   !   now have moved the logic out of ModelOperator3D but keeping this temporarily
   !   for historic reasons, until we can revisit [AK]
  integer   :: counter,j,k,status
     allocate (E0_from_file(nTx_nPol), STAT=status)
     counter = 0
     do j=1,eAll_larg%nTx
        !  hard to imagine anything but 1 polarization per transmitter here!
        do k = 1,eAll_larg%solns(1)%nPol
           counter = counter+1
           E0_from_file(counter)=eAll_larg%solns(j)%pol(k)
        end do
     end do

end subroutine copyE0fromFile

   !**********************************************************************
   subroutine initSolver(iTx,sigma,grid,e0,e,comb)
   !   Initializes forward solver for transmitter iTx.
   !     Idea is to call this before calling fwdSolve or sensSolve,
   !     in particular before the first solution for each transmitter
   !     (frequency).  If called for the first time (in a program run,
   !     or after a call to exitSolver), full initialization
   !     (after deallocation/cleanup if required) is performed.
   !
   !   iTx defines transmitter: for 2D MT, this provides info about
   !       frequency and TE/TM mode; for 3D frequency and number
   !       of polarizations
   !
   !   This now does all setup (including matrix factorization) for
   !     the appropriate mode/frequency
   !   NOTE: e and comb are optional calling arguments;
   !     both should be present if one is

   integer, intent(in)				:: iTx
   type(modelParam_t),intent(in), target		:: sigma
   type(grid_t), intent(in), target         :: grid
   !  following structures are initialized
   !	solution vector for forward problem
   type(solnVector_t), intent(inout)			:: e0
   !	solution vector for sensitivity
   type(solnVector_t), intent(inout), optional	:: e
   !	forcing for sensitivity
   type(rhsVector_t), intent(inout), optional		:: comb

   !  local variables
   integer		:: IER,k,ix,iy,iz,counter,nzAir
   character*80 :: gridType,paramtype
   logical		:: initForSens,sigmaNotCurrent

   type (modelParam_t) :: sigmaTemp
   type(timer_t) :: timeDipole
   real :: timeD
   

	

   initForSens = present(comb)

   !  allocate for background solution
   call create_solnVector(grid,iTx,e0)

   if(initForSens) then
      !  allocate for sensitivity solution, RHS - same for all TX types
      !  assuming here that we don't use sparse storage ... we could!
      call create_solnVector(grid,iTx,e)
      comb%nonzero_source = .true.
      comb%sparse_source = .false.
      comb%nonzero_bc = .false.
      call create_rhsVector(grid,iTx,comb)
!      do k = 1,comb%nPol
!        comb%b(k)%sparse_Source = .false.
!        comb%b(k)%adj = ''
!        !  using all this information, reallocate storage for each polarization
!        call create_RHS(grid,iTx,comb%b(k))
!      enddo
   endif


   if(.NOT.modelDataInitialized) then
   !   Initialize modelData, setup model operators
      call ModelDataInit(grid)
      call ModelOperatorSetup()
      modelDataInitialized = .true.
	  
	!=======================================================================
      !Mar. 13, 2011================================= Use for CSEM Calculation
      !=======================================================================
	  !   CALL setPrimaryCond(sigma, .true.) ! This subroutine created by Aihua
		    ! This subroutine can create the primary model by averaging Sigma
			! or read it from file. If the second argument is true, this routine
			! will read Primary model from file.
			
	  !====================================================================
	  !===================== Create 1D Layered Model ======================
	  !====================================================================
	  	!nlay1D    ! Number of layers
	  	!sig1D => ! (S/m) Layer conductivities 
		!zlay1D => ! (m)   Depth to top of each layer, first layer ignored  

		

   endif
   if (txDict(iTx)%Tx_type=='CSEM' .or. txDict(iTx)%Tx_type=='DC') then	
		xTx1D = txDict(iTx)%xyzTx(1)
		yTx1D = txDict(iTx)%xyzTx(2)   
		Call set1DModel(sigma,xTx1D,yTx1D)
   end if
   

!    the following needs work ... want to avoid reinitializing
!     operator coefficients when conductivity does not change;
!     need to have a way to reset sigmaNotCurrent to false when
!     conductivity changes (one idea: assign a random number whenever
!     a conductivity parameter is modified (by any of the routines in
!     module ModelSpace); store this in the modelOperator module (which
!     is where updateCond sits) and have updateCond compare the random key
!     with what is stored)
!  if(sigmaNotCurrent) then
       call updateCond(sigma)
!      sigmaNotCurrent = .false.
!   endif

      if (txDict(iTx)%Tx_type=='DC') then
	      call init_input_for_DC(grid,sigma,iTx)
          !write(6,*) 'init_input_for_DC'
          call setAnomCond(sigma)
          call create_cvector(grid,V_p,EDGE)
     end if	
   
   if (txDict(iTx)%Tx_type=='SFF') then
      ! compute sigma-sigma1D for the source... NOT PHYSICAL!
      Call linComb_modelParam(ONE,sigma,MinusONE,sigmaPrimary,sigmaTemp)
      ! sigmaTemp is the anomalous conductivity, map it onto edges
      Call ModelParamToEdge(sigmaTemp,condAnomaly)
      write(0,*) 'DEBUG size condAnomaly ',condAnomaly%nx,condAnomaly%ny,condAnomaly%nz
   end if
 
if (.NOT. initForSens) then

   if (txDict(iTx)%Tx_type=='CSEM') then	
					   !=======================================================================
					   !============================================== Use for CSEM Calculation
					   !=======================================================================
						  n1D = (grid%Nx)*(grid%Ny+1)*(grid%Nz+1)
						  n1D = n1D + (grid%Nx+1)*(grid%Ny)*(grid%Nz+1)
						  n1D = n1D + (grid%Nx+1)*(grid%Ny+1)*(grid%Nz)

						  if (allocated (x1D)) then  
							 Deallocate(x1D, y1D, z1D)
							 Deallocate(ex1D,ey1D,jz1D)
							 Deallocate(bx1D,by1D,bz1D)
						  end if
						  
						  allocate(x1D(n1D), y1D(n1D), z1D(n1D))
						  allocate (ex1D(n1D),ey1D(n1D),jz1D(n1D))
						  allocate (bx1D(n1D),by1D(n1D),bz1D(n1D))
						
						  call create_cvector(grid,E_p,EDGE)
						  
						
						  !====================================================================
						  ! Create position vector that the primary field has to be calculated
						  !====================================================================
						  counter = 1
						  ! E-field corresponing to these nodes is Ex
						  Do iz = 1,grid%Nz+1 !Edge Z
							  Do iy = 1,grid%Ny+1 !Edge Y
								  Do ix = 1,grid%Nx !Center X
									x1D(counter) = grid%xCenter(ix)
									y1D(counter) = grid%yEdge(iy)
									z1D(counter) = grid%zEdge(iz)
									counter = counter + 1
									End Do
							  End Do
						  End Do
						  
						  ! E-field corresponing to these nodes is Ey
						  Do iz = 1,grid%Nz+1 !Edge Z
							  Do iy = 1,grid%Ny !Center y
								  Do ix = 1,grid%Nx+1 !Edge x
									  x1D(counter) = grid%xEdge(ix)
									  y1D(counter) = grid%yCenter(iy)
									  z1D(counter) = grid%zEdge(iz)
									  counter = counter + 1
								  End Do
							  End Do
						  End Do
						  
						  ! E-field corresponing to these nodes is Ez
						  Do iz = 1,grid%Nz !Center Z
							  Do iy = 1,grid%Ny+1 !Edge y
								  Do ix = 1,grid%Nx+1 !Edge x
									  x1D(counter) = grid%xEdge(ix)
									  y1D(counter) = grid%yEdge(iy)
									  z1D(counter) = grid%zCenter(iz)
									  counter = counter + 1
								  End Do
							  End Do
						  End Do
						  
						  
						  
						  !====================================================================
						  !================== Create Anomalous Conductivity ===================
						  !====================================================================
						  call setAnomCond(sigma) ! This subroutine created by Aihua
							 ! This subroutine is used for creating the anomalous conductivity
							 ! by subtract current conductivity with primary conductivity
						
						  !====================================================================
						  !============ Create or Load Transmitter Parameters =================
						  !====================================================================
						  
							xTx1D = txDict(iTx)%xyzTx(1)
							yTx1D = txDict(iTx)%xyzTx(2)
							zTx1D = txDict(iTx)%xyzTx(3)
							ftx1D = 1.0d0/txDict(iTx)%PERIOD
							sdm1D = txDict(iTx)%Moment           ! (Am), dipole moment. Normalize to unit source moment
							azimuthTx1D = txDict(iTx)%azimuthTx ! (degrees) 
							dipTx1D     = txDict(iTx)%dipTx
							
							HTmethod1D      = 'kk_ht_201'    ! Use 201 point HT digital filters.
							outputdomain1D  = 'spatial'      ! Assume spatial domain comps
							lbcomp          = .false.        ! This is changed to true if magnetics in data file
							lUseSpline1D    = .true.         ! Use spline interpolation for faster 1D computations
							linversion      = .false.        ! Compute derivatives with respect to sigma(layers)
							
							phaseConvention = 'lag'          ! The usual default is lag, where phase becomes larger 
												 !    positive values with increasing range.
						call reset_time(timeDipole)
						call comp_dipole1D               ! Calculate E-Field by Key's code
						timeD = elapsed_time(timeDipole)
						write(*,*) timeD	
	

	

						!====================================================================
						  !================ Map E-field back to modular format=================
						  !====================================================================
						  counter = 1
						  ! E-field corresponing to these nodes is Ex
						   Do iz = 1,grid%Nz+1 !Edge Z
							  Do iy = 1,grid%Ny+1 !Edge Y
								  Do ix = 1,grid%Nx !Center X	  		  	  
									  E_p%x(ix,iy,iz) = ex1D(counter)				  
									  counter = counter + 1
								  End Do
							  End Do
						  End Do
							  
						  ! E-field corresponing to these nodes is Ey
						  Do iz = 1,grid%Nz+1 !Edge Z
							  Do iy = 1,grid%Ny !Center y
								  Do ix = 1,grid%Nx+1 !Edge x	  				  
									  E_p%y(ix,iy,iz) = ey1D(counter)
									  counter = counter + 1
								  End Do
							  End Do
						  End Do
						  
						  ! E-field corresponing to these nodes is Ez
						  Do iz = 1,grid%Nz !Center Z
							  Do iy = 1,grid%Ny+1 !Edge y
								  Do ix = 1,grid%Nx+1 !Edge x
									  E_p%z(ix,iy,iz) = jz1D(counter)
									  counter = counter + 1
								  End Do
							  End Do
						  End Do
		


					
						
    end if
end if	
	  
   ! This needs to be called before solving for a different frequency
   !!!!!!!  BUT AFTER UPDATECOND !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call UpdateFreq(txDict(iTx)%omega)

   end subroutine initSolver

   !**********************************************************************
   subroutine exitSolver(e0,e,comb)
   !   deallocates b0, comb, e0, e and solver arrays
   type(solnVector_t), intent(inout), optional  :: e0
   type(solnVector_t), intent(inout), optional	::e
   type(rhsVector_t), intent(inout), optional	::comb

   ! local variables
   logical			:: initForSens
   character(10)    :: txType = 'MT'
 
   initForSens = present(comb)

   !if(present(e0)) then
   !   if(e0%allocated) then
   !     txType = txDict(e0%tx)%tx_type
   !     if (txType =='DC') then
   !      call de_ini_private_data_DC
   !     end if
   !   endif
   !endif

   if(present(e0)) then
      call deall_solnVector(e0)
   endif

   if(initForSens) then
      call deall_rhsVector(comb)
      call deall_solnVector(e)
   endif

   if(modelDataInitialized) then
      ! cleanup/deallocation routines for model operators
      call ModelDataCleanUp() ! FWD/modelOperator3D.f90
      call ModelOperatorCleanUp() ! FWD/EMsolve3D.f90
      modelDataInitialized = .false.
   endif

   end subroutine exitSolver


!**********************************************************************
! Sets boundary conditions in RHS b0, and the initial conditions in e0,
! for a transmitter index iTx. Typically, this is done for all polarizations.
! However in the MPI case, e0%nPol is artificially set to 1 to ensure optimal
! load distribution, so need to be careful with the indexing here.
! In all cases, this is meant to be called just before fwdSolve ONLY.
! The difference with fwdSolve is that here b0 is computed; in fwdSolve
! it is input only.
!
! Boundary conditions would have already been initialized from file
! into bAll for option BC_FROM_RHS_FILE. Otherwise, need to be computed.
! At present, we are not using this routine to set up the initial value
! of e0, but we might do so in the future.
!
! We are doing all this in a separate routine because we want this to be
! a high-level function that can be called from the upper level.
!
! A. Kelbert, 24 May 2018; last edited 18 Aug 2021
  Subroutine fwdSetup(iTx,e0,b0)

    !  Input mode, period
    integer, intent(in)     :: iTx
    ! Output electric field first guess (for iterative solver)
    type(solnVector_t), intent(inout) :: e0
    ! Output boundary conditions
    type(rhsVector_t), intent(inout)  :: b0

    ! local
    type(cboundary)     :: BC
    integer             :: iMode,j
    real(kind=prec)     :: omega
    complex(kind=prec)  :: i_omega_mu

    ! local variables for TIDE
   integer		 :: ios,istat
   character*80 :: file_name,comment
   character*2  :: tidal_component
   logical		 :: exists
   type(sparsevecc) :: jInt

    ! For RHS vector, the major differences between MT and CSEM are
    ! (1) Boundary condition for MT isn't zeros while CSEM is zeros
    ! (2) CSEM has a source term while that of MT is zeros
    ! Initialize the RHS vector; should we always clean it up on input?
    if (.not. b0%allocated) then
      select case (txDict(iTx)%Tx_type)
      case ('CSEM','DC','SFF')
        b0%nonzero_Source = .true.
        b0%sparse_Source = .false.
        b0%nonzero_BC = .false.
      case ('MT')
        b0%nonzero_Source = .false.
        b0%sparse_Source = .false.
        b0%nonzero_BC = .true.
      case('TIDE')
        b0%nonzero_Source = .true.
        b0%sparse_Source = .false.
        b0%nonzero_BC = .true.
      case default
        write(0,*) node_info,'Unknown FWD problem type',trim(txDict(iTx)%Tx_type),'; unable to initialize RHS'
      end select
      call create_rhsVector(e0%grid,iTx,b0)
    end if


    ! careful here with imode indexing. In MPI, we trick the code into thinking
    ! that there is only one mode for each processor. So instead of using plain indexing
    ! to determine the mode, we use Pol_index integer variable.
    do j = 1,e0%nPol
        iMode = e0%Pol_index(j)

        select case (txDict(iTx)%Tx_type)

            case ('DC')
                ! do nothing for now: all in fwdSolve; will clean up later

            case ('CSEM')
                !  set period, complete setup of 3D EM equation system
                omega = txDict(iTx)%omega
                i_omega_mu = cmplx(0.,ISIGN*MU_0*omega,kind=prec)
                ! Now finish up the computation of the general b0%s = - ISIGN * i\omega\mu_0 j
                call diagMult(condAnomaly,E_P,b0%b(j)%s)
                call scMult(-i_omega_mu,b0%b(j)%s,b0%b(j)%s)

            case ('SFF')
                ! this is currently implemented only for 1 mode - check for this...
                !if (iMode .ne. 1) then
                !  write(0,*) 'ERROR: SFF only implemented for one mode at present. Exiting...'
                !  stop
                !end if
                ! we've read eAllPrimary from EM soln file already
                E_P = eAllPrimary%solns(iTx)%pol(iMode)
                !  set period, complete setup of 3D EM equation system
                omega = txDict(iTx)%omega
                i_omega_mu = cmplx(0.,ISIGN*MU_0*omega,kind=prec)
                ! Now finish up the computation of the general b0%s = - ISIGN * i\omega\mu_0 j
                call diagMult(condAnomaly,E_P,b0%b(j)%s)
                call scMult(-i_omega_mu,b0%b(j)%s,b0%b(j)%s)

            case ('MT')
                if (BC_FROM_RHS_FILE) then
                    ! in this case, we've read bAll from RHS file already
                    write (*,'(a12,a29,a12,i4,a15,i2)') node_info, 'Setting the BC from RHS file ', &
                        ' for period ',iTx,' & mode # ',iMode
                    BC = bAll%combs(iTx)%b(iMode)%bc

                elseif (BC_FROM_E0_FILE) then
                    ! TEMPORARY, TO REPLICATE TIDES - WILL FIX THIS LATER
                    ! we are going to make a huge assumption here: nPol == 1 always for this case
                    !  and of course transmitters are in same order always
                    write (*,'(a12,a28,a12,i4,a15,i2)') node_info, 'Setting the BC from E0 file ', &
                        ' for period ',iTx,' & mode # ',iMode
                    e0%pol(j) = E0_from_file(iTx)
                    call getBC(e0%pol(j),BC)
                    !   do we now need to set boundary edges of E0 == 0?

                elseif (NESTED_BC) then
                    ! The BC are already computed from a larger grid for all transmitters and modes and stored in BC_from_file.
                    ! Overwrite BC with BC_from_file.
                    ! Note [NM]: Right now we are using the same period layout for both grid.
                    ! This why, it is enough to know the period and mode index to pick up the BC from BC_from_file vector.
                    write (*,'(a12,a35,a12,i4,a15,i2)') node_info, 'Setting the BC from nested E0 file ', &
                        ' for period ',iTx,' & mode # ',iMode
                    BC = BC_from_file((iTx*2)-(2-iMode))

                elseif (COMPUTE_BC) then
                    ! For e0 and b0, use the same fake polarization index j for MPI modeling context
                    write (*,'(a12,a28,a12,i4,a15,i2)') node_info, 'Computing the BC internally ', &
                        ' for period ',iTx,' & mode # ',iMode
                    BC = b0%b(j)%bc
                    call ComputeBC(iTx,iMode,e0%pol(j),BC)

                end if
                ! store the BC in b0 and set up the forward problem - use fake indexing in MPI
                b0%b(j)%adj = 'FWD'
                b0%b(j)%bc = BC

            case ('TIDE') ! in the future, may use the BC options from MT block above

               file_name = trim(txDict(iTx)%id)//'.source'
               inquire(FILE=file_name,EXIST=exists)
               if (exists) then
                  write(*,*) node_info,'Reading source - i \omega \mu \sigma_E (v x B) from interior source file: ',trim(file_name)
                  open(ioREAD,file=file_name,status='unknown',form='formatted',iostat=ios)
                  read(ioREAD,'(a35)',iostat=istat) comment
                  read(ioREAD,'(a2)',iostat=istat) tidal_component
                  if (tidal_component .ne. trim(txDict(iTx)%id)) then
                     write(0,*) node_info,'Warning: tidal component ',tidal_component,' is read from file ',trim(file_name)
                  end if
                  call read_sparsevecc(ioREAD,jInt)
                  close(ioREAD)
               end if

               ! Assume that the source - i \omega \mu \sigma_E (v x B) in on the edges
               if(jInt%allocated) then
                  write(0,*) node_info,'Using interior forcing to compute the RHS for the FWD problem'
                  call add_scvector(ISIGN*C_ONE,jInt,b0%b(j)%s)
                  call deall_sparsevecc(jInt)
               end if

            case default
                write(0,*) node_info,'Unknown FWD problem type',trim(txDict(iTx)%Tx_type),'; unable to compute RHS'
        end select
    end do

  end subroutine fwdSetup

   !**********************************************************************
   subroutine fwdSolve(iTx,e0,b0)

   !  driver for 3d forward solver; sets up for transmitter iTx, returns
   !   solution in e0 ; rhs vector (b0) is generated locally--i.e.
   !   boundary conditions are set internally (NOTE: could use transmitter
   !   dictionary to indireclty provide information about boundary
   !    conditions.  Presently we set BC using WS approach.

   integer, intent(in)		:: iTx
   type(solnVector_t), intent(inout)	:: e0
   type(rhsVector_t), intent(in)        :: b0

   ! local variables
   real(kind=prec)	:: omega,term
   integer			:: IER,iMode
   !complex(kind=prec)	:: i_omega_mu
   integer :: iz, ix, iy, counter,k,j,i
   !DC parameters
   character (len=80)              	:: Desc = ''
   

   omega = txDict(iTx)%omega
   !  set period, complete setup of 3D EM equation system
   !i_omega_mu = cmplx(0.,ISIGN*MU_0*omega,kind=prec)
   !write(55,*)omega,period,2.0d0*PI/Period,i_omega_mu 
   !
   !   the rhs from mult of primary field and wave number difference
   !   due to the disturbance of conductivity
   !

   if (txDict(iTx)%Tx_type=='DC') then 
     AdjFWD_DC='FWD_DC'
     !call compute_Pri_potential
     !call put_v_in_V_P(e0)
     !call diagMult(condAnomaly,V_P,b0%s)
     !call extract_RHS_for_FWD_DC(b0,e0)
       call RHS_DC_FWD
	   call FWDsolve3D_DC
       call put_v_in_e(e0)
       !call de_ini_private_data_DC
   elseif (txDict(iTx)%Tx_type=='CSEM') then 
 
   ! Now finish up the computation of the general b0%s = - ISIGN * i\omega\mu_0 j
   !call diagMult(condAnomaly,E_P,b0%s)
   !call scMult(-i_omega_mu,b0%s,b0%s)
   !   call forward solver, compute secondary field
	write(*,'(a12,a3,a25,i3,a14,es15.7)') 'Solving the ','FWD', &
                  ' problem for transmitter ',iTx,' at frequency ',txDict(iTx)%PERIOD
   call zero_solnVector(e0)
   call FWDsolve3D(b0%b(1),omega,e0%pol(1))

   !   add primary field to secondary field
   !e0%pol(1)=E_p
   call add(E_p,e0%pol(1),e0%pol(1))
	  !term=1.0/10.0 ! txDict(iTx)%Moment  
      !call scMult(term,e0%pol(1),e0%pol(1))

   elseif (txDict(iTx)%Tx_type=='SFF') then 
 
      ! General b0%s = - ISIGN * i\omega\mu_0 (sigma-sigma1d) E1D already computed
      do iMode = 1,e0%nPol
         ! Extract primary solution again...
         E_P = eAllPrimary%solns(iTx)%pol(iMode)
		   ! call forward solver, compute secondary field
         ! set the starting solution to zero
		   ! NOTE that in the MPI parallelization, e0 may only contain a single mode;
		   ! mode number is determined by Pol_index, NOT by its order index in e0
		   ! ... but b0 uses the same fake indexing as e0
		   write (*,'(a12,a12,a3,a20,i4,a2,es13.6,a15,i2)') node_info, 'Solving the ','SFF', &
			   	' problem for period ',iTx,': ',(2*PI)/omega,' secs & mode # ',e0%Pol_index(iMode)
         call zero(e0%pol(iMode))
		   call FWDsolve3D(b0%b(iMode),omega,e0%pol(iMode))
		   write (6,*)node_info,'FINISHED solve, nPol',e0%nPol
         ! now add primary field to secondary field
         call add(E_p,e0%pol(1),e0%pol(1))
   	enddo

   elseif ((txDict(iTx)%Tx_type=='MT') .or. (txDict(iTx)%Tx_type=='TIDE')) then
   !call fwdSetup(iTx,e0,b0)
   	do iMode = 1,e0%nPol
		   ! compute boundary conditions for polarization iMode
		   !   uses cell conductivity already set by updateCond
		   !call setBound(iTx,e0%Pol_index(iMode),e0%pol(imode),b0%bc)
		   ! NOTE that in the MPI parallelization, e0 may only contain a single mode;
		   ! mode number is determined by Pol_index, NOT by its order index in e0
		   ! ... but b0 uses the same fake indexing as e0
		   write (*,'(a12,a12,a3,a20,i4,a2,es13.6,a15,i2)') node_info, 'Solving the ','FWD', &
			   	' problem for period ',iTx,': ',(2*PI)/omega,' secs & mode # ',e0%Pol_index(iMode)
		   call FWDsolve3D(b0%b(iMode),omega,e0%pol(iMode))
		   write (6,*)node_info,'FINISHED solve, nPol',e0%nPol
   	enddo
   else
    write(0,*) node_info,'Unknown FWD problem type',trim(txDict(iTx)%Tx_type),'; unable to run fwdSolve'
   endif

   ! update pointer to the transmitter in solnVector
   e0%tx = iTx

   end subroutine fwdSolve

   !**********************************************************************
   subroutine sensSolve(iTx,FWDorADJ,e,comb)
   !   Uses forcing input from comb, which must be set before calling
   !    solves forward or adjoint problem, depending on comb%ADJ
   !  NOTE that this routine DOES NOT call UpdateFreq routine to complete
   !   initialization of solver for a particular frequency.
   !  This final initialization step must (at present) be done by
   !    calling fwdSolve before calling this routine.

   integer, intent(in)          	:: iTx
   character*3, intent(in)		:: FWDorADJ
   type(solnVector_t), intent(inout)		:: e
   type(rhsVector_t), intent(inout)		:: comb

   ! local variables
   integer      			:: IER,iMode
   real(kind=prec) 		:: omega, period

!  zero starting solution, solve for all modes
   call zero_solnVector(e)
   
   if (txDict(iTx)%Tx_type=='MT' .or. txDict(iTx)%Tx_type=='CSEM' ) then 
   	omega = txDict(iTx)%omega
   	period = txDict(iTx)%period
   	do iMode = 1,e%nPol
      		comb%b(e%Pol_index(iMode))%adj = FWDorADJ
      		write (*,'(a12,a12,a3,a20,i4,a2,es13.6,a15,i2)') node_info,'Solving the ',FWDorADJ, &
				' problem for period ',iTx,': ',(2*PI)/omega,' secs & mode # ',e%Pol_index(iMode)
      		call FWDsolve3d(comb%b(e%Pol_index(iMode)),omega,e%pol(imode))
   	enddo
   elseif (txDict(iTx)%Tx_type=='DC') then
      AdjFWD_DC='Adj_DC'
     call extract_RHS_for_Adjoint_DC (comb)
	 call FWDsolve3D_DC
     call put_v_in_e(e)
     !call de_ini_private_data_DC
   else
    write(0,*) node_info,'Unknown FWD problem type',trim(txDict(iTx)%Tx_type),'; unable to run sensSolve'
   endif

   ! update pointer to the transmitter in solnVector
   e%tx = iTx

   end subroutine sensSolve

!=======================================================================
!Mar. 13, 2011================================= Use for CSEM Calculation
!=======================================================================
! I add two subroutines that were created by Aihua
! (1) set1DModel
! (2) setAnomCond

!==========================================================================
!=============================================================== set1DModel
!==========================================================================
subroutine set1DModel(sigma,xTx1D,yTx1D,FromFile)

   !   this is a private routine, used to extract layer averages from
   !   a 3D conductivity parameter (sigma) and set up
   !   (1) nlay1D    ! Number of layers
   !   (2) sig1D => ! (S/m) Layer conductivities 
   !   (3) zlay1D => ! (m)   Depth to top of each layer, first layer ignored the 1D Model  z_P, sigma_P
  
   type(modelParam_t),intent(in), target		:: sigma
   real(kind=prec),intent(in)                           :: xTx1D,yTx1D 
   logical, intent(in), optional                        :: FromFile

   !   local variables ... this is an easy, but not necessarily most efficient
   !   way to get an average background layered conductivity ...
   !    could add routines to modelParameter module to do this more directly
   type(rscalar)	::	 sigmaCell 
   character(len=80)    ::       paramtype,modelFile
   type(rscalar)        ::       model
   type(modelParam_t)   ::       aModel,Anomalous_model,sigma_temp,sigma_temp1
   type(rvector)::  cond


   integer	:: nzEarth,Nz,nzAir,i,j,k,ixTx,iyTx
   real(kind=prec)	:: wt,vAir,asigma
   character(len=256)   ::       PrimaryFile
   
   IF( PRESENT(FromFile) ) THEN
	write(*,*) "Under Develop ^^"
   Else

	!   first define conductivity on cells  
   	!   (extract into variable which is public)
   		call modelParamToCell(sigma, sigmaCell, paramtype)
   		nlay1D = sigmaCell%nz
   		nzEarth = sigmaCell%grid%nzEarth
   		nzAir = sigmaCell%grid%nzAir
          
		 ixTx= minNode(xTx1D, sigmaCell%grid%xEdge)  
		 iyTx= minNode(yTx1D, sigmaCell%grid%yEdge)  

   		!   for layer boundaries use z-edges of 3D grid
   		  	 if(allocated(zlay1D)) then
			    Deallocate(zlay1D, sig1D)
		     end if
			 
			 allocate(zlay1D(nlay1D))
			 allocate(sig1D(nlay1D))
			  do k=1,nlay1D
			  zlay1D(k) = sigmaCell%grid%zEdge(k)
			  end do
			 



   		! For create sig1D, we divide this process into two parts (1) for air layers and 
   		!    (2) for earth layers
   		! For air layer, sig1D equal to air layer conductivity
   		! For earth layer, The Geometric mean is be used to create sig1D

   		sig1D(1:nzAir) = sigmaCell%v(1,1,1:nzAir)
   		do k = nzAir+1,nlay1D
		    sig1D(k) = R_ZERO
			wt = R_ZERO
      			do i = 1,sigmaCell%grid%Nx
         			do j = 1,sigmaCell%grid%Ny
             				wt = wt + sigmaCell%grid%dx(i)*sigmaCell%grid%dy(j)
             				sig1D(k) = sig1D(k) + log(sigmaCell%v(i,j,k))* &
						    sigmaCell%grid%dx(i)*sigmaCell%grid%dy(j)
         			end do
      			end do           
      			sig1D(k) = sig1D(k)/wt
      			sig1D(k) = sigmaCell%v(1,1,k) !exp(sig1D(k)) !sigmaCell%v(1,1,k) !exp(sig1D(k)) 
				!write(22,*)k,1.0/sig1D(k),1.0/sigmaCell%v(1,1,k)
   		end do
   		call getValue_modelParam(sigma,paramType,model,vAir)

   END IF
   
   



	model%v=R_ZERO
   ! Put the background (Primary) "condNomaly" conductivities in ModEM model format
   do k = 1,nzEarth
      asigma = sig1D(k+nzAir)
      if( trim(ParamType) == LOGE) then
        asigma = log(asigma)
      elseif (trim(ParamType) == LOG_10) then
        asigma = log10(asigma)
      end if
      do i = 1,sigmaCell%grid%Nx
         do j = 1,sigmaCell%grid%Ny	    
            model%v(i,j,k) = asigma
         end do
     end do
   end do   
   
!  do k = 1,nzEarth
!     asigma = sig1D(k+nzAir)
!     if( trim(ParamType) == LOGE) asigma = log(asigma)
!     do i = 1,sigmaCell%grid%Nx
!        do j = 1,sigmaCell%grid%Ny
!          if ( k == 1) then
! 		     model%v(i,j,k)= log(sigmaCell%v(i,j,k+nzAir)) !R_ZERO !log(sigmaCell%v(i,j,k+nzAir))
!		   else	 
!            model%v(i,j,k) = asigma !sigmaCell%v(i,j,k+nzAir)-asigma
!		   end if	 
!        end do
!    end do
!  end do   
   
   call copy_modelParam(amodel,sigma)   
   call setType_modelParam(amodel,paramType)
   call setValue_modelParam(amodel,paramType,model,vAir)   
   call ModelParamToEdge(amodel,condNomaly)
   

	  


	  
  !   clean up
   call deall_modelParam(amodel)
   call deall_rscalar(model)
   call deall_rscalar(sigmaCell)
end subroutine set1DModel
   
!==========================================================================
!============================================================== setAnomCond
!==========================================================================
subroutine setAnomCond(sigma)
   !   This is a private routine that sets anomalous conductivity
   !    in module variable condAnomaly using input model parameter sigma,
   !     and layered background conductivity (already set in module 
   !     variables z_P and sigma_P by a call to setPrimaryCond)

   type(modelParam_t),intent(in), target		:: sigma
   
   type(rvector)::  cond
   real (kind = kind(0.0d0)) :: temp
   type(modelParam_t)                   		:: out_sigma,sigma_temp
   character(len=80)    ::       paramtype,modelFile
   type(rscalar)	::	 sigmaCell 
   integer   k,i,j,ix,iy,iz
   real(kind=prec)	:: vAir

   
  
	  
   !   map conductivity onto edges
   Call ModelParamToEdge(sigma,cond)
   
   CondAnomaly = subtract_rvector_f(cond , condNomaly )     

	  
   
   call deall_rvector(cond)
   
end subroutine setAnomCond

!###########################################################
subroutine init_input_for_DC(grid,sigma,iTx)
   type(grid_t), intent(in), target         :: grid
   type(modelParam_t),intent(in), target		:: sigma
   integer 										:: iTx
   
    type(rscalar)	::	 sigmaCell 
    integer         :: ix,iy,iz,nzAir
	character*80    ::   paramtype
	
	call modelParamToCell(sigma, sigmaCell, paramtype)
    nzAir = sigmaCell%grid%nzAir
	!imax---> in X-direction
    !jmax---> in y-direction
    !kmax---> in z-direction   
    !However, in Klaus's code the grid is 90 Deg. rotetted. In our corredinate system, this means:
 	!imax---> in y-direction
    !jmax---> in x-direction
    !kmax---> in z-direction 
    
    imax=grid%Ny+1   ! Number of edags
	jmax=grid%Nx+1   ! Number of edags    
    kmax=sigmaCell%nz-nzAir+1 ! Number of edags in the Earth

	! The same is for the Tx position:
    
          rs1(1,1) = txDict(iTx)%xyzTx(2)
	      rs1(1,2) = txDict(iTx)%xyzTx(1)
	      rs1(1,3) = txDict(iTx)%xyzTx(3)
          
		   if (allocated (x)) then  
		    deallocate(x,y,z)
           end if
		   allocate(x(0:imax+1),y(0:jmax+1),z(0:kmax+1))
           
		      Do ix = 1,imax
		        x(ix) = grid%yEdge(ix)
              end do
              
		      Do iy = 1,jmax
		        y(iy) = grid%xEdge(iy)
              end do
              
			  Do iz = 1,kmax
		        z(iz) = grid%zEdge(iz+nzAir)
              end do
              
		   if (allocated (sigmac)) then
		     deallocate(sigmac)
           end if
           
		    allocate(sigmac(0:imax,0:jmax,0:kmax))
			
			


		    do iz = 1,kmax-1
         			do iy = 1,jmax-1
                       do ix = 1,imax-1
                         sigmac(ix,iy,iz) =sigmaCell%v(iy,ix,nzAir+iz)
         			end do
      			end do           
            end do
            
            call ini_private_data_DC
		
end  subroutine init_input_for_DC
!###############################################################
subroutine put_v_in_V_P(e)
   type(solnVector_t), intent(inout)	:: e
   
   integer ix,iy,iz,nzAir,nx,ny,nz
   
   nzAir = e%grid%nzAir
   Nz=e%grid%nz-nzAir
   Nx=e%grid%nx
   Ny=e%grid%ny
   
   !open(10,file=AdjFWD_DC)
   do iz=1, Nz+1
	   do iy=1, Ny
            do ix=1, Nx+1
     !        write(10,'(f10.2)')v(iy,ix,iz)
	        V_p%y(ix,iy,iz+nzAir)=dcmplx(v(iy,ix,iz),R_ZERO)
		 end do
       end do
   end do
    ! close(10)
end subroutine put_v_in_V_P
!###############################################################
subroutine put_v_in_e(e)
   type(solnVector_t), intent(inout)	:: e
   
   integer ix,iy,iz,nzAir,nx,ny,nz
   
   nzAir = e%grid%nzAir
   Nz=e%grid%nz-nzAir
   Nx=e%grid%nx
   Ny=e%grid%ny
   
   !open(10,file=AdjFWD_DC)

   do iz=1, Nz+1
	   do iy=1, Ny+1
            do ix=1, Nx
	        e%pol(1)%x(ix,iy,iz+nzAir)=dcmplx(v(iy,ix,iz),R_ZERO)
		 end do
       end do
   end do
 
   do iz=1, Nz+1
	   do iy=1, Ny
            do ix=1, Nx+1
	        e%pol(1)%y(ix,iy,iz+nzAir)=dcmplx(v(iy,ix,iz),R_ZERO)
		 end do
       end do
   end do
   
   do iz=1, Nz
	   do iy=1, Ny+1
            do ix=1, Nx+1
	        e%pol(1)%z(ix,iy,iz+nzAir)=dcmplx(v(iy,ix,iz),R_ZERO)
		 end do
       end do
   end do
    ! close(10)
end subroutine put_v_in_e
!###############################################################
subroutine extract_RHS_for_FWD_DC (b0,e)
  type(RHS_t), intent(in)		:: b0
  type(solnVector_t), intent(inout)	:: e
   integer ix,iy,iz,nzAir,counter,nx,ny,nz
   !comb%b(1)%s%y=R_Zero
   
   nzAir = e%grid%nzAir
   Nz=e%grid%nz-nzAir
   Nx=e%grid%nx
   Ny=e%grid%ny
   
   b_RHS_DC=R_Zero

   
   !   open(10,file=AdjFWD_DC)
   do iz=1, Nz+1
	   do iy=1, Ny
            do ix=1, Nx+1
             b_RHS_DC(iy,ix,iz)= (real(b0%s%y(ix,iy,iz+nzAir)))
         end do
       end do
   end do
!100   close(10)
end subroutine extract_RHS_for_FWD_DC  
!##############################################################
subroutine take_one_component(e0)
  type(solnVector_t), intent(inout)	:: e0
    type(solnVector_t)	:: e_temp
    e_temp=e0
     call zero_solnVector(e0)
     e0%pol(1)%y=real(e_temp%pol(1)%y)
     

      call deall_solnVector(e_temp) 
end subroutine take_one_component
!###############################################################
subroutine take_one_component_of_comb(comb,itx)
  type(rhsVector_t), intent(inout)		:: comb
  integer,intent(in)		:: itx
   integer ix,iy,iz,nzAir,counter,nx,ny,nz
   type(rhsVector_t)		:: comb_temp
   
            call create_rhsVector(comb%b(1)%s%grid,iTx,comb_temp)
            !comb_temp%b(1)=(comb%b(1))
      
            !call deall_rhsVector (comb)
     !call zero_RHS(comb%b(1))
     !comb%b(1)=comb_temp%b(1)
     !call create_rhsVector(comb%b(1)%s%grid,iTx,comb)
    nzAir = comb%b(1)%s%grid%nzAir
   Nz=comb%b(1)%s%grid%nz-nzAir
   Nx=comb%b(1)%s%grid%nx
   Ny=comb%b(1)%s%grid%ny
  

      
     !comb%b(1)%s%y=real(comb_temp%b(1)%s%y)
          call deall_rhsVector(comb_temp)
          
end subroutine take_one_component_of_comb
!###############################################################
subroutine extract_RHS_for_Adjoint_DC (comb)
  type(rhsVector_t), intent(in)		:: comb
   integer ix,iy,iz,nzAir,counter,nx,ny,nz
   !comb%b(1)%s%y=R_Zero
   
   nzAir = comb%b(1)%s%grid%nzAir
   Nz=comb%b(1)%s%grid%nz-nzAir
   Nx=comb%b(1)%s%grid%nx
   Ny=comb%b(1)%s%grid%ny
   
   b_RHS_DC=R_Zero

   
   !   open(10,file=AdjFWD_DC)
   do iz=1, Nz
	   do iy=1, Ny+1
            do ix=1, Nx+1
             b_RHS_DC(iy,ix,iz)= -(real(comb%b(1)%s%z(ix,iy,iz+nzAir)))
            ! write(10,*) iz,iy,ix,b_RHS_DC(iy,ix,iz)
         end do
       end do
   end do
!100   close(10)
end subroutine extract_RHS_for_Adjoint_DC    
 !**********************************************************************
  ! uses nestedEM module to extract the boundary conditions directly from
  ! a full EMsolnMTX vector on a larger (and coarser) grid

  subroutine Interpolate_BC_from_E_soln(eAll_larg,Larg_Grid,grid)

  type(grid_t)  ,intent(in)                 ::  Larg_Grid
  type(solnVectorMTX_t),intent(in)          ::  eAll_larg
  type(grid_t),intent(in)                   ::  Grid

    ! local variables needed for nesting calculations
  integer                      :: status,iMode,iTx,counter

  Call setup_BC_from_file(Grid,eAll_larg%nTx,eAll_larg%solns(1)%nPol)

  ! For now extract the BC from eAll_larg
    counter=0
    do iTx = 1, eAll_larg%nTx
       do iMode=1, eAll_larg%solns(iTx)%nPol
         counter=counter+1
         Call compute_BC_from_file(Larg_Grid,eAll_larg%solns(iTx)%pol(iMode),Grid,counter)
       end do
    end do

  end subroutine Interpolate_BC_from_E_soln


end module ForwardSolver
