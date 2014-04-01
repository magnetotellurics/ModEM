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

!    keep data structures used only by
!    routines in this module private
!   rhs data structures for solving forward, sensitivity probs
type(RHS_t), save, private			:: b0

!=======================================================================
!Mar. 13, 2011============================== Use only for MT Calculation
!=======================================================================

!! Used for interpolating the BC from a larger grid.
  type(grid_t),save,public               	 ::	 Larg_Grid   
  type(solnVectorMTX_t),save,public          ::  eAll_larg
  integer              ,save,public          ::   nTx_nPol
!!  initialization routines (call Fwd version if no sensitivities are
!!     are calculated).  Note that these routines are set up to
!!    automatically manage memory and to figure out which initialization
!!    (or reinitialization) steps are required (e.g., when the frequency
!!    changes from the previous solver call, appropriate solver
!!    coefficients are updated, matrices factored, etc.).  This
!!    functionality needs to be maintained in implementations for new
!!    problems!

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
subroutine Interpolate_BC(grid)
  type(grid_t), intent(in)        :: grid
  integer                                  :: iTx,ix,iy,iz

  


	 b0%nonzero_Source = .false.
     b0%nonzero_bc = .true.
	 
	 
     b0%adj = 'FWD'
     b0%sparse_Source = .false.
     iTx=1
     call create_RHS(grid,iTx,b0)
	!In case of interpolating the BC from eAll_larg   
	! If eAll_ larg solution is already allocated, then use that to interpolate the BC from it
     
	  if (eAll_larg%allocated) then
           write(15,*) ' Start interploating',grid%nx,grid%ny,grid%nz,b0%bc%nx,b0%bc%ny,b0%bc%nz
		call Interpolate_BC_from_E_soln (eAll_larg,Larg_Grid,Grid,b0)
        !Once we are ready from eAll_larg, deallocate it, and keep track, that BC_from_file are already Initialized.
        call deall(eAll_larg)
        call deall_grid(Larg_Grid)
         BC_from_file_Initialized=.true.
      end if
        write(15,*) ' End interploating',BC_from_file(1)%yXMin(10,11)
  


        
      
end subroutine Interpolate_BC
!**********************************************************************
subroutine ini_BC_from_file(grid)
  type(grid_t), intent(in)        :: grid    
  integer                                  :: iTx,j,status
  
	 b0%nonzero_Source = .false.
     b0%nonzero_bc = .true.
     b0%adj = 'FWD'
     b0%sparse_Source = .false.
     iTx=1
     call create_RHS(grid,iTx,b0)
     
     
     allocate (BC_from_file(nTx_nPol), STAT=status)
     do j=1,nTx_nPol
       BC_from_file(j)=b0%bc
     end do
     BC_from_file_Initialized=.true.
     

end subroutine ini_BC_from_file     
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
   character*80 :: gridType,paramtype,file_name
   logical		:: initForSens,sigmaNotCurrent

   type(timer_t) :: timeDipole
   real :: timeD
   

	

   initForSens = present(comb)
   
   !=======================================================================
   !Mar. 13, 2011============================== Use only for MT Calculation
   !=======================================================================

!	!  allocate for scratch rhsVector structure for background, sensitivity
 if (txDict(iTx)%Tx_type=='CSEM' .or. txDict(iTx)%Tx_type=='DC') then
	 b0%nonzero_Source = .true.
     b0%nonzero_bc = .false.
 elseif(txDict(iTx)%Tx_type=='MT') then
	 b0%nonzero_Source = .false.
     b0%nonzero_bc = .true.
end if
	 
	 
     b0%adj = 'FWD'
     b0%sparse_Source = .false.
     call create_RHS(grid,iTx,b0)

   	!!In case of interpolating the BC from eAll_larg   
	!! If eAll_ larg solution is already allocated, then use that to interpolate the BC from it
	!  if (eAll_larg%allocated) then
	!	call Interpolate_BC_from_E_soln (eAll_larg,Larg_Grid,grid,b0%bc)
    !    !Once we are ready from eAll_larg, deallocate it, and keep track, that BC_from_file are already Initialized.
    !   call deall(eAll_larg)
    !    call deall_grid(Larg_Grid)
    !   BC_from_file_Initialized=.true.
    !  end if
      !write(6,*) BC_from_file_Initialized,b0%s%allocated,itx,txDict(iTx)%Tx_type
      if (BC_from_file_Initialized) then
        b0%bc%read_E_from_file=.true.
      end if

   !=======================================================================
   !Mar. 13, 2011================================= Use for CSEM Calculation
   !=======================================================================  
   ! For RHS vector, the major differences between MT and CSEM are 
   ! (1) Boundary condition for MT isn't zeros while CSEM is zeros
   ! (2) CSEM has a source term while that of MT is zeros  

   


        
   !  allocate for background solution
   call create_solnVector(grid,iTx,e0)

   if(initForSens) then
      !  allocate for sensitivity solution, RHS
      call create_solnVector(grid,iTx,e)
      call create_rhsVector(grid,iTx,comb)
      do k = 1,comb%nPol
	   if (txDict(iTx)%Tx_type=='CSEM' .or. txDict(iTx)%Tx_type=='DC') then
        comb%b(k)%nonzero_source = .true.
        comb%b(k)%nonzero_bc = .false.
	   elseif(txDict(iTx)%Tx_type=='MT') then
        comb%b(k)%nonzero_source = .true.
        comb%b(k)%nonzero_bc = .false.
       end if		
        !  assuming here that we don't use sparse storage ... we could!
        comb%b(k)%sparse_Source = .false.
        comb%b(k)%adj = ''
        !  using all this information, reallocate storage for each polarization
        call create_RHS(grid,iTx,comb%b(k))
      enddo
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

   initForSens = present(comb)

   call deall_RHS(b0)
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

   if (txDict(e0%tx)%tx_type =='DC') then
     call de_ini_private_data_DC
  end if
   
   
   end subroutine exitSolver

   !**********************************************************************
   subroutine fwdSolve(iTx,e0)

   !  driver for 3d forward solver; sets up for transmitter iTx, returns
   !   solution in e0 ; rhs vector (b0) is generated locally--i.e.
   !   boundary conditions are set internally (NOTE: could use transmitter
   !   dictionary to indireclty provide information about boundary
   !    conditions.  Presently we set BC using WS approach.
   !  NOTE that this routine calls UpdateFreq routine to complete
   !   initialization of solver for a particular frequency.

   integer, intent(in)		:: iTx
   type(solnVector_t), intent(inout)	:: e0

   ! local variables
   real(kind=prec)	:: period, omega,term
   integer			:: IER,iMode
   complex(kind=prec)	:: i_omega_mu
   integer :: iz, ix, iy, counter,k,j,i
   !DC parameters
   character (len=80)              	:: Desc = ''
   

   omega = txDict(iTx)%omega
   period = txDict(iTx)%period
   !  set period, complete setup of 3D EM equation system
   i_omega_mu = cmplx(0.,-1.0d0*ISIGN*MU_0*omega,kind=prec)
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
 
   call diagMult(condAnomaly,E_P,b0%s)
   call scMult(i_omega_mu,b0%s,b0%s)
   !   call forward solver, compute secondary field
	write(*,'(a12,a3,a25,i3,a14,es12.7)') 'Solving the ','FWD', &
                  ' problem for transmitter ',iTx,' at frequency ',txDict(iTx)%PERIOD
   call zero_solnVector(e0)
   call FWDsolve3D(b0,omega,e0%pol(1))

   !   add primary field to secondary field
   !e0%pol(1)=E_p
   call add(E_p,e0%pol(1),e0%pol(1))
	  !term=1.0/10.0 ! txDict(iTx)%Moment  
      !call scMult(term,e0%pol(1),e0%pol(1))
		  
elseif (txDict(iTx)%Tx_type=='MT') then
   do iMode = 1,e0%nPol
      ! compute boundary conditions for polarization iMode
      !   uses cell conductivity already set by updateCond
      call SetBound(e0%Pol_index(iMode),period,e0%pol(imode),b0%bc,iTx)
      write (*,'(a12,a12,a3,a20,i4,a2,es12.6,a15,i2)') node_info, 'Solving the ','FWD', &
				' problem for period ',iTx,': ',(2*PI)/omega,' secs & mode # ',e0%Pol_index(iMode)
      call FWDsolve3D(b0,omega,e0%pol(imode))
      write (6,*)node_info,'FINISHED solve, nPol',e0%nPol
   enddo
end if
   
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
      write (*,'(a12,a12,a3,a20,i4,a2,es12.6,a15,i2)') node_info,'Solving the ',FWDorADJ, &
				' problem for period ',iTx,': ',(2*PI)/omega,' secs & mode # ',e%Pol_index(iMode)
      call FWDsolve3d(comb%b(e%Pol_index(iMode)),omega,e%pol(imode))
   enddo
elseif (txDict(iTx)%Tx_type=='DC') then
      AdjFWD_DC='Adj_DC'
     call extract_RHS_for_Adjoint_DC (comb)
	 call FWDsolve3D_DC
     call put_v_in_e(e)
     !call de_ini_private_data_DC
end if

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
      if( trim(ParamType) == LOGE) asigma = log(asigma)
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

  subroutine Interpolate_BC_from_E_soln(eAll_larg,Larg_Grid,grid,b0)

  type(grid_t)  ,intent(in)                 ::  Larg_Grid
  type(solnVectorMTX_t),intent(in)          ::  eAll_larg
  type(grid_t),intent(in)                   ::  Grid
  type(RHS_t),intent(in)                    ::  b0

    ! local variables needed for nesting calculations
  integer                      :: status,iMode,iTx,counter

  Call setup_BC_from_file(Grid,b0%bc,eAll_larg%nTx,eAll_larg%solns(1)%nPol)

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
