module SensComp
!  higher level module that calls routines in DataSens, SolverSens
!  and ForwardSolver to implement Jacobian computations, including
!  1) computation of predicted data functionals;
!  2) multiplication by J and J^T;
!  3) evaluation of the full sensitivity matrix J.
!
!  Before any of these routines may be called, the transmitter (txDict),
!  data type (typeDict) and receiver (rxDict) dictionaries must be
!  created and initialized. These are heavily used by Level II routines
!  in DataFunc, SolverSens and ForwardSolver, inherited by this module;
!  "pointers" to dictionary entries are attached to data vector d.

  use math_constants
  use utilities
  use SensMatrix
  use DataSens
  use SolverSens
  use ForwardSolver

  implicit none

  public 	:: calcJ,   Jrows
  public	:: Jmult,   Jmult_TX
  public	:: JmultT,  JmultT_TX
  public	:: fwdPred, fwdPred_TX
  public	:: setGrid, cleanUp

  ! local timer
  type (timer_t), save, private         :: timer

  ! numerical discretization used to compute the EM solution
  !  (may be different from the grid stored in model parameter)
  type (grid_t), target, save, private     :: grid

  ! temporary EM fields, that are saved for efficiency - to avoid
  !  memory allocation & deallocation for each transmitter
  type(solnVector_t), save, private		:: e,e0
  type(rhsVector_t) , save, private		:: comb


Contains

  !**********************************************************************
  subroutine Jrows(iTx,iDt,iRx,sigma0,emsoln,Jreal,Jimag)
   !  Given background model parameter sigma0 and background
   !  solution vector emsoln, calculate a row of sensitivity matrix
   !  for specified transmitter, data type and receiver,
   !  for all data type components.
   !  Generic form that will work for any problem.
   !   returns a model parameter for each real data type component.
   !
   !   indices into transmitter, receiver & data type dictionaries
   integer, intent(in)				              :: iTx, iRx, iDt
   !   sigma0 is background conductivity parameter
   type(modelParam_t), intent(in)	          :: sigma0
   !   to use this routine, need to first supply background EM soln
   type(solnVector_t), intent(in)  	        :: emsoln
   !   Jreal is the output array of data sensitivities,
   !   one for each component in the data type.  Each sensitivity
   !    is an element of type modelParam, an abstract
   !    data type that defines the unknown conductivity
   !   NOTE: Jreal and Jimag both exist regardless of whether the data
   !     are real or complex, since J itself is complex
   type(modelParam_t), pointer   	          :: Jreal(:), Jimag(:)

   !  local variables
   integer 		:: istat,ii,nFunc,nComp,iFunc
   type(sparseVector_t), pointer	:: L(:)
   type(modelParam_t), pointer    :: Qreal(:),Qimag(:)
   logical      :: Qzero

   ! the vectors Jreal, Jimag have to be allocated
   if((.not. associated(Jreal)) .or. (.not. associated(Jimag))) then
      call errStop('output rows of J have to be allocated in Jrows')
   endif
   nFunc = size(Jreal)
   if(size(Jimag) .ne. nFunc) then
      call errStop('outputs of incompatible size in Jrows')
   endif

   ! allocate and initialize sensitivity values
   do iFunc = 1,nFunc
      ! this makes a copy of modelParam, then zeroes it
      Jreal(iFunc) = sigma0
      call zero(Jreal(iFunc))
      Jimag(iFunc) = sigma0
      call zero(Jimag(iFunc))
   enddo

   !  manage any necessary initialization for this transmitter
   call initSolver(iTx,sigma0,grid,e0,e,comb)

   !  store the provided solnVector in e0
   e0 = emsoln

   allocate(L(nFunc),STAT=istat)
   allocate(Qreal(nFunc),STAT=istat)
   allocate(Qimag(nFunc),STAT=istat)
   
	  do iFunc=1,nFunc
		  call create_sparseVector(e0%grid,iTx,L(iFunc))
	  end do
   ! compute linearized data functional(s) : L
   call Lrows(e0,sigma0,iDt,iRx,L)

   ! compute linearized data functional(s) : Q
   call Qrows(e0,sigma0,iDt,iRx,Qzero,Qreal,Qimag)

   ! loop over functionals  (e.g., for 2D TE/TM impedances nFunc = 1)
   do iFunc = 1,nFunc

      ! solve transpose problem for each of nFunc functionals
      call zero_rhsVector(comb)
      call add_sparseVrhsV(C_ONE,L(iFunc),comb)

      call sensSolve(iTx,TRN,e,comb)

      ! multiply by P^T and add the rows of Q
      call PmultT(e0,sigma0,e,Jreal(iFunc),Jimag(iFunc))
      if (.not. Qzero) then
        call scMultAdd(ONE,Qreal(iFunc),Jreal(iFunc))
        call scMultAdd(ONE,Qimag(iFunc),Jimag(iFunc))
      endif

      ! deallocate temporary vectors
      call deall_sparseVector(L(iFunc))
      call deall_modelParam(Qreal(iFunc))
      call deall_modelParam(Qimag(iFunc))

   enddo  ! iFunc

   !  deallocate local arrays
   deallocate(L,STAT=istat)
   deallocate(Qreal,STAT=istat)
   deallocate(Qimag,STAT=istat)

  end subroutine Jrows

  !**********************************************************************
  subroutine calcJ(d,sigma0,sens)
   !  Calculate sensitivity matrix for data in d;
   !  Generic form that will work for any problem.
   !
   !   d is the input data vector, here just used to identify
   !     receiver, data type and transmitter combinations
   !     to compute sensitivities for
   type(dataVectorMTX_t), intent(in)	  :: d
   !   sigma0 is background conductivity parameter
   type(modelParam_t), intent(in)	      :: sigma0
   !   sens is the output structure of data sensitivities,
   !     stored as model parameters,
   !     one for each element in the data array.
   type(sensMatrix_t), pointer 		      :: sens(:)

   !  local variables
   type(modelParam_t), pointer   :: Jreal(:),Jimag(:)
   type(dataBlock_t)             :: dTemplate
   logical      :: isComplex
   integer 		  :: i,j,nTx,k,nSite,nTotal,ii,iTx, &
				iDT,nfunc,ncomp,iRx,iFunc,iComp,istat

   ! nTotal is number of real data
   nTotal = countData(d)
   write(0,'(a32,i6,a15)') 'Computing the sensitivities for ',nTotal,' data values...'

   ! now, allocate for sensitivity values, if necessary
   if(.not. associated(sens)) then
      call create_sensMatrixMTX(d, sigma0, sens)
   endif

   ! nTx is number of transmitters
   nTx = d%nTx

   ! loop over frequencies, computing all sensitivities for
   !   one frequency
   do j = 1,nTx

      ! identify the transmitter
      iTx = d%d(j)%tx

      !  manage any necessary initilization for this transmitter
      call initSolver(iTx,sigma0,grid,e0,e,comb)

      !  solve forward problem; result is stored in e0
      call fwdSolve(iTx,e0)

      ! now loop over data types
      do i = 1,d%d(j)%nDt

        ! make a local copy of the dataBlock for this transmitter and dataType
        dTemplate = d%d(j)%data(i)

        ! get index into dataType dictionary for this dataVec
        iDT = dTemplate%dataType

        ! get the dimensions of this data vector
        nSite = dTemplate%nSite
        nComp = dTemplate%nComp

        ! keep the user informed
        write(0,'(a35,i4,a12,i4,a4,i4,a6)') &
        	'Computing the sensitivities for tx ',iTx,' & dataType ',iDt,' at ',nSite,' sites'

        ! allocate the rows of Jacobian operator
        isComplex = dTemplate%isComplex
		    if(isComplex) then
		       !  data are complex; one sensitivity calculation can be
		       !   used for both real and imaginary parts
		       if(mod(nComp,2).ne.0) then
		         call errStop('for complex data # of components must be even in calcJ')
		       endif
		       nFunc = nComp/2
		    else
		       !  data are treated as real: full sensitivity computation is required
		       !   for each component
		       nFunc = nComp
		    endif
		    allocate(Jreal(nFunc),STAT=istat)
        allocate(Jimag(nFunc),STAT=istat)

        ! loop over sites, computing sensitivity for all components for each site
        do k = 1,nSite

           iRx = dTemplate%rx(k)

           ! compute the sensitivities for these transmitter, data type & receiver
           call Jrows(iTx,iDt,iRx,sigma0,e0,Jreal,Jimag)

           ! store in the full sensitivity matrix
           ii = 1
           do iFunc = 1,nFunc
              if(isComplex) then
                 sens(j)%v(i)%dm(ii,k)   = Jreal(iFunc)
                 sens(j)%v(i)%dm(ii+1,k) = Jimag(iFunc)
                 ii = ii + 2
              else
                 ! for real data, throw away the imaginary part
                 sens(j)%v(i)%dm(ii,k)   = Jreal(iFunc)
                 ii = ii + 1
              endif
           enddo

        enddo  ! sites

        ! deallocate temporary vectors
		    do iFunc = 1,nFunc
		       call deall_modelParam(Jreal(iFunc))
		       call deall_modelParam(Jimag(iFunc))
		    enddo
		    deallocate(Jreal, STAT=istat)
		    deallocate(Jimag, STAT=istat)
		    call deall_dataBlock(dTemplate)

      enddo  ! dataType's
   enddo  ! tx

  end subroutine calcJ

  !**********************************************************************
  subroutine Jmult_TX(dsigma,sigma0,d,emsoln)

   !  Calculate product of sensitivity matrix and a model parameter
   !    for one transmitter (but possibly for multiple data types)
   !
   !  First need to set up transmitter and receiver dictionaries;
   !  indicies into dictionaries are attached to data vector d .
   !
   !  If optional input parameter e1 is present, it must contain
   !   solutions for this transmitter for conductivity sigma0
   !   (This can be used to avoid recomputing forward solutions)
   !
   !   sigma0 is background conductivity parameter
   type(modelParam_t), intent(in)	:: sigma0
   !   dsigma is the input conductivity parameter perturbation
   type(modelParam_t), intent(in)	:: dsigma
   !   d is the computed (output) data vector, also used to identify
   !     receiver transmitter pairs for various computations
   type(dataVector_t), intent(inout)		:: d
   type(solnVector_t), intent(inout), optional	:: emsoln

   !  local variables
   type(dataVector_t) :: d1,d2
   integer 		:: i,j,iTx,iDT
   logical		:: savedSolns

	  savedSolns = present(emsoln)
	  if(savedSolns) then
	     if(d%tx .ne. emsoln%tx) then
	        call errStop('EM soln and d do not agree in Jmult')
	     endif
	  endif

	  !   get index into transmitter dictionary for this dataVector
	  iTx = d%tx

	  !  manage any necessary initilization for this transmitter
	  call initSolver(iTx,sigma0,grid,e0,e,comb)

	  !  initialize the temporary data vectors
	  d1 = d
	  d2 = d

	  if(savedSolns) then
	     ! e0 = emsoln
	     call copy_solnVector(e0,emsoln)
	  else
	     ! solve forward problem; result is stored in e0
	     call fwdSolve(iTx,e0)
	  endif

	  !  compute rhs (stored in comb) for forward sensitivity
	  !  calculation, using conductivity perturbations and
	  !  background soln:
	  call Pmult(e0,sigma0,dsigma,comb)

	  ! solve forward problem with source in comb
	  call sensSolve(iTx,FWD,e,comb)

	  ! multiply e by operator L to obtain perturbation in the data
	  call Lmult(e0,sigma0,e,d1)

	  ! multiply dsigma by operator Q to obtain perturbation in the data
	  call Qmult(e0,sigma0,dsigma,d2)

	  ! add the two together to compute the total sensitivity
	  call linComb_dataVector(ONE,d1,ONE,d2,d)

	  ! clean up
	  call deall_dataVector(d1)
	  call deall_dataVector(d2)

  end subroutine Jmult_TX


  !**********************************************************************
  subroutine Jmult(dsigma,sigma0,d,eAll)

   !  Calculate product of sensitivity matrix and a model parameter
   !    for all transmitters in a datavector (i.e., multiple dataVec objects)
   !
   !  First need to set up transmitter and receiver dictionaries;
   !  indicies into dictionaries are attached to data vector d .
   !
   !  If optional input parameter eAll is present, it must contain
   !   solutions for all transmitters for conductivity sigma0
   !   (This can be used to avoid recomputing forward solutions)
   !
   !   sigma0 is background conductivity parameter
   type(modelParam_t), intent(in)	:: sigma0
   !   dsigma is the input conductivity parameter perturbation
   type(modelParam_t), intent(in)	:: dsigma
   !   d is the computed (output) data vector, also used to identify
   !     receiver transmitter pairs for various computations
   type(dataVectorMTX_t), intent(inout)		:: d
   type(solnVectorMTX_t), intent(inout), optional	:: eAll

   !  local variables
   integer 		:: j
   logical		:: savedSolns

   savedSolns = present(eAll)
   if(savedSolns) then
      if(d%nTx .ne. eAll%nTx) then
         call errStop('dimensions of eAll and d do not agree in Jmult_MTX')
      endif
   endif

   ! loop over frequencies : solve forward system twice,
   !    to compute background, perturbation solutions
   !    apply data functionals (linearized about background soln)
   !    to perturbation solution
   do j = 1,d%nTx

    if(savedSolns) then
   	   ! compute d = J x m for a single transmitter
   	   call Jmult_TX(dsigma,sigma0,d%d(j),eAll%solns(j))
    else
       ! do not pass the EM soln to Jmult - it will be computed
   	   call Jmult_TX(dsigma,sigma0,d%d(j))
    endif

   enddo  ! tx


  end subroutine Jmult


  !**********************************************************************
  subroutine JmultT_TX(sigma0,d,dsigma,emsoln)

   !  Transpose of Jmult multiplied by data vector d for one transmitter;
   !   output is a single conductivity parameter in dsigma
   !
   !   sigma0 is background conductivity parameter
   !
   !  If optional input parameter e0 is present, it must contain
   !   solutions for this transmitter for conductivity sigma0
   !   (This can be used to avoid recomputing forward solutions,
   !    e.g., in a CG solution scheme)

   type(modelParam_t), intent(in)	          :: sigma0
   !   d is the computed (output) data vector, also used to identify
   !     receiver transmitter pairs for various computations
   type(dataVector_t), intent(in)		        :: d
   !   dsigma is the output conductivity parameter
   type(modelParam_t), intent(inout)  	    :: dsigma
   type(solnVector_t), intent(in), optional	:: emsoln

   !  local variables
   type(modelParam_t)	:: Qcomb
   integer 		:: i,j,iTx,iDT
   logical		:: savedSolns

	  savedSolns = present(emsoln)
	  if(savedSolns) then
	     if(d%tx .ne. emsoln%tx) then
	        write(0,*) 'd%tx, emsoln%tx: ',d%tx,emsoln%tx
	        call errStop('solution and data vectors do not agree in JmultT')
	     endif
	  endif

    ! initialize local variable and output
    dsigma = sigma0
    call zero(dsigma)
    Qcomb = dsigma

	  !   get index into transmitter dictionary for this dataVector
	  iTx = d%tx

	  !  manage any necessary initilization for this transmitter
	  call initSolver(iTx,sigma0,grid,e0,e,comb)

	  if(savedSolns) then
	     ! e0 = e1
	     call copy_solnVector(e0,emsoln)
	  else
	     ! solve forward problem; result is stored in e0
	     call fwdSolve(iTx,e0)
	  endif

    ! set up comb using linearized data functionals
    call LmultT(e0,sigma0,d,comb)
    ! solve transpose problem with source in comb
    call sensSolve(iTx,TRN,e,comb)

    ! map from edges to conductivity space ... result P^T e
    !  NOTE: here we throw away imaginary part, even for complex
    !     data (conceivably might want to save this in some cases!)
    call PmultT(e0,sigma0,e,dsigma)

    !  ... for dataVectors with data functionals depending on conductivity
    !   parameter also compute analogous comb in parameter space ...
    !  throw away imaginary part here
    call QmultT(e0,sigma0,d,Qcomb)

    !  add Qcomb to P^T e
    call scMultAdd(ONE,Qcomb,dsigma)

    !  clean up
    call deall_modelParam(Qcomb)

  end subroutine JmultT_TX


  !**********************************************************************
  subroutine JmultT(sigma0,d,dsigma,eAll)

   !  Transpose of Jmult multiplied by data vector d for all transmitters;
   !   output is a single conductivity parameter in dsigma
   !
   !   sigma0 is background conductivity parameter
   !
   !  If optional input parameter eAll is present, it must contain
   !   solutions for all transmitters for conductivity sigma0
   !   IN THE PROPER ORDER (at present) !!!!
   !   (This can be used to avoid recomputing forward solutions,
   !    e.g., in a CG solution scheme)
   !
   !  First need to set up transmitter and receiver dictionaries;
   !  "pointers" to dictionary entries are attached to multi-transmitter
   !   data vector d

   type(modelParam_t), intent(in)	:: sigma0
   !   d is the computed (output) data vector, also used to identify
   !     receiver transmitter pairs for various computations
   type(dataVectorMTX_t), intent(in)		:: d
   !   dsigma is the output conductivity parameter
   type(modelParam_t), intent(out)  	:: dsigma
   type(solnVectorMTX_t), intent(in), optional	:: eAll

   !  local variables
   type(modelParam_t)	:: sigmaTemp
   integer 		:: j
   logical		:: savedSolns

   savedSolns = present(eAll)
   if(savedSolns) then
      if(d%nTx .ne. eAll%nTx) then
         call errStop('dimensions of eAll and d do not agree in JmultT_MTX')
      endif
   endif

   dsigma = sigma0
   call zero(dsigma)

   ! loop over transmitters
   do j = 1,d%nTx

    if(savedSolns) then
   	   ! compute sigmaTemp = JT x d for a single transmitter
   	   call JmultT_TX(sigma0,d%d(j),sigmaTemp,eAll%solns(j))
    else
       ! do not pass the EM soln to JmultT - it will be computed
   	   call JmultT_TX(sigma0,d%d(j),sigmaTemp)
    endif
   	! ... add to dsigma
   	call linComb_modelParam(ONE,dsigma,ONE,sigmaTemp,dsigma)

   enddo  ! tx

   !  clean up
   call deall_modelParam(sigmaTemp)

  end subroutine JmultT


  !**********************************************************************
  subroutine fwdPred_TX(sigma,d,emsoln)

   !  Calculate predicted data for single-transmitter dataVector object d
   !    and for conductivity parameter sigma
   !  Also returns the EM solution e for this transmitter
   !
   !  First need to set up transmitter, receiver, dataType dictionaries;
   !  "pointers" to dictionaries are attached to data vector d .
   !
   !   sigma is input conductivity parameter
   type(modelParam_t), intent(in)	:: sigma
   !   d is the computed (output) data vector, also used to identify
   !     receiver/transmitter
   type(dataVector_t), intent(inout)	:: d
   !  structure containing array of solution vectors (should be
   !   allocated before calling)
   type(solnVector_t), intent(inout)	:: emsoln

   ! local variables
   integer				:: iTx,iDt,i,j

      if(.not.d%allocated) then
         call errStop('data vector not allocated on input to fwdPred')
      end if

      ! get index into transmitter dictionary for this dataVector
      iTx = d%tx

      !  do any necessary initialization for transmitter iTx
      call initSolver(iTx,sigma,grid,emsoln)

      ! compute forward solution
      call fwdSolve(iTx,emsoln)

      ! cycle over data types
      do i = 1,d%nDt

         ! set errorBar=.false. since predicted data do not have
         ! well-defined error bars (important for the inversion)
         d%data(i)%errorBar = .false.

         iDt = d%data(i)%dataType

         ! apply data functionals - loop over sites
		     do j = 1,d%data(i)%nSite

		        ! output is a real vector: complex values come in pairs
		        call dataResp(emsoln,sigma,iDt,d%data(i)%rx(j),d%data(i)%value(:,j))


		     enddo

      enddo

  end subroutine fwdPred_TX


  !**********************************************************************
  subroutine fwdPred(sigma,d,eAll)

   !  Calculate predicted data for dataVectorMTX object d
   !    and for conductivity parameter sigma
   !  Optionally returns array of EM solutions eAll, one for
   !    each transmitter (if eAll is present)
   !
   !  First need to set up transmitter, receiver, dataType dictionaries;
   !  "pointers" to dictionaries are attached to multi-transmitter
   !   data vector d .
   !
   !   sigma is input conductivity parameter
   type(modelParam_t), intent(in)	:: sigma
   !   d is the computed (output) data vector, also used to identify
   !     receiver/transmitter
   type(dataVectorMTX_t), intent(inout)	:: d
   !  structure containing array of solution vectors (should be
   !   allocated before calling)
   type(solnVectorMTX_t), intent(inout), optional	:: eAll

   ! local variables
   integer				:: j

   if(.not.d%allocated) then
      call errStop('data vector not allocated on input to fwdPred')
   end if

   if(present(eAll)) then
      if(.not. eAll%allocated) then
         call create_solnVectorMTX(d%nTx,eAll)
      else if(d%nTx .ne. eAll%nTx) then
         call errStop('dimensions of eAll and d do not agree in fwdPred')
      endif
   endif

   ! loop over transmitters: solve forward system for each,
   !    apply (non-linear) data functionals
   do j = 1,d%nTx

      call fwdPred_TX(sigma,d%d(j),e0)

      if(present(eAll)) then
         call copy_solnVector(eAll%solns(j),e0)
      endif

   enddo

  end subroutine fwdPred


  !**********************************************************************
  subroutine setGrid(newgrid)

   !  Use to set and/or update the numerical grid, that is then used
   !   all computations in this module;
   !   This is not a pointer target.
   !  Might also have to run exitSolver at this point, if we are updating
   !   the grid during an inversion; that restarts the ForwardSolver module.

   type(grid_t), intent(inout)     :: newgrid

   grid = newgrid

  end subroutine setGrid


  !**********************************************************************
  subroutine cleanUp()

   ! Subroutine to deallocate all memory stored in this module

   call exitSolver(e0,e,comb)
   call deall_grid(grid)

  end subroutine cleanUp

end module SensComp
