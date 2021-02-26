module DataSens
!  higher level module that calls routines in DataFunc to
!   --> evaluate non-linear and linear data functionals for all sites
!       for a dataVector object
!   --> construct the comb for forcing the adjoint model for the linear
!       combination of data functionals corresponding to a dataVector object
!   This module is generic, and should work for a broad range of different
!     problems (in particular, 2D and 3D MT).  Much of the code in these
!     routines is for handling two general efficiency issues that arise
!     with frequency domain EM data: Even at a single site there can be multiple
!     components that can be evaluated from a single solnVector object: e.g.,
!     in 2D MT, TE mode impedance and Tipper; or for 3D 4 components of
!     impedance tensor; and data can be complex or real.  These modules
!     handle these different cases efficiently, so that these details
!     are hidden from manipulations in the sensitivity calculation module
!   This module also handles the sparse vectors Q, which define the
!     direct sensitivity of the data functionals to model parameters.

!*****************************************************************************
  use math_constants
  use utilities
  use DataSpace
  use DataFunc

implicit none

 public 		:: Lmult,  LmultT, Qmult, QmultT

Contains

  subroutine Lmult(e0,Sigma0,ef,d)
  ! given the background model parameter (Sigma0) and both
  ! measured and background electric field solutions (ef,e0)
  ! evaluate linearized data functionals for all sites and data types
  ! found in a dataVector object.
  !  implements d = L x ef.

  !  electric field solutions are stored as type solnVector
  type (solnVector_t), intent(in)			:: ef,e0
  ! d provides indices into receiver and data type dictionaries on
  ! input.  Predicted impedances are computed using
  ! these and the input electric field solutions
  type (dataVector_t), intent(inout)      :: d
  !  background model parameter
  type (modelParam_t), intent(in)	        :: Sigma0

  !  local variables
  logical               	:: Conj_Case = .false.
  integer               	:: iSite, ncomp, iTx, &
					nFunc, iDt, iRx, iComp, iFunc, j
	logical                 :: isComplex, exists
  complex(kind=prec)	    :: Z
  type(sparseVector_t), pointer	:: Lz(:)

  !  first check consistency of e0, ef, d
  !             ... all should point to same transmitter!
  if((d%tx.ne.e0%tx) .or. (d%tx.ne.ef%tx) .or. (e0%tx.ne.ef%tx)) then
     call errStop('transmitter incompatability in Lmult')
  endif

  iTx = d%tx

  do j = 1,d%ndt
    iDt = d%data(j)%dataType
	  ncomp = d%data(j)%nComp
	  isComplex = d%data(j)%isComplex
	  if(isComplex) then
	     !  data are complex; one sensitivity calculation can be
	     !   used for both real and imaginary parts
	     if(mod(ncomp,2).ne.0) then
	        call errStop('for complex data # of components must be even in Lmult')
	     endif
	     nFunc = ncomp/2
	  else
	     !  data are treated as real: full sensitivity computation is required
	     !   for each component
	     nFunc = ncomp
	  endif

	  ! allocate space for sparse rows of L
	  allocate(Lz(nFunc))
	  do iFunc=1,nFunc
		  call create_sparseVector(e0%grid,iTx,Lz(iFunc))
	  end do

    ! make sure that predicted data do not have error bars
    d%data(j)%errorBar = .false.

	  !  loop over sites
	  do iSite = 1,d%data(j)%nSite
	     ! compute sparse vector representations of linearized
	     ! data functionals for transfer function elements
	     ! at one site, for one data type.
	     !   Lrows returns one Lz for each of nFunc functionals
	     iRx = d%data(j)%rx(iSite)
	     call Lrows(e0,Sigma0,iDt,iRx,Lz)
	     iComp = 1
	     do iFunc  = 1, nFunc
	        exists = d%data(j)%exist(iComp,iSite)
	        if(exists) then
	           Z = dotProd_sparseVsolnV(Lz(iFunc),ef,Conj_Case)
	        else
	           Z = C_ZERO
	        endif
	        if(isComplex) then
	           d%data(j)%value(iComp,iSite) = real(Z)
	           iComp = iComp + 1
	           d%data(j)%value(iComp,iSite) = dimag(Z)
	           iComp = iComp + 1
	        else
	           d%data(j)%value(iComp,iSite) = real(Z)
	           iComp = iComp + 1
	        endif
	     enddo ! iFunc
	  enddo ! iSite
	  !  deallocate local arrays
	  do iFunc = 1,nFunc
	     call deall_sparseVector(Lz(iFunc))
	  enddo
	  deallocate(Lz)
  enddo ! j

  end subroutine Lmult

!*****************************************************************************
  subroutine LmultT(e0,Sigma0,d,comb)

  ! given background solution vector e0
  ! and a dataVector object d (element of data space containing data
  ! for one frequency, one or more sites and data types) compute adjoint
  ! of measurement operator: i.e., the comb constructed from the scaled
  ! superposition of data kernals (scaled by conjugate of data values ...
  !  e.g., if d contains residuals, this can be used to set up for
  !  gradient calculation).
  !  NOTE: we are only supporting full storage sources in comb;
  !    the elements of comb should be allocated and zeroed before calling.
  !  implements comb = L^T x d.
  !   multiplies each row of L by respective data component and sums up.

  ! background model parameter
  type (modelParam_t), intent(in)  :: Sigma0
  !  electric field solutions are stored as type solnVector
  type (solnVector_t), intent(in)     :: e0
  ! d provides indices into receiver dictionary on
  ! input.  Predicted impedances are computed using
  ! these and the input electric field solutions
  type (dataVector_t), intent(in)               :: d
  ! Output
  type (rhsVector_t), intent(inout)          		:: comb

  !  local variables
  integer               		:: iSite, iTx, iDt, iRx, nComp, &
					 nFunc, iComp, iFunc, j
  logical                   :: isComplex, exists
  complex(kind=prec)		    :: Z
  type(sparseVector_t), pointer    	:: Lz(:)

  iTx = d%tx

  call zero_rhsVector(comb)

  do j = 1,d%ndt
	  iDt = d%data(j)%dataType
	  ncomp = d%data(j)%nComp
	  isComplex = d%data(j)%isComplex
	  if(isComplex) then
	     !  data are complex; one sensitivity calculation can be
	     !   used for both real and imaginary parts
	     if(mod(ncomp,2).ne.0) then
	        call errStop('for complex data # of components must be even in LmultT')
	     endif
	     nFunc = ncomp/2
	  else
	     !  data are treated as real: full sensitivity computation is required
	     !   for each component
	     nFunc = ncomp
	  endif

	  ! allocate space for sparse rows of L and zero the output
	  allocate(Lz(nFunc))
	  do iFunc=1,nFunc
	      call create_sparseVector(e0%grid,iTx,Lz(iFunc))
	  end do

	  !  loop over sites
	  do iSite = 1,d%data(j)%nSite
	     ! compute sparse vector representations of linearized
	     ! data functionals for transfer function elements
	     ! at one site, for one data type.
	     !   Lrows returns one Lz for each of nFunc functionals
	     iRx = d%data(j)%rx(iSite)
	     call Lrows(e0,Sigma0,iDt,iRx,Lz)
	     iComp = 1
	     do iFunc  = 1, nFunc
	        exists = d%data(j)%exist(iComp,iSite)
	        if(isComplex) then
	           !  move real data from dataVector into complex conjugate of TF (impedance)
	           !  multiply this by data kernel for complex impedance ...
	           !      (take real part in parameter space)
	           Z = cmplx(d%data(j)%value(iComp,iSite),-d%data(j)%value(iComp+1,iSite),8)
	           iComp = iComp+2
	        else
	           Z = cmplx(d%data(j)%value(iComp,iSite),0.0,8)
	           iComp = iComp+1
	        endif
	        if(exists) then
	           call add_sparseVrhsV(Z,Lz(iFunc),comb)
	        endif
	     enddo ! iFunc
	  enddo ! iSite
	  !  deallocate local arrays
	  do iFunc = 1,nFunc
	     call deall_sparseVector(Lz(iFunc))
	  enddo
	  deallocate(Lz)
  enddo ! j

  end subroutine LmultT

!****************************************************************************

  subroutine Qmult(e0,Sigma0,dSigma,d)
  ! given the background model parameter (Sigma0) and the
  ! background electric field solution (e0)
  ! evaluate linearized data functionals with respect to model parameters
  ! for all sites and all data types in a dataVector object.
  !  implements d = Q x dSigma.

  ! background electric field solution
  type (solnVector_t), intent(in)           :: e0
  ! d provides indices into receiver dictionary on
  ! input.  Predicted impedances are computed using
  ! these and the input electric field solutions
  type (dataVector_t), intent(inout)        :: d
  !  background model parameter
  type (modelParam_t), intent(in)           :: Sigma0
  ! conductivity perturbation dSigma to multiply
  type (modelParam_t), intent(in)           :: dSigma

  !  local variables
  integer                   :: iSite, ncomp, iTx, &
                    nFunc, iDt, iRx, iComp, iFunc, j, istat
  logical                   :: isComplex, zeroQ, exists
  real(kind=prec)           :: Zreal, Zimag
  type (modelParam_t), pointer  :: sigmaQreal(:), sigmaQimag(:)

  !  first check consistency of e0 and d
  !             ... they should point to same transmitter!
  if(d%tx.ne.e0%tx) then
     call errStop('transmitter incompatability in Qmult')
  endif


  iTx = d%tx

  do j = 1,d%ndt
	  iDt = d%data(j)%dataType
	  ncomp = d%data(j)%nComp
	  isComplex = d%data(j)%isComplex
	  if(isComplex) then
	     !  data are complex; one sensitivity calculation can be
	     !   used for both real and imaginary parts
	     if(mod(ncomp,2).ne.0) then
	        call errStop('for complex data # of components must be even in Qmult')
	     endif
	     nFunc = ncomp/2
	  else
	     !  data are treated as real: full sensitivity computation is required
	     !   for each component
	     nFunc = ncomp
	  endif

	  ! allocate space for full rows of Q
	  allocate(sigmaQreal(nFunc),STAT=istat)
    allocate(sigmaQimag(nFunc),STAT=istat)

    ! make sure that predicted data do not have error bars
    d%data(j)%errorBar = .false.

	  !  loop over sites
	  do iSite = 1,d%data(j)%nSite
 	     ! compute sparse vector representations of linearized
	     ! data functionals for transfer function
	     ! elements at one site and for one data type
	     !   Qrows returns complex sigmaQ for each of nFunc functionals
	     iRx = d%data(j)%rx(iSite)
	     call Qrows(e0,Sigma0,iDt,iRx,zeroQ,sigmaQreal,sigmaQimag)
	     iComp = 1
	     do iFunc  = 1, nFunc
            exists = d%data(j)%exist(iComp,iSite)
	        if(isComplex) then
	           if (zeroQ .or. (.not. exists)) then
	               Zreal = R_ZERO
	           else
	               Zreal = dotProd_modelParam(sigmaQreal(iFunc),dSigma)
	           endif
	           d%data(j)%value(iComp,iSite) = Zreal
	           iComp = iComp + 1
               if (zeroQ .or. (.not. exists)) then
                   Zimag = R_ZERO
               else
                   Zimag = dotProd_modelParam(sigmaQimag(iFunc),dSigma)
               endif
	           d%data(j)%value(iComp,iSite) = Zimag
	           iComp = iComp + 1
	        else
               if (zeroQ .or. (.not. exists)) then
                   Zreal = R_ZERO
               else
	               Zreal = dotProd_modelParam(sigmaQreal(iFunc),dSigma)
	           endif
	           d%data(j)%value(iComp,iSite) = Zreal
	           iComp = iComp + 1
	        endif
	     enddo ! iFunc
	  enddo ! iSite
	  !  deallocate local arrays
	  do iFunc = 1, nFunc
        call deall_modelParam(sigmaQreal(iFunc))
        call deall_modelParam(sigmaQimag(iFunc))
	  enddo
	  deallocate(sigmaQreal,STAT=istat)
	  deallocate(sigmaQimag,STAT=istat)
  enddo ! j

  end subroutine Qmult

!*****************************************************************************
  subroutine QmultT(e0,Sigma0,d,QcombReal,QcombImag)

  ! given background solution vector e0
  ! and a dataVector object d (element of data space containing data
  ! for one frequency, one or more sites and data types) compute adjoint
  ! of measurement operator: i.e., the comb constructed from the scaled
  ! superposition of data kernals (scaled by conjugate of data values ...
  !  e.g., if d contains residuals, this can be used to set up for
  !  gradient calculation).
  ! Returns conductivity parameter data structure corresponding
  !   to part of sensitivity due to derivative of data functionals
  !   with respect to model parameters
  !  implements Qcomb = Q^T x d.
  !   multiplies each row of Q by respective data component and sums up.

  ! background model parameter
  type (modelParam_t), intent(in)  :: Sigma0
  !  electric field solutions are stored as type solnVector
  type (solnVector_t), intent(in)     :: e0
  ! d provides indices into receiver dictionary on
  ! input.  Predicted impedances are computed using
  ! these and the input electric field solutions
  type (dataVector_t), intent(in)      :: d
  ! Output
  type (modelParam_t), intent(inout)            :: QcombReal
  type (modelParam_t), intent(inout), optional  :: QcombImag

  !  local variables
  integer                   :: iSite, ncomp, iTx, &
                    nFunc, iDt, iRx, iComp, iFunc, j, istat
  logical                   :: isComplex, zeroQ, exists
  real(kind=prec)           :: Zreal, Zimag
  type (modelParam_t)       :: sigmaTemp
  type (modelParam_t), pointer  :: sigmaQreal(:), sigmaQimag(:)

  !  first check consistency of e0 and d
  !             ... they should point to same transmitter!
  if(d%tx.ne.e0%tx) then
     call errStop('transmitter incompatability in QmultT')
  endif

  iTx = d%tx

  ! initialize local variable and output
  sigmaTemp = Sigma0
  call zero(sigmaTemp)
  QcombReal = sigmaTemp
  if(present(QcombImag)) then
     QcombImag = sigmaTemp
  endif

  do j = 1,d%ndt
	  iDt = d%data(j)%dataType
	  ncomp = d%data(j)%nComp
	  isComplex = d%data(j)%isComplex
	  if(isComplex) then
	     !  data are complex; one sensitivity calculation can be
	     !   used for both real and imaginary parts
	     if(mod(ncomp,2).ne.0) then
	        call errStop('for complex data # of components must be even in QmultT')
	     endif
	     nFunc = ncomp/2
	  else
	     !  data are treated as real: full sensitivity computation is required
	     !   for each component
	     nFunc = ncomp
	  endif

	  ! allocate space for full rows of Q
	  allocate(sigmaQreal(nFunc),STAT=istat)
    allocate(sigmaQimag(nFunc),STAT=istat)

	  !  loop over sites
	  do iSite = 1,d%data(j)%nSite
	     ! compute sparse vector representations of linearized
	     ! data functionals for transfer function
	     ! elements at one site and for one data type
	     !   Qrows returns one sigmaQ for each of nFunc functionals
	     iRx = d%data(j)%rx(iSite)
	     call Qrows(e0,Sigma0,iDt,iRx,zeroQ,sigmaQreal,sigmaQimag)
	     iComp = 1
	     do iFunc  = 1, nFunc
	        exists = d%data(j)%exist(iComp,iSite)
	        if(isComplex) then
		       ! implement the real variant of Q^T conj(d), where conj(d) is component-wise
		       ! complex conjugate... the data are stored as real, so don't use complex
		       ! multiplication here
		       if ((.not. zeroQ) .and. exists) then
		           Zreal = d%data(j)%value(iComp,iSite)
		           Zimag = d%data(j)%value(iComp+1,iSite)
		           ! Re( Q^T conj(d) ) = Re(Q^T) Re(d) + Im(Q^T) Im(d)
		           call linComb(Zreal,sigmaQreal(iFunc),Zimag,sigmaQimag(iFunc),sigmaTemp)
		           call scMultAdd(ONE,sigmaTemp,QcombReal)
		           ! Im( Q^T conj(d) ) = Im(Q^T) Re(d) - Re(Q^T) Im(d)
		           if(present(QcombImag)) then
		               call linComb(Zreal,sigmaQimag(iFunc),-Zimag,sigmaQreal(iFunc),sigmaTemp)
	                 call scMultAdd(ONE,sigmaTemp,QcombImag)
	               endif
               endif
	           iComp = iComp + 2
	        else
               if ((.not. zeroQ) .and. exists) then
	               Zreal = d%data(j)%value(iComp,iSite)
	               call scMultAdd(Zreal,sigmaQreal(iFunc),QcombReal)
	           endif
	           iComp = iComp + 1
	        endif
	     enddo ! iFunc
	  enddo ! iSite
	  !  deallocate these arrays for each data type
	  do iFunc = 1, nFunc
	    call deall_modelParam(sigmaQreal(iFunc))
	    call deall_modelParam(sigmaQimag(iFunc))
	  enddo
	  deallocate(sigmaQreal,STAT=istat)
	  deallocate(sigmaQimag,STAT=istat)
  enddo ! j
  ! finally, deallocate sigmaTemp
  call deall_modelParam(sigmaTemp)

  end subroutine QmultT

!****************************************************************************

end module DataSens
