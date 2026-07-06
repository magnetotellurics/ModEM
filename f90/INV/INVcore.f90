module INVcore

!use math_constants
!use utilities
use senscomp
use dataio
use DistortionParam
use DistortGradient
use DataSens, only: set_dataSens_distortion_ptr

#ifdef MPI
  use Main_MPI
  use Sub_MPI
#endif

   ! inherits datasens,  dataspace, dataFunc, SolnSpace,
   !            modelspace, soln2d

implicit none

public          :: printf, func, gradient
public          :: CdInvMult, CmSqrtMult
public          :: func_dist, gradient_dist, psi_C
public          :: reComputedHat_dist


Contains

!**********************************************************************
   subroutine printf(comment,lambda,alpha,f,mNorm,rms,logfile)

   ! Compute the full penalty functional F
   ! Also output the predicted data and the EM solution
   ! that can be used for evaluating the gradient
   ! Assuming that the model norm is already scaled by Nmodel

    character(*), intent(in)               :: comment
    real(kind=prec), intent(in)  :: lambda, alpha, f, mNorm, rms
    character(*), intent(in), optional	:: logfile
    integer  :: io_unit, ios
    logical  :: opened

    if (present(logfile)) then
    	io_unit = ioLog
    	inquire(file=logfile,opened=opened)
    	if (.not. opened) then
    		open (unit=ioLog,file=logfile,status='unknown',position='append',iostat=ios)
    	end if
    else
    	io_unit = 6
    end if

	write(io_unit,'(a10)',advance='no') trim(comment)//':'
	write(io_unit,'(a3,es12.6)',advance='no') ' f=',f
	write(io_unit,'(a4,es12.6)',advance='no') ' m2=',mNorm
	write(io_unit,'(a5,f11.6)',advance='no') ' rms=',rms
	write(io_unit,'(a8,es12.6)',advance='no') ' lambda=',lambda
	write(io_unit,'(a7,es12.6)') ' alpha=',alpha

	! flush(io_unit): this has the effect of flushing the buffer
	if (present(logfile)) then
		close(io_unit)
		open (unit=ioLog,file=logfile,status='old',position='append',iostat=ios)
	end if

   end subroutine printf


!**********************************************************************
   subroutine func(lambda,d,m0,mHat,F,mNorm,dHat,eAll,RMS)

   ! Compute the full penalty functional F
   ! Also output the predicted data and the EM solution
   ! that can be used for evaluating the gradient

   real(kind=prec), intent(in)  :: lambda
   type(dataVectorMTX_t), intent(in)              :: d
   type(modelParam_t), intent(in)           :: m0
   type(modelParam_t), intent(in)           :: mHat
   real(kind=prec), intent(out) :: F, mNorm
   type(dataVectorMTX_t), optional, intent(inout)   :: dHat
   type(solnVectorMTX_t), optional, intent(inout) :: eAll
   real(kind=prec), optional, intent(out) :: RMS

   !  local variables
   type(dataVectorMTX_t)    :: res,Nres
   type(modelParam_t) :: m,JTd
   real(kind=prec) :: SS
   integer :: Ndata, Nmodel

   ! compute the smoothed model parameter vector
   call CmSqrtMult(mHat,m)

   ! overwriting input with output
   call linComb(ONE,m,ONE,m0,m)

   ! initialize dHat
   dHat = d

   !  compute predicted data for current model parameter m
   !   also sets up forward solutions for all transmitters in eAll
   !   (which is created on the fly if it doesn't exist)
#ifdef MPI
      call Master_Job_fwdPred(m,dHat,eAll)
#else
      call fwdPred(m,dHat,eAll)
#endif

!	call write_Z_ascii(fidWrite,cfile,nPer,periods,modes, &
!			nSites,sites,allData)

   ! initialize res
   res = d

   ! compute residual: res = d-dHat
   call linComb(ONE,d,MinusONE,dHat,res)

   ! normalize residuals, compute sum of squares
   call CdInvMult(res,Nres)
   SS = dotProd(res,Nres)
   Ndata = countData(res)

   ! compute the model norm
   mNorm = dotProd(mHat,mHat)
   Nmodel = countModelParam(mHat)

   ! penalty functional = sum of squares + scaled model norm
   F = SS/Ndata + (lambda * mNorm/Nmodel)

   ! scale mNorm for output
   mNorm = mNorm/Nmodel

   ! if required, compute the Root Mean Squared misfit
   if (present(RMS)) then
       RMS = sqrt(SS/Ndata)
   end if

   call deall_dataVectorMTX(res)
   call deall_dataVectorMTX(Nres)
   call deall_modelParam(m)
   call deall_modelParam(JTd)

   end subroutine func

!**********************************************************************
   subroutine gradient(lambda,d,m0,mHat,grad,dHat,eAll)

   !  Computes the gradient of the penalty functional,
   !  using EM solution (eAll) and the predicted data (dHat)
   !  Here, mHat denotes the non-regularized model parameter that
   !  is normally referred to as \tilde{m} = C_m^{-1/2}(m - m_0),
   !  and the gradient is computed with respect to \tilde{m}.
   !  Before calling this routine, the forward solver must be run:
   !  call CmSqrtMult(mHat,m)
   !  call linComb(ONE,m,ONE,m0,m)
   !  call fwdPred(m,dHat,eAll)

   real(kind=prec), intent(in)  :: lambda
   type(dataVectorMTX_t), intent(in)              :: d
   type(modelParam_t), intent(in)           :: m0
   type(modelParam_t), intent(in)           :: mHat
   type(modelParam_t), intent(inout)          :: grad
   type(dataVectorMTX_t), intent(inout)              :: dHat
   type(solnVectorMTX_t), intent(inout)            :: eAll

   !  local variables
   real(kind=prec)       :: Ndata,Nmodel
   type(dataVectorMTX_t)    :: res
   type(modelParam_t) :: m,JTd,CmJTd

   ! integer :: j, Ny, NzEarth

   ! compute the smoothed model parameter vector
   call CmSqrtMult(mHat,m)

   ! overwriting the input with output
   call linComb(ONE,m,ONE,m0,m)

   ! initialize res
   res = d

   ! compute residual: res = (d-dHat)/Ndata
   call linComb(ONE,d,MinusONE,dHat,res)

   ! multiply by J^T
   call CdInvMult(res)
#ifdef MPI
        call Master_job_JmultT(m,res,JTd,eAll)
#else
#ifdef MPI
    call Master_job_JmultT(m,res,JTd,eAll)
#else
    call JmultT(m,res,JTd,eAll)
#endif
#endif

   call CmSqrtMult(JTd,CmJTd)

   ! initialize grad
   grad = m

   ! compute the number of data and model parameters for scaling
   Ndata = countData(res)
   Nmodel = countModelParam(mHat)

   ! multiply by 2 (to be consistent with the formula)
   ! and add the gradient of the model norm
   call linComb(MinusTWO/Ndata,CmJTd,TWO*lambda/Nmodel,mHat,grad)

   call deall_dataVectorMTX(res)
   call deall_modelParam(m)
   call deall_modelParam(JTd)
   call deall_modelParam(CmJTd)
   !call deall(eAll)
   end subroutine gradient

!**********************************************************************
   subroutine CdInvMult(d_in,d_out)

   ! Divides by the data covariance C_d, which is a diagonal
   ! operator. Divides by the variances (squared error bars)
   ! and scales by the number of data (degrees of freedom).

   type(dataVectorMTX_t), intent(inout)           :: d_in
   type(dataVectorMTX_t), optional, intent(out)   :: d_out
   type(dataVectorMTX_t)                          :: d
   !integer                                :: Ndata

    d = d_in

    ! divide each data component by its variance
    call normalize_dataVectorMTX(d,2)

    ! divide by the number of data
    !Ndata = countData(d)
    !d = scMult(ONE/Ndata,d)

   	if (present(d_out)) then
   		d_out = d
   	else
   	    d_in = d
   	end if

   	call deall(d)

   end subroutine CdInvMult


!**********************************************************************
   subroutine CmSqrtMult(m_in,m_out)

   ! Multiplies by the square root of the model covariance,
   ! which is viewed as a smoothing operator. Intended
   ! to be used to compute m = C_m^{1/2} \tilde{m} + m_0.
   ! For efficiency, CmSqrt is a saved, private variable inside
   ! the modelParam module. Before this routine can be called,
   ! it has to be initialized by calling create_CmSqrt(m).
   ! Now that multBy_CmSqrt routine exists in modelParam,
   ! this routine is no longer needed. Leaving it here for now,
   ! to minimize changes.

   type(modelParam_t), intent(in)              :: m_in
   type(modelParam_t), intent(inout)             :: m_out

	! apply the operator Cm^(1/2) here
	! m_out = m_in
	m_out = multBy_CmSqrt(m_in)

   end subroutine CmSqrtMult

!**********************************************************************
   function psi_C(distC) result(psi)
   type(distortionParam_t), intent(in) :: distC
   real(kind=prec) :: psi
   integer :: i
   psi = R_ZERO
   do i = 1, distC%nSites
      psi = psi + (distC%C(1,1,i) - ONE)**2 &
                + distC%C(1,2,i)**2 &
                + distC%C(2,1,i)**2 &
                + (distC%C(2,2,i) - ONE)**2
   end do
   psi = psi * 0.5_prec
   end function psi_C

!**********************************************************************
   subroutine func_dist(lambda, nu, d, m0, mHat, distC, F, mNorm, distNorm, dHat, eAll, RMS)
   ! Compute the full penalty functional F, including distortion regularization
   ! Also output the predicted data 

   real(kind=prec), intent(in)              :: lambda, nu
   type(dataVectorMTX_t), intent(in)        :: d
   type(modelParam_t), intent(in)           :: m0
   type(modelParam_t), intent(in)           :: mHat
   type(distortionParam_t), intent(in)      :: distC
   real(kind=prec), intent(out)             :: F, mNorm, distNorm
   type(dataVectorMTX_t), optional, intent(inout) :: dHat
   type(solnVectorMTX_t), optional, intent(inout) :: eAll
   real(kind=prec), optional, intent(out)         :: RMS

   type(dataVectorMTX_t)                    :: res,Nres
   type(modelParam_t)                       :: m, JTd
   real(kind=prec)                          :: SS
   integer                                  :: Ndata, Nmodel

   call CmSqrtMult(mHat,m)
   call linComb(ONE,m,ONE,m0,m)

   dHat = d
   call set_distortionParam_ptr(distC)
   call set_dataSens_distortion_ptr(distC)

#ifdef MPI
   call Master_Job_fwdPred(m,dHat,eAll)
#else
   call fwdPred(m,dHat,eAll)
#endif

   res = d
   call linComb(ONE,d,MinusONE,dHat,res)
   call CdInvMult(res,Nres)
   SS = dotProd(res,Nres)
   Ndata = countData(res)

   mNorm = dotProd(mHat,mHat)
   Nmodel = countModelParam(mHat)

   distNorm = psi_C(distC)

   F = SS/Ndata + lambda * mNorm/Nmodel + nu * distNorm
   mNorm = mNorm/Nmodel

   if (present(RMS)) RMS = sqrt(SS/Ndata)

   call deall_dataVectorMTX(res)
   call deall_dataVectorMTX(Nres)
   call deall_modelParam(m)
   call deall_modelParam(JTd)

   end subroutine func_dist

!**********************************************************************
   subroutine gradient_dist(lambda, nu, d, m0, mHat, distC, gradM, gradC, dHat, eAll)
   ! calculate the gradient of the penalty functional with respect to both 
   ! model parameters and distortion parameters
   ! note that full eAll is required for the distortion gradient calculation, 
   ! so this routine should be called after func_dist

   real(kind=prec), intent(in)             :: lambda, nu
   type(dataVectorMTX_t), intent(in)       :: d
   type(modelParam_t), intent(in)          :: m0
   type(modelParam_t), intent(in)          :: mHat
   type(distortionParam_t), intent(in)     :: distC
   type(modelParam_t), intent(inout)       :: gradM
   type(distortionParam_t), intent(inout)  :: gradC
   type(dataVectorMTX_t), intent(inout)    :: dHat
   type(solnVectorMTX_t), intent(inout)    :: eAll

   type(dataVectorMTX_t)                   :: res
   type(modelParam_t)                      :: m, JTd, CmJTd
   real(kind=prec)                         :: Ndata, Nmodel

     call CmSqrtMult(mHat,m)
     call linComb(ONE,m,ONE,m0,m)

     call set_distortionParam_ptr(distC)
     call set_dataSens_distortion_ptr(distC)

     res = d
     call linComb(ONE,d,MinusONE,dHat,res)
     call CdInvMult(res)
#ifdef MPI
     call Master_job_JmultT(m,res,JTd,eAll)
#else
     call JmultT(m,res,JTd,eAll)
#endif
      call CmSqrtMult(JTd,CmJTd)

      gradM = m
      Ndata = countData(res)
      Nmodel = countModelParam(mHat)
      call linComb(MinusTWO/Ndata,CmJTd,TWO*lambda/Nmodel,mHat,gradM)

      call compDistGrad(d, eAll, m, distC, nu, Ndata, gradC)

     call deall_dataVectorMTX(res)
     call deall_modelParam(m)
     call deall_modelParam(JTd)
     call deall_modelParam(CmJTd)

   end subroutine gradient_dist

!**********************************************************************
   subroutine reComputedHat_dist(d, eAll, sigma, distC, dHat, RMS)
   
   ! recompute the  predicted data and RMS after C-update, without 
   ! re-running the full forward solve, using locked EM fields (eAll) 
   ! Useful when only distortion (distC) changes, (conductivity model 
   ! sigma is fixed)
   ! currently used in NLCG_DIST inversion method only.
   ! put it here as a separate subroutine for future use in other 
   ! inversion methods

   type(dataVectorMTX_t), intent(in)  :: d      ! observed data
   type(solnVectorMTX_t), intent(in)  :: eAll   ! locked EM field solutions
   type(modelParam_t), intent(in)     :: sigma  ! conductivity model (smoothed)
   type(distortionParam_t), intent(in):: distC  ! distortion matrices
   type(dataVectorMTX_t), intent(out) :: dHat   ! updated predicted data
   real(kind=prec), intent(out)       :: RMS    ! output misfit RMS

   ! local variables
   type(dataVectorMTX_t) :: res, Nres
   real(kind=prec) :: SS
   integer :: Ndata
   integer :: j, i, iSite, iRx, iTx, iDt, nSites, iDistSite
   real(kind=prec) :: C_site(2,2)

   ! Set up distortion for data sensitivity calculations
   call set_distortionParam_ptr(distC)
   call set_dataSens_distortion_ptr(distC)

   ! Initialize predicted data to observed (will be overwritten)
   dHat = d

   ! Iterate through data blocks and recompute predictions using frozen eAll
   ! This follows the pattern of compDistGrad
   do j = 1, d%nTx
      iTx = d%d(j)%tx
      do i = 1, d%d(j)%nDt
         ! note that predicted data should not carry observational error bars.
         dHat%d(j)%data(i)%errorBar = .false.

         iDt = d%d(j)%data(i)%dataType
         
         ! Only process Full_Impedance_Dist data blocks
         if (iDt == Full_Impedance_Dist) then
            nSites = d%d(j)%data(i)%nSite
            
            do iSite = 1, nSites
               iRx = d%d(j)%data(i)%rx(iSite)
               
               ! Map receiver to distortion site
               iDistSite = find_dist_site(distC, iRx)
               if (iDistSite < 1) cycle
               
               ! Extract 2x2 distortion matrix for this site
               C_site(1,1) = distC%C(1,1,iDistSite)
               C_site(1,2) = distC%C(1,2,iDistSite)
               C_site(2,1) = distC%C(2,1,iDistSite)
               C_site(2,2) = distC%C(2,2,iDistSite)
               
               ! Call dataResp_dist with frozen eAll and site-specific C matrix
               ! This applies C*Z multiplication internally
               call dataResp_dist(eAll%solns(j), sigma, C_site, iDt, iRx, &
                                  dHat%d(j)%data(i)%value(:,iSite), &
                                  d%d(j)%data(i)%orient(iSite))
            end do
         end if
      end do
   end do

   ! Compute residual and misfit RMS
   res = d
   call linComb(ONE, d, MinusONE, dHat, res)
   call CdInvMult(res, Nres)
   SS = dotProd(res, Nres)
   Ndata = countData(res)

   RMS = sqrt(SS / Ndata)

   call deall_dataVectorMTX(res)
   call deall_dataVectorMTX(Nres)

   end subroutine reComputedHat_dist

end module INVcore
