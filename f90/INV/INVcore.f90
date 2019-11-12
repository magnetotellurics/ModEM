module INVcore

!use math_constants
!use utilities
use senscomp
use dataio

#ifdef MPI
  use Main_MPI
  use Sub_MPI
#endif

   ! inherits datasens,  dataspace, dataFunc, SolnSpace,
   !            modelspace, soln2d

implicit none

public          :: printf 
public          :: func, gradient
public          :: func2, gradient2 ! in untransformed modelspace
public          :: CdInvMult
public          :: CmInvMult, CmSqrtMult ! used only in untransformed modelspace


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
        call JmultT(m,res,JTd,eAll)
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
   subroutine func2(lambda,d,m0,m,F,mNorm,dHat,eAll,RMS)

   ! Compute the full penalty functional F 
   ! for the "untransformed" model space and transformed data space
   ! Also output the predicted data and the EM solution
   ! that can be used for evaluating the gradient

   real(kind=prec), intent(in)                      :: lambda
   type(dataVectorMTX_t), intent(in)                :: d
   type(modelParam_t), intent(in)                   :: m0
   type(modelParam_t), intent(in)                   :: m
   real(kind=prec), intent(out)                     :: F, mNorm
   type(dataVectorMTX_t), optional, intent(inout)   :: dHat
   type(solnVectorMTX_t), optional, intent(inout)   :: eAll
   real(kind=prec), optional, intent(out)           :: RMS

   !  local variables
   type(dataVectorMTX_t)                            :: res,Nres
   real(kind=prec)                                  :: SS
   type(modelParam_t)                               :: mHat,m_minus_m0
   integer                                          :: Ndata, Nmodel

   ! compute the smoothed model parameter vector
   ! no longer useful in "untransformed" model domain
   ! call CmSqrtMult(mHat,m)

   ! overwriting input with output
   ! call linComb(ONE,m,ONE,m0,m)

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

   ! initialize res
   res = d

   ! compute residual: 
   ! res = d - f(m)
   call linComb(ONE,d,MinusONE,dHat,res)
   ! and in the transformed space: 
   ! \tilde(res) = C_d^(-1)*res
   call CdInvMult(res,Nres)
   ! we should have used a "CdSqrtInvMult" if we have one!
   ! ss = \tilde{res}^T*\tilde{res}, residue norm
   SS = dotProd(res,Nres)
   Ndata = countData(res)

   ! compute the model norm, still using the mHat convention
   call linComb(ONE, m, MinusONE, m0, m_minus_m0)
   call CmInvMult(m_minus_m0,mHat)
   
   ! mNorm = (m-m0)^T*C_m^(-1)*(m-m0), model norm
   mNorm = dotProd(m_minus_m0,mHat)
   Nmodel = countModelParam(mHat)

   ! scale mNorm for output
   mNorm = mNorm/Nmodel

   ! penalty functional = sum of squares + scaled model norm
   F = SS/Ndata + (lambda * mNorm)

   ! if required, compute the Root Mean Squared misfit
   if (present(RMS)) then
       RMS = sqrt(SS/Ndata)
   end if

   call deall_dataVectorMTX(res)
   call deall_dataVectorMTX(Nres)
   call deall_modelParam(mHat)
   call deall_modelParam(m_minus_m0)
   ! call deall_modelParam(JTd)

   end subroutine func2

!**********************************************************************
   subroutine gradient2(lambda,d,m0,m,grad,dHat,eAll)

   !  Computes the gradient of the penalty functional in the "untransformed"
   !  model space and tranformed data space
   !  using EM solution (eAll) and the predicted data (dHat)
   !  NOTE: the gradient here is computed with respect to m (as in normal 
   !  inversions), instead of \tilde{m} = C_m^{-1/2}(m - m_0),
   !  which is in Gary's tranformed function space
   !  Before calling this routine, the following subroutines should be called 
   !  call CmSqrtMult(mHat,m) ! no longer needed
   !  call linComb(ONE,m,ONE,m0,m)  ! no longer needed
   !  call fwdPred(m,dHat,eAll)

   real(kind=prec), intent(in)                      :: lambda
   type(dataVectorMTX_t), intent(in)                :: d
   type(modelParam_t), intent(in)                   :: m0
   type(modelParam_t), intent(in)                   :: m
   type(modelParam_t), intent(inout)                :: grad
   type(dataVectorMTX_t), intent(inout)             :: dHat
   type(solnVectorMTX_t), intent(inout)             :: eAll

   !  local variables
   real(kind=prec)                                  :: Ndata,Nmodel
   type(dataVectorMTX_t)                            :: res
   type(modelParam_t)                               :: mHat,JTd,CmJTd
   type(modelParam_t)                               :: m_minus_m0

   ! integer :: j, Ny, NzEarth

   ! compute the smoothed model parameter vector
   ! no longer useful in "normal" model domain
   ! call CmSqrtMult(mHat,m)
   ! overwriting the input with output
   ! call linComb(ONE,m,ONE,m0,m)

   ! initialize res
   res = d

   ! compute residual: res = (d-dHat)/Ndata
   call linComb(ONE,d,MinusONE,dHat,res)
   !
   ! now the first part of the gradient w.r.t untransformed model
   ! apply for the *full* data covariance
   call CdInvMult(res)
   ! multiply by J^T: 
#ifdef MPI
        call Master_job_JmultT(m,res,JTd,eAll)
#else
        call JmultT(m,res,JTd,eAll)
#endif
   ! apply for the *half* model covariance
   ! call CmSqrtMult(JTd,CmJTd)

   ! compute the number of data and model parameters for scaling
   Ndata = countData(res)

   ! now the second part of the gradient w.r.t untransfomred model
   ! mHat is now used to store the second part of the gradient -
   ! although it is sort of playing the same role as in the transformed 
   ! function space!

   m_minus_m0 = m
   call linComb(ONE, m, MinusONE, m0, m_minus_m0)
   ! apply for the whole model covariance
   call CmInvMult(m_minus_m0, mHat)
   Nmodel = countModelParam(mHat)

   grad = m
   ! multiply by 2 (to be consistent with the formula)
   ! well, think of (x^2)' = 2x...
   ! and add the gradient of the model norm
   ! grad = -2 J^T C_d^(-1) r + 2 lambda C_m^(-1) m 
   ! in Gary's roughed domain, it was like:
   ! grad = -2 C_m^(0.5) J^T C_d^(-1) r + 2 lambda C_m^(-0.5) m 
   call linComb(MinusTWO/Ndata,JTd,TWO*lambda/Nmodel,mHat,grad)
   !
   call deall_dataVectorMTX(res)
   call deall_modelParam(mHat)
   call deall_modelParam(m_minus_m0)
   call deall_modelParam(JTd)
   call deall_modelParam(CmJTd) 
   !call deall(eAll)
   end subroutine gradient2

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
   type(modelParam_t), intent(out)             :: m_out

	! apply the operator Cm^(1/2) here
	! m_out = m_in
	m_out = multBy_CmSqrt(m_in)

   end subroutine CmSqrtMult

!**********************************************************************
   subroutine CmInvMult(m_in,m_out)

   ! Multiplies by the inverse of the full model covariance,
   ! which is viewed as a roughing operator. Intended
   ! to be used to compute C_m (m - m_0).
   ! For efficiency, CmSqrt is a saved, private variable inside
   ! the modelParam module. Before this routine can be called,
   ! it has to be initialized by calling create_CmSqrt(m).
   ! Now that multBy_CmSqrt routine exists in modelParam,
   ! this routine is no longer needed. Leaving it here for now,
   ! to minimize changes.

   type(modelParam_t), intent(in)              :: m_in
   type(modelParam_t), intent(out)             :: m_out

   ! apply the operator Cm^(-1) here
   m_out = multBy_CmInv(m_in)

   end subroutine CmInvMult

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

end module INVcore
