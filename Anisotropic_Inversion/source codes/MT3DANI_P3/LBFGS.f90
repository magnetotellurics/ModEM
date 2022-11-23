module mLBFGS

!use math_constants
!use utilities
use senscomp
use dataio

#ifdef MPI
  use MPI_main
  use MPI_sub
#endif

   ! inherits datasens,  dataspace, dataFunc, SolnSpace,
   !            modelspace, soln2d

implicit none


!-----------------------------------------------------------------------------wkp添加lbfgs全局变量

integer, private      :: lp,mp
real(8), private      :: gtol,stpmin,stpmax

private :: lbfgs,lb1,daxpy,mcsrch,mcstep
!private :: ddot    !ddot 使用的是wsBLAS.f90的 代替 lbfgs里的ddot

!-----------------------------------------------------------------------------wkp添加lbfgs全局变量

public  :: LBFGSsolver

! iteration control for the LBFGS solver is initialized once
! and saved in the module to be used by most subroutines

  type  :: LBFGSiterControl_t
     ! maximum number of iterations in one call to iterative solver
     integer					:: maxIter
     ! convergence criteria: return from solver if rms < rmsTol
     real (kind=prec)	:: rmsTol
     ! the condition to identify when the inversion stalls
     real (kind=prec)   :: fdiffTol
     ! initial value of lambda (will not override the LBFGS input argument)
     real (kind=prec)   :: lambda
     ! exit if lambda < lambdaTol approx. 1e-4
     real (kind=prec)   :: lambdaTol
     ! set lambda_i = lambda_{i-1}/k when the inversion stalls
     real (kind=prec)   :: k
     ! the factor that ensures sufficient decrease in the line search
     real (kind=prec)   :: c
     ! restart CG every nCGmax iterations to ensure conjugacy
     integer                    :: nCGmax
     ! restart CG if orthogonality is lost (not necessarily needed)
     ! real (kind=prec)   :: delta ! 0.5
     ! the starting step for the line search
     real (kind=prec)   :: alpha_1
     ! if alpha_{i+1} < alpha_i * k_{alpha}, set alpha_{i+1} = alpha_i/2
     ! real (kind=prec)   :: alpha_k ! 0.1
     ! if alpha_{i+1} - alpha_i < tol_{alpha}, set alpha_{i+1} = alpha_i/2
     ! real (kind=prec)   :: alpha_tol ! 1.0e-2
     ! maximum initial delta mHat (overrides alpha_1)
     real (kind=prec)   :: startdm
     ! optional relaxation parameter (Renormalised Steepest Descent algorithm)
     real (kind=prec)   :: gamma
     ! model and data output file name
     character(80)              :: fname
  end type LBFGSiterControl_t

  type(LBFGSiterControl_t), private, save :: iterControl

Contains

!**********************************************************************
   subroutine set_LBFGSiterControl(iterControl)
   !  at present this is a simple routine to setup for
   !  iteration control; allowing parameters to be read from
   !   a file would be desireable.

   type(LBFGSiterControl_t), intent(inout)	:: iterControl

     ! maximum number of iterations in one call to iterative solver
     iterControl%maxIter = 600
     ! convergence criteria: return from solver if rms < rmsTol
     iterControl%rmsTol  = 1.05
     ! inversion stalls when abs(rms - rmsPrev) < fdiffTol (2e-3 works well)
     iterControl%fdiffTol = 2.0e-3
     ! initial value of lambda (will not override the LBFGS input argument)
     iterControl%lambda = 1.
     ! exit if lambda < lambdaTol approx. 1e-4
     iterControl%lambdaTol = 1.0e-8
     ! set lambda_i = lambda_{i-1}/k when the inversion stalls
     iterControl%k = 10.
     ! the factor that ensures sufficient decrease in the line search >=1e-4
     iterControl%c = 1.0e-4
     ! restart CG every nCGmax iterations to ensure conjugacy
     iterControl%nCGmax = 8
     ! the starting step for the line search
     iterControl%alpha_1 = 20.
     ! maximum initial delta mHat (overrides alpha_1)
     iterControl%startdm = 20.
     ! optional relaxation parameter (Renormalised Steepest Descent algorithm)
     iterControl%gamma = 0.99
     ! model and data output file name
     iterControl%fname = 'Modular'

   end subroutine set_LBFGSiterControl


   ! ***************************************************************************
   ! * read_LBFGSiterControl reads the inverse solver configuration from file

   subroutine read_LBFGSiterControl(iterControl,rFile,fileExists)

	type(LBFGSiterControl_t), intent(inout)	:: iterControl
    character(*), intent(in)		        :: rFile
	logical, intent(out), optional          :: fileExists
    integer									:: ios
	logical                             	:: exists
	character(80)							:: string

    ! Initialize inverse solver configuration

    call set_LBFGSiterControl(iterControl)

    inquire(FILE=rFile,EXIST=exists)
    if (present(fileExists)) then
       fileExists = exists
    end if

    if (.not. exists) then
       return
    else
       write(*,*) 'Reading inverse configuration from file ',trim(rFile)
    end if

    open (unit=ioInvCtrl,file=rFile,status='old',iostat=ios)

    if(ios/=0) then
       write(0,*) 'Error opening file: ', rFile
    end if

    ! This is the list of options specified in the startup file

    read (ioInvCtrl,'(a36,a80)') string,iterControl%fname      ! 反演每次迭代模型及相应输出名
    if (output_level > 2) then
       write (*,*)
       write (*,'(a36,a80)') string,iterControl%fname
    end if
    iterControl%fname = adjustl(iterControl%fname)
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%lambda   ! 初始 正则化因子λ ？
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%lambda
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%k        ! 把当前λ缩小 
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%k
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%startdm
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%startdm
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%fdiffTol  ! 前后两次拟合差 差值 , 用于反应步长过小
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%fdiffTol
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%rmsTol     ! 反演拟合差最小值
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%rmsTol
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%lambdaTol  ! 正则化因子容许最小值
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%lambdaTol
    end if
    read (ioInvCtrl,'(a36,i4)') string,iterControl%maxIter       ! 最大迭代次数
    if (output_level > 2) then
       write (*,'(a36,i4)') string,iterControl%maxIter
       write (*,*)
    end if

    close(ioInvCtrl)
   end subroutine read_LBFGSiterControl


!**********************************************************************
   subroutine printfLBFGS(comment,lambda,alpha,f,mNorm,rms,logfile)

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
!	write(io_unit,'(a5,f11.6)',advance='no') ' rms=',rms
    write(io_unit,'(a5,es15.6)',advance='no') ' rms=',rms
    

	write(io_unit,'(a8,es12.6)',advance='no') ' lambda=',lambda
	write(io_unit,'(a7,es12.6)') ' alpha=',alpha

	! flush(io_unit): this has the effect of flushing the buffer
	if (present(logfile)) then
		close(io_unit)
		open (unit=ioLog,file=logfile,status='old',position='append',iostat=ios)
	end if

   end subroutine printfLBFGS


!**********************************************************************
   subroutine funcLBFGS(lambda,d,m0,mHat,F,mNorm,dHat,eAll,RMS)

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

   ! compute the smoothed model parameter vector  这个 m是 local variables  m= Cm^(1/2) * mHat, 而 mHat = m(~)=Cm^(-1/2)(m-m0)
   call CmSqrtMultLBFGS(mHat,m)

   ! overwriting input with output   这一步应该是 m= Cm^(1/2) * mHat + m0 ,目的是为了恢复正在应该参与正演的模型参数
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

   ! compute residual: res = d-dHat    ! 观测数据与正演响应 残差 ? MinusONE= -1 定义在 math_constants.f90
   call linComb(ONE,d,MinusONE,dHat,res)

   ! normalize residuals, compute sum of squares
   call CdInvMultLBFGS(res,Nres)
   SS = dotProd(res,Nres)
   Ndata = countData(res)

   ! compute the model norm    对于第一次反演正演的时候 ,mHat = m(~)=Cm^(-1/2)(m-m0),所以第一次mHat在理论上是0 所以在Main.f90 dsigma 直接给0
   mNorm = dotProd(mHat,mHat)
   Nmodel = countModelParam(mHat)

   ! penalty functional = sum of squares + scaled model norm  目标函数值 F
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

   end subroutine funcLBFGS

!**********************************************************************
   subroutine gradientLBFGS(lambda,d,m0,mHat,grad,dHat,eAll)

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
   type(dataVectorMTX_t), intent(in)              :: d       ! 输入的观测数据 dataVectorMTX_t结构体定义在DataSpace.f90
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

   ! compute the smoothed model parameter vector  m= Cm^(1/2) * mHat
   call CmSqrtMultLBFGS(mHat,m)

   ! overwriting the input with output            m= Cm^(1/2) * mHat + m0
   call linComb(ONE,m,ONE,m0,m)

   ! initialize res
   res = d

   ! compute residual: res = (d-dHat)/Ndata
   call linComb(ONE,d,MinusONE,dHat,res)

   ! multiply by J^T  
   call CdInvMultLBFGS(res)                ! res 残差 除以 自己的 协方差 

#ifdef MPI
        call Master_job_JmultT(m,res,JTd,eAll)
#else
        call JmultT(m,res,JTd,eAll)
#endif


   call CmSqrtMultLBFGS(JTd,CmJTd)         ! CmJTd = Cm^(1/2) * JTd

   ! initialize grad
   grad = m

   ! compute the number of data and model parameters for scaling
   Ndata = countData(res)
   Nmodel = countModelParam(mHat)

   ! multiply by 2 (to be consistent with the formula)
   ! and add the gradient of the model norm  模型部分 的 梯度
   call linComb(MinusTWO/Ndata,CmJTd,TWO*lambda/Nmodel,mHat,grad)

   call deall_dataVectorMTX(res)
   call deall_modelParam(m)
   call deall_modelParam(JTd)
   call deall_modelParam(CmJTd)
   !call deall(eAll)
   end subroutine gradientLBFGS
   

!**********************************************************************
   subroutine CdInvMultLBFGS(d_in,d_out)

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

   end subroutine CdInvMultLBFGS


!**********************************************************************
   subroutine CmSqrtMultLBFGS(m_in,m_out)

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

   end subroutine CmSqrtMultLBFGS

!**********************************************************************
   subroutine LBFGSsolver(d,lambda,m0,m,fname)

   ! computes inverse solution minimizing penalty functional
   !   for fixed value of regularization parameter, using
   !   a variant of non-linear conjugate gradient search.
   !   Various flavours of the algorithm and of the line search
   !   can be called from this routine
   !
   !  Note about the starting model:
   !  The starting model has to be in the smoothed model space,
   !  i.e. of the form m = C_m^{1/2} \tilde{m} + m_0.
   !  In order to compute \tilde{m} from the starting model,
   !  C_m^{-1/2} has to be implemented. To avoid this issue
   !  altogether, we are always starting from the prior,
   !  with \tilde{m} = 0. However, in general we could also
   !  start with the result of a previous search.

   !  d is data; on output it contains the responses for the inverse model
   type(dataVectorMTX_t), intent(inout)		   :: d
   !  lambda is regularization parameter
   real(kind=prec), intent(inout)  :: lambda
   !  m0 is prior model parameter
   type(modelParam_t), intent(in)		       :: m0
   !  m is solution parameter ... on input m contains starting guess
   type(modelParam_t), intent(inout)	       :: m
   !  flavor is a string that specifies the algorithm to use
   character(*), intent(in), optional        :: fname
   !  initial step size in the line search direction in model units
   real(kind=prec)                           :: startdm
   !  flavor is a string that specifies the algorithm to use
   character(80)                           :: flavor = 'Cubic'

   !  local variables
   type(dataVectorMTX_t)			:: dHat, res
   type(modelParam_t)			:: mHat, m_minus_m0, grad, g, h, gPrev
   !type(LBFGSiterControl_t)			:: iterControl
   real(kind=prec)		:: value, valuePrev, rms, rmsPrev, alpha, beta, gnorm, mNorm, Nmodel
   real(kind=prec)      :: grad_dot_h, g_dot_g, g_dot_gPrev, gPrev_dot_gPrev, g_dot_h
   integer				:: iter, nCG, nLS, nfunc, ios
   logical              :: ok
   character(3)         :: iterChar
   character(100)       :: mFile, mHatFile, gradFile, dataFile, resFile, logFile
   type(solnVectorMTX_t)      :: eAll
   !----前面有些原来NLCG用的结构体，可以不需要了,暂时还没注释前后的无用结构体



   !----wkp add for calling subroutine lbfgs-----------------------------------------------------------------lbfgs相关标量
   !x(ndim),g(ndim),diag(ndim),w(nwork)
   real(8), dimension(:), allocatable ::x_lbfgs(:),g_lbfgs(:),diag_lbfgs(:),w_lbfgs(:),g0_lbfgs(:),g1_lbfgs(:)
   integer :: ndim_lbfgs , nwork_lbfgs
   integer :: msave_lbfgs
   double precision f_lbfgs,eps,xtol,f0_lbfgs,f1_lbfgs      
   integer iprint(2),iflag,icall,m_lbfgs    
   logical diagco
   ! 考虑地形以及不参与反演区域的时候也是这样的吗?包括后面的x_lbfgs的模型及梯度转换??????????????????????????????????????
   ndim_lbfgs  =  nd_lbfgs(m0)   ! 这里少一维,此函数定义在ModelSpace最后(wkp add),完成m0%Ny * m0%NzEarth
   msave_lbfgs = 7
   nwork_lbfgs = ndim_lbfgs*(2*msave_lbfgs+1)+2*msave_lbfgs
   allocate ( x_lbfgs(ndim_lbfgs),g_lbfgs(ndim_lbfgs),diag_lbfgs(ndim_lbfgs),w_lbfgs(nwork_lbfgs) )
   allocate (g0_lbfgs(ndim_lbfgs),g1_lbfgs(ndim_lbfgs))

   lp=6
   mp=6
   gtol=9.0d-01 
   stpmin=1.0d-20
   stpmax=1.0d+20
   !本模块开头把这5个变量定义为了private(wkp add for lbfgs)


      !n=100
      m_lbfgs=5
      iprint(1)= 1
      iprint(2)= 0

      diagco= .false.
      eps= 1.0d-5
      xtol= 1.0d-16
      icall=0
      iflag=0

   !write(*,*) ndim_lbfgs
  !----wkp add for calling subroutine lbfgs-----------------------------------------------------------------lbfgs相关标量

   if (present(fname)) then
      call read_LBFGSiterControl(iterControl,fname,ok)   ! 读取反演用的控制文件 iterControl 是本模块定义的私有 LBFGSiterControl_t 结构体 fname是控制文件名
      if (ok) then                                      ! ok 指示 fname 是否存在
         lambda = iterControl%lambda
      end if
   else
      call set_LBFGSiterControl(iterControl)             ! 如果 反演控制文件不存在 就在这里使用程序的默认值 控制
   end if

   ! initialize the output to log file                  ! 创建反演记录文件
!   logFile = trim(iterControl%fname)//'_LBFGS.log'
!   open (unit=ioLog,file=logFile,status='unknown',position='append',iostat=ios)

   ! initialize the line search
   alpha = iterControl%alpha_1
   startdm = iterControl%startdm

   write(*,'(a41,es8.1)') 'The initial damping parameter lambda is ',lambda
   write(*,'(a55,f12.6)') 'The initial line search step size (in model units) is ',startdm

!   write(ioLog,'(a41,es8.1)') 'The initial damping parameter lambda is ',lambda
!   write(ioLog,'(a55,f12.6)') 'The initial line search step size (in model units) is ',startdm


   ! starting from the prior hardcoded by setting mHat = 0 and m = m0
   ! m = m0
   ! mHat = m0
   ! call zero(mHat)

!----------------------------------------------------------------------------------------------lbfgs调用开始

   ! starting model contains the rough deviations from the prior
   mHat = m

 2016   continue


   !  compute the penalty functional and predicted data
   call funcLBFGS(lambda,d,m0,mHat,value,mNorm,dHat,eAll,rms) 
   
   call get_x_g(mHat,ndim_lbfgs,x_lbfgs)   ! 定义在ModelSpace里,带地形以及不参与反演区域也这样？
   f_lbfgs = value

   ! compute gradient of the full penalty functional
   call gradientLBFGS(lambda,d,m0,mHat,grad,dHat,eAll)

   call get_x_g(grad,ndim_lbfgs,g_lbfgs)   ! 定义在ModelSpace里,带地形以及不参与反演区域也这样？
  

!   if (output_level > 3) then
!     gradFile = trim(iterControl%fname)//'_LBFGS_'//iterChar//'.grt'
!     call write_modelParam(grad,trim(gradFile))
!   end if

!------例子的目标函数f 及 梯度计算
!      f= 0.d0
!      do 3016 j=1,n,2
!        t1= 1.d0-x(j)
!        t2= 1.d1*(x(j+1)-x(j)**2)
!        g(j+1)= 2.d1*t2
!        g(j)= -2.d0*(x(j)*g(j+1)+t1)
!        f= f+t1**2+t2**2
! 3016   continue

 !    n 变量总数, m 需要保存前面多少次的梯度
 !    x 解向量  , f 目标函数值  , g 梯度
 !    call lbfgs(n,m,x,f,g,diagco,diag,iprint,eps,xtol,w,iflag)
 !    wkp为调用顺利,有以下变量名改变
 !    ndim_lbfgs = n 变量总数 
 !    m_lbfgs 需要保存前面多少次的梯度
 !    x_lbfgs 解向量
 !    f_lbfgs 目标函数值
 !    g_lbfgs 梯度
      call lbfgs(ndim_lbfgs,m_lbfgs,x_lbfgs,f_lbfgs,g_lbfgs,diagco,diag_lbfgs,  &
                  iprint,eps,xtol,w_lbfgs,iflag,m0,mHat,rms,lambda,f0_lbfgs,f1_lbfgs,g0_lbfgs,g1_lbfgs)

      if(iflag.le.0) go to 5016

      
      call x_to_mHat(x_lbfgs,ndim_lbfgs,mHat)  ! 定义在ModelSpace里,带地形以及不参与反演区域也这样？
      
      icall=icall + 1
 
!c     we allow at most 2000 evaluations of f and g
      if(icall.gt.2000) go to 5016  
      go to 2016
 5016   continue

 !----------------------------------------------------------------------------------------------lbfgs调用结束
 deallocate ( x_lbfgs,g_lbfgs,diag_lbfgs,w_lbfgs,g0_lbfgs,g1_lbfgs )

 !----------------------------------------------------------------------------------------------lbfgs调用结束 



   ! cleaning up
   call deall_dataVectorMTX(dHat)
   call deall_dataVectorMTX(res)
   call deall_modelParam(mHat)
   call deall_modelParam(m_minus_m0)
   call deall_modelParam(grad)
   call deall_modelParam(g)
   call deall_modelParam(h)
   call deall_modelParam(gPrev)
   call deall_solnVectorMTX(eAll)

   end subroutine LBFGSsolver


!--------------------------------------LBFGS-------------------------------------------------start!
!c     ----------------------------------------------------------------------
!c     this file contains the lbfgs algorithm and supporting routines
!c
!c     ****************
!c     lbfgs subroutine
!c     ****************
!c


!      block data lb2     !给全局变量赋初值
!      integer lp,mp
!      double precision gtol,stpmin,stpmax
!      common /lb3/mp,lp,gtol,stpmin,stpmax
!      data mp,lp,gtol,stpmin,stpmax/6,6,9.0d-01,1.0d-20,1.0d+20/
!      end


      subroutine lbfgs(n,m,x,f,g,diagco,diag,iprint,eps,xtol,w,iflag,m0,mHat,rms,lambda,f0,f1,g0,g1)

      implicit none
!c
      integer n,m,iprint(2),iflag
      double precision x(n),g(n),diag(n),w(n*(2*m+1)+2*m)
      double precision f,eps,xtol
      logical diagco
      double precision f0,f1,g0(n),g1(n)   ! 用于线搜索
     

      double precision rms,lambda,rmsPrev,fprev,gnormPrev  !用于更新λ
      type(modelParam_t), intent(in)		       :: m0
      type(modelParam_t), intent(inout)		   :: mHat
!c
!c        limited memory bfgs method for large scale optimization
!c                          jorge nocedal
!c                        *** july 1990 ***
!c
!c 
!c     this subroutine solves the unconstrained minimization problem
!c 
!c                      min f(x),    x= (x1,x2,...,xn),
!c
!c      using the limited memory bfgs method. the routine is especially
!c      effective on problems involving a large number of variables. in
!c      a typical iteration of this method an approximation hk to the
!c      inverse of the hessian is obtained by applying m bfgs updates to
!c      a diagonal matrix hk0, using information from the previous m steps.
!c      the user specifies the number m, which determines the amount of
!c      storage required by the routine. the user may also provide the
!c      diagonal matrices hk0 if not satisfied with the default choice.
!c      the algorithm is described in "on the limited memory bfgs method
!c      for large scale optimization", by d. liu and j. nocedal,
!c      mathematical programming b 45 (1989) 503-528.
!c 
!c      the user is required to calculate the function value f and its
!c      gradient g. in order to allow the user complete control over
!c      these computations, reverse  communication is used. the routine
!c      must be called repeatedly under the control of the parameter
!c      iflag. 
!c
!c      the steplength is determined at each iteration by means of the
!c      line search routine mcvsrch, which is a slight modification of
!c      the routine csrch written by more' and thuente.
!c 
!c      the calling statement is 
!c 
!c          call lbfgs(n,m,x,f,g,diagco,diag,iprint,eps,xtol,w,iflag)
!c 
!c      where
!c 
!c     n       is an integer variable that must be set by the user to the
!c             number of variables. it is not altered by the routine.
!c             restriction: n>0.
!c 
!c     m       is an integer variable that must be set by the user to
!c             the number of corrections used in the bfgs update. it
!c             is not altered by the routine. values of m less than 3 are
!c             not recommended; large values of m will result in excessive
!c             computing time. 3<= m <=7 is recommended. restriction: m>0.
!c 
!c     x       is a double precision array of length n. on initial entry
!c             it must be set by the user to the values of the initial
!c             estimate of the solution vector. on exit with iflag=0, it
!c             contains the values of the variables at the best point
!c             found (usually a solution).
!c 
!c     f       is a double precision variable. before initial entry and on
!c             a re-entry with iflag=1, it must be set by the user to
!c             contain the value of the function f at the point x.
!c 
!c     g       is a double precision array of length n. before initial
!c             entry and on a re-entry with iflag=1, it must be set by
!c             the user to contain the components of the gradient g at
!c             the point x.
!c 
!c     diagco  is a logical variable that must be set to .true. if the
!c             user  wishes to provide the diagonal matrix hk0 at each
!c             iteration. otherwise it should be set to .false., in which
!c             case  lbfgs will use a default value described below. if
!c             diagco is set to .true. the routine will return at each
!c             iteration of the algorithm with iflag=2, and the diagonal
!c              matrix hk0  must be provided in the array diag.
!c 
!c 
!c     diag    is a double precision array of length n. if diagco=.true.,
!c             then on initial entry or on re-entry with iflag=2, diag
!c             it must be set by the user to contain the values of the 
!c             diagonal matrix hk0.  restriction: all elements of diag
!c             must be positive.
!c 
!c     iprint  is an integer array of length two which must be set by the
!c             user.
!c 
!c             iprint(1) specifies the frequency of the output:
!c                iprint(1) < 0 : no output is generated,
!c                iprint(1) = 0 : output only at first and last iteration,
!c                iprint(1) > 0 : output every iprint(1) iterations.
!c 
!c             iprint(2) specifies the type of output generated:
!c                iprint(2) = 0 : iteration count, number of function 
!c                                evaluations, function value, norm of the
!c                                gradient, and steplength,
!c                iprint(2) = 1 : same as iprint(2)=0, plus vector of
!c                                variables and  gradient vector at the
!c                                initial point,
!c                iprint(2) = 2 : same as iprint(2)=1, plus vector of
!c                                variables,
!c                iprint(2) = 3 : same as iprint(2)=2, plus gradient vector.
!c 
!c 
!c     eps     is a positive double precision variable that must be set by
!c             the user, and determines the accuracy with which the solution
!c             is to be found. the subroutine terminates when
!c
!c                         ||g|| < eps max(1,||x||),
!c
!c             where ||.|| denotes the euclidean norm.
!c 
!c     xtol    is a  positive double precision variable that must be set by
!c             the user to an estimate of the machine precision (e.g.
!c             10**(-16) on a sun station 3/60). the line search routine will
!c             terminate if the relative width of the interval of uncertainty
!c             is less than xtol.
!c 
!c     w       is a double precision array of length n(2m+1)+2m used as
!c             workspace for lbfgs. this array must not be altered by the
!c             user.
!c 
!c     iflag   is an integer variable that must be set to 0 on initial entry
!c             to the subroutine. a return with iflag<0 indicates an error,
!c             and iflag=0 indicates that the routine has terminated without
!c             detecting errors. on a return with iflag=1, the user must
!c             evaluate the function f and gradient g. on a return with
!c             iflag=2, the user must provide the diagonal matrix hk0.
!c 
!c             the following negative values of iflag, detecting an error,
!c             are possible:
!c 
!c              iflag=-1  the line search routine mcsrch failed. the
!c                        parameter info provides more detailed information
!c                        (see also the documentation of mcsrch):
!c
!c                       info = 0  improper input parameters.
!c
!c                       info = 2  relative width of the interval of
!c                                 uncertainty is at most xtol.
!c
!c                       info = 3  more than 20 function evaluations were
!c                                 required at the present iteration.
!c
!c                       info = 4  the step is too small.
!c
!c                       info = 5  the step is too large.
!c
!c                       info = 6  rounding errors prevent further progress. 
!c                                 there may not be a step which satisfies
!c                                 the sufficient decrease and curvature
!c                                 conditions. tolerances may be too small.
!c
!c 
!c              iflag=-2  the i-th diagonal element of the diagonal inverse
!c                        hessian approximation, given in diag, is not
!c                        positive.
!c           
!c              iflag=-3  improper input parameters for lbfgs (n or m are
!c                        not positive).
!c 
!c
!c
!c    on the driver:
!c
!c    the program that calls lbfgs must contain the declaration:
!c
!c                       external lb2
!c
!c    lb2 is a block data that defines the default values of several
!c    parameters described in the common section. 
!c
!c 
!c 
!c    common:
!c 
!c     the subroutine contains one common area, which the user may wish to
!c    reference:
!c 
!         common /lb3/mp,lp,gtol,stpmin,stpmax
!c 
!c    mp  is an integer variable with default value 6. it is used as the
!c        unit number for the printing of the monitoring information
!c        controlled by iprint.
!c 
!c    lp  is an integer variable with default value 6. it is used as the
!c        unit number for the printing of error messages. this printing
!c        may be suppressed by setting lp to a non-positive value.
!c 
!c    gtol is a double precision variable with default value 0.9, which
!c        controls the accuracy of the line search routine mcsrch. if the
!c        function and gradient evaluations are inexpensive with respect
!c        to the cost of the iteration (which is sometimes the case when
!c        solving very large problems) it may be advantageous to set gtol
!c        to a small value. a typical small value is 0.1.  restriction:
!c        gtol should be greater than 1.d-04.
!c 
!c    stpmin and stpmax are non-negative double precision variables which
!c        specify lower and uper bounds for the step in the line search.
!c        their default values are 1.d-20 and 1.d+20, respectively. these
!c        values need not be modified unless the exponents are too large
!c        for the machine being used, or unless the problem is extremely
!c        badly scaled (in which case the exponents should be increased).
!c 
!c
!c  machine dependencies
!c
!c        the only variables that are machine-dependent are xtol,
!c        stpmin and stpmax.
!c 
!c
!c  general information
!c 
!c    other routines called directly:  daxpy, ddot, lb1, mcsrch
!c 
!c    input/output  :  no input; diagnostic messages on unit mp and
!c                     error messages on unit lp.
!c 
!c 
!c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!c
      double precision gtol,one,zero,gnorm,ddot,stp1,ftol,stpmin,  &
                       stpmax,stp,ys,yy,sq,yr,beta,xnorm
      integer mp,lp,iter,nfun,point,ispt,iypt,maxfev,info,          &
              bound,npt,cp,i,nfev,inmc,iycn,iscn
      logical finish
!c
      save
      data one,zero/1.0d+0,0.0d+0/
!c
!c     initialize
!c     ----------
!c


      if(iflag.eq.0) go to 10
      go to (172,100) iflag           !每一次从主函数得到新的f和g时,iflag=1,此句的意思是iflag为1时跳到172,iflag为2时跳到100
  10  iter= 0
      if(n.le.0.or.m.le.0) go to 196  !返回 iflag=-3, n or m 有负(不正确) ,直接结束
      if(gtol.le.1.d-04) then
        if(lp.gt.0) write(lp,245)
        gtol=9.d-01
      endif
      nfun= 1
      point= 0
      finish= .false.
      if(diagco) then                !主函数输入的diagco = .fasle.,来自主函数,全程为假
         do 30 i=1,n
 30      if (diag(i).le.zero) go to 195
      else                           !diagco为假, 则对角线diag有值
         do 40 i=1,n
 40      diag(i)= 1.0d0               !第一次迭代的拟牛顿海森矩阵,对角线1
      endif

!c
!c     the work vector w is divided as follows:
!c     ---------------------------------------
!c     the first n locations are used to store the gradient and
!c         other temporary information.
!c     locations (n+1)...(n+m) store the scalars rho.
!c     locations (n+m+1)...(n+2m) store the numbers alpha used
!c         in the formula that computes h*g.
!c     locations (n+2m+1)...(n+2m+nm) store the last m search
!c         steps.
!c     locations (n+2m+nm+1)...(n+2m+2nm) store the last m
!c         gradient differences.
!c
!c     the search steps and gradient differences are stored in a
!c     circular order controlled by the parameter point.
!c
      ispt= n+2*m
      iypt= ispt+n*m     
      do 50 i=1,n
 50   w(ispt+i)= -g(i)*diag(i)       !w=-gh,牛顿方向
      gnorm= dsqrt(ddot(n,g,1,g,1))  !gt*g
      stp1=  iterControl%startdm / gnorm     !采用缩放因子,第一次迭代的初始尝试步长 ( iterControl%startdm * gnorm ) /gnorm 
!c
!c     parameters for line search routine
!c     
      ftol= 1.0d-4
      maxfev= 4

!c     lb1输出iter,nfun,f,gnorm,stp
      if(iprint(1).ge.0) call lb1(iprint,iter,nfun,     &
                          gnorm,n,m,x,f,g,stp,finish,m0,mHat,rms,lambda)
!c
!c    --------------------
!c     main iteration loop
!c    --------------------
!c
 80   iter= iter+1

      rmsPrev= rms                    !用于后面的λ更新判断
      fPrev = f
      gnormPrev=gnorm

      info=0
      bound=iter-1                !iter < = m 时, bound=iter-1 (?)
      if(iter.eq.1) go to 165    !第一次迭代不需要用下面的部分计算-h*g,故直接跳到165后面开始线搜索
      if (iter .gt. m)bound=m
!c
         ys= ddot(n,w(iypt+npt+1),1,w(ispt+npt+1),1)
      if(.not.diagco) then
         yy= ddot(n,w(iypt+npt+1),1,w(iypt+npt+1),1)
         do 90 i=1,n
   90    diag(i)= ys/yy
      else
         iflag=2
         return
      endif
 100  continue
      if(diagco) then
        do 110 i=1,n
 110    if (diag(i).le.zero) go to 195
      endif

!c
!c     compute -h*g using the formula given in: nocedal, j. 1980,
!c     "updating quasi-newton matrices with limited storage",
!c     mathematics of computation, vol.24, no.151, pp. 773-782.
!c     ---------------------------------------------------------
!c
      cp= point
      if (point.eq.0) cp=m
      w(n+cp)= one/ys
      do 112 i=1,n
 112  w(i)= -g(i)                                !注意这里的初值,w=-g
      cp= point                                  !当前计算次数
      do 125 i= 1,bound                          !前面有边界bound计算说明
         cp=cp-1                                 !往回－1,这里cp指示下面的数组
         if (cp.eq. -1)cp=m-1                    !防止cp=1
         sq= ddot(n,w(ispt+cp*n+1),1,w,1)
         inmc=n+m+cp+1
         iycn=iypt+cp*n
         w(inmc)= w(n+cp+1)*sq
         call daxpy(n,-w(inmc),w(iycn+1),1,w,1)
 125  continue
!c
      do 130 i=1,n
 130  w(i)=diag(i)*w(i)                          !前面有关于diag的计算,即H_k(0)
!c
      do 145 i=1,bound
         yr= ddot(n,w(iypt+cp*n+1),1,w,1)
         beta= w(n+cp+1)*yr
         inmc=n+m+cp+1
         beta= w(inmc)-beta
         iscn=ispt+cp*n
         call daxpy(n,beta,w(iscn+1),1,w,1)
         cp=cp+1
         if (cp.eq.m)cp=0
 145  continue


!c
!c     store the new search direction
!c     ------------------------------
!c
       do 160 i=1,n
 160   w(ispt+point*n+i)= w(i)   !w 就是 -hg
!c
!c     obtain the one-dimensional minimizer of the function 
!c     by using the line search routine mcsrch
!c     ----------------------------------------------------
 165  nfev=0
      stp=one                   !初始步长为one=1 
      if (iter.eq.1) stp=stp1   !第一次迭代时的搜索步长
      do 170 i=1,n
 170  w(i)=g(i)                 !w 就是 -hg
 172  continue


 !----线搜索
      
      call mcsrch(iter,n,x,f,fprev,g,w(ispt+point*n+1),stp,ftol,     &
                  xtol,maxfev,info,nfev,diag,f0,f1,g0,g1)


      if (info .eq. -1) then     ! info = -1 返回主函数,计算目标函数及梯度
        iflag=1
        return
      endif
      if (info .ne. 1) go to 190 ! info==1的情况除了满足wolfe条件,还可能是线搜索找不到足够下降的步长
      nfun= nfun + nfev
!!c 放到后面去
!!c     compute the new step and gradient change 
!!c     -----------------------------------------
!!c
!      npt=point*n
!      do 175 i=1,n
!      w(ispt+npt+i)= stp*w(ispt+npt+i)
! 175  w(iypt+npt+i)= g(i)-w(i)
!      point=point+1
!      if (point.eq.m)point=0
!c
!c     termination test
!c     ----------------
!c
      gnorm= dsqrt(ddot(n,g,1,g,1))  
      xnorm= dsqrt(ddot(n,x,1,x,1))
      xnorm= dmax1(1.0d0,xnorm)
      !if (gnorm/xnorm .le. eps) finish=.true.
      !if(rms.le.1.0) finish=.true.
!-----通过外部最小rms和最大迭代次数控制反演
      if( (rms.lt.iterControl%rmsTol).or.(iter.ge.iterControl%maxIter) ) then
          write(*,*) '迭代停止条件满足'
          write(*,*) rms,iterControl%rmsTol
          write(*,*) iter,iterControl%maxIter
          finish=.true.
      endif

      if ( ( abs(rmsPrev - rms) < iterControl%fdiffTol  ) ) then

      write(*,'(a55)') 'Unable to get out of a local minimum. Exiting...'
      finish=.true.
      

      endif

!c
      if(iprint(1).ge.0) call lb1(iprint,iter,nfun,     &
                     gnorm,n,m,x,f,g,stp,finish,m0,mHat,rms,lambda)
      if (finish) then
         iflag=0
         return 
      endif

      write(9151,801) iter , gnormPrev , gnorm , abs(gnormPrev-gnorm), fPrev,f ,abs( fPrev-f ), lambda
!-----------添加正则化因子更新部分---------------------------------------------------
	  ! if alpha is too small, we are not making progress: update lambda      
!        if ( ( abs(rmsPrev - rms) < iterControl%fdiffTol  )             &   
!               .or. ( abs( gnormPrev - gnorm ) < 0.01     )             &               
!               .or. ( abs( fPrev-f ) < 0.01               )             ) then
        if ( ( abs(rmsPrev - rms)   <  iterControl%fdiffTol  )           ) then 
      		! update lambda, penalty functional and gradient
             call update_lambda(lambda,x,f,g,n)
      		! update alpha
      		gnorm = dsqrt(ddot(n,g,1,g,1)) !sqrt(dotProd(grad,grad)) ! 更新梯度(因为λ缩小了)            
!            write(*,'(a34,es12.6)') 'The norm of the last gradient is ',gnorm
!            write(ioLog,'(a34,es12.6)') 'The norm of the last gradient is ',gnorm
            !alpha = min(iterControl%alpha_1,startdm/gnorm)
      		!alpha = min(ONE,startdm)/gnorm
            !stp1 = min(ONE,iterControl%startdm)/gnorm
            ! iter = 0
!            write(555,es12.6)'The value of line search step alpha updated to ',step1
!      		write(*,'(a48,es12.6)') 'The value of line search step alpha updated to ',alpha
!            write(ioLog,'(a48,es12.6)') 'The value of line search step alpha updated to ',alpha
      		! g = - grad
!			call linComb(MinusONE,grad,R_ZERO,grad,g)
			! check that lambda is still at a reasonable value
			if (lambda < iterControl%lambdaTol) then
				write(*,'(a55)') 'Unable to get out of a local minimum. Exiting...'
                 iflag=0
                 return
!                write(ioLog,'(a55)') 'Unable to get out of a local minimum. Exiting...'
!				! multiply by C^{1/2} and add m_0
!                call CmSqrtMult(mHat,m_minus_m0)
!                call linComb(ONE,m_minus_m0,ONE,m0,m)
!                d = dHat
				return
			end if
	  	! restart
!			write(*,'(a55)') 'Restarting LBFGS with the damping parameter updated'
!			call printf('to',lambda,alpha,value,mNorm,rms)
!			write(ioLog,'(a55)') 'Restarting LBFGS with the damping parameter updated'
!			call printf('to',lambda,alpha,value,mNorm,rms,logFile)
!	  	h = g
!	  	nCG = 0
!	  	cycle
	  end if
!----------------------------------------------------------------------------------      
write(9152,801) iter , gnormPrev , gnorm , abs(gnormPrev-gnorm), fPrev,f ,abs( fPrev-f ) , lambda

!c     compute the new step and gradient change  
!c     -----------------------------------------
!c
      npt=point*n
      do 175 i=1,n
      w(ispt+npt+i)= stp*w(ispt+npt+i)
 175  w(iypt+npt+i)= g(i)-w(i)
      point=point+1
      if (point.eq.m)point=0


      go to 80
!c
!c     ------------------------------------------------------------
!c     end of main iteration loop. error exits.
!c     ------------------------------------------------------------
!c
 190  iflag=-1
      if(lp.gt.0) write(lp,200) info
      write(*,200) info
      return
 195  iflag=-2
      if(lp.gt.0) write(lp,235) i
      write(*,235) i
      return
 196  iflag= -3
      if(lp.gt.0) write(lp,240)
      write(*,240)
!c
!c     formats
!c     -------
!c
801   format(1(i4,1x),5x,7(1pd12.5,2x))

 200  format(/' iflag= -1 ',/' line search failed. see'                          &
               ' documentation of routine mcsrch',/' error return'               &
               ' of line search: info= ',i2,/                                    &
               ' possible causes: function or gradient are incorrect',/,         &
               ' or incorrect tolerances')
 235  format(/' iflag= -2',/' the',i5,'-th diagonal element of the',/,         &
            ' inverse hessian approximation is not positive')
 240  format(/' iflag= -3',/' improper input parameters (n or m',             &
            ' are not positive)')
 245  format(/'  gtol is less than or equal to 1.d-04',                     &
            / ' it has been reset to 9.d-01')
      return
      end  subroutine lbfgs
!c
!c     last line of subroutine lbfgs
!c
!c


      subroutine update_lambda(lambda,x,f,g,n)

      real*8  , intent(inout) :: lambda,f
      integer , intent(in)    :: n
      real*8  , intent(in)    :: x(n)
      real*8  , intent(inout) :: g(n)

      double precision ddot
      real*8 :: SS,mNorm
      integer :: i

      mNorm = ddot(n,x,1,x,1)
      SS = f - (lambda*mNorm) / n

      do i= 1, n  !将模型部分的梯度项抽离

         g(i) = g(i) + (MinusTWO*lambda/n) * x(i)

      enddo

      lambda = lambda/iterControl%k     !更新λ

      f = SS + ( lambda * mNorm / n )

      do i= 1, n  !将模型部分的梯度项重新加上

         g(i) = g(i) + (TWO*lambda/n) * x(i)

      enddo
      
      end subroutine update_lambda


      subroutine lb1(iprint,iter,nfun,                     &
                           gnorm,n,m,x,f,g,stp,finish,m0,mHat,rms,lambda)
!c
!c     -------------------------------------------------------------
!c     this routine prints monitoring information. the frequency and
!c     amount of output are controlled by iprint.
!c     -------------------------------------------------------------
!c
      !use lb23 注释放在本模块前
      
      implicit none
!----wkp add for implicit none
      integer :: i

      integer iprint(2),iter,nfun,n,m   !,lp,mp
      double precision x(n),g(n),f,gnorm,stp  !,gtol,stpmin,stpmax
      logical finish
      !common /lb3/mp,lp,gtol,stpmin,stpmax


      double precision rms,lambda
      type(modelParam_t), intent(in)		       :: m0

      type(modelParam_t), intent(inout)		   :: mHat
      type(modelParam_t)    m_minus_m0,mout
      character(3)         :: iterChar
      character(100)       :: mFile
!c
      if (iter.eq.0)then
           write(mp,10)
           write(mp,20) n,m
           write(mp,30)f,gnorm
                          
                          write(555,80)iter,nfun,f,gnorm,stp,rms,lambda
                          call x_to_mHat(x,n,mHat)
                          ! write out the intermediate model solution and responses
                          call CmSqrtMultLBFGS(mHat,m_minus_m0)
   	                      call linComb(ONE,m_minus_m0,ONE,m0,mout)
   	                      write(iterChar,'(i3.3)') iter
   	                      !if (output_level > 1) then
   	                      mFile = trim(iterControl%fname)//'_lbfgs_'//iterChar//'.rho'
                          call write_modelParam(mout,trim(mFile))

                 if (iprint(2).ge.1)then
                     write(mp,40)
                     write(mp,50) (x(i),i=1,n)
                     write(mp,60)
                     write(mp,50) (g(i),i=1,n)
                  endif
           write(mp,10)
           write(mp,70)
      else
          if ((iprint(1).eq.0).and.(iter.ne.1.and..not.finish))return
              if (iprint(1).ne.0)then
                   if(mod(iter-1,iprint(1)).eq.0.or.finish)then
                         if(iprint(2).gt.1.and.iter.gt.1) write(mp,70)
                         write(mp,80)iter,nfun,f,gnorm,stp,rms,lambda

                         write(555,80)iter,nfun,f,gnorm,stp,rms,lambda

                         call x_to_mHat(x,n,mHat)
                          ! write out the intermediate model solution and responses
                          call CmSqrtMultLBFGS(mHat,m_minus_m0)
   	                      call linComb(ONE,m_minus_m0,ONE,m0,mout)
   	                      write(iterChar,'(i3.3)') iter
   	                      !if (output_level > 1) then
   	                        mFile = trim(iterControl%fname)//'_LBFGS_'//iterChar//'.rho'
                            call write_modelParam(mout,trim(mFile))
                          !end if
!   	                      if (output_level > 2) then
!   	                        mHatFile = trim(iterControl%fname)//'_LBFGS_'//iterChar//'.prm'
!                            call write_modelParam(mHat,trim(mHatFile))
!                          end if
!   	                      if (output_level > 2) then
!   	                        dataFile = trim(iterControl%fname)//'_LBFGS_'//iterChar//'.dat'
!                            call write_dataVectorMTX(dHat,trim(dataFile))
!                          end if

                   else
                         return
                   endif
              else
                   if( iprint(2).gt.1.and.finish) write(mp,70)
                   write(mp,80)iter,nfun,f,gnorm,stp
              endif
              if (iprint(2).eq.2.or.iprint(2).eq.3)then
                    if (finish)then
                        write(mp,90)
                    else
                        write(mp,40)
                    endif
                      write(mp,50)(x(i),i=1,n)
                  if (iprint(2).eq.3)then
                      write(mp,60)
                      write(mp,50)(g(i),i=1,n)
                  endif
              endif
            if (finish) write(mp,100)
      endif
!c
 10   format('*************************************************') 
 20   format('  n=',i5,'   number of corrections=',i2,           &
             /,  '       initial values')
 30   format(' f= ',1pd10.3,'   gnorm= ',1pd10.3)
 40   format(' vector x= ')
 50   format(6(2x,1pd10.3))
 60   format(' gradient vector g= ')
 70   format(/'   i   nfn',4x,'func',8x,'gnorm',7x,'steplength',4x,'rms',8x,'lambda'/)
 80   format(2(i4,1x),5x,5(1pd12.5,2x))
 90   format(' final point x= ')
 100  format(/' the minimization terminated without detecting errors.',   &
             /' iflag = 0')


     !call deall_modelParam(mHat) 
     call deall_modelParam(m_minus_m0)
     call deall_modelParam(mout)
!c
      return
      end subroutine lb1
!c     ******
!c
!c
!c   ----------------------------------------------------------
!c     data 我把这个放到前面去了
!c   ----------------------------------------------------------
!c
!      block data lb2
!      integer lp,mp
!      double precision gtol,stpmin,stpmax
!      common /lb3/mp,lp,gtol,stpmin,stpmax
!      data mp,lp,gtol,stpmin,stpmax/6,6,9.0d-01,1.0d-20,1.0d+20/
!      end
!c
!c
!c   ----------------------------------------------------------
!c
      subroutine daxpy(n,da,dx,incx,dy,incy)
!c
!c     constant times a vector plus a vector.
!c     uses unrolled loops for increments equal to one.
!c     jack dongarra, linpack, 3/11/78.
!c
      implicit none

      double precision dx(1),dy(1),da
      integer i,incx,incy,ix,iy,m,mp1,n
!c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
!c
!c        code for unequal increments or equal increments
!c          not equal to 1
!c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!c
!c        code for both increments equal to 1
!c
!c
!c        clean-up loop
!c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end subroutine daxpy
!c
!c
!c   ----------------------------------------------------------
!c   使用的是wsBLAS.f90的 代替 lbfgs里的ddot
!!      double precision function ddot(n,dx,incx,dy,incy)
!!!c
!!!c     forms the dot product of two vectors.
!!!c     uses unrolled loops for increments equal to one.
!!!c     jack dongarra, linpack, 3/11/78.
!!!c
!!      implicit none
!!
!!      double precision dx(1),dy(1),dtemp
!!      integer i,incx,incy,ix,iy,m,mp1,n
!!!c
!!      ddot = 0.0d0
!!      dtemp = 0.0d0
!!      if(n.le.0)return
!!      if(incx.eq.1.and.incy.eq.1)go to 20
!!!c
!!!c        code for unequal increments or equal increments
!!!c          not equal to 1
!!!c
!!      ix = 1
!!      iy = 1
!!      if(incx.lt.0)ix = (-n+1)*incx + 1
!!      if(incy.lt.0)iy = (-n+1)*incy + 1
!!      do 10 i = 1,n
!!        dtemp = dtemp + dx(ix)*dy(iy)
!!        ix = ix + incx
!!        iy = iy + incy
!!   10 continue
!!      ddot = dtemp
!!      return
!!!c
!!!c        code for both increments equal to 1
!!!c
!!!c
!!!c        clean-up loop
!!!c
!!   20 m = mod(n,5)
!!      if( m .eq. 0 ) go to 40
!!      do 30 i = 1,m
!!        dtemp = dtemp + dx(i)*dy(i)
!!   30 continue
!!      if( n .lt. 5 ) go to 60
!!   40 mp1 = m + 1
!!      do 50 i = mp1,n,5
!!        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +         &
!!         dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
!!   50 continue
!!   60 ddot = dtemp
!!      return
!!      end function ddot
!!!c    ------------------------------------------------------------------
!c
!c     **************************
!c     line search routine mcsrch
!c     **************************
!c
      subroutine mcsrch(iter,n,x,f,fprev,g,s,stp,ftol,xtol,maxfev,info,nfev,wa,f0,f1,g0,g1)

      !use lb23  注释放在本模块前

      implicit none

      !integer :: mp,lp

      integer n,maxfev,info,nfev,iter
      double precision f,stp,ftol,xtol,fprev  !,gtol,stpmin,stpmax
      double precision x(n),g(n),s(n),wa(n)

      double precision  f0,f1,g0(n),g1(n)  ! 增加 为了避免无限循环
      !common /lb3/mp,lp,gtol,stpmin,stpmax
      save
!c
!c                     subroutine mcsrch
!c                
!c     a slight modification of the subroutine csrch of more' and thuente.
!c     the changes are to allow reverse communication, and do not affect
!c     the performance of the routine. 
!c
!c     the purpose of mcsrch is to find a step which satisfies
!c     a sufficient decrease condition and a curvature condition.
!c
!c     at each stage the subroutine updates an interval of
!c     uncertainty with endpoints stx and sty. the interval of
!c     uncertainty is initially chosen so that it contains a
!c     minimizer of the modified function
!c
!c          f(x+stp*s) - f(x) - ftol*stp*(gradf(x)'s).
!c
!c     if a step is obtained for which the modified function
!c     has a nonpositive function value and nonnegative derivative,
!c     then the interval of uncertainty is chosen so that it
!c     contains a minimizer of f(x+stp*s).
!c
!c     the algorithm is designed to find a step which satisfies
!c     the sufficient decrease condition  充分下降条件
!c
!c           f(x+stp*s) .le. f(x) + ftol*stp*(gradf(x)'s),
!c
!c     and the curvature condition        强曲率条件
!c
!c           abs(gradf(x+stp*s)'s)) .le. gtol*abs(gradf(x)'s).
!c
!c     if ftol is less than gtol and if, for example, the function
!c     is bounded below, then there is always a step which satisfies
!c     both conditions. if no step can be found which satisfies both
!c     conditions, then the algorithm usually stops when rounding
!c     errors prevent further progress. in this case stp only
!c     satisfies the sufficient decrease condition.
!c
!c     the subroutine statement is
!c
!c        subroutine mcsrch(n,x,f,g,s,stp,ftol,xtol, maxfev,info,nfev,wa)
!c     where
!c
!c       n is a positive integer input variable set to the number
!c         of variables.
!c
!c       x is an array of length n. on input it must contain the
!c         base point for the line search. on output it contains
!c         x + stp*s.
!c
!c       f is a variable. on input it must contain the value of f
!c         at x. on output it contains the value of f at x + stp*s.
!c
!c       g is an array of length n. on input it must contain the
!c         gradient of f at x. on output it contains the gradient
!c         of f at x + stp*s.
!c
!c       s is an input array of length n which specifies the
!c         search direction.
!c
!c       stp is a nonnegative variable. on input stp contains an
!c         initial estimate of a satisfactory step. on output
!c         stp contains the final estimate.
!c
!c       ftol and gtol are nonnegative input variables. (in this reverse
!c         communication implementation gtol is defined in a common
!c         statement.) termination occurs when the sufficient decrease
!c         condition and the directional derivative condition are
!c         satisfied.
!c
!c       xtol is a nonnegative input variable. termination occurs
!c         when the relative width of the interval of uncertainty
!c         is at most xtol.
!c
!c       stpmin and stpmax are nonnegative input variables which
!c         specify lower and upper bounds for the step. (in this reverse
!c         communication implementatin they are defined in a common
!c         statement).
!c
!c       maxfev is a positive integer input variable. termination
!c         occurs when the number of calls to fcn is at least
!c         maxfev by the end of an iteration.
!c
!c       info is an integer output variable set as follows:
!c
!c         info = 0  improper input parameters.
!c
!c         info =-1  a return is made to compute the function and gradient.
!c
!c         info = 1  the sufficient decrease condition and the
!c                   directional derivative condition hold.
!c
!c         info = 2  relative width of the interval of uncertainty
!c                   is at most xtol.
!c
!c         info = 3  number of calls to fcn has reached maxfev.
!c
!c         info = 4  the step is at the lower bound stpmin.
!c
!c         info = 5  the step is at the upper bound stpmax.
!c
!c         info = 6  rounding errors prevent further progress.
!c                   there may not be a step which satisfies the
!c                   sufficient decrease and curvature conditions.
!c                   tolerances may be too small.
!c
!c       nfev is an integer output variable set to the number of
!c         calls to fcn.
!c
!c       wa is a work array of length n.
!c
!c     subprograms called
!c
!c       mcstep
!c
!c       fortran-supplied...abs,max,min
!c
!c     argonne national laboratory. minpack project. june 1983
!c     jorge j. more', david j. thuente
!c
!c     **********
      integer infoc,j
      logical brackt,stage1
      double precision dg,dgm,dginit,dgtest,dgx,dgxm,dgy,dgym,    &
             finit,ftest1,fm,fx,fxm,fy,fym,p5,p66,stx,sty,          &
             stmin,stmax,width,width1,xtrapf,zero
      data p5,p66,xtrapf,zero /0.5d0,0.66d0,4.0d0,0.0d0/
      if(info.eq.-1) go to 45
      infoc = 1


!c
!c     check the input parameters for errors.
!c
      if (n .le. 0 .or. stp .le. zero .or. ftol .lt. zero .or.           &
         gtol .lt. zero .or. xtol .lt. zero .or. stpmin .lt. zero       &
         .or. stpmax .lt. stpmin .or. maxfev .le. 0) return


!c
!c     compute the initial gradient in the search direction
!c     and check that s is a descent direction.
!c     判断输入的方向s(即之前的w=-hg)是否为下降方向
      dginit = zero
      do 10 j = 1, n
         dginit = dginit + g(j)*s(j) !s是之前的-hg,若s为下降方向,梯度g与s相乘应<0
   10    continue
      if (dginit .ge. zero) then
         write(lp,15)
   15    format(/'  the search direction is not a descent direction')
         return
         endif


!c
!c     initialize local variables.
!c
      brackt = .false.
      stage1 = .true.
      nfev = 0
      finit = f               !带进线搜索的目标函数
      dgtest = ftol*dginit    !ftol= 1.0d-4,dgtest=c1*gkt*dk
      width = stpmax - stpmin
      width1 = width/p5
      do 20 j = 1, n          !x变量赋予wa
         wa(j) = x(j)
   20    continue


!c
!c     the variables stx, fx, dgx contain the values of the step,
!c     function, and directional derivative at the best step.
!c     the variables sty, fy, dgy contain the value of the step,
!c     function, and derivative at the other endpoint of
!c     the interval of uncertainty.
!c     the variables stp, f, dg contain the values of the step,
!c     function, and derivative at the current step.
!c
      stx = zero
      fx  = finit
      dgx = dginit
      sty = zero
      fy  = finit
      dgy = dginit
!c
!c     start of iteration.
!c
   30 continue


!c
!c        set the minimum and maximum steps to correspond
!c        to the present interval of uncertainty.
!c
         if (brackt) then      !第一次线搜索时,brackt为假
            stmin = min(stx,sty)
            stmax = max(stx,sty)
         else
            stmin = stx         !第一次线搜索,stmin=stx=0,stmax有所改变
            stmax = stp + xtrapf*(stp - stx)
            end if
!c
!c        force the step to be within the bounds stpmax and stpmin.
!c        步长只能在 步长区间内  stpmin<=stp>=stpmax
         stp = max(stp,stpmin)
         stp = min(stp,stpmax)
!c
!c        if an unusual termination is to occur then let
!c        stp be the lowest point obtained so far.
!c        出现下面任何情况,都使步长为步长区间的最小
         if ((brackt .and. (stp .le. stmin .or. stp .ge. stmax))            &
            .or. nfev .ge. maxfev-1 .or. infoc .eq. 0                       &
            .or. (brackt .and. stmax-stmin .le. xtol*stmax)) stp = stx



!c
!c        evaluate the function and gradient at stp
!c        and compute the directional derivative.
!c        we return to main program to obtain f and g.
!c
         do 40 j = 1, n
            x(j) = wa(j) + stp*s(j)  !修正x,stp为步长,s为修正方向
   40       continue
         info=-1                     !-1会让lbfgs直接跳到线搜索
         return                     !并且得到新的f,g后,前面会指示直接跳到45行来
!--------返回主函数

!c
   45    info=0                      !info重置为0  
         nfev = nfev + 1
         dg = zero
         do 50 j = 1, n
            dg = dg + g(j)*s(j)      !s是之前的-hg,若s为下降方向,梯度g与s相乘应<0
   50       continue
         ftest1 = finit + stp*dgtest !ftest1=f(前一个目标函数)+c1*αk*gkt*dk


!c
!c        test for convergence.
!c
         if ((brackt .and. (stp .le. stmin .or. stp .ge. stmax))           &
             .or. infoc .eq. 0) info = 6
         if (stp .eq. stpmax .and.                                         &
             f .le. ftest1 .and. dg .le. dgtest) info = 5                  
         if (stp .eq. stpmin .and.                                         &
            (f .gt. ftest1 .or. dg .ge. dgtest)) info = 4                
         if (nfev .ge. maxfev) info = 3
         if (brackt .and. stmax-stmin .le. xtol*stmax) info = 2
!        强wolfe条件同时满足:充分下降条件及强曲率条件
!         if (iter==1 .and. nfev <= 2) then  去掉这一步 因为我发现效率其实是在降低
!             info=0       
!         else
             if (f .le. ftest1 .and. abs(dg) .le. gtol*(-dginit)) then 
                 info = 1
             endif
!        endif
!80   format(2(i4,1x),5x,5(1pd12.5,2x)) 
write(691,890) iter,nfev,fprev,f,stp

         if(nfev==1) then 
            f0=f
            g0=g
         endif

         if(nfev==2) then
            f1=f
            g1=g
         endif

         if(nfev==2 .and. info /= 1) then 

            if(f0<f1) then

               f=f0
               g=g0
               info=1

            else

               f=f1
               g=g1
               info=1 

            endif

         endif



!!--------注释 for EPS 线搜索最大次数改为4
!             if ( abs( fprev - f ) < iterControl%fdiffTol ) then   !这么做是防止线搜索过多
!                 info = 1
!             endif



!c
!c        check for termination.
!c
         if (info .ne. 0) return
!c
!c        in the first stage we seek a step for which the modified
!c        function has a nonpositive value and nonnegative derivative.
!c        ftol = 1.0d-4,gtol = 9.0d-01 = 0.9 
!c        下面不是强wolfe条件,看清楚,满足这个条件,stage1 = .false.
         if (stage1 .and. f .le. ftest1 .and.                               &
             dg .ge. min(ftol,gtol)*dginit) stage1 = .false.  
!c
!c        a modified function is used to predict the step only if
!c        we have not obtained a step for which the modified
!c        function has a nonpositive function value and nonnegative
!c        derivative, and if a lower function value has been
!c        obtained but the decrease is not sufficient.
!c        第一个if执行条件,stage1为真,表示不满足“非强wolfe条件”
         if (stage1 .and. f .le. fx .and. f .gt. ftest1) then
!c
!c           define the modified function and derivative values.
!c
            fm = f - stp*dgtest
            fxm = fx - stx*dgtest
            fym = fy - sty*dgtest
            dgm = dg - dgtest
            dgxm = dgx - dgtest
            dgym = dgy - dgtest
!c
!c           call cstep to update the interval of uncertainty
!c           and to compute the new step.
!c           更新步长区间,计算新步长
            call mcstep(iter,stx,fxm,dgxm,sty,fym,dgym,stp,fm,dgm,   &
                       brackt,stmin,stmax,infoc)
!c
!c           reset the function and gradient values for f.
!c
            fx = fxm + stx*dgtest
            fy = fym + sty*dgtest
            dgx = dgxm + dgtest
            dgy = dgym + dgtest
         else
!c
!c           call mcstep to update the interval of uncertainty
!c           and to compute the new step.
!c
            call mcstep(iter,stx,fx,dgx,sty,fy,dgy,stp,f,dg,      &
                       brackt,stmin,stmax,infoc)
            end if
!c
!c        force a sufficient decrease in the size of the
!c        interval of uncertainty.
!c
         if (brackt) then
            if (abs(sty-stx) .ge. p66*width1)                 &
               stp = stx + p5*(sty - stx)
            width1 = width
            width = abs(sty-stx)
            end if
!c
!c        end of iteration.
!c
         go to 30
!c
!c     last line of subroutine mcsrch.
!c

890   format(2(i4,1x),5x,3(1pd12.5,2x)) 

      end subroutine mcsrch

      subroutine mcstep(iter,stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,     &
                       stpmin,stpmax,info)

      implicit none

      integer info,iter
      double precision stx,fx,dx,sty,fy,dy,stp,fp,dp,stpmin,stpmax
      logical brackt,bound
!c
!c     subroutine mcstep
!c
!c     the purpose of mcstep is to compute a safeguarded step for
!c     a linesearch and to update an interval of uncertainty for
!c     a minimizer of the function.
!c
!c     the parameter stx contains the step with the least function
!c     value. the parameter stp contains the current step. it is
!c     assumed that the derivative at stx is negative in the
!c     direction of the step. if brackt is set true then a
!c     minimizer has been bracketed in an interval of uncertainty
!c     with endpoints stx and sty.
!c
!c     the subroutine statement is
!c
!c       subroutine mcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
!c                        stpmin,stpmax,info)
!c
!c     where
!c
!c       stx, fx, and dx are variables which specify the step,
!c         the function, and the derivative at the best step obtained
!c         so far. the derivative must be negative in the direction
!c         of the step, that is, dx and stp-stx must have opposite
!c         signs. on output these parameters are updated appropriately.
!c
!c       sty, fy, and dy are variables which specify the step,
!c         the function, and the derivative at the other endpoint of
!c         the interval of uncertainty. on output these parameters are
!c         updated appropriately.
!c
!c       stp, fp, and dp are variables which specify the step,
!c         the function, and the derivative at the current step.
!c         if brackt is set true then on input stp must be
!c         between stx and sty. on output stp is set to the new step.
!c
!c       brackt is a logical variable which specifies if a minimizer
!c         has been bracketed. if the minimizer has not been bracketed
!c         then on input brackt must be set false. if the minimizer
!c         is bracketed then on output brackt is set true.
!c
!c       stpmin and stpmax are input variables which specify lower
!c         and upper bounds for the step.
!c
!c       info is an integer output variable set as follows:
!c         if info = 1,2,3,4,5, then the step has been computed
!c         according to one of the five cases below. otherwise
!c         info = 0, and this indicates improper input parameters.
!c
!c     subprograms called
!c
!c       fortran-supplied ... abs,max,min,sqrt
!c
!c     argonne national laboratory. minpack project. june 1983
!c     jorge j. more', david j. thuente
!c
      double precision gamma,p,q,r,s,sgnd,stpc,stpf,stpq,theta
      info = 0
!c
!c     check the input parameters for errors.
!c
      if ((brackt .and. (stp .le. min(stx,sty) .or.              &
          stp .ge. max(stx,sty))) .or.                           &
          dx*(stp-stx) .ge. 0.0 .or. stpmax .lt. stpmin) return
!c
!c     determine if the derivatives have opposite sign.
!c
      sgnd = dp*(dx/abs(dx))
!c
!c     first case. a higher function value.
!c     the minimum is bracketed. if the cubic step is closer
!c     to stx than the quadratic step, the cubic step is taken,
!c     else the average of the cubic and quadratic steps is taken.
!c
      if (fp .gt. fx) then
         info = 1
         bound = .true.
         theta = 3*(fx - fp)/(stp - stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))
         gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
         if (stp .lt. stx) gamma = -gamma
         p = (gamma - dx) + theta
         q = ((gamma - dx) + gamma) + dp
         r = p/q
         stpc = stx + r*(stp - stx)
         stpq = stx + ((dx/((fx-fp)/(stp-stx)+dx))/2)*(stp - stx)
         if (abs(stpc-stx) .lt. abs(stpq-stx)) then
            stpf = stpc
         else
           stpf = stpc + (stpq - stpc)/2
           end if
         brackt = .true.
!c
!c     second case. a lower function value and derivatives of
!c     opposite sign. the minimum is bracketed. if the cubic
!c     step is closer to stx than the quadratic (secant) step,
!c     the cubic step is taken, else the quadratic step is taken.
!c
      else if (sgnd .lt. 0.0) then
         info = 2
         bound = .false.
         theta = 3*(fx - fp)/(stp - stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))
         gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
         if (stp .gt. stx) gamma = -gamma
         p = (gamma - dp) + theta
         q = ((gamma - dp) + gamma) + dx
         r = p/q
         stpc = stp + r*(stx - stp)
         stpq = stp + (dp/(dp-dx))*(stx - stp)
         if (abs(stpc-stp) .gt. abs(stpq-stp)) then
            stpf = stpc
         else
            stpf = stpq
            end if
         brackt = .true.
!c
!c     third case. a lower function value, derivatives of the
!c     same sign, and the magnitude of the derivative decreases.
!c     the cubic step is only used if the cubic tends to infinity
!c     in the direction of the step or if the minimum of the cubic
!c     is beyond stp. otherwise the cubic step is defined to be
!c     either stpmin or stpmax. the quadratic (secant) step is also
!c     computed and if the minimum is bracketed then the the step
!c     closest to stx is taken, else the step farthest away is taken.
!c
      else if (abs(dp) .lt. abs(dx)) then
         info = 3
         bound = .true.
         theta = 3*(fx - fp)/(stp - stx) + dx + dp
         s = max(abs(theta),abs(dx),abs(dp))
!c
!c        the case gamma = 0 only arises if the cubic does not tend
!c        to infinity in the direction of the step.
!c
         gamma = s*sqrt(max(0.0d0,(theta/s)**2 - (dx/s)*(dp/s)))
         if (stp .gt. stx) gamma = -gamma
         p = (gamma - dp) + theta
         q = (gamma + (dx - dp)) + gamma
         r = p/q
         if (r .lt. 0.0 .and. gamma .ne. 0.0) then
            stpc = stp + r*(stx - stp)
         else if (stp .gt. stx) then
            stpc = stpmax
         else
            stpc = stpmin
            end if
         stpq = stp + (dp/(dp-dx))*(stx - stp)
         if (brackt) then
            if (abs(stp-stpc) .lt. abs(stp-stpq)) then
               stpf = stpc
            else
               stpf = stpq
               end if
         else
            if (abs(stp-stpc) .gt. abs(stp-stpq)) then
               stpf = stpc
            else
               stpf = stpq
               end if
            end if
!c
!c     fourth case. a lower function value, derivatives of the
!c     same sign, and the magnitude of the derivative does
!c     not decrease. if the minimum is not bracketed, the step
!c     is either stpmin or stpmax, else the cubic step is taken.
!c
      else
         info = 4
         bound = .false.
         if (brackt) then
            theta = 3*(fp - fy)/(sty - stp) + dy + dp
            s = max(abs(theta),abs(dy),abs(dp))
            gamma = s*sqrt((theta/s)**2 - (dy/s)*(dp/s))
            if (stp .gt. sty) gamma = -gamma
            p = (gamma - dp) + theta
            q = ((gamma - dp) + gamma) + dy
            r = p/q
            stpc = stp + r*(sty - stp)
            stpf = stpc
         else if (stp .gt. stx) then
            stpf = stpmax
         else
            stpf = stpmin
            end if
         end if
!c
!c     update the interval of uncertainty. this update does not
!c     depend on the new step or the case analysis above.
!c
      if (fp .gt. fx) then
         sty = stp
         fy = fp
         dy = dp
      else
         if (sgnd .lt. 0.0) then
            sty = stx
            fy = fx
            dy = dx
            end if
         stx = stp
         fx = fp
         dx = dp
         end if
!c
!c     compute the new step and safeguard it.
!c
!     if(iter==1)then
!        stpmax = 1.0
!     endif

      stpf = min(stpmax,stpf)
      stpf = max(stpmin,stpf)
      stp = stpf
      if (brackt .and. bound) then
         if (sty .gt. stx) then
            stp = min(stx+0.66*(sty-stx),stp)
         else
            stp = max(stx+0.66*(sty-stx),stp)
            end if
         end if
      return
!c
!c     last line of subroutine mcstep.
!c
      end subroutine mcstep
!--------------------------------------------------------------------------------------------finish!


end module mLBFGS
