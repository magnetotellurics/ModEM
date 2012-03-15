! *****************************************************************************
module dataMisfit
  ! Module containing the definitions and subroutines related to the penalty
  ! functional only. Only use it if we want to compute the data misfit.

  use responses
  use functionals
  use model_operators
  use UserData
  use SolnSpace
  use DataSpace
  use dataTypes
  use DataFunc
  implicit none

  ! ***************************************************************************
  ! * type sensitivity_t contains the full information about the data sensitivities
  ! * with respect to the original model parameters and to
  ! * each cell resistivity for all frequencies
  ! * used in the Jacobian calculations only
  type :: sensitivity_t

    complex(8), pointer, dimension(:,:,:,:) :: da  !(nfreq,nfunc,nobs,nvar)
    real(8), pointer, dimension(:,:,:,:)    :: da_real  !(nfreq,nfunc,nobs,nvar)
    real(8), pointer, dimension(:,:,:,:)    :: da_imag  !(nfreq,nfunc,nobs,nvar)
    type (rscalar)                          :: drho_real
    type (rscalar)                          :: drho_imag
    !real(8), pointer, dimension(:,:,:)     :: drho_real  !(nx,ny,nz) - single freq, all func
    !real(8), pointer, dimension(:,:,:)     :: drho_imag  !(nx,ny,nz) - single freq, all func

  end type sensitivity_t


  ! ***************************************************************************
  ! * type misfit_t contains the full information about the misfit and its'
  ! * partial derivatives with respect to the original model parameters and to
  ! * each cell resistivity for all frequencies
  ! * used in the Jacobian and derivative calculations
  type :: misfit_t

    character(80)                           :: name
    real(8)                   :: damping
    real(8), pointer, dimension(:,:)        :: value  !(nfreq,nfunc)
    integer, pointer, dimension(:,:)        :: ndat   !(nfreq,nfunc)
    real(8), pointer, dimension(:)          :: weight !nfunc
    real(8), pointer, dimension(:,:,:)      :: dRda   !(nfreq,nfunc,ncoeff)
    !real(8), pointer, dimension(:,:,:,:)   :: dRda   !(nfreq,nfunc,nlayer,nparam)

  end type misfit_t

  type (misfit_t), save                            :: misfit
  type (sensitivity_t), save                       :: sens

  ! ***************************************************************************
  ! * storing data misfit for each freq and func, dim(nfreq,nfunc)
  ! * total number of observations per data type & misfit weights, dim(nfunc)
  !real(8), allocatable, dimension(:,:), save       :: misfit
  integer, allocatable, dimension(:,:), save        :: ndat
  !real(8), allocatable, dimension(:), save     :: weight
  real(8), allocatable, dimension(:), save      :: misfitValue
  real(8), allocatable, dimension(:,:), save        :: dmisfitValue

  !type (misfitInfo_t),dimension(:,:),allocatable,save   :: misfitInfo ! (nfreq,nfunc)

  ! ***************************************************************************
  ! * type sensMatrix_t stores the full sensitivity matrix or the
  ! * partial derivatives with respect to the original model parameters
  ! * for all frequencies and all data functionals

  type :: sensMatrix_tmp_t

    type(modelParam_t), pointer, dimension(:,:)  :: dm    !(nfreq,nfunc)
    logical                                      :: allocated=.false.

  end type sensMatrix_tmp_t


Contains

  ! ***************************************************************************
  ! * calcResponses is the subroutine to compute the values of data functionals

  subroutine calcResponses(H,dat,psi,m0)

	!uses: TFList,obsList

	type (dataVector_t), intent(in)			:: dat
	type (dataVector_t), intent(inout)		:: psi
	type (cvector), intent(in)				:: H
	type (modelParam_t), intent(in)         :: m0
	complex(8)                              :: Resp
	integer									:: i,j,k,itype,iobs,nSite

    psi = dat
	i = dat%tx

	do j=1,dat%ndt
	  itype = dat%data(j)%dataType
	  nSite = dat%data(j)%nSite
!	  do k=1,nSite
!	      iobs = psi%data(j)%rx(k)
!	      call dataResp(H,m0,itype,iobs,psi%data(j)%value(:,k))
!          psi%data(j)%error = 0.0d0
!          psi%data(j)%errorBar = .FALSE.
!	  end do
	  select case ( TFList%info(itype)%name )
	  case ('C')
		do k=1,nSite
		  iobs = psi%data(j)%rx(k)
		  Resp = dataResp_rx(C_ratio,obsList%info(iobs),H)
		  psi%data(j)%value(1,k) = dreal(Resp)
          psi%data(j)%value(2,k) = dimag(Resp)
		  psi%data(j)%error(:,k) = 0.0d0
		end do
        psi%data(j)%errorBar = .FALSE.
	  case ('D')
		do k=1,nSite
          iobs = psi%data(j)%rx(k)
          Resp = dataResp_rx(D_ratio,obsList%info(iobs),H)
          psi%data(j)%value(1,k) = dreal(Resp)
          psi%data(j)%value(2,k) = dimag(Resp)
          psi%data(j)%error(:,k) = 0.0d0
		end do
        psi%data(j)%errorBar = .FALSE.
	  case default
		write(0,*) 'Warning: no transfer functions specified to compute'
		stop
	  end select
	  do k=1,nSite
        iobs = psi%data(j)%rx(k)
		if (.not.obsList%info(iobs)%defined) then
		  !psi%data(j)%exists(k) = .FALSE.
		end if
	  end do
	end do

  end subroutine calcResponses  ! calcResponses


  ! ***************************************************************************
  ! * calcResiduals is the subroutine to compute the data residuals from the
  ! * data and the corresponding responses obtained by first computing the fields.

  subroutine calcResiduals(dat,psi,res,weighted)

	type (dataVector_t), intent(in)				:: dat,psi
	type (dataVector_t), intent(inout)			:: res
	logical, intent(in), optional				:: weighted
	integer										:: i,j,k,itype,iobs,nSite
	real(8)									    :: value(2)
	real(8)										:: error

    res = dat
    i = dat%tx

	! Cycle over functional types and the observatories for this frequency
    do j=1,dat%ndt
      itype = dat%data(j)%dataType
      res%data(j)%value = dat%data(j)%value - psi%data(j)%value
	end do

	! If instead weighted responses are required, divide by error squared
	if (present(weighted)) then
	  if (weighted) then
		do j=1,dat%ndt
		  nSite = res%data(j)%nSite
		  do k=1,nSite
			value = res%data(j)%value(:,k)
			error = res%data(j)%error(1,k)
			res%data(j)%value(:,k) = value/(error**2)
			res%data(j)%error(:,k) = ONE
		  end do
		end do
	  end if
	end if

  end subroutine calcResiduals  ! calcResiduals


  ! ***************************************************************************
  ! * calcMisfit performs the misfit computations per frequency

  subroutine calcMisfit(res,misfit)

	type (misfit_t), intent(inout)			:: misfit
	type (dataVector_t), intent(in)			:: res
	real(8)										:: val,err
	integer										:: iTx,iDt,j,k,icomp,nSite

	iTx = res%tx

	!call OutputResiduals(freq,res,TFList,obsList)
	do j=1,res%ndt
       iDt = res%data(j)%dataType

	   do k=1,res%data(j)%nSite
		!if (.not.res%data(j)%exists(k)) then
		!  cycle
		!end if
		do icomp = 1,res%data(j)%nComp

		val = res%data(j)%value(icomp,k)
		err = res%data(j)%error(icomp,k)

		! Calculate the misfit
		select case (trim(misfit%name))
		case ('Mean Squared')
		  misfit%value(iTx,iDt) = &
			 misfit%value(iTx,iDt) + (val/err)**2
		case default
		  write(0,*) 'Warning: (compute_misfit) unknown penalty functional ',trim(misfit%name)
		  return
		end select
		end do ! components
	  end do ! sites

	  !write(0,'(a40,2i3,i6,g15.7)') 'iTx,iDt,ndat,misfit = ',iTx,iDt,&
	  !				misfit%ndat(iTx,iDt),misfit%value(iTx,iDt)/(2*misfit%ndat(iTx,iDt))
	end do

  end subroutine calcMisfit	! calcMisfit


  ! ***************************************************************************
  ! * misfitSumUp

  subroutine misfitSumUp(p_input,misfit,misfitValue,dmisfitValue)

	! this is temporary: output to files will be moved to output.f90
	use iotypes
	use iospec
	! uses: TFList
	type (misfit_t), intent(inout)			  :: misfit
    type (modelParam_t), intent(in)           :: p_input
	type (modelParam_t)							  :: param,dparam,weighted_norm
  	real(8), dimension(:), intent(out)			  :: misfitValue  !nfunc
	real(8), dimension(:,:), intent(out),optional :: dmisfitValue !nfunc,ncoeff
	type (modelCoeff_t), dimension(p_input%nc)	  :: norm,grad,point
	character(80)                    :: comment
	real(8)										  :: v1,v2
	integer										  :: j,l

	misfitValue = 0.0d0
	write(0,*) 'Penalty functional used: ',trim(misfit%name)
	do j=1,TFList%n
	  misfitValue(j) = sum(misfit%value(:,j))/(2*sum(misfit%ndat(:,j)))
	  write(0,*) 'Total data misfit for ',&
				  trim(TFList%info(j)%name),' response = ',misfitValue(j)
	end do

	v1=dot_product(misfit%weight,misfitValue)
	v2=dotProd_modelParam(p_input,p_input)
	v2=v2*misfit%damping

	write(0,*) 'Weighted model norm = ',v2

	write(0,*) 'Total data misfit   = ',v1+v2

	do j=1,TFList%n
	  misfitValue(j) = misfit%weight(j)*misfitValue(j)+v2
	end do

	!---------------------------------------------------------------------------
	! Important output to file; this will be a separate subroutine in output.f90
	call initFileWrite(cUserDef%fn_misfit, ioWRITE)

	write(ioWRITE,'(a10)',advance='no') 'MISFIT = '
	write (ioWRITE,'(g17.9)') v1+v2

	close(ioWRITE)
	!---------------------------------------------------------------------------

	! Output derivative summary
    param = multBy_CmSqrt(p_input)
	dparam = param

	if(.not.present(dmisfitValue)) then
	  return
	end if

	do j=1,TFList%n
		write(comment,*) 'Total derivative for ',trim(TFList%info(j)%name),' response = '
!	  if (output_level>2) then
!	  write(0,*) 'Total derivative for ',&
!				  trim(TFList%info(j)%name),' response = '
!	  end if
	  do l=1,param%nc
		dmisfitValue(j,l) = sum(misfit%dRda(:,j,l))/(2*sum(misfit%ndat(:,j)))
	  end do
	  call fillParamValues_modelParam(dparam,dmisfitValue(j,:))
		call print_modelParam(dparam,output_level,trim(comment))
		call getCoeffArray_modelParam(dparam,grad)
!	  do l=1,param%nc
!		grad(l)=getCoeff_modelParam(dparam,l)
!		if (output_level>1) then
!		  write (0,'(2i8,g17.9)') grad(l)%L%num,grad(l)%F%num,grad(l)%value
!	        end if
!	  end do

	end do

	!if (output_level>1) then
	!  write(0,*) 'Weighted model norm = '
	!end if

	!do l=1,param%nc
	!  norm(l)=getCoeff_modelParam(p_input,l)
	!  if (output_level>1) then
	!    write (0,'(2i8,g17.9)') norm(l)%L%num,norm(l)%F%num,2.0d0*misfit%damping*norm(l)%value
	!  end if
	!end do
	weighted_norm = p_input
	call getCoeffArray_modelParam(p_input,norm)
	call scMult_modelParam(TWO*misfit%damping,p_input,weighted_norm)
	call print_modelParam(weighted_norm,output_level,"Weighted model norm derivative = ")

	! Compute total gradient including model norm term
	do l=1,param%nc
	  grad(l)%value=dot_product(misfit%weight,dmisfitValue(:,l))
	end do
	!dparam = dparam + weighted_norm

	norm%value = 2.0d0 * misfit%damping * norm%value
	grad%value = grad%value + norm%value
	call fillParam_modelParam(dparam,grad)


	call print_modelParam(dparam,output_level,"Total derivative = ")
!	do l=1,param%nc
!	  write (0,'(2i8,g17.9)') grad(l)%L%num,grad(l)%F%num,grad(l)%value
!	end do

	do j=1,TFList%n
	  dmisfitValue(j,:) = misfit%weight(j)*dmisfitValue(j,:)+norm%value
	end do

	call getCoeffArray_modelParam(p_input,point)

	!---------------------------------------------------------------------------
	! Important output to file; this will be a separate subroutine in output.f90
	call initFileWrite(cUserDef%fn_gradient, ioWRITE)

	write(ioWRITE,'(a15)',advance='no') 'DERIVATIVE = '
	do l=1,param%nc
	  if(grad(l)%frozen) then
		cycle
	  end if
	  write (ioWRITE,'(g17.9)',advance='no') grad(l)%value
	end do

	close(ioWRITE)
	!---------------------------------------------------------------------------
	! Important output to file; this will be a separate subroutine in output.f90
	call initFileWrite(cUserDef%fn_point, ioWRITE)

	write(ioWRITE,'(i1)') 1
	do l=1,param%nc
	  if(point(l)%frozen) then
		cycle
	  end if
	  write (ioWRITE,'(g17.9)',advance='no') point(l)%value
	end do

	close(ioWRITE)
	!---------------------------------------------------------------------------
	! Free memory
	call deall_modelParam(dparam)
	call deall_modelParam(weighted_norm)

  end subroutine misfitSumUp  ! misfitSumUp


end module dataMisfit
