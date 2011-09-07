! *****************************************************************************
module model_operators
	! This module contains the elementary operators for computations involving
	! model parametrization structures. 'Y' stands for spherical harmonics;
	! however, many of the subroutines here would be identical for a different
	! kind of layered parametrization.
	! For all parametrization <-> vector interactions, vectors are constructed
	! from those coefficients in the parametrization, whose values 'exist'.

  use math_constants
  use file_units
  use utilities
  use modeldef
#ifdef MPI
  use MPI_declaration
#endif
  implicit none

  ! * BOP
  ! equal to
  interface assignment (=)
     MODULE PROCEDURE copy_modelParam
	 MODULE PROCEDURE fillParam_modelParam
  end interface

  interface zero
     MODULE PROCEDURE zero_modelParam
  end interface

  interface iszero
     MODULE PROCEDURE iszero_modelParam
  end interface

  interface deall
     MODULE PROCEDURE deall_modelParam
  end interface

  interface countModelParam
     MODULE PROCEDURE count_modelParam
  end interface

!  INTERFACE OPERATOR (+)
!     MODULE PROCEDURE add_modelParam
!  END INTERFACE
!
!  INTERFACE OPERATOR (-)
!     MODULE PROCEDURE subtract_modelParam
!  END INTERFACE
!
!  ! multiplication
!  INTERFACE OPERATOR (*)
!     MODULE PROCEDURE mult_modelParam
!     MODULE PROCEDURE scMult_modelParam
!  END INTERFACE
!
!  ! division
!  INTERFACE OPERATOR (/)
!     MODULE PROCEDURE scDiv_modelParam
!  END INTERFACE

  INTERFACE dotProd
     MODULE PROCEDURE dotProd_modelParam
     MODULE PROCEDURE dotProdVec_modelParam
  END INTERFACE

  INTERFACE linComb
     MODULE PROCEDURE linComb_modelParam
  END INTERFACE

  INTERFACE scMult
     MODULE PROCEDURE scMult_modelParam
  END INTERFACE

  INTERFACE scMultAdd
     MODULE PROCEDURE scMultAdd_modelParam
  END INTERFACE

  INTERFACE compare
     MODULE PROCEDURE compare_modelParam
     MODULE PROCEDURE compareLayers_modelParam
  END INTERFACE

  INTERFACE getDegree
     MODULE PROCEDURE getDegree_modelParam
     MODULE PROCEDURE getLayerDegree_modelParam
  END INTERFACE

  INTERFACE getRadial
     MODULE PROCEDURE getRadial_modelParam
  END INTERFACE
  ! * EOP

  public			:: create_modelParam, deall_modelParam, setup_modelParam, copy_modelParam
  public			:: fillParam_modelParam, fillParamValues_modelParam, random_modelParam
  public			:: compare_modelParam, compareLayers_modelParam, adjustLayers_modelParam
  public			:: zero_modelParam, getParamValues_modelParam
  public			:: setLayer_modelParam, setCoeffValue_modelParam, setCrust_modelParam
  public            :: getCoeffValue_modelParam, getCoeff_modelParam, getCoeffArray_modelParam
  public			:: getDegree_modelParam, getLayerDegree_modelParam
  public			:: getRadial_modelParam
  public			:: add_modelParam, subtract_modelParam, linComb_modelParam
  public			:: mult_modelParam, dotProd_modelParam, dotProdVec_modelParam
  public			:: scMult_modelParam, scMultAdd_modelParam, scDiv_modelParam
  public			:: print_modelParam, write_modelParam, read_modelParam
  public			:: smoothV_modelParam, smoothH_modelParam
  public			:: multBy_CmSqrt, multBy_Cm
  public            :: count_modelParam

Contains

#ifdef MPI
#include "model_operators_MPI.inc"
#endif


  ! **********************************************************************
  ! * BOP
  subroutine deall_modelParam(P)

    implicit none
    type (modelParam_t)               :: P
    ! * EOP

    integer                           :: status

    if(P%allocated) then
	  deallocate(P%F,STAT=status)
	  deallocate(P%L,STAT=status)
	  deallocate(P%c,STAT=status)
	  deallocate(P%crust%cond, STAT=status)
	  P%crust%allocated = .false.
    end if

    if(P%rho%allocated) then
      call deall_rscalar(P%rho)
    end if

	P%nF = 0
	P%nL = 0
	P%nc = 0
	P%smoothed  = .FALSE.
	P%zeroValued = .FALSE.
	P%allocated = .FALSE.

  end subroutine deall_modelParam


  ! **********************************************************************
  ! * Test whether two parametrizations have the same basic definitions
  ! * Returns 1 if parametrizations are compatible, zero otherwise
  ! * BOP
  function compare_modelParam(P2,P1) result (status)

    implicit none
    type (modelParam_t), intent(in)                  :: P1
    type (modelParam_t), intent(in)                  :: P2
	logical                                       :: status
    ! * EOP

	integer                                        :: j

	status = .FALSE.

	if(.not.(P1%allocated .and. P2%allocated)) then
		call warning('(compare_modelParam) parametrization not allocated yet')
		return
	end if

    if ((trim(P1%type) .ne. 'harmonic') .or. (trim(P2%type) .ne. 'harmonic')) then
        call warning('(compare_modelParam) no comparison will be performed for this parametrization type')
        return
    end if

	! compare functional construction
	if (P1%nL /= P2%nL) then
		call warning('(compare_modelParam) the two models have a different number of layers')
		return
	end if

	! compare layer construction
	if (P1%nF /= P2%nF) then
		call warning('(compare_modelParam) the two models have a different number of functionals')
		return
	end if

	! compare parameter construction
	if (P1%nc /= P2%nc) then
		call warning('(compare_modelParam) the two models have a different number of coefficients')
		return
	end if

	do j = 1,P1%nL
		status = compareLayers_modelParam(P1%L(j),P2%L(j))
		if (status .eqv. .FALSE.) then
			write(0,*) '(compare_modelParam) incompatible layer',j
			return
		end if
	end do

  end function compare_modelParam


  ! **********************************************************************
  ! * Test whether two layer structures have the same basic definitions
  ! * Returns 1 if layers are identical, zero otherwise
  ! * BOP
  function compareLayers_modelParam(L1,L2) result (status)

    implicit none
    type (modelLayer_t), intent(in)      :: L1,L2
    logical                            :: status
    ! * EOP

	status = .FALSE.

	if (L1%if_log .neqv. L2%if_log) then
		call warning('(compareLayers) layers have different structure')
		return
	end if

	if (clean(L1%depth)  /= clean(L2%depth)) then
		call warning('(compareLayers) layers have different depths')
		return
	end if

	if (clean(L1%width)  /= clean(L2%width)) then
		call warning('(compareLayers) layers have different widths')
		return
	end if

	if (clean(L1%lbound) /= clean(L2%lbound)) then
		call warning('(compareLayers) layers have different lower bounds')
		return
	end if

	if (clean(L1%ubound) /= clean(L2%ubound)) then
		call warning('(compareLayers) layers have different upper bounds')
		return
	end if

	if (clean(L1%alpha)  /= clean(L2%alpha)) then
		call warning('(compareLayers) layers have different horizontal regularization')
		return
	end if

	if (clean(L1%beta)   /= clean(L2%beta)) then
		call warning('(compareLayers) layers have different vertical regularization')
		return
	end if

	status = .TRUE.
	return

  end function compareLayers_modelParam


  ! *************************************************************************
  ! * Adjust the layer boundaries in the model parameter to match the grid
  ! * (for now, only used to make sure parametrization covers all of the grid
  ! * vertically, to the core-mantle boundary)
  ! * BOP
  subroutine adjustLayers_modelParam(P,r)

    implicit none
    type (modelParam_t), intent(inout)              :: P
    real(8), dimension(:), intent(in)               :: r
    ! * EOP
    integer											:: n
    real(8)											:: CMB ! core-mantle boundary

    n = size(r)
    CMB = r(n)

    if (trim(P%type) .ne. 'harmonic') then
        call warning('(adjustLayers) no layers to adjust for this parametrization type')
        return
    end if

    ! Test to make sure no grid is defined outside the layered region
	if (CMB < P%L(P%nL)%lbound) then
	  P%L(P%nL)%lbound = CMB
	end if

  end subroutine adjustLayers_modelParam


  ! **********************************************************************
  ! * BOP
  subroutine create_modelParam(P,nL,degree)

    implicit none
    type (modelParam_t), intent(inout)               :: P
    integer, intent(in)                            :: nL,degree
    ! * EOP

    character(80)                                  :: code
    integer                                        :: status
    integer                                        :: l,m,i,j,nV

    if (P%allocated) then
       call deall_modelParam(P)
    end if

    ! create spherical harmonic functionals
    nV=0
    do l=0,degree
      nV = nV + (2*l+1)
    end do
	allocate(P%F(nV),STAT=status)

	P%nF = nV
	P%F(:)%name='Y'

	i=1
	P%F(i)%num=i
	P%F(i)%l=0
	P%F(i)%m=0
	do l=1,degree
	  i=i+1
	  P%F(i)%num=i
	  P%F(i)%l=l
	  P%F(i)%m=0
	  do m=1,l
		i=i+1
  		P%F(i)%num=i
		P%F(i)%l=l
		P%F(i)%m=m
		i=i+1
		P%F(i)%num=i
		P%F(i)%l=l
		P%F(i)%m=-m
	  end do
	end do


	! create layer construction
	allocate(P%L(nL),STAT=status)
	P%nL = nL
	do i=1,nL
	  P%L(i)%num = i
	  P%L(i)%lbound = R_ZERO
	  P%L(i)%ubound = R_ZERO
	  P%L(i)%width  = R_ZERO
	  P%L(i)%depth  = R_ZERO
	  P%L(i)%alpha  = 0.0d0
	  P%L(i)%beta   = 1.0d0
	  P%L(i)%if_log =.FALSE.
	end do


	! create parameter construction
	allocate(P%c(nL,nV),STAT=status)

	do j=1,nL
	  do i=1,nV
		P%c(j,i)%L => P%L(j)
		P%c(j,i)%F => P%F(i)
		write(code,'(2i5)') P%F(i)%l,P%F(i)%m
		P%c(j,i)%code   = trim(code)
		P%c(j,i)%value  = R_ZERO
		P%c(j,i)%min    = R_ZERO
		P%c(j,i)%max    = R_ZERO
		P%c(j,i)%frozen =.TRUE.
		P%c(j,i)%exists =.FALSE.
		! initialize pointers
        if (j==1 .and. j==P%nL) then
          nullify(P%c(j,i)%prev)
          nullify(P%c(j,i)%next)
		else if (j==1) then
		  nullify(P%c(j,i)%prev)
		  P%c(j,i)%next => P%c(j+1,i)
		else if (j==P%nL) then
		  P%c(j,i)%prev => P%c(j-1,i)
		  nullify(P%c(j,i)%next)
		else
		  P%c(j,i)%prev => P%c(j-1,i)
		  P%c(j,i)%next => P%c(j+1,i)
		end if
	  end do
	end do

	P%nc = 0
	P%smoothed = .FALSE.
	P%zeroValued = .FALSE.
	P%allocated = .TRUE.

  end subroutine create_modelParam


  ! **********************************************************************
  ! * Can be used to create a new zero-valued parametrization from another
  ! * parametrization that have been already initialized
  ! * BOP
  subroutine setup_modelParam(P,L,F)

    implicit none
    type (modelParam_t), intent(inout)               :: P
	type (modelLayer_t), dimension(:), intent(in)    :: L
	type (modelFunc_t), dimension(:), intent(in)	   :: F
    ! * EOP

	if (P%allocated) then
	  call deall_modelParam(P)
	end if

	call create_modelParam(P,size(L),maxval(F%l))

	if (P%nF /= size(F)) then
	   call deall_modelParam(P)
	   call errStop('(setup_modelParam) spherical harmonics initialization incorrect')
	end if

	P%L = L
	P%F = F

  end subroutine setup_modelParam


  ! **********************************************************************
  ! * Creates a random perturbation and stores in a model parametrization;
  ! * extensively used for testing
  subroutine random_modelParam(P,eps)

    implicit none
    type (modelParam_t), intent(inout)               :: P
    real(8), intent(in), optional                    :: eps
    ! local
    real(8), allocatable :: dm(:)
    integer              :: istat


    if (.not. P%allocated) then
      call errStop('model parameter not allocated in random_modelParam')
    end if

    if (trim(P%type) .ne. 'harmonic') then
        call warning('(random_modelParam) model will not be randomised for this parametrization type')
        return
    end if

    ! this includes 'frozen' variables, but we're just testing
    allocate(dm(P%nc),STAT=istat)
    call random_number(dm)
    if (present(eps)) then
        dm = dm * eps
    else
        dm = dm * 0.05
    end if

    ! initialize and create the perturbation
    call zero_modelParam(P)
    call fillParamValues_modelParam(P,dm)
    deallocate(dm,STAT=istat)

  end subroutine random_modelParam


  ! **********************************************************************
	! * Insert new coefficients into an existing parametrization structure.
	! * Only replaces
  ! * BOP
  subroutine fillParam_modelParam(P,c)

    implicit none
    type (modelParam_t), intent(inout)               :: P
	type (modelCoeff_t), dimension(:), intent(in)    :: c
    ! * EOP

	integer										   :: i,j,iL,iV

    if(.not.P%allocated) then
       call errStop('(fillParam_modelParam) parametrization not allocated yet')
	end if

	do i=1,size(c)
	  iL = c(i)%L%num
	  iV = c(i)%F%num
	  if (P%c(iL,iV)%code /= c(i)%code) then
        call errStop('(fillParam_modelParam) parametrization setup error')
	  else if(.not.P%c(iL,iV)%exists) then
        call warning('(fillParam_modelParam) unable to replace coefficient - exists=.FALSE.')
		cycle
	  end if
	  P%c(iL,iV) = c(i)
	end do

	P%nc = count(P%c%exists)
	P%zeroValued = .false.


  end subroutine fillParam_modelParam


  ! **********************************************************************
  ! * BOP
  subroutine fillParamValues_modelParam(P,v)

    implicit none
    type (modelParam_t), intent(inout)               :: P
    real(8), dimension(:), intent(in)			   :: v
    ! * EOP

	integer										   :: i,j,k

    if(.not.P%allocated) then
       call errStop('(fillParamValues_modelParam) parametrization not allocated yet')
	end if

	if(count(P%c%exists) /= size(v)) then
       call errStop('(fillParamValues_modelParam) wrong number of input coefficients')
	end if

	k=0
	do j=1,P%nL
	  do i=1,P%nF
		if(.not.P%c(j,i)%exists) then
		  cycle
		end if
		k=k+1
		P%c(j,i)%value=v(k)
	  end do
	end do

	P%zeroValued = .false.

  end subroutine fillParamValues_modelParam

  ! **********************************************************************
  ! * P1 is the 1D (radial) part of P.
  ! * Always use a new variable for P1; should not overwrite P!!!
  function getRadial_modelParam(P) result (P1)

    implicit none
    type (modelParam_t), intent(in)				   :: P
    type (modelParam_t)							   :: P1
    integer							   			   :: lmax
    ! * EOP

	integer										   :: i,j,k

    if(.not.P%allocated) then
       call errStop('(getRadial_modelParam) parametrization not allocated yet')
	end if

    P1 = P

    if(trim(P%type) .ne. 'harmonic') then
       call warning('(getRadial_modelParam) unable to obtain the radial structure for this model type')
       return
    end if

	do j=1,P%nL
	  do i=1,P%nF
		if((P%F(i)%l == 0) .and. (P%F(i)%m == 0)) then
		  cycle
		end if
		P1%c(j,i)%value=R_ZERO
	  end do
	end do

    P1%crust%allocated = .FALSE.
    P1%zeroValued = P%zeroValued
	P1%temporary = .TRUE.

  end function getRadial_modelParam

  ! **********************************************************************
  ! * BOP
  function getDegree_modelParam(P) result (lmax)

    implicit none
    type (modelParam_t), intent(in)				   :: P
    integer							   			   :: lmax
    ! * EOP

	integer										   :: i,j,k

    if(.not.P%allocated) then
       call errStop('(getDegree_modelParam) parametrization not allocated yet')
	end if

	lmax = 0

	do i=1,P%nF
		if (P%F(i)%l > lmax) then
			lmax = P%F(i)%l
		end if
	end do

  end function getDegree_modelParam

  ! **********************************************************************
  ! * BOP
  function getLayerDegree_modelParam(P,iLayer) result (lmax)

    implicit none
    type (modelParam_t), intent(in)				   :: P
    integer, intent(in)				   			   :: iLayer
    integer							   			   :: lmax
    ! * EOP

	integer										   :: i,j,k

    if(.not.P%allocated) then
       call errStop('(getLayerDegree_modelParam) parametrization not allocated yet')
	end if

	lmax = 0

	do i=1,P%nF
		if ((.not.P%c(iLayer,i)%frozen) .and. (P%c(iLayer,i)%F%l > lmax)) then
			lmax = P%F(i)%l
		end if
	end do

  end function getLayerDegree_modelParam

  ! **********************************************************************
  ! * BOP
  subroutine getParamValues_modelParam(P,v)

    implicit none
    type (modelParam_t), intent(in)				   :: P
    real(8), dimension(:), intent(out)			   :: v
    ! * EOP

	integer										   :: i,j,k

    if(.not.P%allocated) then
       call errStop('(getParamValues_modelParam) parametrization not allocated yet')
	end if

	!do j=1,P%nL
	!  do i=1,P%nF
	!	if(P%c(j,i)%exists) then
	!	  print *, j,i,P%c(j,i)%value
	!	end if
	!  end do
	!end do


	if(count(P%c%exists) /= size(v)) then
       write(0,*) '(getParamValues_modelParam) wrong number of input coefficients: ',size(v)
       stop
	end if

	k=0
	do j=1,P%nL
	  do i=1,P%nF
		if(.not.P%c(j,i)%exists) then
		  cycle
		end if
		k=k+1
		v(k)=P%c(j,i)%value
	  end do
	end do

  end subroutine getParamValues_modelParam


  ! **********************************************************************
  ! * BOP
  subroutine setLayer_modelParam(P,iL,upperb,lowerb,alpha,beta,if_log)

    implicit none
    type (modelParam_t), intent(inout)               :: P
	integer, intent(in)							   :: iL
    real(8), intent(in)							   :: upperb,lowerb
    real(8), intent(in)							   :: alpha,beta
	logical, intent(in)							   :: if_log
    ! * EOP

	real(8)										   :: depth, width
    integer                                        :: status

    if(.not.P%allocated) then
       call errStop('(setLayer_modelParam) parametrization not allocated yet')
	end if

    depth = EARTH_R - lowerb
	width  = upperb - lowerb
	if (width <= R_ZERO) then
       write(0,*) '(setLayer_modelParam) incorrect depth of layer ',iL
       stop
	end if

	P%L(iL)%depth  = depth
	P%L(iL)%width  = width
	P%L(iL)%lbound = lowerb
	P%L(iL)%ubound = upperb
	P%L(iL)%alpha  = alpha
	P%L(iL)%beta   = beta
	P%L(iL)%if_log = if_log


  end subroutine setLayer_modelParam


  ! **********************************************************************
  ! * BOP
  subroutine setCoeffValue_modelParam(P,iL,l,m,v,min,max,frozen)

    implicit none
    type (modelParam_t), intent(inout)               :: P
	integer, intent(in)							   :: iL,l,m
    real(8), intent(in)							   :: v,min,max
	logical, intent(in)							   :: frozen
    ! * EOP

	integer										   :: iV

    if(.not.P%allocated) then
       call errStop('(setCoeffValue_modelParam) parametrization not allocated yet')
	end if

	iV = iY(l,m)
	if (P%c(iL,iV)%exists) then
       write(0,*) '(setCoeffValue_modelParam) coefficient ',iL,iV,' already exists, overwrite...'
	end if

	if ((P%c(iL,iV)%F%l /= l).or.(P%c(iL,iV)%F%m /= m)) then
       call errStop('(setCoeffValue_modelParam) major error in parametrization setup')
	end if

	P%c(iL,iV)%value = v
	P%c(iL,iV)%min = min
	P%c(iL,iV)%max = max
	P%c(iL,iV)%frozen = frozen

	P%nc = P%nc + 1
	P%zeroValued = .false.

	P%c(iL,iV)%exists = .TRUE.

  end subroutine setCoeffValue_modelParam


  ! **********************************************************************
  ! * BOP
  subroutine getCoeffValue_modelParam(P,iL,l,m,v,min,max,frozen)

    implicit none
    type (modelParam_t), intent(in)				   :: P
	integer, intent(in)							   :: iL,l,m
    real(8), intent(out)                           :: v
    real(8), intent(out), optional			       :: min,max
	logical, intent(out), optional			       :: frozen
    ! * EOP

	integer										   :: iV

    if(.not.P%allocated) then
       call errStop('(getCoeffValue_modelParam) parametrization not allocated yet')
	end if

	iV = iY(l,m)
	if (.not.P%c(iL,iV)%exists) then
	   call warning('(getCoeffValue_modelParam) required coefficient does not exist')
	end if

	v = P%c(iL,iV)%value

	if (present(min) .and. present(max)) then
	    min = P%c(iL,iV)%min
	    max = P%c(iL,iV)%max
	end if

    if (present(frozen)) then
	    frozen = P%c(iL,iV)%frozen
	end if

  end subroutine getCoeffValue_modelParam


  ! **********************************************************************
  ! * BOP
  function getCoeff_modelParam(P,n) result (c)

    implicit none
    type (modelParam_t), intent(in)				   :: P
	integer, intent(in)							   :: n
	type (modelCoeff_t)							   :: c
    ! * EOP

	integer										   :: i,j,k

    if(.not.P%allocated) then
       call errStop('(getCoeff_modelParam) parametrization not allocated yet')
	end if

	k=0
	search: do j=1,P%nL
	  do i=1,P%nF
		if(P%c(j,i)%exists) then
		  k=k+1
		  if(k==n) then
			exit search
		  end if
		end if
	  end do
	end do search

	c=P%c(j,i)

  end function getCoeff_modelParam

  ! **********************************************************************
  ! * BOP
  subroutine getCoeffArray_modelParam(P,c)

    implicit none
    type (modelParam_t), intent(in)				   :: P
	type (modelCoeff_t), dimension(:), intent(out)   :: c
    ! * EOP

	integer										   :: i,j,k

    if(.not.P%allocated) then
       call errStop('(getCoeffArray_modelParam) parametrization not allocated yet')
	end if

	if(size(c) /= P%nc) then
	   call errStop('(getCoeffArray_modelParam) coefficient array is of a wrong size')
	end if


	k=0
	search: do j=1,P%nL
	  do i=1,P%nF
		if(P%c(j,i)%exists) then
		  k=k+1
			c(k)=P%c(j,i)
		end if
	  end do
	end do search

  end subroutine getCoeffArray_modelParam

  ! **********************************************************************
  ! * count_modelParam: counts the number of variable model parameters
  ! * BOP
  function count_modelParam(P) result (N)

    implicit none
    type (modelParam_t), intent(in)                :: P
    integer                                        :: N
    ! * EOP

    if(.not.P%allocated) then
       call errStop('(count_modelParam) parametrization not allocated yet')
    end if

    N = count(P%c%exists)

  end function count_modelParam

  ! **********************************************************************
  ! * BOP
  subroutine copy_modelParam(P2,P1)

    implicit none
    type (modelParam_t), intent(in)                  :: P1
    type (modelParam_t), intent(inout)               :: P2
    ! * EOP

    integer                                        :: status
	integer										   :: i,j
	integer                                        :: nx1,ny1,nx2,ny2

    ! if output is allocated and of different size, reallocate
	if (P2%allocated) then
	  if (.not. compare(P2,P1)) then
	    call deall_modelParam(P2)
	  end if
	end if

	! copy type
	P2%type = P1%type

	! copy functional construction
	P2%nF = P1%nF
	if (.not. associated(P2%F)) allocate(P2%F(P2%nF),STAT=status)
	P2%F = P1%F

	! copy layer construction
	P2%nL = P1%nL
	if (.not. associated(P2%L)) allocate(P2%L(P2%nL),STAT=status)
	P2%L = P1%L

	! copy parameter construction
	P2%nc = P1%nc
	if (.not. associated(P2%c)) allocate(P2%c(P2%nL,P2%nF),STAT=status)
	P2%c = P1%c

	! copy values of a grid, if exist
	if (P1%rho%allocated) then
	    P2%rho = P1%rho
	end if

	! copy crust, if exists
	if (P1%crust%allocated) then
		nx1 = size(P1%crust%cond,1)
		ny1 = size(P1%crust%cond,2)
		if (associated(P2%crust%cond)) then
			nx2 = size(P1%crust%cond,1)
			ny2 = size(P1%crust%cond,2)
			if ((nx1 .ne. nx2) .or. (ny1 .ne. ny2)) then
	   			deallocate(P2%crust%cond,STAT=status)
	   			allocate(P2%crust%cond(nx1,ny1),STAT=status)
	   		end if
	   	else
	   		allocate(P2%crust%cond(nx1,ny1),STAT=status)
		end if
		P2%crust%cond = P1%crust%cond
		P2%crust%allocated = .true.
	else
		P2%crust%allocated = .false.
	end if

	do j=1,P2%nL
	  do i=1,P2%nF
		P2%c(j,i)%L => P2%L(j)
		P2%c(j,i)%F => P2%F(i)
	  end do
	end do

    P2%smoothed  = P1%smoothed
    P2%zeroValued = P1%zeroValued
	P2%allocated = P1%allocated

	if (P1%temporary) then
		call deall_modelParam(P1)
	end if

  end subroutine copy_modelParam


  ! **********************************************************************
  ! * A function is more logical here, but the subroutine is much less
  ! * error-prone, since m = zero(m) might create trouble.
  ! * BOP
  subroutine zero_modelParam(P)

    implicit none
    type (modelParam_t), intent(inout)	:: P
    ! * EOP

    if(.not.(P%allocated)) then
    	call errStop('(zero_modelParam) input parametrization not allocated yet')
	end if

	if (P%rho%allocated) then
	    call zero_rscalar(P%rho)
	end if

	! Set all values to zero
	P%c%value = R_ZERO

	! Update the logical
	P%zeroValued = .true.

	! All logicals except allocated and zeroValued should be false again
	P%smoothed = .false.

  end subroutine zero_modelParam

  ! **********************************************************************
  ! * Get zeroValued logical from model parameter
  ! * BOP
  logical function iszero_modelParam(P) result (zeroValued)

    implicit none
    type (modelParam_t), intent(inout)  :: P
    ! * EOP

    zeroValued = P%zeroValued

  end function iszero_modelParam

  ! **********************************************************************
  ! * We can add or subtract values with exists==.FALSE. (they equal zero)
  ! * BOP
  function add_modelParam(P1,P2) result(P)

    implicit none
    type (modelParam_t), intent(in)		:: P1
    type (modelParam_t), intent(in)		:: P2
    type (modelParam_t)					:: P
    ! * EOP

    integer					:: i,j

	if (.not. compare(P1,P2)) then
		call errStop('(add_modelParam) input parametrization structures incompatible')
	end if

	P = P1

    if (P1%rho%allocated .and. P2%rho%allocated) then
        call add_rscalar(P1%rho,P2%rho,P%rho)
        P%temporary = .true.
        return
    end if

	do j=1,P%nL
	  do i=1,P%nF
	    if((P1%c(j,i)%F%l==P2%c(j,i)%F%l).and.(P1%c(j,i)%F%m==P2%c(j,i)%F%m)) then
		  P%c(j,i)%value = P1%c(j,i)%value + P2%c(j,i)%value
	    else
	      call errStop('(add_modelParam) parametrization not set up correctly')
	    end if
	  end do
	end do

	P%temporary = .true.
	P%zeroValued = P1%zeroValued .and. P2%zeroValued

  end function add_modelParam


  ! **********************************************************************
  ! * We can add or subtract values with exists==.FALSE. (they equal zero)
  ! * BOP
  function subtract_modelParam(P1,P2) result(P)

    implicit none
    type (modelParam_t), intent(in)		:: P1
    type (modelParam_t), intent(in)		:: P2
    type (modelParam_t)					:: P
    ! * EOP

    integer					:: i,j

	if (.not. compare(P1,P2)) then
		call errStop('(subtract_modelParam) input parametrization structures incompatible')
	end if

	P = P1

    if (P1%rho%allocated .and. P2%rho%allocated) then
        call subtract_rscalar(P1%rho,P2%rho,P%rho)
        P%temporary = .true.
        return
    end if

	do j=1,P%nL
	  do i=1,P%nF
	    if((P1%c(j,i)%F%l==P2%c(j,i)%F%l).and.(P1%c(j,i)%F%m==P2%c(j,i)%F%m)) then
		  P%c(j,i)%value = P1%c(j,i)%value - P2%c(j,i)%value
	    else
	      call errStop('(subtract_modelParam) parametrization not set up correctly')
	    end if
	  end do
	end do

	P%temporary = .true.
    P%zeroValued = P1%zeroValued .and. P2%zeroValued

  end function subtract_modelParam


  ! **********************************************************************
  ! * We can multiply values with exists==.FALSE. (they equal zero)
  ! * BOP
  function mult_modelParam(P1,P2) result(P)

    implicit none
    type (modelParam_t), intent(in)		:: P1
    type (modelParam_t), intent(in)		:: P2
    type (modelParam_t)					:: P
    ! * EOP

    integer					:: i,j

	if (.not. compare(P1,P2)) then
		call errStop('(mult_modelParam) input parametrization structures incompatible')
	end if

	P = P1

    if ((trim(P1%type) .ne. 'harmonic') .or. (trim(P2%type) .ne. 'harmonic')) then
        call warning('(mult_modelParam) no multiplication will be performed for this parametrization type')
        return
    end if

	do j=1,P%nL
	  do i=1,P%nF
	    if((P1%c(j,i)%F%l==P2%c(j,i)%F%l).and.(P1%c(j,i)%F%m==P2%c(j,i)%F%m)) then
		  P%c(j,i)%value = P1%c(j,i)%value * P2%c(j,i)%value
	    else
	      call errStop('(mult_modelParam) parametrization not set up correctly')
	    end if
	  end do
	end do

	P%temporary = .true.
    P%zeroValued = P1%zeroValued .and. P2%zeroValued

  end function mult_modelParam

  ! **********************************************************************
  ! * We can multiply values with exists==.FALSE. (they equal zero)
  ! * BOP
  function dotProd_modelParam(P1,P2) result(r)

    implicit none
    type (modelParam_t), intent(in)		:: P1
    type (modelParam_t), intent(in)		:: P2
    real(8)								:: r
    ! * EOP

    integer					:: i,j

	if (.not. compare(P1,P2)) then
		call errStop('(dotProd_modelParam) parametrization structures incompatible')
	end if

    r = R_ZERO

    if ((trim(P1%type) .ne. 'harmonic') .or. (trim(P2%type) .ne. 'harmonic')) then
        call warning('(dotProd_modelParam) no multiplication will be performed for this parametrization type')
        return
    end if

    if (P1%zeroValued .and. P2%zeroValued) then
        return
    end if

	do j=1,P2%nL
	  do i=1,P2%nF
	    if((P1%c(j,i)%F%l==P2%c(j,i)%F%l).and.(P1%c(j,i)%F%m==P2%c(j,i)%F%m)) then
		  if(.not.P1%c(j,i)%frozen) then
			r = r + P1%c(j,i)%value * P2%c(j,i)%value
		  end if
	    else
	      call errStop('(dotProd_modelParam) parametrization not set up correctly')
	    end if
	  end do
	end do

  end function dotProd_modelParam


  ! **********************************************************************
  ! * If we multiply by a vector we only use values with exists==.TRUE.
  ! * BOP
  function dotProdVec_modelParam(v,P2) result(r)

    implicit none
    type (modelParam_t), intent(in)		:: P2
    real(8), dimension(:), intent(in)	:: v
    real(8)								:: r
    ! * EOP

    integer					:: i,j,k

    if(.not.P2%allocated) then
    	call errStop('(dotProdVec_modelParam) parametrization not allocated yet')
	end if

	if(count(P2%c%exists) /= size(v)) then
		call errStop('(dotProdVec_modelParam) wrong number of input coefficients')
	end if

	if (trim(P2%type) .ne. 'harmonic') then
        call warning('(dotProdVec_modelParam) no multiplication will be performed for this parametrization type')
        return
    end if

    r = R_ZERO
	k = 0

	do j=1,P2%nL
	  do i=1,P2%nF
	    if(P2%c(j,i)%exists) then
		  if(P2%c(j,i)%frozen) then
		  else
			r = r + v(k) * P2%c(j,i)%value
		  end if
		  k = k+1
	    else
		  cycle
	    end if
	  end do
	end do

  end function dotProdVec_modelParam

  ! **********************************************************************
  ! * BOP
  subroutine scMult_modelParam(v,P1,P)

    implicit none
    type (modelParam_t), intent(in)		:: P1
    real(8), intent(in)					:: v
    type (modelParam_t), intent(inout)	:: P
    ! * EOP

	if (.not. P1%allocated) then
		call errStop('(scMult_modelParam) input parametrization not allocated yet')
	else if (.not. P%allocated) then
		call errStop('(scMult_modelParam) output structure has to be allocated before calling')
	else if (.not. compare(P1,P)) then
		call errStop('(scMult_modelParam) output parametrization is incompatible')
	end if

    if (P1%rho%allocated) then
        call scMult_rscalar(v,P1%rho,P%rho)
    end if

	P%c%value = v * P1%c%value
	P%zeroValued = P1%zeroValued
	P%smoothed = P1%smoothed


  end subroutine scMult_modelParam

  ! **********************************************************************
  ! * BOP
  subroutine scMultAdd_modelParam(v,P1,P)

    implicit none
    type (modelParam_t), intent(in)     :: P1
    real(8), intent(in)                 :: v
    type (modelParam_t), intent(inout)  :: P
    ! * EOP

    if (.not. P1%allocated) then
        call errStop('(scMult_modelParam) input parametrization not allocated yet')
    else if (.not. P%allocated) then
        call errStop('(scMult_modelParam) output structure has to be allocated before calling')
    else if (.not. compare(P1,P)) then
        call errStop('(scMult_modelParam) output parametrization is incompatible')
    end if

    if (P1%rho%allocated) then
        call scMultAdd_rscalar(v,P1%rho,P%rho)
    end if

    P%c%value = P%c%value + v * P1%c%value
    P%zeroValued = P%zeroValued .and. P1%zeroValued
    P%smoothed = P1%smoothed


  end subroutine scMultAdd_modelParam

  ! **********************************************************************
  ! * BOP
  function scDiv_modelParam(P1,v) result(P)

    implicit none
    type (modelParam_t), intent(in)		:: P1
    real(8), intent(in)					:: v
    type (modelParam_t)					:: P
    ! * EOP

	if (.not. P1%allocated) then
		call errStop('(scDiv_modelParam) input parametrization not allocated yet')
	end if

	P = P1

	if(v==R_ZERO) then
	   write(0,*) 'Error: (scDiv_modelParam) division by zero'
	   return
	end if

	P%c%value = P1%c%value / v
	P%zeroValued = P1%zeroValued
    P%smoothed = P1%smoothed
	P%temporary = .true.

  end function scDiv_modelParam

  ! **********************************************************************
  ! * We can add or subtract values with exists==.FALSE. (they equal zero)
  ! *
  ! * Output has to be allocated and of the right size before calling.
  ! * This is necessary to ensure we don't modify the inputs in a call
  ! * such as m1 = a*m1 + b*m2 (thus allowing to overwrite input with output).
  ! * BOP
  subroutine linComb_modelParam(r1,P1,r2,P2,P)

     !  forms the linear combination of model parameters
     !    P = r1*P1 + r2*P2
     !  where r1 and r2 are real constants and P1 and P2
     !   are model parameters
     !   output P may overwrite P1 or P2

    implicit none
    type (modelParam_t), intent(in)				   :: P1
    type (modelParam_t), intent(in)				   :: P2
    real(8), intent(in)							   :: r1
    real(8), intent(in)							   :: r2
    type (modelParam_t), intent(inout)             :: P
    ! * EOP

	integer										   :: i,j

	if (P1%rho%allocated .and. P2%rho%allocated) then
	    call linComb_rscalar(r1,P1%rho,r2,P2%rho,P%rho)
	    return
	end if

	if (.not. compare(P1,P2)) then
		call errStop('(linComb_modelParam) input parametrization structures incompatible')
	else if (.not. P%allocated) then
		call errStop('(linComb_modelParam) output structure has to be allocated before calling')
	else if (.not. compare(P1,P)) then
		call errStop('(linComb_modelParam) output parametrization is incompatible')
	end if

	do j=1,P%nL
	  do i=1,P%nF
	    if((P1%c(j,i)%F%l==P2%c(j,i)%F%l).and.(P1%c(j,i)%F%m==P2%c(j,i)%F%m)) then
		  P%c(j,i)%value = r1 * P1%c(j,i)%value + r2 * P2%c(j,i)%value
	    else
	      call errStop('(linComb_modelParam) parametrization not set up correctly')
	    end if
	  end do
	end do

    P%zeroValued = P1%zeroValued .and. P2%zeroValued
	P%smoothed = P1%smoothed .or. P2%smoothed

  end subroutine linComb_modelParam


  ! **********************************************************************
  ! * BOP
  subroutine smoothH_modelParam(P)

    implicit none
    type (modelParam_t), intent(inout)               :: P
    ! * EOP

	integer										   :: i,j

    if(.not.P%allocated) then
       call warning('(smoothH_modelParam) parametrization not allocated yet')
	   return
	end if

	do j=1,P%nL
	  do i=1,P%nF
			P%c(j,i)%value = ((P%F(i)%l+1)**(-(P%L(j)%alpha)/2)) * P%c(j,i)%value
	  end do
	end do

  end subroutine smoothH_modelParam


  ! **********************************************************************
  ! * Values with exists==.FALSE. have been initialized to zero and have
  ! * a smoothing effect: neighbouring parameters will have a smaller value
  ! * BOP
  subroutine smoothV_modelParam(P)

    implicit none
    type (modelParam_t), intent(inout)               :: P
    ! * EOP

	integer										   :: i,j
	real(8),dimension(P%nL,P%nF)				   :: v,v_prev,v_next

    if(.not.P%allocated) then
       call warning('(smoothV_modelParam) parametrization not allocated yet')
	   return
	end if


	do j=1,P%nL
	  do i=1,P%nF

			if(associated(P%c(j,i)%prev).and.associated(P%c(j,i)%next)) then
			v_prev(j,i) = P%c(j,i)%prev%value * (1-P%L(j-1)%beta)/2
			v(j,i)      = P%c(j,i)%value * P%L(j)%beta
			v_next(j,i) = P%c(j,i)%next%value * (1-P%L(j)%beta)/2
			else if (associated(P%c(j,i)%next)) then
			v_prev(j,i) = R_ZERO
			v(j,i)      = P%c(j,i)%value * P%L(j)%beta
			v_next(j,i) = P%c(j,i)%next%value * (1-P%L(j)%beta)/2
			else if (associated(P%c(j,i)%prev)) then
			v_prev(j,i) = P%c(j,i)%prev%value * (1-P%L(j-1)%beta)/2
			v(j,i)      = P%c(j,i)%value * P%L(j)%beta
			v_next(j,i) = R_ZERO
			else
			v_prev(j,i) = R_ZERO
			v(j,i)      = P%c(j,i)%value
			v_next(j,i) = R_ZERO
			end if

	  end do
	end do


	P%c%value = v_prev + v + v_next


  end subroutine smoothV_modelParam


  ! **********************************************************************
  ! * Preconditioning operator C_p^{1/2}
  ! * One of the few occasions where it is illegal to overwrite the input
  ! * i.e. no calls like m = multBy_CmSqrt(m) allowed
  ! * BOP
  function multBy_CmSqrt(P1) result (P)

    implicit none
    type (modelParam_t), intent(in)				   :: P1
    type (modelParam_t)							   :: P
    ! * EOP

	if (.not. P1%allocated) then
		call errStop('(multBy_CmSqrt) input parametrization not allocated yet')
	else if (P%allocated) then
		call errStop('(multBy_CmSqrt) output structure cannot to be allocated before calling')
	end if

	P = P1

	if (trim(P%type) .ne. 'harmonic') then
        call warning('(multBy_CmSqrt) no smoothing will be performed for this parametrization type')
        P%smoothed = .false.
        P%temporary = .true.
        return
    end if

	! make the smoothing operator symmetric
	if (.not. P%zeroValued) then
	   call smoothV_modelParam(P)
	   call smoothH_modelParam(P)
	   call smoothV_modelParam(P)
	end if

	P%smoothed = .true.
	P%temporary = .true.

  end function multBy_CmSqrt


  ! **********************************************************************
  ! * Preconditioning operator C_p
  ! * One of the few occasions where it is illegal to overwrite the input
  ! * i.e. no calls like m = multBy_CmSqrt(m) allowed
  ! * BOP
  function multBy_Cm(P1) result (P)

    implicit none
    type (modelParam_t), intent(in)				   :: P1
    type (modelParam_t)							   :: P2
    type (modelParam_t)							   :: P
    ! * EOP

	if (.not. P1%allocated) then
		call errStop('(multBy_Cm) input parametrization not allocated yet')
	else if (P%allocated) then
		call errStop('(multBy_Cm) output structure cannot be allocated before calling')
	end if

	! apply operator C_p^{1/2} twice
	P2 = multBy_CmSqrt(P1)
	P  = multBy_CmSqrt(P2)

	call deall_modelParam(P2)

	P%temporary = .true.

  end function multBy_Cm

  ! **********************************************************************
  ! * BOP
  subroutine setCrust_modelParam(crust,P)

    type (modelShell_t), intent(in)          :: crust
    type (modelParam_t), intent(inout)       :: P
    ! * EOP
    integer                                  :: nx,ny,status

    if(.not.P%allocated) then
       call errStop('(setCrust_modelParam) parametrization not allocated yet')
	end if

	if(.not.crust%allocated) then
	   ! crust not allocated intentionally; do nothing
	   P%crust%allocated = .false.
	   return
	end if

	if(.not. associated(crust%cond)) then
	   call errStop('(setCrust_modelParam) crust not allocated yet')
	end if

	if(associated(P%crust%cond)) then
	   deallocate(P%crust%cond, STAT=status)
	end if

    nx = size(crust%cond,1)
    ny = size(crust%cond,2)
	allocate(P%crust%cond(nx,ny), STAT=status)

    P%crust%cond = crust%cond
    P%crust%allocated = .true.

  end subroutine setCrust_modelParam

  ! **********************************************************************
  ! * BOP
  subroutine print_modelParam(P,verbose,comment)

    implicit none
    type (modelParam_t), intent(in)         :: P
	integer, intent(in)                   	:: verbose
	character(*),intent(in),optional		:: comment
    ! * EOP

	type (modelCoeff_t)						:: coeff
	integer									:: i,j

    if(.not.P%allocated) then
       call warning('(print_modelParam) parametrization not allocated yet')
       return
    else if(trim(P%type) .ne. 'harmonic') then
       call warning('(print_modelParam) skipped printing this model parametrization type')
	   return
	end if

	write(0,*)
	if(present(comment)) then
		write(0,*) node_info,comment
	end if

	if (verbose>0) then
  		write(0,'(a12,a50,i4)') node_info,'Number of layers in script: ',P%nL
		do j=1,P%nL
    		write(0,'(a12,a46,i2,a2,i4)') node_info,'Number of coefficients in layer ',j,': ',count(.not.P%c(j,:)%frozen)
		end do
  		write(0,'(a12,a50,i4)') node_info,'Number of variable parameters in script: ',count(.not.P%c%frozen)
		write(0,*)
	end if

	do j=1,P%nL
    	write(0,'(a12,a46,i2,a2,g15.7)') node_info,'Degree and order zero coefficient in layer ',j,': ',P%c(j,1)%value
		end do
	write(0,*)

	if (verbose>0) then
		do j=1,P%nL
    		write(0,'(a12,a46,i2,a2,g10.5)') node_info,'Horizontal regularization in layer ',j,': ',P%L(j)%alpha
		end do
		write(0,*)
	end if

	if (verbose>0) then
		do j=1,P%nL
    		write(0,'(a12,a46,i2,a2,g10.5)') node_info,'Vertical regularization in layer ',j,': ',P%L(j)%beta
		end do
		write(0,*)
	end if

	if (verbose>3) then
		do i=1,P%nc
	  		coeff=getCoeff_modelParam(P,i)
	  		write (0,'(a12,2i8,g17.9)') node_info,coeff%L%num,coeff%F%num,coeff%value
		end do
	end if

  end subroutine print_modelParam


  ! **********************************************************************
  ! * This is an essential subroutine that is used to output the model
  ! * solution in e.g. NLCG. This solution is then used as the prior/initial
  ! * model, or to compute the fields.
  ! * We output both the final "smooth" model m = Cm^{1/2} mhat + m0,
  ! * and the "rough" model mhat, which is the result of the inverse search.
  ! * The smooth model m is best used to save and to plot the model; while
  ! * the rough model mhat is best used in conjunction with the prior to
  ! * restart the inversion, or compute the responses / derivatives.
  ! * Therefore, to avoid confusion we DO NOT write the regularization
  ! * parameters for a model if the smoothing has already been applied.
  ! * Regularization parameters are only needed to smooth a rough model,
  ! * so we leave them out for m (and leave them in for mhat).
  subroutine write_modelParam(P,cfile)

    implicit none
    type (modelParam_t), intent(in)         :: P
	character(*), intent(in)				:: cfile

	integer									:: lmax,i,j,k,istat
    character(6)							:: if_log_char,if_var_char

    if(.not.P%allocated) then
       call warning('(write_modelParam) parametrization not allocated yet')
       return
    else if(trim(P%type) .ne. 'harmonic') then
       call warning('(write_modelParam) skipped writing this model parametrization type')
       return
    end if

	open(unit=ioPrm, file=cfile, status='unknown', iostat=istat)

	lmax = getDegree(P)

	write(ioPrm,'(a24,i2,a8,i2)') 'Format: harmonic layers ',P%nL,' degree ',lmax
	write(ioPrm,*)

	do j=1,P%nL
		if (P%L(j)%if_log) then
			if_log_char = 'log'
		else
			if_log_char = 'linear'
		end if
		lmax = getDegree(P,j)
		write(ioPrm,'(a6)',advance='no') if_log_char
		write(ioPrm,'(a8,i2)',advance='no') ' degree ',lmax
		write(ioPrm,'(a7,g10.5)',advance='no') ' layer ',P%L(j)%depth
		if (P%smoothed) then
          ! OUTPUT SMOOTH MODEL SO DO NOT WRITE THE REGULARISATION PARAMETERS...
		  ! INSTEAD, START NEW LINE:
		  write(ioPrm,*)
		else
		  write(ioPrm,'(a5,2g10.5)') ' reg ',P%L(j)%alpha, P%L(j)%beta
		end if
		write(ioPrm,*) '  l   m   value  	  min     	max'
		do i=1,P%nF
			if (.not.P%c(j,i)%exists) then
				cycle
			end if
			write(ioPrm,'(2i4,g15.7)',advance='no') P%F(i)%l,P%F(i)%m,P%c(j,i)%value
			if (P%c(j,i)%frozen) then
				if_var_char = 'const'
			else
				if_var_char = 'range'
			end if
			write(ioPrm,'(2g15.7,a6)') P%c(j,i)%min,P%c(j,i)%max,if_var_char
		end do
    	write(ioPrm,*)
	end do
	write(ioPrm,*)

	close(ioPrm)

  end subroutine write_modelParam


  ! **********************************************************************
  ! * Use to output a vector of model parameters, normally for the full
  ! * sensitivity matrix (NOTE ioSens HARDCODED).
  ! * Assume that the vector consists of parameters of the same type
  ! * (size, etc); could be verified on input.
  ! * Assume also that the file is already open; will be closed in the
  ! * calling subroutine.
  subroutine writeVec_modelParam(np,P,header,cfile)

    implicit none
    integer, intent(in)                     :: np
    type (modelParam_t), intent(in)         :: P(np)
    character(*), intent(in)                :: header, cfile
    ! * EOP

    integer                                 :: lmax,i,j,k,istat
    logical                                 :: opened
    character(6)                            :: if_log_char,if_var_char

    inquire(file=cfile, opened=opened)
    if (.not. opened) then
        open(unit=ioSens, file=cfile, status='unknown', iostat=istat)
    endif

    lmax = getDegree(P(1))

    write(ioSens,'(a24,i2,a8,i2)') 'Format: harmonic layers ',P(1)%nL,' degree ',lmax
    write(ioSens,*)

    do j=1,P(1)%nL
        if (P(1)%L(j)%if_log) then
            if_log_char = 'log'
        else
            if_log_char = 'linear'
        end if
        lmax = getDegree(P(1),j)
        write(ioSens,'(a6)',advance='no') if_log_char
        write(ioSens,'(a8,i2)',advance='no') ' degree ',lmax
        write(ioSens,'(a7,g10.5)',advance='no') ' layer ',P(1)%L(j)%depth
        if (P(1)%smoothed) then
          ! OUTPUT SMOOTH MODEL SO DO NOT WRITE THE REGULARISATION PARAMETERS...
          ! INSTEAD, START NEW LINE:
          write(ioPrm,*)
        else
          write(ioPrm,'(a5,2g10.5)') ' reg ',P(1)%L(j)%alpha, P(1)%L(j)%beta
        end if
        write(ioSens,*) '  l   m   value       min       max'
        do i=1,P(1)%nF
            if (.not.P(1)%c(j,i)%exists) then
                cycle
            end if
            write(ioSens,'(2i4)',advance='no') P(1)%F(i)%l,P(1)%F(i)%m
            do k=1,size(P)
                write(ioSens,'(g15.7)',advance='no') P(k)%c(j,i)%value
            end do
            if (P(1)%c(j,i)%frozen) then
                if_var_char = 'const'
            else
                if_var_char = 'range'
            end if
            write(ioSens,'(2g15.7,a6)') P(1)%c(j,i)%min,P(1)%c(j,i)%max,if_var_char
        end do
        write(ioSens,*)
    end do
    write(ioSens,*)

    !... do not close!
    !close(ioSens)

  end subroutine writeVec_modelParam

  ! ***************************************************************************
  ! * read_modelParam reads the parametrization info to store in the derived data
  ! * type variables (modelCoeff_t, modelLayer_t) in the special case when the
  ! * parametrization is in terms of spherical harmonics of various degree/order
  ! * per layer, each of the parameters being either variable (keyword 'range') or
  ! * constant (keyword 'const'). We keep this information internally in such
  ! * a format, that generalizations of this case would be easy to implement.
  ! * In the parametrization type variables, each layer has the same (maximum)
  ! * degree and order of spherical harmonics, and hence identical numbers of
  ! * coefficients to store. Out of these coefficients, those that have been
  ! * defined in the script bear logical keyword 'exists'. Out of those, variable
  ! * coefficients have logical 'frozen==.FALSE.'
  ! * Obviously, the information on the range is not required for the forward
  ! * solver to operate. This is provided for the inversion, which will share
  ! * the same input format for now.

  subroutine read_modelParam(P,cfile)

    character(*), intent(in)            :: cfile
    type (modelParam_t), intent(inout)  :: P
    integer                             :: ilayer,i,j,k,n,l,m
    integer                             :: nF,nL
    integer                             :: sum,sum0,degree
    integer                             :: ios,istat
    real(8)                             :: upperb,lowerb,width,depth,alpha,beta
    character(6)                        :: if_log_char,if_var_char
    logical                             :: if_log, if_fixed, exists
    character(80)                       :: prmname, string
    real(8)                             :: v,min,max

    lowerb = EARTH_R
    depth = 0.0d0
    width = 0.0d0
    sum = 0
    sum0 = 0

    inquire(FILE=trim(cfile),EXIST=exists)
    if(exists) then
      open(ioPrm,file=cfile,status='old',form='formatted',iostat=ios)
      write(6,*) node_info,'Reading from the parametrization file ',trim(cfile)
    else
      write(0,*) node_info,'Error: (read_modelParam) input file does not exist'
      stop
    end if

    read(ioPrm,'(a8,a80)') string,prmname

    if (index(prmname,'harmonic')==0) then
       write(0, *) node_info,'Error: (read_modelParam) not a spherical harmonic parametrization'
       stop
    else
      i = index(prmname,'layers')
      read(prmname(i+7:len(prmname)),'(i2)') nL
      i = index(prmname,'degree')
      read(prmname(i+7:len(prmname)),'(i2)') degree
      write(6,*) node_info,trim(prmname)
    end if


    call create_modelParam(P,nL,degree)

    !write(6,'(a50,i3)') 'Number of layers in script: ',P%nL


    ! Unwind spherical harmonic parametrization

    do n=1,nL

       read(ioPrm, *)

       upperb = lowerb

       ! If you wish to read the line in sections, use advance='no' specifier
     !  read(ioPrm,'(a3,1x,a6,1x,i2,a6,g15.7)',iostat=ios) &
            !if_log_char,string,degree,string,depth

            !read(ioPrm,'(a3)',iostat=ios,advance='no') if_log_char

            read(ioPrm,'(a80)',iostat=ios) string
            i = index(string,'degree')
            j = index(string,'layer')
            k = index(string,'reg')


            if (k==0) then
                    ! no regularisation specified for this layer
                    alpha = 0.0d0
                    beta  = 1.0d0
                    k = len(string)
            else
                    read(string(k+4:len(string)),*) alpha,beta
            end if

            read(string(1:i-1),*) if_log_char
            read(string(i+7:j),*) degree
            read(string(j+6:k),*) depth

       read(ioPrm,'(a80)',iostat=ios) string

       if (if_log_char == 'log') then  ! log means log_{10}
          if_log = .TRUE.
       else
          if_log = .FALSE.
       end if
       lowerb = EARTH_R - depth

       call setLayer_modelParam(P,n,upperb,lowerb,alpha,beta,if_log)

       sum0=sum0+sum
       sum=0
       do l=0,degree
          sum = sum + (2*l+1)
       end do
       !write(6,'(a46,i2,a2,i3)') 'Number of coefficients in layer ',n,': ',sum

       do i=1,sum
          read(ioPrm,*,iostat=ios) l,m,v,min,max,if_var_char
          if (if_var_char == 'range') then
            if_fixed = .FALSE.
          else if (if_var_char == 'const') then
            if_fixed = .TRUE.
          else
            write(0, *) 'Error: (read_modelParam) wrong character constant for ',n,l,m,': ',if_var_char
            write(0, *) 'Error: (read_modelParam) if file looks OK, try running dos2unix on it. Exiting...'
            stop
          end if
          call setCoeffValue_modelParam(P,n,l,m,v,min,max,if_fixed)
          !print *,'Values: ',l,m,v,min,max,if_fixed
       end do

    end do

    ! this check is not necessary, but helps to debug parametrization scripts
    !write(6,'(a50,i3)') 'Number of variable parameters in script: ',count(.not.P%c%frozen)
    write(6,*)

    close(ioPrm)

    P%zeroValued = .FALSE.
    P%smoothed = .FALSE.

    return

  end subroutine read_modelParam

end module model_operators
