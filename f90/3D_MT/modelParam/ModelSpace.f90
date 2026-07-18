module ModelSpace

! Defines the modelParam_t and modelCov_t abstract data types.
! Both have private attributes, to force coding of other modules
! to be completely independent of the specific parameterization
! instance/implementation. Contains all methods required for creation,
! algebra, and inner products, including covariances.
!
! Must include the following public operations:
!    create_modelParam, deall_modelParam, zero_modelParam, copy_modelParam,
!    linComb_modelParam, scMult_modelParam, dotProd_modelParam
!
! and the following interfaces:
!    assignment (=), operator (.dot.), operator (*),
!    create_CmSqrt, deall_CmSqrt, multBy_CmSqrt,
!    read_modelParam, write_modelParam.
!
! Also includes conductivity mappings on the grid:
!    ModelParamToCell, ModelParamToEdge, EdgeToModelParam,
!    QtoModelParam, sigC

#ifdef HDF5
  use hdf5
#endif
  use gridcalc
  use file_units
  use math_constants
  use utilities
  use sg_scalar
  use sg_vector
  use sg_sparse_vector
#ifdef SP_FORM
  use spOpTools
#endif
#ifdef MPI
  use Declaration_MPI
#endif

  implicit none


  ! supported model parameter types (conductivity only)
   character(len=80), parameter		:: LOGE = 'LOGE'
   character(len=80), parameter     :: LOG_10 = 'LOG10'
   character(len=80), parameter		:: LINEAR = 'LINEAR'

  type :: modelParam_t
     !  Parameters which define conductivity distribution.
     !   for the present approach, with each cell of the computational grid
     !   allowed a separate (constant) conductivity) a separate type here
     !   is hardly needed ...
     !  The idea is that the conductivity parameter should be treated as an
     !   "abstract data type", defined along with a set of routines for
     !    mapping to the internal representation of conductivity used directly
     !    by the forward solver.
     !  HERE we are starting to implement this idea ... cond_param will be
     !   public, but all attributes will be private.  Only routines in this
     !   module can use the actual attributes that define the specific realization
     !   of the model parameterization
     private
     integer			:: Nx,Ny,NzEarth
     type(rscalar)		:: cellCond
     type(grid_t),pointer    :: grid
     real (kind=prec)   :: AirCond
     logical			:: allocated = .false.
     !  this logical is set true by zero_modelParam ONLY
     logical            :: zeroValued = .false.
     !  necessary to avoid memory leaks; only true for function outputs
     logical			:: temporary = .false.
     !  another logical that gets set to true every time the model
     !  parameter is modified; used by the ForwardSolver for updateCond
     logical			:: updated = .false.
     !  supported paramType at present: LINEAR, LOG10 and LOGE
     character (len=80)	:: paramType = ''
  end type modelParam_t

  ! Used for HDF5 attributes
  type :: custom_att
    ! Define different attribute types
    character(len=100)          :: att_string
    integer                     :: att_int
    real(8)                     :: att_real

  end type custom_att

  character(len=80), save           :: userParamType = 'LOGE'


interface assignment (=)
   MODULE PROCEDURE copy_modelParam
end interface

interface zero
   MODULE PROCEDURE zero_modelParam
end interface

interface iszero
   MODULE PROCEDURE iszero_modelParam
end interface

interface scMult ! operator (*)
   MODULE PROCEDURE scMult_modelParam
end interface

interface scMultAdd
   MODULE PROCEDURE scMultAdd_modelParam
end interface

interface dotProd
   MODULE PROCEDURE dotProd_modelParam
end interface

interface linComb
   MODULE PROCEDURE linComb_modelParam
end interface

interface deall
   MODULE PROCEDURE deall_modelParam
end interface

interface countModelParam
   MODULE PROCEDURE count_modelParam
end interface

! IO Module Interfaces - These are implemented in modelParamIO (and it's submodules)
interface
    module subroutine setup_modelParamIO(input_ftype, output_ftype)
        implicit none
        character(len=*), intent(in) :: input_ftype
        character(len=*), intent(in) :: output_ftype
    end subroutine setup_modelParamIO

	module subroutine write_modelParam(m, cfile, comment)
        implicit none
        type(modelParam_t), intent(in):: m
        character(len=*), intent(in) :: cfile
        character(len=*), intent(in), optional :: comment
    end subroutine write_modelParam

    module subroutine read_modelParam(grid, airLayers, m, cfile)
        use GridDef, only : grid_t, airLayers_t
        implicit none
        type(grid_t), target, intent(inout)  :: grid
	    type(airLayers_t), intent(inout)	   :: airLayers
        type(modelParam_t), intent(out)	   :: m
        character(*), intent(in)             :: cfile
    end subroutine read_modelParam
end interface

interface
    module subroutine readVec_modelParam(grid,nSigma,sigma,header,cfile)
      implicit none
      type(grid_t), target, intent(inout)   :: grid
      integer, intent(in)		:: nSigma
      type(modelParam_t), intent(inout) 	:: sigma(nSigma)
      character(*), intent(out)		:: header
      character(*), intent(in)		:: cfile
    end subroutine readVec_modelParam

    module subroutine writeVec_modelParam(nSigma,sigma,header,cfile)
        implicit none
        integer, intent(in)		:: nSigma
        character(*), intent(in)		:: header, cfile
        type(modelParam_t), intent(in)	:: sigma(nSigma)
    end subroutine writeVec_modelParam
end interface

! definitions for CmSqrt: must be consistent with the include file below

#include "modelCov/RecursiveAR.hd"
#ifdef SP_FORM
#include "modelCov/BiHelmholtz.hd"
#endif

  ! CovType: 1 = AR (default), 2 = Bi-Helmholtz
  integer, save :: CovType = 1

Contains

! *****************************************************************************
!  routines which define mappings between the "natural" representation
!  of conductivity/resistivity on the model grid and the formal
!  model parameter structure
#include "ModelMap.inc"

!  Backend implementations
#include "modelCov/RecursiveAR.inc"
#ifdef SP_FORM
#include "modelCov/BiHelmholtz.inc"
#endif

!  dummy functions for covtype (always available, not dependent on spOpTools)
subroutine set_CovType(covTypeIn)
  integer, intent(in) :: covTypeIn
  CovType = covTypeIn
end subroutine set_CovType

function get_CovType() result(covTypeOut)
  integer :: covTypeOut
  covTypeOut = CovType
end function get_CovType

!  real functions: route to CovType=1 (AR) or CovType=2 (Bi-Helmholtz) backends

  function multBy_Cm(mhat) result(dm)
    type (modelParam_t), intent(in) :: mhat
    type (modelParam_t)             :: dm
    select case (CovType)
    case (1)
       dm = multBy_Cm_AR(mhat)
#ifdef SP_FORM
    case (2)
       dm = multBy_Cm_BiHelm(mhat)
#endif
    case default
       call errStop('Unknown CovType in multBy_Cm')
    end select
  end function multBy_Cm

  function multBy_CmSqrt(mhat) result(dm)
    type (modelParam_t), intent(in) :: mhat
    type (modelParam_t)             :: dm
    select case (CovType)
    case (1)
       dm = multBy_CmSqrt_AR(mhat)
#ifdef SP_FORM
    case (2)
       dm = multBy_CmSqrt_BiHelm(mhat)
#endif
    case default
       call errStop('Unknown CovType in multBy_CmSqrt')
    end select
  end function multBy_CmSqrt

  function multBy_CmSqrtInv(dm) result(mhat)
    type (modelParam_t), intent(in) :: dm
    type (modelParam_t)             :: mhat
    select case (CovType)
    case (1)
       mhat = multBy_CmSqrtInv_AR(dm)
#ifdef SP_FORM
    case (2)
       mhat = multBy_CmSqrtInv_BiHelm(dm)
#endif
    case default
       call errStop('Unknown CovType in multBy_CmSqrtInv')
    end select
  end function multBy_CmSqrtInv

  subroutine create_CmSqrt(m, cfile)
    type (modelParam_t), intent(in) :: m
    character(*), intent(in), optional :: cfile
    select case (CovType)
    case (1)
       call create_CmSqrt_AR(m, cfile)
#ifdef SP_FORM
    case (2)
       call create_CmSqrt_BiHelm(m, cfile)
#endif
    case default
       call errStop('Unknown CovType in create_CmSqrt')
    end select
  end subroutine create_CmSqrt

  subroutine deall_CmSqrt()
    select case (CovType)
    case (1)
       call deall_CmSqrt_AR()
#ifdef SP_FORM
    case (2)
       call deall_CmSqrt_BiHelm()
#endif
    case default
       call errStop('Unknown CovType in deall_CmSqrt')
    end select
  end subroutine deall_CmSqrt

!  MPI model parameter, if needed
#ifdef MPI
#include "ModelParam_MPI.inc"
#endif

!**********************************************************************
!
   !  create_modelParam allocates and initializes arrays for
   !   conductivity parameter structure
   !   Pass grid of type grid_t to set array sizes
   ! optional arguments v and vAir are assumed to be consistent
   ! with paramType
   subroutine create_modelParam(grid,paramtype,m,v,vAir)

     implicit none
     type (grid_t), intent(in), target	:: grid
     character(80), intent(in)			    :: paramtype
     type (modelParam_t), intent(inout)		:: m
     type (rscalar), intent(in), optional   :: v
     real (kind=prec), intent(in), optional :: vAir

     if(m%allocated) then
        call deall_modelParam(m)
     endif

     call create_rscalar(grid,m%cellCond,CELL_EARTH)
     m%Nx = grid%Nx
     m%Ny = grid%Ny
     m%NzEarth = grid%NzEarth
     m%grid => grid
     m%paramType = paramtype
     m%allocated = .true.

     if(present(v)) then
        m%cellCond = v
     endif

     if(present(vAir)) then
        m%AirCond = vAir
     else
		if(paramtype .eq. LOGE) then
		   m%AirCond = log(SIGMA_AIR)
		elseif(paramtype .eq. LOG_10) then
		   m%AirCond = log10(SIGMA_AIR)
		else
		   m%AirCond = SIGMA_AIR
		endif
     endif

     m%updated = .true.
     m%zeroValued = .false.

   end subroutine create_modelParam

  !************************************************************
   subroutine deall_modelParam(m)
     implicit none
     type (modelParam_t)   :: m

     if(m%allocated) then
        call deall_rscalar(m%cellCond)
        nullify(m%grid)
        m%allocated = .false.
        m%zeroValued = .false.
        m%paramType = ''
     endif

   end subroutine deall_modelParam

	!**********************************************************************
  subroutine getType_modelParam(m,paramType)
      type(modelParam_t), intent(in)    :: m
      character(*), intent(out)		      :: paramType
	  paramType=trim(m%paramType)
 end subroutine getType_modelParam

   !**********************************************************************
   ! Converts the input model parameter structure to paramType, by
   ! comparing paramType with m%paramType and performing the necessary
   ! computations if the two strings differ; assumes that m is allocated.
   subroutine setType_modelParam(m,paramType)

     type(modelParam_t), intent(inout)    :: m
     character(*), intent(in)		      :: paramType

     if(.not.(m%allocated)) then
        call errstop('modelParam must be allocated before calling setType_modelParam')
     endif

	 if(trim(paramType) .eq. trim(m%paramType)) then
	    ! we are done
	 elseif(m%paramType == LINEAR) then
	    ! convert to log
	    if(paramType == LOGE) then
	        m%cellCond%v = log(m%cellCond%v)
	        m%AirCond=log(m%AirCond)
	    else if(paramType == LOG_10) then
            m%cellCond%v = log10(m%cellCond%v)
            m%AirCond=log10(m%AirCond)
	    endif
	 elseif(paramType == LINEAR) then
	    ! convert from log to linear
	    if(m%paramType == LOGE) then
	        m%cellCond%v = exp(m%cellCond%v)
	        m%AirCond=exp(m%AirCond)
	    else if(m%paramType == LOG_10) then
            m%cellCond%v = exp(m%cellCond%v * log(10.))
            m%AirCond=exp(m%AirCond * log(10.))
        endif
     elseif ((m%paramType == LOGE) .and. (paramType == LOG_10)) then
        ! convert from natural log to log10
        m%cellCond%v = m%cellCond%v / log(10.)
        m%AirCond=m%AirCond / log(10.)
     elseif ((m%paramType == LOG_10) .and. (paramType == LOGE)) then
        ! convert from log10 to natural log
        m%cellCond%v = m%cellCond%v * log(10.)
        m%AirCond=m%AirCond * log(10.)
	 else
        call errstop('unknown paramType in setType_modelParam')
     endif

     m%paramType = paramType
     return

   end subroutine setType_modelParam

!**********************************************************************
   function dotProd_modelParam(m1,m2) result(r)

   !   dot product of two model space parameter objects
   !    here implemented for 3D numerical grid
     real(kind=prec)          :: r
     type(modelParam_t), intent(in)       :: m1,m2

     ! local variables
     integer    :: j,k

     if((m1%Ny .ne. m2%Ny).or. (m1%Nx .ne. m2%Nx).or. &
		 (m1%NzEarth .ne. m2%NzEarth)) then
        write(6,*) 'Nx, Ny, Nz ',m1%Nx, m2%Nx,    &
			m1%Ny,m2%Ny,m1%NzEarth,m2%NzEarth
        call errStop('size of m1, m2 incompatable in dotProd_modelParam')
     endif

     ! if one of the model parameters is zero valued, no need to compute r
     r = R_ZERO
     if(.not. m1%zeroValued .and. .not. m2%zeroValued) then
        r = dotProd_rscalar_f(m1%cellCond, m2%cellCond)
     endif

   end function dotProd_modelParam

!**********************************************************************

   subroutine zero_modelParam(m)

     !  zeros a model space object

     type(modelParam_t), intent(inout) 		:: m

     call zero_rscalar(m%cellCond)

     m%updated = .true.
     m%zeroValued = .true.

   end subroutine zero_modelParam


!**********************************************************************

   logical function iszero_modelParam(m) result (zeroValued)

     type(modelParam_t), intent(in)      :: m

     zeroValued = m%zeroValued

   end function iszero_modelParam

!**********************************************************************
   subroutine random_modelParam(m,eps)

     !  generated a random model parameter perturbation [0,1) in log
     !  space; otherwise, exp of that in linear space

     type(modelParam_t), intent(inout)  :: m
     real(kind=prec), intent(in), optional :: eps

     if (.not. m%allocated) then
       call errStop('model parameter not allocated in random_modelParam')
     end if

     if (present(eps)) then
        call random_rscalar(m%cellCond,eps)
     else
        call random_rscalar(m%cellCond)
     endif

     if(m%paramType == LINEAR) then
        m%cellCond%v = exp(m%cellCond%v)
     endif

     m%updated = .true.
     m%zeroValued = .false.

   end subroutine random_modelParam

!**********************************************************************

   subroutine linComb_modelParam(a1,m1,a2,m2,m)

     !  forms the linear combination of model parameters
     !    m = a1*m1 + a2*m2
     !  where a1 and a2 are real constants and m1 and m2
     !   are model parameters
     !   output m may overwrite m1 or m2

     real(kind=prec),intent(in)       :: a1,a2
     type(modelParam_t), intent(in)		:: m1,m2
     type(modelParam_t), intent(inout)		:: m

     if((m1%Ny .ne. m2%Ny).or. (m1%Nx .ne. m2%Nx) .or. &
		(m1%NzEarth .ne. m2%NzEarth)) then
        call errStop('size of m1, m2 incompatable in linComb_modelParam')
     endif
     if(m1%paramType .ne. m2%paramType) then
        call errStop('paramType incompatable in linComb_modelParam')
     endif

     ! make sure m is allocated, and of same size; if not same
     !   size, deallocate and allocate as correct size; otherwise
     !   do nothing (use m as input ... allowing m to overwrite m1)
     if(m%allocated) then
        if((m1%Ny .ne. m2%Ny).or. (m1%Nx .ne. m2%Nx) .or. &
		(m1%NzEarth .ne. m2%NzEarth)) then
           call deall_modelParam(m)
           call create_modelParam(m1%grid,m1%paramtype,m)
        endif
     else
        call create_modelParam(m1%grid,m1%paramtype,m)
     endif
     m%cellCond%v = a1*m1%cellCond%v + a2*m2%cellCond%v
     !  we are implicitly assuming that if m1 and m2 are the same
     !   size other parameters are of the same type
     m%AirCond = m1%AirCond
     m%zeroValued = m1%zeroValued .and. m2%zeroValued
     m%updated = .true.

   end subroutine linComb_modelParam

  ! **********************************************************************
    subroutine scMult_modelParam(a,mIn,mOut)
  !  computes mOut = a * mIn for modelParam object m and real scalar a

    real (kind=prec), intent(in)				:: a
    type(modelParam_t), intent(in)	            :: mIn
    type(modelParam_t), intent(inout)           :: mOut

    ! check to see that input m is allocated
    if(.not.mIn%allocated) then
       call errStop('input not allocated on call to scMult_modelParam')
    endif

    ! if mIn is zero valued, no need to compute mOut
    if (mIn%zeroValued) then
        mOut = mIn
    else
    	call linComb_modelParam(R_ZERO,mIn,a,mIn,mOut)
    endif

  end subroutine scMult_modelParam

  ! **********************************************************************
    subroutine scMultAdd_modelParam(a,mIn,mOut)
  !  computes mOut = a * mIn + mOut for modelParam object m and real scalar a

    real (kind=prec), intent(in)                :: a
    type(modelParam_t), intent(in)              :: mIn
    type(modelParam_t), intent(inout)           :: mOut

    ! check to see that input m is allocated
    if(.not.(mIn%allocated .and. mOut%allocated)) then
       call errStop('input not allocated on call to scMultAdd_modelParam')
    endif

    ! if mIn is zero valued, no need to compute mOut
    if (.not. mIn%zeroValued) then
        call linComb_modelParam(a,mIn,ONE,mOut,mOut)
    endif

  end subroutine scMultAdd_modelParam

   !**********************************************************************
   ! Sets cell conductivities in a model parameter object m
   !
   ! the values of v and vAir are determined by paramType;
   ! if different from that of the model parameter, returns
   ! an error. To avoid this error, first convert m to the
   ! required paramType using setType_modelParam.
   subroutine setValue_modelParam(m,paramType,v,vAir)

     type(modelParam_t), intent(inout)    :: m
     character(80), intent(in)		      :: paramType
     type(rscalar), intent(in)		      :: v
     real(kind=prec), intent(in), optional :: vAir

     if(.not.(m%allocated)) then
        call errstop('output modelParam must be allocated before calling setValue_modelParam')
     endif

     !  error checking
     if((m%Ny .ne. v%Ny).or. (m%Nx .ne. v%Nx) .or. (m%NzEarth .ne. v%Nz)) then
        call errstop('modelParam/rscalar dimensions disagree in setValue_modelParam')
     else if(paramType .ne. m%paramType) then
        call errstop('paramTypes not consistent in setValue_modelParam')
     endif

	 ! set values
	 m%cellCond = v
	 if(present(vAir)) then
	    m%AirCond=vAir
	 endif
     m%updated = .true.
     m%zeroValued = .false.

   end subroutine setValue_modelParam

   !**********************************************************************
   ! Gets cell conductivities from a model parameter object m
   !
   ! Extracts the values of v and vAir;
   ! values that are extracted are converted to paramType.
   ! Different from ModelParamToCell in that the value
   ! that gets extracted is exactly what is stored in the modelParam:
   ! it does not contain any air layers. This is needed for BC_x0_WS.
   subroutine getValue_modelParam(m,paramType,v,vAir)

     type(modelParam_t), intent(in)       :: m
     character(80), intent(inout)		  :: paramType
     type(rscalar), intent(out)		      :: v
     real(kind=prec), intent(out), optional :: vAir
     ! local variable
     type(modelParam_t)                   :: mTemp

     if(.not.(m%allocated)) then
        call errstop('input modelParam must be allocated before calling getValue_modelParam')
     endif

     if(trim(paramType) .eq. '') then
     	paramType = m%paramType
     endif

     if (v%allocated) then
        call deall_rscalar(v)
     endif

     ! create a temporary copy
     mTemp = m

     ! convert model to the required type
     call setType_modelParam(mTemp,paramType)

	 ! set values
	 v = mTemp%cellCond
	 if(present(vAir)) then
	    vAir = mTemp%AirCond
	 endif

	 ! deallocate temporary model parameter
	 call deall_modelParam(mTemp)

   end subroutine getValue_modelParam

   !**********************************************************************
   subroutine copy_modelParam(mOut,mIn)

     type(modelParam_t), intent(in)       :: mIn
     type(modelParam_t), intent(inout)    :: mOut

     ! if mOut is allocated, check to see if if is of same size as
     !   mIn; if not, deallocate and reallocate as correct size; otherwise
     !   use as input
     if(mOut%allocated) then
        if((mOut%Ny .ne. mIn%Ny).or. (mOut%Nx .ne. mIn%Nx) .or. &
		(mOut%NzEarth .ne. mIn%NzEarth)) then
           call deall_modelParam(mOut)
           call create_modelParam(mIn%grid,mIn%paramtype,mOut)
        endif
     else
        call create_modelParam(mIn%grid,mIn%paramtype,mOut)
     endif
     mOut%cellCond%v = mIn%cellCond%v
     mOut%AirCond = mIn%AirCond
     mOut%zeroValued = mIn%zeroValued

     ! no need to set updated to TRUE - this just makes a copy.
     ! clean up the memory if this is used with an "=" sign.
     if(mIn%temporary) then
     	call deall_modelParam(mIn)
     endif

   end subroutine copy_modelParam

   !************************************************************************
   !  count_modelParam counts the number of variable model parameters
   function count_modelParam(cond) result (N)

     implicit none
     type (modelParam_t), intent(in)      :: cond
     integer                              :: N

     if (.not.cond%allocated) then
        call errStop('Model parameter not allocated in count_modelParam')
     end if

     N = cond%Nx * cond%Ny * cond%NzEarth

   end function count_modelParam

  !**********************************************************************
  !  extracts the grid from a modelParam object; this is only needed since
  !  the attributes are private (in F2003, can declare everything private
  !  while grid and allocated attributes could be public)

  subroutine getGrid_modelParam(grid,mIn)

    type (grid_t), intent(out)          :: grid
    type (modelParam_t), intent(in)     :: mIn

    if (.not. mIn%allocated) then
       call warning('model vector not allocated on call to getGrid_modelParam')
    endif

    grid = mIn%grid

  end subroutine getGrid_modelParam

  !**********************************************************************
  !  sets modelParam%updated to FALSE; this is only needed since
  !  the attributes are private (in F2003, can declare everything private
  !  while grid and allocated attributes could be public)
  !  used when the model parameter is no longer considered "new",
  !  e.g. after the updateModelData routine in ForwardSolver.

  subroutine setValueUpdated_modelParam(m)

    type (modelParam_t), intent(inout)		:: m

    m%updated = .false.

  end subroutine setValueUpdated_modelParam

  !**********************************************************************
  !  checks whether a modelParam is updated; this is only needed since
  !  the attributes are private (in F2003, can declare everything private
  !  while grid and allocated attributes could be public)

  subroutine getValueUpdated_modelParam(m,updated)

    type (modelParam_t), intent(in)		:: m
    logical, intent(out)				:: updated

    updated = m%updated

  end subroutine getValueUpdated_modelParam

!**********************************************************************
!  Detects the filetype (CovType) from line 17 of a covariance file.
!  Reads the first non-blank line after the 16-line header and parses:
!    - 4 ints (Nx Ny NzEarth FileType): new format, FileType must be 1 or 2
!    - 3 ints (Nx Ny NzEarth):          legacy format, defaults to FileType=1
!    - anything else:                   error
!
!  Both read_CmSqrt_AR and read_CmSqrt_BiHelm call this routine.

  subroutine detect_cov_filetype(fid, Nx, Ny, NzEarth, fileType)

    integer, intent(in)                          :: fid
    integer, intent(out)                         :: Nx, Ny, NzEarth, fileType
    character(len=256)                           :: line
    integer                                      :: istat

    ! skip blank lines after header
    do
       read(fid, '(A)', iostat=istat) line
       if (istat /= 0) call errStop('Error reading covariance file after header')
       line = adjustl(line)
       if (line /= '') exit
    end do

    ! try 4 ints (new format: Nx Ny NzEarth FileType)
    read(line, *, iostat=istat) Nx, Ny, NzEarth, fileType
    if (istat == 0 .and. fileType >= 1 .and. fileType <= 2) return

    ! try 3 ints (legacy: Nx Ny NzEarth, default to CovType=1)
    read(line, *, iostat=istat) Nx, Ny, NzEarth
    if (istat == 0) then
       fileType = 1
       return
    end if

    call errStop('Invalid grid dimensions / CovType in covariance file')

  end subroutine detect_cov_filetype

end module ModelSpace
