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
! Also includes condictivity mappings on the grid:
!    CellToNode, NodeToCell, ModelParamToCell,
!    rhoC, CellToEdge, EdgeToCell, QtoModelParam

use file_units
use math_constants
use utilities
use emfield
#ifdef MPI
  use MPI_declaration
#endif

implicit none

  ! supported model parameter types (conductivity only)
   character(len=80), parameter		:: LOGE = 'LOGE'
   character(len=80), parameter		:: LINEAR = 'LINEAR'

type ::  modelParam_t
   !  modelParam is derived data type used to store parameters that
   !   determine conductivity/resistivity model
   !   This specific implementation is based on blocks
   !    allows for linear/log conductivity/resistivity
   private
   integer   				:: Ny = 0
   integer   				:: NzEarth = 0
   type(grid_t), pointer  		:: grid
   real (kind=prec), allocatable, dimension(:,:) :: v
   real (kind=prec)		:: AirCond = SIGMA_AIR
   logical           			:: allocated = .false.
     !  this logical is set true by zero_modelParam ONLY
     logical                    :: zeroValued = .false.
   !   necessary to avoid memory leaks; only true for function outputs
   logical						:: temporary = .false.
   !   parameter types supported:
   !   LINEAR = linear conductivity of each grid Earth cell is specified
   !   LOGE = natural log of conductivity in Earth cells specified
   character*80      			:: paramType = ''
end type modelParam_t

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

!  I/O interfaces

interface write_modelParam
   MODULE PROCEDURE write_modelParam_mackie
end interface

interface read_modelParam
   MODULE PROCEDURE read_modelParam_mackie
end interface

interface writeVec_modelParam
   MODULE PROCEDURE writeVec_modelParam_binary
end interface

interface readVec_modelParam
   MODULE PROCEDURE readVec_modelParam_binary
end interface

! definitions for CmSqrt: must be consistent with the include file below

#include "modelCov/Diffusion.hd"

Contains

! *****************************************************************************
!  conductivity mappings and adjoints for 2D MT modeling and inversion code
!  routines that are public
#include "ModelMap.inc"

!  The included file must contain subroutines create_CmSqrt, deall_CmSqrt, multBy...
#include "modelCov/Diffusion.inc"

!  I/O choices
#include "modelParamIO/Binary.inc"
#include "modelParamIO/Mackie.inc"

!  MPI model parameter, if needed
#ifdef MPI
#include "ModelParam_MPI.inc"
#endif
!************************************************************************
   !  allocateEarthCond allocates and initializes arrays for
   !   Earth-cell conductivity structure;
   !   Pass grid of type grid_t to set array sizes
   subroutine create_modelParam(grid,paramtype,cond,value,airCond)

     implicit none
     type (grid_t), intent(in), target   	:: grid
     character(*), intent(in)   		    :: paramtype
     type (modelParam_t), intent(inout)   	:: cond
     real (kind=prec), intent(in), optional :: value(:,:)
     real (kind=prec), intent(in), optional :: airCond
     !  local variables
     integer ::       Nz,Ny,Nza,NzEarth,istat

	 if (cond%allocated) then
        call warning('Model parameter already allocated in create_modelParam')
        return
	 end if

     Nz = grid%Nz
     Nza = grid%Nza
     NzEarth = Nz-Nza
     Ny = grid%Ny
     cond%Ny = Ny
     cond%NzEarth = NzEarth
     cond%grid => grid
     cond%paramType = paramtype

     allocate(cond%v(Ny,NzEarth), STAT=istat)
     
     if (present(value)) then
        if((size(value,1)==Ny).and.(size(value,2)==NzEarth)) then
           cond%v = value
        else
           call errStop('Wrong conductivity array size in create_modelParam')
        end if
     else
        cond%v = R_ZERO
     end if
     cond%allocated = .true.
     cond%zeroValued = .false.

     if (present(airCond)) then
        cond%AirCond = airCond
     end if

   end subroutine create_modelParam

  !************************************************************
   subroutine deall_modelParam(cond)
     implicit none
     type (modelParam_t)   :: cond
     ! local
     integer		:: istat

     if(cond%allocated) then
        deallocate(cond%v, STAT=istat)
        !nullify(cond%v)
        nullify(cond%grid)
        cond%allocated = .false.
        cond%zeroValued = .false.
        cond%paramType = ''
     endif

   end subroutine deall_modelParam

!**********************************************************************
   function dotProd_modelParam(m1,m2) result(r)

   !   dot product of two model space parameter objects
   !    here implemented for 2D blocks
     real(kind=prec)		:: r
     type(modelParam_t), intent(in)	:: m1,m2

     ! local variables
     integer	:: j,k

     if((m1%Ny .ne. m2%Ny).or. (m1%NzEarth .ne. m2%NzEarth)) then
        write(6,*) 'Ny ',m1%Ny,m2%Ny,m1%NzEarth,m2%NzEarth
        call errStop('size of m1, m2 incompatible in dotProd_modelParam')
     endif

     r = R_ZERO

     ! if m1 and m2 are zero valued, no need to compute r
     if (m1%zeroValued .and. m2%zeroValued) then
        return
     endif

     do k = 1,m1%NzEarth
        do j = 1,m1%Ny
           r = r + m1%v(j,k)*m2%v(j,k)
        enddo
     enddo

   end function dotProd_modelParam

!**********************************************************************
   function maxNorm_modelParam(m) result(r)

   !   Max norm = max(|m_ij|)
   !    here implemented for 2D blocks
     real(kind=prec)		:: r
     type(modelParam_t), intent(in)	:: m

     ! local variables
     integer	:: j,k

     if(.not.(m%allocated)) then
        call errStop('m not allocated in maxNorm_modelParam')
     endif

     r = max(abs(minval(m%v)),abs(maxval(m%v)))

   end function maxNorm_modelParam

!**********************************************************************
   subroutine zero_modelParam(m)

     !  zeros a model space object

     type(modelParam_t), intent(inout)	:: m

     m%v = R_ZERO
     m%zeroValued = .true.

   end subroutine zero_modelParam

!**********************************************************************

   logical function iszero_modelParam(m) result (zeroValued)

     type(modelParam_t), intent(in)      :: m

     zeroValued = m%zeroValued

   end function iszero_modelParam

!**********************************************************************
   subroutine random_modelParam(m,eps)

     !  zeros a model space object

     type(modelParam_t), intent(inout)  :: m
     real(kind=prec), intent(in), optional :: eps

     call random_number(m%v)
     if (present(eps)) then
        m%v = m%v * eps
     else
        m%v = m%v * 0.05
     endif
     m%zeroValued = .false.

   end subroutine random_modelParam

!**********************************************************************
   subroutine linComb_modelParam(a1,m1,a2,m2,m)

     !  forms the linear combination of model parameters
     !    m = a1*m1 + a2*m2
     !  where a1 and a2 are real constants and m1 and m2
     !   are model parameters
     !   output m may overwrite m1 or m2

     real(kind=prec),intent(in)		:: a1,a2
     type(modelParam_t), intent(in)	:: m1,m2
     type(modelParam_t), intent(inout)	:: m

     if((m1%Ny .ne. m2%Ny).or. (m1%NzEarth .ne. m2%NzEarth)) then
        call errStop('size of m1, m2 incompatable in linComb_modelParam')
     endif
     if(m1%paramType .ne. m2%paramType) then
        call errStop('paramType incompatable in linComb_modelParam')
     endif

     ! make sure m is allocated, and of same size; if not same
     !   size, deallocate and allocate as correct size; otherwise
     !   do nothing (use m as input ... allowing m to overwrite m1)
     if(m%allocated) then
        if((m%Ny .ne. m1%Ny).or. (m%NzEarth .ne. m1%NzEarth)) then
           call deall_modelParam(m)
           call create_modelParam(m1%grid,m1%paramtype,m)
        endif
     else
        call create_modelParam(m1%grid,m1%paramtype,m)
     endif
     m%v = a1*m1%v + a2*m2%v
     ! m is zero is both m1 and m2 are (ignoring the constants)
     m%zeroValued = m1%zeroValued .and. m2%zeroValued
     !  we are implicitly assuming that if m1 and m2 are the same
     !   size other parameters are of the same type
     m%AirCond = m1%AirCond

   end subroutine linComb_modelParam

  ! **********************************************************************
    subroutine scMult_modelParam(a,mIn,mOut)
  !  computes mOut = a * mIn for modelParam object m and real scalar a

    real (kind=prec), intent(in)			:: a
    type(modelParam_t), intent(in)	        :: mIn
    type(modelParam_t), intent(inout)       :: mOut

    ! check to see that input m is allocated
    if(.not.mIn%allocated) then
       call errStop('input not allocated on call to scMult_modelParam')
    endif

    ! if mIn is zero, no need to compute mOut
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
   subroutine copy_modelParam(mOut,mIn)

     type(modelParam_t), intent(in)	:: mIn
     type(modelParam_t), intent(inout)	:: mOut

     ! if mOut is allocated, check to see if if is of same size as
     !   mIn; if not, deallocate and reallocate as correct size; otherwise
     !   use as input
     if(mOut%allocated) then
        if((mOut%Ny .ne. mIn%Ny).or. (mOut%NzEarth .ne. mIn%NzEarth)) then
           call deall_modelParam(mOut)
           call create_modelParam(mIn%grid,mIn%paramtype,mOut)
        endif
     else
        call deall_modelParam(mOut)
        call create_modelParam(mIn%grid,mIn%paramtype,mOut)
     endif
     mOut%v = mIn%v
     mOut%AirCond = mIn%AirCond
     mOut%zeroValued = mIn%zeroValued

     if(mIn%temporary) then
     	call deall_modelParam(mIn)
     endif

   end subroutine copy_modelParam

   !************************************************************************
   !  count_modelParam counts the number of variable model parameters
   function count_modelParam(cond) result (N)

     implicit none
     type (modelParam_t), intent(in)   	  :: cond
     integer                              :: N

     if (.not.cond%allocated) then
        call errStop('Model parameter not allocated in count_modelParam')
     end if

     N = cond%Ny * cond%NzEarth

   end function count_modelParam

   !************************************************************************
   !  getValue_modelParam extracts information for a modelParam_t variable;
   !  note that the output array must already be of correct size on input.
   !  Therefore, grid is assumed to be known before calling this subroutine.
   subroutine getValue_modelParam(cond,value,paramtype,airCond)

     implicit none
     type (modelParam_t), intent(in)   	        :: cond
     real (kind=prec), intent(inout)    :: value(:,:)
     character*80, intent(out), optional   		:: paramtype
     real (kind=prec), intent(out), optional :: airCond
     !  local variables
     integer ::       Nz,Ny,Nza,NzEarth

     if (.not.cond%allocated) then
        call errStop('Model parameter not allocated in getValue_modelParam')
     end if

     Ny = cond%Ny
     NzEarth = cond%NzEarth

     if((size(value,1)==Ny).and.(size(value,2)==NzEarth)) then
        value = cond%v
     else
        call errStop('Wrong conductivity array size in getValue_modelParam')
     end if

     if (present(paramtype)) then
        paramtype = cond%paramType
     end if

     if (present(airCond)) then
        airCond = cond%AirCond
     end if

   end subroutine getValue_modelParam

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

end module ModelSpace
