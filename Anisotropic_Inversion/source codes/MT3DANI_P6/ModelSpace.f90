module ModelSpace ! 已改为电阻率形式的模型参数化方式

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

  use gridcalc
  use file_units
  use math_constants
  use utilities
  use sg_scalar
  use sg_vector
  use sg_sparse_vector
  
  use anisConstraint !!!
  
#ifdef MPI
  use MPI_declaration
#endif

  implicit none

  ! supported model parameter types (conductivity only)
   character(len=80), parameter		:: LOGE = 'LOGE'
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
     
     ! modified to array
     type(rscalar)		:: cellRho(6)
     
     type(grid_t),pointer    :: grid
     real (kind=prec)   :: AirRho
     logical			:: allocated = .false.
     !  this logical is set true by zero_modelParam ONLY
     logical            :: zeroValued = .false.
     !  necessary to avoid memory leaks; only true for function outputs
     logical			:: temporary = .false.
     !  another logical that gets set to true every time the model
     !  parameter is modified; used by the ForwardSolver for updateCond
     logical			:: updated = .false.
     !  supported paramType at present: LINEAR and LOGE
     character (len=80)	:: paramType = ''
     ! logical variables for controlling inversion parameters
     logical           :: fixed1 = .false.
     logical           :: fixed2 = .false.
     logical           :: fixed3 = .false.
     logical           :: fixed4 = .false.
     logical           :: fixed5 = .false.
     logical           :: fixed6 = .false.
     ! logical  variables for isotropy or special anisotropy
     logical           :: equal12 = .false.
     logical           :: equal13 = .false.
     logical           :: equal23 = .false.
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
   MODULE PROCEDURE write_modelParam_WS
end interface

interface read_modelParam
   MODULE PROCEDURE read_modelParam_WS
end interface

!interface writeVec_modelParam
!   MODULE PROCEDURE writeVec_modelParam_binary
!end interface
!
!interface readVec_modelParam
!   MODULE PROCEDURE readVec_modelParam_binary
!end interface

! definitions for CmSqrt: must be consistent with the include file below

include "modelCov/RecursiveAR.hd"
!#include "modelCov/Diffusion.hd"
Contains

! *****************************************************************************
!  routines which define mappings between the "natural" representation
!  of conductivity/resistivity on the model grid and the formal
!  model parameter structure
include "ModelMap.inc"

!  The included file must contain subroutines create_CmSqrt, deall_CmSqrt, multBy...
include "modelCov/RecursiveAR.inc"
!#include "modelCov/Diffusion.inc"
!  I/O choices
!#include "modelParamIO/Binary.inc"
!#include "modelParamIO/Mackie.inc"
include "modelParamIO/WS.inc"

!  MPI model parameter, if needed
#ifdef MPI
include "ModelParam_MPI.inc"
#endif

!**********************************************************************
!
   !  create_modelParam allocates and initializes arrays for
   !   conductivity parameter structure
   !   Pass grid of type grid_t to set array sizes
   ! optional arguments v and vAir are assumed to be consistent
   ! with paramType
   subroutine create_modelParam(grid,paramtype,m,v,flag,vAir)

     implicit none
     type (grid_t), intent(in), target	:: grid
     character(80), intent(in)			    :: paramtype
     type (modelParam_t), intent(inout)		:: m
     integer, intent(in),optional :: flag(9)
     type (rscalar), intent(in), optional   :: v(6)
     real (kind=prec), intent(in), optional :: vAir
     ! local
     integer :: ia

     if(m%allocated) then
        call deall_modelParam(m)
     endif
     
     do ia = 1,6
       call create_rscalar(grid,m%cellRho(ia),CELL_EARTH)
     enddo
     m%Nx = grid%Nx
     m%Ny = grid%Ny
     m%NzEarth = grid%NzEarth
     m%grid => grid  ! 在读取模型时创建了网格
     m%paramType = paramtype
     m%allocated = .true.

     if(present(v)) then
       do ia = 1,6
         m%cellRho(ia) = v(ia)
       enddo
     endif

     if(present(vAir)) then
       m%AirRho = vAir
     else
		   if(paramtype .eq. LOGE) then
		     m%AirRho = -log(SIGMA_AIR)
		   else
		     m%AirRho = ONE/SIGMA_AIR
		   endif
     endif

     ! these logical variables are used to control sepecific inversions
     if(present(flag)) then
       if(flag(1).eq.1) then
         m%fixed1 = .true.
       endif
       if(flag(2).eq.1) then
         m%fixed2 = .true.
       endif
       if(flag(3).eq.1) then
         m%fixed3 = .true.
       endif
       if(flag(4).eq.1) then
         m%fixed4 = .true.
       endif
       if(flag(5).eq.1) then
         m%fixed5 = .true.
       endif
       if(flag(6).eq.1) then
         m%fixed6 = .true.
       endif 
       if(flag(7).eq.1) then
         m%equal12 = .true.
       endif  
       if(flag(8).eq.1) then
         m%equal13 = .true.
       endif 
       if(flag(9).eq.1) then
         m%equal23 = .true.
       endif                                                     
     end if
     
     m%updated = .true.
     m%zeroValued = .false.

   end subroutine create_modelParam

  !************************************************************
   subroutine deall_modelParam(m)
     implicit none
     type (modelParam_t)   :: m
     ! local
     integer :: ia

     if(m%allocated) then
        do ia = 1,6
          call deall_rscalar(m%cellRho(ia))
        enddo
        nullify(m%grid)
        m%allocated = .false.
        m%zeroValued = .false.
        m%paramType = ''
        m%fixed1 = .false.
        m%fixed2 = .false.
        m%fixed3 = .false.
        m%fixed4 = .false.
        m%fixed5 = .false.
        m%fixed6 = .false.
        m%equal12 = .false.
        m%equal13 = .false.
        m%equal23 = .false.
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
     integer :: ia

     if(.not.(m%allocated)) then
        call errstop('modelParam must be allocated before calling setType_modelParam')
     endif

	   if(trim(paramType) .eq. trim(m%paramType)) then
	      ! we are done
	   else if((paramType == LOGE) .and. (m%paramType == LINEAR)) then
	      ! convert to log
	      do ia = 1,3
	        m%cellRho(ia)%v = log(m%cellRho(ia)%v)
	      enddo
	      m%AirRho=log(m%AirRho)
	   else if((paramType == LINEAR) .and. (m%paramType == LOGE)) then
	      ! convert to cell
	      do ia = 1,3
	        m%cellRho(ia)%v = exp(m%cellRho(ia)%v)
	      enddo
	      m%AirRho=exp(m%AirRho)
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
     real(kind=prec)          :: r,r1,r2,r3,r4,r5,r6
     type(modelParam_t), intent(in)       :: m1,m2
     integer :: ia
     logical :: equal12,equal13,equal23
     logical :: fixed1,fixed2,fixed3,fixed4,fixed5,fixed6

     ! local variables
     integer    :: j,k

     if((m1%Ny .ne. m2%Ny).or. (m1%Nx .ne. m2%Nx).or. &
		 (m1%NzEarth .ne. m2%NzEarth)) then
        write(6,*) 'Nx, Ny, Nz ',m1%Nx, m2%Nx,    &
			m1%Ny,m2%Ny,m1%NzEarth,m2%NzEarth
        call errStop('size of m1, m2 incompatable in dotProd_modelParam')
     endif

     fixed1 = m1%fixed1 .or. m2%fixed1
     fixed2 = m1%fixed2 .or. m2%fixed2
     fixed3 = m1%fixed3 .or. m2%fixed3
     fixed4 = m1%fixed4 .or. m2%fixed4
     fixed5 = m1%fixed5 .or. m2%fixed5
     fixed6 = m1%fixed6 .or. m2%fixed6
     equal12 = m1%equal12 .or. m2%equal12
     equal13 = m1%equal13 .or. m2%equal13
     equal23 = m1%equal23 .or. m2%equal23
     ! if one of the model parameters is zero valued, no need to compute r
     r = R_ZERO     
     if(.not. m1%zeroValued .and. .not. m2%zeroValued) then

        ! needed to be refined
        if((equal12.and.equal23) .or. (equal12.and.equal13)  &
          .or. (equal13.and.equal23)) then
          ! isotropy, 10/10/10
          r1 = dotProd_rscalar_f(m1%cellRho(1), m2%cellRho(1))
          r2 = R_ZERO
          r3 = R_ZERO
          r4 = R_ZERO
          r5 = R_ZERO
          r6 = R_ZERO
        elseif(equal12.and.(.not.equal23).and.(.not.equal13)) then
          ! Vertical anisotropy or dipping anisotropy, S angle fails, 10/10/100
          r1 = dotProd_rscalar_f(m1%cellRho(1), m2%cellRho(1))
          r2 = R_ZERO
          r3 = dotProd_rscalar_f(m1%cellRho(3), m2%cellRho(3))
          r4 = R_ZERO
          r5 = dotProd_rscalar_f(m1%cellRho(5), m2%cellRho(5))
          r6 = R_ZERO
        elseif(equal23.and.(.not.equal12).and.(.not.equal13)) then      
          ! biaxial anisotropy or horizontal anisotropy, 10/100/100
          r1 = dotProd_rscalar_f(m1%cellRho(1), m2%cellRho(1))
          r2 = dotProd_rscalar_f(m1%cellRho(2), m2%cellRho(2))
          r3 = R_ZERO
          r4 = dotProd_rscalar_f(m1%cellRho(4), m2%cellRho(4))
          r5 = R_ZERO
          r6 = R_ZERO
        elseif(equal13.and.(.not.equal12).and.(.not.equal23)) then      
          ! biaxial anisotropy or horizontal anisotropy or dipping anisotropy, 100/10/100
          r1 = dotProd_rscalar_f(m1%cellRho(1), m2%cellRho(1))
          r2 = dotProd_rscalar_f(m1%cellRho(2), m2%cellRho(2))
          r3 = R_ZERO
          r4 = dotProd_rscalar_f(m1%cellRho(4), m2%cellRho(4))
          r5 = dotProd_rscalar_f(m1%cellRho(5), m2%cellRho(5))
          r6 = R_ZERO
        else
          r1 = dotProd_rscalar_f(m1%cellRho(1), m2%cellRho(1))
          r2 = dotProd_rscalar_f(m1%cellRho(2), m2%cellRho(2))
          r3 = dotProd_rscalar_f(m1%cellRho(3), m2%cellRho(3))
          r4 = dotProd_rscalar_f(m1%cellRho(4), m2%cellRho(4))
          r5 = dotProd_rscalar_f(m1%cellRho(5), m2%cellRho(5))
          r6 = dotProd_rscalar_f(m1%cellRho(6), m2%cellRho(6))
        endif
        
        if(fixed1) then
          r1 = R_ZERO
        endif
        if(fixed2) then
          r2 = R_ZERO
        endif
        if(fixed3) then
          r3 = R_ZERO
        endif
        if(fixed4) then
          r4 = R_ZERO
        endif
        if(fixed5) then
          r5 = R_ZERO
        endif
        if(fixed6) then
          r6 = R_ZERO
        endif
          
        r = r1 + r2 + r3 + r4 + r5 + r6

     endif

   end function dotProd_modelParam


!**********************************************************************
   subroutine dotProd_mNorm(m1,m2,mNorm,TotalNorm)

   !   dot product of two model space parameter objects
   !    here implemented for 3D numerical grid
     real(kind=prec)          :: TotalNorm,mNorm(6)
     type(modelParam_t), intent(in)       :: m1,m2
     integer :: ia
     logical :: equal12,equal13,equal23
     logical :: fixed1,fixed2,fixed3,fixed4,fixed5,fixed6

     ! local variables
     integer    :: j,k

     if((m1%Ny .ne. m2%Ny).or. (m1%Nx .ne. m2%Nx).or. &
		 (m1%NzEarth .ne. m2%NzEarth)) then
        write(6,*) 'Nx, Ny, Nz ',m1%Nx, m2%Nx,    &
			m1%Ny,m2%Ny,m1%NzEarth,m2%NzEarth
        call errStop('size of m1, m2 incompatable in dotProd_modelParam')
     endif

     fixed1 = m1%fixed1 .or. m2%fixed1
     fixed2 = m1%fixed2 .or. m2%fixed2
     fixed3 = m1%fixed3 .or. m2%fixed3
     fixed4 = m1%fixed4 .or. m2%fixed4
     fixed5 = m1%fixed5 .or. m2%fixed5
     fixed6 = m1%fixed6 .or. m2%fixed6
     equal12 = m1%equal12 .or. m2%equal12
     equal13 = m1%equal13 .or. m2%equal13
     equal23 = m1%equal23 .or. m2%equal23
     ! if one of the model parameters is zero valued, no need to compute r 
     TotalNorm = R_ZERO   
     if(.not. m1%zeroValued .and. .not. m2%zeroValued) then

        ! needed to be refined
        if((equal12.and.equal23) .or. (equal12.and.equal13)  &
          .or. (equal13.and.equal23)) then
          ! isotropy, 10/10/10
          mNorm(1) = dotProd_rscalar_f(m1%cellRho(1), m2%cellRho(1))
          mNorm(2) = R_ZERO
          mNorm(3) = R_ZERO
          mNorm(4) = R_ZERO
          mNorm(5) = R_ZERO
          mNorm(6) = R_ZERO
        elseif(equal12.and.(.not.equal23).and.(.not.equal13)) then
          ! Vertical anisotropy or dipping anisotropy, S angle fails, 10/10/100
          mNorm(1) = dotProd_rscalar_f(m1%cellRho(1), m2%cellRho(1))
          mNorm(2) = R_ZERO
          mNorm(3) = dotProd_rscalar_f(m1%cellRho(3), m2%cellRho(3))
          mNorm(4) = R_ZERO
          mNorm(5) = dotProd_rscalar_f(m1%cellRho(5), m2%cellRho(5))
          mNorm(6) = R_ZERO
        elseif(equal23.and.(.not.equal12).and.(.not.equal13)) then      
          ! biaxial anisotropy or horizontal anisotropy, 10/100/100
          mNorm(1) = dotProd_rscalar_f(m1%cellRho(1), m2%cellRho(1))
          mNorm(2) = dotProd_rscalar_f(m1%cellRho(2), m2%cellRho(2))
          mNorm(3) = R_ZERO
          mNorm(4) = dotProd_rscalar_f(m1%cellRho(4), m2%cellRho(4))
          mNorm(5) = R_ZERO
          mNorm(6) = R_ZERO
        elseif(equal13.and.(.not.equal12).and.(.not.equal23)) then      
          ! biaxial anisotropy or horizontal anisotropy or dipping anisotropy, 100/10/100
          mNorm(1) = dotProd_rscalar_f(m1%cellRho(1), m2%cellRho(1))
          mNorm(2) = dotProd_rscalar_f(m1%cellRho(2), m2%cellRho(2))
          mNorm(3) = R_ZERO
          mNorm(4) = dotProd_rscalar_f(m1%cellRho(4), m2%cellRho(4))
          mNorm(5) = dotProd_rscalar_f(m1%cellRho(5), m2%cellRho(5))
          mNorm(6) = R_ZERO
        else
          mNorm(1) = dotProd_rscalar_f(m1%cellRho(1), m2%cellRho(1))
          mNorm(2) = dotProd_rscalar_f(m1%cellRho(2), m2%cellRho(2))
          mNorm(3) = dotProd_rscalar_f(m1%cellRho(3), m2%cellRho(3))
          mNorm(4) = dotProd_rscalar_f(m1%cellRho(4), m2%cellRho(4))
          mNorm(5) = dotProd_rscalar_f(m1%cellRho(5), m2%cellRho(5))
          mNorm(6) = dotProd_rscalar_f(m1%cellRho(6), m2%cellRho(6))
        endif
        
        if(fixed1) then
          mNorm(1) = R_ZERO
        endif
        if(fixed2) then
          mNorm(2) = R_ZERO
        endif
        if(fixed3) then
          mNorm(3) = R_ZERO
        endif
        if(fixed4) then
          mNorm(4) = R_ZERO
        endif
        if(fixed5) then
          mNorm(5) = R_ZERO
        endif
        if(fixed6) then
          mNorm(6) = R_ZERO
        endif
        
        do ia = 1,6
          TotalNorm = TotalNorm + mNorm(ia)
        enddo
     endif

   end subroutine dotProd_mNorm


!**********************************************************************

   subroutine zero_modelParam(m)

     !  zeros a model space object

     type(modelParam_t), intent(inout) 		:: m
     integer :: ia

     do ia = 1,6
       call zero_rscalar(m%cellRho(ia))
     enddo

     m%updated = .true.
     m%zeroValued = .true.

   end subroutine zero_modelParam


!**********************************************************************

   logical function iszero_modelParam(m) result (zeroValued)

     type(modelParam_t), intent(in)      :: m

     zeroValued = m%zeroValued

   end function iszero_modelParam

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
     integer :: ia

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
     do ia = 1,6
       m%cellRho(ia)%v = a1*m1%cellRho(ia)%v + a2*m2%cellRho(ia)%v
     enddo
     !  we are implicitly assuming that if m1 and m2 are the same
     !   size other parameters are of the same type
     m%AirRho = m1%AirRho
     m%zeroValued = m1%zeroValued .and. m2%zeroValued
     m%updated = .true.


     m%fixed1 = m1%fixed1.or.m2%fixed1
     m%fixed2 = m1%fixed2.or.m2%fixed2
     m%fixed3 = m1%fixed3.or.m2%fixed3
     m%fixed4 = m1%fixed4.or.m2%fixed4
     m%fixed5 = m1%fixed5.or.m2%fixed5
     m%fixed6 = m1%fixed6.or.m2%fixed6
     m%equal12 = m1%equal12.or.m2%equal12
     m%equal13 = m1%equal13.or.m2%equal13
     m%equal23 = m1%equal23.or.m2%equal23
        
   end subroutine linComb_modelParam


!**********************************************************************

   subroutine linComb_gradient(a1,m1,a2,lambda,m2,m) 

     !  forms the linear combination of model parameters
     !    m = a1*m1 + a2*m2
     !  where a1 and a2 are real constants and m1 and m2
     !   are model parameters
     !   output m may overwrite m1 or m2

     real(kind=prec),intent(in)       :: a1,a2
     real(kind=prec),intent(in)       :: lambda(8)
     type(modelParam_t), intent(in)		:: m1,m2
     type(modelParam_t), intent(inout)		:: m
     integer :: ia

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
     do ia = 1,6
       m%cellRho(ia)%v = a1*m1%cellRho(ia)%v + a2*lambda(ia)*m2%cellRho(ia)%v
     enddo
     !  we are implicitly assuming that if m1 and m2 are the same
     !   size other parameters are of the same type
     m%AirRho = m1%AirRho
     m%zeroValued = m1%zeroValued .and. m2%zeroValued
     m%updated = .true.


     m%fixed1 = m1%fixed1.or.m2%fixed1
     m%fixed2 = m1%fixed2.or.m2%fixed2
     m%fixed3 = m1%fixed3.or.m2%fixed3
     m%fixed4 = m1%fixed4.or.m2%fixed4
     m%fixed5 = m1%fixed5.or.m2%fixed5
     m%fixed6 = m1%fixed6.or.m2%fixed6
     m%equal12 = m1%equal12.or.m2%equal12
     m%equal13 = m1%equal13.or.m2%equal13
     m%equal23 = m1%equal23.or.m2%equal23
        
   end subroutine linComb_gradient
   
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
     type(rscalar), intent(in)		      :: v(6)
     real(kind=prec), intent(in), optional :: vAir
     integer :: ia

     if(.not.(m%allocated)) then
        call errstop('output modelParam must be allocated before calling setValue_modelParam')
     endif

     !  error checking     
     if((m%Ny .ne. v(1)%Ny).or. (m%Nx .ne. v(1)%Nx) .or. (m%NzEarth .ne. v(1)%Nz)) then
        call errstop('modelParam/rscalar dimensions disagree in setValue_modelParam')
     else if(paramType .ne. m%paramType) then
        call errstop('paramTypes not consistent in setValue_modelParam')
     endif

	   ! set values
	   do ia = 1,6
	     m%cellRho(ia) = v(ia)
	   enddo
	   if(present(vAir)) then
	     m%AirRho=vAir
	   endif
     m%updated = .true.
     m%zeroValued = .false.

   end subroutine setValue_modelParam

   !**********************************************************************
   ! Gets cell parameters from a model parameter object m
   !
   ! Extracts the values of v and vAir;
   ! values that are extracted are converted to paramType.
   ! Different from ModelParamToCell in that the value
   ! that gets extracted is exactly what is stored in the modelParam:
   ! it does not contain any air layers. This is needed for BC_x0_WS.
   subroutine getValue_modelParam(m,paramType,v,vAir)

     type(modelParam_t), intent(in)       :: m
     character(80), intent(inout)		  :: paramType
     type(rscalar), intent(out)		      :: v(6)
     real(kind=prec), intent(out), optional :: vAir
     ! local variable
     integer :: ia
     type(modelParam_t)                   :: mTemp

     if(.not.(m%allocated)) then
        call errstop('input modelParam must be allocated before calling getValue_modelParam')
     endif

     if(trim(paramType) .eq. '') then
     	 paramType = m%paramType
     endif

     do ia = 1,6
       if (v(ia)%allocated) then
         call deall_rscalar(v(ia))
       endif
     enddo

     ! create a temporary copy
     mTemp = m

     ! convert model to the required type
     call setType_modelParam(mTemp,paramType)

	   ! set values
	   do ia = 1,6
	     v(ia) = mTemp%cellRho(ia)
	   enddo
	   if(present(vAir)) then
	     vAir = mTemp%AirRho
	   endif

	   ! deallocate temporary model parameter
	   call deall_modelParam(mTemp)

   end subroutine getValue_modelParam


   !**********************************************************************
   ! Gets axial resistivities from a model parameter object m
   subroutine getResValue_modelParam(m,paramType,v,vAir)

     type(modelParam_t), intent(in)       :: m
     character(80), intent(inout)		  :: paramType
     type(rscalar), intent(out)		      :: v(3)
     real(kind=prec), intent(out), optional :: vAir
     ! local variable
     integer :: ia
     type(modelParam_t)                   :: mTemp

     if(.not.(m%allocated)) then
        call errstop('input modelParam must be allocated before calling getValue_modelParam')
     endif

     if(trim(paramType) .eq. '') then
     	 paramType = m%paramType
     endif

     do ia = 1,3
       if (v(ia)%allocated) then
         call deall_rscalar(v(ia))
       endif
     enddo

     ! create a temporary copy
     mTemp = m

     ! convert model to the required type
     call setType_modelParam(mTemp,paramType)

	   ! set values
	   do ia = 1,3
	     v(ia) = mTemp%cellRho(ia)
	   enddo
	   if(present(vAir)) then
	     vAir = mTemp%AirRho
	   endif

	   ! deallocate temporary model parameter
	   call deall_modelParam(mTemp)

   end subroutine getResValue_modelParam


   !**********************************************************************
   subroutine copy_modelParam(mOut,mIn)

     type(modelParam_t), intent(in)       :: mIn
     type(modelParam_t), intent(inout)    :: mOut
     
     integer :: ia

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
     
     do ia = 1,6
       mOut%cellRho(ia)%v = mIn%cellRho(ia)%v
     enddo
     mOut%AirRho = mIn%AirRho
     mOut%zeroValued = mIn%zeroValued

     mOut%fixed1 = mIn%fixed1
     mOut%fixed2 = mIn%fixed2
     mOut%fixed3 = mIn%fixed3
     mOut%fixed4 = mIn%fixed4
     mOut%fixed5 = mIn%fixed5
     mOut%fixed6 = mIn%fixed6
     mOut%equal12 = mIn%equal12
     mOut%equal13 = mIn%equal13
     mOut%equal23 = mIn%equal23
     
     ! no need to set updated to TRUE - this just makes a copy.
     ! clean up the memory if this is used with an "=" sign.
     if(mIn%temporary) then
     	call deall_modelParam(mIn)
     endif

   end subroutine copy_modelParam

   !************************************************************************
   !  count_modelParam counts the number of variable model parameters
   function count_modelParam(rho) result (N)

     implicit none
     type (modelParam_t), intent(in)      :: rho
     integer                              :: N

     if (.not.rho%allocated) then
        call errStop('Model parameter not allocated in count_modelParam')
     end if

     N = rho%Nx * rho%Ny * rho%NzEarth

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
  


  subroutine getValueFlag_modelParam(m,flag)

    type (modelParam_t), intent(in)		:: m
    integer, intent(out)				:: flag(9)

    if(m%fixed1) then
      flag(1) = 1
    else
      flag(1) = 0
    endif
    
    if(m%fixed2) then
      flag(2) = 1
    else
      flag(2) = 0
    endif
    
    if(m%fixed3) then
      flag(3) = 1
    else
      flag(3) = 0
    endif
    
    if(m%fixed4) then
      flag(4) = 1
    else
      flag(4) = 0
    endif
    
    if(m%fixed5) then
      flag(5) = 1
    else
      flag(5) = 0
    endif
    
    if(m%fixed6) then
      flag(6) = 1
    else
      flag(6) = 0
    endif
    
    if(m%equal12) then
      flag(7) = 1
    else
      flag(7) = 0
    endif
    
    if(m%equal13) then
      flag(8) = 1
    else
      flag(8) = 0
    endif
    
    if(m%equal23) then
      flag(9) = 1
    else
      flag(9) = 0
    endif                               
    
  end subroutine getValueFlag_modelParam
    


  subroutine getQrow_modelparam(iPar,nFunc,mx,isComplex,Qreal,Qimag,Resp)
      implicit none
      integer :: nFunc
      integer :: mx(3),iPar
      logical :: isComplex
      type(modelParam_t), intent(inout)     :: Qreal(:), Qimag(:)
      real(kind=prec),intent(in)	:: Resp(:)
      ! local variables
      integer iFunc, iComp

      iComp = 0
      do iFunc  = 1, nFunc
	      if(isComplex) then
	         iComp = iComp + 1
	         Qreal(iFunc)%cellRho(iPar)%v(mx(1),mx(2),mx(3)) = Resp(iComp)
	         iComp = iComp + 1
	         Qimag(iFunc)%cellRho(iPar)%v(mx(1),mx(2),mx(3)) = Resp(iComp)
	         Qreal(iFunc)%zeroValued = .false.
	         Qimag(iFunc)%zeroValued = .false.
	      else
	         iComp = iComp + 1
	         Qreal(iFunc)%cellRho(iPar)%v(mx(1),mx(2),mx(3)) = Resp(iComp)
	         Qimag(iFunc)%cellRho(iPar)%v(mx(1),mx(2),mx(3)) = R_ZERO
	         Qreal(iFunc)%zeroValued = .false.
	         Qimag(iFunc)%zeroValued = .true.
	      endif
	    enddo	   
	  
  end subroutine getQrow_modelparam
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine AnalyzeAniModel(m0,mReal)
    implicit none
    type(modelParam_t), intent(in)     :: m0
    type(modelParam_t), intent(inout)     :: mReal 
    
    if((m0%equal12.and.m0%equal23) .or. (m0%equal12.and.m0%equal13)  &
       .or. (m0%equal13.and.m0%equal23)) then
       ! isotropy, 10/10/10
       call add_rscalar(mReal%cellRho(1), mReal%cellRho(2), mReal%cellRho(1))
       call add_rscalar(mReal%cellRho(1), mReal%cellRho(3), mReal%cellRho(1)) 
       mReal%cellRho(2)%v = mReal%cellRho(1)%v
       mReal%cellRho(3)%v = mReal%cellRho(1)%v
    elseif(m0%equal12.and.(.not.m0%equal23).and.(.not.m0%equal13)) then
       ! Vertical anisotropy or dipping anisotropy, S angle fails, 10/10/100
       call add_rscalar(mReal%cellRho(1), mReal%cellRho(2), mReal%cellRho(1)) 
       mReal%cellRho(2)%v = mReal%cellRho(1)%v
    elseif(m0%equal23.and.(.not.m0%equal12).and.(.not.m0%equal13)) then      
       ! biaxial anisotropy or horizontal anisotropy, 10/100/100
       call add_rscalar(mReal%cellRho(2), mReal%cellRho(3), mReal%cellRho(2))  
       mReal%cellRho(3)%v = mReal%cellRho(2)%v
    elseif(m0%equal13.and.(.not.m0%equal12).and.(.not.m0%equal23)) then      
       ! biaxial anisotropy or horizontal anisotropy, 100/10/100
       call add_rscalar(mReal%cellRho(1), mReal%cellRho(3), mReal%cellRho(1)) 
       mReal%cellRho(3)%v = mReal%cellRho(1)%v           
    endif
    
    ! if fixed, then the corresponding model parameter do not updata, mHat%cellRho(ia) = 0
    if(m0%fixed1) then
      call zero_rscalar(mReal%cellRho(1))
    endif
    if(m0%fixed2) then
      call zero_rscalar(mReal%cellRho(2))
    endif
    if(m0%fixed3) then
      call zero_rscalar(mReal%cellRho(3))
    endif
    if(m0%fixed4) then
      call zero_rscalar(mReal%cellRho(4))
    endif
    if(m0%fixed5) then
      call zero_rscalar(mReal%cellRho(5))
    endif
    if(m0%fixed6) then
      call zero_rscalar(mReal%cellRho(6))
    endif
      
  end subroutine AnalyzeAniModel 
  
    !---------------------------------------------------------------------------------------------add for lbfgs begining
    function nd_lbfgs(m0_lbfgs) result(ndim_lbfgs)    ! for lbfgs wkp add. 2016.4.6, KONG-2018-10-24

    type(modelParam_t), intent(in)	  :: m0_lbfgs
    integer :: ndim_lbfgs

    ndim_lbfgs = 6 * m0_lbfgs%Nx * m0_lbfgs%Ny * m0_lbfgs%NzEarth !少一维

    end function nd_lbfgs

    subroutine get_x_g(m_in,ndim_lbfgs,m_out)

    type(modelParam_t), intent(in)	  :: m_in
    integer , intent(in)   :: ndim_lbfgs
    real*8	, intent(inout)  :: m_out(ndim_lbfgs)

    integer :: ix,iy,iz,cont,ia

    cont = 1
    do ia = 1,6
      do ix=1,m_in%Nx
         do iy=1,m_in%Ny
            do iz=1,m_in%NzEarth
                m_out(cont) = m_in%cellRho(ia)%v(ix,iy,iz)
                cont = cont + 1
                !write(*,*) cont,ndim_lbfgs
            enddo
         enddo
      enddo
    enddo
    
    end subroutine get_x_g


    subroutine x_to_mHat(m_in,ndim_lbfgs,m_out)

    integer , intent(in)  :: ndim_lbfgs
    real*8	, intent(in)  :: m_in(ndim_lbfgs)   
    type(modelParam_t), intent(inout)	  :: m_out      

    integer :: ix,iy,iz,cont,ia
    
    !  write(*,*) "aaaaa",m_out%Ny,m_out%NzEarth
    cont = 1
    do ia = 1,6
      do ix=1,m_out%Nx
         do iy=1,m_out%Ny
            do iz=1,m_out%NzEarth
                m_out%cellRho(ia)%v(ix,iy,iz) = m_in(cont)
                cont = cont + 1
            enddo
         enddo
      enddo
    enddo
    end subroutine x_to_mHat
    !---------------------------------------------------------------------------------------------add for lbfgs ending
    
    ! added by Kong, 2018-12-19
    subroutine cal_anisCsNorm(m,aniNorm,lambdaANI)
      implicit none
      type(modelParam_t),intent(in) :: m
      real(kind=prec),intent(inout) :: aniNorm
      type(DampingFactorANI),intent(in) :: lambdaANI
      ! locals
      integer nx,ny,nzEarth
      
      nx = m%Nx
      ny = m%Ny
      nzEarth = m%NzEarth 
      
      call anisCsNorm(nx,ny,nzEarth,lambdaANI,m%cellRho(1)%v,m%cellRho(2)%v,m%cellRho(3)%v,aniNorm)
      
    end subroutine cal_anisCsNorm
    

    subroutine cal_anisCsVec(m,aniVec,lambdaANI)
      implicit none
      type(modelParam_t),intent(in) :: m
      type(modelParam_t),intent(inout) :: aniVec 
      type(DampingFactorANI),intent(in) :: lambdaANI
      ! locals
      integer nx,ny,nzEarth
      
      nx = m%Nx
      ny = m%Ny
      nzEarth = m%NzEarth
      
      call zero_modelParam(aniVec)
      
      call anisCsVec(nx,ny,nzEarth,lambdaANI,m%cellRho(1)%v,m%cellRho(2)%v,m%cellRho(3)%v, &
              aniVec%cellRho(1)%v,aniVec%cellRho(2)%v,aniVec%cellRho(3)%v)
      
    end subroutine cal_anisCsVec
        
end module ModelSpace
