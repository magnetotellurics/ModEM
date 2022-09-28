!**
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
!
! VTI - Anisotropic case
!*

module ModelSpace
  
  use gridcalc
  use file_units
  use math_constants
  use utilities
  use sg_scalar
  use sg_vector
  use sg_sparse_vector
#ifdef MPI
  use Declaration_MPI
#endif
  
  implicit none
  
  ! supported model parameter types (conductivity only)
  character(len = 80), parameter :: LOGE   = 'LOGE'
  character(len = 80), parameter :: LOG_10 = 'LOG10'
  character(len = 80), parameter :: LINEAR = 'LINEAR'

  
  type :: modelParam_t
     !**
     ! Parameters which define conductivity distribution.
     !
     ! For the present approach, with each cell of the computational grid
     ! allowed a separate (constant) conductivity) a separate type here
     ! is hardly needed ...
     !
     ! The idea is that the conductivity parameter should be treated as an
     ! "abstract data type", defined along with a set of routines for
     ! mapping to the internal representation of conductivity used directly
     ! by the forward solver.
     !
     !*
     
     
     integer               :: Nx, Ny, NzEarth
     type(rscalar)         :: cellCond_h    ! Horizontal conductivity
     type(rscalar)         :: cellCond_v    ! Vertical conductivity
     type(grid_t), pointer :: grid
     real (kind = prec)    :: AirCond
     
     logical :: allocated = .false.
     logical :: isVTI = .false.
     
     ! This logical is set true by zero_modelParam ONLY.
     logical :: zeroValued = .false.

     ! Necessary to avoid memory leaks; only true for function outputs.
     logical :: temporary = .false.

     ! Another logical that gets set to true every time the model
     ! parameter is modified; used by the ForwardSolver for updateCond.
     
     logical :: updated = .false.

     ! Supported paramType at present: LINEAR, LOG10 and LOGE
     character (len = 80) :: paramType = ''
     
  end type modelParam_t
  
  character(len = 80), save :: userParamType = 'LOGE'
  
  interface assignment (=)
     MODULE PROCEDURE copy_modelParam
  end interface assignment (=)
  
  interface zero
     MODULE PROCEDURE zero_modelParam
  end interface zero
  
  interface iszero
     MODULE PROCEDURE iszero_modelParam
  end interface iszero
  
  interface scMult ! operator (*)
     MODULE PROCEDURE scMult_modelParam
  end interface scMult
  
  interface scMultAdd
     MODULE PROCEDURE scMultAdd_modelParam
  end interface scMultAdd
  
  interface dotProd
     MODULE PROCEDURE dotProd_modelParam
  end interface dotProd
  
  interface linComb
     MODULE PROCEDURE linComb_modelParam
  end interface linComb
  
  interface deall
     MODULE PROCEDURE deall_modelParam
  end interface deall
  
  interface countModelParam
     MODULE PROCEDURE count_modelParam
  end interface countModelParam
  
  !  I/O interfaces
  
  interface write_modelParam
     MODULE PROCEDURE write_modelParam_WS
  end interface write_modelParam
  
  interface read_modelParam
     MODULE PROCEDURE read_modelParam_WS
  end interface read_modelParam
  
  interface writeVec_modelParam
     MODULE PROCEDURE writeVec_modelParam_binary
  end interface writeVec_modelParam
  
  interface readVec_modelParam
     MODULE PROCEDURE readVec_modelParam_binary
  end interface readVec_modelParam
  
  ! definitions for CmSqrt: must be consistent with the include file below
  
#include "modelCov/RecursiveAR.hd"
  !#include "modelCov/Diffusion.hd"
contains
  
  !**
  ! Routines which define mappings between the "natural" representation
  ! of conductivity/resistivity on the model grid and the formal
  ! model parameter structure.
#include "ModelMap.inc"
  
  ! The included file must contain subroutines
  ! create_CmSqrt, deall_CmSqrt, multBy...
#include "modelCov/RecursiveAR.inc"
  !#include "modelCov/Diffusion.inc"
  !  I/O choices
#include "modelParamIO/Binary.inc"
#include "modelParamIO/Mackie.inc"
#include "modelParamIO/WS.inc"
  
  !  MPI model parameter, if needed
#ifdef MPI
#include "ModelParam_MPI.inc"
#endif

  !**
  ! create_modelParam allocates and initializes arrays for
  ! conductivity parameter structure
  ! Pass grid of type grid_t to set array sizes
  ! optional arguments v and vAir are assumed to be consistent
  ! with paramType.
  subroutine create_modelParam(grid, paramtype, m, v_h, v_v, vAir)
    implicit none
    ! Arguments
    type (grid_t)      , intent(in), target   :: grid
    character(80)      , intent(in)           :: paramtype
    type (modelParam_t), intent(inout)        :: m
    type (rscalar)     , intent(in), optional :: v_h, v_v
    real (kind = prec) , intent(in), optional :: vAir
    !
    !***********************
    ! Executable statements
    !***********************
    if (m%allocated) then
       call deall_modelParam(m)
    end if

    call create_rscalar(grid, m%cellCond_h, CELL_EARTH)
       
    m%Nx = grid%Nx
    m%Ny = grid%Ny
    m%NzEarth = grid%NzEarth
    m%grid => grid
    m%paramType = paramtype
    m%allocated = .true.
    
    if (present(v_h)) then
       m%cellCond_h = v_h
    end if

    if (present(v_v)) then
       call create_rscalar(grid, m%cellCond_v, CELL_EARTH)
       m%cellCond_v = v_v
       m%isVTI = .true.
    else if (grid%isVTI) then
       call create_rscalar(grid, m%cellCond_v, CELL_EARTH)
       m%isVTI = .true.
    else
       m%isVTI = .false.
    end if

    if (present(vAir)) then
       m%AirCond = vAir
    else
       if (paramtype .eq. LOGE) then
          m%AirCond = log(SIGMA_AIR)
       else if(paramtype .eq. LOG_10) then
          m%AirCond = log10(SIGMA_AIR)
       else
          m%AirCond = SIGMA_AIR
       end if
    end if
    
    m%updated = .true.
    m%zeroValued = .false.
    
  end subroutine create_modelParam
  
  !**
  !
  !*
  subroutine deall_modelParam(m)
    implicit none
    ! Arguments
    type (modelParam_t) :: m
    !
    !***********************
    ! Executable statements
    !***********************
    !
    if (m%allocated) then
       if (m%isVTI) then
          call deall_rscalar(m%cellCond_v)
       end if

       call deall_rscalar(m%cellCond_h)
       
       nullify(m%grid)
       
       m%allocated = .false.
       m%zeroValued = .false.
       m%paramType = ''
    end if
    
  end subroutine deall_modelParam
  
  !**
  !
  !*
  subroutine getType_modelParam(m, paramType)
    ! Arguemtns
    type(modelParam_t), intent(in)  :: m
    character(*)      , intent(out) :: paramType

    paramType = trim(m%paramType)
  end subroutine getType_modelParam
  
  !**
  ! Converts the input model parameter structure to paramType, by
  ! comparing paramType with m%paramType and performing the necessary
  ! computations if the two strings differ; assumes that m is allocated.
  subroutine setType_modelParam(m, paramType)
    implicit none
    ! Arguments
    type(modelParam_t), intent(inout) :: m
    character(*)      , intent(in)    :: paramType
    !
    !***********************
    ! Executable statements
    !***********************
    !
    if (.not.(m%allocated)) then
       call errstop('modelParam must be allocated before calling'//&
            'setType_modelParam')
    end if
    
    if (trim(paramType) .eq. trim(m%paramType)) then
       ! we are done
    else if(m%paramType == LINEAR) then
       ! convert to log
       if (paramType == LOGE) then
          if (m%isVTI) then
             m%cellCond_v%v = log(m%cellCond_v%v)
          end if
          m%cellCond_h%v = log(m%cellCond_h%v)
          m%AirCond = log(m%AirCond)

       else if(paramType == LOG_10) then
          if (m%isVTI) then
             m%cellCond_v%v = log10(m%cellCond_v%v)
          end if
          m%cellCond_h%v = log10(m%cellCond_h%v)
          m%AirCond = log10(m%AirCond)
       end if
       
    else if(paramType == LINEAR) then
       ! convert from log to linear
       if (m%paramType == LOGE) then
          if (m%isVTI) then
             m%cellCond_v%v = exp(m%cellCond_v%v)
          end if
          m%cellCond_h%v = exp(m%cellCond_h%v)
          m%AirCond = exp(m%AirCond)
          
       else if(m%paramType == LOG_10) then
          if (m%isVTI) then
             m%cellCond_v%v = exp(m%cellCond_v%v * log(10.))
          end if
          m%cellCond_h%v = exp(m%cellCond_h%v * log(10.))
          m%AirCond = exp(m%AirCond * log(10.))
       end if
       
    else if ((m%paramType == LOGE) .and. (paramType == LOG_10)) then
       ! convert from natural log to log10
       if (m%isVTI) then
          m%cellCond_v%v = m%cellCond_v%v / log(10.)
       end if
       m%cellCond_h%v = m%cellCond_h%v / log(10.)
       m%AirCond = m%AirCond / log(10.)
       
    else if ((m%paramType == LOG_10) .and. (paramType == LOGE)) then
       ! convert from log10 to natural log
       if (m%isVTI) then
          m%cellCond_v%v = m%cellCond_v%v * log(10.)
       end if
       m%cellCond_h%v = m%cellCond_h%v * log(10.)
       m%AirCond = m%AirCond * log(10.)
       
    else
       call errstop('unknown paramType in setType_modelParam')
    endif
    
    m%paramType = paramType
    
    return
    
  end subroutine setType_modelParam
  
  !**
  ! dot product of two model space parameter objects
  ! here implemented for 3D numerical grid.
  !*
  function dotProd_modelParam(m1, m2) result(r)
    implicit none
    ! Arguments
    real(kind = prec) :: r
    type(modelParam_t), intent(in) :: m1, m2
    ! Local variables
    integer :: j, k
    !
    !***********************    
    ! Executable statements
    !***********************
    !
    if ((m1%Ny .ne. m2%Ny).or. (m1%Nx .ne. m2%Nx).or. &
         (m1%NzEarth .ne. m2%NzEarth)) then
       write(6,*) 'Nx, Ny, Nz ',m1%Nx, m2%Nx,    &
            m1%Ny,m2%Ny,m1%NzEarth,m2%NzEarth
       call errStop('size of m1, m2 incompatable in dotProd_modelParam')
    end if
    
    ! if one of the model parameters is zero valued, no need to compute r
    r = R_ZERO
    if (.not. m1%zeroValued .and. .not. m2%zeroValued) then
       ! TODO - Anisotropic conductivity
       r = dotProd_rscalar_f(m1%cellCond_h, m2%cellCond_h)
    end if
    
  end function dotProd_modelParam
  
  !**
  !  zeros a model space object
  !*
  subroutine zero_modelParam(m)
    implicit none
    ! Arguments 
    type(modelParam_t), intent(inout) :: m
    !
    !***********************
    ! Executable statements
    !***********************
    !
    if (m%isVTI) then
       call zero_rscalar(m%cellCond_v)
    end if
    call zero_rscalar(m%cellCond_h)
    
    m%updated = .true.
    m%zeroValued = .true.
    
  end subroutine zero_modelParam
    
  !**
  !
  !*
  logical function iszero_modelParam(m) result (zeroValued)
    implicit none
    ! Arguments
    type(modelParam_t), intent(in) :: m
    !
    !***********************
    ! Executable statements
    !***********************
    !
    zeroValued = m%zeroValued
    
  end function iszero_modelParam
  
  !**
  ! 
  ! generated a random model parameter perturbation [0,1) in log
  ! space; otherwise, exp of that in linear space.
  !*.
  subroutine random_modelParam(m, eps)
    implicit none
    ! Arguments
    type(modelParam_t), intent(inout)  :: m
    real(kind = prec) , intent(in), optional :: eps
    !
    !***********************
    ! Executable statements
    !***********************
    !
    if (.not. m%allocated) then
       call errStop('Model parameter not allocated in random_modelParam.')
    end if
    
    if (present(eps)) then
       if (m%isVTI) then
          call random_rscalar(m%cellCond_v, eps)
       end if
       call random_rscalar(m%cellCond_h, eps)       
    else
       if (m%isVTI) then
          call random_rscalar(m%cellCond_v)
       end if
       call random_rscalar(m%cellCond_h)
    end if
    
    if (m%paramType == LINEAR) then
       if (m%isVTI) then
          m%cellCond_v%v = exp(m%cellCond_v%v)
       end if
       m%cellCond_h%v = exp(m%cellCond_h%v)
    end if
    
    m%updated = .true.
    m%zeroValued = .false.
    
  end subroutine random_modelParam
  
  !**
  ! Forms the linear combination of model parameters
  !     m = a1*m1 + a2*m2
  ! where a1 and a2 are real constants and m1 and m2
  ! are model parameters output m may overwrite m1 or m2.
  !*
  subroutine linComb_modelParam(a1, m1, a2, m2, m)
    implicit none
    ! Arguments
    real(kind = prec) , intent(in)    :: a1, a2
    type(modelParam_t), intent(in)    :: m1, m2
    type(modelParam_t), intent(inout) :: m
    !
    !***********************
    ! Executable statements
    !***********************
    !
    if ((m1%Ny .ne. m2%Ny).or. (m1%Nx .ne. m2%Nx) .or. &
         (m1%NzEarth .ne. m2%NzEarth)) then
       call errStop('size of m1, m2 incompatable in linComb_modelParam.')
    end if

    if (m1%paramType .ne. m2%paramType) then
       call errStop('paramType incompatable in linComb_modelParam.')
    end if
    
    ! Make sure m is allocated, and of same size; if not same
    ! size, deallocate and allocate as correct size; otherwise
    ! do nothing (use m as input ... allowing m to overwrite m1)
    if (m%allocated) then
       if ((m1%Ny .ne. m2%Ny).or. (m1%Nx .ne. m2%Nx) .or. &
            (m1%NzEarth .ne. m2%NzEarth)) then
          call deall_modelParam(m)
          call create_modelParam(m1%grid, m1%paramtype, m)
       end if
    else
       call create_modelParam(m1%grid, m1%paramtype, m)
    end if
    
    if (m%isVTI) then
       m%cellCond_v%v = a1*m1%cellCond_v%v + a2*m2%cellCond_v%v
    end if
    m%cellCond_h%v = a1*m1%cellCond_h%v + a2*m2%cellCond_h%v
    
    ! We are implicitly assuming that if m1 and m2 are the same
    ! size other parameters are of the same type.
    m%AirCond = m1%AirCond
    m%zeroValued = m1%zeroValued .and. m2%zeroValued
    m%updated = .true.
    
  end subroutine linComb_modelParam
  
  !**
  ! Computes mOut = a * mIn for modelParam object
  ! 'm' and real scalar 'a'.
  subroutine scMult_modelParam(a, mIn, mOut)
    implicit none
    ! Arguments
    real (kind = prec), intent(in)    :: a
    type(modelParam_t), intent(in)    :: mIn
    type(modelParam_t), intent(inout) :: mOut
    !
    !***********************    
    ! Executable statements
    !***********************
    !
    
    ! check to see that input m is allocated
    if (.not.mIn%allocated) then
       call errStop('IUnput not allocated on call to scMult_modelParam.')
    end if
    
    ! if mIn is zero valued, no need to compute mOut.
    if (mIn%zeroValued) then
       mOut = mIn
    else
       call linComb_modelParam(R_ZERO, mIn, a,mIn, mOut)
    end if
    
  end subroutine scMult_modelParam
  
  !**
  ! Computes mOut = a * mIn + mOut for modelParam
  ! object m and real scalar 'a'.
  !*
  subroutine scMultAdd_modelParam(a, mIn, mOut)
    implicit none
    ! Arguments
    real (kind = prec), intent(in)    :: a
    type(modelParam_t), intent(in)    :: mIn
    type(modelParam_t), intent(inout) :: mOut
    !
    !***********************
    ! Executable statements
    !***********************
    !
    
    ! Ceck to see that input m is allocated
    if (.not.(mIn%allocated .and. mOut%allocated)) then
       call errStop('input not allocated on call to scMultAdd_modelParam.')
    end if
    
    ! if mIn is zero valued, no need to compute mOut
    if (.not. mIn%zeroValued) then
       call linComb_modelParam(a, mIn, ONE, mOut, mOut)
    end if
    
  end subroutine scMultAdd_modelParam
  
  !**
  ! Sets cell conductivities in a model parameter object m
  ! the values of v and vAir are determined by paramType;
  ! if different from that of the model parameter, returns
  ! an error. To avoid this error, first convert m to the
  ! required paramType using setType_modelParam.
  !*
  subroutine setValue_modelParam(m, paramType, v_h, vAir, v_v)
    implicit none
    ! Arguments
    type(modelParam_t), intent(inout)        :: m
    character(80)     , intent(in)           :: paramType
    type(rscalar)     , intent(in)           :: v_h
    real(kind = prec) , intent(in), optional :: vAir
    type(rscalar)     , intent(in), optional :: v_v

    !
    !***********************
    ! Executable statements
    !***********************
    !
    if (.not.(m%allocated)) then
       call errstop('Output modelParam must be allocated'//&
            'before calling setValue_modelParam.')
    end if
    
    ! Error checking
    if ((m%Ny .ne. v_h%Ny) .or. (m%Nx .ne. v_h%Nx) .or.&
         (m%NzEarth .ne. v_h%Nz)) then
       call errstop('modelParam/rscalar dimensions disagree'//&
            ' in setValue_modelParam.')
    else if(paramType .ne. m%paramType) then
       call errstop('paramTypes not consistent in setValue_modelParam.')
    end if
    
    ! Set values
    m%cellCond_h = v_h

    if (present(v_v)) then
       if (m%isVTI) then
          m%cellCond_v = v_v
       else
          m%cellCond_v = m%cellCond_h
       end if
    end if
    
    if(present(vAir)) then
       m%AirCond = vAir
    end if
    
    m%updated = .true.
    m%zeroValued = .false.
    
  end subroutine setValue_modelParam
  
  !**
  ! Gets cell conductivities from a model parameter object m
  !
  ! Extracts the values of v and vAir;
  ! values that are extracted are converted to paramType.
  ! Different from ModelParamToCell in that the value
  ! that gets extracted is exactly what is stored in the modelParam:
  ! it does not contain any air layers. This is needed for BC_x0_WS.
  subroutine getValue_modelParam(m, paramType, v_h, vAir, v_v)
    implicit none
    ! Arguments
    type(modelParam_t), intent(in)    :: m
    character(80)     , intent(inout) :: paramType
    type(rscalar)     , intent(out)   :: v_h
    real(kind = prec)  , intent(out), optional :: vAir
    type(rscalar)     , intent(out), optional :: v_v

    ! Local variables
    type(modelParam_t) :: mTemp
    !
    !***********************
    ! Executable statements
    !***********************
    !
    if (.not.(m%allocated)) then
       call errstop('Input modelParam must be allocated'//&
            'before calling getValue_modelParam.')
    end if
    
    if (trim(paramType) .eq. '') then
       paramType = m%paramType
    end if
    
    if (v_h%allocated) then
       call deall_rscalar(v_h)
    end if

    if (present(v_v)) then
       if (v_v%allocated) then
          call deall_rscalar(v_v)
       end if
    end if
    
    ! create a temporary copy
    mTemp = m
    
    ! convert model to the required type
    call setType_modelParam(mTemp, paramType)
    
    ! set values
    v_h = mTemp%cellCond_h
    
    if (present(v_v)) then
       if (m%isVTI) then
          v_v = mTemp%cellCond_v
       else
          v_v = v_h
       end if
    end if
    
    if (present(vAir)) then
       vAir = mTemp%AirCond
    end if
    
    ! deallocate temporary model parameter
    call deall_modelParam(mTemp)
    
  end subroutine getValue_modelParam
  
  !**
  ! if mOut is allocated, check to see if if is of same size as
  ! mIn; if not, deallocate and reallocate as correct size; otherwise
  ! use as input.
  !*
  subroutine copy_modelParam(mOut, mIn)
    implicit none
    ! Arguments
    type(modelParam_t), intent(in)    :: mIn
    type(modelParam_t), intent(inout) :: mOut
    !
    !***********************
    ! Executable statements
    !***********************
    !    
    if (mOut%allocated) then
       if ((mOut%Ny .ne. mIn%Ny).or. (mOut%Nx .ne. mIn%Nx) .or. &
            (mOut%NzEarth .ne. mIn%NzEarth)) then
          call deall_modelParam(mOut)
          call create_modelParam(mIn%grid, mIn%paramtype, mOut)
       end if
    else
       call create_modelParam(mIn%grid, mIn%paramtype, mOut)
    end if
    
    mOut%cellCond_h%v = mIn%cellCond_h%v

    if (mIn%isVTI.neqv.mOut%isVTI) then
       write(*, *) 'copy_modelParam: Both mIn and mOut need to be VTI'
       STOP
    end if

    if (mIn%isVTI.and.mOut%isVTI) then
       mOut%cellCond_v%v = mIn%cellCond_v%v
    end if

    mOut%AirCond = mIn%AirCond
    mOut%zeroValued = mIn%zeroValued

    ! no need to set updated to TRUE - this just makes a copy.
    ! clean up the memory if this is used with an "=" sign.
    if (mIn%temporary) then
       call deall_modelParam(mIn)
    endif
    
  end subroutine copy_modelParam
  
  !**
  ! count_modelParam counts the number of variable model parameters.
  !*
  function count_modelParam(cond) result (N)
    implicit none
    ! Arguments
    type (modelParam_t), intent(in) :: cond
    ! Local variables
    integer :: N
    !
    !***********************
    ! Executable statements
    !***********************
    !
    if (.not.cond%allocated) then
       call errStop('Model parameter not allocated in count_modelParam.')
    end if
    
    N = cond%Nx * cond%Ny * cond%NzEarth
    
  end function count_modelParam
  
  !**
  ! Extracts the grid from a modelParam object; this is only needed since
  ! the attributes are private (in F2003, can declare everything private
  ! while grid and allocated attributes could be public).
  !*
  subroutine getGrid_modelParam(grid, mIn)
    implicit none
    ! Arguments
    type (grid_t)      , intent(out) :: grid
    type (modelParam_t), intent(in)  :: mIn
    !
    !***********************    
    ! Executabl estatements
    !***********************
    !
    if (.not. mIn%allocated) then
       call warning('Model vector not allocated on call to getGrid_modelParam.')
    end if
    
    grid = mIn%grid
    
  end subroutine getGrid_modelParam
  
  !**
  ! Sets modelParam%updated to FALSE; this is only needed since
  ! the attributes are private (in F2003, can declare everything private
  ! while grid and allocated attributes could be public)
  ! used when the model parameter is no longer considered "new",
  ! e.g. after the updateModelData routine in ForwardSolver.
  !*
  subroutine setValueUpdated_modelParam(m)
    implicit none
    ! Arguments
    type (modelParam_t), intent(inout) :: m
    !
    !***********************
    ! Executable statements
    !***********************
    !
    m%updated = .false.
    
  end subroutine setValueUpdated_modelParam
  
  !**
  ! Checks whether a modelParam is updated; this is only needed since
  ! the attributes are private (in F2003, can declare everything private
  ! while grid and allocated attributes could be public).
  !*
  subroutine getValueUpdated_modelParam(m, updated)
    implicit none
    ! Arguments
    type (modelParam_t), intent(in)  :: m
    logical            , intent(out) :: updated
    !
    !***********************
    ! Executable statements
    !***********************
    !
    updated = m%updated
    
  end subroutine getValueUpdated_modelParam
  
end module ModelSpace

