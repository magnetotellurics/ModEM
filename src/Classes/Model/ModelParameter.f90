module ModelParameter
  use Constants
  use rScalar
  use rVector

  public :: LOGE, LOG_10, LINEAR
  
  !**
  ! Supported model parameter types (conductivity only).
  !*
  character (len = 4), parameter :: LOGE   = 'LOGE'
  character (len = 5), parameter :: LINEAR = 'LINEAR'
  character (len = 6), parameter :: LOG_10 = 'LOG_10'

  type, abstract :: ModelParameter_t
     integer             :: mKey
     character(len = 80) :: paramType   = ''     
     real(kind = prec)   :: airCond     = 1E-7_prec
     logical             :: zeroValued  = .false.
     logical             :: isAllocated = .false.
     
     procedure(iface_SigMap), pointer, nopass :: SigMap_ptr => SigMap_Linear
     
   contains
     procedure(iface_Length) , deferred, public :: Length     
     
     procedure(iface_Zeros)   , deferred, public :: Zeros
     procedure(iface_CopyFrom), deferred, public :: CopyFrom

     ! Model mapping methods
     procedure(iface_PDEmapping)  , deferred, public :: PDEmapping
     procedure(iface_dPDEmapping) , deferred, public :: dPDEmapping
     procedure(iface_dPDEmappingT), deferred, public :: dPDEmappingT
     
     procedure, public :: SigMap
     procedure, public :: SetSigMap
     
     procedure(iface_SetType), deferred, public :: SetType
     procedure, public :: GetType
     
  end type ModelParameter_t

  abstract interface
     
     pure function iface_SigMap(x, p_job) result(y)
       import :: prec
       real(kind = prec), intent(in) :: x
       character(*)     , intent(in), optional :: p_job
       real(kind = prec) :: y
     end function iface_SigMap
     
     !**
     ! Length
     ! This just returns total number of model parameters
     !*
     function iface_Length(self) result(nParam)
       import :: ModelParameter_t
       class(ModelParameter_t), intent(in) :: self
       integer :: nParam
     end function iface_Length
     
     !**
     ! Zeros
     ! Set model parameter to zero.
     !*
     subroutine iface_Zeros(self)
       import :: ModelParameter_t
       class(ModelParameter_t), intent(inout) :: self
     end subroutine iface_Zeros
     
     !**
     ! Copy model parameter.
     !*
     subroutine iface_CopyFrom(self, rhs)
       import :: ModelParameter_t       
       class(ModelParameter_t), intent(inout) :: self
       class(ModelParameter_t), intent(in)    :: rhs
     end subroutine iface_CopyFrom

     !**
     ! SetType
     !*
     subroutine iface_SetType(self, paramType)
       import :: ModelParameter_t
       class(ModelParameter_t), intent(inout) :: self
       character(*)           , intent(in)    :: paramType
     end subroutine iface_SetType

     !**
     ! These are what interfaces would be if ModelMap
     ! is a distinct object;
     ! If these methods were in an extension of
     ! ModelParameter object then m (or m0) would
     ! be omitted (self would be the model parameter already).
     !*
     function iface_PDEmapping(self) result(eVec)
       import :: ModelParameter_t, rVector_t
       class(ModelParameter_t), intent(inout)  :: self
       class(rVector_t), allocatable :: eVec
     end function iface_PDEmapping
     
     !**
     !
     !*
     function iface_dPDEmapping(self, dm) result(eVec)
       import :: ModelParameter_t, rVector_t
       class(ModelParameter_t), intent(in)  :: self
       class(ModelParameter_t), intent(in)  :: dm
       class(rVector_t), allocatable :: eVec
     end function iface_dPDEmapping
     
     !**
     ! NOTE: For transpose (adjoint) dm is output`
     !       and eVec is input.
     !*
     function iface_dPDEmappingT(self, eVec) result(dm)
       import :: ModelParameter_t, rVector_t
       class(ModelParameter_t), intent(in)  :: self
       class(rVector_t)       , intent(in)  :: eVec
       class(ModelParameter_t), allocatable :: dm
     end function iface_dPDEmappingT
     
  end interface

contains
  
  elemental function SigMap(self, x, job) result(y)
    class(ModelParameter_t), intent(in) :: self
    real(kind = prec), intent(in) :: x
    character(*)     , intent(in), optional :: job
    ! Local variables
    real(kind = prec) :: y
    
    y = self%Sigmap_ptr(x)
    
  end function SigMap
  
  subroutine SetSigMap(self, paramType)
    class(ModelParameter_t) :: self
    character(*), intent(in) :: paramType
    
    select case(paramType)
    case('LOGE')
       self%SigMap_ptr => SigMap_Log       
    case('LINEAR')
       self%SigMap_ptr => SigMap_Linear
    endselect
    
  end subroutine SetSigMap
  
  pure function SigMap_Linear(x, job) result(y)
    real(kind = prec), intent(in) :: x
    character(*)     , intent(in), optional :: job
    real(kind = prec) :: y
    
    y = x
  end function SigMap_Linear
  
  pure function SigMap_Log(x, p_job) result(y)
    real(kind = prec), intent(in) :: x
    character(*)     , intent(in), optional :: p_job
    real(kind = prec) :: y
    ! Local variables
    character(30) :: job
    
    if (.not.present(p_job)) then
       job = 'forward'
    else
       job = p_job
    end if
  end function SigMap_Log

  function GetType(self) result(pType)
    ! Arguments
    class(ModelParameter_t), intent(in) :: self
    ! Local variables
    character(len = 80) :: pType

    pType = trim(self%paramType)
  end function GetType
  
end module ModelParameter
