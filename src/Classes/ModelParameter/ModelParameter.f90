module ModelParameter
   !
   use Constants
   use rScalar
   use rVector
   use Grid2D
   use ModelParameter1D
   use ModelParameter2D
   use Grid
   use MetricElements
   !
   type, abstract :: ModelParameter_t
       !
	   !
       ! Pointer to the original grid
       class( Grid_t ), pointer :: grid

       ! Pointer to metric elements -- useful for model mappings
       !    provides Viedge, Vcell
       class( MetricElements_t ), pointer :: metric
	   !
	   integer             :: mKey
       character(len = 80) :: paramType    = ''       
       real(kind = prec)   :: airCond       = 1E-7_prec
       logical             :: zeroValued   = .false.
       logical             :: isAllocated = .false.
       !
       procedure( iface_SigMap ), pointer, nopass :: SigMap_ptr => SigMap_Linear
       !
    contains
	   !
       procedure( interface_length_model_parameter ) , deferred, public :: length
       !
       procedure( interface_zeros_model_parameter )    , deferred, public :: zeros
       procedure( interface_copy_from_model_parameter ), deferred, public :: copyFrom
       !
       ! Model mapping methods
       procedure( iface_PDEmapping ), deferred, public   :: PDEmapping
       procedure( iface_dPDEmapping ), deferred, public  :: dPDEmapping
       procedure( iface_dPDEmappingT ), deferred, public :: dPDEmappingT
	   !
	   procedure( interface_slice_1d_model_parameter ) , deferred, public :: Slice1D
	   procedure( interface_slice_2d_model_parameter ) , deferred, public :: Slice2D
	   !
	   procedure( interface_avg_model_1d_model_parameter ) , deferred, public :: AvgModel1D
       procedure, public :: SigMap
       procedure, public :: SetSigMap
       !
       procedure(iface_SetType), deferred, public :: SetType
       procedure, public :: GetType
	   procedure, public :: setMetric => setMetricModelOperator
       !
   end type ModelParameter_t

   abstract interface
       !
       function interface_slice_1d_model_parameter( self, ix, iy ) result( m1D )
			import :: ModelParameter_t, ModelParameter1D_t
			class( ModelParameter_t ), intent(in) :: self
			integer, intent(in)                   :: ix, iy
			type( ModelParameter1D_t )            :: m1D 
       end function interface_slice_1d_model_parameter
	   !
       function interface_avg_model_1d_model_parameter( self ) result( m1D )
			import :: ModelParameter_t, ModelParameter1D_t
			class( ModelParameter_t ), intent(in) :: self
			type( ModelParameter1D_t )            :: m1D 
       end function interface_avg_model_1d_model_parameter
	   !
       function interface_slice_2d_model_parameter( self, axis, j ) result( m2D )
			import :: ModelParameter_t, ModelParameter2D_t
			class( ModelParameter_t ), intent(in) :: self
			integer, intent( in )                 :: axis
			integer, intent(in)                   :: j
			type( ModelParameter2D_t )            :: m2D 
       end function interface_slice_2d_model_parameter
	   !
	   !
	   !
       pure function iface_SigMap( x, p_job ) result( y )
          import :: prec
          real( kind=prec ), intent(in)      :: x
          character(*), intent(in), optional :: p_job
          real( kind=prec )                  :: y
       end function iface_SigMap
       !**
       ! length
       ! This just returns total number of model parameters
       !*
       function interface_length_model_parameter( self ) result( nParam )
          import :: ModelParameter_t
          class( ModelParameter_t ), intent(in) :: self
          integer                               :: nParam
       end function interface_length_model_parameter
       !**
       ! zeros
       ! Set model parameter to zero.
       !*
       subroutine interface_zeros_model_parameter(self)
          import :: ModelParameter_t
          class( ModelParameter_t ), intent(inout) :: self
       end subroutine interface_zeros_model_parameter
       !**
       ! Copy model parameter.
       !*
       subroutine interface_copy_from_model_parameter(self, rhs)
          import :: ModelParameter_t          
          class( ModelParameter_t ), intent( inout ) :: self
          class( ModelParameter_t ), intent( in )    :: rhs
       end subroutine interface_copy_from_model_parameter
       !**
       ! SetType
       !*
       subroutine iface_SetType(self, paramType)
          import :: ModelParameter_t
          class( ModelParameter_t ), intent( inout ) :: self
          character(*), intent( in )                 :: paramType
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
          class(ModelParameter_t), intent(inout)   :: self
          class(rVector_t), allocatable :: eVec
       end function iface_PDEmapping
       !**
       !
       !*
       function iface_dPDEmapping(self, dm) result(eVec)
          import :: ModelParameter_t, rVector_t
          class(ModelParameter_t), intent(in)   :: self
          class(ModelParameter_t), intent(in)   :: dm
          class(rVector_t), allocatable :: eVec
       end function iface_dPDEmapping
       !**
       ! NOTE: For transpose (adjoint) dm is output`
       !          and eVec is input.
       !*
       function iface_dPDEmappingT(self, eVec) result(dm)
          import :: ModelParameter_t, rVector_t
          class(ModelParameter_t), intent(in)   :: self
          class(rVector_t)          , intent(in)   :: eVec
          class(ModelParameter_t), allocatable :: dm
       end function iface_dPDEmappingT
       
   end interface

contains
   !
   elemental function SigMap(self, x, job) result(y)
      class(ModelParameter_t), intent(in) :: self
      real(kind = prec), intent(in) :: x
      character(*)       , intent(in), optional :: job
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
      character(*)       , intent(in), optional :: job
      real(kind = prec) :: y
      
      y = x
   end function SigMap_Linear
   
   pure function SigMap_Log(x, p_job) result(y)
      real(kind = prec), intent(in) :: x
      character(*)       , intent(in), optional :: p_job
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
   
   !
   subroutine setMetricModelOperator( self, metric )
      !
      class( ModelParameter_t ), intent( inout )       :: self
      class( MetricElements_t ), target, intent( in )  :: metric
      !
      self%metric => metric
      !
   end subroutine setMetricModelOperator
   !
   
end module ModelParameter
