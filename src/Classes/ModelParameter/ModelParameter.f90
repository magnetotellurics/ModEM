!
!> Abstract Base class to define a ModelParameter
!
module ModelParameter
    !
    use Constants
    use Vector
    use Grid2D
    use ModelParameter1D
    use ModelParameter2D
    use Grid
    use MetricElements
    !
    !> Global name for e_solution file
    character(:), allocatable :: dsigma_file_name
    !
    type, abstract :: ModelParameter_t
        !
        class( MetricElements_t ), pointer :: metric
        !
        integer :: mKey(8)
        real( kind=prec ) :: air_cond
        character(:), allocatable :: param_type
        logical :: zero_valued, is_allocated, is_vti
        !
        procedure( interface_sigmap_model_parameter ), pointer, nopass :: SigMap_ptr
        !
        contains
            !
            procedure, public :: init => initializeModelParameter
            !
            procedure( interface_set_type_model_parameter ), deferred, public :: SetType
            !
            procedure, public :: setMetric => setMetricModelParameter
            procedure, public :: SigMap    => SigMapModelParameter
            procedure, public :: SetSigMap => SetSigMapModelParameter
            !
            !> Interfaces
            procedure( interface_zeros_model_parameter ), deferred, public :: zeros
            procedure( interface_copy_from_model_parameter ), deferred, public :: copyFrom
            !
            procedure( interface_count_model_parameter ), deferred, public :: countModel
            !
            procedure( interface_lin_comb_parameter ), deferred, public :: linComb
            !
            procedure( interface_dot_product_model_parameter ), deferred, public :: dotProd
            generic :: operator(.dot.) => dotProd
            !
            procedure( interface_pdemapping_model_parameter ), deferred, public :: PDEmapping
            procedure( interface_dpdemapping_model_parameter ), deferred, public :: dPDEmapping
            procedure( interface_dpdemapping_t_model_parameter ), deferred, public :: dPDEmappingT
            !
            procedure( interface_slice_1d_model_parameter ), deferred, public :: Slice1D
            procedure( interface_slice_2d_model_parameter ), deferred, public :: Slice2D
            !
            procedure( interface_add_model_parameter ), deferred, public :: add
            !
            procedure( interface_cond_rmsd_model_parameter ), deferred, public :: rmsd
            !
            procedure( interface_write_model_parameter ), deferred, public :: write
            procedure( interface_print_model_parameter ), deferred, public :: print
            !
    end type ModelParameter_t
    !
    abstract interface
        !
        !> No interface function briefing
        function interface_slice_1d_model_parameter( self, ix, iy ) result( model_param_1D )
            import :: ModelParameter_t, ModelParameter1D_t
            class( ModelParameter_t ), intent( in ) :: self
            integer, intent( in ) :: ix, iy
            type( ModelParameter1D_t ) ::  model_param_1D 
        end function interface_slice_1d_model_parameter
        !
        !> No interface function briefing
        function interface_avg_model_1d_model_parameter( self ) result( model_param_1D )
            import :: ModelParameter_t, ModelParameter1D_t
            class( ModelParameter_t ), intent( in ) :: self
            type( ModelParameter1D_t ) :: model_param_1D
        end function interface_avg_model_1d_model_parameter
        !
        !> No interface function briefing
        function interface_slice_2d_model_parameter( self, axis, j ) result( m2D )
            import :: ModelParameter_t, ModelParameter2D_t
            class( ModelParameter_t ), intent( in ) :: self
            integer, intent( in ) :: axis, j
            type( ModelParameter2D_t ) :: m2D 
        end function interface_slice_2d_model_parameter
        !
        !> No interface function briefing
        subroutine interface_add_model_parameter( self, other )
            import :: ModelParameter_t
            class( ModelParameter_t ), intent( inout ) :: self
            class( ModelParameter_t ), intent( in ) :: other
        end subroutine interface_add_model_parameter
        !
        !> No interface function briefing
        function interface_cond_rmsd_model_parameter( self, other ) result( rmsd )
            import :: ModelParameter_t, prec
            class( ModelParameter_t ), intent( in ) :: self, other
            complex( kind=prec ) :: rmsd
        end function interface_cond_rmsd_model_parameter
        !
        !> No interface briefing
        pure function interface_sigmap_model_parameter( x, p_job ) result( y )
            import :: prec
            real( kind=prec ), intent( in ) :: x
            character(*), intent( in ), optional :: p_job
            real( kind=prec ) :: y
        end function interface_sigmap_model_parameter
        !
        !> No interface subroutine briefing
        subroutine interface_zeros_model_parameter(self)
            import :: ModelParameter_t
            class( ModelParameter_t ), intent( inout ) :: self
        end subroutine interface_zeros_model_parameter
        !
        !> No interface subroutine briefing
        subroutine interface_copy_from_model_parameter( self, rhs )
            import :: ModelParameter_t            
            class( ModelParameter_t ), intent( inout ) :: self
            class( ModelParameter_t ), intent( in ) :: rhs
        end subroutine interface_copy_from_model_parameter
        !
        !> No interface function briefing
        function interface_count_model_parameter( self ) result( counter )
            import :: ModelParameter_t
            class( ModelParameter_t ), intent( in ) :: self
            integer :: counter
        end function interface_count_model_parameter
        !
        !> No interface subroutine briefing
        subroutine interface_lin_comb_parameter( self, a1, a2, rhs )
            import :: ModelParameter_t, prec
            class( ModelParameter_t ), intent( inout ) :: self
            real( kind=prec ), intent( in ) :: a1,a2
            class( ModelParameter_t ), intent( in ) :: rhs
        end subroutine interface_lin_comb_parameter
        !
        !> No interface function briefing
        function interface_dot_product_model_parameter( self, rhs ) result( rvalue )
            import :: ModelParameter_t, prec
            class( ModelParameter_t ), intent( in ) :: self, rhs
            real( kind=prec ) :: rvalue
        end function interface_dot_product_model_parameter
        !
        !> No interface subroutine briefing
        subroutine interface_set_type_model_parameter( self, param_type )
            import :: ModelParameter_t
            class( ModelParameter_t ), intent( inout ) :: self
            character(:), allocatable, intent( in ) :: param_type
        end subroutine interface_set_type_model_parameter
        !
        !> No interface subroutine briefing
        subroutine interface_pdemapping_model_parameter( self, eVec )
            import :: ModelParameter_t, Vector_t
            class( ModelParameter_t ), intent( in ) :: self
            class( Vector_t ), intent( inout ) :: eVec
        end subroutine interface_pdemapping_model_parameter
        !
        !> No interface subroutine briefing
        subroutine interface_dpdemapping_model_parameter( self, dsigma, eVec )
            import :: ModelParameter_t, Vector_t
            class( ModelParameter_t ), intent( in ) :: self, dsigma
            class( Vector_t ), intent( inout ) :: eVec
        end subroutine interface_dpdemapping_model_parameter
        !
        !> No interface function briefing
        subroutine interface_dpdemapping_t_model_parameter( self, eVec, dsigma )
            import :: ModelParameter_t, Vector_t
            class( ModelParameter_t ), intent( in ) :: self
            class( Vector_t ), intent( in ) :: eVec
            class( ModelParameter_t ), allocatable, intent( out ) :: dsigma
        end subroutine interface_dpdemapping_t_model_parameter
        !
        !> No interface subroutine briefing
        subroutine interface_write_model_parameter( self, comment )
            import :: ModelParameter_t, Vector_t
            class( ModelParameter_t ), intent( in ) :: self
            character(*), intent( in ), optional :: comment
        end subroutine interface_write_model_parameter
        !
        !> No interface subroutine briefing
        subroutine interface_print_model_parameter( self )
            import :: ModelParameter_t
            class( ModelParameter_t ), intent( in ) :: self
        end subroutine interface_print_model_parameter
        ! !
    end interface
    !
contains
    !
    !> No subroutine briefing
    subroutine setMetricModelParameter( self, metric )
        implicit none
        !
        class( ModelParameter_t ), intent( inout ) :: self
        class( MetricElements_t ), target, intent( in ) :: metric
        !
        self%metric => metric
        !
    end subroutine setMetricModelParameter
    !
    !> No procedure briefing
    elemental function SigMapModelParameter( self, x, job ) result( y )
        implicit none
        !
        class( ModelParameter_t), intent( in ) :: self
        real( kind=prec ), intent( in ) :: x
        character(*), intent( in ), optional :: job
        !
        real( kind=prec ) :: y
        !
        y = self%Sigmap_ptr( x )
        !
    end function SigMapModelParameter
    !
    !> No subroutine briefing
    subroutine SetSigMapModelParameter( self, param_type )
        implicit none
        !
        class( ModelParameter_t ) :: self
        character(*), intent( in ) :: param_type
        !
        select case( param_type )
            case( LOGE )
                self%SigMap_ptr => SigMap_Log
            case( LINEAR )
                self%SigMap_ptr => SigMap_Linear
        end select
        !
    end subroutine SetSigMapModelParameter
    !
    !> No procedure briefing
    pure function SigMap_Linear( x, job ) result( y )
        implicit none
        !
        real( kind=prec ), intent( in ) :: x
        character(*), intent( in ), optional :: job
        real( kind=prec ) :: y
        !
        y = x
        !
    end function SigMap_Linear
    !
    !> No procedure briefing
    pure function SigMap_Log( x, p_job ) result( y )
        implicit none
        !
        real( kind=prec ), intent( in ) :: x
        character(*), intent( in ), optional :: p_job
        real( kind=prec ) :: y
        !
        character(30) :: job
        !
        if(.NOT.present( p_job ) ) then
            job = FORWARD
        else
            job = p_job
        endif
        !
        select case( job )
           case ( FORWARD )
               y = exp( x )
           case ( DERIV )
               y = exp( x )
           case ( INVERSE )
               y = log( x )
        end select
        !
    end function SigMap_Log
    !
    !> No subroutine briefing
    subroutine initializeModelParameter( self )
        implicit none
        !
        class( ModelParameter_t ), intent( inout ) :: self
        !
        call date_and_time( values=self%mKey )
        !
        self%metric => null()
        !
        self%SigMap_ptr => null()
        !
        self%param_type   = ""
        self%air_cond     = SIGMA_AIR
        self%zero_valued  = .FALSE.
        self%is_allocated = .FALSE.
        self%is_vti       = .FALSE.
        !
    end subroutine initializeModelParameter
    !
end module ModelParameter
