!
!> Abstract Base class to define a ModelParameter
!
module ModelParameter
    !
    use Utilities
    use Vector
    use Grid2D
    use ModelParameter1D
    use ModelParameter2D
    use MetricElements
    use Scalar
    !
    type, abstract :: ModelParameter_t
        !
        class( MetricElements_t ), pointer :: metric
        !
        integer :: anisotropic_level, mKey(8)
        !
        real( kind=prec ) :: air_cond
        !
        character(:), allocatable :: param_type
        !
        logical :: is_allocated
        !
        procedure( interface_sigmap_model_parameter ), pointer, nopass :: sigMap_ptr
        !
        contains
            !
            procedure, public :: baseInit => initialize_ModelParameter
            !
            procedure, public :: setMetric => setMetric_ModelParameter
            procedure, public :: sigMap => sigMap_ModelParameter
            procedure, public :: setSigMap => setSigMap_ModelParameter
            !
            !> Interfaces
            procedure( interface_set_type_model_parameter ), deferred, public :: setType
            !
            procedure( interface_get_one_cond_model_parameter ), deferred, public :: getOneCond
            procedure( interface_get_all_cond_model_parameter ), deferred, public :: getAllCond
            generic :: getCond => getOneCond, getAllCond
            !
            procedure( interface_set_one_cond_model_parameter ), deferred, public :: setOneCond
            procedure( interface_set_all_cond_model_parameter ), deferred, public :: setAllCond
            generic :: setCond => setOneCond, setAllCond
            !
            procedure( interface_zeros_model_parameter ), deferred, public :: zeros
            !
            procedure( interface_copy_from_model_parameter ), deferred, public :: copyFrom
            generic :: assignment(=) => copyFrom
            !
            procedure( interface_count_model_parameter ), deferred, public :: countModel
            !
            procedure( interface_lin_comb_model_model_parameter ), deferred, public :: linComb
            !
            procedure( interface_dot_product_model_parameter ), deferred, public :: dotProd
            !
            procedure( interface_pdemapping_model_parameter ), deferred, public :: PDEmapping
            procedure( interface_dpdemapping_model_parameter ), deferred, public :: dPDEmapping
            procedure( interface_dpdemapping_t_model_parameter ), deferred, public :: dPDEmapping_T
            !
            procedure( interface_slice_1d_model_parameter ), deferred, public :: slice1D
            procedure( interface_slice_2d_model_parameter ), deferred, public :: slice2D
            !
            procedure( interface_write_model_parameter ), deferred, public :: write
            procedure( interface_print_model_parameter ), deferred, public :: print
            !
    end type ModelParameter_t
    !
    abstract interface
        !
        !> No interface function briefing
        !
        function interface_slice_1d_model_parameter( self, ix, iy ) result( model_param_1D )
            import :: ModelParameter_t, ModelParameter1D_t
            !
            class( ModelParameter_t ), intent( inout ) :: self
            integer, intent( in ) :: ix, iy
            !
            type( ModelParameter1D_t ) ::  model_param_1D 
            !
        end function interface_slice_1d_model_parameter
        !
        !> No interface function briefing
        !
        function interface_avg_model_1d_model_parameter( self ) result( model_param_1D )
            import :: ModelParameter_t, ModelParameter1D_t
            !
            class( ModelParameter_t ), intent( inout ) :: self
            !
            type( ModelParameter1D_t ) :: model_param_1D
            !
        end function interface_avg_model_1d_model_parameter
        !
        !> No interface function briefing
        !
        function interface_slice_2d_model_parameter( self, axis, j ) result( m2D )
            import :: ModelParameter_t, ModelParameter2D_t
            !
            class( ModelParameter_t ), intent( inout ) :: self
            integer, intent( in ) :: axis, j
            !
            type( ModelParameter2D_t ) :: m2D
            !
        end function interface_slice_2d_model_parameter
        !
        !> No interface briefing
        !
        pure function interface_sigmap_model_parameter( x, p_job ) result( y )
            import :: prec
            !
            real( kind=prec ), intent( in ) :: x
            character(*), intent( in ), optional :: p_job
            !
            real( kind=prec ) :: y
            !
        end function interface_sigmap_model_parameter
        !
        !> No interface subroutine briefing
        !
        function interface_get_one_cond_model_parameter( self, i_cond ) result( cell_cond )
            import :: ModelParameter_t, Scalar_t
            !
            class( ModelParameter_t ), intent( in ) :: self
            integer, intent( in ) :: i_cond
            !
            class( Scalar_t ), allocatable :: cell_cond
            !
        end function interface_get_one_cond_model_parameter
        !
        !> No interface subroutine briefing
        !
        function interface_get_all_cond_model_parameter( self ) result( cell_cond )
            import :: ModelParameter_t, GenScalar_t
            !
            class( ModelParameter_t ), intent( in ) :: self
            !
            type( GenScalar_t ), allocatable, dimension(:) :: cell_cond
            !
        end function interface_get_all_cond_model_parameter
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_one_cond_model_parameter( self, cell_cond, i_cond )
            import :: ModelParameter_t, Scalar_t
            !
            class( ModelParameter_t ), intent( inout ) :: self
            class( Scalar_t ), intent( in ) :: cell_cond
            integer, intent( in ) :: i_cond
            !
        end subroutine interface_set_one_cond_model_parameter
        !
        !> No interface subroutine briefing
        subroutine interface_set_all_cond_model_parameter( self, cell_cond )
            import :: ModelParameter_t, GenScalar_t
            !
            class( ModelParameter_t ), intent( inout ) :: self
            type( GenScalar_t ), allocatable, dimension(:), intent( in ) :: cell_cond
            !
        end subroutine interface_set_all_cond_model_parameter
        !
        !> No interface subroutine briefing
        subroutine interface_zeros_model_parameter( self )
            import :: ModelParameter_t
            !
            class( ModelParameter_t ), intent( inout ) :: self
            !
        end subroutine interface_zeros_model_parameter
        !
        !> No interface subroutine briefing
        subroutine interface_copy_from_model_parameter( self, rhs )
            import :: ModelParameter_t
            !
            class( ModelParameter_t ), intent( inout ) :: self
            class( ModelParameter_t ), intent( in ) :: rhs
            !
        end subroutine interface_copy_from_model_parameter
        !
        !> No interface function briefing
        function interface_count_model_parameter( self ) result( counter )
            import :: ModelParameter_t
            !
            class( ModelParameter_t ), intent( in ) :: self
            !
            integer :: counter
            !
        end function interface_count_model_parameter
        !
        !> No interface subroutine briefing
        subroutine interface_lin_comb_model_model_parameter( self, a1, a2, rhs )
            import :: ModelParameter_t, prec
            !
            class( ModelParameter_t ), intent( inout ) :: self
            real( kind=prec ), intent( in ) :: a1, a2
            class( ModelParameter_t ), intent( inout ) :: rhs
            !
        end subroutine interface_lin_comb_model_model_parameter
        !
        !> No interface subroutine briefing
        !
        subroutine interface_lin_comb_scalar_model_parameter( self, a1, a2, rhs )
            import :: ModelParameter_t, prec, Scalar_t
            !
            class( ModelParameter_t ), intent( inout ) :: self
            real( kind=prec ), intent( in ) :: a1, a2
            class( Scalar_t ), intent( in ) :: rhs
            !
        end subroutine interface_lin_comb_scalar_model_parameter
        !
        !> No interface function briefing
        !
        function interface_dot_product_model_parameter( self, rhs ) result( rvalue )
            import :: ModelParameter_t, prec
            !
            class( ModelParameter_t ), intent( inout ) :: self, rhs
            real( kind=prec ) :: rvalue
            !
        end function interface_dot_product_model_parameter
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_type_model_parameter( self, param_type )
            import :: ModelParameter_t
            !
            class( ModelParameter_t ), intent( inout ) :: self
            character(:), allocatable, intent( in ) :: param_type
            !
        end subroutine interface_set_type_model_parameter
        !
        !> No interface subroutine briefing
        !
        subroutine interface_pdemapping_model_parameter( self, e_vec )
            import :: ModelParameter_t, Vector_t
            !
            class( ModelParameter_t ), intent( in ) :: self
            class( Vector_t ), intent( inout ) :: e_vec
            !
        end subroutine interface_pdemapping_model_parameter
        !
        !> No interface subroutine briefing
        !
        subroutine interface_dpdemapping_model_parameter( self, dsigma, e_vec )
            import :: ModelParameter_t, Vector_t
            !
            class( ModelParameter_t ), intent( inout ) :: self
            class( ModelParameter_t ), intent( in ) :: dsigma
            class( Vector_t ), intent( inout ) :: e_vec
            !
        end subroutine interface_dpdemapping_model_parameter
        !
        !> No interface function briefing
        !
        subroutine interface_dpdemapping_t_model_parameter( self, e_vec, dsigma )
            import :: ModelParameter_t, Vector_t
            !
            class( ModelParameter_t ), intent( inout ) :: self
            class( Vector_t ), intent( in ) :: e_vec
            class( ModelParameter_t ), allocatable, intent( out ) :: dsigma
            !
        end subroutine interface_dpdemapping_t_model_parameter
        !
        !> No interface subroutine briefing
        !
        subroutine interface_write_model_parameter( self, file_name, comment )
            import :: ModelParameter_t, Vector_t
            !
            class( ModelParameter_t ), intent( inout ) :: self
            character(*), intent( in ) :: file_name
            character(*), intent( in ), optional :: comment
            !
        end subroutine interface_write_model_parameter
        !
        !> No interface subroutine briefing
        !
        subroutine interface_print_model_parameter( self )
            import :: ModelParameter_t
            !
            class( ModelParameter_t ), intent( inout ) :: self
            !
        end subroutine interface_print_model_parameter
        ! !
    end interface
    !
contains
    !
    !> No subroutine briefing
    !
    subroutine setMetric_ModelParameter( self, metric )
        implicit none
        !
        class( ModelParameter_t ), intent( inout ) :: self
        class( MetricElements_t ), target, intent( in ) :: metric
        !
        self%metric => metric
        !
    end subroutine setMetric_ModelParameter
    !
    !> No procedure briefing
    !
    elemental function sigMap_ModelParameter( self, x, job ) result( y )
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
    end function sigMap_ModelParameter
    !
    !> No subroutine briefing
    !
    subroutine setSigMap_ModelParameter( self, param_type )
        implicit none
        !
        class( ModelParameter_t ) :: self
        character(*), intent( in ) :: param_type
        !
        select case( param_type )
            case( LOGE )
                self%sigMap_ptr => sigMap_Log
            case( LINEAR )
                self%sigMap_ptr => sigMap_Linear
        end select
        !
    end subroutine setSigMap_ModelParameter
    !
    !> No procedure briefing
    !
    pure function sigMap_Linear( x, job ) result( y )
        implicit none
        !
        real( kind=prec ), intent( in ) :: x
        character(*), intent( in ), optional :: job
        real( kind=prec ) :: y
        !
        y = x
        !
    end function sigMap_Linear
    !
    !> No procedure briefing
    !
    pure function sigMap_Log( x, p_job ) result( y )
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
           case( FORWARD )
               y = exp( x )
           case( DERIV )
               y = exp( x )
           case( INVERSE )
               y = log( x )
        end select
        !
    end function sigMap_Log
    !
    !> No subroutine briefing
    !
    subroutine initialize_ModelParameter( self )
        implicit none
        !
        class( ModelParameter_t ), intent( inout ) :: self
        !
        self%metric => null()
        !
        call date_and_time( values=self%mKey )
        !
        self%air_cond = SIGMA_AIR
        !
        self%param_type = ""
        !
        self%is_allocated = .FALSE.
        !
        self%sigMap_ptr => null()
        !
    end subroutine initialize_ModelParameter

end module ModelParameter
