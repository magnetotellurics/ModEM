!
!> Abstract Base class to define a ModelOperator
!
module ModelOperator
    !
    use MetricElements
    use ModelParameter
    !
    character(:), allocatable :: model_operator_type
    character( len=15 ), parameter :: MODELOP_MF = "MatrixFreeField"
    character( len=17 ), parameter :: MODELOP_SP = "SparseMatrixField"
    character( len=19 ), parameter :: MODELOP_SP2 = "SparseMatrixFieldV2"
    !
    type, abstract :: ModelOperator_t
        !
        class( MetricElements_t ), allocatable :: metric
        !
        integer :: mKey(8)
        !
        logical :: eqset, is_allocated
        !
        contains 
            !
            procedure( interface_set_equations_model_operator ), deferred, public :: setEquations
            procedure( interface_set_cond_model_operator ), deferred, public :: setCond
            procedure( interface_amult_model_operator ), deferred, public :: amult
            procedure( interface_multaib_model_operator ), deferred, public :: multAib
            procedure( interface_multcurl_t_model_operator ), deferred, public :: multCurlT
            !
            procedure( interface_divcor_setup_model_operator ), deferred, public :: divCorSetUp
            procedure( interface_divc_grad_model_operator ), deferred, public :: divCGrad
            procedure( interface_divc_model_operator ), deferred, public :: divC
            procedure( interface_grad_model_operator ), deferred, public :: grad
            procedure( interface_div_model_operator ), deferred, public :: div
            !
            procedure( interface_print_model_operator ), deferred, public :: print
            !
            procedure, public :: baseInit => baseInit_ModelOperator
            procedure, public :: baseDealloc => baseDealloc_ModelOperator
            !
    end type ModelOperator_t
    !
    abstract interface
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_equations_model_operator( self )
            import :: ModelOperator_t
            !
            class( ModelOperator_t ), intent( inout ) :: self
        end subroutine interface_set_equations_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_cond_model_operator( self, sigma ) 
            import :: ModelOperator_t, ModelParameter_t
            !
            class( ModelOperator_t ), intent( inout ) :: self
            class( ModelParameter_t ), intent( in ) :: sigma
        end subroutine interface_set_cond_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_amult_model_operator( self, omega, in_e, out_e, p_adjoint )
            import :: ModelOperator_t, prec, Vector_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            real( kind=prec ), intent( in ), optional :: omega
            class( Vector_t ), intent( inout ) :: in_e
            class( Vector_t ), intent( inout ) :: out_e
            logical, intent( in ), optional :: p_adjoint
        end subroutine interface_amult_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_multaib_model_operator( self, in_e, out_e )
            import :: ModelOperator_t, Vector_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            class( Vector_t ), intent( inout ) :: in_e
            class( Vector_t ), intent( inout ) :: out_e
        end subroutine interface_multaib_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_multcurl_t_model_operator( self, in_e, out_e )
            import :: ModelOperator_t, Vector_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            class( Vector_t ), intent( inout ) :: in_e
            class( Vector_t ), allocatable, intent( inout ) :: out_e
        end subroutine interface_multcurl_t_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_divcor_setup_model_operator( self )
            import :: ModelOperator_t
            class( ModelOperator_t ), intent( inout ) :: self
        end subroutine interface_divcor_setup_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_divc_grad_model_operator( self, in_phi, out_phi )
            import :: ModelOperator_t, Scalar_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            class( Scalar_t ), intent( inout ) :: in_phi
            class( Scalar_t ), intent( inout ) :: out_phi
        end subroutine interface_divc_grad_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_divc_model_operator( self, in_e, out_phi )
            import :: ModelOperator_t, Vector_t, Scalar_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            class( Vector_t ), intent( inout ) :: in_e
            class( Scalar_t ), intent( inout ) :: out_phi
        end subroutine interface_divc_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_grad_model_operator( self, in_phi, out_e )
            import :: ModelOperator_t, Scalar_t, Vector_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            class( Scalar_t ), intent( inout ) :: in_phi
            class( Vector_t ), intent( inout ) :: out_e
        end subroutine interface_grad_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_div_model_operator( self, in_e, out_phi )
            import :: ModelOperator_t, Vector_t, Scalar_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            class( Vector_t ), intent( inout ) :: in_e
            class( Scalar_t ), intent( inout ) :: out_phi
        end subroutine interface_div_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_print_model_operator( self )
            import :: ModelOperator_t
            !
            class( ModelOperator_t ), intent( in ) :: self
        end subroutine interface_print_model_operator
        !
    end interface
    !
contains
    !
    !> No subroutine briefing
    subroutine baseInit_ModelOperator( self )
        implicit none
        !
        class( ModelOperator_t ), intent( inout ) :: self
        !
        self%eqset = .FALSE.
        !
        self%is_allocated = .FALSE.
        !
        call date_and_time( values=self%mKey )
        !
    end subroutine baseInit_ModelOperator
    !
    !> No subroutine briefing
    subroutine baseDealloc_ModelOperator( self )
        implicit none
        !
        class( ModelOperator_t ), intent( inout ) :: self
        !
        if( allocated( self%metric ) ) deallocate( self%metric )
        !
    end subroutine baseDealloc_ModelOperator
    !
end module ModelOperator
