!
!> Abstract base class to define a ModelOperator
!
module ModelOperator
    !
    use MetricElements
    use ModelParameter
    !
    character(:), allocatable :: model_operator_type
    character( len=11 ), parameter :: MODELOP_MF = "Matrix Free"
    character( len=13 ), parameter :: MODELOP_SP = "Sparse Matrix"
    character( len=16 ), parameter :: MODELOP_SP2 = "Sparse Matrix V2"
    !
    type, abstract :: ModelOperator_t
        !
        class( MetricElements_t ), pointer :: metric
        !
        integer :: mKey(8)
        !
        logical :: eqset, is_allocated
        !
        contains 
            !
            !> Abstract Interfaces
            !
            !> Setup
            procedure( interface_create_model_operator ), deferred, public :: create
            procedure( interface_set_equations_model_operator ), deferred, public :: setEquations
            procedure( interface_set_cond_model_operator ), deferred, public :: setCond
            !
            procedure( interface_dealloc_operator ), deferred, public :: dealloc
            !
            procedure( interface_divcor_setup_model_operator ), deferred, public :: divCorSetUp
            !
            !> Operations
            procedure( interface_amult_model_operator ), deferred, public :: amult
            procedure( interface_multaib_model_operator ), deferred, public :: multAib
            !
            procedure( interface_multcurl_t_model_operator ), deferred, public :: multCurlT
            !
            procedure( interface_div_model_operator ), deferred, public :: div
            procedure( interface_divc_model_operator ), deferred, public :: divC
            procedure( interface_divc_grad_model_operator ), deferred, public :: divCGrad
            !
            procedure( interface_grad_model_operator ), deferred, public :: grad
            !
            !> Miscellaneous
            procedure( interface_print_model_operator ), deferred, public :: print
            !
            !> Base procedures
            procedure, public :: baseInit => baseInit_ModelOperator
            procedure, public :: baseDealloc => baseDealloc_ModelOperator
            !
    end type ModelOperator_t
    !
    !> Public Global Generic ModelOperator object
    class( ModelOperator_t ), allocatable :: model_operator
    !
    public :: Maxwell
    !
    abstract interface
        !
        !> No interface subroutine briefing
        !
        subroutine interface_create_model_operator( self, grid )
            import :: ModelOperator_t, Grid_t
            !
            class( ModelOperator_t ), intent( inout ) :: self
            class( Grid_t ), target, intent( in ) :: grid
            !
        end subroutine interface_create_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_equations_model_operator( self )
            import :: ModelOperator_t
            !
            class( ModelOperator_t ), intent( inout ) :: self
            !
        end subroutine interface_set_equations_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_set_cond_model_operator( self, sigma, omega_in ) 
            import :: ModelOperator_t, ModelParameter_t, prec
            !
            class( ModelOperator_t ), intent( inout ) :: self
            class( ModelParameter_t ), intent( in ) :: sigma
            real( kind=prec ), intent( in ) :: omega_in
            !
        end subroutine interface_set_cond_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_dealloc_operator( self )
            import :: ModelOperator_t
            !
            class( ModelOperator_t ), intent( inout ) :: self
            !
        end subroutine interface_dealloc_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_divcor_setup_model_operator( self )
            import :: ModelOperator_t
            !
            class( ModelOperator_t ), intent( inout ) :: self
            !
        end subroutine interface_divcor_setup_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_amult_model_operator( self, in_e, out_e, omega, adjoint )
            import :: ModelOperator_t, prec, Vector_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            class( Vector_t ), intent( in ) :: in_e
            class( Vector_t ), intent( inout ) :: out_e
            real( kind=prec ), intent( in ) :: omega
            logical, intent( in ) :: adjoint
            !
        end subroutine interface_amult_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_multaib_model_operator( self, in_e, out_e )
            import :: ModelOperator_t, Vector_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            class( Vector_t ), intent( in ) :: in_e
            class( Vector_t ), intent( inout ) :: out_e
            !
        end subroutine interface_multaib_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_multcurl_t_model_operator( self, in_b, out_e )
            import :: ModelOperator_t, Vector_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            class( Vector_t ), intent( inout ) :: in_b
            class( Vector_t ), allocatable, intent( out ) :: out_e
            !
        end subroutine interface_multcurl_t_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_div_model_operator( self, in_e, out_phi )
            import :: ModelOperator_t, Vector_t, Scalar_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            class( Vector_t ), intent( in ) :: in_e
            class( Scalar_t ), intent( inout ) :: out_phi
            !
        end subroutine interface_div_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_divc_model_operator( self, in_e, out_phi )
            import :: ModelOperator_t, Vector_t, Scalar_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            class( Vector_t ), intent( in) :: in_e
            class( Scalar_t ), intent( inout ) :: out_phi
            !
        end subroutine interface_divc_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_divc_grad_model_operator( self, in_phi, out_phi )
            import :: ModelOperator_t, Scalar_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            class( Scalar_t ), intent( in ) :: in_phi
            class( Scalar_t ), intent( inout ) :: out_phi
            !
        end subroutine interface_divc_grad_model_operator
        !
        !> No interface subroutine briefing
        !
        subroutine interface_grad_model_operator( self, in_phi, out_e )
            import :: ModelOperator_t, Scalar_t, Vector_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            class( Scalar_t ), intent( in ) :: in_phi
            class( Vector_t ), intent( inout ) :: out_e
            !
        end subroutine interface_grad_model_operator
        !!
        !> No interface subroutine briefing
        !
        subroutine interface_print_model_operator( self )
            import :: ModelOperator_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            !
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
        if( associated( self%metric ) ) deallocate( self%metric )
        !
    end subroutine baseDealloc_ModelOperator
    !
end module ModelOperator
!