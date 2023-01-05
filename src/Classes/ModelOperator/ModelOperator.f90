!
!> Abstract Base class to define a ModelOperator
!
module ModelOperator
    !
    use Constants
    use MetricElements
    use ModelParameter
    !
    type, abstract :: ModelOperator_t
        !
        class( MetricElements_t ), allocatable :: metric
        !
        logical :: is_allocated
        !
        contains 
            !
            procedure, public :: dealloc => deallocateModelOperator
            procedure, public :: init    => initializeModelOperator
            !
            procedure( interface_set_equations_model_operator ), deferred, public :: setEquations
            procedure( interface_set_cond_model_operator ), deferred, public :: setCond
            procedure( interface_amult_model_operator ), deferred, public :: amult
            procedure( interface_multaib_model_operator ), deferred, public :: multAib
            procedure( interface_multcurl_t_model_operator ), deferred, public :: multCurlT
            !
            procedure( interface_divcor_setup_model_operator ), deferred, public :: divCorSetup
            procedure( interface_divc_grad_model_operator ), deferred, public :: divCgrad
            procedure( interface_divc_model_operator ), deferred, public :: divC
            procedure( interface_grad_model_operator ), deferred, public :: grad
            procedure( interface_div_model_operator ), deferred, public :: div
			!
            procedure( interface_adj_bc_model_operator ), deferred, public :: AdjtBC
            
            !
            procedure( interface_print_model_operator ), deferred, public :: print
            !
    end type ModelOperator_t
    !
    abstract interface
        !
        !> No interface subroutine briefing
        subroutine interface_set_equations_model_operator( self )
            import :: ModelOperator_t
            !
            class( ModelOperator_t ), intent( inout ) :: self
        end subroutine interface_set_equations_model_operator
        !
        !> No interface subroutine briefing
        subroutine interface_set_cond_model_operator( self, sigma ) 
            import :: ModelOperator_t, ModelParameter_t
            !
            class( ModelOperator_t ), intent( inout ) :: self
            class( ModelParameter_t ), intent( in ) :: sigma
        end subroutine interface_set_cond_model_operator
        !
        !> No interface subroutine briefing
        subroutine interface_multaib_model_operator( self, bdry, outE )
            import :: ModelOperator_t, Vector_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            class( Vector_t ), intent( in ) :: bdry
            class( Vector_t ), intent( inout ) :: outE
        end subroutine interface_multaib_model_operator
        !
        !> No interface subroutine briefing
        subroutine interface_multcurl_t_model_operator( self, inH, outE )
            import :: ModelOperator_t, Vector_t
            !
            class( ModelOperator_t ) , intent( in ) :: self
            class( Vector_t )            , intent( inout ) :: inH
            class( Vector_t ), allocatable, intent( inout ) :: outE
        end subroutine interface_multcurl_t_model_operator
        !
        !> No interface subroutine briefing
        subroutine interface_amult_model_operator( self, omega, inE, outE, p_adjoint )
            import :: ModelOperator_t, Vector_t, prec
            !
            class( ModelOperator_t ), intent( in ) :: self
            real( kind = prec ), intent( in ), optional :: omega
            class( Vector_t ), intent( in ) :: inE
            class( Vector_t ), intent( inout ) :: outE
            logical, intent( in ), optional :: p_adjoint
        end subroutine interface_amult_model_operator
        !
        !> No interface subroutine briefing
        subroutine interface_divcor_setup_model_operator( self )
            import :: ModelOperator_t
            class( ModelOperator_t ) , intent( inout ) :: self
        end subroutine interface_divcor_setup_model_operator
        !
        !> No interface subroutine briefing
        subroutine interface_divc_grad_model_operator( self, inPhi, outPhi )
            import :: ModelOperator_t, Scalar_t
            !
            class( ModelOperator_t ) , intent( in ) :: self
            class( Scalar_t )            , intent( in ) :: inPhi
            class( Scalar_t )            , intent( inout ) :: outPhi
        end subroutine interface_divc_grad_model_operator
        !
        !> No interface subroutine briefing
        subroutine interface_divc_model_operator( self, inE, outPhi )
            import :: ModelOperator_t, Vector_t, Scalar_t
            !
            class( ModelOperator_t ) , intent( in ) :: self
            class( Vector_t )            , intent( in ) :: inE
            class( Scalar_t )            , intent( inout ) :: outPhi
        end subroutine interface_divc_model_operator
        !
        !> No interface subroutine briefing
        subroutine interface_grad_model_operator( self, inPhi, outE )
            import :: ModelOperator_t, Vector_t, Scalar_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            class( Scalar_t ), intent( in ) :: inPhi
            class( Vector_t ), intent( inout ) :: outE
        end subroutine interface_grad_model_operator
        !
        !> No interface subroutine briefing
        subroutine interface_div_model_operator( self, inE, outPhi )
            import :: ModelOperator_t, Vector_t, Scalar_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            class( Vector_t ), intent( in ) :: inE
            class( Scalar_t ), intent( inout ) :: outPhi
        end subroutine interface_div_model_operator
        !
        !> No interface subroutine briefing
        subroutine interface_adj_bc_model_operator( self, eIn, BC )
            import :: ModelOperator_t, Vector_t
            !
            class( ModelOperator_t ), intent( in ) :: self
            class( Vector_t ), intent( in ) :: eIn
            class( Vector_t ), intent( inout ) :: BC
        end subroutine interface_adj_bc_model_operator
        !
        !
        !> No interface subroutine briefing
        subroutine interface_print_model_operator( self )
            import :: ModelOperator_t
            !
            class( ModelOperator_t ), intent( in ) :: self
        end subroutine interface_print_model_operator
    end interface
    !
contains
    !
    !> No subroutine briefing
    subroutine initializeModelOperator( self )
        implicit none
        !
        class( ModelOperator_t ), intent( inout ) :: self
        !
        self%is_allocated = .FALSE.
        !
    end subroutine initializeModelOperator
    !
    !> No subroutine briefing
    subroutine deallocateModelOperator( self )
        implicit none
        !
        class( ModelOperator_t ), intent( inout ) :: self
        !
        if( allocated( self%metric ) ) deallocate( self%metric )
        !
    end subroutine deallocateModelOperator
    !
end module ModelOperator
