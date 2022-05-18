module ModelOperator
    !
    use Constants
    use cVector
    use cScalar
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
            !
            !procedure( iface_UpdateFrequency ), deferred, public :: UpdateFrequency
            procedure( interface_set_equations_model_operator ), deferred, public :: setEquations
            procedure( interface_set_cond_model_operator ), deferred, public        :: setCond
            procedure( interface_amult_model_operator ), deferred, public            :: amult
            procedure( interface_multaib_model_operator ), deferred, public         :: multAib
            procedure( interface_multcurl_t_model_operator ), deferred, public     :: multCurlT
            !     following procedures are generally used for divergence correction
            !     and might in some cases (e.g., model operator for "SP2" case)
            !     only be implemented as dummy procedures
            procedure( interface_divcor_setup_model_operator ), deferred, public :: divCorSetup
            procedure( interface_divc_grad_model_operator ), deferred, public :: divCgrad
            procedure( interface_divc_model_operator ), deferred, public        :: divC
            procedure( interface_grad_model_operator ), deferred, public        :: grad
            procedure( interface_div_model_operator ), deferred, public         :: div
            !     these will be coded to return cScalar/cVector of type appropriate
            !         for specific ModelOperator implementation
            !        I see no need for real versions -- but we can add if needed!
            !
            !procedure( interface_create_scalar_model_operator ), deferred, public :: createScalar
            !procedure( interface_create_vector_model_operator ), deferred, public :: createVector
            !
            procedure( interface_print_model_operator ), deferred, public :: print
            !
    end type ModelOperator_t
    
    abstract interface
         !**
         ! setEquations
         !*
         subroutine interface_set_equations_model_operator( self )
             import :: ModelOperator_t
             !
             class( ModelOperator_t ), intent( inout ) :: self
         end subroutine interface_set_equations_model_operator
         !**
         ! setCond
         subroutine interface_set_cond_model_operator( self, ModPar ) 
             import :: ModelOperator_t, ModelParameter_t
             !
             class( ModelOperator_t ), intent( inout )  :: self
             class( ModelParameter_t ), intent( in ) :: ModPar
         end subroutine interface_set_cond_model_operator
         !**
         ! multAib
         !*
         subroutine interface_multaib_model_operator( self, bdry, outE )
             import :: ModelOperator_t, cVector_t
             !
             class( ModelOperator_t ), intent( in ) :: self
             class( cVector_t ), intent( in )         :: bdry
             class( cVector_t ), intent( inout )     :: outE
         end subroutine interface_multaib_model_operator
         !**
         ! multCurlT
         !*
         subroutine interface_multcurl_t_model_operator( self, inH, outE )
             import :: ModelOperator_t, cVector_t
             !
             class( ModelOperator_t ) , intent( in )             :: self
             class( cVector_t )             , intent( inout )     :: inH
             class( cVector_t ), allocatable, intent( inout ) :: outE
         end subroutine interface_multcurl_t_model_operator
         !
         subroutine interface_amult_model_operator( self, omega, x, y, p_adjt )
             import :: ModelOperator_t, cVector_t, prec
             !
             class( ModelOperator_t ), intent( in )         :: self
             real( kind = prec ), intent( in ), optional  :: omega
             class( cVector_t ), intent( in )                 :: x
             class( cVector_t ), intent( inout )             :: y
             logical, intent( in ), optional                  :: p_adjt
         end subroutine interface_amult_model_operator
         !
         ! these are for divergence correction, might be dummies in some cases
         subroutine interface_divcor_setup_model_operator( self )
             import :: ModelOperator_t
             class( ModelOperator_t ) , intent( inout )         :: self
         end subroutine interface_divcor_setup_model_operator
         !
         subroutine interface_divc_grad_model_operator( self, inPhi, outPhi )
             import :: ModelOperator_t, cScalar_t
             !
             class( ModelOperator_t ) , intent( in )         :: self
             class( cScalar_t )             , intent( in )     :: inPhi
             class( cScalar_t )             , intent( inout ) :: outPhi
         end subroutine interface_divc_grad_model_operator
         !
         subroutine interface_divc_model_operator( self, inE, outPhi )
             import :: ModelOperator_t, cVector_t, cScalar_t
             !
             class( ModelOperator_t ) , intent( in ) :: self
             class( cVector_t )             , intent( in ) :: inE
             class( cScalar_t )             , intent( inout ) :: outPhi
         end subroutine interface_divc_model_operator
         !
         subroutine interface_grad_model_operator( self, inPhi, outE )
             import :: ModelOperator_t, cVector_t, cScalar_t
             !
             class( ModelOperator_t ), intent( in ) :: self
             class( cScalar_t ), intent( in )         :: inPhi
             class( cVector_t ), intent( inout )     :: outE
         end subroutine interface_grad_model_operator
         !
         subroutine interface_div_model_operator(self, inE, outPhi)
             import :: ModelOperator_t, cVector_t, cScalar_t
             !
             class( ModelOperator_t ), intent( in ) :: self
             class( cVector_t ), intent( in )         :: inE
             class( cScalar_t ), intent( inout )     :: outPhi
         end subroutine interface_div_model_operator
         !
         function interface_create_vector_model_operator( self, gridType ) result( cVec )
             import :: ModelOperator_t, cVector_t
             !
             class( ModelOperator_t ) , intent( in )     :: self
             character( len=80 ), intent( in ), optional :: gridType
             class( cVector_t ), allocatable             :: cVec
         end function interface_create_vector_model_operator
         !
         function interface_create_scalar_model_operator( self, gridType ) result( cSclr )
             import :: ModelOperator_t, cScalar_t
             !
             class( ModelOperator_t ), intent( in )      :: self
             character( len=80 ), intent( in ), optional :: gridType
             class( cScalar_t ), allocatable             :: cSclr
         end function interface_create_scalar_model_operator
         !
        subroutine interface_print_model_operator( self )
             import :: ModelOperator_t
             !
             class( ModelOperator_t ), intent( in ) :: self
         end subroutine interface_print_model_operator
    end interface
    !
contains
    !
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
