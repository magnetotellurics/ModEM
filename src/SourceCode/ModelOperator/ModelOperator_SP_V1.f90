!
!> Derived class to define a ModelOperator
!> with basic operations for Sparse Matrices - Version 1
!
module ModelOperator_SP_V1
    !
    use ModelOperator_SP
    !
    type, extends( ModelOperator_SP_t ) :: ModelOperator_SP_V1_t
        !
        ! No derived properties
        !
        contains
            !
            final :: ModelOperator_SP_V1_dtor
            !
    end type ModelOperator_SP_V1_t
    !
    interface ModelOperator_SP_V1_t
        module procedure ModelOperator_SP_V1_ctor
    end interface ModelOperator_SP_V1_t
    !
contains
    !
    !> No subroutine briefing
    !
    function ModelOperator_SP_V1_ctor( grid ) result( self )
        implicit none
        !
        class( Grid_t ), target, intent( in ) :: grid
        !
        type( ModelOperator_SP_V1_t ) :: self
        !
        call self%baseInit
        !
        call self%create( grid )
        !
    end function ModelOperator_SP_V1_ctor
    !
    !> ModelOperator_SP_V1 destructor
    subroutine ModelOperator_SP_V1_dtor( self )
        implicit none
        !
        type( ModelOperator_SP_V1_t ), intent( inout ) :: self
        !
        call self%baseDealloc
        !
        call self%dealloc
        !
    end subroutine ModelOperator_SP_V1_dtor
    !
    !> No subroutine briefing
    !
    subroutine print_ModelOperator_SP_V1( self )
        implicit none
        !
        class( ModelOperator_SP_V1_t ), intent( in ) :: self
        !
        call errStop( "print_ModelOperator_SP_V1 not implemented" )
        !
    end subroutine print_ModelOperator_SP_V1
    !
end module ModelOperator_SP_V1
