!
!> Derived class to define a SourceInteriorForce
!
module SourceInteriorForce
    !
    use Constants
    use cVector3D_SG
    use Source
    use ModelOperator
    use ModelParameter1D
    use Forward1D
    !
    type, extends( Source_t ) :: SourceInteriorForce_t
        !
        !> No derived properties
        !
        contains
            !
            procedure, public :: createE => createE_SourceInteriorForce
            !
            procedure, public :: createRHS => createRHS_SourceInteriorForce
            !
    end type SourceInteriorForce_t
    !
    interface SourceInteriorForce_t
        module procedure SourceInteriorForce_ctor
    end interface SourceInteriorForce_t
    !
contains
    !
    !> SourceInteriorForce constructor
    function SourceInteriorForce_ctor( model_operator, sigma, period, for_transpose ) result( self )
        implicit none
        !
        class( ModelOperator_t ), target, intent( in ) :: model_operator
        class( ModelParameter_t ), target, intent( in ) :: sigma
        real( kind=prec ), intent( in ) :: period
        logical, optional, intent( in ) :: for_transpose
        !
        type( SourceInteriorForce_t ) :: self
        !
        !write( *, * ) "Constructor SourceInteriorForce_t"
        !
        call self%init
        !
        self%model_operator => model_operator
        !
        self%sigma => sigma
        !
        self%period = period
        !
        self%calc_sens = .TRUE.
        !
        if( present( for_transpose ) ) then
            !
            self%for_transpose = for_transpose
        else
            self%for_transpose = .FALSE.
            !
        endif
        !
        self%non_zero_source = .TRUE.
        !
        self%non_zero_bc = .TRUE.
        !
    end function SourceInteriorForce_ctor
    !
    !> Set self%E from Forward Modeling 1D
    subroutine createE_SourceInteriorForce( self )
        implicit none
        !
        class( SourceInteriorForce_t ), intent( inout ) :: self
        !
        stop "Error: Dummy createE_SourceInteriorForce, not to be implemented"
        !
    end subroutine createE_SourceInteriorForce
    !
    !> Set RHS from self%E
    subroutine createRHS_SourceInteriorForce( self )
        implicit none
        !
        class( SourceInteriorForce_t ), intent( inout ) :: self
        !
        integer :: pol
        !
        !> RHS = E
        if( allocated( self%rhs ) ) deallocate( self%rhs )
        allocate( self%rhs, source = self%E )
        !
        do pol = 1, size( self%rhs )
            !
            if( self%for_transpose ) then
                !
                !> E = E / DIV
                call self%E( pol )%div( self%model_operator%metric%VEdge )
                !
            else
                !> RHS = E * V_E
                call self%rhs( pol )%mult( self%model_operator%metric%VEdge )
                !
            endif
            !
        enddo
        !
    end subroutine createRHS_SourceInteriorForce
    !
end module SourceInteriorForce
