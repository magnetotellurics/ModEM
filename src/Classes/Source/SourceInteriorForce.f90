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
            procedure, public :: createE => createESourceInteriorForce
            procedure, public :: createRHS => createRHSSourceInteriorForce
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
    function SourceInteriorForce_ctor( model_operator, sigma, period, trans ) result( self )
        implicit none
        !
        class( ModelOperator_t ), target, intent( in ) :: model_operator
        class( ModelParameter_t ), target, intent( in ) :: sigma
        real( kind=prec ), intent( in ) :: period
        logical, optional, intent( in ) :: trans
        !
        type( SourceInteriorForce_t ) :: self
        !
        !write( *, * ) "Constructor SourceInteriorForce_t"
        !
        call self%init()
        !
        self%model_operator => model_operator
        !
        self%sigma => sigma
        !
        self%period = period
        !
        if( present( trans ) ) then
            !
            self%trans = trans
        else
            self%trans = .FALSE.
            !
        endif
        !
        self%sens = .TRUE.
        !
        self%non_zero_source = .TRUE.
        !
    end function SourceInteriorForce_ctor
    !
    !> Set self%E from Forward Modeling 1D
    subroutine createESourceInteriorForce( self )
        implicit none
        !
        class( SourceInteriorForce_t ), intent( inout ) :: self
        !
        stop "Error: Dummy createESourceInteriorForce, not to be implemented"
        !
    end subroutine createESourceInteriorForce
    !
    !> Set RHS from self%E
    subroutine createRHSSourceInteriorForce( self )
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
        if( self%trans ) then
            !
        else
            !
            !> RHS = E * V_E
            do pol = 1, size( self%rhs )
                !
                call self%rhs( pol )%mult( self%model_operator%metric%VEdge )
                !
            enddo
            !
        endif
        !
    end subroutine createRHSSourceInteriorForce
    !
end module SourceInteriorForce
