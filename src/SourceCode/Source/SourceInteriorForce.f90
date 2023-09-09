!
!> Derived class to define a SourceInteriorForce
!> Based on a pre-determined Eletric Field (E)
!
module SourceInteriorForce
    !
    use Source
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
    !
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
        call self%baseInit
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
    !> Dummy subroutine
    !> not to be implemented for this Source type
    !
    subroutine createE_SourceInteriorForce( self )
        implicit none
        !
        class( SourceInteriorForce_t ), intent( inout ) :: self
        !
        call errStop( "createE_SourceInteriorForce not to be implemented" )
        !
    end subroutine createE_SourceInteriorForce
    !
    !> Build the proper Source RHS from its E
    !
    subroutine createRHS_SourceInteriorForce( self )
        implicit none
        !
        class( SourceInteriorForce_t ), intent( inout ) :: self
        !
        integer :: pol
        type( cVector3D_MR_t ) :: temp_vec_mr
        !
        !> RHS = E
        allocate( self%rhs( size( self%E ) ) )
        !
        do pol = 1, size( self%rhs )
            !
            !> Check if grid is MR 
            !> RHS calculated as MR vector
            select type( grid => self%model_operator%metric%grid )
                !
                class is( Grid3D_SG_t )
                    !
                    allocate( self%rhs( pol )%v, source = self%E( pol ) )
                    !
                class is( Grid3D_MR_t )
                    !
                    temp_vec_mr = cVector3D_MR_t( grid, self%E( pol )%grid_type )
                    !
                    call temp_vec_mr%fromSG( self%E( pol ) )
                    !
                    allocate( self%rhs( pol )%v, source = temp_vec_mr )
                    !
                class default
                    call errStop( "createRHS_SourceInteriorForce > model_operator must be SP V1 or V2" )
                !
            end select
            !
            if( self%for_transpose ) then
                !
                !> E = E / V_E
                call self%E( pol )%div( self%model_operator%metric%v_edge )
                !
            else
                !> RHS = RHS * V_E
                call self%rhs( pol )%v%mult( self%model_operator%metric%v_edge )
                !
            endif
            !
        enddo
        !
    end subroutine createRHS_SourceInteriorForce
    !
end module SourceInteriorForce
