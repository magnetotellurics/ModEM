!
!> Derived class to define a SourceAdjoint
!> Based on a pre-determined Eletric Field (E)
!
module SourceAdjoint
    !
    use Source
    !
    type, extends( Source_t ) :: SourceAdjoint_t
        !
        !> No derived properties
        !
        contains
            !
            procedure, public :: createE => createE_SourceAdjoint
            !
            procedure, public :: createRHS => createRHS_SourceAdjoint
            !
    end type SourceAdjoint_t
    !
    interface SourceAdjoint_t
        module procedure SourceAdjoint_ctor
    end interface SourceAdjoint_t
    !
contains
    !
    !> SourceAdjoint constructor
    !
    function SourceAdjoint_ctor( model_operator, sigma, period, for_transpose ) result( self )
        implicit none
        !
        class( ModelOperator_t ), target, intent( in ) :: model_operator
        class( ModelParameter_t ), target, intent( in ) :: sigma
        real( kind=prec ), intent( in ) :: period
        logical, optional, intent( in ) :: for_transpose
        !
        type( SourceAdjoint_t ) :: self
        !
        !write( *, * ) "Constructor SourceAdjoint_t"
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
    end function SourceAdjoint_ctor
    !
    !> Dummy subroutine
    !> not to be implemented for this Source type
    !
    subroutine createE_SourceAdjoint( self )
        implicit none
        !
        class( SourceAdjoint_t ), intent( inout ) :: self
        !
        call errStop( "createE_SourceAdjoint not to be implemented" )
        !
    end subroutine createE_SourceAdjoint
    !
    !> Build the proper Source RHS from its E
    !
    subroutine createRHS_SourceAdjoint( self )
        implicit none
        !
        class( SourceAdjoint_t ), intent( inout ) :: self
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
                    call errStop( "createRHS_SourceAdjoint > model_operator must be SP V1 or V2" )
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
    end subroutine createRHS_SourceAdjoint
    !
end module SourceAdjoint
