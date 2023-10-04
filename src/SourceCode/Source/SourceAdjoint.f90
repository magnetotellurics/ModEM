!
!> Derived class to define a SourceAdjoint
!> Based on a pre-determined Eletric Field (E)
!
module SourceAdjoint
    !
    use Source
    use ModelOperator_SP_V2
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
        complex( kind=prec ) :: iOmegaMuInv
        complex( kind=prec ), allocatable, dimension(:) :: v_edge_v, rhs_v, rhs_v_int, in_e, in_e_int, out_e
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
                !> Check if model_operator is SP2
                select type( model_operator => self%model_operator )
                    !
                    class is( ModelOperator_SP_V2_t )
                        !
                        iOmegaMuInv = isign / ( cmplx( 0.0, real( 2.0 * ( PI / self%period ) * MU_0, kind=prec ), kind=prec ) )
                        !
                        in_e = self%E(1)%getArray()
                        !
                        v_edge_v = self%model_operator%metric%v_edge%getArray()
                        !
                        in_e_int = in_e( self%E(1)%indInterior() ) / v_edge_v( self%model_operator%metric%v_edge%indInterior() )
                        !
                        out_e = in_e_int
                        out_e = C_ZERO
                        !
                        call RMATxCVEC( model_operator%GDii, in_e_int, out_e )
                        !
                        out_e = out_e * iOmegaMuInv * v_edge_v( self%model_operator%metric%v_edge%indInterior() )
                        !
                        rhs_v = self%rhs(1)%v%getArray()
                        !
                        rhs_v_int = ( rhs_v( self%rhs(1)%v%indInterior() ) + out_e )
                        !
                        rhs_v = C_ZERO
                        rhs_v( self%rhs(1)%v%indInterior() ) = rhs_v_int
                        !
                        call self%rhs(1)%v%setArray( rhs_v )
                        !
                end select
                !
                !> E = E / V_E
                call self%E( pol )%div( self%model_operator%metric%v_edge )
                !
            else
                !
                !> Check if model_operator is SP2
                select type( model_operator => self%model_operator )
                    !
                    class is( ModelOperator_SP_V2_t )
                        !
                        !RHS = RHS + (i omega mu)^-1 GDii* E(EDGEi)
                        iOmegaMuInv = isign / ( cmplx( 0.0, real( 2.0 * ( PI / self%period ) * MU_0, kind=prec ), kind=prec ) )
                        !
                        in_e = self%E(1)%getArray()
                        !
                        in_e_int = in_e( self%E(1)%indInterior() )
                        !
                        out_e = in_e_int
                        out_e = C_ZERO
                        !
                        call RMATxCVEC( model_operator%GDii, in_e_int, out_e )
                        !
                        v_edge_v = self%model_operator%metric%v_edge%getArray()
                        !
                        out_e = out_e * iOmegaMuInv
                        !
                        rhs_v = self%rhs(1)%v%getArray()
                        !
                        rhs_v_int = ( rhs_v( self%rhs(1)%v%indInterior() ) + out_e )
                        !
                        rhs_v = C_ZERO
                        rhs_v( self%rhs(1)%v%indInterior() ) = rhs_v_int
                        !
                        call self%rhs(1)%v%setArray( rhs_v )
                        !
                end select
                !
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
