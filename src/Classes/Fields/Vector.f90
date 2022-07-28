module Vector
    !
    use Constants
    use Scalar
    !
    type, abstract, extends( Field_t ) :: Vector_t
        !
        integer, dimension(3) :: NdX, NdY, NdZ, Nxyz
        !
    contains
        !
        procedure( interface_boundary_vector ), deferred :: boundary
        procedure( interface_interiovector ), deferred :: interior
        !
        procedure( interface_vector_mult_by_field ), deferred, public :: multByField
        procedure( interface_vector_mult_by_value ), deferred, public :: multByValue
        generic :: mult => multByField, multByValue
        !
        procedure( interface_vector_div_by_field ), deferred, public :: divByField
        procedure( interface_vector_div_by_value ), deferred, public :: divByValue
        generic :: div => divByField, divByValue
        !
        procedure( interface_dot_product_vector ), deferred, public :: dotProd
        generic :: operator(.dot.) => dotProd
        !
        procedure( interface_diag_mult_vector ), deferred, public :: diagMult
        !
        procedure( interface_lin_combs_func_vector ), deferred, public  :: linCombS
        procedure( interface_scmultadds_func_vector ), deferred, public :: scMultAddS
        procedure( interface_interp_func_vector ), deferred, public     :: interpFunc
        !
        procedure( interface_sum_edges_vector ), deferred, public :: sumEdges
        procedure( interface_sum_cells_vector ), deferred, public :: sumCells
        !
    end type Vector_t
    !
    abstract interface
        !
        ! Boundary operations
        subroutine interface_boundary_vector( self, boundary )
            import :: Vector_t
            class( Vector_t ), intent( in ) :: self
            class( Vector_t ), allocatable  :: boundary
        end subroutine interface_boundary_vector
        !
        subroutine interface_interiovector( self, interior )
            import :: Vector_t
            class( Vector_t ), intent( in ) :: self
            class( Vector_t ), allocatable  :: interior
        end subroutine interface_interiovector
        !
        ! Arithmetic/algebraic operations
        subroutine interface_vector_mult_by_field( self, rhs )
            import :: Field_t, Vector_t
            class( Vector_t ), intent( inout ) :: self
            class( Field_t ), intent( in )     :: rhs
        end subroutine interface_vector_mult_by_field
        !
        subroutine interface_vector_mult_by_value( self, cvalue )
            import :: Vector_t, prec
            class( Vector_t ), intent( inout ) :: self
            complex( kind=prec ), intent( in ) :: cvalue
        end subroutine interface_vector_mult_by_value
        !
        subroutine interface_vector_div_by_field( self, rhs )
            import :: Field_t, Vector_t
            class( Vector_t ), intent( inout ) :: self
            class( Field_t ), intent( in )     :: rhs
        end subroutine interface_vector_div_by_field
        !
        subroutine interface_vector_div_by_value( self, cvalue )
            import :: Vector_t, prec
            class( Vector_t ), intent( inout ) :: self
            complex( kind=prec ), intent( in ) :: cvalue
        end subroutine interface_vector_div_by_value
        !
        ! Miscellaneous
        function interface_dot_product_vector( self, rhs ) result( cvalue )
            import :: Vector_t, prec
            class( Vector_t ), intent( in ) :: self, rhs
            complex( kind=prec )            :: cvalue
        end function interface_dot_product_vector
        !
        function interface_diag_mult_vector( self, rhs ) result ( diag_mult )
            import :: Vector_t
            class( Vector_t ), intent( in ) :: self, rhs
            class( Vector_t ), allocatable  :: diag_mult
        end function interface_diag_mult_vector
        !
        subroutine interface_lin_combs_func_vector( self, rhs, c1, c2 )
            import :: Vector_t, prec
            class( Vector_t ), intent( inout ) :: self
            class( Vector_t ), intent( in )    :: rhs
            complex( kind=prec ), intent( in ) :: c1, c2
        end subroutine interface_lin_combs_func_vector
        !
        subroutine interface_scmultadds_func_vector( self, rhs, cvalue )
            import :: Vector_t, prec
            class( Vector_t ), intent( in )    :: self
            class( Vector_t ), intent( inout ) :: rhs
            complex( kind=prec ), intent( in )  :: cvalue
        end subroutine interface_scmultadds_func_vector
        !
        subroutine interface_interp_func_vector( self, location, xyz, interp )
            import :: Vector_t, prec
            class( Vector_t ), intent( in )    :: self
            real( kind=prec ), intent( in )    :: location(3)
            character, intent( in )            :: xyz
            class( Vector_t ), intent( inout ) :: interp
        end subroutine interface_interp_func_vector
        !
        function interface_sum_edges_vector( self, interior_only ) result( cell_obj )
            import :: Vector_t, Scalar_t
            class( Vector_t ), intent( in ) :: self
            logical, optional, intent( in ) :: interior_only
            class( Scalar_t ), allocatable  :: cell_obj
        end function interface_sum_edges_vector
        !
        subroutine interface_sum_cells_vector( self, E_in, ptype )
            import :: Vector_t, Scalar_t
            class( Vector_t ), intent( inout )   :: self
            class( Scalar_t ), intent( in )      :: E_in
            character(*), intent( in ), optional :: ptype
        end subroutine interface_sum_cells_vector
        !
    end interface
    !
end module Vector
