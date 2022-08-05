module Scalar
    !
    use Constants
    use Field
    !
    type, abstract, extends( Field_t ) :: Scalar_t
        !
        integer, dimension(3) :: NdV
        !
        integer :: Nxyz
        !
    contains
        !
        procedure( interface_scalar_mult_by_field ), deferred, public :: multByField
        procedure( interface_scalar_mult_by_value ), deferred, public :: multByValue
        generic :: mult => multByField, multByValue
        !
        procedure( interface_scalar_div_by_field ), deferred, public :: divByField
        procedure( interface_scalar_div_by_value ), deferred, public :: divByValue
        generic :: div => divByField, divByValue
        !
        procedure( interface_lin_combs_scalar ), deferred, public    :: linCombS
        !
        procedure( interface_sc_mult_adds_scalar ), deferred, public :: scMultAddS
        !
        procedure( interface_dot_product_scalar ), deferred, public :: dotProd
        generic :: operator(.dot.) => dotProd
        !
    end type Scalar_t
    !
    abstract interface
        !
        ! Arithmetic/algebraic operations
        subroutine interface_scalar_mult_by_field( self, rhs )
            import :: Field_t, Scalar_t
            class( Scalar_t ), intent( inout ) :: self
            class( Field_t ), intent( in )     :: rhs
        end subroutine interface_scalar_mult_by_field
        !
        subroutine interface_scalar_mult_by_value( self, cvalue )
            import :: Scalar_t, prec
            class( Scalar_t ), intent( inout )  :: self
            complex( kind=prec ), intent( in )  :: cvalue
        end subroutine interface_scalar_mult_by_value
        !
        subroutine interface_scalar_div_by_field( self, rhs )
            import :: Field_t, Scalar_t
            class( Scalar_t ), intent( inout ) :: self
            class( Field_t ), intent( in )     :: rhs
        end subroutine interface_scalar_div_by_field
        !
        subroutine interface_scalar_div_by_value( self, cvalue )
            import :: Scalar_t, prec
            class( Scalar_t ), intent( inout ) :: self
            complex( kind=prec ), intent( in ) :: cvalue
        end subroutine interface_scalar_div_by_value
        !
        subroutine interface_lin_combs_scalar( self, rhs, c1, c2 )
            import :: Scalar_t, prec
            class( Scalar_t ), intent( inout )  :: self
            class( Scalar_t ), intent( in )     :: rhs
            complex( kind=prec ), intent( in ) :: c1, c2
        end subroutine interface_lin_combs_scalar
        !
        subroutine interface_sc_mult_adds_scalar( self, rhs, cvalue )
            import :: Scalar_t, prec
            class( Scalar_t ), intent( in )    :: self
            class( Scalar_t ), intent( inout ) :: rhs
            complex( kind=prec ), intent( in ) :: cvalue
        end subroutine interface_sc_mult_adds_scalar
        !
        function interface_dot_product_scalar( self, rhs ) result( cvalue )
            import :: Scalar_t, prec
            class( Scalar_t ), intent( in ) :: self, rhs
            complex( kind=prec )            :: cvalue
        end function interface_dot_product_scalar
        !
    end interface
    !
end module Scalar
