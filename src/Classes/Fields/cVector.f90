module cVector
    !
    use Constants
    use rScalar
    use rVector
    !
    type, abstract :: cVector_t
        !
        logical :: is_allocated
        !
    contains
        !
        procedure, public :: init => initializeCVector
        !**
        ! Input/Output
        !*
        procedure( interface_read_c_vector ) , deferred, public :: read
        procedure( interface_write_c_vector ), deferred, public :: write
        !**
        ! boundary operations
        !*
        procedure( interface_set_all_boundary_c_vector ), deferred :: setAllboundary
        procedure( interface_set_one_boundary_c_vector ), deferred :: setOneboundary
        procedure( interface_set_all_interior_c_vector ), deferred :: setAllinterior
        procedure( interface_int_bdry_indices_c_vector ), deferred :: intBdryIndices
        procedure( interface_boundary_c_vector ), deferred         :: boundary
        procedure( interface_interior_c_vector ), deferred         :: interior
        !**
        ! Data access
        !*
        procedure( interface_length_c_vector )    , deferred :: length
        procedure( interface_get_array_c_vector ), deferred :: getArray
        procedure( interface_set_array_c_vector ), deferred :: setArray
        !**
        ! Arithmetic/algebraic operations
        procedure( interface_zeros_c_vector ), deferred, public :: zeros
        !
        procedure( interface_add1_c_vector ) , deferred, public :: add1
        generic :: add => add1
        generic :: operator(+) => add1
        !
        procedure( interface_sub1_c_vector ) , deferred, public :: sub1
        generic :: sub => sub1
        generic :: operator(-) => sub1
        !
        procedure( interface_mult1_c_vector ), deferred, public :: mult1
        procedure( interface_mult2_c_vector ), deferred, public :: mult2
        procedure( interface_mult3_c_vector ), deferred, public :: mult3
        procedure( interface_mult4_c_vector ), deferred, public :: mult4
		!
        generic :: mult => mult1, mult2, mult3, mult4
        generic :: operator(*) => mult1, mult2, mult3, mult4
        !
        procedure( interface_mults1_c_vector ), deferred, public :: mults1
        procedure( interface_mults3_c_vector ), deferred, public :: mults3
        generic :: mults => mults1, mults3
        !
        procedure( interface_div1_c_vector ) , deferred, public :: div1
        procedure( interface_div2_c_vector ) , deferred, public :: div2
        generic :: div => div1, div2
        generic :: operator(/) => div1, div2
        !
        procedure( interface_divs2_c_vector ) , deferred, public :: divs2
        generic :: divs =>  divs2
        !
        procedure( interface_linCombS_c_vector ) , deferred, public :: linCombS
        procedure( interface_scMultAddS_c_vector ) , deferred, public :: scMultAddS
        !
        procedure( interface_dot_product_c_vector ), deferred, public :: dotProd
        generic :: operator(.dot.) => dotProd
        !
        procedure( interface_copy_from_c_vector ), deferred, public :: copyFrom
        generic :: assignment(=) => copyFrom
        !
        procedure( interface_print_c_vector ), deferred, public :: print
        !**
        ! Miscellaneous
        !**
        procedure( interface_interp_func_c_vector ), deferred :: interpFunc
        !
    end type cVector_t
    !
    abstract interface
        !**
        ! Input/Output
        !*

        !**
        ! Read
        !*
        subroutine interface_read_c_vector( self, fid )
            import :: cVector_t
            class( cVector_t ), intent( inout ) :: self
            integer, intent( in )               :: fid
        end subroutine interface_read_c_vector
        !**
        ! write
        !*
        subroutine interface_write_c_vector( self, fid )
            import :: cVector_t
            class( cVector_t ), intent( in ) :: self
            integer, intent( in )            :: fid
        end subroutine interface_write_c_vector
        !**
        ! boundary operations
        !*
        !**
        ! setAllboundary
        !*
        subroutine interface_set_all_boundary_c_vector( self, c_in )
            import :: cVector_t, prec
            class( cVector_t ), intent( inout ) :: self
            complex( kind=prec ), intent( in ) :: c_in
        end subroutine interface_set_all_boundary_c_vector
        !**
        ! setOneboundary
        !*
        subroutine interface_set_one_boundary_c_vector( self, bdry, c, int_only )
            import :: cVector_t, prec
            class( cVector_t ), intent( inout ) :: self
            character(*), intent( in )          :: bdry
            complex( kind=prec ), intent( in )  :: c
            logical, optional, intent( in )     :: int_only
        end subroutine interface_set_one_boundary_c_vector
        !**
        ! setAllinterior
        !*
        subroutine interface_set_all_interior_c_vector( self, c_in )
            import :: cVector_t, prec
            class( cVector_t ), intent( inout ) :: self
            complex( kind=prec ), intent( in )  :: c_in
        end subroutine interface_set_all_interior_c_vector
        !**
        ! intBdryIndices
        !*
        subroutine interface_int_bdry_indices_c_vector( self, ind_i, ind_b )
            import :: cVector_t
            class( cVector_t ), intent( in )  :: self
            integer, allocatable, intent(out) :: ind_i(:), ind_b(:)
        end subroutine interface_int_bdry_indices_c_vector
        !**
        ! boundary
        !*
        function interface_boundary_c_vector( self ) result( E )
            import :: cVector_t
            class( cVector_t ), intent( in ) :: self
            class( cVector_t ), allocatable  :: E
        end function interface_boundary_c_vector
        !**
        ! interior
        !*
        function interface_interior_c_vector( self ) result( E )
            import :: cVector_t
            class( cVector_t ), intent( in ) :: self
            class( cVector_t ), allocatable  :: E
        end function interface_interior_c_vector
        !**
        ! Data access
        !*
        !**
        ! length
        !*
        function interface_length_c_vector( self ) result( n )
            import :: cVector_t
            class( cVector_t ), intent( in ) :: self
            integer                          :: n
        end function interface_length_c_vector
        !**
        ! getArray
        !*
        subroutine interface_get_array_c_vector( self, v )
            import :: cVector_t, prec
            class( cVector_t ), intent( in )                 :: self
            complex( kind=prec ), allocatable, intent( out ) :: v(:)
        end subroutine interface_get_array_c_vector
        !**
        ! setArray
        !*
        subroutine interface_set_array_c_vector( self, v )
            import :: cVector_t, prec
            class( cVector_t ), intent( inout ) :: self
            complex( kind=prec ), intent( in )  :: v(:)
        end subroutine interface_set_array_c_vector
        !**
        ! Arithmetic/algebraic operations
        !*
        !**
        ! zeros
        !*
        subroutine interface_zeros_c_vector( self )
            import :: cVector_t
            class( cVector_t ), intent( inout ) :: self
        end subroutine interface_zeros_c_vector
        !
        function interface_add1_c_vector( lhs, rhs ) result( Eout )
            import :: cVector_t
            class( cVector_t ), intent( in ) :: lhs, rhs
            class( cVector_t ), allocatable  :: Eout
        end function interface_add1_c_vector
        !
        function interface_sub1_c_vector( lhs, rhs ) result( Eout )
            import :: cVector_t
            class( cVector_t ), intent( in ) :: lhs, rhs
            class( cVector_t ), allocatable  :: Eout
        end function interface_sub1_c_vector
        !
        function interface_mult1_c_vector( lhs, rhs ) result( Eout )
            import :: cVector_t
            class( cVector_t ), intent( in ) :: lhs, rhs
            class( cVector_t ), allocatable  :: Eout
        end function interface_mult1_c_vector
        !
        function interface_mult2_c_vector( self, c ) result( Eout )
            import :: cVector_t, prec
            class( cVector_t ), intent( in )   :: self
            complex( kind=prec ), intent( in ) :: c
            class( cVector_t ), allocatable    :: Eout
        end function interface_mult2_c_vector
        !
        function interface_mult3_c_vector( lhs, rhs ) result( Eout )
            import :: cVector_t, rVector_t
            class( cVector_t ), intent( in ) :: lhs
            class( rVector_t ), intent( in ) :: rhs
            class( cVector_t ), allocatable  :: Eout
        end function interface_mult3_c_vector
        !
        function interface_mult4_c_vector(lhs, rhs) result(Eout)
            import :: cVector_t, rScalar_t
            class( cVector_t ), intent( in ) :: lhs
            class( rScalar_t ), intent( in ) :: rhs
            class( cVector_t ), allocatable  :: Eout
            !
        end function interface_mult4_c_vector
        subroutine interface_mults3_c_vector( lhs, rhs )
            !    subroutine version that overwrites lhs with lhs*rhs
            import :: cVector_t, rVector_t
            class( cVector_t ), intent( inout ) :: lhs
            class( rVector_t ), intent( in )      :: rhs
            class( cVector_t ), allocatable     :: Eout
        end subroutine interface_mults3_c_vector
        !
        subroutine interface_mults1_c_vector( self, c )
            import :: cVector_t, prec
            class( cVector_t ), intent( inout ) :: self
            complex( kind=prec ), intent( in )  :: c
        end subroutine interface_mults1_c_vector
        !
        function interface_div1_c_vector( lhs, rhs ) result( Eout )
            import :: cVector_t
            class( cVector_t ), intent( in ) :: lhs
            class( cVector_t ), intent( in ) :: rhs
            class( cVector_t ), allocatable  :: Eout
        end function interface_div1_c_vector
        !
        function interface_div2_c_vector( lhs, rhs ) result( Eout )
            import :: cVector_t, rVector_t
            class( cVector_t ), intent( in ) :: lhs
            class(rVector_t), intent( in )   :: rhs
            class( cVector_t ), allocatable  :: Eout
        end function interface_div2_c_vector
        !
        subroutine interface_divs2_c_vector( lhs, rhs )
            import :: cVector_t, rVector_t
            class( cVector_t ), intent( inout ) :: lhs
            class( rVector_t ), intent( in )    :: rhs
            class( cVector_t ), allocatable     :: Eout
        end subroutine interface_divs2_c_vector
        !
        subroutine interface_linCombS_c_vector( lhs, rhs, c1, c2 )
            import :: cVector_t,  prec
            class( cVector_t ), intent( inout ) :: lhs
            class( cVector_t ), intent( in )    :: rhs
            complex( kind=prec ), intent( in )  :: c1
            complex( kind=prec ), intent( in )  :: c2
        end subroutine interface_linCombS_c_vector
        !
        subroutine interface_scMultAddS_c_vector( lhs, rhs, c )
              import :: cVector_t, prec
              class( cVector_t ), intent( in )    :: lhs
              class( cVector_t ), intent( inout ) :: rhs
              complex( kind=prec ), intent( in )  :: c
        end subroutine interface_scMultAddS_c_vector
        !
        function interface_dot_product_c_vector( lhs, rhs ) result( r )
            import :: cVector_t, prec
            class( cVector_t ), intent( in ) :: lhs, rhs
            complex( kind=prec ) :: r
        end function interface_dot_product_c_vector
        !
        subroutine interface_copy_from_c_vector( self, rhs )
            import :: cVector_t
            class( cVector_t ), intent( inout ) :: self
            class( cVector_t ), intent( in )    :: rhs
        end subroutine interface_copy_from_c_vector
        !**
        ! Create a Vector object containing weights needed for
        ! interpolation of xyz component of obj1 to location.
        !*
        subroutine interface_interp_func_c_vector( self, location, xyz, E )
            import :: cVector_t, prec
            class( cVector_t ), intent( in )              :: self
            real( kind=prec ), intent( in )               :: location(3)
            character, intent( in )                       :: xyz
            class( cVector_t ) , intent(out), allocatable :: E
        end subroutine interface_interp_func_c_vector
        !
        subroutine interface_print_c_vector( self, io_unit, title, append )
            import :: cVector_t
            class( cVector_t ) , intent( in )    :: self
            integer, intent( in ), optional      :: io_unit
            character(*), intent( in ), optional :: title
            logical, intent( in ), optional      :: append
        end subroutine interface_print_c_vector
        !
    end interface
    !
contains
    !
    subroutine initializeCVector( self )
        implicit none
        !
        class( cVector_t ), intent( inout ) :: self
        !
        self%is_allocated = .FALSE.
        !
    end subroutine initializeCVector
    !
end module cVector
