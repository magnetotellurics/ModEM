module cScalar
    !
    use Constants
    use Grid
    use rScalar 
    !
    type, abstract :: cScalar_t
        !
        logical :: is_allocated
        !
     contains
        !
        procedure, public :: init => initializeCScalar
        !
        !**
        ! Input/Output
        !*
        procedure( interface_read_c_scalar ) , deferred, public :: read
        procedure( interface_write_c_scalar ), deferred, public :: write
        !**
        ! Boundary operations
        !*
        procedure( interface_setAllBoundary_c_scalar ), deferred, public :: setAllBoundary
        procedure( interface_setOneBoundary_c_scalar ), deferred, public :: setOneBoundary
        procedure( interface_setAllInterior_c_scalar ), deferred, public :: setAllInterior
        procedure( interface_intBdryIndices_c_scalar ), deferred, public :: intBdryIndices
        !**
        ! Data access
        !*
        procedure( interface_length_c_scalar ), deferred, public     :: length
        procedure( interface_get_array_c_scalar ), deferred, public :: getArray
        procedure( interface_set_array_c_scalar ), deferred, public :: setArray
        
        !**
        ! Arithmetic/algebraic operations
        !*
        procedure( interface_zeros_c_scalar ), deferred, public :: zeros
        procedure( interface_add1_c_scalar ), deferred, public :: add1
        generic :: add => add1
        generic :: operator(+) => add1
        
        procedure( interface_sub1_c_scalar ), deferred, public :: sub1
        generic :: sub => sub1
        generic :: operator(-) => sub1
        
        procedure( interface_mult1_c_scalar ), deferred, public :: mult1
        procedure( interface_mult2_c_scalar ), deferred, public :: mult2
        procedure( interface_mult3_c_scalar ), deferred, public :: mult3
        generic :: mult => mult1, mult2, mult3
        generic :: operator(*) => mult1, mult2, mult3

        procedure( interface_mults1_c_scalar ), deferred, public :: mults1
        procedure( interface_mults2_c_scalar ), deferred, public :: mults2
        procedure( interface_mults3_c_scalar ), deferred, public :: mults3
        generic :: mults => mults1, mults2, mults3 

        procedure( interface_div1_c_scalar ), deferred, public :: div1
        generic :: div => div1
        generic :: operator(/) => div1

        procedure( interface_div3S_c_scalar ), deferred, public :: divs3
        generic :: divS => divs3

        procedure( interface_linCombS_c_scalar), deferred, public :: linCombS
        procedure( interface_scMultAddS_c_scalar), deferred, public :: scMultAddS
        
        procedure( interface_dot_product_c_scalar ), deferred, public :: dotProd
        generic :: operator(.dot.) => dotProd

        !**
        ! Miscellaneous
        !*
        procedure( interface_is_compatible1_c_scalar ), deferred, public :: isCompatible1
        procedure( interface_is_compatible2_c_scalar ), deferred, public :: isCompatible2
        procedure( interface_copy_from_c_scalar ), deferred, public :: copyFrom
        generic :: assignment(=) => copyFrom
        
    end type cScalar_t
    
    abstract interface
        
        !
        !************************************************
        ! Input/Output
        !************************************************
        !
        subroutine interface_read_c_scalar(self, fid, ftype)
            import :: cScalar_t
            class(cScalar_t)           , intent(inout) :: self
            integer                          , intent(in)        :: fid
            character(*), optional, intent(in)        :: ftype        
        end subroutine interface_read_c_scalar
        
        subroutine interface_write_c_scalar(self, fid, ftype)
            import :: cScalar_t
            class(cScalar_t), intent(in) :: self
            integer                , intent(in) :: fid
            character(*), optional, intent(in) :: ftype
        end subroutine interface_write_c_scalar
        
        !
        !************************************************
        ! Boundary operations
        !************************************************
        !

        !**
        ! setAllBoundary Set all boundary edges to specific constant.
        !
        ! All interior edges are set to zero (Overwrites previous values).
        !
        ! INPUTS:
        !            self    Reference to the instance of scalar_real_t class.
        !            c_in    (defaults to 1.0) Constant value to set boundary edges.
        !*
        subroutine interface_setAllBoundary_c_scalar(self, c_in)
            use Constants
            import :: cScalar_t
            class(cScalar_t)        , intent(inout) :: self
            complex(kind = prec), intent(in)        :: c_in
        end subroutine interface_setAllBoundary_c_scalar

        !**
        ! setOneBoundary Set all edges (faces) on specified boundary to a constant.
        !
        ! Possible values for bdry : 'x1','x2','y1',y2','z1','z2'
        !
        ! Note: If optional argument intOnly = true, only set edges on
        ! interior of face (this is meaningless for type = 'face')
        ! for edges allow only one vector component to be set (now
        !    only for z1, z2 case.
        !
        ! INPUTS:
        !            self    Reference to the instance of scalar_real_t class.
        !            bdry    Boundary to be set.
        !                  c    Constant value.
        !     int_only    See documentation note above.
        !*
        subroutine interface_setOneBoundary_c_scalar(self, bdry, c, int_only)
            use Constants
            import :: cScalar_t
            class(cScalar_t)        , intent(inout) :: self
            character(*)               , intent(in) :: bdry
            complex(kind = prec), intent(in) :: c
            logical            , intent(in), optional :: int_only
        end subroutine interface_setOneBoundary_c_scalar

        !**
        ! SetallInterior
        !*
        subroutine interface_setAllInterior_c_scalar(self, c_in)
            use Constants
            import :: cScalar_t
            class(cScalar_t)        , intent(inout) :: self
            complex(kind = prec), intent(in)        :: c_in
        end subroutine interface_setAllInterior_c_scalar
        
        !**
        ! intBdryIndices
        !*
        subroutine interface_intBdryIndices_c_scalar(self, ind_i, ind_b)
            use Constants
            import :: cScalar_t
            class(cScalar_t)        , intent(in)    :: self
            integer, allocatable, intent(out) :: ind_i(:), ind_b(:)
        end subroutine interface_intBdryIndices_c_scalar

        !
        !****************************************
        ! Data access
        !****************************************
        !

        !**
        ! length
        ! Returns the total number of elements in this vector.
        !
        ! INPUTS:
        !            self    Reference to the instance of scalar_real_t class.
        !
        ! OUTPUTS:
        !                  n    The total number of elements (edges/faces) in this vector.
        !*
        function interface_length_c_scalar(self) result(n)
            import :: cScalar_t
            class(cScalar_t), intent(in) :: self
            integer :: n
        end function interface_length_c_scalar
        
        !**
        ! getArray
        ! Return 1D array with all elements in this vector.
        !
        ! Follows Matlab/Fortran ordering for multidimensional arrays.
        !
        ! INPUTS:
        !            self    Reference to the instance of scalar_real_t class.
        !
        ! OUTPUTS:
        !                  v    1D array with vector elements.
        !*
        subroutine interface_get_array_c_scalar(self, v)
            use Constants
            import :: cScalar_t
            class(cScalar_t)                              , intent(in)    :: self
            complex(kind = prec), allocatable, intent(out) :: v(:)
        end subroutine interface_get_array_c_scalar
        
        !**
        ! setArray
        ! Set vector components from 1D array.
        !
        ! Follows Matlab/Fortran ordering for multidimensional arrays.
        !
        ! INPUTS:
        !            self    Reference to the instance of scalar_real_t class.
        !                  v    1D array with vector elements.
        !
        !*
        subroutine interface_set_array_c_scalar(self, v)
            use Constants
            import :: cScalar_t
            class(cScalar_t)        , intent(inout) :: self
            complex(kind = prec), intent(in)        :: v(:)
        end subroutine interface_set_array_c_scalar
        
        !**
        ! Arithmetic operations
        !*

        !**
        ! zeros
        !*
        subroutine interface_zeros_c_scalar(self)
            import :: cScalar_t
            class(cScalar_t), intent(inout) :: self
        end subroutine interface_zeros_c_scalar

        function interface_add1_c_scalar(lhs, rhs) result(Eout)
            import :: cScalar_t
            class(cScalar_t), intent(in)    :: lhs, rhs
            class(cScalar_t), allocatable :: Eout
        end function interface_add1_c_scalar

        function interface_sub1_c_scalar(lhs, rhs) result(Eout)
            import :: cScalar_t
            class(cScalar_t), intent(in) :: lhs, rhs
            class(cScalar_t), allocatable :: Eout
        end function interface_sub1_c_scalar

        function interface_mult1_c_scalar(lhs, rhs) result(Eout)
            import :: cScalar_t
            class(cScalar_t), intent(in) :: lhs, rhs
            class(cScalar_t), allocatable :: Eout
        end function interface_mult1_c_scalar

        function interface_mult2_c_scalar(lhs, rhs) result(Eout)
            import :: cScalar_t, prec
            class(cScalar_t), intent(in)        :: lhs
            complex( kind=prec ), intent(in) :: rhs
            class(cScalar_t), allocatable        :: Eout
        end function interface_mult2_c_scalar
        
        function interface_mult3_c_scalar(lhs, rhs) result(Eout)
            import :: cScalar_t, rScalar_t
            class(cScalar_t), intent(in)    :: lhs
            class(rScalar_t), intent(in)    :: rhs
            class(cScalar_t), allocatable :: Eout
        end function interface_mult3_c_scalar

        subroutine interface_mults1_c_scalar(lhs, rhs)
            !    subroutine version
            import :: cScalar_t
            class(cScalar_t), intent(inout)  :: lhs
            class(cScalar_t), intent(in)  :: rhs
        end subroutine interface_mults1_c_scalar

        subroutine interface_mults2_c_scalar(lhs, c)
            !    subroutine version
            import :: cScalar_t, prec
            class(cScalar_t), intent(inout)  :: lhs
            complex( kind=prec ), intent(in) :: c
        end subroutine interface_mults2_c_scalar

        subroutine interface_mults3_c_scalar(lhs, rhs)
            !    subroutine version
            import :: cScalar_t, rScalar_t
            class(cScalar_t), intent(inout)  :: lhs
            class(rScalar_t), intent(in)  :: rhs
        end subroutine interface_mults3_c_scalar
  
        function interface_div1_c_scalar(lhs, rhs) result(Eout)
            import :: cScalar_t
            class(cScalar_t), intent(in) :: lhs, rhs
            class(cScalar_t), allocatable :: Eout
        end function interface_div1_c_scalar

        subroutine interface_div3S_c_scalar(lhs, rhs)
            ! subroutine (overwrite) version)
            import :: cScalar_t, rscalar_t
            class(cScalar_t), intent(inout) :: lhs
            class(rScalar_t), intent(in) :: rhs
        end subroutine interface_div3S_c_scalar

        subroutine interface_linCombS_c_scalar(lhs,rhs,c1,c2)
            import :: cScalar_t, prec
            class(cScalar_t), intent(inout) :: lhs
            class(cScalar_t), intent(in) :: rhs
            complex(kind=prec), intent(in)  :: c1, c2
        end subroutine interface_linCombS_c_scalar

        subroutine interface_scMultAddS_c_scalar(lhs,rhs,c)
            import :: cScalar_t, prec
            class(cScalar_t), intent(in) :: lhs
            class(cScalar_t), intent(inout)  :: rhs
            complex(kind=prec), intent(in)    :: c
        end subroutine interface_scMultAddS_c_scalar

        function interface_dot_product_c_scalar(lhs, rhs) result(r)
            import :: cScalar_t, prec
            class(cScalar_t), intent(in) :: lhs, rhs
            real(kind = prec) :: r
        end function interface_dot_product_c_scalar

        function interface_is_compatible1_c_scalar(self, rhs) result( is_compatible )
            import :: cScalar_t
            class(cScalar_t), intent(in) :: self
            class(cScalar_t), intent(in) :: rhs
            logical :: is_compatible
        end function interface_is_compatible1_c_scalar
        
        function interface_is_compatible2_c_scalar(self, rhs) result( is_compatible )
            import :: cScalar_t, rScalar_t
            class(cScalar_t), intent(in) :: self
            class(rScalar_t), intent(in) :: rhs
            logical :: is_compatible
        end function interface_is_compatible2_c_scalar

        subroutine interface_copy_from_c_scalar(self, rhs)
            import :: cScalar_t
            class(cScalar_t), intent(inout) :: self
            class(cScalar_t), intent(in)    :: rhs
        end subroutine interface_copy_from_c_scalar

    end interface
    !
contains
    !
    subroutine initializeCScalar( self )
        implicit none
        !
        class( cScalar_t ), intent( inout ) :: self
        !
        self%is_allocated = .FALSE.
        !
    end subroutine initializeCScalar
    !
end module cScalar
