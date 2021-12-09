module rVector
   !
   use Constants
   use rScalar
   !
   type, abstract :: rVector_t
       
       logical :: isAllocated = .false.
       !
    contains
       
        !**
        ! Input/Output
        !*
        procedure( interface_read_r_vector ) , deferred, public :: read
        procedure( interface_write_r_vector ), deferred, public :: write

        !**
        ! Boundary operations
        !*
        procedure( interface_set_all_boundary_r_vector ), deferred :: setAllBoundary
        procedure( interface_set_one_boundary_r_vector ), deferred :: setOneBoundary
        procedure( interface_set_all_interior_r_vector ), deferred :: setAllInterior
        procedure( interface_int_bdry_indices_r_vector ), deferred :: intBdryIndices
		procedure( interface_boundary_r_vector )         , deferred :: boundary
        procedure( interface_interior_r_vector )         , deferred :: interior

        !**
        ! Data access
        !*
        procedure( interface_length_r_vector )   , deferred :: length
        procedure( interface_get_array_r_vector ), deferred :: getArray
        procedure( interface_set_array_r_vector ), deferred :: setArray

        !**
        ! Arithmetic/algebraic operations
        procedure( interface_zeros_r_vector ), deferred, public :: zeros

        procedure( interface_add1_r_vector ) , deferred, public :: add1
        generic :: add => add1
        generic :: operator(+) => add1

        procedure( interface_sub1_r_vector ) , deferred, public :: sub1
        generic :: sub => sub1
        generic :: operator(-) => sub1

        procedure( interface_mult1_r_vector ), deferred, public :: mult1
        procedure( interface_mult2_r_vector ), deferred, public, pass(self) :: mult2   
        generic :: mult => mult1, mult2
        generic :: operator(*) => mult1, mult2

        procedure( interface_div1_r_vector ) , deferred, public :: div1
        generic :: div => div1
        generic :: operator(/) => div1

        !   subroutine versions -- all overwrite first argument
        !     divide rVector by rVector
        procedure( interface_divs1_r_vector ) , deferred, public :: divS1
        generic :: divS => divS1
        !     multiply rVector by rVector
        procedure( interface_mults1_r_vector ) , deferred, public :: multS1
        !     multiply rVector by real
        procedure( interface_mults2_r_vector ) , deferred, public :: multS2
        generic :: multS => multS1, multS2


        procedure( interface_dot_product_r_vector ), deferred, public :: dotProd
        generic :: operator(.dot.) => dotProd

        procedure( interface_diag_mult_r_vector ), deferred, public :: diagMult

        procedure( interface_copy_from_r_vector ), deferred, public :: copyFrom
        generic :: assignment(=) => copyFrom

        !**
        ! Miscellaneous
        !*
        procedure( interface_is_compatible_r_vector ), deferred, public :: isCompatible
        procedure( interface_interp_func_r_vector ), deferred, public :: InterpFunc

        procedure( interface_sum_edges_r_vector ), deferred, public :: SumEdges
        procedure( interface_sum_cells_r_vector ), deferred, public :: SumCells
		
   end type rVector_t
   
   abstract interface
       !**
       ! Input/Output
       !*

       !**
       ! Read
       !*
       subroutine interface_read_r_vector(self, fid)
          import :: rVector_t
          class(rVector_t), intent(inout) :: self
          integer             , intent(in)      :: fid
       end subroutine interface_read_r_vector

       !**
       ! write
       !*
       subroutine interface_write_r_vector(self, fid)
          import :: rVector_t
          class(rVector_t), intent(in) :: self
          integer             , intent(in) :: fid
       end subroutine interface_write_r_vector
       
       !**
       ! Boundary operations
       !*

       !**
       ! setAllBoundary
       !*
       subroutine interface_set_all_boundary_r_vector(self, c_in)
          import :: rVector_t, prec
          class(rVector_t) , intent(inout) :: self
          real(kind = prec), intent(in)      :: c_in
       end subroutine interface_set_all_boundary_r_vector

       !**
       ! setOneBoundary
       !*
       subroutine interface_set_one_boundary_r_vector(self, bdry, c, int_only)
          import :: rVector_t, prec
          class(rVector_t) , intent(inout) :: self
          character(*)       , intent(in)      :: bdry
          real(kind = prec), intent(in)      :: c
          logical, optional, intent(in)      :: int_only
       end subroutine interface_set_one_boundary_r_vector

       !**
       ! setAllInterior
       !*
       subroutine interface_set_all_interior_r_vector(self, c_in)
          import :: rVector_t, prec
          class(rVector_t) , intent(inout) :: self
          real(kind = prec), intent(in)      :: c_in
       end subroutine interface_set_all_interior_r_vector

       !**
       ! intBdryIndices
       !*
       subroutine interface_int_bdry_indices_r_vector(self, ind_i, ind_b)
          import :: rVector_t
          class(rVector_t)      , intent(in)   :: self
          integer, allocatable, intent(out) :: ind_i(:), ind_b(:)
       end subroutine interface_int_bdry_indices_r_vector

       !**
       ! boundary
       !*
       function interface_boundary_r_vector(self) result(E)
          import :: rVector_t
          class(rVector_t), intent(in) :: self
          class(rVector_t), allocatable :: E
       end function interface_boundary_r_vector

       !**
       ! interior
       !*
       function interface_interior_r_vector(self) result(E)
          import :: rVector_t
          class(rVector_t), intent(in) :: self
          class(rVector_t), allocatable :: E
       end function interface_interior_r_vector
       !**
       ! Data access
       !*

       !**
       ! length
       !*
       function interface_length_r_vector(self) result(n)
          import :: rVector_t
          class(rVector_t), intent(in) :: self
          integer :: n
       end function interface_length_r_vector

       !**
       ! getArray
       !*
       subroutine interface_get_array_r_vector(self, v)
          import :: rVector_t, prec
          class(rVector_t), intent(in)   :: self
          real(kind = prec), allocatable, intent(out) :: v(:)
       end subroutine interface_get_array_r_vector

       !**
       ! setArray
       !*
       subroutine interface_set_array_r_vector(self, v)
          import :: rVector_t, prec
          class(rVector_t) , intent(inout) :: self
          real(kind = prec), intent(in)      :: v(:)
       end subroutine interface_set_array_r_vector
       
       !**
       ! Arithmetic/algebraic operations
       !*
       
       !**
       ! zeros
       !*
       subroutine interface_zeros_r_vector(self)
          import :: rVector_t
          class(rVector_t), intent(inout) :: self
       end subroutine interface_zeros_r_vector

       function interface_add1_r_vector(lhs, rhs) result(Eout)
          import :: rVector_t
          class(rVector_t), intent(in)   :: lhs, rhs
          class(rVector_t), allocatable :: Eout
       end function interface_add1_r_vector

       function interface_sub1_r_vector(lhs, rhs) result(Eout)
          import :: rVector_t
          class(rVector_t), intent(in) :: lhs, rhs
          class(rVector_t), allocatable :: Eout
       end function interface_sub1_r_vector

       function interface_mult1_r_vector(lhs, rhs) result(Eout)
          import :: rVector_t
          class(rVector_t), intent(in) :: lhs, rhs
          class(rVector_t), allocatable :: Eout
       end function interface_mult1_r_vector
       
       function interface_mult2_r_vector(c, self) result(Eout)
          import :: rVector_t, prec
          real(kind = prec), intent(in) :: c
          class(rVector_t) , intent(in) :: self
          class(rVector_t), allocatable :: Eout
       end function interface_mult2_r_vector
       
       function interface_div1_r_vector(lhs, rhs) result(Eout)
          import :: rVector_t
          class(rVector_t), intent(in) :: lhs, rhs
          class(rVector_t), allocatable :: Eout
       end function interface_div1_r_vector

       subroutine interface_divs1_r_vector( lhs, rhs )
          import :: rVector_t
          class(rVector_t), intent(inout) :: lhs
          class(rVector_t), intent(in)    :: rhs
       end subroutine interface_divs1_r_vector

       subroutine interface_mults1_r_vector(lhs,rhs)
          import :: rVector_t
          class(rVector_t), intent(inout) :: lhs
          class(rVector_t), intent(in)    :: rhs
       end subroutine interface_mults1_r_vector

       subroutine interface_mults2_r_vector(lhs,r)
          import :: rVector_t, prec
          class(rVector_t), intent(inout) :: lhs
          real(kind=prec), intent(in)    :: r
       end subroutine interface_mults2_r_vector
       !   end subroutine versions

       function interface_dot_product_r_vector(lhs, rhs) result(r)
          import :: rVector_t, prec
          class(rVector_t), intent(in) :: lhs, rhs
          real(kind = prec) :: r
       end function interface_dot_product_r_vector
       !
      function interface_diag_mult_r_vector( self, rhs ) result ( diag_mult )
            import :: rVector_t
            class(rVector_t), intent(in)   :: self
            class(rVector_t), intent(in)   :: rhs
            class(rVector_t), allocatable :: diag_mult
      end function interface_diag_mult_r_vector

       function interface_is_compatible_r_vector(self, rhs) result(status)
          import :: rVector_t
          class(rVector_t), intent(in) :: self
          class(rVector_t), intent(in) :: rhs
          logical :: status
       end function interface_is_compatible_r_vector

       subroutine interface_copy_from_r_vector(self, rhs)
          import :: rVector_t
          class(rVector_t), intent(inout) :: self
          class(rVector_t), intent(in)      :: rhs
       end subroutine interface_copy_from_r_vector

       !**
       ! Create a Vector object containing weights needed for
       ! interpolation of xyz component of obj1 to location.
       !*
       function interface_interp_func_r_vector(self, location, xyz) result(E)
          import :: rVector_t, prec
          class(rVector_t)   , intent(in) :: self
          real(kind = prec), intent(in) :: location(3)
          character            , intent(in) :: xyz
          class(rVector_t), allocatable :: E
       end function interface_interp_func_r_vector

       !**
       ! SumEdges
       ! Sum all edges (or faces) that bound a cell, return as
       ! rScalar object.
       ! If interior only is true, only sum over interior edges (faces).
       function interface_sum_edges_r_vector(self, InteriorOnly) result(cellObj)
          import :: rVector_t, rScalar_t
          class(rVector_t) , intent(in) :: self
          logical, optional, intent(in) :: InteriorOnly
          class(rScalar_t), allocatable :: cellObj
       end function interface_sum_edges_r_vector

       !**
       ! SumCells
       ! Sum scalar object (cell type only) onto rVector object;
       ! default is to sum onto all bounding interior edges;
       ! boundary edges are presently zero;
       ! faces case is coded now.
       !*
       subroutine interface_sum_cells_r_vector(self, E_in, ptype)
          import :: rVector_t, rScalar_t
          class(rVector_t), intent(inout) :: self
          class(rScalar_t), intent(in)      :: E_in
          character(*)      , intent(in), optional :: ptype
       end subroutine interface_sum_cells_r_vector
       
   end interface
   
end module rVector
