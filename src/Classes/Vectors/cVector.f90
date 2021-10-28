module cVector
  use Constants
  use rScalar
  use rVector
  
  type, abstract :: cVector_t
     
     logical :: isAllocated = .false.
          
   contains
     
     !**
     ! Input/Output
     !*
     procedure(iface_Read) , deferred, public :: read
     procedure(iface_Write), deferred, public :: write
     
     !**
     ! boundary operations
     !*
     procedure(iface_setAllboundary), deferred :: setAllboundary
     procedure(iface_setOneboundary), deferred :: setOneboundary
     procedure(iface_setAllinterior), deferred :: setAllinterior
     procedure(iface_intBdryIndices), deferred :: intBdryIndices
     procedure(iface_boundary)      , deferred :: boundary
     procedure(iface_interior)      , deferred :: interior

     !**
     ! Data access
     !*
     procedure(iface_length)  , deferred :: length
     procedure(iface_getArray), deferred :: getArray
     procedure(iface_setArray), deferred :: setArray
     
     !**
     ! Arithmetic/algebraic operations
     procedure(iface_zeros_cVector), deferred, public :: zeros

     procedure(iface_add1) , deferred, public :: add1
     generic :: add => add1
     generic :: operator(+) => add1
     
     procedure(iface_sub1) , deferred, public :: sub1
     generic :: sub => sub1
     generic :: operator(-) => sub1
     
     procedure(iface_mult1), deferred, public :: mult1
     procedure(iface_mult2), deferred, public, pass(self) :: mult2
     procedure(iface_mult3), deferred, public :: mult3     
     generic :: mult => mult1, mult2, mult3
     generic :: operator(*) => mult1, mult2, mult3

     procedure(iface_mults1), deferred, public :: mults1
     generic :: mults => mults1
     
     procedure(iface_div1) , deferred, public :: div1
     procedure(iface_div2) , deferred, public :: div2     
     generic :: div => div1, div2
     generic :: operator(/) => div1, div2
     
     procedure(iface_dotProd), deferred, public :: dotProd
     generic :: operator(.dot.) => dotProd

     procedure(iface_copyFrom_cVector), deferred, public :: copyFrom
     generic :: assignment(=) => copyFrom
     
     !**
     ! Miscellaneous
     !*     
     procedure(iface_interpFunc), deferred :: interpFunc

  end type cVector_t
  
  abstract interface
     !**
     ! Input/Output
     !*

     !**
     ! Read
     !*
     subroutine iface_Read(self, fid)
       import :: cVector_t
       class(cVector_t), intent(inout) :: self
       integer         , intent(in)    :: fid
     end subroutine iface_Read

     !**
     ! write
     !*
     subroutine iface_Write(self, fid)
       import :: cVector_t
       class(cVector_t), intent(in) :: self
       integer         , intent(in) :: fid
     end subroutine iface_Write
     
     !**
     ! boundary operations
     !*

     !**
     ! setAllboundary
     !*
     subroutine iface_setAllboundary(self, c_in)
       import :: cVector_t, prec
       class(cVector_t)    , intent(inout) :: self
       complex(kind = prec), intent(in)    :: c_in
     end subroutine iface_setAllboundary

     !**
     ! setOneboundary
     !*
     subroutine iface_setOneboundary(self, bdry, c, int_only)
       import :: cVector_t, prec
       class(cVector_t)    , intent(inout) :: self
       character(*)        , intent(in)    :: bdry
       complex(kind = prec), intent(in)    :: c
       logical, optional, intent(in)    :: int_only
     end subroutine iface_setOneboundary

     !**
     ! setAllinterior
     !*
     subroutine iface_setAllinterior(self, c_in)
       import :: cVector_t, prec
       class(cVector_t)    , intent(inout) :: self
       complex(kind = prec), intent(in)    :: c_in
     end subroutine iface_setAllinterior

     !**
     ! intBdryIndices
     !*
     subroutine iface_intBdryIndices(self, ind_i, ind_b)
       import :: cVector_t
       class(cVector_t)    , intent(in)  :: self
       integer, allocatable, intent(out) :: ind_i(:), ind_b(:)
     end subroutine Iface_intBdryIndices
     
     !**
     ! boundary
     !*
     function iface_boundary(self) result(E)
       import :: cVector_t
       class(cVector_t), intent(in) :: self
       class(cVector_t), allocatable :: E
     end function Iface_boundary

     !**
     ! interior
     !*
     function iface_interior(self) result(E)
       import :: cVector_t
       class(cVector_t), intent(in) :: self
       class(cVector_t), allocatable :: E
     end function Iface_interior


     !**
     ! Data access
     !*

     !**
     ! length
     !*
     function iface_length(self) result(n)
       import :: cVector_t
       class(cVector_t), intent(in) :: self
       integer :: n
     end function iface_length

     !**
     ! getArray
     !*
     subroutine iface_getArray(self, v)
       import :: cVector_t, prec
       class(cVector_t), intent(in)  :: self
       complex(kind = prec), allocatable, intent(out) :: v(:)
     end subroutine Iface_getArray

     !**
     ! setArray
     !*
     subroutine iface_setArray(self, v)
       import :: cVector_t, prec
       class(cVector_t)    , intent(inout) :: self
       complex(kind = prec), intent(in)    :: v(:)
     end subroutine iface_setArray
     
     !**
     ! Arithmetic/algebraic operations
     !*
     
     !**
     ! zeros
     !*
     subroutine iface_zeros_cVector(self)
       import :: cVector_t
       class(cVector_t), intent(inout) :: self
     end subroutine iface_zeros_cVector

     function iface_add1(lhs, rhs) result(Eout)
       import :: cVector_t
       class(cVector_t), intent(in)  :: lhs, rhs
       class(cVector_t), allocatable :: Eout
     end function iface_add1

     function iface_sub1(lhs, rhs) result(Eout)
       import :: cVector_t
       class(cVector_t), intent(in) :: lhs, rhs
       class(cVector_t), allocatable :: Eout
     end function iface_sub1

     function iface_mult1(lhs, rhs) result(Eout)
       import :: cVector_t
       class(cVector_t), intent(in) :: lhs, rhs
       class(cVector_t), allocatable :: Eout
     end function iface_mult1

     function iface_mult2(c, self) result(Eout)
       import :: cVector_t, prec
       complex(kind = prec), intent(in) :: c
       class(cVector_t)    , intent(in) :: self
       class(cVector_t), allocatable :: Eout
     end function iface_mult2

     function iface_mult3(lhs, rhs) result(Eout)
       import :: cVector_t, rVector_t
       class(cVector_t), intent(in) :: lhs
       class(rVector_t), intent(in) :: rhs       
       class(cVector_t), allocatable :: Eout
     end function iface_mult3
     
     subroutine iface_mults1(self, c)
       import :: cVector_t, prec
       class(cVector_t)    , intent(inout) :: self
       complex(kind = prec), intent(in) :: c       
     end subroutine iface_mults1
     
     function iface_div1(lhs, rhs) result(Eout)
       import :: cVector_t
       class(cVector_t), intent(in) :: lhs
       class(cVector_t), intent(in) :: rhs
       class(cVector_t), allocatable :: Eout
     end function iface_div1

     function iface_div2(lhs, rhs) result(Eout)
       import :: cVector_t, rVector_t
       class(cVector_t), intent(in) :: lhs
       class(rVector_t), intent(in) :: rhs
       class(cVector_t), allocatable :: Eout
     end function iface_div2
     
     function iface_dotProd(lhs, rhs) result(r)
       import :: cVector_t, prec
       class(cVector_t), intent(in) :: lhs, rhs
       complex(kind = prec) :: r
     end function iface_dotProd

     subroutine iface_copyFrom_cVector(self, rhs)
       import :: cVector_t
       class(cVector_t), intent(inout) :: self
       class(cVector_t), intent(in)    :: rhs
     end subroutine iface_copyFrom_cVector
     
     !**
     ! Create a Vector object containing weights needed for
     ! interpolation of xyz component of obj1 to location.
     !*
     subroutine iface_interpFunc(self, location, xyz, E)
       import :: cVector_t, prec
       class(cVector_t) , intent(in)  :: self
       real(kind = prec), intent(in)  :: location(3)
       character        , intent(in)  :: xyz
       class(cVector_t) , intent(out), allocatable :: E
     end subroutine iface_interpFunc
     
  end interface
  
end module cVector
