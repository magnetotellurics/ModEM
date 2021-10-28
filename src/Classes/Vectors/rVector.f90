module rVector
  use Constants
  use rScalar
  
  implicit none

  private
  
  public :: rVector_t
  
  type, abstract :: rVector_t
     
     logical :: isAllocated = .false.
          
   contains
     
     !**
     ! Input/Output
     !*
     procedure(iface_Read) , deferred, public :: read
     procedure(iface_Write), deferred, public :: write
     
     !**
     ! Boundary operations
     !*
     procedure(iface_setAllBoundary), deferred :: setAllBoundary
     procedure(iface_setOneBoundary), deferred :: setOneBoundary
     procedure(iface_setAllInterior), deferred :: setAllInterior
     procedure(iface_intBdryIndices), deferred :: intBdryIndices
          
     !**
     ! Data access
     !*
     procedure(iface_length)  , deferred :: length
     procedure(iface_getArray), deferred :: getArray
     procedure(iface_setArray), deferred :: setArray
     
     !**
     ! Arithmetic/algebraic operations
     procedure(iface_zeros), deferred, public :: zeros

     procedure(iface_add1) , deferred, public :: add1
     generic :: add => add1
     generic :: operator(+) => add1
     
     procedure(iface_sub1) , deferred, public :: sub1
     generic :: sub => sub1
     generic :: operator(-) => sub1
     
     procedure(iface_mult1), deferred, public :: mult1
     procedure(iface_mult2), deferred, public, pass(self) :: mult2  
     generic :: mult => mult1, mult2
     generic :: operator(*) => mult1, mult2
     
     procedure(iface_div1) , deferred, public :: div1
     generic :: div => div1
     generic :: operator(/) => div1
     
     procedure(iface_dotProd), deferred, public :: dotProd
     generic :: operator(.dot.) => dotProd
     
     procedure(iface_diagMult), deferred, public :: diagMult

     procedure(iface_copyFrom), deferred, public :: copyFrom
     generic :: assignment(=) => copyFrom

     !**
     ! Miscellaneous
     !*
     procedure(iface_isCompatible), deferred, public :: isCompatible
     procedure(iface_interpFunc), deferred, public :: InterpFunc

     procedure(iface_SumEdges), deferred, public :: SumEdges     
     procedure(iface_SumCells), deferred, public :: SumCells
  end type rVector_t
  
  abstract interface
     !**
     ! Input/Output
     !*

     !**
     ! Read
     !*
     subroutine iface_Read(self, fid)
       import :: rVector_t
       class(rVector_t), intent(inout) :: self
       integer         , intent(in)    :: fid
     end subroutine iface_Read

     !**
     ! write
     !*
     subroutine iface_Write(self, fid)
       import :: rVector_t
       class(rVector_t), intent(in) :: self
       integer         , intent(in) :: fid
     end subroutine iface_Write
     
     !**
     ! Boundary operations
     !*

     !**
     ! setAllBoundary
     !*
     subroutine iface_setAllBoundary(self, c_in)
       import :: rVector_t, prec
       class(rVector_t) , intent(inout) :: self
       real(kind = prec), intent(in)    :: c_in
     end subroutine iface_setAllBoundary

     !**
     ! setOneBoundary
     !*
     subroutine iface_setOneBoundary(self, bdry, c, int_only)
       import :: rVector_t, prec
       class(rVector_t) , intent(inout) :: self
       character(*)     , intent(in)    :: bdry
       real(kind = prec), intent(in)    :: c
       logical, optional, intent(in)    :: int_only
     end subroutine iface_setOneBoundary

     !**
     ! setAllInterior
     !*
     subroutine iface_setAllInterior(self, c_in)
       import :: rVector_t, prec
       class(rVector_t) , intent(inout) :: self
       real(kind = prec), intent(in)    :: c_in
     end subroutine iface_setAllInterior

     !**
     ! intBdryIndices
     !*
     subroutine iface_intBdryIndices(self, ind_i, ind_b)
       import :: rVector_t
       class(rVector_t)    , intent(in)  :: self
       integer, allocatable, intent(out) :: ind_i(:), ind_b(:)
     end subroutine Iface_intBdryIndices

     !**
     ! Data access
     !*

     !**
     ! length
     !*
     function iface_length(self) result(n)
       import :: rVector_t
       class(rVector_t), intent(in) :: self
       integer :: n
     end function iface_length

     !**
     ! getArray
     !*
     subroutine iface_getArray(self, v)
       import :: rVector_t, prec
       class(rVector_t), intent(in)  :: self
       real(kind = prec), allocatable, intent(out) :: v(:)
     end subroutine Iface_getArray

     !**
     ! setArray
     !*
     subroutine iface_setArray(self, v)
       import :: rVector_t, prec
       class(rVector_t) , intent(inout) :: self
       real(kind = prec), intent(in)    :: v(:)
     end subroutine iface_setArray
     
     !**
     ! Arithmetic/algebraic operations
     !*
     
     !**
     ! zeros
     !*
     subroutine iface_zeros(self)
       import :: rVector_t
       class(rVector_t), intent(inout) :: self
     end subroutine iface_zeros

     function iface_add1(lhs, rhs) result(Eout)
       import :: rVector_t
       class(rVector_t), intent(in)  :: lhs, rhs
       class(rVector_t), allocatable :: Eout
     end function iface_add1

     function iface_sub1(lhs, rhs) result(Eout)
       import :: rVector_t
       class(rVector_t), intent(in) :: lhs, rhs
       class(rVector_t), allocatable :: Eout
     end function iface_sub1

     function iface_mult1(lhs, rhs) result(Eout)
       import :: rVector_t
       class(rVector_t), intent(in) :: lhs, rhs
       class(rVector_t), allocatable :: Eout
     end function iface_mult1
     
     function iface_mult2(c, self) result(Eout)
       import :: rVector_t, prec
       real(kind = prec), intent(in) :: c
       class(rVector_t) , intent(in) :: self
       class(rVector_t), allocatable :: Eout
     end function iface_mult2
     
     function iface_div1(lhs, rhs) result(Eout)
       import :: rVector_t
       class(rVector_t), intent(in) :: lhs, rhs
       class(rVector_t), allocatable :: Eout
     end function iface_div1

     function iface_dotProd(lhs, rhs) result(r)
       import :: rVector_t, prec
       class(rVector_t), intent(in) :: lhs, rhs
       real(kind = prec) :: r
     end function iface_dotProd
     !
    function iface_diagMult( self, rhs ) result ( diag_mult )
        import :: rVector_t
        class(rVector_t), intent(in)  :: self
        class(rVector_t), intent(in)  :: rhs
        class(rVector_t), allocatable :: diag_mult
    end function iface_diagMult

     function iface_isCompatible(self, rhs) result(status)
       import :: rVector_t
       class(rVector_t), intent(in) :: self
       class(rVector_t), intent(in) :: rhs
       logical :: status
     end function iface_isCompatible

     subroutine iface_copyFrom(self, rhs)
       import :: rVector_t
       class(rVector_t), intent(inout) :: self
       class(rVector_t), intent(in)    :: rhs
     end subroutine iface_copyFrom

     !**
     ! Create a Vector object containing weights needed for
     ! interpolation of xyz component of obj1 to location.
     !*
     function iface_InterpFunc(self, location, xyz) result(E)
       import :: rVector_t, prec
       class(rVector_t)  , intent(in) :: self
       real(kind = prec), intent(in) :: location(3)
       character        , intent(in) :: xyz
       class(rVector_t), allocatable :: E
     end function iface_InterpFunc

     !**
     ! SumEdges
     ! Sum all edges (or faces) that bound a cell, return as
     ! rScalar object.
     ! If interior only is true, only sum over interior edges (faces).
     function iface_SumEdges(self, InteriorOnly) result(cellObj)
       import :: rVector_t, rScalar_t
       class(rVector_t) , intent(in) :: self
       logical, optional, intent(in) :: InteriorOnly
       class(rScalar_t), allocatable :: cellObj
     end function iface_SumEdges

     !**
     ! SumCells
     ! Sum scalar object (cell type only) onto rVector object;
     ! default is to sum onto all bounding interior edges;
     ! boundary edges are presently zero;
     ! faces case is coded now.
     !*
     subroutine iface_SumCells(self, E_in, ptype)
       import :: rVector_t, rScalar_t
       class(rVector_t), intent(inout) :: self
       class(rScalar_t), intent(in)    :: E_in
       character(*)    , intent(in), optional :: ptype
     end subroutine iface_SumCells
  end interface
  
end module rVector
