module rScalar
  use Constants
  use Grid
  use cScalar
  
  implicit none
  
  private
  
  public :: rScalar_t
  
  type, abstract :: rScalar_t
     
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
     procedure(iface_setAllBoundary), deferred, public :: setAllBoundary
     procedure(iface_setOneBoundary), deferred, public :: setOneBoundary
     procedure(iface_setAllInterior), deferred, public :: setAllInterior
     procedure(iface_intBdryIndices), deferred, public :: intBdryIndices
          
     !**
     ! Data access
     !*
     procedure(iface_length)  , deferred, public :: length
     procedure(iface_getArray), deferred, public :: getArray
     procedure(iface_setArray), deferred, public :: setArray
     
     !**
     ! Arithmetic/algebraic operations
     !*
     procedure(iface_zeros) , deferred, public :: zeros
     procedure(iface_add1) , deferred, public :: add1
     generic :: add => add1
     generic :: operator(+) => add1
     
     procedure(iface_sub1) , deferred, public :: sub1
     generic :: sub => sub1
     generic :: operator(-) => sub1
     
     procedure(iface_mult1), deferred, public :: mult1
	 procedure(iface_mult2), deferred, public :: mult2
     generic :: mult => mult1, mult2
     generic :: operator(*) => mult1, mult2
     
     procedure(iface_div1) , deferred, public :: div1
     generic :: div => div1
     generic :: operator(/) => div1
     
     procedure(iface_dotProd), deferred, public :: dotProd
     generic :: operator(.dot.) => dotProd

     !**
     ! Miscellaneous
     !*
     procedure(iface_isCompatible), deferred, public :: isCompatible
     procedure(iface_copyFrom), deferred, public :: copyFrom
     generic :: assignment(=) => copyFrom
     
  end type rScalar_t
  
  abstract interface
     
     !
     !************************************************
     ! Input/Output
     !************************************************
     !
     subroutine iface_Read(self, fid, ftype)
       import :: rScalar_t
       class(rScalar_t)      , intent(inout) :: self
       integer               , intent(in)    :: fid
       character(*), optional, intent(in)    :: ftype    
     end subroutine iface_Read
     
     subroutine iface_Write(self, fid, ftype)
       import :: rScalar_t
       class(rScalar_t), intent(in) :: self
       integer         , intent(in) :: fid
       character(*), optional, intent(in) :: ftype
     end subroutine iface_Write
     
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
     !       self  Reference to the instance of scalar_real_t class.
     !       c_in  (defaults to 1.0) Constant value to set boundary edges.
     !*
     subroutine iface_setAllBoundary(self, c_in)
       use Constants
       import :: rScalar_t
       class(rScalar_t) , intent(inout) :: self
       real(kind = prec), intent(in)    :: c_in
     end subroutine iface_setAllBoundary

     !**
     ! setOneBoundary Set all edges (faces) on specified boundary to a constant.
     !
     ! Possible values for bdry : 'x1','x2','y1',y2','z1','z2'
     !
     ! Note: If optional argument intOnly = true, only set edges on
     ! interior of face (this is meaningless for type = 'face')
     ! for edges allow only one vector component to be set (now
     !  only for z1, z2 case.
     !
     ! INPUTS:
     !       self  Reference to the instance of scalar_real_t class.
     !       bdry  Boundary to be set.
     !          c  Constant value.
     !   int_only  See documentation note above.
     !*
     subroutine iface_setOneBoundary(self, bdry, c, int_only)
       use Constants
       import :: rScalar_t
       class(rScalar_t) , intent(inout) :: self
       character(*)     , intent(in) :: bdry
       real(kind = prec), intent(in) :: c
       logical       , intent(in), optional :: int_only
     end subroutine iface_setOneBoundary

     !**
     ! SetallInterior
     !*
     subroutine iface_setAllInterior(self, c_in)
       use Constants
       import :: rScalar_t
       class(rScalar_t) , intent(inout) :: self
       real(kind = prec), intent(in)    :: c_in
     end subroutine iface_setAllInterior
     
     !**
     ! intBdryIndices
     !*
     subroutine iface_intBdryIndices(self, ind_i, ind_b)
       use Constants
       import :: rScalar_t
       class(rScalar_t)    , intent(in)  :: self
       integer, allocatable, intent(out) :: ind_i(:), ind_b(:)
     end subroutine iface_intBdryIndices

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
     !       self  Reference to the instance of scalar_real_t class.
     !
     ! OUTPUTS:
     !          n  The total number of elements (edges/faces) in this vector.
     !*
     function iface_length(self) result(n)
       import :: rScalar_t
       class(rScalar_t), intent(in) :: self
       integer :: n
     end function iface_length
     
     !**
     ! getArray
     ! Return 1D array with all elements in this vector.
     !
     ! Follows Matlab/Fortran ordering for multidimensional arrays.
     !
     ! INPUTS:
     !       self  Reference to the instance of scalar_real_t class.
     !
     ! OUTPUTS:
     !          v  1D array with vector elements.
     !*
     subroutine iface_getArray(self, v)
       use Constants
       import :: rScalar_t
       class(rScalar_t)              , intent(in)  :: self
       real(kind = prec), allocatable, intent(out) :: v(:)
     end subroutine iface_getArray
     
     !**
     ! setArray
     ! Set vector components from 1D array.
     !
     ! Follows Matlab/Fortran ordering for multidimensional arrays.
     !
     ! INPUTS:
     !       self  Reference to the instance of scalar_real_t class.
     !          v  1D array with vector elements.
     !
     !*
     subroutine iface_setArray(self, v)
       use Constants
       import :: rScalar_t
       class(rScalar_t)  , intent(inout) :: self
       real(kind = prec), intent(in)    :: v(:)
     end subroutine iface_setArray
     
     !**
     ! Arithmetic operations
     !*

     !**
     ! zeros
     !*
     subroutine iface_zeros(self)
       import :: rScalar_t
       class(rScalar_t), intent(inout) :: self
     end subroutine iface_zeros

     function iface_add1(lhs, rhs) result(Eout)
       import :: rScalar_t
       class(rScalar_t), intent(in)  :: lhs, rhs
       class(rScalar_t), allocatable :: Eout
     end function iface_add1

     function iface_sub1(lhs, rhs) result(Eout)
       import :: rScalar_t
       class(rScalar_t), intent(in) :: lhs, rhs
       class(rScalar_t), allocatable :: Eout
     end function iface_sub1

     function iface_mult1(lhs, rhs) result(Eout)
       import :: rScalar_t
       class(rScalar_t), intent(in) :: lhs, rhs
       class(rScalar_t), allocatable :: Eout
     end function iface_mult1
	 
	 function iface_mult2(lhs, rhs) result(Eout)
       import :: rScalar_t, cScalar_t
       class(rScalar_t), intent(in)  :: lhs
	   class(cScalar_t), intent(in)  :: rhs
       class(cScalar_t), allocatable :: Eout
     end function iface_mult2
     
     function iface_div1(lhs, rhs) result(Eout)
       import :: rScalar_t
       class(rScalar_t), intent(in) :: lhs, rhs
       class(rScalar_t), allocatable :: Eout
     end function iface_div1

     function iface_dotProd(lhs, rhs) result(r)
       import :: rScalar_t, prec
       class(rScalar_t), intent(in) :: lhs, rhs
       real(kind = prec) :: r
     end function iface_dotProd

     function iface_isCompatible(self, rhs) result(status)
       import :: rScalar_t
       class(rScalar_t), intent(in) :: self
       class(rScalar_t), intent(in) :: rhs
       logical :: status
     end function iface_isCompatible

     subroutine iface_copyFrom(self, rhs)
       import :: rScalar_t
       class(rScalar_t), intent(inout) :: self
       class(rScalar_t), intent(in)    :: rhs
     end subroutine iface_copyFrom


  end interface
  
end module rScalar
