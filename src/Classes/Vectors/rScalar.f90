module rScalar
  use Constants
  use Grid
  
  implicit none
  
  private
  
  public :: rScalar_t
  
  type, abstract :: rScalar_t
     
     logical :: isAllocated = .false.
     
   contains
     
     !**
     ! Input/Output
     !*
     procedure(iface_Read) , deferred, public :: Read
     procedure(iface_Write), deferred, public :: Write
     
     !**
     ! Boundary operations
     !*
     procedure(iface_SetAllBoundary), deferred, public :: SetAllBoundary
     procedure(iface_SetOneBoundary), deferred, public :: SetOneBoundary
     procedure(iface_SetAllInterior), deferred, public :: SetAllInterior
     procedure(iface_IntBdryIndices), deferred, public :: IntBdryIndices
          
     !**
     ! Data access
     !*
     procedure(iface_Length)  , deferred, public :: Length
     procedure(iface_GetArray), deferred, public :: GetArray
     procedure(iface_SetArray), deferred, public :: SetArray
     
     !**
     ! Arithmetic/algebraic operations
     !*
     procedure(iface_Zeros) , deferred, public :: Zeros
     procedure(iface_Add_1) , deferred, public :: Add_1
     generic :: Add => Add_1
     generic :: operator(+) => Add_1
     
     procedure(iface_Sub_1) , deferred, public :: Sub_1
     generic :: Sub => Sub_1
     generic :: operator(-) => Sub_1
     
     procedure(iface_Mult_1), deferred, public :: Mult_1
     generic :: Mult => Mult_1
     generic :: operator(*) => Mult_1
     
     procedure(iface_Div_1) , deferred, public :: div_1
     generic :: Div => Div_1
     generic :: operator(/) => Div_1
     
     procedure(iface_dotProd), deferred, public :: dotProd
     generic :: operator(.dot.) => dotProd

     !**
     ! Miscellaneous
     !*
     procedure(iface_isCompatible), deferred, public :: isCompatible
     procedure(iface_CopyFrom), deferred, public :: CopyFrom
     generic :: assignment(=) => CopyFrom
     
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
     ! SetAllBoundary Set all boundary edges to specific constant.
     !
     ! All interior edges are set to zero (Overwrites previous values).
     !
     ! INPUTS:
     !       self  Reference to the instance of scalar_real_t class.
     !       c_in  (defaults to 1.0) Constant value to set boundary edges.
     !*
     subroutine iface_SetAllBoundary(self, c_in)
       use Constants
       import :: rScalar_t
       class(rScalar_t) , intent(inout) :: self
       real(kind = prec), intent(in)    :: c_in
     end subroutine iface_SetAllBoundary

     !**
     ! SetOneBoundary Set all edges (faces) on specified boundary to a constant.
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
     subroutine iface_SetOneBoundary(self, bdry, c, int_only)
       use Constants
       import :: rScalar_t
       class(rScalar_t) , intent(inout) :: self
       character(*)     , intent(in) :: bdry
       real(kind = prec), intent(in) :: c
       logical       , intent(in), optional :: int_only
     end subroutine iface_SetOneBoundary

     !**
     ! SetallInterior
     !*
     subroutine iface_SetAllInterior(self, c_in)
       use Constants
       import :: rScalar_t
       class(rScalar_t) , intent(inout) :: self
       real(kind = prec), intent(in)    :: c_in
     end subroutine iface_SetAllInterior
     
     !**
     ! IntBdryIndices
     !*
     subroutine iface_IntBdryIndices(self, ind_i, ind_b)
       use Constants
       import :: rScalar_t
       class(rScalar_t)    , intent(in)  :: self
       integer, allocatable, intent(out) :: ind_i(:), ind_b(:)
     end subroutine iface_IntBdryIndices

     !
     !****************************************
     ! Data access
     !****************************************
     !

     !**
     ! Length
     ! Returns the total number of elements in this vector.
     !
     ! INPUTS:
     !       self  Reference to the instance of scalar_real_t class.
     !
     ! OUTPUTS:
     !          n  The total number of elements (edges/faces) in this vector.
     !*
     function iface_Length(self) result(n)
       import :: rScalar_t
       class(rScalar_t), intent(in) :: self
       integer :: n
     end function iface_Length
     
     !**
     ! GetArray
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
     subroutine iface_GetArray(self, v)
       use Constants
       import :: rScalar_t
       class(rScalar_t)              , intent(in)  :: self
       real(kind = prec), allocatable, intent(out) :: v(:)
     end subroutine iface_GetArray
     
     !**
     ! SetArray
     ! Set vector components from 1D array.
     !
     ! Follows Matlab/Fortran ordering for multidimensional arrays.
     !
     ! INPUTS:
     !       self  Reference to the instance of scalar_real_t class.
     !          v  1D array with vector elements.
     !
     !*
     subroutine iface_SetArray(self, v)
       use Constants
       import :: rScalar_t
       class(rScalar_t)  , intent(inout) :: self
       real(kind = prec), intent(in)    :: v(:)
     end subroutine iface_SetArray
     
     !**
     ! Arithmetic operations
     !*

     !**
     ! Zeros
     !*
     subroutine iface_Zeros(self)
       import :: rScalar_t
       class(rScalar_t), intent(inout) :: self
     end subroutine iface_Zeros

     function iface_Add_1(lhs, rhs) result(Eout)
       import :: rScalar_t
       class(rScalar_t), intent(in)  :: lhs, rhs
       class(rScalar_t), allocatable :: Eout
     end function iface_Add_1

     function iface_Sub_1(lhs, rhs) result(Eout)
       import :: rScalar_t
       class(rScalar_t), intent(in) :: lhs, rhs
       class(rScalar_t), allocatable :: Eout
     end function iface_Sub_1

     function iface_Mult_1(lhs, rhs) result(Eout)
       import :: rScalar_t
       class(rScalar_t), intent(in) :: lhs, rhs
       class(rScalar_t), allocatable :: Eout
     end function iface_Mult_1
     
     function iface_Div_1(lhs, rhs) result(Eout)
       import :: rScalar_t
       class(rScalar_t), intent(in) :: lhs, rhs
       class(rScalar_t), allocatable :: Eout
     end function iface_Div_1

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

     subroutine iface_CopyFrom(self, rhs)
       import :: rScalar_t
       class(rScalar_t), intent(inout) :: self
       class(rScalar_t), intent(in)    :: rhs
     end subroutine iface_CopyFrom


  end interface
  
end module rScalar
