module cVector
  use Constants
  use rScalar
  use rVector
  
  implicit none

  private

  public :: cVector_t
  
  type, abstract :: cVector_t
     
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
     procedure(iface_SetAllBoundary), deferred :: SetAllBoundary
     procedure(iface_SetOneBoundary), deferred :: SetOneBoundary
     procedure(iface_SetAllInterior), deferred :: SetAllInterior
     procedure(iface_IntBdryIndices), deferred :: IntBdryIndices
     procedure(iface_Boundary)      , deferred :: Boundary
     procedure(iface_Interior)      , deferred :: Interior

     !**
     ! Data access
     !*
     procedure(iface_Length)  , deferred :: Length
     procedure(iface_GetArray), deferred :: GetArray
     procedure(iface_SetArray), deferred :: SetArray
     
     !**
     ! Arithmetic/algebraic operations
     procedure(iface_Zeros), deferred, public :: Zeros

     procedure(iface_Add_1) , deferred, public :: Add_1
     generic :: Add => Add_1
     generic :: operator(+) => Add_1
     
     procedure(iface_Sub_1) , deferred, public :: Sub_1
     generic :: Sub => Sub_1
     generic :: operator(-) => Sub_1
     
     procedure(iface_Mult_1), deferred, public :: Mult_1
     procedure(iface_Mult_2), deferred, public, pass(self) :: Mult_2
     procedure(iface_Mult_3), deferred, public :: Mult_3     
     generic :: Mult => Mult_1, Mult_2, Mult_3
     generic :: operator(*) => Mult_1, Mult_2, Mult_3

     procedure(iface_Mult_s_1), deferred, public :: Mult_s_1
     generic :: Mult_s => Mult_s_1
     
     procedure(iface_Div_1) , deferred, public :: div_1
     procedure(iface_Div_2) , deferred, public :: div_2     
     generic :: Div => Div_1, Div_2
     generic :: operator(/) => Div_1, Div_2
     
     procedure(iface_dotProd), deferred, public :: dotProd
     generic :: operator(.dot.) => dotProd

     procedure(iface_CopyFrom), deferred, public :: CopyFrom
     generic :: assignment(=) => CopyFrom
     
     !**
     ! Miscellaneous
     !*     
     procedure(iface_interpFunc), deferred :: InterpFunc

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
     ! Boundary operations
     !*

     !**
     ! SetAllBoundary
     !*
     subroutine iface_SetAllBoundary(self, c_in)
       import :: cVector_t, prec
       class(cVector_t)    , intent(inout) :: self
       complex(kind = prec), intent(in)    :: c_in
     end subroutine iface_SetAllBoundary

     !**
     ! SetOneBoundary
     !*
     subroutine iface_SetOneBoundary(self, bdry, c, int_only)
       import :: cVector_t, prec
       class(cVector_t)    , intent(inout) :: self
       character(*)        , intent(in)    :: bdry
       complex(kind = prec), intent(in)    :: c
       logical, optional, intent(in)    :: int_only
     end subroutine iface_SetOneBoundary

     !**
     ! SetAllInterior
     !*
     subroutine iface_SetAllInterior(self, c_in)
       import :: cVector_t, prec
       class(cVector_t)    , intent(inout) :: self
       complex(kind = prec), intent(in)    :: c_in
     end subroutine iface_SetAllInterior

     !**
     ! IntBdryIndices
     !*
     subroutine iface_IntBdryIndices(self, ind_i, ind_b)
       import :: cVector_t
       class(cVector_t)    , intent(in)  :: self
       integer, allocatable, intent(out) :: ind_i(:), ind_b(:)
     end subroutine Iface_IntBdryIndices
     
     !**
     ! Boundary
     !*
     function iface_Boundary(self) result(E)
       import :: cVector_t
       class(cVector_t), intent(in) :: self
       class(cVector_t), allocatable :: E
     end function Iface_Boundary

     !**
     ! Interior
     !*
     function iface_Interior(self) result(E)
       import :: cVector_t
       class(cVector_t), intent(in) :: self
       class(cVector_t), allocatable :: E
     end function Iface_Interior


     !**
     ! Data access
     !*

     !**
     ! Length
     !*
     function iface_Length(self) result(n)
       import :: cVector_t
       class(cVector_t), intent(in) :: self
       integer :: n
     end function iface_Length

     !**
     ! GetArray
     !*
     subroutine iface_GetArray(self, v)
       import :: cVector_t, prec
       class(cVector_t), intent(in)  :: self
       complex(kind = prec), allocatable, intent(out) :: v(:)
     end subroutine Iface_GetArray

     !**
     ! SetArray
     !*
     subroutine iface_SetArray(self, v)
       import :: cVector_t, prec
       class(cVector_t)    , intent(inout) :: self
       complex(kind = prec), intent(in)    :: v(:)
     end subroutine iface_SetArray
     
     !**
     ! Arithmetic/algebraic operations
     !*
     
     !**
     ! Zeros
     !*
     subroutine iface_Zeros(self)
       import :: cVector_t
       class(cVector_t), intent(inout) :: self
     end subroutine iface_Zeros

     function iface_Add_1(lhs, rhs) result(Eout)
       import :: cVector_t
       class(cVector_t), intent(in)  :: lhs, rhs
       class(cVector_t), allocatable :: Eout
     end function iface_Add_1

     function iface_Sub_1(lhs, rhs) result(Eout)
       import :: cVector_t
       class(cVector_t), intent(in) :: lhs, rhs
       class(cVector_t), allocatable :: Eout
     end function iface_Sub_1

     function iface_Mult_1(lhs, rhs) result(Eout)
       import :: cVector_t
       class(cVector_t), intent(in) :: lhs, rhs
       class(cVector_t), allocatable :: Eout
     end function iface_Mult_1

     function iface_Mult_2(c, self) result(Eout)
       import :: cVector_t, prec
       complex(kind = prec), intent(in) :: c
       class(cVector_t)    , intent(in) :: self
       class(cVector_t), allocatable :: Eout
     end function iface_Mult_2

     function iface_Mult_3(lhs, rhs) result(Eout)
       import :: cVector_t, rVector_t
       class(cVector_t), intent(in) :: lhs
       class(rVector_t), intent(in) :: rhs       
       class(cVector_t), allocatable :: Eout
     end function iface_Mult_3
     
     subroutine iface_Mult_s_1(self, c)
       import :: cVector_t, prec
       class(cVector_t)    , intent(inout) :: self
       complex(kind = prec), intent(in) :: c       
     end subroutine iface_Mult_s_1
     
     function iface_Div_1(lhs, rhs) result(Eout)
       import :: cVector_t
       class(cVector_t), intent(in) :: lhs
       class(cVector_t), intent(in) :: rhs
       class(cVector_t), allocatable :: Eout
     end function iface_Div_1

     function iface_Div_2(lhs, rhs) result(Eout)
       import :: cVector_t, rVector_t
       class(cVector_t), intent(in) :: lhs
       class(rVector_t), intent(in) :: rhs
       class(cVector_t), allocatable :: Eout
     end function iface_Div_2
     
     function iface_dotProd(lhs, rhs) result(r)
       import :: cVector_t, prec
       class(cVector_t), intent(in) :: lhs, rhs
       complex(kind = prec) :: r
     end function iface_dotProd

     subroutine iface_CopyFrom(self, rhs)
       import :: cVector_t
       class(cVector_t), intent(inout) :: self
       class(cVector_t), intent(in)    :: rhs
     end subroutine iface_CopyFrom
     
     !**
     ! Create a Vector object containing weights needed for
     ! interpolation of xyz component of obj1 to location.
     !*
     subroutine iface_InterpFunc(self, location, xyz, E)
       import :: cVector_t, prec
       class(cVector_t) , intent(in)  :: self
       real(kind = prec), intent(in)  :: location(3)
       character        , intent(in)  :: xyz
       class(cVector_t) , intent(out), allocatable :: E
     end subroutine iface_InterpFunc
     
  end interface
  
end module cVector
