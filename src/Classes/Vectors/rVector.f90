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
     procedure(iface_Read) , deferred, public :: Read
     procedure(iface_Write), deferred, public :: Write
     
     !**
     ! Boundary operations
     !*
     procedure(iface_SetAllBoundary), deferred :: SetAllBoundary
     procedure(iface_SetOneBoundary), deferred :: SetOneBoundary
     procedure(iface_SetAllInterior), deferred :: SetAllInterior
     procedure(iface_IntBdryIndices), deferred :: IntBdryIndices
          
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
     generic :: Mult => Mult_1, Mult_2
     generic :: operator(*) => Mult_1, Mult_2
     
     procedure(iface_Div_1) , deferred, public :: div_1
     generic :: Div => Div_1
     generic :: operator(/) => Div_1
     
     procedure(iface_dotProd), deferred, public :: dotProd
     generic :: operator(.dot.) => dotProd
     
     procedure(iface_diagMult), deferred, public :: diagMult

     procedure(iface_CopyFrom), deferred, public :: CopyFrom
     generic :: assignment(=) => CopyFrom

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
     ! SetAllBoundary
     !*
     subroutine iface_SetAllBoundary(self, c_in)
       import :: rVector_t, prec
       class(rVector_t) , intent(inout) :: self
       real(kind = prec), intent(in)    :: c_in
     end subroutine iface_SetAllBoundary

     !**
     ! SetOneBoundary
     !*
     subroutine iface_SetOneBoundary(self, bdry, c, int_only)
       import :: rVector_t, prec
       class(rVector_t) , intent(inout) :: self
       character(*)     , intent(in)    :: bdry
       real(kind = prec), intent(in)    :: c
       logical, optional, intent(in)    :: int_only
     end subroutine iface_SetOneBoundary

     !**
     ! SetAllInterior
     !*
     subroutine iface_SetAllInterior(self, c_in)
       import :: rVector_t, prec
       class(rVector_t) , intent(inout) :: self
       real(kind = prec), intent(in)    :: c_in
     end subroutine iface_SetAllInterior

     !**
     ! IntBdryIndices
     !*
     subroutine iface_IntBdryIndices(self, ind_i, ind_b)
       import :: rVector_t
       class(rVector_t)    , intent(in)  :: self
       integer, allocatable, intent(out) :: ind_i(:), ind_b(:)
     end subroutine Iface_IntBdryIndices

     !**
     ! Data access
     !*

     !**
     ! Length
     !*
     function iface_Length(self) result(n)
       import :: rVector_t
       class(rVector_t), intent(in) :: self
       integer :: n
     end function iface_Length

     !**
     ! GetArray
     !*
     subroutine iface_GetArray(self, v)
       import :: rVector_t, prec
       class(rVector_t), intent(in)  :: self
       real(kind = prec), allocatable, intent(out) :: v(:)
     end subroutine Iface_GetArray

     !**
     ! SetArray
     !*
     subroutine iface_SetArray(self, v)
       import :: rVector_t, prec
       class(rVector_t) , intent(inout) :: self
       real(kind = prec), intent(in)    :: v(:)
     end subroutine iface_SetArray
     
     !**
     ! Arithmetic/algebraic operations
     !*
     
     !**
     ! Zeros
     !*
     subroutine iface_Zeros(self)
       import :: rVector_t
       class(rVector_t), intent(inout) :: self
     end subroutine iface_Zeros

     function iface_Add_1(lhs, rhs) result(Eout)
       import :: rVector_t
       class(rVector_t), intent(in)  :: lhs, rhs
       class(rVector_t), allocatable :: Eout
     end function iface_Add_1

     function iface_Sub_1(lhs, rhs) result(Eout)
       import :: rVector_t
       class(rVector_t), intent(in) :: lhs, rhs
       class(rVector_t), allocatable :: Eout
     end function iface_Sub_1

     function iface_Mult_1(lhs, rhs) result(Eout)
       import :: rVector_t
       class(rVector_t), intent(in) :: lhs, rhs
       class(rVector_t), allocatable :: Eout
     end function iface_Mult_1
     
     function iface_Mult_2(c, self) result(Eout)
       import :: rVector_t, prec
       real(kind = prec), intent(in) :: c
       class(rVector_t) , intent(in) :: self
       class(rVector_t), allocatable :: Eout
     end function iface_Mult_2
     
     function iface_Div_1(lhs, rhs) result(Eout)
       import :: rVector_t
       class(rVector_t), intent(in) :: lhs, rhs
       class(rVector_t), allocatable :: Eout
     end function iface_Div_1

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

     subroutine iface_CopyFrom(self, rhs)
       import :: rVector_t
       class(rVector_t), intent(inout) :: self
       class(rVector_t), intent(in)    :: rhs
     end subroutine iface_CopyFrom

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
