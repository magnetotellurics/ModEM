module ModelOperator_File
   !
   use Constants
   use Grid3D_SG
   use rScalar3D_SG
   use cScalar3D_SG
   use cVector3D_SG
   use rVector
   use rVector3D_SG   
   use MetricElements_CSG
   use ModelParameter
   use ModelOperator
   !
   type, extends( ModelOperator_MF_t ) :: ModelOperator_File_t
       !   These are additional variables ...
       integer  :: n
       !    note the full A matrix has a complex diagonal -- will
       !    only write out the real part, leaving the frequncy dependent
       !     imaginary part of the diagonal to be handled internally
       real(kind=prec), allocatable, dimension(:,:) :: A
       character(:), allocatable :: Matrix_File_Name
    contains
       !
       final :: ModelOperator_File_dtor
       !
       procedure, public :: Amult
       !    this last procedure is unique to this class extension
       procedure, public :: SetCoefficientMatrix
       !
   end type ModelOperator_File_t
   
   interface ModelOperator_File_t
	   module procedure ModelOperator_File_ctor
   end interface ModelOperator_File_t
   
contains
   
   !
   function ModelOperator_File_ctor( grid,fName ) result( self )
      implicit none
      !
      class( Grid3D_SG_t ), target, intent( in ) :: grid
      character(:), intent(in)  :: fName    !  I am assuming that the creator
                      !   for this extension can have an extra input argument
      !
      type( ModelOperator_File_t ) :: self
      !
      !write(*,*) "Constructor ModelOperator_MF"
      !
      ! Instantiation of the specific object MetricElements
      allocate( self%metric, source = MetricElements_CSG_t( grid ) )
      !
      call self%create( grid )
      call self%SetCoefficientMatrix(fName)
      !
   end function ModelOperator_File_ctor
   !
   ! ModelOperator_MF destructor
   subroutine ModelOperator_File_dtor( self )
      implicit none
      !
      type( ModelOperator_File_t ), intent( inout ) :: self
      !
      !write(*,*) "Destructor ModelOperator_MF_t"
      !
      call self%deallocate()
     
   end subroutine ModelOperator_File_dtor
   !
   ! ***
   !
   subroutine SetCoefficientMatrix(self,fName)
      type( ModelOperator_File_t ),intent(inout) :: self
      character(:), intent(in)  :: fName   

      !   local variables
      integer :: inUnit, n

      !    not sure we need a separate routine for this
      !    just do directly in ctor?
      self%Matrix_File_Name = fName
      !   not sure how you want to handle input unit --- just making a 
      !    number up here
      inUnit = 77
      open(inUnit,file = cfile, form = 'unformatted',action = 'read')
      !    first sequential record is matrix size
      read(inUnit) n
      allocate(self%A(n,n))
      !    second sequential record is real part of coefficient matrix
      read(inUnit) self%A
      close(inUnit)      

   end subroutine SetCoefficientMatrix
   !**
   ! SetEquations  -- will remain as it was, so that
   !          all of the other routines in ModelOperator still work
   !*
   !**
   ! AMult -- this is the only overloaded procedure
   ! This does the matrix-vector multiply (A+iwB)inE = outE
   ! parameterized by input frequency omega, but using real matrix self%A
   !    loaded from file
   !*
   subroutine Amult(self, omega, x, y, p_adjt)
      class(ModelOperator_File_t), intent(in)        :: self
      real( kind = prec ), intent(in), optional      :: omega
      !   do these have to be abstract????
      class(cVector_t)         , intent(in)                :: x
      class(cVector_t)         , intent(inout)            :: y
      logical                      , intent(in), optional :: p_adjt

      ! Local variables
      integer :: ix, iy, iz
      complex(kind=prec), allocatable, dimension(:) :: xVec, yVec
      complex( kind=prec ) :: c
      logical :: adjt

      if ( present( p_adjt ) ) then
          adjt = p_adjt
      else
          adjt = .false.
      end if
      !
      if (adjt) then
          c = -C_ONE * omega * ISIGN * MU_0
      else
          c = C_ONE * omega * ISIGN * MU_0
      end if
      !
      select type(x)
      class is(cVector3D_SG_t)
          if(.not.y%isAllocated) then
               write(*,*) 'ERROR: Amult in   ModelOperator_MF'
               write(*,*) 'output vector y not allocated'
               STOP
          endif
          select type(y)
          class is(cVector3D_SG_t)

             !   convert input cVector to column format 
             call x%getArrayCVector3D_SG(xVec)
             !   allocate for product A*xVec
             allocate(yVec(n))
             !   direct matrix-vector multiplication`
             yVec = A*xVec
             !   add in imaginary diagonal part of operator
             yVec = yVec + c*xVec
             !   convert result back to cVector`
             call y%setArrayCVector3D_SG(yVec)
             ! Finally multiply by VEdge (in place)
             call y%mults3(self%metric%Vedge)

          end select
      class default
          write(*, *) 'ERROR:ModelOperator_MF::Amult:'
          write(*, *) '         Incompatible input [x]. Exiting.'
          
          STOP
      end select

   end subroutine Amult
   
end module ModelOperator_File
