module ModelOperator_File
   !
   use ModelOperator_MF
   !
   type, extends( ModelOperator_MF_t ) :: ModelOperator_File_t
       !   These are additional variables ...
       integer  :: n = 0
       !    note the full A matrix has a complex diagonal -- will
       !    only write out the real part, leaving the frequncy dependent
       !     imaginary part of the diagonal to be handled internally
       real( kind=prec ), allocatable, dimension(:,:) :: A
       character(:), allocatable                      :: Matrix_File_Name
       !
    contains
       !
       final :: ModelOperator_File_dtor
       !
       procedure, public :: amult
       !    this last procedure is unique to this class extension
       procedure, public :: setCoefficientMatrix
       !
       procedure, public :: print
       !
   end type ModelOperator_File_t
   !
   interface ModelOperator_File_t
       module procedure ModelOperator_File_ctor
   end interface ModelOperator_File_t
   !
contains
   !
   function ModelOperator_File_ctor( grid, fname ) result( self )
      implicit none
      !
      class( Grid3D_SG_t ), target, intent( in ) :: grid
      character(:), allocatable, intent( in )    :: fname ! I am assuming that the creator
                                                          ! for this extension can have an extra input argument
      !
      type( ModelOperator_File_t ) :: self
      !
      !write( *, * ) "Constructor ModelOperator_File"
      !
      ! Instantiation of the specific object MetricElements
      allocate( self%metric, source = MetricElements_CSG_t( grid ) )
      !
      call self%create( grid )
      call self%setCoefficientMatrix( fname )
      !
   end function ModelOperator_File_ctor
   !
   ! ModelOperator_MF destructor
   subroutine ModelOperator_File_dtor( self )
      implicit none
      !
      type( ModelOperator_File_t ), intent( inout ) :: self
      !
      !write( *, * ) "Destructor ModelOperator_File"
      !
      call self%deallocate()
     
   end subroutine ModelOperator_File_dtor
   !
   ! ***
   !
   subroutine setCoefficientMatrix( self, fname )
      implicit none
      !
      class( ModelOperator_File_t ), intent( inout ) :: self
      character(:), allocatable, intent( in )        :: fname
      !
      !   local variables
      integer :: inUnit
      !
      !    not sure we need a separate routine for this
      !    just do directly in ctor?
      self%Matrix_File_Name = fname
      !   not sure how you want to handle input unit --- just making a 
      !    number up here
      inUnit = 77
      open( inUnit, file = fname, form = 'unformatted', action = 'read' )
      !    first sequential record is matrix size
      read( inUnit ) self%n
      allocate( self%A( self%n, self%n ) )
      !    second sequential record is real part of coefficient matrix
      read( inUnit ) self%A
      !
      close( inUnit )
      !
   end subroutine setCoefficientMatrix
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
   subroutine Amult( self, omega, x, y, p_adjt )
      implicit none
      !
      class( ModelOperator_File_t ), intent( in ) :: self
      real( kind=prec ), intent( in ), optional   :: omega
      !   do these have to be abstract????
      class( cVector_t ), intent( in )            :: x
      class( cVector_t ), intent( inout )         :: y
      logical, intent( in ), optional             :: p_adjt
      !
      ! Local variables
      integer :: i, j
      complex( kind=prec ), allocatable, dimension(:) :: xVec, yVec
      real( kind=prec ), allocatable, dimension(:) :: Sigma, Vedge
      complex( kind=prec ) :: c
      logical :: adjt
      !
      if ( present( p_adjt ) ) then
          adjt = p_adjt
      else
          adjt = .false.
      end if
      !
      if ( adjt ) then
          c = -ONE_I * omega * ISIGN * MU_0
      else
          c = ONE_I * omega * ISIGN * MU_0
      end if
      !call self%print()
      write(*,*)  'omega = ',omega
      !
      select type( x )
      class is( cVector3D_SG_t )
          !
          select type( y )
          class is( cVector3D_SG_t )
             !
             if( .not. y%isAllocated ) then
                write( *, * ) 'ERROR: amult in   ModelOperator_MF'
                stop       'output vector y not allocated'
             endif
             write(*,*) 'grid dimensions on input to Amult'
             write(*,*) x%grid%nx, x%grid%ny, x%grid%nz, x%grid%nzAir
             call x%print(664,'input x-vector')
             !
             !   convert input cVector to column format 
             call x%getArray( xVec )
             !
            !write( 1111, * ) "xVec:"
            ! do i = 1, self%n
            !    if( abs(xVec( i )) .gt. R_TINY ) write( 1111, * ) xVec( i )
            ! enddo
             !
             !   allocate for product A*xVec
             allocate( yVec( self%n ) )
             !
             do i = 1, self%n
                yVec( i ) = C_ZERO
                do j = 1, self%n
                   yVec( i ) = yVec( i ) + complex( self%A( i, j ), 0.0 ) * xVec( j )
                enddo
             enddo
             !
             !   add in imaginary diagonal part of operator
             call self%Sigma_E%getArray(Sigma)
             call self%Metric%Vedge%getArray(Vedge)
             !call self%sigma_E%print(665,'edge conductivity')
             yVec = yVec + c * sigma * Vedge * xVec
             !yVec = yVec + c * xVec
             !
             !   convert result back to cVector`
             call y%setArray( yVec )
             !
          end select
          !
      class default
          write( *, * ) 'ERROR:ModelOperator_MF::amult:'
          stop        '         Incompatible input [x]. Exiting.'
          !
      end select
   end subroutine amult
   !
   subroutine print( self )
      implicit none
      !
      ! Arguments
      class( ModelOperator_File_t ), intent(in) :: self
      integer :: i, j
      !
      write( 1111, * ) "Matrix A, n=", self%n
      write( 1111, * ) "i", "j", "non_zero_value"
      do i = 1, self%n
         do j = 1, self%n
            if( self%A( i, j ) /= 0 ) then
               write( 1111, * ) i, j, ":[", self%A( i, j ), "]"
            endif
         enddo
      enddo
      !
   end subroutine print
   !
end module ModelOperator_File
