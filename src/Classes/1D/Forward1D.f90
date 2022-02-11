!   Forward  modeling class for 1D, used for creating boundary  data for 
!    for 3D forward problem.   I am creating a single explicit claass; as
!    as we develop alternative 1D Fwd instances (e.g., for CSEM, anisotropy)
!    we might create an abstract base class and make this an extension, or
!    just extend this class, overloading setup and solve routines
module Forward1D
   use Constants
   use Grid1D
   use ModelParameter1D

  implicit none
   public  :: Forward1D_t

   type  :: Forward1D_t
      class(ModelParameter1D_t), pointer  :: m => null()
      !    note that the model parameter object contains two grids --
      !    one that defines the numerical grid (which solutions will be
      !    computed on, and one that defines the model parameter grid
      real(kind=prec)                     :: omega

   contains

      procedure, public  :: SetCond
      procedure, public  :: SetFreq
      procedure, public  :: Solve
      procedure, public  :: GridedSolution

   end type Forward1D_t

   interface Forward1D_t
      module procedure Forward1D_ctor
   end interface Forward1D_t

 contains

   !**
   ! Class constructor
   !*
   function Forward1D_ctor(m) result(fwd1d)
      !    pass the model parameter, create and return object
      class(ModelParameter1D_t), target , intent(in) :: m
      ! Local variables
      type(Forward1D_t) :: fwd1d

      fwd1d%m => m
   end function Forward1D_ctor
   !
   !**********
   !
   subroutine SetCond(self,m)
      !   set pointer to possibly different model parameter object
      class(Forward1D_t), intent(inout) :: self
      type(ModelParameter1D_t), target , intent(in) :: m

      self%m => m
   end subroutine SetCond
   !
   !**********
   !
   subroutine SetFreq(self,omega)
      !   set frequency for  to possibly different model parameter object
      class(Forward1D_t), intent(inout) :: self
      real(kind=prec), intent(in) :: omega

      self%omega = omega

   end subroutine SetFreq
   !
   !**********
   !
   subroutine Solve( self, E1D )
      !   solve 1D equations for model parameter m, frequncy omega, returning
      !   solution in complex array E1D on grid defined by m%grid (pointer to
      !   numerical grid)
      class(Forward1D_t), intent(in) :: self
      complex(kind=prec), intent(out), dimension(:) :: E1D

      !  local variables
      integer ::  ku,kl,nlayer,n,nz,lda,info,i,j,nRHS
      integer, allocatable, dimension(:)  :: ipiv
      complex(kind=prec), allocatable, dimension(:,:)  :: A
      complex(kind=prec), allocatable, dimension(:) :: beta
      real(kind=prec), allocatable, dimension(:) :: h
      complex(kind=prec), allocatable, dimension(:) :: b
      character*1 :: TRANS
      !integer *  !   for debugging
      character*80  debugFile

      !   temporary code for debugging
      !* = 1
      !debugFile = 'Test/Solver1D_Coefficient_Matrix'
      !open(*,file = debugFile, form = 'formatted')
      
      !  For now I am just hard-coding for this case: analytic layered solution
      !  This subroutine will set up equations and solve
      !  Might break this into pieces, to better reuse code

      !  The system matrix is banded with two subdiagonals, and two superdiagonals
      TRANS = 'N'   !  character string that controls ZGBTRS:
                    !  N = solve base system Ax = b
                    !  T = solve transposed system A'x = b
                    !  C = solve Hermitian transpose A^* = b   
      ku = 2
      kl = 2
      lda = 2*kl+ku+1    !   storage required for banded matrix
                 !   factorization using ZGBTRF
      nlayer = self%m%ParamGrid%nz 
      nz = self%m%grid%nz
      n = 2*nlayer - 1 !   number of unknowns to solve for
      !   NOTE: we assume that the last layer (which will have finite thickness
      !    in the ParamGrid) actually extends to infinity -- thus this thickness
      !    is not really used in computations, only the bottom-layer conductivity

      !   allocate arrays for computations: coefficient matrix, rhs, solution
      !    E interpolated to numerical grid
      allocate(A(lda,n))
      allocate(b(n))
      allocate(ipiv(n))
      allocate(beta(nlayer))
      allocate(h(nlayer))
      !allocate(E1D(self%m%grid%nz))   !  probably should be allocated before
                                      !   calling
      !   compute complex beta coefcients, for each layer in ParamGrid
      !   I am going to assume that paramType for 1D model parameter is always
      !    LINEAR -- need to make sure that this is true when creating
      !     or could add check/transformation here ...
      do i = 1,nlayer
          beta(i) = sqrt(ISIGN*ONE_I*MU_0*self%omega*self%m%CellCond(i))
          h(i) = self%m%ParamGrid%dz(i)    !   layer thickness -- last is not used
      enddo
         
      !   create cofficient matrix -- layer by layer -- first and last are unique
      A = C_ZERO
      do i = 2,nlayer-1
         A(lda-3,2*i-1) = C_ONE
         A(lda-4,2*i) = C_ONE
         A(lda-2,2*i-1) = -beta(i)
         A(lda-3,2*i) = beta(i)
         A(lda-1,2*i-1) = -exp(-beta(i)*h(i))
         A(lda-2,2*i) = -exp(beta(i)*h(i))
         A(lda,2*i-1) = beta(i)*exp(-beta(i)*h(i))
         A(lda-1,2*i) = -beta(i)*exp(beta(i)*h(i))
      enddo

      !   special case for first layer
      A(lda-2,1) = C_ONE
      A(lda-3,2) = C_ONE
      A(lda-1,1) = -exp(-beta(1)*h(1))
      A(lda-2,2) = -exp(beta(1)*h(1))
      A(lda,1) = beta(1)*exp(-beta(1)*h(1))
      A(lda-1,2) = -beta(1)*exp(beta(1)*h(1))

      !  special case for last layer
      A(lda-3,2*nlayer-1) = C_ONE
      A(lda-2,2*nlayer-1) = -beta(i)

      !   some debugging output
      !write(*,*) 'nlayer,ku,kl,lda,n'
      !write(*,*) nlayer,ku,kl,lda,n
      !write(*,*) 'omega',self%omega
      !write(*,*) 'Coefficient Matrix'
      !do i = 1,n
         !write(*,'(7(e11.3,e11.3,1x))') (A(j,i),j=1,lda)
      !enddo

      !   create rhs
      b = C_ZERO
      b(1) = C_ONE
      nRHS = ONE

      !   solve equations using lapack
      !    LU decomposition of complex banded matrix
      call zgbtrf(n,n,ku,kl,A,lda,ipiv,info)
      call zgbtrs(TRANS,n,kl,ku,nRHS,A,lda,ipiv,b,n,info)

      !write(*,*) 'solution coefficients'
      !write(*,*) b
      !close(*)

      !  calculate E1D on model grid
      call self%GridedSolution( b, beta, h, E1D )

      !   deallocate local arrays (is this needed???)
      deallocate(A)
      deallocate(b)
      deallocate(beta)
      deallocate(h)
      deallocate(ipiv)
   end subroutine Solve
   !
   !*********
   !
   subroutine GridedSolution( self, b, beta, h, E1D )
      !   given coefficients b, complex wavenumbers beta, layer thickness h
      !   compute electric field solution E1D evaluated at grid levels
      !   specfied in z , and normalized so that E1D(1) = 1.0
      class(Forward1D_t), intent(in) :: self
      complex(kind=prec), intent(in), dimension(:) :: beta,b
      real(kind=prec), intent(in), dimension(:) :: h
      complex(kind=prec), dimension(:) :: E1D
      ! local variables
      integer :: n,nlayer,nz,i,j,layerInd
      real(kind=prec), allocatable, dimension(:)  :: interfaceDepth, z
      complex(kind=prec)  :: d, u,scaleFactor
      real(kind=prec)  :: zLayer

      n = size(b)
      nlayer = size(beta)
      !
      z = self%m%grid%zEdge
      nz = size(z)
      write(*,*) 'length of E1D = ', nz
      !
      allocate(interfaceDepth(nlayer))
      interfaceDepth = R_ZERO
      do i = 2,nlayer
         interfaceDepth(i) = interfaceDepth(i-1)+h(i-1)
      enddo
      !   assign layer index for each depth z(i)
      do i = 1,nz
          layerInd = 0
          do j = nlayer,1,-1
             if(interfaceDepth(j).lt.z(i)) then 
                !   this should be deepest interface above z(i)
                layerInd = j
                exit
             endif
          enddo
          if(layerInd.gt. 0) then
              zLayer = z(i)-interFaceDepth(layerInd)
          else
              zLayer = z(i)
          endif
          if(layerInd .gt.0) then
             ! downgoing and upgoing coeffcients for this layer
             d = b(2*layerInd-1)
             if(layerInd .lt. nlayer) then
                u = b(2*layerInd)
             else
                u = C_ZERO
             endif
             E1D(i) =  d*exp(-beta(layerInd)*zLayer)+ &
               u*exp(beta(layerInd)*zLayer)
          else  !  this z is in the air -- user linear expression
             E1D(i) = C_ONE + beta(1)*(b(2)-b(1))*z(i)
          endif
       enddo
       scaleFactor = C_ONE/E1D(1)
       E1D(:) = E1D(:)*scaleFactor


       !if( allocated( interfaceDepth) ) deallocate( interfaceDepth )

   end subroutine GridedSolution

end Module Forward1D
