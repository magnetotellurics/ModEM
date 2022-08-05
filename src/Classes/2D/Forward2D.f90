!    Forward  modeling class for 2D, used for creating boundary  data for 
!     for 3D forward problem.    I am creating a single explicit claass; as
!     as we develop alternative 2D Fwd instances
!     we might create an abstract base class and make this an extension, or
!     just extend this class, overloading setup and solve routines
module Forward2D
    !
	use Constants
    use Grid2D
    use ModelParameter2D
    use Esoln2DTM
    !
    type  :: Forward2D_t
        class(ModelParameter2D_t), pointer  :: m
        !     note that the model parameter object contains two grids --
        !     one that defines the numerical grid (which solutions will be
        !     computed on, and one that defines the model parameter grid
        real( kind=prec )                            :: omega
        !     I am going code this for the 2D TM, E field case -- and use
        !     property names approprate for this case.  The same variables
        !     could be used for other cases, but names won"t be intuitive
        !    Following are variables to save curl-curl operator in termns of
        !     non-zero diagonals
        integer :: nDiag
        integer, allocatable, dimension(:) :: Diag
        real( kind=prec ), allocatable, dimension(:,:) :: CC
        !
        contains
            !
            procedure, public  :: createCC => createCCForward2D
            procedure, public  :: setCond => setCondForward2D
            procedure, public  :: setFreq => setFreqForward2D
            procedure, public  :: solve => solveForward2D
            !
    end type Forward2D_t

    interface Forward2D_t
        module procedure Forward2D_ctor
    end interface Forward2D_t

 contains

    !**
    ! Class constructor
    !*
    function Forward2D_ctor(m) result( self )
        !     pass the model parameter, create and return object
        class(ModelParameter2D_t), target , intent(in) :: m
        ! Local variables
        type(Forward2D_t) :: self

        self%nDiag = 9
        self%m => m     ! Set poinbter to model parameter

        call self%createCC()    !  create curl-curl operator
        !
        call self%setCond( m )
        !
    end function Forward2D_ctor
    !
    !**********
    !
    subroutine createCCForward2D( self )
        class(Forward2D_t), intent(inout) :: self
        ! local variables
        integer :: nYedge, nZedge, Ny, Nz, j, k, jk, nEdge
        real( kind=prec ), allocatable, dimension(:) :: dY, dZ, delY, delZ

        !     copy some things from the grid to make code easier to read
        Ny = self%m%grid%ny
        Nz = self%m%grid%nz
        allocate(dY(Ny))
        allocate(dZ(Nz))
        allocate(delY(Ny+1))
        allocate(delZ(Nz+1))
        dY = self%m%grid%dY
        delY = self%m%grid%delY
        dZ  = self%m%grid%dZ
        delZ = self%m%grid%delZ

        call self%m%grid%NumberOfEdges( nYedge, nZedge )
        nEdge = nYedge+nZedge
        allocate(self%Diag(self%nDiag))
        allocate(self%CC(self%nDiag,nEdge))
        !    set diagonal indices
        self%Diag(1) = 0
        self%Diag(2) = -1
        self%diag(3) = 1
        self%diag(4) = -Nz
        self%diag(5) = -Nz-1
        self%diag(6) = Nz
        self%diag(7) = Nz+1
        self%diag(8) = -2*Nz-1
        self%diag(9) = 2*Nz+1
        !    set values on corresponding diagonals for CC operator
        !    Y edges
        jk = -1
        self%CC = R_ZERO  !    initialize to zero
        self%CC(1,:) = ONE  !  except for main diagonal
        !    jk represents the edge numbering -- column # in banded storage is jk+diag!
        jk = -1
        do j=1,Ny
            jk = jk+Nz+2     !    Nz to skip all vertical edges, 2 to skip top and bottom bdry
            do k = 2,Nz !  loop limits are set to only include interior edges
                jk = jk+1
                self%CC(1,jk) = ONE/(dZ(k)*delZ(k))+ONE/(dZ(k-1)*delZ(k))
                self%CC(2,jk+self%diag(2)) = -ONE/(dZ(k-1)*delZ(k))  
                self%CC(3,jk+self%diag(3)) = -ONE/(dZ(k)*delZ(k))  
                self%CC(4,jk+self%diag(4)) =  -ONE/(dY(j)*delZ(k))
                self%CC(5,jk+self%diag(5)) =  ONE/(dY(j)*delZ(k))
                self%CC(6,jk+self%diag(6)) =  -ONE/(dY(j)*delZ(k))
                self%CC(7,jk+self%diag(7)) =  ONE/(dY(j)*delZ(k))
            enddo
        enddo
        !    Z edges
        jk = 2*Nz+1  !    skip first z column + y column
        do j=2,Ny  !  loop limits are set to only include interior edges
            do k = 1,Nz
                jk = jk+1
                self%CC(1,jk) = ONE/(dY(j)*delY(j))+ONE/(dY(j-1)*delY(j))
                self%CC(4,jk+self%diag(4)) = -ONE/(dZ(k)*delY(j))  
                self%CC(5,jk+self%diag(5)) = ONE/(dZ(k)*delY(j))  
                self%CC(6,jk+self%diag(6)) = -ONE/(dZ(k)*delY(j))  
                self%CC(7,jk+self%diag(7)) = ONE/(dZ(k)*delY(j))  
                self%CC(8,jk+self%diag(8)) =  -ONE/(dY(j-1)*delY(j))
                self%CC(9,jk+self%diag(9)) =  -ONE/(dY(j)*delY(j))
            enddo
            jk = jk + Nz+1    ! skip next y column
        enddo
        !
        deallocate(dY)
        deallocate(delY)
        deallocate(dZ)
        deallocate(delZ)
        !
    end subroutine createCCForward2D
    !**********
    !
    subroutine setCondForward2D(self,m)
        !    set pointer to possibly different model parameter object
        class(Forward2D_t), intent(inout) :: self
        type(ModelParameter2D_t), target , intent(in) :: m

        self%m => m
    end subroutine setCondForward2D
    !
    !**********
    !
    subroutine setFreqForward2D(self,omega)
        !    set frequency for  to possibly different model parameter object
        class(Forward2D_t), intent(inout) :: self
        real( kind=prec ), intent(in) :: omega

        self%omega = omega

    end subroutine setFreqForward2D
    !
    !**********
    !
    subroutine solveForward2D( self, E2D )
        implicit none
        !    solve 2D equations for model parameter m, frequncy omega
        !    E2D: on input contains boundary data on edges;
        !          on output contains solution, including input boundary data
        class(Forward2D_t), intent(in)    :: self
        type(Esoln2DTM_t), intent(inout) :: E2D

        !  local variables
        integer ::  ku, kl, nYedge, nZedge, lda, n, info, i, j, k, nRHS, ny, nz, iRow
        !    arrays to allocate for lapack solver
        integer, allocatable, dimension(:) :: ipiv
        complex( kind=prec ), dimension(:,:), allocatable :: A
        complex( kind=prec ), dimension(:), allocatable    :: b
        complex:: cFac
        real( kind=prec ), allocatable, dimension(:) :: sigmaEdge
        character*1 :: TRANS
        integer debugUnit, error_code
        character*80  debugFile

        if( .NOT. E2D%is_allocated) then
            stop "Error in Forward2D%Solve: Esoln2D object not allocated"
        endif
        !
        !    temporary code for debugging
        !debugUnit = 1
        !debugFile = "Test/Solver2D_test"
        !open(debugUnit,file = debugFile, form = "formatted")
        
        !  For now I am just hard-coding for this case: analytic layered solution
        !  This subroutine will set up equations and solve
        !  Might break this into pieces, to better reuse code

        !  The system matrix is banded with two subdiagonals, and two superdiagonals
        TRANS = "N"    !  character string that controls ZGBTRS:
                          !  N = solve base system Ax = b
                          !  T = solve transposed system A"x = b
                          !  C = solve Hermitian transpose A^* = b    
        ku = maxval(self%Diag)
        kl = -minval(self%Diag)
        lda = 2*kl+ku+1     !    storage required for banded matrix
                      !    factorization using ZGBTRF
        print*,"ku,kl,lda",ku,kl,lda

        call self%m%grid%NumberOfEdges( nYedge, nZedge )
        nRHS = 1
        n = nYedge+nZedge !    number of unknowns to solve for
        Ny = self%m%grid%ny
        Nz = self%m%grid%nz
        print*,"nYedge,nZedge,n",nYedge,nZedge,n
        !    NOTE: we assume that the last layer (which will have finite thickness
        !     in the ParamGrid) actually extends to infinity -- thus this thickness
        !     is not really used in computations, only the bottom-layer conductivity

        !    allocate arrays for computations: coefficient matrix, rhs, solution
        !     E interpolated to numerical grid
        allocate( A( lda, n ), stat=error_code )
        !
        allocate(b(n))
        allocate( sigmaEdge(n), stat=error_code )
        !
        allocate( ipiv(n), stat=error_code )
        !
        !
        !    copy curl-curl coefficients stored in CC_y and CC_z into band-storage array
        A = C_ZERO
        !  copy CC into appropriate rows of A
        do i = 1,self%nDiag
            ! kl + ku + 1 is row index for main diagona
            iRow = kl+ku+1-self%Diag(i)
            A(iRow,:) = self%CC(i,:)
        enddo
        !    add 1i*omega*mu*sigma
        sigmaEdge = self%m%PDEmapping()    !  make this routine return zero on boundary edges
        iRow = kl+ku+1     !    row for diagonal in band matrix storage
        cFac = ISIGN*ONE_I*mu_0*self%omega
        A(iRow,:) = A(iRow,:) + cFac*sigmaEdge(:)
        !    set up RHS -- lets use a procedure in E2D
        b = E2D%GetBoundary()     !    this just extracts what is on the boundary
                                     ! in the Ez, Ey arrays
        !    solve equations using lapack
        !     LU decomposition of complex banded matrix
        call zgbtrf(n,n,ku,kl,A,lda,ipiv,info)
        call zgbtrs(TRANS,n,kl,ku,nRHS,A,lda,ipiv,b,n,info)
        !
        !    put the solution into the E2D object]
        call E2D%SetArray(b)
        !
        !    deallocate local arrays (is this needed???)
        deallocate(A)
        deallocate(b)
        deallocate(ipiv)
        deallocate(sigmaEdge)
        !
    end subroutine solveForward2D

end Module Forward2D
