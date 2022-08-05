Module Esoln2DTM
    !
    use Constants
    use Grid2D
    !
    type :: Esoln2DTM_t
      !
      complex( kind=prec ), allocatable, dimension(:,:) :: Ey, Ez
      class( Grid2D_t ), pointer :: grid
      logical :: is_allocated
      !
     contains
         !
         final :: Esoln2DTM_dtor
         !
         procedure, public :: getArray => getArrayEsoln2DTM
         procedure, public :: setArray => setArrayEsoln2DTM
         procedure, public :: getBoundary => getBoundaryEsoln2DTM
         procedure, public :: nEdges => nEdgesEsoln2DTM
         !
    end type Esoln2DTM_t
    !
    interface Esoln2DTM_t
     module procedure Esoln2DTM_ctor
    end interface Esoln2DTM_t
	!
	contains
	!**
	! Class constructor
	!*
	function Esoln2DTM_ctor( grid ) result( self )
		!    pass the model grid, create and return object
		class( Grid2D_t ), target, intent(in) :: grid
		type( Esoln2DTM_t ) :: self
		!    local variables
		integer :: ny, nz, status
		!
		self%is_allocated = .FALSE.
		ny = grid%Ny
		nz = grid%Nz
		self%grid => grid    ! Set pointer to grid

		allocate(self%Ey(ny,nz+1),STAT=status)
		if(status /= 0) then
			write(*,*) "Error:Esoln2DTM_ctor: allocation failed"
		endif
		!
		allocate(self%Ez(ny+1,nz),STAT=status)
		!
		if(status /= 0) then
			write(*,*) "Error:Esoln2DTM_ctor: allocation failed"
		endif
		self%is_allocated = .true.
		self%Ey = C_ZERO
		self%Ez = C_ZERO
	end function Esoln2DTM_ctor
	!
	!     Class destructor
	!
	subroutine Esoln2DTM_dtor(self)
		type(Esoln2DTM_T)    ::    self

		if(self%is_allocated) then
			deallocate(self%Ey)
			deallocate(self%Ez)
		endif
		self%grid => null()
	end subroutine Esoln2DTM_dtor
	!
	!*************************
	!
	function nEdgesEsoln2DTM(self) result(n)
		class(Esoln2DTM_t), intent(in)    :: self
		integer :: n
		!    local variables
		integer    :: nYedge, nZedge

		call self%grid%NumberOfEdges(nYedge,nZedge)

		n = nYedge+nZedge
	end function nEdgesEsoln2DTM
	!
	!*************************
	!
    function GetArrayEsoln2DTM(self) result(x)
		! puts Ey and Ez into a column vector, with interleaved
		!    vertical grid columns
		implicit none
		class(Esoln2DTM_t), intent(in) :: self
		complex(kind=prec),allocatable    :: x(:)
		! local variables
		integer :: j,k,jk,ny,nz, nEdge
		!
		if(.NOT.self%is_allocated) then
			stop "Error in GetArray: Esoln2DTM object not allocated"
		endif

		allocate(x(self%nEdges())) 

		ny = self%grid%Ny
		nz = self%grid%Nz
		!    insert Ez, Ey into a column vector, interleaving vertical columns of edges
		!    maybe more efficient to use array operations?    
		jk = 0
		!    first Ez
		do j = 1,ny+1
			do k = 1,nz
				jk = jk+1
				x(jk) = self%Ez(j,k)
			enddo
			jk = jk+Nz+1
		enddo
		!    then Ey
		jk = 0
		do j = 1,ny
			jk = jk+Nz
			do k = 1,nz+1
				jk = jk+1
				x(jk) = self%Ey(j,k)
			enddo
		enddo
    end function getArrayEsoln2DTM
    !
    !*************************
    !
    subroutine setArrayEsoln2DTM(self,x)
		! extracts Ey and Ez and puts into a column vector, with interleaved
		!    vertical grid columns
		implicit none
		class(Esoln2DTM_t),intent(inout) :: self
		complex(kind=prec)    :: x(:)
		! local variables
		integer :: j,k,ny,nz,jk
		!
		if(.NOT.self%is_allocated) then
			stop "Error in GetArray: Esoln2DTM object not allocated"
		endif
		!
		ny = self%grid%Ny
		nz = self%grid%Nz
		!    insert Ez, Ey, using    interleaved vertical columns of edges
		!    in input vector x
		jk = 0
		!    first Ez
		do j = 1,ny+1
			do k = 1,nz
				jk = jk+1
				self%Ez(j,k) = x(jk)
			enddo
			jk = jk+nz+1
		enddo
		!    then Ey
		jk = 0
		do j = 1,ny
			jk = jk+nz
			do k = 1,nz+1
				jk = jk+1
				self%Ey(j,k) = x(jk)
			enddo
		enddo
		!
    end subroutine setArrayEsoln2DTM
    !
    !*************************
    !
    function getBoundaryEsoln2DTM(self) result(x)
		!    comparable to GetArray, but only boundary values are inserted into
		!    the output vector
		implicit none
		class(Esoln2DTM_t), intent(in) :: self
		complex(kind=prec), allocatable    :: x(:)
		! local variables
		integer :: j,k,ny,nz,jkLeft,jkRight,jkTop,jkBottom, nEdge
		!
		!
		if(.NOT.self%is_allocated) then
			stop "Error in GetArray: Esoln2DTM object not allocated"
		endif

		allocate(x(self%nEdges()))
		x = C_ZERO

		ny = self%grid%Ny
		nz = self%grid%Nz
		!    first Ez on left/right boundaries
		jkLeft = 0
		jkRight = ny*(2*nz+1)
        do k = 1,nz
			jkLeft = jkLeft+1
			jkRight = jkRight+1
			x(jkLeft) = self%Ez(1,k)
			x(jkRight) = self%Ez(ny+1,k)
        enddo
        jkTop = nz+1
        jkBottom = 2*nz+1
        do j = 1,ny
			x(jkTop) = self%Ey(j,1)
			x(jkBottom) = self%Ey(j,nz+1)
			jkTop = jkTop+2*nz+1
			jkBottom = jkBottom+2*nz+1
        enddo
        !    then Ey on top/bottom boundaries
    end function getBoundaryEsoln2DTM
    !
end module Esoln2DTM
