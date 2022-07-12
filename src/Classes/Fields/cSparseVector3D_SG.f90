module cSparseVector3D_SG  
    !
    use Constants
    use cVector3D_SG
    !
    type :: cSparsevector3D_SG_t
        !
        character( len=4 ) :: gridType
        !
        integer  :: nCoeff
        !
        integer , allocatable, dimension(:) :: i, j, k, xyz
        !
        complex ( kind=prec ), allocatable, dimension(:) :: c
        !
        logical :: is_allocated
        !
        contains
            !
            final :: cSparsevector3D_SG_dtor
            !
    end type cSparsevector3D_SG_t
    !
    ! Constructors for Scalar3d_csg_real_t
    interface cSparsevector3D_SG_t
        module procedure cSparsevector3D_SG_ctor
    end interface cSparsevector3D_SG_t
    !
contains
    function cSparsevector3D_SG_ctor() result ( self )
        implicit none
        !
        type( cSparsevector3D_SG_t ) :: self
        !
        !write(*,*) "Constructor cSparsevector3D_SG"
        !
        self%gridType = ""
        self%nCoeff = 0
        self%is_allocated = .FALSE.
        !
    end function cSparsevector3D_SG_ctor
    !
    subroutine cSparsevector3D_SG_dtor( self )
        implicit none
        !
        type( cSparsevector3D_SG_t ), intent( inout ) :: self
        !
        !write(*,*) "Destructor cSparsevector3D_SG_t:"
        !
        self%gridType = ""
        self%nCoeff = 0
        self%is_allocated = .FALSE.
        !
        if( allocated( self%i ) ) deallocate( self%i )
        if( allocated( self%j ) ) deallocate( self%j )
        if( allocated( self%k ) ) deallocate( self%k )
        !
        if( allocated( self%xyz ) ) deallocate( self%xyz )
        !
        if( allocated( self%c ) ) deallocate( self%c )
        !
    end subroutine cSparsevector3D_SG_dtor
    !
    function dotProdSparse( self, cvector ) result( cvalue )
        implicit none
        !
        type( cSparsevector3D_SG_t ), intent( in ) :: self
        type( cVector3D_SG_t ), intent( in )       :: cvector
        !
        complex( kind=prec )                       :: cvalue
        !
        integer :: i, xi, yi, zi
        !
        cvalue = C_ZERO
        !
        if( .NOT. self%is_allocated ) then
            stop "SELF not is_allocated yet for dotProdSparse"
        endif
        !
        if( .NOT. cvector%is_allocated ) then
            stop "RHS not is_allocated yet for dotProdSparse"
        endif
        !
        if ( self%gridType /= cvector%gridType ) then
            stop "dotProdSparse: not compatible usage for dotProdSparse"
        endif
        !
        ! sum over  non-zero terms in sparse vector (conjugate sparse)
        ! (need to check xyz the component)
        ! Remember, xyz = 1,2,3 refers to x, y or z components
        do i = 1, self%nCoeff
            !
            ! generic test for both edge and face (all the components)
            if( ( self%i(i) .LE. cvector%grid%nx + 1 ) .OR. &
                ( self%j(i) .LE. cvector%grid%ny + 1 ) .OR. &
                ( self%k(i) .LE. cvector%grid%nz + 1 ) ) then
                !
                ! dealing with x-components
                if( self%xyz(i) == 1 ) then
                    xi = self%i(i)
                    yi = self%j(i)
                    zi = self%k(i)
                    cvalue = cvalue + conjg( self%c(i) ) * cvector%x( xi, yi, zi )
                !
                ! dealing with y-component
                else if( self%xyz(i) == 2 ) then
                    xi = self%i(i)
                    yi = self%j(i)
                    zi = self%k(i)
                    cvalue = cvalue + conjg( self%c(i) ) * cvector%y( xi, yi, zi )
                !
                ! dealing with z-component
                else if( self%xyz(i) == 3 ) then
                    xi = self%i(i)
                    yi = self%j(i)
                    zi = self%k(i)
                    cvalue = cvalue + conjg( self%c(i) ) * cvector%z( xi, yi, zi )
                end if
            !
            else
                stop "IJK out of bounds for dotProdSparse"
            !
            endif
            !
        enddo
        !
    end function dotProdSparse
    !
    subroutine getFullVector( self, cvector )
        implicit none
        !
        type( cSparsevector3D_SG_t ), intent( in ) :: self
        type( cVector3D_SG_t ), intent( inout )    :: cvector
        !
		!
    end subroutine getFullVector
    !
    subroutine full2Sparse( self, cvector )
        implicit none
        !
        type( cSparsevector3D_SG_t ), intent( inout ) :: self
        type( cVector3D_SG_t ), intent( in )          :: cvector
        !
        integer, allocatable, dimension(:,:,:)  :: Ix, Jx, Kx, XYZ1
        integer, allocatable, dimension(:,:,:)  :: Iy, Jy, Ky, XYZ2
        integer, allocatable, dimension(:,:,:)  :: Iz, Jz, Kz, XYZ3
        !
        integer :: i, j, k, Nx, Ny, Nz
        !
        Ix = cvector%x
        Jx = cvector%x
        Kx = cvector%x
        !
        Iy = cvector%y
        Jy = cvector%y
        Ky = cvector%y
        !
        Iz = cvector%z
        Jz = cvector%z
        Kz = cvector%z
        !
        XYZ1 = cvector%x
        XYZ2 = cvector%y
        XYZ3 = cvector%z
        !
        ! X component of the cvector%x
        do i = 1, size( cvector%x, 1 )
            Ix(i,:,:) = i
        enddo
        !
        do j = 1, size( cvector%x, 2 )
            Jx(:,j,:) = j
        enddo
        !
        do  k= 1, size( cvector%x, 3 )
            kx(:,:,k) = k
        enddo
        !
        XYZ1 = 1
        !
        ! Y component of the cvector%y
        do i = 1, size( cvector%y, 1 )
            Iy(i,:,:) = i
        enddo
        !
        do j = 1, size( cvector%y, 2 )
            Jy(:,j,:) = j
        enddo
        !
        do k = 1, size( cvector%y, 3 )
            ky(:,:,k)=k
        enddo
        !
        XYZ2 = 2
        !
        ! Z component of the cvector%z
        do i = 1, size( cvector%z, 1 )
            Iz(i,:,:) = i
        enddo
        !
        do j = 1, size( cvector%z, 2 )
            Jz(:,j,:) = j
        enddo
        !
        do k = 1, size( cvector%z, 3 )
            kz(:,:,k) = k
        enddo
        !
        XYZ3 = 3
        !
        ! Get indexes of Non-Zero coefficients
        if( allocated( self%i ) ) deallocate( self%i )
        allocate( self%i, source = (/ pack(Ix,cvector%x /= 0), pack(Iy,cvector%y /= 0), pack(Iz,cvector%z /= 0) /) )
        !
        if( allocated( self%j ) ) deallocate( self%j )
        allocate( self%j, source = (/ pack(Jx,cvector%x /= 0), pack(Jy,cvector%y /= 0), pack(Jz,cvector%z /= 0) /) )
        !
        if( allocated( self%k ) ) deallocate( self%k )
        allocate( self%k, source = (/ pack(Kx,cvector%x /= 0), pack(Ky,cvector%y /= 0), pack(Kz,cvector%z /= 0) /) )
        !
        ! Get Values of Non-Zero coefficients
        if( allocated( self%c ) ) deallocate( self%c )
        allocate( self%c, source = (/ pack(cvector%x,cvector%x /= 0), pack(cvector%y,cvector%y /= 0), pack(cvector%z,cvector%z /= 0) /) )
        ! Get Components
        if( allocated( self%xyz ) ) deallocate( self%xyz )
        allocate( self%xyz, source = (/ pack(XYZ1,cvector%x /= 0), pack(XYZ2,cvector%y /= 0), pack(XYZ3,cvector%z /= 0) /) )
        !
        ! Get number of Non-Zero coefficients
        self%nCoeff = size( self%c )
        !
        ! Get gridType
        self%gridType = cvector%gridType
        !
        self%is_allocated = .TRUE.
        !
        deallocate( Ix, Jx, Kx, XYZ1 )
        deallocate( Iy, Jy, Ky, XYZ2 )
        deallocate( Iz, Jz, Kz, XYZ3 )
        !
    end subroutine full2Sparse
    !
end module cSparseVector3D_SG  
