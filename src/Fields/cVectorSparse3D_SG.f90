!
!> Lone class to define a CVectorSparse3D_SG
!
module cVectorSparse3D_SG
    !
    use Field
    use cVector3D_SG
    !
    type, extends( Field_t ) :: cVectorSparse3D_SG_t
        !
        integer :: nCoeff
        !
        integer, allocatable, dimension(:) :: i, j, k, xyz
        !
        complex( kind=prec ), allocatable, dimension(:) :: c
        !
        contains
            !
            final :: CVectorSparse3D_SG_dtor
            !
            !> Boundary operations
            procedure, public :: setAllBoundary => setAllBoundaryCVectorSparse3D_SG
            procedure, public :: setOneBoundary => setOneBoundaryCVectorSparse3D_SG
            procedure, public :: setAllInterior => setAllInteriorCVectorSparse3D_SG
            procedure, public :: intBdryIndices => intBdryIndicesCVectorSparse3D_SG
            !
            !> Dimensioning operations
            procedure, public :: length => lengthCVectorSparse3D_SG
            procedure, public :: setVecComponents => setVecComponentsCVectorSparse3D_SG
            !
            !> Arithmetic/algebraic unary operations
            procedure, public :: zeros => zerosCVectorSparse3D_SG
            procedure, public :: sumEdges => sumEdgesCVectorSparse3D_SG
            procedure, public :: avgCells => avgCellsCVectorSparse3D_SG
            procedure, public :: conjugate => conjugateCVectorSparse3D_SG
            !
            !> Arithmetic/algebraic binary operations
            procedure, public :: add => addCVectorSparse3D_SG
            !
            procedure, public :: linComb => linCombCVectorSparse3D_SG
            !
            procedure, public :: subValue => subValueCVectorSparse3D_SG
            procedure, public :: subField => subFieldCVectorSparse3D_SG
            !
            procedure, public :: multByField => multByFieldCVectorSparse3D_SG
            procedure, public :: multByComplex => multByComplexVectorCVectorSparse3D_SG
            procedure, public :: multByReal => multByRealVectorCVectorSparse3D_SG
            !
            procedure, public :: diagMult => diagMultCVectorSparse3D_SG
            !
            procedure, public :: multAdd => multAddCVectorSparse3D_SG
            !
            procedure, public :: dotProd => dotProdCVectorSparse3D_SG
            !
            procedure, public :: divByField => divByFieldCVectorSparse3D_SG
            procedure, public :: divByValue => divByValueCVectorSparse3D_SG
            !
            procedure, public :: interpFunc => interpFuncCVectorSparse3D_SG
            !
            !> Miscellaneous
            procedure, public :: getReal => getRealCVectorSparse3D_SG
            procedure, public :: getArray => getArrayCVectorSparse3D_SG
            procedure, public :: setArray => setArrayCVectorSparse3D_SG
            procedure, public :: switchStoreState => switchStoreStateCVectorSparse3D_SG
            procedure, public :: copyFrom => copyFromCVectorSparse3D_SG
            !
            !> I/O operations
            procedure, public :: read => readCVectorSparse3D_SG
            procedure, public :: write => writeCVectorSparse3D_SG
            procedure, public :: print => printCVectorSparse3D_SG
			!
            procedure, public :: setInteriorMask => setInteriorMaskCVectorSparse3D_SG
            !
            !> Module routines
            procedure, public :: fromFullVector => fromFullVectorCVectorSparse3D_SG
            procedure, public :: getFullVector => getFullVectorCVectorSparse3D_SG
            !
            procedure, public :: reall => reallocateCVectorSparse3D_SG
            !
            procedure, private :: deallocateCVectorSparse3D_SG
            !
    end type cVectorSparse3D_SG_t
    !
    interface cVectorSparse3D_SG_t
        module procedure CVectorSparse3D_SG_ctor
    end interface cVectorSparse3D_SG_t
    !
contains
    !
    !> No subroutine briefing
    !
    function CVectorSparse3D_SG_ctor( nCoeff, grid_type ) result( self )
        implicit none
        !
        integer, intent( in ) :: nCoeff
        character( len=4 ), intent( in ) :: grid_type
        !
        type( cVectorSparse3D_SG_t ) :: self
        !
        integer                    :: status
        !
        !write( *, * ) "Constructor CVectorSparse3D_SG"
        !
        call self%init
        !
        ! the old baggage is out of the door
        if(self%is_allocated) then
            deallocate(self%i, STAT = status)
            deallocate(self%j, STAT = status)
            deallocate(self%k, STAT = status)
            deallocate(self%xyz, STAT = status)
            deallocate(self%c, STAT = status)
            self%grid_type = ''
            self%is_allocated = .FALSE.
        endif
        !
        self%is_allocated = .TRUE.
        allocate(self%i(nCoeff),STAT=status)
        self%is_allocated = self%is_allocated .AND. (status .eq. 0)
        allocate(self%j(nCoeff),STAT=status)
        self%is_allocated = self%is_allocated .AND. (status .eq. 0)
        allocate(self%k(nCoeff),STAT=status)
        self%is_allocated = self%is_allocated .AND. (status .eq. 0)
        allocate(self%xyz(nCoeff),STAT=status)
        self%is_allocated = self%is_allocated .AND. (status .eq. 0)
        allocate(self%c(nCoeff),STAT=status)
        self%is_allocated = self%is_allocated .AND. (status .eq. 0)
        !
        self%nCoeff = nCoeff
        self%i = 0
        self%j = 0
        self%k = 0
        self%xyz = 0
        self%c = C_ZERO
        self%grid_type = grid_type
        !
    end function CVectorSparse3D_SG_ctor
    !
    !> No subroutine briefing
    !
    subroutine CVectorSparse3D_SG_dtor( self )
        implicit none
        !
        type( cVectorSparse3D_SG_t ), intent( inout ) :: self
        !
        !write( *, * ) "Destructor cVectorSparse3D_SG_t:"
        !
        self%grid_type = ""
        self%nCoeff = 0
        self%is_allocated = .FALSE.
        !
        call self%deallocateCVectorSparse3D_SG()
        !
    end subroutine CVectorSparse3D_SG_dtor
    !
    !> No subroutine briefing
    !
    subroutine deallocateCVectorSparse3D_SG( self )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        !
        if( allocated( self%i ) ) deallocate( self%i )
        if( allocated( self%j ) ) deallocate( self%j )
        if( allocated( self%k ) ) deallocate( self%k )
        !
        if( allocated( self%xyz ) ) deallocate( self%xyz )
        !
        if( allocated( self%c ) ) deallocate( self%c )
        !
    end subroutine deallocateCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setAllBoundaryCVectorSparse3D_SG( self, cvalue )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        stop "Error: setAllBoundaryCVectorSparse3D_SG not implemented yet!"
        !
    end subroutine setAllBoundaryCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setAllInteriorCVectorSparse3D_SG( self, cvalue )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        stop "Error: setAllInteriorCVectorSparse3D_SG not implemented yet!"
        !
    end subroutine setAllInteriorCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setOneBoundaryCVectorSparse3D_SG( self, bdry, cvalue, int_only )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        character(*), intent( in ) :: bdry
        complex( kind=prec ), intent( in ) :: cvalue
        logical, intent( in ), optional :: int_only
        logical :: int_only_p
        !
        stop "Error: setOneBoundaryCVectorSparse3D_SG not implemented yet!"
        !
    end subroutine setOneBoundaryCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine intBdryIndicesCVectorSparse3D_SG( self, ind_i, ind_b )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        integer, allocatable, intent( out ) :: ind_i(:), ind_b(:)
        !
        stop "Error: intBdryIndicesCVectorSparse3D_SG not implemented yet!"
        !
    end subroutine intBdryIndicesCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    function dotProdCVectorSparse3D_SG( self, rhs ) result( cvalue )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( in ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        complex( kind=prec ) :: cvalue
        !
        integer :: i, xi, yi, zi
        !
        cvalue = C_ZERO
        !
        if( .NOT. self%is_allocated ) then
            stop "Self not is_allocated yet for dotProdSparse"
        endif
        !
        if( .NOT. rhs%is_allocated ) then
            stop "rhs not is_allocated yet for dotProdSparse"
        endif
        !
        if( self%grid_type /= rhs%grid_type ) then
            stop "dotProdSparse: not compatible usage for dotProdSparse"
        endif
        !
        !> sum over  non-zero terms in sparse vector (conjugate sparse)
        !> (need to check xyz the component)
        !> Remember, xyz = 1,2,3 refers to x, y or z components
        !
        select type( rhs )
            !
            class is( cVector3D_SG_t )
                !
                do i = 1, self%nCoeff
                    !
                    !> generic test for both edge and face (all the components)
                    if( ( self%i(i) .LE. rhs%grid%nx + 1 ) .OR. &
                        ( self%j(i) .LE. rhs%grid%ny + 1 ) .OR. &
                        ( self%k(i) .LE. rhs%grid%nz + 1 ) ) then
                        !
                        !> dealing with x-components
                        if( self%xyz(i) == 1 ) then
                            xi = self%i(i)
                            yi = self%j(i)
                            zi = self%k(i)
                            cvalue = cvalue + conjg( self%c(i) ) * rhs%x( xi, yi, zi )
                        !
                        !> dealing with y-component
                        elseif( self%xyz(i) == 2 ) then
                            xi = self%i(i)
                            yi = self%j(i)
                            zi = self%k(i)
                            cvalue = cvalue + conjg( self%c(i) ) * rhs%y( xi, yi, zi )
                        !
                        !> dealing with z-component
                        elseif( self%xyz(i) == 3 ) then
                            xi = self%i(i)
                            yi = self%j(i)
                            zi = self%k(i)
                            cvalue = cvalue + conjg( self%c(i) ) * rhs%z( xi, yi, zi )
                        endif
                    !
                    else
                        stop "IJK out of bounds for dotProdSparse"
                    !
                    endif
                    !
                enddo
                !
            class default
                stop "Error: dotProdCVector3D_SG: undefined rhs"
            !
        end select
        !
    end function dotProdCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multByComplexVectorCVectorSparse3D_SG( self, cvalue )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        !
        complex( kind=prec ), intent( in ) :: cvalue
        !
        integer :: ii
        !
        do ii = 1, size( self%xyz )
            self%c(ii) = self%c(ii) * cvalue
        enddo
        !
        self%is_allocated = .TRUE.
        !
    end subroutine multByComplexVectorCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multByRealVectorCVectorSparse3D_SG( self, rvalue )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        !
        real( kind=prec ), intent( in ) :: rvalue
        !
        integer :: ii
        !
        do ii = 1, size( self%xyz )
            self%c(ii) = self%c(ii) * rvalue
        enddo
        !
        self%is_allocated = .TRUE.
        !
    end subroutine multByRealVectorCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    function getFullVectorCVectorSparse3D_SG( self ) result( cvector )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( in ) :: self
        !
        type( cVector3D_SG_t ) :: cvector
        !
        integer :: ii
        !
        cvector = cVector3D_SG_t( self%grid, self%grid_type )
        !
        call cvector%zeros
        !
        do ii = 1, size( self%xyz )
            if( self%xyz(ii) == 1 ) then
                cvector%x( self%i(ii), self%j(ii), self%k(ii) ) = self%c(ii)
            elseif( self%xyz(ii) == 2 ) then
                cvector%y( self%i(ii), self%j(ii), self%k(ii) ) = self%c(ii)
            elseif( self%xyz(ii) == 3 ) then
                cvector%z( self%i(ii), self%j(ii), self%k(ii) ) = self%c(ii)
            endif
        enddo
        !
    end function getFullVectorCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine fromFullVectorCVectorSparse3D_SG( self, cvector )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        class( Vector_t ), intent( in ) :: cvector
        !
        integer, allocatable, dimension(:,:,:) :: Ix, Jx, Kx, XYZ1
        integer, allocatable, dimension(:,:,:) :: Iy, Jy, Ky, XYZ2
        integer, allocatable, dimension(:,:,:) :: Iz, Jz, Kz, XYZ3
        !
        integer :: i, j, k, Nx, Ny, Nz
        !
        select type( cvector )
            !
            class is( cVector3D_SG_t )
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
                !> X component of the cvector%x
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
                !> Y component of the cvector%y
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
                !> Z component of the cvector%z
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
                !> Get indexes of Non-Zero coefficients
                self%i = (/ pack(Ix,cvector%x /= 0), pack(Iy,cvector%y /= 0), pack(Iz,cvector%z /= 0) /)
                !
                self%j = (/ pack(Jx,cvector%x /= 0), pack(Jy,cvector%y /= 0), pack(Jz,cvector%z /= 0) /)
                !
                self%k = (/ pack(Kx,cvector%x /= 0), pack(Ky,cvector%y /= 0), pack(Kz,cvector%z /= 0) /)
                !
                !> Get Values of Non-Zero coefficients
                self%c = (/ pack(cvector%x,cvector%x /= 0), pack(cvector%y,cvector%y /= 0), pack(cvector%z,cvector%z /= 0) /)
                !> Get Components
                self%xyz = (/ pack(XYZ1,cvector%x /= 0), pack(XYZ2,cvector%y /= 0), pack(XYZ3,cvector%z /= 0) /)
                !
                !> Set number of Non-Zero coefficients
                self%nCoeff = size( self%c )
                !
                !> Set grid
                self%grid => cvector%grid
                !
                !> Set grid_type
                self%grid_type = trim( cvector%grid_type )
                !
                self%is_allocated = .TRUE.
                !
                deallocate( Ix, Jx, Kx, XYZ1 )
                deallocate( Iy, Jy, Ky, XYZ2 )
                deallocate( Iz, Jz, Kz, XYZ3 )
                !
        end select
        !
    end subroutine fromFullVectorCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    function lengthCVectorSparse3D_SG( self ) result( n )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( in ) :: self
        !
        integer :: n
        !
        stop "Error: lengthCVectorSparse3D_SG not implemented yet!"
        !
    end function lengthCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    function getArrayCVectorSparse3D_SG( self ) result( array )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( in ) :: self
        !
        complex( kind=prec ), allocatable, dimension(:) :: array
        !
        stop "Error: getArrayCVectorSparse3D_SG not implemented yet!"
        !
    end function getArrayCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setArrayCVectorSparse3D_SG( self, array )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), dimension(:), intent( in ) :: array
        !
        stop "Error: setArrayCVectorSparse3D_SG not implemented yet!"
        !
    end subroutine setArrayCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setVecComponentsCVectorSparse3D_SG( self, xyz, &
            &                              xmin, xstep, xmax, &
            &                              ymin, ystep, ymax, &
            &                              zmin, zstep, zmax, rvalue )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        character, intent( in ) :: xyz
        integer, intent( in ) :: xmin, xstep, xmax
        integer, intent( in ) :: ymin, ystep, ymax
        integer, intent( in ) :: zmin, zstep, zmax
        real( kind=prec ), intent ( in ) :: rvalue
        !
        stop "Error: setVecComponentsCVectorSparse3D_SG not implemented yet!"
        !
    end subroutine setVecComponentsCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine zerosCVectorSparse3D_SG( self )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        !
        self%c = R_ZERO
        !
    end subroutine zerosCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine addCVectorSparse3D_SG( self, rhs )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        stop "Error: addCVectorSparse3D_SG not implemented yet!"
        !
    end subroutine addCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine subValueCVectorSparse3D_SG( self, cvalue )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        stop "Error: subValueCVectorSparse3D_SG not implemented yet!"
        !
    end subroutine subValueCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine subFieldCVectorSparse3D_SG( self, rhs )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        stop "Error: subFieldCVectorSparse3D_SG not implemented yet!"
        !
    end subroutine subFieldCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multByFieldCVectorSparse3D_SG( self, rhs )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        stop "Error: multByFieldCVectorSparse3D_SG not implemented yet!"
        !
    end subroutine multByFieldCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine divByFieldCVectorSparse3D_SG( self, rhs )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        stop "Error: divByFieldCVectorSparse3D_SG not implemented yet!"
        !
    end subroutine divByFieldCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine divByValueCVectorSparse3D_SG( self, cvalue )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        !
        stop "Error: divByValueCVectorSparse3D_SG not implemented yet!"
        !
    end subroutine divByValueCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    function diagMultCVectorSparse3D_SG( self, rhs ) result( diag_mult )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        class( Vector_t ), intent( in ) :: rhs
        !
        class( Vector_t ), allocatable :: diag_mult
        !
        stop "Error: diagMultCVectorSparse3D_SG not implemented yet!"
        !
    end function diagMultCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine sumEdgesCVectorSparse3D_SG( self, cell_obj, interior_only )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        class( Field_t ), allocatable, intent( inout ) :: cell_obj
        logical, optional, intent( in ) :: interior_only
        !
        stop "Error: sumEdgesCVectorSparse3D_SG not implemented yet!"
        !
    end subroutine sumEdgesCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine avgCellsCVectorSparse3D_SG( self, E_in, ptype )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: E_in
        character(*), intent( in ), optional :: ptype
        !
        stop "Error: avgCellsCVectorSparse3D_SG not implemented yet!"
        !
    end subroutine avgCellsCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine conjugateCVectorSparse3D_SG( self )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        !
        stop "Error: conjugateCVectorSparse3D_SG not implemented yet!"
        !
    end subroutine conjugateCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine linCombCVectorSparse3D_SG( self, rhs, c1, c2 )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        complex( kind=prec ), intent( in ) :: c1, c2
        !
        stop "Error: linCombCVectorSparse3D_SG not implemented yet!"
        !
    end subroutine linCombCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine multAddCVectorSparse3D_SG( self, cvalue, rhs )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        complex( kind=prec ), intent( in ) :: cvalue
        class( Field_t ), intent( in ) :: rhs
        !
        stop "Error: multAddCVectorSparse3D_SG not implemented yet!"
        !
    end subroutine multAddCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine interpFuncCVectorSparse3D_SG( self, location, xyz, interp )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( in ) :: self
        real( kind=prec ), intent( in ) :: location(3)
        character, intent( in ) :: xyz
        class( Vector_t ), allocatable, intent( inout ) :: interp
        !
        stop "Error: interpFuncCVectorSparse3D_SG not implemented yet!"
        !
    end subroutine interpFuncCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine getRealCVectorSparse3D_SG( self, r_field )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( in ) :: self
        class( Field_t ), allocatable, intent( out ) :: r_field
        !
        stop "Error: getRealCVectorSparse3D_SG not implemented yet!"
        !
    end subroutine getRealCVectorSparse3D_SG
    !
    subroutine switchStoreStateCVectorSparse3D_SG( self )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        !
        stop "Error: switchStoreStateCVectorSparse3D_SG not implemented yet!"
        !
    end subroutine switchStoreStateCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine copyFromCVectorSparse3D_SG( self, rhs )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        class( Field_t ), intent( in ) :: rhs
        !
        if( .NOT. rhs%is_allocated ) then
            stop "Error: copyFromCVectorSparse3D_SG > rhs not allocated"
        endif
        !
        self%grid => rhs%grid
        self%grid_type = rhs%grid_type
        self%nx = rhs%nx
        self%ny = rhs%ny
        self%nz = rhs%nz
        self%is_allocated = .TRUE.
        !
        select type( rhs )
            !
            class is( cVectorSparse3D_SG_t )
                !
                self%nCoeff = rhs%nCoeff
                !
                self%i = rhs%i
                self%j = rhs%j
                self%k = rhs%k
                self%xyz = rhs%xyz
                !
                self%c = rhs%c
                !
            class default
                stop "Error: copyFromCVectorSparse3D_SG > Incompatible rhs"
            !
        end select
        !
    end subroutine copyFromCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine readCVectorSparse3D_SG( self, funit, ftype )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        integer, intent( in ) :: funit
        character(:), allocatable, intent( in ), optional :: ftype
        !
        stop "Error: readCVectorSparse3D_SG not implemented yet!"
        !
    end subroutine readCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine writeCVectorSparse3D_SG( self, funit, ftype )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        integer, intent( in ) :: funit
        character(:), allocatable, intent( in ), optional :: ftype
        !
        stop "Error: setAllBoundaryCVectorSparse3D_SG not implemented yet!"
        !
    end subroutine writeCVectorSparse3D_SG
    !
    !> Reallocates an object of type sparsevecc. The object has to already be
    !> allocated. If allocated and shorter than nCoeff, more memory is
    !> allocated at the end and the contents are preserved.
    !> If allocated and longer than nCoeff, truncates to the first nCoeff values.
    !> This is useful when we need to store the information somewhere, but do
    !> not yet know the final length of the vector. Once it is fully read and
    !> the number of coefficients is known, use this routine to truncate
    !> to the correct length preserving all the values already stored.
    !
    subroutine reallocateCVectorSparse3D_SG( self, nCoeff )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        integer, intent( in ) :: nCoeff
        !
        type( cVectorSparse3D_SG_t ) :: tempLC
        integer :: n, status
        !
        ! the old baggage is out of the door
        if( .NOT. self%is_allocated ) then
            stop "Error: reallocateCVectorSparse3D_SG > Not self%is_allocated"
        endif
        !
        tempLC = self
        !
        if( tempLC%nCoeff .EQ. nCoeff ) then
            ! do nothing
        else
            call self%deallocateCVectorSparse3D_SG()
            self%is_allocated = .true.
            allocate(self%i(nCoeff),STAT=status)
            self%is_allocated = self%is_allocated .AND. (status .eq. 0)
            allocate(self%j(nCoeff),STAT=status)
            self%is_allocated = self%is_allocated .AND. (status .eq. 0)
            allocate(self%k(nCoeff),STAT=status)
            self%is_allocated = self%is_allocated .AND. (status .eq. 0)
            allocate(self%xyz(nCoeff),STAT=status)
            self%is_allocated = self%is_allocated .AND. (status .eq. 0)
            allocate(self%c(nCoeff),STAT=status)
            self%is_allocated = self%is_allocated .AND. (status .eq. 0)
            self%grid_type = tempLC%grid_type
            self%nCoeff = nCoeff
        endif
        !
        if( tempLC%nCoeff > nCoeff ) then
            ! new vector will be shorter
            do n = 1, nCoeff
                self%i(n) = tempLC%i(n)
                self%j(n) = tempLC%j(n)
                self%k(n) = tempLC%k(n)
                self%xyz(n) = tempLC%xyz(n)
                self%c(n) = tempLC%c(n)
            enddo
            !
        elseif( tempLC%nCoeff < nCoeff ) then
            ! new vector will be longer; copy the old values
            do n = 1, tempLC%nCoeff
                self%i(n) = tempLC%i(n)
                self%j(n) = tempLC%j(n)
                self%k(n) = tempLC%k(n)
                self%xyz(n) = tempLC%xyz(n)
                self%c(n) = tempLC%c(n)
            enddo
            ! ... then pad with zeroes
            do n = tempLC%nCoeff + 1, nCoeff
                self%i(n) = 0
                self%j(n) = 0
                self%k(n) = 0
                self%xyz(n) = 0
                self%c(n) = C_ZERO
            enddo
            !
        endif
        !
    end subroutine reallocateCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine printCVectorSparse3D_SG( self, io_unit, title, append )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        integer, intent( in ), optional :: io_unit
        character(*), intent( in ), optional :: title
        logical, intent( in ), optional :: append
        !
        integer :: ii, funit
        !
        if( present( io_unit ) ) then
            funit = io_unit
        else
            funit = 0
        endif
        !
        if( present( title ) ) write( funit, * ) title
        !
        do ii = 1, size( self%xyz )
            if( self%xyz(ii) == 1 ) then
                write( funit, * ) "x(i,j,k):[", self%i(ii), self%j(ii), self%k(ii), "]=", self%c(ii)
            elseif( self%xyz(ii) == 2 ) then
                write( funit, * ) "y(i,j,k):[", self%i(ii), self%j(ii), self%k(ii), "]=", self%c(ii)
            elseif( self%xyz(ii) == 3 ) then
                write( funit, * ) "z(i,j,k):[", self%i(ii), self%j(ii), self%k(ii), "]=", self%c(ii)
            endif
        enddo
       !
    end subroutine printCVectorSparse3D_SG
    !
    !> No subroutine briefing
    !
    subroutine setInteriorMaskCVectorSparse3D_SG( self )
        implicit none
        !
        class( cVectorSparse3D_SG_t ), intent( inout ) :: self
        !
	end subroutine setInteriorMaskCVectorSparse3D_SG
	!
end module cVectorSparse3D_SG  
