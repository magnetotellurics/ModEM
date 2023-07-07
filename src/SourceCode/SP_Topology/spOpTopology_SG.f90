!
!> This file is part of the ModEM modeling and inversion package.
!> 
!> LICENSING information
!
!> Copyright (C) 2020 ModEM research group.
!> Contact: http://
!
!> GNU General Public License Usage
!> This file may be used under the terms of the GNU
!> General Public License version 3.0 as published by the Free Software
!> Foundation and appearing in the file LICENSE.GPL included in the
!> packaging of this file.  Please review the following information to
!> ensure the GNU General Public License version 3.0 requirements will be
!> met: http://www.gnu.org/copyleft/gpl.html.
!> 
!> SUMMARY
!> 
!> No Module briefing
!
module spOpTopology_SG
    !
    use spOpTopology
    use Grid
    !
    type, extends( spOpTopology_t ) :: spOpTopology_SG_t
        !
        class( Grid_t ), pointer :: grid
        !
        contains
            !
            procedure, public :: curl => curlSpOpTopology_SG
            !
            procedure, public :: grad => gradSpOpTopology_SG
            !
    end type spOpTopology_SG_t
    !
    interface spOpTopology_SG_t
        module procedure spOpTopology_SG_ctor
    end interface spOpTopology_SG_t
    !
contains
    !
    function spOpTopology_SG_ctor( grid ) result( self )
        implicit none
        !
        class( Grid_t ), pointer, intent( in ) :: grid
        !
        type( spOpTopology_SG_t ) :: self
        !
        self%grid => grid
        !
    end function spOpTopology_SG_ctor
    !
    subroutine curlSpOpTopology_SG( self, T ) 
        implicit none
        !
        class( spOpTopology_SG_t ), intent( inout ) :: self
        type( spMatCSR_Real ), intent( inout ) :: T
        !  
        integer :: nXedge, nYedge, nZedge
        integer :: nXface, nYface, nZface
        integer :: ii, jj, n, m, nz
        integer, dimension (:), allocatable :: IndVec, I, J, K
        !
        call self%grid%numberOfEdges( nXedge, nYedge, nZedge )
        !
        call self%grid%numberOfFaces( nXface, nYface, nZface )
        !
        m = nXface + nYface + nZface
        n = nXedge + nYedge + nZedge
        nz = 4 * m
        !
        call create_spMatCSR( m, n, nz, T )
        !
        deallocate( T%row )
        allocate( T%row( T%nRow + 1 ) )
        !
        deallocate( T%col )
        allocate( T%col( 4 * T%nRow ) )
        !
        deallocate( T%val )
        allocate( T%val( 4 * T%nRow ) )
        !
        T%row(1) = 1
        !
        do ii = 1, T%nRow
            T%row( ii +  1) = 4 * ii + 1
        end do
        !
        ! xFaces
        allocate(IndVec(nXface))
        allocate(I(nXface))
        allocate(J(nXface))
        allocate(K(nXface))
        !
        do ii = 1, nXface
            IndVec(ii) = ii
        end do
        !
        call self%grid%gridIndex( XFACE, IndVec, I, J, K )
        call self%grid%vectorIndex( YEDGE, I, J, K, IndVec )
        !
        do ii = 1, nXface
            T%col( 4 * ii - 3 ) = IndVec(ii) + nXedge
            T%val( 4 * ii - 3 ) = 1
        end do
        !
        call self%grid%vectorIndex( ZEDGE, I, J, K, IndVec )
        !
        do ii = 1, nXface
            T%col( 4 * ii - 1 ) = IndVec(ii) + nXedge + nYedge
            T%val( 4 * ii - 1 ) = -1
        end do
        !
        K = K + 1
        !
        call self%grid%vectorIndex( YEDGE, I, J, K, IndVec )
        !
        do ii = 1, nXface
            T%col(4*ii - 2) = IndVec(ii) + nXedge
            T%val(4*ii - 2) = -1
        end do
        !
        K = K - 1
        J = J + 1
        !
        call self%grid%vectorIndex( ZEDGE, I, J, K, IndVec )
        !
        do ii = 1, nXface
            T%col(4*ii) = IndVec(ii) + nXedge + nYedge
            T%val(4*ii) = 1
        end do
        !
        deallocate( IndVec, I, J, K )
        !
        ! yFaces
        allocate( IndVec( nYface ) )
        allocate( I( nYface ) )
        allocate( J( nYface ) )
        allocate( K( nYface ) )
        !
        do ii = 1, nYface
            IndVec(ii) = ii
        end do
        !
        call self%grid%gridIndex( YFACE, IndVec, I, J, K )
        call self%grid%vectorIndex( XEDGE, I, J, K, IndVec )
        !
        do ii = 1, nYface
            jj = ii + nXface
            T%col(4*jj - 3) = IndVec(ii)
            T%val(4*jj - 3) = -1
        end do
        !
        call self%grid%vectorIndex( ZEDGE, I, J, K, IndVec )
        !
        do ii = 1, nYface
            jj = ii + nXface
            T%col(4*jj - 1) = IndVec(ii) + nXedge + nYedge
            T%val(4*jj - 1) = 1
        end do
        !
        K = K + 1
        !
        call self%grid%vectorIndex( XEDGE, I, J, K, IndVec )
        !
        do ii = 1, nYface
            jj = ii + nXface
            T%col(4*jj - 2) = IndVec(ii)
            T%val(4*jj - 2) = 1
        end do
        !
        K = K - 1
        I = I + 1
        !
        call self%grid%vectorIndex( ZEDGE, I, J, K, IndVec )
        !
        do ii = 1, nYface
            jj = ii + nXface
            T%col(4*jj) = IndVec(ii) + nXedge + nYedge
            T%val(4*jj) = -1
        end do
        !
        deallocate( IndVec, I, J, K )
        !
        ! zFaces
        allocate(IndVec(nZface))
        allocate(I(nZface))
        allocate(J(nZface))
        allocate(K(nZface))
        !
        do ii = 1, nZface
            IndVec(ii) = ii
        end do
        !
        call self%grid%gridIndex( ZFACE, IndVec, I, J, K )
        call self%grid%vectorIndex( XEDGE, I, J, K, IndVec )
        !
        do ii = 1, nZface
            jj = ii + nXface + nYface
            T%col(4*jj - 3) = IndVec(ii)
            T%val(4*jj - 3) = 1
        end do
        !
        call self%grid%vectorIndex( YEDGE, I, J, K, IndVec )
        !
        do ii = 1, nZface
            jj = ii + nXface + nYface
            T%col(4*jj - 1) = IndVec(ii) + nXedge
            T%val(4*jj - 1) = -1
        end do
        !
        J = J + 1
        !
        call self%grid%vectorIndex( XEDGE, I, J, K, IndVec )
        !
        do ii = 1, nZface
            jj = ii + nXface + nYface
            T%col(4*jj - 2) = IndVec(ii)
            T%val(4*jj - 2) = -1
        end do
        !
        I = I + 1
        J = J - 1
        !
        call self%grid%vectorIndex( YEDGE, I, J, K, IndVec )
        !
        do ii = 1, nZface
            jj = ii + nXface + nYface
            T%col(4*jj) = IndVec(ii) + nXedge
            T%val(4*jj) = 1
        end do
        !
        deallocate( IndVec, I, J, K )
        !
    end subroutine curlSpOpTopology_SG

    subroutine gradSpOpTopology_SG( self, G )
        implicit none
        !
        class( spOpTopology_SG_t ), intent( inout ) :: self
        type( spMatCSR_Real ), intent( inout ) :: G
        !
        integer :: nXedge, nYedge, nZedge, nNodes
        integer :: ii, jj, n, m, nx, ny, nz
        integer, dimension (:), allocatable :: IndVec, I, J, K
        !
        call self%grid%numberOfEdges( nXedge, nYedge, nZedge )
        call self%grid%limits( NODE, nx, ny, nz )
        !
        nNodes = nx * ny * nz
        m = nXedge + nYedge + nZedge
        n = nNodes
        nz = 2 * m
        !
        call create_spMatCSR( m, n, nz, G )
        !
        G%row(1) = 1
        !
        do ii = 1,G%nRow
            G%row(ii + 1) = 2*ii + 1
            G%val(2*ii - 1) = -1
            G%val(2*ii) = 1
        end do
        !
        ! xedges
        allocate(IndVec(nXedge))
        allocate(I(nXedge))
        allocate(J(nXedge))
        allocate(K(nXedge))
        !
        do ii = 1, nXedge
            IndVec(ii) = ii
        end do
        !
        call self%grid%gridIndex( XEDGE, IndVec, I, J, K )
        call self%grid%vectorIndex( NODE, I, J, K, IndVec )
        !
        do ii = 1, nXedge
            G%col(2*ii - 1) = IndVec(ii)
        end do
        !
        I = I + 1
        !
        call self%grid%vectorIndex( NODE, I, J, K, IndVec )
        !
        do ii = 1, nXedge
            G%col(2*ii) = IndVec(ii)
        end do
        !
        deallocate( IndVec, I, J, K )
        !
        ! yedges
        allocate(IndVec(nYedge))
        allocate(I(nYedge))
        allocate(J(nYedge))
        allocate(K(nYedge))
        !
        do ii = 1, nYedge
        IndVec(ii) = ii
        end do
        !
        call self%grid%gridIndex( YEDGE, IndVec, I, J, K )
        call self%grid%vectorIndex( NODE, I, J, K, IndVec )
        !
        do ii = 1, nYedge
        jj = ii + nXedge
        G%col(2*jj - 1) = IndVec(ii)
        end do
        !
        J = J + 1
        !
        call self%grid%vectorIndex(NODE, I, J, K, IndVec)
        !
        do ii = 1, nYedge
        jj = ii + nXedge
        G%col(2*jj) = IndVec(ii)
        end do
        !
        deallocate( IndVec, I, J, K )
        !
        ! zedges
        allocate(IndVec(nZedge))
        allocate(I(nZedge))
        allocate(J(nZedge))
        allocate(K(nZedge))
        !
        do ii = 1, nZedge
            IndVec(ii) = ii
        end do
        !
        call self%grid%gridIndex( ZEDGE, IndVec, I, J, K )
        call self%grid%vectorIndex( NODE, I, J, K, IndVec )
        !
        do ii = 1, nZedge
            jj = ii + nXedge + nYedge
            G%col(2*jj - 1) = IndVec(ii)
        end do
        !
        K = K + 1
        !
        call self%grid%vectorIndex(NODE, I, J, K, IndVec)
        !
        do ii = 1, nZedge
            jj = ii + nXedge + nYedge
            G%col(2*jj) = IndVec(ii)
        end do
        !
        deallocate( IndVec, I, J, K )
        !
    end subroutine gradSpOpTopology_SG

end module spOpTopology_SG
