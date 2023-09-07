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
module SpOpTopology_SG
    !
    use SpOpTopology
    use Grid3D_SG
    !
    type, extends( SpOpTopology_t ) :: SpOpTopology_SG_t
        !
        type( Grid3D_SG_t ), pointer :: grid
        !
        contains
            !
            procedure, public :: curl => curl_SpOpTopology_SG
            !
            procedure, public :: grad => grad_SpOpTopology_SG
            !
    end type SpOpTopology_SG_t
    !
    interface SpOpTopology_SG_t
        module procedure SpOpTopology_SG_ctor
    end interface SpOpTopology_SG_t
    !
contains
    !
    function SpOpTopology_SG_ctor( grid ) result( self )
        implicit none
        !
        type( Grid3D_SG_t ), pointer, intent( in ) :: grid
        !
        type( SpOpTopology_SG_t ) :: self
        !
        self%grid => grid
        !
    end function SpOpTopology_SG_ctor
    !
    subroutine curl_SpOpTopology_SG( self, curl ) 
        implicit none
        !
        class( SpOpTopology_SG_t ), intent( in ) :: self
        type( spMatCSR_Real ), intent( inout ) :: curl
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
        call create_spMatCSR( m, n, nz, curl )
        !
        deallocate( curl%row )
        allocate( curl%row( curl%nRow + 1 ) )
        !
        deallocate( curl%col )
        allocate( curl%col( 4 * curl%nRow ) )
        !
        deallocate( curl%val )
        allocate( curl%val( 4 * curl%nRow ) )
        !
        curl%row(1) = 1
        !
        do ii = 1, curl%nRow
            curl%row( ii +  1) = 4 * ii + 1
        enddo
        !
        ! xFaces
        allocate(IndVec(nXface))
        allocate(I(nXface))
        allocate(J(nXface))
        allocate(K(nXface))
        !
        do ii = 1, nXface
            IndVec(ii) = ii
        enddo
        !
        call self%grid%gridIndex( XFACE, IndVec, I, J, K )
        call self%grid%vectorIndex( YEDGE, I, J, K, IndVec )
        !
        do ii = 1, nXface
            curl%col( 4 * ii - 3 ) = IndVec(ii) + nXedge
            curl%val( 4 * ii - 3 ) = 1
        enddo
        !
        call self%grid%vectorIndex( ZEDGE, I, J, K, IndVec )
        !
        do ii = 1, nXface
            curl%col( 4 * ii - 1 ) = IndVec(ii) + nXedge + nYedge
            curl%val( 4 * ii - 1 ) = -1
        enddo
        !
        K = K + 1
        !
        call self%grid%vectorIndex( YEDGE, I, J, K, IndVec )
        !
        do ii = 1, nXface
            curl%col(4*ii - 2) = IndVec(ii) + nXedge
            curl%val(4*ii - 2) = -1
        enddo
        !
        K = K - 1
        J = J + 1
        !
        call self%grid%vectorIndex( ZEDGE, I, J, K, IndVec )
        !
        do ii = 1, nXface
            curl%col(4*ii) = IndVec(ii) + nXedge + nYedge
            curl%val(4*ii) = 1
        enddo
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
        enddo
        !
        call self%grid%gridIndex( YFACE, IndVec, I, J, K )
        call self%grid%vectorIndex( XEDGE, I, J, K, IndVec )
        !
        do ii = 1, nYface
            jj = ii + nXface
            curl%col(4*jj - 3) = IndVec(ii)
            curl%val(4*jj - 3) = -1
        enddo
        !
        call self%grid%vectorIndex( ZEDGE, I, J, K, IndVec )
        !
        do ii = 1, nYface
            jj = ii + nXface
            curl%col(4*jj - 1) = IndVec(ii) + nXedge + nYedge
            curl%val(4*jj - 1) = 1
        enddo
        !
        K = K + 1
        !
        call self%grid%vectorIndex( XEDGE, I, J, K, IndVec )
        !
        do ii = 1, nYface
            jj = ii + nXface
            curl%col(4*jj - 2) = IndVec(ii)
            curl%val(4*jj - 2) = 1
        enddo
        !
        K = K - 1
        I = I + 1
        !
        call self%grid%vectorIndex( ZEDGE, I, J, K, IndVec )
        !
        do ii = 1, nYface
            jj = ii + nXface
            curl%col(4*jj) = IndVec(ii) + nXedge + nYedge
            curl%val(4*jj) = -1
        enddo
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
        enddo
        !
        call self%grid%gridIndex( ZFACE, IndVec, I, J, K )
        call self%grid%vectorIndex( XEDGE, I, J, K, IndVec )
        !
        do ii = 1, nZface
            jj = ii + nXface + nYface
            curl%col(4*jj - 3) = IndVec(ii)
            curl%val(4*jj - 3) = 1
        enddo
        !
        call self%grid%vectorIndex( YEDGE, I, J, K, IndVec )
        !
        do ii = 1, nZface
            jj = ii + nXface + nYface
            curl%col(4*jj - 1) = IndVec(ii) + nXedge
            curl%val(4*jj - 1) = -1
        enddo
        !
        J = J + 1
        !
        call self%grid%vectorIndex( XEDGE, I, J, K, IndVec )
        !
        do ii = 1, nZface
            jj = ii + nXface + nYface
            curl%col(4*jj - 2) = IndVec(ii)
            curl%val(4*jj - 2) = -1
        enddo
        !
        I = I + 1
        J = J - 1
        !
        call self%grid%vectorIndex( YEDGE, I, J, K, IndVec )
        !
        do ii = 1, nZface
            jj = ii + nXface + nYface
            curl%col(4*jj) = IndVec(ii) + nXedge
            curl%val(4*jj) = 1
        enddo
        !
        deallocate( IndVec, I, J, K )
        !
    end subroutine curl_SpOpTopology_SG

    subroutine grad_SpOpTopology_SG( self, grad )
        implicit none
        !
        class( SpOpTopology_SG_t ), intent( in ) :: self
        type( spMatCSR_Real ), intent( inout ) :: grad
        !
        integer :: nXedge, nYedge, nZedge, nNodes
        integer :: ii, jj, n, m, nx, ny, nz
        integer, dimension (:), allocatable :: IndVec, I, J, K
        !
        call self%grid%numberOfEdges( nXedge, nYedge, nZedge )
        !
        call self%grid%setLimits( NODE, nx, ny, nz )
        !
        nNodes = nx * ny * nz
        m = nXedge + nYedge + nZedge
        n = nNodes
        nz = 2 * m
        !
        call create_spMatCSR( m, n, nz, grad )
        !
        grad%row(1) = 1
        !
        do ii = 1,grad%nRow
            grad%row(ii + 1) = 2*ii + 1
            grad%val(2*ii - 1) = -1
            grad%val(2*ii) = 1
        enddo
        !
        ! xedges
        allocate(IndVec(nXedge))
        allocate(I(nXedge))
        allocate(J(nXedge))
        allocate(K(nXedge))
        !
        do ii = 1, nXedge
            IndVec(ii) = ii
        enddo
        !
        call self%grid%gridIndex( XEDGE, IndVec, I, J, K )
        call self%grid%vectorIndex( NODE, I, J, K, IndVec )
        !
        do ii = 1, nXedge
            grad%col(2*ii - 1) = IndVec(ii)
        enddo
        !
        I = I + 1
        !
        call self%grid%vectorIndex( NODE, I, J, K, IndVec )
        !
        do ii = 1, nXedge
            grad%col(2*ii) = IndVec(ii)
        enddo
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
        enddo
        !
        call self%grid%gridIndex( YEDGE, IndVec, I, J, K )
        call self%grid%vectorIndex( NODE, I, J, K, IndVec )
        !
        do ii = 1, nYedge
        jj = ii + nXedge
        grad%col(2*jj - 1) = IndVec(ii)
        enddo
        !
        J = J + 1
        !
        call self%grid%vectorIndex(NODE, I, J, K, IndVec)
        !
        do ii = 1, nYedge
        jj = ii + nXedge
        grad%col(2*jj) = IndVec(ii)
        enddo
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
        enddo
        !
        call self%grid%gridIndex( ZEDGE, IndVec, I, J, K )
        call self%grid%vectorIndex( NODE, I, J, K, IndVec )
        !
        do ii = 1, nZedge
            jj = ii + nXedge + nYedge
            grad%col(2*jj - 1) = IndVec(ii)
        enddo
        !
        K = K + 1
        !
        call self%grid%vectorIndex(NODE, I, J, K, IndVec)
        !
        do ii = 1, nZedge
            jj = ii + nXedge + nYedge
            grad%col(2*jj) = IndVec(ii)
        enddo
        !
        deallocate( IndVec, I, J, K )
        !
    end subroutine grad_SpOpTopology_SG

end module SpOpTopology_SG
