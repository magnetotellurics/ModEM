module spOpTopology_SG

use  spOpTools
use  vectranslate

     implicit none
     save
     ! private
     public
!    not sure what this module really represents ... if just topology
!    do we want to keep the grid here?
!     type(grid_t), private, target   ::      mGrid   ! The model grid

!    Sparse grad topology : maps from all nodes to all edges
     type(spMatCSR_real)     :: G

!    Sparse curl topology : maps from all edges to all faces
     type(spMatCSR_real)     :: T

!    this could just be topology, or we could put metric elements here
!    also ... then we would also have actual curl-curl + diagonal (A)

Contains
   subroutine setCurlTopology(grid)
   !  this maps from all edges to all faces; rows of matrix correspond to  faces
   !     enumerated with all xFaces, then yFaces, then zFaces
   !   Columns follow a similar pattern for edges
      type (grid_t), intent(in)            :: grid 
      integer    :: nXedge,nYedge,nZedge,nXface,nYface,nZface,ii,jj,n,m,nz
      integer, allocatable, dimension(:)    :: IndVec,I,J,K

     
      call nEdges(grid,nXedge,nYedge,nZedge)
      call nFaces(grid,nXface,nYface,nZface)
      m = nXface+nYface+nZface
      n = nXedge+nYedge+nZedge
      nz = 4*m
      call create_spMatCSR(m,n,nz,T)
!     write(0,*) 'm,n,nz for Curl', m,n,nz
      deallocate(T%row)
      allocate(T%row(T%nRow+1))
      deallocate(T%col)
      allocate(T%col(4*T%nRow))
      deallocate(T%val)
      allocate(T%val(4*T%nRow))
      T%row(1) = 1
      do ii = 1,T%nRow
         T%row(ii+1) = 4*ii+1
         ! probably better to set values as we go ...
      enddo

     !   xFaces
      allocate(IndVec(nXface))
      allocate(I(nXface))
      allocate(J(nXface))
      allocate(K(nXface))
      do ii=1,nXface
         IndVec(ii) = ii
      enddo
      call gridIndex('XFACE',grid,IndVec,I,J,K)
      call vectorIndex('YEDGE',grid,I,J,K,IndVec)
      do ii = 1,nXface
         T%col(4*ii-3) = IndVec(ii)+nXedge
         T%val(4*ii-3) = 1
      enddo
      call vectorIndex('ZEDGE',grid,I,J,K,IndVec)
      do ii = 1,nXface
         T%col(4*ii-1) = IndVec(ii)+nXedge+nYedge
         T%val(4*ii-1) = -1
      enddo
      K = K+1
      call vectorIndex('YEDGE',grid,I,J,K,IndVec)
      do ii = 1,nXface
         T%col(4*ii-2) = IndVec(ii)+nXedge
         T%val(4*ii-2) = -1
      enddo
      K = K-1
      J = J+1
      call vectorIndex('ZEDGE',grid,I,J,K,IndVec)
      do ii = 1,nXface
         T%col(4*ii) = IndVec(ii)+nXedge+nYedge
         T%val(4*ii) = 1
      enddo
      deallocate(IndVec)
      deallocate(I)
      deallocate(J)
      deallocate(K)

     !   yFaces
      allocate(IndVec(nYface))
      allocate(I(nYface))
      allocate(J(nYface))
      allocate(K(nYface))
      do ii=1,nYface
         IndVec(ii) = ii
      enddo
      call gridIndex('YFACE',grid,IndVec,I,J,K)
      call vectorIndex('XEDGE',grid,I,J,K,IndVec)
      do ii = 1,nYface
         jj = ii + nXface
         T%col(4*jj-3) = IndVec(ii)
         T%val(4*jj-3) = -1
      enddo
      call vectorIndex('ZEDGE',grid,I,J,K,IndVec)
      do ii = 1,nYface
         jj = ii + nXface
         T%col(4*jj-1) = IndVec(ii)+nXedge+nYedge
         T%val(4*jj-1) = 1
      enddo
      K = K+1
      call vectorIndex('XEDGE',grid,I,J,K,IndVec)
      do ii = 1,nYface
         jj = ii + nXface
         T%col(4*jj-2) = IndVec(ii)
         T%val(4*jj-2) = 1
      enddo
      K = K-1
      I = I+1
      call vectorIndex('ZEDGE',grid,I,J,K,IndVec)
      do ii = 1,nYface
         jj = ii + nXface
         T%col(4*jj) = IndVec(ii)+nXedge+nYedge
         T%val(4*jj) = -1
      enddo
      deallocate(IndVec)
      deallocate(I)
      deallocate(J)
      deallocate(K)

     !   zFaces
      allocate(IndVec(nZface))
      allocate(I(nZface))
      allocate(J(nZface))
      allocate(K(nZface))
      do ii=1,nZface
         IndVec(ii) = ii
      enddo
      call gridIndex('ZFACE',grid,IndVec,I,J,K)
      call vectorIndex('XEDGE',grid,I,J,K,IndVec)
      do ii = 1,nZface
         jj = ii + nXface + nYface
         T%col(4*jj-3) = IndVec(ii)
         T%val(4*jj-3) = 1
      enddo
      call vectorIndex('YEDGE',grid,I,J,K,IndVec)
      do ii = 1,nZface
         jj = ii + nXface + nYface
         T%col(4*jj-1) = IndVec(ii)+nXedge
         T%val(4*jj-1) = -1
      enddo
      J = J+1
      call vectorIndex('XEDGE',grid,I,J,K,IndVec)
      do ii = 1,nZface
         jj = ii + nXface + nYface
         T%col(4*jj-2) = IndVec(ii)
         T%val(4*jj-2) = -1
      enddo
      I = I+1
      J = J-1
      call vectorIndex('YEDGE',grid,I,J,K,IndVec)
      do ii = 1,nZface
         jj = ii + nXface + nYface
         T%col(4*jj) = IndVec(ii)+nXedge
         T%val(4*jj) = 1
      enddo
      deallocate(IndVec)
      deallocate(I)
      deallocate(J)
      deallocate(K)
      return
   end subroutine
!**********************************************************************
   subroutine setGradTopology(grid)
   !  this maps from all nodes to all edges; rows of matrix correspond to  edges
   !     enumerated with all xEdges, then yEdges, then zEdges
      type (grid_t), intent(in)            :: grid 
      integer    :: nXedge,nYedge,nZedge,nNodes,ii,jj,n,m,nz,nx,ny
      integer, allocatable, dimension(:)    :: IndVec,I,J,K
 
      call nEdges(grid,nXedge,nYedge,nZedge)
      call setLimits(CORNER,grid,nx,ny,nz)
!     write(0,*) 'nx,ny,nz',nz,ny,nz
      nNodes = nx*ny*nz
      m = nXedge+nYedge+nZedge
      n = nNodes
      nz = 2*m
      call create_spMatCSR(m,n,nz,G)
!     write(0,*) 'm,n,nz for G', m,n,nz
!     write(0,*) '# edges',nXedge,nYedge,nZedge
      G%row(1) = 1 
      do ii = 1,G%nRow
         G%row(ii+1) = 2*ii + 1
         G%val(2*ii-1) = -1
         G%val(2*ii) = 1
      enddo
      !   xedges
      allocate(IndVec(nXedge))
      allocate(I(nXedge))
      allocate(J(nXedge))
      allocate(K(nXedge))
      do ii=1,nXedge
         IndVec(ii) = ii
      enddo
      call gridIndex('XEDGE',grid,IndVec,I,J,K)
      call vectorIndex(CORNER,grid,I,J,K,IndVec)
      do ii = 1,nXedge
         G%col(2*ii-1) = IndVec(ii)
      enddo
      I = I+1
      call vectorIndex(CORNER,grid,I,J,K,IndVec)
      do ii = 1,nXedge
         G%col(2*ii) = IndVec(ii)
      enddo
      deallocate(IndVec)
      deallocate(I)
      deallocate(J)
      deallocate(K)
      !   yedges
      allocate(IndVec(nYedge))
      allocate(I(nYedge))
      allocate(J(nYedge))
      allocate(K(nYedge))
      do ii=1,nYedge
         IndVec(ii) = ii
      enddo
      call gridIndex('YEDGE',grid,IndVec,I,J,K)
      call vectorIndex(CORNER,grid,I,J,K,IndVec)
      do ii = 1,nYedge
         jj = ii+nXedge
         G%col(2*jj-1) = IndVec(ii)
      enddo
      J = J+1
      call vectorIndex(CORNER,grid,I,J,K,IndVec)
      do ii = 1,nYedge
         jj = ii+nXedge
         G%col(2*jj) = IndVec(ii)
      enddo
      deallocate(IndVec)
      deallocate(I)
      deallocate(J)
      deallocate(K)
      !   zedges
      allocate(IndVec(nZedge))
      allocate(I(nZedge))
      allocate(J(nZedge))
      allocate(K(nZedge))
      do ii=1,nZedge
         IndVec(ii) = ii
      enddo
      call gridIndex('ZEDGE',grid,IndVec,I,J,K)
      call vectorIndex(CORNER,grid,I,J,K,IndVec)
      do ii = 1,nZedge
         jj = ii+nXedge+nYedge
         G%col(2*jj-1) = IndVec(ii)
      enddo
      K = K+1
      call vectorIndex(CORNER,grid,I,J,K,IndVec)
      do ii = 1,nZedge
         jj = ii+nXedge+nYedge
         G%col(2*jj) = IndVec(ii)
      enddo
      deallocate(IndVec)
      deallocate(I)
      deallocate(J)
      deallocate(K)
   end subroutine
end module
