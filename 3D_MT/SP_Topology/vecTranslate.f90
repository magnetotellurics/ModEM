! *****************************************************************************
module vecTranslate

  use utilities
  use sg_vector
  use sg_scalar
  implicit none

! Possible node types:  possibly this (and subroutines 
!    setLimits, gridIndex, vectorIndex) should be added to module gridDef
! such as cvector, cscalar, rvector, rscalar, sparsevecc.
!   These names are already defined in gridDef:
!   character(len=80), parameter            :: CENTER = 'CELL'
!     character(len=80), parameter          :: CORNER = 'NODE'
  character(len=80), parameter          :: XFACE = 'XFACE'
  character(len=80), parameter          :: XEDGE = 'XEDGE'
  character(len=80), parameter          :: YFACE = 'YFACE'
  character(len=80), parameter          :: YEDGE = 'YEDGE'
  character(len=80), parameter          :: ZFACE = 'ZFACE'
  character(len=80), parameter          :: ZEDGE = 'ZEDGE'

  INTERFACE getVector
     module procedure getRvector
     module procedure getCvector
  END INTERFACE
  INTERFACE getScalar
     module procedure getRscalar
     module procedure getCscalar
  END INTERFACE
  INTERFACE setVector
     module procedure setRvector
     module procedure setCvector
  END INTERFACE
  INTERFACE setScalar
     module procedure setRscalar
     module procedure setCscalar
  END INTERFACE

Contains

  subroutine getCvector(E,v)
  !    convert an input Cvector E to a 1-D complex array,
  !      following standard staggered grid ordering

    implicit none
    type (cvector), intent(in)            :: E
    complex (kind=prec), pointer, dimension(:), intent(inout)    :: v
    integer nVec(3),nVecT,id(1),i1,i2

    if(.not.E%allocated) then
       call errStop('Input vector not created on entry to getCvector')
    endif
    nVec(1) = size(E%x)
    nVec(2) = size(E%y)
    nVec(3) = size(E%z)
    nVecT = nVec(1)+nVec(2)+nVec(3)
    if(associated(v).and.nVecT.ne.size(v)) then
       deallocate(v)
    endif
    if(.not.associated(v)) then
       allocate(v(nVecT))
    endif
   !   now that we know v is allocated, an of proper size
   !     just copy contents of E into v
    id(1) = nVec(1)
    i1 = 1
    i2 = nVec(1)
    v(i1:i2) = reshape(E%x,id)
    id(1) = nVec(2)
    i1 = i2+1
    i2 = i2+nVec(2)
    v(i1:i2) = reshape(E%y,id)
    id(1) = nVec(3)
    i1 = i2+1
    i2 = i2+nVec(3)
    v(i1:i2) = reshape(E%z,id)

  end subroutine
!*************************************************************
  subroutine setCvector(v,E)
  !    copy contents of v into an already created and allocated
  !        Cvector

    implicit none
    type (cvector), intent(inout)            :: E
    complex (kind=prec), dimension(:), intent(in)    :: v
    integer nVec(3,3),nVecT,id(3),i,i1,i2
    if(.not.E%allocated) then
       call errStop('Output not created before call to setCvector')
    endif
    do i = 1,3 
       nVec(1,i) = size(E%x,i)
       nVec(2,i) = size(E%y,i)
       nVec(3,i) = size(E%z,i)
    enddo
    nVect = 0
    do i =1,3
       nVecT = nVecT + nVec(i,1)*nVec(i,2)*nVec(i,3)
    enddo
    if(nVecT.ne.size(v)) then
       call errStop('Input vector of incorrect size in setCvector')
    endif
   !     copy contents of v into E
    i1 = 1
    i2 = nVec(1,1)*nVec(1,2)*nVec(1,3)
    id = nVec(1,:)
    E%x = reshape(v(i1:i2),id)
    i1 = i2+1
    i2 = i2+nVec(2,1)*nVec(2,2)*nVec(2,3)
    id = nVec(2,:)
    E%y = reshape(v(i1:i2),id)
    i1 = i2+1
    i2 = i2+nVec(3,1)*nVec(3,2)*nVec(3,3)
    id = nVec(3,:)
    E%z = reshape(v(i1:i2),id)
  end subroutine
!**********************************************************************
  subroutine getCscalar(E,v)
  !    convert an input Cscalar E to a 1-D complex array,
  !      following standard staggered grid ordering

    implicit none
    type (cscalar), intent(in)            :: E
    complex (kind=prec), pointer, dimension(:), intent(inout)    :: v
    integer nVecT,id(1)

    if(.not.E%allocated) then
       call errStop('Input vector not created on entry to getCscalar')
    endif
    nVecT = size(E%v)
    if(associated(v).and.nVecT.ne.size(v)) then
       deallocate(v)
    endif
    if(.not.associated(v)) then
       allocate(v(nVecT))
    endif
   !   now that we know v is allocated, an of proper size
   !     just copy contents of E into v
    id(1) = nVecT
    v = reshape(E%v,id)

  end subroutine
!*************************************************************
  subroutine setCscalar(v,E)
  !    copy contents of v into an already created and allocated
  !        Cvector

    implicit none
    type (cscalar), intent(inout)            :: E
    complex (kind=prec), dimension(:), intent(in)    :: v
    integer nVec(3),nVecT,i
    if(.not.E%allocated) then
       call errStop('Output not created before call to setCscalar')
    endif
    do i = 1,3 
       nVec(i) = size(E%v,i)
    enddo
    nVecT = nVec(1)*nVec(2)*nVec(3)

    if(nVecT.ne.size(v)) then
       call errStop('Input vector of incorrect size in setCvector')
    endif
   !     copy contents of v into E
    E%v = reshape(v,nVec)
  end subroutine
  !**********************************************************************
  subroutine getRvector(E,v)
  !    convert an input Rvector E to a 1-D real array,
  !      following standard staggered grid ordering

    implicit none
    type (rvector), intent(in)            :: E
    real (kind=prec), pointer, dimension(:), intent(inout)    :: v
    integer nVec(3),nVecT,id(1),i1,i2

    if(.not.E%allocated) then
       call errStop('Input vector not created on entry to getRvector')
    endif
    nVec(1) = size(E%x)
    nVec(2) = size(E%y)
    nVec(3) = size(E%z)
    nVecT = nVec(1)+nVec(2)+nVec(3)
    if(associated(v)) then
        if(nVect.ne.size(v)) then
           deallocate(v)
        endif
    endif
    if(.not.associated(v)) then
       allocate(v(nVecT))
    endif
    if (nVec(1)>100 .or. nVec(2)>100) then
        write(0,*) 'NOTE: Reshape of large arrays is a common cause of runtime failure at stacksize limit.'
        write(0,*) 'NOTE: Raise the stacksize limit with -heap-arrays compiler flag or ulimit -s unlimited.'
    endif
   !   now that we know v is allocated, an of proper size
   !     just copy contents of E into v
    id(1) = nVec(1)
    i1 = 1
    i2 = nVec(1)
    v(i1:i2) = reshape(E%x,id)
    id(1) = nVec(2)
    i1 = i2+1
    i2 = i2+nVec(2)
    v(i1:i2) = reshape(E%y,id)
    id(1) = nVec(3)
    i1 = i2+1
    i2 = i2+nVec(3)
    v(i1:i2) = reshape(E%z,id)

  end subroutine
!*************************************************************
  subroutine setRvector(v,E)
  !    copy contents of v into an already created and allocated
  !        Rvector

    implicit none
    type (rvector), intent(inout)            :: E
    real (kind=prec), dimension(:), intent(in)    :: v
    integer nVec(3,3),nVecT,id(3),i,i1,i2
    if(.not.E%allocated) then
       call errStop('Output not created before call to setRvector')
    endif
    do i = 1,3 
       nVec(1,i) = size(E%x,i)
       nVec(2,i) = size(E%y,i)
       nVec(3,i) = size(E%z,i)
    enddo
    nVect = 0
    do i =1,3
       nVecT = nVecT + nVec(i,1)*nVec(i,2)*nVec(i,3)
    enddo
    if(nVecT.ne.size(v)) then
       call errStop('Input vector of incorrect size in setRvector')
    endif
   !     copy contents of v into E
    i1 = 1
    i2 = nVec(1,1)*nVec(1,2)*nVec(1,3)
    id = nVec(1,:)
    E%x = reshape(v(i1:i2),id)
    i1 = i2+1
    i2 = i2+nVec(2,1)*nVec(2,2)*nVec(2,3)
    id = nVec(2,:)
    E%y = reshape(v(i1:i2),id)
    i1 = i2+1
    i2 = i2+nVec(3,1)*nVec(3,2)*nVec(3,3)
    id = nVec(3,:)
    E%z = reshape(v(i1:i2),id)
  end subroutine
!**********************************************************************
  subroutine getRscalar(E,v)
  !    convert an input Rscalar E to a 1-D real array,
  !      following standard staggered grid ordering

    implicit none
    type (rscalar), intent(in)            :: E
    real (kind=prec), pointer, dimension(:), intent(inout)    :: v
    integer nVecT,id(1)

    if(.not.E%allocated) then
       call errStop('Input vector not created on entry to getRscalar')
    endif
    nVecT = size(E%v)
    if(associated(v).and.nVecT.ne.size(v)) then
       deallocate(v)
    endif
    if(.not.associated(v)) then
       allocate(v(nVecT))
    endif
   !   now that we know v is allocated, an of proper size
   !     just copy contents of E into v
    id(1) = nVecT
    v = reshape(E%v,id)

  end subroutine
!*************************************************************
  subroutine setRscalar(v,E)
  !    copy contents of v into an already created and allocated
  !        Rscalar
    implicit none
    type (rscalar), intent(inout)            :: E
    real (kind=prec), dimension(:), intent(in)    :: v
    integer nVec(3),nVecT,i
    if(.not.E%allocated) then
       call errStop('Output not created before call to setRscalar')
    endif
    do i = 1,3 
       nVec(i) = size(E%v,i)
    enddo
    nVecT = nVec(1)*nVec(2)*nVec(3)

    if(nVecT.ne.size(v)) then
       call errStop('Input vector of incorrect size in setRvector')
    endif
   !     copy contents of v into E
    E%v = reshape(v,nVec)
  end subroutine
!*************************************************************
!   following subroutines depend only on grid, but are used for
!   converting between ModEM data structures/matrix-free operators
!   and simple column vectors/sparse matrices
  subroutine setLimits(nodeType,grid,nx,ny,nz)
     character(*), intent(in)    :: nodeType
     type(grid_t), intent(in)    :: grid
     integer, intent(out)        :: nx,ny,nz

     selectcase(nodeType)
        case(CENTER)
           nx = grid%nx
           ny = grid%ny
           nz = grid%nz
        case(CORNER)
           nx = grid%nx+1
           ny = grid%ny+1
           nz = grid%nz+1
        case(XEDGE)
           nx = grid%nx
           ny = grid%ny+1
           nz = grid%nz+1
        case(XFACE)
           nx = grid%nx+1
           ny = grid%ny
           nz = grid%nz
        case(YEDGE)
           nx = grid%nx+1
           ny = grid%ny
           nz = grid%nz+1
        case(YFACE)
           nx = grid%nx
           ny = grid%ny+1
           nz = grid%nz
        case(ZEDGE)
           nx = grid%nx+1
           ny = grid%ny+1
           nz = grid%nz
        case(ZFACE)
           nx = grid%nx
           ny = grid%ny
           nz = grid%nz+1
      end select
      return
   end subroutine
!*************************************************************
   subroutine nEdges(grid,nXedge,nYedge,nZedge)
      type(grid_t), intent(in)    :: grid
      integer, intent(out)  ::  nXedge,nYedge,nZedge
      integer nx,ny,nz

      call setLimits(XEDGE,grid,nx,ny,nz)
      nXedge = nx*ny*nz
      call setLimits(YEDGE,grid,nx,ny,nz)
      nYedge = nx*ny*nz
      call setLimits(ZEDGE,grid,nx,ny,nz)
      nZedge = nx*ny*nz
      return
   end subroutine
!*************************************************************
   subroutine nFaces(grid,nXface,nYface,nZface)
      type(grid_t), intent(in)    :: grid
      integer,intent(out) ::  nXface,nYface,nZface
      integer nx,ny,nz

      call setLimits(XFACE,grid,nx,ny,nz)
      nXface = nx*ny*nz
      call setLimits(YFACE,grid,nx,ny,nz)
      nYface = nx*ny*nz
      call setLimits(ZFACE,grid,nx,ny,nz)
      nZface = nx*ny*nz
      return
   end subroutine
!*************************************************************
   subroutine gridIndex(nodeType,grid,IndVec,I,J,K)
   !   based on matlab method of same name in class TGrid3D
   !    IndVec is the index within the list of nodes of a fixed type
   !     e.g., among the list of y-Faces.   An offset needs to be
   !     added to get index in list of all faces (for example)
      character(*), intent(in)    :: nodeType
      type(grid_t), intent(in)    :: grid
      integer, dimension(:), intent(in)        :: IndVec
      integer, dimension(:), intent(out)        :: I,J,K
      integer            ::   nx,ny,nz,nxy,nVec,ii
      real(4)           ::   rNxy,rNx
 
      call setLimits(nodeType,grid,nx,ny,nz)
      nVec = size(IndVec)
      if(nVec.ne.size(I)) then
         call errStop('size of IndVec and I do not agree')
      endif
      if(nVec.ne.size(J)) then
         call errStop('size of IndVec and J do not agree')
      endif
      if(nVec.ne.size(K)) then
         call errStop('size of IndVec and K do not agree')
      endif
      rNxy = float(nx*ny)
      rNx = float(nx)
      do ii = 1,nVec
         I(ii) = mod(IndVec(ii),nx)
         J(ii) = mod(ceiling(float(IndVec(ii))/rNx),ny)
         K(ii) = ceiling(float(IndVec(ii))/rNxy)
      enddo
      where(I.eq.0) I = nx
      where(J.eq.0) J = ny
      where(K.eq.0) K = nz
 
   end subroutine
!*************************************************************
   subroutine vectorIndex(nodeType,grid,I,J,K,IndVec)
   !   based on matlab method of same name in class TGrid3D
   !    returned array IndVec gives numbering of nodes within
   !      the list for nodeType; need to add an offset for position
   !        in full list of all faces or edges (not nodes and cells)
      character(*), intent(in)    :: nodeType
      type(grid_t), intent(in)    :: grid
      integer, dimension(:), intent(out)        :: IndVec
      integer, dimension(:), intent(in)        :: I,J,K
      integer            ::   nx,ny,nz,nxy,nVec,ii
 
      call setLimits(nodeType,grid,nx,ny,nz)
      nVec = size(IndVec)
      if(nVec.ne.size(I)) then
         call errStop('size of IndVec and I do not agree')
      endif
      if(nVec.ne.size(J)) then
         call errStop('size of IndVec and J do not agree')
      endif
      if(nVec.ne.size(K)) then
         call errStop('size of IndVec and K do not agree')
      endif
      nxy = nx*ny
      !   IndVec = (K-1)*nxy+(J-1)*nx+I
      do ii = 1,nVec
         IndVec(ii) = (K(ii)-1)*nxy+(J(ii)-1)*nx + I(ii)
      enddo
   end subroutine
!*************************************************************
   subroutine boundaryIndex(gridType,grid,INDb,INDi)
   !   for a given type find indicies for boundary and 
   !      interior nodes
      character(*), intent(in)    :: gridType
      type(grid_t), intent(in)    :: grid
      integer,intent(inout),allocatable, dimension(:) :: INDb,INDi

      integer   :: nVec(3),nVecT,nBdry,nb,ni,i
      real(kind=prec),pointer   :: temp(:)
      type(rvector)     :: E
      type(rscalar)     :: phi

!     write(0,*) gridType
!     write(0,*) 'grid.nx,ny,nz', grid%nx, grid%ny, grid%nz, grid%nzAir
      selectcase(gridType)
         case(EDGE)
            call create_rvector(grid,E,EDGE)
            nVec(1) = size(E%x)
            nVec(2) = size(E%y)
            nVec(3) = size(E%z)
            nVecT = nVec(1)+nVec(2)+nVec(3)
            !write(0,*) nVec, nVecT
            allocate(temp(nVecT))
            E%x(:,1,:) = 1
            E%x(:,E%ny+1,:) = 1
            E%x(:,:,1) = 1
            E%x(:,:,E%nz+1) = 1
            E%y(1,:,:) = 1
            E%y(E%nx+1,:,:) = 1
            E%y(:,:,1) = 1
            E%y(:,:,E%nz+1) = 1
            E%z(1,:,:) = 1
            E%z(E%nx+1,:,:) = 1
            E%z(:,1,:) = 1
            E%z(:,E%ny+1,:) = 1
            call getVector(E,temp) 
            !write(0,*) 'getVector complete'
            call deall_rvector(E)
        case(FACE)
            call create_rvector(grid,E,FACE)
            nVec(1) = size(E%x)
            nVec(2) = size(E%y)
            nVec(3) = size(E%z)
            nVecT = nVec(1)+nVec(2)+nVec(3)
            allocate(temp(nVecT))
            E%x(1,:,:) = 1
            E%x(E%nx+1,:,:) = 1
            E%y(:,1,:) = 1
            E%y(:,E%ny+1,:) = 1
            E%z(:,:,1) = 1
            E%z(:,:,E%nz+1) = 1
            call getVector(E,temp) 
            call deall_rvector(E)
        case(CORNER)
            call create_rscalar(grid,phi,CORNER)
            nVecT = size(phi%v)
            allocate(temp(nVecT))
            phi%v(1,:,:) = 1
            phi%v(phi%nx+1,:,:) = 1
            phi%v(:,1,:) = 1
            phi%v(:,phi%ny+1,:) = 1
            phi%v(:,:,1) = 1
            phi%v(:,:,phi%nz+1) = 1
            call getRscalar(phi,temp) 
            call deall_rscalar(phi)
        end select 
        nBdry = 0
        do i = 1,nVecT
           nBdry = nBdry+nint(temp(i))
        enddo
        if(allocated(INDi)) then
           deallocate(INDi)
        endif
        allocate(INDi(nVecT-nBdry))
        if(allocated(INDb)) then
           deallocate(INDb)
        endif
        allocate(INDb(nBdry))
        nb = 0
        ni = 0
        do i = 1,nVecT
           if(nint(temp(i)).eq.1) then
              nb = nb+1
              INDb(nb) = i
           else
              ni = ni+1
              INDi(ni) = i
           endif
        enddo
        deallocate(temp)
        return 
    end subroutine
end module
