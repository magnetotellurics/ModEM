! *****************************************************************************
! Initializes and does basic calculations for the grid. Computes basic derived
! grid parameters that are used repeatedly by other routines, such as the unit
! lengths and volumes, and other grid details.
! Important note: the spherical staggered-grid elements borrowed from the global
! code assume that the dual grid is shifted by + 1/2 element length in all directions.
! This works ok if all vectors on the dual grid contain no boundary elements;
! however if we want to represent the boundaries in this convention (incl. North
! pole for the spherical code) we are left with zero index for these values.
! Instead, we assume - 1/2 element length shifts for the dual grid throughout the code.
! We call the length and area elements carefully to match these conventions.
module GridCalc

  use sg_vector
  use sg_scalar
  use sg_spherical
  use elements
  implicit none

  save

  !!!!!!!>>>>>>>>> block of precomputed grid elements (public)
  !!!!!!!>>>>>>>>> initialized in ModelDataInit to accommodate
  !!!!!!!>>>>>>>>> for potential on-the-fly grid modification
  type(rvector), public     :: V_E, l_E, S_F ! defined on primary grid
  type(rvector), public     :: V_F, l_F, S_E ! defined on dual grid
  type(rscalar), public     :: V_N, V_C ! nodes and cells

  public      :: create_gridElements, deall_gridElements
  public      :: EdgeVolume, FaceVolume, NodeVolume, CellVolume
  public      :: EdgeLength, DualEdgeLength
  public      :: FaceArea,   DualFaceArea
  public      :: Cell2Edge,  Edge2Cell
  public      :: Cell2Node

  !!!!!!!>>>>>>>>> global parameter key to enable switching
  !!!!!!!>>>>>>>>> between cartesian and spherical grids
  !!!!!!!>>>>>>>>> (spherical can be global or regional)
  character (len=80), public, parameter  :: gridCoords = SPHERICAL

Contains


  ! *************************************************************************
  ! * EdgeVolume creates volume elements centered around the edges of
  ! * the grid, and stores them as real vectors with gridType=EDGE.
  ! *
  ! * A diagonal matrix multiplication of the edge volume with the difference
  ! * equations enables us to make a symmetrical matrix. Remember,
  ! * the electrical fields are defined on the center of the edges, therefore,
  ! * the edge volume is centered about the electrical field measurement.

  subroutine EdgeVolume(grid, V_E, l_E, S_E)

      implicit none
      type (grid_t), intent(in)           :: grid     ! input model
      type (rvector), intent(inout)       :: V_E       ! edge volume
      type (rvector), intent(in), optional:: l_E, S_E ! optional inputs
      ! local
      type (rvector)                      :: length, area
      logical                             :: compute_elements

      write(*,*) 'Computing edge volume elements...'

      if (.not. V_E%allocated) then
        call create_rvector(grid, V_E, EDGE)
      else
        ! Checks that the size and type are the same
        if (((grid%nx == V_E%nx).and.(grid%ny == V_E%ny).and.(grid%nz == V_E%nz)) &
            .and. (V_E%gridType == EDGE)) then
            ! vector adequately allocated
        else
            ! reallocate
            call deall_rvector(V_E)
            call create_rvector(grid, V_E, EDGE)
        end if
      end if

      compute_elements = .true.
      if ((present(l_E)) .and. (present(S_E))) then
        if ((l_E%allocated) .and. (S_E%allocated)) then
            compute_elements = .false.
        endif
      endif

      ! create length and surface elements vectors
      if (compute_elements) then
        call EdgeLength(grid,length)
        call DualFaceArea(grid,area)
      else
        length = l_E
        area = S_E
      end if

      ! compute volume elements
      call diagMult_rvector(length,area,V_E)

      call deall_rvector(length)
      call deall_rvector(area)

  end subroutine EdgeVolume  ! EdgeVolume

  ! *************************************************************************
  ! * FaceVolume creates volume elements centered around the edges of
  ! * the dual grid, and stores them as real vectors with gridType=FACE.
  ! *
  ! * A diagonal matrix multiplication by the face volumes is part of the
  ! * natural adjoint of the curl operator.

  subroutine FaceVolume(grid, V_F, l_F, S_F)

      implicit none
      type (grid_t), intent(in)           :: grid     ! input model
      type (rvector), intent(inout)       :: V_F      ! face volume
      type (rvector), intent(in), optional:: l_F, S_F ! optional inputs
      ! local
      type (rvector)                      :: length, area
      logical                             :: compute_elements

      write(*,*) 'Computing face volume elements...'

      if (.not. V_F%allocated) then
        call create_rvector(grid, V_F, FACE)
      else
        ! Checks that the size and type are the same
        if (((grid%nx == V_F%nx).and.(grid%ny == V_F%ny).and.(grid%nz == V_F%nz)) &
            .and. (V_F%gridType == FACE)) then
            ! vector adequately allocated
        else
            ! reallocate
            call deall_rvector(V_F)
            call create_rvector(grid, V_F, FACE)
        end if
      end if

      compute_elements = .true.
      if ((present(l_F)) .and. (present(S_F))) then
        if ((l_F%allocated) .and. (S_F%allocated)) then
            compute_elements = .false.
        endif
      endif

      ! create length and surface elements vectors
      if (compute_elements) then
        call DualEdgeLength(grid,length)
        call FaceArea(grid,area)
      else
        length = l_F
        area = S_F
      end if

      ! compute volume elements
      call diagMult_rvector(length,area,V_F)

      call deall_rvector(length)
      call deall_rvector(area)

  end subroutine FaceVolume  ! FaceVolume

  ! *************************************************************************
  ! * NodeVolume creates volume elements centered around the corners (nodes)
  ! * of the grid, and stores them as real scalars with gridType=CORNER.

  subroutine NodeVolume(grid, V_N)

      type (grid_t), intent(in)          :: grid  ! input grid
      type (rscalar), intent(inout)      :: V_N   ! node/corner volume as output
      type (rscalar)                     :: gV_N  ! temporary vector in global coord system
      ! local variables
      integer                            :: nx,ny,nz
      real(8),dimension(:),allocatable   :: x,y,z
      real(8)                            :: volume,sijk2,zm,zp,ym,yp
      integer                            :: i, j, k, istat

      write(*,*) 'Computing node volume elements...'

      if (.not. V_N%allocated) then
        call create_rscalar(grid, V_N, CORNER)
      else
        ! Checks that the size and type are the same
        if (((grid%nx == V_N%nx).and.(grid%ny == V_N%ny).and.(grid%nz == V_N%nz)) &
            .and. (V_N%gridType == CORNER)) then
            ! vector adequately allocated
        else
            ! reallocate
            call deall_rscalar(V_N)
            call create_rscalar(grid, V_N, CORNER)
        end if
      end if

      ! flip the x & y directions to use the conventions in the spherical code
      call coordFlip_rscalar(gV_N, V_N, Spherical)

      nx = grid%ny !grid%nx
      ny = grid%nx !grid%ny
      nz = grid%nz

      ! need to be very careful here:
      ! x historically denotes *interval* in longitude;
      ! y denotes co-latitude from North to South pole;
      ! z is +ve down but denotes radii from the Earth's center
      ! moreover, all elements assume that dual grid is shifted
      ! in the positive direction - careful with indices
      allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
      x= grid%dy*pi/180.                !x = grid%x
      y= (90.- grid%xEdge)*pi/180.      !y = grid%y
      z = EARTH_R * KM2M - grid%zEdge           !z = grid%z


      ! node volumes are only using the internal corner nodes
      ! but need to be careful on a global grid
      do k=1,nz+1
        do j=1,ny+1
            do i=1,nx+1
               call volume_vijk2(i-1,i,j-1,k-1,x,y,z,volume)
               gV_N%v(i, j, k) = volume
            enddo
        enddo
      enddo

      deallocate(x,y,z)

      ! flip the coordinate directions back to MT code notations
      call coordFlip_rscalar(V_N, gV_N, Cartesian)

  end subroutine NodeVolume

  ! *************************************************************************
  ! * CellVolumes may prove useful in mappings from cells to faces and back.

  subroutine CellVolume(grid, V_C)

      implicit none
      type (grid_t), intent(in)          :: grid     ! input model
      type (rscalar), intent(inout)      :: V_C       ! cell volume
      type (rscalar)                     :: gV_C  ! temporary vector in global coord system
      ! local variables
      integer                            :: nx,ny,nz
      real(8),dimension(:),allocatable   :: x,y,z
      real(8)                            :: volume
      integer                            :: i, j, k, istat

      write(*,*) 'Computing cell volume elements...'

      if (.not. V_C%allocated) then
        call create_rscalar(grid, V_C, CENTER)
      else
        ! Checks that the size and type are the same
        if (((grid%nx == V_C%nx).and.(grid%ny == V_C%ny).and.(grid%nz == V_C%nz)) &
            .and. (V_C%gridType == CENTER)) then
            ! vector adequately allocated
        else
            ! reallocate
            call deall_rscalar(V_C)
            call create_rscalar(grid, V_C, CENTER)
        end if
      end if

      ! flip the x & y directions to use the conventions in the spherical code
      call coordFlip_rscalar(gV_C, V_C, Spherical)

      nx = grid%ny !grid%nx
      ny = grid%nx !grid%ny
      nz = grid%nz

      ! need to be very careful here:
      ! x historically denotes *interval* in longitude;
      ! y denotes co-latitude from North to South pole;
      ! z is +ve down but denotes radii from the Earth's center
      ! moreover, all elements assume that dual grid is shifted
      ! in the positive direction - careful with indices
      allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
      x= grid%dy*pi/180.                !x = grid%x
      y= (90.- grid%xEdge)*pi/180.      !y = grid%y
      z = EARTH_R * KM2M - grid%zEdge           !z = grid%z


      ! cell volumes are computed in all cells
      ! & need to be careful on a global grid
      do k=1,nz
        do j=1,ny
            do i=1,nx
               call volume_vijk(i,j,k,x,y,z,volume)
               gV_C%v(i, j, k) = volume
            enddo
        enddo
      enddo

      deallocate(x,y,z)

      ! flip the coordinate directions back to MT code notations
      call coordFlip_rscalar(V_C, gV_C, Cartesian)

  end subroutine CellVolume  ! CellVolume

  ! *************************************************************************
  ! * EdgeLength creates line elements defined on edges of the primary grid.
  ! * Edge length elements are defined on interior and boundary edges.

  subroutine EdgeLength(grid,l_E)

      type(grid_t), intent(in)      :: grid
      type(rvector), intent(inout)  :: l_E
      type(rvector)                 :: gl_E  ! temporary vector in global coord system
      ! local variables
      integer                   :: nx,ny,nz
      real(8),dimension(:),allocatable      :: x,y,z
      real(8)                   :: xlen,ylen,zlen
      integer                   :: ic,i,j,k,ii,istat

      write(*,*) 'Computing edge length elements...'

      call create_rvector(grid, l_E, EDGE)

      ! flip the x & y directions to use the conventions in the spherical code
      call coordFlip_rvector(gl_E, l_E, Spherical)

      nx = grid%ny !grid%nx
      ny = grid%nx !grid%ny
      nz = grid%nz

      ! need to be very careful here:
      ! x historically denotes *interval* in longitude;
      ! y denotes co-latitude from North to South pole;
      ! z is +ve down but denotes radii from the Earth's center
      ! moreover, all elements assume that dual grid is shifted
      ! in the positive direction - careful with indices
      allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
      x= grid%dy*pi/180.                !x = grid%x
      y= (90.- grid%xEdge)*pi/180.      !y = grid%y
      z = EARTH_R * KM2M - grid%zEdge           !z = grid%z

      do i=1,nx
        do k=1,nz+1
          do j=1,ny+1
            call leng_xijk(i,j,k,x,y,z,xlen)
            gl_E%x(i,j,k)=xlen
          end do
        end do
      end do

      do j=1,ny
        do k=1,nz+1
          call leng_yijk(j,k,y,z,ylen)
          do i=1,nx+1
            gl_E%y(i,j,k)=ylen
          end do
        end do
      end do

      do k=1,nz
        call leng_zijk(k,z,zlen)
        do j=1,ny+1
          do i=1,nx+1
            gl_E%z(i,j,k)=zlen
          end do
        end do
      end do

      deallocate(x,y,z)

      ! flip the coordinate directions back to MT code notations
      call coordFlip_rvector(l_E, gl_E, Cartesian)

  end subroutine EdgeLength

  ! *************************************************************************
  ! * DualEdgeLength creates line elements defined on faces of the primary grid.
  ! * These line elements are perpendicular to a face with center coinciding
  ! * with the face center. They correspond to the edges of the dual grid.
  ! *
  ! * Face length elements are defined on interior faces only.

  subroutine DualEdgeLength(grid,l_F)

      type(grid_t), intent(in)      :: grid
      type(rvector), intent(inout)  :: l_F
      type(rvector)                 :: gl_F  ! temporary vector in global coord system
      ! local variables
      integer                   :: nx,ny,nz
      real(8),dimension(:),allocatable      :: x,y,z
      real(8)                   :: xlen,ylen,zlen
      integer                   :: ic,i,j,k,ii,istat

      write(*,*) 'Computing face length elements...'

      call create_rvector(grid, l_F, FACE)

      ! flip the x & y directions to use the conventions in the spherical code
      call coordFlip_rvector(gl_F, l_F, Spherical)

      nx = grid%ny !grid%nx
      ny = grid%nx !grid%ny
      nz = grid%nz

      ! need to be very careful here:
      ! x historically denotes *interval* in longitude;
      ! y denotes co-latitude from North to South pole;
      ! z is +ve down but denotes radii from the Earth's center
      ! moreover, all elements assume that dual grid is shifted
      ! in the positive direction - careful with indices
      allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
      x= grid%dy*pi/180.                !x = grid%x
      y= (90.- grid%xEdge)*pi/180.      !y = grid%y
      z = EARTH_R * KM2M - grid%zEdge           !z = grid%z


      ! edges of the dual grid in the longitudinal direction depend
      ! on latitude, and also on longitude for non-uniform grids
      do k=1,nz
        do j=1,ny
          do i=1,nx+1
            call leng_xijk2(i-1,i,j,k,x,y,z,xlen)
            gl_F%x(i,j,k)=xlen
          end do
        end do
      end do

      ! edges of the dual grid in the longitudinal direction
      ! are all the same for a chosen longitude and depth;
      ! they are undefined (zero) at the poles
      do k=1,nz
        do i=1,nx
          do j=1,ny+1
            call leng_yijk2(j-1,k,y,z,ylen)
            gl_F%y(i,j,k)=ylen
          end do
        end do
      end do

      ! vertical edges of the dual grid are defined at the poles
      ! with the usual formula, but not on the upper and lower
      ! boundaries of the domain
      do j=1,ny
        do i=1,nx
          do k=1,nz+1
            call leng_zijk2(k-1,z,zlen)
            gl_F%z(i,j,k) = zlen
          end do
        end do
      end do

      deallocate(x,y,z)

      ! flip the coordinate directions back to MT code notations
      call coordFlip_rvector(l_F, gl_F, Cartesian)

      return
      end subroutine DualEdgeLength

  ! ***************************************************************************
  ! * DualFaceArea: surface area elements perpendicular to the edges of the primary
  ! * grid, with indices matching the primary grid edges. These correspond to
  ! * faces of the dual grid.
  ! *
  ! * Edge surface area elements are defined on interior edges only.

  subroutine DualFaceArea(grid,S_E)

      type(grid_t), intent(in)      :: grid
      type(rvector), intent(inout)  :: S_E
      type(rvector)                 :: gS_E  ! temporary vector in global coord system
      ! local variables
      integer                   :: nx,ny,nz
      real(8),dimension(:),allocatable      :: x,y,z
      real(8)                   :: sijk2,sjki2,skij2
      real(8)                   :: zm,ym,yp,zp
      integer                   :: ic,i,j,k,ii,istat

      write(*,*) 'Computing edge area elements...'

      call create_rvector(grid, S_E, EDGE)

      ! flip the x & y directions to use the conventions in the spherical code
      call coordFlip_rvector(gS_E, S_E, Spherical)

      nx = grid%ny !grid%nx
      ny = grid%nx !grid%ny
      nz = grid%nz

      ! need to be very careful here:
      ! x historically denotes *interval* in longitude;
      ! y denotes co-latitude from North to South pole;
      ! z is +ve down but denotes radii from the Earth's center
      ! moreover, all elements assume that dual grid is shifted
      ! in the positive direction - careful with indices
      allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
      x= grid%dy*pi/180.                !x = grid%x
      y= (90.- grid%xEdge)*pi/180.      !y = grid%y
      z = EARTH_R * KM2M - grid%zEdge           !z = grid%z


      ! There must be paddle wheels of faces at j=1 and j=ny+1
      ! that are used to compute Hr at the poles; these are also
      ! defined for the dual grid. They likely won't be used;
      ! instead we'll want to define the North and South cap volumes.
      ! Still, need to be careful there.
      do k=1,nz+1
        do j=1,ny+1
          call area_sijk2(j-1,k-1,y,z,sijk2)
          do i=1,nx
            gS_E%x(i,j,k) = sijk2
          end do
        end do
      end do

      do k=1,nz+1
        do j=1,ny
          do i=1,nx+1
            call area_sjki2(i-1,i,j,k-1,x,y,z,sjki2)
            gS_E%y(i,j,k) = sjki2
          end do
        end do
      end do

      ! The uppermost face of the dual grid is already within the domain
      do k=1,nz
        do j=1,ny+1
          do i=1,nx+1
            call area_skij2(i-1,i,j-1,k,x,y,z,skij2)
            gS_E%z(i,j,k) = skij2
          end do
        end do
      end do

      deallocate(x,y,z)

      ! flip the coordinate directions back to MT code notations
      call coordFlip_rvector(S_E, gS_E, Cartesian)

      return
      end subroutine DualFaceArea

  ! *************************************************************************
  ! * FaceArea computes surface area elements on faces of the primary grid.
  ! * Face surface area elements are defined on interior and boundary faces.

  subroutine FaceArea(grid,S_F)

      type(grid_t), intent(in)      :: grid
      type(rvector), intent(inout)  :: S_F
      type(rvector)                 :: gS_F  ! temporary vector in global coord system
      ! local variables
      integer                   :: nx,ny,nz
      real(8),dimension(:),allocatable      :: x,y,z
      real(8)                   :: sijk,sjki,skij
      integer                   :: i,j,k,istat

      write(*,*) 'Computing face area elements...'

      call create_rvector(grid, S_F, FACE)

      ! flip the x & y directions to use the conventions in the spherical code
      call coordFlip_rvector(gS_F, S_F, Spherical)

      nx = grid%ny !grid%nx
      ny = grid%nx !grid%ny
      nz = grid%nz

      ! need to be very careful here:
      ! x historically denotes *interval* in longitude;
      ! y denotes co-latitude from North to South pole;
      ! z is +ve down but denotes radii from the Earth's center
      ! moreover, all elements assume that dual grid is shifted
      ! in the positive direction - careful with indices
      allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
      x= grid%dy*pi/180.                !x = grid%x
      y= (90.- grid%xEdge)*pi/180.      !y = grid%y
      z = EARTH_R * KM2M - grid%zEdge           !z = grid%z


      ! surface areas perpendicular to the dual edges
      ! in the longitudinal direction are all the same
      ! for a chosen latitude and depth;
      ! j=1 is North pole, j=ny is South pole, both are defined;
      ! this also takes care of zero longitude
      do k=1,nz
        do j=1,ny
          call area_sijk(j,k,y,z,sijk)
          do i=1,nx+1
            gS_F%x(i,j,k) = sijk
          end do
        end do
      end do

      ! surface areas perpendicular to the dual edges
      ! in the latitudinal direction depend on longitude
      ! for non-uniform grids
      do k=1,nz
        do i=1,nx
          do j=1,ny+1
            call area_sjki(i,j,k,x,y,z,sjki)
            gS_F%y(i,j,k) = sjki
          end do
        end do
      end do

      ! horizontal surface areas are defined at the poles
      ! and on upper and lower boundaries too; the usual
      ! formula may be used at the poles
      do j=1,ny
        do i=1,nx
          do k=1,nz+1
            call area_skij(i,j,k,x,y,z,skij)
            gS_F%z(i,j,k) = skij
          end do
        end do
      end do

      deallocate(x,y,z)

      ! flip the coordinate directions back to MT code notations
      call coordFlip_rvector(S_F, gS_F, Cartesian)

      return
      end subroutine FaceArea

  ! *************************************************************************
  ! * UNWEIGHTED MAPPING OPERATORS THAT ALSO COMPUTE BOUNDARY EDGES (TESTED)
  ! * SAME FOR SPHERICAL AND CARTESIAN COORDINATES
  ! * (NOTE: most likely, this will also work for multigrid: we would use
  ! * separate subroutines to map from egdes to multigrid edges and back)
  ! *************************************************************************

  ! *************************************************************************
  ! * Cell2Edge will be used by forward model mappings

  subroutine Cell2Edge(grid,C,E)

      type(grid_t), intent(in)      :: grid
      type(rscalar), intent(in)     :: C
      type(rvector), intent(out)    :: E
      ! local variables
      integer                   :: ix,iy,iz

      call create_rvector(grid, E, EDGE)

      ! for x-components inside the domain
      do ix = 1,grid%nx
         do iy = 2,grid%ny
            do iz = 2,grid%nz

               E%x(ix, iy, iz) = (C%v(ix, iy-1, iz-1) + C%v(ix, iy, iz-1) + &
                                  C%v(ix, iy-1, iz) + C%v(ix, iy, iz))/4.0d0

            enddo
         enddo
      enddo

      ! for y-components inside the domain
      do ix = 2,grid%nx
         do iy = 1,grid%ny
            do iz = 2,grid%nz

               E%y(ix, iy, iz) = (C%v(ix-1, iy, iz-1) + C%v(ix, iy, iz-1) + &
                                  C%v(ix-1, iy, iz) + C%v(ix, iy, iz))/4.0d0

            enddo
         enddo
      enddo

      ! for z-components inside the domain
      do ix = 2,grid%nx
         do iy = 2,grid%ny
            do iz = 1,grid%nz

               E%z(ix, iy, iz) = (C%v(ix-1, iy-1, iz) + C%v(ix-1, iy, iz) + &
                                  C%v(ix, iy-1, iz) + C%v(ix, iy, iz))/4.0d0

            enddo
         enddo
      enddo

      ! upper boundary
      iz = 1
      do iy = 1,grid%ny
         do ix = 2,grid%nx
            E%y(ix, iy, iz) = (C%v(ix-1, iy, iz) + C%v(ix, iy, iz))/2.0d0
         enddo
      enddo
      do ix = 1,grid%nx
         do iy = 2,grid%ny
            E%x(ix, iy, iz) = (C%v(ix, iy-1, iz) + C%v(ix, iy, iz))/2.0d0
         enddo
      enddo

      ! lower boundary
      iz = grid%nz+1
      do iy = 1,grid%ny
         do ix = 2,grid%nx
            E%y(ix, iy, iz) = (C%v(ix-1, iy, iz-1) + C%v(ix, iy, iz-1))/2.0d0
         enddo
      enddo
      do ix = 1,grid%nx
         do iy = 2,grid%ny
            E%x(ix, iy, iz) = (C%v(ix, iy-1, iz-1) + C%v(ix, iy, iz-1))/2.0d0
         enddo
      enddo

      ! side boundaries: left x-side
      ix = 1
      do iy = 1,grid%ny
         do iz = 2,grid%nz
            E%y(ix, iy, iz) = (C%v(ix, iy, iz-1) + C%v(ix, iy, iz))/2.0d0
         enddo
      enddo
      do iz = 1,grid%nz
        do iy = 2,grid%ny
            E%z(ix, iy, iz) = (C%v(ix, iy-1, iz) + C%v(ix, iy, iz))/2.0d0
        enddo
      enddo

      ! side boundaries: right x-side
      ix = grid%nx+1
      do iy = 1,grid%ny
         do iz = 2,grid%nz
            E%y(ix, iy, iz) = (C%v(ix-1, iy, iz-1) + C%v(ix-1, iy, iz))/2.0d0
         enddo
      enddo
      do iz = 1,grid%nz
         do iy = 2,grid%ny
            E%z(ix, iy, iz) = (C%v(ix-1, iy-1, iz) + C%v(ix-1, iy, iz))/2.0d0
         enddo
      enddo

      ! side boundaries: left y-side
      iy = 1
      do ix = 1,grid%nx
         do iz = 2,grid%nz
            E%x(ix, iy, iz) = (C%v(ix, iy, iz-1) + C%v(ix, iy, iz))/2.0d0
         enddo
      enddo
      do iz = 1,grid%nz
         do ix = 2,grid%nx
            E%z(ix, iy, iz) = (C%v(ix-1, iy, iz) + C%v(ix, iy, iz))/2.0d0
         enddo
      enddo

      ! side boundaries: right y-side
      iy = grid%ny+1
      do ix = 1,grid%nx
         do iz = 2,grid%nz
            E%x(ix, iy, iz) = (C%v(ix, iy-1, iz-1) + C%v(ix, iy-1, iz))/2.0d0
         enddo
      enddo
      do iz = 1,grid%nz
         do ix = 2,grid%nx
            E%z(ix, iy, iz) = (C%v(ix-1, iy-1, iz) + C%v(ix, iy-1, iz))/2.0d0
         enddo
      enddo

      ! z-component corners
      do iz = 1,grid%nz
         ix = 1; iy = 1
         E%z(ix, iy, iz) = C%v(ix, iy, iz)
         ix = 1; iy = grid%ny+1
         E%z(ix, iy, iz) = C%v(ix, iy-1, iz)
         ix = grid%nx+1; iy = 1
         E%z(ix, iy, iz) = C%v(ix-1, iy, iz)
         ix = grid%nx+1; iy = grid%ny+1
         E%z(ix, iy, iz) = C%v(ix-1, iy-1, iz)
      enddo

      ! y-component corners
      do iy = 1,grid%ny
         ix = 1; iz = 1
         E%y(ix, iy, iz) = C%v(ix, iy, iz)
         ix = 1; iz = grid%nz+1
         E%y(ix, iy, iz) = C%v(ix, iy, iz-1)
         ix = grid%nx+1; iz = 1
         E%y(ix, iy, iz) = C%v(ix-1, iy, iz)
         ix = grid%nx+1; iz = grid%nz+1
         E%y(ix, iy, iz) = C%v(ix-1, iy, iz-1)
      enddo

      ! x-component corners
      do ix = 1,grid%nx
         iy = 1; iz = 1
         E%x(ix, iy, iz) = C%v(ix, iy, iz)
         iy = 1; iz = grid%nz+1
         E%x(ix, iy, iz) = C%v(ix, iy, iz-1)
         iy = grid%ny+1; iz = 1
         E%x(ix, iy, iz) = C%v(ix, iy-1, iz)
         iy = grid%ny+1; iz = grid%nz+1
         E%x(ix, iy, iz) = C%v(ix, iy-1, iz-1)
      enddo

  end subroutine Cell2Edge

  ! *************************************************************************
  ! * Cell2Node will be used by modified system equation with "node based
  ! * scaling factor
  ! * used by ModelParamToNode (currently only employed in SP2)
  ! * essentially this maps cell based scalar value onto Nodes.  

  subroutine Cell2Node(grid,C,N)
      type(grid_t), intent(in)      :: grid
      type(rscalar), intent(in)     :: C
      type(rscalar), intent(out)    :: N
      ! local variables
      integer                   :: ix,iy,iz
      call create_rscalar(grid, N, CORNER) ! node-based
      ! non-boundary nodes
      do ix = 2,grid%nx
         do iy = 2,grid%ny
            do iz = 2,grid%nz
            ! i-1,j-1,k-1
               N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy-1,iz-1)/8.0d0
            ! i-1,j-1,k
               N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy-1, iz)/8.0d0
            ! i-1,j,k-1
               N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy, iz-1)/8.0d0
            ! i-1,j,k
               N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy, iz)/8.0d0
            ! i,j-1,k-1
               N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy-1, iz-1)/8.0d0
            ! i,j-1,k
               N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy-1, iz)/8.0d0
            ! i,j,k-1
               N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy, iz-1)/8.0d0
            ! i,j,k
               N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy, iz)/8.0d0
            enddo
         enddo
      enddo

      ! left boundary
      iy = 1
      do ix = 2,grid%nx
          do iz = 2,grid%nz
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy, iz-1)/4.0d0
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy, iz-1)/4.0d0
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy, iz)/4.0d0
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy, iz)/4.0d0
          enddo
      enddo

      ! right boundary
      iy = grid%ny+1
      do ix = 2,grid%nx
          do iz = 2,grid%nz
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy-1, iz-1)/4.0d0
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy-1, iz-1)/4.0d0
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy-1, iz)/4.0d0
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy-1, iz)/4.0d0
          enddo
      enddo


      ! front boundary
      ix = 1
      do iy = 2,grid%ny
          do iz = 2,grid%nz
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy-1, iz-1)/4.0d0
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy, iz-1)/4.0d0
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy-1, iz)/4.0d0
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy, iz)/4.0d0
          enddo
      enddo
      ! back boundary
      ix = grid%nx+1
      do iy = 2,grid%ny
          do iz = 2,grid%nz
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy-1, iz-1)/4.0d0
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy, iz-1)/4.0d0
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy-1, iz)/4.0d0
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy, iz)/4.0d0
          enddo
      enddo

      ! upper boundary
      iz = 1
      do iy = 2,grid%ny
          do ix = 2,grid%nx
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy-1, iz)/4.0d0
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy-1, iz)/4.0d0
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy, iz)/4.0d0
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy, iz)/4.0d0
          enddo
      enddo

      ! lower boundary
      iz = grid%nz+1
      do iy = 2,grid%ny
          do ix = 2,grid%nx
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy-1, iz-1)/4.0d0
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy-1, iz-1)/4.0d0
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy, iz-1)/4.0d0
            N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy, iz-1)/4.0d0
          enddo
      enddo
      ! 12 edges of the model domain
      ! average of two cells
      ! x edges
      do ix = 2, grid%nx
          ! upper, left
          iy = 1
          iz = 1
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy, iz)/2.0d0
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy, iz)/2.0d0
          ! lower, left
          iz = grid%nz+1
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy, iz-1)/2.0d0
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy, iz-1)/2.0d0
          ! upper, right
          iy = grid%ny+1
          iz = 1
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy-1, iz)/2.0d0
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy-1, iz)/2.0d0
          ! lower, right
          iz = grid%nz+1
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy-1, iz-1)/2.0d0
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy-1, iz-1)/2.0d0
      enddo
      ! y edges
      do iy = 2, grid%ny
          ! upper, front
          ix = 1
          iz = 1
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy-1, iz)/2.0d0
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy, iz)/2.0d0
          ! lower, front
          iz = grid%nz+1
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy-1, iz-1)/2.0d0
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy, iz-1)/2.0d0
          ! upper, back
          ix = grid%nx+1
          iz = 1
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy-1, iz)/2.0d0
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy, iz)/2.0d0
          ! lower, back
          iz = grid%nz+1
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy-1, iz-1)/2.0d0
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy, iz-1)/2.0d0
      enddo
      ! z edges
      do iz = 2, grid%nz
          ! front, left
          ix = 1
          iy = 1
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy, iz-1)/2.0d0
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy, iz)/2.0d0
          ! front, right
          iy = grid%ny+1
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy-1, iz-1)/2.0d0
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix, iy-1, iz)/2.0d0
          ! back, left
          ix = grid%nx+1
          iy = 1
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy, iz-1)/2.0d0
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy, iz)/2.0d0
          ! back, right
          iy = grid%ny+1
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy-1, iz-1)/2.0d0
          N%v(ix, iy, iz) = N%v(ix, iy, iz) + C%v(ix-1, iy-1, iz)/2.0d0
      enddo
      ! 8 corners of the model domain
      ! front, left, upper
      N%v(1,1,1) = C%v(1,1,1)
      ! back, left, upper
      N%v(grid%nx+1,1,1) = C%v(grid%nx,1,1)
      ! front, right, upper
      N%v(1,grid%ny+1,1) = C%v(1,grid%ny,1)
      ! back, right, upper
      N%v(grid%nx+1,grid%ny+1,1) = C%v(grid%nx,grid%ny,1)
      ! front, left, lower
      N%v(1,1,grid%nz+1) = C%v(1,1,grid%nz)
      ! back, left, lower
      N%v(grid%nx+1,1,grid%nz+1) = C%v(grid%nx,1,grid%nz)
      ! front, right, lower
      N%v(1,grid%ny+1,grid%nz+1) = C%v(1,grid%ny,grid%nz)
      ! back, right, upper
      N%v(grid%nx+1,grid%ny+1,grid%nz+1) = C%v(grid%nx,grid%ny,grid%nz)
      !
  end subroutine Cell2Node
  ! *************************************************************************
  ! * Edge2Cell will be used by adjoint model mappings

  subroutine Edge2Cell(grid,E,C)

      type(grid_t), intent(in)      :: grid
      type(rvector), intent(in)     :: E
      type(rscalar), intent(out)    :: C
      ! local variables
      integer                   :: ix,iy,iz

      call create_rscalar(grid, C, CENTER)

      ! for x-components inside the domain
      do ix = 1,grid%nx
         do iy = 2,grid%ny
            do iz = 2,grid%nz

               C%v(ix, iy-1, iz-1) = C%v(ix, iy-1, iz-1) + E%x(ix, iy, iz)/4.0d0
               C%v(ix, iy, iz-1) = C%v(ix, iy, iz-1) + E%x(ix, iy, iz)/4.0d0
               C%v(ix, iy-1, iz) = C%v(ix, iy-1, iz) + E%x(ix, iy, iz)/4.0d0
               C%v(ix, iy, iz) = C%v(ix, iy, iz) + E%x(ix, iy, iz)/4.0d0

            enddo
         enddo
      enddo

      ! for y-components inside the domain
      do ix = 2,grid%nx
         do iy = 1,grid%ny
            do iz = 2,grid%nz

               C%v(ix-1, iy, iz-1) = C%v(ix-1, iy, iz-1) + E%y(ix, iy, iz)/4.0d0
               C%v(ix, iy, iz-1) = C%v(ix, iy, iz-1) + E%y(ix, iy, iz)/4.0d0
               C%v(ix-1, iy, iz) = C%v(ix-1, iy, iz) + E%y(ix, iy, iz)/4.0d0
               C%v(ix, iy, iz) = C%v(ix, iy, iz) + E%y(ix, iy, iz)/4.0d0

            enddo
         enddo
      enddo

      ! for z-components inside the domain
      do ix = 2,grid%nx
         do iy = 2,grid%ny
            do iz = 1,grid%nz

               C%v(ix-1, iy-1, iz) = C%v(ix-1, iy-1, iz) + E%z(ix, iy, iz)/4.0d0
               C%v(ix-1, iy, iz) = C%v(ix-1, iy, iz) + E%z(ix, iy, iz)/4.0d0
               C%v(ix, iy-1, iz) = C%v(ix, iy-1, iz) + E%z(ix, iy, iz)/4.0d0
               C%v(ix, iy, iz) = C%v(ix, iy, iz) + E%z(ix, iy, iz)/4.0d0

            enddo
         enddo
      enddo

      ! upper boundary
      iz = 1
      do iy = 1,grid%ny
         do ix = 2,grid%nx
            C%v(ix-1, iy, iz) = C%v(ix-1, iy, iz) + E%y(ix, iy, iz)/2.0d0
            C%v(ix, iy, iz) = C%v(ix, iy, iz) + E%y(ix, iy, iz)/2.0d0
         enddo
      enddo
      do ix = 1,grid%nx
         do iy = 2,grid%ny
            C%v(ix, iy-1, iz) = C%v(ix, iy-1, iz) + E%x(ix, iy, iz)/2.0d0
            C%v(ix, iy, iz) = C%v(ix, iy, iz) + E%x(ix, iy, iz)/2.0d0
         enddo
      enddo

      ! lower boundary
      iz = grid%nz+1
      do iy = 1,grid%ny
         do ix = 2,grid%nx
            C%v(ix-1, iy, iz-1) = C%v(ix-1, iy, iz-1) + E%y(ix, iy, iz)/2.0d0
            C%v(ix, iy, iz-1) = C%v(ix, iy, iz-1) + E%y(ix, iy, iz)/2.0d0
         enddo
      enddo
      do ix = 1,grid%nx
         do iy = 2,grid%ny
            C%v(ix, iy-1, iz-1) = C%v(ix, iy-1, iz-1) + E%x(ix, iy, iz)/2.0d0
            C%v(ix, iy, iz-1) = C%v(ix, iy, iz-1) + E%x(ix, iy, iz)/2.0d0
         enddo
      enddo

      ! side boundaries: left x-side
      ix = 1
      do iy = 1,grid%ny
         do iz = 2,grid%nz
            C%v(ix, iy, iz-1) = C%v(ix, iy, iz-1) + E%y(ix, iy, iz)/2.0d0
            C%v(ix, iy, iz) = C%v(ix, iy, iz) + E%y(ix, iy, iz)/2.0d0
         enddo
      enddo
      do iz = 1,grid%nz
         do iy = 2,grid%ny
            C%v(ix, iy-1, iz) = C%v(ix, iy-1, iz) + E%z(ix, iy, iz)/2.0d0
            C%v(ix, iy, iz) = C%v(ix, iy, iz) + E%z(ix, iy, iz)/2.0d0
         enddo
      enddo

      ! side boundaries: right x-side
      ix = grid%nx+1
      do iy = 1,grid%ny
         do iz = 2,grid%nz
            C%v(ix-1, iy, iz-1) = C%v(ix-1, iy, iz-1) + E%y(ix, iy, iz)/2.0d0
            C%v(ix-1, iy, iz) = C%v(ix-1, iy, iz) + E%y(ix, iy, iz)/2.0d0
         enddo
      enddo
      do iz = 1,grid%nz
         do iy = 2,grid%ny
            C%v(ix-1, iy-1, iz) = C%v(ix-1, iy-1, iz) + E%z(ix, iy, iz)/2.0d0
            C%v(ix-1, iy, iz) = C%v(ix-1, iy, iz) + E%z(ix, iy, iz)/2.0d0
         enddo
      enddo

      ! side boundaries: left y-side
      iy = 1
      do ix = 1,grid%nx
         do iz = 2,grid%nz
            C%v(ix, iy, iz-1) = C%v(ix, iy, iz-1) + E%x(ix, iy, iz)/2.0d0
            C%v(ix, iy, iz) = C%v(ix, iy, iz) + E%x(ix, iy, iz)/2.0d0
         enddo
      enddo
      do iz = 1,grid%nz
         do ix = 2,grid%nx
            C%v(ix-1, iy, iz) = C%v(ix-1, iy, iz) + E%z(ix, iy, iz)/2.0d0
            C%v(ix, iy, iz) = C%v(ix, iy, iz) + E%z(ix, iy, iz)/2.0d0
         enddo
      enddo

      ! side boundaries: right y-side
      iy = grid%ny+1
      do ix = 1,grid%nx
         do iz = 2,grid%nz
            C%v(ix, iy-1, iz-1) = C%v(ix, iy-1, iz-1) + E%x(ix, iy, iz)/2.0d0
            C%v(ix, iy-1, iz) = C%v(ix, iy-1, iz) + E%x(ix, iy, iz)/2.0d0
         enddo
      enddo
      do iz = 1,grid%nz
         do ix = 2,grid%nx
            C%v(ix-1, iy-1, iz) = C%v(ix-1, iy-1, iz) + E%z(ix, iy, iz)/2.0d0
            C%v(ix, iy-1, iz) = C%v(ix, iy-1, iz) + E%z(ix, iy, iz)/2.0d0
         enddo
      enddo

      ! z-component corners
      do iz = 1,grid%nz
         ix = 1; iy = 1
         C%v(ix, iy, iz) = C%v(ix, iy, iz) + E%z(ix, iy, iz)
         ix = 1; iy = grid%ny+1
         C%v(ix, iy-1, iz) = C%v(ix, iy-1, iz) + E%z(ix, iy, iz)
         ix = grid%nx+1; iy = 1
         C%v(ix-1, iy, iz) = C%v(ix-1, iy, iz) + E%z(ix, iy, iz)
         ix = grid%nx+1; iy = grid%ny+1
         C%v(ix-1, iy-1, iz) = C%v(ix-1, iy-1, iz) + E%z(ix, iy, iz)
      enddo

      ! y-component corners
      do iy = 1,grid%ny
         ix = 1; iz = 1
         C%v(ix, iy, iz) = C%v(ix, iy, iz) + E%y(ix, iy, iz)
         ix = 1; iz = grid%nz+1
         C%v(ix, iy, iz-1) = C%v(ix, iy, iz-1) + E%y(ix, iy, iz)
         ix = grid%nx+1; iz = 1
         C%v(ix-1, iy, iz) = C%v(ix-1, iy, iz) + E%y(ix, iy, iz)
         ix = grid%nx+1; iz = grid%nz+1
         C%v(ix-1, iy, iz-1) = C%v(ix-1, iy, iz-1) + E%y(ix, iy, iz)
      enddo

      ! x-component corners
      do ix = 1,grid%nx
         iy = 1; iz = 1
         C%v(ix, iy, iz) = C%v(ix, iy, iz) + E%x(ix, iy, iz)
         iy = 1; iz = grid%nz+1
         C%v(ix, iy, iz-1) = C%v(ix, iy, iz-1) + E%x(ix, iy, iz)
         iy = grid%ny+1; iz = 1
         C%v(ix, iy-1, iz) = C%v(ix, iy-1, iz) + E%x(ix, iy, iz)
         iy = grid%ny+1; iz = grid%nz+1
         C%v(ix, iy-1, iz-1) = C%v(ix, iy-1, iz-1) + E%x(ix, iy, iz)
      enddo


  end subroutine Edge2Cell

  ! *************************************************************************
  ! * create_gridElements(grid) is the master call to initialize all arrays

  subroutine create_gridElements(grid)
      implicit none
      type (grid_t), intent(in)           :: grid     ! input model

      ! First deallocate, if already exist
      call deall_gridElements()

      ! Now, initialize all grid elements using the input grid
      call FaceArea(grid, S_F)
      call DualFaceArea(grid, S_E)
      call EdgeLength(grid, l_E)
      call DualEdgeLength(grid, l_F)
      call NodeVolume(grid, V_N)
      call EdgeVolume(grid, V_E, l_E, S_E)
      call CellVolume(grid, V_C)

  end subroutine create_gridElements

  ! *************************************************************************
  ! * deall_gridElements() is the master call to deallocate all saved arrays
  subroutine deall_gridElements()

      call deall_rvector(S_F)
      call deall_rvector(S_E)
      call deall_rvector(l_E)
      call deall_rvector(l_F)
      call deall_rvector(V_E)
      call deall_rscalar(V_N)
      call deall_rscalar(V_C)

   end subroutine deall_gridElements

end module GridCalc
