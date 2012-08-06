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
  use elements
  implicit none

  public      :: EdgeVolume, FaceVolume, NodeVolume, CellVolume
  public      :: EdgeLength, FaceLength
  public      :: EdgeArea,   FaceArea

Contains

  ! *************************************************************************
  ! * EdgeVolume creates volume elements centered around the edges of
  ! * the grid, and stores them as real vectors with gridType=EDGE.
  ! *
  ! * A diagonal matrix multiplication of the edge volume with the difference
  ! * equations enables us to make a symmetrical matrix. Remember,
  ! * the electrical fields are defined on the center of the edges, therefore,
  ! * the edge volume is centered about the electrical field measurement.

  subroutine EdgeVolume(grid, V_E)

      implicit none
      type (grid_t), intent(in)           :: grid     ! input model
      type (rvector), intent(inout)       :: V_E       ! edge volume
      type (rvector), intent(in), optional:: l_E, S_E ! optional inputs
      ! local
      type (rvector)                      :: length, area
      logical                             :: compute_elements

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
        call FaceLength(grid,length)
        call FaceArea(grid,area)
      else
        length = l_E
        area = S_E
      end if

      ! compute volume elements
      call diagMult_rvector(l_E,S_E,V_E)

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
        call FaceLength(grid,length)
        call FaceArea(grid,area)
      else
        length = l_F
        area = S_F
      end if

      ! compute volume elements
      call diagMult_rvector(l_F,S_F,V_F)

      call deall_rvector(length)
      call deall_rvector(area)

  end subroutine FaceVolume  ! FaceVolume

  ! *************************************************************************
  ! * NodeVolume creates volume elements centered around the corners (nodes)
  ! * of the grid, and stores them as real scalars with gridType=CORNER.

  subroutine NodeVolume(grid, V_N)

      type (grid_t), intent(in)          :: grid  ! input grid
      type (rscalar), intent(inout)      :: V_N   ! node/corner volume as output
      ! local variables
      integer                            :: nx,ny,nz
      real(8),dimension(:),allocatable   :: x,y,z
      real(8)                            :: volume,sijk2,zm,zp,ym,yp
      integer                            :: i, j, k

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

      nx = grid%nx
      ny = grid%ny
      nz = grid%nz

      ! need to be very careful here:
      ! x historically denotes *interval* in longitude;
      ! y denotes co-latitude from North to South pole;
      ! z is +ve down but denotes radii from the Earth's center
      ! moreover, all elements assume that dual grid is shifted
      ! in the positive direction - careful with indices
      allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
      x = grid%x
      y = grid%y
      z = grid%z

      ! node volumes are only using the internal corner nodes
      ! but need to be careful on a global grid
      do k=2,nz
        ! North pole cap
        j=1
        zm=(z(k-1)+z(k))/2.0d0
        zp=(z(k)+z(k+1))/2.0d0
        yp=(y(j+1))/2.0d0
        sijk2=( (zm**2)-(zp**2) )*( yp )/2.d0
        V_N%v(:,j,k) = 2.0d0*PI*sijk2
        do j=2,ny
            ! Zero longitude
            i=1
            call volume_vijk2(nx,1,j-1,k-1,x,y,z,volume)
            V_N%v(i, j, k) = volume
            do i=2,nx
               call volume_vijk2(i-1,i,j-1,k-1,x,y,z,volume)
               V_N%v(i, j, k) = volume
            enddo
            ! Zero longitude
            i=nx+1
            call volume_vijk2(nx,1,j-1,k-1,x,y,z,volume)
            V_N%v(i, j, k) = volume
        enddo
        ! South pole cap
        j=ny+1
        zm=(z(k-1)+z(k))/2.0d0
        zp=(z(k)+z(k+1))/2.0d0
        ym=(y(j)-y(j-1))/2.0d0
        sijk2=( (zm**2)-(zp**2) )*( ym )/2.d0
        V_N%v(:,j,k) = 2.0d0*PI*sijk2
      enddo

      deallocate(x,y,z)

      call validate_rscalar(V_N,.true.)

  end subroutine NodeVolume

  ! *************************************************************************
  ! * CellVolumes may prove useful in mappings from cells to faces and back.

  subroutine CellVolume(grid, V_C)

      implicit none
      type (grid_t), intent(in)          :: grid     ! input model
      type (rvector), intent(inout)      :: V_C       ! cell volume
      ! local variables
      integer                            :: nx,ny,nz
      real(8),dimension(:),allocatable   :: x,y,z
      real(8)                            :: volume
      integer                            :: i, j, k

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

      nx = grid%nx
      ny = grid%ny
      nz = grid%nz

      ! need to be very careful here:
      ! x historically denotes *interval* in longitude;
      ! y denotes co-latitude from North to South pole;
      ! z is +ve down but denotes radii from the Earth's center
      ! moreover, all elements assume that dual grid is shifted
      ! in the positive direction - careful with indices
      allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
      x = grid%x
      y = grid%y
      z = grid%z

      ! cell volumes are computed in all cells
      ! & need to be careful on a global grid
      do k=1,nz
        do j=1,ny
            do i=1,nx
               call volume_vijk(i,j,k,x,y,z,volume)
               V_C%v(i, j, k) = volume
            enddo
        enddo
      enddo

      deallocate(x,y,z)

      call validate_rscalar(V_C,.true.)

  end subroutine CellVolume  ! CellVolume

  ! *************************************************************************
  ! * EdgeLength creates line elements defined on edges of the primary grid.
  ! * Edge length elements are defined on interior and boundary edges.

  subroutine EdgeLength(grid,l_F)

      type(grid_t), intent(in)      :: grid
      type(rvector), intent(inout)  :: l_E
      ! local variables
      integer                   :: nx,ny,nz
      real(8),dimension(:),allocatable      :: x,y,z
      real(8)                   :: xlen,ylen,zlen
      integer                   :: ic,i,j,k,ii,istat

      call create_rvector(grid, l_E, EDGE)

      nx = grid%nx
      ny = grid%ny
      nz = grid%nz

      ! need to be very careful here:
      ! x historically denotes *interval* in longitude;
      ! y denotes co-latitude from North to South pole;
      ! z is +ve down but denotes radii from the Earth's center
      ! moreover, all elements assume that dual grid is shifted
      ! in the positive direction - careful with indices
      allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
      x = grid%x
      y = grid%y
      z = grid%z

      do i=1,nx
        do k=1,nz+1
          ! North pole
          j = 1
          l_E%x(i,j,k) = 0.0d0
          do j=2,ny
            call leng_xijk(i,j,k,x,y,z,xlen)
            l_E%x(i,j,k)=xlen
          end do
          ! South pole
          j = ny+1
          l_E%x(i,j,k) = 0.0d0
        end do
      end do

      do j=1,ny
        do k=1,nz+1
          call leng_yijk(j,k,y,z,ylen)
          do i=1,nx+1
            l_E%y(i,j,k)=ylen
          end do
        end do
      end do

      do k=1,nz
        call leng_zijk(k,z,zlen)
        do j=1,ny+1
          do i=1,nx+1
            l_E%z(i,j,k)=zlen
          end do
        end do
      end do

      deallocate(x,y,z)

      call validate_rvector(l_E,.true.)

  end subroutine EdgeLength

  ! *************************************************************************
  ! * FaceLength creates line elements defined on faces of the primary grid.
  ! * These line elements are perpendicular to a face with center coinciding
  ! * with the face center. They correspond to the edges of the dual grid.
  ! *
  ! * Face length elements are defined on interior faces only.

  subroutine FaceLength(grid,l_F)

      type(grid_t), intent(in)      :: grid
      type(rvector), intent(inout)  :: l_F
      ! local variables
      integer                   :: nx,ny,nz
      real(8),dimension(:),allocatable      :: x,y,z
      real(8)                   :: xlen,ylen,zlen
      integer                   :: ic,i,j,k,ii,istat

      call create_rvector(grid, l_F, FACE)

      nx = grid%nx
      ny = grid%ny
      nz = grid%nz

      ! need to be very careful here:
      ! x historically denotes *interval* in longitude;
      ! y denotes co-latitude from North to South pole;
      ! z is +ve down but denotes radii from the Earth's center
      ! moreover, all elements assume that dual grid is shifted
      ! in the positive direction - careful with indices
      allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
      x = grid%x
      y = grid%y
      z = grid%z

      ! edges of the dual grid in the longitudinal direction depend
      ! on latitude, and also on longitude for non-uniform grids
      do k=1,nz
        do j=1,ny
          ! Zero longitude
          i=1
          call leng_xijk2(nx,1,j,k,x,y,z,xlen)
          l_F%x(i,j,k)=xlen
          do i=2,nx+1
            call leng_xijk2(i-1,i,j,k,x,y,z,xlen)
            l_F%x(i,j,k)=xlen
          end do
        end do
      end do

      ! edges of the dual grid in the longitudinal direction
      ! are all the same for a chosen longitude and depth;
      ! they are undefined (zero) at the poles
      do k=1,nz
        do i=1,nx
          ! North pole
          j=1
          l_F%y(i,j,k)=0.0d0
          do j=2,ny
            call leng_yijk2(j-1,k,y,z,ylen)
            l_F%y(i,j,k)=ylen
          end do
          ! South pole
          j=ny+1
          l_F%y(i,j,k)=0.0d0
        end do
      end do

      ! vertical edges of the dual grid are defined at the poles
      ! with the usual formula, but not on the upper and lower
      ! boundaries of the domain
      do j=1,ny
        do i=1,nx
          ! Upper boundary
          k=1
          l_F%z(i,j,k) = 0.0d0
          do k=2,nz
            call leng_zijk2(k-1,z,zlen)
            l_F%z(i,j,k) = zlen
          end do
          ! Lower boundary
          k=nz+1
          l_F%z(i,j,k) = 0.0d0
        end do
      end do

      deallocate(x,y,z)

      call validate_rvector(l_F,.true.)

      return
      end subroutine FaceLength

  ! ***************************************************************************
  ! * EdgeArea: surface area elements perpendicular to the edges of the primary
  ! * grid, with indices matching the primary grid edges. These correspond to
  ! * faces of the dual grid.
  ! *
  ! * Edge surface area elements are defined on interior edges only.

  subroutine EdgeArea(grid,S_E)

      type(grid_t), intent(in)      :: grid
      type(rvector), intent(inout)  :: S_E
      ! local variables
      integer                   :: nx,ny,nz
      real(8),dimension(:),allocatable      :: x,y,z
      real(8)                   :: sijk2,sjki2,skij2
      real(8)                   :: zm,ym,yp,zp
      integer                   :: ic,i,j,k,ii,istat

      call create_rvector(grid, S_E, EDGE)

      nx = grid%nx
      ny = grid%ny
      nz = grid%nz

      ! need to be very careful here:
      ! x historically denotes *interval* in longitude;
      ! y denotes co-latitude from North to South pole;
      ! z is +ve down but denotes radii from the Earth's center
      ! moreover, all elements assume that dual grid is shifted
      ! in the positive direction - careful with indices
      ! Upper and lower boundaries not defined in this routine.
      allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
      x = grid%x
      y = grid%y
      z = grid%z

      ! There must be paddle wheels of faces at j=1 and j=ny+1
      ! that are used to compute Hr at the poles; these are also
      ! defined for the dual grid. They likely won't be used;
      ! instead we'll want to define the North and South cap volumes.
      ! Still, need to be careful there.
      do k=2,nz
        ! North pole cap "paddle wheel"
        j=1
        zm=(z(k-1)+z(k))/2.0d0
        zp=(z(k)+z(k+1))/2.0d0
        yp=(y(j+1))/2.0d0
        sijk2=( (zm**2)-(zp**2) )*( yp )/2.d0
        S_E%x(:,j,k) = sijk2
        do j=2,ny
          call area_sijk2(j-1,k-1,y,z,sijk2)
          do i=1,nx
            S_E%x(i,j,k) = sijk2
          end do
        end do
        ! South pole cap "paddle wheel"
        j=ny+1
        zm=(z(k-1)+z(k))/2.0d0
        zp=(z(k)+z(k+1))/2.0d0
        ym=(y(j)-y(j-1))/2.0d0
        sijk2=( (zm**2)-(zp**2) )*( ym )/2.d0
        S_E%x(:,j,k) = sijk2
      end do

      do k=2,nz
        do j=1,ny
          ! Zero longitude
          i = 1
          call area_sjki2(nx,i,j,k-1,x,y,z,sjki2)
          S_E%y(i,j,k) = sjki2
          do i=2,nx
            call area_sjki2(i-1,i,j,k-1,x,y,z,sjki2)
            S_E%y(i,j,k) = sjki2
          end do
          ! Zero longitude
          i = nx+1
          S_E%y(i,j,k) = S_E%y(1,j,k)
        end do
      end do

      ! The uppermost face of the dual grid is already within the domain
      do k=1,nz
        ! North pole cap
        j=1
        zp=(z(k)+z(k+1))/2.0d0
        yp=(y(j)+y(j+1))/2.0d0
        skij2=2.0d0*pi*(zp**2)*(1.0d0-dcos(yp))
        S_E%z(:,j,k) = skij2
        do j=2,ny
          ! Zero longitude
          call area_skij2(nx,1,j-1,k,x,y,z,skij2)
          S_E%z(i,j,k) = skij2
          do i=2,nx
            call area_skij2(i-1,i,j-1,k,x,y,z,skij2)
            S_E%z(i,j,k) = skij2
          end do
        end do
        ! South pole cap
        j=ny+1
        zp=(z(k)+z(k+1))/2.0d0
        ym=(y(j)+y(j-1))/2.0d0
        skij2=2.0d0*pi*(zp**2)*(dcos(ym)+1.0d0)
        S_E%z(:,j,k) = skij2
      end do

      deallocate(x,y,z)

      call validate_rvector(S_E,.true.)

      return
      end subroutine EdgeArea

  ! *************************************************************************
  ! * FaceArea computes surface area elements on faces of the primary grid.
  ! * Face surface area elements are defined on interior and boundary faces.

  subroutine FaceArea(grid,S_F)

      type(grid_t), intent(in)      :: grid
      type(rvector), intent(inout)  :: S_F
      ! local variables
      integer                   :: nx,ny,nz
      real(8),dimension(:),allocatable      :: x,y,z
      real(8)                   :: sijk,sjki,skij
      integer                   :: i,j,k

      call create_rvector(grid, S_F, FACE)

      nx = grid%nx
      ny = grid%ny
      nz = grid%nz

      ! need to be very careful here:
      ! x historically denotes *interval* in longitude;
      ! y denotes co-latitude from North to South pole;
      ! z is +ve down but denotes radii from the Earth's center
      ! moreover, all elements assume that dual grid is shifted
      ! in the positive direction - careful with indices
      allocate(x(nx),y(ny+1),z(nz+1),STAT=istat)
      x = grid%x
      y = grid%y
      z = grid%z


      ! surface areas perpendicular to the dual edges
      ! in the longitudinal direction are all the same
      ! for a chosen latitude and depth;
      ! j=1 is North pole, j=ny is South pole, both are defined;
      ! this also takes care of zero longitude
      do k=1,nz
        do j=1,ny
          call area_sijk(j,k,y,z,sijk)
          do i=1,nx+1
            S_F%x(i,j,k) = sijk
          end do
        end do
      end do

      ! surface areas perpendicular to the dual edges
      ! in the latitudinal direction depend on longitude
      ! for non-uniform grids
      do k=1,nz
        do i=1,nx
          ! North pole
          j=1
          S_F%y(i,j,k) = 0.0d0
          do j=2,ny
            call area_sjki(i,j,k,x,y,z,sjki)
            S_F%y(i,j,k) = sjki
          end do
          ! South pole
          j=ny+1
          S_F%y(i,j,k) = 0.0d0
        end do
      end do

      ! horizontal surface areas are defined at the poles
      ! and on upper and lower boundaries too; the usual
      ! formula may be used at the poles
      do j=1,ny
        do i=1,nx
          do k=1,nz+1
            call area_skij(i,j,k,x,y,z,skij)
            S_F%z(i,j,k) = skij
          end do
        end do
      end do

      deallocate(x,y,z)

      call validate_rvector(S_F,.true.)

      return
      end subroutine FaceArea

end module GridCalc
