! *****************************************************************************
! initializes and does basic calculations for the grid. Computes basic derived
! grid parameters that are used repeatedly by other routines, such as the unit
! lengths and volumes, and other grid details
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

  subroutine EdgeVolume(grid, V_E, l_E, S_E)

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
        call EdgeLength(grid,length)
        call EdgeArea(grid,area)
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
      call diagMult_rvector(length,area,V_F)

      call deall_rvector(length)
      call deall_rvector(area)

  end subroutine FaceVolume  ! FaceVolume

  ! *************************************************************************
  ! * NodeVolume creates volume elements centered around the corners (nodes)
  ! * of the grid, and stores them as real scalars with gridType=CORNER.
  ! * IMPORTANT: Boundary volume elements are also defined.

  subroutine NodeVolume(grid, V_N)

      type (grid_t), intent(in)          :: grid  ! input grid
      type (rscalar), intent(inout)      :: V_N   ! node/corner volume as output
      ! local variables
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

      ! center volume is only using the internal corner nodes
      do i = 1, grid%nx+1
         do j = 1, grid%ny+1
            do k = 1, grid%nz+1

               ! note that we are multiplying
               ! using the distances with corner of a cell as a center
               V_N%v(i, j, k) = grid%delX(i)*grid%delY(j)*grid%delZ(k)

            enddo
         enddo
      enddo

  end subroutine NodeVolume

  ! *************************************************************************
  ! * CellVolumes may prove useful in mappings from cells to faces and back.

  subroutine CellVolume(grid, V_C)

      implicit none
      type (grid_t), intent(in)          :: grid     ! input model
      type (rscalar), intent(inout)      :: V_C       ! cell volume
      ! local variables
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

      ! cell volume is only using the internal corner nodes
      do i = 1, grid%nx
         do j = 1, grid%ny
            do k = 1, grid%nz

               V_C%v(i, j, k) = grid%dx(i)*grid%dy(j)*grid%dz(k)

            enddo
         enddo
      enddo

  end subroutine CellVolume  ! CellVolume

  ! *************************************************************************
  ! * EdgeLength creates line elements defined on edges of the primary grid.
  ! * Edge length elements are defined on interior and boundary edges.

  subroutine EdgeLength(grid,l_E)

      type(grid_t), intent(in)      :: grid
      type(rvector), intent(inout)  :: l_E
      ! local variables
      integer                   :: ix,iy,iz

      call create_rvector(grid, l_E, EDGE)

      ! x-component edge length elements
      do ix = 1,grid%nx
         l_E%x(ix, :, :) = grid%dx(ix)
      enddo

      ! y-component edge length elements
      do iy = 1,grid%ny
         l_E%y(:, iy, :) = grid%dy(iy)
      enddo

      ! z-component edge length elements
      do iz = 1,grid%nz
         l_E%z(:, :, iz) = grid%dz(iz)
      enddo

  end subroutine EdgeLength

  ! *************************************************************************
  ! * FaceLength creates line elements defined on faces of the primary grid.
  ! * These line elements are perpendicular to a face with center coinciding
  ! * with the face center. They correspond to the edges of the dual grid.
  ! *
  ! * Face length elements are correctly defined on the boundary faces as the
  ! * perpendicular distances from the boundary faces to primary cell centers.
  ! * For a cartesian grid, these are 1/2 of the corresponding edge lengths.
  ! * These have been dealt with correctly as delX, delY, delZ in the grid_setup;
  ! * however, we intent to replace these completely with this subroutine as
  ! * the more general version. Therefore, recomputing them here.

  subroutine FaceLength(grid,l_F)

      type(grid_t), intent(in)      :: grid
      type(rvector), intent(inout)  :: l_F
      ! local variables
      integer                   :: ix,iy,iz

      call create_rvector(grid, l_F, FACE)

      ! x-component face length elements
      ! are the dual grid edge lengths
      l_F%x(1,:,:) = grid%dx(1)
      do ix = 2,grid%nx
         l_F%x(ix, :, :) = grid%dx(ix-1) + grid%dx(ix)
      enddo
      l_F%x(grid%nx+1,:,:) = grid%dx(grid%nx)
      l_F%x = l_F%x/2.0d0

      ! y-component face length elements
      ! are the dual grid edge lengths
      l_F%y(:,1,:) = grid%dy(1)
      do iy = 2,grid%ny
         l_F%y(:, iy, :) = grid%dy(iy-1) + grid%dy(iy)
      enddo
      l_F%y(:,grid%ny+1,:) = grid%dy(grid%ny)
      l_F%y = l_F%y/2.0d0

      ! z-component face length elements
      ! are the dual grid edge lengths
      l_F%z(:,:,1) = grid%dz(1)
      do iz = 2,grid%nz
         l_F%z(:, :, iz) = grid%dz(iz-1) + grid%dz(iz)
      enddo
      l_F%z(:,:,grid%nz+1) = grid%dz(grid%nz)
      l_F%z = l_F%z/2.0d0

  end subroutine FaceLength

  ! ***************************************************************************
  ! * EdgeArea: surface area elements perpendicular to the edges of the primary
  ! * grid, with indices matching the primary grid edges. These correspond to
  ! * faces of the dual grid.
  ! *
  ! * Edge surface area elements are correctly defined on the boundary edge as
  ! * the interior area element perpendicular to the edge.
  ! * For a cartesian grid, these use 1/2's of the corresponding dual edge lengths.
  ! * These have been dealt with correctly as delX, delY, delZ in the grid_setup;
  ! * however, we intent to replace these completely with this subroutine as
  ! * the more general version. This might be recoded in the future to first
  ! * compute the dual grid edge lengths.

  subroutine EdgeArea(grid,S_E)

      type(grid_t), intent(in)      :: grid
      type(rvector), intent(inout)  :: S_E
      ! local variables
      integer                   :: ix,iy,iz

      call create_rvector(grid, S_E, EDGE)

      ! edge areas are made for all the edges
      ! for x-components
      do ix = 1,grid%nx
         do iy = 1,grid%ny+1
            do iz = 1,grid%nz+1

               ! S_E%x values are centered within dx.
               S_E%x(ix, iy, iz) = grid%delY(iy)*grid%delZ(iz)

            enddo
         enddo
      enddo

      ! edge areas are made for all the edges
      ! for y-components
      do ix = 1,grid%nx+1
         do iy = 1,grid%ny
            do iz = 1,grid%nz+1

               ! S_E%y values are centered within dy.
               S_E%y(ix, iy, iz) = grid%delX(ix)*grid%delZ(iz)

            enddo
         enddo
      enddo

      ! edge areas are made for all the edges
      ! for z-components
      do ix = 1,grid%nx+1
         do iy = 1,grid%ny+1
            do iz = 1,grid%nz

               ! S_E%z values are centered within dz.
               S_E%z(ix, iy, iz) = grid%delX(ix)*grid%delY(iy)

            enddo
         enddo
      enddo

      end subroutine EdgeArea

  ! *************************************************************************
  ! * FaceArea computes surface area elements on faces of the primary grid.
  ! * Face surface area elements are defined on interior and boundary faces.

  subroutine FaceArea(grid,S_F)

      type(grid_t), intent(in)      :: grid
      type(rvector), intent(inout)  :: S_F
      ! local variables
      integer                   :: ix,iy,iz

      call create_rvector(grid, S_F, FACE)

      ! face areas are made for all the faces
      ! for x-components
      do ix = 1,grid%nx+1
         do iy = 1,grid%ny
            do iz = 1,grid%nz

               S_F%x(ix, iy, iz) = grid%dy(iy)*grid%dz(iz)

            enddo
         enddo
      enddo

      ! face areas are made for all the faces
      ! for y-components
      do ix = 1,grid%nx
         do iy = 1,grid%ny+1
            do iz = 1,grid%nz

               S_F%y(ix, iy, iz) = grid%dx(ix)*grid%dz(iz)

            enddo
         enddo
      enddo

      ! face areas are made for all the faces
      ! for z-components
      do ix = 1,grid%nx
         do iy = 1,grid%ny
            do iz = 1,grid%nz+1

               S_F%z(ix, iy, iz) = grid%dx(ix)*grid%dy(iy)

            enddo
         enddo
      enddo

      end subroutine FaceArea

  ! *************************************************************************
  ! * UNWEIGHTED MAPPING OPERATORS THAT LEAVE ALL BOUNDARY EDGES SET TO ZERO
  ! * (NOTE: most likely, this will also work for multigrid: we would use
  ! * separate subroutines to map from egdes to multigrid edges and back)
  ! *************************************************************************

  ! *************************************************************************
  ! * Cell2Edge will be used by forward model mappings
  ! * might need to revisit boundary edges

  subroutine Cell2Edge(grid,C,E)

      type(grid_t), intent(in)      :: grid
      type(rscalar), intent(in)     :: C
      type(rvector), intent(out)    :: E
      ! local variables
      integer                   :: ix,iy,iz

      call create_rvector(grid, E, EDGE)

      ! for x-components
      do ix = 1,grid%nx
         do iy = 2,grid%ny
            ! special case of upper boundary
            iz = 1
            E%x(ix, iy, iz) = (C%v(ix, iy-1, iz) + C%v(ix, iy, iz))/2.0d0
            ! inside the Earth
            do iz = 2,grid%nz

               E%x(ix, iy, iz) = (C%v(ix, iy-1, iz-1) + C%v(ix, iy, iz-1) + &
                                  C%v(ix, iy-1, iz) + C%v(ix, iy, iz))/4.0d0

            enddo
            ! special case of lower boundary
            iz = grid%nz+1
            E%x(ix, iy, iz) = (C%v(ix, iy-1, iz-1) + C%v(ix, iy, iz-1))/2.0d0
         enddo
      enddo

      ! for y-components
      do ix = 2,grid%nx
         do iy = 1,grid%ny
            ! special case of upper boundary
            iz = 1
            E%y(ix, iy, iz) = (C%v(ix-1, iy, iz) + C%v(ix, iy, iz))/2.0d0
            ! inside the Earth
            do iz = 2,grid%nz

               E%y(ix, iy, iz) = (C%v(ix-1, iy, iz-1) + C%v(ix, iy, iz-1) + &
                                  C%v(ix-1, iy, iz) + C%v(ix, iy, iz))/4.0d0

            enddo
            ! special case of lower boundary
            iz = grid%nz+1
            E%y(ix, iy, iz) = (C%v(ix-1, iy, iz-1) + C%v(ix, iy, iz-1))/2.0d0
         enddo
      enddo

      ! for z-components
      do ix = 2,grid%nx
         do iy = 2,grid%ny
            do iz = 1,grid%nz

               E%z(ix, iy, iz) = (C%v(ix-1, iy-1, iz) + C%v(ix-1, iy, iz) + &
                                  C%v(ix, iy-1, iz) + C%v(ix, iy, iz))/4.0d0

            enddo
         enddo
      enddo

  end subroutine Cell2Edge

  ! *************************************************************************
  ! * Edge2Cell will be used by adjoint model mappings

  subroutine Edge2Cell(grid,E,C)

      type(grid_t), intent(in)      :: grid
      type(rvector), intent(in)     :: E
      type(rscalar), intent(out)    :: C
      ! local variables
      integer                   :: ix,iy,iz

      call create_rscalar(grid, C, CENTER)

      ! for x-components
      do ix = 1,grid%nx
         do iy = 2,grid%ny
            ! special case of upper boundary
            iz = 1
            C%v(ix, iy-1, iz) = C%v(ix, iy-1, iz) + E%x(ix, iy, iz)/2.0d0
            C%v(ix, iy, iz) = C%v(ix, iy, iz) + E%x(ix, iy, iz)/2.0d0
            ! inside the Earth
            do iz = 2,grid%nz

               C%v(ix, iy-1, iz-1) = C%v(ix, iy-1, iz-1) + E%x(ix, iy, iz)/4.0d0
               C%v(ix, iy, iz-1) = C%v(ix, iy, iz-1) + E%x(ix, iy, iz)/4.0d0
               C%v(ix, iy-1, iz) = C%v(ix, iy-1, iz) + E%x(ix, iy, iz)/4.0d0
               C%v(ix, iy, iz) = C%v(ix, iy, iz) + E%x(ix, iy, iz)/4.0d0

            enddo
            ! special case of lower boundary
            iz = grid%nz+1
            C%v(ix, iy-1, iz-1) = C%v(ix, iy-1, iz-1) + E%x(ix, iy, iz)/2.0d0
            C%v(ix, iy, iz-1) = C%v(ix, iy, iz-1) + E%x(ix, iy, iz)/2.0d0
         enddo
      enddo

      ! for y-components
      do ix = 2,grid%nx
         do iy = 1,grid%ny
            ! special case of upper boundary
            iz = 1
            C%v(ix-1, iy, iz) = C%v(ix-1, iy, iz) + E%y(ix, iy, iz)/2.0d0
            C%v(ix, iy, iz) = C%v(ix, iy, iz) + E%y(ix, iy, iz)/2.0d0
            ! inside the Earth
            do iz = 2,grid%nz

               C%v(ix-1, iy, iz-1) = C%v(ix-1, iy, iz-1) + E%y(ix, iy, iz)/4.0d0
               C%v(ix, iy, iz-1) = C%v(ix, iy, iz-1) + E%y(ix, iy, iz)/4.0d0
               C%v(ix-1, iy, iz) = C%v(ix-1, iy, iz) + E%y(ix, iy, iz)/4.0d0
               C%v(ix, iy, iz) = C%v(ix, iy, iz) + E%y(ix, iy, iz)/4.0d0

            enddo
            ! special case of lower boundary
            iz = grid%nz+1
            C%v(ix-1, iy, iz-1) = C%v(ix-1, iy, iz-1) + E%y(ix, iy, iz)/2.0d0
            C%v(ix, iy, iz-1) = C%v(ix, iy, iz-1) + E%y(ix, iy, iz)/2.0d0
         enddo
      enddo

      ! for z-components
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

  end subroutine Edge2Cell

end module GridCalc
