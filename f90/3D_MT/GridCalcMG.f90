! *****************************************************************************
! initializes and does basic calculations for the grid. Computes basic derived
! grid parameters that are used repeatedly by other routines, such as the unit
! lengths and volumes, and other grid details
module GridCalcMG

  use sg_vector_mg
  use sg_scalar_mg
  use gridcalc
  implicit none
  !private

  public      :: EdgeVolume, NodeVolume
  !public      :: EdgeVolume, FaceVolume, NodeVolume, CellVolume
  !public      :: EdgeLength, FaceLength
  !public      :: EdgeArea,   FaceArea

  ! Interfaces allow for function overloading (polymorphism)
  ! for derived (extended) grid types
  interface EdgeVolume
    module procedure EdgeVolumeFDMG
  end interface

!  interface FaceVolume
!    module procedure FaceVolumeFDMG
!  end interface

  interface NodeVolume
    module procedure NodeVolumeFDMG
  end interface

!  interface CellVolume
!    module procedure CellVolumeFDMG
!  end interface
!
!  interface EdgeLength
!    module procedure EdgeLengthFDMG
!  end interface
!
!  interface FaceLength
!    module procedure FaceLengthFDMG
!  end interface
!
!  interface EdgeArea
!    module procedure EdgeAreaFDMG
!  end interface
!
!  interface FaceArea
!    module procedure FaceAreaFDMG
!  end interface

Contains

  ! *************************************************************************
  ! * EdgeVolume creates volume elements centered around the edges of
  ! * the subgrids, and stores them as real vectors with gridType=EDGE.
  ! *
  ! * A diagonal matrix multiplication of the edge volume with the difference
  ! * equations enables us to make a symmetrical matrix. Remember,
  ! * the electrical fields are defined on the center of the edges, therefore,
  ! * the edge volume is centered about the electrical field measurement.
  subroutine EdgeVolumeFDMG(mgrid, V_E)

  implicit none
    type (grid_t), intent(in)  :: mgrid            ! input MULTIGRID
    type (rvector_mg), intent(inout)  :: V_E       ! edge volume elements
    ! local variables
    integer  :: ix, iy, iz, izCum, imgrid
    integer  :: nx, ny, nz, nzCum

    ! Checks whether the size is the same
    if (mgrid%mgridSize == V_E%mgridSize) then
       if (V_E%gridType == EDGE) then

         do imgrid = 1, mgrid%mgridSize  ! global loop on subgrid
            call EdgeVolumeFD(mgrid%gridArray(imgrid),V_E%rvArray(imgrid))
         enddo ! global loop on subgrid

       else
         write (0, *) 'EdgeVolume: not compatible usage for existing data types'
       end if


    else
       write(0, *) 'Error-grid size and edge volume are not the same size (mgridSzie)'
    endif

  end subroutine EdgeVolumeFDMG  ! EdgeVolume Multigrid case

!  subroutine EdgeVolumeFDMG(inGr, eV)
!
!  implicit none
!    type (grid_t), intent(in)  :: inGr            ! input MULTIGRID
!    type (rvector_mg), intent(inout)  :: eV       ! edge volume
!    ! local variables
!    integer  :: ix, iy, iz, izCum, imgrid
!    integer  :: nx, ny, nz, nzCum
!
!    ! Checks whether the size is the same
!    if (inGr%mgridSize == eV%mgridSize) then
!       if (eV%gridType == EDGE) then
!
!         do imgrid = 1, inGr%mgridSize  ! global loop on subgrid
!
!            call EdgeVolumeFD(inGr%gridArray(imgrid),eV%rvArray(imgrid))
!
!           nx = inGr%gridArray(imgrid)%nx  ! nx current subgrid
!           ny = inGr%gridArray(imgrid)%ny  ! ny current subgrid
!           nz = inGr%gridArray(imgrid)%nz  ! nz=nzEarth+nzAir current
!
!           if ((nx == eV%rvArray(imgrid)%nx).and. &  ! compares nx ny nz (if...)
!              (ny == eV%rvArray(imgrid)%ny).and. &
!              (nz == eV%rvArray(imgrid)%nz)) then
!
!              ! edge volume are made for all the edges
!              ! for x-components
!              do ix = 1, nx
!                do iy = 1, ny+1
!                  do iz = 1, nz+1
!                    ! eV%x values are centered within dx.
!                    eV%rvArray(imgrid)%x(ix, iy, iz) = inGr%gridArray(imgrid)%dx(ix)*inGr%gridArray(imgrid)%delY(iy) &
!                                                                                                 *inGr%gridArray(imgrid)%delZ(iz)
!                 enddo
!              enddo
!            enddo
!
!            ! edge volume are made for all the edges
!            ! for y-components
!            do ix = 1, nx+1
!               do iy = 1, ny
!                  do iz = 1, nz+1
!                     ! eV%y values are centered within dy.
!                     eV%rvArray(imgrid)%y(ix, iy, iz) = inGr%gridArray(imgrid)%delX(ix)*inGr%gridArray(imgrid)%dy(iy) &
!                                                                                                  *inGr%gridArray(imgrid)%delZ(iz)
!                  enddo
!               enddo
!            enddo
!
!            ! edge volume are made for all the edges
!            ! for z-components
!            do ix = 1, nx+1
!               do iy = 1,ny+1
!                  do iz = 1, nz
!                     ! eV%z values are centered within dz.
!                     eV%rvArray(imgrid)%z(ix, iy, iz) = inGr%gridArray(imgrid)%delX(ix)*inGr%gridArray(imgrid)%delY(iy) &
!                                                                                               *inGr%gridArray(imgrid)%dz(iz)
!                  enddo
!               enddo
!            enddo
!
!
!           else ! compares nx ny nz (if...)
!            write (0, *) 'Error-grid size and edge volume are not the same size (nx not= ny ...)'
!           end if  ! compares nx ny nz (if...)
!
!         enddo ! global loop on subgrid
!
!       else
!         write (0, *) 'EdgeVolume: not compatible usage for existing data types'
!       end if
!
!
!    else
!       write(0, *) 'Error-grid size and edge volume are not the same size (mgridSzie)'
!    endif
!
!  end subroutine EdgeVolumeFDMG  ! EdgeVolume Multigrid case

  ! ******************************************************************************************
  ! * NodeVolume creates volume elements centered around the corners of
  ! * the grid, and stores them as real scalars with gridType=CORNER.

  subroutine NodeVolumeFDMG(mgrid, V_N)

  implicit none
    type (grid_t), intent(in)  :: mgrid                 ! input multigrid
    type (rscalar_mg), intent(inout)  :: V_N            ! node volume elements
    integer  :: ix, iy, iz
    integer  :: imgrid
    integer  :: nx, ny, nz


    if (V_N%gridType == CORNER) then

      if (V_N%mgridSize == mgrid%mgridSize) then

        do imgrid = 1, mgrid%mgridSize
            call NodeVolumeFD(mgrid%gridArray(imgrid),V_N%rscArray(imgrid))
        enddo

      else
       print *, 'NodeVolume: error mgridSize'
      end if

    else
       print *, 'NodeVolumeMG: not compatible usage for existing data types'
    end if

  end subroutine NodeVolumeFDMG


!  subroutine NodeVolumeFDMG(mgrid, cV)
!
!  implicit none
!    type (grid_t), intent(in)  :: mgrid                 ! input multigrid
!    type (rscalar_mg), intent(inout)  :: cV             ! center volume as output
!    integer  :: ix, iy, iz
!    integer  :: imgrid
!    integer  :: nx, ny, nz
!
!
!    if (cV%gridType == CORNER) then
!
!      if (cV%mgridSize == mgrid%mgridSize) then
!
!        do imgrid = 1, mgrid%mgridSize
!          nx = mgrid%gridArray(imgrid)%nx
!          ny = mgrid%gridArray(imgrid)%ny
!          nz = mgrid%gridArray(imgrid)%nz
!
!          !Checks whether the size is the same
!          if ((nx == cV%rscArray(imgrid)%nx).and.&
!              (ny == cV%rscArray(imgrid)%ny).and.&
!              (nz == cV%rscArray(imgrid)%nz)) then
!            ! center volume is only using the internal corner nodes
!
!            if(imgrid == 1) then
!              do ix = 2, nx
!                 do iy = 2, ny
!                   do iz = 2, nz+1
!
!                     ! note that we are multiplying
!                     ! using the distances with corner of a cell as a center
!                     cV%rscArray(imgrid)%v(ix, iy, iz) = mgrid%gridArray(imgrid)%delX(ix)*mgrid%gridArray(imgrid)%delY(iy) &
!                                                       *mgrid%gridArray(imgrid)%delZ(iz)
!                    enddo
!                 enddo
!              enddo
!            else if(imgrid == mgrid%mgridSize)then
!              do ix = 2, nx
!                 do iy = 2, ny
!                     do iz = 1, nz
!
!                     ! note that we are multiplying
!                     ! using the distances with corner of a cell as a center
!                     cV%rscArray(imgrid)%v(ix, iy, iz) = mgrid%gridArray(imgrid)%delX(ix)*mgrid%gridArray(imgrid)%delY(iy) &
!                                                       *mgrid%gridArray(imgrid)%delZ(iz)
!                    enddo
!                 enddo
!              enddo
!            else
!              do ix = 2, nx
!                 do iy = 2, ny
!                    do iz = 1, nz+1
!
!                     ! note that we are multiplying
!                     ! using the distances with corner of a cell as a center
!                     cV%rscArray(imgrid)%v(ix, iy, iz) = mgrid%gridArray(imgrid)%delX(ix)*mgrid%gridArray(imgrid)%delY(iy) &
!                                                       *mgrid%gridArray(imgrid)%delZ(iz)
!                    enddo
!                 enddo
!              enddo
!            endif
!          else
!            print *, 'Error-grid size and center volume are not the same size in NodeVolumeMG'
!          endif
!
!        enddo
!
!      else
!       print *, 'NodeVolume: error mgridSize'
!      end if
!
!    else
!       print *, 'NodeVolumeMG: not compatible usage for existing data types'
!    end if
!
!  end subroutine NodeVolumeFDMG
 ! *************************************************************************
end module GridCalcMG
