! *****************************************************************************
! initializes and does basic calculations for the grid. Computes basic derived
! grid parameters that are used repeatedly by other routines, such as the unit
! lengths and volumes, and other grid details
module GridCalc

  use sg_vector_mg
  use sg_scalar_mg
  implicit none

  public  :: EdgeVolume, CornerVolume

Contains

  ! *************************************************************************
  ! * EdgeVolume creates volume elements centered around the edges of
  ! * the subgrids, and stores them as real vectors with gridType=EDGE.
  ! *
  ! * A diagonal matrix multiplication of the edge volume with the difference
  ! * equations enables us to make a symmetrical matrix. Remember,
  ! * the electrical fields are defined on the center of the edges, therefore,
  ! * the edge volume is centered about the electrical field measurement.

  subroutine EdgeVolume(inGr, eV)

  implicit none
    type (grid_t), intent(in)  :: inGr            ! input MULTIGRID
    type (rvector_mg), intent(inout)  :: eV       ! edge volume
    ! local variables
    integer  :: ix, iy, iz, izCum, imgrid
    integer  :: nx, ny, nz, nzCum

    ! Checks whether the size is the same
    if (inGr%mgridSize == eV%mgridSize) then
       if (eV%gridType == EDGE) then

         do imgrid = 1, inGr%mgridSize  ! global loop on subgrid

           nx = inGr%gridArray(imgrid)%nx  ! nx current subgrid
           ny = inGr%gridArray(imgrid)%ny  ! ny current subgrid
           nz = inGr%gridArray(imgrid)%nz  ! nz=nzEarth+nzAir current

           if ((nx == eV%rvArray(imgrid)%nx).and. &  ! compares nx ny nz (if...)
              (ny == eV%rvArray(imgrid)%ny).and. &
              (nz == eV%rvArray(imgrid)%nz)) then

              ! edge volume are made for all the edges
              ! for x-components
              do ix = 1, nx
                do iy = 1, ny+1
                  do iz = 1, nz+1
                    ! eV%x values are centered within dx.
                    eV%rvArray(imgrid)%x(ix, iy, iz) = inGr%gridArray(imgrid)%dx(ix)*inGr%gridArray(imgrid)%delY(iy) &
                                                                                                 *inGr%gridArray(imgrid)%delZ(iz)
                 enddo
              enddo
            enddo

            ! edge volume are made for all the edges
            ! for y-components
            do ix = 1, nx+1
               do iy = 1, ny
                  do iz = 1, nz+1
                     ! eV%y values are centered within dy.
                     eV%rvArray(imgrid)%y(ix, iy, iz) = inGr%gridArray(imgrid)%delX(ix)*inGr%gridArray(imgrid)%dy(iy) &
                                                                                                  *inGr%gridArray(imgrid)%delZ(iz)
                  enddo
               enddo
            enddo

            ! edge volume are made for all the edges
            ! for z-components
            do ix = 1, nx+1
               do iy = 1,ny+1
                  do iz = 1, nz
                     ! eV%z values are centered within dz.
                     eV%rvArray(imgrid)%z(ix, iy, iz) = inGr%gridArray(imgrid)%delX(ix)*inGr%gridArray(imgrid)%delY(iy) &
                                                                                               *inGr%gridArray(imgrid)%dz(iz)
                  enddo
               enddo
            enddo


           else ! compares nx ny nz (if...)
            write (0, *) 'Error-grid size and edge volume are not the same size (nx not= ny ...)'
           end if  ! compares nx ny nz (if...)

         enddo ! global loop on subgrid

       else
         write (0, *) 'EdgeVolume: not compatible usage for existing data types'
       end if


    else
       write(0, *) 'Error-grid size and edge volume are not the same size (mgridSzie)'
    endif

  end subroutine EdgeVolume  ! EdgeVolume Multigrid case

  ! ******************************************************************************************
  ! * CornerVolume creates volume elements centered around the corners of
  ! * the grid, and stores them as real scalars with gridType=CORNER.

  subroutine CornerVolume(mgrid, cV)

  implicit none
    type (grid_t), intent(in)  :: mgrid                 ! input multigrid
    type (rscalar_mg), intent(inout)  :: cV             ! center volume as output
    integer  :: ix, iy, iz
    integer  :: imgrid
    integer  :: nx, ny, nz


    if (cV%gridType == CORNER) then

      if (cV%mgridSize == mgrid%mgridSize) then

        do imgrid = 1, mgrid%mgridSize
          nx = mgrid%gridArray(imgrid)%nx
          ny = mgrid%gridArray(imgrid)%ny
          nz = mgrid%gridArray(imgrid)%nz

          !Checks whether the size is the same
          if ((nx == cV%rscArray(imgrid)%nx).and.&
              (ny == cV%rscArray(imgrid)%ny).and.&
              (nz == cV%rscArray(imgrid)%nz)) then
            ! center volume is only using the internal corner nodes

            if(imgrid == 1) then
              do ix = 2, nx
                 do iy = 2, ny
                   do iz = 2, nz+1

                     ! note that we are multiplying
                     ! using the distances with corner of a cell as a center
                     cV%rscArray(imgrid)%v(ix, iy, iz) = mgrid%gridArray(imgrid)%delX(ix)*mgrid%gridArray(imgrid)%delY(iy) &
                                                       *mgrid%gridArray(imgrid)%delZ(iz)
                    enddo
                 enddo
              enddo
            else if(imgrid == mgrid%mgridSize)then
              do ix = 2, nx
                 do iy = 2, ny
                     do iz = 1, nz

                     ! note that we are multiplying
                     ! using the distances with corner of a cell as a center
                     cV%rscArray(imgrid)%v(ix, iy, iz) = mgrid%gridArray(imgrid)%delX(ix)*mgrid%gridArray(imgrid)%delY(iy) &
                                                       *mgrid%gridArray(imgrid)%delZ(iz)
                    enddo
                 enddo
              enddo
            else
              do ix = 2, nx
                 do iy = 2, ny
                    do iz = 1, nz+1

                     ! note that we are multiplying
                     ! using the distances with corner of a cell as a center
                     cV%rscArray(imgrid)%v(ix, iy, iz) = mgrid%gridArray(imgrid)%delX(ix)*mgrid%gridArray(imgrid)%delY(iy) &
                                                       *mgrid%gridArray(imgrid)%delZ(iz)
                    enddo
                 enddo
              enddo
            endif
          else
            print *, 'Error-grid size and center volume are not the same size; CorverVolume;'
          endif

        enddo

      else
       print *, 'CornerVolume: error mgridSize'
      end if

    else
       print *, 'CornerVolumeMG: not compatible usage for existing data types'
    end if

  end subroutine CornerVolume
 ! *************************************************************************
end module GridCalc
