! *****************************************************************************
! initializes and does basic calculations for the grid. Computes basic derived
! grid parameters that are used repeatedly by other routines, such as the unit
! lengths and volumes, and other grid details
module GridCalc

  use sg_vector
  use sg_scalar
  use elements
  implicit none

  public      :: EdgeVolume, CornerVolume

Contains

  ! *************************************************************************
  ! * EdgeVolume creates volume elements centered around the edges of
  ! * the grid, and stores them as real vectors with gridType=EDGE.
  ! *
  ! * A diagonal matrix multiplication of the edge volume with the difference
  ! * equations enables us to make a symmetrical matrix. Remember,
  ! * the electrical fields are defined on the center of the edges, therefore,
  ! * the edge volume is centered about the electrical field measurement.

  subroutine EdgeVolume(inGr, eV)

    implicit none
    type (grid_t), intent(in)           :: inGr     ! input model
    type (rvector), intent(inout)         :: eV       ! edge volume
    integer                               :: ix, iy, iz

    ! Checks whether the size is the same
    if ((inGr%nx == eV%nx).and.&
         (inGr%ny == eV%ny).and.&
         (inGr%nz == eV%nz)) then

       if (eV%gridType == EDGE) then

          ! edge volume are made for all the edges
          ! for x-components
          do ix = 1,inGr%nx
             do iy = 1,inGr%ny+1
                do iz = 1,inGr%nz+1

                   ! eV%x values are centered within dx.
                   eV%x(ix, iy, iz) = inGr%dx(ix)*inGr%delY(iy)*inGr%delZ(iz)

                enddo
             enddo
          enddo

          ! edge volume are made for all the edges
          ! for y-components
          do ix = 1,inGr%nx+1
             do iy = 1,inGr%ny
                do iz = 1,inGr%nz+1

                   ! eV%y values are centered within dy.
                   eV%y(ix, iy, iz) = inGr%delX(ix)*inGr%dy(iy)*inGr%delZ(iz)

                enddo
             enddo
          enddo

          ! edge volume are made for all the edges
          ! for z-components
          do ix = 1,inGr%nx+1
             do iy = 1,inGr%ny+1
                do iz = 1,inGr%nz

                   ! eV%z values are centered within dz.
                   eV%z(ix, iy, iz) = inGr%delX(ix)*inGr%delY(iy)*inGr%dz(iz)

                enddo
             enddo
          enddo

       else
          write (0, *) 'EdgeVolume: not compatible usage for existing data types'
       end if


    else
       write(0, *) 'Error-grid size and edge volume are not the same size'
    endif

  end subroutine EdgeVolume  ! EdgeVolume

  ! *************************************************************************
  ! * CornerVolume creates volume elements centered around the corners of
  ! * the grid, and stores them as real scalars with gridType=CORNER.

  subroutine CornerVolume(inGr, cV)

    type (grid_t), intent(in)        :: inGr  ! input grid
    type (rscalar), intent(inout)      :: cV    ! center volume as output
    integer                            :: ix, iy, iz

    ! Checks whether the size is the same
    if ((inGr%nx == cV%nx).and.&
         (inGr%ny == cV%ny).and.&
         (inGr%nz == cV%nz)) then

       if (cV%gridType == CORNER) then

          ! center volume is only using the internal corner nodes
          do ix = 2, inGr%nx
             do iy = 2, inGr%ny
                do iz = 2, inGr%nz

                   ! note that we are multiplying
                   ! using the distances with corner of a cell as a center
                   cV%v(ix, iy, iz) = inGr%delX(ix)*inGr%delY(iy)*inGr%delZ(iz)

                enddo
             enddo
          enddo

       else
          write (0, *) 'CornerVolume: not compatible usage for existing data types'
       end if

    else
       write(0, *) 'Error-grid size and center volume are not the same size'
    endif

  end subroutine CornerVolume

end module GridCalc
