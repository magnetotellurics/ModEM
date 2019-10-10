! *****************************************************************************
module sg_spherical
  ! This module contains additional subroutines required to adapt the general
  ! SG_Basics scheme specifically to the Staggered Spherical grid. In these
  ! subroutines we assume that our model domain is enclosed between two spheres
  ! and has an upper and lower boundaries at k=1 and k=nz+1, respectively.
  ! The vectors defined in module sg_vector are reconsidered with this domain
  ! in mind. In particular, the North and South poles and the surrounding vector
  ! components, and the zero meridian has to be viewed as special cases and
  ! dealt with separately to make physical sense.

  use math_constants
  use sg_vector
  use sg_scalar
  use utilities
  implicit none

  ! Generic interfaces are done through subroutines
  ! creates edge/ face nodes
  INTERFACE validate
     module procedure validate_rvector
     module procedure validate_cvector
  END INTERFACE

  public		::   validate_rvector, validate_cvector

Contains
  ! VALIDATE GRID_edge/ face VECTORS
  ! * subroutine validate_rvector(igrid, E, gridType)
  ! * subroutine validate_cvector(igrid, E, gridType)


  ! ***************************************************************************
  ! validate_rvector corrects the vector if necessary to make physical sense
  ! when used in the context of a Spherical Staggered grid. Draw pictures to
  ! make sure everything makes sense.
  ! gridType is a character string to describe intended usage

  subroutine validate_rvector(E,verbose)

    implicit none
	logical, intent(in), optional		:: verbose
    type (rvector), intent(inout)        :: E
    type (rvector)						 :: inE,diffE
    character (len=80)					:: gridType

	logical								:: identical
    integer                            :: nx,ny,nz
	integer								:: i,j,k


    if(.not.E%allocated) then
       ! output an error and exit
		 write(0,*) 'Error: (validate_cvector) E not allocated'
		 stop
    end if

    ! Grid dimensions
    nx = E%nx
    ny = E%ny
    nz = E%nz

    ! gridType
    gridType = E%gridType

	! save the input vector
	inE = E

	! Assuming the vector has been correctly allocated, we now check
	! if it makes physical sense in the context of spherical grid.
	! If a value is undefined but included in the vector, we check that
	! it has been set to zero. This is used extensively in the program
	! to simplify special cases. If a value is repetitious, we check
	! that the relevant entries are identical.
    if (gridType == EDGE) then
	   ! North pole (j=1):
	   ! phi component: undefined
	   ! theta component: defined
	   ! r component: 1 edge
	   E%x(:,1,:) = R_ZERO
	   do k=1,nz
		do i=2,nx
		  E%z(i,1,k) = E%z(1,1,k)
		end do
	   end do

	   ! South pole (j=ny+1):
	   ! phi component: undefined
	   ! theta component: not included
	   ! r component: 1 edge
	   E%x(:,ny+1,:) = R_ZERO
	   do k=1,nz
		do i=2,nx
		  E%z(i,ny+1,k) = E%z(1,ny+1,k)
		end do
	   end do

	   ! Zero meridian (i=nx+1):
	   ! phi component: not included
	   ! theta component: duplicated
	   ! r component: duplicated
	   E%y(nx+1,:,:) = E%y(1,:,:)
	   E%z(nx+1,:,:) = E%z(1,:,:)


    else if (gridType == FACE) then
	   ! North pole (j=1):
	   ! phi component: defined
	   ! theta component: undefined
	   ! r component: defined
	   E%y(:,1,:) = R_ZERO

	   ! South pole (j=ny+1):
	   ! phi component: not included
	   ! theta component: undefined
	   ! r component: not included
	   E%y(:,ny+1,:) = R_ZERO

	   ! Zero meridian (i=nx+1):
	   ! phi component: duplicated
	   ! theta component: not included
	   ! r component: not included
	   E%x(nx+1,:,:) = E%x(1,:,:)

    else
       write (0, *) 'not a known tag'
    end if

	if (present(verbose)) then

	  ! Compare the resultant vector to the original. If not identical,
	  ! output a warning. diffE = E - inE
	  diffE = E
	  call linComb_rvector(ONE,E,MinusONE,inE,diffE)

	  if (gridType == EDGE) then
		! North pole
		identical = (sum(abs(diffE%x(:,1,:))) == R_ZERO)
		if (.not.identical) then
			write(0,*) 'Warning: (validate_cvector) North pole x-components corrected'
		end if
		identical = (sum(abs(diffE%z(:,1,:))) == R_ZERO)
		if (.not.identical) then
			write(0,*) 'Warning: (validate_cvector) North pole z-components corrected'
		end if
		! South pole
		identical = (sum(abs(diffE%x(:,ny+1,:))) == R_ZERO)
		if (.not.identical) then
			write(0,*) 'Warning: (validate_cvector) South pole x-components corrected'
		end if
		identical = (sum(abs(diffE%z(:,ny+1,:))) == R_ZERO)
		if (.not.identical) then
			write(0,*) 'Warning: (validate_cvector) South pole z-components corrected'
		end if
		! Zero meridian
		identical = (sum(abs(diffE%y(nx+1,:,:))) == R_ZERO)
		if (.not.identical) then
			write(0,*) 'Warning: (validate_cvector) Zero meridian y-components corrected'
		end if
		identical = (sum(abs(diffE%z(nx+1,:,:))) == R_ZERO)
		if (.not.identical) then
			write(0,*) 'Warning: (validate_cvector) Zero meridian z-components corrected'
		end if

	  else if (gridType == FACE) then
		! North pole
		identical = (sum(abs(diffE%y(:,1,:))) == R_ZERO)
		if (.not.identical) then
			write(0,*) 'Warning: (validate_cvector) North pole y-components corrected'
		end if
		! South pole
		identical = (sum(abs(diffE%y(:,ny+1,:))) == R_ZERO)
		if (.not.identical) then
			write(0,*) 'Warning: (validate_cvector) South pole y-components corrected'
		end if
		! Zero meridian
		identical = (sum(abs(diffE%x(nx+1,:,:))) == R_ZERO)
		if (.not.identical) then
			write(0,*) 'Warning: (validate_cvector) Zero meridian x-components corrected'
		end if
	  end if

	end if


  end subroutine validate_rvector  ! validate_rvector


  ! ***************************************************************************
  ! validate_cvector corrects the vector if necessary to make physical sense
  ! when used in the context of a Spherical Staggered grid. Draw pictures to
  ! make sure everything makes sense.
  ! gridType is a character string to describe intended usage
  subroutine validate_cvector(E,verbose)

    implicit none
	logical, intent(in), optional		  :: verbose
    type (cvector), intent(inout)         :: E
	type (cvector)						  :: inE,diffE
    character (len=80)					  :: gridType

	logical								:: identical
    integer                            :: nx,ny,nz
	integer								:: i,j,k


    if(.not.E%allocated) then
       ! output an error and exit
		 write(0,*) 'Error: (validate_cvector) E not allocated'
		 stop
    end if

    ! Grid dimensions
    nx = E%nx
    ny = E%ny
    nz = E%nz

    ! gridType
    gridType = E%gridType

	! save the input vector
	inE = E

	! Assuming the vector has been correctly allocated, we now check
	! if it makes physical sense in the context of spherical grid.
	! If a value is undefined but included in the vector, we check that
	! it has been set to zero. This is used extensively in the program
	! to simplify special cases. If a value is repetitious, we check
	! that the relevant entries are identical.
    if (gridType == EDGE) then
	   ! North pole (j=1):
	   ! phi component: undefined
	   ! theta component: defined
	   ! r component: 1 edge
	   E%x(:,1,:) = C_ZERO
	   do k=1,nz
		do i=2,nx
		  E%z(i,1,k) = E%z(1,1,k)
		end do
	   end do

	   ! South pole (j=ny+1):
	   ! phi component: undefined
	   ! theta component: not included
	   ! r component: 1 edge
	   E%x(:,ny+1,:) = C_ZERO
	   do k=1,nz
		do i=2,nx
		  E%z(i,ny+1,k) = E%z(1,ny+1,k)
		end do
	   end do

	   ! Zero meridian (i=nx+1):
	   ! phi component: not included
	   ! theta component: duplicated
	   ! r component: duplicated
	   E%y(nx+1,:,:) = E%y(1,:,:)
	   E%z(nx+1,:,:) = E%z(1,:,:)


    else if (gridType == FACE) then
	   ! North pole (j=1):
	   ! phi component: defined
	   ! theta component: undefined
	   ! r component: defined
	   E%y(:,1,:) = C_ZERO

	   ! South pole (j=ny+1):
	   ! phi component: not included
	   ! theta component: undefined
	   ! r component: not included
	   E%y(:,ny+1,:) = C_ZERO

	   ! Zero meridian (i=nx+1):
	   ! phi component: duplicated
	   ! theta component: not included
	   ! r component: not included
	   E%x(nx+1,:,:) = E%x(1,:,:)

    else
       write (0, *) 'not a known tag'
    end if


	if (present(verbose)) then

	  ! Compare the resultant vector to the original. If not identical,
	  ! output a warning. diffE = E - inE
	  diffE = E
	  call linComb_cvector(C_ONE,E,C_MinusONE,inE,diffE)

	  if (gridType == EDGE) then
		! North pole
		identical = (sum(abs(diffE%x(:,1,:))) == R_ZERO)
		if (.not.identical) then
			write(0,*) 'Warning: (validate_cvector) North pole x-components corrected'
		end if
		identical = (sum(abs(diffE%z(:,1,:))) == R_ZERO)
		if (.not.identical) then
			write(0,*) 'Warning: (validate_cvector) North pole z-components corrected'
		end if
		! South pole
		identical = (sum(abs(diffE%x(:,ny+1,:))) == R_ZERO)
		if (.not.identical) then
			write(0,*) 'Warning: (validate_cvector) South pole x-components corrected'
		end if
		identical = (sum(abs(diffE%z(:,ny+1,:))) == R_ZERO)
		if (.not.identical) then
			write(0,*) 'Warning: (validate_cvector) South pole z-components corrected'
		end if
		! Zero meridian
		identical = (sum(abs(diffE%y(nx+1,:,:))) == R_ZERO)
		if (.not.identical) then
			write(0,*) 'Warning: (validate_cvector) Zero meridian y-components corrected'
		end if
		identical = (sum(abs(diffE%z(nx+1,:,:))) == R_ZERO)
		if (.not.identical) then
			write(0,*) 'Warning: (validate_cvector) Zero meridian z-components corrected'
		end if

	  else if (gridType == FACE) then
		! North pole
		identical = (sum(abs(diffE%y(:,1,:))) == R_ZERO)
		if (.not.identical) then
			write(0,*) 'Warning: (validate_cvector) North pole y-components corrected'
		end if
		! South pole
		identical = (sum(abs(diffE%y(:,ny+1,:))) == R_ZERO)
		if (.not.identical) then
			write(0,*) 'Warning: (validate_cvector) South pole y-components corrected'
		end if
		! Zero meridian
		identical = (sum(abs(diffE%x(nx+1,:,:))) == R_ZERO)
		if (.not.identical) then
			write(0,*) 'Warning: (validate_cvector) Zero meridian x-components corrected'
		end if
	  end if

	end if

	call deall_cvector(inE)
	call deall_cvector(diffE)

  end subroutine validate_cvector  ! validate_cvector


  !****************************************************************************
  ! coordFlip_rvector makes an exact copy of derived data type
  ! rvector, except the coordinate systems are flipped;   NOTE: first argument is output
  ! used to switch between the MT coordinate system North-East-Down
  ! and the global code coordinate system Earth-South-Down.
  ! NOTE ALSO that this does NOTHING to the underlying grid.
  ! If we only use this to compute grid elements in GridCalc module,
  ! as we do now, we can save ourselves some work and not reverse the latitude order:
  ! it will still match the grid but they will be decreasing rather than increasing.
  ! NOT OK TO OVERWRITE THE INPUT WITH OUTPUT because the "spherical" vector is illegal -
  ! the grid doesn't match it's shape, so we cannot use copy statements.
  subroutine coordFlip_rvector(E, E1, newCoords, verbose)

    implicit none
    type (rvector), intent(in)       :: E1
    type (rvector), intent(inout)    :: E
    character(*), intent(in)         :: newCoords
    logical, intent(in), optional    :: verbose
    integer                      :: nx,ny,nz,status,k,nlat,nz1(3)

    ! check to see if RHS (E1) is active (allocated)
    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for copy_rvector'
    else

       if(E%allocated) then
          call deall_rvector(E)
       end if

       ! flip x & y coordinate dimensions
       nx = E1%ny
       ny = E1%nx
       nz = E1%nz

        ! allocate memory for x,y,z
        ! E%allocated will be true if all allocations succeed
        E%allocated = .true.
        if (E1%gridType == EDGE) then
           ! For spherical problem:
           ! 1) E%x(:,1,:) and E%x(:,ny+1,:) are undefined,
           ! 2) E%y(nx+1,:,:) and E%z(nx+1,:,:) are repetitios,
           ! 3) E%x(-1,0,ny+1,ny+2,:,:) will be needed for interpolation.
           allocate(E%x(nx,ny+1,nz+1), STAT=status)
           E%allocated = E%allocated .and. (status .EQ. 0)
           allocate(E%y(nx+1,ny,nz+1), STAT=status)
           E%allocated = E%allocated .and. (status .EQ. 0)
           allocate(E%z(nx+1,ny+1,nz), STAT=status)
           E%allocated = E%allocated .and. (status .EQ. 0)
           nz1(1)=E1%nz+1;nz1(2)=E1%nz+1;nz1(3)=E1%nz
        else if (E1%gridType == FACE) then
           ! For spherical problem:
           ! 1) E%y(:,1,:) and E%y(:,ny+1,:) are undefined,
           ! 2) E%x(nx+1,:,:) is repetitios and equals E%x(1,:,:),
           ! 3) E%z(:,:,1) and E%z(:,:,nz+1) are undefined.
           allocate(E%x(nx+1,ny,nz), STAT=status)
           E%allocated = E%allocated .and. (status .EQ. 0)
           allocate(E%y(nx,ny+1,nz), STAT=status)
           E%allocated = E%allocated .and. (status .EQ. 0)
           allocate(E%z(nx,ny,nz+1), STAT=status)
           E%allocated = E%allocated .and. (status .EQ. 0)
           nz1(1)=E1%nz;nz1(2)=E1%nz;nz1(3)=E1%nz+1
        else
           write (0, *) 'not a known tag in coordFlip_rvector'
        end if

        if (E%allocated) then
           E%x = R_ZERO
           E%y = R_ZERO
           E%z = R_ZERO
        else
            write (0, *) 'Warning: unable to allocate rvector - invalid grid supplied'
        end if

          if (trim(newCoords) == trim(Cartesian)) then

             ! switch horizontal components and the direction of latitude
             E%nx = E1%ny
             E%ny = E1%nx
             E%nz = E1%nz
             !nlat = size(E1%x,2)
             !E%x = E1%y !(:, nlat:1:-1, :) ! co-latitude South->latitude North
             !nlat = size(E1%y,2)
             !E%y = E1%x !(:, nlat:1:-1, :) ! longitude East
             !nlat = size(E1%z,2)
             !E%z = E1%z !(:, nlat:1:-1, :) ! Down
             do k=1,nz1(1)
                E%x(:,:,k) = transpose(E1%y(:,:,k))
             enddo
             do k=1,nz1(2)
                E%y(:,:,k) = transpose(E1%x(:,:,k))
             enddo
             do k=1,nz1(3)
                E%z(:,:,k) = transpose(E1%z(:,:,k))
             end do
             E%gridType = E1%gridType
             E%grid => E1%grid

          elseif (trim(newCoords) == trim(Spherical)) then

             ! switch horizontal components and the direction of latitude
             E%nx = E1%ny
             E%ny = E1%nx
             E%nz = E1%nz
             !nlat = size(E1%x,1)
             !E%x = E1%y !(nlat:1:-1, :, :) ! longitude East
             !nlat = size(E1%y,1)
             !E%y = E1%x !(nlat:1:-1, :, :) ! latitude North->co-latitude South
             !nlat = size(E1%z,1)
             !E%z = E1%z !(nlat:1:-1, :, :) ! Down
             do k=1,nz1(1)
                E%x(:,:,k) = transpose(E1%y(:,:,k))
             enddo
             do k=1,nz1(2)
                E%y(:,:,k) = transpose(E1%x(:,:,k))
             enddo
             do k=1,nz1(3)
                E%z(:,:,k) = transpose(E1%z(:,:,k))
             enddo
             E%gridType = E1%gridType
             E%grid => E1%grid

          else
             write (0, *) 'not compatible usage for coordFlip_rvector'
          end if

    end if

    ! if the input was a temporary function output, deallocate
    if (E1%temporary) then
        call deall_rvector(E1)
    end if

    if (present(verbose)) then
        write(*,*) '(coordFlip_rvector) successfully flipped the vector to ',trim(newCoords)
        write(*,*) 'without changing the underlying grid...'
        write(*,*) 'this renders the output vector illegal. Do not use except as a quick fix.'
    end if

  end subroutine coordFlip_rvector

  !****************************************************************************
  ! coordFlip_rscalar makes an exact copy of derived data type
  ! rscalar, except the coordinate systems are flipped;   NOTE: first argument is output
  ! used to switch between the MT coordinate system North-East-Down
  ! and the global code coordinate system Earth-South-Down.
  ! NOTE ALSO that this does NOTHING to the underlying grid.
  ! If we only use this to compute grid elements in GridCalc module,
  ! as we do now, we can save ourselves some work and not reverse the latitude order:
  ! it will still match the grid but they will be decreasing rather than increasing.
  ! NOT OK TO OVERWRITE THE INPUT WITH OUTPUT because the "spherical" vector is illegal -
  ! the grid doesn't match it's shape, so we cannot use copy statements.
  subroutine coordFlip_rscalar(E, E1, newCoords, verbose)

    implicit none
    type (rscalar), intent(in)       :: E1
    type (rscalar), intent(inout)    :: E
    character(*), intent(in)         :: newCoords
    logical, intent(in), optional    :: verbose
    integer                      :: nx,ny,nz,status,k,nlat,i1

    ! check to see if RHS (E1) is active (allocated)
    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for copy_rvector'
    else

       if(E%allocated) then
          call deall_rscalar(E)
       end if

       ! flip x & y coordinate dimensions
       nx = E1%ny
       ny = E1%nx
       nz = E1%nz
       i1=0
        ! allocate memory for v
        ! E%allocated will be true if all allocations succeed
        if (E1%gridType == CENTER) then
           allocate(E%v(nx,ny,nz), STAT=status)
        else if (E1%gridType == CORNER) then
           allocate(E%v(nx+1,ny+1,nz+1), STAT=status)
           i1=1
        else if(E1%gridType == CELL_EARTH) then
           nz = E1%grid%nzEarth
           allocate(E%v(nx,ny,nz), STAT=status)
        else
           write (0, *) 'gridType == ',trim(E1%gridType),' undefined in coordFlip_rscalar'
        end if
        E%allocated = status .EQ. 0

        if (E%allocated) then
           E%v = R_ZERO
        end if

          if (trim(newCoords) == trim(Cartesian)) then

             ! switch horizontal components and the direction of latitude
             nlat = size(E1%v,2)
             E%nx = E1%ny
             E%ny = E1%nx
             E%nz = E1%nz
             !E%v = E1%v !(:, nlat:1:-1, :) ! co-latitude South->latitude North
             do k=1,E1%nz+i1
                E%v(:,:,k) = transpose(E1%v(:,:,k))
             end do
             E%gridType = E1%gridType
             E%grid => E1%grid

          elseif (trim(newCoords) == trim(Spherical)) then

             ! switch horizontal components and the direction of latitude
             nlat = size(E1%v,1)
             E%nx = E1%ny
             E%ny = E1%nx
             E%nz = E1%nz
             !E%v = E1%v !(E1%nx+1:1:-1, :, :) ! latitude North->co-latitude South
             do k=1,E1%nz+i1
                E%v(:,:,k) = transpose(E1%v(:,:,k))
             end do
             E%gridType = E1%gridType
             E%grid => E1%grid

          else
             write (0, *) 'not compatible usage for coordFlip_rscalar'
          end if

    end if

    ! if the input was a temporary function output, deallocate
    if (E1%temporary) then
        call deall_rscalar(E1)
    end if

    if (present(verbose)) then
        write(*,*) '(coordFlip_rscalar) successfully flipped the vector to ',trim(newCoords)
        write(*,*) 'without changing the underlying grid...'
        write(*,*) 'this renders the output vector illegal. Do not use except as a quick fix.'
    end if

  end subroutine coordFlip_rscalar

end module sg_spherical	! sg_spherical
