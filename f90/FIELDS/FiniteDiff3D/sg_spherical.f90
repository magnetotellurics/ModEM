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


end module sg_spherical	! sg_spherical
