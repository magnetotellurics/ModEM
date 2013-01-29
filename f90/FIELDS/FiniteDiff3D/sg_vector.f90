! *****************************************************************************
module sg_vector
  ! This module creates data types for vector fields defined on
  ! edge/ face nodes of a staggered Cartesian grid, along with basic
  ! algebraic (vector space) operations. Example of these basic operations
  ! are allocation, deallocation, intialization, copying, and algebriac
  ! operations (linear combinations, scalar products, dot products).
  ! Belongs to SG_Basics class: staggered cartesian grid, data
  ! types defined on this grid, and operations defined on these data types. Not
  ! specific to EM problem, no dependency on outside (from other classes) modules.

  use math_constants		! math/ physics constants
  use utilities
  use griddef
  implicit none

  ! Important - overloading the '=' assignment
  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE copy_cvector
     MODULE PROCEDURE copy_rvector
  END INTERFACE

  ! Generic interfaces are done through subroutines
  ! creates edge/ face nodes
  INTERFACE create
     module procedure create_rvector
     module procedure create_cvector
  END INTERFACE

  ! deallocates the edge/ face nodes
  INTERFACE deall
     module procedure deall_rvector
     module procedure deall_cvector
  END INTERFACE

  ! set the values to zero
  INTERFACE zero
     module procedure zero_rvector
     module procedure zero_cvector
  END INTERFACE

  ! compatibility test
  INTERFACE compare
     module procedure compare_rvector_f
     module procedure compare_cvector_f
  END INTERFACE

  ! scalar value multiplies the edge/ face nodes
  INTERFACE scMult
     module procedure scMult_rvector
     module procedure scMult_cvector
     module procedure scMultReal_cvector
  END INTERFACE

  INTERFACE scMultAdd
     module procedure scMultAdd_cvector
  END INTERFACE

  INTERFACE linComb
     module procedure linComb_rvector
     module procedure linComb_cvector
  END INTERFACE

  ! adds the edge/ face nodes
  INTERFACE add
     module procedure add_rvector
     module procedure add_cvector
  END INTERFACE

  ! subtracts the edge/ face nodes
  INTERFACE subtract
     module procedure subtract_rvector
     module procedure subtract_cvector
  END INTERFACE

  ! pointwise vector (two vector data types) multiplication of edge/ face
  ! nodes
  ! and pointwise real-complex (mixed) multiplication of edge/ face nodes
  ! Both are vector data types
  INTERFACE diagMult
     module procedure diagMult_rvector
     module procedure diagMult_cvector
     module procedure diagMult_rcvector
     module procedure diagMult_crvector
  END INTERFACE

  ! pointwise real-complex (mixed) division of edge/ face nodes
  ! Both are vector data types
  INTERFACE diagDiv
     module procedure diagDiv_rvector
     module procedure diagDiv_rcvector
     module procedure diagDiv_crvector
  END INTERFACE

  INTERFACE dotProd
     MODULE PROCEDURE dotProd_rvector_f
     MODULE PROCEDURE dotProd_cvector_f
  END INTERFACE

  INTERFACE dotProd_noConj
     MODULE PROCEDURE dotProd_noConj_cvector_f
  END INTERFACE

  ! overload some intrinsic functions for complex numbers
  INTERFACE conjg
     MODULE PROCEDURE conjg_cvector_f
  END INTERFACE

  INTERFACE cmplx
     MODULE PROCEDURE cmplx_rvector_f
  END INTERFACE

  INTERFACE real
     MODULE PROCEDURE real_cvector_f
  END INTERFACE

  INTERFACE imag
     MODULE PROCEDURE imag_cvector_f
  END INTERFACE

  ! Interface operators are done through functions
!  INTERFACE OPERATOR (+)
!     MODULE PROCEDURE add_rvector_f
!     MODULE PROCEDURE add_cvector_f
!  END INTERFACE
!
!  INTERFACE OPERATOR (-)
!     MODULE PROCEDURE subtract_rvector_f
!     MODULE PROCEDURE subtract_cvector_f
!  END INTERFACE
!
!  INTERFACE OPERATOR (*)
!     MODULE PROCEDURE scMult_rvector_f
!     MODULE PROCEDURE scMult_cvector_f
!     MODULE PROCEDURE scMultReal_cvector_f
!     MODULE PROCEDURE diagMult_rvector_f
!     MODULE PROCEDURE diagMult_cvector_f
!     MODULE PROCEDURE diagMult_rcvector_f
!     MODULE PROCEDURE diagMult_crvector_f
!  END INTERFACE
!
!  INTERFACE OPERATOR (/)
!     MODULE PROCEDURE diagDiv_rcvector_f
!     MODULE PROCEDURE diagDiv_crvector_f
!  END INTERFACE

  public		::   create_rvector,  create_cvector, &
       deall_rvector, deall_cvector, &
       write_rvector, write_cvector, &
       read_rvector, read_cvector, &
       copy_rvector, copy_cvector, &
       zero_rvector, zero_cvector, &
       compare_rvector_f, compare_cvector_f, &
       scMult_rvector, scMult_cvector, scMultReal_cvector, &
       scMult_rvector_f, scMult_cvector_f, scMultReal_cvector_f, &
       add_rvector, add_cvector, &
       add_rvector_f, add_cvector_f, &
       subtract_rvector, subtract_cvector, &
       subtract_rvector_f, subtract_cvector_f, &
       diagMult_rvector, diagMult_cvector, &
       diagMult_rvector_f, diagMult_cvector_f, &
       diagMult_rcvector, diagMult_crvector, &
       diagMult_rcvector_f, diagMult_crvector_f, &
       diagDiv_rcvector, diagDiv_crvector, &
       diagDiv_rcvector_f, diagDiv_crvector_f, &
       dotProd_rvector_f, dotProd_cvector_f, &
       linComb_cvector, linComb_rvector, scMultAdd_cvector, conjg_cvector_f, &
       cmplx_rvector_f, real_cvector_f, imag_cvector_f


  ! ***************************************************************************
  ! type vector defines vector for either edge or face in a staggered grid as
  ! a complex field
  type :: cvector

     ! store the intention of the use in a character string defined
     ! in GridDef as a parameter: EDGE or FACE are two possibilities
     character (len=80)	                             :: gridType

     ! Typical usage:  electrical fields on cell edges of
     ! staggered grid
     ! For example, in an EDGE, the dimensions would be
     ! x: edge nodes in x-direction: dimension Nx, Ny+1, Nz+1
     ! y: edge nodes in y-direction: dimension Nx+1, Ny, Nz+1
     ! z: edge nodes in z-direction: dimension Nx+1, Ny+1, Nz
     ! Note that the arrays are defined through dynamic memory allocation
     complex(kind=prec), pointer, dimension(:,:,:)  :: x, y, z

     ! Grid Dimensions:
     ! nx is grid dimension (number of cells) in the x-direction
     ! ny is grid dimension (number of cells) in the y-direction
     ! nz is grid dimension (number of cells) in the z-direction:
     integer                                          :: nx = 0, ny = 0, nz = 0

     ! allocated:  .true.  x, y, z arrays have been allocated
     logical		                              :: allocated = .false.

     ! temporary:  .true. for function outputs only; necessary to avoid memory leaks
     ! (probably will not be needed in the future when compilers will support
     ! ISO/IEC 15581 - the "allocatable array extension")
     logical										:: temporary = .false.

     ! pointer to parent grid
     type (grid_t), pointer                             :: grid

  end type cvector

  ! ***************************************************************************
  ! type vector defines vector for either edge or face in a staggered grid as
  ! a real field
  type :: rvector

     ! store the intention of the use in a character string defined
     ! in GridDef as a parameter: EDGE or FACE are two possibilities
     character (len=80)	                              :: gridType

     ! Typical usage:  conductivity averaged on cell edges of
     ! staggered grid
     ! x: edge nodes in x-direction: dimension Nx, Ny+1, Nz+1
     ! y: edge nodes in y-direction: dimension Nx+1, Ny, Nz+1
     ! z: edge nodes in z-direction: dimension Nx+1, Ny+1, Nz
     ! Note that the arrays are defined through dynamic memory allocation
     real(kind=prec), pointer, dimension(:,:,:)      :: x, y, z

     ! Grid Dimensions:
     ! nx is grid dimension (number of cells) in the x-direction
     ! ny is grid dimension (number of cells) in the y-direction
     ! nz is grid dimension (number of cells) in the z-direction:
     integer                                          :: nx = 0, ny = 0, nz = 0


     ! allocated:  .true.  x, y, z arrays have been allocated
     logical		                              :: allocated = .false.

     ! temporary:  .true. for function outputs only; necessary to avoid memory leaks
     ! (probably will not be needed in the future when compilers will support
     ! ISO/IEC 15581 - the "allocatable array extension")
     logical										:: temporary = .false.

     ! pointer to parent grid
     type (grid_t), pointer                             :: grid

  end type rvector

Contains
  ! CREATE GRID_edge/ face VECTORS
  ! * subroutine create_rvector(igrid, E, gridType)
  ! * subroutine create_cvector(igrid, E, gridType)

  ! DEALLOCATE GRID_edge/ face VECTORS
  ! * subroutine deall_rvector(E)
  ! * subroutine deall_cvector(E)

  ! COPY GRID_edge/ face VECTORS :  (=)
  ! * subroutine copy_rvector(E2,E1)
  ! * subroutine copy_cvector(E2,E1)

  ! ZERO GRID_edge/ face VECTORS
  ! * subroutine zero_rvector(E)
  ! * subroutine zero_cvector(E)

  ! SCALAR MULTIPLICATION:
  ! * subroutine scMult_rvector(c, E1, E2)
  ! * subroutine scMult_cvector(c, E1, E2)
  ! * subroutine scMultReal_cvector(c, E1, E2)

  ! SCALAR MULTIPLICATION: FUNCTION VERSION (c*E)
  ! * function scMult_rvector_f(c, E1) result(E2)
  ! * function scMult_cvector_f(c, E1) result(E2)
  ! * function scMultReal_cvector_f(c, E1) result(E2)

  ! VECTOR SUM:
  ! * subroutine add_rvector(E1, E2, E3)
  ! * subroutine add_cvector(E1, E2, E3)

  ! VECTOR SUM: FUNCTION VERSION (E3 = E1+E2)
  ! * function add_rvector_f(E1, E2) result(E3)
  ! * function add_cvector_f(E1, E2) result(E3)

  ! VECTOR SUBTRACT:
  ! * subroutine subtract_rvector(E1, E2, E3)
  ! * subroutine subtract_cvector(E1, E2, E3)

  ! VECTOR SUBTRACT: FUNCTION VERSION (E3 = E1-E2)
  ! * function subtract_rvector_f(E1, E2) result(E3)
  ! * function subtract_cvector_f(E1, E2) result(E3)

  ! POINTWISE MULTIPLICATION OF VECTORS:
  ! * subroutine diagMult_rvector(E1, E2, E3)
  ! * subroutine diagMult_cvector(E1, E2, E3)

  ! POINTWISE MULTIPLICATION OF VECTORS: FUNCTION VERSION (c*E)
  ! * function diagMult_rvector_f(E1, E2) result(E3)
  ! * function diagMult_cvector_f(E1, E2) result(E3)

  ! POINTWISE MULTIPLICATION OF VECTORS and SCALARS:
  ! one combination is real and complex vectors and another is complex
  ! and real vectors as the sequence for inputs
  ! * subroutine diagMult_crvector(E1, E2, E3)
  ! * subroutine diagMult_rcvector(E1, E2, E3)

  ! POINTWISE MULTIPLICATION OF VECTORS and SCALARS: FUNCTION VERSION (c*E)
  ! one combination is real and complex vectors and another is complex
  ! and real vectors as the sequence for inputs
  ! * function diagMult_crvector_f(E1, E2) result(E3)
  ! * function diagMult_rcvector_f(E1, E2) result(E3)

  ! POINTWISE DIVISION OF VECTORS and SCALARS:
  ! one combination is real and complex vectors and another is complex
  ! and real vectors as the sequence for inputs
  ! * subroutine diagDiv_crvector(E1, E2, E3)
  ! * subroutine diagDiv_rcvector(E1, E2, E3)

  ! POINTWISE DIVISION OF VECTORS and SCALARS: FUNCTION VERSION (c*E)
  ! one combination is real and complex vectors and another is complex
  ! and real vectors as the sequence for inputs
  ! * function diagDiv_crvector_f(E1, E2) result(E3)
  ! * function diagDiv_rcvector_f(E1, E2) result(E3)

  ! VECTOR DOT PRODUCT
  ! * function dotProd_rvector(E1, E2) result(r)
  ! * function dotProd_cvector(E1, E2) result(c)

  ! COMPLEX LINEAR COMBINATION ... formally redundent with above
  ! only doing ones that we are sure to want
  ! * subroutine linCom_cvector(inc1, E1, inc2, E2, E3)
  ! * subroutine scMultAdd_V_node(c, E1, E2)  (E2 = E2+c*E1)

  ! COMBINE REAL VECTORS TO PRODUCE A COMPLEX VECTOR, AND THE CONVERSE
  ! * function conjg_cvector(E1) result (E2)
  ! * function cmplx_rvector(E1, E2) result (E3)
  ! * function real_cvector(E1) result (E2)
  ! * function imag_cvector(E1) result (E2)

  ! The algebraic routines expect all input and output
  ! variables to be of the correct type, already allocated,
  ! and of the correct size.


  !****************************************************************************
  ! create_rvector creates variable of derived type rvector,
  ! using grid definition in structure "grid" ;
  ! allocates memory in x,y,z component arrays
  ! gridType is a character string to describe intended usage
  subroutine create_rvector(igrid, E, gridType)

    implicit none
    type(grid_t), target, intent(in)    :: igrid
    ! the grid for which an edge/ face node field is being initialized
    type (rvector), intent(inout)      :: E

    integer                            :: status,nx,ny,nz

    character (len=80), intent(in)     :: gridType

	! First deallocate anything, that's allocated
	call deall_rvector(E)

    ! Set pointer
    E%grid => igrid

    ! Grid dimensions
    nx = igrid%nx
    ny = igrid%ny
    nz = igrid%nz
    E%nx = nx
    E%ny = ny
    E%nz = nz

    ! gridType
    E%gridType = gridType

    ! allocate memory for x,y,z
    ! E%allocated will be true if all allocations succeed
    E%allocated = .true.
    if (E%gridType == EDGE) then
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
    else if (E%gridType == FACE) then
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
    else
       write (0, *) 'not a known tag'
    end if

    if (E%allocated) then
       E%x = R_ZERO
       E%y = R_ZERO
       E%z = R_ZERO
    else
        write (0, *) 'Warning: unable to allocate rvector - invalid grid supplied'
    end if

  end subroutine create_rvector  ! create_rvector


  !****************************************************************************
  ! create_cvector creates variable of derived type cvector,
  ! using grid definition in structure "grid" ;
  ! allocates memory in x,y,z component arrays
  ! gridType is a character string to describe intended usage
  subroutine create_cvector(igrid, E, gridType)

    implicit none
    type(grid_t), target, intent(in)     :: igrid
    ! the grid for which an edge/ face node field is being initialized
    type (cvector), intent(inout)       :: E

    integer                             :: status,nx,ny,nz

    character (len=80), intent(in)      :: gridType

	! First deallocate anything, that's allocated
    call deall_cvector(E)

    ! Set pointer
    E%grid => igrid

    ! Grid dimensions
    nx = igrid%nx
    ny = igrid%ny
    nz = igrid%nz
    E%nx = nx
    E%ny = ny
    E%nz = nz
    ! print *, 'nx, ny, nz', nx, ny, nz

    ! gridType
    E%gridType = gridType

    ! allocate memory for x,y,z ;
    ! E%allocated will be true if all allocations succeed
    E%allocated = .true.
    if (E%gridType == EDGE) then
	   ! For spherical problem:
	   ! 1) E%x(:,1,:) and E%x(:,ny+1,:) are undefined,
	   ! 2) E%y(nx+1,:,:) and E%z(nx+1,:,:) are repetitious,
	   ! 3) E%z(:,1,k) and E%z(:,ny+1,k) indep. of i for a fixed k,
	   ! 4) E%x(-1,0,ny+1,ny+2,:,:) will be needed for interpolation.
       allocate(E%x(nx,ny+1,nz+1), STAT=status)
       E%allocated = E%allocated .and. (status .EQ. 0)
       allocate(E%y(nx+1,ny,nz+1), STAT=status)
       E%allocated = E%allocated .and. (status .EQ. 0)
       allocate(E%z(nx+1,ny+1,nz), STAT=status)
       E%allocated = E%allocated .and. (status .EQ. 0)
    else if (E%gridType == FACE) then
	   ! For spherical problem:
	   ! 1) E%y(:,1,:) and E%y(:,ny+1,:) are undefined,
	   ! 2) E%x(nx+1,:,:) is repetitious and equals E%x(1,:,:).
       allocate(E%x(nx+1,ny,nz), STAT=status)
       E%allocated = E%allocated .and. (status .EQ. 0)
       allocate(E%y(nx,ny+1,nz), STAT=status)
       E%allocated = E%allocated .and. (status .EQ. 0)
       allocate(E%z(nx,ny,nz+1), STAT=status)
       E%allocated = E%allocated .and. (status .EQ. 0)
    else
       write (0, *) 'not a known tag'
    end if

    if (E%allocated) then
       E%x = C_ZERO
       E%y = C_ZERO
       E%z = C_ZERO
    else
        write (0, *) 'Warning: unable to allocate cvector - invalid grid supplied'
    end if

  end subroutine create_cvector  ! create_cvector


  !****************************************************************************
  ! deall_rvector destoys variable of derived type rvector,
  ! deallocating memory
  subroutine deall_rvector(E)

    implicit none
    type (rvector)  :: E
    integer	    :: status

    ! deallocate memory for x,y,z
    if (E%allocated) then 
    deallocate(E%x, STAT=status)
    deallocate(E%y, STAT=status)
    deallocate(E%z, STAT=status)
    end if
	!if(associated(E%x)) deallocate(E%x, STAT=status)
	!if(associated(E%y)) deallocate(E%y, STAT=status)
	!if(associated(E%z)) deallocate(E%z, STAT=status)
    if(associated(E%grid)) nullify(E%grid)

    E%nx = 0
    E%ny = 0
    E%nz = 0
    E%gridType = ''
    E%allocated = .false.

  end subroutine deall_rvector  ! deall_rvector


  !****************************************************************************
  ! deall_cvector destoys variable of derived type cvector,
  ! deallocating memory
  subroutine deall_cvector(E)

    implicit none
    type (cvector)  :: E
    integer	    :: status

    ! deallocate memory for x,y,z
    if (E%allocated) then 
    deallocate(E%x, STAT=status)
    deallocate(E%y, STAT=status)
    deallocate(E%z, STAT=status)
    end if
    

	!if(associated(E%x)) deallocate(E%x, STAT=status)
	!if(associated(E%y)) deallocate(E%y, STAT=status)
	!if(associated(E%z)) deallocate(E%z, STAT=status)
    if(associated(E%grid)) nullify(E%grid)

    E%nx = 0
    E%ny = 0
    E%nz = 0
    E%gridType = ''
    E%allocated = .false.

  end subroutine deall_cvector  ! deall_cvector

  !****************************************************************************
  ! write_rvector writes an rvector in a simple ASCII format; rvector has
  ! to exist and be allocated before calling this routine, and the file unit
  ! must already be available for writing.
  subroutine write_rvector(fid, E)

      integer,        intent(in)		:: fid
      type (rvector), intent(in)		:: E

      !  local variables
      integer 		                    :: Nx, Ny, Nz
      integer                           :: i, j, k, istat
      real (kind(E%x)), allocatable, dimension(:,:,:)  :: x, y, z

      if(.not. E%allocated) then
         write(0, *) 'rvector must be allocated before call to write_rvector'
         return
      endif

      write(fid,'(3i5,a10)',iostat=istat) E%nx,E%ny,E%nz,trim(E%gridType)

	  Nx = E%nx
	  Ny = E%ny
	  Nz = E%nz

      allocate(x(Nx+1,Ny+1,Nz+1),y(Nx+1,Ny+1,Nz+1),z(Nx+1,Ny+1,Nz+1),STAT=istat)
      x = R_ZERO
      y = R_ZERO
      z = R_ZERO

      if (E%gridType == EDGE) then
	   ! For spherical problem:
	   ! 1) E%x(:,1,:) and E%x(:,ny+1,:) are undefined,
	   ! 2) E%y(nx+1,:,:) and E%z(nx+1,:,:) are repetitios,
	   ! 3) E%x(-1,0,ny+1,ny+2,:,:) will be needed for interpolation.
	   x(1:Nx,:,:) = E%x
	   y(:,1:Ny,:) = E%y
	   z(:,:,1:Nz) = E%z
      else if (E%gridType == FACE) then
	   ! For spherical problem:
	   ! 1) E%y(:,1,:) and E%y(:,ny+1,:) are undefined,
	   ! 2) E%x(nx+1,:,:) is repetitios and equals E%x(1,:,:).
	   x(:,1:Ny,1:Nz) = E%x
	   y(1:Nx,:,1:Nz) = E%y
	   z(1:Nx,1:Ny,:) = E%z
      else
       write (0, *) 'not a known tag'
      end if

      do k = 1,Nz+1
      	do j = 1,Ny+1
      		do i = 1,Nx+1
      			write(fid,'(3i5,3es14.6)',iostat=istat) i,j,k,x(i,j,k),y(i,j,k),z(i,j,k)
      		end do
      	end do
      end do

	  deallocate(x,y,z,STAT=istat)

  end subroutine write_rvector

  !****************************************************************************
  ! write_cvector writes a cvector in a simple ASCII format; cvector has
  ! to exist and be allocated before calling this routine, and the file unit
  ! must already be available for writing.
  subroutine write_cvector(fid, E)

      integer,        intent(in)		:: fid
      type (cvector), intent(in)		:: E

      !  local variables
      integer 		                    :: Nx, Ny, Nz
      integer                           :: i, j, k, istat
      complex (kind(E%x)), allocatable, dimension(:,:,:)  :: x, y, z

      if(.not. E%allocated) then
         write(0, *) 'cvector must be allocated before call to write_cvector'
         return
      endif

      write(fid,'(3i5,a10)',iostat=istat) E%nx,E%ny,E%nz,trim(E%gridType)

	  Nx = E%nx
	  Ny = E%ny
	  Nz = E%nz

      allocate(x(Nx+1,Ny+1,Nz+1),y(Nx+1,Ny+1,Nz+1),z(Nx+1,Ny+1,Nz+1),STAT=istat)
      x = C_ZERO
      y = C_ZERO
      z = C_ZERO

      if (E%gridType == EDGE) then
	   ! For spherical problem:
	   ! 1) E%x(:,1,:) and E%x(:,ny+1,:) are undefined,
	   ! 2) E%y(nx+1,:,:) and E%z(nx+1,:,:) are repetitios,
	   ! 3) E%x(-1,0,ny+1,ny+2,:,:) will be needed for interpolation.
	   x(1:Nx,:,:) = E%x
	   y(:,1:Ny,:) = E%y
	   z(:,:,1:Nz) = E%z
      else if (E%gridType == FACE) then
	   ! For spherical problem:
	   ! 1) E%y(:,1,:) and E%y(:,ny+1,:) are undefined,
	   ! 2) E%x(nx+1,:,:) is repetitios and equals E%x(1,:,:).
	   x(:,1:Ny,1:Nz) = E%x
	   y(1:Nx,:,1:Nz) = E%y
	   z(1:Nx,1:Ny,:) = E%z
      else
       write (0, *) 'not a known tag'
      end if

      do k = 1,Nz+1
      	do j = 1,Ny+1
      		do i = 1,Nx+1
      			write(fid,'(3i5,6es14.6)',iostat=istat) i,j,k,x(i,j,k),y(i,j,k),z(i,j,k)
      		end do
      	end do
      end do

	  deallocate(x,y,z,STAT=istat)

  end subroutine write_cvector

  !****************************************************************************
  ! read_rvector reads an rvector in a simple ASCII format; rvector must match
  ! the input grid; file unit must already be available for reading.
  subroutine read_rvector(fid, E, grid)

      integer,        intent(in)		:: fid
      type (rvector), intent(inout)		:: E
      type (grid_t), target, intent(in) :: grid

      !  local variables
      integer 		                    :: Nx, Ny, Nz
      character(80)						:: gridType
      integer                           :: i, j, k, ii, jj, kk, istat
      real (kind(E%x)), allocatable, dimension(:,:,:)  :: x, y, z

      read(fid,*,iostat=istat) Nx,Ny,Nz,gridType

      if ((Nx .ne. grid%nx) .or. (Ny .ne. grid%ny) .or. (Nz .ne. grid%nz)) then
         write(0, *) 'rvector size does not match the grid in read_rvector'
         return
      endif

      if(.not. E%allocated) then
      	call create_rvector(grid,E,gridType)
      else
      	E%nx = Nx
      	E%ny = Ny
      	E%nz = Nz
      	E%gridType = gridType
      	E%grid => grid
      endif

      allocate(x(Nx+1,Ny+1,Nz+1),y(Nx+1,Ny+1,Nz+1),z(Nx+1,Ny+1,Nz+1),STAT=istat)
      x = R_ZERO
      y = R_ZERO
      z = R_ZERO

      do k = 1,Nz+1
      	do j = 1,Ny+1
      		do i = 1,Nx+1
      			read(fid,*,iostat=istat) ii,jj,kk,x(i,j,k),y(i,j,k),z(i,j,k)
      		end do
      	end do
      end do

      if (E%gridType == EDGE) then
	   ! For spherical problem:
	   ! 1) E%x(:,1,:) and E%x(:,ny+1,:) are undefined,
	   ! 2) E%y(nx+1,:,:) and E%z(nx+1,:,:) are repetitios,
	   ! 3) E%x(-1,0,ny+1,ny+2,:,:) will be needed for interpolation.
	   E%x = x(1:Nx,:,:)
	   E%y = y(:,1:Ny,:)
	   E%z = z(:,:,1:Nz)
      else if (E%gridType == FACE) then
	   ! For spherical problem:
	   ! 1) E%y(:,1,:) and E%y(:,ny+1,:) are undefined,
	   ! 2) E%x(nx+1,:,:) is repetitios and equals E%x(1,:,:).
	   E%x = x(:,1:Ny,1:Nz)
	   E%y = y(1:Nx,:,1:Nz)
	   E%z = z(1:Nx,1:Ny,:)
      else
       write (0, *) 'not a known tag'
      end if

	  deallocate(x,y,z,STAT=istat)

  end subroutine read_rvector

  !****************************************************************************
  ! read_cvector reads a cvector in a simple ASCII format; cvector must match
  ! the input grid; file unit must already be available for reading.
  subroutine read_cvector(fid, E, grid)

      integer,        intent(in)		:: fid
      type (cvector), intent(inout)		:: E
      type (grid_t), target, intent(in) :: grid

      !  local variables
      integer 		                    :: Nx, Ny, Nz
      character(80)						:: gridType
      integer                           :: i, j, k, ii, jj, kk, istat
      real (kind(E%x))                  :: xr, xi, yr, yi, zr, zi
      complex (kind(E%x)), allocatable, dimension(:,:,:)  :: x, y, z

      read(fid,*,iostat=istat) Nx,Ny,Nz,gridType

      if ((Nx .ne. grid%nx) .or. (Ny .ne. grid%ny) .or. (Nz .ne. grid%nz)) then
         write(0, *) 'cvector size does not match the grid in read_cvector'
         return
      endif

      if(.not. E%allocated) then
      	call create_cvector(grid,E,gridType)
      else
      	E%nx = Nx
      	E%ny = Ny
      	E%nz = Nz
      	E%gridType = gridType
      	E%grid => grid
      endif

      allocate(x(Nx+1,Ny+1,Nz+1),y(Nx+1,Ny+1,Nz+1),z(Nx+1,Ny+1,Nz+1),STAT=istat)
      x = C_ZERO
      y = C_ZERO
      z = C_ZERO

      do k = 1,Nz+1
      	do j = 1,Ny+1
      		do i = 1,Nx+1
      		    ! makes the acceptable formatting more flexible
      			read(fid,*,iostat=istat) ii,jj,kk,xr,xi,yr,yi,zr,zi
      			x(i,j,k) = cmplx(xr,xi)
      			y(i,j,k) = cmplx(yr,yi)
      			z(i,j,k) = cmplx(zr,zi)
      		end do
      	end do
      end do

      if (E%gridType == EDGE) then
	   ! For spherical problem:
	   ! 1) E%x(:,1,:) and E%x(:,ny+1,:) are undefined,
	   ! 2) E%y(nx+1,:,:) and E%z(nx+1,:,:) are repetitios,
	   ! 3) E%x(-1,0,ny+1,ny+2,:,:) will be needed for interpolation.
	   E%x = x(1:Nx,:,:)
	   E%y = y(:,1:Ny,:)
	   E%z = z(:,:,1:Nz)
      else if (E%gridType == FACE) then
	   ! For spherical problem:
	   ! 1) E%y(:,1,:) and E%y(:,ny+1,:) are undefined,
	   ! 2) E%x(nx+1,:,:) is repetitios and equals E%x(1,:,:).
	   E%x = x(:,1:Ny,1:Nz)
	   E%y = y(1:Nx,:,1:Nz)
	   E%z = z(1:Nx,1:Ny,:)
      else
       write (0, *) 'not a known tag'
      end if

	  deallocate(x,y,z,STAT=istat)

  end subroutine read_cvector

  !****************************************************************************
  ! copy_rvector makes an exact copy of derived data type
  ! rvector;   NOTE: first argument is output
  subroutine copy_rvector(E2, E1)

    implicit none
    type (rvector), intent(in)       :: E1
    type (rvector), intent(inout)    :: E2
    integer	                     :: status

    ! check to see if RHS (E1) is active (allocated)
    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for copy_rvector'
    else

       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! just copy components
             E2%x = E1%x
             E2%y = E1%y
             E2%z = E1%z
             E2%gridType = E1%gridType
             E2%grid => E1%grid

          else
             write (0, *) 'not compatible usage for copy_rvector'
          end if

       else

          if(E2%allocated) then
             ! first deallocate memory for x,y,z
             deallocate(E2%x, E2%y, E2%z,STAT=status)
          end if

          !  then allocate E2 as correct size ...
          Call create_rvector(E1%grid, E2, E1%gridType)
          !   .... and copy E1
          E2%x = E1%x
          E2%y = E1%y
          E2%z = E1%z
          E2%gridType = E1%gridType
          E2%grid => E1%grid

       end if

    end if

    ! if the input was a temporary function output, deallocate
    if (E1%temporary) then
    	call deall_rvector(E1)
    end if

  end subroutine copy_rvector  ! copy_rvector


  !****************************************************************************
  ! copy_cvector makes an exact copy of derived data type
  ! cvector;
  subroutine copy_cvector(E2, E1)

    implicit none
    type (cvector), intent(in)            :: E1
    type (cvector), intent(inout)         :: E2
    integer                               :: status

    ! check to see if RHS (E1) is active (allocated)
    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for copy_cvector'
    else

       if((E2%nx == E1%nx).and.(E2%ny == E1%ny).and.(E2%nz == E1%nz)) then

          if  (E1%gridType == E2%gridType) then

             ! just copy components
             E2%x = E1%x
             E2%y = E1%y
             E2%z = E1%z
             E2%gridType = E1%gridType
             E2%grid => E1%grid

          else
             write (0, *) 'not compatible usage for copy_cvector'
          end if

       else

          if(E2%allocated) then
             ! first deallocate memory for x,y,z
             deallocate(E2%x, E2%y, E2%z,STAT=status)
          end if

          !  then allocate E2 as correct size ...
          Call create_cvector(E1%grid, E2, E1%gridType)
          !   .... and copy E1
          E2%x = E1%x
          E2%y = E1%y
          E2%z = E1%z
          E2%gridType = E1%gridType
          E2%grid => E1%grid

       end if

    end if

    ! if the input was a temporary function output, deallocate
    if (E1%temporary) then
    	call deall_cvector(E1)
    end if

  end subroutine copy_cvector  ! copy_cvector


  ! ***************************************************************************
  ! zero_rvector zeros variable of derived data type
  ! rvector;
  subroutine zero_rvector(E)

    implicit none
    type (rvector), intent(inout)   :: E

    ! check to see if E is active (allocated)
    if(.not.E%allocated) then
       write(0,*) 'Error in zero_rvector: E not allocated'
    else

       E%x = R_ZERO
       E%y = R_ZERO
       E%z = R_ZERO

    end if

  end subroutine zero_rvector


  ! ***************************************************************************
  ! zero_cvector zeros variable of derived data type
  ! cvector;
  subroutine zero_cvector(E)

    implicit none
    type (cvector), intent(inout) :: E

    ! check to see if E is active (allocated)
    if(.not.E%allocated) then
       write(0,*) 'Error in zero_cvector: E not allocated'
    else

       E%x = C_ZERO
       E%y = C_ZERO
       E%z = C_ZERO

    end if

  end subroutine zero_cvector ! zero_cvector


  ! **********************************************************************
  ! * Creates a random perturbation in rvector - used for testing

  subroutine random_rvector(E,eps)

    implicit none
    type (rvector), intent(inout)                    :: E
    real(8), intent(in), optional                    :: eps

    if (.not. E%allocated) then
      call warning('cvector not allocated in random_rvector')
      return
    end if

    call zero_rvector(E)

    call random_number(E%x)
    call random_number(E%y)
    call random_number(E%z)
    if (present(eps)) then
        E%x = E%x * eps
        E%y = E%y * eps
        E%z = E%z * eps
    else
        E%x = E%x * 0.05
        E%y = E%y * 0.05
        E%z = E%z * 0.05
    end if

  end subroutine random_rvector


  ! **********************************************************************
  ! * Creates a random perturbation in cvector - used for testing

  subroutine random_cvector(E,eps)

    implicit none
    type (cvector), intent(inout)                    :: E
    real(8), intent(in), optional                    :: eps
    ! local
    real (kind(E%x)), allocatable, dimension(:,:,:)  :: xr,xi,yr,yi,zr,zi
    integer              :: Nx,Ny,Nz,istat

    if (.not. E%allocated) then
      call warning('cvector not allocated in random_cvector')
      return
    end if

    call zero_cvector(E)

    ! make some random vectors
    Nx = E%nx
    Ny = E%ny
    Nz = E%nz
    allocate(xr(Nx+1,Ny+1,Nz+1),xi(Nx+1,Ny+1,Nz+1),STAT=istat)
    allocate(yr(Nx+1,Ny+1,Nz+1),yi(Nx+1,Ny+1,Nz+1),STAT=istat)
    allocate(zr(Nx+1,Ny+1,Nz+1),zi(Nx+1,Ny+1,Nz+1),STAT=istat)
    call random_number(xr)
    call random_number(xi)
    call random_number(yr)
    call random_number(yi)
    call random_number(zr)
    call random_number(zi)
    if (present(eps)) then
        xr = xr * eps
        xi = xi * eps
        yr = yr * eps
        yi = yi * eps
        zr = zr * eps
        zi = zi * eps
    else
        xr = xr * 0.05
        xi = xi * 0.05
        yr = yr * 0.05
        yi = yi * 0.05
        zr = zr * 0.05
        zi = zi * 0.05
    end if

	if (E%gridType == EDGE) then
	   E%x = cmplx(xr(1:Nx,:,:),xi(1:Nx,:,:))
	   E%y = cmplx(yr(:,1:Ny,:),yi(:,1:Ny,:))
	   E%z = cmplx(zr(:,:,1:Nz),zi(:,:,1:Nz))
	else if (E%gridType == FACE) then
	   E%x = cmplx(xr(:,1:Ny,1:Nz),xi(:,1:Ny,1:Nz))
	   E%y = cmplx(yr(1:Nx,:,1:Nz),yi(1:Nx,:,1:Nz))
	   E%z = cmplx(zr(1:Nx,1:Ny,:),zi(1:Nx,1:Ny,:))
	else
	   write (0, *) 'not a known tag'
	end if

    deallocate(xr,xi,yr,yi,zr,zi,STAT=istat)

  end subroutine random_cvector

  ! **********************************************************************
  ! * check two vectors for compatibility for linear operator purposes

  function compare_cvector_f(E1,E2) result (status)

    type (cvector), intent(in)                  :: E1
    type (cvector), intent(in)                  :: E2
    logical                                     :: status

	status = .FALSE.

    ! Check whether all the vector nodes are of the same size
    if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then
       if (E1%gridType == E2%gridType) then
          status = .TRUE.
       end if
    end if


  end function compare_cvector_f

  ! **********************************************************************
  ! * check two vectors for compatibility for linear operator purposes

  function compare_rvector_f(E1,E2) result (status)

    type (rvector), intent(in)                  :: E1
    type (rvector), intent(in)                  :: E2
    logical                                     :: status

	status = .FALSE.

    ! Check whether all the vector nodes are of the same size
    if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then
       if (E1%gridType == E2%gridType) then
          status = .TRUE.
       end if
    end if


  end function compare_rvector_f

  ! ***************************************************************************
  ! scMult_cvector multiplies vector stored as derived data type
  ! cvector with a complex scalar; subroutine version
  ! E2 can overwrite E1
  subroutine scMult_cvector(c, E1, E2)

    implicit none
    complex(kind=prec), intent(in)                      :: c
    ! a complex scalar to be multiplied with
    type (cvector), intent(in)                       :: E1
    type (cvector), intent(inout)                    :: E2

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMult_cvector'
       return
    endif

    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated yet for scMult_cvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! complex scalar multiplication for x,y,z-components
             E2%x = E1%x * c
             E2%y = E1%y * c
             E2%z = E1%z * c

          else
             write (0, *) 'not compatible usage for scMult_cvector'
          end if

       else
          write(0, *) 'Error:scMult_cvector: vectors not same size'

       end if
    end if

  end subroutine scMult_cvector ! scMult_cvector


  ! ***************************************************************************
  ! scMult_cvector_f multiplies vector stored as derived data type
  ! cvector with a complex scalar; function version
  function scMult_cvector_f(c, E1) result(E2)

    implicit none
    complex(kind=prec), intent(in)                      :: c
    ! a complex scalar to be multiplied with
    type (cvector), intent(in)                       :: E1
    type (cvector)                                   :: E2

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMult_cvector_f'
       return
    endif

    ! In function version, appropriate data types need to be created
    Call create_cvector(E1%grid, E2, E1%gridType)
    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated for scMult_cvector_f'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! complex scalar multiplication for x,y,z-components
             E2%x = E1%x * c
             E2%y = E1%y * c
             E2%z = E1%z * c

          else
             write (0, *) 'not compatible usage for scMult_cvector_f'
          end if

       else

          write(0, *) 'Error:scMult_cvector_f: vectors not same size'

       end if
    end if

    E2%temporary = .true.

  end function scMult_cvector_f ! scMult_cvector_f


  ! ***************************************************************************
  ! scMultReal_cvector multiplies vector stored as derived data type
  ! cvector with a real scalar; subroutine version
  ! E2 can overwrite E1
  subroutine scMultReal_cvector(c, E1, E2)

    implicit none
    real (kind=prec), intent(in)                         :: c
    ! a real scalar to be multiplied with
    type (cvector), intent(in)                       :: E1
    type (cvector), intent(inout)                    :: E2

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMultReal_cvector'
       return
    endif

    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated yet for scMultReal_cvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! complex scalar multiplication for x,y,z-components
             E2%x = E1%x * c
             E2%y = E1%y * c
             E2%z = E1%z * c

          else
             write (0, *) 'not compatible usage for scMultReal_cvector'
          end if

       else
          write(0, *) 'Error:scMultReal_cvector: vectors not same size'

       end if
    end if

  end subroutine scMultReal_cvector ! scMultReal_cvector


  !****************************************************************************
  ! scMult_cvector_f multiplies vector stored as derived data type
  ! cvector with a real scalar; function version
  function scMultReal_cvector_f(c, E1) result(E2)

    implicit none
    real (kind=prec), intent(in)			     :: c
    ! a real scalar to be multiplied with
    type (cvector), intent(in)                       :: E1
    type (cvector)                                   :: E2

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMultReal_cvector_f'
       return
    endif

    ! In function version, appropriate data types need to be created
    Call create_cvector(E1%grid, E2, E1%gridType)
    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated for scMultReal_cvector_f'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! complex scalar multiplication for x,y,z-components
             E2%x = E1%x * c
             E2%y = E1%y * c
             E2%z = E1%z * c

          else
             write (0, *) 'not compatible usage for scMultReal_cvector_f'
          end if

       else

          write(0, *) 'Error:scMultReal_cvector_f: vectors not same size'

       end if
    end if

    E2%temporary = .true.

  end function scMultReal_cvector_f ! scMultReal_cvector_f


  ! ***************************************************************************
  ! scMult_rvector multiplies vector stored as derived data type
  ! rvector with a real scalar; subroutine version
  ! E2 can overwrite E1
  subroutine scMult_rvector(c, E1, E2)

    implicit none
    real (kind=prec), intent(in)                         :: c
    ! a real scalar to be multiplied with
    type (rvector), intent(in)                       :: E1
    type (rvector), intent(inout)                    :: E2

   if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMult_rvector'
       return
    endif

    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated for scMult_rvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! real scalar multiplication for x,y,z-components
             E2%x = E1%x * c
             E2%y = E1%y * c
             E2%z = E1%z * c

          else
             write (0, *) 'not compatible usage for scMult_rvector'
          end if

       else

          write(0, *) 'Error:scMult_rvector: vectors not same size'

       end if
    end if

  end subroutine scMult_rvector ! scMult_rvector


  !****************************************************************************
  ! scMult_rvector_f multiplies vector stored as derived data type
  ! rvector with a real scalar; function version
  function scMult_rvector_f(c, E1) result(E2)

    implicit none
    real (kind=prec), intent(in)                         :: c
    ! a complex scalar to be multiplied with
    type (rvector), intent(in)                       :: E1
    type (rvector)                                   :: E2

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMult_rvector_f'
       return
    endif

    ! In function version, appropriate data types need to be created
    Call create_rvector(E1%grid, E2, E1%gridType)
    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated for scMult_rvector_f'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! real scalar multiplication for x,y,z-components
             E2%x = E1%x * c
             E2%y = E1%y * c
             E2%z = E1%z * c

          else
             write (0, *) 'not compatible usage for scMult_rvector_f'
          end if

       else

          write(0, *) 'Error:scMult_rvector_f: vectors not same size'

       end if
    end if

    E2%temporary = .true.

  end function scMult_rvector_f ! scMult_rvector_f


  !****************************************************************************
  ! add_rvector adds vectors stored as derived data type
  ! rvector with ; subroutine version
  ! E3 can overwrite E1 and E2
  subroutine add_rvector(E1, E2, E3)

    implicit none
    type (rvector), intent(in)               :: E1, E2
    type (rvector), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for add_rvector'
       return
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for add_rvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! add x,y,z-components
             E3%x = E1%x + E2%x
             E3%y = E1%y + E2%y
             E3%z = E1%z + E2%z

          else
             write (0, *) 'not compatible usage for add_rvector'
          end if

       else

          write(0, *) 'Error:add_rvector: vectors not same size'

       end if
    end if

  end subroutine add_rvector ! add_rvector


  !****************************************************************************
  ! add_rvector_f adds vectors stored as derived data type
  ! rvector with ; function version
  function add_rvector_f(E1, E2) result(E3)

    implicit none
    type (rvector), intent(in)               :: E1, E2
    type (rvector)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for add_rvector_f'
       return
    endif

    ! In function version, appropriate data types need to be created
    Call create_rvector(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for add_rvector_f'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! add x,y,z-components
             E3%x = E1%x + E2%x
             E3%y = E1%y + E2%y
             E3%z = E1%z + E2%z

          else
             write (0, *) 'not compatible usage for add_rvector_f'
          end if

       else

          write(0, *) 'Error:add_rvector_f: vectors not same size'

       end if
    end if

    E3%temporary = .true.

  end function add_rvector_f ! add_rvector_f


  !****************************************************************************
  ! add_cvector adds vectors stored as derived data type
  ! cvector with ; subroutine version
  ! E3 can overwrite E1 and E2
  subroutine add_cvector(E1, E2, E3)

    implicit none
    type (cvector), intent(in)               :: E1, E2
    type (cvector), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for add_cvector'
       return
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS not allocated for add_cvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! add x,y,z-components
             E3%x = E1%x + E2%x
             E3%y = E1%y + E2%y
             E3%z = E1%z + E2%z

          else
             write (0, *) 'not compatible usage for add_cvector'
          end if

       else

          write(0, *) 'Error:add_cvector: vectors not same size'

       end if
    end if

  end subroutine add_cvector ! add_cvector


  !****************************************************************************
  ! add_cvector_f adds vectors stored as derived data type
  ! cvector with ; function version
  function add_cvector_f(E1, E2) result(E3)

    implicit none
    type (cvector), intent(in)               :: E1, E2
    type (cvector)                           :: E3

     if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for add_cvector_f'
       return
    endif

    ! In function version, appropriate data types need to be created
    Call create_cvector(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for add_cvector_f'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! add x,y,z-components
             E3%x = E1%x + E2%x
             E3%y = E1%y + E2%y
             E3%z = E1%z + E2%z

          else
             write (0, *) 'not compatible usage for add_cvector_f'
          end if

       else

          write(0, *) 'Error:add_cvector: vectors not same size'

       end if
    end if

    E3%temporary = .true.

  end function add_cvector_f ! add_cvector_f


  !****************************************************************************
  ! subtract_rvector subtracts vectors stored as derived data type rvector with
  ! ; subroutine version
  ! E3 can overwrite E1 and E2
  subroutine subtract_rvector(E1, E2, E3)

    implicit none
    type (rvector), intent(in)               :: E1, E2
    type (rvector), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for subtract_rvector'
       return
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for subtract_rvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! subtract x,y,z-components
             E3%x = E1%x - E2%x
             E3%y = E1%y - E2%y
             E3%z = E1%z - E2%z

          else
             write (0, *) 'not compatible usage for subtract_rvector'
          end if

       else

          write(0, *) 'Error: subtract_rvector: vectors not same size'

       end if
    end if

  end subroutine subtract_rvector ! subtract_rvector


  !****************************************************************************
  ! subtract_rvector_f subtracts vectors stored as derived data type
  ! rvector with ; function version
  function subtract_rvector_f(E1, E2) result(E3)

    implicit none
    type (rvector), intent(in)               :: E1, E2
    type (rvector)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for add_rvector_f'
       return
    endif

    ! In function version, appropriate data types need to be created
    Call create_rvector(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for add_rvector_f'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! subtract x,y,z-components
             E3%x = E1%x - E2%x
             E3%y = E1%y - E2%y
             E3%z = E1%z - E2%z

          else
             write (0, *) 'not compatible usage for subtract_rvector_f'
          end if

       else

          write(0, *) 'Error:subtract_rvector_f: vectors not same size'

       end if
    end if

    E3%temporary = .true.

  end function subtract_rvector_f ! subtract_rvector_f


  !****************************************************************************
  ! subtract_cvector subtracts vectors stored as derived data type
  ! cvector with ; subroutine version
  ! E3 can overwrite E1 and E2
  subroutine subtract_cvector(E1, E2, E3)

    implicit none
    type (cvector), intent(in)               :: E1, E2
    type (cvector), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for subtract_cvector'
       return
    endif

    ! check to see if LHS (E2) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS not allocated for subtract_cvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! subtract x,y,z-components
             E3%x = E1%x - E2%x
             E3%y = E1%y - E2%y
             E3%z = E1%z - E2%z

          else
             write (0, *) 'not compatible usage for subtract_cvector'
          end if

       else

          write(0, *) 'Error:subtract_cvector: vectors not same size'

       end if
    end if

  end subroutine subtract_cvector ! subtract_cvector


  !****************************************************************************
  ! subtract_cvector_f subtracts vectors stored as derived data type
  ! cvector with ; function version
  function subtract_cvector_f(E1, E2) result(E3)

    implicit none
    type (cvector), intent(in)               :: E1, E2
    type (cvector)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for subtract_cvector_f'
       return
    endif

    ! In function version, appropriate data types need to be created
    Call create_cvector(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for subtract_cvector_f'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! subtract x,y,z-components
             E3%x = E1%x - E2%x
             E3%y = E1%y - E2%y
             E3%z = E1%z - E2%z

          else
             write (0, *) 'not compatible usage for subtract_cvector_f'
          end if

       else

          write(0, *) 'Error: subtract_cvector_f: vectors not same size'

       end if
    end if

    E3%temporary = .true.

  end function subtract_cvector_f ! subtract_cvector_f


  !****************************************************************************
  ! diagMult_rvector multiplies two vectors E1, E2 stored as derived data
  ! type rvector pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagMult_rvector(E1, E2, E3)

    implicit none
    type (rvector), intent(in)               :: E1, E2
    type (rvector), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_rvector'
       return
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_rvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for x,y,z-components
             E3%x = E1%x * E2%x
             E3%y = E1%y * E2%y
             E3%z = E1%z * E2%z

          else
             write (0, *) 'not compatible usage for diagMult_rvector'
          end if

       else

          write(0, *) 'Error:diagMult_rvector: vectors not same size'

       end if
    end if

  end subroutine diagMult_rvector ! diagMult_rvector


  !****************************************************************************
  ! diagMult_rvector_f multiplies two vectors E1, E2 stored as derived
  ! data type rvector pointwise; function version
  function diagMult_rvector_f(E1, E2) result(E3)

    implicit none
    type (rvector), intent(in)               :: E1, E2
    type (rvector)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_rvector_f'
       return
    endif

    ! In function version, appropriate data types need to be created
    Call create_rvector(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_rvector_f'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for x,y,z-components
             E3%x = E1%x * E2%x
             E3%y = E1%y * E2%y
             E3%z = E1%z * E2%z

          else
             write (0, *) 'not compatible usage for diagMult_rvector_f'
          end if

       else

          write(0, *) 'Error:diagMult_rvector_f: vectors not same size'

       end if
    end if

    E3%temporary = .true.

  end function diagMult_rvector_f ! diagMult_rvector_f


  !****************************************************************************
  ! diagMult_cvector multiplies two vectors E1, E2 stored as derived data
  ! type cvector pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagMult_cvector(E1, E2, E3)

    implicit none
    type (cvector), intent(in)               :: E1, E2
    type (cvector), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_cvector'
       return
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diaMult_cvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for x,y,z-components
             E3%x = E1%x * E2%x
             E3%y = E1%y * E2%y
             E3%z = E1%z * E2%z

          else
             write (0, *) 'not compatible usage for diagMult_cvector'
          end if

       else

          write(0, *) 'Error:diagMult_cvector: vectors not same size'

       end if
    end if

  end subroutine diagMult_cvector ! diagMult_cvector


  !****************************************************************************
  ! diagMult_cvector_f multiplies two vectors E1, E2 stored as derived
  ! data  type cvector pointwise; function version
  function diagMult_cvector_f(E1, E2) result(E3)

    implicit none
    type (cvector), intent(in)               :: E1, E2
    type (cvector)                           :: E3

   if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_cvector_f'
       return
    endif

    ! In function version, appropriate data types need to be created
    Call create_cvector(E1%grid, E3, E1%gridType)
    ! check to see if RHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'RHS was not allocated for diagMult_cvector_f'
    else

       ! Check whether both vectors are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for x,y,z-components
             E3%x = E1%x * E2%x
             E3%y = E1%y * E2%y
             E3%z = E1%z * E2%z

          else
             write (0, *) 'not compatible usage for diagMult_cvector_f'
          end if

       else

          write(0, *) 'Error:diagMult_cvector_f: vectors not same size'

       end if
    end if

    E3%temporary = .true.

  end function diagMult_cvector_f ! diagMult_cvector_f


  !****************************************************************************
  ! diagMult_crvector multiplies complex vector E1 with scalar vector E2
  ! stored as derived type cvector pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagMult_crvector(E1, E2, E3)

    implicit none
    type (cvector), intent(in)               :: E1
    type (rvector), intent(in)               :: E2
    type (cvector), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_crvector'
       return
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_crvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for x,y,z-components
             E3%x = E1%x * E2%x
             E3%y = E1%y * E2%y
             E3%z = E1%z * E2%z

          else
             write (0, *) 'not compatible usage for diagMult_crvector'
          end if

       else

          write(0, *) 'Error:diagMult_crvector: vectors not same size'

       end if
    end if

  end subroutine diagMult_crvector ! diagMult_crvector


  !****************************************************************************
  ! diagMult_crvector_f multiplies complex vector E1 with real vector
  ! E2 stored as derived data type cvector pointwise; function version
  function diagMult_crvector_f(E1, E2) result(E3)

    implicit none
    type (cvector), intent(in)               :: E1
    type (rvector), intent(in)               :: E2
    type (cvector)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_crvector_f'
       return
    endif

    ! In function version, appropriate data types need to be created
    Call create_cvector(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_crvector_f'
    else

       ! Check whether both vectors are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for x,y,z-components
             E3%x = E1%x * E2%x
             E3%y = E1%y * E2%y
             E3%z = E1%z * E2%z

          else
             write (0, *) 'not compatible usage for diagMult_crvector_f'
          end if

       else

          write(0, *) 'Error:diagMult_Node_MixedCR_f: vectors not same size'

       end if
    end if

    E3%temporary = .true.

  end function diagMult_crvector_f ! diagMult_crvector_f


  !****************************************************************************
  ! diagMult_rcvector multiplies real vector E1 with complex vector E2
  ! stored as derived type cvector pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagMult_rcvector(E1, E2, E3)

    implicit none
    type (rvector), intent(in)               :: E1
    type (cvector), intent(in)               :: E2
    type (cvector), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_rcvector'
       return
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_rcvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for x,y,z-components
             E3%x = E1%x * E2%x
             E3%y = E1%y * E2%y
             E3%z = E1%z * E2%z

          else
             write (0, *) 'not compatible usage for diagMult_rcvector'
          end if

       else

          write(0, *) 'Error:diagMult_rcvector: vectors not same size'

       end if
    end if

  end subroutine diagMult_rcvector ! diagMult_rcvector


  !****************************************************************************
  ! diagMult_rcvector_f multiplies real vector E1 with complex vector E2
  ! stored as derived data type cvector pointwise; function version
  function diagMult_rcvector_f(E1, E2) result(E3)

    implicit none
    type (rvector), intent(in)               :: E1
    type (cvector), intent(in)               :: E2
    type (cvector)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_rcvector_f'
       return
    endif

    ! In function version, appropriate data types need to be created
    Call create_cvector(E1%grid, E3, E1%gridType)
    ! check to see if RHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'RHS was not allocated for diagMult_rcvector_f'
    else

       ! Check whether both vectors are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for x,y,z-components
             E3%x = E1%x * E2%x
             E3%y = E1%y * E2%y
             E3%z = E1%z * E2%z

          else
             write (0, *) 'not compatible usage for diagMult_rcvector_f'
          end if

       else

          write(0, *) 'Error:diagMult_rcvector_f: vectors not same size'

       end if
    end if

    E3%temporary = .true.

  end function diagMult_rcvector_f ! diagMult_rcvector_f

  !****************************************************************************
  ! diagDiv_rvector divides real vector E1 with real vector E2
  ! stored as derived type rvector pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagDiv_rvector(E1, E2, E3)

    implicit none
    type (rvector), intent(in)               :: E1
    type (rvector), intent(in)               :: E2
    type (rvector), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagDiv_rvector'
       return
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagDiv_rvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise division for x,y,z-components
             E3%x = E1%x / E2%x
             E3%y = E1%y / E2%y
             E3%z = E1%z / E2%z

          else
             write (0, *) 'not compatible usage for diagDiv_rvector'
          end if

       else

          write(0, *) 'Error:diagDiv_rvector: vectors not same size'

       end if
    end if

  end subroutine diagDiv_rvector ! diagDiv_rvector

  !****************************************************************************
  ! diagDiv_crvector divides complex vector E1 with real vector E2
  ! stored as derived type cvector pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagDiv_crvector(E1, E2, E3)

    implicit none
    type (cvector), intent(in)               :: E1
    type (rvector), intent(in)               :: E2
    type (cvector), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagDiv_crvector'
       return
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagDiv_crvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise division for x,y,z-components
             E3%x = E1%x / E2%x
             E3%y = E1%y / E2%y
             E3%z = E1%z / E2%z

          else
             write (0, *) 'not compatible usage for diagDiv_crvector'
          end if

       else

          write(0, *) 'Error:diagDiv_crvector: vectors not same size'

       end if
    end if

  end subroutine diagDiv_crvector ! diagDiv_crvector


  !****************************************************************************
  ! diagDiv_crvector_f divides complex vector E1 with real vector
  ! E2 stored as derived data type cvector pointwise; function version
  function diagDiv_crvector_f(E1, E2) result(E3)

    implicit none
    type (cvector), intent(in)               :: E1
    type (rvector), intent(in)               :: E2
    type (cvector)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagDiv_crvector_f'
       return
    endif

    ! In function version, appropriate data types need to be created
    Call create_cvector(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagDiv_crvector_f'
    else

       ! Check whether both vectors are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise division for x,y,z-components
             E3%x = E1%x / E2%x
             E3%y = E1%y / E2%y
             E3%z = E1%z / E2%z

          else
             write (0, *) 'not compatible usage for diagDicrvector_f'
          end if

       else

          write(0, *) 'Error:diagDicrvector_f: vectors not same size'

       end if
    end if

    E3%temporary = .true.

  end function diagDiv_crvector_f ! diagDiv_crvector_f


  !****************************************************************************
  ! diagDiv_rcvector divides real vector E1 with complex vector E2
  ! stored as derived type cvector pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagDiv_rcvector(E1, E2, E3)

    implicit none
    type (rvector), intent(in)               :: E1
    type (cvector), intent(in)               :: E2
    type (cvector), intent(inout)            :: E3

       if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagDiv_rcvector'
       return
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagDiv_rcvector'
    else

       ! Check whether all the vector nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise divides for x,y,z-components
             E3%x = E2%x / E1%x
             E3%y = E2%y / E1%y
             E3%z = E2%z / E1%z

          else
             write (0, *) 'not compatible usage for diagDiv_rcvector'
          end if

       else

          write(0, *) 'Error:diagDiv_rcvector: vectors not same size'

       end if
    end if

  end subroutine diagDiv_rcvector ! diagDiv_rcvector


  !****************************************************************************
  ! diagDiv_rcvector_f divides real vector E1 with complex vector E2
  ! stored as derived data type cvector pointwise; function version
  function diagDiv_rcvector_f(E1, E2) result(E3)

    implicit none
    type (rvector), intent(in)               :: E1
    type (cvector), intent(in)               :: E2
    type (cvector)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagDiv_rcvector_f'
       return
    endif

    ! In function version, appropriate data types need to be created
    Call create_cvector(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagDiv_rcvector_f'
    else

       ! Check whether both vectors are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise divides for x,y,z-components
             E3%x = E2%x / E1%x
             E3%y = E2%y / E1%y
             E3%z = E2%z / E1%z

          else
             write (0, *) 'not compatible usage for diagDiv_rcvector_f'
          end if

       else

          write(0, *) 'Error:diagDiv_rcvector_f: vectors not same size'

       end if
    end if

    E3%temporary = .true.

  end function diagDiv_rcvector_f ! diagDiv_rcvector_f


  ! ***************************************************************************
  ! dotProd_rvector computes dot product of two vectors stored
  ! as derived data type rvector, returning a real number
  function dotProd_rvector_f(E1, E2) result(r)

    implicit none
    type (rvector), intent(in)   :: E1, E2
	type (rvector)				 :: E3
    real(kind=prec)		         :: r
	integer						 :: nx,ny,nz

    r = R_ZERO

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'input vectors not allocated yet for dotProd_rvector_f'
       return
    endif

	E3 = E1
	nx = E3%nx
	ny = E3%ny
	nz = E3%nz

	if (E3%grid%coords == Spherical) then
	  ! needs special treatment
	  if (E3%gridType == EDGE) then
		! delete duplicate edges
		E3%y(nx+1,:,:) = R_ZERO
		E3%z(nx+1,:,:) = R_ZERO
		E3%z(2:nx,1,:) = R_ZERO
		E3%z(2:nx,ny+1,:) = R_ZERO
	  else if (E3%gridType == FACE) then
		! delete duplicate faces
		E3%x(nx+1,:,:) = R_ZERO
	  else
		write(0,*) 'unknown gridType ',trim(E3%gridType),' in dotProd_rvector_f'
	  end if
	else if (E3%grid%coords == Cartesian) then
	  ! do nothing
	else
	  write(0,*) 'unknown coordinate system ',trim(E3%grid%coords),' in dotProd_rvector_f'
	  return
	end if

    ! Check whether both input vectors are of the same size
    if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

       if ((E1%gridType == E2%gridType)) then

          r = r + sum(E3%x * E2%x)
          r = r + sum(E3%y * E2%y)
          r = r + sum(E3%z * E2%z)

       else
          write (0,*) 'not compatible input vectors in dotProd_rvector_f'
       end if

    else

       write(0,*) 'vectors not the same size in dotProd_rvector_f'

    end if

    call deall_rvector(E3)

  end function dotProd_rvector_f  ! dotProd_rvector_f


  ! ***************************************************************************
  ! dotProd_cvector computes dot product of two vectors stored
  ! as derived data type cvector, returning a complex number
  function dotProd_cvector_f(E1, E2) result(c)

    implicit none
    type (cvector), intent(in)       :: E1, E2
	type (cvector)					 :: E3
    complex(kind=prec)		         :: c
	integer						     :: nx,ny,nz

    c = R_ZERO

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for dotProdCC'
       return
    endif

 	E3 = E1
	nx = E3%nx
	ny = E3%ny
	nz = E3%nz

	if (E3%grid%coords == Spherical) then
	  ! needs special treatment
	  if (E3%gridType == EDGE) then
		! delete duplicate edges
		E3%y(nx+1,:,:) = C_ZERO
		E3%z(nx+1,:,:) = C_ZERO
		E3%z(2:nx,1,:) = C_ZERO
		E3%z(2:nx,ny+1,:) = C_ZERO
	  else if (E3%gridType == FACE) then
		! delete duplicate faces
		E3%x(nx+1,:,:) = C_ZERO
	  else
		write(0,*) 'unknown gridType ',trim(E3%gridType),' in dotProd_cvector_f'
	  end if
	else if (E3%grid%coords == Cartesian) then
	  ! do nothing
	else
	  write(0,*) 'unknown coordinate system ',trim(E3%grid%coords),' in dotProd_cvector_f'
	  return
	end if

   ! Check whether both input vectors are of the same size
    if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

       if ((E1%gridType == E2%gridType)) then

          c = c + sum(conjg(E3%x) * E2%x)
          c = c + sum(conjg(E3%y) * E2%y)
          c = c + sum(conjg(E3%z) * E2%z)

       else
          write (0, *) 'not compatible input vectors in dotProd_cvector_f'
       end if

    else

       write(0, *) 'vectors not the same size in dotProd_cvector_f'

    end if

    call deall_cvector(E3)

  end function dotProd_cvector_f ! dotProd_cvector


  ! ***************************************************************************
  ! dotProd_noConj_cvector computes dot product of two vectors stored
  ! as derived data type cvector, returning a complex number
  function dotProd_noConj_cvector_f(E1, E2) result(c)

    implicit none
    type (cvector), intent(in)       :: E1, E2
	type (cvector)					 :: E3
    complex(kind=prec)		         :: c
	integer						     :: nx,ny,nz

    c = R_ZERO

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for dotProdCC'
       return
    endif

 	E3 = E1
	nx = E3%nx
	ny = E3%ny
	nz = E3%nz

	if (E3%grid%coords == Spherical) then
	  ! needs special treatment
	  if (E3%gridType == EDGE) then
		! delete duplicate edges
		E3%y(nx+1,:,:) = C_ZERO
		E3%z(nx+1,:,:) = C_ZERO
		E3%z(2:nx,1,:) = C_ZERO
		E3%z(2:nx,ny+1,:) = C_ZERO
	  else if (E3%gridType == FACE) then
		! delete duplicate faces
		E3%x(nx+1,:,:) = C_ZERO
	  else
		write(0,*) 'unknown gridType ',trim(E3%gridType),' in dotProd_noConj_cvector_f'
	  end if
	else if (E3%grid%coords == Cartesian) then
	  ! do nothing
	else
	  write(0,*) 'unknown coordinate system ',trim(E3%grid%coords),' in dotProd_noConj_cvector_f'
	  return
	end if

   ! Check whether both input vectors are of the same size
    if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

       if ((E1%gridType == E2%gridType)) then

          c = c + sum(E3%x * E2%x)
          c = c + sum(E3%y * E2%y)
          c = c + sum(E3%z * E2%z)

       else
          write (0, *) 'not compatible input vectors in dotProd_noConj_cvector_f'
       end if

    else

       write(0, *) 'vectors not the same size in dotProd_noConj_cvector_f'

    end if

    call deall_cvector(E3)

  end function dotProd_noConj_cvector_f ! dotProd__noConj_cvector_f


  !****************************************************************************
  ! linComb_cvector computes linear combination of two vectors
  ! stored as derived data type cvector; subroutine, not a function
  ! both input vectors must have the same dimension
  subroutine linComb_cvector(inc1, E1, inc2, E2, E3)

    implicit none
    !   input vectors
    type (cvector), intent(in)             :: E1, E2
    !  input complex scalars
    complex (kind=prec), intent(in)           :: inc1, inc2
    type (cvector), intent(inout)          :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       call warning('inputs not allocated yet for linComb_cvector')
       return
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       call warning('output has to be allocated before call to linComb_cvector')

    elseif (compare(E1,E2) .and. compare(E1,E3)) then
        ! form linear combinatoin
        E3%x = inc1*E1%x + inc2*E2%x
        E3%y = inc1*E1%y + inc2*E2%y
        E3%z = inc1*E1%z + inc2*E2%z

    else
        call warning('not compatible usage for linComb_cvector')
    end if

  end subroutine linComb_cvector ! linComb_cvector

  !****************************************************************************
  ! linComb_rvector computes linear combination of two vectors
  ! stored as derived data type cvector; subroutine, not a function
  ! both input vectors must have the same dimension
  subroutine linComb_rvector(inc1, E1, inc2, E2, E3)

    implicit none
    !   input vectors
    type (rvector), intent(in)             :: E1, E2
    !  input real scalars
    real (kind=prec), intent(in)           :: inc1, inc2
    type (rvector), intent(inout)          :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       call warning('inputs not allocated yet for linComb_rvector')
       return
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       call warning('output has to be allocated before call to linComb_rvector')
       return

    elseif (compare(E1,E2) .and. compare(E1,E3)) then
        ! form linear combinatoin
        E3%x = inc1*E1%x + inc2*E2%x
        E3%y = inc1*E1%y + inc2*E2%y
        E3%z = inc1*E1%z + inc2*E2%z

    else
        call warning('not compatible usage for linComb_rvector')
        return
    end if

  end subroutine linComb_rvector ! linComb_rvector

  !****************************************************************************
  ! scMultadd_cvector multiplies vector E1 stored as derived data type
  ! cvector with a complex scalar c, adding result to output vector E2
  subroutine scMultAdd_cvector(c, E1, E2)

    implicit none
    complex(kind=prec), intent(in)                   :: c
    ! a complex scalar to be multiplied with
    type (cvector), intent(in)                       :: E1
    type (cvector)                                   :: E2

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMultAdd_cvector'
       return
    endif

    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated for scMultAdd_cvector'
    else

       ! Check whether both vectors are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if ((E1%gridType == E2%gridType)) then

             ! complex scalar multiplication for x,y,z-components
             E2%x = E2%x + E1%x * c
             E2%y = E2%y + E1%y * c
             E2%z = E2%z + E1%z * c

          else
             write (0, *) 'not compatible usage for scMultAdd_cvector'
          end if

       else

          write(0, *) 'Error:scMultAdd_cvector: vectors not same size'

       end if
    end if

  end subroutine scMultAdd_cvector ! scMultAdd_cvector


  ! ***************************************************************************
  function conjg_cvector_f(E1) result (E2)
	! conjg_cvector_f computes a conjugate of a derived data type cvector
	! A.K.
    implicit none
    type (cvector), intent(in)            :: E1
    type (cvector)                        :: E2

    integer                               :: status

    ! check to see if RHS (E1) is active (allocated)
    if(.not.E1%allocated) then
       write(0,*) 'input not allocated yet for conjg_cvector_f'
    else

       if((E2%nx == E1%nx).and.(E2%ny == E1%ny).and.(E2%nz == E1%nz)) then

          if  (E1%gridType == E2%gridType) then

             ! just conjugate components
             E2%x = conjg(E1%x)
             E2%y = conjg(E1%y)
             E2%z = conjg(E1%z)
             E2%gridType = E1%gridType

          else
             write (0, *) 'not compatible usage for conjg_cvector_f'
          end if

       else

          if(E2%allocated) then
             ! first deallocate memory for x,y,z
             deallocate(E2%x, E2%y, E2%z,STAT=status)
          end if

          !  then allocate E2 as correct size ...
          Call create_cvector(E1%grid, E2, E1%gridType)
          !   .... and conjugate E1
          E2%x = conjg(E1%x)
          E2%y = conjg(E1%y)
          E2%z = conjg(E1%z)
          E2%gridType = E1%gridType

       end if

    end if

    E2%temporary = .true.

  end function conjg_cvector_f  ! conjg_cvector_f


  ! ***************************************************************************
  function cmplx_rvector_f(E1, E2) result (E3)
  ! inputs two real vectors, merges them as real1 + imag(real2), a complex
  ! vector
    implicit none
    type (rvector), intent(in)            :: E1
    type (rvector), intent(in)            :: E2
    type (cvector)                        :: E3

    integer                               :: status

    ! check to see if RHS (E1 and E2) are active (allocated)
    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for cmplx_rvector_f'
    else

       if((E3%nx == E1%nx).and.(E3%ny == E1%ny).and.(E3%nz == E1% nz).and.&
        (E3%nx == E2%nx).and.(E3%ny == E2%ny).and.(E3%nz == E2%nz))  then

          if  ((E1%gridType == E2%gridType).and.(E1%gridType == E3% gridType)) then

             ! create a complex pair
             E3%x = cmplx(E1%x, E2%x, prec)
             E3%y = cmplx(E1%y, E2%y, prec)
             E3%z = cmplx(E1%z, E2%z, prec)
             E3%gridType = E1%gridType

          else
             write (0, *) 'not compatible usage for cmplx_rvector_f'
          end if

       else

          if(E3%allocated) then
             ! first deallocate memory for x,y,z
             deallocate(E3%x, E3%y, E3%z,STAT=status)
          end if

          !  then allocate E3 as correct size ...
          Call create_cvector(E1%grid, E3, E1%gridType)
          !   .... and create a complex pair
          E3%x = cmplx(E1%x, E2%x, prec)
          E3%y = cmplx(E1%y, E2%y, prec)
          E3%z = cmplx(E1%z, E2%z, prec)
          E3%gridType = E1%gridType

       end if

    end if

    E3%temporary = .true.

  end function cmplx_rvector_f  ! cmplx_rvector_f


  ! ***************************************************************************
  ! real_cvector_f copies the real part of the derived data type cvector variable;
  ! to produce a derived data type rvector
  function real_cvector_f(E1) result (E2)

    implicit none
    type (cvector), intent(in)            :: E1
    type (rvector)                        :: E2

    integer                               :: status

    ! check to see if RHS (E1) is active (allocated)
    if(.not.E1%allocated) then
       write(0,*) 'input not allocated yet for real_cvector_f'
    else

      ! we know nothing about E2 ... deallocate just in case
      Call deall_rvector(E2)
      !  then allocate E2 as correct size ...
      Call create_rvector(E1%grid, E2, E1%gridType)
      !   .... and copy E1
      E2%x = real(E1%x)
      E2%y = real(E1%y)
      E2%z = real(E1%z)
      E2%gridType = E1%gridType

    end if

    E2%temporary = .true.

  end function real_cvector_f  ! real_cvector_f


  ! ***************************************************************************
  ! imag_cvector_f copies the imag part of the derived data type cvector variable;
  ! to produce a derived data type rvector
  function imag_cvector_f(E1) result (E2)

    implicit none
    type (cvector), intent(in)            :: E1
    type (rvector)                        :: E2

    integer                               :: status

    ! check to see if RHS (E1) is active (allocated)
    if(.not.E1%allocated) then
       write(0,*) 'input not allocated yet for imag_cvector_f'
    else

      ! we know nothing about E2 ... deallocate just in case
      Call deall_rvector(E2)
      !  then allocate E2 as correct size ...
      Call create_rvector(E1%grid, E2, E1%gridType)
      !   .... and copy E1
      E2%x = aimag(E1%x)
      E2%y = aimag(E1%y)
      E2%z = aimag(E1%z)
      E2%gridType = E1%gridType

    end if

    E2%temporary = .true.

  end function imag_cvector_f  ! imag_cvector_f

end module sg_vector ! sg_vector
