! *****************************************************************************
module sg_scalar
  ! This module creates data types for scalar fields defined on
  ! scalar (center or corner) nodes of a staggered Cartesian grid, along with
  ! basic algebraic (scalar) operations. Such basic operations are: allocation,
  ! deallocation, initialization, copying, and algebraic operations (linear
  ! combinations, scalar products, dot products).
  ! Belongs to SG_Basics class: staggered cartesian grid, data
  ! types defined on this grid, and operations defined on these data types. Not
  ! specific to EM problem, no dependency on outside (from other classes) modules.

  use math_constants		! math/ physics constants
  use utilities             ! for error and warning messages
  use griddef

  implicit none

  ! Generic interfaces are done through subroutines
  ! creates scalar (center or corner) nodes
  INTERFACE create
     module procedure create_rscalar
     module procedure create_cscalar
     module procedure create_iscalar
  END INTERFACE

  ! deallocates the scalar (center or corner) nodes
  INTERFACE deall
     module procedure deall_rscalar
     module procedure deall_cscalar
     module procedure deall_iscalar
  END INTERFACE

  ! reads the scalar (center or corner) nodes from fid
  INTERFACE read
     module procedure read_rscalar
     module procedure read_cscalar
     module procedure read_iscalar
  END INTERFACE

  ! writes the scalar (center or corner) nodes to fid
  INTERFACE write
     module procedure write_rscalar
     module procedure write_cscalar
     module procedure write_iscalar
  END INTERFACE

  ! a real scalar multiplies the scalar (center or corner) nodes
  INTERFACE scMult
     module procedure scMult_rscalar
     module procedure scMult_cscalar
  END INTERFACE

  INTERFACE scMultAdd
     module procedure scMultAdd_rscalar
     module procedure scMultAdd_cscalar
  END INTERFACE

  INTERFACE linComb
     MODULE PROCEDURE linComb_rscalar
     MODULE PROCEDURE linComb_cscalar
  END INTERFACE

  ! adds the scalar (center or corner) nodes
  INTERFACE add
     module procedure add_rscalar
     module procedure add_cscalar
  END INTERFACE

  ! subtracts the scalar (center or corner) nodes
  INTERFACE subtract
     module procedure subtract_rscalar
     module procedure subtract_cscalar
  END INTERFACE

  ! pointwise scalar multiplication of scalar (center or corner) nodes and
  ! pointwise complex-real (mixed) multiplication of scalar (center or corner)
  ! nodes
  ! Both are scalar data types
  INTERFACE diagMult
     module procedure diagMult_rscalar
     module procedure diagMult_cscalar
     module procedure diagMult_rcscalar
     module procedure diagMult_crscalar
  END INTERFACE

  ! zeros the scalar (center or corner) nodes
  INTERFACE zero
     module procedure zero_rscalar
     module procedure zero_cscalar
  END INTERFACE

  ! computes the dot product
  INTERFACE dotProd
     MODULE PROCEDURE dotProd_rscalar_f
     MODULE PROCEDURE dotProd_cscalar_f
  END INTERFACE


  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE copy_cscalar
     MODULE PROCEDURE copy_rscalar
     MODULE PROCEDURE copy_iscalar
  END INTERFACE

  ! Interface operators are done through functions
  INTERFACE OPERATOR (+)
     MODULE PROCEDURE add_rscalar_f
     MODULE PROCEDURE add_cscalar_f
  END INTERFACE

  INTERFACE OPERATOR (-)
     MODULE PROCEDURE subtract_rscalar_f
     MODULE PROCEDURE subtract_cscalar_f
  END INTERFACE

  INTERFACE OPERATOR (*)
     MODULE PROCEDURE scMult_rscalar_f
     MODULE PROCEDURE scMult_cscalar_f
     MODULE PROCEDURE diagMult_rscalar_f
     MODULE PROCEDURE diagMult_cscalar_f
     MODULE PROCEDURE diagMult_rcscalar_f
     MODULE PROCEDURE diagMult_crscalar_f
  END INTERFACE


  public		::   create_rscalar,  create_cscalar, create_iscalar, &
       deall_rscalar, deall_cscalar, deall_iscalar, &
       read_rscalar, read_cscalar, read_iscalar, &
       write_rscalar, write_cscalar, write_iscalar, &
       copy_rscalar, copy_cscalar, &
       zero_rscalar, zero_cscalar, &
       scMult_rscalar, scMult_cscalar, &
       scMult_rscalar_f, scMult_cscalar_f, &
       add_rscalar, add_cscalar, &
       add_rscalar_f, add_cscalar_f, &
       subtract_rscalar, subtract_cscalar, &
       subtract_rscalar_f, subtract_cscalar_f, &
       diagMult_rscalar, diagMult_cscalar, &
       diagMult_rscalar_f, diagMult_cscalar_f, &
       diagMult_rcscalar, diagMult_crscalar, &
       diagMult_rcscalar_f, diagMult_crscalar_f, &
       dotProd_rscalar_f, dotProd_cscalar_f, &
       linComb_rscalar, linComb_cscalar, scMultAdd_cscalar

  ! ***************************************************************************
  ! type cscalar defines scalar for either edge or face in a staggered grid as
  ! a complex field
  type :: cscalar

     ! store the intention of the use in a character string defined
     ! in GridDef as a parameter: CENTER, CORNER, CELL_EARTH
     character (len=80)	                               :: gridType

     ! Note that the arrays are defined through dynamic memory allocation
     complex(kind=prec), pointer, dimension(:,:,:)    :: v

     ! Grid Dimensions:
     ! nx is grid dimension (number of cells) in the x-direction
     ! ny is grid dimension (number of cells) in the y-direction
     ! nz is grid dimension (number of cells) in the z-direction:
     integer                                           :: nx = 0, ny = 0, nz = 0

     ! allocated:  .true.  v array has been allocated
     logical		                               :: allocated = .false.

     ! temporary:  .true. for function outputs only; necessary to avoid memory leaks
     ! (probably will not be needed in the future when compilers will support
     ! ISO/IEC 15581 - the "allocatable array extension")
     logical										:: temporary = .false.

     ! pointer to parent grid
     type (grid_t), pointer                              :: grid

  end type cscalar

  ! ***************************************************************************
  ! type rscalar defines scalar for either edge or face in a staggered grid as
  ! a real field
  type :: rscalar

     ! store the intention of the use in a character string defined
     ! in GridDef as a parameter: CENTER, CORNER, CELL_EARTH
     character (len=80)	                           	:: gridType

     ! Typical usage:  conductivity averaged on centers of
     ! staggered grid
     ! v: dimension Nx, Ny, Nz
     ! Note that the arrays are defined through dynamic memory allocation
     real(kind=prec), pointer, dimension(:,:,:)        :: v

     ! Grid Dimensions:
     ! nx is grid dimension (number of cells) in the x-direction
     ! ny is grid dimension (number of cells) in the y-direction
     ! nz is grid dimension (number of cells) in the z-direction:
     integer                                            :: nx = 0, ny = 0, nz =0

     ! allocated:  .true.  v array has been allocated
     logical		                                :: allocated = .false.

     ! temporary:  .true. for function outputs only; necessary to avoid memory leaks
     ! (probably will not be needed in the future when compilers will support
     ! ISO/IEC 15581 - the "allocatable array extension")
     logical										:: temporary = .false.

     ! pointer to parent grid
     type (grid_t), pointer                               :: grid

  end type rscalar

  ! ***************************************************************************
  ! type iscalar defines scalar for either edge or face in a staggered grid as
  ! an integer field; useful for mask arrays, so only create and deall needed
  type :: iscalar

     ! store the intention of the use in a character string defined
     ! in GridDef as a parameter: CENTER, CORNER, CELL_EARTH
     character (len=80)	                           	:: gridType

     ! Typical usage:  conductivity averaged on centers of
     ! staggered grid
     ! v: dimension Nx, Ny, Nz
     ! Note that the arrays are defined through dynamic memory allocation
     integer(kind=prec), pointer, dimension(:,:,:)        :: v

     ! Grid Dimensions:
     ! nx is grid dimension (number of cells) in the x-direction
     ! ny is grid dimension (number of cells) in the y-direction
     ! nz is grid dimension (number of cells) in the z-direction:
     integer                                            :: nx = 0, ny = 0, nz =0

     ! allocated:  .true.  v array has been allocated
     logical		                                :: allocated = .false.

     ! temporary:  .true. for function outputs only; necessary to avoid memory leaks
     ! (probably will not be needed in the future when compilers will support
     ! ISO/IEC 15581 - the "allocatable array extension")
     logical										:: temporary = .false.

     ! pointer to parent grid
     type (grid_t), pointer                               :: grid

  end type iscalar

Contains
  ! CREATE GRID_scalar
  ! * subroutine create_rscalar(igrid, E, gridType)
  ! * subroutine create_cscalar(igrid, E, gridType)
  ! * subroutine create_iscalar(igrid, E, gridType)

  ! DEALLOCATE GRID_scalar
  ! * subroutine deall_rscalar(E)
  ! * subroutine deall_cscalar(E)
  ! * subroutine deall_iscalar(E)

  ! READ GRID_scalar
  ! * subroutine read_rscalar(fid, E)
  ! * subroutine read_cscalar(fid, E)
  ! * subroutine read_iscalar(fid, E)

  ! WRITE GRID_scalar
  ! * subroutine write_rscalar(fid, E)
  ! * subroutine write_cscalar(fid, E)
  ! * subroutine write_iscalar(fid, E)

  ! COPY GRID_scalar:  (=)
  ! * subroutine copy_rscalar(E2,E1)
  ! * subroutine copy_cscalar(E2,E1)

  ! ZERO GRID_scalar
  ! * subroutine zero_rscalar(E)
  ! * subroutine zero_cscalar(E)

  ! SCALAR MULTIPLICATION:
  ! * subroutine scMult_rscalar(c, E1, E2)
  ! * subroutine scMult_cscalar(c, E1, E2)

  ! SCALAR MULTIPLICATION: FUNCTION VERSION (c*E)
  ! * function scMult_rscalar_f(c, E1) result(E2)
  ! * function scMult_cscalar_f(c, E1) result(E2)

  ! scalar SUM:
  ! * subroutine add_rscalar(E1, E2, E3)
  ! * subroutine add_cscalar(E1, E2, E3)

  ! scalar SUM: FUNCTION VERSION (E3 = E1+E2)
  ! * function add_rscalar_f(E1, E2) result(E3)
  ! * function add_cscalar_f(E1, E2) result(E3)

  ! scalar SUBTRACT:
  ! * subroutine subtract_rscalar(E1, E2, E3)
  ! * subroutine subtract_cscalar(E1, E2, E3)

  ! scalar SUBTRACT: FUNCTION VERSION (E3 = E1-E2)
  ! * function subtract_rscalar_f(E1, E2) result(E3)
  ! * function subtract_cscalar_f(E1, E2) result(E3)

  ! POINTWISE MULTIPLICATION OF scalars:
  ! * subroutine diagMult_rscalar(E1, E2, E3)
  ! * subroutine diagMult_cscalar(E1, E2, E3)

  ! POINTWISE MULTIPLICATION OF scalars: FUNCTION VERSION (c*E)
  ! * function diagMult_rscalar_f(E1, E2) result(E3)
  ! * function diagMult_cscalar_f(E1, E2) result(E3)

  ! POINTWISE MULTIPLICATION OF scalars and scalars:
  ! one combination is real and complex scalars and another is complex
  ! and real scalars as the sequence for inputs
  ! * subroutine diagMult_crscalar(E1, E2, E3)
  ! * subroutine diagMult_rcscalar(E1, E2, E3)

  ! POINTWISE MULTIPLICATION OF scalarS and SCALARS: FUNCTION VERSION (c*E)
  ! one combination is real and complex scalars and another is complex
  ! and real scalars as the sequence for inputs
  ! * function diagMult_crscalar_f(E1, E2) result(E3)
  ! * function diagMult_rcscalar_f(E1, E2) result(E3)

  ! scalar DOT PRODUCT
  ! * function dotProd_rscalar(E1, E2) result(r)
  ! * function dotProd_cscalar(E1, E2) result(c)

  ! COMPLEX LINEAR COMBINATION ... formally redundent with above
  ! only doing ones that we are sure to want
  ! * subroutine linCom_cscalar(inc1, E1, inc2, E2, E3)
  ! * subroutine scMultAdd_S_node(c, E1, E2)  (E2 = E2+c*E1)

  ! The algebraic routines expect all input and output
  ! variables to be of the correct type, already allocated,
  ! and of the correct size.


  !****************************************************************************
  ! create_rscalar creates variable of derived type rscalar,
  ! using grid definition in structure "grid" ;
  ! allocates memory in v component array
  ! gridType is a character string to describe intended usage
  subroutine create_rscalar(igrid, E, gridType)

    implicit none
    type(grid_t), target, intent(in)    :: igrid
    ! the grid for which an scalar (center or corner) node field is being
    ! initialized
    type (rscalar), intent(inout)      :: E

    integer                            :: status,nx,ny,nz,nzEarth

    character (len=80)                 :: gridType

    if(E%allocated) then
       ! first deallocate memory for v
       deallocate(E%v, STAT=status)
    end if

    ! Set pointer
    E%grid => igrid

    ! Grid dimensions
    nx = igrid%nx
    ny = igrid%ny
    nz = igrid%nz
    nzEarth = igrid%nzEarth
    E%nx = nx
    E%ny = ny
    E%nz = nz
    E%gridType = gridType

    ! allocate memory for v
    ! E%allocated will be true if all allocations succeed
    if (E%gridType == CENTER) then
       allocate(E%v(nx,ny,nz), STAT=status)
    else if (E%gridType == CORNER) then
       allocate(E%v(nx+1,ny+1,nz+1), STAT=status)
    else if(E%gridType == CELL_EARTH) then
       E%nz = nzEarth
       allocate(E%v(nx,ny,nzEarth), STAT=status)
    else
       write (0, *) 'gridType == ',trim(E%gridType),' undefined in create_rscalar'
    end if
    E%allocated = status .EQ. 0

    if (E%allocated) then
       E%v = R_ZERO
    end if

  end subroutine create_rscalar

  ! ***************************************************************************
  ! create_cscalar creates variable of derived type cscalar,
  ! using grid definition in structure "grid" ;
  ! allocates memory in v component array
  ! gridType is a character string to describe intended usage
  subroutine create_cscalar(igrid, E, gridType)

    implicit none
    type (grid_t), target, intent(in)     :: igrid
    ! the grid for which an scalar (center or corner) node field is being
    ! initialized
    type (cscalar), intent(inout)       :: E

    integer                             :: status,nx,ny,nz

    character (len=80)                  :: gridType

    if(E%allocated) then
       ! first deallocate memory for v
       deallocate(E%v,STAT=status)
    end if

    ! Set pointer
    E%grid => igrid

    ! Grid dimensions
    nx = igrid%nx
    ny = igrid%ny
    nz = igrid%nz
    E%nx = nx
    E%ny = ny
    E%nz = nz
    E%gridType = gridType

    ! allocate memory for v ;
    ! E%allocated will be true if all allocations succeed
    E%allocated = .true.
    if (E%gridType == CENTER) then
       allocate(E%v(nx,ny,nz), STAT=status)
    else if (E%gridType == CORNER) then
       allocate(E%v(nx+1,ny+1,nz+1), STAT=status)
    else
       write (0, *) 'gridType == ',trim(E%gridType),' undefined in create_cscalar'
    end if
    E%allocated = E%allocated .and. (status .EQ. 0)

    if (E%allocated) then
       E%v = C_ZERO
    end if
    ! print *, 'E%allocated', E%allocated

  end subroutine create_cscalar

  !****************************************************************************
  ! create_iscalar creates variable of derived type iscalar,
  ! using grid definition in structure "grid" ;
  ! allocates memory in v component array
  ! gridType is a character string to describe intended usage; but iscalar
  ! is only useful for gridType == CELL_EARTH
  subroutine create_iscalar(igrid, E, gridType)

    implicit none
    type (grid_t), target, intent(in)    :: igrid
    type (iscalar),  intent(inout)         :: E

    integer                            :: status,nx,ny,nz,nzEarth

    character (len=80)                 :: gridType

    if(E%allocated) then
       ! first deallocate memory for v
       deallocate(E%v, STAT=status)
    end if

    ! Set pointer
    E%grid => igrid

    ! Grid dimensions
    nx = igrid%nx
    ny = igrid%ny
    nz = igrid%nz
    nzEarth = igrid%nzEarth
    E%nx = nx
    E%ny = ny
    E%nz = nz
    E%gridType = gridType

    ! allocate memory for v
    ! E%allocated will be true if all allocations succeed
    if(E%gridType == CELL_EARTH) then
       E%nz = nzEarth
       allocate(E%v(nx,ny,nzEarth), STAT=status)
    else
       write (0, *) 'gridType == ',trim(E%gridType),' undefined in create_iscalar'
    end if
    E%allocated = status .EQ. 0

    if (E%allocated) then
       E%v = 0
    end if

  end subroutine create_iscalar

  !****************************************************************************
  ! deall_rscalar destoys variable of derived type rscalar,
  ! deallocating memory
  subroutine deall_rscalar(E)

    implicit none
    type (rscalar)  :: E
    integer	    :: status

    if(E%allocated) then
       ! deallocate memory for v
       deallocate(E%v,STAT=status)
    end if

    E%nx = 0
    E%ny = 0
    E%nz = 0
    E%gridType = ''
    E%allocated = .false.

  end subroutine deall_rscalar

  !****************************************************************************
  ! deall_cscalar destoys variable of derived type cscalar,
  ! deallocating memory
  subroutine deall_cscalar(E)

    implicit none
    type (cscalar)  :: E
    integer	    :: status

    ! deallocate memory for v
    if(E%allocated) then
       ! deallocate memory for v
       deallocate(E%v,STAT=status)
    end if

    E%nx = 0
    E%ny = 0
    E%nz = 0
    E%gridType = ''
    E%allocated = .false.

  end subroutine deall_cscalar

  !****************************************************************************
  ! deall_iscalar destoys variable of derived type iscalar,
  ! deallocating memory
  subroutine deall_iscalar(E)

    implicit none
    type (iscalar)  :: E
    integer	    :: status

    ! deallocate memory for v
    if(E%allocated) then
       ! deallocate memory for v
       deallocate(E%v,STAT=status)
    end if

    E%nx = 0
    E%ny = 0
    E%nz = 0
    E%gridType = ''
    E%allocated = .false.

  end subroutine deall_iscalar

  !****************************************************************************
  ! read_rscalar reads in an rscalar in a simple format; rscalar has
  ! to exist and be allocated before calling this routine, and the file unit
  ! must already be available for reading; also need to be careful that
  ! dimensions in the input file match those of the rscalar;
  ! assuming the "intuitive geographic" file convention by reading
  ! x from the bottom up (x points North), y left to right (y points East)
  ! optional file type may be 'ascii' or 'binary', shortcuts allowed;
  ! defaults to 'ascii'. use binary for better accuracy.
  subroutine read_rscalar(fid, E, ftype)

      integer,        intent(in)		:: fid
      type (rscalar), intent(inout)		:: E
      character(*), optional, intent(in):: ftype

      !  local variables
      integer 		                    :: Nx, Ny, Nz
      integer                           :: i, j, k, k1, k2, istat
      real (kind(E%v)), allocatable     :: temp(:)
      logical                           :: ok, hasname, binary
      character(80)                     :: fname, isbinary

      if(.not. E%allocated) then
         write(0, *) 'rscalar must be allocated before call to read_rscalar'
         stop
      endif

      if (.not. present(ftype)) then
         binary = .false.
      elseif (index(ftype,'b')>0) then
         binary = .true.
      else
         binary = .false.
      endif

      inquire(fid, opened=ok, named=hasname, name=fname, unformatted=isbinary)

      ! check that the file is unformatted if binary, formatted if ascii
      if ((index(isbinary,'yes')>0 .or. index(isbinary,'YES')>0) .and. .not. binary) then
         write(0,*) 'Warning: Unable to read rscalar from unformatted file ',trim(fname)
      elseif ((index(isbinary,'no')>0 .or. index(isbinary,'NO')>0) .and. binary) then
         write(0,*) 'Warning: Unable to read rscalar from formatted file ',trim(fname)
      endif

      if (binary) then
         ! read binary from unformatted files and exit
         read(fid) E%nx,E%ny,E%nz,E%gridType
         read(fid) E%v
         return
      end if

	  Nx = size(E%v,1)
	  Ny = size(E%v,2)
	  Nz = size(E%v,3)

	  allocate(temp(Ny),STAT=istat)
	  i = 1
	  do
	    read(fid,*,iostat=istat) k1, k2 ! block numbers have to be there to read
	    if (istat /= 0) exit
		if ((k1 < 0) .or. (k2 > Nz)) then
	    	write(0, *) 'Error reading the ',i,'th block in read_rscalar'
	    	stop
	    else if (k1 > k2) then
	    	write(0, *) 'Warning: block ',i,' in read_rscalar will be ignored'
	    end if
	    do j = Nx,1,-1
	    	read(fid,*,iostat=istat) temp
	    	if (istat /= 0) then
	    	  	write(0, *) 'Error reading the ',j,'th row in ',i,'th block in read_rscalar'
	    	  	stop
	    	end if
	    	do k = k1,k2
	    	  	E%v(j,:,k) = temp
	    	end do
	    end do
	    if (k == Nz) exit
	    i = i+1
      end do
	  deallocate(temp,STAT=istat)

  end subroutine read_rscalar

  !****************************************************************************
  ! read_cscalar reads in an cscalar in a simple format; cscalar has
  ! to exist and be allocated before calling this routine, and the file unit
  ! must already be available for reading; also need to be careful that
  ! dimensions in the input file match those of the cscalar. Each complex
  ! value must be of the form (v1, v2);
  ! assuming the "intuitive geographic" file convention by reading
  ! x from the bottom up (x points North), y left to right (y points East)
  ! optional file type may be 'ascii' or 'binary', shortcuts allowed;
  ! defaults to 'ascii'. use binary for better accuracy.
  subroutine read_cscalar(fid, E, ftype)

      integer,        intent(in)		:: fid
      type (cscalar), intent(inout)		:: E
      character(*), optional, intent(in):: ftype

      !  local variables
      integer 		                    :: Nx, Ny, Nz
      integer                           :: i, j, k, k1, k2, istat
      complex (kind(E%v)), allocatable  :: temp(:)
      logical                           :: ok, hasname, binary
      character(80)                     :: fname, isbinary

      if(.not. E%allocated) then
         write(0, *) 'cscalar must be allocated before call to read_cscalar'
         stop
      endif

      if (.not. present(ftype)) then
         binary = .false.
      elseif (index(ftype,'b')>0) then
         binary = .true.
      else
         binary = .false.
      endif

      inquire(fid, opened=ok, named=hasname, name=fname, unformatted=isbinary)

      ! check that the file is unformatted if binary, formatted if ascii
      if ((index(isbinary,'yes')>0 .or. index(isbinary,'YES')>0) .and. .not. binary) then
         write(0,*) 'Warning: Unable to read cscalar from unformatted file ',trim(fname)
      elseif ((index(isbinary,'no')>0 .or. index(isbinary,'NO')>0) .and. binary) then
         write(0,*) 'Warning: Unable to read cscalar from formatted file ',trim(fname)
      endif

      if (binary) then
         ! read binary from unformatted files and exit
         read(fid) E%nx,E%ny,E%nz,E%gridType
         read(fid) E%v
         return
      end if

	  Nx = size(E%v,1)
	  Ny = size(E%v,2)
	  Nz = size(E%v,3)

	  allocate(temp(Ny),STAT=istat)
	  i = 1
	  do
	    read(fid,*,iostat=istat) k1, k2 ! block numbers have to be there to read
	    if (istat /= 0) exit
		if ((k1 < 0) .or. (k2 > Nz)) then
	    	write(0, *) 'Error reading the ',i,'th block in read_cscalar'
	    	stop
	    else if (k1 > k2) then
	    	write(0, *) 'Warning: block ',i,' in read_cscalar will be ignored'
	    end if
	    do j = Nx,1,-1
	    	read(fid,*,iostat=istat) temp
	    	if (istat /= 0) then
	    	  	write(0, *) 'Error reading the ',j,'th row in ',i,'th block in read_cscalar'
	    	  	stop
	    	end if
	    	do k = k1,k2
	    	  	E%v(j,:,k) = temp
	    	end do
	    end do
	    if (k == Nz) exit
	    i = i+1
      end do
	  deallocate(temp,STAT=istat)

  end subroutine read_cscalar

  !****************************************************************************
  ! read_iscalar reads in an iscalar in a simple ASCII format; iscalar has
  ! to exist and be allocated before calling this routine, and the file unit
  ! must already be available for reading; also need to be careful that
  ! dimensions in the input file match those of the iscalar;
  ! assuming the "intuitive geographic" file convention by reading
  ! x from the bottom up (x points North), y left to right (y points East)
  ! ... conforms to the model format of Weerachai Siripunvaraporn
  ! optional file type may be 'ascii' or 'binary', shortcuts allowed;
  ! defaults to 'ascii'.
  subroutine read_iscalar(fid, E, ftype)

      integer,        intent(in)		:: fid
      type (iscalar), intent(inout)		:: E
      character(*), optional, intent(in):: ftype

      !  local variables
      integer 		                    :: Nx, Ny, Nz
      integer                           :: i, j, k, k1, k2, istat
      integer, allocatable              :: temp(:)
      logical                           :: ok, hasname, binary
      character(80)                     :: fname, isbinary

      if(.not. E%allocated) then
         write(0, *) 'iscalar must be allocated before call to read_iscalar'
         stop
      endif

      if (.not. present(ftype)) then
         binary = .false.
      elseif (index(ftype,'b')>0) then
         binary = .true.
      else
         binary = .false.
      endif

      inquire(fid, opened=ok, named=hasname, name=fname, unformatted=isbinary)

      ! check that the file is unformatted if binary, formatted if ascii
      if ((index(isbinary,'yes')>0 .or. index(isbinary,'YES')>0) .and. .not. binary) then
         write(0,*) 'Warning: Unable to read iscalar from unformatted file ',trim(fname)
      elseif ((index(isbinary,'no')>0 .or. index(isbinary,'NO')>0) .and. binary) then
         write(0,*) 'Warning: Unable to read iscalar from formatted file ',trim(fname)
      endif

      if (binary) then
         ! read binary from unformatted files and exit
         read(fid) E%nx,E%ny,E%nz,E%gridType
         read(fid) E%v
         return
      end if

	  Nx = size(E%v,1)
	  Ny = size(E%v,2)
	  Nz = size(E%v,3)

	  allocate(temp(Ny),STAT=istat)
	  i = 1
	  do
	    read(fid,*,iostat=istat) k1, k2 ! block numbers have to be there to read
	    if (istat /= 0) exit
	    if ((k1 < 0) .or. (k2 > Nz)) then
	    	write(0, *) 'Error reading the ',i,'th block in read_iscalar'
	    	stop
	    else if (k1 > k2) then
	    	write(0, *) 'Warning: block ',i,' in read_iscalar will be ignored'
	    end if
	    do j = Nx,1,-1
	    	read(fid,*,iostat=istat) temp
	    	if (istat /= 0) then
	    	  	write(0, *) 'Error reading the ',j,'th row in ',i,'th block in read_iscalar'
	    	  	stop
	    	end if
	    	do k = k1,k2
	    	  	E%v(j,:,k) = temp
	    	end do
	    end do
	    if (k == Nz) exit
	    i = i+1
      end do
	  deallocate(temp,STAT=istat)

  end subroutine read_iscalar

  !****************************************************************************
  ! write_rscalar writes an iscalar in a simple ASCII format; rscalar has
  ! to exist and be allocated before calling this routine, and the file unit
  ! must already be available for writing;
  ! assuming the "intuitive geographic" file convention by writing
  ! x from the bottom up (x points North), y left to right (y points East)
  ! optional file type may be 'ascii' or 'binary', shortcuts allowed;
  ! defaults to 'ascii'. use binary for better accuracy.
  subroutine write_rscalar(fid, E, ftype)

      integer,        intent(in)		:: fid
      type (rscalar), intent(in)		:: E
      character(*), optional, intent(in):: ftype

      !  local variables
      integer 		                    :: Nx, Ny, Nz
      integer                           :: i, j, k, k1, k2, istat
      real (kind(E%v)), allocatable     :: temp(:,:)
      logical                           :: ok, hasname, binary
      character(80)                     :: fname, isbinary

      if(.not. E%allocated) then
         write(0, *) 'rscalar must be allocated before call to write_rscalar'
         stop
      endif

      if (.not. present(ftype)) then
         binary = .false.
      elseif (index(ftype,'b')>0) then
         binary = .true.
      else
         binary = .false.
      endif

      inquire(fid, opened=ok, named=hasname, name=fname, unformatted=isbinary)

      ! check that the file is unformatted if binary, formatted if ascii
      if ((index(isbinary,'yes')>0 .or. index(isbinary,'YES')>0) .and. .not. binary) then
         write(0,*) 'Warning: Unable to write rscalar to unformatted file ',trim(fname)
      elseif ((index(isbinary,'no')>0 .or. index(isbinary,'NO')>0) .and. binary) then
         write(0,*) 'Warning: Unable to write rscalar to formatted file ',trim(fname)
      endif

      ! write binary to unformatted files
      if (binary) then
         write(fid) E%nx,E%ny,E%nz,E%gridType
         write(fid) E%v
         return
      end if

      Nx = size(E%v,1)
	  Ny = size(E%v,2)
	  Nz = size(E%v,3)

	  allocate(temp(Nx,Ny),STAT=istat)
      k1 = 1
	  do
	  	k2 = Nz
	    do k = k1,Nz-1
	    	temp = abs(E%v(:,:,k+1) - E%v(:,:,k))
	    	if (maxval(temp) > TOL6) then
	    		k2 = k
	    		exit
	    	end if
	    end do
	    write(fid,'(2i5)',iostat=istat) k1, k2
	    if (istat /= 0) then
	    	write(0, *) 'Error writing to file in write_rscalar'
	    	stop
	    end if
	    temp = E%v(:,:,k1)
	    do i = Nx,1,-1
	    	do j = 1,Ny
	    		write(fid,'(es13.5)',iostat=istat,advance='no') E%v(i,j,k1)
	    	end do
	    	write(fid,*)
	    end do
	    k1 = k2+1
	    if (k1 > Nz) exit
      end do
	  deallocate(temp,STAT=istat)

  end subroutine write_rscalar

  !****************************************************************************
  ! write_cscalar writes an iscalar in a simple format; cscalar has
  ! to exist and be allocated before calling this routine, and the file unit
  ! must already be available for writing;
  ! assuming the "intuitive geographic" file convention by writing
  ! x from the bottom up (x points North), y left to right (y points East)
  ! optional file type may be 'ascii' or 'binary', shortcuts allowed;
  ! defaults to 'ascii'. use binary for better accuracy.
  subroutine write_cscalar(fid, E, ftype)

      integer,        intent(in)		:: fid
      type (cscalar), intent(in)		:: E
      character(*), optional, intent(in):: ftype

      !  local variables
      integer 		                    :: Nx, Ny, Nz
      integer                           :: i, j, k, k1, k2, istat
      complex (kind(E%v)), allocatable  :: temp(:,:)
      logical                           :: ok, hasname, binary
      character(80)                     :: fname, isbinary

      if(.not. E%allocated) then
         write(0, *) 'cscalar must be allocated before call to write_cscalar'
         stop
      endif

      if (.not. present(ftype)) then
         binary = .false.
      elseif (index(ftype,'b')>0) then
         binary = .true.
      else
         binary = .false.
      endif

      inquire(fid, opened=ok, named=hasname, name=fname, unformatted=isbinary)

      ! check that the file is unformatted if binary, formatted if ascii
      if ((index(isbinary,'yes')>0 .or. index(isbinary,'YES')>0) .and. .not. binary) then
         write(0,*) 'Warning: Unable to write cscalar to unformatted file ',trim(fname)
      elseif ((index(isbinary,'no')>0 .or. index(isbinary,'NO')>0) .and. binary) then
         write(0,*) 'Warning: Unable to write cscalar to formatted file ',trim(fname)
      endif

      ! write binary to unformatted files
      if (binary) then
         write(fid) E%nx,E%ny,E%nz,E%gridType
         write(fid) E%v
         return
      end if

      Nx = size(E%v,1)
	  Ny = size(E%v,2)
	  Nz = size(E%v,3)

	  allocate(temp(Nx,Ny),STAT=istat)
      k1 = 1
	  do
	  	k2 = Nz
	    do k = k1,Nz-1
	    	temp = abs(E%v(:,:,k+1) - E%v(:,:,k))
	    	if ((maxval(real(temp)) > TOL6) .or. (maxval(imag(temp)) > TOL6)) then
	    		k2 = k
	    		exit
	    	end if
	    end do
	    write(fid,'(2i5)',iostat=istat) k1, k2
	    if (istat /= 0) then
	    	write(0, *) 'Error writing to file in write_cscalar'
	    	stop
	    end if
	    temp = E%v(:,:,k1)
	    do i = Nx,1,-1
	    	do j = 1,Ny
	    		write(fid,'(es13.5)',iostat=istat,advance='no') E%v(i,j,k1)
	    	end do
	    	write(fid,*)
	    end do
	    k1 = k2+1
	    if (k1 > Nz) exit
      end do
	  deallocate(temp,STAT=istat)

  end subroutine write_cscalar

  !****************************************************************************
  ! write_iscalar writes an iscalar in a simple format; iscalar has
  ! to exist and be allocated before calling this routine, and the file unit
  ! must already be available for writing;
  ! assuming the "intuitive geographic" file convention by writing
  ! x from the bottom up (x points North), y left to right (y points East)
  ! ... ascii conforms to the model format of Weerachai Siripunvaraporn
  ! optional file type may be 'ascii' or 'binary', shortcuts allowed;
  ! defaults to 'ascii'.
  subroutine write_iscalar(fid, E, ftype)

      integer,        intent(in)		:: fid
      type (iscalar), intent(in)		:: E
      character(*), optional, intent(in):: ftype

      !  local variables
      integer 		                    :: Nx, Ny, Nz
      integer                           :: i, j, k, k1, k2, istat
      integer, allocatable              :: temp(:,:)
      character(10)                     :: fmt
      logical                           :: ok, hasname, binary
      character(80)                     :: fname, isbinary

      if(.not. E%allocated) then
         write(0, *) 'iscalar must be allocated before call to write_iscalar'
         stop
      endif

      if (.not. present(ftype)) then
         binary = .false.
      elseif (index(ftype,'b')>0) then
         binary = .true.
      else
         binary = .false.
      endif

      inquire(fid, opened=ok, named=hasname, name=fname, unformatted=isbinary)

      ! check that the file is unformatted if binary, formatted if ascii
      if ((index(isbinary,'yes')>0 .or. index(isbinary,'YES')>0) .and. .not. binary) then
         write(0,*) 'Warning: Unable to write iscalar to unformatted file ',trim(fname)
      elseif ((index(isbinary,'no')>0 .or. index(isbinary,'NO')>0) .and. binary) then
         write(0,*) 'Warning: Unable to write iscalar to formatted file ',trim(fname)
      endif

      ! write binary to unformatted files
      if (binary) then
         ! (AK) CORRECT:
         !write(fid) E%nx,E%ny,E%nz,trim(E%gridType)
         ! (AK) DEBUG & need version compatibility; output gibberish header:
         write(fid) R_ZERO,E%nx,E%ny,'*******testing******'
         ! (AK) END DEBUG
         write(fid) E%v
         return
      end if

	  Nx = size(E%v,1)
	  Ny = size(E%v,2)
	  Nz = size(E%v,3)

	  allocate(temp(Nx,Ny),STAT=istat)
      k1 = 1
	  do
	  	k2 = Nz
	    do k = k1,Nz-1
	    	temp = abs(E%v(:,:,k+1) - E%v(:,:,k))
	    	if (maxval(temp) > 0) then
	    		k2 = k
	    		exit
	    	end if
	    end do
	    if (Nz < 1000) then
	    	fmt = '(2i4)'
	    else
	    	fmt = '(2i6)'
	    end if
	    write(fid,fmt,iostat=istat) k1, k2
	    if (istat /= 0) then
	    	write(0, *) 'Error writing to file in write_iscalar'
	    	stop
	    end if
	    temp = E%v(:,:,k1)
   		if (maxval(temp) <= 9) then
   			fmt = '(i2)'
   		else
   			fmt = '(i6)'
   		end if
	    do i = Nx,1,-1
	    	do j = 1,Ny
	    		write(fid,fmt,iostat=istat,advance='no') E%v(i,j,k1)
	    	end do
	    	write(fid,*)
	    end do
	    k1 = k2+1
	    if (k1 > Nz) exit
      end do
	  deallocate(temp,STAT=istat)

  end subroutine write_iscalar

  !****************************************************************************
  ! copy_rscalar makes an exact copy of derived data type
  ! rscalar;   NOTE: first argument is output
  subroutine copy_rscalar(E2, E1)

    implicit none
    type (rscalar), intent(in)       :: E1
    type (rscalar), intent(inout)    :: E2
    integer	                     :: status

    ! check to see if RHS (E1) is active (allocated)
    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for copy_rscalar'
    else

       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! just copy components
             E2%v = E1%v
             E2%gridType = E1%gridType

          else
             write (0, *) 'not compatible usage in copy_rscalar'
          end if

       else

          if(E2%allocated) then
             ! first deallocate memory for v
             deallocate(E2%v,STAT=status)
          end if

          !  then allocate E2 as correct size ...
          Call create_rscalar(E1%grid, E2, E1%gridType)
          !   .... and copy E1
          ! for the large grids we need to code this explicitly(compiler problems)

          E2%v = E1%v

          E2%gridType = E1%gridType

       end if

    end if

    if(E1%temporary) then
    	call deall_rscalar(E1)
    end if

  end subroutine copy_rscalar  ! copy_rscalar


  !****************************************************************************
  ! copy_cscalar makes an exact copy of derived data type
  ! cscalar; and NOTE: E2 is the output
  subroutine copy_cscalar(E2, E1)

    implicit none
    type (cscalar), intent(in)            :: E1
    type (cscalar), intent(inout)         :: E2
    integer                               :: status

    ! check to see if RHS (E1) is active (allocated)
    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for copy_cscalar'
    else

       if((E2%nx == E1%nx).and.(E2%ny == E1%ny).and.(E2%nz == E1%nz)) then

          if (E1%gridType == E2%gridType) then

             ! just copy components
             E2%v = E1%v
             E2%gridType = E1%gridType

          else
             write (0, *) 'not compatible usage in copy_cscalar'
          end if

       else

          if(E2%allocated) then
             ! first deallocate memory for v
             deallocate(E2%v, STAT=status)
          end if

          !  then allocate E2 as correct size ...
          Call create(E1%grid, E2, E1%gridType) ! create_cscalar
          !   .... and copy E1
          E2%v = E1%v
          E2%gridType = E1%gridType

       end if

    end if

    if(E1%temporary) then
    	call deall(E1) ! deall_cscalar
    end if

  end subroutine copy_cscalar


  !****************************************************************************
  ! copy_iscalar makes an exact copy of derived data type
  ! rscalar;   NOTE: first argument is output
  subroutine copy_iscalar(E2, E1)

    implicit none
    type (iscalar), intent(in)       :: E1
    type (iscalar), intent(inout)    :: E2
    integer	                     :: status

    ! check to see if RHS (E1) is active (allocated)
    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for copy_iscalar'
    else

       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! just copy components
             E2%v = E1%v
             E2%gridType = E1%gridType

          else
             write (0, *) 'not compatible usage in copy_iscalar'
          end if

       else

          if(E2%allocated) then
             ! first deallocate memory for v
             deallocate(E2%v,STAT=status)
          end if

          !  then allocate E2 as correct size ...
          Call create(E1%grid, E2, E1%gridType) ! create_iscalar
          !   .... and copy E1
          E2%v = E1%v
          E2%gridType = E1%gridType

       end if

    end if

    if(E1%temporary) then
    	call deall(E1) !deall_iscalar
    end if

  end subroutine copy_iscalar

  !****************************************************************************
  ! zero_rscalar zeros variable of derived data type
  ! rscalar;
  subroutine zero_rscalar(E)

    implicit none
    type (rscalar), intent(inout)   :: E

    ! check to see if E is active (allocated)
    if(.not.E%allocated) then
       write(0,*) 'Error in zero_rscalar: E not allocated'
    else

       E%v = R_ZERO

    end if

  end subroutine zero_rscalar


  !****************************************************************************
  ! zero_cscalar zeros variable of derived data type
  ! cscalar
  subroutine zero_cscalar(E)

    implicit none
    type (cscalar), intent(inout) :: E

    ! check to see if E is active (allocated)
    if(.not.E%allocated) then
       write(0,*) 'Error in zero_cscalar: E not allocated'
    else

       E%v = C_ZERO

    end if

  end subroutine zero_cscalar

  !****************************************************************************
  ! scMult_cscalar multiplies scalar stored as devired data type
  ! cscalar with a complex scalar; subroutine version
  ! E2 can overwrite E1
  subroutine scMult_cscalar(c, E1, E2)

    implicit none
    complex(kind=prec), intent(in)                      :: c
    ! a complex scalar to be multiplied with
    type (cscalar), intent(in)                       :: E1
    type (cscalar), intent(inout)                    :: E2

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMult_cscalar'
       stop
    endif

    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated yet for scMult_cscalar'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! complex scalar multiplication for v-component

             E2%v = E1%v * c

          else
             write (0, *) 'not compatible usage in scMult_cscalar'
          end if

       else
          write(0, *) 'Error:scMult_cscalar: scalars not same size'

       end if
    end if

  end subroutine scMult_cscalar

  !****************************************************************************
  ! scMult_cscalar_f multiplies scalar stored as devired data type
  ! cscalar with a complex scalar; function version
  function scMult_cscalar_f(c, E1) result(E2)

    implicit none
    complex(kind=prec), intent(in)                      :: c
    ! a complex scalar to be multiplied with
    type (cscalar), intent(in)                       :: E1
    type (cscalar)                                   :: E2

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMult_cscalar_f'
       stop
    endif

    ! In function version, appropriate data types need to be created
    Call create_cscalar(E1%grid, E2, E1%gridType)
    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated for scMult_cscalar_f'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! complex scalar multiplication for v-component

             E2%v = E1%v * c

          else
             write (0, *) 'not compatible usage in scMult_cscalar_f'
          end if

       else

          write(0, *) 'Error:scMult_cscalar_f: scalars not same size'

       end if
    end if

    E2%temporary = .true.

  end function scMult_cscalar_f

  ! ***************************************************************************
  ! scMult_rscalar multiplies scalar stored as devired data type
  ! rscalar with a real scalar; subroutine version
  ! E2 can overwrite E1
  subroutine scMult_rscalar(c, E1, E2)

    implicit none
    real (kind=prec), intent(in)                         :: c
    ! a real scalar to be multiplied with
    type (rscalar), intent(in)                       :: E1
    type (rscalar), intent(inout)                    :: E2

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMult_rscalar'
       stop
    endif

    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated for scMult_rscalar'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! real scalar multiplication for v-component
             E2%v = E1%v * c

          else
             write (0, *) 'not compatible usage for scMult_rscalar'
          end if

       else

          write(0, *) 'Error:scMult_rscalar: scalars not same size'

       end if
    end if

  end subroutine scMult_rscalar

  !****************************************************************************
  ! scMult_rscalar_f multiplies scalar stored as devired data type
  ! rscalar with a real scalar; function version
  function scMult_rscalar_f(c, E1) result(E2)

    implicit none
    real (kind=prec), intent(in)                         :: c
    ! a complex scalar to be multiplied with
    type (rscalar), intent(in)                       :: E1
    type (rscalar)                                   :: E2

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMult_rscalar_f'
       stop
    endif

    ! In function version, appropriate data types need to be created
    Call create_rscalar(E1%grid, E2, E1%gridType)
    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated for scMult_rscalar_f'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! real scalar multiplication for v-component
             E2%v = E1%v * c

          else
             write (0, *) 'not compatible usage for scMult_rscalar_f'
          end if

       else

          write(0, *) 'Error:scMult_rscalar_f: scalars not same size'

       end if
    end if

    E2%temporary = .true.

  end function scMult_rscalar_f

  !****************************************************************************
  ! scMultadd_rscalar multiplies scalar E1 stored as derived data type
  ! rscalar with a real scalar r, adding result to output scalar E2
  subroutine scMultAdd_rscalar(r, E1, E2)

    implicit none
    real(kind=prec), intent(in)                      :: r
    ! a real scalar to be multiplied with
    type (rscalar), intent(in)                       :: E1
    type (rscalar), intent(inout)                    :: E2

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMultAdd_rscalar'
       stop
    endif

    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated for scMultAdd_rscalar'
    else

       ! Check whether both scalars are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! complex scalar multiplication for v-component
             E2%v = E2%v + E1%v * r

          else
             write (0, *) 'not compatible usage for scMultAdd_rscalar'
          end if

       else

          write(0, *) 'Error:scMultAdd_rscalar: scalars not same size'

       end if
    end if

  end subroutine scMultAdd_rscalar

  !****************************************************************************
  ! add_rscalar adds scalars stored as devired data type
  ! rscalar with ; subroutine version
  ! E3 can overwrite E1 and E2
  subroutine add_rscalar(E1, E2, E3)

    implicit none
    type (rscalar), intent(in)               :: E1, E2
    type (rscalar), intent(inout)            :: E3

   if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for add_rscalar'
       stop
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for add_rscalar'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! add v-component
             E3%v = E1%v + E2%v

          else
             write (0, *) 'not compatible usage for add_rscalar'
          end if

       else

          write(0, *) 'Error:add_rscalar: scalars not same size'

       end if
    end if

  end subroutine add_rscalar

  !****************************************************************************
  ! add_rscalar_f adds scalars stored as devired data type
  ! rscalar with ; function version
  function add_rscalar_f(E1, E2) result(E3)

    implicit none
    type (rscalar), intent(in)               :: E1, E2
    type (rscalar)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for add_rscalar_f'
       stop
    endif

    ! In function version, appropriate data types need to be created
    Call create_rscalar(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for add_rscalar_f'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! add v-component
             E3%v = E1%v + E2%v

          else
             write (0, *) 'not compatible usage for add_rscalar_f'
          end if

       else

          write(0, *) 'Error:add_rscalar_f: scalars not same size'

       end if
    end if

    E3%temporary = .true.

  end function add_rscalar_f


  !****************************************************************************
  ! add_cscalar adds scalars stored as devired data type
  ! cscalar with ; subroutine version
  ! E3 can overwrite E1 and E2
  subroutine add_cscalar(E1, E2, E3)

    implicit none
    type (cscalar), intent(in)               :: E1, E2
    type (cscalar), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for add_cscalar'
       stop
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS not allocated for add_cscalar'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! add v-component
             E3%v = E1%v + E2%v

          else
             write (0, *) 'not compatible usage for add_cscalar'
          end if

       else

          write(0, *) 'Error:add_cscalar: scalars not same size'

       end if
    end if

  end subroutine add_cscalar

  !****************************************************************************
  ! add_cscalar_f adds scalars stored as devired data type
  ! cscalar with ; function version
  function add_cscalar_f(E1, E2) result(E3)

    implicit none
    type (cscalar), intent(in)               :: E1, E2
    type (cscalar)                           :: E3

     if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for add_cscalar_f'
       stop
    endif

    ! In function version, appropriate data types need to be created
    Call create_cscalar(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for add_cscalar_f'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! add v-component
             E3%v = E1%v + E2%v

          else
             write (0, *) 'not compatible usage for add_cscalar_f'
          end if

       else

          write(0, *) 'Error:add_cscalar_f: scalars not same size'

       end if
    end if

    E3%temporary = .true.

  end function add_cscalar_f

  !****************************************************************************
  ! subtract_rscalar subtracts scalars stored as devired data type
  ! rscalar with ; subroutine version
  ! E3 can overwrite E1 and E2
  subroutine subtract_rscalar(E1, E2, E3)

    implicit none
    type (rscalar), intent(in)               :: E1, E2
    type (rscalar), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for subtract_rscalar'
       stop
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for subtract_rscalar'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! subtract v-component
             E3%v = E1%v - E2%v

          else
             write (0, *) 'not compatible usage for subtract_rscalar'
          end if

       else

          write(0, *) 'Error:add_rscalar: scalars not same size'

       end if
    end if

  end subroutine subtract_rscalar

  !****************************************************************************
  ! subtract_rscalar_f subtracts scalars stored as devired data type
  ! rscalar with ; function version
  function subtract_rscalar_f(E1, E2) result(E3)

    implicit none
    type (rscalar), intent(in)               :: E1, E2
    type (rscalar)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for add_rscalar_f'
       stop
    endif

    ! In function version, appropriate data types need to be created
    Call create_rscalar(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for add_rscalar_f'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! subtract v-component
             E3%v = E1%v - E2%v

          else
             write (0, *) 'not compatible usage for subtract_rscalar_f'
          end if

       else

          write(0, *) 'Error:subtract_rscalar_f: scalars not same size'

       end if
    end if

    E3%temporary = .true.

  end function subtract_rscalar_f

  !****************************************************************************
  ! subtract_cscalar subtracts scalars stored as devired data type
  ! cscalar with ; subroutine version
  ! E3 can overwrite E1 and E2
  subroutine subtract_cscalar(E1, E2, E3)

    implicit none
    type (cscalar), intent(in)               :: E1, E2
    type (cscalar), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for subtract_cscalar'
       stop
    endif

    ! check to see if LHS (E2) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS not allocated for subtract_cscalar'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! subtract v-component
             E3%v = E1%v - E2%v

          else
             write (0, *) 'not compatible usage for subtract_cscalar'
          end if

       else

          write(0, *) 'Error:subtract_cscalar: scalars not same size'

       end if
    end if

  end subroutine subtract_cscalar

  !****************************************************************************
  ! subtract_cscalar_f subtracts scalars stored as devired data type
  ! cscalar with ; function version
  function subtract_cscalar_f(E1, E2) result(E3)

    implicit none
    type (cscalar), intent(in)               :: E1, E2
    type (cscalar)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for subtract_cscalar_f'
       stop
    endif

    ! In function version, appropriate data types need to be created
    Call create_cscalar(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for subtract_cscalar_f'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! subtract v-component
             E3%v = E1%v - E2%v

          else
             write (0, *) 'not compatible usage for subtract_cscalar_f'
          end if

       else

          write(0, *) 'Error:subtract_cscalar_f: scalars not same size'

       end if
    end if

    E3%temporary = .true.

  end function subtract_cscalar_f

  !****************************************************************************
  ! diagMult_rscalar multiplies two scalars E1, E2 stored as devired data
  ! type rscalar pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagMult_rscalar(E1, E2, E3)

    implicit none
    type (rscalar), intent(in)               :: E1, E2
    type (rscalar), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_rscalar'
       stop
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_rscalar'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for v-component

             E3%v = E1%v * E2%v

          else
             write (0, *) 'not compatible usage for diagMult_rscalar'
          end if

       else

          write(0, *) 'Error:diagMult_rscalar: scalars not same size'

       end if
    end if

  end subroutine diagMult_rscalar

  !****************************************************************************
  ! diagMult_rscalar_f multiplies two scalars E1, E2 stored as devired
  ! data type rscalar pointwise; function version
  function diagMult_rscalar_f(E1, E2) result(E3)

    implicit none
    type (rscalar), intent(in)               :: E1, E2
    type (rscalar)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_rscalar_f'
       stop
    endif

    ! In function version, appropriate data types need to be created
    Call create_rscalar(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_rscalar_f'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for v-component
             E3%v = E1%v * E2%v

          else
             write (0, *) 'not compatible usage for diagMult_rscalar_f'
          end if

       else

          write(0, *) 'Error:diagMult_rscalar_f: scalars not same size'

       end if
    end if

    E3%temporary = .true.

  end function diagMult_rscalar_f

  !****************************************************************************
  ! diagDiv_rscalar divides scalar E1 by scalar E2 to get scalar E3
  ! type rscalar pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagDiv_rscalar(E1, E2, E3)

    implicit none
    type (rscalar), intent(in)               :: E1, E2
    type (rscalar), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagDiv_rscalar'
       stop
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagDiv_rscalar'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise division for v-component; NO ZERO TESTING FOR EFFICIENCY
             E3%v = E1%v / E2%v

          else
             write (0, *) 'not compatible usage for diagDiv_rscalar'
          end if

       else

          write(0, *) 'Error:diagDiv_rscalar: scalars not same size'

       end if
    end if

  end subroutine diagDiv_rscalar ! diagDiv_rscalar
  !****************************************************************************
  ! diagDiv_crscalar divides scalar E1 by scalar E2 to get scalar E3
  ! type cscalar pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagDiv_crscalar(E1, E2, E3)

    implicit none
    type (cscalar), intent(in)               :: E1
    type (rscalar), intent(in)               :: E2
    type (cscalar), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagDiv_rscalar'
       stop
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagDiv_rscalar'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise division for v-component; NO ZERO TESTING FOR EFFICIENCY
             E3%v = E1%v / E2%v

          else
             write (0, *) 'not compatible usage for diagDiv_rscalar'
          end if

       else

          write(0, *) 'Error:diagDiv_rscalar: scalars not same size'

       end if
    end if

  end subroutine diagDiv_crscalar ! diagDiv_cscalar


  !****************************************************************************
  ! diagMult_cscalar multiplies two scalars E1, E2 stored as devired data
  ! type cscalar pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagMult_cscalar(E1, E2, E3)

    implicit none
    type (cscalar), intent(in)               :: E1, E2
    type (cscalar), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_cscalar'
       stop
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diaMult_cscalar'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for v-component
             E3%v = E1%v * E2%v

          else
             write (0, *) 'not compatible usage for diagMult_cscalar'
          end if

       else

          write(0, *) 'Error:diagMult_cscalar: scalars not same size'

       end if
    end if

  end subroutine diagMult_cscalar

  !****************************************************************************
  ! diagMult_cscalar_f multiplies two scalars E1, E2 stored as devired
  ! data  type cscalar pointwise; function version
  function diagMult_cscalar_f(E1, E2) result(E3)

    implicit none
    type (cscalar), intent(in)               :: E1, E2
    type (cscalar)                           :: E3

   if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_cscalar_f'
       stop
    endif

    ! In function version, appropriate data types need to be created
    Call create_cscalar(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_cscalar_f'
    else

       ! Check whether both scalars are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for v-component
             E3%v = E1%v * E2%v

          else
             write (0, *) 'not compatible usage for diagMult_cscalar_f'
          end if

       else

          write(0, *) 'Error:diagMult_cscalar_f: scalars not same size'

       end if
    end if

    E3%temporary = .true.

  end function diagMult_cscalar_f

  !****************************************************************************
  ! diagMult_crscalar multiplies scalar E1 with scalar E2 stored as
  ! devired type cscalar pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagMult_crscalar(E1, E2, E3)

    implicit none
    type (cscalar), intent(in)               :: E1
    type (rscalar), intent(in)               :: E2
    type (cscalar), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_crscalar'
       stop
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_crscalar'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for v-component
             E3%v = E1%v * E2%v

          else
             write (0, *) 'not compatible usage for diagMult_crscalar'
          end if

       else

          write(0, *) 'Error:diagMult_crscalar: scalars not same size'

       end if
    end if

  end subroutine diagMult_crscalar

  !****************************************************************************
  ! diagMult_crscalar_f multiplies scalar E1 with scalar E2 stored as
  ! derived data type cscalar pointwise; function version
  function diagMult_crscalar_f(E1, E2) result(E3)

    implicit none
    type (cscalar), intent(in)               :: E1
    type (rscalar), intent(in)               :: E2
    type (cscalar)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_crscalar_f'
       stop
    endif

    ! In function version, appropriate data types need to be created
    Call create_cscalar(E1%grid, E3, E1%gridType)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_crscalar_f'
    else

       ! Check whether both scalars are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for v-component
             E3%v = E1%v * E2%v

          else
             write (0, *) 'not compatible usage for diagMult_crscalar_f'
          end if

       else

          write(0, *) 'Error:diagMult_crscalar_f: scalars not same size'

       end if
    end if

    E3%temporary = .true.

  end function diagMult_crscalar_f

  !****************************************************************************
  ! diagMult_rcscalar multiplies scalar E1 with scalar E2 stored as
  ! derived type cscalar pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagMult_rcscalar(E1, E2, E3)

    implicit none
    type (rscalar), intent(in)               :: E1
    type (cscalar), intent(in)               :: E2
    type (cscalar), intent(inout)            :: E3


   if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_rcscalar'
       stop
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_rcscalar'
    else

       ! Check whether all the scalar nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for v-component

                   E3%v = E1%v * E2%v

          else
             write (0, *) 'not compatible usage for diagMult_rcscalar'
          end if

       else

          write(0, *) 'Error:diagMult_rcscalar: scalars not same size'

       end if
    end if


  end subroutine diagMult_rcscalar

  !****************************************************************************
  ! diagMult_rcscalar_f multiplies scalar E1 with scalar E2 stored as
  ! derived data type cscalar pointwise; function version
  function diagMult_rcscalar_f(E1, E2) result(E3)

    implicit none
    type (rscalar), intent(in)               :: E1
    type (cscalar), intent(in)               :: E2
    type (cscalar)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_rcscalar_f'
       stop
    endif

    ! In function version, appropriate data types need to be created
    Call create_cscalar(E1%grid, E3, E1%gridType)
    ! check to see if RHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_rcscalar_f'
    else

       ! Check whether both scalars are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and. &
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! pointwise multiplication for v-component
             E3%v = E1%v * E2%v

          else
             write (0, *) 'not compatible usage for diagMult_rcscalar_f'
          end if

       else

          write(0, *) 'Error:diagMult_rcscalar_f: scalars not same size'

       end if
    end if

    E3%temporary = .true.

  end function diagMult_rcscalar_f

  !****************************************************************************
  ! dotProd_rscalar computes dot product of two vecors stored
  ! as derived data type rscalar, returning a real number
  function dotProd_rscalar_f(E1, E2) result(r)

    implicit none
    type (rscalar), intent(in)   :: E1, E2
    real (kind=prec)		     :: r

    r = R_ZERO

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for dotProd_rscalar'
       stop
    endif

    ! Check whether both input scalars are of the same size
    if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

       if (E1%gridType == E2%gridType) then

          r = r + sum(E1%v * E2%v)

       else
          write (0, *) 'not compatible usage for dotProd_rscalar'
       end if

    else

       write(0, *) 'Error:dotProd_rscalar: scalars not same size'

    end if

  end function dotProd_rscalar_f

  !****************************************************************************
  ! dotProd_cscalar computes dot product of two vecors stored
  ! as derived data type cscalar, returning a complex number
  function dotProd_cscalar_f(E1, E2) result(c)

    implicit none
    type (cscalar), intent(in)       :: E1, E2
    complex(kind=prec)		     :: c

    c = C_ZERO

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for dotProd_cscalar'
       stop
    endif

    ! Check whether both input scalars are of the same size
    if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

       if(E1%gridType == E2%gridType) then

          c = c + sum(conjg(E1%v) * E2%v)

       else
          write (0, *) 'not compatible usage for dotProd_cscalar'
       end if

    else

       write(0, *) 'Error:dotProd_cscalar: scalars not same size'

    end if

  end function dotProd_cscalar_f

  !****************************************************************************
  ! linComb_rscalar computes linear combination of two scalars
  ! stored as derived data type rscalar; subroutine, not a function
  ! both input scalars must have the same dimension
  subroutine linComb_rscalar(inc1, E1, inc2, E2, E3)

    implicit none
    !   input scalars
    type (rscalar), intent(in)             :: E1, E2
    !  input complex scalars
    real (kind=prec), intent(in)           :: inc1, inc2
    type (rscalar), intent(inout)          :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for linComb_rscalar'
       stop
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for linComb_rscalar'
    else

       ! Check whether all scalars are of the same size
       if ((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! form linear combinatoin
             E3%v = inc1*E1%v + inc2*E2%v

          else
             write (0, *) 'not compatible usage for linComb_rscalar'
          end if

       else

          write(0, *) 'Error:linComb_rscalar:  scalars not same size'

       end if
    end if

  end subroutine linComb_rscalar

  !****************************************************************************
  ! linComb_cscalar computes linear combination of two scalars
  ! stored as derived data type cscalar; subroutine, not a function
  ! both input scalars must have the same dimension
  subroutine linComb_cscalar(inc1, E1, inc2, E2, E3)

    implicit none
    !   input scalars
    type (cscalar), intent(in)             :: E1, E2
    !  input complex scalars
    complex (kind=prec), intent(in)           :: inc1, inc2
    type (cscalar), intent(inout)          :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for linComb_cscalar'
       stop
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for linComb_cscalar'
    else

       ! Check whether all scalars are of the same size
       if ((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then

          if ((E1%gridType == E2%gridType).and.(E1%gridType == E3%gridType)) then

             ! form linear combinatoin
             E3%v = inc1*E1%v + inc2*E2%v

          else
             write (0, *) 'not compatible usage for linComb_cscalar'
          end if

       else

          write(0, *) 'Error:linComb_cscalar:  scalars not same size'

       end if
    end if

  end subroutine linComb_cscalar

  !****************************************************************************
  ! scMultadd_cscalar multiplies scalar E1 stored as derived data type
  ! cscalar with a complex scalar c, adding result to output scalar E2
  subroutine scMultAdd_cscalar(c, E1, E2)

    implicit none
    complex(kind=prec), intent(in)                      :: c
    ! a complex scalar to be multiplied with
    type (cscalar), intent(in)                       :: E1
    type (cscalar), intent(inout)                    :: E2

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMultAdd_cscalar'
       stop
    endif

    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated for scMultAdd_cscalar'
    else

       ! Check whether both scalars are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then

          if (E1%gridType == E2%gridType) then

             ! complex scalar multiplication for v-component
             E2%v = E2%v + E1%v * c

          else
             write (0, *) 'not compatible usage for scMultAdd_cscalar'
          end if

       else

          write(0, *) 'Error:scMultAdd_cscalar: scalars not same size'

       end if
    end if

  end subroutine scMultAdd_cscalar

  ! **********************************************************************
  ! * Creates a random perturbation in rscalar - used for testing
  ! * (developed by Lana Erofeeva & tested by Anna Kelbert)
  subroutine random_rscalar(E,eps)

    implicit none
    type (rscalar), intent(inout)                    :: E
    real(8), intent(in), optional                    :: eps
    ! local
    real (kind(E%v)), allocatable, dimension(:,:,:)  :: zr
    integer              :: Nx,Ny,Nz,istat

    if (.not. E%allocated) then
      call warning('rscalar not allocated in random_rscalar')
      return
    end if

    call zero_rscalar(E)

    ! make some random vectors
    Nx = E%nx
    Ny = E%ny
    Nz = E%nz
    allocate(zr(Nx+1,Ny+1,Nz+1),STAT=istat)
    call random_number(zr)
    if (present(eps)) then
        zr = zr * eps
    else
        zr = zr * 0.05
    end if
    if ((E%gridType == CENTER) .or. (E%gridType == CELL_EARTH)) then
       E%v = zr(1:Nx,1:Ny,1:Nz)
    else if (E%gridType == CORNER) then
       E%v = zr
    else
       write (0, *) 'gridType == ',trim(E%gridType),' undefined in random_rscalar'
    end if
    deallocate(zr,STAT=istat)

  end subroutine random_rscalar

  ! **********************************************************************
  ! * Creates a random perturbation in cscalar - used for testing
  ! * (developed by Lana Erofeeva & tested by Anna Kelbert)
  subroutine random_cscalar(E,eps)

    implicit none
    type (cscalar), intent(inout)                    :: E
    real(8), intent(in), optional                    :: eps
    ! local
    real (kind(E%v)), allocatable, dimension(:,:,:)  :: zr,zi
    integer              :: Nx,Ny,Nz,istat

    if (.not. E%allocated) then
      call warning('cscalar not allocated in random_cscalar')
      return
    end if

    call zero_cscalar(E)

    ! make some random vectors
    Nx = E%nx
    Ny = E%ny
    Nz = E%nz
    allocate(zr(Nx+1,Ny+1,Nz+1),zi(Nx+1,Ny+1,Nz+1),STAT=istat)
    call random_number(zr)
    call random_number(zi)
    if (present(eps)) then
        zr = zr * eps
        zi = zi * eps
    else
        zr = zr * 0.05
        zi = zi * 0.05
    end if
    if (E%gridType == CENTER) then
	   E%v = cmplx(zr(1:Nx,1:Ny,1:Nz),zi(1:Nx,1:Ny,1:Nz))
    else if (E%gridType == CORNER) then
	   E%v = cmplx(zr,zi)
    else
       write (0, *) 'gridType == ',trim(E%gridType),' undefined in random_cscalar'
    end if
    deallocate(zr,zi,STAT=istat)

  end subroutine random_cscalar

end module sg_scalar
