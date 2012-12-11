! *****************************************************************************
module sg_boundary
  ! This module creates data types for 2D vector fields defined on
  ! boundary nodes of a staggered Cartesian grid, along with basic
  ! algebraic (2D vector space) operations. Such operations include
  ! allocation, deallocation, initialization, copying, algebriac operations
  ! (linear combinations, scalar products, dot products) and other operations
  ! like working between boundary conditions and complex vector data
  ! structures. Belongs to SG_Basics class: staggered cartesian grid, data
  ! types defined on this grid, and operations defined on these data types. Not
  ! specific to EM problem, no dependency on outside (from other classes) modules.

  use math_constants
  use griddef
  use sg_vector
  use sg_vector_mg
  implicit none

  ! Generic interfaces are done through subroutines
  ! creates boundary nodes
  INTERFACE create
     module procedure create_cboundary
     module procedure create_cboundary_mg
  END INTERFACE

  ! deallocates the boundary nodes
  INTERFACE deall
     module procedure deall_cboundary
     module procedure deall_cboundary_mg
  END INTERFACE

  ! scalar multiplies the boundary nodes
  INTERFACE scMult
     module procedure scMult_cboundary
     module procedure scMult_cboundary_mg
  END INTERFACE

  INTERFACE linComb
     MODULE PROCEDURE linComb_cboundary
     module procedure linComb_cboundary_mg
  END INTERFACE

  ! adds the boundary nodes
  INTERFACE add
     module procedure add_cboundary
  END INTERFACE

  ! pointwise vector multiplication of boundary nodes
  INTERFACE diagMult
     module procedure diagMult_cboundary
  END INTERFACE

  ! zeros the boundary nodes
  INTERFACE zero
     module procedure zero_cboundary
  END INTERFACE

  INTERFACE dotProd
     MODULE PROCEDURE dotProd_cboundary
  END INTERFACE

  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE copy_cboundary
     module procedure b2bmg
  END INTERFACE

  ! Interface operators are done through functions
  INTERFACE OPERATOR (+)
     MODULE PROCEDURE add_cboundary_f
  END INTERFACE

  INTERFACE OPERATOR (*)
     MODULE PROCEDURE scMult_cboundary_f
     MODULE PROCEDURE diagMult_cboundary_f
  END INTERFACE

  ! gets boundary conditions from the vector field
  INTERFACE getBC
     module procedure copy_cbvector
  END INTERFACE

  ! sets boundary nodes in the vector field
  INTERFACE setBC
     !  module procedure copy_bcvector
     module procedure b2cmg
     module procedure bmg2cmg
  END INTERFACE

  public                :: copy_cbvector, copy_bcvector,b2cmg,bmg2cmg,b2bmg
  public				:: create_cboundary, create_cboundary_mg,   &
       deall_cboundary,deall_cboundary_mg, copy_cboundary, zero_cboundary, &
       scMult_cboundary,scMult_cboundary_mg, scMult_cboundary_f, add_cboundary, &
       add_cboundary_f, diagMult_cboundary, diagMult_cboundary_f, &
       dotProd_cboundary, &
       ! the two routines below are not part of any interface
       linComb_cboundary,linComb_cboundary_mg, scMultAdd_cboundary


  ! ***************************************************************************
  ! type cboundary defines complex 2D vector for boundary edge in a staggered
  ! grid
  type :: cboundary

     ! the boundary data are stored in two dimensional data structures
     ! xYMax are the boundary conditions for x component on the maximum Y end
     ! zYMax are the boundary conditions for z component on the maximum Y end
     ! xYMin are the boundary conditions for x component on the minimum Y end
     !   .... etc.
     ! the dimensions are:
     ! the dimensions for the YMax and YMin face is nx, nz+1 for x component for
     ! and nx+1, nz for z component
     ! the dimensions for the XMax and XMin face is ny, nz+1 for y component for
     ! and ny+1, nz for z component
     ! the dimensions for the ZMax and ZMin face is nx, ny+1 for x component for
     !  and nx+1, ny for z component

     complex (kind=prec), pointer, dimension(:,:)    :: xYMax, zYMax
     complex (kind=prec), pointer, dimension(:,:)    :: xYMin, zYMin
     complex (kind=prec), pointer, dimension(:,:)    :: yXMax, zXMax
     complex (kind=prec), pointer, dimension(:,:)    :: yXMin, zXMin
     complex (kind=prec), pointer, dimension(:,:)    :: xZMin, yZMin
     complex (kind=prec), pointer, dimension(:,:)    :: xZMax, yZMax

     ! Grid Dimensions:
     ! nx is grid dimension (number of cells) in the x-direction
     ! ny is grid dimension (number of cells) in the y-direction
     ! nz is grid dimension (number of cells) in the z-direction:
     integer                                          :: nx = 0, ny = 0, nz = 0

     ! allocated:  .true.  x, y, z arrays have been allocated
     logical		                              :: allocated = .false.
     ! a Logical variable used to determined if we are going to read an electrical field solution from a larger grid.
     logical                                      :: read_E_from_file=.false.
     ! pointer to parent grid
     type (grid_t), pointer                             :: grid

! *******************************************************************************************************************
  end type cboundary

  ! extension of the original cboundary type for the multi-grid
  ! needed for AdjBC
  type, extends(cboundary)            :: cboundary_mg
    ! bcArray contains BC in the sub-grid
    type(cboundary), pointer          :: bcArray(:)
    ! number of the sub-grids
    integer                           :: mgridSize
    ! coarseness
    integer, allocatable               :: coarseness(:)

  end type cboundary_mg
! *******************************************************************************************************************

Contains
  ! CREATE GRID_EDGE VECTORS
  ! * subroutine create_cboundary(igrid, E)

  ! DEALLOCATE GRID_EDGE VECTORS
  ! * subroutine deall_cboundary(E)

  ! COPY GRID_EDGE VECTORS :  (=)
  ! * subroutine copy_cboundary(E2,E1)

  ! ZERO GRID_EDGE VECTORS
  ! * subroutine zero_cboundary(E)

  ! SCALAR MULTIPLICATION :
  ! * subroutine scMult_cboundary(c, E1, E2)

  ! SCALAR MULTIPLICATION : FUNCTION VERSION (c*E)
  ! * function scMult_cboundary_f(c, E1) result(E2)

  ! VECTOR SUM :
  ! * subroutine add_cboundary(E1, E2, E3)

  ! VECTOR SUM : FUNCTION VERSION (E3 = E1+E2)
  ! * function add_cboundary_f(E1, E2) result(E3)

  ! POINTWISE MULTIPLICATION OF VECTORS:
  ! * subroutine diagMult_cboundary(E1, E2, E3)

  ! POINTWISE MULTIPLICATION OF VECTORS: FUNCTION VERSION (c*E)
  ! * function diagMult_cboundary_f(E1, E2) result(E3)

  ! VECTOR DOT PRODUCT
  ! * function dotProd_cboundary(E1, E2) result(c)

  ! COMPLEX LINEAR COMBINATION ... formally redundent with above
  ! only doing ones that we are sure to want
  ! * subroutine linCom_cboundary(inc1, E1, inc2, E2, E3)
  ! * subroutine scMultAdd_cboundary(c, E1, E2)  (E2 = E2+c*E1)

  ! The algebraic routines expect all input and output
  ! variables to be of the correct type, already allocated,
  ! and of the correct size.


  !****************************************************************************
  ! create_cboundary creates variable of derived type cboundary,
  ! using grid definition in structure "grid" ;
  ! allocates memory in different 2D faces of boundaries component arrays

  subroutine create_cboundary(igrid, E)

    implicit none
    type(grid_t), target, intent(in)     :: igrid
    ! the grid for which boundary edge fields are being initialized
    type (cboundary), intent(inout)     :: E
    integer                             :: status,nx,ny,nz

    if(E%allocated) then
       ! first deallocate memory for different 2D faces of boundaries
       deallocate(E%xYMax, E%zYMax, E%xYMin, E%zYMin, E%yXMax, &
            E%zXMax, E%yXMin, E%zXMin, E%xZMin, E%yZMin, &
            E%xZMax, E%yZMax, STAT=status)
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

    ! allocate memory for different 2D faces of boundaries;
    ! E%allocated will be true if all allocations succeed
    E%allocated = .true.
    allocate(E%xYMax(nx,nz+1), STAT=status)
    E%allocated = E%allocated .and.(status .EQ. 0)
    allocate(E%zYMax(nx+1,nz), STAT=status)
    E%allocated = E%allocated .and.(status .EQ. 0)
    allocate(E%xYMin(nx,nz+1), STAT=status)
    E%allocated = E%allocated .and.(status .EQ. 0)
    allocate(E%zYMin(nx+1,nz), STAT=status)
    E%allocated = E%allocated .and.(status .EQ. 0)
    allocate(E%yXMax(ny,nz+1), STAT=status)
    E%allocated = E%allocated .and.(status .EQ. 0)
    allocate(E%zXMax(ny+1,nz), STAT=status)
    E%allocated = E%allocated .and.(status .EQ. 0)
    allocate(E%yXMin(ny,nz+1), STAT=status)
    E%allocated = E%allocated .and.(status .EQ. 0)
    allocate(E%zXMin(ny+1,nz), STAT=status)
    E%allocated = E%allocated .and.(status .EQ. 0)
    allocate(E%xZMin(nx,ny+1), STAT=status)
    E%allocated = E%allocated .and.(status .EQ. 0)
    allocate(E%yZMin(nx+1,ny), STAT=status)
    E%allocated = E%allocated .and.(status .EQ. 0)
    allocate(E%xZMax(nx,ny+1), STAT=status)
    E%allocated = E%allocated .and.(status .EQ. 0)
    allocate(E%yZMax(nx+1,ny), STAT=status)
    E%allocated = E%allocated .and.(status .EQ. 0)

    if (E%allocated) then
       E%xYMax = C_ZERO
       E%zYMax = C_ZERO
       E%xYMin = C_ZERO
       E%zYMin = C_ZERO
       E%yXMax = C_ZERO
       E%zXMax = C_ZERO
       E%yXMin = C_ZERO
       E%zXMin = C_ZERO
       E%xZMin = C_ZERO
       E%yZMin = C_ZERO
       E%xZMax = C_ZERO
       E%yZMax = C_ZERO
    end if
    ! print *, 'E%allocated', E%allocated

  end subroutine create_cboundary  ! create_cboundary
  ! ***************************************************************************
  ! create_cboundary_mg creates variable of derived type cboundary_mg,
  ! using original procedure create_cboundary

  subroutine create_cboundary_mg(mgrid, E)

    implicit none
    type(grid_t), target, intent(in)     :: mgrid
    ! the grid for which boundary edge fields are being initialized
    type (cboundary_mg), intent(inout)   :: E
    !local
    integer                              :: imgrid

    ! first deallocate everything
    if(E%allocated)then
      call deall(E)
    endif

    ! allocate BC for sub-grids
    allocate(E%bcArray(mgrid%mgridSize))
    allocate(E%coarseness(mgrid%mgridSize))
    E%mgridSize = mgrid%mgridSize
    E%coarseness = mgrid%coarseness
    do imgrid = 1, E%mgridSize ! Loop over sub-grids
      call create(mgrid%gridArray(imgrid),E%bcArray(imgrid))
    enddo

    E%allocated = .true.

  end subroutine create_cboundary_mg
  !****************************************************************************
  ! deall_cboundary destoys variable of derived type cboundary,
  ! deallocating memory
  subroutine deall_cboundary(E)

    implicit none
    type (cboundary)  :: E
    integer	    :: status

    ! deallocate memory for different 2D faces of boundaries
    if(E%allocated) then
       deallocate(E%xYMax, E%zYMax, E%xYMin, E%zYMin, &
            E%yXMax, E%zXMax, E%yXMin, E%zXMin, &
            E%xZMin, E%yZMin, E%xZMax, E%yZMax, STAT=status)
    end if

    if(associated(E%grid)) nullify(E%grid)
    E%nx = 0
    E%ny = 0
    E%nz = 0
    E%allocated = .false.

  end subroutine deall_cboundary  ! deall_cboundary

! ****************************************************************************
  ! deall_cboundary_mg destroys variable of derived type cboundary_mg,
  ! deallocating memory
  subroutine deall_cboundary_mg(E)

    implicit none
    type (cboundary_mg)  :: E
    integer     :: status, imgrid

    do imgrid = 1, E%mgridSize
      ! deallocate memory for different 2D faces of boundaries
      if(E%allocated) then
        deallocate(E%bcArray(imgrid)%xYMax, E%bcArray(imgrid)%zYMax, E%bcArray(imgrid)%xYMin, E%bcArray(imgrid)%zYMin, &
            E%bcArray(imgrid)%yXMax, E%bcArray(imgrid)%zXMax, E%bcArray(imgrid)%yXMin, E%bcArray(imgrid)%zXMin, &
            E%bcArray(imgrid)%xZMin, E%bcArray(imgrid)%yZMin, E%bcArray(imgrid)%xZMax, E%bcArray(imgrid)%yZMax, STAT=status)
      end if

      if(associated(E%bcArray(imgrid)%grid)) nullify(E%bcArray(imgrid)%grid)

      E%bcArray(imgrid)%nx = 0
      E%bcArray(imgrid)%ny = 0
      E%bcArray(imgrid)%nz = 0
    enddo

    if(E%allocated) then
       deallocate(E%bcArray, STAT=status)
       deallocate(E%coarseness, STAT=status)
    end if

    E%mgridSize = 0
    E%allocated = .false.

  end subroutine deall_cboundary_mg  ! deall_cboundary

  !****************************************************************************
  ! copy_E nodeC makes an exact copy of derived data type
  ! cboundary;
  subroutine copy_cboundary(E2, E1)

    implicit none
    type (cboundary), intent(in)            :: E1
    type (cboundary), intent(inout)         :: E2
    integer                               :: status

    ! check to see if RHS (E1) is active (allocated)
    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for copy_cboundary'
    else

       if((E2%nx == E1%nx).and.(E2%ny == E1%ny).and.(E2%nz == E1%nz)) then
         ! just copy components
          E2%xYMax = E1%xYMax
          E2%zYMax = E1%zYMax
          E2%xYMin = E1%xYMin
          E2%zYMin = E1%zYMin
          E2%yXMax = E1%yXMax
          E2%zXMax = E1%zXMax
          E2%yXMin = E1%yXMin
          E2%zXMin = E1%zXMin
          E2%xZMin = E1%xZMin
          E2%yZMin = E1%yZMin
          E2%xZMax = E1%xZMax
          E2%yZMax = E1%yZMax
       else
          if(E2%allocated) then
             ! first deallocate memory for different 2D faces of boundaries
             deallocate(E2%xYMax, E2%zYMax, E2%xYMin, E2%zYMin, E2%yXMax, &
                  E2%zXMax, E2%yXMin, E2%zXMin, E2%xZMin, E2%yZMin, &
                  E2%xZMax, E2%yZMax, STAT=status)
          end if
          !  then allocate E2 as correct size ...
          Call create(E1%grid, E2) !create_cboundary
          !   .... and copy E1
          E2%xYMax = E1%xYMax
          E2%zYMax = E1%zYMax
          E2%xYMin = E1%xYMin
          E2%zYMin = E1%zYMin
          E2%yXMax = E1%yXMax
          E2%zXMax = E1%zXMax
          E2%yXMin = E1%yXMin
          E2%zXMin = E1%zXMin
          E2%xZMin = E1%xZMin
          E2%yZMin = E1%yZMin
          E2%xZMax = E1%xZMax
          E2%yZMax = E1%yZMax
          E2%read_E_from_file = E1%read_E_from_file

       end if

    end if

  end subroutine copy_cboundary  ! copy_cboundary


  !****************************************************************************
  ! zero_E nodeC zeros variable of derived data type
  ! cboundary;

  subroutine zero_cboundary(E)

    implicit none
    type (cboundary), intent(inout) :: E

    ! check to see if E is active (allocated)
    if(.not.E%allocated) then
       write(0,*) 'Error in zero_cboundary: E not allocated'
    else

       E%xYMax = C_ZERO
       E%zYMax = C_ZERO
       E%xYMin = C_ZERO
       E%zYMin = C_ZERO
       E%yXMax = C_ZERO
       E%zXMax = C_ZERO
       E%yXMin = C_ZERO
       E%zXMin = C_ZERO
       E%xZMin = C_ZERO
       E%yZMin = C_ZERO
       E%xZMax = C_ZERO
       E%yZMax = C_ZERO

    end if

  end subroutine zero_cboundary ! zero_cboundary


  !****************************************************************************
  ! scMult_cboundary multiplies vector stored as derived data type
  ! cboundary with a complex scalar; subroutine version
  ! E2 can overwrite E1
  subroutine scMult_cboundary(c, E1, E2)

    implicit none
    complex(kind=prec), intent(in)                      :: c
    ! a complex scalar to be multiplied with
    type (cboundary), intent(in)                       :: E1
    type (cboundary), intent(inout)                    :: E2

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for ScMult_cboundary'
       stop
    endif

    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated yet for ScMult_cboundary'
    else

       ! Check whether all the boundary nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then
             ! complex scalar multiplication for different 2D faces of
             ! boundaries components
          E2%xYMax = E1%xYMax * c
          E2%zYMax = E1%zYMax * c
          E2%xYMin = E1%xYMin * c
          E2%zYMin = E1%zYMin * c
          E2%yXMax = E1%yXMax * c
          E2%zXMax = E1%zXMax * c
          E2%yXMin = E1%yXMin * c
          E2%zXMin = E1%zXMin * c
          E2%xZMin = E1%xZMin * c
          E2%yZMin = E1%yZMin * c
          E2%xZMax = E1%xZMax * c
          E2%yZMax = E1%yZMax * c
       else
          write(0, *) 'Error:scMult_cboundary: vectors not same size'
       end if
    end if

  end subroutine scMult_cboundary

  !****************************************************************************
  ! scMult_cboundary multiplies vector stored as derived data type
  ! cboundary with a complex scalar; subroutine version
  ! E2 can overwrite E1
  subroutine scMult_cboundary_mg(c, E1, E2)

    implicit none
    complex(kind=prec), intent(in)                      :: c
    ! a complex scalar to be multiplied with
    type (cboundary_mg), intent(in)                     :: E1
    type (cboundary_mg), intent(inout)                  :: E2
    ! local
    integer                                             :: imgrid

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for ScMult_cboundary'
       stop
    endif

    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated yet for ScMult_cboundary_mg'
    else

      if(E1%mgridSize == E2%mgridSize) then

        do imgrid = 1, E1%mgridSize

          ! Check whether all the boundary nodes are of the same size
          if((E1%bcArray(imgrid)%nx == E2%bcArray(imgrid)%nx).and.(E1%bcArray(imgrid)%ny == E2%bcArray(imgrid)%ny)&
                .and.(E1%bcArray(imgrid)%nz == E2%bcArray(imgrid)%nz)) then
             ! complex scalar multiplication for different 2D faces of
             ! boundaries components
             E2%bcArray(imgrid)%xYMax = E1%bcArray(imgrid)%xYMax * c
             E2%bcArray(imgrid)%zYMax = E1%bcArray(imgrid)%zYMax * c
             E2%bcArray(imgrid)%xYMin = E1%bcArray(imgrid)%xYMin * c
             E2%bcArray(imgrid)%zYMin = E1%bcArray(imgrid)%zYMin * c
             E2%bcArray(imgrid)%yXMax = E1%bcArray(imgrid)%yXMax * c
             E2%bcArray(imgrid)%zXMax = E1%bcArray(imgrid)%zXMax * c
             E2%bcArray(imgrid)%yXMin = E1%bcArray(imgrid)%yXMin * c
             E2%bcArray(imgrid)%zXMin = E1%bcArray(imgrid)%zXMin * c
             E2%bcArray(imgrid)%xZMin = E1%bcArray(imgrid)%xZMin * c
             E2%bcArray(imgrid)%yZMin = E1%bcArray(imgrid)%yZMin * c
             E2%bcArray(imgrid)%xZMax = E1%bcArray(imgrid)%xZMax * c
             E2%bcArray(imgrid)%yZMax = E1%bcArray(imgrid)%yZMax * c
           else
             print*, 'Error:scMult_cboundary_mg: vectors not same size'
          end if
        enddo
      else
        print*, 'Error:scMult_cboundary_mg: vectors not same size'
      endif
    end if

  end subroutine scMult_cboundary_mg

  !****************************************************************************
  ! scMult_cboundary_f multiplies vector stored as derived data type
  ! cboundary with a complex scalar; function version
  function scMult_cboundary_f(c, E1) result(E2)

    implicit none
    complex(kind=prec), intent(in)                      :: c
    ! a complex scalar to be multiplied with
    type (cboundary), intent(in)                       :: E1
    type (cboundary)                                   :: E2

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for ScMult_cboundary_f'
       stop
    endif

    ! In function version, appropriate data types need to be created
    Call create_cboundary(E1%grid, E2)
    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated for scMult_cboundary_f'
    else

       ! Check whether all the boundary nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then
          ! complex scalar multiplication for different 2D faces of
          ! boundaries components
          E2%xYMax = E1%xYMax * c
          E2%zYMax = E1%zYMax * c
          E2%xYMin = E1%xYMin * c
          E2%zYMin = E1%zYMin * c
          E2%yXMax = E1%yXMax * c
          E2%zXMax = E1%zXMax * c
          E2%yXMin = E1%yXMin * c
          E2%zXMin = E1%zXMin * c
          E2%xZMin = E1%xZMin * c
          E2%yZMin = E1%yZMin * c
          E2%xZMax = E1%xZMax * c
          E2%yZMax = E1%yZMax * c
       else
          write(0, *) 'Error:scMult_cboundary_f: vectors not same size'
       end if
    end if

  end function scMult_cboundary_f

  !****************************************************************************
  ! add_cboundary adds vectors stored as derived data type
  ! cboundary with ; subroutine version
  ! E3 can overwrite E1 and E2
  subroutine add_cboundary(E1, E2, E3)

    implicit none
    type (cboundary), intent(in)               :: E1, E2
    type (cboundary), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for add_cboundary'
       stop
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS not allocated for add_cboundary'
    else

       ! Check whether all the boundary nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
          (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then
          ! add different 2D faces of boundaries components
          E3%xYMax = E1%xYMax + E2%xYMax
          E3%zYMax = E1%zYMax + E2%zYMax
          E3%xYMin = E1%xYMin + E2%xYMin
          E3%zYMin = E1%zYMin + E2%zYMin
          E3%yXMax = E1%yXMax + E2%yXMax
          E3%zXMax = E1%zXMax + E2%zXMax
          E3%yXMin = E1%yXMin + E2%yXMin
          E3%zXMin = E1%zXMin + E2%zXMin
          E3%xZMin = E1%xZMin + E2%xZMin
          E3%yZMin = E1%yZMin + E2%yZMin
          E3%xZMax = E1%xZMax + E2%xZmax
          E3%yZMax = E1%yZMax + E2%yZmax
       else
          write(0, *) 'Error:add_cboundary: vectors not same size'
       end if
    end if

  end subroutine add_cboundary ! add_cboundary


  !****************************************************************************
  ! add_cboundary_f adds vectors stored as derived data type
  ! cboundary with ; function version
  function add_cboundary_f(E1, E2) result(E3)

    implicit none
    type (cboundary), intent(in)               :: E1, E2
    type (cboundary)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for add_cboundary_f'
       stop
    endif

    ! In function version, appropriate data types need to be created
    Call create(E1%grid, E3) !create_cboundary
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for add_cboundary_f'
    else

       ! Check whether all the boundary nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then
          ! add different 2D faces of boundaries components
          E3%xYMax = E1%xYMax + E2%xYMax
          E3%zYMax = E1%zYMax + E2%zYMax
          E3%xYMin = E1%xYMin + E2%xYMin
          E3%zYMin = E1%zYMin + E2%zYMin
          E3%yXMax = E1%yXMax + E2%yXMax
          E3%zXMax = E1%zXMax + E2%zXMax
          E3%yXMin = E1%yXMin + E2%yXMin
          E3%zXMin = E1%zXMin + E2%zXMin
          E3%xZMin = E1%xZMin + E2%xZMin
          E3%yZMin = E1%yZMin + E2%yZMin
          E3%xZMax = E1%xZMax + E2%xZmax
          E3%yZMax = E1%yZMax + E2%yZmax
       else
          write(0, *) 'Error:add_cboundary_f: vectors not same size'
       end if
    end if
  end function add_cboundary_f


  !****************************************************************************
  ! diagMult_cboundary multiplies two vectors E1, E2 stored as derived data
  ! type cboundary pointwise; subroutine version
  ! E3 can overwrite E1 or E2
  subroutine diagMult_cboundary(E1, E2, E3)

    implicit none
    type (cboundary), intent(in)               :: E1, E2
    type (cboundary), intent(inout)            :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_cboundary'
       stop
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_cboundary'
    else
       ! Check whether all the boundary nodes are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
         (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then
          ! pointwise multiplication for different 2D faces of
          ! boundaries components
          E3%xYMax = E1%xYMax * E2%xYMax
          E3%zYMax = E1%zYMax * E2%zYMax
          E3%xYMin = E1%xYMin * E2%xYMin
          E3%zYMin = E1%zYMin * E2%zYMin
          E3%yXMax = E1%yXMax * E2%yXMax
          E3%zXMax = E1%zXMax * E2%zXMax
          E3%yXMin = E1%yXMin * E2%yXMin
          E3%zXMin = E1%zXMin * E2%zXMin
          E3%xZMin = E1%xZMin * E2%xZMin
          E3%yZMin = E1%yZMin * E2%yZMin
          E3%xZMax = E1%xZMax * E2%xZmax
          E3%yZMax = E1%yZMax * E2%yZmax
       else
          write(0, *) 'Error:diagMult_cboundary: vectors not same size'
       end if
    end if
  end subroutine diagMult_cboundary


  !****************************************************************************
  ! diagMult_cboundary_f multiplies two vectors E1, E2 stored as derived
  ! data  type cboundary pointwise; function version
  function diagMult_cboundary_f(E1, E2) result(E3)

    implicit none
    type (cboundary), intent(in)               :: E1, E2
    type (cboundary)                           :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_cboundary_f'
       stop
    endif

    ! In function version, appropriate data types need to be created
    Call create_cboundary(E1%grid, E3)
    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for diagMult_cboundary_f'
    else
       ! Check whether both vectors are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
            (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then
          ! pointwise multiplication for different 2D faces of
          ! boundaries components
          E3%xYMax = E1%xYMax * E2%xYMax
          E3%zYMax = E1%zYMax * E2%zYMax
          E3%xYMin = E1%xYMin * E2%xYMin
          E3%zYMin = E1%zYMin * E2%zYMin
          E3%yXMax = E1%yXMax * E2%yXMax
          E3%zXMax = E1%zXMax * E2%zXMax
          E3%yXMin = E1%yXMin * E2%yXMin
          E3%zXMin = E1%zXMin * E2%zXMin
          E3%xZMin = E1%xZMin * E2%xZMin
          E3%yZMin = E1%yZMin * E2%yZMin
          E3%xZMax = E1%xZMax * E2%xZmax
          E3%yZMax = E1%yZMax * E2%yZmax
       else
          write(0, *) 'Error:diagMult_cboundary_f: vectors not same size'
       end if
    end if
  end function diagMult_cboundary_f

  !****************************************************************************
  ! dotProd_cboundary computes dot product of two vecors stored
  ! as derived data type cboundary, returning a complex number
  function dotProd_cboundary(E1, E2) result(c)

    implicit none
    type (cboundary), intent(in)       :: E1, E2
    complex(kind=prec)		       :: c

    c = C_ZERO

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for dotProd_cboundary'
       stop
    endif

    ! Check whether both input vectors are of the same size
    if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then
       c = c + sum(conjg(E1%xYMax) * E2%xYMax)
       c = c + sum(conjg(E1%zYMax) * E2%zYMax)
       c = c + sum(conjg(E1%xYMin) * E2%xYMin)
       c = c + sum(conjg(E1%zYMin) * E2%zYMin)
       c = c + sum(conjg(E1%yXMax) * E2%yXMax)
       c = c + sum(conjg(E1%zXMax) * E2%zXMax)
       c = c + sum(conjg(E1%yXMin) * E2%yXMin)
       c = c + sum(conjg(E1%zXMin) * E2%zXMin)
       c = c + sum(conjg(E1%xZMin) * E2%xZMin)
       c = c + sum(conjg(E1%yZMin) * E2%yZMin)
       c = c + sum(conjg(E1%xZMax) * E2%xZmax)
       c = c + sum(conjg(E1%yZMax) * E2%yZmax)
    else
       write(0, *) 'Error:dotProd_cboundary: vectors not same size'
    end if
  end function dotProd_cboundary

  !****************************************************************************
  ! linComb_cboundary computes linear combination of two vectors
  ! stored as derived data type cboundary; subroutine, not a function
  ! both input vectors must have the same dimension
  subroutine linComb_cboundary(inc1, E1, inc2, E2, E3)

    implicit none
    !   input vectors
    type (cboundary), intent(in)             :: E1, E2
    !  input complex scalars
    complex (kind=8), intent(in)             :: inc1, inc2
    type (cboundary), intent(inout)          :: E3

    if((.not.E1%allocated).or.(.not.E2%allocated)) then
       write(0,*) 'RHS not allocated yet for linComb_cboundary'
       stop
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.E3%allocated) then
       write(0,*) 'LHS was not allocated for linComb_cboundary'
    else

       ! Check whether all vectors are of the same size
       if ((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz).and.&
         (E1%nx == E3%nx).and.(E1%ny == E3%ny).and.(E1%nz == E3%nz)) then
          ! form linear combinatoin
          E3%xYMax = inc1*E1%xYMax + inc2*E2%xYMax
          E3%zYMax = inc1*E1%zYMax + inc2*E2%zYMax
          E3%xYMin = inc1*E1%xYMin + inc2*E2%xYMin
          E3%zYMin = inc1*E1%zYMin + inc2*E2%zYMin
          E3%yXMax = inc1*E1%yXMax + inc2*E2%yXMax
          E3%zXMax = inc1*E1%zXMax + inc2*E2%zXMax
          E3%yXMin = inc1*E1%yXMin + inc2*E2%yXMin
          E3%zXMin = inc1*E1%zXMin + inc2*E2%zXMin
          E3%xZMin = inc1*E1%xZMin + inc2*E2%xZMin
          E3%yZMin = inc1*E1%yZMin + inc2*E2%yZMin
          E3%xZMax = inc1*E1%xZMax + inc2*E2%xZmax
          E3%yZMax = inc1*E1%yZMax + inc2*E2%yZmax
       else
          write(0, *) 'Error:linComb_cboundary:  vectors not same size'
       end if
    end if

  end subroutine linComb_cboundary
  !****************************************************************************
  ! linComb_cboundary computes linear combination of two vectors
  ! stored as derived data type cboundary; subroutine, not a function
  ! both input vectors must have the same dimension
  subroutine linComb_cboundary_mg(inc1, E1, inc2, E2, E3)

    implicit none
    !   input vectors
    type (cboundary_mg), intent(in)             :: E1
    type (cboundary), intent(in)                :: E2
    !  input complex scalars
    complex (kind=8), intent(in)                :: inc1, inc2
    type (cboundary_mg), intent(inout)          :: E3
    ! local
    type(cboundary_mg)                          :: temp
    integer                                     :: imgrid

    call create(E1%grid,temp)
    temp = E2   ! b2bmg

    if(E1%mgridSize == temp%mgridSize .and. E1%mgridSize == E3%mgridSize) then
      do imgrid = 1, E1%mgridsize
        call linComb_cboundary(inc1,E1%bcArray(imgrid), inc2, temp%bcArray(imgrid), E3%bcArray(imgrid))
      enddo
    else
      print*, 'Error:linComb_cboundary_mg:  vectors not same size'
    end if

    call deall(temp)

  end subroutine linComb_cboundary_mg

  !****************************************************************************
  ! scMultadd_cboundary multiplies vector E1 stored as derived data type
  ! cboundary with a complex scalar c, adding result to output vector E2
  subroutine scMultAdd_cboundary(c, E1, E2)

    implicit none
    complex(kind=prec), intent(in)                        :: c
    ! a complex scalar to be multiplied with
    type (cboundary), intent(in)                       :: E1
    type (cboundary), intent(inout)                    :: E2

    if(.not.E1%allocated) then
       write(0,*) 'RHS not allocated yet for scMultAdd_cboundary'
       stop
    endif

    ! check to see if LHS (E2) is active (allocated)
    if(.not.E2%allocated) then
       write(0,*) 'LHS was not allocated for scMultAdd_cboundary'
    else

       ! Check whether both vectors are of the same size
       if((E1%nx == E2%nx).and.(E1%ny == E2%ny).and.(E1%nz == E2%nz)) then
             ! complex scalar multiplication for different 2D faces of
             ! boundaries components
          E2%xYMax = E2%xYMax + E1%xYMax * c
          E2%zYMax = E2%zYMax + E1%zYMax * c
          E2%xYMin = E2%xYMin + E1%xYMin * c
          E2%zYMin = E2%zYMin + E1%zYMin * c
          E2%yXMax = E2%yXMax + E1%yXMax * c
          E2%zXMax = E2%zXMax + E1%zXMax * c
          E2%yXMin = E2%yXMin + E1%yXMin * c
          E2%zXMin = E2%zXMin + E1%zXMin * c
          E2%xZMin = E2%xZMin + E1%xZMin * c
          E2%yZMin = E2%yZMin + E1%yZMin * c
          E2%xZMax = E2%xZMax + E1%xZmax * c
          E2%yZMax = E2%yZMax + E1%yZmax * c
       else
          write(0, *) 'Error:scMultAddy_cboundary: vectors not same size'
       end if
    end if

  end subroutine scMultAdd_cboundary

  ! ***************************************************************************
  ! * copy_cbvector reads the reads a complex vector used for EDGE as an input and
  ! * extracts the boundary condition values and writes them into boundary
  ! * conditions derived data types as an output for a specific frequency/
  ! * period and mode.
  subroutine copy_cbvector(inE, outBC)

    implicit none
    type (cvector),target, intent(in)              :: inE
    ! the electrical field as an input
    type (cboundary), intent(inout)                :: outBC
    ! boundary conditions as an output
    integer                      :: ix, iy, iz        ! dummy integers
    integer                      :: YMax, YMin        ! ends for BC extraction
    integer                      :: XMax, XMin        ! ends for BC extraction
    integer                      :: ZMin, ZMax        ! ends for BC extraction

    if (.not.inE%allocated) then
      write(0,*) 'inE in copy_cbvector not allocated yet'
      stop
    end if

    if (.not.outBC%allocated) then
      write(0,*) 'outBC in copy_cbvector not allocated yet'
      stop
    end if

    ! Check whether the bounds are the same
    if ((inE%nx == outBC%nx).and.&
         (inE%ny == outBC%ny).and.&
         (inE%nz == outBC%nz)) then

       if (inE%gridType == EDGE) then

          YMax = inE%ny+1
          YMin = 1
          XMax = inE%nx+1
          XMin = 1
          ZMin = 1
          ZMax = inE%nz+1

          ! extracting boundary conditions from the electrical field
          ! extracting boundary values for the YMin and YMax face
          ! x-component
          do ix = 1, inE%nx
             do iz = 1, inE%nz+1
                ! Extracting the YMin face
                outBC%xYMin(ix, iz) = inE%x(ix, YMin, iz)
                ! Extracting the YMax face
                outBC%xYMax(ix, iz) = inE%x(ix, YMax, iz)
             enddo ! iz
          enddo    ! ix
          ! z-component
          do ix = 1, inE%nx+1
             do iz = 1, inE%nz
                ! Extracting the YMin face
                outBC%zYMin(ix, iz) = inE%z(ix, YMin, iz)
                ! Extracting the YMax face
                outBC%zYMax(ix, iz) = inE%z(ix, YMax, iz)
             enddo  ! iz
          enddo     ! ix

          ! extracting boundary values for the XMax and XMin face
          ! y-component
          do iy = 1, inE%ny
             do iz = 1, inE%nz+1
                ! Extracting the XMax face
                outBC%yXMax(iy, iz) = inE%y(XMax, iy, iz)
                ! Extracting the XMin face
                outBC%yXMin(iy, iz) = inE%y(XMin, iy, iz)
             enddo  ! iy
          enddo     ! iz
          ! z-component
          do iy = 1, inE%ny+1
             do iz = 1, inE%nz
                ! Extracting the XMax face
                outBC%zXMax(iy, iz) = inE%z(XMax, iy, iz)
                ! Extracting the XMin face
                outBC%zXMin(iy, iz) = inE%z(XMin, iy, iz)
             enddo  ! iy
          enddo     ! iz

          ! extracting boundary values for the ZMin and ZMax face
          ! x-component
          do ix = 1, inE%nx
             do iy = 1, inE%ny+1
                ! Extracting the ZMin face
                outBC%xZMin(ix, iy) = inE%x(ix, iy, ZMin)
                ! Extracting the ZMax face
                outBC%xZMax(ix, iy) = inE%x(ix, iy, ZMax)
             enddo   ! iy
          enddo      ! ix
          ! y-component
          do ix = 1, inE%nx+1
             do iy = 1, inE%ny
                ! Extracting the ZMin face
                outBC%yZMin(ix, iy) = inE%y(ix, iy, ZMin)
                ! Extracting the ZMax face
                outBC%yZMax(ix, iy) = inE%y(ix, iy, ZMax)
             enddo   ! ix
          enddo      ! iy

       else
          write (0, *) 'not compatible usage for existing data types'
       end if

    else
       write(0, *) 'Error-complex vector and BC not same size'
    endif

  end subroutine copy_cbvector ! copy_cbvector

  ! ***************************************************************************
  ! * copy_bcvector reads the boundary conditions as an input and writes them to the
  ! * complex vector used as EDGE on the designated boundaries for a
  ! * specific frequency/ period and mode.
  subroutine copy_bcvector(inBC, outE)

    implicit none
    type (cboundary), intent(in)     :: inBC
    ! boundary conditions as an input
    type (cvector), intent(inout)    :: outE
    ! the electrical field as an output
    integer                      :: ix, iy, iz       ! dummy integers
    integer                      :: YMax, YMin       ! ends for BC writing
    integer                      :: XMax, XMin       ! ends for BC writing
    integer                      :: ZMin, ZMax       ! ends for BC writing

    if (.not.outE%allocated) then
      write(0,*) 'outE in copy_bcvector not allocated yet'
      stop
    end if

    if (.not.inBC%allocated) then
      write(0,*) 'inBC in copy_bcvector not allocated yet'
      stop
    end if

    ! Check whether the bounds are the same
    if ((inBC%nx == outE%nx).and.&
         (inBC%ny == outE%ny).and.&
         (inBC%nz == outE%nz)) then

       if (outE%gridType == EDGE) then

          YMax = outE%ny+1
          YMin = 1
          XMax = outE%nx+1
          XMin = 1
          ZMin = 1
          ZMax = outE%nz+1

          ! passing boundary condition values to the complete description
          ! of the vector field by  passing boundary values for the YMin
          ! and YMax face
          ! x-component
          do ix = 1, outE%nx
             do iz = 1, outE%nz+1
                ! Passing the YMin face
                outE%x(ix, YMin, iz) = inBC%xYMin(ix, iz)
                ! Passing the YMax face
                outE%x(ix, YMax, iz) = inBC%xYMax(ix, iz)
             enddo ! iz
          enddo    ! ix
          ! z-component
          do ix = 1, outE%nx
             do iz = 1, outE%nz
                ! Passing the YMin face
                outE%z(ix, YMin, iz) = inBC%zYMin(ix, iz)
                ! Passing the YMax face
                outE%z(ix, YMax, iz) = inBC%zYMax(ix, iz)
             enddo  ! iz
          enddo     ! ix

          ! passing boundary values for the XMax and XMin face
          ! y-component
          do iy = 1, outE%ny
             do iz = 1, outE%nz+1
                ! Passing the XMax face
                outE%y(XMax, iy, iz) = inBC%yXMax(iy, iz)
                ! Passing the XMin face
                outE%y(XMin, iy, iz) = inBC%yXMin(iy, iz)
             enddo  ! iy
          enddo     ! iz
          ! z-component
          do iy = 1, outE%ny+1
             do iz = 1, outE%nz
                ! Passing the XMax face
                outE%z(XMax, iy, iz) = inBC%zXMax(iy, iz)
                ! Passing the XMin face
                outE%z(XMin, iy, iz) = inBC%zXMin(iy, iz)
             enddo  ! iy
          enddo     ! iz

          ! passing boundary values for the ZMin and ZMax face
          ! x-component
          do ix = 1, outE%nx
             do iy = 1, outE%ny+1
                ! Passing the ZMin face
                outE%x(ix, iy, ZMin) = inBC%xZMin(ix, iy)
                ! Passing the ZMax face
                outE%x(ix, iy, ZMax) = inBC%xZMax(ix, iy)
             enddo   ! iy
          enddo      ! ix
          ! y-component
          do ix = 1, outE%nx+1
             do iy = 1, outE%ny
                ! Passing the ZMin face
                outE%y(ix, iy, ZMin) = inBC%yZMin(ix, iy)
                ! Passing the ZMax face
                outE%y(ix, iy, ZMax) = inBC%yZMax(ix, iy)
             enddo   ! ix
          enddo      ! iy

       else
          write (0, *) 'copy_bcvector: not compatible usage for existing data types'
       end if

    else
       write(0, *) 'Error-complex vector and BC not same size'
    endif

  end subroutine copy_bcvector

  ! ***************************************************************************
  ! * bmg2cmg reads the boundary conditions defined on the multi-grid as an input
  !   and writes them to the complex vector on the multi-grid
  ! created by Cherevatova(Nov,2012)

  subroutine bmg2cmg(inBC, outE)

    implicit none
    type (cboundary_mg), intent(in)     :: inBC
    ! boundary conditions as an input
    type (cvector_mg), intent(inout)    :: outE
    ! the electrical field as an output
    integer                      :: ix, iy, iz, imgrid       ! dummy integers
    integer                      :: YMax, YMin               ! ends for BC writing
    integer                      :: XMax, XMin               ! ends for BC writing
    integer                      :: ZMin, ZMax               ! ends for BC writing

    if (.not.outE%allocated) then
      print*,  'outE in bmg2cmg not allocated yet'
      stop
    end if

    if (.not.inBC%allocated) then
      print*,  'inBC in copy_bcvector not allocated yet'
      stop
    end if

    ! Check whether the bounds are the same
    if(inBC%mgridSize == outE%mgridSize) then

      if (outE%gridType == EDGE) then
        XMin = 1
        YMin = 1
        do imgrid = 1, inBC%mgridSize  ! Global loop over sub-grids

          if ((inBC%bcArray(imgrid)%nx == outE%cvArray(imgrid)%nx).and. &
             (inBC%bcArray(imgrid)%ny == outE%cvArray(imgrid)%ny).and.&
             (inBC%bcArray(imgrid)%nz == outE%cvArray(imgrid)%nz)) then

              YMax = outE%cvArray(imgrid)%ny+1
              XMax = outE%cvArray(imgrid)%nx+1

              ! passing boundary condition values to the complete description
              ! of the vector field by  passing boundary values for the YMin
              ! and YMax face
              ! x-component
              do ix = 1, outE%cvArray(imgrid)%nx
                 do iz = 1, outE%cvArray(imgrid)%nz+1
                    ! Passing the YMin face
                    outE%cvArray(imgrid)%x(ix, YMin, iz) = inBC%bcArray(imgrid)%xYMin(ix, iz)
                    ! Passing the YMax face
                    outE%cvArray(imgrid)%x(ix, YMax, iz) = inBC%bcArray(imgrid)%xYMax(ix, iz)
                 enddo ! iz
              enddo    ! ix
              ! z-component
              do ix = 1, outE%cvArray(imgrid)%nx
                 do iz = 1, outE%cvArray(imgrid)%nz
                    ! Passing the YMin face
                    outE%cvArray(imgrid)%z(ix, YMin, iz) = inBC%bcArray(imgrid)%zYMin(ix, iz)
                    ! Passing the YMax face
                    outE%cvArray(imgrid)%z(ix, YMax, iz) = inBC%bcArray(imgrid)%zYMax(ix, iz)
                 enddo  ! iz
              enddo     ! ix

              ! passing boundary values for the XMax and XMin face
              ! y-component
              do iy = 1, outE%cvArray(imgrid)%ny
                 do iz = 1, outE%cvArray(imgrid)%nz+1
                    ! Passing the XMax face
                    outE%cvArray(imgrid)%y(XMax, iy, iz) = inBC%bcArray(imgrid)%yXMax(iy, iz)
                    ! Passing the XMin face
                    outE%cvArray(imgrid)%y(XMin, iy, iz) = inBC%bcArray(imgrid)%yXMin(iy, iz)
                 enddo  ! iy
              enddo     ! iz
              ! z-component
              do iy = 1, outE%cvArray(imgrid)%ny+1
                 do iz = 1, outE%cvArray(imgrid)%nz
                    ! Passing the XMax face
                    outE%cvArray(imgrid)%z(XMax, iy, iz) = inBC%bcArray(imgrid)%zXMax(iy, iz)
                    ! Passing the XMin face
                    outE%cvArray(imgrid)%z(XMin, iy, iz) = inBC%bcArray(imgrid)%zXMin(iy, iz)
                 enddo  ! iy
              enddo     ! iz
           else
             print*, 'bmg2cmg: error-complex vector and BC not same size'
           end if
        enddo  ! Global loop over sub-grids

        ZMin = 1
        ZMax = outE%cvArray(outE%mgridSize)%nz+1
        ! passing boundary values for the ZMin and ZMax face
        ! x-component
        do ix = 1, outE%cvArray(1)%nx
          do iy = 1, outE%cvArray(1)%ny+1
           ! Passing the ZMin face
            outE%cvArray(1)%x(ix, iy, ZMin) = inBC%bcArray(1)%xZMin(ix, iy)
          enddo
        enddo
        ! y-component
        do ix = 1, outE%cvArray(1)%nx+1
          do iy = 1, outE%cvArray(1)%ny
           ! Passing the ZMin face
           outE%cvArray(1)%y(ix, iy, ZMin) = inBC%bcArray(1)%yZMin(ix, iy)
          enddo
        enddo
        ! x-component
        do ix = 1, outE%cvArray(outE%mgridSize)%nx
          do iy = 1, outE%cvArray(outE%mgridSize)%ny+1
            !  Passing the ZMax face
            outE%cvArray(outE%mgridSize)%x(ix, iy, ZMax) = inBC%bcArray(outE%mgridSize)%xZMax(ix, iy)
          enddo   ! iy
        enddo      ! ix
        ! y-component
        do ix = 1, outE%cvArray(outE%mgridSize)%nx+1
          do iy = 1, outE%cvArray(outE%mgridSize)%ny
           ! Passing the ZMax face
           outE%cvArray(outE%mgridSize)%y(ix, iy, ZMax) = inBC%bcArray(outE%mgridSize)%yZMax(ix, iy)
          enddo   ! ix
        enddo      ! iy

      else
        print*, 'bmg2cmg: not compatible usage for existing data types'
       end if

    else
       print*, 'Error-complex vector and BC not same size'
    endif

  end subroutine bmg2cmg

  ! ***************************************************************************
  ! * copy_bcvector reads the boundary conditions as an input and writes them to the
  ! * complex vector used as EDGE on the designated boundaries for a
  ! * specific frequency/ period and mode.
  ! * editied for the case of cvector_mg !!! Multi-grid
  subroutine b2cmg(inBC, outE)

    implicit none
    ! boundary conditions as an input
    type (cboundary), intent(in)  :: inBC
    ! the electrical field as an output
    type (cvector_mg), intent(inout)    :: outE

    ! local
    type(cboundary_mg)            :: temp

    if (.not.outE%allocated) then
      print *, 'outE in copy_bcvector_mg not allocated yet'
      stop
    end if

    if (.not.inBC%allocated) then
      print *, 'inBC in copy_bcvector_mg not allocated yet'
      stop
    end if

    call create(outE%grid,temp)
    temp=inBC
    call bmg2cmg(temp,outE)
    ! deallocate temp array
    call deall(temp)

end subroutine b2cmg

  ! ***************************************************************************
  ! *copy cboundary to cboundary_mg
  subroutine b2bmg(outBC, inBC)
  ! created by Cherevatova (Nov, 2012)

    implicit none
    ! boundary conditions as an input
    type (cboundary), intent(in)          :: inBC
    ! averaged boundary conditions for the multi-grid as an input
    type (cboundary_mg), intent(inout)    :: outBC
    ! local
    integer   :: ix, iy, iz, izv, ic, imgrid    ! dummy integers
    integer   :: YMax, YMin       ! ends for BC writing
    integer   :: XMax, XMin       ! ends for BC writing
    integer   :: ZMin, ZMax       ! ends for BC writing
    integer   :: ccoeff_current, nzCum
    integer  :: errAll

    if (.not.outBC%allocated) then
      call create(inBC%grid,outBC)
    end if

    if (.not.inBC%allocated) then
      print *, 'inBC in b2mg not allocated yet'
      stop
    end if

      YMin = 1
      XMin = 1
      ZMin = 1
      ZMax = outBC%bcarray(outBC%mgridSize)%nz+1
      nzCum = 0
     ! Check whether the bounds are the same
     do imgrid = 1, outBC%mgridSize

        YMax = outBC%bcArray(imgrid)%ny+1
        XMax = outBC%bcArray(imgrid)%nx+1
        ccoeff_current = 2**outBC%coarseness(imgrid)

       ! passing boundary condition values to the complete description
       ! of the vector field by  passing boundary values for the YMin
       ! and YMax face
       ! x-component
       do iz = 1, outBC%bcArray(imgrid)%nz+1
          izv = iz + nzCum
          do ix = 1, outBC%bcArray(imgrid)%nx
            do ic = 1, ccoeff_current
            ! Passing the YMin face
              outBC%bcArray(imgrid)%xYMin(ix,iz)  = (outBC%bcArray(imgrid)%xYMin(ix,iz) + &
                              inBC%xYMin(ccoeff_current*(ix-1)+ic, izv))*inBC%grid%dx(ccoeff_current*(ix-1)+ic)/&
                              inBC%grid%gridArray(imgrid)%dx(ix)
            ! Passing the YMax face
              outBC%bcArray(imgrid)%xYMax(ix,iz) = (outBC%bcArray(imgrid)%xYMax(ix,iz) + &
                              inBC%xYMax(ccoeff_current*(ix-1)+ic, izv))*inBC%grid%dx(ccoeff_current*(ix-1)+ic)/&
                              inBC%grid%gridArray(imgrid)%dx(ix)
            enddo  ! ic
          enddo  ! ix
        enddo    ! iz

        ! z-component
        do iz = 1, outBC%bcArray(imgrid)%nz
           izv = iz + nzCum
          do ix = 1, outBC%bcArray(imgrid)%nx+1
            ! Passing the YMin face
            outBC%bcArray(imgrid)%zYmin(ix,iz) = inBC%zYMin(ccoeff_current*(ix-1)+1,izv)
            ! Passing the YMax face
            outBC%bcArray(imgrid)%zYMax(ix,iz) = inBC%zYMax(ccoeff_current*(ix-1)+1,izv)
           enddo  ! ix
        enddo     ! iz

       ! passing boundary values for the XMax and XMin face
       ! y-component
       do iz = 1, outBC%bcArray(imgrid)%nz+1
          izv = iz + nzCum
         do iy = 1, outBC%bcArray(imgrid)%ny
           do ic = 1, ccoeff_current
           ! Passing the XMin face
              outBC%bcArray(imgrid)%yXMin(iy, iz) = (outBC%bcArray(imgrid)%yXMin(iy, iz) + &
                                  inBC%yXMin(ccoeff_current*(iy-1)+ic, izv))*inBC%grid%dy(ccoeff_current*(iy-1)+ic)/&
                                  inBC%grid%gridArray(imgrid)%dy(iy)
           ! Passing the XMax face
              outBC%bcArray(imgrid)%yXMax(iy, iz) = (outBC%bcArray(imgrid)%yXMax(iy, iz) + &
                                  inBC%yXMax(ccoeff_current*(iy-1)+ic, izv))*inBC%grid%dy(ccoeff_current*(iy-1)+ic)/&
                                  inBC%grid%gridArray(imgrid)%dy(iy)
           enddo ! ic
         enddo  ! iy
        enddo   ! iz

       ! z-component
       do iz = 1, outBC%bcArray(imgrid)%nz
          izv = iz + nzCum
         do iy = 1, outBC%bcArray(imgrid)%ny+1
           ! Passing the XMin face
           outBC%bcArray(imgrid)%zXMin(iy, iz) = inBC%zXMin(ccoeff_current*(iy-1)+1, izv)
           ! Passing the XMax face
           outBC%bcArray(imgrid)%zXMax(iy, iz) = inBC%zXMax(ccoeff_current*(iy-1)+1, izv)
         enddo  ! iy
       enddo  ! iz
      nzCum = outBC%bcArray(imgrid)%nz + nzCum
     enddo ! Global loop over sub-grids

     ! passing boundary values for the ZMin and ZMax face
     ! x-component
     do iy = 1, outBC%bcArray(1)%ny+1
       do ix = 1, outBC%bcArray(1)%nx
         do ic =1, 2**outBC%coarseness(1)
            ! Passing the ZMin face
            outBC%bcArray(1)%xZMin(ix, iy) = (outBC%bcArray(1)%xZMin(ix, iy) + &
                                             inBC%xZMin(2**outBC%coarseness(1)*(ix-1)+ic, 2**outBC%coarseness(1)*(iy-1)+1))*&
                                             inBC%grid%dx(2**outBC%coarseness(1)*(ix-1)+ic)/&
                                             inBC%grid%gridArray(1)%dx(ix)
         enddo ! ic
       enddo
     enddo
     do iy = 1, outBC%bcArray(outBC%mgridSize)%ny+1
       do ix = 1, outBC%bcArray(outBC%mgridSize)%nx
            ! Passing the ZMax face
         do ic =1, 2**outBC%coarseness(outBC%mgridSize)
              outBC%bcArray(outBC%mgridSize)%xZMax(ix, iy) = (outBC%bcArray(outBC%mgridSize)%xZMax(ix, iy)&
                                                   + inBC%xZMax(2**outBC%coarseness(outBC%mgridSize)*(ix-1)+ic, 2**outBC%coarseness(outBC%mgridSize)*(iy-1)+1))*&
                                                   inBC%grid%dx(2**outBC%coarseness(outBC%mgridSize)*(ix-1)+ic)/&
                                                   inBC%grid%gridArray(outBC%mgridSize)%dx(ix)
         enddo   ! ic
       enddo  ! ix
     enddo  !iy

     ! y-component
     do ix = 1, outBC%bcArray(1)%nx+1
        do iy = 1, outBC%bcArray(1)%ny
          do ic = 1, 2**outBC%coarseness(1)
            outBC%bcArray(1)%yZMin(ix, iy) =  (outBC%bcArray(1)%yZMin(ix, iy) + inBC%yZMin(2**outBC%coarseness(1)*(ix-1)+1, 2**outBC%coarseness(1)*(iy-1)+ic))*&
                                                      inBC%grid%dy(2**outBC%coarseness(1)*(iy-1)+ic)/&
                                                      inBC%grid%gridArray(1)%dy(iy)
          enddo ! ic
        enddo
     enddo
     do ix = 1, outBC%bcArray(outBC%mgridSize)%nx+1
        do iy = 1, outBC%bcArray(outBC%mgridSize)%ny
          do ic = 1, 2**outBC%coarseness(outBC%mgridSize)
            outBC%bcArray(outBC%mgridSize)%yZMax(ix, iy) = (outBC%bcArray(outBC%mgridSize)%yZMax(ix, iy)+ inBC%yZMax(2**outBC%coarseness(outBC%mgridSize)*(ix-1)+1, 2**outBC%coarseness(outBC%mgridSize)*(iy-1)+ic))*&
                                                     inBC%grid%dy(2**outBC%coarseness(outBC%mgridSize)*(iy-1)+ic)/&
                                                     inBC%grid%gridArray(outBC%mgridSize)%dy(iy)
          enddo ! ic
        enddo  ! iy
      enddo  ! ix

    end subroutine b2bmg

end module sg_boundary ! sg_boundary
