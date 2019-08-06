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
  !
  ! Dimensions:
  !    E%xYMax(nx,nz+1)
  !    E%zYMax(nx+1,nz)
  !    E%xYMin(nx,nz+1)
  !    E%zYMin(nx+1,nz)
  !    E%yXMax(ny,nz+1)
  !    E%zXMax(ny+1,nz)
  !    E%yXMin(ny,nz+1)
  !    E%zXMin(ny+1,nz)
  !    E%xZMin(nx,ny+1)
  !    E%yZMin(nx+1,ny)
  !    E%xZMax(nx,ny+1)
  !    E%yZMax(nx+1,ny)

  use math_constants
  use griddef
  use sg_vector
  use sg_sparse_vector
  implicit none

  ! Generic interfaces are done through subroutines
  ! creates boundary nodes
  INTERFACE create
     module procedure create_cboundary
  END INTERFACE

  ! deallocates the boundary nodes
  INTERFACE deall
     module procedure deall_cboundary
  END INTERFACE

  ! scalar multiplies the boundary nodes
  INTERFACE scMult
     module procedure scMult_cboundary
  END INTERFACE

  INTERFACE linComb
     MODULE PROCEDURE linComb_cboundary
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
     module procedure copy_sparseveccb
  END INTERFACE

  ! sets boundary nodes in the vector field
  INTERFACE setBC
     module procedure copy_bcvector
     module procedure copy_bsparsevecc
  END INTERFACE

  public                :: copy_cbvector, copy_bcvector
  public                :: copy_sparseveccb, copy_bsparsevecc
  public				:: create_cboundary, &
       deall_cboundary, copy_cboundary, zero_cboundary, &
       scMult_cboundary, scMult_cboundary_f, add_cboundary, &
       add_cboundary_f, diagMult_cboundary, diagMult_cboundary_f, &
       dotProd_cboundary, &
       ! the two routines below are not part of any interface
       linComb_cboundary, scMultAdd_cboundary, random_cboundary, &
       ! and these are for reading and writing
       write_cboundary, read_cboundary


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

     ! pointer to parent grid
     type (grid_t), pointer                             :: grid

  end type cboundary


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

    E%nx = 0
    E%ny = 0
    E%nz = 0
    E%allocated = .false.

  end subroutine deall_cboundary  ! deall_cboundary


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
          Call create_cboundary(E1%grid, E2)
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

  ! **********************************************************************
  ! * Creates a random perturbation in cboundary;  AK 2018-04-25

  subroutine random_cboundary(E,eps)

    implicit none
    type (cboundary), intent(inout)                  :: E
    real(8), intent(in), optional                    :: eps
    ! local
    real (kind(E%xYMax)), allocatable, dimension(:,:) :: xr,xi
    integer              :: nx,ny,nz,nmax,istat

    if (.not. E%allocated) then
      call warning('cboundary not allocated in random_cboundary')
      return
    end if

    call zero_cboundary(E)

    ! make some random vectors
    nx = E%nx
    ny = E%ny
    nz = E%nz
    nmax = max(max(nx+1,ny+1),nz+1)

    allocate(xr(nmax,nmax),xi(nmax,nmax),STAT=istat)

    call random_number(xr)
    call random_number(xi)
    E%xYMax = cmplx(xr(1:nx,1:nz+1),xi(1:nx,1:nz+1))

    call random_number(xr)
    call random_number(xi)
    E%zYMax = cmplx(xr(1:nx+1,1:nz),xi(1:nx+1,1:nz))

    call random_number(xr)
    call random_number(xi)
    E%xYMin = cmplx(xr(1:nx,1:nz+1),xi(1:nx,1:nz+1))

    call random_number(xr)
    call random_number(xi)
    E%zYMin = cmplx(xr(1:nx+1,1:nz),xi(1:nx+1,1:nz))

    call random_number(xr)
    call random_number(xi)
    E%yXMax = cmplx(xr(1:ny,1:nz+1),xi(1:ny,1:nz+1))

    call random_number(xr)
    call random_number(xi)
    E%zXMax = cmplx(xr(1:ny+1,1:nz),xi(1:ny+1,1:nz))

    call random_number(xr)
    call random_number(xi)
    E%yXMin = cmplx(xr(1:ny,1:nz+1),xi(1:ny,1:nz+1))

    call random_number(xr)
    call random_number(xi)
    E%zXMin = cmplx(xr(1:ny+1,1:nz),xi(1:ny+1,1:nz))

    call random_number(xr)
    call random_number(xi)
    E%xZMin = cmplx(xr(1:nx,1:ny+1),xi(1:nx,1:ny+1))

    call random_number(xr)
    call random_number(xi)
    E%yZMin = cmplx(xr(1:nx+1,1:ny),xi(1:nx+1,1:ny))

    call random_number(xr)
    call random_number(xi)
    E%xZMax = cmplx(xr(1:nx,1:ny+1),xi(1:nx,1:ny+1))

    call random_number(xr)
    call random_number(xi)
    E%yZMax = cmplx(xr(1:nx+1,1:ny),xi(1:nx+1,1:ny))

    deallocate(xr,xi,STAT=istat)

    if (present(eps)) then
       E%xYMax = E%xYMax * eps
       E%zYMax = E%zYMax * eps
       E%xYMin = E%xYMin * eps
       E%zYMin = E%zYMin * eps
       E%yXMax = E%yXMax * eps
       E%zXMax = E%zXMax * eps
       E%yXMin = E%yXMin * eps
       E%zXMin = E%zXMin * eps
       E%xZMin = E%xZMin * eps
       E%yZMin = E%yZMin * eps
       E%xZMax = E%xZMax * eps
       E%yZMax = E%yZMax * eps
    else
       E%xYMax = E%xYMax * 0.05
       E%zYMax = E%zYMax * 0.05
       E%xYMin = E%xYMin * 0.05
       E%zYMin = E%zYMin * 0.05
       E%yXMax = E%yXMax * 0.05
       E%zXMax = E%zXMax * 0.05
       E%yXMin = E%yXMin * 0.05
       E%zXMin = E%zXMin * 0.05
       E%xZMin = E%xZMin * 0.05
       E%yZMin = E%yZMin * 0.05
       E%xZMax = E%xZMax * 0.05
       E%yZMax = E%yZMax * 0.05
    end if


  end subroutine random_cboundary

  !******************************************************************************
  ! write_cboundary writes a BC vector  in a simple ASCII format; vector has
  ! to exist and be allocated before calling this routine, and the file unit
  ! must already be available for writing; AK 2018-04-25
  subroutine write_cboundary(fid, E)

      integer,          intent(in)      :: fid
      type (cboundary), intent(in)      :: E

      !  local variables
      type (sparsevecc)                 :: SV
      integer                           :: ii, istat

      if(.not. E%allocated) then
         write(0, *) 'BC vector must be allocated before call to write_cboundary'
         return
      endif

      ! write as sparse vector in ascii
      call copy_bsparsevecc(E,SV)
      call write_sparsevecc(fid,SV)
      call deall_sparsevecc(SV)

!      do ii = 1,SV%nCoeff
!         write(fid,'(4i5,2es14.6)',iostat=istat) SV%i(ii),SV%j(ii),SV%k(ii),SV%xyz(ii),real(SV%c(ii)),aimag(SV%c(ii))
!      end do
!
!      write(fid,'(a11)',iostat=istat) 'xYMin xYMax'
!      do k = 1,E%nz+1
!          do i = 1,E%nx
!              write(fid,'(3i5,6es14.6)',iostat=istat) i,k,E%xYMin(i,k),E%xYMax(i,k)
!          end do
!      end do
!
!      write(fid,'(a11)',iostat=istat) 'zYMin zYMax'
!      do k = 1,E%nz
!          do i = 1,E%nx+1
!              write(fid,'(3i5,6es14.6)',iostat=istat) i,k,E%zYMin(i,k),E%zYMax(i,k)
!          end do
!      end do
!
!      write(fid,'(a11)',iostat=istat) 'yXMin yXMax'
!      do k = 1,E%nz+1
!          do j = 1,E%ny
!              write(fid,'(3i5,6es14.6)',iostat=istat) j,k,E%yXMin(j,k),E%yXMax(j,k)
!          end do
!      end do


  end subroutine write_cboundary

  !**********************************************************************************
  ! read_cboundary reads a BC vector in a simple ASCII format; vector must match
  ! the input grid; file unit must already be available for reading; AK 2018-04-25
  subroutine read_cboundary(fid, E, grid)

      integer,        intent(in)          :: fid
      type (cboundary), intent(inout)     :: E
      type (grid_t), intent(in), optional :: grid

      !  local variables
      type (sparsevecc)                 :: SV
      integer                           :: ii, istat

      if (.not. E%allocated) then
        if (present(grid)) then
            call create_cboundary(grid,E)
        else
            write(0, *) 'BC vector not allocated: unable to read cboundary from file'
            stop
        end if
      end if

      ! write as sparse vector in ascii
      call read_sparsevecc(fid,SV)

      ! if grid is available, make sure the two are consistent
      if (present(grid)) then
        if ((maxval(SV%i) > grid%nx+1) .or. (maxval(SV%j) > grid%ny+1) .or. (maxval(SV%k) > grid%nz+1)) then
            write(0, *) 'BC vector size does not match the grid in read_cboundary'
            write(0, *) 'NX: ',grid%nx,maxval(SV%i)
            write(0, *) 'NY: ',grid%ny,maxval(SV%j)
            write(0, *) 'NZ: ',grid%nz,maxval(SV%k)
            return
        endif
      endif

      call copy_sparseveccb(SV,E)
      call deall_sparsevecc(SV)

  end subroutine read_cboundary

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
    Call create_cboundary(E1%grid, E3)
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
          do ix = 1, outE%nx+1
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
  ! * copy_sparsevectorb converts a compatible sparse vector directly into BC
  ! * avoiding the sometimes unnecessary conversion to cvector. BC needs to be
  ! * pre-allocated on input to this routine, and sparse vector needs to
  ! * store the BC for this operation to make sense. AK 2018-04-25
  subroutine copy_sparseveccb(inSV, outBC)

    implicit none
    type (sparsevecc), intent(in)                  :: inSV
    type (cboundary), intent(inout)                :: outBC
    integer                      :: ix, iy, iz, ic, ib, ibprev  ! dummy integers
    integer                      :: nx, ny, nz, nCoeff, duplCoeff   ! dimensions
    integer                      :: YMax, YMin       ! ends for BC writing
    integer                      :: XMax, XMin       ! ends for BC writing
    integer                      :: ZMin, ZMax       ! ends for BC writing


    if (.not.inSV%allocated) then
      write(0,*) 'inSV in copy_sparseveccb not allocated yet'
      stop
    end if

    if (.not.outBC%allocated) then
      write(0,*) 'outBC in copy_sparseveccb not allocated yet'
      stop
    else
        nx = outBC%nx
        ny = outBC%ny
        nz = outBC%nz
        nCoeff = 2*(nx*(nz+1)) + 2*((nx+1)*nz) + 2*(ny*(nz+1)) &
               + 2*((ny+1)*nz) + 2*(nx*(ny+1)) + 2*((nx+1)*ny)
        duplCoeff = 4*nx + 4*ny + 4*nz ! duplicate cube edges
        if (inSV%nCoeff < nCoeff - duplCoeff) then
            write(0,*) 'Error-complex sparse vector does not have complete BC info in copy_sparseveccb'
            write(0,*) 'Values in BC sparse vector: ',inSV%nCoeff
            write(0,*) 'Expected number of values : ',nCoeff,' less ',duplCoeff,' duplicate edges'
            stop
        elseif (inSV%gridType .ne. EDGE) then
            write(0,*) 'Unable to copy sparsevecc of type ',trim(inSV%gridType),' to BC vector'
            write(0,*) 'Input BC sparse vector is not defined on grid edges'
            stop
        end if
    end if

    YMax = outBC%ny+1
    YMin = 1
    XMax = outBC%nx+1
    XMin = 1
    ZMax = outBC%nz+1
    ZMin = 1

    ib = 0
    ibprev = 0

    do ic = 1,inSV%nCoeff
        ix = inSV%i(ic)
        iy = inSV%j(ic)
        iz = inSV%k(ic)
        ibprev = ib
        ! passing boundary values for the YMin and YMax sides
        if ((iy == YMin) .and. (inSV%xyz(ic) == 1)) then ! x-component YMin face
            outBC%xYMin(ix, iz) = inSV%c(ic)
            ib = ib+1
        elseif ((iy == YMax) .and. (inSV%xyz(ic) == 1)) then ! x-component YMax face
            outBC%xYMax(ix, iz) = inSV%c(ic)
            ib = ib+1
        elseif ((iy == YMin) .and. (inSV%xyz(ic) == 3)) then ! z-component YMin face
            outBC%zYMin(ix, iz) = inSV%c(ic)
            ib = ib+1
        elseif ((iy == YMax) .and. (inSV%xyz(ic) == 3)) then ! z-component YMax face
            outBC%zYMax(ix, iz) = inSV%c(ic)
            ib = ib+1
        end if
        ! passing boundary values for the XMin and XMax sides
        if ((ix == XMin) .and. (inSV%xyz(ic) == 2)) then ! y-component XMin face
            outBC%yXMin(iy, iz) = inSV%c(ic)
            ib = ib+1
        elseif ((ix == XMax) .and. (inSV%xyz(ic) == 2)) then ! y-component XMax face
            outBC%yXMax(iy, iz) = inSV%c(ic)
            ib = ib+1
        elseif ((ix == XMin) .and. (inSV%xyz(ic) == 3)) then ! z-component XMin face
            outBC%zXMin(iy, iz) = inSV%c(ic)
            ib = ib+1
        elseif ((ix == XMax) .and. (inSV%xyz(ic) == 3)) then ! z-component XMax face
            outBC%zXMax(iy, iz) = inSV%c(ic)
            ib = ib+1
        end if
        ! passing boundary values for the ZMin and ZMax sides
        if ((iz == ZMin) .and. (inSV%xyz(ic) == 1)) then ! x-component ZMin face
            outBC%xZMin(ix, iy) = inSV%c(ic)
            ib = ib+1
        elseif ((iz == ZMax) .and. (inSV%xyz(ic) == 1)) then ! x-component ZMax face
            outBC%xZMax(ix, iy) = inSV%c(ic)
            ib = ib+1
        elseif ((iz == ZMin) .and. (inSV%xyz(ic) == 2)) then ! y-component ZMin face
            outBC%yZMin(ix, iy) = inSV%c(ic)
            ib = ib+1
        elseif ((iz == ZMax) .and. (inSV%xyz(ic) == 2)) then ! y-component ZMax face
            outBC%yZMax(ix, iy) = inSV%c(ic)
            ib = ib+1
        end if
    end do

    if (ib .ne. nCoeff) then
        write(0,*) 'Warning: BC copy from complex sparse vector might have failed in copy_sparseveccb'
        write(0,*) 'copied ',ib,' values to fill in ',nCoeff,' values on BC edges...'
    end if

end subroutine copy_sparseveccb

  ! ***************************************************************************
  ! * copy_bsparsevector reads the boundary conditions as an input and writes
  ! * them to a complex sparse vector; AK 2018-04-25
  subroutine copy_bsparsevecc(inBC, outSV)

    implicit none
    type (cboundary), intent(in)     :: inBC
    ! boundary conditions as an input
    type (sparsevecc), intent(inout)    :: outSV
    ! the electrical field as an output
    integer                      :: ix, iy, iz, ic       ! dummy integers
    integer                      :: nx, ny, nz, nCoeff   ! dimensions
    integer                      :: YMax, YMin       ! ends for BC writing
    integer                      :: XMax, XMin       ! ends for BC writing
    integer                      :: ZMin, ZMax       ! ends for BC writing

    if (.not.inBC%allocated) then
      write(0,*) 'inBC in copy_bsparsevecc not allocated yet'
      stop
    end if

    nx = inBC%nx
    ny = inBC%ny
    nz = inBC%nz

    nCoeff = 2*(nx*(nz+1)) + 2*((nx+1)*nz) + 2*(ny*(nz+1)) &
           + 2*((ny+1)*nz) + 2*(nx*(ny+1)) + 2*((nx+1)*ny) &
           - (4*nx + 4*ny + 4*nz) ! no need to store duplicate edges

    call create_sparsevecc(nCoeff,outSV,EDGE)

    outSV%gridType = EDGE
    ic = 1

    YMax = inBC%ny+1
    YMin = 1
    XMax = inBC%nx+1
    XMin = 1
    ZMax = inBC%nz+1
    ZMin = 1

    ! NOTE: mind the 12 cube edges. Unless something special is done there,
    ! x-components on x-edges at iy=1 iz=1, iy=YMax iz=1, iy=1 iz=ZMax, iy=YMax iz=ZMax get written twice
    ! y-components on y-edges at ix=1 iz=1, ix=XMax iz=1, ix=1 iz=ZMax, ix=XMax iz=ZMax get written twice
    ! z-components on z-edges at iy=1 ix=1, iy=YMax ix=1, iy=1 ix=XMax, iy=YMax ix=XMax get written twice

    ! passing boundary condition values to the complete description
    ! of the vector field. To avoid duplication at cube edges, wrap around
    ! the cube and ignore the last set of edges for each side

    ! Going around the cube, clockwise:
    ! x-component - ZMin,YMax,ZMax,YMin
    ! passing boundary values for the ZMin face
    ! x-component at all y-edges except last
    do ix = 1, nx
        do iy = 1, ny
            ! Passing the ZMin face
            outSV%i(ic) = ix
            outSV%j(ic) = iy
            outSV%k(ic) = ZMin
            outSV%xyz(ic) = 1
            outSV%c(ic) = inBC%xZMin(ix, iy)
            ic = ic+1
        enddo   ! iy
    enddo      ! ix
    ! passing boundary values for the YMax face
    ! x-component at all y-edges except last
    do ix = 1, nx
        do iz = 1, nz
            ! Passing the YMax face
            outSV%i(ic) = ix
            outSV%j(ic) = YMax
            outSV%k(ic) = iz
            outSV%xyz(ic) = 1
            outSV%c(ic) = inBC%xYMax(ix, iz)
            ic = ic+1
        enddo ! iz
    enddo    ! ix
    ! passing boundary values for the ZMax face
    ! x-component at all y-edges except first
    do ix = 1, nx
        do iy = 2, ny+1
            ! Passing the ZMax face
            outSV%i(ic) = ix
            outSV%j(ic) = iy
            outSV%k(ic) = ZMax
            outSV%xyz(ic) = 1
            outSV%c(ic) = inBC%xZMax(ix, iy)
            ic = ic+1
        enddo   ! iy
    enddo      ! ix
    ! passing boundary values for the YMin face
    ! x-component at all z-edges except first
    do ix = 1, nx
        do iz = 2, nz+1
            ! Passing the YMin face
            outSV%i(ic) = ix
            outSV%j(ic) = YMin
            outSV%k(ic) = iz
            outSV%xyz(ic) = 1
            outSV%c(ic) = inBC%xYMin(ix, iz)
            ic = ic+1
        enddo ! iz
    enddo    ! ix

    ! Going around the cube, clockwise:
    ! y-component - ZMin,XMax,ZMax,XMin
    ! passing boundary values for the ZMin face
    ! y-component at all x-edges except last
    do ix = 1, nx
        do iy = 1, ny
            ! Passing the ZMin face
            outSV%i(ic) = ix
            outSV%j(ic) = iy
            outSV%k(ic) = ZMin
            outSV%xyz(ic) = 2
            outSV%c(ic) = inBC%yZMin(ix, iy)
            ic = ic+1
        enddo   ! ix
    enddo      ! iy
    ! passing boundary values for the XMax face
    ! y-component at all z-edges except last
    do iy = 1, ny
        do iz = 1, nz
            ! Passing the XMax face
            outSV%i(ic) = XMax
            outSV%j(ic) = iy
            outSV%k(ic) = iz
            outSV%xyz(ic) = 2
            outSV%c(ic) = inBC%yXMax(iy, iz)
            ic = ic+1
        enddo  ! iy
    enddo     ! iz
    ! passing boundary values for the ZMax face
    ! y-component at all x-edges except first
    do ix = 2, nx+1
        do iy = 1, ny
            ! Passing the ZMax face
            outSV%i(ic) = ix
            outSV%j(ic) = iy
            outSV%k(ic) = ZMax
            outSV%xyz(ic) = 2
            outSV%c(ic) = inBC%yZMax(ix, iy)
            ic = ic+1
        enddo   ! ix
    enddo      ! iy
    ! passing boundary values for the XMin face
    ! y-component at all z-edges except first
    do iy = 1, ny
        do iz = 2, nz+1
            ! Passing the XMin face
            outSV%i(ic) = XMin
            outSV%j(ic) = iy
            outSV%k(ic) = iz
            outSV%xyz(ic) = 2
            outSV%c(ic) = inBC%yXMin(iy, iz)
            ic = ic+1
        enddo  ! iy
    enddo     ! iz

    ! Going around the cube, clockwise:
    ! z-component = XMin,YMax,XMax,YMin
    ! passing boundary values for the XMin face
    ! z-component at all y-edges except last
    do iy = 1, ny
        do iz = 1, nz
            ! Passing the XMin face
            outSV%i(ic) = XMin
            outSV%j(ic) = iy
            outSV%k(ic) = iz
            outSV%xyz(ic) = 3
            outSV%c(ic) = inBC%zXMin(iy, iz)
            ic = ic+1
        enddo  ! iy
    enddo     ! iz
    ! passing boundary values for the YMax face
    ! z-component at all x-edges except last
    do ix = 1, nx
        do iz = 1, nz
            ! Passing the YMax face
            outSV%i(ic) = ix
            outSV%j(ic) = YMax
            outSV%k(ic) = iz
            outSV%xyz(ic) = 3
            outSV%c(ic) = inBC%zYMax(ix, iz)
            ic = ic+1
        enddo  ! iz
    enddo     ! ix
    ! passing boundary values for the XMax face
    ! z-component at all y-edges except first
    do iy = 2, ny+1
        do iz = 1, nz
            ! Passing the XMax face
            outSV%i(ic) = XMax
            outSV%j(ic) = iy
            outSV%k(ic) = iz
            outSV%xyz(ic) = 3
            outSV%c(ic) = inBC%zXMax(iy, iz)
            ic = ic+1
        enddo  ! iy
    enddo     ! iz
    ! passing boundary values for the YMin face
    ! z-component at all x-edges except first
    do ix = 2, nx+1
        do iz = 1, nz
            ! Passing the YMin face
            outSV%i(ic) = ix
            outSV%j(ic) = YMin
            outSV%k(ic) = iz
            outSV%xyz(ic) = 3
            outSV%c(ic) = inBC%zYMin(ix, iz)
            ic = ic+1
        enddo  ! iz
    enddo     ! ix

    outSV%allocated = .true.
    outSV%temporary = .false.

end subroutine copy_bsparsevecc

end module sg_boundary ! sg_boundary
