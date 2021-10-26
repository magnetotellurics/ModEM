!**
! SUMMARY
!
! Standard cartesian grid vectors.
!*
module cVector3D_SG
  use Constants
  use MatUtils
  use Grid
  use Grid3D_SG
  use rScalar
  use rScalar3D_SG
  use rVector
  use rVector3D_SG
  use cVector
  
  implicit none
  
  private
  
  public :: cVector3D_SG_t
  
  !**
  !*
  type, extends(cVector_t) :: cVector3D_SG_t
     !**
     ! pointer to parent grid
     class(Grid3D_SG_t), pointer :: grid => null()
     
     !**
     ! Store the intention of the use in a character
     ! string defined in GridDef as a parameter: EDGE
     ! or FACE are two possibilities.
     !*
     character(len = 80) :: gridType = ""

     !**
     ! Grid Dimensions:
     ! nx is grid dimension (number of cells) in the x-direction
     ! ny is grid dimension (number of cells) in the y-direction
     ! nz is grid dimension (number of cells) in the z-direction:
     !*
     integer :: nx = 0, ny = 0, nz = 0

     integer, dimension(3) :: NdX = 0, NdY = 0, NdZ = 0
     integer, dimension(3) :: Nxyz = 0

     !**
     ! Typical usage:  electrical fields on cell edges of
     ! staggered grid
     ! For example, in an EDGE, the dimensions would be
     ! x: edge nodes in x-direction: dimension Nx, Ny+1, Nz+1
     ! y: edge nodes in y-direction: dimension Nx+1, Ny, Nz+1
     ! z: edge nodes in z-direction: dimension Nx+1, Ny+1, Nz
     ! Note that the arrays are defined through dynamic
     ! memory allocation
     !*
     complex(kind = prec), allocatable, dimension(:, :, :) :: x, y, z

   contains
     
     !**
     ! Initialization and finalization
     !*
     final :: cVector3D_SG_dtor
     
     !**
     ! Input/Output
     !*
     procedure, public :: Read
     procedure, public :: Write
     
     !**
     ! Boundary operations
     !*
     procedure, public :: SetAllBoundary
     procedure, public :: SetOneBoundary
     procedure, public :: SetAllInterior
     procedure, public :: IntBdryIndices
     procedure, public :: Boundary
     procedure, public :: Interior

     !**
     ! Data access
     !*
     procedure, public :: Length
     procedure, public :: GetArray
     procedure, public :: SetArray
     procedure, public :: SetVecComponents
          
     !**
     ! Arithmetic operations
     !*

     procedure, public :: Zeros   
     procedure, public :: Add_1   
     procedure, public :: Sub_1   
     procedure, public :: Mult_1
     procedure, public, pass(self) :: Mult_2
     procedure, public :: Mult_3
     procedure, public :: Mult_s_1
     procedure, public :: Div_1
     procedure, public :: Div_2     
     procedure, public :: dotProd

     !**
     ! Miscellaneous
     !*
     procedure, public :: isCompatible_1
     procedure, public :: isCompatible_2
     generic :: isCompatible => isCompatible_1, isCompatible_2
     
     procedure, public :: CopyFrom         
     procedure, public :: interpFunc
     
  end type cVector3D_SG_t
  
  interface cVector3D_SG_t
     module procedure cVector3D_SG_ctor
  end interface cVector3D_SG_t
  
contains
  
  !
  !************************************************
  ! Initialization and finalization
  !************************************************
  !
  
  !**
  ! Constructor for REAL vectors.
  !
  ! Arguments
  !   igrid      Underlying grid.
  !   gridtype  Definied in GridDef.f90
  !
  !*
  function cVector3D_SG_ctor(igrid, gridType) result (E)
    implicit none
    ! Arguments
    class(Grid3D_SG_t), target, intent(in) :: igrid
    character(*)              , intent(in) :: gridType
    ! Local variables
    integer :: nx, ny, nz, nz_earth, nzAir
    integer :: status
    type(cVector3D_SG_t) :: E
    
    E%grid => igrid
    
    ! Grid dimensions
    nx = igrid%nx
    ny = igrid%ny
    nz = igrid%nz
    
    E%nx = nx
    E%ny = ny
    E%nz = nz
    
    E%gridType = trim(gridType)
    E%isAllocated = .false.
    
    if (E%gridType == EDGE) then
       allocate(E%x(nx, ny + 1, nz + 1), STAT = status)
       E%isAllocated = status.eq.0
       
       allocate(E%y(nx + 1, ny, nz + 1), STAT = status)
       E%isAllocated = E%isAllocated.and.(status.eq.0)
       
       allocate(E%z(nx + 1, ny + 1, nz), STAT = status)
       E%isAllocated = E%isAllocated.and.(status.eq.0)

       E%NdX = (/E%nx, E%ny + 1, E%nz + 1/)
       E%NdY = (/E%nx + 1, E%ny, E%nz + 1/)
       E%NdZ = (/E%nx + 1, E%ny + 1, E%nz/)         
       
    else if (E%gridType == FACE) then
       allocate(E%x(nx + 1, ny, nz), STAT = status)
       E%isAllocated = status.eq.0
       
       allocate(E%y(nx, ny + 1, nz), STAT = status)
       E%isAllocated = E%isAllocated.and.(status.eq.0)
       
       allocate(E%z(nx, ny, nz + 1), STAT = status)
       E%isAllocated = E%isAllocated.and.(status.eq.0)
       
       E%NdX = (/E%nx + 1, E%ny, E%nz/)
       E%NdY = (/E%nx, E%ny + 1, E%nz/)
       E%NdZ = (/E%nx, E%ny, E%nz + 1/)
       
    else
       write(*, *) 'ERROR:cVector3D_SG_ctor:'
       write(*, *) '      Only EDGE or FACE types allowed. Exiting.'
       
       STOP       
    end if
    
    if (E%isAllocated) then
       E%x = R_ZERO
       E%y = R_ZERO
       E%z = R_ZERO      
    else
       write(*, *) 'ERROR:cVector3D_SG_ctor:'
       write(*, *) '      Unable to allocate vector. Exiting.'
       
       STOP
    end if
    
    E%Nxyz = (/product(E%NdX), &
         product(E%NdY), &
         product(E%NdZ)/)
    
  end function cVector3D_SG_ctor
  
  subroutine cVector3D_SG_dtor(E)
    implicit none
    ! Arguments
    type(cVector3D_SG_t), intent(inout) :: E
    ! Local variables
    integer :: status
    
    if (E%isAllocated) then 
       deallocate(E%x, STAT = status)
       deallocate(E%y, STAT = status)
       deallocate(E%z, STAT = status)
    end if
    
    nullify(E%grid)
    
    E%nx = 0
    E%ny = 0
    E%nz = 0
    
    E%gridType = ''
    E%isAllocated = .false.        
    
  end subroutine cVector3D_SG_dtor
  
  !
  !************************************************
  ! Input/Output
  !************************************************
  !

  !**
  ! Read
  !*
  subroutine Read(self, fid)
    implicit none
    ! Arguments
    class(cVector3D_SG_t), intent(inout) :: self
    integer              , intent(in)    :: fid
    ! Local variables
    integer :: Nx, Ny, Nz
    character(80) :: gridType
    integer :: i, j, k, ii, jj, kk, istat
    real(kind = prec), allocatable, dimension(:, :, :) :: x, y, z    
    logical :: ok, hasname, binary = .true.
    character(80) :: fname, isbinary
    
    inquire(fid, opened = ok, named = hasname, name = fname, unformatted = isbinary)
    
    ! Check that the file is unformatted if binary, formatted if ascii.
    if ((index(isbinary, 'yes') > 0 .or.index(isbinary, 'YES') > 0) &
         .and. .not.binary ) then
       write(0, *) 'ERROR:cVector3D_SG_t::Read: '
       write(0, *) '      Unable to read vector from unformatted file. ', &
            trim(fname), '. Exiting.'
       
       STOP
    else if ((index(isbinary, 'no') > 0 .or.index(isbinary, 'NO') > 0) &
         .and.binary) then
       write(0, *) 'ERROR:cVector3D_SG_t::Read: '
       write(0, *) '      Unable to read vector from formatted file ', &
            trim(fname), '. Exiting.'
       
       STOP
    end if
    
    read(fid) Nx, Ny, Nz, gridType
    
    if (.not.self%isAllocated) then
       write(*, *) 'ERROR:cVector3D_SG_t::Read: '
       write(*, *) '      Vector must be allocated before reading from ', &
            trim(fname), '. Exiting.'
       
       STOP
    else if (self%gridType.ne.gridType) then
       write (*, *) 'ERROR:cVector3D_SG_t::Read: '
       write (*, *) '      Vector must be of type ', gridType, &
            &       '      before reading from ', trim (fname), '. Exiting.'

       STOP
    else if ((self%nx.ne.Nx).or. &
         (self%ny.ne.Ny).or.(self%nz.ne.Nz)) then
       write(*, *) 'ERROR:cVector3D_SG_t::Read: '
       write(*, *) '      Wrong size of vector on input from ', &
            trim (fname), '. Exiting.'
       
       STOP
    end if
    
    read(fid) self%x
    read(fid) self%y
    read(fid) self%z          
    
  end subroutine Read
  
  !**
  ! Write
  !*
  subroutine Write(self, fid)
    implicit none
    ! Arguments
    class(cVector3D_SG_t), intent(in) :: self
    integer              , intent(in) :: fid
    ! Local variables
    integer :: Nx, Ny, Nz
    integer :: i, j, k, istat
    real(kind = prec), allocatable, dimension(:, :, :) :: x, y, z    
    logical :: ok, hasname, binary = .true.
    character(80) :: fname, isbinary
    
    if (.not.self%isAllocated) then
       write(*, *) 'ERROR:cVector3D_SG::Write: '
       write(*, *) '      Vector not allocated. Exiting.'

       STOP
    end if
    
    inquire(fid, opened = ok, named = hasname, name = fname, unformatted = isbinary)
    
    ! Check that the file is unformatted if binary, formatted if ascii.
    if ((index(isbinary, 'yes') > 0.or.index(isbinary, 'YES') > 0) &
         .and. .not.binary) then
       write (0, *) 'ERROR:cVector3D_SG::Write: '
       write (0, *) '      Unable to write vector to unformatted file. ', &
            trim(fname), '. Exiting.'
       
       STOP       
    else if ((index(isbinary,'no') > 0.or.index(isbinary,'NO') > 0) &
         .and.binary) then
       write (0, *) 'ERROR:cVector3D_SG::Write: '
       write (0, *) '      Unable to write vector to formatted file. ', &
            trim(fname), '. Exiting.'
       
       STOP
    end if
    
    write(fid) self%nx, self%ny, self%nz, self%gridType
    write(fid) self%x
    write(fid) self%y
    write(fid) self%z          
    
  end subroutine Write
  
  !
  !************************************************
  ! Boundary operations
  !************************************************
  !

  !** 
  ! SetAllBoundary
  !
  !*
  subroutine SetAllBoundary(self, c_in)
    implicit none
    ! Arguments
    class(cVector3D_SG_t), intent(inout) :: self
    complex(kind = prec)    , intent(in) :: c_in
    
    select case(self%gridType)
    case(EDGE)
       self%x(:, (/1, self%NdX(2)/), :) = c_in
       self%x(:, :, (/1, self%NdX(3)/)) = c_in
       self%y((/1, self%NdY(1)/), :, :) = c_in
       self%y(:, :, (/1, self%NdY(3)/)) = c_in
       self%z(:, (/1, self%NdZ(2)/), :) = c_in
       self%z((/1, self%NdZ(1)/), :, :) = c_in
       
    case(FACE)
       self%x((/1, self%NdX(1)/), :, :) = c_in
       self%y(:, (/1, self%NdY(2)/), :) = c_in
       self%z(:, :, (/1, self%NdZ(3)/)) = c_in

    case default
       write(*, *) 'ERROR:cVector3D_SG::SetAllBoundary: '
       write(*, *) '      Invalid grid type. Exiting.'
       
       STOP
    end select
    
  end subroutine SetAllBoundary
  
  !** 
  ! SetAllInterior
  !
  !*
  subroutine SetAllInterior(self, c_in)
    implicit none
    ! Arguments
    class(cVector3D_SG_t), intent(inout) :: self
    complex(kind = prec) , intent(in)    :: c_in
    
    select case(self%gridType)
    case(EDGE)
       self%x(:, 2:self%NdX(2)-1, :) = c_in
       self%x(:, :, 2:self%NdX(3)-1) = c_in
       self%y(2:self%NdY(1)-1, :, :) = c_in
       self%y(:, :, 2:self%NdY(3)-1) = c_in
       self%z(:, 2:self%NdZ(2)-1, :) = c_in
       self%z(2:self%NdZ(1)-1, :, :) = c_in
       
    case(FACE)
       self%x(2:self%NdX(1)-1, :, :) = c_in
       self%y(:, 2:self%NdY(2)-1, :) = c_in
       self%z(:, :, 2:self%NdZ(3)-1) = c_in
       
    case default
       write(*, *) 'ERROR:cVector3D_SG::SetAllInterior: '
       write(*, *) '      Invalid grid type. Exiting.'
       
       STOP
    end select
    
  end subroutine SetAllInterior

  !**
  ! SetOneBoundary
  !
  !*
  subroutine SetOneBoundary(self, bdry, c, int_only)
    implicit none
    ! Arguments
    class(cVector3D_SG_t), intent(inout) :: self
    character(*)         , intent(in)    :: bdry
    complex(kind = prec) , intent(in)    :: c
    logical, optional    , intent(in)    :: int_only
    ! Local variables
    logical :: int_only_p
    
    if (.not.present(int_only)) then
       int_only_p = .false.
    else
       int_only_p = int_only
    end if
    
    select case(self%gridType)
    case(EDGE)
       if (int_only_p) then
          select case (bdry)
          case('x1')
             self%z(1, 2:self%NdZ(2)-1, :) = c
             self%y(1, :, 2:self%NdY(3)-1) = c             
          case('x2')
             self%z(self%NdZ(1), 2:self%NdZ(2)-1, :) = c
             self%y(self%NdY(1), :, 2:self%NdY(3)-1) = c             
          case('y1')
             self%z(2:self%NdZ(1)-1, 1, :) = c
             self%x(:, 1, 2:self%NdX(3)-1) = c             
          case('y2')
             self%z(2:self%NdZ(1)-1, self%NdZ(2), :) = c
             self%x(:, self%NdX(2), 2:self%NdX(3)-1) = c             
          case('z1')
             self%x(:, 2:self%NdX(2)-1, 1) = c
             self%y(2:self%NdY(1)-1, :, 1) = c             
          case('z2')
             self%x(:, 2:self%NdX(2)-1, self%NdX(3)) = c
             self%y(2:self%NdY(1)-1, :, self%NdY(3)) = c             
          case('z1_x')
             self%x(:, 2:self%NdX(2)-1, 1) = c             
          case('z2_x')
             self%x(:, 2:self%NdX(2)-1, self%NdX(3)) = c             
          case('z1_y')
             self%y(2:self%NdY(1)-1, :, 1) = c             
          case('z2_y')
             self%y(2:self%NdY(1)-1, :, self%NdY(3)) = c             
          end select          
       else
          select case (bdry)
          case('x1')
             self%z(1, :, :) = c
             self%y(1, :, :) = c             
          case('x2')
             self%z(self%NdZ(1), :, :) = c
             self%y(self%NdY(1), :, :) = c             
          case('y1')
             self%z(:, 1, :) = c
             self%x(:, 1, :) = c             
          case('y2')
             self%z(:, self%NdZ(2), :) = c
             self%x(:, self%NdX(2), :) = c             
          case('z1')
             self%x(:, :, 1) = c
             self%y(:, :, 1) = c             
          case('z2')
             self%x(:, :, self%NdX(3)) = c
             self%y(:, :, self%NdY(3)) = c             
          case('z1_x')
             self%x(:, :, 1) = c             
          case('z2_x')
             self%x(:, :, self%NdX(3)) = c   
          case('z1_y')
             self%y(:, :, 1) = c  
          case('z2_y')
             self%y(:, :, self%NdY(3)) = c  
          end select
       end if
    case(FACE)
       select case (bdry)
       case('x1')
          self%x(1, :, :) = c          
       case('x2')
          self%x(self%NdX(1), :, :) = c          
       case('y1')
          self%y(:, 1, :) = c          
       case('y2')
          self%y(:, self%NdY(2), :) = c          
       case('z1')
          self%z(:, :, 1) = c          
       case('z2')
          self%z(:, :, self%NdZ(3)) = c          
       end select
    case default
       write(*, *) 'ERROR:cVector3D_SG::SetOneBoundary: '
       write(*, *) '      Invalid grid type. Exiting.'
       
       STOP
    end select
    
  end subroutine SetOneBoundary

  !**
  ! IntBdryIndices
  !*
  subroutine IntBdryIndices(self, ind_i, ind_b)
    implicit none
    ! Arguments
    class(cVector3D_SG_t), intent(in)  :: self
    integer, allocatable, intent(out) :: ind_i(:), ind_b(:)
    ! Local variables
    integer :: nVec(3), nVecT, nBdry, nb, ni, i
    complex(kind = prec), allocatable :: temp(:)
    type(cVector3D_SG_t) :: E
    
    if (self%isAllocated) then
       select type(grid => self%grid)
         class is(Grid3D_SG_t)
           E = cVector3D_SG_t(grid, self%gridType)
       end select
    else
       write(*, *) 'ERROR:cVector3D_SG_t:IntBdryIndices'
       write(*, *) '      Not allocated. Exiting.'
       
       STOP
    end if
    
    select case(self%gridType)
    case(EDGE)
       nVec(1) = size(E%x)
       nVec(2) = size(E%y)
       nVec(3) = size(E%z)
       nVecT = nVec(1) + nVec(2) + nVec(3)
       
       allocate(temp(nVecT))
       
       E%x(:, 1, :) = 1
       E%x(:, E%ny + 1, :) = 1
       E%x(:, :, 1) = 1
       E%x(:, :, E%nz + 1) = 1
       E%y(1, :, :) = 1
       E%y(E%nx + 1, :, :) = 1
       E%y(:, :, 1) = 1
       E%y(:, :, E%nz + 1) = 1
       E%z(1, :, :) = 1
       E%z(E%nx + 1, :, :) = 1
       E%z(:, 1, :) = 1
       E%z(:, E%ny + 1, :) = 1

       call E%GetArray(temp)       

    case(FACE)
       nVec(1) = size(E%x)
       nVec(2) = size(E%y)
       nVec(3) = size(E%z)
       nVecT = nVec(1) + nVec(2) + nVec(3)
       
       allocate(temp(nVecT))
       
       E%x(1, :, :) = 1
       E%x(E%nx + 1, :, :) = 1
       E%y(:, 1, :) = 1
       E%y(:, E%ny + 1, :) = 1
       E%z(:, :, 1) = 1
       E%z(:, :, E%nz + 1) = 1
       
       call E%GetArray(temp) 
       
    end select
    
    nBdry = 0
    do i = 1, nVecT
       nBdry = nBdry + nint(real(temp(i)))
    end do
    
    if (allocated(ind_i)) then
       deallocate(ind_i)
    end if
    
    allocate(ind_i(nVecT - nBdry))
    
    if (allocated(ind_b)) then
       deallocate(ind_b)
    end if
    
    allocate(ind_b(nBdry))
    
    nb = 0
    ni = 0
    do i = 1, nVecT
       if (nint(real(temp(i))).eq.1) then
          nb = nb + 1
          ind_b(nb) = i
       else
          ni = ni + 1
          ind_i(ni) = i
       end if
    end do
    
  end subroutine IntBdryIndices
  
  !**
  ! Boundary
  ! Returns a copy of this vector with all interior elements ser to zero.
  !*
  function Boundary(self) result(E)
    class(cVector3D_SG_t), intent(in) :: self
    ! Local variables
    class(cVector_t), allocatable :: E
    
    allocate(E, source = self)
    call E%SetAllInterior(C_ZERO)
    
  end function Boundary
  
  !**
  ! Interior
  ! Returns a copy of this vector with all boundary elements ser to zero.
  !*
  function Interior(self) result(E)
    class(cVector3D_SG_t), intent(in) :: self
    ! Local variables
    class(cVector_t), allocatable :: E

    allocate(E, source = self)
    call E%SetAllBoundary(C_ZERO)
    
  end function Interior
  
  !
  !************************************************
  ! Data access
  !************************************************
  !

  !**
  ! Length
  !
  !*
  function Length(self) result(n)
    implicit none
    ! Arguments
    class(cVector3D_SG_t), intent(in) :: self
    ! Local variables
    integer :: n
    
    n = self%Nxyz(1) + self%Nxyz(2) + self%Nxyz(3)
    
  end function Length

  !**
  ! GetArray
  !
  !*
  subroutine GetArray(self, v)
    implicit none
    ! Arguments
    class(cVector3D_SG_t), intent(in)  :: self
    complex(kind = prec), allocatable, intent(out) :: v(:)
    allocate (v(self%Length ()))
    v = (/reshape(self%x, (/self%Nxyz(1), 1/)), &
         reshape(self%y, (/self%Nxyz(2), 1/)), &
         reshape(self%z, (/self%Nxyz(3), 1/))/)
    
  end subroutine GetArray
  
  !**
  ! SetArray
  !
  !*
  subroutine SetArray(self, v)
    implicit none
    ! Arguments
    class(cVector3D_SG_t), intent(inout) :: self
    complex(kind = prec) , intent(in)    :: v(:)
    ! Local variables
    integer :: i1, i2
    
    ! Ex
    i1 = 1; i2 = self%Nxyz(1)
    self%x = reshape(v(i1:i2), self%NdX)
    
    ! Ey
    i1 = i2 + 1; i2 = i2 + self%Nxyz(2)
    self%y = reshape(v(i1:i2), self%NdY)
    
    ! Ez
    i1 = i2 + 1; i2 = i2 + self%Nxyz(3)
    self%z = reshape(v(i1:i2), self%NdZ)
    
  end subroutine SetArray
  
  !**
  ! SetVecComponents
  !
  !*
  subroutine SetVecComponents(self, xyz, &
       &                      xmin, xstep, xmax, &
       &                      ymin, ystep, ymax, &
       &                      zmin, zstep, zmax, c)
    implicit none
    ! Arguments
    class(cVector3D_SG_t), intent(inout) :: self
    character           , intent(in)    :: xyz
    integer, intent(in) :: xmin, xstep, xmax
    integer, intent(in) :: ymin, ystep, ymax
    integer, intent(in) :: zmin, zstep, zmax    
    real(kind = prec), intent (in) :: c
    ! Local variables
    integer :: x1, x2
    integer :: y1, y2
    integer :: z1, z2

    x1 = xmin; x2 = xmax
    y1 = ymin; y2 = ymax
    z1 = zmin; z2 = zmax

    select case(xyz)
    case('x')
       if (xmin == 0) x1 = self%NdX(1)
       if (xmax <= 0) x2 = self%NdX(1) + xmax
       
       if (ymin == 0) y1 = self%NdX(2)
       if (ymax <= 0) y2 = self%NdX(2) + ymax
       
       if (zmin == 0) z1 = self%NdX(3)
       if (zmax <= 0) z2 = self%NdX(3) + zmax
       
       self%x(x1:x2:xstep, y1:y2:ystep, z1:z2:zstep) = c

    case('y')       
       if (xmin == 0) x1 = self%NdY(1)
       if (xmax <= 0) x2 = self%NdY(1) + xmax

       if (ymin == 0) y1 = self%NdY(2)
       if (ymax <= 0) y2 = self%NdY(2) + ymax

       if (zmin == 0) z1 = self%NdY(3)
       if (zmax <= 0) z2 = self%NdY(3) + zmax

       self%y(x1:x2:xstep, y1:y2:ystep, z1:z2:zstep) = c

    case('z')       
       if (xmin == 0) x1 = self%NdZ(1)
       if (xmax <= 0) x2 = self%NdZ(1) + xmax

       if (ymin == 0) y1 = self%NdZ(2)
       if (ymax <= 0) y2 = self%NdZ(2) + ymax

       if (zmin == 0) z1 = self%NdZ(3)
       if (zmax <= 0) z2 = self%NdZ(3) + zmax

       self%z(x1:x2:xstep, y1:y2:ystep, z1:z2:zstep) = c

    case default
       write(*, *) 'ERROR:cVector3D_SG::SetVecComponents: '
       write(*, *) '      Invalid "xyz" argument. Exiting.'
       
       STOP
    end select
    
  end subroutine SetVecComponents

  !
  !************************************************
  ! Arithmetic/algebraic operations
  !************************************************
  !

  !**
  ! Zeros
  !*
  subroutine Zeros(self)
    implicit none
    class(cVector3D_SG_t), intent(inout) :: self

    self%x = R_ZERO
    self%y = R_ZERO
    self%z = R_ZERO
    
  end subroutine Zeros

  !**
  ! Add_1
  !*
  function Add_1(lhs, rhs) result(Eout)
    class(cVector3D_SG_t), intent(in) :: lhs
    class(cVector_t)     , intent(in) :: rhs
    class(cVector_t), allocatable :: Eout

    if (lhs%isCompatible(rhs)) then
       allocate(Eout, source = cVector3D_SG_t(lhs%grid, lhs%gridType))
       select type(Eout)
       class is(cVector3D_SG_t)
          select type(rhs)
          class is(cVector3D_SG_t)
             Eout%x = lhs%x + rhs%x
             Eout%y = lhs%y + rhs%y
             Eout%z = lhs%z + rhs%z             
          end select
       end select
    else
       write(*, *) 'ERROR:cVector3D_SG::Add_1'
       write(*, *) '  Incompatible inputs. Exiting.'

       STOP
    end if
  end function Add_1

  !**
  ! Sub_1
  !*
  function Sub_1(lhs, rhs) result(Eout)
    class(cVector3D_SG_t), intent(in) :: lhs
    class(cVector_t)     , intent(in) :: rhs
    class(cVector_t), allocatable :: Eout

    if (lhs%isCompatible(rhs)) then
       allocate(Eout, source = cVector3D_SG_t(lhs%grid, lhs%gridType))
       select type(Eout)
       class is(cVector3D_SG_t)
          select type(rhs)
          class is(cVector3D_SG_t)
             Eout%x = lhs%x - rhs%x
             Eout%y = lhs%y - rhs%y
             Eout%z = lhs%z - rhs%z             
          end select
       end select
    else
       write(*, *) 'ERROR:cVector3D_SG::Sub_1'
       write(*, *) '  Incompatible inputs. Exiting.'
       
       STOP
    end if
  end function Sub_1

  !**
  ! Mult_1
  !*
  function Mult_1(lhs, rhs) result(Eout)
    class(cVector3D_SG_t), intent(in) :: lhs
    class(cVector_t)     , intent(in) :: rhs
    class(cVector_t), allocatable :: Eout

    if (lhs%isCompatible(rhs)) then
       allocate(Eout, source = cVector3D_SG_t(lhs%grid, lhs%gridType))
       select type(Eout)
       class is(cVector3D_SG_t)
          select type(rhs)
          class is(cVector3D_SG_t)
             Eout%x = lhs%x * rhs%x
             Eout%y = lhs%y * rhs%y
             Eout%z = lhs%z * rhs%z             
          end select
       end select
    else
       write(*, *) 'ERROR:cVector3D_SG::mult_1'
       write(*, *) '  Incompatible inputs. Exiting.'

       STOP
    end if
  end function Mult_1

  !**
  ! Mult_2
  !*
  function Mult_2(c, self) result(Eout)
    complex(kind = prec) , intent(in) :: c
    class(cVector3D_SG_t), intent(in) :: self
    class(cVector_t), allocatable :: Eout

    allocate(Eout, source = cVector3D_SG_t(self%grid, self%gridType))
    select type(Eout)
    class is(cVector3D_SG_t)
       Eout%x = c * self%x
       Eout%y = c * self%y
       Eout%z = c * self%z             
    end select
    
  end function Mult_2

  !**
  ! Mult_3
  !*
  function Mult_3(lhs, rhs) result(Eout)
    class(cVector3D_SG_t), intent(in) :: lhs
    class(rVector_t)     , intent(in) :: rhs
    class(cVector_t), allocatable :: Eout

    if (lhs%isCompatible(rhs)) then
       allocate(Eout, source = cVector3D_SG_t(lhs%grid, lhs%gridType))
       select type(Eout)
       class is(cVector3D_SG_t)
          select type(rhs)
          class is(rVector3D_SG_t)
             Eout%x = lhs%x * rhs%x
             Eout%y = lhs%y * rhs%y
             Eout%z = lhs%z * rhs%z             
          end select
       end select
    else
       write(*, *) 'ERROR:cVector3D_SG::Mult_3'
       write(*, *) '  Incompatible inputs. Exiting.'

       STOP
    end if
  end function Mult_3

  !**
  ! Mult_s_1
  ! Subroutine version of Mult_2
  !*
  subroutine Mult_s_1(self, c)
    class(cVector3D_SG_t), intent(inout) :: self
    complex(kind = prec) , intent(in) :: c

    self%x = c * self%x
    self%y = c * self%y
    self%z = c * self%z             
  end subroutine Mult_s_1

  !**
  ! Div_1
  !*
  function Div_1(lhs, rhs) result(Eout)
    class(cVector3D_SG_t), intent(in) :: lhs
    class(cVector_t)     , intent(in) :: rhs
    class(cVector_t), allocatable :: Eout

    if (lhs%isCompatible(rhs)) then
       allocate(Eout, source = cVector3D_SG_t(lhs%grid, lhs%gridType))
       select type(Eout)
       class is(cVector3D_SG_t)
          select type(rhs)
          class is(cVector3D_SG_t)
             Eout%x = lhs%x / rhs%x
             Eout%y = lhs%y / rhs%y
             Eout%z = lhs%z / rhs%z             
          end select
       end select
    else
       write(*, *) 'ERROR:cVector3D_SG::Div_1'
       write(*, *) '  Incompatible inputs. Exiting.'

       STOP
    end if
  end function Div_1
  
  !**
  ! Div_2
  !*
  function Div_2(lhs, rhs) result(Eout)
    class(cVector3D_SG_t), intent(in) :: lhs
    class(rVector_t)     , intent(in) :: rhs
    class(cVector_t), allocatable :: Eout

    if (lhs%isCompatible(rhs)) then
       allocate(Eout, source = cVector3D_SG_t(lhs%grid, lhs%gridType))
       select type(Eout)
       class is(cVector3D_SG_t)
          select type(rhs)
          class is(rVector3D_SG_t)
             Eout%x = lhs%x / rhs%x
             Eout%y = lhs%y / rhs%y
             Eout%z = lhs%z / rhs%z             
          end select
       end select
    else
       write(*, *) 'ERROR:cVector3D_SG::Div_2'
       write(*, *) '  Incompatible inputs. Exiting.'
       
       STOP
    end if
  end function Div_2

  function dotProd(lhs, rhs) result(c)
    ! Arguments
    class(cVector3D_SG_t), intent(in) :: lhs
    class(cVector_t)     , intent(in) :: rhs
    ! Local variables
    complex(kind = prec) :: c
    integer :: nx, ny, nz
        
    c = C_ZERO
    
    if ((.not.lhs%isAllocated).or.(.not.rhs%isAllocated)) then
       write(*, *) 'ERROR:cVector3D_SG::dotProd: '
       write(*, *) '      Input vectors not allocated. Exiting.'

       STOP
    end if
    
    nx = lhs%nx
    ny = lhs%ny
    nz = lhs%nz
    
    if (lhs%isCompatible(rhs)) then
       select type(rhs)
       class is(cVector3D_SG_t)
          c = c + sum(conjg(lhs%x) * rhs%x)
          c = c + sum(conjg(lhs%y) * rhs%y)
          c = c + sum(conjg(lhs%z) * rhs%z)     
       end select
    else
       write(*, *) 'ERROR:cVector3D_SG:dotProd:'
       write(*, *) '      Incompatible input. Exiting'
       
       STOP
    end if
  end function dotProd
  
  !**
  ! InterpFunc
  ! Creates a Vector object containing weights needed for
  ! interpolation of xyz component of obj1 to location.
  !*
  subroutine InterpFunc(self, location, xyz, E)
    implicit none
    class(cVector3D_SG_t), intent(in) :: self
    real(kind = prec)    , intent(in) :: location(3)
    character            , intent(in) :: xyz
    class(cVector_t)     , intent(out), allocatable :: E
    ! Local variables
    real(kind = prec), allocatable, dimension(:) :: xC, yC, zC
    integer :: ix, iy, iz, i
    real(kind = prec) :: wx, wy, wz
    logical, dimension(:), allocatable :: tmp
    class(Grid3D_SG_t), pointer :: grid
    
    grid => self%grid
    
    select case(self%gridType)
    case(EDGE)
       allocate(E, source = cVector3D_SG_t(grid, EDGE))
       select case(xyz)
       case('x')
          allocate(xC(size(grid%delX)))
          allocate(yC(size(grid%dy + 1)))
          allocate(zC(size(grid%dz + 1)))
          
          xC = CumSum(grid%delX)
          yC = CumSum([0._prec, grid%dy])
          zC = CumSum([0._prec, grid%dz])
       case('y')
          allocate(xC(size(grid%dx + 1)))
          allocate(yC(size(grid%delY)))
          allocate(zC(size(grid%dz)))
          
          xC = CumSum([0._prec, grid%dx])
          yC = CumSum([grid%delY])
          zC = CumSum([0._prec, grid%dz])
       case('z')
          allocate(xC(size(grid%dx + 1)))
          allocate(yC(size(grid%dy + 1)))
          allocate(zC(size(grid%delZ)))
          
          xC = CumSum([0._prec, grid%dx])
          yC = CumSum([0._prec, grid%dy])
          zC = CumSum([grid%delZ])
       end select
              
    case(FACE)
       allocate(E, source = cVector3D_SG_t(self%grid, FACE))
       select case(xyz)
       case('x')
          allocate(xC(size(grid%dx + 1)))
          allocate(yC(size(grid%delY)))
          allocate(zC(size(grid%delZ)))
          
          xC = CumSum([0._prec, grid%dx])
          yC = CumSum([grid%delY])
          zC = CumSum([grid%delZ])
       case('y')
          allocate(xC(size(grid%delX)))
          allocate(yC(size(grid%dy + 1)))
          allocate(zC(size(grid%delZ)))
          
          xC = CumSum([grid%delX])
          yC = CumSum([0._prec, grid%dy])
          zC = CumSum([grid%delZ])
       case('z')
          allocate(xC(size(grid%delX)))
          allocate(yC(size(grid%delY)))
          allocate(zC(size(grid%dz + 1)))
          
          xC = CumSum([grid%delX])
          yC = CumSum([grid%delY])
          zC = CumSum([0._prec, grid%dz])
       end select
    end select
    
    xC = xC - grid%ox
    yC = yC - grid%oy
    zC = zC - sum(grid%dz(1:grid%nzAir)) - grid%oz
    
    ! Find indicies of "minimal corner"
    ! Code below should be replaced by intrinsic
    ! 'findloc' if it is avaiable.
    tmp = location(1) > xC
    do i = size(tmp), 1, -1 
       if (tmp(i)) then
          ix = i
          exit
       end if
    end do
    tmp = location(2) > yC
    do i = size(tmp), 1, -1 
       if (tmp(i)) then
          iy = i
          exit
       end if
    end do
    tmp = location(3) > zC
    do i = size(tmp), 1, -1 
       if (tmp(i)) then
          iz = i
          exit
       end if
    end do
    !ix = findloc(location(1) > xC, .true., back = .true., dim = 1)
    !iy = findloc(location(2) > yC, .true., back = .true., dim = 1)
    !iz = findloc(location(3) > zC, .true., back = .true., dim = 1)
    ! Find weights
    wx = (xC(ix + 1) - location(1))/(xC(ix + 1) - xC(ix))
    wy = (yC(iy + 1) - location(2))/(yC(iy + 1) - yC(iy))
    wz = (zC(iz + 1) - location(3))/(zC(iz + 1) - zC(iz))

    select type(E)
    class is(cVector3D_SG_t)
       select case(xyz)
       case('x')
          E%x(ix,iy,iz) = wx*wy*wz
          E%x(ix+1,iy,iz) = (1-wx)*wy*wz
          E%x(ix,iy+1,iz) = wx*(1-wy)*wz
          E%x(ix,iy,iz+1) = wx*wy*(1-wz)
          E%x(ix,iy+1,iz+1) = wx*(1-wy)*(1-wz)
          E%x(ix+1,iy,iz+1) = (1-wx)*wy*(1-wz)
          E%x(ix+1,iy+1,iz) = (1-wx)*(1-wy)*wz
          E%x(ix+1,iy+1,iz+1) = (1-wx)*(1-wy)*(1-wz)
          
       case('y')
          E%y(ix,iy,iz) = wx*wy*wz
          E%y(ix+1,iy,iz) = (1-wx)*wy*wz
          E%y(ix,iy+1,iz) = wx*(1-wy)*wz
          E%y(ix,iy,iz+1) = wx*wy*(1-wz)
          E%y(ix,iy+1,iz+1) = wx*(1-wy)*(1-wz)
          E%y(ix+1,iy,iz+1) = (1-wx)*wy*(1-wz)
          E%y(ix+1,iy+1,iz) = (1-wx)*(1-wy)*wz
          E%y(ix+1,iy+1,iz+1) = (1-wx)*(1-wy)*(1-wz)
          
       case('z')
          E%z(ix,iy,iz) = wx*wy*wz
          E%z(ix+1,iy,iz) = (1-wx)*wy*wz
          E%z(ix,iy+1,iz) = wx*(1-wy)*wz
          E%z(ix,iy,iz+1) = wx*wy*(1-wz)
          E%z(ix,iy+1,iz+1) = wx*(1-wy)*(1-wz)
          E%z(ix+1,iy,iz+1) = (1-wx)*wy*(1-wz)
          E%z(ix+1,iy+1,iz) = (1-wx)*(1-wy)*wz
          E%z(ix+1,iy+1,iz+1) = (1-wx)*(1-wy)*(1-wz)
       end select
    end select
    
  end subroutine InterpFunc

  subroutine CopyFrom(self, rhs)
    ! Arguments
    class(cVector3D_SG_t), intent(inout) :: self
    class(cVector_t)     , intent(in)    :: rhs

    if (.not.rhs%isAllocated) then
       write(*, *) 'ERROR:cVector3D_SG::CopyFrom'
       write(*, *) '  Inputs not allocated. Exiting.'
       
       STOP
    end if

    select type(rhs)
    class is(cVector3D_SG_t)
       if (.not.self%isAllocated) then
          self%grid => rhs%grid
          self%gridType = rhs%gridType
          self%nx = rhs%nx
          self%ny = rhs%ny
          self%nz = rhs%nz
          self%NdX = rhs%NdX
          self%NdY = rhs%NdY
          self%NdZ = rhs%NdZ
          self%Nxyz = rhs%Nxyz
          self%x = rhs%x
          self%y = rhs%y
          self%z = rhs%z
          
          self%isAllocated = .true.
       else
          if (self%isCompatible(rhs)) then
             self%grid => rhs%grid
             self%gridType = rhs%gridType
             self%nx = rhs%nx
             self%ny = rhs%ny
             self%nz = rhs%nz
             self%NdX = rhs%NdX
             self%NdY = rhs%NdY
             self%NdZ = rhs%NdZ
             self%Nxyz = rhs%Nxyz
             self%x = rhs%x
             self%y = rhs%y
             self%z = rhs%z
             
             self%isAllocated = .true.
          else
             write(*, *) 'ERROR:cVector3D_SG::CopyFrom'
             write(*, *) '  Incompatible input. Exiting.'
             
             STOP
          end if
       end if
    class default
       write(*, *) 'ERROR:cVector3D_SG::CopyFrom:'
       write(*, *) '      Incompatible input type. Exiting.'
       STOP           
    end select
  end subroutine CopyFrom

  function isCompatible_1(self, rhs) result(status)
    ! Arguments
    class(cVector3D_SG_t), intent(in) :: self
    class(cVector_t), intent(in) :: rhs
    logical :: status
    
    status = .false.
    select type(rhs)
    class is(cVector3D_SG_t)
       if((self%nx == rhs%nx).and.(self%ny == rhs%ny).and.(self%nz == rhs%nz)) then 
          if (self%gridType == rhs%gridType) then
             status = .true.
          end if
       end if
    end select
  end function isCompatible_1
  
  function isCompatible_2(self, rhs) result(status)
    ! Arguments
    class(cVector3D_SG_t), intent(in) :: self
    class(rVector_t)     , intent(in) :: rhs
    logical :: status
    
    status = .false.
    select type(rhs)
    class is(rVector3D_SG_t)
       if((self%nx == rhs%nx).and.(self%ny == rhs%ny).and.(self%nz == rhs%nz)) then 
          if (self%gridType == rhs%gridType) then
             status = .true.
          end if
       end if
    end select
  end function IsCompatible_2
  
end module cVector3D_SG
