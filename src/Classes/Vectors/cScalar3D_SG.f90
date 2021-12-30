!**
! SUMMARY
!
! Standard cartesian grid scalars.
!*
module cScalar3D_SG
   use Constants
   use Grid
   use Grid3D_SG
   use cScalar
   use rScalar3D_SG
   !**
   ! Type cScalar3D_SG_t defines real scalar fields defined on
   ! cells/nodes of a regular 3D cartesian grid.
   !*
   type, extends(cScalar_t) :: cScalar3D_SG_t
       !**
       ! Pointer to parent grid
       class(Grid3D_SG_t), pointer :: grid => null()
       
       !**
       ! Store the intention of the use in a character
       ! string defined in GridDef as a parameter: EDGE
       ! or FACE are two possibilities.
       !*
       character(len = 4) :: gridType = ""
       
       !**
       ! Grid Dimensions:
       ! nx is grid dimension (number of cells) in the x-direction
       ! ny is grid dimension (number of cells) in the y-direction
       ! nz is grid dimension (number of cells) in the z-direction:
       !*
       integer :: nx = 0, ny = 0, nz = 0

       integer, dimension(3) :: NdV = 0
       integer :: Nxyz = 0
       
       complex(kind = prec), allocatable, dimension(:, :, :) :: v

    contains
       private
       
       !**
       ! Initialization and finalization
       !*
       final :: cScalar3D_SG_dtor
       
       !**
       ! Input/Output
       !*
       procedure, public :: read => readCScalar3D_SG
       procedure, public :: write => writeCScalar3D_SG
       
       !**
       ! Boundary operations
       !*
       procedure, public :: setAllBoundary => setAllBoundaryCScalar3D_SG
       procedure, public :: setOneBoundary => setOneBoundaryCScalar3D_SG
       procedure, public :: setAllInterior => setAllInteriorCScalar3D_SG
       procedure, public :: intBdryIndices => intBdryIndicesCScalar3D_SG
       !procedure, public :: Boundary
       !procedure, public :: Interior

       !**
       ! Data access
       !*
       procedure, public :: length => lengthCScalar3D_SG
       procedure, public :: getArray => getArrayCScalar3D_SG
       procedure, public :: setArray => setArrayCScalar3D_SG
       procedure, public :: setVecComponents => setVecComponentsCScalar3D_SG
       
       !**
       ! Arithmetic/algebraic operations
       !*
       
       ! Overriden methods
       procedure, public :: zeros => zerosCScalar3D_SG
       procedure, public :: add1 => add1CScalar3D_SG
       procedure, public :: sub1 => sub1CScalar3D_SG
       procedure, public :: mult1 => mult1CScalar3D_SG
       procedure, public :: mult2 => mult2CScalar3D_SG
	   procedure, public :: mult3 => mult3CScalar3D_SG
	   !
       procedure, public :: mults1 => mults1CScalar3D_SG
       procedure, public :: mults2 => mults2CScalar3D_SG
       procedure, public :: mults3 => mults3CScalar3D_SG
	   !
       procedure, public :: div1 => div1CScalar3D_SG
       !
       procedure, public :: divs3 => divs3CScalar3D_SG
       procedure, public :: dotProd => dotProdCScalar3D_SG
       
       !   linear combinations for cScalars
       procedure, public :: linCombS => linCombSCcalar3D_SG
       procedure, public :: scMultAddS =>  SCMultAddCscalar3D_SG

       !**
       ! Miscellaneous
       !*
       procedure, public :: isCompatible1 => isCompatible1CScalar3D_SG
       procedure, public :: isCompatible2 => isCompatible2CScalar3D_SG
       generic :: isCompatible => isCompatible1, isCompatible2

       ! Overriden methods
       procedure, public :: copyFrom => copyFromCScalar3D_SG
       
   end type cScalar3D_SG_t
   
   ! Constructors for Scalar3d_csg_real_t
   interface cScalar3D_SG_t
       module procedure cScalar3D_SG_ctor
   end interface cScalar3D_SG_t

contains
   
   !**
   ! Parameterized constructor.
   !
   ! Arguments
   !    igrid       Underlying grid.
   !    grid_type   Definied in GridDef.f90
   !
   !*
   function cScalar3D_SG_ctor(igrid, gridType) result (E)
      implicit none
      ! Arguments
      class(Grid3D_SG_t), target , intent(in) :: igrid
      character(*)                      , intent(in) :: gridType
      ! Local variables
      integer :: nx, ny, nz, nzAir, nz_earth
      integer :: status
      type(cScalar3D_SG_t) :: E
      
      E%grid => igrid
      E%gridType = gridType
      
      ! Grid dimensions
      call igrid%GetDimensions(nx, ny, nz, nzAir)
      nz_earth = nz - nzAir
      
      E%nx = nx
      E%ny = ny
      E%nz = nz
      
      ! allocate memory for x,y,z ;
      ! E%allocated will be true if all allocations succeed
      E%isAllocated = .true.
            
      if (gridType == CORNER) then          
          allocate(E%v(nx + 1, ny + 1, nz + 1), STAT = status)   
          E%NdV = (/E%nx + 1, E%ny + 1, E%nz + 1/)
          
      else if (gridType == CENTER) then          
          allocate(E%v(nx, ny, nz), STAT = status) 
          E%NdV = (/E%nx, E%ny, E%nz/)
          
      else if (gridType == CELL_EARTH) then
          E%nz = nz_earth
          allocate(E%v(nx, ny, nz_earth), STAT = status)
          E%NdV = (/nx, ny, nz_earth/)
          
      else          
          write(*, *) 'ERROR:cScalar3D_SG_t::ctor:'
          write(*, *) '         Unrecognized grid type: ', gridType, '. Exiting.'
          
          STOP          
      end if
      
      E%isAllocated = E%isAllocated.and.(status .EQ. 0)
      if (E%isAllocated) then
          E%v = R_ZERO          
      else          
          write(*, *) 'ERROR:cScalar3D_SG_t::ctor:'
          write(*, *) '         Unable to allocate rScalar - invalid grid supplied. Exiting.'
          
          STOP          
      end if
      
      E%Nxyz = product(E%NdV)
      
   end function cScalar3D_SG_ctor
   
   !**
   ! Destructor.
   !
   !*
   subroutine cScalar3D_SG_dtor(E)
      implicit none
      ! Arguments
      type(cScalar3D_SG_t) :: E
      ! Local variables
      integer :: status
      
      if (E%isAllocated) then 
          deallocate(E%v, STAT = status)
      end if
      
      E%grid => null()
            
      E%nx = 0
      E%ny = 0
      E%nz = 0
      
      E%gridType = ''
      E%isAllocated = .false.
      
   end subroutine cScalar3D_SG_dtor
   
   !
   !***************
   ! Input/Output
   !***************
   !

   !**
   ! readCScalar3D_SG
   !
   !*
   subroutine readCScalar3D_SG(self, fid, ftype)
      implicit none
      ! Arguments
      class(cScalar3D_SG_t), intent (inout) :: self
      integer                     , intent (in)      :: fid
      character(*)             , intent (in), optional :: ftype
      ! Local variables
      integer :: Nx, Ny, Nz
      character(4) :: gridType
      integer :: i, j, k, k1, k2, istat
      complex(kind = prec), allocatable, dimension (:) :: temp      
      logical :: ok, hasname, binary
      character(80) :: fname, isbinary
      
      if (.not. present (ftype)) then
          binary = .false.
      else if (index (ftype, 'b') > 0) then
          binary = .true.
      else
          binary = .false.
      end if
      
      inquire(fid, opened = ok, named = hasname, &
             name = fname, unformatted = isbinary)
      
      ! check that the file is unformatted if binary, formatted if ascii
      if ((index(isbinary, 'yes') > 0.or.index(isbinary, 'YES') > 0) &
             .and. .not.binary) then          
          write(*, *) 'ERROR:cScalar3D_SG_t::readCScalar3D_SG: '
          write(*, *) '         Unable to read scalar from unformatted file ', &
                  trim(fname), '.Exiting.'

          STOP
      else if ((index(isbinary, 'no') > 0.or.index(isbinary, 'NO') > 0) &
             .and.binary) then
          write(*, *) 'ERROR:cScalar3D_SG_t::readCScalar3D_SG: '
          write(*, *) '         Unable to read scalar from formatted file ', &
                  trim(fname), '. Exiting.'
          
          STOP
      end if
      
      if (binary) then
          ! read binary from unformatted files
          read(fid) self%Nx, self%Ny, self%Nz, gridType
          read(fid) self%v          
      end if
      
      Nx = size(self%v, 1)
      Ny = size(self%v, 2)
      Nz = size(self%v, 3)
      
      allocate(temp(Ny), STAT = istat)
          
      i = 1
      do
          read(fid, *, iostat = istat) k1, k2
          if (istat /= 0) exit
          
          if ((k1 < 0) .or. (k2 > Nz)) then
               write(*, *) 'ERROR:cScalar3D_SG::readCScalar3D_SG: '
               write(*, *) '         While reading the ', i, 'th block. Exiting.'
               
               STOP
          else if (k1 > k2) then
               write(*, *) 'WARNING:cScalar3D_SG::readCScalar3D_SG: '
               write(*, *) '            Block ', i, ' will be ignored.'
          end if
          
          do j = Nx, 1, -1
               read(fid, *, iostat = istat) temp
               
               if (istat /= 0) then
                   write(*, *) 'ERROR:cScalar3D_SG::readCScalar3D_SG: '
                   write(*, *) '         While reading the ', j, 'th row in ', i,'th block. Exiting.'
                   
                   STOP
               end if
               
               do k = k1, k2
                   self%v(j, :, k) = temp
               end do
          end do
          
          if (k == Nz) exit
          
          i = i + 1
          
      end do
   end subroutine readCScalar3D_SG

   !**
   ! writeCScalar3D_SG
   !
   !*
   subroutine writeCScalar3D_SG(self, fid, ftype)
      implicit none
      ! Arguments
      class(cScalar3D_SG_t), intent(in) :: self
      integer                     , intent(in) :: fid
      character(*)             , intent(in), optional :: ftype
      ! Local variables
      integer :: Nx, Ny, Nz
      integer :: i, j, k, k1, k2, istat
      complex(kind = prec), allocatable, dimension(:, :) :: temp
      logical :: ok, hasname, binary
      character(80) :: fname, isbinary

      if (.not.self%isAllocated) then
          write(0, *) 'ERROR:cScalar3D_SG::writeCScalar3D_SG: '
          write(0, *) '         Not allocated. Exiting.'
          
          STOP
      end if
      
      if (.not.present(ftype)) then
          binary = .false.
      else if (index(ftype, 'b') > 0) then
          binary = .true.
      else
          binary = .false.
      end if
      
      inquire(fid, opened = ok, named = hasname, &
             name = fname, unformatted = isbinary)
      
      if ((index(isbinary, 'yes') > 0.or.index(isbinary, 'YES') > 0) &
             .and..not.binary) then          
          write(*, *) 'ERROR:cScalar3D_SG::writeCScalar3D_SG: '
          write(*, *) '         Unable to write vector to unformatted file ', &
                  trim(fname), '. Exiting.'
          
          STOP
      else if ((index(isbinary,'no') > 0.or.index(isbinary,'NO') > 0) &
             .and.binary) then
          write(*, *) 'ERROR:cScalar3D_SG::writeCScalar3D_SG: '
          write(*, *) ' Unable to write vector to formatted file ', &
                  trim(fname), '. Exiting.'
          
          STOP
      end if
      
      if (binary) then
          write(fid) self%nx, self%ny, self%nz, self%gridType
          write(fid) self%v          
          return
      end if
      
      !**
      ! ASCII format
      !*
      write(fid, '(3i5,a10)', iostat = istat) self%nx, &
             self%ny, self%nz, trim(self%gridType)

      Nx = size(self%v, 1)
      Ny = size(self%v, 2)
      Nz = size(self%v, 3)
          
      allocate(temp(Nx, Ny), STAT = istat)
          
      k1 = 1
      do
          k2 = Nz
          do k = k1, Nz - 1
               temp = abs(self%v(:, :, k + 1) - self%v(:, :, k))
               if (maxval(real(temp)) > TOL6) then
                   k2 = k
                   exit
               end if
          end do
          
          write(fid, '(2i5)', iostat = istat) k1, k2
          
          if (istat /= 0) then
               write(*, *) 'ERROR:cScalar3D_SG::writeCScalar3D_SG: '
               write(*, *) '         Failed while writing to file. Exiting.'
               
               STOP
          end if
          
          temp = self%v(:, :, k1)
          
          do i = Nx, 1, -1
               do j = 1, Ny
                   write(fid, '(es13.5)', iostat = istat, &
                           advance = 'no') self%v(i, j, k1)
               end do
               write(fid, *)
          end do
          
          k1 = k2 + 1
          
          if (k1 > Nz) exit
      end do
      
   end subroutine writeCScalar3D_SG
   
   !
   !************************
   ! Boundary operations
   !************************
   !
   
   !**
   ! setAllBoundaryCScalar3D_SG
   !
   !*
   subroutine setAllBoundaryCScalar3D_SG(self, c_in)
      implicit none
      ! Arguments
      class(cScalar3D_SG_t), intent(inout) :: self
      complex(kind = prec) , intent(in)      :: c_in
      ! Local variables
      character(80) :: gridType
      
      gridType = self%gridType
      select case(gridType)
      case (CORNER) 
          self%v((/1, self%NdV(1)/), :, :) = c_in
          self%v(:, (/1, self%NdV(2)/), :) = c_in
          self%v(:, :, (/1, self%NdV(3)/)) = c_in
          
      case default
          write (*, *) 'ERROR:cScalar3D_SG_t::setAllBoundaryCScalar3D_SG: '
          write (*, *) '         Grid type not recognized. Exiting. '
          
          STOP          
      end select
      
   end subroutine setAllBoundaryCScalar3D_SG

   !**
   ! setOneBoundaryCScalar3D_SG
   !
   !*
   subroutine setOneBoundaryCScalar3D_SG(self, bdry, c, int_only)
      implicit none
      ! Arguments
      class(cScalar3D_SG_t), intent(inout) :: self
      character(*)             , intent(in) :: bdry
      complex(kind = prec) , intent(in) :: c
      logical                     , intent(in), optional :: int_only
      ! Local variables
      logical :: int_only_p
      character(80) :: gridType
      
      if (.not.present (int_only)) then
          int_only_p = .false.
      else 
          int_only_p = int_only
      end if
      
      gridType = self%gridType
      select case(gridType)
      case (CORNER)          
          if (int_only_p) then               
               select case (bdry)
               case('x1')
                   self%v(1, 2:self%NdV(2)-1, 2:self%NdV(3)-1) = c 
               case('x2')
                   self%v(self%NdV(1), 2:self%NdV(2)-1, 2:self%NdV(3)-1) = c                   
               case('y1')
                   self%v(2:self%NdV(1)-1, 1, 2:self%NdV(3)-1) = c
               case('y2')
                   self%v(2:self%NdV(1)-1, self%NdV(2), 2:self%NdV(3)-1) = c                   
               case('z1')
                   self%v(2:self%NdV(1)-1, 2:self%NdV(2)-1, 1) = c                   
               case('z2')
                   self%v(2:self%NdV(1)-1, 2:self%NdV(2)-1, self%NdV(3)) = c
               end select
          else
               select case(bdry)
               case('x1')
                   self%v(1, :, :) = c                   
               case('x2')
                   self%v(self%NdV(1), :, :) = c
               case('y1')
                   self%v(:, 1, :) = c                   
               case('y2')
                   self%v(:, self%NdV(2), :) = c
               case('z1')
                   self%v(:, :, 1) = c                   
               case('z2')
                   self%v(:, :, self%NdV(3)) = c
               end select
          end if
          
      case(FACE)
          select case(bdry)
          case('x1')
               self%v(1, :, :) = c
          case('x2')
               self%v(self%NdV(1), :, :) = c
          case('y1')
               self%v(:, 1, :) = c
          case('y2')
               self%v(:, self%NdV(2), :) = c
          case('z1')
               self%v(:, :, 1) = c
          case('z2')
               self%v(:, :, self%NdV(3)) = c
          end select
          
      case default          
          write (*, *) 'ERROR:cScalar3D_SG_t::setOneBoundaryCScalar3D_SG: '
          write (*, *) '         Invalid grid type. Exiting.'
          
          STOP          
      end select
      
   end subroutine setOneBoundaryCScalar3D_SG

   !**
   ! setAllInteriorCScalar3D_SG
   !*
   subroutine setAllInteriorCScalar3D_SG(self, c_in)
      implicit none
      ! Arguments
      class(cScalar3D_SG_t), intent(inout) :: self
      complex(kind = prec) , intent(in)      :: c_in
   end subroutine setAllInteriorCScalar3D_SG
   
   !**
   ! intBdryIndicesCScalar3D_SG
   !
   !*
   subroutine intBdryIndicesCScalar3D_SG(self, ind_i, ind_b)
      implicit none
      ! Arguments
      class(cScalar3D_SG_t) , intent(in)   :: self
      integer, allocatable, intent(out) :: ind_i(:), ind_b(:)
      ! Local variables
      integer :: nVec(3), nVecT, nBdry, nb, ni, i
      complex(kind = prec), allocatable :: temp(:)
      type(cScalar3D_SG_t) :: phi
      character(80) :: gridType

      if (self%isAllocated) then
	    select type(grid => self%grid)
               class is(Grid3D_SG_t)
                   phi = cScalar3D_SG_t(grid, self%gridType)
	    end select
      else
          write(*, *) 'ERROR:cScalar3D_SG::intBdryIndicesCScalar3D_SG:'
          write(*, *) '         Not allocated. Exiting.'

          STOP
      end if
      
      gridType = self%gridType
      select case(gridType)
      case(CORNER)
          nVecT = size(phi%v)
               
          allocate(temp(nVecT))
               
          phi%v(1,:,:) = 1
          phi%v(phi%nx+1,:,:) = 1
          phi%v(:,1,:) = 1
          phi%v(:,phi%ny+1,:) = 1
          phi%v(:,:,1) = 1
          phi%v(:,:,phi%nz+1) = 1
          
          call phi%getArray(temp)

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
               nb = nb+1
               ind_b(nb) = i
          else
               ni = ni+1
               ind_i(ni) = i
          end if
      end do
   end subroutine intBdryIndicesCScalar3D_SG

   ! !**
   ! ! Boundary
   ! ! Returns a copy of this Scalar with
   ! ! all interior elements ser to zero.
   ! !*
   ! function Boundary(self) result(E)
   !    class(cScalar3D_SG_t), intent(in) :: self
   !    ! Local variables
   !    type(cScalar3D_SG_t) :: E
      
   !    E = self
   !    call E%setAllInteriorCScalar3D_SG(R_ZERO)
   ! end function Boundary
   
   ! !**
   ! ! Interior
   ! ! Returns a copy of this vector with
   ! ! all boundary elements ser to zero.
   ! !*
   ! function Interior(self) result(E)
   !    class(cScalar3D_SG_t), intent(in) :: self
   !    ! Local variables
   !    type(cScalar3D_SG_t) :: E
      
   !    E = self
   !    call E%setAllBoundaryCScalar3D_SG(R_ZERO)
      
   ! end function Interior
   
   !
   !***************
   ! Data access
   !***************
   !
   
   !**
   ! lengthCScalar3D_SG
   !
   !*
   function lengthCScalar3D_SG(self) result(n)
      implicit none
      ! Arguments
      class(cScalar3D_SG_t), intent(in) :: self
      ! Local variables
      integer :: n
      
      n = self%Nxyz
      
   end function lengthCScalar3D_SG

   !**
   ! getArrayCScalar3D_SG
   !
   !*
   subroutine getArrayCScalar3D_SG(self, v)
      implicit none
      ! Arguments
      class(cScalar3D_SG_t)                  , intent(in)   :: self
      complex(kind = prec), allocatable, intent(out) :: v(:)
      !
      allocate(v(self%length()))      
      v = (/reshape(self%v, (/self%Nxyz, 1/))/)
      
   end subroutine getArrayCScalar3D_SG
   
   !**
   ! setArrayCScalar3D_SG
   !
   !*
   subroutine setArrayCScalar3D_SG(self, v)
      implicit none
      ! Arguments
      class(cScalar3D_SG_t), intent(inout) :: self
      complex(kind = prec) , intent(in)      :: v(:)
      
      self%v = reshape(v, (/self%NdV(1), self%NdV(2), self%NdV(3)/))
   end subroutine setArrayCScalar3D_SG
   
   !
   !************************************************
   ! Arithmetic operations
   !************************************************
   !
   
   !**
   ! zerosCScalar3D_SG
   ! zerosCScalar3D_SG array.
   !
   !*
   subroutine zerosCScalar3D_SG(self)
      implicit none
      ! Arguments
      class(cScalar3D_SG_t), intent(inout) :: self
      
      if (.not.self%isAllocated) then
          write (*, *) 'ERROR:cScalar3D_SG_t::zerosCScalar3D_SG: '
          write (*, *) '         Not allocated. Exiting.'
          
          STOP
      end if
      
      self%v = C_ZERO

   end subroutine zerosCScalar3D_SG

   !**
   ! add1CScalar3D_SG
   !*
   function add1CScalar3D_SG(lhs, rhs) result(Eout)
      class(cScalar3D_SG_t), intent(in) :: lhs
      class(cScalar_t)       , intent(in) :: rhs
      class(cScalar_t), allocatable :: Eout

      !if (lhs%isCompatible(rhs)) then
          allocate(Eout, source = cScalar3D_SG_t(lhs%grid, lhs%gridType))
          select type(Eout)
          class is(cScalar3D_SG_t)
               select type(rhs)
               class is(cScalar3D_SG_t)
                   Eout%v = lhs%v + rhs%v
               end select
          end select
      !else
          !write(*, *) 'ERROR:cScalar3D_SG::add1CScalar3D_SG'
          !write(*, *) '   Incompatible inputs. Exiting.'
          
          !STOP
      !end if
   end function add1CScalar3D_SG
   
   !**
   ! sub1CScalar3D_SG
   !*
   function sub1CScalar3D_SG(lhs, rhs) result(Eout)
      class(cScalar3D_SG_t), intent(in)   :: lhs
      class(cScalar_t)       , intent(in)   :: rhs
      class(cScalar_t), allocatable :: Eout

      if (lhs%isCompatible(rhs)) then
          allocate(Eout, source = cScalar3D_SG_t(lhs%grid, lhs%gridType))
          select type(Eout)
          class is(cScalar3D_SG_t)
               select type(rhs)
               class is(cScalar3D_SG_t)
                   Eout%v = lhs%v - rhs%v
               end select
          end select
      else
          write(*, *) 'ERROR:cScalar3D_SG::sub1CScalar3D_SG'
          write(*, *) '   Incompatible inputs. Exiting.'
          
          STOP
      end if
   end function sub1CScalar3D_SG
   
   !**
   ! mult1CScalar3D_SG
   !*
   function mult1CScalar3D_SG(lhs, rhs) result(Eout)
      class(cScalar3D_SG_t), intent(in)   :: lhs
      class(cScalar_t)       , intent(in)   :: rhs
      class(cScalar_t), allocatable :: Eout

      if (lhs%isCompatible(rhs)) then
          allocate(Eout, source = cScalar3D_SG_t(lhs%grid, lhs%gridType))
          select type(Eout)
          class is(cScalar3D_SG_t)
               select type(rhs)
               class is(cScalar3D_SG_t)
                   Eout%v = lhs%v * rhs%v
               end select
          end select
      else
          write(*, *) 'ERROR:cScalar3D_SG::mult1CScalar3D_SG'
          write(*, *) '   Incompatible inputs. Exiting.'
          
          STOP
      end if
   end function mult1CScalar3D_SG
   
   !**
   ! mult2CScalar3D_SG
   !*
   function mult2CScalar3D_SG(lhs, rhs) result(Eout)
      class(cScalar3D_SG_t), intent(in)   :: lhs
      complex( kind=prec ), intent(in) :: rhs
      class(cScalar_t), allocatable :: Eout

          allocate(Eout, source = cScalar3D_SG_t(lhs%grid, lhs%gridType))
          select type(Eout)
          class is(cScalar3D_SG_t)
               !
               Eout%v = lhs%v * rhs
               !
          end select
   end function mult2CScalar3D_SG
   
   !**
   ! mult3CScalar3D_SG
   !*
   function mult3CScalar3D_SG(lhs, rhs) result(Eout)
      class(cScalar3D_SG_t), intent(in)   :: lhs
      class(rScalar_t)       , intent(in)   :: rhs
      class(cScalar_t), allocatable :: Eout

      !if (lhs%isCompatibleRScalar3D_SG(rhs)) then
          allocate(Eout, source = cScalar3D_SG_t(lhs%grid, lhs%gridType))
          select type(Eout)
          class is(cScalar3D_SG_t)
               select type(rhs)
               class is(rScalar3D_SG_t)
                   Eout%v = lhs%v * rhs%v
               end select
          end select
      !else
          !write(*, *) 'ERROR:rScalar3D_SG::mult3CScalar3D_SG'
          !write(*, *) '   Incompatible inputs. Exiting.'

          !STOP
      !end if
   end function mult3CScalar3D_SG
   !**
   ! mults2CScalar3D_SG
   !*
   !  subroutine versions -- these overwrite input
   !**
   ! mults1CScalar3D_SG
   !*
   subroutine mults1CScalar3D_SG(lhs, rhs)
      class(cScalar3D_SG_t), intent(inout)   :: lhs
      class(cScalar_t)       , intent(in)   :: rhs

      if (lhs%isCompatible(rhs)) then
         select type(rhs)
           class is(cScalar3D_SG_t)
               lhs%v = lhs%v * rhs%v
         end select
      else
          write(*, *) 'ERROR:cScalar3D_SG::mults1CScalar3D_SG'
          write(*, *) '   Incompatible inputs. Exiting.'
          STOP
      end if
   end subroutine mults1CScalar3D_SG
   !**
   !   mults2CScalar3D_SG
   !**
   subroutine mults2CScalar3D_SG(lhs, c)
      class(cScalar3D_SG_t), intent(inout)   :: lhs
      complex( kind=prec ), intent(in) :: c

         lhs%v = lhs%v * c

   end subroutine mults2CScalar3D_SG

   !**
   !
   !*
   subroutine mults3CScalar3D_SG(lhs, rhs)
      class(cScalar3D_SG_t), intent(inout)   :: lhs
      class(rScalar_t)       , intent(in)   :: rhs

      if (lhs%isCompatible(rhs)) then
         select type(rhs)
           class is(rScalar3D_SG_t)
               lhs%v = lhs%v * rhs%v
         end select
      else
          write(*, *) 'ERROR:cScalar3D_SG::mults3CScalar3D_SG'
          write(*, *) '   Incompatible inputs. Exiting.'
          STOP
      end if
   end subroutine mults3CScalar3D_SG
   !**
   !
   !*
   subroutine divs3CScalar3D_SG(lhs, rhs)
      class(cScalar3D_SG_t), intent(inout)   :: lhs
      class(rScalar_t)       , intent(in)   :: rhs

      if (lhs%isCompatible(rhs)) then
         select type(rhs)
           class is(rScalar3D_SG_t)
               lhs%v = lhs%v / rhs%v
         end select
      else
          write(*, *) 'ERROR:cScalar3D_SG::mults3CScalar3D_SG'
          write(*, *) '   Incompatible inputs. Exiting.'
          STOP
      end if
   end subroutine divs3CScalar3D_SG
   !**
   ! div1CScalar3D_SG
   !*
   function div1CScalar3D_SG(lhs, rhs) result(Eout)
      class(cScalar3D_SG_t), intent(in)   :: lhs
      class(cScalar_t)       , intent(in)   :: rhs
      class(cScalar_t), allocatable :: Eout

      if (lhs%isCompatible(rhs)) then
          allocate(Eout, source = cScalar3D_SG_t(lhs%grid, lhs%gridType))
          select type(Eout)
          class is(cScalar3D_SG_t)
               select type(rhs)
               class is(cScalar3D_SG_t)
                   Eout%v = lhs%v / rhs%v
               end select
          end select
      else
          write(*, *) 'ERROR:cScalar3D_SG::div1CScalar3D_SG'
          write(*, *) '   Incompatible inputs. Exiting.'
          
          STOP
      end if
   end function div1CScalar3D_SG
   !**
   ! linCombSCScalar3D_SG
   !*
   subroutine linCombSCScalar3D_SG(lhs, rhs, c1, c2)
      class(cScalar3D_SG_t), intent(inout) :: lhs
      class(cScalar_t)       , intent(in) :: rhs
      complex(kind=prec), intent(in)  :: c1, c2

      !  linear combination, in place: lhs = c1*lhs+c2*rhs
       if (lhs%isCompatible(rhs)) then
           select type(rhs)
           class is(cScalar3D_SG_t)
              lhs%v = c1*lhs%v + c2*rhs%v
           end select
       else
           write(*, *) 'ERROR:cScalar3D_SG::linCombScScalar3D_SG'
           write(*, *) '   Incompatible inputs. Exiting.'
           STOP
       end if
    end subroutine linCombSCScalar3D_SG
   !**
   ! scMultAddSCScalar3D_SG
   !*
   subroutine scMultAddSCScalar3D_SG(lhs, rhs, c)
      !   returns rhs = rhs + c*lhs
      class(cScalar3D_SG_t), intent(in) :: lhs
      class(cScalar_t)       , intent(inout) :: rhs
      complex(kind=prec), intent(in)  :: c

       if (lhs%isCompatible(rhs)) then
           select type(rhs)
           class is(cScalar3D_SG_t)
              rhs%v = rhs%v + c*lhs%v
           end select
       else
           write(*, *) 'ERROR:cScalar3D_SG:scMultAddScScalar3D_SG'
           write(*, *) '   Incompatible inputs. Exiting.'
           STOP
       end if
    end subroutine scMultAddSCScalar3D_SG
    !
    !*********
    !
    function dotProdCScalar3D_SG(lhs, rhs) result(c)
      class(cScalar3D_SG_t), intent(in) :: lhs
      class(cScalar_t)       , intent(in) :: rhs
      real(kind = prec) :: c

      if (lhs%isCompatible(rhs)) then
          select type(rhs)
          class is(cScalar3D_SG_t)
               c = sum(conjg(lhs%v) * rhs%v)                         
          end select
      else
          write(*, *) 'ERROR:cScalar3D_SG::dotProdCScalar3D_SG'
          write(*, *) '   Incompatible inputs. Exiting.'
          
          STOP
      end if
      
   end function dotProdCScalar3D_SG
   
   !
   !********************
   ! Miscellaneous
   !********************
   !

   !**
   ! setVecComponentsCScalar3D_SG
   !*
   subroutine setVecComponentsCScalar3D_SG(self, xyz, &
          &       xmin, xstep, xmax, &
          &       ymin, ystep, ymax, &
          &       zmin, zstep, zmax, c)
      implicit none
      ! Arguments
      class(cScalar3D_SG_t), intent(inout) :: self
      character                  , intent(in)      :: xyz
      integer                     , intent(in)      :: xmin, xstep, xmax
      integer                     , intent(in)      :: ymin, ystep, ymax
      integer                     , intent(in)      :: zmin, zstep, zmax      
      complex(kind = prec), intent(in) :: c
      ! Local variables
      integer :: x1, x2
      integer :: y1, y2
      integer :: z1, z2
      character(len = 20) :: xstr
      
      x1 = xmin; x2 = xmax
      y1 = ymin; y2 = ymax
      z1 = zmin; z2 = zmax
      
      if (xmin == 0) x1 = self%NdV(1)
      if (xmax <= 0) x2 = self%NdV(1) + xmax
      
      if (ymin == 0) y1 = self%NdV(2)
      if (ymax <= 0) y2 = self%NdV(2) + ymax
      
      if (zmin == 0) z1 = self%NdV(3)
      if (zmax <= 0) z2 = self%NdV(3) + zmax
      
      self%v(x1:x2:xstep, y1:y2:ystep, z1:z2:zstep) = c
      
   end subroutine setVecComponentsCScalar3D_SG

    !**
   ! linCombSCcalar3D_SG
   !*
   subroutine linCombSCcalar3D_SG(lhs, rhs, c1, c2)
     !   lhs = c1*lhs+c2*rhs
      class(cScalar3D_SG_t), intent(inout) :: lhs
      class(cScalar_t)       , intent(in) :: rhs
      complex(kind=prec), intent(in)  :: c1, c2

      !  linear combination, in place: lhs = c1*lhs+c2*rhs
      if (lhs%isCompatible(rhs)) then
          select type(rhs)
          class is(cScalar3D_SG_t)
             lhs%v = c1*lhs%v + c2*rhs%v
          end select
      else
          write(*, *) 'ERROR:cScalar3D_SG_t::linCombS'
          write(*, *) '   Incompatible inputs. Exiting.'

          STOP
      end if
   end subroutine linCombSCcalar3D_SG
   !*
   ! SCMultAddCscalar3D_SG
   !*
   subroutine SCMultAddCscalar3D_SG(lhs, rhs, c)
     !   returns rhs = rhs + c*lhs
     class(cScalar3D_SG_t), intent(in) :: lhs
     class(cScalar_t)       , intent(inout) :: rhs
     complex(kind=prec), intent(in)  :: c
       
    if (lhs%isCompatible(rhs)) then
       select type(rhs)
       class is(cScalar3D_SG_t) 
          rhs%v = rhs%v + c*lhs%v 
       end select
    else
      write(*, *) 'ERROR:cVector3D_SG::scMultAdd'
      write(*, *) '   Incompatible inputs. Exiting.'
      STOP
   end if
   end subroutine SCMultAddCscalar3D_SG
	
   function isCompatible1Cscalar3D_SG(self, rhs) result(status)
      ! Arguments
      class(cScalar3D_SG_t), intent(in) :: self
      class(cScalar_t), intent(in) :: rhs
      logical :: status
      
      status = .false.
      select type(rhs)
         class is(cScalar3D_SG_t)
            if((self%nx == rhs%nx).and.(self%ny == rhs%ny).and.(self%nz == rhs%nz)) then 
               if (self%gridType == rhs%gridType) then
                  status = .true.
               end if
            end if
        end select
   end function isCompatible1CScalar3D_SG
   !
   function isCompatible2Cscalar3D_SG(self, rhs) result(status)
      ! Arguments
      class(cScalar3D_SG_t), intent(in) :: self
      class(rScalar_t), intent(in) :: rhs
      logical :: status
      
      status = .false.
      select type(rhs)
         class is(rScalar3D_SG_t)
            if((self%nx == rhs%nx).and.(self%ny == rhs%ny).and.(self%nz == rhs%nz)) then 
               if (self%gridType == rhs%gridType) then
                    status = .true.
               end if
            end if
      end select
   end function isCompatible2CScalar3D_SG

   subroutine copyFromCScalar3D_SG(self, rhs)
      ! Arguments
      class(cScalar3D_SG_t), intent(inout) :: self
      class(cScalar_t)       , intent(in)      :: rhs

      if (.not.rhs%isAllocated) then
          write(*, *) 'ERROR:cScalar3D_SG::copyFromCScalar3D_SG'
          write(*, *) '   Input not allocated. Exiting.'
          
          STOP
      end if

      select type(rhs)
      class is(cScalar3D_SG_t)
          if (.not.self%isAllocated) then
               self%grid => rhs%grid
               self%gridType = rhs%gridType
               self%nx = rhs%nx
               self%ny = rhs%ny
               self%nz = rhs%nz
               self%NdV = rhs%NdV
               self%Nxyz = rhs%Nxyz
               self%v = rhs%v
               self%isAllocated = .true.          
          else
               if (self%isCompatible(rhs)) then
                   self%grid => rhs%grid
                   self%gridType = rhs%gridType
                   self%nx = rhs%nx
                   self%ny = rhs%ny
                   self%nz = rhs%nz
                   self%NdV = rhs%NdV
                   self%Nxyz = rhs%Nxyz
                   self%v = rhs%v
                   self%isAllocated = .true.
               else
                   write(*, *) 'ERROR:cScalar3D_SG::copyFromCScalar3D_SG'
                   write(*, *) '   Incompatible input. Exiting.'
                   
                   STOP
               end if
          end if
      class default
          write(*, *) 'ERROR:cScalar3D_SG::copyFromCScalar3D_SG:'
          write(*, *) '         Incompatible input type. Exiting.'
          STOP 
      end select
   end subroutine copyFromCScalar3D_SG

   subroutine printCScalar3D_SG( self, io_unit, title )
      implicit none
      !
      ! Arguments
      class( cScalar3D_SG_t ), intent( in ) :: self
      integer, intent( in ), optional       :: io_unit
      character(*), intent( in ), optional       :: title
      !
      integer :: ix, iy, iz,fid
      !
      if( present( io_unit ) ) then
         fid = io_unit
      else
         fid = 6   !   usually this will work to write to standard output
      endif
      if(present(title)) then
        write(fid,*) title
      end if

      write( fid, * ) self%nx, self%ny, self%nz
      !
      !
      write(fid,*) 'scalar field'
      do ix = 1, self%nx
          do iy = 1, self%ny
              do iz = 1, self%nz
                  if( self%v( ix, iy, iz ) /= 0 ) then
                     write(fid,*) ix,iy,iz, ":[", self%v( ix, iy, iz ), "]"
                  endif
              enddo
          enddo
      enddo

   end subroutine printCScalar3D_SG


end module cScalar3D_SG
