!**
! SUMMARY
!
! Standard cartesian grid scalars.
!*
module rScalar3D_SG
   !
   use Constants
   use Grid
   use Grid3D_SG
   use rScalar
   !
   !**
   ! Type rScalar3D_SG_t defines real scalar fields defined on
   ! cells/nodes of a regular 3D cartesian grid.
   !*
   type, extends(rScalar_t) :: rScalar3D_SG_t
       !**
       ! Pointer to parent grid
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

       integer, dimension(3) :: NdV = 0
       integer :: Nxyz = 0
       
       real(kind = prec), allocatable, dimension(:, :, :) :: v

    contains
       private
       
       !**
       ! Initialization and finalization
       !*
       final :: rScalar3D_SG_dtor
       
       !**
       ! Input/Output
       !*
       procedure, public :: read   => readRScalar3D_SG
       procedure, public :: write => writeRScalar3D_SG
       
       !**
       ! Boundary operations
       !*
       procedure, public :: setAllBoundary => setAllBoundaryRScalar3D_SG
       procedure, public :: setOneBoundary => setOneBoundaryRScalar3D_SG
       procedure, public :: setAllInterior => setAllInteriorRScalar3D_SG
       procedure, public :: intBdryIndices => intBdryIndicesRScalar3D_SG
       !procedure, public :: Boundary
       !procedure, public :: Interior

       !**
       ! Data access
       !*
       procedure, public :: length => lengthRScalar3D_SG
       procedure, public :: getArray => getArrayRScalar3D_SG
       procedure, public :: setArray => setArrayRScalar3D_SG
       procedure, public :: setVecComponents => setVecComponentsRScalar3D_SG
       
       !**
       ! Arithmetic/algebraic operations
       !*
       
       procedure, public :: zeros => zerosRScalar3D_SG
       procedure, public :: add1 => add1RScalar3D_SG
       procedure, public :: sub1 => sub1RScalar3D_SG
       procedure, public :: mult1 => mult1RScalar3D_SG
       procedure, public :: div1 => div1RScalar3D_SG
       procedure, public :: dotProd => dotProdRScalar3D_SG
       
       !**
       ! Arithmetic operations, subroutine versions: first argument is overwritrten
       !*
       !procedure, public :: divS1 => divS1RScalar3D_SG
       !procedure, public :: multS1 => multS1RScalar3D_SG
       !procedure, public :: multS2 => multS2RScalar3D_SG

       !**
       ! Miscellaneous
       !*

       procedure, public :: isCompatible1 => isCompatible1RScalar3D_SG
       procedure, public :: copyFrom => copyFromCompatibleRScalar3D_SG
               
   end type rScalar3D_SG_t
   
   ! Constructors for Scalar3d_csg_real_t
   interface rScalar3D_SG_t
       module procedure rScalar3D_SG_ctor
   end interface rScalar3D_SG_t

contains
   
   !**
   ! Parameterized constructor.
   !
   ! Arguments
   !    igrid       Underlying grid.
   !    grid_type   Definied in GridDef.f90
   !
   !*
   function rScalar3D_SG_ctor(igrid, gridType) result (E)
      implicit none
      ! Arguments
      class(Grid3D_SG_t), target , intent(in) :: igrid
      character(*)                      , intent(in) :: gridType
      ! Local variables
      integer :: nx, ny, nz, nzAir, nz_earth
      integer :: status
      type(rScalar3D_SG_t) :: E
      
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
          write(*, *) 'ERROR:rScalar3D_SG::ctor:'
          write(*, *) '         Unrecognized grid type: ', gridType, '. Exiting.'
          
          STOP          
      end if
      
      E%isAllocated = E%isAllocated.and.(status .EQ. 0)
      if (E%isAllocated) then
          E%v = R_ZERO          
      else          
          write(*, *) 'ERROR:rScalar3D_SG::ctor'
          write(*, *) '         Unable to allocate - invalid grid supplied. Exiting.'
          
          STOP          
      end if
      
      E%Nxyz = product(E%NdV)
      
   end function rScalar3D_SG_ctor
   
   !**
   ! Destructor.
   !
   !*
   subroutine rScalar3D_SG_dtor(E)
      implicit none
      ! Arguments
      type(rScalar3D_SG_t) :: E
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
      
   end subroutine rScalar3D_SG_dtor
   
   !
   !***************
   ! Input/Output
   !***************
   !

   !**
   ! readRScalar3D_SG
   !
   !*
   subroutine readRScalar3D_SG(self, fid, ftype)
      implicit none
      ! Arguments
      class(rScalar3D_SG_t), intent (inout) :: self
      integer                     , intent (in)      :: fid
      character(*)             , intent (in), optional :: ftype
      ! Local variables
      integer :: Nx, Ny, Nz
      character(80) :: gridType
      integer :: i, j, k, k1, k2, istat
      real(kind = prec), allocatable, dimension (:) :: temp      
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
          write(*, *) 'ERROR:rScalar3D_SG::readRScalar3D_SG: '
          write(*, *) '         Unable to read vector from unformatted file ', &
                  trim (fname), '. Exiting.'
          
          STOP
          
      else if ((index(isbinary, 'no') > 0.or.index(isbinary, 'NO') > 0) &
             .and.binary) then
          write(*, *) 'ERROR:rScalar3D_SG::readRScalar3D_SG: '
          write(*, *) '         Unable to read vector from formatted file ', &
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
               write(*, *) 'ERROR:rScalar3D_SG::readRScalar3D_SG: '
               write(*, *) '         While reading the ', i, 'th block. Exiting.'               

               STOP
          else if (k1 > k2) then
               write(*, *) 'WARNING:rScalar3D_SG_t::readRScalar3D_SG: '
               write(*, *) '            Block ', i, ' will be ignored.'
          end if
          
          do j = Nx, 1, -1
               read(fid, *, iostat = istat) temp
               
               if (istat /= 0) then
                   write(*, *) 'ERROR:rScalar3D_SG_t::readRScalar3D_SG: '
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
   end subroutine readRScalar3D_SG

   !**
   ! writeRScalar3D_SG
   !
   !*
   subroutine writeRScalar3D_SG(self, fid, ftype)
      implicit none
      ! Arguments
      class(rScalar3D_SG_t), intent(in) :: self
      integer                     , intent(in) :: fid
      character(*)             , intent(in), optional :: ftype
      ! Local variables
      integer :: Nx, Ny, Nz
      integer :: i, j, k, k1, k2, istat
      real(kind = prec), allocatable, dimension(:, :) :: temp
      logical :: ok, hasname, binary
      character(80) :: fname, isbinary, gridType

      if (.not.self%isAllocated) then
          write(0, *) 'ERROR:rScalar3D_SG::writeRScalar3D_SG: '
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
          write(*, *) 'ERROR:rScalar3D_SG::writeRScalar3D_SG:'
          write(*, *) '         Unable to write vector to unformatted file ', &
                  trim(fname), '. Exiting.'          

          STOP          
      else if ((index(isbinary,'no') > 0.or.index(isbinary,'NO') > 0) &
             .and.binary) then
          write(*, *) 'ERROR:rScalar3D_SG_t::writeRScalar3D_SG: '
          write(*, *) ' Unable to write vector to formatted file ', &
                  trim(fname), '. Exiting.'
          
          STOP
      end if
      
      gridType = self%gridType

      if (binary) then
          write(fid) self%nx, self%ny, self%nz, gridType
          write(fid) self%v          
          return
      end if
      
      !**
      ! ASCII format
      !*
      write(fid, '(3i5,a10)', iostat = istat) self%nx, &
             self%ny, self%nz, trim(gridType)

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
               write(*, *) 'ERROR:rScalar3D_SG::writeRScalar3D_SG:'
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
      
   end subroutine writeRScalar3D_SG
   
   !
   !************************
   ! Boundary operations
   !************************
   !
   
   !**
   ! setAllBoundaryRScalar3D_SG
   !
   !*
   subroutine setAllBoundaryRScalar3D_SG(self, c_in)
      implicit none
      ! Arguments
      class(rScalar3D_SG_t), intent(inout) :: self
      real(kind = prec)      , intent(in)      :: c_in
      ! Local variables
      character(80) :: gridType
      
      gridType = self%gridType
      select case(gridType)
      case (CORNER) 
          self%v((/1, self%NdV(1)/), :, :) = c_in
          self%v(:, (/1, self%NdV(2)/), :) = c_in
          self%v(:, :, (/1, self%NdV(3)/)) = c_in
          
      case default
          write (*, *) 'ERROR:rScalar3D_SG_t::setAllBoundaryRScalar3D_SG:'
          write (*, *) '         Grid type not recognized. Exiting.'
          
          STOP          
      end select
      
   end subroutine setAllBoundaryRScalar3D_SG

   !**
   ! setOneBoundaryRScalar3D_SG
   !
   !*
   subroutine setOneBoundaryRScalar3D_SG(self, bdry, c, int_only)
      implicit none
      ! Arguments
      class(rScalar3D_SG_t), intent(inout) :: self
      character(*)             , intent(in) :: bdry
      real(kind = prec)      , intent(in) :: c
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
          write (*, *) 'ERROR:rScalar3D_SG_t::setOneBoundaryRScalar3D_SG:'
          write (*, *) '         Invalid grid type. Exiting.'
          
          STOP          
      end select
      
   end subroutine setOneBoundaryRScalar3D_SG

   !**
   ! setAllInteriorRScalar3D_SG
   !*
   subroutine setAllInteriorRScalar3D_SG(self, c_in)
      implicit none
      ! Arguments
      class(rScalar3D_SG_t), intent(inout) :: self
      real(kind = prec)      , intent(in)      :: c_in
   end subroutine setAllInteriorRScalar3D_SG
   
   !**
   ! intBdryIndicesRScalar3D_SG
   !
   !*
   subroutine intBdryIndicesRScalar3D_SG(self, ind_i, ind_b)
      implicit none
      ! Arguments
      class(rScalar3D_SG_t) , intent(in)   :: self
      integer, allocatable, intent(out) :: ind_i(:), ind_b(:)
      ! Local variables
      integer :: nVec(3), nVecT, nBdry, nb, ni, i
      real(kind = prec), allocatable :: temp(:)
      type(rScalar3D_SG_t) :: phi
      character(80) :: gridType

      if (self%isAllocated) then
	    select type(grid => self%grid)
               class is(Grid3D_SG_t)
                   phi = rScalar3D_SG_t(grid, self%gridType)
	    end select
      else
          write(*, *) 'ERROR:rScalar3D_SG::intBdryIndicesRScalar3D_SG:'
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
          nBdry = nBdry + nint(temp(i))
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
          if (nint(temp(i)).eq.1) then
               nb = nb+1
               ind_b(nb) = i
          else
               ni = ni+1
               ind_i(ni) = i
          end if
      end do
   end subroutine intBdryIndicesRScalar3D_SG

   ! !**
   ! ! Boundary
   ! ! Returns a copy of this Scalar with
   ! ! all interior elements ser to zero.
   ! !*
   ! function Boundary(self) result(E)
   !    class(rScalar3D_SG_t), intent(in) :: self
   !    ! Local variables
   !    type(rScalar3D_SG_t) :: E
      
   !    E = self
   !    call E%setAllInteriorRScalar3D_SG(R_ZERO)
   ! end function Boundary
   
   ! !**
   ! ! Interior
   ! ! Returns a copy of this vector with
   ! ! all boundary elements ser to zero.
   ! !*
   ! function Interior(self) result(E)
   !    class(rScalar3D_SG_t), intent(in) :: self
   !    ! Local variables
   !    type(rScalar3D_SG_t) :: E
      
   !    E = self
   !    call E%setAllBoundaryRScalar3D_SG(R_ZERO)
      
   ! end function Interior
   
   !
   !***************
   ! Data access
   !***************
   !
   
   !**
   ! lengthRScalar3D_SG
   !
   !*
   function lengthRScalar3D_SG(self) result(n)
      implicit none
      ! Arguments
      class(rScalar3D_SG_t), intent(in) :: self
      ! Local variables
      integer :: n
      
      n = self%Nxyz
      
   end function lengthRScalar3D_SG

   !**
   ! getArrayRScalar3D_SG
   !
   !*
   subroutine getArrayRScalar3D_SG(self, v)
      implicit none
      ! Arguments
      class(rScalar3D_SG_t)             , intent(in)   :: self
      real(kind = prec), allocatable, intent(out) :: v(:)
      !
      allocate(v(self%length()))      
      v = (/reshape(self%v, (/self%Nxyz, 1/))/)
      
   end subroutine getArrayRScalar3D_SG
   
   !**
   ! setArrayRScalar3D_SG
   !
   !*
   subroutine setArrayRScalar3D_SG(self, v)
      implicit none
      ! Arguments
      class(rScalar3D_SG_t), intent(inout) :: self
      real(kind = prec)      , intent(in)      :: v(:)
      
      self%v = reshape(v, (/self%NdV(1), self%NdV(2), self%NdV(3)/))
   end subroutine setArrayRScalar3D_SG

   !
   !************************************************
   ! Arithmetic operations
   !************************************************
   !
   
   !**
   ! zerosRScalar3D_SG
   ! zerosRScalar3D_SG array.
   !
   !*
   subroutine zerosRScalar3D_SG(self)
      implicit none
      ! Arguments
      class(rScalar3D_SG_t), intent(inout) :: self
      
      if (.not.self%isAllocated) then
          write (*, *) 'ERROR:rScalar3D_SG::zerosRScalar3D_SG:'
          write (*, *) '         Not allocated. Exiting.'
          
          STOP
      end if
      
      self%v = R_ZERO
   end subroutine zerosRScalar3D_SG

   !**
   ! add1RScalar3D_SG
   !*
   function add1RScalar3D_SG(lhs, rhs) result(Eout)
      class(rScalar3D_SG_t), intent(in) :: lhs
      class(rScalar_t)       , intent(in) :: rhs
      class(rScalar_t), allocatable :: Eout

      if (lhs%isCompatible(rhs)) then
          allocate(Eout, source = rScalar3D_SG_t(lhs%grid, lhs%gridType))
          select type(Eout)
          class is(rScalar3D_SG_t)
               select type(rhs)
               class is(rScalar3D_SG_t)
                   Eout%v = lhs%v + rhs%v
               end select
          end select
      else
          write(*, *) 'ERROR:rScalar3D_SG::add1RScalar3D_SG:'
          write(*, *) '   Incompatible inputs. Exiting.'
          
          STOP
      end if
   end function add1RScalar3D_SG
   
   !**
   ! sub1RScalar3D_SG
   !*
   function sub1RScalar3D_SG(lhs, rhs) result(Eout)
      class(rScalar3D_SG_t), intent(in)   :: lhs
      class(rScalar_t)       , intent(in)   :: rhs
      class(rScalar_t), allocatable :: Eout

      if (lhs%isCompatible(rhs)) then
          allocate(Eout, source = rScalar3D_SG_t(lhs%grid, lhs%gridType))
          select type(Eout)
          class is(rScalar3D_SG_t)
               select type(rhs)
               class is(rScalar3D_SG_t)
                   Eout%v = lhs%v - rhs%v
               end select
          end select
      else
          write(*, *) 'ERROR:rScalar3D_SG::sub1RScalar3D_SG'
          write(*, *) '   Incompatible inputs. Exiting.'

          STOP
      end if
   end function sub1RScalar3D_SG
   
   !**
   ! mult1RScalar3D_SG
   !*
   function mult1RScalar3D_SG(lhs, rhs) result(Eout)
      class(rScalar3D_SG_t), intent(in)   :: lhs
      class(rScalar_t)       , intent(in)   :: rhs
      class(rScalar_t), allocatable :: Eout

      if (lhs%isCompatible(rhs)) then
          allocate(Eout, source = rScalar3D_SG_t(lhs%grid, lhs%gridType))
          select type(Eout)
          class is(rScalar3D_SG_t)
               select type(rhs)
               class is(rScalar3D_SG_t)
                   Eout%v = lhs%v * rhs%v
               end select
          end select
      else
          write(*, *) 'ERROR:rScalar3D_SG::mult1RScalar3D_SG'
          write(*, *) '   Incompatible inputs. Exiting.'

          STOP
      end if
   end function mult1RScalar3D_SG
   
   !**
   ! div1RScalar3D_SG
   !*
   function div1RScalar3D_SG(lhs, rhs) result(Eout)
      class(rScalar3D_SG_t), intent(in)   :: lhs
      class(rScalar_t)       , intent(in)   :: rhs
      class(rScalar_t), allocatable :: Eout

      if (lhs%isCompatible(rhs)) then
          allocate(Eout, source = rScalar3D_SG_t(lhs%grid, lhs%gridType))
          select type(Eout)
          class is(rScalar3D_SG_t)
               select type(rhs)
               class is(rScalar3D_SG_t)
                   Eout%v = lhs%v / rhs%v
               end select
          end select
      else
          write(*, *) 'ERROR:rScalar3D_SG::div1RScalar3D_SG'
          write(*, *) '   Incompatible inputs. Exiting.'

          STOP
      end if
   end function div1RScalar3D_SG

   function dotProdRScalar3D_SG(lhs, rhs) result(c)
      class(rScalar3D_SG_t), intent(in) :: lhs
      class(rScalar_t)       , intent(in) :: rhs
      real(kind = prec) :: c

      if (lhs%isCompatible(rhs)) then
          select type(rhs)
          class is(rScalar3D_SG_t)
               c = sum(lhs%v * rhs%v)
          end select
      else
          write(*, *) 'ERROR:rScalar3D_SG::dotProdRScalar3D_SG'
          write(*, *) '   Incompatible inputs. Exiting.'

          STOP
      end if
      
   end function dotProdRScalar3D_SG
   
   !
   !********************
   ! Miscellaneous
   !********************
   !

   !**
   ! setVecComponentsRScalar3D_SG
   !*
   subroutine setVecComponentsRScalar3D_SG(self, xyz, &
          &       xmin, xstep, xmax, &
          &       ymin, ystep, ymax, &
          &       zmin, zstep, zmax, c)
      implicit none
      ! Arguments
      class(rScalar3D_SG_t), intent(inout) :: self
      character                  , intent(in)      :: xyz
      integer                     , intent(in)      :: xmin, xstep, xmax
      integer                     , intent(in)      :: ymin, ystep, ymax
      integer                     , intent(in)      :: zmin, zstep, zmax      
      real(kind = prec), intent(in) :: c
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
      
   end subroutine setVecComponentsRScalar3D_SG

   function isCompatible1RScalar3D_SG(self, rhs) result(status)
      ! Arguments
      class(rScalar3D_SG_t), intent(in) :: self
      class(rScalar_t), intent(in) :: rhs
      logical :: status
      
      status = .false.
      !if (same_type_as(self, rhs)) then

          select type(rhs)
          class is(rScalar3D_SG_t)

               if((self%nx == rhs%nx).and.(self%ny == rhs%ny).and.(self%nz == rhs%nz)) then 

                   if (self%gridType == rhs%gridType) then

                        status = .true.
                   end if
               end if
          end select
      !end if
   end function isCompatible1RScalar3D_SG

   subroutine copyFromCompatibleRScalar3D_SG(self, rhs)
      ! Arguments
      class(rScalar3D_SG_t), intent(inout) :: self
      class(rScalar_t)       , intent(in)      :: rhs

      if (.not.rhs%isAllocated) then
          write(*, *) 'ERROR:rScalar3D_SG::copyFromCompatibleRScalar3D_SG'
          write(*, *) '   Input not allocated. Exiting.'
          
          STOP
      end if

      select type(rhs)
      class is(rScalar3D_SG_t)
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
                   write(*, *) 'ERROR:rScalar3D_SG::copyFromCompatibleRScalar3D_SG'
                   write(*, *) '   Incompatible input. Exiting.'
                   
                   STOP
               end if
          end if
      class default
          write(*, *) 'ERROR:rScalar3D_SG::copyFromCompatibleRScalar3D_SG:'
          write(*, *) '         Incompatible input type. Exiting.'
          STOP 
      end select
   end subroutine copyFromCompatibleRScalar3D_SG
   
end module rScalar3D_SG
