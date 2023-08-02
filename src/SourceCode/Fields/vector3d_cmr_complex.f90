!**
!* This file is part of the ModEM modeling and inversion package.
!
! LICENSING information
!
!* Copyright (C) 2020 ModEM research group.
!* Contact: http://
!
!* GNU General Public License Usage
!* This file may be used under the terms of the GNU
!* General Public License version 3.0 as published by the Free Software
!* Foundation and appearing in the file LICENSE.GPL included in the
!* packaging of this file.  Please review the following information to
!* ensure the GNU General Public License version 3.0 requirements will be
!* met: http://www.gnu.org/copyleft/gpl.html.
!
! SUMMARY
!
! This module specializes the abstract Vector3d_complex class for
! complex vector fields on a multi-resolution staggered grid.
!
!*
module Vector3d_cmr_complex
  use GridDef, only : grid_t
  use GridCalc, only : edgeLength
  use sg_vector, only : cvector, create_cvector
  use sg_vector, only : rvector, create_rvector
  use vecTranslate, only : getRVector, setRVector, setlimits
  use vecTranslate, only : getCVector, setCVector

  use MR_constants
  
  use Grid3d_cmr
  use Vector3d_csg_complex

  private 

  public :: &
       assignment (=)
  
  interface assignment (=)
     module procedure operator_assignment
  end interface assignment (=)

  !**
  ! Type vector defines vector for either edge or 
  ! face in a multi-resolution staggered grid as
  ! a complex field.
  !*
  type, public :: Vector3d_cmr_complex_t
     !**
     ! store the intention of the use in a character
     ! string defined in GridDef as a parameter: EDGE
     ! or FACE are two possibilities.
     !*
     character (len = 80) :: grid_type

     ! allocated:  .TRUE.  
     logical :: is_allocated = .false.

     !**
     ! pointer to parent grid
     class (Grid3d_cmr_t), pointer :: grid_ptr

     class (Vector3d_csg_complex_t), allocatable :: sub_vectors(:)

     ! List of all active (interior + boundary) indices
     integer, dimension (:), allocatable :: ind_active
     integer, dimension (:), allocatable :: ind_interior
     integer, dimension (:), allocatable :: ind_boundaries

   contains

     !**
     ! Initialization and finalization
     !*
     procedure, public :: Initialize => Vector3d_cmr_complex_t_initialize
     final :: Vector3d_cmr_complex_t_dtor

     !**
     ! Input/Output
     !*
     procedure, public :: Write => Write_
     procedure, public :: Read  => Read_

     !**
     ! Boundary operations
     !*
     procedure, public :: Set_all_boundary => Set_all_boundary_
     procedure, public :: Set_one_boundary => Set_one_boundary_
     procedure, public :: Get_int_bdry_indices => Get_int_bdry_indices_

     procedure, public :: Set_active_int_boundary => Set_active_int_boundary_

     !**
     ! Data access
     !*
     procedure, public :: length    => length_
     procedure, public :: Get_array => Get_array_
     procedure, public :: Set_array => Set_array_
     procedure, public :: Copy_from => Copy_from_

     procedure, public :: length_full => length_full_
     procedure, public :: Get_full    => Get_full_
     procedure, public :: Set_full    => Set_full_

     !**
     ! Arithmetic operations
     !*
     procedure, public :: Zero        => Zero_
     procedure, public :: Add         => Add_
     procedure, public :: Subtract    => Subtract_
     procedure, public :: Multiply_by_scalar => Multiply_by_scalar_
     procedure, public :: Multiply_by => Multiply_by_
     procedure, public :: Divide_by   => Divide_by_

     !**
     ! Algebraic operations
     !*
     procedure, public :: Dot_product_with => Dot_product_with_

     !**
     ! Miscellaneous
     !*
     procedure, public :: Interpolate => Interpolate_
     procedure, public :: Compare_with => Compare_with_
     procedure, public :: Is_zero => Is_zero_
     
     procedure, public :: cvector_to_mr => cvector_to_mr_
     procedure, public :: cvector_to_mr_b => cvector_to_mr_b_

     procedure, public :: mr_to_cvector => mr_to_cvector_

  end type Vector3d_cmr_complex_t

  ! Constructors for TCVector3D_CSG
  interface Vector3d_cmr_complex_t
     module procedure Vector3d_cmr_complex_t_ctor_default
     module procedure Vector3d_cmr_complex_t_ctor_copy
     module procedure Vector3d_cmr_complex_t_ctor1
  end interface Vector3d_cmr_complex_t

contains

  !
  !************************************************
  ! Initialization and finalization
  !************************************************
  !

  !**
  ! Vector3d_cmr_complex_t default constructor.
  !
  !*
  function Vector3d_cmr_complex_t_ctor_default () result (E)
    implicit none
    ! Local variables
    type (Vector3d_cmr_complex_t) :: E
    !
    !***********************
    ! Executable statements
    !***********************
    !
    E%grid_ptr => null ()
    E%is_allocated = .false.
  end function Vector3d_cmr_complex_t_ctor_default

  !**
  ! 
  ! Copy constructor
  !
  ! Arguments
  !
  !*
  function Vector3d_cmr_complex_t_ctor_copy (E_in) result (E)
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent (in) :: E_in
    ! Local variables
    type (Vector3d_cmr_complex_t) :: E
    integer :: i
    !
    !***********************
    ! Executable statements
    !***********************
    !
    call E%Initialize (E_in%grid_ptr, E_in%grid_type)
    
    E%ind_active = E_in%ind_active
    E%ind_interior = E_in%ind_interior
    E%ind_boundaries = E_in%ind_boundaries
    
    do i = 1, E_in%grid_ptr%ngrids
       E%sub_vectors(i) = E_in%sub_vectors(i)
    end do
 
  end function Vector3d_cmr_complex_t_ctor_copy

  !**
  ! Vector3d_cmr_complex_t parameterized constructor
  !
  ! Arguments
  !   igrid     Underlying grid.
  !   grid_type  Definied in GridDef.f90
  !
  !*
  function Vector3d_cmr_complex_t_ctor1 (igrid, grid_type) result (E)
    implicit none
    ! Arguments
    class (Grid3d_cmr_t), target, intent (in) :: igrid
    character (*)              , intent (in) :: grid_type
    ! Local variables
    type (Vector3d_cmr_complex_t) :: E
    !
    !***********************
    ! Executable statements
    !***********************
    !
    call E%Initialize (igrid, grid_type)

  end function Vector3d_cmr_complex_t_ctor1

  !**
  ! INITIALIZE
  !
  !*
  subroutine Vector3d_cmr_complex_t_initialize (self, igrid, grid_type) 
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent (in out) :: self
    class (Grid3d_cmr_t), pointer, intent (in)     :: igrid
    character (*)                , intent (in)     :: grid_type
    ! Local variables
    integer :: status
    integer :: i
    !
    !***********************
    ! Executable statements
    !***********************
    !    
    self%grid_ptr => igrid
    self%grid_type = grid_type

    self%is_allocated = .TRUE.
    allocate (self%sub_vectors(self%grid_ptr%ngrids), STAT = status)
    self%is_allocated = self%is_allocated .AND. (status.EQ.0)

    do i = 1, self%grid_ptr%ngrids
       self%sub_vectors(i) = Vector3d_csg_complex_t (&
            igrid%sub_grids(i), grid_type)
    end do

    call self%set_active_int_boundary ()
    call self%Zero ()

  end subroutine Vector3d_cmr_complex_t_initialize

  !**
  ! Vector3d_cmr_complex_t destructor
  !
  !*
  subroutine Vector3d_cmr_complex_t_dtor (self)
    implicit none
    ! Arguments
    type (Vector3d_cmr_complex_t), intent (in out) :: self
    !
    character(80) :: ctmp

    ctmp = self%grid_type
  end subroutine Vector3d_cmr_complex_t_dtor

  !
  !************************************************
  ! Input/Output
  !************************************************
  !
  !**
  ! WRITE
  !
  !*
  subroutine Write_ (self, fid, ftype)
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent (in) :: self
    integer                      , intent (in) :: fid
    character (*), optional      , intent (in) :: ftype
     integer :: itmp
    character(80) :: ctmp

    ctmp = self%grid_type
    ctmp = ftype    
    itmp = fid

    STOP 'vector3d_cmr_complex.f90:Write: Not implemented.'
  end subroutine Write_

  !**
  ! READ
  !
  !*
  subroutine Read_ (self, fid, ftype)
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent (in out) :: self
    integer                      , intent (in)     :: fid
    character (*)                , intent (in), optional :: ftype
    ! Local variables
    integer :: itmp
    character(80) :: ctmp

    ctmp = self%grid_type
    ctmp = ftype    
    itmp = fid

    STOP 'vector3d_cmr_complex.f90:Read: Not implemented.'
  end subroutine Read_

  !
  !************************************************
  ! Boundary operations
  !************************************************
  !

  !**
  ! SET_ALL_BOUNDARY
  !
  !*
  subroutine Set_all_boundary_ (self, c_in)
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent (in out)       :: self
    complex (kind=prec)        , intent (in), optional :: c_in
    ! Local variables
    complex(kind=prec) :: rtmp
    character(80) :: ctmp
    !
    !***********************
    ! Executable statements
    !***********************
    !
    ctmp = self%grid_type
    rtmp = c_in
    
    STOP 'vector3d_cmr_complex.f90:Set_all_boundary: Not implemented.'
  end subroutine Set_all_boundary_

  !**
  ! SET_ONE_BOUNDARY
  !
  !*
  subroutine Set_one_boundary_ (self, bdry, c, int_only)
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent (in out) :: self
    character (*)                , intent (in)     :: bdry
    complex (kind=prec)        , intent (in)     :: c
    logical                      , intent (in), optional :: int_only
    ! Local variables
    character(80) ::ctmp
    complex (kind=prec) :: rtmp

    if (int_only) then
       ctmp = self%grid_type
    end if
    
    ctmp = bdry
    rtmp = c
    
    STOP 'vector3d_cmr_complex.f90:Set_one_boundary: Not implemented.'
    
  end subroutine Set_one_boundary_

  !**
  ! SET_ACTIVE_INT_BOUNDARY
  !*
  subroutine Set_active_int_boundary_ (self, xy_in) 
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t)  :: self
    logical, intent (in), optional :: xy_in
    ! Local variables
    logical :: xy, int_only
    integer :: i, k
    integer :: n_full, n_active, n_interior, n_boundaries
    complex (kind=prec), dimension (:), allocatable :: v_1, v_2
    !
    !***********************
    ! Executable statements
    !***********************
    !
    if (.not.present (xy_in)) then
       xy = .false.
    else
       xy = xy_in
    end if

    ! Loop over subgrids, setting boundary edges to one,
    ! interior to  zero
    do k = 1, self%grid_ptr%ngrids
       call self%sub_vectors(k)%set_all_boundary (cmplx(1._prec, 0._prec, prec))
    end do

    ! Loop over interfaces: set redundant interface edges to 2
    select case (self%grid_type)
    case ('EDGE')
       int_only = .TRUE.

    case ('FACE')
       int_only = .false.

    case ('NODE')
       int_only = .TRUE.

    case default
       write(*, *) 'ERROR:Vector3d_cmr_complex_t:set_active_int_boundary: '
       write(*, *) '      Invalid grid type option!'

       STOP
    end select

    do k = 2, self%grid_ptr%ngrids
       if (self%grid_ptr%Coarseness(k - 1, 1) < self%grid_ptr%Coarseness(k, 1)) then
          ! upper grid is finer: grid k-1 interface nodes are
          ! not active; also reset interior part of interface
          ! edges to 0
          if (xy) then
             call self%sub_vectors(k-1)%set_one_boundary('z2_x', &
                  cmplx(-1._prec, 0._prec, prec))
             call self%sub_vectors(k-1)%set_one_boundary('z2_y', &
                  cmplx(-10._prec, 0._prec, prec))
          else
             call self%sub_vectors(k-1)%set_one_boundary('z2', &
                  cmplx(-1._prec, 0._prec, prec))
          end if
          call self%sub_vectors(k)%set_one_boundary('z1', &
               cmplx(0._prec, 0._prec, prec), int_only)
       else
          if (xy) then
             call self%sub_vectors(k)%set_one_boundary('z1_x', &
                  cmplx(-1.0_prec, 0._prec, prec))
             call self%sub_vectors(k)%set_one_boundary('z1_y', &
                  cmplx(-10.0_prec, 0._prec, prec))
          else                
             call self%sub_vectors(k)%set_one_boundary('z1', &
                  cmplx(-1.0_prec, 0._prec, prec))
          end if
          call self%sub_vectors(k-1)%set_one_boundary('z2', &
               cmplx(0._prec, 0._prec, prec), int_only)
       end if
    end do

    !******                                        ******
    ! ***  Set active, interior, and boundary edges. ***
    !******                                        ******

    call self%Get_full(v_1)

    n_full = size (v_1)

    ! Start - Set indices of active edges
    !
    n_active = 0
    do k = 1, n_full
       if (real (v_1(k)) >= 0) then
          n_active = n_active + 1
       end if
    end do

    allocate (self%ind_active(n_active))

    i = 0
    do k = 1, n_full
       if (real (v_1(k)) >= 0) then
          i = i + 1
          self%ind_active(i) = k
       end if
    end do
    !
    ! End - Set indices of active edges

    ! Start - Set indices of interior edges
    !
    n_interior = 0
    do k = 1, n_full
       if (real (v_1(k)) == 0) then
          n_interior = n_interior + 1
       end if
    end do

    allocate (v_2(n_active))
    v_2 = v_1(self%ind_active)

    allocate (self%ind_interior(n_interior))

    i = 0
    do k = 1, n_active
       if (real (v_2(k)) == 0) then
          i = i + 1
          self%ind_interior(i) = k               
       end if
    end do
    !
    ! End - Set indices of interior edges

    ! Start - Set indices of boundary edges
    !
    n_boundaries = 0
    do k = 1, n_active
       if (real (v_2(k)) == 1) then
          n_boundaries = n_boundaries + 1
       end if
    end do

    allocate (self%ind_boundaries(n_boundaries)) 

    i = 0
    do k = 1, n_active
       if (real(v_2(k)) == 1) then
          i = i + 1
          self%ind_boundaries(i) = k             
       end if
    end do
    !
    ! End - Set indices of boundary edges

  end subroutine Set_active_int_boundary_

  subroutine Get_int_bdry_indices_ (self, ind_i, ind_b)
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent (in)  :: self
    integer, allocatable         , intent (out) :: ind_i(:), ind_b(:)
    ! Local variables
    integer :: m, n
    !
    !***********************
    ! Executable statements
    !***********************
    !
    m = size (self%ind_interior)
    n = size (self%ind_boundaries)

    allocate (ind_i(m))
    allocate (ind_b(n))

    ind_i = self%ind_interior
    ind_b = self%ind_boundaries
  end subroutine Get_int_bdry_indices_

  !
  !************************************************
  ! Data access
  !************************************************

  !**
  ! length
  !
  !*
  function length_ (self) result (n)
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent (in) :: self
    ! Local variables
    integer :: n
    !
    !***********************
    ! Executable statements
    !***********************
    !
    n = size (self%ind_active)

  end function length_

  !**
  ! GET_ARRAY
  !*
  subroutine Get_array_ (self, v)
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t)     , intent (in)  :: self
    complex (kind=prec), allocatable, intent (out) :: v(:)
    ! Local variables
    complex (kind=prec), allocatable :: v_full(:)
    !
    !***********************
    ! Executable statements
    !***********************
    !
    allocate (v(size(self%ind_active)))
    call self%Get_full (v_full)
    v = v_full(self%ind_active)

  end subroutine Get_array_

  !**
  ! SET_ARRAY
  !
  !*
  subroutine Set_array_ (self, v)
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent (in out) :: self
    complex (kind=prec)         , intent (in)     :: v(:)
    ! Local variables
    complex (kind=prec), allocatable :: vFull(:)
    !
    !***********************
    ! Executable statements
    !***********************
    !
    allocate (vFull(self%length_full ()))
    vFull = 0.0
    vFull(self%ind_active) = v
    call self%set_full (vFull)
  end subroutine Set_array_

  !**
  ! COPY_FROM
  !
  !*
  subroutine Copy_from_ (self, rhs)
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent (in out) :: self
    class (Vector3d_cmr_complex_t), intent (in)     :: rhs
    ! Local variables
    integer :: i, nvecs
    !
    !***********************
    ! Executable statements
    !***********************
    !
    if (.not.rhs%is_allocated) then
       write(*, *) 'ERROR:Vector3d_cmr_complex_t::copy_from: '
       write(*, *) '      Input vector not allocated.'

       STOP
    end if

    select type (rhs)
    class is (Vector3d_cmr_complex_t)  
       if (.not.self%is_allocated) then
          call self%Initialize(rhs%grid_ptr, rhs%grid_type)
       else
          if (self%grid_type /= rhs%grid_type) then
             write(*, *) 'ERROR:Vector3d_cmr_complex_t::copy_from: '
             write(*, *) '      Grid types not compatible.'

             STOP       
          end if
       end if

       nvecs = size (self%sub_vectors)
       do i = 1, nvecs
          call self%sub_vectors(i)%copy_from (rhs%sub_vectors(i))
       end do

       class default
       write(*, *) 'ERROR:Vector3d_cmr_complex_t::copy_from: '
       write(*, *) '      Runtime incompatible class for output vector "rhs".'

       STOP       
    end select

  end subroutine Copy_from_

  !**
  ! GET_FULL Creates standard vector (1-D array) of E fields
  ! for all sub-vectors,  INCLUDING redundant interface edges
  !
  !*
  subroutine Get_full_ (self, v)
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent (in) :: self
    ! Local variables
    complex (kind=prec), allocatable, intent (out) :: v(:)
    integer :: n, i1, i2, k
    complex (kind=prec), allocatable :: v_temp(:)
    !
    !***********************
    ! Executable statements
    !***********************
    !
    n = self%length_full ()
    allocate (v(n))

    v = C_ZERO

    i1 = 1
    i2 = 0

    do k = 1, self%grid_ptr%ngrids
       n = self%sub_vectors(k)%length()         
       i2 = i2 + n
       call self%sub_vectors(k)%Get_array (v_temp)
       v(i1:i2) = v_temp
       i1 = i1 + n
    end do

  end subroutine Get_full_

  subroutine Set_full_ (self, v)
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent (in out) :: self
    complex (kind=prec)        , intent (in)     :: v(:)
    ! Local variables
    integer :: i1, i2, k, n
    !
    !***********************
    ! Executable statements
    !***********************
    !
    i1 = 1; i2 = 0;
    do k = 1, self%grid_ptr%ngrids
       n = self%sub_vectors(k)%length ()
       i2 = i2 + n
       call self%sub_vectors(k)%set_array (v(i1:i2))
       i1 = i1 + n
    end do
  end subroutine Set_full_

  !**
  ! length_FULL Find total number of edges in the
  ! multiresolution object, including redundant edges.
  !
  !*
  function length_full_ (self) result (n)
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t) :: self
    ! Local variables
    integer :: k
    integer :: n
    !
    !***********************
    ! Executable statements
    !***********************
    !
    n = 0
    do k = 1, self%grid_ptr%ngrids
       n = n + self%sub_vectors(k)%length()
    end do

  end function length_full_

  ! !**
  ! ! FIND_FULL
  ! !*
  ! function Find_full_ (self, c) result (I)
  !   implicit none
  !   ! Arguments
  !   class (Vector3d_cmr_complex_t), intent (in) :: self
  !   complex (kind=prec)        , intent (in) :: c
  !   ! Local variables
  !   integer, allocatable :: I(:)
  !   complex (kind=prec), allocatable :: v(:)
  !   integer :: n, n_I, k
  !   !
  !   !***********************
  !   ! Executable statements
  !   !***********************
  !   !
  !   n = self%length_full()
  !   call self%Get_full (v)

  !   n_I = 0
  !   do k = 1, n
  !      if (v(k) == c) n_I = n_I + 1
  !   end do

  !   allocate (I(n_I))

  !   n_I = 0
  !   do k = 1, n
  !      if (v(k) == c) then
  !         n_I = n_I + 1
  !         I(n_I) = k
  !      end if
  !   end do

  !   !**
  !   ! Clean up

  !   !*
  ! end function Find_full_

  !**
  ! FIND
  !*
  ! subroutine Find_ (self, I, c) 
  !   implicit none 
  !   ! Arguments
  !   class (Vector3d_cmr_complex_t), intent (in)  :: self
  !   integer, allocatable         , intent (out) :: I(:)
  !   complex (kind=prec)        , intent (in)  :: c
  !   ! Local variables
  !   complex (kind=prec), allocatable :: v(:)
  !   integer :: n, n_I, k
  !   real (kind=prec), parameter :: TOL = 1E-5
  !   !
  !   !***********************
  !   ! Executable statements
  !   !***********************
  !   !
  !   n = self%length()
  !   allocate (v(n))
  !   call self%Get_array (v)

  !   n_I = 0
  !   do k = 1, n
  !      if (abs (v(k) - c)/abs (c) <= TOL) n_I = n_I + 1
  !   end do

  !   allocate (I(n_I))

  !   n_I = 0
  !   do k = 1, n
  !      if (abs (v(k) - c)/abs (c) <= TOL) then
  !         n_I = n_I + 1
  !         I(n_I) = k
  !      end if
  !   end do

  ! end subroutine Find_

  !
  !************************************************
  ! Arithmetic operations
  !************************************************
  !

  !**
  ! ZERO zeros variable of derived data type
  ! cvector;
  !*
  subroutine Zero_ (self)
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent (in out) :: self
    ! Local variables
    integer :: i
    !
    !***********************
    ! Executable statements
    !***********************
    !
    do i = 1, self%grid_ptr%NGrids
       call self%sub_vectors(i)%Zero ()
    end do

  end subroutine Zero_

  !**
  ! ADD
  !*
  subroutine Add_ (self, rhs) 
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent (in out) :: self
    class (Vector3d_cmr_complex_t), intent (in)     :: rhs
    ! Local variables
    character(80) :: ctmp
    !
    !***********************
    ! Executable statements
    !***********************
    !
    ctmp = self%grid_type
    ctmp = rhs%grid_type
    
    STOP 'Vector3d_cmr_complex::Add: Not implemented.'
  end subroutine Add_

  !**
  ! SUBTRACT 
  !*
  subroutine Subtract_ (self, rhs) 
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent (in out) :: self
    class (Vector3d_cmr_complex_t), intent (in)     :: rhs
    ! Local variables
    character(80) :: ctmp
    !
    !***********************
    ! Executable statements
    !***********************
    !
    ctmp = self%grid_type
    ctmp = rhs%grid_type
    
    STOP 'Vector3d_cmr_complex::Subtract: Not implemented.'

  end subroutine Subtract_

  !**
  ! MULTIPLY_BY
  !
  ! Multiplies two vectors E1, E2 stored as 
  ! derived data type cvector pointwise; 
  ! subroutine version
  ! E3 can overwrite E1 or E2
  !*
  subroutine Multiply_by_ (self, rhs) 
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent (in out) :: self
    class (Vector3d_cmr_complex_t), intent (in)     :: rhs
    ! Local variables
    character(80) :: ctmp
    !
    !***********************
    ! Executable statements
    !***********************
    !
    ctmp = self%grid_type
    ctmp = rhs%grid_type
    
    STOP 'Vector3d_cmr_complex::Multiply_by: Not implemented.'
  end subroutine Multiply_by_

  !**
  ! MULTIPLY_BY_SCALAR
  !
  ! Multiplies vector stored as derived data 
  ! type cvector with a real scalar; subroutine version
  ! E2 can overwrite E1
  !*
  subroutine Multiply_by_scalar_ (self, c) 
    implicit none
    ! Arguments
    !   a real scalar to be multiplied with
    class (Vector3d_cmr_complex_t), intent (in out) :: self
    complex (kind=prec)        , intent (in)     :: c    
    ! Local variables
    character(80) :: ctmp
    complex (kind=prec) :: rtmp
    !
    !***********************
    ! Executable statements
    !***********************
    !
    ctmp = self%grid_type
    rtmp = c
    
    STOP 'Vector3d_cmr_complex_t::Multiply_by_scalar: Not implemented.'
  end subroutine Multiply_by_scalar_

  subroutine Divide_by_ (self, rhs) 
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent (in out) :: self
    class (Vector3d_cmr_complex_t), intent (in)     :: rhs
    ! Local variables
    character(80) :: ctmp
    !
    !***********************
    ! Executable statements
    !***********************
    !
    ctmp = self%grid_type
    ctmp = rhs%grid_type
    
    STOP 'Vector3d_cmr_complex_t::divide_by: Not implemented.'

  end subroutine Divide_by_

  !
  !************************************************
  ! Algebraic operations
  !************************************************
  !

  !**
  ! DOT_PRODUCT_WITH
  !
  ! Computes dot product of two vectors stored
  ! as derived data type cvector, returning a real number
  !*
  function Dot_product_with_ (self, rhs) result (c)
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent (in) :: self
    class (Vector3d_cmr_complex_t), intent (in) :: rhs
    ! Local variables
    complex (kind=prec) :: c
    character(80) :: ctmp
    !
    !***********************
    ! Executable statements
    !***********************
    !
    ctmp = self%grid_type
    ctmp = rhs%grid_type
    c = 0
    
    STOP 'Vector3d_cmr_complex_t::divide_by: Not implemented.'

  end function Dot_product_with_

  !
  !************************************************
  ! Miscellaneous
  !************************************************
  !
  function Interpolate_ (self, component, rx) result (ret)
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent (in) :: self
    integer                       , intent (in) :: component
    real (kind=prec)            , intent (in) :: rx(3)
    ! Local variables
    complex (kind=prec) :: ret
    character(80) :: ctmp
    real (kind=prec) :: rtmp
    integer :: itmp
    !
    !***********************
    ! Executable statements
    !***********************
    !
    itmp = component
    ctmp = self%grid_type
    rtmp = rx(1)
    ret = 0
    
    STOP 'vector3d_cmr_complex.f90:Interpolate: Not implemented.' 
  end function Interpolate_

  !**
  ! COMPARE_WITH
  !
  ! Check two vectors for compatibility for linear operator
  ! purposes.
  !*
  function Compare_with_ (self, rhs) result (status)
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent (in) :: self
    class (Vector3d_cmr_complex_t), intent (in) :: rhs
    ! Local variables
    logical :: status
    character(80) :: ctmp
    !
    !***********************
    ! Executable statements
    !***********************
    !
    status = .FALSE.

    ctmp = self%grid_type
    ctmp = rhs%grid_type

    STOP 'vector3d_cmr_complex.f90:Compare_with: Not implemented.' 
  end function Compare_with_

  !**
  ! IS_ZERO
  !
  !*
  function Is_zero_ (self) result (ret)
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent (in) :: self
    ! Local variables
    logical :: ret
    !
    !***********************
    ! Executable statements
    !***********************
    !
    ret = self%is_allocated

  end function Is_zero_

  subroutine operator_assignment (lhs, rhs)
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent (in out) :: lhs
    class (Vector3d_cmr_complex_t), intent (in)     :: rhs
    ! Local variables
    
    !
    !***********************
    ! Executable statements
    !***********************
    !
    call lhs%Copy_from (rhs)
    
  end subroutine operator_assignment



  !**
  ! MR_TO_CVECTOR
  !
  ! Convert an MR TVector to a full SG, filling in the full
  ! fine grid.
  ! Converts MR Vector object to SG
  ! copying from variable resolution subgrids to completely fill in the
  ! underlying fine grid.
  !*
  subroutine mr_to_cvector_(self, sg_v, sg_grid)
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent(in)    :: self
    type (cvector)                , intent(inout) :: sg_v
    type (grid_t), pointer        , intent (in)   :: sg_grid
    ! Local variables
    type (cvector) :: temp
    integer :: x_nx, x_ny, x_nz
    integer :: y_nx, y_ny, y_nz
    integer :: z_nx, z_ny, z_nz
    integer :: last, Cs, i1, i2, i, k
    real (kind=prec) w1, w2
    !
    !***********************
    ! Executable statements
    !***********************
    !
    call create_cvector(sg_grid, temp, self%grid_type)

    x_nx = size(temp%x, 1)
    x_ny = size(temp%x, 2)
    x_nz = size(temp%x, 3)
    
    y_nx = size(temp%y, 1)
    y_ny = size(temp%y, 2)
    y_nz = size(temp%y, 3)

    z_nx = size(temp%z, 1)
    z_ny = size(temp%z, 2)
    z_nz = size(temp%z, 3)

    sg_v%x = 0; sg_v%y = 0; sg_v%z = 0
    
    select case (self%grid_type)
    case (EDGE)
       do k = 1, self%grid_ptr%NGrids
          Cs = 2**self%grid_ptr%Coarseness(k, 1)
          i1 = self%grid_ptr%Coarseness(k, 3)
          i2 = self%grid_ptr%Coarseness(k, 4)
          ! Copy  x and y components in x and y directions
          ! edges that aligned with subgrid edge.
          do i = 1, Cs 
             temp%x(i:x_nx:Cs, 1:x_ny:Cs, i1:i2+1) = self%sub_vectors(k)%x
             temp%y(1:y_nx:Cs, i:y_ny:Cs, i1:i2+1) = self%sub_vectors(k)%y
             
             w1 = 1. - (i - 1.)/Cs
             w2 = 1. - w1
             
             if (i == 1) then
                temp%z(1:z_nx:Cs, 1:z_ny:Cs, i1:i2) = self%sub_vectors(k)%z
             else
                last = size(self%sub_vectors(k)%z, 1)
                temp%z(i:z_nx:Cs, 1:z_ny:Cs, i1:i2) = &
                     self%sub_vectors(k)%z(1:last-1, :, :) * w1 +&
                     self%sub_vectors(k)%z(2:last, :, :) * w2
             end if
          end do
          ! edges that subdivide the subgrid
          ! interpolate  in y and x direxctions
          ! copy/interpolate x in y direction
          ! copy x and y along x and y directions
          ! respectively
          do i = 2, Cs
             w1 = 1. - (i - 1.)/Cs
             w2 = 1. - w1

             temp%x(:, i:x_ny:Cs, i1:i2+1) = temp%x(:, 1:x_ny-Cs:Cs, i1:i2+1)*w1 + &
                  temp%x(:, Cs+1:x_ny:Cs, i1:i2+1)*w2
             
             temp%y(i:y_nx:Cs, :, i1:i2+1) = temp%y(1:y_nx-Cs:Cs, :, i1:i2+1)*w1 + &
                  temp%y(Cs+1:y_nx:Cs, :, i1:i2+1)*w2
             
             ! added by zhhq, 2017
             temp%z(:, i:z_ny:Cs, i1:i2) = temp%z(:, 1:z_ny-Cs:Cs, i1:i2)*w1 + &
                  temp%z(:, Cs+1:z_ny:Cs, i1:i2) * w2
             ! temp.z(i:Cs:end,i:Cs:end,i1:i2) = temp.z(:,1:Cs:end-Cs,i1:i2)*w1+ ...
             ! temp.z(:,Cs+1:Cs:end,i1:i2)*w2;
             ! added by zhhq, 2017
          end do
          
          sg_v%x(:, :, i1:i2+1) = sg_v%x(:, :, i1:i2+1) + temp%x(:, :, i1:i2+1)
          sg_v%y(:, :, i1:i2+1) = sg_v%y(:, :, i1:i2+1) + temp%y(:, :, i1:i2+1)
          sg_v%z(:, :, i1:i2)   = temp%z(:, :, i1:i2)
          
       end do

    case default
       STOP 'ERROR:mr_to_cvector: Unrecognized grid type.'
    end select
    
  end subroutine mr_to_cvector_

  subroutine cvector_to_mr_ (self, sg_v)
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent(inout) :: self
    type (cvector)                , intent(in)    :: sg_v
    ! Local variables
    type (grid_t), pointer :: grid
    
    complex (kind=prec), pointer, dimension(:) :: tempe
    type(rvector) :: templ_r
    real (kind=prec), pointer, dimension(:) :: templ
    complex (kind=prec), pointer, dimension(:) :: temple
    complex (kind=prec), allocatable, dimension(:,:,:) :: lengthx, lengthy
    integer :: sx1, sx2, sx3, sy1, sy2, sy3, s1, s2
    integer :: Cs, i1, i2

    type (cvector) :: tempEL
    type (rvector) :: temp_L
        
    integer :: k, i
    !
    !***********************
    ! Executable statements
    !***********************
    !
    
    grid => sg_v%grid
    
    call edgeLength (grid, templ_r)
    templ => null ()
    call getRVector (templ_r, templ)
    
    tempe => null ()
    call getCVector(sg_v, tempe)

    allocate(temple(size(templ)))
    temple = templ*tempe
    
    call create_rvector (grid, temp_L, EDGE)
    call setRVector(templ, temp_L)

    call create_cvector (grid, tempEL, EDGE)
    call setCVector(temple, tempEL)

    do k = 1, self%grid_ptr%NGrids
       Cs = 2**self%grid_ptr%Coarseness(k, 1)
       i1 = self%grid_ptr%Coarseness(k, 3)
       i2 = self%grid_ptr%Coarseness(k, 4)

       sx1 = size(self%sub_vectors(k)%x, 1)
       sx2 = size(self%sub_vectors(k)%x, 2)
       sx3 = size(self%sub_vectors(k)%x, 3)
       allocate (lengthx(sx1, sx2, sx3))
       lengthx = 0.0
       
       sy1 = size(self%sub_vectors(k)%y, 1)
       sy2 = size(self%sub_vectors(k)%y, 2)
       sy3 = size(self%sub_vectors(k)%y, 3)
       allocate (lengthy(sy1, sy2, sy3))
       lengthy = 0.0

       do i = 1, Cs
          s1 = size(tempEL%x, 1)
          s2 = size(tempEL%x, 2)
          self%sub_vectors(k)%x = self%sub_vectors(k)%x + &
               tempEL%x(i:s1:Cs, 1:s2:Cs, i1:i2+1)

          s1 = size(tempEL%y, 1)
          s2 = size(tempEL%y, 2)
          self%sub_vectors(k)%y = self%sub_vectors(k)%y + &
               tempEL%y(1:s1:Cs,i:s2:Cs, i1:i2+1)

          s1 = size(temp_L%x, 1)
          s2 = size(temp_L%x, 2)          
          lengthx = lengthx + temp_L%x(i:s1:Cs, 1:s2:Cs, i1:i2+1);
          s1 = size(temp_L%y, 1)
          s2 = size(temp_L%y, 2)
          lengthy = lengthy + temp_L%y(1:s1:Cs, i:s2:Cs, i1:i2+1);
       end do
       
       self%sub_vectors(k)%x = self%sub_vectors(k)%x/lengthx
       self%sub_vectors(k)%y = self%sub_vectors(k)%y/lengthy
       s1 = size(sg_v%z, 1)
       s2 = size(sg_v%z, 2)
       self%sub_vectors(k)%z = sg_v%z(1:s1:Cs, 1:s2:Cs, i1:i2)

       deallocate (lengthx, lengthy)
    end do    
    
  end subroutine cvector_to_mr_

  !**
  ! CVECTOR_TO_MR
  !
  ! Converts SG Vector object to MR
  ! by averaging
  ! used only for e0
  !*
  subroutine cvector_to_mr_b_(self, sg_v)
    implicit none
    ! Arguments
    class (Vector3d_cmr_complex_t), intent(inout) :: self
    type (cvector)                , intent(in)    :: sg_v
    ! Local varaibles
    type (grid_t), pointer :: grid

    integer :: x_nx, x_ny, x_nz
    integer :: y_nx, y_ny, y_nz
    integer :: z_nx, z_ny, z_nz
    integer :: last, Cs, i1, i2, i, k

    !
    !***********************
    ! Executable statements
    !***********************
    !
    grid => sg_v%grid

    x_nx = size(sg_v%x, 1)
    x_ny = size(sg_v%x, 2)
    x_nz = size(sg_v%x, 3)
    
    y_nx = size(sg_v%y, 1)
    y_ny = size(sg_v%y, 2)
    y_nz = size(sg_v%y, 3)

    z_nx = size(sg_v%z, 1)
    z_ny = size(sg_v%z, 2)
    z_nz = size(sg_v%z, 3)
    
    do k = 1, self%grid_ptr%NGrids
       Cs = 2**self%grid_ptr%Coarseness(k, 1)
       i1 = self%grid_ptr%Coarseness(k, 3)
       i2 = self%grid_ptr%Coarseness(k, 4)
       
       do i = 1, Cs
          last = size(grid%Dx)
          self%sub_vectors(k)%x = self%sub_vectors(k)%x + &
               sg_v%x(i:x_nx:Cs, 1:x_ny:Cs, i1:i2+1) *    &
               rep_mat(grid%Dx(i:last:Cs), &
               1, &
               self%grid_ptr%sub_grids(k)%Ny + 1, &
               self%grid_ptr%sub_grids(k)%Nz + 1, .false.)
          
          last = size(grid%Dy)
          self%sub_vectors(k)%y = self%sub_vectors(k)%y + &
               sg_v%y(1:y_nx:Cs, i:y_ny:Cs, i1:i2+1) *  & 
               rep_mat(grid%Dy(i:last:Cs), &
               self%grid_ptr%sub_grids(k)%Nx + 1, &
               1, &
               self%grid_ptr%sub_grids(k)%Nz + 1, .TRUE.)
       end do
       
       self%sub_vectors(k)%x = self%sub_vectors(k)%x / &
            rep_mat(self%grid_ptr%sub_grids(k)%Dx, &
            1, &
            self%grid_ptr%sub_grids(k)%Ny + 1, &
            self%grid_ptr%sub_grids(k)%Nz + 1, .false.)
       
       self%sub_vectors(k)%y = self%sub_vectors(k)%y / &
            rep_mat(self%grid_ptr%sub_grids(k)%Dy, &
            self%grid_ptr%sub_grids(k)%Nx + 1, &
            1, &
            self%grid_ptr%sub_grids(k)%Nz + 1, .TRUE.)
       
       self%sub_vectors(k)%z = sg_v%z(1:z_nx:Cs, 1:z_ny:Cs, i1:i2)
    end do
    
  end subroutine cvector_to_mr_b_
  
  function rep_mat(m_in, nx, ny, nz, transp) result(m_out)
    implicit none
    ! Arguments
    real (kind=prec), dimension(:), intent(in) :: m_in
    integer, intent(in) :: nx, ny, nz
    logical, intent(in) :: transp
    ! Local variables
    real (kind=prec), dimension(:,:,:), allocatable :: m_out
    integer :: i, i1, i2, n_in
    !
    !***********************
    ! Executable statements
    !***********************
    !
    n_in = size(m_in)

    if (transp) then
       allocate(m_out(nx, n_in*ny, nz))

       !*
       ! Copy along 2nd dimension.
       i1 = 1; i2 = n_in
       do i = 1, ny
          m_out(1, i1:i2, 1) = m_in
          i1 = i2 + 1
          i2 = i2 + n_in
       end do

       !*
       ! Copy along 1st dimension.
       do i = 1, nx
          m_out(i, :, 1) = m_out(1, :, 1)
       end do

       !*
       ! Copy along 3rd dimension
       do i = 1, nz
          m_out(:, :, i) = m_out(:, :, 1)
       end do
       
    else
       allocate(m_out(n_in*nx, ny, nz))

       !*
       ! Copy along 1st dimension.
       i1 = 1; i2 = n_in
       do i = 1, nx
          m_out(i1:i2, 1, 1) = m_in
          i1 = i2 + 1
          i2 = i2 + n_in
       end do

       !*
       ! Copy along 2nd dimension.
       do i = 1, ny
          m_out(:, i, 1) = m_out(:, 1, 1)
       end do

       !*
       ! Copy along 3rd dimension
       do i = 1, nz
          m_out(:, :, i) = m_out(:, :, 1)
       end do
       
    end if

  end function rep_mat
  
end module Vector3d_cmr_complex
