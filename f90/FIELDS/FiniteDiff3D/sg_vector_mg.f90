module sg_vector_mg ! MULTIGRID field vectors

  ! This module creates(allocates) cvectors/rvectors on staggered MULTIgrid
  ! Deallocates  and does some basic algebraic operations.
  ! For algebraic operations, new routines use the old ones (from sg_vector)
  ! and do loop on subgrids

  ! Derived cvector_mg and rvector_mg data types are defined in sg_vector
  ! module (after cvector, rvector respectively)

use math_constants        ! math/ physics constants
use utilities
use griddef
use sg_vector
implicit none

  ! Important - overloading the '=' assignment
  interface assignment (=)
    module procedure copy_cvector_mg
    module procedure c2mg
    module procedure mg2c
  end interface

  ! Generic interfaces are done through subroutines
  ! creates edge/ face nodes
  interface create
    module procedure create_rvector_mg
    module procedure create_cvector_mg
  end interface

  ! deallocates the edge/ face nodes
  interface deall
    module procedure deall_rvector_mg
    module procedure deall_cvector_mg
  end interface

  ! set the values to zero
  interface zero
    module procedure zero_cvector_mg
  end interface

  ! scalar value multiplies the edge/ face nodes
  interface scMult
    module procedure scMult_cvector_mg
    module procedure scMultReal_cvector_mg
  end interface

  INTERFACE scMultAdd
     module procedure scMultAdd_cvector_mg
  END INTERFACE

  interface linComb
    module procedure linComb_cvector_mg
  end interface

! adds the edge/ face nodes
interface add
  module procedure add_cvector_mg
  module procedure add_rvector_mg
end interface

! pointwise vector (two vector data types) multiplication of edge/ face
! nodes
! and pointwise real-complex (mixed) multiplication of edge/ face nodes
! Both are vector data types

! IS IT REALLY NEEDED?
interface diagMult
  module procedure diagMult_cvector_mg
  module procedure diagMult_rvector_mg
  module procedure diagMult_crvector_mg
  module procedure diagMult_rcvector_mg
end interface

! subtracts the edge/ face nodes
interface subtract
  module procedure subtract_rvector_mg
  module procedure subtract_cvector_mg
end interface

interface dotProd
  module procedure dotProd_cvector_mg_f
end interface

interface dotProd_noConj
  module procedure dotProd_noConj_cvector_mg_f
end interface

interface  diagDiv
  module procedure diagDiv_rcvector_mg
  module procedure diagDiv_crvector_mg
end interface

interface UpdateZ
  module procedure updateZ_cvector
  module procedure updateZ_rvector
end interface

public  :: create_rvector_mg, create_cvector_mg, &
           deall_rvector_mg, deall_cvector_mg, &
           copy_cvector_mg, c2mg, mg2c, &
           zero_cvector_mg, &
           subtract_rvector_mg, subtract_cvector_mg, &
           scMult_cvector_mg, &
           add_cvector_mg,add_rvector_mg, &
           diagMult_cvector_mg, diagMult_crvector_mg, &
           diagMult_rvector_mg, diagMult_rcvector_mg, &
           linComb_cvector_mg, &
           dotProd_cvector_mg_f,  dotProd_noConj_cvector_mg_f, &
           diagDiv_rcvector_mg, diagDiv_crvector_mg, &
           updateZ_cvector, updateZ_rvector

! *******************************************************************************
! this type defines cvectors on MULTIgrid

type :: cvector_mg
   ! number of subgrids
   integer  :: mgridSize
   ! cvectors for subgrids
   type(cvector), pointer :: cvArray(:)
   ! coarseness
   integer, allocatable  :: coarseness(:)
   ! allocated:  .true.  x, y, z arrays have been allocated
   logical  :: allocated = .false.
   ! store the intention of the use in a character string defined
   ! in GridDef as a parameter: EDGE or FACE are two possibilities
   character (len=80)   :: gridType
   ! pointer to parent grid
   type (grid_t), pointer  :: grid
end type cvector_mg

! ******************************************************************************
! this type defines multigrids rvector
type rvector_mg
  ! number of subgrids
  integer  :: mgridSize
  ! store the intention of the use in a character string defined
  ! in GridDef as a parameter: EDGE or FACE are two possibilities
  character (len=80)  :: gridType
  ! rvectors for subgrids
  type(rvector), pointer :: rvArray(:)
  ! coarseness
  integer, allocatable  :: coarseness(:)
  ! allocated:  .true.  x, y, z arrays have been allocated
  logical  :: allocated = .false.

end type rvector_mg


Contains

  ! *******************************************************************************
  ! allocates multigrids rvector
  subroutine create_rvector_mg (mgrid, e, gridType)

  implicit none
    type(grid_t), target, intent(in)  :: mgrid
    ! the grid for which an edge/ face node field is being initialized
    type (rvector_mg), intent(inout)  :: e
    character (len=80), intent(in)    :: gridType
    !local
    integer  :: status, imgrid


    ! First deallocate anything, that's allocated
    call deall_rvector_mg(e)

   ! allocate memory for rvector arrays
    allocate(e%rvArray(mgrid%mgridSize), STAT = status)
    ! allocate rvectors
    do imgrid = 1, mgrid%mgridSize
      call create_rvector(mgrid%gridArray(imgrid),e%rvArray(imgrid),gridType)
    enddo
    ! allocate memory for coarseness array
    allocate(e%coarseness(mgrid%mgridSize), STAT = status)
    ! number of subgrids
    e%mgridSize = mgrid%mgridSize
    ! gridType
    e%gridType = gridType
    ! e%allocated will be true if all allocations succeed
    e%allocated = .true.

    ! copy coarseness from multigrid
    e%coarseness = mgrid%coarseness

  end subroutine create_rvector_mg

  ! ******************************************************************************
  ! allocates cvector for multigrid
  subroutine create_cvector_mg (mgrid, e, gridType)

  implicit none
    type(grid_t), target, intent(in)  :: mgrid
    ! the grid for which an edge/ face node field is being initialized
    type (cvector_mg), intent(inout) :: e
    character (len=80), intent(in) :: gridType
    !local
    integer  :: status, imgrid

    ! First deallocate anything, that's allocated
    call deall_cvector_mg(e)


    ! Set pointer
    e%grid => mgrid
    ! allocate memory for cvector arrays
    allocate(e%cvArray(mgrid%mgridSize), STAT = status)
    ! allocate cvectors
    do imgrid = 1, mgrid%mgridSize
      call create_cvector(mgrid%gridArray(imgrid),e%cvArray(imgrid),gridType)
    enddo
    ! allocate memory for coarseness array
    allocate(e%coarseness(mgrid%mgridSize), STAT = status)
    ! number of subgrids
    e%mgridSize = mgrid%mgridSize
    ! gridType
    e%gridType = gridType
    ! e%allocated will be true if all allocations succeed
    e%allocated = .true.

    ! copy coarseness from multigrid
    e%coarseness = mgrid%coarseness

  end subroutine create_cvector_mg

! *****************************************************************************
  subroutine deall_rvector_mg(e)

  implicit none
    type (rvector_mg)  :: e
    integer  :: status, imgrid


      if (e%allocated) then
    ! deallocate cvectors of subgrids
      do imgrid = 1, e%mgridSize
        call deall_rvector(e%rvArray(imgrid))
      enddo
      deallocate(e%rvArray, STAT = status)
      end if

      e%mgridSize = 0
      E%gridType = ''
      E%allocated = .false.

end  subroutine deall_rvector_mg

! *****************************************************************************
  subroutine deall_cvector_mg(e)

  implicit none
    type (cvector_mg)  :: e
    integer  :: status, imgrid


      if (e%allocated) then
    ! deallocate cvectors of subgrids
      do imgrid = 1, e%mgridSize
        call deall_cvector(e%cvArray(imgrid))
      enddo
      deallocate(e%coarseness, STAT=status)
      deallocate(e%cvArray, STAT = status)
      end if

      if(associated(e%grid)) nullify(e%grid)


      e%mgridSize = 0
      E%gridType = ''
      E%allocated = .false.

end  subroutine deall_cvector_mg

  ! *****************************************************************************
  ! copy cvector_mg to cvector_mg
  ! 20.04.2012
  subroutine copy_cvector_mg(e2,e1)

  ! first argument is output
  implicit none
    type (cvector_mg), intent(in)  :: e1
    type (cvector_mg), intent(inout)  :: e2
    !local
    integer  :: imgrid

    if(.not.e1%allocated) then
      print *, 'RHS not allocated yet for copy_cvector_mg'
    endif

    if(.not.e2%allocated) then
      call create(e1%grid, e2, e1%gridType)
    endif

     if (e1%gridType == e2%gridType.and.e1%mgridSize == e2%mgridSize) then
        do imgrid = 1, e1%mgridSize
          if(e2%cvArray(imgrid)%nx ==  e1%cvArray(imgrid)%nx.and. &
             e2%cvArray(imgrid)%ny ==  e1%cvArray(imgrid)%ny.and. &
             e2%cvArray(imgrid)%nz ==  e1%cvArray(imgrid)%nz) then

             e2%cvArray(imgrid)%x = e1%cvArray(imgrid)%x
             e2%cvArray(imgrid)%y = e1%cvArray(imgrid)%y
             e2%cvArray(imgrid)%z = e1%cvArray(imgrid)%z
             e2%cvArray(imgrid)%gridType = e1%cvArray(imgrid)%gridType
             e2%cvArray(imgrid)%grid => e1%cvArray(imgrid)%grid
          else
            print *, 'e1 and e2 are not the same size; copy_cvector_mg'
          endif
        enddo
        e2%gridType = e1%gridType
        e2%grid => e1%grid
        e2%coarseness = e1%coarseness
      else
        print *, 'not compatible usage for copy_cvector'
      endif


end subroutine copy_cvector_mg
  ! *****************************************************************************
  ! averages or copies fields from cvector to cvector_mg

  ! 12.07.2012
  subroutine c2mg(e2,e1)

  implicit none

    type(cvector), intent(in)  :: e1       ! input cvector
    type(cvector_mg), intent(inout)  ::e2  ! output cvector

    ! local
    complex (kind=prec),allocatable, dimension(:,:,:)    :: tempE2
    integer  :: imgrid,ifine, ix,iy,iz,izv,ic
    integer  :: nx,ny,nz,nzCum,ccoeff_current
    integer  :: errAll

    if (.not.e1%allocated)then
      print *, 'Error c2mg; e1 (cvector) is not allocated'
    endif

    if (.not.e2%allocated)then
      print *, 'Error c2mg; e2 (cvector_mg) is not allocated'
    endif

    ! allocate temp array
    allocate(tempE2(e1%grid%nx+1,e1%grid%ny+1,e1%grid%nz+1), STAT= errAll)

    ! which gridType and check
    if (e1%gridType == EDGE .and. e1%gridType == e2%gridType) then
      nzCum = 0
      do imgrid = 1, e2%mgridSize  ! Global loop on sub-grids
        nx = e2%cvArray(imgrid)%nx
        ny = e2%cvArray(imgrid)%ny
        nz = e2%cvArray(imgrid)%nz
        ccoeff_current = 2**e2%coarseness(imgrid)
        ! re-count x component
        tempE2 = C_ZERO
        do iz = 1, nz+1
           izv = iz + nzCum
          do iy = 1, ny+1
            do ix =1, nx
              do ic = 1, ccoeff_current
              ! we may average the fields, but it does not work properly
              !  tempE2(ix,iy,iz) = tempE2(ix,iy,iz) + e1%x(ccoeff_current*(ix-1)+ic,ccoeff_current*(iy-1)+1,izv)*e1%grid%dx(ccoeff_current*(ix-1)+ic)
              ! we just copy them
                tempE2(ix,iy,iz) =  e1%x(ccoeff_current*(ix-1)+ic,ccoeff_current*(iy-1)+1,izv)
              enddo ! ic
              !  tempE2(ix,iy,iz) = tempE2(ix,iy,iz)/e2%grid%gridArray(imgrid)%dx(ix)
                e2%cvArray(imgrid)%x(ix,iy,iz) = tempE2(ix,iy,iz)
            enddo  !ix
          enddo !iy
        ! re-count y component
        tempE2 = C_ZERO
          do ix = 1, nx+1
            do iy = 1, ny
              do ic = 1, ccoeff_current
               ! tempE2(ix,iy,iz) = tempE2(ix,iy,iz) +  e1%y(ccoeff_current*(ix-1)+1,ccoeff_current*(iy-1)+ic,izv)*e1%grid%dy(ccoeff_current*(iy-1)+ic)
                tempE2(ix,iy,iz) = e1%y(ccoeff_current*(ix-1)+1,ccoeff_current*(iy-1)+ic,izv)
              enddo !ic
               ! tempE2(ix,iy,iz) = tempE2(ix,iy,iz)/e2%grid%gridArray(imgrid)%dy(iy)
                e2%cvArray(imgrid)%y(ix,iy,iz) = tempE2(ix,iy,iz)
            enddo !iy
          enddo !ix
        enddo ! iz
        ! re-count z component
        do iy =1, ny+1
          do ix= 1, nx+1
            do iz = 1, nz
               izv = iz + nzCum
              e2%cvArray(imgrid)%z(ix,iy,iz) = e1%z(ccoeff_current*(ix-1)+1,ccoeff_current*(ix-1)+1,izv)
           enddo !iz
          enddo !ix
        enddo !iy
        nzCum =  nzCum + nz
      enddo  ! Global loop over sub-grids
    else if(e1%gridType == FACE .and. e1%gridType == e2%gridType) then
      nzCum = 0
      do imgrid = 1, e2%mgridSize  ! Global loop on sub-grids
        nx = e2%cvArray(imgrid)%nx
        ny = e2%cvArray(imgrid)%ny
        nz = e2%cvArray(imgrid)%nz
        ccoeff_current = 2**e2%coarseness(imgrid)
        ! re-count x component
        tempE2 = C_ZERO
        do iz = 1, nz
           izv = iz + nzCum
          do iy = 1, ny
            do ix =1, nx+1
              do ic = 1, ccoeff_current
              ! we may average the fields, but it does not work properly
              !  tempE2(ix,iy,iz) = tempE2(ix,iy,iz) + e1%x(ccoeff_current*(ix-1)+ic,ccoeff_current*(iy-1)+1,izv)*e1%grid%dx(ccoeff_current*(ix-1)+ic)
              ! we just copy them
                tempE2(ix,iy,iz) =  e1%x(ccoeff_current*(ix-1)+ic,ccoeff_current*(iy-1)+1,izv)
              enddo ! ic
              !  tempE2(ix,iy,iz) = tempE2(ix,iy,iz)/e2%grid%gridArray(imgrid)%dx(ix)
                e2%cvArray(imgrid)%x(ix,iy,iz) = tempE2(ix,iy,iz)
            enddo  !ix
          enddo !iy
        ! re-count y component
        tempE2 = C_ZERO
          do ix = 1, nx
            do iy = 1, ny+1
              do ic = 1, ccoeff_current
               ! tempE2(ix,iy,iz) = tempE2(ix,iy,iz) +  e1%y(ccoeff_current*(ix-1)+1,ccoeff_current*(iy-1)+ic,izv)*e1%grid%dy(ccoeff_current*(iy-1)+ic)
                tempE2(ix,iy,iz) = e1%y(ccoeff_current*(ix-1)+1,ccoeff_current*(iy-1)+ic,izv)
              enddo !ic
               ! tempE2(ix,iy,iz) = tempE2(ix,iy,iz)/e2%grid%gridArray(imgrid)%dy(iy)
                e2%cvArray(imgrid)%y(ix,iy,iz) = tempE2(ix,iy,iz)
            enddo !iy
          enddo !ix
        enddo ! iz
        ! re-count z component
        do iy =1, ny
          do ix= 1, nx
            do iz = 1, nz+1
               izv = iz + nzCum
              e2%cvArray(imgrid)%z(ix,iy,iz) = e1%z(ccoeff_current*(ix-1)+1,ccoeff_current*(ix-1)+1,izv)
           enddo !iz
          enddo !ix
        enddo !iy
        nzCum =  nzCum + nz
      enddo  ! Global loop over sub-grids
    else
      print *, 'Error c2mg; cvector and cvector_mg not are the same gridType'
    endif

    ! deallocate temp array
    deallocate(tempE2, STAT=errAll)

  end subroutine c2mg

  ! *****************************************************************************
  ! convert cvector_mg to cvector

  subroutine mg2c(e2, e1)

  implicit none
    type(cvector_mg), intent(in)  :: e1  ! cvector_mg in
    type(cvector), intent(inout)  :: e2  ! cvector out

    ! local
    integer :: nzCum, ccoeff_current, nx,ny,nz
    integer  :: imgrid,ix,iy,iz,izv, ic, iyy

    if (.not.e1%allocated)then
      print *, 'Error mg2c; e1 (cvector_mg) is not allocated'
    endif

    if (.not.e2%allocated)then
      print *, 'Error mg2c; e2 (cvector) is not allocated'
    endif

    ! check gridType
    if (e1%gridType == e2%gridType) then

        nzCum = 0
        do imgrid = 1, e1%mgridSize  ! Global loop on sub-grids
           nx = e1%cvArray(imgrid)%nx
           ny = e1%cvArray(imgrid)%ny
           nz = e1%cvArray(imgrid)%nz
           ccoeff_current = 2**e1%coarseness(imgrid)
           do iz = 2, nz
              izv = iz + nzCum
              do iy = 1, ny
                 do ix =1, nx
                    do ic = 1, ccoeff_current
                       e2%x(ccoeff_current*(ix-1)+ic,iy,izv) = e1%cvArray(imgrid)%x(ix,iy,iz)
                       e2%y(ccoeff_current*(ix-1)+ic,iy,izv) = e1%cvArray(imgrid)%y(ix,iy,iz)
                       e2%z(ccoeff_current*(ix-1)+ic,iy,izv) = e1%cvArray(imgrid)%z(ix,iy,iz)
                       e2%x(ix,ccoeff_current*(iy-1)+ic,izv) = e1%cvArray(imgrid)%x(ix,iy,iz)
                       e2%y(ix,ccoeff_current*(iy-1)+ic,izv) = e1%cvArray(imgrid)%y(ix,iy,iz)
                       e2%z(ix,ccoeff_current*(iy-1)+ic,izv) = e1%cvArray(imgrid)%z(ix,iy,iz)
                     enddo
                  enddo
               enddo
            enddo
        nzCum = nzCum + nz
        enddo   ! Global loop on sub-grids
    else
      print *, 'Error mg2c; cvector and cvector_mg not are the same gridType'
    endif

  end subroutine mg2c

! *****************************************************************************
  subroutine zero_cvector_mg(e)

  implicit none
  type(cvector_mg), intent(inout)  :: e
  ! local
  integer  ::imgrid

    do imgrid = 1, e%mgridSize
      call zero_cvector(e%cvArray(imgrid))
    enddo

  end subroutine zero_cvector_mg

! *****************************************************************************

subroutine scMult_cvector_mg(c, e1, e2)

implicit none
! a complex scalar to be multiplied with
complex(kind=prec), intent(in)  :: c
type (cvector_mg), intent(in)   :: e1
type (cvector_mg), intent(inout)  :: e2
!local
integer  :: imgrid

! check whether e1 and e2 have the same number of subgrids
  if (e1%mgridSize == e2%mgridSize) then

    do imgrid = 1, e1%mgridSize
      call scMult_cvector(c, e1%cvArray(imgrid),e2%cvArray(imgrid))
    enddo

  else
    write (0, *) 'Error:scMult_cvector_mg: vectors not same subgrids size'
  endif

end subroutine scMult_cvector_mg
! *****************************************************************************
 subroutine scMultReal_cvector_mg(c, e1, e2)

! MULTIGRId case

! NOT DEBUGED

 implicit none
  real (kind=prec), intent(in)   :: c
  ! a real scalar to be multiplied with
  type (cvector_mg), intent(in)  :: e1
  type (cvector_mg), intent(inout)  :: e2

  integer  :: imgrid

    if(.not.e1%allocated) then
      print *, 'RHS not allocated yet for scMultReal_cvector_mg'
      return
    endif

    ! check to see if LHS (E2) is active (allocated)
    if(.not.e2%allocated) then
       print *,'LHS was not allocated yet for scMultReal_cvector_mg'
    else

      if ((e1%gridType == e2%gridType).and.(e1%mgridSize == e2%mgridSize)) then

        do imgrid = 1, e1%mgridSize
         ! Check whether all the vector nodes are of the same size
         if((e1%cvArray(imgrid)%nx == e2%cvArray(imgrid)%nx).and.(e1%cvArray(imgrid)%ny == e2%cvArray(imgrid)%ny) &
                                           .and.(e1%cvArray(imgrid)%nz == e2%cvArray(imgrid)%nz)) then
             ! complex scalar multiplication for x,y,z-components
             e2%cvArray(imgrid)%x = e1%cvArray(imgrid)%x * c
             e2%cvArray(imgrid)%y = e1%cvArray(imgrid)%y * c
             e2%cvArray(imgrid)%z = e1%cvArray(imgrid)%z * c

         else
           print *, 'Error:scMultReal_cvector_mg: vectors not same size'
         endif
        enddo

      else
        print *, 'not compatible usage for scMultReal_cvector'
      end if
    endif

  end subroutine scMultReal_cvector_mg ! scMultReal_cvector_mg

! *****************************************************************************
subroutine linComb_cvector_mg(inc1, e1, inc2, e2, e3)

implicit none
!   input vectors
type (cvector_mg), intent(in)  :: e1, e2
!  input complex scalars
complex (kind=prec), intent(in)  :: inc1, inc2
! lin comp cvector
type (cvector_mg), intent(inout)          :: e3
!local
integer  :: imgrid

! check whether e1 and e2 have the same number of subgrids
  if (e1%mgridSize == e2%mgridSize) then
    e3%mgridSize = e1%mgridSize
    do imgrid = 1, e1%mgridSize
      call linComb_cvector(inc1, e1%cvArray(imgrid), inc2, &
                        e2%cvArray(imgrid), e3%cvArray(imgrid))
    enddo
   else

    write (0, *) 'Error:linComb_cvector_mg: vectors not same subgrids size'
  endif

end subroutine linComb_cvector_mg

! *****************************************************************************************
subroutine add_cvector_mg(e1, e2, e3)

implicit none
type (cvector_mg), intent(in)  :: e1, e2
type (cvector_mg), intent(inout)  :: e3
!local
integer  :: imgrid

! check whether e1 and e2 have the same number of subgrids
  if ((e1%mgridSize == e2%mgridSize).and.(e1%mgridSize == e3%mgridSize).and. &
                                   (e2%mgridSize == e3%mgridSize)) then

    do imgrid = 1, e1%mgridSize
      call add_cvector(e1%cvArray(imgrid), e2%cvArray(imgrid), e3%cvArray(imgrid))
    enddo

  else
     write(0, *) 'Error:add_cvector_mg: vectors not same subgrids size'

  endif

end subroutine add_cvector_mg
! *****************************************************************************************
subroutine add_rvector_mg(e1, e2, e3)

implicit none
type (rvector_mg), intent(in)  :: e1, e2
type (rvector_mg), intent(inout)  :: e3
!local
integer  :: imgrid

! check whether e1 and e2 have the same number of subgrids
  if ((e1%mgridSize == e2%mgridSize).and.(e1%mgridSize == e3%mgridSize).and. &
                                   (e2%mgridSize == e3%mgridSize)) then

    do imgrid = 1, e1%mgridSize
      call add_rvector(e1%rvArray(imgrid), e2%rvArray(imgrid), e3%rvArray(imgrid))
    enddo

  else
     write(0, *) 'Error:add_rvector_mg: vectors not same subgrids size'

  endif

end subroutine add_rvector_mg
! **********************************************************************************************
subroutine diagMult_cvector_mg(e1, e2, e3)

implicit none
type (cvector_mg), intent(in)  :: e1, e2
type (cvector_mg), intent(inout)  :: e3
!local
integer  :: imgrid

! check whether e1 and e2 have the same number of subgrids
  if ((e1%mgridSize == e2%mgridSize).and.(e1%mgridSize == e3%mgridSize).and. &
                                   (e2%mgridSize == e3%mgridSize)) then

    do imgrid = 1, e1%mgridSize
      call diagMult_cvector(e1%cvArray(imgrid), e2%cvArray(imgrid), e3%cvArray(imgrid))
    enddo

  else
     write(0, *) 'Error::diagMult_cvector_mg vectors not same subgrids size'

  endif

end subroutine diagMult_cvector_mg
! *****************************************************************************************
subroutine diagMult_rvector_mg(e1, e2, e3)

implicit none
type (rvector_mg), intent(in)  :: e1, e2
type (rvector_mg), intent(inout)  :: e3
!local
integer  :: imgrid

! check whether e1 and e2 have the same number of subgrids
  if ((e1%mgridSize == e2%mgridSize).and.(e1%mgridSize == e3%mgridSize).and. &
                                   (e2%mgridSize == e3%mgridSize)) then

    do imgrid = 1, e1%mgridSize
      call diagMult_rvector(e1%rvArray(imgrid), e2%rvArray(imgrid), e3%rvArray(imgrid))
    enddo

  else
     write(0, *) 'Error::diagMult_rvector_mg vectors not same subgrids size'

  endif

end subroutine diagMult_rvector_mg
! *****************************************************************************************
subroutine diagMult_crvector_mg(e1, e2, e3)

implicit none
type (cvector_mg), intent(in)  :: e1
type (rvector_mg), intent(in)  :: e2
type (cvector_mg), intent(inout)  :: e3
!local
integer  :: imgrid

! check whether e1 and e2 have the same number of subgrids
  if ((e1%mgridSize == e2%mgridSize).and.(e1%mgridSize == e3%mgridSize).and. &
                                   (e2%mgridSize == e3%mgridSize)) then

    do imgrid = 1, e1%mgridSize
      call diagMult_crvector(e1%cvArray(imgrid), e2%rvArray(imgrid), e3%cvArray(imgrid))
    enddo
  else
     write(0, *) 'Error::diagMult_crvector_mg vectors not same subgrids size'

  endif

end subroutine diagMult_crvector_mg
    ! **********************************************************************************************
    subroutine diagMult_rcvector_mg(e1, e2, e3)

    implicit none
    type (rvector_mg), intent(in)  :: e1
    type (cvector_mg), intent(in)  :: e2
    type (cvector_mg), intent(inout)  :: e3
    ! local
    integer  :: imgrid

    ! check whether e1 and e2 have the same number of subgrids
      if ((e1%mgridSize == e2%mgridSize).and.(e1%mgridSize == e3%mgridSize).and. &
                                       (e2%mgridSize == e3%mgridSize)) then

        do imgrid = 1, e1%mgridSize
          call diagMult_rcvector(e1%rvArray(imgrid), e2%cvArray(imgrid), e3%cvArray(imgrid))
        enddo

      else
         write(0, *) 'Error::diagMult_rcvector vectors not same subgrids size'

      endif

    end subroutine diagMult_rcvector_mg
! ********************************************************************************************88
subroutine subtract_rvector_mg(e1, e2, e3)

implicit none
type (rvector_mg), intent(in)  :: e1, e2
type (rvector_mg), intent(inout)  :: e3
! local variables
integer  :: imgrid

! check whether e1 and e2 have the same number of subgrids
  if ((e1%mgridSize == e2%mgridSize).and.(e1%mgridSize == e3%mgridSize).and. &
                                   (e2%mgridSize == e3%mgridSize)) then
    do imgrid = 1, e1%mgridSize
       call subtract_rvector(e1%rvArray(imgrid), e2%rvArray(imgrid), e3%rvArray(imgrid))
     enddo

  else

     write(0, *) 'Error::subtract_rvector_mg vectors not same subgrids size'

  endif

end subroutine subtract_rvector_mg
! ********************************************************************************************88
subroutine subtract_cvector_mg(e1, e2, e3)

implicit none
type (cvector_mg), intent(in)  :: e1, e2
type (cvector_mg), intent(inout)  :: e3
! local variables
integer  :: imgrid

! check whether e1 and e2 have the same number of subgrids
  if ((e1%mgridSize == e2%mgridSize).and.(e1%mgridSize == e3%mgridSize).and. &
                                   (e2%mgridSize == e3%mgridSize)) then
    do imgrid = 1, e1%mgridSize
       call subtract_cvector(e1%cvArray(imgrid), e2%cvArray(imgrid), e3%cvArray(imgrid))
     enddo

  else

     write(0, *) 'Error::subtract_cvector_mg vectors not same subgrids size'

  endif

end subroutine subtract_cvector_mg

  ! ***************************************************************************
  ! dotProd_cvector_mg computes dot product of two vectors stored
  ! as derived data type cvector_mg, returning a complex number
  function dotProd_cvector_mg_f(e1, e2) result(c)

  implicit none

    type (cvector_mg), intent(in)  :: e1, e2
    complex(kind=prec)  :: c
    ! local
    type (cvector_mg) :: e3
    integer  :: imgrid
    complex (kind=prec),allocatable  :: ctemp(:)

    c = R_ZERO

    if((.not.e1%allocated).or.(.not.e2%allocated)) then
      write(0,*) 'RHS not allocated yet for dotProdCC multi-grid'
      return
    endif

    e3 = e1

   ! check whether e1 and e2 have the same number of subgrids
   if (e1%mgridSize == e2%mgridSize) then
     allocate(ctemp(e1%mgridSize))

     do imgrid = 1, e1%mgridSize
       ! delete duplicate edges
       e3%cvarray(imgrid)%x(:,:,e1%cvarray(imgrid)%nz+1) = C_ZERO
       e3%cvarray(imgrid)%y(:,:,e1%cvarray(imgrid)%nz+1) = C_ZERO

       ctemp(imgrid) = dotProd_cvector_f(e3%cvArray(imgrid), e2%cvArray(imgrid))

     enddo

  else
    write(0,*) 'Error :: dotProd_cvector_mg_f: e1 and e2 not the same subgrids'
  end if

  c = sum(ctemp)

  deallocate(ctemp)
  call deall_cvector_mg(e3)

end function dotProd_cvector_mg_f

! **************************************************************************8
! dotProd_noConj_cvector_mg computes dot product of two vectors stored
! as derived data type cvector_mg, returning a complex number
function dotProd_noConj_cvector_mg_f(e1, e2) result(c)

implicit none
type (cvector_mg), intent(in)  :: e1, e2
complex(kind=prec)  :: c
! local
integer  :: imgrid
complex (kind=prec),allocatable  :: ctemp(:)

! check whether e1 and e2 have the same number of subgrids
 if (e1%mgridSize == e2%mgridSize) then
   allocate(ctemp(e1%mgridSize))
   do imgrid = 1, e1%mgridSize
     ctemp(imgrid) = dotProd_cvector_f(e1%cvArray(imgrid), e2%cvArray(imgrid))
   enddo
 else
    write(0,*) 'Error :: dotProd_noConj_cvector_mg_f: e1 and e2 not the same subgrids'
 end if
 c = sum(ctemp)

  deallocate(ctemp)

end function dotProd_noConj_cvector_mg_f
! ***************************************************************************8
subroutine diagDiv_crvector_mg(e1,e2,e3)! may be it is not needed

 implicit none
 type (cvector_mg), intent(in)  :: e1
 type (rvector_mg), intent(in)  :: e2
 type (cvector_mg), intent(inout)  :: e3

 integer   :: imgrid

 if((.not.e1%allocated).or.(.not.e2%allocated)) then
   print*, 'RHS not allocated yet for diagDiv_crvector_mg'
   return
 endif

 ! check to see if LHS (e3) is active (allocated)
 if(.not.e3%allocated) then
   print *, 'LHS was not allocated for diagDiv_crvector_mg'
 else
   do imgrid = 1, e1%mgridSize ! loop on subgrids
     Call diagDiv_crvector(e1%cvArray(imgrid), e2%rvArray(imgrid), e3%cvArray(imgrid))

   enddo
 endif

end subroutine diagDiv_crvector_mg
! ******************************************************************************
 subroutine diagDiv_rcvector_mg(e1, e2, e3) ! may be it is not needed

 implicit none
 type (rvector_mg), intent(in)  :: e1
 type (cvector_mg), intent(in)  :: e2
 type (cvector_mg), intent(inout)  :: e3

 integer   :: imgrid

 if((.not.e1%allocated).or.(.not.e2%allocated)) then
   print*, 'RHS not allocated yet for diagDiv_rcvector_mg'
   return
 endif

 ! check to see if LHS (e3) is active (allocated)
 if(.not.e3%allocated) then
   print *, 'LHS was not allocated for diagDiv_rcvector_mg'
 else
   do imgrid = 1, e1%mgridSize ! loop on subgrids
     Call diagDiv_rcvector(e1%rvArray(imgrid), e2%cvArray(imgrid), e3%cvArray(imgrid))
   enddo
 endif

end subroutine diagDiv_rcvector_mg
! *****************************************************************************
 subroutine scMultAdd_cvector_mg(c, e1, e2)
! Multigrid case

! NOt DEBUGED

    implicit none
    complex(kind=prec), intent(in)  :: c
    ! a complex scalar to be multiplied with
    type (cvector_mg), intent(in)   :: e1
    type (cvector_mg)  :: e2

! local variables
   integer  :: imgrid

    if(.not.e1%allocated) then
       print *, 'RHS not allocated yet for scMultAdd_cvector_mg'
       return
    endif

    ! check to see if LHS (e2) is active (allocated)
    if(.not.e2%allocated) then
       print *, 'LHS was not allocated for scMultAdd_cvector_mg'
    else

      if ((e1%gridType == e2%gridType)) then
        do imgrid = 1, e1%mgridSize  ! loop on subgrids
          ! Check whether both vectors are of the same size
          if((e1%cvArray(imgrid)%nx == e2%cvArray(imgrid)%nx).and.(e1%cvArray(imgrid)%ny == e2%cvArray(imgrid)%ny) &
                                                  .and.(e1%cvArray(imgrid)%nz == e2%cvArray(imgrid)%nz)) then

             ! complex scalar multiplication for x,y,z-components
             e2%cvArray(imgrid)%x = e2%cvArray(imgrid)%x + e1%cvArray(imgrid)%x * c
             e2%cvArray(imgrid)%y = e2%cvArray(imgrid)%y + e1%cvArray(imgrid)%y * c
             e2%cvArray(imgrid)%z = e2%cvArray(imgrid)%z + e1%cvArray(imgrid)%z * c
       else
         print *, 'Error:scMultAdd_cvector_mg: vectors not same size'

       end if
        enddo ! loop on subgrids
      else
        print *, 'not compatible usage for scMultAdd_cvector'
      end if

    end if

  end subroutine scMultAdd_cvector_mg ! scMultAdd_cvector MULTIGRID case
! ******************************************************************************************************
  subroutine updateZ_cvector(outE, whichZ,imgrid)
  ! created 24.04.2012
  ! modified 27.04.2012
  ! modified 30.05.2012

  ! Updates first z layer in cvector_mg (E fields)
  implicit none

  type(cvector_mg), intent(inout)  :: outE
  integer, intent(in)  :: imgrid
  character(len=10), intent(in)  :: whichZ


  ! local variables
  integer  :: ix, iy, iz
  character(len=10),parameter  :: first = 'first'
  character(len=10),parameter  :: last = 'last'

    if(whichZ.eq.first)then
      if(imgrid == 1) then
        ! do nothing
        ! this is upper boundary
        ! leave ...%z(:,:,1) as it is
        return
      endif

        if(outE%coarseness(imgrid).lt.outE%coarseness(imgrid-1))then
          ! interface : coarser grid to finer
          ! copy fields from the previous sub-grid (nz+1 layer) to the current sub-grid (1 layer)
          do iy = 1, outE%cvArray(imgrid-1)%ny
            do ix = 1, outE%cvArray(imgrid-1)%nx
              ! copy Ex odd values
              outE%cvArray(imgrid)%x(2*ix-1,2*iy-1,1) = outE%cvArray(imgrid-1)%x(ix,iy,outE%cvArray(imgrid-1)%nz+1)
              ! copy Ex even values
              outE%cvArray(imgrid)%x(2*ix,2*iy,1) = outE%cvArray(imgrid)%x(2*ix-1,2*iy-1,1)
              ! copy Ex odd values
              outE%cvArray(imgrid)%y(2*ix-1,2*iy-1,1) = outE%cvArray(imgrid-1)%y(ix,iy,outE%cvArray(imgrid-1)%nz+1)
              ! copy Ex even values
              outE%cvArray(imgrid)%y(2*ix,2*iy,1) = outE%cvArray(imgrid)%y(2*ix-1,2*iy-1,1)
            enddo
          enddo

       else if (outE%coarseness(imgrid).gt.outE%coarseness(imgrid-1)) then
         ! interface : finer grid to coarser

         do iy = 1,  outE%cvArray(imgrid)%ny
           do ix = 1,  outE%cvArray(imgrid)%nx
             outE%cvArray(imgrid)%x(ix,iy,1) =  outE%cvArray(imgrid-1)%x(2*ix-1,2*iy-1,outE%cvArray(imgrid-1)%nz+1)
             !outE%cvArray(imgrid-1)%x(2*ix,2*iy,outE%cvArray(imgrid-1)%nz+1) = &
             !  outE%cvArray(imgrid-1)%x(2*ix-1,2*iy-1,outE%cvArray(imgrid-1)%nz+1)

             outE%cvArray(imgrid)%y(ix,iy,1) =  outE%cvArray(imgrid-1)%y(2*ix-1,2*iy-1,outE%cvArray(imgrid-1)%nz+1)
             !outE%cvArray(imgrid-1)%y(2*ix,2*iy,outE%cvArray(imgrid-1)%nz+1) = &
             !  outE%cvArray(imgrid-1)%y(2*ix-1,2*iy-1,outE%cvArray(imgrid-1)%nz+1)
           enddo
         enddo

       else if (outE%coarseness(imgrid) .eq. outE%coarseness(imgrid-1)) then
        ! coarseness does not change
        ! copy fields

         do iy = 1,  outE%cvArray(imgrid)%ny
           do ix = 1,  outE%cvArray(imgrid)%nx
             outE%cvArray(imgrid)%x(ix,iy,1) = outE%cvArray(imgrid-1)%x(ix,iy,outE%cvArray(imgrid-1)%nz+1)
             outE%cvArray(imgrid)%y(ix,iy,1) = outE%cvArray(imgrid-1)%y(ix,iy,outE%cvArray(imgrid-1)%nz+1)
           enddo
         enddo
      endif

    else if(whichZ.eq.last)then

      if(imgrid == outE%mgridSize) then
        ! do nothing
        ! this is lower boundary
        ! leave ...%z(:,:,nz+1) as it is
      return
      endif
        if(outE%coarseness(imgrid).lt.outE%coarseness(imgrid+1))then
          ! interface : finer to coarser
          ! copy fields from the previous subgrid (1 layer) to the current subgrid (nz+2 layer)
          do iy = 1, outE%cvArray(imgrid+1)%ny
            do ix = 1, outE%cvArray(imgrid+1)%nx
              outE%cvArray(imgrid)%x(2*ix-1,2*iy-1,outE%cvArray(imgrid)%nz+1) = outE%cvArray(imgrid+1)%x(ix,iy,1)
              outE%cvArray(imgrid)%x(2*ix,2*iy,outE%cvArray(imgrid)%nz+1) = outE%cvArray(imgrid)%x(2*ix-1,2*iy-1,outE%cvArray(imgrid)%nz+1)

              outE%cvArray(imgrid)%y(2*ix-1,2*iy-1,outE%cvArray(imgrid)%nz+1) = outE%cvArray(imgrid+1)%y(ix,iy,1)
              outE%cvArray(imgrid)%y(2*ix,2*iy,outE%cvArray(imgrid)%nz+1) = outE%cvArray(imgrid)%y(2*ix-1,2*iy-1,outE%cvArray(imgrid)%nz+1)
            enddo
          enddo

        else if (outE%coarseness(imgrid).gt.outE%coarseness(imgrid+1)) then
          ! interface : coarser to finer
          do iy = 1, outE%cvArray(imgrid)%ny
            do ix = 1, outE%cvArray(imgrid)%nx
              outE%cvArray(imgrid)%x(ix,iy,outE%cvArray(imgrid)%nz+1) = outE%cvArray(imgrid+1)%x(2*ix-1,2*iy-1,1)
              !outE%cvArray(imgrid+1)%x(2*ix,2*iy,1) = outE%cvArray(imgrid+1)%x(2*ix,2*iy,1)

              outE%cvArray(imgrid)%y(ix,iy,outE%cvArray(imgrid)%nz+1) = outE%cvArray(imgrid+1)%y(2*ix-1,2*iy-1,1)
              !outE%cvArray(imgrid+1)%y(2*ix,2*iy,1) = outE%cvArray(imgrid+1)%y(2*ix,2*iy,1)
            enddo
          enddo

        else if (outE%coarseness(imgrid) .eq. outE%coarseness(imgrid+1)) then
          ! coarseness does not change
          ! copy fields

      do iy = 1, outE%cvArray(imgrid)%ny
         do ix = 1, outE%cvArray(imgrid)%nx
           outE%cvArray(imgrid)%x(ix,iy,outE%cvArray(imgrid)%nz+1) = outE%cvArray(imgrid+1)%x(ix,iy,1)
           outE%cvArray(imgrid)%y(ix,iy,outE%cvArray(imgrid)%nz+1) = outE%cvArray(imgrid+1)%y(ix,iy,1)
         enddo
      enddo
        endif
    endif

end subroutine updateZ_cvector

! ***********************************************************************************************************
  subroutine updateZ_rvector(outE, imgrid)
  ! created 25.04.2012
  !modified 25.05.2012

  ! Updates first z layer in rvector_mg
  ! starting from the second subgrid
  ! in the first subgris z1 == 0 (upper boundary)
  ! in the next subgrids first layer must be equal to the last layer in the previous subgrid
  ! e(igrid)(:,:,1) == E(igrid-1)(:,:,nz)
  implicit none

  type(rvector_mg), intent(inout)  :: outE
  integer, intent(in)  :: imgrid

  ! local variables
  integer  :: ix, iy, iz
  integer  :: nx, ny, nz


    if(imgrid == 1) then
    ! do nothing
    ! this is upper boundary
    ! z1 == 0 as it is
      return
    endif

      nz = outE%rvArray(imgrid)%nz
    if(outE%coarseness(imgrid).lt.outE%coarseness(imgrid-1))then

      ! interface : coarser grid to finer
      ! copy fields from the prioues subgrid (nz layer) to the current subgrid (1 layer)
      do iy = 1, outE%rvArray(imgrid+1)%ny
        do ix = 1, outE%rvArray(imgrid+1)%nx
           outE%rvArray(imgrid)%x(ix*2-1,iy*2-1,1) = &
             outE%rvArray(imgrid-1)%x(ix,iy,outE%rvArray(imgrid-1)%nz+1)
           outE%rvArray(imgrid)%x(ix*2,iy*2,1) = outE%rvArray(imgrid)%x(ix*2-1,iy*2-1,1)

           outE%rvArray(imgrid)%y(ix*2-1,iy*2-1,1) = &
             outE%rvArray(imgrid-1)%y(ix,iy,outE%rvArray(imgrid-1)%nz+1)
           outE%rvArray(imgrid)%y(ix*2,iy*2,1) = outE%rvArray(imgrid)%y(ix*2-1,iy*2-1,1)
        enddo
      enddo

    else if (outE%coarseness(imgrid).gt.outE%coarseness(imgrid-1)) then

      ! interface : finer grid to coarser
      ! now copies fileds

      ! probably we need averaging here !!!!!!!!!

      do iy = 1, outE%rvArray(imgrid)%ny
        do ix = 1, outE%rvArray(imgrid)%nx
          outE%rvArray(imgrid)%x(ix,iy,1) = &
            outE%rvArray(imgrid-1)%x(2*(ix-1)+1,2*(iy-1)+1,outE%rvArray(imgrid-1)%nz+1)
          outE%rvArray(imgrid-1)%x(2*(ix-1)+2,2*(iy-1)+2,outE%rvArray(imgrid-1)%nz+1) = &
            outE%rvArray(imgrid-1)%x(2*(ix-1)+1,2*(iy-1)+1,outE%rvArray(imgrid-1)%nz+1)

          outE%rvArray(imgrid)%y(ix,iy,1) = &
            outE%rvArray(imgrid-1)%y(2*(ix-1)+1,2*(iy-1)+1,outE%rvArray(imgrid-1)%nz+1)
          outE%rvArray(imgrid-1)%y(2*(ix-1)+2,2*(iy-1)+2,outE%rvArray(imgrid-1)%nz+1) = &
            outE%rvArray(imgrid-1)%y(2*(ix-1)+1,2*(iy-1)+1,outE%rvArray(imgrid-1)%nz+1)
        enddo
      enddo

    else if (outE%coarseness(imgrid) .eq. outE%coarseness(imgrid-1)) then

    ! coarseness does not change
    ! copy fields

      do iy = 1, outE%rvArray(imgrid)%ny
         do ix = 1, outE%rvArray(imgrid)%nx
           outE%rvArray(imgrid)%x(ix,iy,1) = outE%rvArray(imgrid-1)%x(ix,iy,outE%rvArray(imgrid-1)%nz+1)
           outE%rvArray(imgrid)%y(ix,iy,1) = outE%rvArray(imgrid-1)%y(ix,iy,outE%rvArray(imgrid-1)%nz+1)
         enddo
      enddo

    endif

end subroutine updateZ_rvector
! *********************************************************************************************************


end module sg_vector_mg
