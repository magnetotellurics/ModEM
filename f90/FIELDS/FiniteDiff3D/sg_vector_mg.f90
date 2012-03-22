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
 ! module procedure scMultReal_cvector ! do we realy need it
end interface

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

public  :: create_rvector_mg, create_cvector_mg, &
           deall_rvector_mg, deall_cvector_mg, &
           copy_cvector_mg, &
           zero_cvector_mg, &
           subtract_rvector_mg, subtract_cvector_mg, &
           scMult_cvector_mg, &
           add_cvector_mg,add_rvector_mg, &
           diagMult_cvector_mg, diagMult_crvector_mg, &
           diagMult_rvector_mg, diagMult_rcvector_mg, &
           linComb_cvector_mg, &
           dotProd_cvector_mg_f,  dotProd_noConj_cvector_mg_f

! *******************************************************************************
! this type defines cvectors on MULTIgrid
type :: cvector_mg
! number of subgrids
    integer  :: mgridSize
! cvectors for subgrids
    type(cvector), pointer  :: cvector_mg(:)

end type cvector_mg

! ******************************************************************************
! this type defines rvectors on MULTIgrid
type rvector_mg
! number of subgrids
    integer  :: mgridSize
! rvectors for subgrids
    type(rvector), pointer :: rvArray(:)
! allocated:  .true.  x, y, z arrays have been allocated
    logical  :: allocated = .false.
end type rvector_mg

! *******************************************************************************
Contains

! allocates rvector for MULTIgrid
subroutine create_rvector_mg (mgrid, e, gridType)

implicit none
type(grid_t), intent(in)     :: mgrid
! the grid for which an edge/ face node field is being initialized
type (rvector_mg), intent(inout)       :: e
character (len=80), intent(in)      :: gridType
!local
integer  :: imgrid

  e%mgridSize = mgrid%mgridSize

  allocate(e%rvArray(e%mgridSize))

  do imgrid = 1, e%mgridSize
    call create_rvector(mgrid%gridArray(imgrid),e%rvArray(imgrid),gridType)
  enddo

end subroutine create_rvector_mg

! ******************************************************************************

! allocates cvector for multigrid
subroutine create_cvector_mg (mgrid, e, gridType)

implicit none
type(grid_t), intent(in)     :: mgrid
! the grid for which an edge/ face node field is being initialized
type (cvector_mg), intent(inout)  :: e
character (len=80), intent(in)      :: gridType
!local
integer  :: imgrid

  e%mgridSize = mgrid%mgridSize

  allocate(e%cvector_mg(e%mgridSize))

  do imgrid = 1, e%mgridSize
    call create_cvector(mgrid%gridArray(imgrid),e%cvector_mg(imgrid),gridType)
  enddo

end subroutine create_cvector_mg

! *****************************************************************************
subroutine deall_rvector_mg(e)

implicit none
type (rvector_mg)  :: e
integer  :: imgrid

  do imgrid = 1, e%mgridSize
    call deall_rvector(e%rvArray(imgrid))
  enddo

  e%mgridSize = 0

end  subroutine deall_rvector_mg

! *****************************************************************************
subroutine deall_cvector_mg(e)

implicit none
type (cvector_mg)  :: e
integer  :: imgrid

  do imgrid = 1, e%mgridSize
    call deall_cvector(e%cvector_mg(imgrid))
  enddo

  e%mgridSize = 0

end  subroutine deall_cvector_mg

! *****************************************************************************
subroutine copy_cvector_mg(e2,e1)
! first argument is output
implicit none
type (cvector_mg), intent(in)  :: e1
type (cvector_mg), intent(inout)  :: e2
!local
integer  :: imgrid

  e2%mgridSize = e1%mgridSize
  e2%cvector_mg = e1%cvector_mg

  do imgrid = 1, e1%mgridSize
    call copy_cvector(e2%cvector_mg(imgrid),e1%cvector_mg(imgrid))
  enddo

end subroutine copy_cvector_mg

! *****************************************************************************
subroutine zero_cvector_mg(e)

implicit none
type(cvector_mg), intent(inout)  :: e
! local
integer  ::imgrid

  do imgrid = 1, e%mgridSize
    call zero_cvector(e%cvector_mg(imgrid))
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
      call scMult_cvector(c, e1%cvector_mg(imgrid),e2%cvector_mg(imgrid))
    enddo

  else
    write (0, *) 'Error:scMult_cvector_mg: vectors not same subgrids size'
  endif

end subroutine scMult_cvector_mg

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
      call linComb_cvector (inc1, e1%cvector_mg(imgrid), inc2, &
                        e2%cvector_mg(imgrid), e3%cvector_mg(imgrid))
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
      call add_cvector(e1%cvector_mg(imgrid), e2%cvector_mg(imgrid), e3%cvector_mg(imgrid))
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
      call diagMult_cvector(e1%cvector_mg(imgrid), e2%cvector_mg(imgrid), e3%cvector_mg(imgrid))
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
      call diagMult_crvector(e1%cvector_mg(imgrid), e2%rvArray(imgrid), e3%cvector_mg(imgrid))
    enddo

  else
     write(0, *) 'Error::diagMult_crvector vectors not same subgrids size'

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
      call diagMult_rcvector(e1%rvArray(imgrid), e2%cvector_mg(imgrid), e3%cvector_mg(imgrid))
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
       call subtract_cvector(e1%cvector_mg(imgrid), e2%cvector_mg(imgrid), e3%cvector_mg(imgrid))
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
integer  :: imgrid
complex (kind=prec),allocatable  :: ctemp(:)

! check whether e1 and e2 have the same number of subgrids
 if (e1%mgridSize == e2%mgridSize) then
   allocate(ctemp(e1%mgridSize))
   do imgrid = 1, e1%mgridSize
     ctemp(imgrid) = dotProd_cvector_f(e1%cvector_mg(imgrid), e2%cvector_mg(imgrid))
   enddo
 else
    write(0,*) 'Error :: dotProd_cvector_mg_f: e1 and e2 not the same subgrids'
 end if
 c = sum(ctemp)

  deallocate(ctemp)

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
     ctemp(imgrid) = dotProd_cvector_f(e1%cvector_mg(imgrid), e2%cvector_mg(imgrid))
   enddo
 else
    write(0,*) 'Error :: dotProd_noConj_cvector_mg_f: e1 and e2 not the same subgrids'
 end if
 c = sum(ctemp)

  deallocate(ctemp)

end function dotProd_noConj_cvector_mg_f


end module sg_vector_mg
