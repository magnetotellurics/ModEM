! *****************************************************************************
module sg_scalar_mg ! for MULTIgrid
  ! This module creates (allocates) rscalar, cscalar on staggered MULTIgrid
  ! Deallocates  and does some basic algebraic operations.
  ! For algebraic operations, new routines use the old ones (from sg_scalar)
  ! and do loop on subgrids

  ! Derived cscalar_mg and rscalar_mg data types are defined in sg_scalar
  ! module (after cscalar, rscalar respectively)

  use math_constants		! math/ physics constants
  use griddef
  use sg_scalar
  implicit none

  ! Generic interfaces are done through subroutines
  ! creates scalar (center or corner) nodes
interface create
  module procedure create_rscalar_mg
  module procedure create_cscalar_mg
end interface

  ! deallocates the scalar (center or corner) nodes
interface deall
  module procedure deall_rscalar_mg
  module procedure deall_cscalar_mg
end interface

interface assignment (=)
   module procedure copy_cscalar_mg
   module procedure copy_rscalar_mg
end interface

! pointwise scalar multiplication of scalar (center or corner) nodes
interface diagMult
     module procedure diagMult_rscalar_mg
end interface

! zeros the scalar (center or corner) nodes
interface zero
  module procedure zero_rscalar_mg
end interface

 ! computes the dot product
interface dotProd
  module procedure dotProd_rscalar_mg_f
end interface

public  ::    create_rscalar_mg,  create_cscalar_mg, &
              deall_rscalar_mg, deall_cscalar_mg, &
              copy_rscalar_mg, copy_cscalar_mg, &
              diagMult_rscalar_mg, &
              zero_rscalar_mg, dotProd_rscalar_mg_f

! **********************************************************************************
type  :: rscalar_mg

! number of subgrids
  integer  :: mgridSize
! rscalar for subgrids
  type(rscalar), pointer  :: rscArray(:)
! allocated: true if rscArray allocated
  logical  :: allocated = .false.

end type   rscalar_mg

! **********************************************************************************
type  :: cscalar_mg

! number of subgrids
  integer  :: mgridSize
! cscalar for subgrids
  type(cscalar), pointer  :: cscalar_mg(:)

end type   cscalar_mg

! ***********************************************************************************
Contains

!****************************************************************************
! create_rscalar_mg creates variable of derived type rscalar_mg
! in MULTIgrid. Uses create_rscalar
! allocates memory in v component array
! gridType is a character string to describe intended usage
subroutine create_rscalar_mg(igrid, e, gridType)

implicit none
type(grid_t), target, intent(in)  :: igrid
! the MULTIgrid for which an scalar (center or corner) node field is being
! initialized
type (rscalar_mg), intent(inout)  :: e
character (len=80)  :: gridType
! local
integer  :: imgrid

  e%mgridSize = igrid%mgridSize

  allocate(e%rscArray(e%mgridSize))

  do imgrid = 1, e%mgridSize

    call create_rscalar(igrid, e%rscArray(imgrid), gridType)

  enddo

end subroutine create_rscalar_mg

! ***************************************************************************
! create_cscalar_mg creates variable of derived type cscalar_mg
! in MULTIgrid. Uses create_cscalar
! allocates memory in v component array
! gridType is a character string to describe intended usage
subroutine create_cscalar_mg(igrid, e, gridType)

implicit none
type (grid_t), target, intent(in)   :: igrid
! the grid for which an scalar (center or corner) node field is being
! initialized
type (cscalar_mg), intent(inout)  :: e
character (len=80)  :: gridType
! local
integer  :: imgrid

  e%mgridSize = igrid%mgridSize

  allocate(e%cscalar_mg(e%mgridSize))

  do imgrid = 1, e%mgridSize

    call create_cscalar(igrid, e%cscalar_mg(imgrid), gridType)

  enddo

end subroutine create_cscalar_mg
! **************************************************************************************
! deall_rscalar_mg destoys variable of derived type rscalar_mg,
! deallocating memory
subroutine deall_rscalar_mg(e)

implicit none
type (rscalar_mg)  :: e
!local
integer  :: imgrid

  do imgrid = 1, e%mgridSize
    call deall_rscalar(e%rscArray(imgrid))
  enddo

  e%mgridSize = 0


end subroutine deall_rscalar_mg
! ****************************************************************************************
! deall_cscalar_mg destoys variable of derived type cscalar_mg,
! deallocating memory
subroutine deall_cscalar_mg(e)

implicit none
type (cscalar_mg)  :: e
!local
integer  :: imgrid

  do imgrid = 1, e%mgridSize
    call deall_cscalar(e%cscalar_mg(imgrid))
  enddo

  e%mgridSize = 0

end subroutine deall_cscalar_mg
! ***************************************************************************************
! copy_rscalar_mg makes an exact copy of derived data type
! rscalar_mg;   NOTE: first argument is output
subroutine copy_rscalar_mg(e2, e1)

implicit none
type (rscalar_mg), intent(in)  :: e1
type (rscalar_mg), intent(inout)   :: e2
! local
integer  :: imgrid

! check whether e1 and e2 have the same number of subgrids
  if (e1%mgridSize == e2%mgridSize) then

    do imgrid = 1, e1%mgridSize

      call copy_rscalar(e2%rscArray(imgrid), e1%rscArray(imgrid))

    enddo

  else

    write(0,*) 'Error :: copy_rscalar_mg: e1 and e2 not the same subgrids'

  end if

end subroutine copy_rscalar_mg

! ************************************************************************************
! copy_cscalar_mg makes an exact copy of derived data type
! cscalar_mg;   NOTE: first argument is output
subroutine copy_cscalar_mg(e2, e1)

implicit none
type (cscalar_mg), intent(in)  :: e1
type (cscalar_mg), intent(inout)   :: e2
! local
integer  :: imgrid

! check whether e1 and e2 have the same number of subgrids
  if (e1%mgridSize == e2%mgridSize) then

    do imgrid = 1, e1%mgridSize

      call copy_cscalar(e2%cscalar_mg(imgrid), e1%cscalar_mg(imgrid))

    enddo

  else

    write(0,*) 'Error :: copy_cscalar_mg: e1 and e2 not the same subgrids'

  end if

end subroutine copy_cscalar_mg
! ****************************************************************************************88
! diagMult_rscalar_mg multiplies two scalars e1, e2 stored as devired data
! type rscalar_mg pointwise;
subroutine diagMult_rscalar_mg(e1, e2, e3)

implicit none
type (rscalar_mg), intent(in)  :: e1, e2
type (rscalar_mg), intent(inout)  :: e3
! local
integer  :: imgrid

! check whether e1 and e2 have the same number of subgrids
  if ((e1%mgridSize == e2%mgridSize).and.(e1%mgridSize == e3%mgridSize).and. &
                                    (e2%mgridSize == e3%mgridSize))  then

    do imgrid = 1, e1%mgridSize

      call diagMult_rscalar(e1%rscArray(imgrid), e2%rscArray(imgrid), &
                                              e3%rscArray(imgrid))

    enddo

  else

    write(0,*) 'Error :: diagMult_rscalar_mg : e1,e2,e3 not the same subgrids'

  end if

end subroutine diagMult_rscalar_mg
! ****************************************************************************************
! zero_rscalar_mg zeros variable of derived data type rscalar_mg
subroutine zero_rscalar_mg(e)

implicit none
type (rscalar_mg), intent(inout)  :: e
!local
integer  :: imgrid

  do imgrid = 1, e%mgridSize
    call zero_rscalar(e%rscArray(imgrid))
  enddo

end subroutine zero_rscalar_mg
! ****************************************************************************************88
! dotProd_rscalar_mg function computes dot product of two vecors stored
! as derived data type rscalar_mg, returning a real number
function dotProd_rscalar_mg_f(e1, e2) result(r)

implicit none
type (rscalar_mg), intent(inout)  :: e1, e2
real (kind=prec)  :: r
!local
integer  :: imgrid
real (kind=prec),allocatable  :: rtemp(:)

 ! check whether e1 and e2 have the same number of subgrids
  if (e1%mgridSize == e2%mgridSize) then
    allocate(rtemp(e1%mgridSize))
    do imgrid = 1, e1%mgridSize
      rtemp(imgrid) = dotProd_rscalar_f(e1%rscArray(imgrid), e2%rscArray(imgrid))
    enddo
  else
    write(0,*) 'Error :: dotProd_rscalar_mg_f: e1 and e2 not the same subgrids'
  end if
   r = sum(rtemp)

  deallocate(rtemp)

end function dotProd_rscalar_mg_f

end module sg_scalar_mg
