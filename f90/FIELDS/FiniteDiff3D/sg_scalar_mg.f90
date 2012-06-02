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
     module procedure diagMult_rcscalar_mg
end interface

! zeros the scalar (center or corner) nodes
interface zero
  module procedure zero_rscalar_mg
end interface

 ! computes the dot product
interface dotProd
  module procedure dotProd_rscalar_mg_f
  module procedure dotProd_cscalar_mg_f
end interface

interface scMult
  module procedure scMult_cscalar_mg
end interface

interface scMultAdd
  module procedure scMultAdd_cscalar_mg
endinterface

! subtracts the scalar (center or corner) nodes
interface subtract
     module procedure subtract_cscalar_mg
end interface

interface linComb
  module procedure linComb_cscalar_mg
end interface

interface UpdateZ
  module procedure updateZ_RS
  module procedure updateZ_CS
end interface

public  ::    create_rscalar_mg,  create_cscalar_mg, &
              deall_rscalar_mg, deall_cscalar_mg, &
              copy_rscalar_mg, copy_cscalar_mg, &
              diagMult_rscalar_mg, diagMult_rcscalar_mg, scMult_cscalar_mg, scMultAdd_cscalar_mg, &
              zero_rscalar_mg, dotProd_rscalar_mg_f,subtract_cscalar_mg,linComb_cscalar_mg, &
              updateZ_RS, updateZ_CS

! **********************************************************************************
type  :: rscalar_mg

  ! number of subgrids
  integer  :: mgridSize
  ! rscalar for subgrids
  type(rscalar), pointer  :: rscArray(:)
  ! coarseness
  integer, allocatable  :: coarseness(:)
  ! allocated: true if rscArray allocated
  logical  :: allocated = .false.
  ! in GridDef as a parameter: CENTER, CORNER, CELL_EARTH
  character (len=80)  :: gridType

end type   rscalar_mg

! **********************************************************************************
type  :: cscalar_mg

  ! number of subgrids
  integer  :: mgridSize
  ! cscalar for subgrids
  type(cscalar), pointer  :: csArray(:)
  ! coarseness
  integer, allocatable  :: coarseness(:)
  ! allocated: true if csArray allocated
  logical  :: allocated = .false.
  ! store the intention of the use in a character string defined
  ! in GridDef as a parameter: CENTER, CORNER, CELL_EARTH
  character (len=80)  :: gridType
  ! pointer to parent grid
  type (grid_t), pointer  :: grid

end type   cscalar_mg


Contains

  !****************************************************************************
  ! create_rscalar_mg creates variable of derived type rscalar_mg on multigrid. Uses create_rscalar
  ! allocates memory in v component array
  ! gridType is a character string to describe intended usage
  subroutine create_rscalar_mg(igrid, e, gridType)

  implicit none
    type(grid_t), target, intent(in)  :: igrid
    type (rscalar_mg), intent(inout)  :: e
    character (len=80)  :: gridType
    ! local
    integer  :: status, imgrid

    ! First deallocate anything, that's allocated
    call deall_rscalar_mg(e)

    e%mgridSize = igrid%mgridSize
    ! gridType
    e%gridType = gridType
    allocate(e%rscArray(e%mgridSize), STAT = status)

    do imgrid = 1, igrid%mgridSize
      call create_rscalar(igrid%gridArray(imgrid), e%rscArray(imgrid), gridType)
    enddo
    ! allocate memory for coarseness array
    allocate(e%coarseness(igrid%mgridSize), STAT = status)
    e%allocated = .true.

    ! copy coarseness from multigrid
    e%coarseness = igrid%coarseness

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
    integer  :: status, imgrid

    if(e%allocated) then
      call deall_cscalar_mg(e)
    end if
    ! Set pointer
    e%grid => igrid
    e%mgridSize = igrid%mgridSize
    e%gridType = gridType

    allocate(e%csArray(e%mgridSize), STAT = status)

    do imgrid = 1, e%mgridSize
      call create_cscalar(igrid%gridArray(imgrid), e%csArray(imgrid), gridType)
    enddo
    ! allocate memory for coarseness array
    allocate(e%coarseness(igrid%mgridSize), STAT = status)
    e%allocated = .true.

    ! copy coarseness from multigrid
    e%coarseness = igrid%coarseness

end subroutine create_cscalar_mg

  ! **************************************************************************************
  ! deall_rscalar_mg destoys variable of derived type rscalar_mg,
  ! deallocating memory
  subroutine deall_rscalar_mg(e)

  implicit none
    type (rscalar_mg)  :: e
    !local
    integer  :: status, imgrid

    do imgrid = 1, e%mgridSize
      call deall_rscalar(e%rscArray(imgrid))
    enddo
    ! deallocate memory for csArray
    deallocate(e%rscArray, STAT = status)
    e%mgridSize = 0
    e%gridType = ''
    e%allocated = .false.

end subroutine deall_rscalar_mg

  ! ****************************************************************************************
  ! deall_cscalar_mg destoys variable of derived type cscalar_mg,
  ! deallocating memory
  subroutine deall_cscalar_mg(e)

  implicit none
    type (cscalar_mg)  :: e
    !local
    integer  :: status, imgrid

    do imgrid = 1, e%mgridSize
      call deall_cscalar(e%csArray(imgrid))
    enddo
    ! deallocate memory for csArray
    deallocate(e%csArray, STAT = status)
    if(associated(e%grid)) nullify(e%grid)

    e%mgridSize = 0
    e%gridType = ''
    e%allocated = .false.

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

    if(.not.e2%allocated) then
      call create_cscalar_mg(e1%grid, e2, e1%gridType)
    endif

! check whether e1 and e2 have the same number of subgrids
  if (e1%mgridSize == e2%mgridSize) then

    do imgrid = 1, e1%mgridSize

      call copy_cscalar(e2%csArray(imgrid), e1%csArray(imgrid))

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

  !****************************************************************************
  ! diagMult_rcscalar multiplies scalar E1 with scalar E2 stored as
  ! derived type cscalar pointwise; subroutine version
  ! E3 can overwrite E1 or E2

  subroutine diagMult_rcscalar_mg(e1, e2, e3)

    implicit none
    type (rscalar_mg), intent(in)  :: e1
    type (cscalar_mg), intent(in)   :: e2
    type (cscalar_mg), intent(inout)   :: e3

    integer  :: imgrid

   if((.not.e1%allocated).or.(.not.e2%allocated)) then
       write(0,*) 'RHS not allocated yet for diagMult_rcscalar'
       stop
    endif

    ! check to see if LHS (E3) is active (allocated)
    if(.not.e3%allocated) then
      print *,'LHS was not allocated for diagMult_rcscalar'
    else
      !check mgridSize
      if ((e1%mgridSize == e2%mgridSize).and.(e1%mgridSize == e3%mgridSize))then

        do imgrid = 1, e1%mgridSize
          call diagMult_rcscalar(e1%rscArray(imgrid), e2%csArray(imgrid), e3%csArray(imgrid))
        enddo
      else
        print *, 'Error mgridSize in diagMult_rcscalar_mg'
      endif
    endif

  end subroutine diagMult_rcscalar_mg ! diagMult_rcscalar

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
! ***************************************************************************
 function dotProd_cscalar_mg_f(e1, e2) result(c)

    implicit none
    type (cscalar_mg), intent(in)   :: e1, e2
    complex(kind=prec)           :: c

    ! local
    complex(kind=prec), allocatable  :: ctemp(:)
    integer  :: imgrid

    c = C_ZERO

    if((.not.e1%allocated).or.(.not.e2%allocated)) then
      print *, 'RHS not allocated yet for dotProd_cscalar_mg'
      stop
    endif

    if (e1%mgridSize == e2%mgridSize) then
      allocate(ctemp(e1%mgridSize))

      if(e1%gridType == e2%gridType) then

        do imgrid = 1, e1%mgridSize
          ! Check whether both input scalars are of the same size
          if((e1%csArray(imgrid)%nx == e2%csArray(imgrid)%nx).and.(e1%csArray(imgrid)%ny == e2%csArray(imgrid)%ny) &
                                    .and.(e1%csArray(imgrid)%nz == e2%csArray(imgrid)%nz)) then

          c = c+dotProd_cscalar_f(e1%csArray(imgrid), e2%csArray(imgrid))

          else
            print *, 'Error:dotProd_cscalar: scalars not same size'
          end if
        enddo

      else
        print *, 'not compatible usage for dotProd_cscalar'
      end if
    else
      print *, 'Error:dotProd_cscalar_mg: scalars not same ,mgridsize'
    endif


  deallocate(ctemp)
  end function dotProd_cscalar_mg_f ! dotProd_cscalar_mg
!****************************************************************************
! scMult_cscalar multiplies scalar stored as devired data type
! cscalar with a complex scalar; subroutine version
! E2 can overwrite E1

! MULTIGRID case

! NOT DEBUGED

subroutine scMult_cscalar_mg(c, e1, e2)

implicit none
  complex(kind=prec), intent(in)  :: c
  ! a complex scalar to be multiplied with
  type (cscalar_mg), intent(in)   :: e1
  type (cscalar_mg), intent(inout)  :: e2
  ! local
  integer  :: imgrid

  if(.not.e1%allocated) then
    print *, 'RHS not allocated yet for scMult_cscalar_mg'
    stop
  endif

  ! check to see if LHS (E2) is active (allocated)
  if(.not.e2%allocated) then
    print *, 'LHS was not allocated yet for scMult_cscalar'
  else

    do imgrid = 1, e1%mgridSize

      call scMult_cscalar(c,e1%csArray(imgrid), e2%csArray(imgrid))

    enddo
  endif
end subroutine scMult_cscalar_mg

  !****************************************************************************
  ! scMultadd_cscalar multiplies scalar E1 stored as derived data type
  ! cscalar with a complex scalar c, adding result to output scalar E2
  ! MULTIGRID CASE

  ! NOT DEBUGGED

  subroutine scMultAdd_cscalar_mg(c, e1, e2)

    implicit none
    complex(kind=prec), intent(in)  :: c
    ! a complex scalar to be multiplied with
    type (cscalar_mg), intent(in)   :: e1
    type (cscalar_mg), intent(inout)  :: e2

    integer  :: imgrid

    if(.not.e1%allocated) then
      print *, 'RHS not allocated yet for scMultAdd_cscalar'
      stop
    endif

    ! check to see if LHS (E2) is active (allocated)
    if(.not.e2%allocated) then
      print *, 'LHS was not allocated for scMultAdd_cscalar'
    else
      if ((e1%gridType == e2%gridType).and.(e1%mgridSize == e2%mgridSize)) then
        do imgrid = 1, e1%mgridSize
         ! Check whether both scalars are of the same size
         if((e1%csArray(imgrid)%nx == e2%csArray(imgrid)%nx).and.(e1%csArray(imgrid)%ny == e2%csArray(imgrid)%ny) &
                                  .and.(e1%csArray(imgrid)%nz == e2%csArray(imgrid)%nz)) then
             ! complex scalar multiplication for v-component
             e2%csArray(imgrid)%v = e2%csArray(imgrid)%v + e1%csArray(imgrid)%v * c
         else
           print *, 'Error:scMultAdd_cscalar: scalars not same size'
         end if
        enddo
      else
        print *, 'not compatible usage for scMultAdd_cscalar'
      endif
    end if

  end subroutine scMultAdd_cscalar_mg ! scMultAdd_cscalar



! *************************************************************************************
! NOT DEBUGED
  subroutine subtract_cscalar_mg(e1, e2, e3)

    implicit none
    type (cscalar_mg), intent(in) :: e1, e2
    type (cscalar_mg), intent(inout)  :: e3

    ! local
    integer  :: nx, ny, nz
    integer  :: imgrid

    if((.not.e1%allocated).or.(.not.e2%allocated)) then
       print *, 'RHS not allocated yet for subtract_cscalar_mg'
       stop
    endif

    ! check to see if LHS (E2) is active (allocated)
    if(.not.e3%allocated) then
       print *, 'LHS not allocated for subtract_cscalar_mg'
    else

      if ((e1%gridType == e2%gridType).and.(e1%gridType == e3%gridType).and.(e1% mgridSize== e2%mgridSize) &
                                    .and.(e1%mgridSize == e3%mgridSize)) then

        do imgrid = 1, e1%mgridSize
          nx = e1%csArray(imgrid)%nx
          ny = e1%csArray(imgrid)%ny
          nz = e1%csArray(imgrid)%nz
         ! Check whether all the scalar nodes are of the same size
          if((nx == e2%csArray(imgrid)%nx).and.(ny == e2%csArray(imgrid)%ny).and.(nz == e2%csArray(imgrid)%nz).and.&
              (nx == e3%csArray(imgrid)%nx).and.(ny == e3%csArray(imgrid)%ny).and.(nz == e3%csArray(imgrid)%nz)) then
            ! subtract v-component
            e3%csArray(imgrid)%v = e1%csArray(imgrid)%v - e2%csArray(imgrid)%v
          else
            print *, 'Error:subtract_cscalar_mg: scalars not same size'
          end if
        enddo
       else
         print *, 'not compatible usage for subtract_cscalar_mg'
       end if
    end if

  end subroutine subtract_cscalar_mg ! subtract_cscalar_mg
! ******************************************************************************************************************
! MULTIGRID case

! NOT DEBUGGED

subroutine linComb_cscalar_mg(inc1, e1, inc2, e2, e3)

  implicit none
  !   input scalars
  type (cscalar_mg), intent(in)  :: e1, e2
  !  input complex scalars
  complex (kind=prec), intent(in)  :: inc1, inc2
  type (cscalar_mg), intent(inout)  :: e3

  integer  :: imgrid

  if((.not.e1%allocated).or.(.not.e2%allocated)) then
    print *,'RHS not allocated yet for linComb_cscalar_mg'
    stop
  endif

  ! check to see if LHS (E3) is active (allocated)
  if(.not.e3%allocated) then
     print *, 'LHS was not allocated for linComb_cscalar'
  else

    if ((e1%gridType == e2%gridType).and.(e1%gridType == e3%gridType).and. &
        (e1%mgridSize == e2%mgridSize).and.(e1%mgridSize == e3%mgridSize)) then

      do imgrid = 1, e1%mgridSize

       ! Check whether all scalars are of the same size
       if ((e1%csArray(imgrid)%nx == e2%csArray(imgrid)%nx).and.(e1%csArray(imgrid)%ny == e2%csArray(imgrid)%ny).and. &
           (e1%csArray(imgrid)%nz == e2%csArray(imgrid)%nz).and.&
            (e1%csArray(imgrid)%nx == e3%csArray(imgrid)%nx).and.(e1%csArray(imgrid)%ny == e3%csArray(imgrid)%ny).and. &
            (e1%csArray(imgrid)%nz == e3%csArray(imgrid)%nz)) then

             ! form linear combinatoin
             e3%csArray(imgrid)%v = inc1*e1%csArray(imgrid)%v + inc2*e2%csArray(imgrid)%v

        else
          print *, 'Error:linComb_cscalar:  scalars not same size'
        end if

      enddo
    else
      print *, 'not compatible usage for linComb_cscalar'
    end if

  end if

  end subroutine linComb_cscalar_mg ! linComb_cscalar_mg

! ***********************************************************************************************************
  subroutine updateZ_RS(outSc, imgrid)
  ! created 25.04.2012
  ! modified 25.05.2012

  implicit none

    type (rscalar_mg), intent(inout)   :: outSc
    integer, intent(in)  :: imgrid

    ! local variables
    integer  :: ix, iy, iz
    integer  :: nx, ny, nz

      nx = outSc%rscArray(imgrid)%nx
      ny = outSc%rscArray(imgrid)%ny
      nz = outSc%rscArray(imgrid)%nz

      if(imgrid == outSc%mgridSize) then
      ! the last multigrid layer
      ! 0 layer must be already filled in
      ! nz+1 should be zero
      ! nothing to do
      return
      endif

      ! the basic strategy is to fill in the interface layer with values from the finer grid!
      if (outSc%coarseness(imgrid).lt.outSc%coarseness(imgrid+1))then
        ! interface layer: finer -> coarser
        ! Update 1 of Div taking valies from nz of the previous sub-grid
        ! decided just to copy values from finer grid
        do iy = 1, ny/2
          do ix = 1, nx/2
            outSc%rscArray(imgrid+1)%v(ix,iy,1) = outSc%rscArray(imgrid)%v(2*(ix-1)+1, 2*(iy-1)+1,nz+1)
          enddo
        enddo

      else if (outSc%coarseness(imgrid).gt.outSc%coarseness(imgrid+1)) then
       ! interface : finer grid to coarser
       ! vice versa
       do iy = 1, ny
         do ix = 1, nx
           outSc%rscArray(imgrid+1)%v(2*(ix-1)+1, 2*(iy-1)+1,1) = outSc%rscArray(imgrid)%v(ix,iy,nz+1)
           outSc%rscArray(imgrid+1)%v(2*(ix-1)+2, 2*(iy-1)+2,1) = outSc%rscArray(imgrid)%v(2*(ix-1)+1, 2*(iy-1)+1,1)
         enddo
       enddo

     else if (outSc%coarseness(imgrid).eq.outSc%coarseness(imgrid+1)) then
       do ix = 1, nx
         do iy = 1, ny
           outSc%rscarray(imgrid+1)%v(ix,iy,1) = outSc%rscarray(imgrid)%v(ix,iy,nz)
         enddo
       enddo
     endif

end subroutine updateZ_RS
! **************************************************************************************************************
! Updates first and last z layer in Div result
  subroutine UpdateZ_CS(outSc, imgrid)
  ! created 24.05.2012

  implicit none

    type (cscalar_mg), intent(inout)   :: outSc
    integer, intent(in)  :: imgrid

    ! local variables
    integer  :: ix, iy, iz
    integer  :: nx, ny, nz

      nx = outSc%csArray(imgrid)%nx
      ny = outSc%csArray(imgrid)%ny
      nz = outSc%csArray(imgrid)%nz

      if(imgrid == outSc%mgridSize) then
      ! the last multigrid layer
      ! 0 layer must be already filled in
      ! nz+1 should be zero
      ! nothing to do
      return
      endif

      ! the basic strategy is to fill in the interface layer with values from the finer grid!
      if (outSc%coarseness(imgrid).lt.outSc%coarseness(imgrid+1))then
        ! interface layer: finer -> coarser
        ! Update 1 of Div taking valies from nz of the previous sub-grid
        ! decided just to copy values from finer grid
        do iy = 1, ny/2
          do ix = 1, nx/2
            outSc%csArray(imgrid+1)%v(ix,iy,1) = outSc%csArray(imgrid)%v(2*(ix-1)+1, 2*(iy-1)+1,nz+1)
          enddo
        enddo

      else if (outSc%coarseness(imgrid).gt.outSc%coarseness(imgrid+1)) then
       ! interface : finer grid to coarser
       ! vice versa
       do iy = 1, ny
         do ix = 1, nx
           outSc%csArray(imgrid+1)%v(2*(ix-1)+1, 2*(iy-1)+1,1) = outSc%csArray(imgrid)%v(ix,iy,nz+1)
           outSc%csArray(imgrid+1)%v(2*(ix-1)+2, 2*(iy-1)+2,1) = outSc%csArray(imgrid)%v(2*(ix-1)+1, 2*(iy-1)+1,1)
         enddo
       enddo

     else if (outSc%coarseness(imgrid).eq.outSc%coarseness(imgrid+1)) then
       do ix = 1, nx
         do iy = 1, ny
           outSc%csarray(imgrid+1)%v(ix,iy,1) = outSc%csarray(imgrid)%v(ix,iy,nz+1)
         enddo
       enddo
     endif

  end  subroutine UpdateZ_CS



! *********************************************************************************************************

end module sg_scalar_mg
