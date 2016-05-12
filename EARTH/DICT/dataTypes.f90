! *****************************************************************************
module dataTypes
  ! This module contains the data type dictionary (typeDict) for EARTH

  use math_constants
  use utilities
  use iotypes

  implicit none

  public			:: initTF, deall_TFList

  ! ***************************************************************************
  ! * type functional_t contains the definitions of a data functional used in
  ! * the expression for the penalty functional; we will further specify the rules
  ! * by which they are computed, perhaps by pointing at a function
  type :: functional_t

    ! isComplex is also stored in the dataVector; the sole purpose of the dataType
    ! should be to identify the mapping from the EM solution to the data vector.
    logical                                 :: isComplex = .true.
    ! Here specify rules used to compute this data functional
    character(80)                           :: name
    ! nComp is number of *complex* EM components observed
    integer                                 :: nComp
    ! The weight of this data type in the calculation of misfit
    real(8)                                 :: w

  end type functional_t

  ! ***************************************************************************
  ! * contains the list of transfer functions
  type :: TF_List

    integer                                     :: n
    type (functional_t), pointer, dimension(:)  :: info !nfunc

  end type TF_List

  ! data type dictionary must be public; some attributes are referenced
  !   by top-level inversion routines
  type (TF_List), save, public                  :: TFList


Contains

  ! ***************************************************************************
  ! * initTF reads the file fn_func that contains the information about the
  ! * number and the types of data functionals to use
  subroutine initTF(cUserDef,TFList)

    implicit none
    type (userdef_control), intent(in)                           :: cUserDef
    type (TF_List), intent(out)                             :: TFList
    integer                                                 :: num,i,ios

    open(ioDT,file=cUserDef%fn_func,status='old',form='formatted',iostat=ios)

    write(6,*) node_info,'Reading from the transfer functions file ',trim(cUserDef%fn_func)
    read(ioDT,'(a)') label

    read(ioDT,*) num
    allocate(TFList%info(num))
    do i=1,num

      read(ioDT,*) TFList%info(i)%name,TFList%info(i)%nComp,TFList%info(i)%w

    end do

    close(ioDT)

    TFList%n = num

  end subroutine initTF ! initTF

! **************************************************************************
! Cleans up and deletes type dictionary at end of program execution
  subroutine deall_TFList()

	integer     :: j, istat

	if (associated(TFList%info)) then
        deallocate(TFList%info,STAT=istat)
    end if

  end subroutine deall_TFList

end module dataTypes
