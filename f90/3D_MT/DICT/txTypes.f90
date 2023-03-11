! *****************************************************************************
module txTypes
  ! This module contains the transmitter types, i.e., the different types of
  ! forward problems that are supported.
  ! It links to both transmitter dictionary (txDict) and data type dictionary (typeDict).
  ! It contains bookkeeping rather than computational details. Will only be used
  ! for bookkeeping purposes (e.g., input and output routines); not in ForwardSolver.
  ! Right now, the program does NOT allow for extensions to non-EM problems such as
  ! gravity and seismic. To make this possible, a substantial rethinking would need
  ! to be required, allowing for custom, problem-dependent penalty functional terms.

  use dataTypes

  implicit none
  private

  type :: transmitterType_t

     logical, pointer, dimension(:) :: definedTypes

  end type transmitterType_t

  ! transmitter type dictionary will only be used for bookkeeping;
  ! currently in I/O but possibly later in inversion
  type (transmitterType_t), pointer, save, public, dimension(:) :: txTypeDict

  integer, parameter, public   :: MT = 1
  integer, parameter, public   :: CSEM = 2
  integer, parameter, public   :: SFF = 3
  integer, parameter, public   :: TIDE = 4
  integer, parameter, public   :: GLOBAL = 5

!**************************************************************************
! Initializes and sets up transmitter type dictionary
  subroutine setup_txTypeDict()

     integer     :: istat

     allocate(txTypeDict(5),STAT=istat)
     do iTxt=1,size(txTypeDict)
        allocate(txTypeDict(iTxt)%definedTypes(nDt),STAT=istat)
        txTypeDict(iTxt)%definedTypes = .false.
     end do


     ! MT
     txTypeDict(MT)%definedType(Full_Impedance) = .true.

     allocate(txTypeDict(MT)%dataTypes(6),STAT=istat)
     txTypeDict(MT)%dataTypes(1) = Full_Impedance
     txTypeDict(MT)%dataTypes(2) = Off_Diagonal_Impedance
     txTypeDict(MT)%dataTypes(3) = Full_Vertical_Components
     txTypeDict(MT)%dataTypes(4) = Full_Interstation_TF
     txTypeDict(MT)%dataTypes(5) = Off_Diagonal_Rho_Phase
     txTypeDict(MT)%dataTypes(6) = Phase_Tensor

     allocate(txTypeDict(SFF)%dataTypes(6),STAT=istat)
     txTypeDict(SFF)%dataTypes(1) = Ex_Field
     txTypeDict(SFF)%dataTypes(2) = Ey_Field
     txTypeDict(SFF)%dataTypes(3) = Bx_Field
     txTypeDict(SFF)%dataTypes(4) = By_Field
     txTypeDict(SFF)%dataTypes(5) = Bz_Field
     txTypeDict(SFF)%dataTypes(6) = Full_Impedance

     allocate(txTypeDict(CSEM)%dataTypes(5),STAT=istat)
     txTypeDict(CSEM)%dataTypes(1) = Ex_Field
     txTypeDict(CSEM)%dataTypes(2) = Ey_Field
     txTypeDict(CSEM)%dataTypes(3) = Bx_Field
     txTypeDict(CSEM)%dataTypes(4) = By_Field
     txTypeDict(CSEM)%dataTypes(5) = Bz_Field

     allocate(txTypeDict(TIDE)%dataTypes(5),STAT=istat)
     txTypeDict(TIDE)%dataTypes(1) = Ex_Field
     txTypeDict(TIDE)%dataTypes(2) = Ey_Field
     txTypeDict(TIDE)%dataTypes(3) = Bx_Field
     txTypeDict(TIDE)%dataTypes(4) = By_Field
     txTypeDict(TIDE)%dataTypes(5) = Bz_Field

     allocate(txTypeDict(GLOBAL)%dataTypes(6),STAT=istat)
     txTypeDict(GLOBAL)%dataTypes(1) = Bx_Field
     txTypeDict(GLOBAL)%dataTypes(2) = By_Field
     txTypeDict(GLOBAL)%dataTypes(3) = Bz_Field
     txTypeDict(GLOBAL)%dataTypes(4) = C_Response
     txTypeDict(GLOBAL)%dataTypes(5) = D_Response
     txTypeDict(GLOBAL)%dataTypes(6) = Q_Response

  end subroutine setup_txTypeDict

! **************************************************************************
! Cleans up and deletes transmitter type dictionary at end of program execution
  subroutine deall_txTypeDict()

    integer     :: j, istat

    if (associated(txTypeDict)) then

       do j = 1,size(txTypeDict)
          if (associated(txTypeDict(j)%dataTypes)) then
             deallocate(txTypeDict(j)%dataTypes,STAT=istat)
          end if
       end do

       deallocate(txTypeDict,STAT=istat)

    end if

  end subroutine deall_txTypeDict

!**********************************************************************
! Figures out the data type from its name

  function txType(typeName) result (txType)

    character(*), intent(in)    :: typeName
    integer                     :: txType

    select case (trim(typeName))

       case('MT')
          txType = MT

       case('SFF')
          txType = SFF

       case('CSEM')
          txType = CSEM

       case('TIDE')
          txType = TIDE

       case('GLOBAL')
          txType = GLOBAL

       case default
          call errStop('Unknown transmitter type:'//trim(typeName))

    end select

  end function txType



end module txTypes
