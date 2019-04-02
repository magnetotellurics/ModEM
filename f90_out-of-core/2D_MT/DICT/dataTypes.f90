! *****************************************************************************
module dataTypes
  ! This module contains the data type dictionary (typeDict) for 2D MT

  use math_constants
  use utilities

  implicit none

  type :: dataType

     !  stores information about the "data type"
     !   The following two attributes must be defined for all
     !    data types; these are accessed and used by the top-level
     !    inversion routines.
     logical                    :: isComplex = .false.
     !    Other attributes might be different (different number,
     !        different names, types, etc.) for  different applications.
     ! character(2)                :: mode = ''! = 'TE' or 'TM'
     character(80)               :: name = ''
     ! the units of the data type
     character(80)              :: units
     ! the number of components in the data type
     integer           			:: nComp
     !  could add rxDictNumber to keep track of reciever dictionary
     !  number used for this dataType (only 1 receiver dictionary now,
     !   so this is omitted)
     ! id(nComp) are the text comments that describe the components; these are
     ! initialized, but can be overwritten by the info from the data file.
     character(15), pointer, dimension(:) :: id

     ! these lists contain the indices into the data vector for each data type;
     ! they make it possible to sort the data by receiver for output.
     ! no data denoted by zero index; dimensions (nTx) and (nTx,nRx).
     integer, pointer, dimension(:)   :: tx_index
     integer, pointer, dimension(:)   :: dt_index
     integer, pointer, dimension(:,:) :: rx_index

  end type dataType

  ! data type dictionary must be public; some attributes are referenced
  !   by top-level inversion routines
  type (dataType), pointer, save, public, dimension(:) :: typeDict

  ! add data types here ... this all needs work!
  integer, parameter    :: TE_Impedance = 1
  integer, parameter    :: TM_Impedance = 2
  integer, parameter    :: Tzy_Impedance = 3
  integer, parameter    :: Rho_Phs_TM = 4  
Contains

!**************************************************************************
! Initializes and sets up data type dictionary
  subroutine setup_typeDict()

  	 integer     :: istat

     allocate(typeDict(4),STAT=istat)
     
     typeDict(TE_Impedance)%name = 'TE_Impedance'
     typeDict(TE_Impedance)%isComplex = .true.
     typeDict(TE_Impedance)%units     = '[V/m]/[T]'
     typeDict(TE_Impedance)%nComp     = 2
     allocate(typeDict(TE_Impedance)%id(1),STAT=istat)
     typeDict(TE_Impedance)%id(1)     = 'TE'
     
     typeDict(TM_Impedance)%name = 'TM_Impedance'
     typeDict(TM_Impedance)%isComplex = .true.
     typeDict(TM_Impedance)%units     = '[V/m]/[T]'
     typeDict(TM_Impedance)%nComp     = 2
     allocate(typeDict(TM_Impedance)%id(1),STAT=istat)
     typeDict(TM_Impedance)%id(1)     = 'TM'
     
     typeDict(Tzy_Impedance)%name = 'Tzy_Impedance'
     typeDict(Tzy_Impedance)%isComplex = .true.
     typeDict(Tzy_Impedance)%units     = '[]'
     typeDict(Tzy_Impedance)%nComp     = 2
     allocate(typeDict(Tzy_Impedance)%id(1),STAT=istat)
     typeDict(Tzy_Impedance)%id(1)     = 'Ty'
     
     typeDict(Rho_Phs_TM)%name = 'Rho_Phs_TM'
     typeDict(Rho_Phs_TM)%isComplex = .false.
     typeDict(Rho_Phs_TM)%units     = '[]'
     typeDict(Rho_Phs_TM)%nComp     = 2
     allocate(typeDict(Rho_Phs_TM)%id(2),STAT=istat)
     typeDict(Rho_Phs_TM)%id(1)     = 'Rho_TM'
     typeDict(Rho_Phs_TM)%id(2)     = 'Phs_TM'   
     
  end subroutine setup_typeDict

! **************************************************************************
! Cleans up and deletes type dictionary at end of program execution
  subroutine deall_typeDict()

	integer     :: j,istat

	if (associated(typeDict)) then

	   do j = 1,size(typeDict)
	      if (associated(typeDict(j)%id)) then
	         deallocate(typeDict(j)%id,STAT=istat)
	      end if
          if (associated(typeDict(j)%tx_index)) then
             deallocate(typeDict(j)%tx_index,STAT=istat)
          end if
          if (associated(typeDict(j)%dt_index)) then
             deallocate(typeDict(j)%dt_index,STAT=istat)
          end if
          if (associated(typeDict(j)%rx_index)) then
             deallocate(typeDict(j)%rx_index,STAT=istat)
          end if
	   end do

       deallocate(typeDict,STAT=istat)

    end if

  end subroutine deall_typeDict

!**********************************************************************
! Computes the value by which the data must be multiplied to convert
! from the old units to the new units.
! The units may be any of the following.
! 1) SI units for E/B: [V/m]/[T] (used in ModEM code)
! 2) practical units for E/B: [mV/km]/[nT]
! 3) SI units for E/H: [V/m]/[A/m] = Ohm

  function ImpUnits(oldUnits,newUnits) result (SI_factor)

	character(*), intent(in)    :: oldUnits, newUnits
	real(kind=prec)             :: SI_factor
	! local
	real(kind=prec)             :: factor1, factor2

    
    if ((index(oldUnits,'[]')>0) .or. (index(newUnits,'[]')>0)) then
	   SI_factor = ONE
	   return
    end if
        
	! first convert the old units to [V/m]/[T]
	if (index(oldUnits,'[V/m]/[T]')>0) then
	   ! SI units for E/B
	   factor1 = ONE
	else if (index(oldUnits,'[mV/km]/[nT]')>0) then
	   ! practical units for E/B
	   factor1 = ONE * 1000.0
	else if ((index(oldUnits,'[V/m]/[A/m]')>0) .or. (index(oldUnits,'Ohm')>0)) then
	   ! SI units for E/H
	   factor1 = ONE * 1000.0 * 10000.0/(4*PI) ! approx. 796000.0
	else
	   call errStop('Unknown input units in ImpUnits: '//trim(oldUnits))
	end if

	! now convert [V/m]/[T] to the new units
	if (index(newUnits,'[V/m]/[T]')>0) then
	   ! SI units for E/B
	   factor2 = ONE
	else if (index(newUnits,'[mV/km]/[nT]')>0) then
	   ! practical units for E/B
	   factor2 = ONE / (1000.0)
	else if ((index(newUnits,'[V/m]/[A/m]')>0) .or. (index(newUnits,'Ohm')>0)) then
	   ! SI units for E/H
	   factor2 = ONE / (1000.0 * 10000.0/(4*PI))
	else
	   call errStop('Unknown output units in ImpUnits: '//trim(newUnits))
	end if

	SI_factor = factor1 * factor2

  end function ImpUnits

!**********************************************************************
! Figures out the data type from the mode name

  function ImpType(mode) result (dataType)

    character(*), intent(in)    :: mode
    integer                     :: dataType

    
     select case (trim(mode))
         case('TE_Impedance')
            dataType = TE_Impedance
         case('TM_Impedance')  
             dataType = TM_Impedance
         case('Tzy_Impedance')
             dataType = Tzy_Impedance
         case('Rho_Phs_TM')
             dataType = Rho_Phs_TM
         case default
             call errStop('Unknown data type:'//trim(mode))             
     end select  

  end function ImpType


!**********************************************************************
! Figures out the component index from its name for any data type

  function ImpComp(compid,dataType) result (icomp)

    character(*), intent(in)    :: compid
    integer,      intent(in)    :: dataType
    integer                     :: icomp
    ! local
    integer                     :: i, ncomp

    ncomp = typeDict(dataType)%ncomp
    if (typeDict(dataType)%isComplex) then
        ncomp = ncomp/2
    end if
    icomp = 0

    do i = 1,ncomp
        if (index(typeDict(dataType)%id(i),trim(compid))>0) then
            icomp = i
            exit
        end if
    end do

    if (icomp == 0) then
        call errStop('Problem locating the component '//trim(compid)//' in data type '//trim(typeDict(dataType)%name))
    end if

  end function ImpComp


end module dataTypes
