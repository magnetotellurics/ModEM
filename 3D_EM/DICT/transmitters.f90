! *****************************************************************************
module transmitters
  ! This module contains the general EM transmitter dictionary (txDict)
  !
  ! Currently defined are the following problems:
  !	MT      2D and 3D magnetotelluric modeling with plane-wave sources
  !	CSEM    3D controlled source EM
  !     SFF     Secondary field formulation used with any EM primary fields
  !	TIDE    3D EM modeling with tidal sources
  !     GLOBAL  3D EM global with spherical coordinate source representation 
  !     DC      Direct current - old code that is not maintained
  !
  ! Not all of these problems are fully implemented or included in this specific
  ! version of the code. Also, not all of these problems are currently working
  ! in the inversion mode or included in joint inversion.
  ! However, we shall maintain the complete transmitter dictionaries here
  ! to streamline code maintenance.
  !
  ! A. Kelbert, Nov 16, 2022

  use math_constants

  implicit none

  public			:: setup_txDict, update_txDict, deall_txDict

  type :: transmitter_t

     ! defines the kind of transmitter: MT, CSEM, SFF, TIDE, GLOBAL, DC
     character(10)		        :: tx_type=''
     ! attributes common for all transmitter types:
     integer				:: nPol !while setting up the Tx, nPol=2 for MT and 1 for CSEM
     ! angular frequency (radians/sec), and for convenience period (s)
     real(kind=prec)            :: omega = R_ZERO
     real(kind=prec)            :: period = R_ZERO
     ! index number to frequency/ period in solution file
     integer                    :: iPer

!######################################################	 		  
! CSEM details
     ! Specific Dipole Type (Electric or Magnetic)
     character(8)		:: Dipole
     !   location of transmitter, relative to grid 
     real(kind=prec)            :: xyzTx(3)
	 ! Source azimuth from x axis (positive clockwise)
     real(kind=prec)            :: azimuthTx ! (degrees) 
     ! Vertical dip angle of source along azimuthTx, positive down 
     real(kind=prec)            :: dipTx ! (degrees) 
     ! Source dipole moment
     real(kind=prec)            :: moment ! (A.m) for electric, (A.m^2) for magnetic

!######################################################
! Tidal details
    ! in some cases (e.g., tides), might want to give the transmitter a name
    character(20)              :: id = ''
    ! ocean tides also have amplitude which might be useful
    real(kind=prec)            :: amplitude = R_ZERO
    ! internal source for this transmitter, stored as a sparse vector on the grid
    !   this is supported for some rare circumstances (e.g., tides);
    !   doesn't exist for MT problem and should be ignored by most users
    !type(sparsevecc)          :: jInt
    ! for now, hard code the name in the ForwardSolver and read it there.
    !   this is very crude but may just do for our purposes.
    !character(120)            :: fn_intsource = ''		  

  end type transmitter_t


   ! NOTE: could have multiple transmitter dictionaries, each of which
   !    could constist of elements of different types; nothing about
   !    the dictionary or the elements that it consists of is used
   !    in higher level routines
   ! In the future, the plan is to use submodules as soon as these are
   ! universally supported, and have separate submodules to set up each
   ! of the transmitter types (MT, CSEM, tidal etc)
   ! Then the master transmitter dictionary will do the bookkeeping.
   ! e.g., transmitter dictionary txDict for 3D-CSEM data will be an array of
   ! type VMDtx (one element  for each transmitter)
   ! type MTtx (for magnetotellurics)
   ! type TIDEtx (for tidal source)
   type (transmitter_t), pointer, save, public, dimension (:)   :: txDict

  ! transmitter types; correspond to index iTxt in the data vectors
  !  these will be heavily used in inversion routines
  integer, parameter   :: MT = 1
  integer, parameter   :: CSEM = 2
  integer, parameter   :: SFF = 3
  integer, parameter   :: TIDE = 4
  integer, parameter   :: GLOBAL = 5
  integer, parameter   :: DC = 6

Contains

!**********************************************************************
! Initializes and sets up transmitter dictionary for MT,
!  This is just a simple example of a routine for setting up the TX
!   dictionary; In this example we assume that there are nPer periods
!   for either TE or TM modes, or for both.

  subroutine setup_txDict(nTx,Periods,nPol)

     integer, intent(in)            :: nTx
     real(kind=prec), intent(in)    :: Periods(nTx)
     integer, intent(in), optional	:: nPol

     ! local variables
     integer                     :: iTx,istat

     if (.not. associated(txDict)) then
    	allocate(txDict(nTx),STAT=istat)
     end if

     do iTx = 1, nTx
        txDict(iTx)%period = Periods(iTx)
        txDict(iTx)%omega = (2*PI)/ txDict(iTx)%period
        if (present(nPol)) then
        	txDict(iTx)%nPol = nPol
        endif
     enddo

  end subroutine setup_txDict

!**********************************************************************
! Updates the transmitter dictionary with a new source
! Returns the index of the new element.
! This is not efficient; but this would only be used a few times, with
! a small number of values, so convenience is much more of an issue here
!
  function update_txDict(aTx) result (iTx)

     type(transmitter_t)                :: aTx
     integer                            :: iTx
     ! local
     type(transmitter_t), pointer, dimension(:)  :: temp
     integer                            :: nTx, istat

     ! Create a transmitter for this period
     nTx = size(txDict)

     aTx%iPer   = nTx + 1

     ! If txDict doesn't yet exist, create it
     if(.not. associated(txDict)) then
     	allocate(txDict(1),STAT=istat)
     	txDict(1) = aTx
     	iTx = 1
     	return
     end if

     ! If this period isn't new, do nothing
     do iTx = 1,nTx
     	if ( compare_tx(aTx,txDict(iTx) ) ) then
     	  return
     	end if
     end do

     ! If the period really is new, append to the end of the dictionary
     allocate(temp(nTx+1),STAT=istat)
     temp(1:nTx) = txDict
     temp(nTx+1) = aTx
     deallocate(txDict,STAT=istat)
     allocate(txDict(nTx+1),STAT=istat)
     txDict = temp
     deallocate(temp,STAT=istat)
     iTx = nTx+1

  end function update_txDict

!**********************************************************************
! Writes the transmitter dictionary to screen. Useful for debugging.

  subroutine print_txDict()

     ! local variables
     integer                     :: iTx

     if (.not. associated(txDict)) then
        return

     end if


     write(*,*) 'Transmitter dictionary:'
     do iTx = 1, size(txDict)
        write(*,*) iTx,txDict(iTx)%period,txDict(iTx)%nPol
     enddo

  end subroutine print_txDict

! **************************************************************************
! Cleans up and deletes transmitter dictionary at end of program execution

  subroutine deall_txDict()

    integer     :: istat

    if (associated(txDict)) then
       deallocate(txDict,STAT=istat)
    end if

  end subroutine deall_txDict

! **************************************************************************
! Used to compare two transmitters for updating the dictionary

  function compare_tx(Txa,Txb) result (YESNO)

    type(transmitter_t), intent(in):: Txa
    type(transmitter_t), intent(in):: Txb
    logical                  YESNO

    YESNO = .false.
    if (trim(Txa%Tx_type) .eq. 'DC') then
      if( ABS(Txa%xyzTx(1) - Txb%xyzTx(1)) < TOL6 .AND.   &
          ABS(Txa%xyzTx(2) - Txb%xyzTx(2)) < TOL6 .AND.   &
          ABS(Txa%xyzTx(3) - Txb%xyzTx(3)) < TOL6 )  then
          YESNO = .true.
      end if
    elseif (trim(Txa%Tx_type) .eq. 'CSEM') then
      if( ABS(Txa%period - Txb%period) < TOL6  .AND.      &
          ABS(Txa%xyzTx(1) - Txb%xyzTx(1)) < TOL6 .AND.   &
          ABS(Txa%xyzTx(2) - Txb%xyzTx(2)) < TOL6 .AND.   &
          ABS(Txa%xyzTx(3) - Txb%xyzTx(3)) < TOL6 .AND.   &
          ABS(Txa%moment - Txb%moment) < TOL6 .AND. &
          ABS(Txa%azimuthTx - Txb%azimuthTx) < TOL6 .AND. &
          ABS(Txa%dipTx - Txb%dipTx) < TOL6 ) then
          if (Txa%Dipole .Eq. Txb%Dipole) then
              YESNO = .true.
          end if
      end if
    elseif ((trim(Txa%Tx_type) .eq. 'MT') .or. (trim(Txa%Tx_type) .eq. 'SFF')) then
      if(ABS(Txa%period - Txb%period) < TOL6  .and. Txa%nPol == Txb%nPol) then
        YESNO = .true.
      end if
    elseif ((trim(Txa%Tx_type) .eq. 'TIDE') .or. (trim(Txa%Tx_type) .eq. 'GLOBAL')) then
      if (trim(Txa%id) .eq. trim(Txb%id)) then
        YESNO = .true.
      end if
    else
        write(0,*) 'Unknown transmitter type ',trim(Txa%Tx_type)
   end if
 
  end function compare_tx

! **************************************************************************
! Used to extract tx_type character name from transmitter type index iTxt
!
  function tx_type_name(iTxt) result (tx_type)

    integer, intent(in)                 :: iTxt
    character(10)                       :: tx_type

    select case (iTxt)
       case(MT)
          tx_type = 'MT'
       case(CSEM)
          tx_type = 'CSEM'
       case(SFF)
          tx_type = 'SFF'
       case(TIDE)
          tx_type = 'TIDE'
       case(GLOBAL)
          tx_type = 'GLOBAL'
       case(DC)
          tx_type = 'DC'
       case default
          write(0,*) 'Unknown transmitter type #',iTxt
    end select

  end function tx_type_name

! **************************************************************************
! Used to extract transmitter type index iTxt from transmitter type name.
! All this is only needed because key-value lists aren't allowed in Fortran!
! In the future, we should stick to the transmitter integer indicator
! and keep the name for input/output only. The integer is all that
! the data vector should ever know of.
!
  function tx_type_index(tx_type) result (iTxt)

    character(*), intent(in)            :: tx_type
    integer                             :: iTxt

    select case (trim(adjustl(tx_type)))
       case('MT')
          iTxt = MT
       case('CSEM')
          iTxt = CSEM
       case('SFF')
          iTxt = SFF
       case('TIDE')
          iTxt = TIDE
       case('GLOBAL')
          iTxt = GLOBAL
       case('DC')
          iTxt = DC
       case default
          write(0,*) 'Unknown transmitter type: ',trim(tx_type)
    end select

  end function tx_type_index

end module transmitters
