! *****************************************************************************
module transmitters
  ! This module contains the transmitter dictionary (txDict) for 3D MT and CSEM

  use math_constants

  implicit none

  public			:: setup_txDict, update_txDict, deall_txDict

type :: transmitter_t
!defines what kind of transmitter: MT or CSEM        
           character(10)		    :: tx_type=''
		   
! both MT and CSEM have common attributes:         
		 integer					:: nPol !while setting up the Tx, nPol=2 for MT and 1 for CSEM
		 ! angular frequency (radians/sec), and for convenience period (s)
		 real(kind=prec)            :: omega = R_ZERO
		 real(kind=prec)            :: period = R_ZERO
		 ! index number to frequency/ period in solution file
		 integer                    :: iPer
		  
!######################################################	 		  
! CSEM detalies
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

end type transmitter_t 





   ! transmitter dictionary txDict for 3D-CSEM data will be an array of
   ! type VMDtx (one element  for each transmitter)
   !
   ! NOTE: could have multiple transmitter dictionaries, each of which
   !    could constist of elements of different types; nothing about
   !    the dictionary or the elements that it consists of is used
   !    in higher level routines
   type (transmitter_t), pointer, save, public, dimension (:)   :: txDict


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
! Updates the transmitter dictionary for CSEM with a new source
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
     	if ( IsTxExist(aTx,txDict(iTx) ) ) then
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
!
! **************************************************************************
!
  function IsTxExist(Txa,Txb) result (YESNO)
  type(transmitter_t), intent(in):: Txa
  type(transmitter_t), intent(in):: Txb
  logical                  YESNO

  YESNO = .FALSE.
if (Txa%Tx_type=='DC') then
  if( ABS(Txa%xyzTx(1) - Txb%xyzTx(1)) < TOL6 .AND.   &
      ABS(Txa%xyzTx(2) - Txb%xyzTx(2)) < TOL6 .AND.   &
      ABS(Txa%xyzTx(3) - Txb%xyzTx(3)) < TOL6 )  then
	  YESNO = .true.
  end if
elseif (Txa%Tx_type=='CSEM') then
  if( ABS(Txa%Period - Txb%period) < TOL6  .AND.      &
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
elseif (Txa%Tx_type=='MT') then
      if(ABS(Txa%Period - Txb%period) < TOL6  .and. Txa%nPol == Txb%nPol) then
    YESNO = .true.
  end if
end if
 

  return
  end function IsTxExist
  
!  
! **************************************************************************
! Cleans up and deletes transmitter dictionary at end of program execution

  subroutine deall_txDict()

	integer     :: istat

    if (associated(txDict)) then
       deallocate(txDict,STAT=istat)
    end if

  end subroutine deall_txDict

end module transmitters
