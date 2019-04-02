! *****************************************************************************
module transmitters
  ! This module contains the transmitter dictionary (txDict) for 3D MT

  use math_constants

  implicit none

  public			:: setup_txDict,  update_txDict, deall_txDict

  type :: MTtx
     !  An MT source is defined by frequency and boundary conditions
     !   at present there does not seem to be much need for BC info ... add
     !    if needed.  Other sorts of EM data may have more
     !    complex tx descriptions
     ! required attribute - number of polarizations
     integer					:: nPol = 1 ! = 1 for 2D MT
     character(2)               :: mode = ''! = 'TE' or 'TM'
     ! angular frequency (radians/sec), and for convenience period (s)
     real(kind=prec)            :: omega = R_ZERO
     real(kind=prec)            :: period = R_ZERO
     ! index number to frequency/ period in solution file
     integer                    :: iPer
     character(2)               :: Tx_type = 'MT'
   end type MTtx

   ! transmitter dictionary for MT data will be an array of
   ! type mt_forcing (one element  for each frequency)
   ! Perhaps this should be moved to ForwardSolver module (and be private
   !    to this module)
   type (MTtx), pointer, save, public, dimension (:)   :: txDict

Contains

!**********************************************************************
! Initializes and sets up transmitter dictionary for MT,
!  This is just a simple example of a routine for setting up the TX
!   dictionary; In this example we assume that there are nPer periods
!   for either TE or TM modes, or for both.

  subroutine setup_txDict(nTx,Periods,modes)

     integer, intent(in)         :: nTx
     real*8, intent(in)          :: periods(nTx)
     character*2, intent(in)     :: modes(nTx)

     ! local variables
     integer                     :: iTx,istat

	 if (.not. associated(txDict)) then
    	allocate(txDict(nTx),STAT=istat)
     end if

     do iTx = 1, nTx
        txDict(iTx)%period = Periods(iTx)
        txDict(iTx)%omega = (2*PI)/ txDict(iTx)%period
        txDict(iTx)%mode = modes(iTx)
        txDict(iTx)%iPer = iTx
        txDict(iTx)%Tx_type = 'MT'
     enddo

  end subroutine setup_txDict

!**********************************************************************
! Updates the transmitter dictionary for MT with a new period (in secs)
! Returns the index of the new element.
! This is not efficient; but this would only be used a few times, with
! a small number of values, so convenience is much more of an issue here!

  function update_txDict(Period,Mode) result (iTx)

     real(kind=prec), intent(in)        :: Period
     character(2), intent(in)           :: Mode
     integer                            :: iTx
     ! local
     type(MTtx)                         :: new
     type(MTtx), pointer, dimension(:)  :: temp
     integer                            :: nTx, istat

     ! Create a transmitter for this period
     new%period = Period
     new%omega  = (2*PI)/Period
     new%mode   = Mode
     new%iPer   = nTx + 1

     ! If txDict doesn't yet exist, create it
     if(.not. associated(txDict)) then
     	allocate(txDict(1),STAT=istat)
     	txDict(1) = new
     	iTx = 1
     	return
     end if

     nTx = size(txDict)

     ! If this transmitter isn't new, do nothing
     do iTx = 1,nTx
     	if ((abs(Period - txDict(iTx)%period) < TOL6) .and. (Mode .eq. txDict(iTx)%mode)) then
     		return
     	end if
     end do

     ! If the period really is new, append to the end of the dictionary
     allocate(temp(nTx+1),STAT=istat)
     temp(1:nTx) = txDict
     temp(nTx+1) = new
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
        write(*,*) iTx,txDict(iTx)%period,txDict(iTx)%mode,txDict(iTx)%nPol
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

end module transmitters
