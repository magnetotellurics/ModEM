! *****************************************************************************
module transmitters
  ! This module contains the transmitter dictionary (txDict) for 3D MT

  use math_constants

  implicit none

  public			:: setup_txDict, update_txDict, deall_txDict

  type :: MTtx
     !  An MT source is defined by frequency and boundary conditions
     !   at present there does not seem to be much need for BC info ... add
     !    if needed.  Other sorts of EM data may have more
     !    complex tx descriptions
     ! required attribute - number of polarizations
     ! this could be different for e.g. active source 3D applications
     integer					:: nPol = 2 ! = 2 for 3D MT
     ! angular frequency (radians/sec), and for convenience period (s)
     real(kind=prec)            :: omega = R_ZERO
     real(kind=prec)            :: period = R_ZERO
     ! index number to frequency/ period in solution file
     integer                    :: iPer
   end type MTtx

   ! transmitter dictionary txDict for 3D-MT data will be an array of
   ! type MTtx (one element  for each frequency)
   ! Perhaps this should be moved to ForwardSolver module (and be private
   !    to that module?)
   ! NOTE: could have multiple transmitter dictionaries, each of which
   !    could constist of elements of different types; nothing about
   !    the dictionary or the elements that it consists of is used
   !    in higher level routines
   type (MTtx), pointer, save, public, dimension (:)   :: txDict


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
! Updates the transmitter dictionary for MT with a new period (in secs)
! Returns the index of the new element.
! This is not efficient; but this would only be used a few times, with
! a small number of values, so convenience is much more of an issue here!

  function update_txDict(Period,nPol) result (iTx)

     real(kind=prec), intent(in)        :: Period
     integer, intent(in), optional		:: nPol
     integer                            :: iTx
     ! local
     type(MTtx)                         :: new
     type(MTtx), pointer, dimension(:)  :: temp
     integer                            :: nTx, istat,i
     logical							:: new_Tx


     ! Create a transmitter for this period
     new%period = Period
     new%omega  = (2*PI)/Period
     if (present(nPol)) then
     	new%nPol = nPol
     else
     	new%nPol = 2
     end if
     new%iPer   = nTx + 1

     ! If txDict doesn't yet exist, create it
     if(.not. associated(txDict)) then
     	allocate(txDict(1),STAT=istat)
     	txDict(1) = new
     	iTx = 1
	    new_Tx = .true.
     	return
     end if


     nTx = size(txDict)
       ! If this period isn't new, do nothing
     do i = 1,nTx
     	if ((abs(Period - txDict(i)%period) .lt. TOL6) .and. (new%nPol == txDict(i)%nPol)) then
     	itx=i
	    new_Tx=.false.
     	return
     	end if
     end do

     ! If the period really is new, append to the end of the dictionary
     new_Tx = .true.
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

end module transmitters
