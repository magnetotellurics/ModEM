! *****************************************************************************
module transmitters
  ! This module contains the transmitter dictionary (txDict) for 3D MT

  use math_constants

  implicit none

  public			:: setup_txDict, update_txDict, deall_txDict

  type :: MTtx
     ! defines the kind of transmitter: here, MT only. Needed for compatibility
     ! with more general versions of the code
     character(10)		:: tx_type='MT'
     !  An MT source is defined by frequency and boundary conditions
     !   at present there does not seem to be much need for BC info ... add
     !    if needed.  Other sorts of EM data may have more
     !    complex tx descriptions
     ! required attribute - number of polarizations
     ! this could be different for e.g. active source 3D applications
     integer					:: nPol=2 !while setting up the Tx, nPol=2 for MT and 1 for CSEM
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

  ! transmitter types; correspond to index iTxt in the data vectors
  !  these will be heavily used in inversion routines
  integer, parameter   :: MT = 1
  integer, parameter   :: CSEM = 2

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

  function update_txDict(aTx) result (iTx)

     type(MTtx)              :: aTx
     integer                 :: iTx
	 
     ! local
     type(MTtx), pointer, dimension(:)  :: temp
     integer                            :: nTx, istat,i
     logical							:: new_Tx


     ! Create a transmitter for this period
     nTx = size(txDict)
     aTx%iPer=nTx+1
	 
     ! If txDict doesn't yet exist, create it
     if(.not. associated(txDict)) then
     	allocate(txDict(1),STAT=istat)
     	txDict(1) = aTx
     	iTx = 1
	    new_Tx = .true.
     	return
     end if


     nTx = size(txDict)
     ! If this period isn't new, do nothing
     do iTx = 1,nTx
     	if (IsTxExist(aTx,txDict(iTx) ) ) then
     	  return
     	end if
     end do

     ! If the period really is new, append to the end of the dictionary
     new_Tx = .true.
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
! A function to compare between two transmitter and return YESNO=True if Txa=Txb
  function IsTxExist(Txa,Txb) result (YESNO)
  type(MTtx), intent(in):: Txa
  type(MTtx), intent(in):: Txb
  logical                  YESNO

  YESNO = .FALSE.
if (Txa%Tx_type=='CSEM') then
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
!**********************************************************************
! Writes the transmitter dictionary to a file -- needed to associate a sensitivity
!   (i.e., in full sensitivity matrix J, or for Jmult_MTX computation) with correct
!   data vector elements.   This assumes that iounit is opened for formated io 
!     -- file connection and  closing are done by calling routine

  subroutine write_txDict_asc(iounit)


     ! local variables
     integer                     :: iounit,iTx,nMT,nCSEM

     if (.not. associated(txDict)) then
        return
     end if

     nMT = 0
     nCSEM = 0
     !   this is only coded for MT + CSEM -- could be generalized
     do iTx = 1, size(txDict)
        if (txDict(iTx)%tx_type .eq. 'MT') then
           nMT = nMT + 1
        else
           nCSEM = nCSEM + 1
        endif
     enddo
     write(iounit,'(i6,a20)') nMT,'     MT Transmitters'
     do iTx = 1, size(txDict)
        if (txDict(iTx)%tx_type .eq. 'MT') then
           write(iounit,'(i6,es12.6,i2)') iTx,txDict(iTx)%period,txDict(iTx)%nPol
        endif
     enddo
     write(iounit,'(i6,a15)') nCSEM,'CSEM Transmitters'
     do iTx = 1, size(txDict)
        if (txDict(iTx)%tx_type .eq. 'CSEM') then
           write(iounit,'(i6,2x,es12.6,2x,i2,2x,a8,3f11.2,3f8.2)') iTx,txDict(iTx)%period, &
             txDict(iTx)%nPol,txDict(iTx)%Dipole, txDict(iTx)%xyzTx, & 
             txDict(iTx)%azimuthTx, txDict(iTx)%dipTx,txDict(iTx)%moment
        endif
     enddo

  end subroutine write_txDict_asc
!**********************************************************************
! Writes the transmitter dictionary to a file -- needed to associate a sensitivity
!   (i.e., in full sensitivity matrix J, or for Jmult_MTX computation) with correct
!   data vector elements.   This version assumes that iounit is opened for sequential unformated io 
!     -- file connection and  closing are done by calling routine

  subroutine write_txDict_bin(iounit)


     ! local variables
     integer                     :: iounit,iTx,nMT,nCSEM,nTx
     character(len=80) header

     if (.not. associated(txDict)) then
        return
     end if

     nMT = 0
     nCSEM = 0
     nTx = size(txDict)
     !   this is only coded for MT + CSEM -- could be generalized
     do iTx = 1, nTx 
        if (txDict(iTx)%tx_type .eq. 'MT') then
           nMT = nMT + 1
        else
           nCSEM = nCSEM + 1
        endif
     enddo
     write(iounit) nTx,nMT,nCSEM 
     header = 'Transmitter Dictionary: MT' 
     write(iounit) header
     do iTx = 1, nTx
        if (txDict(iTx)%tx_type .eq. 'MT') then
           write(iounit) iTx,txDict(iTx)%period,txDict(iTx)%nPol
        endif
     enddo
     header = 'Transmitter Dictionary: CSEM' 
     write(iounit) header
     do iTx = 1, nTx
        if (txDict(iTx)%tx_type .eq. 'CSEM') then
           write(iounit) iTx,txDict(iTx)%period,txDict(iTx)%nPol,txDict(iTx)%Dipole, &
             txDict(iTx)%xyzTx,txDict(iTx)%azimuthTx,txDict(iTx)%dipTx,txDict(iTx)%moment
        endif
     enddo
  end subroutine write_txDict_bin

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

    type(MTtx), intent(in):: Txa
    type(MTtx), intent(in):: Txb
    logical                  YESNO

    YESNO = .false.
    if (trim(Txa%Tx_type) .eq. 'MT') then
      if(ABS(Txa%period - Txb%period) < TOL6  .and. Txa%nPol == Txb%nPol) then
        YESNO = .true.
      end if
    else
        write(0,*) 'Unknown transmitter type #',trim(Txa%Tx_type)
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
       case default
          write(0,*) 'Unknown transmitter type: ',trim(tx_type)
    end select

  end function tx_type_index

end module transmitters
