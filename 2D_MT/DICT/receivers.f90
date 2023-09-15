! *****************************************************************************
module receivers
  ! This module contains the receiver dictionary (rxDict) for 2D MT

  use math_constants

  implicit none

  public			:: setup_rxDict, update_rxDict, deall_rxDict

  type :: MTrx
     ! x gives location of EM measurements
     !  multiple receiver dictionaries can be defined, and
     !   different dictionaries can be used for different data types
     !  Additonal elements of MTrx data type can be added to
     !   accomodate additional data types
     real(kind=prec)					::  x(2)
     ! site ID used for input/output and for searching through the list
     character(50)                      ::  id=''
  end type MTrx

  ! receiver dictionary for 2D MT data will be an array of
  !  type MTrx (one element of the array for each site)
  !  Two components of MTrx%x are position along the profile,
  !     and vertical position (generally on the surface)
  !  Note that the receiver dictionary is only used inside the
  !    data functional module
  type (MTrx), pointer, save, public, dimension(:) :: rxDict


Contains

  ! **************************************************************************
  ! Initializes and sets up receiver dictionary
  ! Now the reciever dictionary only contains the location of the point obs
  subroutine setup_rxDict(nSites,siteLocations,siteIDs)
    !  siteLocatins(nSites,2) is array of measurement locations (x,z)
    !   corresponding to grid  (normally z = 0 for flat Earth surface)

    integer, intent(in)	 		:: nSites
    real(kind=prec), intent(in)	:: siteLocations(nSites,2)
    character(*), intent(in), optional  :: siteIDs(nSites)

    !  local variables
    integer      :: i,istat
	character(3) :: id

	if (.not. associated(rxDict)) then
    	allocate(rxDict(nSites),STAT=istat)
    end if

    do i = 1,nSites
    	rxDict(i)%x = siteLocations(i,:)
		if (present(siteIDs)) then
			rxDict(i)%id = siteIDs(i)
		else
			write(id,'(i6.6)') i
			rxDict(i)%id = id
		end if
    end do

  end subroutine setup_rxDict

!**********************************************************************
! Updates the receiver dictionary with a new location and site ID.
! Returns the index of the new element.
! This is not efficient; but this would only be used a few times, with
! a small number of values, so convenience is much more of an issue here!

  function update_rxDict(loc,id) result (iRx)

     character(*), intent(in)           :: id
     real(kind=prec), intent(in)        :: loc(2)
     integer                            :: iRx
     ! local
     type(MTrx)                         :: new
     type(MTrx), pointer, dimension(:)  :: temp
     integer                            :: nRx, istat

     ! Create a receiver for this location
     new%id = id
     new%x  = loc

     ! If rxDict doesn't yet exist, create it
     if(.not. associated(rxDict)) then
     	allocate(rxDict(1),STAT=istat)
     	rxDict(1) = new
     	iRx = 1
     	return
     end if

     nRx = size(rxDict)

     ! If this site isn't new, do nothing
     do iRx = 1,nRx
     	if (new%id .eq. rxDict(iRx)%id) then
     		return
     	end if
     end do

     ! If the site really is new, append to the end of the dictionary
     allocate(temp(nRx+1),STAT=istat)
     temp(1:nRx) = rxDict
     temp(nRx+1) = new
     deallocate(rxDict,STAT=istat)
     allocate(rxDict(nRx+1),STAT=istat)
     rxDict = temp
     deallocate(temp,STAT=istat)
     iRx = nRx+1

  end function update_rxDict

!**********************************************************************
! Writes the receiver dictionary to screen. Useful for debugging.

  subroutine print_rxDict()

     ! local variables
     integer                     :: iRx

     if (.not. associated(rxDict)) then
        return
     end if

     write(*,*) 'Receiver dictionary:'
     do iRx = 1, size(rxDict)
        write(*,'(i6,a50,2f15.6)') iRx,trim(rxDict(iRx)%id),rxDict(iRx)%x
     enddo

  end subroutine print_rxDict

!**********************************************************************
! Writes the receiver dictionary to a file -- needed to associate rows in sensitivity
!  sensitivity matrix J with correct data vector elements.   This assumes that iounit
!  is opened for formated io-- file connection and  closing are done by calling routine

  subroutine write_rxDict_bin(iounit)

     ! local variables
     integer                     :: iounit, iRx, nRx
     character(len=80)          :: header

     if (.not. associated(rxDict)) then
        return
     end if

     nRx = size(rxDict)
     header = 'Receiver Dictionary'
     write(iounit) header
     write(iounit) nRx
     do iRx = 1, size(rxDict)
        write(iounit) iRx,rxDict(iRx)%x
        !   hard to read variable length records from a binary file,
        !   unless written as individual sequential records; hence this
        write(iounit) trim(rxDict(iRx)%id)
     enddo

  end subroutine write_rxDict_bin

  ! **************************************************************************
  ! Cleans up and deletes receiver dictionary at end of program execution
  subroutine deall_rxDict()

	integer     :: istat

    if (associated(rxDict)) then
       deallocate(rxDict,STAT=istat)
    end if

  end subroutine deall_rxDict


end module receivers
