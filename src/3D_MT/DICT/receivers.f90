! *****************************************************************************
module receivers
  ! This module contains the receiver dictionary (rxDict) for 3D MT

  use math_constants

  implicit none

  public			:: setup_rxDict, update_rxDict, deall_rxDict

  !  multiple receiver dictionaries can be defined, as
  !   different sorts of meta-data may be required for different data
  !   types; e.g., intersite TFs would require two locations.
  !   Note that the location must be specified relative to the same
  !     coordinate system used for the grid, with consistent units:
  !     in 3D MT, the origin of the grid is given relative to
  !     a refernce origin, with everything in meters.  The data site
  !     locations should also be given in meters, relative to this
  !     same physical origin.  All 3 components are required, to support
  !     observations anywhere within the solution domain.
  !  Additonal data types, to be used as elements of additional
  !    dictionaries can be added to accomodate additional data types
  type :: MTrx
     ! x gives location of EM measurements;
  	 ! x(1) points North, x(2) points East, x(3) points down

  	 ! NM: add addtional vector to store the location of a reference station.
     ! Same as x: r(1) points North, r(2) points East, r(3) points down
     real (kind=prec)                   ::  x(3)
     real (kind=prec)                   ::  r(3)
	 real (kind=prec)                   ::  Rx_Azi
     ! site ID used for input/output and for searching through the list
     character(50)                      ::  id=''
     character(50)                      ::  id_ref=''
  end type MTrx

  ! receiver dictionary for 3D MT data will be an array of
  !  type MTrx (one element of the array for each site)
  !  Note that receiver dictionaries are only used inside the
  !  data functional module, and can thus be private
  type (MTrx), pointer, save, public, dimension(:) :: rxDict


Contains

  ! **************************************************************************
  ! Initializes and sets up receiver dictionary
  ! The reciever dictionary contains sparse vectors required
  ! for magnetic and electric field vector evaluation
  subroutine setup_rxDict(nSites,siteLocations,siteIDs)

    integer, intent(in)	 				:: nSites
    real(kind=prec), intent(in)			:: siteLocations(nSites,3)
    character(*), intent(in), optional  :: siteIDs(nSites)


    !  local variables
    integer		 :: i,istat
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
! NM: modified to include referance site info.

function update_rxDict(loc,id,Rx_azi,loc_ref,id_ref) result (iRx)

     character(*), intent(in)            :: id
     real(kind=prec), intent(in)         :: loc(3)
     real(kind=prec):: lat,lon
	 real(kind=prec),intent(in),optional :: Rx_azi
     real(kind=prec),intent(in),optional :: loc_ref(3)
     character(*), intent(in),optional   :: id_ref
     integer                             :: iRx
     ! local
     type(MTrx)                          :: new
     type(MTrx), pointer, dimension(:)   :: temp
     integer                             :: nRx, istat,i
     logical							 :: new_Rx

     ! Create a receiver for this location
     new%id = id
     new%x  = loc

     if (present(loc_ref)) then
     	new%r  = loc_ref
     	new%id_ref=id_ref
     end if

     if (present(Rx_azi)) then
     	new%Rx_azi  = Rx_azi
     end if	 
	 

     ! If rxDict doesn't yet exist, create it
     if(.not. associated(rxDict)) then
     	allocate(rxDict(1),STAT=istat)
     	rxDict(1) = new
     	iRx = 1
     	new_Rx = .true.
     	return
     end if

     nRx = size(rxDict)

     ! If this site isn't new, do nothing, unless a new ref. site is found in case of Inter-Stations TF.
     do i = 1,nRx
     	if (new%id .eq. rxDict(i)%id) then
           if (present(loc_ref)) then
           !Check if the this site is associated with the same Ref. site. If not, then continue and append another site to the Rx dictionary. 
              if (new%id_ref .eq. rxDict(i)%id_ref) then 
     	         rxDict(i)%r=loc_ref
     	         rxDict(i)%id_ref=id_ref
                 iRx=i
                 new_Rx = .false.
                 return
              end if
           elseif (present(Rx_azi)) then
		   !Check if the this site has same azimuth as what we have already in the Dic 
               if (new%Rx_azi .eq. rxDict(i)%Rx_azi) then 
                 iRx=i
                 new_Rx = .false.
                 return
              end if              		   
           else    
     	    iRx=i
            new_Rx = .false.
     		return
           end if 
     	end if
     end do

     ! If the site really is new, append to the end of the dictionary
     new_Rx = .true.
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
        write(*,'(i6,a50,4f15.6,a50)') iRx,trim(rxDict(iRx)%id),rxDict(iRx)%x,rxDict(iRx)%Rx_azi,trim(rxDict(iRx)%id_ref)
     enddo

  end subroutine print_rxDict

  ! **************************************************************************
  ! Cleans up and deletes receiver dictionary at end of program execution
  subroutine deall_rxDict()

	integer     :: istat

    if (associated(rxDict)) then
       deallocate(rxDict,STAT=istat)
    end if

  end subroutine deall_rxDict

end module receivers
