! *****************************************************************************
module receivers
  ! This module contains the receiver dictionary (rxDict) for 3D MT

  use math_constants
  use utilities

  implicit none

  public			:: setup_rxDict, update_rxDict, deall_rxDict
  private      :: count_rx

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
     ! Rx_azi is only well-defined for CSEM; for MT, it is not used since we're supporting
     ! the possibility that each field component has its own orientation (stored in DataSpace).
     real (kind=prec)                   ::  Rx_Azi
     ! site ID used for input/output and for searching through the list
     character(50)                      ::  id=''
     character(50)                      ::  id_ref=''
     ! for efficiency, we may initialize with padding, then trim to only leave what's been read
     logical                            ::  defined = .false.
     integer :: iRx = 0
  end type MTrx

  ! receiver dictionary for 3D MT data will be an array of
  !  type MTrx (one element of the array for each site)
  !  Note that receiver dictionaries are only used inside the
  !  data functional module, and can thus be private
  type (MTrx), pointer, save, public, dimension(:) :: rxDict

  type :: MTrx_hash_map_t
      type (MTrx), pointer :: ptr => null()
  end type MTrx_hash_map_t

  type (MTrx_hash_map_t), dimension(:), pointer :: rx_map => null()

  ! for efficiency, we now initialize with MAX_NRX and trim after reading the file
  integer, parameter                    :: MAX_NRX = 200000
  integer, public :: current_rx_count = 0


Contains


  ! **************************************************************************
  ! Initializes and sets up receiver dictionary
  ! The reciever dictionary contains sparse vectors required
  ! for magnetic and electric field vector evaluation
  subroutine setup_rxDict(nSites, siteLocations,siteIDs)

    integer, intent(in) :: nSites
    real(kind=prec), intent(in), optional       :: siteLocations(nSites,3)
    character(*), intent(in), optional          :: siteIDs(nSites)

    !  local variables
    integer      :: i,istat
    character(3) :: id

    if (.not. associated(rxDict)) then
      allocate(rxDict(nSites),STAT=istat)
      allocate(rx_map(nSites),stat=istat)
    end if

    do i = 1, nSites
        rx_map(i) % ptr => null()
    end do

    if (present(siteLocations)) then
      do i = 1,nSites
         rxDict(i)%x = siteLocations(i,:)
         if (present(siteIDs)) then
            rxDict(i)%id = siteIDs(i)
         else
            write(id,'(i6.6)') i
            rxDict(i)%id = id
         end if
      end do
    end if

    current_rx_count = 0

  end subroutine setup_rxDict

  function hash_fnv_1a(key) result(hash)
      
      use iso_fortran_env, only: int32, int64

      implicit none

      character(len=*), intent(in) :: key
      integer(int64), parameter :: FNV_OFFSET = 1469598103934665603_int64
      integer(int64), parameter :: FNV_PRIME  = 1099511628211_int64
      integer(int64) :: hash
      integer :: i, n

      hash = FNV_OFFSET
      n = len_trim(key)
      do i = 1, n
         hash = ieor(hash, int(iachar(key(i:i)), int64))
         hash = hash * FNV_PRIME
      end do

  end function hash_fnv_1a


  subroutine calculate_rx_index(code, total_rxs, idx)

      implicit none

      character(len=*), intent(in) :: code
      integer, intent(in) :: total_rxs
      integer, intent(out) :: idx
      integer :: hash_key

      hash_key = hash_fnv_1a(trim(code))
      idx = modulo(hash_key, total_rxs) + 1

  end subroutine calculate_rx_index

!**********************************************************************
! Updates the receiver dictionary with a new location and site ID.
! Returns the index of the new element.
! This is not efficient; but this would only be used a few times, with
! a small number of values, so convenience is much more of an issue here!
! NM: modified to include referance site info.

function update_rxDict(loc,id,Rx_azi,loc_ref,id_ref,code, total_rxs) result (iRx)

     character(*), intent(in)            :: id
     real(kind=prec), intent(in)         :: loc(3)
     real(kind=prec),intent(in),optional :: Rx_azi
     real(kind=prec),intent(in),optional :: loc_ref(3)
     character(*), intent(in),optional   :: id_ref
     character(*), intent(in),optional   :: code
     integer, intent(in), optional       :: total_rxs
     integer                             :: iRx
     ! local
     type(MTrx), target                  :: new
     type(MTrx), pointer, dimension(:)   :: temp
     integer                             :: nRx, istat,i
     logical							 :: new_Rx
     integer :: hash_idx

     ! Create a receiver for this location
     new%id = id
     new%x  = loc

     call calculate_rx_index(code, total_rxs, hash_idx)

     if (present(loc_ref)) then
     	new%r  = loc_ref
     	new%id_ref=id_ref
     end if

     if (present(Rx_azi)) then
     	new%Rx_azi  = Rx_azi
     end if	 

     ! If rxDict doesn't yet exist, create it
     if(.not. associated(rxDict)) then
        write(0,*) " we are allocating rxDict"
     	allocate(rxDict(1),STAT=istat)
     	rxDict(1) = new
     	iRx = 1
     	new_Rx = .true.
     	return
     end if

     !nRx = count_rx()

     if (associated(rx_map(hash_idx) % ptr)) then
        if (present(loc_ref)) then
           if (new%id_ref .eq. rx_map(hash_idx) % ptr % id_ref) then 
              rx_map(hash_idx) % ptr %r=loc_ref
              rx_map(hash_idx) % ptr %id_ref=id_ref
              iRx= rx_map(hash_idx) % ptr % iRx
              new_Rx = .false.
              return
           end if   
        elseif (present(Rx_azi)) then
           ! Check if the this site has same azimuth as what we have already in the Dic -
           ! this is only well-defined for CSEM; for MT, Rx_azi is not used since we're supporting
           ! the possibility that each field component has its own orientation (stored in DataSpace).
           if (new%Rx_azi .eq. rx_map(hash_idx) % ptr %Rx_azi) then 
              iRx= rx_map(hash_idx) % ptr % iRx
              new_Rx = .false.
              return
           end if              		   
        else    
         iRx= rx_map(hash_idx) % ptr % iRx
         new_Rx = .false.
         return
        end if 
     end if

     ! If the site really is new, append to the end of the dictionary
     new_Rx = .true.
     new%defined = .true.
     new % iRx = current_rx_count
     iRx = current_rx_count
     if (current_rx_count < size(rxDict)) then
         current_rx_count = current_rx_count + 1
         iRx = current_rx_count
         rxDict(iRx) = new
     else
         write(0,*) 'ERROR: The number of receivers (sites) in file has exceeded the hard-coded maximum', MAX_NRX
         write(0,*) 'ERROR: Please edit MAX_NRX in receivers.f90 and recompile'
         call ModEM_abort()
     end if

     !write(0,*) iRx, trim(code), hash_idx
     rx_map(hash_idx) % ptr => rxDict(iRx)

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

!**********************************************************************
! Writes the receiver dictionary to a file -- needed to associate rows in sensitivity
!  sensitivity matrix J with correct data vector elements.   This assumes that iounit
!  is opened for formated io-- file connection and  closing are done by calling routine

  subroutine write_rxDict_asc(iounit)

     ! local variables
     integer                     :: iounit, iRx, nRx

     if (.not. associated(rxDict)) then
        return
     end if

     !nRx = count_rx()
     write(iounit,*) nRx, '   Receivers'
     do iRx = 1, nRx
        write(iounit,'(i6,2x,a20,4f12.3,a20)') iRx,trim(rxDict(iRx)%id),rxDict(iRx)%x,&
          rxDict(iRx)%Rx_azi,trim(rxDict(iRx)%id_ref)
     enddo

  end subroutine write_rxDict_asc

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

     nRx = count_rx()
     header = 'Receiver Dictionary'
     write(iounit) header
     write(iounit) nRx
     do iRx = 1, nRx
        write(iounit) iRx,rxDict(iRx)%x,rxDict(iRx)%Rx_azi
        !   hard to read variable length records from a binary file,
        !   unless written as individual sequential records; hence this
        write(iounit) trim(rxDict(iRx)%id)
     enddo

  end subroutine write_rxDict_bin

  ! **************************************************************************
  ! Finds the end of the defined part of rxDict
  integer function count_rx()  result (nRx)

    !  local variables
    integer :: N, iRx

    if (.not. associated(rxDict)) then
      nRx = 0
      return
    end if

    N = size(rxDict)

    ! return when the first undefined value is encountered
    do iRx = 1,N
      if (.not. rxDict(iRx)%defined) then
         nRx = iRx-1
         return
      end if
    end do

    ! if they are all defined, nRx = size(rxDict)
    nRx = N

  end function count_rx


  ! **************************************************************************
  ! Trims the tail end of the pre-allocated rxDict
  subroutine trim_rxDict(nRx)

     integer, intent(in) :: nRx
     type(MTrx), pointer, dimension(:)   :: temp
     integer                             :: istat

    !nRx = count_rx()

    if (associated(rxDict)) then
       allocate(temp(nRx),STAT=istat)
       temp = rxDict(1:nRx)
       deallocate(rxDict,STAT=istat)
       allocate(rxDict(nRx),STAT=istat)
       rxDict = temp
       deallocate(temp,STAT=istat)
    end if

    write(0,*) 'rxDict: ', size(rxDict)

  end subroutine trim_rxDict


  ! **************************************************************************
  ! Cleans up and deletes receiver dictionary at end of program execution
  subroutine deall_rxDict()

	integer     :: istat

    if (associated(rxDict)) then
       deallocate(rxDict,STAT=istat)
    end if

  end subroutine deall_rxDict

end module receivers
