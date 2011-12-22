! *****************************************************************************
module transmitters
  ! This module contains the transmitter dictionary (txDict) for EARTH

  use math_constants
  use utilities
  use sg_sparse_vector
  use modelspace
  use iotypes

  implicit none

  public			:: initFreq, getFreq, deall_freqList

  ! ***************************************************************************
  ! * type transmitter_t contains the information about the frequencies; and also
  ! * the info on data functionals required to calculate for each frequency
  type :: transmitter_t

    ! the code may be needed to indicate frequency number in the full data set
    character(80)                           :: code
    ! frequency value
    real(8)                                 :: value
    ! period in days used for input/output only
    real(8)                                 :: period ! (in days)
    ! use / don't use secondary field formulation
    logical                                 :: secondaryField = .false.
    ! external source for this transmitter, stored as spherical harmonics
    integer                                 :: degree
    complex(8), pointer, dimension(:)       :: jExt
    ! internal source for this transmitter, stored as a sparse vector on the grid
    type(sparsevecc)                        :: jInt
    ! order index of this frequency that is used for the output
    integer                                 :: i
    ! nMode is number of "modes" for transmitter (e.g., 2 for MT)
    integer                                 :: nMode

  end type transmitter_t

  ! ***************************************************************************
  ! * storing the frequency dimensions and values in ascending order
  type :: Freq_List

    integer                                     :: n
    type (transmitter_t), pointer, dimension(:) :: info !nfreq

  end type Freq_List

  ! transmitter dictionary
  type (Freq_List), save, public                :: freqList

  ! read in the external source for all frequencies, and save it here
  type (modelParam_t), save, private            :: source,source_imag


Contains

  ! ***************************************************************************
  ! * initFreq reads the file fn_period (currently periods in days, could be
  ! * modified in the future) to store the freq or angular frequencies omega
  subroutine initFreq(cUserDef,myfreq)

    ! Not sure yet, whether freq, period or omega will be required
    implicit none
    type (userdef_control), intent(in)                           :: cUserDef
    type (Freq_List), intent(out)                           :: myfreq
    integer                                                 :: num
    integer                                                 :: i,j,ios=0,istat=0
    integer                                                 :: lmax,ncoeff,icoeff
    real(8)                                                 :: tmp
    real(8), dimension(:), allocatable                      :: value,days
    real(8), dimension(:), allocatable                      :: coeff_real,coeff_imag
    complex(8), dimension(:), allocatable                   :: coeff
    character(80)                                           :: basename,code
    complex(8)                                              :: p10(4)
    logical                                                 :: secondaryField

    open(ioTX,file=cUserDef%fn_period,status='old',form='formatted',iostat=ios)

    write(6,*) node_info,'Reading from the periods file ',trim(cUserDef%fn_period)
    read(ioTX,'(a)') label
    ! write(6,*) label

    read(ioTX,*) num
    allocate(value(num),days(num))
    do i=1,num
      read(ioTX,*) days(i) ! reading period in *days*
      value(i)=1/(days(i)*24*60*60) ! turn into freq.
      value(i)=clean(value(i))
    end do

!        read(ioTX,*) value(i) ! reading frequency value
!     days(i)=(1/value(i))/(24*60*60)   ! turn into period
!   end do

    close(ioTX)

    ! sort the values in ascending order
!    do i=1,num
!      do j=i+1,num
!        if(value(i)>value(j)) then
!          tmp=value(i)
!          value(i)=value(j)
!          value(j)=tmp
!          tmp=days(i)
!          days(i)=days(j)
!          days(j)=tmp
!        end if
!      end do
!    end do

    !omega(1:num)=2.0d0*pi*value(1:num)
    allocate(myfreq%info(num))
    myfreq%n=num
    myfreq%info(1:num)%value=value(1:num)
    myfreq%info(1:num)%period=days(1:num)
    do i=1,num
      myfreq%info(i)%i = i
    end do

    deallocate(value)

    basename=cUserDef%fn_period(1:index(cUserDef%fn_period,'.')-1)

    inquire(FILE=trim(basename)//'.codes',EXIST=exists)
    ! If file is not present, set the codes to the frequency numbers
    if (.not.exists) then
       do i=1,num
          write (myfreq%info(i)%code,'(i3.3)') i
       end do
    else
       open(ioTX,file=trim(basename)//'.codes',status='old',form='formatted',iostat=ios)
       read(ioTX,'(a)') label
       read(ioTX,*) num !should be the same as number of frequencies
       do i=1,num
          read(ioTX,*) myfreq%info(i)%code
       end do
       close(ioTX)
    end if

    ! Should we use secondary field formulation?
    !basename=cUserDef%fn_field
    !do i=1,num
    !    code = myfreq%info(i)%code
    !    inquire(FILE=trim(basename)//'_'//trim(code)//'.field',EXIST=exists)
    !    myfreq%info(i)%secondaryField = exists
    !end do

    ! define p10 just in case we need it
    p10(:) = C_ZERO
    p10(2) = C_ONE

    ! use / don't use secondary field formulation
    if ((index(cUserDef%secondary_field,'no') > 0) .or. (index(cUserDef%secondary_field,'0') > 0)) then
        secondaryField = .false.
    else
        secondaryField = .true.
    end if
    do i=1,num
        myfreq%info(i)%secondaryField = secondaryField
    end do

    inquire(FILE=trim(cUserDef%fn_extsource),EXIST=exists)
    if (exists) then
        ! source file contains the complex source for multiple periods (identical to T)
        call read_modelParam(source,cUserDef%fn_extsource,source_imag)
        allocate(coeff_real(source%nc),coeff_imag(source%nc),coeff(source%nc), STAT=istat)
        call getParamValues_modelParam(source,coeff_real)
        call getParamValues_modelParam(source_imag,coeff_imag)
        coeff = dcmplx(coeff_real,coeff_imag)
        lmax = getDegree_modelParam(source)
        ncoeff = (lmax + 1)**2 ! number of SH
        icoeff = 0
        do i=1,num
            myfreq%info(i)%degree = lmax
            allocate(myfreq%info(i)%jExt(ncoeff), STAT=istat)
            myfreq%info(i)%jExt = coeff(icoeff+1:icoeff+ncoeff)
            icoeff = icoeff + ncoeff
        end do
        deallocate(coeff_real,coeff_imag,coeff, STAT=istat)
    else
        ! can't use secondary field formulation: no sources specified, default to P10
        write(6,*) node_info,'Unable to use secondary field formulation: external sources not defined'
        do i=1,num
            myfreq%info(i)%degree = 1
            allocate(myfreq%info(i)%jExt(4), STAT=istat)
            myfreq%info(i)%jExt = p10
            myfreq%info(i)%secondaryField = .false. ! change this if you want to use p10 secondary fields...
        end do
    end if

    return

  end subroutine initFreq  !    initFreq

  ! **************************************************************************
  ! Cleans up and deletes transmitter dictionary at end of program execution

  subroutine deall_freqList()

    integer     :: i,istat

    if (associated(freqList%info)) then
       do i=1,freqList%n
           deallocate(freqList%info(i)%jExt, STAT=istat)
       end do
       deallocate(freqList%info,STAT=istat)
    end if
    call deall_modelParam(source)
    call deall_modelParam(source_imag)

  end subroutine deall_freqList

  ! ***************************************************************************
  ! * Uses frequency (Hz) to locate the transmitter in the list; default zero

  function getFreq(list,freq) result (ifreq)

    type(Freq_List), intent(in)      :: list
    real(8), intent(in)              :: freq
    integer                          :: ifreq
    ! local
    integer                          :: i
    real(8)                          :: const

    const = 1.0d-2

    ifreq = 0
    do i = 1,list%n
      if (abs((list%info(i)%value-freq)/list%info(i)%value) < const) then
        ifreq = i
        exit
      end if
    end do

  end function getFreq


end module transmitters
